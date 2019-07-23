/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.rappsilber.fdr;

import java.awt.GraphicsEnvironment;
import org.rappsilber.fdr.result.FDRResult;
import java.io.FileNotFoundException;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.text.DecimalFormat;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Timer;
import java.util.TimerTask;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import javax.swing.JOptionPane;
import javax.swing.SwingUtilities;
import org.rappsilber.fdr.entities.DBPSM;
import rappsilber.config.DBRunConfig;
import org.rappsilber.fdr.entities.PSM;
import org.rappsilber.fdr.entities.PeptidePair;
import org.rappsilber.fdr.entities.Protein;
import org.rappsilber.fdr.entities.ProteinGroupLink;
import org.rappsilber.fdr.entities.ProteinGroupPair;
import org.rappsilber.fdr.entities.ProteinGroup;
import rappsilber.ms.crosslinker.CrossLinker;
import rappsilber.ms.lookup.peptides.PeptideTree;
import rappsilber.ms.sequence.AminoAcid;
import rappsilber.ms.sequence.Peptide;
import rappsilber.ms.sequence.Sequence;
import rappsilber.ms.sequence.digest.Digestion;
import rappsilber.utils.CountOccurence;
import org.rappsilber.utils.IntArrayList;
import org.rappsilber.utils.RArrayUtils;
import org.rappsilber.utils.Version;
import rappsilber.config.DBConnectionConfig;
import rappsilber.gui.components.db.GetSearch;
import rappsilber.ms.sequence.AminoModification;

/**
 *
 * @author lfischer
 */
public class DBinFDR extends org.rappsilber.fdr.OfflineFDR implements XiInFDR {

    private Connection m_db_connection;
    private long       m_db_last_used;
    private Timer      m_db_autoclosetimer;
    private PreparedStatement updateValidateOverWrite;
    private PreparedStatement updateValidateNonOverWrite;
//    private PreparedStatement updateValidateOverWriteSetFDR;
//    private PreparedStatement updateValidateNonOverWriteSetFDR;
//    private PreparedStatement updateSetFDR;
    private PreparedStatement updateSetConfidence;
//    private PreparedStatement updateValidateOverWriteBatch;
//    private PreparedStatement updateValidateNonOverWriteBatch;
//    private PreparedStatement updateValidateOverWriteSetFDRBatch;
//    private PreparedStatement updateValidateNonOverWriteSetFDRBatch;
//    private PreparedStatement updateSetFDRBatch;
//    private PreparedStatement updateSetConfidenceBatch;
    private IntArrayList m_search_ids = new IntArrayList();
    private String m_connectionString = null;
    private String m_dbuser = null;
    private String m_dbpass = null;
    private HashMap<Long, Sequence> m_proteinSequences = new HashMap<Long, Sequence>();
    private Sequence m_noSequence = new Sequence(new AminoAcid[0]);
    private DBRunConfig m_conf;
    private HashMap<Integer,DBRunConfig> m_configs;
    private boolean m_writePrePostAA = false;
    public static DecimalFormat sixDigits = new DecimalFormat("###0.000000");

    int[] searchIDSetting;
    String filterSetting = "";
    private ArrayList<String> sequenceDBs;

    private GetSearch databaseProvider;
    private static boolean versionSet = false;
    private boolean flagAutoValidated = false;
    private boolean topOnly = true;
    /** should variable modifications (if detected) be marked out */
    private boolean markModifications;
    private boolean markSearchID;
    private boolean markRun;
    private HashSet<String> additionalInfoNames = new HashSet<>();
    
    private String command_line_auto_validate = null;

    static {
        setVersion();
    }

    public static void setVersion() {
        if (!versionSet) {
            versionSet = true;
            Version v = DBinFDR.getXiFDRVersion();
            String prevExt = v.extension;
            v.setExtension("$Rev: 69 $"); 
            //v.setExtension(v.extension +".dev"); 
            if (prevExt != null && !prevExt.isEmpty()) {
                v.setExtension(prevExt + "." + v.extension);
            }
        }
    }

    private Sequence loadSequence(long id) throws SQLException {
        Connection con = getDBConnection();
        Statement st = con.createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
        ResultSet rs = st.executeQuery("SELECT "
                + " sequence, "
                + " CASE "
                + "    WHEN header IS NULL AND accession_number IS NULL THEN name "
                + "    WHEN header IS NULL THEN accession_number ELSE header "
                + " END as header, "
                + " accession_number "
                + " FROM protein where ID = " + id);
        if (rs.next()) {
            String seqStr = rs.getString(1);
            if (seqStr.matches("'.*'")) {
                seqStr = seqStr.substring(1, seqStr.length() - 2);
            }

            Sequence seq = new Sequence(seqStr, getConfig());
            String headerStr = rs.getString(2);
            if (headerStr.matches("'.*'")) {
                headerStr = headerStr.substring(1, headerStr.length() - 2);
            }

            seq.setFastaHeader(rs.getString(2));

            m_proteinSequences.put(id, seq);
            rs.close();
            return seq;
        } else {
            m_proteinSequences.put(id, m_noSequence);
            rs.close();
            return m_noSequence;
        }
    }

    private Sequence getSequence(Protein p) {
        try {
            long id = p.getId();
            Sequence seq = m_proteinSequences.get(id);

            if (seq == null) {
                return loadSequence(id);
            }
            return seq;
        } catch (SQLException ex) {
            Logger.getLogger(DBinFDR.class.getName()).log(Level.SEVERE, null, ex);
        }
        return m_noSequence;
    }

    public String getLinkWindow(PeptidePair pp, int proteingroup, int window) {
//        return getLinkWindow(pp.getLink(), proteingroup, window);
        try {
            int pepLink = 0;
            org.rappsilber.fdr.entities.Peptide pep;
            if (proteingroup == 0) {
                pep = pp.getPeptide1();
                pepLink = pp.getPeptideLinkSite1() - 1;
            } else {
                pep = pp.getPeptide2();
                pepLink = pp.getPeptideLinkSite2() - 1;
            }
            ProteinGroup pg = pep.getProteinGroup();
            HashMap<Protein, IntArrayList> positions = pep.getPositions();

            Sequence pepSeq;
            if (pep.getSequence().matches("X[0-9]+(\\.[0-9]+)?")) {
                pepSeq = new Sequence("X", getConfig());
            } else {
                pepSeq = new Sequence(pep.getSequence(), getConfig());
            }

            int pepLen = pepSeq.length();

            // prepare the return window 
            ArrayList<CountOccurence<AminoAcid>> sequenceWindow = new ArrayList<CountOccurence<AminoAcid>>(2 * window + 1);

            int from = window - pepLink;
            int to = window + (pepLen - pepLink - 1);

            int pepPos = 0;
            for (int i = 0; i <= 2 * window; i++) {
                CountOccurence<AminoAcid> o = new CountOccurence<AminoAcid>();
                sequenceWindow.add(o);
                if (i >= from && i <= to) {
                    o.add(pepSeq.aminoAcidAt(pepPos++));
                } else {
                    o.add(null);
                }
            }

            // add everything from the surounding proteins into the window
            // load in all sequence windows
            int pid = 0;
            for (Protein p : pg) {
                Sequence seq = getSequence(p);

                IntArrayList sites = positions.get(p);
                for (int protPepFrom : sites) {
                    protPepFrom--;
                    int protPepTo = protPepFrom + pepLen;

                    // protein link
                    int lp = pepLink + protPepFrom;

                    int protFrom = lp - window;
                    int protTo = lp + window;

                    int f = Math.max(protFrom, 0);
                    int swp = f - protFrom;

                    for (int aap = f; aap < protPepFrom; aap++) {
                        sequenceWindow.get(swp).add(seq.aminoAcidAt(aap));
                        swp++;
                    }
                    swp = protPepTo - protFrom;
                    int t = Math.min(seq.length() - 1, protTo);
                    for (int aap = protPepTo; aap <= t; aap++) {
                        sequenceWindow.get(swp).add(seq.aminoAcidAt(aap));
                        swp++;
                    }
                    // next sequence
                    pid++;
                }
            }

            StringBuffer ret = new StringBuffer(2 * window + 1);
            for (CountOccurence<AminoAcid> o : sequenceWindow) {
                ArrayList<AminoAcid> so = o.getSortedList(true);
                AminoAcid aa = so.get(0);
                if (aa == null) {
                    if (so.size() == 1) {
                        ret.append(".");
                    } else {
                        ret.append(so.get(1).SequenceID);
                    }
                } else {
                    ret.append(aa.SequenceID);
                }
            }

            return ret.toString();
        } catch (Exception e) {
            return e.getMessage();
        }
    }

    public String getLinkWindow(ProteinGroupLink l, int proteingroup, int window) {
        try {
            ProteinGroup pg = null;
            HashMap<Protein, IntArrayList> positions = null;
            if (proteingroup == 0) {
                pg = l.getProteinGroup1();
                positions = l.getPosition1();
            } else {
                pg = l.getProteinGroup2();
                positions = l.getPosition2();
            }

            ArrayList<CountOccurence<AminoAcid>> sequenceWindow = new ArrayList<CountOccurence<AminoAcid>>(2 * window + 1);
            for (int i = 2 * window; i >= 0; i--) {
                sequenceWindow.add(new CountOccurence<AminoAcid>());
            }

            // load in all sequence windows
            int pid = 0;
            for (Protein p : pg) {
                Sequence seq = getSequence(p);

                IntArrayList sites = positions.get(p);
                for (int lp : sites) {
                    // make it start at 0
                    lp--;

                    int swp = 0;
                    for (int aap = lp - window; aap <= lp + window; aap++) {
                        if (aap < 0 || aap >= seq.length()) {
                            sequenceWindow.get(swp).add(null);
                        } else {
                            sequenceWindow.get(swp).add(seq.baseAminoAcidAt(aap));
                        }
                        swp++;
                    }
                    // next sequence
                    pid++;
                }
            }

            StringBuffer ret = new StringBuffer(2 * window + 1);
            for (CountOccurence<AminoAcid> o : sequenceWindow) {
                ArrayList<AminoAcid> so = o.getSortedList(true);
                AminoAcid aa = so.get(0);
                if (aa == null) {
                    if (so.size() == 1) {
                        ret.append(".");
                    } else {
                        ret.append(so.get(1).SequenceID);
                    }
                } else {
                    ret.append(aa.SequenceID);
                }
            }

            return ret.toString();
        } catch (Exception e) {
            return e.getMessage();
        }

    }

    /**
     * @return the m_db_connection
     */
    public Connection getDBConnection() {
        m_db_last_used = Calendar.getInstance().getTimeInMillis();
        if (m_db_autoclosetimer == null) {
            m_db_autoclosetimer = new Timer("Timer - auto close db", true);
            m_db_autoclosetimer.schedule(new TimerTask() {
                boolean running = false;
                @Override
                public void run() {
                    if (!running) {
                        running = true;
                        Calendar n = Calendar.getInstance();
                        // if the database connection was not requested for 30 minutes close it.
                        if (n.getTimeInMillis() - m_db_last_used > 1800000) {
                            closeConnection();
                        }
                        running = false;
                    }
                }
            }, 60000, 60000);
        }
        try {
            if (!ensureConnection()) {
                Logger.getLogger(this.getClass().getName()).log(Level.INFO, "The database connection seem to be closed and can't be opend now");
                return null;
            }
        } catch (SQLException ex) {
            return null;
        }
        return m_db_connection;
    }

    public void setDBConnection(String connection, String user, String password) throws SQLException {
        m_connectionString = connection;
        m_dbuser = user;
        m_dbpass = password;
        this.m_db_connection = DriverManager.getConnection(connection, user, password);
        setupPreparedStatements();
    }

    /**
     * @param m_db_connection the m_db_connection to set
     */
    public void setDBConnection(Connection connection) throws SQLException {
        this.m_db_connection = connection;
        setupPreparedStatements();
    }

    public void setDatabaseProvider(GetSearch getSearch) throws SQLException {
        databaseProvider = getSearch;
        setDBConnection(databaseProvider.getConnection());
    }

    public synchronized void closeConnection()  {
        try {
            if (m_db_connection!=null) {
                m_db_connection.close();
                m_db_connection=null;
            }
        } catch (SQLException ex) {
            Logger.getLogger(DBinFDR.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    public synchronized boolean ensureConnection() throws SQLException {
        boolean isconnected = false;
//        try {
//            isconnected = m_db_connection.isValid(30);
//            
//        } catch (SQLException ex) {
//            Logger.getLogger(DBinFDR.class.getName()).log(Level.INFO, null, ex);
//        } catch (Exception e) {
        try {
            Statement st = m_db_connection.createStatement();
            ResultSet rs = st.executeQuery("select 1+1");
            rs.close();
            st.close();
            isconnected = true;
        } catch (Exception sex) {
            Logger.getLogger(DBinFDR.class.getName()).log(Level.INFO, "Database connection is closed/ non-functioning. Will try to reopen", sex);
            try {
                m_db_connection.close();
            } catch (Exception e) {
            }

        }
//        }
//        }

        if (!isconnected) {
            Logger.getLogger(this.getClass().getName()).log(Level.WARNING, "Lost connection to database. Will try to reopen");

            try {
                if (m_connectionString != null && m_dbuser != null && m_dbpass != null) {
                    this.m_db_connection = DriverManager.getConnection(m_connectionString, m_dbuser, m_dbpass);
                } else if (databaseProvider != null) {
                    setDBConnection(databaseProvider.getConnection());
                } else {
                    Logger.getLogger(DBinFDR.class.getName()).log(Level.SEVERE, "Currently I don't have a way to reopen the connection");
                    return false;
                }

                isconnected = true;

                setupPreparedStatements();

                if (!GraphicsEnvironment.isHeadless()) {
                    SwingUtilities.invokeLater(new Runnable() {

                        public void run() {
                            JOptionPane.showMessageDialog(null, "Reopened the database connection");
                        }
                    });
                }

            } catch (SQLException ex) {
                Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, "Could not reopen the connection to the database", ex);
                throw ex;
            }
        }
        return isconnected;

    }

    protected void setupPreparedStatements() throws SQLException {
        updateValidateOverWrite = this.m_db_connection.prepareStatement("UPDATE spectrum_match set validated = ? WHERE id = ?");
        updateValidateNonOverWrite = this.m_db_connection.prepareStatement("UPDATE spectrum_match set validated = ? WHERE id = ? AND validated is null");
//        updateValidateNonOverWriteSetFDR = this.m_db_connection.prepareStatement("UPDATE spectrum_match SET validated = CASE WHEN validated is null THEN ? ELSE validated END, set peptide_fdr = ?, set link_fdr = ?, ppi_fdr = ? SET WHERE id = ?");
//        updateValidateOverWriteSetFDR = this.m_db_connection.prepareStatement("UPDATE spectrum_match SET validated = ?, set peptide_fdr = ?, set link_fdr = ?, ppi_fdr = ? SET WHERE id = ?");
//        updateSetFDR = this.m_db_connection.prepareStatement("UPDATE spectrum_match SET peptide_fdr = ?, set link_fdr = ?, ppi_fdr = ? SET WHERE id = ?");
        updateSetConfidence = this.m_db_connection.prepareStatement("UPDATE spectrum_match rescored = ? SET WHERE id = ?");

//        updateValidateOverWriteBatch = this.m_db_connection.prepareStatement("UPDATE spectrum_match set validated = ? WHERE id = any(?)");
//        updateValidateNonOverWriteBatch = this.m_db_connection.prepareStatement("UPDATE spectrum_match set validated = ? WHERE id = any(?) AND validated is null");
//        updateValidateNonOverWriteSetFDRBatch = this.m_db_connection.prepareStatement("UPDATE spectrum_match SET validated = CASE WHEN validated is null THEN ? ELSE validated END, set peptide_fdr = ?, set link_fdr = ?, ppi_fdr = ? SET WHERE id = any(?)");
//        updateValidateOverWriteSetFDRBatch = this.m_db_connection.prepareStatement("UPDATE spectrum_match SET validated = ?, set peptide_fdr = ?, set link_fdr = ?, ppi_fdr = ? SET WHERE id = any(?)");
//        updateSetFDRBatch = this.m_db_connection.prepareStatement("UPDATE spectrum_match SET peptide_fdr = ?, set link_fdr = ?, ppi_fdr = ? SET WHERE id = any(?)");
//        updateSetConfidenceBatch = this.m_db_connection.prepareStatement("UPDATE spectrum_match rescored = ? SET WHERE id = any(?)");
    }

    @Override
    public String getSource() {
        if (m_connectionString != null && !m_connectionString.isEmpty()) {
            return "\n\tDB:" + m_connectionString + "\n\tIDs:" + RArrayUtils.toString(m_search_ids, ";") + (filterSetting == null || filterSetting.isEmpty() ? "" : "\n\tFilter: " + filterSetting);
        } else {
            return "\n\tIDs:" + RArrayUtils.toString(m_search_ids, ",") + (filterSetting == null || filterSetting.isEmpty() ? "" : "\n\tFilter: " + filterSetting);
        }
    }

    public IntArrayList getSearchIDs() {
        return m_search_ids;
    }

    /**
     * @return the m_conf
     */
    @Override
    public DBRunConfig getConfig() {
        return m_conf;
    }

    @Override
    public DBRunConfig getConfig(int searchid) {
        if (m_search_ids.size() == 1 && searchid == m_search_ids.get(0)) {
            return getConfig();
        }
        if (m_configs == null) {
            m_configs = new HashMap<>();
        }
        DBRunConfig ret = m_configs.get(searchid);
        if (ret == null) {
            try {
                ret = new DBRunConfig(getDBConnection());
                ret.readConfig(searchid);
                m_configs.put(searchid, ret);
            } catch (SQLException ex) {
                Logger.getLogger(DBinFDR.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        return ret;
    }

    /**
     * @param m_conf the m_conf to set
     */
    public void setConfig(DBRunConfig m_conf) {
        this.m_conf = m_conf;
    }

    /**
     * @return the writePrePostAA
     */
    public boolean isWritePrePostAA() {
        return m_writePrePostAA;
    }

    /**
     * @param writePrePostAA the writePrePostAA to set
     */
    public void setWritePrePostAA(boolean writePrePostAA) {
        this.m_writePrePostAA = writePrePostAA;
    }

    /**
     * @return the sequenceDBs
     */
    public ArrayList<String> getSequenceDBs() {
        return sequenceDBs;
    }

    public void setMarkSearchID(boolean b) {
        this.markSearchID = b;
    }

    public void setMarkRun(boolean b) {
        this.markRun = b;
    }

    public boolean markSearchID() {
        return this.markSearchID;
    }

    public boolean markRun() {
        return this.markRun;
    }

    private class hashableXiPeptide {

        public rappsilber.ms.sequence.Peptide peptide;
        int hashcode;

        public hashableXiPeptide(Peptide peptide) {
            this.peptide = peptide;
            hashcode = peptide.toString().hashCode();
        }

        public int hashCode() {
            return hashcode;
        }

        public boolean equals(Object o) {
            return ((hashableXiPeptide) o).peptide.equalSequence(peptide);
        }

    }

    public DBinFDR() {
    }

    public DBinFDR(Connection connection) {
        this.m_db_connection = connection;
    }

    public DBinFDR(Connection connection, int[] peptideLengthGroups) {
        super(peptideLengthGroups);
        this.m_db_connection = connection;
    }

    public void readDB() throws SQLException {
        readDB(searchIDSetting, filterSetting, true);
    }

    public void readDB(int searchId, String filter, boolean topOnly) throws SQLException {
        readDB(new int[]{searchId}, filter, topOnly);
    }

    public void readDB(int[] searchIds, String filter, boolean topOnly) throws SQLException {
        this.readDB(searchIds, filter, topOnly, new ArrayList<Long>(), 0);
    }

    public void readDB(int[] searchIds, String filter, boolean topOnly, ArrayList<Long> skip, int tries) throws SQLException {

        class psminfo {

            ArrayList<Integer> protein1IDs = new ArrayList<>();
            ArrayList<Integer> protein2IDs = new ArrayList<>();
            ArrayList<Integer> protein1Positions = new ArrayList<>();
            ArrayList<Integer> protein2Positions = new ArrayList<>();
        };

        if (!ensureConnection()) {
            return;
        }

        setConfig(new DBRunConfig(getDBConnection()));
        getConfig().readConfig(searchIds);
        boolean isTargted = false;

        for (CrossLinker xl : getConfig().getCrossLinker()) {
            if (xl.getName().toLowerCase().contains("targetmodification")) {
                isTargted = true;
            }
        }

        String dbNameQuerry = "Select id,name from search_sequencedb ss inner join sequence_file sf ON ss.search_id in (" + RArrayUtils.toString(searchIds, ",") + ") and ss.seqdb_id = sf.id";

        Statement dbnSt = getDBConnection().createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
        ResultSet rsDBn = dbnSt.executeQuery(dbNameQuerry);
        StringBuilder sb = new StringBuilder();
        sequenceDBs = new ArrayList<String>(1);
        HashSet<Integer> dbIds = new HashSet<Integer>();

        while (rsDBn.next()) {
            int id = rsDBn.getInt(1);
            if (!dbIds.contains(id)) {
                sequenceDBs.add(rsDBn.getString(2));
                dbIds.add(id);
            }
        }

        PeptidePair.ISTARGETED = isTargted;
        boolean xi3db = false;
        try {
            Statement stm = getDBConnection().createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
            stm.setFetchSize(100);
            ResultSet rs = stm.executeQuery("SELECT * FROM spectrum_source limit 0;");
            rs.close();
            xi3db = true;
        } catch (SQLException ex) {
            xi3db = false;
        }

        if (!m_search_ids.isEmpty()) {
            if ((filter == null && !filterSetting.isEmpty()) || (!filter.contentEquals(filterSetting))) {
                filterSetting = filterSetting + "\n Search IDS:" + RArrayUtils.toString(searchIds, ",") +":" + filter;
            }
        } else if (filter != null && !filter.isEmpty()) {
            filterSetting = filter;
        }

        for (int currentsearch = 0 ; currentsearch<searchIds.length;currentsearch++){
            int searchId = searchIds[currentsearch];
            // wee read that one already
            if (m_search_ids.contains(searchId)) {
                continue;
            }

            ArrayList<String> scorenames = new ArrayList<>();
            if (xi3db) {
                // retrive the subscore names
                String scoreNameQuerry = 
                        "SELECT scorenames FROM search WHERE id = " + searchId +";";
                Statement stm = getDBConnection().createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
                stm.setFetchSize(100);
                ResultSet rs = stm.executeQuery(scoreNameQuerry);
                if (rs.next()) {
                    java.sql.Array sn = rs.getArray(1);
                    if (sn != null) {
                        for (String s : (String[])sn.getArray()) {
                            additionalInfoNames.add(s);
                            scorenames.add(s);
                        }
                    }
                }
            }
            


            String matchQuerry;
            matchQuerry
                    = "SELECT * FROM (SELECT sm.id AS psmID, \n"
                    + "p1.sequence AS pepSeq1, \n"
                    + "p2.sequence AS pepSeq2, \n"
                    + "p1.peptide_length as peplen1, \n"
                    + "p2.peptide_length as peplen2, \n"
                    + "mp1.link_position + 1 as site1, \n"
                    + "mp2.link_position + 1 as site2, \n"
                    + "pr1.is_decoy AS isDecoy1, \n"
                    + "pr2.is_decoy AS isDecoy2, \n"
                    + "sm.precursor_charge AS calc_charge, \n"
                    + "sm.score, \n"
                    + "pr1.accession_number AS accession1, \n"
                    + "CASE WHEN pr1.description is not null THEN pr1.description WHEN pr1.name is not null THEN pr1.name ELSE pr1.accession_number END AS  description1, \n"
                    + "pr2.accession_number AS accession2, \n"
                    + "CASE WHEN pr2.description is not null THEN pr2.description WHEN pr2.name is not null THEN pr2.name ELSE pr2.accession_number END AS  description2, \n"
                    + "hp1.peptide_position + 1 AS pepPosition1,  \n"
                    + "hp2.peptide_position + 1 AS pepPosition2, \n"
                    + "CASE WHEN p2.sequence IS NULL THEN 1 ELSE (4.0/5.0+(p1.peptide_length/(p1.peptide_length+p2.peptide_length)))/2 END AS score_ratio \n"
                    + " , s.precursor_charge as exp_charge \n"
                    + " , hp1.display_site as display_site1\n"
                    + " , hp2.display_site as display_site2\n"
                    + " , pr1.id as protein1id\n"
                    + " , pr2.id as protein2id\n"
                    + " , p1.id as peptide1id\n"
                    + " , p2.id as peptide2id\n"
                    + " , v.run_name, v.scan_number \n"
                    + " , pr1.sequence AS protein1sequence \n"
                    + " , pr2.sequence AS protein2sequence\n"
                    + (xi3db  || (filter != null && filter.toLowerCase().contains("scorep1coverage"))
                        ? (xi3db ? ", scorepeptide1matchedconservative" : " , coverage1 ") + " AS scoreP1Coverage \n"
                        : " , 0 AS scoreP1Coverage \n")
                    + (xi3db  || (filter != null && filter.toLowerCase().contains("scorep2coverage"))
                        ? (xi3db ? ", scorepeptide2matchedconservative" : " , coverage2 ") + " AS scoreP2Coverage \n"
                        : " , 0 AS scoreP2Coverage \n")
                    + (xi3db  || (filter != null && filter.toLowerCase().replace("linksitedelta", "").contains("delta"))
                        ? (xi3db ? ", scoredelta " : " , delta ") + " AS deltaScore\n"
                        : " , 0  AS deltaScore \n")
                    + (xi3db  || (filter != null && filter.toLowerCase().contains("linksitedelta"))
                        ? (xi3db ? ", scorelinksitedelta " : " , LinkSiteDelta ") + " AS LinkSiteDelta\n"
                    : " , 0  AS LinkSiteDelta \n")
                    + " , s.precursor_mz AS exp_mz\n"
                    + " , sm.calc_mass\n"
                    + " , sm.precursor_charge as match_charge\n"
                    + " , p1.mass as pep1mass\n"
                    + " , p2.mass as pep2mass\n"
                    + " , sm.search_id\n"
                    + " , sm.spectrum_id\n"
                    + " , pr1.name as protein1name\n"
                    + " , pr2.name as protein2name\n"
                    + " , s.precursor_charge\n"
                    + " , autovalidated\n"
                    + " " + (xi3db ? ", cl.name " : " , v.crosslinker ") + " AS crosslinker \n"
                    + " " + (xi3db ? ", cl.mass " : " , null ") + " AS clmass \n"
                    + " , sm.rank \n"
                    + " , s.scan_index\n"
                    + (xi3db ? " , plf.name as peaklistfile\n":"")
                    + " , " + (xi3db ? " scoremgxrank AS mgxrank\n" : " 0 AS mgxrank\n")
                    + " , " + (xi3db ? " scorecleavclpep1fragmatched AS cleavclpep1fragmatched\n" : " 0 AS cleavclpep1fragmatched\n")
                    + " , " + (xi3db ? " scorecleavclpep2fragmatched AS cleavclpep2fragmatched\n" : " 0 AS cleavclpep2fragmatched\n")
                    + " , " + (xi3db ? " scores AS subscores\n" : " null AS subscores\n")
                    + " , s.precursor_intensity\n"
                    + " , s.elution_time_start as retentiontime\n"
                    + " \n"
                    + "FROM \n"
                    + "  (SELECT * FROM Spectrum_match WHERE Search_id = " + searchId + (topOnly ? " AND dynamic_rank = 't'":"")+ " AND score>0) sm \n"
                    + "  INNER JOIN (\n"
                    + (xi3db ? "SELECT ss.name as run_name, s.scan_number, sm.id as spectrum_match_id FROM (select * from spectrum_match where Search_id = " + searchId + (topOnly ? " AND dynamic_rank = 't'":"")+ " AND score>0) sm inner join spectrum s on sm.spectrum_id = s.id INNER JOIN spectrum_source ss on s.source_id = ss.id\n"
                            : "SELECT run_name, scan_number, spectrum_match_id, crosslinker FROM v_export_materialized WHERE Search_id = " + searchId +  (topOnly ? " AND dynamic_rank = 't'":"")+ " \n")
                    + ") v \n"
                    + "    ON v.spectrum_match_id = sm.id\n"
                    + "    inner join \n"
                    //                    + "   (SELECT * from matched_peptide where search_id  = "+ searchId +" AND match_type =1) mp1 on sm.id = mp1.match_id   LEFT OUTER JOIN \n"
                    //                    + "   (SELECT * from matched_peptide where search_id  = "+ searchId +" AND match_type =2) mp2 on sm.id = mp2.match_id   INNER JOIN \n"
                    + "   matched_peptide mp1 on sm.id = mp1.match_id and mp1.match_type =1  LEFT OUTER JOIN \n"
                    + "   matched_peptide mp2 on sm.id = mp2.match_id AND mp2.match_type = 2  INNER JOIN \n"
                    + "   peptide p1 on mp1.peptide_id = p1.id LEFT OUTER JOIN  \n"
                    + "   peptide p2 on mp2.peptide_id = p2.id INNER JOIN \n"
                    + "   has_protein hp1 ON mp1.peptide_id = hp1.peptide_id LEFT OUTER JOIN  \n"
                    + "   has_protein hp2 ON mp2.peptide_id = hp2.peptide_id INNER JOIN \n"
                    + "   protein pr1 ON hp1.protein_id = pr1.id LEFT OUTER JOIN \n"
                    + "   protein pr2 ON hp2.protein_id = pr2.id\n"
                    + "   INNER JOIN  \n"
                    + "   spectrum s ON sm.spectrum_id = s.id \n"
                    + ((!xi3db) && filter != null && filter.toLowerCase().contains("scorep1coverage")
                    ? "   LEFT OUTER JOIN (SELECT spectrum_match_id, score AS coverage1 FROM spectrum_match_score WHERE  score_id = (select id from score where name ='peptide1 unique matched conservative')) p1c on sm.id = p1c.spectrum_match_id \n"
                    : "")
                    + ((!xi3db) && filter != null && filter.toLowerCase().contains("scorep2coverage")
                    ? "   LEFT OUTER JOIN (SELECT spectrum_match_id, score AS coverage2 FROM spectrum_match_score WHERE  score_id = (select id from score where name ='peptide2 unique matched conservative')) p2c on sm.id = p2c.spectrum_match_id \n"
                    : "")
                    + ((!xi3db) && filter != null && filter.toLowerCase().replace("linksitedelta", "").contains("delta")
                    ? "   LEFT OUTER JOIN (SELECT spectrum_match_id, score AS delta FROM spectrum_match_score WHERE  score_id = (select id from score where name ='delta')) deltascore on sm.id = deltascore.spectrum_match_id \n"
                    : "")
                    + ((!xi3db) && filter != null && filter.toLowerCase().contains("linksitedelta")
                    ? "   LEFT OUTER JOIN (SELECT spectrum_match_id, score AS LinkSiteDelta FROM spectrum_match_score WHERE  score_id = (select id from score where name ='LinkSiteDelta')) LinkSiteDelta on sm.id = LinkSiteDelta.spectrum_match_id \n"
                    : "")
                    + (xi3db ? " LEFT OUTER JOIN peaklistfile plf on s.peaklist_id = plf.id\n":"")
                    + (xi3db ? " LEFT OUTER JOIN crosslinker cl on mp1.crosslinker_id = cl.id \n" : "")
                    //                    + "   LEFT OUTER JOIN spectrum_match_score p2c on sm.id = p2c.spectrum_match_id AND p2c.score_id = (select id from score where name ='peptide1 unique matched conservative') \n"
                    //                    + "   LEFT OUTER JOIN spectrum_match_score sd on sm.id = sd.spectrum_match_id AND sd.score_id = (select id from score where name ='delta') \n"
                    + " WHERE  \n"
                    //  + " (sd.score>0) AND \n"
                    + " (s.precursor_charge = sm.precursor_charge OR sm.precursor_charge < 6) \n"
                    + " ORDER BY sm.score DESC)  i \n"
                    + (filter == null || filter.isEmpty() ? "" : " WHERE ( " + filter + " )");
            //                + " --INNER JOIN  \n"
            //                + " --  spectrum_match_score sms1 on sm.id = sms1.spectrum_match_id AND sms1.score_id = 5 LEFT OUTER JOIN \n"
            //                + " --  spectrum_match_score sms2 on sm.id = sms2.spectrum_match_id AND sms2.score_id = 6 ;"; 

            Logger.getLogger(this.getClass().getName()).log(Level.INFO, "read from db");
            Logger.getLogger(this.getClass().getName()).log(Level.INFO, matchQuerry);
            if (searchIds.length >1)
                Logger.getLogger(this.getClass().getName()).log(Level.INFO, "Search {0} of {1}", new Object[]{currentsearch+1, searchIds.length});
            getDBConnection().setAutoCommit(false);
            Statement stm = getDBConnection().createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
            stm.setFetchSize(100);
            ResultSet rs = stm.executeQuery(matchQuerry);
            Logger.getLogger(this.getClass().getName()).log(Level.INFO, "go through the results");

            int psmIDColumn = rs.findColumn("psmID");
            int pepSeq1Column = rs.findColumn("pepSeq1");
            int pepSeq2Column = rs.findColumn("pepSeq2");
            int peplen1Column = rs.findColumn("peplen1");
            int peplen2Column = rs.findColumn("peplen2");
            int site1Column = rs.findColumn("site1");
            int site2Column = rs.findColumn("site2");
            int isDecoy1Column = rs.findColumn("isDecoy1");
            int isDecoy2Column = rs.findColumn("isDecoy2");
            int calc_chargeColumn = rs.findColumn("calc_charge");
            int scoreColumn = rs.findColumn("score");
            int accession1Column = rs.findColumn("accession1");
            int description1Column = rs.findColumn("description1");
            int accession2Column = rs.findColumn("accession2");
            int description2Column = rs.findColumn("description2");
            int pepPosition1Column = rs.findColumn("pepPosition1");
            int pepPosition2Column = rs.findColumn("pepPosition2");
            int score_ratioColumn = rs.findColumn("score_ratio");
            int exp_chargeColumn = rs.findColumn("exp_charge");
            int display_site1Column = rs.findColumn("display_site1");
            int display_site2Column = rs.findColumn("display_site2");
            int protein1idColumn = rs.findColumn("protein1id");
            int protein2idColumn = rs.findColumn("protein2id");
            int peptide1idColumn = rs.findColumn("peptide1id");
            int peptide2idColumn = rs.findColumn("peptide2id");
            int run_nameColumn = rs.findColumn("run_name");
            int scan_numberColumn = rs.findColumn("scan_number");
            int protein1sequenceColumn = rs.findColumn("protein1sequence");
            int protein2sequenceColumn = rs.findColumn("protein2sequence");
            int scoreP1CoverageColumn = rs.findColumn("scoreP1Coverage");
            int scoreP2CoverageColumn = rs.findColumn("scoreP2Coverage");
            int deltaScoreColumn = rs.findColumn("deltaScore");
            int LinkSiteDeltaColumn = rs.findColumn("LinkSiteDelta");
            int exp_mzColumn = rs.findColumn("exp_mz");
            int calc_massColumn = rs.findColumn("calc_mass");
            int pep1massColumn = rs.findColumn("pep1mass");
            int pep2massColumn = rs.findColumn("pep2mass");
            int search_idColumn = rs.findColumn("search_id");
            int spectrum_idColumn = rs.findColumn("spectrum_id");
            int xlColumn = rs.findColumn("crosslinker");
            int avColumn = rs.findColumn("autovalidated");
            int rankColumn = rs.findColumn("rank");
            int scanIndexColumn = rs.findColumn("scan_index");
            int peaklistfileColumn = rs.findColumn("peaklistfile");
            int mgxrankColumn = rs.findColumn("mgxrank");
            int cleavclpep1fragmatchedColumn = rs.findColumn("cleavclpep1fragmatched");
            int cleavclpep2fragmatchedColumn = rs.findColumn("cleavclpep2fragmatched");
            int subscoresColumn = rs.findColumn("subscores");
            int precIntensityColumn = rs.findColumn("precursor_intensity");
            int retentiontimeColumn = rs.findColumn("retentiontime");
            int clMassColumn = rs.findColumn("clmass");
            Pattern modDetect = Pattern.compile(".*[^A-Z].*");
            HashSet<Double> tmmodcount = new HashSet<>();

            int total = 0;
            try {
                while (rs.next()) {
                    total++;

                    long ipsmID = rs.getLong(psmIDColumn);
                    if (skip.contains(ipsmID)) {
                        continue;
                    }

                    String psmID = Long.toString(ipsmID);
                    String pepSeq1 = rs.getString(pepSeq1Column);
                    String pepSeq2 = rs.getString(pepSeq2Column);
                    int peplen1 = rs.getInt(peplen1Column);
                    int peplen2 = rs.getInt(peplen2Column);
                    int site1 = rs.getInt(site1Column);
                    int site2 = pepSeq2 == null ? -1 : rs.getInt(site2Column);
                    boolean isDecoy1 = rs.getBoolean(isDecoy1Column);
                    boolean isDecoy2 = rs.getBoolean(isDecoy2Column);
                    int charge = rs.getInt(calc_chargeColumn);
                    double score = rs.getDouble(scoreColumn);
                    String accession1 = rs.getString(accession1Column);
                    String description1 = rs.getString(description1Column);
                    String accession2 = rs.getString(accession2Column);
                    String description2 = rs.getString(description2Column);
                    int pepPosition1 = rs.getInt(pepPosition1Column);
                    int pepPosition2 = pepSeq2 == null ? -1 : rs.getInt(pepPosition2Column);
                    //            double coverage1  = rs.getDouble(18);
                    //            double coverage2  = rs.getDouble(19);
                    //            double scoreRatio = coverage1/(coverage1+coverage2);
                    double scoreRatio = rs.getDouble(score_ratioColumn);
                    int spectrum_charge = rs.getInt(exp_chargeColumn);
                    long protein1ID = rs.getLong(protein1idColumn);
                    long protein2ID = rs.getLong(protein2idColumn);
                    long pep1ID = rs.getLong(peptide1idColumn);
                    long pep2ID = rs.getLong(peptide2idColumn);
                    String run = rs.getString(run_nameColumn);
                    String scan = rs.getString(scan_numberColumn);
                    String sequence1 = rs.getString(protein1sequenceColumn);
                    String sequence2 = rs.getString(protein2sequenceColumn);
                    boolean autovalidated = rs.getBoolean(avColumn);
                    int rank = rs.getInt(rankColumn);

                    if (accession1 == null) {
                        accession1 = description1;
                    }

                    if (accession2 == null) {
                        accession2 = description2;
                    }

                    if (sequence1 != null) {
                        sequence1 = sequence1.replaceAll("[^A-Z]", "");
                    }

                    if (sequence2 != null) {
                        sequence2 = sequence2.replaceAll("[^A-Z]", "");
                    }

                    double p1c = rs.getDouble(scoreP1CoverageColumn);
                    double p2c = rs.getDouble(scoreP2CoverageColumn);
                    double pmz = rs.getDouble(exp_mzColumn);
                    double f = 1;
                    if (pepSeq2 != null && !pepSeq2.isEmpty() && p1c + p2c > 0) {
                        ////                    double max = Math.max(p1c,p2c);
                        ////                    double min = Math.min(p1c,p2c);
                        ////                    f = min/(p1c+p2c);
                        ////                    score = score * f;
                        scoreRatio = (p1c) / (p1c + p2c + 1);
                        //                    if (p1c <3 || p2c <3) 
                        //                        continue;
                    }
                    
                    double peptide1score = score*scoreRatio;
                    double peptide2score = score*(1-scoreRatio);
                    

                    String xl = rs.getString(xlColumn);
                    if (xl == null) {
                        xl = "";
                    }
                    double calc_mass = rs.getDouble(calc_massColumn);
                    double pep1mass = rs.getDouble(pep1massColumn);
                    double pep2mass = rs.getDouble(pep2massColumn);
                    int search_id = rs.getInt(search_idColumn);
                    int scan_id = rs.getInt(spectrum_idColumn);
                    if (pepSeq2 != null && pepSeq2.matches("^X-?[0-9\\.]*$")) {
                        double mass = Double.parseDouble(pepSeq2.substring(1));
                        AminoModification am = new AminoModification(pepSeq2, AminoAcid.A, mass-18.0105647);
                        m_conf.addKnownModification(am);
                        getConfig(searchId).addKnownModification(am);
                        tmmodcount.add(((int)mass*10)/10.0);
                    }

                    int mgxrank = rs.getInt(mgxrankColumn);
                    boolean cleavclpep1fragmatched = rs.getBoolean(cleavclpep1fragmatchedColumn);
                    boolean cleavclpep2fragmatched = rs.getBoolean(cleavclpep2fragmatchedColumn);
                    java.sql.Array subscores = rs.getArray(subscoresColumn);

                    
                    
                    DBPSM psm = setUpDBPSM(psmID, run, scan, pep1ID, pep2ID, pepSeq1, pepSeq2, peplen1, peplen2, site1, site2, isDecoy1, isDecoy2, charge, score, protein1ID, accession1, description1, protein2ID, accession2, description2, pepPosition1, pepPosition2, sequence1, sequence2, peptide1score, peptide2score, spectrum_charge, xl, pmz, calc_mass, pep1mass, pep2mass, search_id, scan_id);
                    psm.setRank(rank);
//                    if (subscores != null) {
//                        Float[] scoresf = (Float[])subscores.getArray();
//                        for (int i =0; i<scorenames.size(); i++) {
//                            psm.getOtherInfo().put(scorenames.get(i),scoresf[i]);
//                        }
//                    }
                    psm.getOtherInfo().put("scoreP1Coverage",p1c);
                    psm.getOtherInfo().put("scoreP2Coverage",p2c);
                    psm.getOtherInfo().put("cleavclpep1fragmatched",cleavclpep1fragmatched ? 1:0);
                    psm.getOtherInfo().put("cleavclpep2fragmatched",cleavclpep2fragmatched ? 1:0);
                    psm.getOtherInfo().put("cleavclpepfragmatched",(cleavclpep2fragmatched ? 1:0) + (cleavclpep1fragmatched ? 1:0));
                    psm.getOtherInfo().put("deltaScore",rs.getDouble(deltaScoreColumn));
                    psm.getOtherInfo().put("ScoreDivDelta",rs.getDouble(scoreColumn)/rs.getDouble(deltaScoreColumn));
                    psm.getOtherInfo().put("linkSiteDelta",rs.getDouble(LinkSiteDeltaColumn));
                    psm.getOtherInfo().put("mgxrank",rs.getDouble(mgxrankColumn));
                    psm.getOtherInfo().put("PrecursorIntensity",rs.getDouble(precIntensityColumn));
                    psm.getOtherInfo().put("RetentionTime",rs.getDouble(retentiontimeColumn));
                    psm.setCrosslinkerModMass(rs.getDouble(clMassColumn));

                    if  (autovalidated) {
                        psm.setAutoValidated(autovalidated);
                        if (flagAutoValidated) {
                            psm.setPositiveGrouping("AutoValidated");
                        }
                    }

                    String modLoockup = pepSeq1;
                    if (pepSeq2 != null)
                        modLoockup+=pepSeq2;
                    DBRunConfig psmconfig  =getConfig(psm.getSearchID());
                    Sequence m = new Sequence(modLoockup, psmconfig);
                    for (AminoAcid aa : m) {
                        if (aa instanceof AminoModification) {
                            if (psmconfig.getVariableModifications().contains(aa)) {
                                psm.setHasVarMods(true);
                            }
                            if (psmconfig.getFixedModifications().contains(aa)) {
                                psm.setHasFixedMods(true);
                            }
                        }
                    }
                    if (psm.hasVarMods()) {
                        if (markModifications())
                            psm.addNegativeGrouping("Modified");
                    }

                    if (markRun())
                        psm.getAdditionalFDRGroups().add(psm.getRun());

                    if (markSearchID())
                        psm.getAdditionalFDRGroups().add(""+psm.getSearchID());
                    
                    if (psm.hasVarMods()) {
                        if (markModifications())
                            psm.addNegativeGrouping("Modified");
                    }
                    

                    int scanIndex = rs.getInt(scanIndexColumn);
                    if (!rs.wasNull()) {
                        psm.setFileScanIndex(scanIndex);
                    }
                    String peaklistfile = rs.getString(peaklistfileColumn);
                    if (!rs.wasNull()) {
                        psm.setPeakListName(peaklistfile);
                    }
                    
//                    if (isMultpleSearches) {
//                        psm.addPositiveGrouping("" + searchId);
//                    }
                    if (total % 10000 == 0) {
                        Logger.getLogger(this.getClass().getName()).log(Level.INFO, 
                                "go through the results ({0}) Search {1} of {2}", 
                                new Object[]{allPSMs.size(), currentsearch+1, searchIds.length});
                    }

                }
                m_search_ids.add(searchId);

            } catch (SQLException sex) {
                if (tries < 5) {
                    readDB(searchIds, filter, topOnly, skip, tries + 1);
                }
                Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, "Repeatedly (" + total + ") failed to read from the database giving up now", sex);
                return;
            }
            Logger.getLogger(this.getClass().getName()).log(Level.INFO, "Count results : " + total);
            
            if (tmmodcount.size()>10) {
                PeptidePair.ISTARGETED = false;
            }
        }

    }

    public void readDBSteps(int[] searchIds, String filter, boolean topOnly, ArrayList<Long> skip, int tries) throws SQLException {

        class psminfo {

            ArrayList<Integer> peptide1Positions = new ArrayList<>();
            ArrayList<Integer> peptide2Positions = new ArrayList<>();
            ArrayList<Integer> protein1IDs = new ArrayList<>();
            ArrayList<Integer> protein2IDs = new ArrayList<>();
            int search_id;
            int spectrum_id;
            int psm_id;
            int peptide1_id;
            int peptide2_id;
            String run;
            String scan;
            int exp_charge;
            int match_charge;
            double score;
            double exp_mz;
            double calc_mass;
            boolean autvalidated;
            int crosslinker_id;
            double scoreRatio;
            int linkSite1;
            int linkSite2;
            HashMap<String, Double> filterables = new HashMap<>();
        };

        class pepinfo {

            String sequence;
            double mass;
            short length;

        }

        class protinfo {

            String sequence;
            String name;
            String accession;
            String description;
            boolean isDecoy;
            int length;

            public protinfo(String sequence, String name, String accession, String description, boolean isDecoy, int length) {
                this.sequence = sequence;
                this.name = name;
                this.accession = accession;
                this.description = description;
                this.isDecoy = isDecoy;
                this.length = length;
            }
        }

        HashMap<Integer, protinfo> proteinIDs = new HashMap<>();
        HashMap<Integer, pepinfo> peptideIDs = new HashMap<>();
        HashMap<Integer, String> crosslinker = new HashMap<>();
        ArrayList<psminfo> psms = new ArrayList<>();

        if (!ensureConnection()) {
            return;
        }

        setConfig(new DBRunConfig(getDBConnection()));
        getConfig().readConfig(searchIds);
        boolean isTargted = false;

        for (CrossLinker xl : getConfig().getCrossLinker()) {
            if (xl.getName().toLowerCase().contains("targetmodification")) {
                isTargted = true;
            }
        }

        String dbNameQuerry = "Select id,name from search_sequencedb ss inner join sequence_file sf ON ss.search_id in (" + RArrayUtils.toString(searchIds, ",") + ") and ss.seqdb_id = sf.id";

        Statement dbnSt = getDBConnection().createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
        ResultSet rsDBn = dbnSt.executeQuery(dbNameQuerry);
        StringBuilder sb = new StringBuilder();
        sequenceDBs = new ArrayList<String>(1);
        HashSet<Integer> dbIds = new HashSet<Integer>();

        while (rsDBn.next()) {
            int id = rsDBn.getInt(1);
            if (!dbIds.contains(id)) {
                sequenceDBs.add(rsDBn.getString(2));
                dbIds.add(id);
            }
        }

        PeptidePair.ISTARGETED = isTargted;
        boolean xi3db = false;
        try {
            Statement stm = getDBConnection().createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
            stm.setFetchSize(100);
            ResultSet rs = stm.executeQuery("SELECT * FROM spectrum_source limit 0;");
            rs.close();
            xi3db = true;
        } catch (SQLException ex) {
            xi3db = false;
        }

        if (!m_search_ids.isEmpty()) {
            if ((filter == null && !filterSetting.isEmpty()) || (!filter.contentEquals(filterSetting))) {
                filterSetting = "--- Mixed Filter --- ";
            }
        } else if (filter != null && !filter.isEmpty()) {
            filterSetting = filter;
        }

        boolean multipleSearches = searchIds.length >1;
        for (int searchId : searchIds) {
            // wee read that one already
            if (m_search_ids.contains(searchId)) {
                continue;
            }

            String matchQuerry;
            matchQuerry
                    = "SELECT * FROM (SELECT "
                    + "sm.id AS psmID \n"
                    + " , sm.precursor_charge AS calc_charge \n"
                    + " , sm.score \n"
                    + " , hp.peptide_position + 1 AS pepPosition  \n"
                    + " , s.precursor_charge as exp_charge \n"
                    + " , hp.display_site as display_site\n"
                    + " , hp.protein_id as proteinid\n"
                    + " , mp.peptide_id as peptideid\n"
                    + " , mp.link_position\n"
                    + " , ss.name as run_name, scan_number \n"
                    + " , scorepeptide1matchedconservative AS scoreP1Coverage \n"
                    + " , scorepeptide2matchedconservative AS scoreP2Coverage \n"
                    + " , scoredelta  AS deltaScore\n"
                    + " , scorelinksitedelta AS LinkSiteDelta\n"
                    + " , s.precursor_mz AS exp_mz\n"
                    + " , sm.calc_mass\n"
                    + " , sm.precursor_charge as match_charge\n"
                    + " , sm.search_id\n"
                    + " , sm.spectrum_id\n"
                    + " , s.precursor_charge\n"
                    + " , autovalidated\n"
                    + " , mp.crosslinker_id \n"
                    + " , mp.match_type \n"
                    + " \n"
                    + "FROM \n"
                    + "  (SELECT * FROM Spectrum_match WHERE Search_id = " + searchId +  (topOnly ? " AND dynamic_rank = 't'":"")+ "  AND score>0) sm \n"
                    + "     INNER JOIN \n"
                    + "   matched_peptide mp on sm.id = mp.match_id \n"
                    + "     INNER JOIN \n"
                    + "   has_protein hp ON mp.peptide_id = hp.peptide_id \n"
                    + "     INNER JOIN  \n"
                    + "   spectrum s ON sm.spectrum_id = s.id \n"
                    + "     INNER JOIN  \n"
                    + "   spectrum_source ss ON s.source_id = ss.id \n"
                    + " WHERE  \n"
                    + " (s.precursor_charge = sm.precursor_charge OR sm.precursor_charge < 6) \n"
                    + " ORDER BY sm.id,match_type DESC)  i \n"
                    + (filter == null || filter.isEmpty() ? "" : " WHERE ( " + filter + " )");

            Logger.getLogger(this.getClass().getName()).log(Level.INFO, "read base psm from db");
            Logger.getLogger(this.getClass().getName()).log(Level.INFO, matchQuerry);
            getDBConnection().setAutoCommit(false);
            Statement stm = getDBConnection().createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
            stm.setFetchSize(100);
            ResultSet rs = stm.executeQuery(matchQuerry);
            Logger.getLogger(this.getClass().getName()).log(Level.INFO, "go through the results");

            int psmIDColumn = rs.findColumn("psmID");
            int siteColumn = rs.findColumn("link_position");
            int calc_chargeColumn = rs.findColumn("calc_charge");
            int scoreColumn = rs.findColumn("score");
            int pepPositionColumn = rs.findColumn("pepPosition");
//            int score_ratioColumn = rs.findColumn("score_ratio");
            int exp_chargeColumn = rs.findColumn("exp_charge");
            //int display_siteColumn = rs.findColumn("display_site");
            int proteinidColumn = rs.findColumn("proteinid");
            int peptideidColumn = rs.findColumn("peptideid");
            int run_nameColumn = rs.findColumn("run_name");
            int scan_numberColumn = rs.findColumn("scan_number");
            int scoreP1CoverageColumn = rs.findColumn("scoreP1Coverage");
            int scoreP2CoverageColumn = rs.findColumn("scoreP2Coverage");
            int deltaScoreColumn = rs.findColumn("deltaScore");
            int LinkSiteDeltaColumn = rs.findColumn("LinkSiteDelta");
            int exp_mzColumn = rs.findColumn("exp_mz");
            int calc_massColumn = rs.findColumn("calc_mass");
            int search_idColumn = rs.findColumn("search_id");
            int spectrum_idColumn = rs.findColumn("spectrum_id");
            int xlColumn = rs.findColumn("crosslinker_id");
            int matchTypeColumn = rs.findColumn("match_type");

            int total = 0;
            psminfo psm = null;
            try {
                int lastPSMID = -1;
                while (rs.next()) {
                    total++;
                    if (total % 5000 == 0)
                        Logger.getLogger(this.getClass().getName()).log(Level.INFO, "go through the results (" + total +")");

                    long ipsmID = rs.getLong(psmIDColumn);
                    if (skip.contains(ipsmID)) {
                        continue;
                    }

                    if (ipsmID != lastPSMID) {
                        psm = new psminfo();
                        psms.add(psm);
                        psm.match_charge = rs.getInt(calc_chargeColumn);
                        psm.search_id = rs.getInt(search_idColumn);
                        psm.spectrum_id = rs.getInt(spectrum_idColumn);
                        psm.calc_mass = rs.getInt(calc_massColumn);
                        psm.score = rs.getDouble(scoreColumn);
                        psm.exp_charge = rs.getInt(exp_chargeColumn);
                        psm.run = rs.getString(run_nameColumn);
                        psm.scan = rs.getString(scan_numberColumn);
                        psm.exp_mz = rs.getDouble(exp_mzColumn);
                        psm.filterables.put("delta", rs.getDouble(deltaScoreColumn));
                        psm.filterables.put("LinkSiteDelta", rs.getDouble(LinkSiteDeltaColumn));
                        psm.filterables.put("Peptide1Fragments", rs.getDouble(scoreP1CoverageColumn));
                        psm.filterables.put("Peptide1Fragments", rs.getDouble(scoreP2CoverageColumn));
                        psm.crosslinker_id = rs.getInt(xlColumn);
                        crosslinker.put(psm.crosslinker_id, "");
                    }
                    int match_type = rs.getInt(matchTypeColumn);
                    int peptide = match_type - 1;
                    int proteinID = rs.getInt(proteinidColumn);
                    int peptideID = rs.getInt(peptideidColumn);
                    int pepPos = rs.getInt(pepPositionColumn);

                    if (match_type == 1) {
                        psm.linkSite1 = rs.getInt(siteColumn);
                        if (psm.peptide1_id == 0) {
                            psm.peptide1_id = peptideID;
                            peptideIDs.put(peptideID, null);
                        }
                        psm.peptide1Positions.add(pepPos);
                        psm.protein1IDs.add(proteinID);
                        proteinIDs.put(proteinID,null);
                    } else if (match_type == 2) {
                        psm.linkSite2 = rs.getInt(siteColumn);
                        if (psm.peptide2_id == 0) {
                            psm.peptide2_id = peptideID;
                            peptideIDs.put(peptideID, null);
                        }
                        psm.peptide2Positions.add(pepPos);
                        psm.protein2IDs.add(proteinID);
                        proteinIDs.put(proteinID,null);
                    }
                }
                rs.close();

                Logger.getLogger(this.getClass().getName()).log(Level.INFO, "read proteins from db");
                String proteinQuerry = "SELECT id,name, accession_number, description, sequence, is_decoy, header FROM protein WHERE id in (" + RArrayUtils.toString(proteinIDs.keySet(), ",") + ");";
                rs = getDBConnection().createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY).
                        executeQuery(proteinQuerry);

                while (rs.next()) {
                    String description = rs.getString(4);
                    if (description == null) {
                        description = rs.getString(7);
                    }
                    if (description == null) {
                        description = rs.getString(2);
                    }
                    if (description == null) {
                        description = rs.getString(3);
                    }
                    String sequence = rs.getString(5).replace("[^A-Z]", "");
                    int length = sequence.length();

                    proteinIDs.put(rs.getInt(1), new protinfo(sequence, rs.getString(2), rs.getString(3), rs.getString(4), rs.getBoolean(6), length));
                }
                rs.close();

                Logger.getLogger(this.getClass().getName()).log(Level.INFO, "read peptides from db");
                String peptideQuerry = "SELECT id, sequence, mass, peptide_length FROM peptide WHERE id in (" + RArrayUtils.toString(peptideIDs.keySet(), ",") + ");";
                rs = getDBConnection().createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY).
                        executeQuery(peptideQuerry);

                while (rs.next()) {
                    pepinfo pi = new pepinfo();
                    int id = rs.getInt(1);
                    pi.sequence = rs.getString(2);
                    pi.mass = rs.getDouble(3);
                    pi.length = rs.getShort(4);

                    peptideIDs.put(id, pi);
                }
                rs.close();

                Logger.getLogger(this.getClass().getName()).log(Level.INFO, "read crosslinker from db");
                String crosslinkerQuerry = "SELECT id, name FROM crosslinker WHERE id in (" + RArrayUtils.toString(crosslinker.keySet(), ",") + ");";
                rs = getDBConnection().createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY).
                        executeQuery(proteinQuerry);

                while (rs.next()) {
                    crosslinker.put(rs.getInt(1), rs.getString(2));
                }
                rs.close();

                Logger.getLogger(this.getClass().getName()).log(Level.INFO, "assemble info");

                for (psminfo p : psms) {
                    if (p.protein2IDs.size() > 0) {
                        for (int prot1id = 0; prot1id < p.protein1IDs.size(); prot1id++) {
                            for (int prot2id = 0; prot2id < p.protein2IDs.size(); prot2id++) {
                            }
                        }
                    }
                }
                //setUpDBPSM(psmID, run, scan, pep1ID, pep2ID, pepSeq1, pepSeq2, peplen1, peplen2, site1, site2, isDecoy1, isDecoy2, charge, score, protein1ID, accession1, description1, protein2ID, accession2, description2, pepPosition1, pepPosition2, sequence1, sequence2, scoreRatio, spectrum_charge, xl, pmz, calc_mass, pep1mass, pep2mass, search_id, scan_id);

                m_search_ids.add(searchId);

            } catch (SQLException sex) {
                if (tries < 5) {
                    readDB(searchIds, filter, topOnly, skip, tries + 1);
                }
                Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, "Repeatedly (" + total + ") failed to read from the database giving up now", sex);
                return;
            }
            Logger.getLogger(this.getClass().getName()).log(Level.INFO, "Count results : " + total);

        }

    }

    protected DBPSM setUpDBPSM(String psmID, String run, String scan, long pep1ID, long pep2ID, String pepSeq1, String pepSeq2, int peplen1, int peplen2, int site1, int site2, boolean isDecoy1, boolean isDecoy2, int charge, double score, long protein1ID, String accession1, String description1, long protein2ID, String accession2, String description2, int pepPosition1, int pepPosition2, String sequence1, String sequence2, double peptide1score, double peptide2score, int spectrum_charge, String xl, double pmz, double calc_mass, double pep1mass, double pep2mass, int search_id, int scan_id) {
        DBPSM psm = (DBPSM)addMatch(psmID, run, scan, pep1ID, pep2ID, pepSeq1, pepSeq2, peplen1, peplen2, site1, site2, isDecoy1, isDecoy2, charge, score, protein1ID, accession1, description1, protein2ID, accession2, description2, pepPosition1, pepPosition2, sequence1, sequence2, peptide1score, peptide2score, spectrum_charge == -1?"Unknow Charge" : null, xl);
        if (spectrum_charge == -1) {
            psm.setNegativeGrouping(" UnknownCharge");
        }
        
        psm.setExperimentalMZ(pmz);
        psm.setCalcMass(calc_mass);
        psm.setExpCharge((byte)spectrum_charge);
        Double pepMass = psm.getPeptide1().getMass();
        if (pepMass == null || pepMass <= 0 || Double.isNaN(pepMass)) {
            if (psm.getPeptide1().getSequence().contentEquals(pepSeq1)) {
                pepMass = pep1mass;
            } else {
                pepMass = pep2mass;
            }

            if (pepMass != null && pepMass != 0 && pepMass != Double.NaN) {

                psm.getPeptide1().setMass(pepMass);
            }
        }
        if (sequence2 != null) {
            pepMass = psm.getPeptide2().getMass();
            if (pepMass == null || pepMass <= 0 || Double.isNaN(pepMass)) {
                if (psm.getPeptide2().getSequence().contentEquals(pepSeq1)) {
                    pepMass = pep1mass;
                } else {
                    pepMass = pep2mass;
                }
                if (pepMass != null && pepMass != 0 && pepMass != Double.NaN) {
                    psm.getPeptide2().setMass(pepMass);
                }
            }
        }

        if (Double.isNaN(psm.getPeptide1().getMass()) || (sequence2 != null && Double.isNaN(psm.getPeptide2().getMass()))) {
            Logger.getLogger(this.getClass().getName()).log(Level.WARNING, "No mass for peptide");
        }

//        if (psm.getCalcMass() - psm.getPeptide1().getMass() - psm.getPeptide2().getMass() <0) {
//            Logger.getLogger(this.getClass().getName()).log(Level.WARNING, "wrong mass");
//        }
        psm.setSearchID(search_id);
        psm.setScanID(scan_id);
        return psm;
//                    if (!psm.isLinear()) {
//                        psm.setCrosslinker(rs.getString(xlColumn));
//                    }
    }

//    public void readDB(int[] searchIds, String filter, boolean topOnly, ArrayList<Long> skip, int tries, boolean crosslinksOnly) throws SQLException {
//
//        if (!ensureConnection()) {
//            return;
//        }
//
//        setConfig(new DBRunConfig(getDBConnection()));
//        getConfig().readConfig(searchIds);
//        boolean isTargted = false;
//
//        for (CrossLinker xl : getConfig().getCrossLinker()) {
//            if (xl.getName().toLowerCase().contains("targetmodification")) {
//                isTargted = true;
//            }
//        }
//
//        String dbNameQuerry = "Select id,name from search_sequencedb ss inner join sequence_file sf ON ss.search_id in (" + RArrayUtils.toString(searchIds, ",") + ") and ss.seqdb_id = sf.id";
//
//        Statement dbnSt = getDBConnection().createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
//        ResultSet rsDBn = dbnSt.executeQuery(dbNameQuerry);
//        StringBuilder sb = new StringBuilder();
//        sequenceDBs = new ArrayList<String>(1);
//        HashSet<Integer> dbIds = new HashSet<Integer>();
//
//        while (rsDBn.next()) {
//            int id = rsDBn.getInt(1);
//            if (!dbIds.contains(id)) {
//                sequenceDBs.add(rsDBn.getString(2));
//                dbIds.add(id);
//            }
//        }
//
//        PeptidePair.ISTARGETED = isTargted;
//        boolean xi3db = false;
//        try {
//            Statement stm = getDBConnection().createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
//            stm.setFetchSize(100);
//            ResultSet rs = stm.executeQuery("SELECT * FROM spectrum_source limit 0;");
//            rs.close();
//            xi3db = true;
//        } catch (SQLException ex) {
//            xi3db = false;
//        }
//
//        if (!m_search_ids.isEmpty()) {
//            if ((filter == null && !filterSetting.isEmpty()) || (!filter.contentEquals(filterSetting))) {
//                filterSetting = "--- Mixed Filter --- ";
//            }
//        } else if (filter != null && !filter.isEmpty()) {
//            filterSetting = filter;
//        }
//
//        for (int searchId : searchIds) {
//
//            // wee read that one already
//            if (m_search_ids.contains(searchId)) {
//                continue;
//            }
//
//            Logger.getLogger(this.getClass().getName()).log(Level.INFO, "Read search " + searchId);
//
//            String matchQuerry;
//            matchQuerry
//                    = "SELECT * FROM (SELECT sm.id AS psmID, "
//                    + "p1.sequence AS pepSeq1, "
//                    + "p2.sequence AS pepSeq2, "
//                    + "p1.peptide_length as peplen1, "
//                    + "p2.peptide_length as peplen2, "
//                    + "mp1.link_position + 1 as site1, "
//                    + "mp2.link_position + 1 as site2, "
//                    + "pr1.is_decoy AS isDecoy1, "
//                    + "pr2.is_decoy AS isDecoy2, "
//                    + "sm.precursor_charge AS calc_charge, "
//                    + "sm.score, "
//                    + "pr1.accession_number AS accession1, "
//                    + "CASE WHEN pr1.description is not null THEN pr1.description WHEN pr1.name is not null THEN pr1.name ELSE pr1.accession_number END AS  description1, "
//                    + "pr2.accession_number AS accession2, "
//                    + "CASE WHEN pr2.description is not null THEN pr2.description WHEN pr2.name is not null THEN pr2.name ELSE pr2.accession_number END AS  description2, "
//                    + "hp1.peptide_position + 1 AS pepPosition1,  "
//                    + "hp2.peptide_position + 1 AS pepPosition2, "
//                    + "CASE WHEN p2.sequence IS NULL THEN 1 ELSE (4.0/5.0+(p1.peptide_length/(p1.peptide_length+p2.peptide_length)))/2 END AS score_ratio "
//                    + " , s.precursor_charge as exp_charge "
//                    + " , hp1.display_site as display_site1"
//                    + " , hp2.display_site as display_site2"
//                    + " , pr1.id as protein1id"
//                    + " , pr2.id as protein2id"
//                    + " , p1.id as peptide1id"
//                    + " , p2.id as peptide2id"
//                    + " , v.run_name, v.scan_number "
//                    + " , pr1.sequence AS protein1sequence "
//                    + " , pr2.sequence AS protein2sequence"
//                    + (filter != null && filter.toLowerCase().contains("scorep1coverage")
//                    ? (xi3db ? ", scorepeptide1matchedconservative" : " , coverage1 ") + " AS scoreP1Coverage "
//                    : " , 0 AS scoreP1Coverage ")
//                    + (filter != null && filter.toLowerCase().contains("scorep2coverage")
//                    ? (xi3db ? ", scorepeptide2matchedconservative" : " , coverage2 ") + " AS scoreP2Coverage "
//                    : " , 0 AS scoreP2Coverage ")
//                    + (filter != null && filter.toLowerCase().replace("linksitedelta", "").contains("delta")
//                    ? (xi3db ? ", scoredelta " : " , delta ") + " AS deltaScore"
//                    : " , 0  AS deltaScore")
//                    + (filter != null && filter.toLowerCase().contains("linksitedelta")
//                    ? (xi3db ? ", scorelinksitedelta " : " , LinkSiteDelta ") + " AS LinkSiteDelta"
//                    : " , 0  AS LinkSiteDelta")
//                    + " , s.precursor_mz AS exp_mz"
//                    + " , sm.calc_mass"
//                    + " , sm.precursor_charge as match_charge"
//                    + " , p1.mass as pep1mass"
//                    + " , p2.mass as pep2mass"
//                    + " , sm.search_id"
//                    + " , sm.spectrum_id"
//                    + " , pr1.name as protein1name"
//                    + " , pr2.name as protein2name"
//                    + " , s.precursor_charge"
//                    + " , autovalidated"
//                    + " " + (xi3db ? ", cl.name " : " , v.crosslinker ") + " AS crosslinker "
//                    + " "
//                    + "FROM "
//                    + "  (SELECT * FROM Spectrum_match WHERE Search_id = " + searchId + " AND dynamic_rank = 't' ORDER BY id) sm "
//                    + "  INNER JOIN ("
//                    + (xi3db ? "SELECT ss.name as run_name, s.scan_number, sm.id as spectrum_match_id FROM (select * from spectrum_match where Search_id = " + searchId + " AND dynamic_rank = 't') sm inner join spectrum s on sm.spectrum_id = s.id INNER JOIN spectrum_source ss on s.source_id = ss.id"
//                            : "SELECT run_name, scan_number, spectrum_match_id, crosslinker FROM v_export_materialized WHERE Search_id = " + searchId + " AND dynamic_rank = 't'")
//                    + ") v "
//                    + "    ON v.spectrum_match_id = sm.id"
//                    + "    inner join "
//                    //                    + "   (SELECT * from matched_peptide where search_id  = "+ searchId +" AND match_type =1) mp1 on sm.id = mp1.match_id   LEFT OUTER JOIN "
//                    //                    + "   (SELECT * from matched_peptide where search_id  = "+ searchId +" AND match_type =2) mp2 on sm.id = mp2.match_id   INNER JOIN "
//                    + "   matched_peptide mp1 on sm.id = mp1.match_id and mp1.match_type =1  LEFT OUTER JOIN "
//                    + "   matched_peptide mp2 on sm.id = mp2.match_id AND mp2.match_type = 2  INNER JOIN "
//                    + "   peptide p1 on mp1.peptide_id = p1.id LEFT OUTER JOIN  "
//                    + "   peptide p2 on mp2.peptide_id = p2.id INNER JOIN "
//                    + "   has_protein hp1 ON mp1.peptide_id = hp1.peptide_id LEFT OUTER JOIN  "
//                    + "   has_protein hp2 ON mp2.peptide_id = hp2.peptide_id INNER JOIN "
//                    + "   protein pr1 ON hp1.protein_id = pr1.id LEFT OUTER JOIN "
//                    + "   protein pr2 ON hp2.protein_id = pr2.id"
//                    + "   INNER JOIN  "
//                    + "   spectrum s ON sm.spectrum_id = s.id "
//                    + ((!xi3db) && filter != null && filter.toLowerCase().contains("scorep1coverage")
//                    ? "   LEFT OUTER JOIN (SELECT spectrum_match_id, score AS coverage1 FROM spectrum_match_score WHERE  score_id = (select id from score where name ='peptide1 unique matched conservative')) p1c on sm.id = p1c.spectrum_match_id "
//                    : "")
//                    + ((!xi3db) && filter != null && filter.toLowerCase().contains("scorep2coverage")
//                    ? "   LEFT OUTER JOIN (SELECT spectrum_match_id, score AS coverage2 FROM spectrum_match_score WHERE  score_id = (select id from score where name ='peptide2 unique matched conservative')) p2c on sm.id = p2c.spectrum_match_id "
//                    : "")
//                    + ((!xi3db) && filter != null && filter.toLowerCase().replace("linksitedelta", "").contains("delta")
//                    ? "   LEFT OUTER JOIN (SELECT spectrum_match_id, score AS delta FROM spectrum_match_score WHERE  score_id = (select id from score where name ='delta')) deltascore on sm.id = deltascore.spectrum_match_id "
//                    : "")
//                    + ((!xi3db) && filter != null && filter.toLowerCase().contains("linksitedelta")
//                    ? "   LEFT OUTER JOIN (SELECT spectrum_match_id, score AS LinkSiteDelta FROM spectrum_match_score WHERE  score_id = (select id from score where name ='LinkSiteDelta')) LinkSiteDelta on sm.id = LinkSiteDelta.spectrum_match_id "
//                    : "")
//                    + (xi3db ? " LEFT OUTER JOIN crosslinker cl on mp1.crosslinker_id = cl.id " : "")
//                    //                    + "   LEFT OUTER JOIN spectrum_match_score p2c on sm.id = p2c.spectrum_match_id AND p2c.score_id = (select id from score where name ='peptide1 unique matched conservative') "
//                    //                    + "   LEFT OUTER JOIN spectrum_match_score sd on sm.id = sd.spectrum_match_id AND sd.score_id = (select id from score where name ='delta') "
//                    + " WHERE  "
//                    //  + " (sd.score>0) AND "
//                    + " (s.precursor_charge = sm.precursor_charge OR sm.precursor_charge<6) "
//                    + " ORDER BY sm.score DESC)  i "
//                    + (filter == null || filter.isEmpty() ? "" : " WHERE ( " + filter + " )");
//            //                + " --INNER JOIN  "
//            //                + " --  spectrum_match_score sms1 on sm.id = sms1.spectrum_match_id AND sms1.score_id = 5 LEFT OUTER JOIN "
//            //                + " --  spectrum_match_score sms2 on sm.id = sms2.spectrum_match_id AND sms2.score_id = 6 ;"; 
//
//            Logger.getLogger(this.getClass().getName()).log(Level.INFO, "read from db");
//            Logger.getLogger(this.getClass().getName()).log(Level.INFO, matchQuerry);
//            getDBConnection().setAutoCommit(false);
//            Statement stm = getDBConnection().createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
//            stm.setFetchSize(100);
//            ResultSet rs = stm.executeQuery(matchQuerry);
//            Logger.getLogger(this.getClass().getName()).log(Level.INFO, "go through the results");
//
//            int psmIDColumn = rs.findColumn("psmID");
//            int pepSeq1Column = rs.findColumn("pepSeq1");
//            int pepSeq2Column = rs.findColumn("pepSeq2");
//            int peplen1Column = rs.findColumn("peplen1");
//            int peplen2Column = rs.findColumn("peplen2");
//            int site1Column = rs.findColumn("site1");
//            int site2Column = rs.findColumn("site2");
//            int isDecoy1Column = rs.findColumn("isDecoy1");
//            int isDecoy2Column = rs.findColumn("isDecoy2");
//            int calc_chargeColumn = rs.findColumn("calc_charge");
//            int scoreColumn = rs.findColumn("score");
//            int accession1Column = rs.findColumn("accession1");
//            int description1Column = rs.findColumn("description1");
//            int accession2Column = rs.findColumn("accession2");
//            int description2Column = rs.findColumn("description2");
//            int pepPosition1Column = rs.findColumn("pepPosition1");
//            int pepPosition2Column = rs.findColumn("pepPosition2");
//            int score_ratioColumn = rs.findColumn("score_ratio");
//            int exp_chargeColumn = rs.findColumn("exp_charge");
//            int display_site1Column = rs.findColumn("display_site1");
//            int display_site2Column = rs.findColumn("display_site2");
//            int protein1idColumn = rs.findColumn("protein1id");
//            int protein2idColumn = rs.findColumn("protein2id");
//            int peptide1idColumn = rs.findColumn("peptide1id");
//            int peptide2idColumn = rs.findColumn("peptide2id");
//            int run_nameColumn = rs.findColumn("run_name");
//            int scan_numberColumn = rs.findColumn("scan_number");
//            int protein1sequenceColumn = rs.findColumn("protein1sequence");
//            int protein2sequenceColumn = rs.findColumn("protein2sequence");
//            int scoreP1CoverageColumn = rs.findColumn("scoreP1Coverage");
//            int scoreP2CoverageColumn = rs.findColumn("scoreP2Coverage");
//            int deltaScoreColumn = rs.findColumn("deltaScore");
//            int LinkSiteDeltaColumn = rs.findColumn("LinkSiteDelta");
//            int exp_mzColumn = rs.findColumn("exp_mz");
//            int calc_massColumn = rs.findColumn("calc_mass");
//            int match_chargeColumn = rs.findColumn("match_charge");
//            int pep1massColumn = rs.findColumn("pep1mass");
//            int pep2massColumn = rs.findColumn("pep2mass");
//            int search_idColumn = rs.findColumn("search_id");
//            int spectrum_idColumn = rs.findColumn("spectrum_id");
//            int xlColumn = rs.findColumn("crosslinker");
//
//            int total = 0;
//            try {
//                while (rs.next()) {
//                    total++;
//                    if (total % 100 == 0) {
//                        Logger.getLogger(this.getClass().getName()).log(Level.INFO, "lines read: " + total);
//                    }
//
//                    long ipsmID = rs.getLong(psmIDColumn);
//                    if (skip.contains(ipsmID)) {
//                        continue;
//                    }
//
//                    String psmID = Long.toString(ipsmID);
//                    String pepSeq1 = rs.getString(pepSeq1Column);
//                    String pepSeq2 = rs.getString(pepSeq2Column);
//                    int peplen1 = rs.getInt(peplen1Column);
//                    int peplen2 = rs.getInt(peplen2Column);
//                    int site1 = rs.getInt(site1Column);
//                    int site2 = pepSeq2 == null ? -1 : rs.getInt(site2Column);
//                    boolean isDecoy1 = rs.getBoolean(isDecoy1Column);
//                    boolean isDecoy2 = rs.getBoolean(isDecoy2Column);
//                    int charge = rs.getInt(calc_chargeColumn);
//                    double score = rs.getDouble(scoreColumn);
//                    String accession1 = rs.getString(accession1Column);
//                    String description1 = rs.getString(description1Column);
//                    String accession2 = rs.getString(accession2Column);
//                    String description2 = rs.getString(description2Column);
//                    int pepPosition1 = rs.getInt(pepPosition1Column);
//                    int pepPosition2 = pepSeq2 == null ? -1 : rs.getInt(pepPosition2Column);
//                    //            double coverage1  = rs.getDouble(18);
//                    //            double coverage2  = rs.getDouble(19);
//                    //            double scoreRatio = coverage1/(coverage1+coverage2);
//                    double scoreRatio = rs.getDouble(score_ratioColumn);
//                    int spectrum_charge = rs.getInt(exp_chargeColumn);
//                    long protein1ID = rs.getInt(protein1idColumn);
//                    long protein2ID = rs.getInt(protein2idColumn);
//                    long pep1ID = rs.getInt(peptide1idColumn);
//                    long pep2ID = rs.getInt(peptide2idColumn);
//                    String run = rs.getString(run_nameColumn);
//                    String scan = rs.getString(scan_numberColumn);
//                    String sequence1 = rs.getString(protein1sequenceColumn);
//                    String sequence2 = rs.getString(protein2sequenceColumn);
//
//                    if (accession1 == null) {
//                        accession1 = description1;
//                    }
//
//                    if (accession2 == null) {
//                        accession2 = description2;
//                    }
//
//                    if (sequence1 != null) {
//                        sequence1 = sequence1.replaceAll("[^A-Z]", "");
//                    }
//
//                    if (sequence2 != null) {
//                        sequence2 = sequence2.replaceAll("[^A-Z]", "");
//                    }
//
//                    double p1c = rs.getDouble(scoreP1CoverageColumn);
//                    double p2c = rs.getDouble(scoreP2CoverageColumn);
//                    double pmz = rs.getDouble(exp_mzColumn);
//                    double f = 1;
//                    if (pepSeq2 != null && !pepSeq2.isEmpty() && p1c + p2c > 0) {
//                        ////                    double max = Math.max(p1c,p2c);
//                        ////                    double min = Math.min(p1c,p2c);
//                        ////                    f = min/(p1c+p2c);
//                        ////                    score = score * f;
//                        scoreRatio = (p1c) / (p1c + p2c + 1);
//                        //                    if (p1c <3 || p2c <3) 
//                        //                        continue;
//                    }
//
//                    PSM psm = addMatch(psmID, run, scan, pep1ID, pep2ID, pepSeq1, pepSeq2, peplen1, peplen2, site1, site2, isDecoy1, isDecoy2, charge, score, protein1ID, accession1, description1, protein2ID, accession2, description2, pepPosition1, pepPosition2, sequence1, sequence2, scoreRatio, spectrum_charge == -1, rs.getString(xlColumn));
//                    if (spectrum_charge == -1) {
//                        psm.setFDRGroup(psm.getFDRGroup()+" UnknownCharge");
//                    }
//                    psm.setExperimentalMZ(pmz);
//                    psm.setCalcMass(rs.getDouble(calc_massColumn));
//                    psm.setExpCharge(rs.getByte(match_chargeColumn));
//                    Double pepMass = rs.getDouble(pep1massColumn);
//                    if (pepMass != 0) {
//                        psm.getPeptide1().setMass(pepMass);
//                    }
//                    pepMass = rs.getDouble(pep2massColumn);
//                    if (pepMass != 0) {
//                        psm.getPeptide2().setMass(pepMass);
//                    }
//                    ((DBPSM) psm).setSearchID(rs.getInt(search_idColumn));
//                    ((DBPSM) psm).setScanID(rs.getInt(spectrum_idColumn));
//
////                    if (!psm.isLinear()) {
////                        psm.setCrosslinker(rs.getString(xlColumn));
////                    }
//                }
//                m_search_ids.add(searchId);
//
//            } catch (SQLException sex) {
//                if (tries < 5) {
//                    readDB(searchIds, filter, topOnly, skip, tries + 1);
//                }
//                Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, "Repeatedly (" + total + ") failed to read from the database giving up now", sex);
//                return;
//            }
//            Logger.getLogger(this.getClass().getName()).log(Level.INFO, "Count results : " + total);
//
//        }
//
//    }

    public void getDBSizes() {
        int targetDBsize = 0;
        int decoyDBsize = 0;
        int targetProtDBsize = 0;
        int decoyProtDBsize = 0;
        int targetLinkDBsize = 0;
        int decoyLinkDBsize = 0;
        try {
            for (int id : m_search_ids) {

                Digestion digest = getConfig().getDigestion_method();
                HashSet<hashableXiPeptide> targetPeps = new HashSet<hashableXiPeptide>();
                HashSet<hashableXiPeptide> decoyPeps = new HashSet<hashableXiPeptide>();

                String dbQuerry = "SELECT sequence, is_Decoy FROM protein p inner join (SELECT DISTINCT protein_id FROM "
                        + " spectrum_match sm inner join "
                        + " matched_peptide mp ON sm.search_id = " + id + " AND sm.id = mp.match_id INNER JOIN "
                        + " has_protein hp on mp.peptide_id = hp.peptide_id) i ON p.id = i.protein_id;";

                getDBConnection().setAutoCommit(false);
                Statement stm = getDBConnection().createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
                stm.setFetchSize(100);
                ResultSet rs = stm.executeQuery(dbQuerry);
                Logger.getLogger(this.getClass().getName()).log(Level.INFO, "go through the results");
                int total = 0;
                PeptideTree xltree = new PeptideTree(getConfig().getPrecousorTolerance());
                PeptideTree lineartree = new PeptideTree(getConfig().getPrecousorTolerance());
                while (rs.next()) {

                    Sequence seq = new Sequence(rs.getString(1), getConfig());
                    boolean isdecoy = rs.getBoolean(2);

                    seq.setDecoy(isdecoy);
                    Digestion d = getConfig().getDigestion_method();
                    d.setPeptideLookup(xltree, lineartree);
                    d.setMaxMissCleavages(getConfig().getMaxMissCleavages());

                    seq.digest(getConfig().getDigestion_method(), Double.MAX_VALUE, getConfig().getCrossLinker());
                    ArrayList<rappsilber.ms.sequence.Peptide> peps = digest.digest(seq, getConfig().getCrossLinker());
                    if (isdecoy) {
                        decoyProtDBsize++;
                        seqpos:
                        for (int i = 0; i < seq.length(); i++) {
                            for (CrossLinker cl : getConfig().getCrossLinker()) {
                                if (cl.canCrossLink(seq, i)) {
                                    decoyLinkDBsize++;
                                    continue seqpos;
                                }
                            }
                        }
                    } else {
                        targetProtDBsize++;
                        seqpos:
                        for (int i = 0; i < seq.length(); i++) {
                            for (CrossLinker cl : getConfig().getCrossLinker()) {
                                if (cl.canCrossLink(seq, i)) {
                                    targetLinkDBsize++;
                                    continue seqpos;
                                }
                            }
                        }

                    }
                }
                for (rappsilber.ms.sequence.Peptide peptide : xltree) {
                    if (peptide.isDecoy()) {
                        decoyDBsize++;
                    } else {
                        targetDBsize++;
                    }
                }
            }
            setTargetDBSize(targetDBsize);
            setDecoyDBSize(decoyDBsize);

            setTargetProtDBSize(targetProtDBsize);
            setDecoyProtDBSize(decoyProtDBsize);

            setTargetLinkDBSize(targetLinkDBsize);
            setDecoyLinkDBSize(decoyLinkDBsize);

        } catch (SQLException ex) {
            Logger.getLogger(DBinFDR.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    public void clearDBValidation(String validation) throws SQLException {
        if (!ensureConnection()) {
            return;
        }
        Statement stm = getDBConnection().createStatement();

        if (validation == null || validation.isEmpty()) {
            stm.execute("Update spectrum_match SET validated = null WHERE search_id in (" + RArrayUtils.toString(m_search_ids, ",") + ") AND validated is not null ");

        } else {
            String[] aval = validation.split(",");
            for (int i = 0; i < aval.length; i++) {
                aval[i] = aval[i].trim();
            }

            stm.execute("Update spectrum_match SET validated = null WHERE search_id in (" + RArrayUtils.toString(m_search_ids, ",") + ") AND validated in ('" + RArrayUtils.toString(aval, "','") + "')");
        }
        if (!getDBConnection().getAutoCommit()) {
            getDBConnection().commit();
        }

//        for (int sid : m_search_ids) {
//            stm.execute("Update spectrum_match SET validated = null WHERE search_id = " + sid  + " AND validated in ('" + RArrayUtils.toString(aval, "','") + "')");
//            if (!m_db_connection.getAutoCommit())
//                getDBConnection().commit();
//        }
    }

    public void clearDBFDR() throws SQLException {
        if (!ensureConnection()) {
            return;
        }
        Statement stm = getDBConnection().createStatement();
        for (int sid : m_search_ids) {
            stm.execute("Update spectrum_match SET peptide_fdr = null, link_fdr = null, ppi_fdr = null   WHERE search_id = " + sid);
            if (!m_db_connection.getAutoCommit()) {
                getDBConnection().commit();
            }
        }
    }

    public void writeDBPPIConfidence(FDRResult result) throws SQLException {
        //Collection<PeptidePair> peps = getFDRPeptidePairs();
        int count = 0;
        int maxUpdate = 100;
        if (!ensureConnection()) {
            return;
        }
        boolean autocomit = getDBConnection().getAutoCommit();
        getDBConnection().setAutoCommit(false);
        // also write the fdr-values
        for (PeptidePair pp : result.peptidePairFDR) {
            ProteinGroupPair pgp = pp.getFdrLink().getFdrPPI();
            double confidence = 100 * (1 - pgp.getFDR());
            for (String psmid : pp.getPSMids()) {
                updateSetConfidence.setDouble(1, confidence);
                updateSetConfidence.setLong(2, Long.parseLong(psmid));
                updateSetConfidence.addBatch();
                if ((count = (count + 1) % maxUpdate) == 0) {
                    updateSetConfidence.executeBatch();
                    Logger.getLogger(DBinFDR.class.getName()).log(Level.INFO, count+" updated");
                }
            }
        }
        if (count > 0) {
            updateSetConfidence.executeBatch();
            Logger.getLogger(DBinFDR.class.getName()).log(Level.INFO, count+" updated");
        }
        getDBConnection().commit();
        getDBConnection().setAutoCommit(autocomit);
    }

    public void writeDBPepPairConfidence(FDRResult result) throws SQLException {
        //     Collection<PeptidePair> peps = getFDRPeptidePairs();
        int count = 0;
        int maxUpdate = 100;
        if (!ensureConnection()) {
            return;
        }
        boolean autocomit = getDBConnection().getAutoCommit();
        getDBConnection().setAutoCommit(false);
        // also write the fdr-values
        for (PeptidePair pp : result.peptidePairFDR) {

            double confidence = 100 * (1 - pp.getFDR());
            for (String psmid : pp.getPSMids()) {
                updateSetConfidence.setDouble(1, confidence);
                updateSetConfidence.setLong(2, Long.parseLong(psmid));
                updateSetConfidence.addBatch();
                if ((count = (count + 1) % maxUpdate) == 0) {
                    updateSetConfidence.executeBatch();
                    Logger.getLogger(DBinFDR.class.getName()).log(Level.INFO, count+" updated");
                }
            }
        }
        if (count > 0) {
            updateSetConfidence.executeBatch();
            Logger.getLogger(DBinFDR.class.getName()).log(Level.INFO, count+" updated");
        }
        getDBConnection().commit();
        getDBConnection().setAutoCommit(autocomit);
    }

    public void writeDBLinkConfidence(FDRResult result) throws SQLException {
        //Collection<PeptidePair> peps = getFDRPeptidePairs();
        int count = 0;
        int maxUpdate = 100;
        if (!ensureConnection()) {
            return;
        }
        boolean autocomit = getDBConnection().getAutoCommit();
        getDBConnection().setAutoCommit(false);
        // also write the fdr-values
        for (PeptidePair pp : result.peptidePairFDR) {
            if (!pp.isLinear()) {
                ProteinGroupLink l = pp.getFdrLink();
                double confidence = 100 * (1 - l.getFDR());
                for (String psmid : pp.getPSMids()) {
                    updateSetConfidence.setDouble(1, confidence);
                    updateSetConfidence.setLong(2, Long.parseLong(psmid));
                    updateSetConfidence.addBatch();
                    if ((count = (count + 1) % maxUpdate) == 0) {
                        updateSetConfidence.executeBatch();
                        Logger.getLogger(DBinFDR.class.getName()).log(Level.INFO, count+" updated");
                    }
                }
            }
        }
        if (count > 0) {
            updateSetConfidence.executeBatch();
            Logger.getLogger(DBinFDR.class.getName()).log(Level.INFO, count+" updated");
        }
        getDBConnection().commit();
        getDBConnection().setAutoCommit(autocomit);
    }

    public void writeDB(String validate, boolean overwrite, boolean writeFDR, FDRResult result, boolean within, boolean  between) throws SQLException {
        //Collection<PeptidePair> peps = getFDRPeptidePairs();
        if (!ensureConnection()) {
            return;
        }
        int count = 0;
        int maxUpdate = 100;
        boolean autocomit = getDBConnection().getAutoCommit();
        getDBConnection().setAutoCommit(false);
        if (validate != null) {
            if (overwrite) { // overwrite validation

                // just validate
                for (PSM psm : result.psmFDR) {
                    String psmid = psm.getPsmID();
                    
                    if (psm.isLinear() || (psm.isBetween() && between) || (psm.isInternal() && within)) {
                        updateValidateOverWrite.setString(1, validate);
                        updateValidateOverWrite.setLong(2, Long.parseLong(psmid));
                        updateValidateOverWrite.addBatch();
                        if ((count = (count + 1) % maxUpdate) == 0) {
                            updateValidateOverWrite.executeBatch();
                            Logger.getLogger(DBinFDR.class.getName()).log(Level.INFO, count+" updated");
                    }
                    }
                }
                if (count > 0) {
                    updateValidateOverWrite.executeBatch();
                    Logger.getLogger(DBinFDR.class.getName()).log(Level.INFO, count+" updated");
                }
//                }
            } else {
                // just validate
                for (PSM psm : result.psmFDR) {
                    Long psmid = Long.parseLong(psm.getPsmID());
//                    for (PeptidePair pp : result.peptidePairFDR) {
//                        for (String psmid : pp.getPSMids()) {
                    if (psm.isLinear() || (psm.isBetween() && between) || (psm.isInternal() && within)) {
                        updateValidateNonOverWrite.setString(1, validate);
                        updateValidateNonOverWrite.setLong(2, psmid);
                        updateValidateNonOverWrite.addBatch();
                        if ((count = (count + 1) % maxUpdate) == 0) {
                            updateValidateNonOverWrite.executeBatch();
                            Logger.getLogger(DBinFDR.class.getName()).log(Level.INFO, count+" updated");
                         }
                    }
//                        }
                }
                
                if (count > 0) {
                    updateValidateNonOverWrite.executeBatch();
                    Logger.getLogger(DBinFDR.class.getName()).log(Level.INFO, count+" updated");
                }
//                }

            }
        } else { // do not validate
        }
        getDBConnection().commit();
        getDBConnection().setAutoCommit(autocomit);
    }

    /**
     * adds a psm to the list folds up the scores to peptidespairs links
     * proteinpairs and proteins
     *
     * @param psmID
     * @param pepid1
     * @param pepid2
     * @param peplen1
     * @param peplen2
     * @param site1
     * @param site2
     * @param charge
     * @param score
     * @param proteinId1
     * @param proteinId2
     * @param pepPosition1
     * @param pepPosition2
     * @param scoreRation
     * @return a peptide pair that is supported by the given match
     */
    @Override
    public PSM addMatch(String psmID, org.rappsilber.fdr.entities.Peptide peptide1, org.rappsilber.fdr.entities.Peptide peptide2, int peplen1, int peplen2, int site1, int site2, int charge, double score, Protein proteinId1, Protein proteinId2, int pepPosition1, int pepPosition2, double peptide1score, double peptide2score, String isSpecialCase, String crosslinker, String run, String scan) {
        org.rappsilber.fdr.entities.Peptide npepid1;
        org.rappsilber.fdr.entities.Peptide npepid2;
        int npeplen1;
        int npeplen2;
        byte nsite1;
        byte nsite2;
        Protein nproteinId1;
        Protein nproteinId2;
        int npepPosition1;
        int npepPosition2;
        int protcomp = proteinId1.compareDecoyUnAware(proteinId2);
        int pepcomp = peptide1.compare(peptide2);
        int sitecomp = (site1 - site2);
        double pep1score = peptide1score;
        double pep2score = peptide2score;
        double npep1score;
        double npep2score;

        if (protcomp < 0 || (protcomp == 0 && pepcomp < 0) || (protcomp == 0 && pepcomp == 0 && site1 < site2)) {
            npepid1 = peptide1;
            npepid2 = peptide2;
            npeplen1 = peplen1;
            npeplen2 = peplen2;
            nsite1 = (byte)site1;
            nsite2 = (byte)site2;
            nproteinId1 = proteinId1;
            nproteinId2 = proteinId2;
            npepPosition1 = pepPosition1;
            npepPosition2 = pepPosition2;
            npep1score = pep1score;
            npep2score = pep2score;

        } else {
            npepid1 = peptide2;
            npepid2 = peptide1;
            npeplen1 = peplen2;
            npeplen2 = peplen1;
            nsite1 = (byte)site2;
            nsite2 = (byte)site1;
            nproteinId1 = proteinId2;
            nproteinId2 = proteinId1;
            npepPosition1 = pepPosition2;
            npepPosition2 = pepPosition1;
            npep1score = pep2score;
            npep2score = pep1score;
        }

        if (!PSMScoreHighBetter) {
            score = 10 - (10 * score);
        }

        DBPSM psm = new DBPSM(psmID, npepid1, npepid2, nsite1, nsite2, proteinId1.isDecoy(), proteinId2.isDecoy(), (byte)charge, score, npep1score, npep2score);
        psm.setNegativeGrouping(isSpecialCase);
        psm.setRun(registerRun(run));
        psm.setScan(scan);
        psm.setCrosslinker(crosslinker);

        PSM regpsm = getAllPSMs().register(psm);

        return regpsm;
    }

    public String argList() {
        return "--dbids=A,B,C --connection= " + super.argList() + " --crosslinkonly --Filter --flagModifications";
    }

    public String argDescription() {
        return super.argDescription() + "\n"
                + "--dbids=A,B,C            List of search ids\n"
                + "--connection=            connection string to be used to\n"
                + "                         access the database\n "
                + "--crosslinkonly          only load cross-links\n"
                + "--filter                 a filter to be applied while loading\n"
                + "                         matches from the database\n"
                + "--user=                  db user name (ignore for now)\n"
                + "--password=              db password (ignore for now)\n"
                + "--flagModifications     should modified peptide make their own sub-group\n";

    }

    public String[] parseArgs(String[] argv, FDRSettings settings) {
        ArrayList<String> unknown = new ArrayList<String>();
        String ConnString = null;
        ArrayList<String> filters = new ArrayList<String>();
        String dbUserName = "";
        String dbPassword = "";
        
        argv = super.parseArgs(argv, settings);

        for (String arg : argv) {
            if (arg.toLowerCase().startsWith("--dbids=")) {

                String[] sids = arg.substring(arg.indexOf("=") + 1).trim().split(",");
                searchIDSetting = new int[sids.length];

                for (int i = 0; i < sids.length; i++) {
                    searchIDSetting[i] = Integer.parseInt(sids[i]);
                }

            } else if (arg.toLowerCase().startsWith("--filter=")) {
                filters.add(arg.substring(arg.indexOf("=") + 1).trim());

            } else if (arg.toLowerCase().equals("--crosslinkonly") || arg.equals("-C")) {
                filters.add("p2.accessen_number is not null");

            } else if (arg.toLowerCase().startsWith("--connection=")) {
                ConnString = arg.substring(arg.indexOf("=") + 1).trim();
            } else if (arg.toLowerCase().startsWith("--user=")) {
                dbUserName = arg.substring(arg.indexOf("=") + 1).trim();
            } else if (arg.toLowerCase().startsWith("--password=")) {
                dbPassword = arg.substring(arg.indexOf("=") + 1).trim();
            }else if (arg.toLowerCase().startsWith("--flagmodifications")) {
                setMarkModifications(true);
            } else if (arg.toLowerCase().startsWith("--validate=")) {
                command_line_auto_validate = arg.substring(arg.indexOf("=") + 1).trim();
            } else {
                unknown.add(arg);
            }
            

        }
        if (ConnString == null) {
            Logger.getLogger(DBinFDR.class.getName()).log(Level.SEVERE, "No db-connection defined");
            printUsage();
            System.exit(0);
        } else {
            Exception connectionException = null;
            try {
                setDBConnection(ConnString, dbUserName, dbPassword);
            } catch (SQLException ex) {
                connectionException = ex;
            }
            if (connectionException != null) {
                DBConnectionConfig dbc = new DBConnectionConfig();
                try {    
                    dbc.readConfig();

                    for (DBConnectionConfig.DBServer s: dbc.getServers()) {

                        if (s.name.contentEquals(ConnString)) {
                            connectionException = null;
                            ConnString = s.connectionString;
                            dbUserName = s.user;
                            dbPassword = s.password;

                            try {
                                setDBConnection(ConnString, dbUserName, dbPassword);
                                break;
                            } catch (SQLException ex) {
                                connectionException = ex;
                            }
                        }
                    }
                } catch (FileNotFoundException ex) {
                } catch (ParseException ex) {
                    Logger.getLogger(DBinFDR.class.getName()).log(Level.SEVERE, "Error parsing db config", ex);
                }
                    
            }
            
            if (connectionException != null) {
                Logger.getLogger(DBinFDR.class.getName()).log(Level.SEVERE, null, connectionException);
                printUsage();
                System.exit(0);
            }
        }

        if (searchIDSetting == null) {
            Logger.getLogger(DBinFDR.class.getName()).log(Level.SEVERE, "No search ids defined");
            printUsage();
            System.exit(0);
        }

        filterSetting = RArrayUtils.toString(filters, " AND ");
        String[] ret = new String[unknown.size()];
        ret = unknown.toArray(ret);
        return ret;
    }

    public static void main(String[] argv) throws SQLException, FileNotFoundException {

        DBinFDR ofdr = new DBinFDR();

        FDRSettings settings = new FDRSettingsImpl();
        
        String[] unknowns = ofdr.parseArgs(argv, settings);

        if (unknowns.length > 0) {
            for (String a : unknowns) {
                Logger.getLogger(DBinFDR.class.getName()).log(Level.SEVERE, "Unknown argument : " + a);
            }
            ofdr.printUsage();
            System.exit(1);
        }

        if (ofdr.getCsvOutDirSetting() != null) {
            if (ofdr.getCsvOutBaseSetting() == null) {
                ofdr.setCsvOutBaseSetting("FDR_" + RArrayUtils.toString(ofdr.searchIDSetting, "_"));
            }
        }

        if (ofdr.getCsvOutBaseSetting() != null) {
            if (ofdr.getCsvOutDirSetting() == null) {
                ofdr.setCsvOutDirSetting(".");
            }

            System.out.println("writing results to " + ofdr.getCsvOutDirSetting() + "/" + ofdr.getCsvOutBaseSetting() + "*");
            System.out.flush();
        }

        Logger.getLogger(DBinFDR.class.getName()).log(Level.INFO, "Read datafrom db");
        ofdr.readDB();

        Logger.getLogger(DBinFDR.class.getName()).log(Level.INFO, "Calculate FDR");
        FDRResult result = ofdr.calculateWriteFDR(ofdr.getCsvOutDirSetting(), ofdr.getCsvOutBaseSetting(), ",", settings);

//        Logger.getLogger(DBinFDR.class.getName()).log(Level.INFO, "Result-summary:" + ofdr.summaryString());
//        
//        if (ofdr.getCsvOutBaseSetting() != null) {
//            
//            ofdr.writeFiles(ofdr.getCsvOutDirSetting(), ofdr.getCsvOutDirSetting(), ",");
//        }
        //System.out.println(ofdr.summaryString(result));
        if (ofdr.command_line_auto_validate != null) {
            ofdr.writeDB(ofdr.command_line_auto_validate, true, false, result, true, true);
        }
        System.exit(0);

    }

    @Override
    protected ArrayList<String> getLinkOutputLine(ProteinGroupLink l) {
        ArrayList<String> ret = super.getLinkOutputLine(l);
        HashSet<Integer> searchIDs = new HashSet<Integer>();
        ret.add(getLinkWindow(l, 0, 20));
        ret.add(getLinkWindow(l, 1, 20));
        Collection<PeptidePair> cpeppairs = l.getPeptidePairs();
        HashSet<String> peps = new HashSet<String>();
        for (PeptidePair pp : cpeppairs) {
            for (org.rappsilber.fdr.entities.Peptide p : pp.getPeptides()) {
                if (p.getSequence().matches("(-\\.)?X-?[0-9.]+(\\.-)?")) {
                    peps.add(p.getSequence());
                }
            }
            for (PSM psm : pp.getAllPSMs()) {
                for (PSM upsm : psm.getRepresented()) {
                    searchIDs.add(((DBPSM) upsm).getSearchID());
                }
                //searchIDs.add(((DBPSM) psm).getSearchID());
            }
        }

        ret.add(0, RArrayUtils.toString(searchIDs, ";"));

        if (!peps.isEmpty()) {
            ret.add(RArrayUtils.toString(peps, ";"));
        }

        if (l.getProteinGroup2().accessions().matches("'?OpenModification.*")) {
            ret.add(getLinkWindow(l, 0, 20));
            HashMap<Protein, IntArrayList> pos = l.getPosition1();
            Protein p = pos.keySet().iterator().next();
            ret.add(p.getAccession());
            int pp = pos.get(p).get(0);
            ret.add("" + pp);

        } else if (l.getProteinGroup1().accessions().matches("'?OpenModification.*")) {
            ret.add(getLinkWindow(l, 1, 20));
            HashMap<Protein, IntArrayList> pos = l.getPosition2();
            Protein p = pos.keySet().iterator().next();
            ret.add(p.getAccession());
            int pp = pos.get(p).get(0);
            ret.add("" + pp);
        } else {
            ret.add("");
            ret.add("");
            ret.add("");
        }

        return ret;
    }

    protected ArrayList<String> getLinkOutputHeader() {
        ArrayList<String> ret = super.getLinkOutputHeader();
        ret.add(0, "SearchIDs");
        ret.add("LinkWindow1");
        ret.add("LinkWindow2");
        ret.add("TM_OM_Mods");
        ret.add("TM_OM_ModWindow");
        ret.add("TM_OM_RepresentiveProtein");
        ret.add("TM_OM_RepresentiveModSite");
        return ret;
    }

    @Override
    protected ArrayList<String> getPSMHeader() {
        ArrayList<String> ret = super.getPSMHeader();
//        String header = "PSMID" + seperator + "run" + seperator + "scan" 
//                + seperator + "exp charge"  + seperator + "exp m/z" + seperator + "exp mass" + seperator + "exp fractionalmass" 
//                + seperator + "match charge" + seperator + "match mass" + seperator + "match fractionalmass"   
//                
//                + seperator + "Protein1" + seperator + "Description1" + seperator + "Decoy1" + seperator + "Protein2" + seperator + "Description2" + seperator + "Decoy2" + seperator + "Peptide1" + seperator + "Peptide2" + seperator + "PeptidePosition1" + seperator + "PeptidePosition2" + seperator + "PeptideLength1" + seperator + "PeptideLength2" + seperator + "fromSite" + seperator + "ToSite" + seperator + "Charge" + seperator + "Score" + seperator + "isDecoy" + seperator + "isTT" + seperator + "isTD" + seperator + "isDD" + seperator + "fdrGroup" + seperator + "fdr" + seperator + "" + seperator + "PeptidePairFDR" + seperator + "Protein1FDR" + seperator + "Protein2FDR" + seperator + "LinkFDR" + seperator + "PPIFDR" + seperator + seperator + "peptide pair id" + seperator + "link id" + seperator + "ppi id";
        ret.add(0, "SearchID");
        ret.add(3, "exp charge");
        ret.add(4, "exp m/z");
        ret.add(5, "exp mass");
        ret.add(6, "exp fractionalmass");
        ret.add(7, "match charge");
        ret.add(8, "match mass");
        ret.add(9, "match fractionalmass");
        return ret;
    }

    protected ArrayList<String> getPSMOutputLine(PSM pp) {
        ArrayList<String> ret = super.getPSMOutputLine(pp);
        DBPSM dp = null;
        dp = (DBPSM) pp;
        double mz = pp.getExperimentalMZ();
        int charge = pp.getCharge();
        double mass = (mz - 1.00727646677) * (pp.getExpCharge() == -1 ? charge : pp.getExpCharge());
        double fraction = mass - Math.floor(mass);
        double calcfraction = pp.getCalcMass() - Math.floor(pp.getCalcMass());

        ret.add(0, "" + dp.getSearchID());
        ret.add(3, "" + pp.getExpCharge());
        ret.add(4, sixDigits.format(mz));
        ret.add(5, sixDigits.format(mass));
        ret.add(6, sixDigits.format(fraction));
        ret.add(7, "" + charge);
        ret.add(8, sixDigits.format(pp.getCalcMass()));
        ret.add(9, sixDigits.format(calcfraction));

        return ret;
    }

    @Override
    protected ArrayList<String> getPPIOutputHeader() {
        ArrayList<String> ret = super.getPPIOutputHeader();
        ret.add(0, "SearchID");
        return ret;
    }

    @Override
    protected ArrayList<String> getPPIOutputLine(ProteinGroupPair pgp) {
        HashSet<Integer> searches = new HashSet<Integer>();
        for (PeptidePair pp : pgp.getPeptidePairs()) {
            for (PSM psm : pp.getAllPSMs()) {
                for (PSM upsm : psm.getRepresented()) {
                    searches.add(((DBPSM) upsm).getSearchID());
                }
//                searches.add(((DBPSM)psm).getSearchID());
            }
        }
        ArrayList<String> ret = super.getPPIOutputLine(pgp);
        ret.add(0, RArrayUtils.toString(searches, ";"));
        return ret;
    }

    @Override
    protected ArrayList<String> getProteinGroupOutputHeader() {
        ArrayList<String> ret = super.getProteinGroupOutputHeader();
        ret.add(0, "SearchID");
        return ret;
    }

    @Override
    protected ArrayList<String> getProteinGroupOutput(ProteinGroup pg) {
        HashSet<Integer> searches = new HashSet<Integer>();
        for (PeptidePair pp : pg.getPeptidePairs()) {
            for (PSM psm : pp.getAllPSMs()) {
                for (PSM upsm : psm.getRepresented()) {
                    searches.add(((DBPSM) upsm).getSearchID());
                }
//                searches.add(((DBPSM)psm).getSearchID());
            }
        }
        ArrayList<String> ret = super.getProteinGroupOutput(pg);
        ret.add(0, RArrayUtils.toString(searches, ";"));
        return ret;
    }

    protected String getPeptideSequence(org.rappsilber.fdr.entities.Peptide p) {
        if (!m_writePrePostAA) {
            return p.getSequence();
        }

        ProteinGroup pg = p.getProteinGroup();
        HashMap<org.rappsilber.fdr.entities.Protein, IntArrayList> positions = p.getPositions();
        int len = p.length();
        HashSet<String> N = new HashSet<String>(positions.size());
        HashSet<String> C = new HashSet<String>(positions.size());

        for (org.rappsilber.fdr.entities.Protein prot : positions.keySet()) {
            Sequence seq = getSequence(prot);
            if (seq.length() == 0) {
                N.add("-");
                C.add("-");
            } else {
                for (int pos : positions.get(prot)) {
                    if (pos == 1) {
                        N.add("-");
                    } else {
                        N.add(seq.aminoAcidAt(pos - 2).toString());
                    }
                    if (pos + len > seq.length()) {
                        C.add("-");
                    } else {
                        C.add(seq.aminoAcidAt(pos + len - 1).toString());
                    }
                }
            }
        }

        return RArrayUtils.toString(N, "") + "." + p.getSequence() + "." + RArrayUtils.toString(C, "");
    }

    @Override
    protected ArrayList<String> getXlPepeptideOutputLine(PeptidePair pp) {
        ArrayList<String> sret = super.getXlPepeptideOutputLine(pp);
        ArrayList<String> ret = new ArrayList<String>(sret.size() + 5);
        HashSet<Integer> searches = new HashSet<Integer>();
        for (PSM psm : pp.getAllPSMs()) {
            for (PSM upsm : psm.getRepresented()) {
                searches.add(((DBPSM) upsm).getSearchID());
            }
        }
        ret.add(RArrayUtils.toString(searches, ";"));
        ret.addAll(sret);

        ret.add(getLinkWindow(pp, 0, 20));
        ret.add(getLinkWindow(pp, 1, 20));
        if (pp.getPeptide1().getSequence().matches("X[\\-0-9]+(\\.[0-9]+)?")) {
            ret.add(getLinkWindow(pp, 1, 20));
            ret.add(getLinkWindow(pp.getLink(), 1, 20));
        } else if (pp.getPeptide2().getSequence().matches("X[\\-0-9]+(\\.[0-9]+)?")) {
            ret.add(getLinkWindow(pp, 0, 20));
            ret.add(getLinkWindow(pp.getLink(), 0, 20));
        }
        return ret;
    }

    @Override
    protected ArrayList<String> getXLPepsHeader() {
        ArrayList<String> ret = super.getXLPepsHeader();
        ret.add(0, "SearchIDs");
        ret.add("LinkWindow1");
        ret.add("LinkWindow2");
        ret.add("TM_OM_ModWindow");
        ret.add("TM_OM_ModWindow_noVarMod");
        return ret;
    }

    @Override
    protected ArrayList<String> getLinearPepeptideOutputLine(PeptidePair pp) {
        ArrayList<String> ret = super.getLinearPepeptideOutputLine(pp);
        HashSet<Integer> searches = new HashSet<Integer>();
        for (PSM psm : pp.getAllPSMs()) {
            for (PSM upsm : psm.getRepresented()) {
                searches.add(((DBPSM) upsm).getSearchID());
            }
        }
        ret.add(0, RArrayUtils.toString(searches, ";"));
        return ret;
    }

    @Override
    protected ArrayList<String> getLinearPepsHeader() {
        ArrayList<String> ret = super.getLinearPepsHeader();
        ret.add(0, "SearchIDs");
        return ret;
    }

    public int getFastas(Collection<Integer> searchID, ArrayList<Integer> dbIDs, ArrayList<String> names) {
        int c = 0;
        ArrayList<Integer> i = new ArrayList<Integer>();
        ArrayList<String> n = new ArrayList<String>();
        try {
            Statement s = getDBConnection().createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
            ResultSet rs = s.executeQuery("SELECT  DISTINCT file_name, id FROM search_sequencedb ss INNER JOIN sequence_file sf on ss.seqdb_id = sf.id WHERE search_id in (" + RArrayUtils.toString(m_search_ids, ",") + ");");
            while (rs.next()) {
                i.add(rs.getInt(2));
                n.add(rs.getString(1));
                c++;
            }
        } catch (SQLException ex) {
            Logger.getLogger(DBinFDR.class.getName()).log(Level.SEVERE, null, ex);
            return -1;
        }
        dbIDs.addAll(i);
        names.addAll(n);
        return c;
    }

    public int getFastas(int searchID, ArrayList<Integer> dbIDs, ArrayList<String> names) {
        ArrayList<Integer> id = new ArrayList<Integer>(1);
        id.add(searchID);
        return getFastas(id, dbIDs, names);
    }

    public int getFastas(ArrayList<Integer> dbIDs, ArrayList<String> names) {
        ArrayList<Integer> ids = new ArrayList<Integer>(m_search_ids.size());
        for (Integer i : m_search_ids) {
            ids.add(i);
        }
        return getFastas(ids, dbIDs, names);
    }

    
    /**
     * @return the flagAutoValidated
     */
    public boolean isFlagAutoValidated() {
        return flagAutoValidated;
    }

    /**
     * @param flagAutoValidated the flagAutoValidated to set
     */
    public void setFlagAutoValidated(boolean flagAutoValidated) {
        this.flagAutoValidated = flagAutoValidated;
    }

    /**
     * @return the markModifications
     */
    public boolean markModifications() {
        return markModifications;
    }

    /**
     * @param mm should modification be marked
     */
    public void setMarkModifications(boolean mm) {
        this.markModifications = mm;
    }

    @Override
    public void cleanup() {
        super.cleanup(); //To change body of generated methods, choose Tools | Templates.
        m_db_autoclosetimer.cancel();
    }
    
    
            
    
}
