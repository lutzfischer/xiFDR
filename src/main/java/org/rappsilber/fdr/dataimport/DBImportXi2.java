/*
 * Copyright 2021 Lutz Fischer <lfischer@staffmail.ed.ac.uk>.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.rappsilber.fdr.dataimport;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.UUID;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.rappsilber.fdr.OfflineFDR;
import org.rappsilber.fdr.DBinFDR;
import org.rappsilber.fdr.entities.PSM;
import org.rappsilber.utils.RArrayUtils;

/**
 *
 * @author Lutz Fischer <lfischer@staffmail.ed.ac.uk>
 */
public class DBImportXi2 implements Import{ 
    
    public class Xi2Score {
        int id;
        String name;
        boolean primary_score;
        boolean higher_is_better;
        boolean major_interest;

        public Xi2Score(int id, String name, boolean primary_score, boolean higher_is_better, boolean major_interest) {
            this.id = id;
            this.name = name;
            this.primary_score = primary_score;
            this.higher_is_better = higher_is_better;
            this.major_interest = major_interest;
        }

        @Override
        public String toString() {
            ArrayList<String> qualifiers = new ArrayList<>();
            
            if (primary_score)
                qualifiers.add("*");

            if (higher_is_better)
                qualifiers.add("h>l");
            else
                qualifiers.add("l>h");
            
            if (major_interest)
                qualifiers.add(0,"*");
                qualifiers.add("*");
            
            
            return name + "(" + RArrayUtils.toString(qualifiers, ",") + ")";
        }
        
        
    }
    
    private class XiDBPSM {
        UUID id ;
        UUID search_id;
        Integer pep1_id;
        Integer pep2_id;
        Short site1;
        Short site2;
        short[] link_score_site1;
        short[] link_score_site2;
        float[] link_score;
        Short crosslinker_id;
        Double assumed_prec_mz;
        Short assumed_prec_charge;
        Double calc_mass;
        boolean top_ranking;
        boolean is_decoy;
        boolean is_dd;
        float matchscore;
        float result_score;
        float[] all_scores;
    }

    private class XiDBPeptide {
        Integer id;
        UUID search_id;
        String base_sequence;
        Double mass;
        Integer length;
        Integer[] modification_ids;
        Integer[] modification_position;
        boolean is_decoy;
        ArrayList<Integer> proteinIDs = new ArrayList<Integer>();
        ArrayList<Integer> proteinPositions = new ArrayList<Integer>();
    }
    
    private class XiDBProtein {
        Integer id;
        String accession;
        String name;
        String gen_name;
        String description;
        String sequence;
        String full_header;
        boolean is_decoy;

        @Override
        public int hashCode() {
            // for use here assume proteins are the same if accession and is_decoy are the same
            int result = 17;
            result = 31 * result + accession.hashCode();
            if (is_decoy) result = 31 * result + 1;
            return result;
        }

        @Override
        public boolean equals(Object obj) {
            // for use here assume proteins are the same if accession and is_decoy are the same
            if (this == obj) {
                return true;
            }
            if (obj == null) {
                return false;
            }
            if (getClass() != obj.getClass()) {
                return false;
            }
            final XiDBProtein other = (XiDBProtein) obj;
            if (this.is_decoy != other.is_decoy) {
                return false;
            }
            return this.accession.contentEquals(other.accession);
        }
    }
       
    private UUID resultset_id;
    private HashMap<UUID, Xi2Config> configs = new HashMap<>();
    private Xi2Score new_primary_score;
    private ArrayList<Xi2Score> scores;
    private ArrayList<XiDBPSM> dbPSMs = new ArrayList<>();
    private HashMap<UUID, HashSet<Integer>> search_to_pepID = new HashMap<>();
    private HashMap<UUID, HashMap<Integer, XiDBPeptide>> search_and_id_to_pep = new HashMap<>();
    private HashMap<UUID, HashSet<Integer>> search_to_protID = new HashMap<>();
    private HashMap<UUID, HashMap<Integer, XiDBProtein>> search_and_id_to_protein = new HashMap<>();
    Connection db_connection;
    String connection_string;
    String db_user;
    String db_pass;
    

    /***
     * tests if there is a working database connection aand if not opens up a connection
     * @return true: connected; false not connected
     * @throws SQLException 
     */
    public synchronized boolean ensureConnection() throws SQLException {
        boolean isconnected = false;
        if (db_connection != null) {
            try {
                Statement st = db_connection.createStatement();
                ResultSet rs = st.executeQuery("select 1");
                rs.close();
                st.close();
                isconnected = true;
            } catch (Exception sex) {
                Logger.getLogger(DBinFDR.class.getName()).log(Level.INFO, "Database connection is closed/ non-functioning. Will try to reopen", sex);
                try {
                    db_connection.close();
                } catch (Exception e) {
                }

            }
        }
//        }
//        }

        if (!isconnected) {
            Logger.getLogger(this.getClass().getName()).log(Level.WARNING, "Lost connection to database. Will try to reopen");

            try {
                if (connection_string != null && db_user != null && db_pass != null) {
                    this.db_connection = DriverManager.getConnection(this.connection_string, this.db_user, this.db_pass);
                } else {
                    Logger.getLogger(DBinFDR.class.getName()).log(Level.SEVERE, "Currently I don't have a way to reopen the connection");
                    return false;
                }

                isconnected = true;

                setupPreparedStatements();

            } catch (SQLException ex) {
                Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, 
                        "Could not reopen the connection to the database", ex);
                throw ex;
            }
        }
        return isconnected;
    }
    
    public void readconfig() throws SQLException{
        ensureConnection();
        Statement st = this.db_connection.createStatement(
                ResultSet.TYPE_FORWARD_ONLY, 
                ResultSet.CONCUR_READ_ONLY);
        ResultSet rs = st.executeQuery("SELECT s.id as id, config "
                + "FROM "
                + " Search s"
                + " INNER JOIN "
                + " ResultSearches rss ON rss.search_id = s.id "
                + "WHERE rss.resultset_id = '" + this.resultset_id + "'");
        while (rs.next()) {
            String config = rs.getString("config");
            UUID id = (UUID) rs.getObject("id");
            configs.put(id, new Xi2Config(config));
        }
    }
    
    public void read_score_names() throws SQLException{
        try {
            ensureConnection();
        } catch (SQLException sex) {
            Logger.getLogger(this.getClass().getName()).log(Level.SEVERE,"No connection to DB");
            return;
        }
        String querry = 
                "SELECT "
                + "  score_id, "
                + "  name, "
                + "  primary_score, "
                + "  higher_is_better, "
                + "  major_interest "
                + "FROM ScoreNames WHERE resultset_id = '" +
                this.resultset_id + "' ORDER BY score_id";
        Statement stm = this.db_connection.createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
        ResultSet rs = stm.executeQuery(querry);
        this.scores = new ArrayList<>();
        while (rs.next()) {
            Xi2Score s = new Xi2Score(
                    rs.getInt("score_id"), 
                    rs.getString("name"), 
                    rs.getBoolean("primary_score"), 
                    rs.getBoolean("higher_is_better"), 
                    rs.getBoolean("major_interest"));
            if (s.primary_score)
                this.new_primary_score = s;
            while (s.id >= this.scores.size()) {
                this.scores.add(null);
            }
            this.scores.set(s.id, s);
        }
    }
    

    public void readBasicMatchInfo() throws SQLException {
        
        // querry to retrive basic
        String querry = "SELECT  m.id, m.search_id, m.pep1_id, m.pep2_id, m.site1, m.site2,"
                + " m.link_score_site1, m.link_score_site2, m.link_score,"
                + " m.crosslinker_id, m.assumed_prec_mz, m.assumed_prec_charge,"
                + " m.calc_mass, m.top_ranking, m.is_decoy, m.is_dd, m.score,"
                + " rm.scores, rm.resultset_id "
                + " FROM resultmatch rm INNER JOIN match m ON rm.match_id = match.id "
                + " WHERE rm.resultset_id = " + this.resultset_id;
        
        // make sure we have a working connection to the DB
        ensureConnection();
        
        // make a non-exclusive querry
        Statement st  = this.db_connection.createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
        Logger.getLogger(this.getClass().getName()).log(Level.INFO,"Request CSM data");
        ResultSet rs = st.executeQuery(querry);
        
        Logger.getLogger(this.getClass().getName()).log(Level.INFO,"Start reading CSM data");
        // make sure we are using the right columns later on
        int cPsmID = rs.findColumn("id");
        int cSearchID = rs.findColumn("search_id");
        int cPep1ID = rs.findColumn("pep1_id");
        int cPep2ID = rs.findColumn("pep2_id");
        int cSite1 = rs.findColumn("site1");
        int cSite2 = rs.findColumn("site2");
        int cLinkScoreSite1 = rs.findColumn("link_score_site1");
        int cLinkScoreSite2 = rs.findColumn("link_score_site2");
        int cLinkScore = rs.findColumn("link_score");
        int cXLID = rs.findColumn("crosslinker_id");
        int cAssumedPrecMZ = rs.findColumn("assumed_prec_mz");
        int cAssumedPrecCharge = rs.findColumn("assumed_prec_charge");
        int cCalcMass = rs.findColumn("calc_mass");
        int cTopRanking = rs.findColumn("top_ranking");
        int cIsDecoy = rs.findColumn("is_decoy");
        int cIsDD = rs.findColumn("is_dd");
        int cScore = rs.findColumn("score");
        int cScores  = rs.findColumn("scores ");
        int cResultSetID = rs.findColumn("resultset_id");        
        
        while (rs.next()) {
            XiDBPSM psm = new XiDBPSM();
            // transfer the data into the temporar representation of a [PC]SM
            psm.id = rs.getObject(cPsmID, java.util.UUID.class);
            psm.search_id = rs.getObject(cSearchID, java.util.UUID.class);;
            psm.pep1_id = rs.getInt(cPep1ID);
            psm.pep2_id = rs.getInt(cPep2ID);
            psm.site1 = rs.getShort(cSite1);
            psm.site2 = rs.getShort(cSite2);
            psm.link_score_site1 = (short[])rs.getArray(cLinkScoreSite1).getArray();
            psm.link_score_site2 = (short[])rs.getArray(cLinkScoreSite2).getArray();
            psm.link_score = (float[])rs.getArray(cLinkScoreSite1).getArray();;
            psm.crosslinker_id = rs.getShort(cXLID);
            psm.assumed_prec_mz = rs.getDouble(cAssumedPrecMZ);
            psm.assumed_prec_charge = rs.getShort(cAssumedPrecCharge);
            psm.calc_mass = rs.getDouble(cCalcMass);
            psm.top_ranking = rs.getBoolean(cTopRanking);
            psm.is_decoy = rs.getBoolean(cIsDecoy);
            psm.is_dd  = rs.getBoolean(cIsDD);
            psm.all_scores = (float[]) rs.getArray(cScores).getArray();
            psm.matchscore = rs.getFloat(cScore);
            psm.result_score = psm.all_scores[this.new_primary_score.id];
            this.dbPSMs.add(psm);
            search_to_pepID.get(psm.search_id);
        }
        rs.close();
        st.close();
        
    }

    public void readPeptideInfo() throws SQLException {
        
        for (UUID sid : search_to_pepID.keySet()) {
            HashSet<Integer> pepids = search_to_pepID.get(sid);
            HashSet<Integer> protids = search_to_protID.get(sid);
            HashMap<Integer, XiDBPeptide> peps = search_and_id_to_pep.get(sid);
            String sPepIDs = RArrayUtils.toString(pepids, ",");
            
            String querry = "SELECT id, base_sequence, mass, length, "
                    + "modification_ids, modification_position, is_decoy FROM "
                    + "modifiedpeptide mp WHERE search_id = '" + sid + "' and "
                    + " peptide_id in (" + sPepIDs + ");";
            
            // make a non-exclusive querry
            Statement st  = this.db_connection.createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
            Logger.getLogger(this.getClass().getName()).log(Level.INFO,"Request CSM data");
            ResultSet rs = st.executeQuery(querry);

            Logger.getLogger(this.getClass().getName()).log(Level.INFO,"Start reading CSM data");
            // make sure we are using the right columns later on
            int cID = rs.findColumn("id");
            int cBaseSequence = rs.findColumn("base_sequence");
            int cMass = rs.findColumn("mass");
            int cLength = rs.findColumn("length");
            int cModificationIDs = rs.findColumn("modification_ids");
            int cModificationPosition = rs.findColumn("modification_position");
            int cIsDecoy = rs.findColumn("is_decoy");

            while (rs.next()) {
                XiDBPeptide pep = new XiDBPeptide();
                
                pep.base_sequence = rs.getString(cBaseSequence);
                pep.id = rs.getInt(cID);
                pep.mass = rs.getDouble(cMass);
                pep.length = rs.getInt(cLength);
                pep.modification_ids = (Integer[]) rs.getArray(cModificationIDs).getArray();
                pep.modification_position = (Integer[]) rs.getArray(cModificationIDs).getArray();
                pep.is_decoy = rs.getBoolean(cIsDecoy);
                
                peps.put(pep.id, pep);
            }
            rs.close();
            st.close();

            querry = "SELECT protein_id, mod_pep_id, start FROM "
                    + " modifiedpeptide mp WHERE search_id = '" + sid + "' and "
                    + " peptide_id in (" + sPepIDs + ");";
            
            // make a non-exclusive querry
            st  = this.db_connection.createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
            Logger.getLogger(this.getClass().getName()).log(Level.INFO,"Request CSM data");
            rs = st.executeQuery(querry);

            Logger.getLogger(this.getClass().getName()).log(Level.INFO,"Start reading CSM data");
            // make sure we are using the right columns later on
            int cPepID = rs.findColumn("mod_pep_id");
            int cProteinID= rs.findColumn("protein_id");
            int cStart = rs.findColumn("start");

            while (rs.next()) {
                Integer pepid = rs.getInt(cPepID);
                Integer protid = rs.getInt(cProteinID);
                Integer start = rs.getInt(cStart);
                XiDBPeptide pep = peps.get(pepid);
                pep.proteinIDs.add(protid);
                pep.proteinPositions.add(start);
                protids.add(protid);
            }
        }
    }

    public void readProteinInfo() throws SQLException {
        
        for (UUID sid : search_to_protID.keySet()) {
            HashSet<Integer> protids = search_to_protID.get(sid);
            HashMap<Integer, XiDBProtein> proteins = search_and_id_to_protein.get(sid);
            String sProtIDs = RArrayUtils.toString(protids, ",");
            
            String querry = "SELECT id, sequence, name, "
                    + "description, is_decoy FROM "
                    + " modifiedpeptide mp WHERE search_id = '" + sid + "' and "
                    + " peptide_id in (" + sProtIDs + ");";
            
            // make a non-exclusive querry
            Statement st  = this.db_connection.createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
            Logger.getLogger(this.getClass().getName()).log(Level.INFO,"Request CSM data");
            ResultSet rs = st.executeQuery(querry);

            Logger.getLogger(this.getClass().getName()).log(Level.INFO,"Start reading CSM data");
            // make sure we are using the right columns later on
            int cID = rs.findColumn("id");
            int cSequence = rs.findColumn("sequence");
            int cDescription = rs.findColumn("description");
            int cIsDecoy = rs.findColumn("is_decoy");

            while (rs.next()) {
                XiDBProtein prot = new XiDBProtein();
                
                prot.sequence = rs.getString(cSequence);
                prot.id = rs.getInt(cID);
                prot.description = rs.getString(cDescription);
                prot.is_decoy = rs.getBoolean(cIsDecoy);
                
                proteins.put(prot.id, prot);
            }
            rs.close();
            st.close();
        }
    }
    
    private void setup() throws SQLException {
        
        // querry to retrive basic
        String querry = "SELECT  search_id "
                + " WHERE rm.resultset_id = " + this.resultset_id;
        
        // make sure we have a working connection to the DB
        ensureConnection();
        
        // make a non-exclusive querry
        Statement st  = this.db_connection.createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_READ_ONLY);
        ResultSet rs = st.executeQuery(querry);
        
        while (rs.next()) {
            UUID search_id = rs.getObject(1, java.util.UUID.class);
            search_to_pepID.put(search_id, new HashSet<Integer>());
            search_to_protID.put(search_id, new HashSet<Integer>());
            search_and_id_to_pep.put(search_id, new HashMap<Integer, XiDBPeptide>());
            search_and_id_to_protein.put(search_id, new HashMap<Integer, XiDBProtein>());
        }
        rs.close();
        st.close();

        readconfig();
        read_score_names();
        
    }

    
    @Override
    public ArrayList<PSM> read(OfflineFDR target) {
        
        try {
            
            readBasicMatchInfo();
            readPeptideInfo();
            readProteinInfo();
            
            for (XiDBPSM dbpsm : dbPSMs) {
                
                /*target.addMatch(
                        dbpsm.id.toString(), 
                        dbpsm.pep1_id, 
                        dbpsm.pep2_id, String pepSeq1, String pepSeq2, int peplen1, int peplen2, int site1, int site2, boolean isDecoy1, boolean isDecoy2, int charge, double score, Long protid1, String accession1, String description1, Long protid2, String accession2, String description2, int pepPosition1, int pepPosition2, String Protein1Sequence, String Protein2Sequence, double peptide1score, double peptide2score, String isSpecialCase, String crosslinker, String run, String Scan);
                */
                
            }
            
            return null;
        } catch (SQLException ex) {
            Logger.getLogger(DBImportXi2.class.getName()).log(Level.SEVERE, null, ex);
        }
        return null;
    }

    private void setupPreparedStatements() {
        // none defined at the moment
    }
    
    /**
     * 
     * @param args 
     */
    public static void main(String[] args) throws SQLException {
        DBImportXi2 imp = new DBImportXi2();
        imp.connection_string = "jdbc:postgresql://192.168.121.72:5432/xisearch2";
        imp.db_user = "xisearch2";
        imp.db_pass = "xisearch2";
        imp.resultset_id = UUID.fromString("801a5c73-caa0-4c90-acb7-1608c4ff80bf");
        imp.readconfig();
        imp.setup();
        imp.read_score_names();
        imp.readBasicMatchInfo();
        imp.readPeptideInfo();
        imp.readProteinInfo();
        OfflineFDR fdr = new OfflineFDR() {
            @Override
            public String getSource() {
                throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
            }
        };
        imp.read(fdr);
                
    }
    
}
