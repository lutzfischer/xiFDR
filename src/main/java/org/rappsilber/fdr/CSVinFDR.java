/*
 * Copyright 2015 Lutz Fischer <lfischer at staffmail.ed.ac.uk>.
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
package org.rappsilber.fdr;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.sql.SQLException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.rappsilber.data.csv.CsvParser;
import org.rappsilber.data.csv.condition.CsvCondition;
import org.rappsilber.fdr.entities.PSM;
import org.rappsilber.utils.AutoIncrementValueMap;
import org.rappsilber.utils.RArrayUtils;
import static org.rappsilber.utils.StringUtils.sixDigits;
import org.rappsilber.utils.UpdatableChar;

/**
 *
 * @author lfischer
 */
public class CSVinFDR extends OfflineFDR {
    private ArrayList<String> m_source = new ArrayList<>();
    private ArrayList<CsvCondition> m_filter = new ArrayList<>();
    private ArrayList<String[]> commandLineColumnMapping;
    private Character delimiter;
    private Character quote;
    public static String[][] DEFAULT_COLUMN_MAPPING=new String[][]{
        new String[]{"peptide1","Peptide1"},
        new String[]{"peptide2","Peptide2"},
        new String[]{"peptide length 1","LengthPeptide1"},
        new String[]{"peptide length 2","LengthPeptide2"},
        new String[]{"peptide link 1","Link1"},
        new String[]{"peptide link 2","Link2"},
        new String[]{"is decoy 1","Protein1decoy"},
        new String[]{"is decoy 2","Protein2decoy"},
        new String[]{"precursor charge","PrecoursorCharge"},
        new String[]{"score","match score"},
        new String[]{"accession1","Protein1"},
        new String[]{"accession2","Protein2"},
        new String[]{"description1","Fasta1"},
        new String[]{"description2","Fasta2"},
        new String[]{"peptide position 1","Start1"},
        new String[]{"peptide position 2","Start2"},
        new String[]{"peptide1 score","Pep1Score"},
        new String[]{"peptide2 score","Pep2Score"},
        new String[]{"experimental mz","expMZ"},
        new String[]{"calculated mass","calcMass"},
        new String[]{"negative grouping","NegativeGrouping"},
        new String[]{"positive grouping","PositiveGrouping"},
        new String[]{"info","info"},
    };
    
    public CSVinFDR() {
    }

    public CSVinFDR(int[] peptideLengthGroups) {
        super(peptideLengthGroups);
    }
    

    
    @Override
    protected ArrayList<String> getPSMHeader() {
        ArrayList<String> ret = super.getPSMHeader();
//        String header = "PSMID" + seperator + "run" + seperator + "scan" 
//                + seperator + "exp charge"  + seperator + "exp m/z" + seperator + "exp mass" + seperator + "exp fractionalmass" 
//                + seperator + "match charge" + seperator + "match mass" + seperator + "match fractionalmass"   
//                
//                + seperator + "Protein1" + seperator + "Description1" + seperator + "Decoy1" + seperator + "Protein2" + seperator + "Description2" + seperator + "Decoy2" + seperator + "Peptide1" + seperator + "Peptide2" + seperator + "PeptidePosition1" + seperator + "PeptidePosition2" + seperator + "PeptideLength1" + seperator + "PeptideLength2" + seperator + "fromSite" + seperator + "ToSite" + seperator + "Charge" + seperator + "Score" + seperator + "isDecoy" + seperator + "isTT" + seperator + "isTD" + seperator + "isDD" + seperator + "fdrGroup" + seperator + "fdr" + seperator + "" + seperator + "PeptidePairFDR" + seperator + "Protein1FDR" + seperator + "Protein2FDR" + seperator + "LinkFDR" + seperator + "PPIFDR" + seperator + seperator + "peptide pair id" + seperator + "link id" + seperator + "ppi id";
        ret.add(2, "exp charge");
        ret.add(3, "exp m/z");
        ret.add(4, "exp mass");
        ret.add(5, "exp fractionalmass");
        ret.add(6, "match charge");
        ret.add(7,  "match mass");
        ret.add(8, "match fractionalmass");
        return ret;
    }

    
    protected ArrayList<String> getPSMOutputLine(PSM pp) {
        ArrayList<String> ret = super.getPSMOutputLine(pp);

        double mz = pp.getExperimentalMZ();
        int charge = pp.getExpCharge();
        double mass = (mz-1.00727646677)*charge;
        double fraction = mass-Math.floor(mass);
        double calcfraction = pp.getCalcMass()-Math.floor(pp.getCalcMass());

        
        ret.add(2,""+ charge);
        ret.add(3, sixDigits.format(mz));
        ret.add(4, sixDigits.format(mass));
        ret.add(5, sixDigits.format(fraction));
        ret.add(6,""+ pp.getCharge());
        ret.add(7,  sixDigits.format(pp.getCalcMass()));
        ret.add(8, sixDigits.format(calcfraction));
        

        return ret;
    }    
    
    protected Integer getColumn(CsvParser csv,String column,boolean optional) throws ParseException {
        Integer c = csv.getColumn(column);
        if (c== null && !optional) {
            throw new ParseException("Column " + column + " not found", 0);
        }
        return c;
    }
    
    public void readCSV(File f) throws FileNotFoundException, IOException, ParseException  {
        readCSV(CsvParser.guessCsv(f, 50), null);
    }

    public void readCSV(CsvParser csv, CsvCondition filter) throws FileNotFoundException, IOException, ParseException  {
        if (csv.getInputFile()!= null) {
            m_source.add(csv.getInputFile().getAbsolutePath());
        } else {
            m_source.add(null);
        }
        m_filter.add(filter);
            

        Integer crun = getColumn(csv,"run",true);
        Integer cscan = getColumn(csv,"scan",true);
        Integer cpsmID = getColumn(csv,"psmid",true);
        int cpep1 = getColumn(csv,"peptide1",false);
        int cpep2 = getColumn(csv,"peptide2",false);
        Integer cpep1len = getColumn(csv,"peptide length 1",true);
        Integer cpep2len = getColumn(csv,"peptide length 2",true);
        int cpep1site = getColumn(csv,"peptide link 1",false);
        int cpep2site = getColumn(csv,"peptide link 2",false);
        int cpep1decoy = getColumn(csv,"is decoy 1",false);
        Integer cpep2decoy = getColumn(csv,"is decoy 2",true);
        int cprecZ = getColumn(csv,"precursor charge",false);
        int cscore = getColumn(csv,"score",false);
        int caccession1 = getColumn(csv,"accession1",false);
        int caccession2 = getColumn(csv,"accession2",false);
        Integer cdescription1 = getColumn(csv,"description1",true);
        Integer cdescription2 = getColumn(csv,"description2",true);
        Integer cpeptide_position1 = getColumn(csv,"peptide position 1",true);
        Integer cpeptide_position2 = getColumn(csv,"peptide position 2",true);
        Integer cscoreratio = getColumn(csv,"score ratio",true);
        Integer cPepScore1 = getColumn(csv,"peptide1 score",true);
        Integer cPepScore2 = getColumn(csv,"peptide2 score",true);
        Integer cCrosslinker = getColumn(csv,"crosslinker",true);
        Integer cExpMZ = getColumn(csv,"experimental mz",true);
        Integer cCalcMass = getColumn(csv,"calculated mass",true);
        Integer cInfo = getColumn(csv,"info",true);
        Integer cNegativeGrouping = getColumn(csv,"negative grouping",true);
        Integer cPositiveGrouping = getColumn(csv,"positive grouping",true);
        
        int noID = 0;
        AutoIncrementValueMap<String> PSMIDs = new AutoIncrementValueMap<String>();
        int lineNumber =0;
        String terminalAminoAcids = "(?:[A-Z\\-][a-z]*|\\[(?:[A-Z\\-][a-z]*)+\\])";
        Pattern hasTerminalAminoAcids = Pattern.compile("^" +terminalAminoAcids+"\\.(.*)\\."+terminalAminoAcids + "$");
        Pattern nterminalAminoAcids = Pattern.compile("^" + terminalAminoAcids +"\\.");
        Pattern cterminalAminoAcids = Pattern.compile("\\."+terminalAminoAcids +"$");
        
        int noPSMID =0;
        while (csv.next()) {
            lineNumber++;
            if (filter != null && !filter.fits(csv))
                continue;
            String psmID;

            String pepSeq1 = csv.getValue(cpep1);
            String pepSeq2 = csv.getValue(cpep2);

            // if the sequence looks like K.PEPTIDEK.A assume that the first and 
            // last aminoacid are leading and trailing aminoacids and not really part of the peptide
            Matcher m1=hasTerminalAminoAcids.matcher(pepSeq1);
            Matcher m2=hasTerminalAminoAcids.matcher(pepSeq2);
            if ((pepSeq2 == null || pepSeq2.isEmpty()) && m1.matches()) {
                pepSeq1 = m1.replaceAll("$1");
            }
                
            if (m1.matches() && m2.matches()) {
                pepSeq1 = m1.replaceAll("$1");
                pepSeq2 = m2.replaceAll("$1");
            }

            
            Integer site1 = csv.getInteger(cpep1site,-1);
            Integer site2 = csv.getInteger(cpep2site,-1); //pepSeq2 == null || pepSeq2.trim().isEmpty() ? -1 : csv.getInteger(cpep2site,-1);
            
            // do we have to generate an ID?
            if (cpsmID == null) {
                if (cscan == null || crun == null) {
                    psmID = Integer.toString(noPSMID++);
                } else {
                    String key = "Scan: " + csv.getValue(cscan) + " Run: " + csv.getValue(crun);
                    int c= pepSeq1.compareTo(pepSeq2) ;
                    if (c > 0 || (c==0 && site1 > site2) ) {
                        key=key +" P1_" + csv.getValue(cpep1) + " P2_" + csv.getValue(cpep2) + " " + csv.getInteger(cpep1site) + " " + csv.getInteger(cpep2site);
                    } else {
                        key=key +" P1_" + csv.getValue(cpep2) + " P2_" + csv.getValue(cpep1) + " " + csv.getInteger(cpep2site) + " " + csv.getInteger(cpep1site);;
                    }
                    //psmID = PSMIDs.toIntValue(key);
                    psmID = key;
                }
            }else
                //psmID=csv.getInteger(cpsmID);
                psmID=csv.getValue(cpsmID);
            
            
            // if we have a column for the peptide length take that value
            // otherwise count all capital letters in the sequence and define 
            // this as length 
            int peplen1 = cpep1len == null ? pepSeq1.replaceAll("[^A-Z]", "").length() : csv.getInteger(cpep1len);

            Integer peplen2 = null;
            if (cpep2len == null) 
                if (pepSeq2 == null) {
                    peplen2 = 0; 
                } else {
                    peplen2 = pepSeq2.replaceAll("[^A-Z]", "").length();
                }
            else {
                peplen2 = csv.getInteger(cpep2len, 0);
            }
            
            boolean isDecoy1 = csv.getBool(cpep1decoy,false);
            boolean isDecoy2=  cpep2decoy == null ? false : csv.getBool(cpep2decoy, false);
            int charge = csv.getInteger(cprecZ);
            Double score = csv.getDouble(cscore);
            String saccession1 = csv.getValue(caccession1);
            String sdescription1 = csv.getValue(cdescription1);
            String saccession2 = csv.getValue(caccession2);
            String sdescription2 = csv.getValue(cdescription2);
            String spepPosition1 = csv.getValue(cpeptide_position1);
            String spepPosition2= csv.getValue(cpeptide_position2);
            if (spepPosition2 == null || spepPosition2.trim().isEmpty()) 
                spepPosition2 = "-1";
            
            
            // how to split up the score
            double scoreRatio = csv.getDouble(cscoreratio);
            if (Double.isNaN(scoreRatio)) {
                if (cPepScore1 != null && cPepScore2 != null) {
                    Double s1 = csv.getDouble(cPepScore1);
                    Double s2 = csv.getDouble(cPepScore2, 0.0);
                    
                    if (!Double.isNaN(s2))
                        scoreRatio = s1/(s1+s2);
                    else 
                        scoreRatio = 1;
                    
                    // we don't have a sumed up score but two peptide scores
                    if (Double.isNaN(score)) {
                        if (Double.isNaN(s2))
                            score = s1;
                        else
                            score = s1+s2;
                    }
                } else {
                    scoreRatio = (4.0/5.0+(peplen1/(peplen1+peplen2)))/2;
                }
            }
            
            String[] accessions1 =saccession1.split(";");
            String[] accessions2 =saccession2.split(";");
            String[] descriptions1 =sdescription1.split(";");
            String[] descriptions2 =sdescription2.split(";");
            String[] pepPositions1 =spepPosition1.split(";");
            String[] pepPositions2 =spepPosition2.split(";");
            
            if (!sdescription1.isEmpty()) {
                if (descriptions1.length != accessions1.length)
                    throw new ParseException("Don't know how to handle different numbers of protein accessions and descriptions", lineNumber);
            } else {
                descriptions1 = accessions1;
            }

            if (!sdescription2.isEmpty()) {
                if (descriptions2.length != accessions2.length)
                    throw new ParseException("Don't know how to handle different numbers of protein accessions and descriptions", lineNumber);
            } else {
                descriptions2 = accessions2;
            }
            
            if (accessions1.length ==1) {
                if (pepPositions1.length > 1) {
                    String[] daccessions1 = new String[pepPositions1.length];
                    String[] ddescriptions1 = new String[pepPositions1.length];
                    for (int i = 0; i< pepPositions1.length; i++) {
                        daccessions1[i] = saccession1;
                        ddescriptions1[i] = descriptions1[0];
                    }
                    accessions1 = daccessions1;
                    descriptions1 = ddescriptions1;
                }
            } else if (accessions1.length != pepPositions1.length){
                throw new ParseException("Don't know how to handle different numbers of proteins and peptide positions", lineNumber);
            }

            if (accessions2.length ==1) {
                if (pepPositions2.length > 1) {
                    String[] daccessions2 = new String[pepPositions2.length];
                    String[] ddescriptions2 = new String[pepPositions2.length];
                    for (int i = 0; i< pepPositions1.length; i++) {
                        daccessions2[i] = saccession2;
                        ddescriptions2[i] = descriptions2[0];
                    }
                    accessions2 = daccessions2;
                    descriptions2 = ddescriptions2;
                }
            } else if (accessions2.length != pepPositions2.length){
                throw new ParseException("Don't know how to handle different numbers of proteins and peptide positions", lineNumber);
            }
            
            int[] ipeppos1 = new int[pepPositions1.length];
            for (int i = 0; i<pepPositions1.length; i++) {
                ipeppos1[i] = Integer.parseInt(pepPositions1[i]);
            }
            
            int[] ipeppos2 = new int[pepPositions2.length];
            for (int i = 0; i<pepPositions2.length; i++) {
                ipeppos2[i] = Integer.parseInt(pepPositions2[i]);
            }
            
            
            PSM psm = null;
            for (int p1 = 0; p1< accessions1.length; p1++) {
                for (int p2 = 0; p2< accessions2.length; p2++) {
                    String run = crun == null ? "":csv.getValue(crun);
                    String scan = cscan == null ? "":csv.getValue(cscan);
                    String crosslinker = cCrosslinker == null ? "":csv.getValue(cCrosslinker);
                    String negativeCase = null;
                    String poisitiveCase = null;
                    if (cNegativeGrouping >=0) {
                        negativeCase = csv.getValue(cNegativeGrouping);
                        if (negativeCase.isEmpty())
                            negativeCase=null;
                    }
                    if (cPositiveGrouping >=0) {
                        poisitiveCase = csv.getValue(cPositiveGrouping);
                        if (poisitiveCase.isEmpty())
                            poisitiveCase=null;
                    }
                    
                    psm = addMatch(psmID, pepSeq1, pepSeq2, peplen1, peplen2, 
                            site1, site2, isDecoy1, isDecoy2, charge, score, 
                            accessions1[p1], descriptions1[p1], accessions2[p2],
                            descriptions2[p2], ipeppos1[p1], ipeppos2[p2], 
                            scoreRatio, negativeCase,crosslinker,run,scan);
                    if (poisitiveCase != null) {
                        psm.setPositiveGrouping(poisitiveCase);
                    }
                    Double expMZ = null;
                    if (cExpMZ != null) {
                        expMZ=csv.getDouble(cExpMZ);
                        psm.setExperimentalMZ(expMZ);
                    }

                    // read exp- and calc-mass
                    if (cCalcMass != null) {
                        Double calcMass=csv.getDouble(cCalcMass);
                        
                        psm.setCalcMass(calcMass);
                        
                        if (expMZ != null) {
                            psm.setExpCharge((byte)Math.round(calcMass/expMZ));
                            
                        }
                    }
                    
                    
                    

                    if (cInfo != null) {
                        psm.setInfo(csv.getValue(cInfo));
                    }
//    public PSM          addMatch(String psmID, Integer pepid1, Integer pepid2, String pepSeq1, String pepSeq2, int peplen1, int peplen2, int site1, int site2, boolean isDecoy1, boolean isDecoy2, int charge, double score, Integer protid1, String accession1, String description1, Integer protid2, String accession2, String description2, int pepPosition1, int pepPosition2, String Protein1Sequence, String Protein2Sequence, double scoreRatio, boolean isSpecialCase, String crosslinker, String run, String Scan) {

                }
            }
            
        }
    }

    @Override
    public String getSource() {
        boolean filterFound = false;
        boolean filterSame = true;
        String filter ="";
        String ret = "";
        if (m_filter.get(0) != null) {
            filter = m_filter.get(0).toString();
        }
        // do we have any filter and are they all the same?
        for (CsvCondition c : m_filter) {
            if (c != null) {
                filterFound = true;
                if (!filter.equals(c.toString())) {
                    filterSame = false;
                    break;
                }
            }
        }
        
        if (filterFound) {
            if (filterSame) {
                StringBuilder sb = new StringBuilder(RArrayUtils.toString(m_source, "\n"));
                sb.append("\nfilter:").append(filter);
                ret = sb.toString();
            } else {
                StringBuilder sb = new StringBuilder();
                for (int s = 0; s<m_source.size(); s++) {
                    sb.append(m_source.get(s)).append("\n");
                    sb.append("filter:");
                    if (m_filter.get(s) == null)
                        sb.append("no filter\n");
                    else 
                        sb.append(m_filter.get(s)).append("\n");
                }
                sb.delete(sb.length(), sb.length()+1);
                ret = sb.toString();
            }
        } else {
            if (m_source != null)
                ret = RArrayUtils.toString(m_source, "\n") ;
        }
        return ret;
    }

    public String getSources() {
        if (m_source != null)
            return RArrayUtils.toString(m_source, "\n") ;
        return "";
    }

    public void setSource(String source) {
        m_source=new ArrayList<>();
        m_source.add(source);
    }

    public void addSource(String source) {
        m_source.add(source);
        m_filter.add(null);
    }

    public void addSource(String source,CsvCondition filter) {
        m_source.add(source);
        m_filter.add(filter);
    }
    
    public void addSource(ArrayList<String> source,ArrayList<CsvCondition> filter) {
        m_source.addAll(source);
        m_filter.addAll(filter);
    }
    
    public String execClass() {
        return this.getClass().getName();
    }
    
    public String argList() {
        return super.argList() + " --map=col:name,col:name --delimiter= --quote= csv-file1 csv-file2";
    }
    
    public String argDescription() {
        return super.argDescription() + "\n"
                + "--map=X                  a mapping of column name to expected\n"
                + "                         column names\n"
                + "                         the column are:\n"
                + "                             run (optional)\n"
                + "                             scan (optional)\n"
                + "                             psmid (optional)\n"
                + "                             peptide1\n"
                + "                             peptide2\n"
                + "                             peptide length 1 (optional)\n"
                + "                             peptide length 2 (optional)\n"
                + "                             peptide link 1\n"
                + "                             peptide link 2\n"
                + "                             is decoy 1\n"
                + "                             is decoy 2\n"
                + "                             precursor charge\n"
                + "                             score\n"
                + "                             accession1\n"
                + "                             accession2\n"
                + "                             description1 (optional)\n"
                + "                             description2 (optional)\n"
                + "                             peptide position 1\n"
                + "                             peptide position 2\n"
                + "                             score ratio (optional)\n"
                + "                             peptide1 score (optional)\n"
                + "                             peptide2 score (optional)\n"
                + "                             crosslinker (optional)\n"
                + "                             calculated mass (optional)\n"
                + "                             experimental mz (optional)\n"
                + "                             info (optional)\n"
                + "                         either run and scan numebr must\n"
                + "                         or psmid needs to be specified\n "
                + "                         Example:\n "
                + "                         --map=scan:spectrum,psmid:id\n "
                + "                          maps scan to spectrum and\n"
                + "                          psmid to id as the columns\n"
                + "                          in the CSV\n"
                + "--delimiter              what separates fields in the file\n "
                + "--quote                  how are text fields qoted\n"
                + "                         e.g. each field that contains the\n"
                + "                         delimiter needs to be in quotes\n ";
        
    }
    
        
    @Override
    public String[] parseArgs(String[] argv) {
        ArrayList<String> unknown = new ArrayList<String>();
        
        argv = super.parseArgs(argv);
        
        for (String arg : argv) {
            if (arg.toLowerCase().startsWith("--map=")) {
                String mappings=arg.substring(6);
                String[] mappairs = mappings.split(",");

                if (commandLineColumnMapping == null)
                    commandLineColumnMapping = new ArrayList<>(mappairs.length);
                
                for (String mp : mappairs){ 
                    commandLineColumnMapping.add(mp.split(":"));
                }
            } else if(arg.toLowerCase().startsWith("--quote=")) {
                String quotechar=arg.substring("--quote=".length());
                Character q;
                if (quotechar.contentEquals("\\t")) {
                    quotechar="\t";
                } else if (quotechar.contentEquals("\\s")) {
                    quotechar=" ";
                } else if (quotechar.length() >1) {
                    Logger.getLogger(CSVinFDR.class.getName()).log(Level.SEVERE, "currently only single charachter quotes are supported");
                    System.exit(-1);
                }
                quote = quotechar.charAt(0);
            } else if(arg.toLowerCase().startsWith("--delimiter=")) {
                String delchar=arg.substring("--delimiter=".length());
                Character d;
                if (delchar.contentEquals("\\t")) {
                    delchar="\t";
                } else if (delchar.contentEquals("\\s")) {
                    delchar=" ";
                } else if (delchar.length() >1) {
                    Logger.getLogger(CSVinFDR.class.getName()).log(Level.SEVERE, "currently only single charachter delimiter are supported");
                    System.exit(-1);
                }
                delimiter = delchar.charAt(0);
            } else if(arg.toLowerCase().contentEquals("--help")) {
                printUsage();
                System.exit(0);
            }  else {
               unknown.add(arg);
            }
            
        }        
        String[] ret = new String[unknown.size()];
        ret = unknown.toArray(ret);
        return ret;        
    }
    
    public static void main (String[] argv) throws SQLException, FileNotFoundException {
        
        CSVinFDR ofdr = new CSVinFDR();
        
        String[] files = ofdr.parseArgs(argv);
        
        // assume that everything that was not matched to an argument is a file
        
        if (files.length == 0) {
            ofdr.printUsage();
            System.exit(1);
        }
        
        
        if (ofdr.getCsvOutDirSetting() != null) {
            if (ofdr.getCsvOutBaseSetting() == null) {
                ofdr.setCsvOutBaseSetting("FDR");
            }
        }
        
        if (ofdr.getCsvOutBaseSetting() != null) {
            if (ofdr.getCsvOutDirSetting() == null)
                ofdr.setCsvOutDirSetting(".");
            
            System.out.println("writing results to " + ofdr.getCsvOutDirSetting() + "/" + ofdr.getCsvOutBaseSetting() + "*");
            System.out.flush();
        }

        CsvParser csv = new CsvParser();
        if (ofdr.commandLineColumnMapping != null) {
            for (String[] map : ofdr.commandLineColumnMapping) {
                for (int i = 1; i<map.length;i++) {
                    csv.setAlternative(map[0], map[i]);
                }
            }
        } else {
            for (String[] map : ofdr.DEFAULT_COLUMN_MAPPING) {
                for (int i = 1; i<map.length;i++) {
                    csv.setAlternative(map[0], map[i]);
                }
            }
        }
        
        
        // read in all files
        for (String f : files) {
            Logger.getLogger(CSVinFDR.class.getName()).log(Level.INFO, "setting up csv input");
            
            UpdatableChar delimChar = new UpdatableChar(',');
            UpdatableChar quoteChar = new UpdatableChar('"');
            if (ofdr.delimiter != null) {
                delimChar.value = ofdr.delimiter;
            }
            if (ofdr.quote != null) {
                quoteChar.value = ofdr.quote;
            }
            if (ofdr.delimiter == null || ofdr.quote == null) {
                try {
                    csv.guessDelimQuote(new File(f), 50, delimChar, quoteChar);
                } catch (IOException ex) {
                    Logger.getLogger(CSVinFDR.class.getName()).log(Level.SEVERE, "error while quessing csv-definitions", ex);
                }
            }
            try {
                csv.openFile(new File(f), true);
            } catch (IOException ex) {
                Logger.getLogger(CSVinFDR.class.getName()).log(Level.SEVERE, "Could not read the file", ex);
                System.exit(-1);
            }
            
            Logger.getLogger(CSVinFDR.class.getName()).log(Level.INFO, "Read datafrom CSV");
            try {
                ofdr.readCSV(csv,null);
            } catch (IOException ex) {
                Logger.getLogger(CSVinFDR.class.getName()).log(Level.SEVERE, "Error while reading file: " + f, ex);
                System.exit(-1);
            } catch (ParseException ex) {
                Logger.getLogger(CSVinFDR.class.getName()).log(Level.SEVERE, "Error parsing file: " + f, ex);
                System.exit(-1);
            }
        }
        
        Logger.getLogger(CSVinFDR.class.getName()).log(Level.INFO, "Calculate FDR");
        ofdr.calculateWriteFDR(ofdr.getCsvOutDirSetting(), ofdr.getCsvOutBaseSetting(), ",");



        System.exit(0);

        
    }

    /**
     * @return the commandLineColumnMapping
     */
    public ArrayList<String[]> getCommandLineColumnMapping() {
        return commandLineColumnMapping;
    }

    /**
     * @param commandLineColumnMapping the commandLineColumnMapping to set
     */
    public void setCommandLineColumnMapping(ArrayList<String[]> commandLineColumnMapping) {
        this.commandLineColumnMapping = commandLineColumnMapping;
    }

    /**
     * @return the delimiter
     */
    public Character getDelimiter() {
        return delimiter;
    }

    /**
     * @param delimiter the delimiter to set
     */
    public void setDelimiter(Character delimiter) {
        this.delimiter = delimiter;
    }

    /**
     * @return the quote
     */
    public Character getQuote() {
        return quote;
    }

    /**
     * @param quote the quote to set
     */
    public void setQuote(Character quote) {
        this.quote = quote;
    }

    public CsvCondition getFilter() {
        return m_filter.get(0);
    }

    public Collection<CsvCondition> getFilters() {
        return m_filter;
    }
    
}
