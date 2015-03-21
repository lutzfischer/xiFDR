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
import java.text.ParseException;
import org.rappsilber.data.csv.CsvParser;
import org.rappsilber.fdr.entities.PSM;
import org.rappsilber.utils.AutoIncrementValueMap;
import org.rappsilber.data.csv.ColumnAlternatives;

/**
 *
 * @author lfischer
 */
public class CSVinFDR extends OfflineFDR {
    private String m_source = null;


    public CSVinFDR() {
    }

    public CSVinFDR(int[] peptideLengthGroups) {
        super(peptideLengthGroups);
    }
    
    public void setColumnName(String standartName, String alternative) {
        
    }
    
    

    public void readCSV(File f) throws FileNotFoundException, IOException, ParseException  {
        readCSV(CsvParser.guessCsv(f, 50));
    }

    public void readCSV(CsvParser csv) throws FileNotFoundException, IOException, ParseException  {
        if (csv.getInputFile()!= null)
            m_source = csv.getInputFile().getAbsolutePath();
            
        ColumnAlternatives.setupAlternatives(csv);
        csv.setAlternative("psmid", "id");
        //csv.setAlternative("psmid", "id");
        
//                  sm.id AS psmID, "
//                + "p1.sequence AS pepSeq1, "
//                + "p2.sequence AS pepSeq2, "
//                + "p1.peptide_length as peplen1, "
//                + "p2.peptide_length as peplen2, "
//                + "mp1.link_position as site1, "
//                + "mp2.link_position as site2, "
//                + "pr1.is_decoy AS isDecoy1, "
//                + "pr2.is_decoy AS isDecoy2, "
//                + "sm.precursor_charge AS charge, "
//                + "sm.score, "
//                + "pr1.accession_number AS accession1, "
//                + "pr1.description AS  description1, "
//                + "pr2.accession_number AS accession2, "
//                + "pr2.description AS  description2, "
//                + "hp1.peptide_position AS pepPosition1,  "
//                + "hp2.peptide_position AS pepPosition2, "
//                "CASE WHEN p2.sequence IS NULL THEN 1 ELSE 2.0/3.0 END AS score_ratio "
        Integer crun = csv.getColumn("run");
        Integer cscan = csv.getColumn("scan");
        Integer cpsmID = csv.getColumn("psmid");
        int cpep1 = csv.getColumn("peptide1");
        int cpep2 = csv.getColumn("peptide2");
        Integer cpep1len = csv.getColumn("peptide length 1");
        Integer cpep2len = csv.getColumn("peptide length 2");
        int cpep1site = csv.getColumn("peptide link 1");
        int cpep2site = csv.getColumn("peptide link 2");
        int cpep1decoy = csv.getColumn("is decoy 1");
        Integer cpep2decoy = csv.getColumn("is decoy 2");
        int cprecZ = csv.getColumn("precursor charge");
        int cscore = csv.getColumn("score");
        int caccession1 = csv.getColumn("accession1");
        int caccession2 = csv.getColumn("accession2");
        Integer cdescription1 = csv.getColumn("description1");
        Integer cdescription2 = csv.getColumn("description2");
        Integer cpeptide_position1 = csv.getColumn("peptide_position_1");
        Integer cpeptide_position2 = csv.getColumn("peptide_position_2");
        Integer cscoreratio = csv.getColumn("score ratio");
        Integer cPepScore1 = csv.getColumn("peptide1 score");
        Integer cPepScore2 = csv.getColumn("peptide2 score");
        
        int noID = 0;
        AutoIncrementValueMap<String> PSMIDs = new AutoIncrementValueMap<String>();
        int lineNumber =0;
        
        
        while (csv.next()) {
            lineNumber++;
            String psmID;
            if (cpsmID == null) {
                String key = "Scan: " + csv.getValue(cscan) + " Run: " + csv.getValue(crun);
                //psmID = PSMIDs.toIntValue(key);
                psmID = key;
            }else
                //psmID=csv.getInteger(cpsmID);
                psmID=csv.getValue(cpsmID);
            
  
            String pepSeq1 = csv.getValue(cpep1);
            String pepSeq2 = csv.getValue(cpep2);
            
            // if the sequence looks like K.PEPTIDEK.A assume that the first and 
            // last aminoacid are leading and trailing aminoacids and not really part of the peptide
            if (pepSeq1.matches("[A-Z\\-]+\\..*\\.[A-Z\\-]+")) {
                pepSeq1 = pepSeq1.replaceAll("^[A-Z\\-]+\\.", "").replaceAll("\\.[A-Z\\-]+", "");
            }
            if (pepSeq2.matches("[A-Z\\-]+\\..*\\.[A-Z\\-]+")) {
                pepSeq2 = pepSeq2.replaceAll("^[A-Z\\-]+\\.", "").replaceAll("\\.[A-Z\\-]+", "");
            }
            
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
            
            Integer site1 = csv.getInteger(cpep1site,-1);
            Integer site2 = csv.getInteger(cpep2site,-1); //pepSeq2 == null || pepSeq2.trim().isEmpty() ? -1 : csv.getInteger(cpep2site,-1);
            boolean isDecoy1 = csv.getBool(cpep1decoy);
            boolean isDecoy2=  cpep2decoy == null ? false : csv.getBool(cpep2decoy, false);
            int charge = csv.getInteger(cprecZ);
            Double score = csv.getDouble(cscore);
            String saccession1 = csv.getValue(caccession1);
            String sdescription1 = csv.getValue(cdescription1);
            String saccession2 = csv.getValue(caccession2);
            String sdescription2 = csv.getValue(caccession2);
            String spepPosition1 = csv.getValue(cpeptide_position1);
            String spepPosition2= csv.getValue(cpeptide_position2);
            if (spepPosition2 == null || spepPosition2.trim().isEmpty()) 
                spepPosition2 = "-1";
            
//            double coverage1  = csv.getDouble(18);
//            double coverage2  = csv.getDouble(19);
//            double scoreRatio = coverage1/(coverage1+coverage2);
            
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
            
            for (int p1 = 0; p1< accessions1.length; p1++) {
                for (int p2 = 0; p2< accessions2.length; p2++) {
                    PSM psm = addMatch(psmID, pepSeq1, pepSeq2, peplen1, peplen2, site1, site2, isDecoy1, isDecoy2, charge, score, accessions1[p1], descriptions1[p1], accessions2[p2], descriptions2[p2], ipeppos1[p1], ipeppos2[p2], scoreRatio, false);
                    if (cscan !=null)
                        psm.setScan(csv.getValue(cscan));
                    if (crun !=null)
                        psm.setRun(csv.getValue(crun));
                    
                }
            }
            
        }
    }

    @Override
    public String getSource() {
        if (m_source != null)
            return m_source;
        return "";
    }
}
