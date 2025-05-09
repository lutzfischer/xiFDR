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
package org.rappsilber.fdr.entities;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Objects;
import org.rappsilber.utils.IntArrayList;
import org.rappsilber.utils.RArrayUtils;

/**
 * Represent a single peptide
 * @author lfischer
 */
public class Peptide extends AbstractFDRElement<Peptide>  { //implements Comparable<Peptide>{
    
    /**
     * A dummy peptide used as the "second" peptide for linear matches
     */
    public static Peptide NOPEPTIDE = new Peptide(Integer.MAX_VALUE,"", false, Protein.NOPROTEIN, 0, 0);
    
//    /** the score of the peptide */
//    private double score = Double.NEGATIVE_INFINITY;
    /** Amino-acid sequence of the peptide */
    private String sequence;
    /** is this a decoy peptide*/
    boolean isDecoy;
    /** length of the peptide */
    int length =0;
    /** peptide pairs that this peptide is part of */
    HashSet<PeptidePair> m_peppairs = new HashSet<PeptidePair>();
    /** where in which proteins can this peptide be found */
    private HashMap<Protein,HashSet<Integer>> m_positions = new HashMap<Protein, HashSet<Integer>>();
    //int posSize = 0;
    /** unique ID for this peptide */
    private long id;
    /** assigned FDR for this peptide */
    public double m_fdr = -1;
    /** mass of this peptide */
    public double mass = Double.NaN;
    
    

    /**
     * instantiates a new peptide
     * @param id unique id 
     * @param sequence aa-sequence
     * @param isDecoy is it a decoy peptide
     * @param protein protein it could come from
     * @param pos position in the protein
     * @param length the length of the peptide
     */
    public Peptide(long id, String sequence, boolean isDecoy, Protein protein, int pos, int length) {
        this(id, sequence, isDecoy, new Protein[]{protein}, new int[]{pos}, length);
    }
    
    /**
     * instantiates a new peptide
     * <p>The protein and position array have to have the same size. E.g if a 
     * peptide could come from two sites within a protein, the protein has to be
     * twice in the protein array as well.</p>
     * @param id unique id 
     * @param sequence aa-sequence
     * @param isDecoy is it a decoy peptide
     * @param protein proteins it could come from
     * @param pos positions within the proteins
     * @param length the length of the peptide
     */
    public Peptide(long id, String sequence, boolean isDecoy, Protein[] protein, int[] pos, int length) {
        this.sequence = sequence;
        this.isDecoy = isDecoy;
        this.length = length;
        this.id = id;
                
        if (protein.length == 1) {
            HashSet<Integer> ppos = new HashSet<>();
            ppos.addAll(RArrayUtils.toCollection(pos));
            m_positions.put(protein[0], ppos);
        } else {
            if (pos.length != protein.length) {
                throw new UnsupportedOperationException("Cant handle entries with non-matching peptide positions and proteins");
            }
            for (int p = 0; p < protein.length; p++) {
                HashSet<Integer> ppos = m_positions.get(protein[p]);
                if (ppos == null) {
                    ppos = new HashSet<>();
                    m_positions.put(protein[p], ppos);
                }
                ppos.add(pos[p]);
            }
            
        }
//        posSize = pos.length;
    }

    @Override
    public Site getLinkSite1() {
        return new PeptideSite(this,-1);
    }

    @Override
    public Site getLinkSite2() {
        return null;
    }
    
    @Override
    public int hashCode() {
        return Boolean.hashCode(isDecoy) +  31 * sequence.hashCode();
    }

    /**
     * Tests whether two peptides are equal to each other.
     * <p> here this means it has the same sequence and both are either traget
     * or decoy proteins</p>
     * @param o
     * @return 
     */
    @Override
    public boolean equals(Object o) {
        return ((Peptide) o).id == id || (((Peptide) o).sequence.equals(sequence) && ((Peptide) o).isDecoy == isDecoy);
    }

//    public int compareTo(Peptide o) {
//        return Double.compare(o.getScore(), getScore());
//    }
//
    /**
     * compares two peptide
     * 
     * @param p
     * @return 0 if p.equals(this) returns true otherwise the result of an 
     * integer comparison of the IDs
     */
    public int compare(Peptide p) {
        if (p.equals(this)) {
            return 0;
        }
        if (this.id < p.id) {
            return -1;
        }
        return 1;
    }
    
    
    /**
     * Currently not considered for peptides outside peptide pairs
     * @return 
     */
    @Override
    public double getScore() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    /**
     * Currently not considered for peptides outside peptide pairs
     * @return 
     */
    @Override
    public double getScore(int topN) {
        throw new UnsupportedOperationException("Not supported yet.");
    }    
    /**
     * Currently not considered for peptides outside peptide pairs
     * @return 
     */
    @Override
    public String getFDRGroup() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
//    /**
//     * Currently not considered for peptides outside peptide pairs
//     * @return 
//     */
//    @Override
//    public String getFDRGroupName() {
//        throw new UnsupportedOperationException("Not supported yet.");
//    }

    /**
     * is this a target peptide
     * @return 
     */
    @Override
    public boolean isTT() {
        return !isDecoy;
    }

    /**
     * is this a decoy peptide
     * @return 
     */
    @Override
    public boolean isTD() {
        return isDecoy;
    }
    
    /**
     * single peptides are never considered as decoy decoy - would 
     * screw with the calculation
     * @return 
     */
    @Override
    public boolean isDD() {
        return false;
    }

    /**
     * add the peptide sources from the given peptide to the current instance.
     * @param o 
     */
    @Override
    public void add(Peptide o) {
        if (o == this)
            return;
        m_peppairs.addAll(o.m_peppairs);
        for (Protein p : o.m_positions.keySet()) {
            for (int pos : o.m_positions.get(p)) {
                add(p,pos);
            }
        }
    }

    /**
     * add a new source for the peptide 
     * position was already registered for a protein
     * @param p
     * @param pos 
     */
    public void add(Protein p, int pos) {
        HashSet<Integer> positions = m_positions.get(p);
        if (positions == null) {
            positions = new HashSet<Integer>(1);
            m_positions.put(p, positions);
        }
        positions.add(pos);
//        posSize++;
    }

    /**
     * add a peptidepair to the matches supporting this peptide
     * @param o 
     */
    public void add(PeptidePair o) {
        this.m_peppairs.add(o);
    }
    
    /**
     * Could does this peptide come from the same protein as the given peptide.
     * Proteins and the derived decoy are considered the same.
     * @param pep
     * @return 
     */
    public boolean sameProtein(Peptide pep) {
        // is the same protein in both sources ?
        for (Protein prot : pep.m_positions.keySet()) {
            if (m_positions.containsKey(prot))
                return true;
                    
        }
        // test if the we have an overlap with the decoy complement
        for (Protein prot : pep.m_positions.keySet()) {
            if (m_positions.containsKey(prot.decoyComplement()))
                return true;
        }
        return false;
    }
    
    /**
     * All proteins that this peptide could have come from.
     * @return 
     */
    public ArrayList<Protein> getProteins() {
         ArrayList<Protein> ret =new ArrayList<Protein>(m_positions.size());
         ret.addAll(m_positions.keySet());
         return ret;
    }
    
    /**
     * Returns a ProteinGroup consisting of all proteins where the peptides 
     * could come from
     * @return 
     */
    public ProteinGroup getProteinGroup() {
         ProteinGroup pg = new ProteinGroup(getProteins(), m_peppairs);
         return pg;
    }
    
    
    /**
     * length of the peptide
     * @return 
     */
    public int length() {
        return length;
    }

    /**
     * all matched PeptidePairs containing this peptide.
     * @return 
     */
    public Collection<PeptidePair> getPeptidePairs() {
        return m_peppairs;
    }


    /**
     * returns a collection with only the current instance 
     * @return 
     */
    public Collection<Peptide> getPeptides() {
        ArrayList<Peptide> ret = new ArrayList<Peptide>(1);
        ret.add(this);
        return ret;
    }

    /**
     * Returns the amino-acid sequence of the peptide
     * @return the sequence
     */
    public String getSequence() {
        return sequence;
    }

    /**
     * @return the sequence without any non-capital letter
     */
    public String getUppercaseSequence() {
        return sequence.replaceAll("[^A-Z]", "");
    }
    
    /**
     * A string representation of the peptide
     * @return 
     */
    @Override
    public String toString() {
        if (isDecoy) 
            return "D - " + sequence;
        return "T - " + sequence;
    }
    
    /**
     * is this a decoy peptide
     * @return 
     */
    @Override
    public boolean isDecoy() {
        return isDecoy;
    }    

    /**
     * the FDR represented by the score
     * @param fdr 
     */
    @Override
    public void setFDR(double fdr) {
        m_fdr = fdr;
    }

    /**
     * the FDR represented by the score
     * @return 
     */
    @Override
    public double getFDR() {
        return m_fdr;
    }

    /**
     * remove support information
     */
    public void resetFDR() {
        m_fdr= Double.MAX_VALUE;
        m_peppairs.clear();
        //m_positions = new HashMap<Protein, IntArrayList>();
        setLinkedSupport(1);
    }
    

    /**
     * how many peptidepairs containing this peptide where matched 
     * @return 
     */
    @Override
    public int getPeptidePairCount() {
        return m_peppairs.size();
    }

    /**
     * All positions of the peptide in proteins
     * @return the m_positions
     */
    public HashMap<Protein,HashSet<Integer>> getPositions() {
        return m_positions;
    }

    public void renewPositions() {
        HashMap<Protein, HashSet<Integer>> newpp = new HashMap<>();
        for (Map.Entry<Protein, HashSet<Integer>> e : m_positions.entrySet()) {
            newpp.put(e.getKey(),e.getValue());
        }
        m_positions = newpp;
    }
    
    /**
     * get a string with all accession numbers for the the originating Proteins
     * 
     * @return 
     */
    public String getAccessions() {
        ArrayList<String> accessStrings = new ArrayList<String>();
        for (Protein p : m_positions.keySet()) {
            if (isDecoy) {
                for (int i = m_positions.get(p).size(); i>0 ; i --) {

                    accessStrings.add("DECOY:"+p.getAccession());
                }
            } else {
                for (int i = m_positions.get(p).size(); i>0 ; i --) {

                    accessStrings.add(p.getAccession());
                }
                
            }
        }
        return RArrayUtils.toString(accessStrings, ";");
    }

    /**
     * get a string with all descriptions for the the originating Proteins
     * 
     * @return 
     */
    public String getDescriptions() {
        ArrayList<String> descStrings = new ArrayList<String>();
        for (Protein p : m_positions.keySet()) {
            if (isDecoy) {
                for (int i = m_positions.get(p).size(); i>0 ; i --) {

                    descStrings.add("DECOY:"+p.getDescription());
                }
            } else {
                for (int i = m_positions.get(p).size(); i>0 ; i --) {

                    descStrings.add(p.getDescription());
                }
                
            }
        }
        return RArrayUtils.toString(descStrings, ";");
    }

    /**
     * get a string with all positions of the peptide in the Proteins
     * 
     * @return 
     */
    public String getStringPositions() {
        return getStringPositions(0);
    }

    /**
     * get a string with all positions of the peptide in the Proteins plus the 
     * given offset.
     * <p>The offset can be used to get the linkage positions within proteins 
     * for e.g. peptide-pairs</p>
     * 
     * @param offset
     * @return 
     */
    public String getStringPositions(int offset) {
        StringBuilder ret = new StringBuilder();
        for (Protein p : m_positions.keySet()) {
            for (int i : m_positions.get(p)) {
                ret.append(";");
                ret.append(i+offset);
            }
        }
        return ret.substring(1);
    }
        
    /**
     * as there is only one site - this function always returns the Peptide 
     * instance itself.
     * @param n
     * @return 
     */
    @Override
    public Object getSite(int n) {
        return this;
    }

    /**
     * This is a single peptide and therfore has only one site
     * @return 1
     */
    @Override
    public int getSites() {
        return 1;
    }

    /**
     * @return the unique peptide id
     */
    public long getId() {
        return id;
    }

    /**
     * @param id the peptide-id
     */
    public void setId(long id) {
        this.id = id;
    }

    /**
     * single peptide - therefore always linear
     * @return 
     */
    @Override
    public boolean isLinear() {
        return true;
    }

    /**
     * no linkage - so not internal
     * @return 
     */
    @Override
    public boolean isInternal() {
        return false;
    }
    
    /**
     * no linkage not between either
     * @return 
     */
    @Override
    public boolean isBetween() {
        return false;
    }
    
    /**
     * mass of the peptide
     * @param mass 
     */
    public void setMass(double mass) {
        this.mass = mass;
    }
    /**
     * mass of the peptide
     * @return 
     */
    public double getMass() {
        return this.mass;
    }

    @Override
    public ProteinGroup getProteinGroup1() {
        return getProteinGroup();
    }

    @Override
    public ProteinGroup getProteinGroup2() {
        return ProteinGroup.NOPROTEINGROUP;
    }

    @Override
    public void setFDRGroup(String fdrGroup) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    
//    @Override
//    public boolean hasPositiveGrouping() {
//        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
//    }
//    
//    @Override
//    public void setPositiveGrouping(String av) {
//        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
//    }
//
//    @Override
//    public HashSet<String> getPositiveGrouping() {
//        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
//    }
//
//    @Override
//    public boolean hasNegativeGrouping() {
//        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
//    }
//
//    @Override
//    public void setNegativeGrouping(String v) {
//        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
//    }
//
//    @Override
//    public HashSet<String> getNegativeGrouping() {
//        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
//    }

    public boolean isNonCovalent() {
        return false;
    }
    
}
