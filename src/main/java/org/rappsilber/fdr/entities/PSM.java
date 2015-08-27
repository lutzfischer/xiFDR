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
import java.util.HashSet;
import org.rappsilber.fdr.groups.ProteinGroup;
import org.rappsilber.fdr.utils.AbstractFDRElement;
import org.rappsilber.utils.SelfAddHashSet;


/**
 * represents a single PSM
 */
public class PSM extends AbstractFDRElement<PSM> { 
    /**
     * unique id for the PSM.
     */
    private String psmID;
    /**
     * Run name
     */
    private String run;
    /**
     * scan-number
     */
    private String scan;
    /**
     * the peptide pair represented by this PSM
     */
    private PeptidePair m_peppair;
    
    /**
     * the link represented by this PSM
     */
    private ProteinGroupLink m_link;
    /**
     * the protein pair represented by this PSM
     */
    private ProteinGroupPair m_ppi;

    /** link-site in peptide1 */
    private int pepsite1;
    /** link-site in peptide2 */
    private int pepsite2;
    /** length of  peptide1 */
    private int peplen1;
    /** length of  peptide2 */
    private int peplen2;
    /** positions of peptide1 */
    Peptide peptide1;
    /** positions of peptide2 */
    Peptide peptide2;
    /** overall score of this pair */
    private double score;
    /** top-score per charge state */
    private int charge = 0;
    /**
     * M/Z value of the spectrum
     */
    private double experimentalMZ;
    /**
     * Mass value of the matched peptide(pair)
     */
    private double calcMass;
    /**
     * Charge value of the matched peptide(pair)
     */
    private int match_charge;
    
    /**
     * 
     */
    private double scoreRatio;
    /** hash code for fast access in hashmaps */
    int hashcode;
    /** first peptide comes from a decoy sequence */
    private boolean isDecoy1;
    /** second peptide comes from a decoy sequence */
    private boolean isDecoy2;
    /** it could be a protein internal link */
    private boolean isInternal = false;
    /** is this a linear "peptide pair" - meaning only one peptide */
    private boolean isLinear = false;
    /** is this a loop-linked peptide */
    private boolean isLoop = false;
    /** are all peptides from the target database */
    private boolean isTT = false;
    /** is one peptide from the decoy database */
    private boolean isTD = false;
    /** are all peptides from the decoy database */
    private boolean isDD = false;
    /** an id for the fdr-group of this PSM */
    private int fdrGroup;
    /** what is the calculated FDR for the score of this PSM */
    public double m_fdr = -1;
    /** protein group for the first peptide*/
    public ProteinGroup fdrProteinGroup1 = null;
    /** protein group for the second peptide */
    public ProteinGroup fdrProteinGroup2 = null;
    /** 
     * Is this a special case. E.g. PSMs that where matched with unknown 
     * charge-state have a inherently higher chance to be wrong. Therefore it
     * makes sense to calculate the FDR for these separately
     */
    private boolean specialcase = false;
    
    
    
    /**
     * creates a new instance of a PSM.
     * @param psmID a unique ID for the PSM (e.g. used when a CSV-file contains 
     * several rows for a single PSM that represent each possible origin of that 
     * PSM).
     * @param run run-name
     * @param scan scan-id
     * @param peptide1 first peptide
     * @param peptide2 second peptide
     * @param site1 link-site in the first peptide
     * @param site2 link site in the second peptide
     * @param isDecoy1 is the first peptide a decoy?
     * @param isDecoy2 is the second peptide a decoy?
     * @param charge precursor charge state
     * @param score score for this PSM
     * @param scoreRatio how to split the score between the peptides - this is 
     * only used, if a protein FDR (each single protein) is to be calculated.
     */
    public PSM(String psmID, String run, String scan, Peptide peptide1, Peptide peptide2, int site1, int site2, boolean isDecoy1, boolean isDecoy2, int charge, double score, double scoreRatio) {
        this(psmID, peptide1, peptide2, site1, site2, isDecoy1, isDecoy2, charge, score, scoreRatio);
        this.run = run;
        this.scan = scan;
    }
    
    
    /**
     * creates a new instance of a PSM.
     * @param psmID a unique ID for the PSM (e.g. used when a CSV-file contains 
     * several rows for a single PSM that represent each possible origin of that 
     * PSM).
     * @param peptide1 first peptide
     * @param peptide2 second peptide
     * @param site1 link-site in the first peptide
     * @param site2 link site in the second peptide
     * @param isDecoy1 is the first peptide a decoy?
     * @param isDecoy2 is the second peptide a decoy?
     * @param charge precursor charge state
     * @param score score for this PSM
     * @param scoreRatio how to split the score between the peptides - this is 
     * only used, if a protein FDR (each single protein) is to be calculated.
     */
    public PSM(String psmID, Peptide peptide1, Peptide peptide2, int site1, int site2, boolean isDecoy1, boolean isDecoy2, int charge, double score, double scoreRatio) {
        this.psmID = psmID;
//        this.id1 = id1;
//        this.id2 = id2;
        
        this.pepsite1 = site1;
        this.pepsite2 = site2;
        this.score = score * score;
        this.isDecoy1 = isDecoy1;
        this.isDecoy2 = isDecoy2;
        if (isDecoy1) {
            if (isDecoy2) {
                isDD = true;
            } else {
                isTD = true;
            }
        } else if (isDecoy2) {
            isTD = true;
        } else {
            isTT = true;
        }
        this.charge = charge;
        this.score = score;
        this.scoreRatio = scoreRatio;
        this.peptide1 = peptide1;
        this.peptide2 = peptide2;
        
        this.hashcode = (peptide1.hashCode() + peptide2.hashCode()) % 100000 * 31 + (pepsite1 + pepsite2) % 10;
        
        // first distinction is whether it could be a protein internal link
        this.isLinear = peptide1 == Peptide.NOPEPTIDE || peptide2 == Peptide.NOPEPTIDE;
        this.isLoop = this.isLinear && site1 >= 0 && site2 >= 0;
        this.isInternal = peptide1.sameProtein(peptide2);
        
//        setFDRGroup();
        
    }

    @Override
    public Site getLinkSite1() {
        if (peptide1 == Peptide.NOPEPTIDE)
            return null;
        return new PeptideSite(peptide1,pepsite1);
    }

    @Override
    public Site getLinkSite2() {
        if (peptide2 == Peptide.NOPEPTIDE)
            return null;
        return new PeptideSite(peptide2,pepsite2);
    }

    /**
     * hash code of this PSM
     * @return 
     */
    @Override
    public int hashCode() {
        return hashcode;
    }

    /**
     * tests, whether to PSMs are the same.
     * Used to join up PSM entries that where generated for different source of 
     * the peptides (ambiguous peptides/non prototypic peptides)
     * @param l
     * @return 
     */
    @Override
    public boolean equals(Object l) {
        PSM c = (PSM) l;
//        return this.score == c.score && this.charge == c.charge &&  this.psmID.contentEquals(c.psmID);
        return this.score == c.score && this.charge == c.charge &&  this.psmID.contentEquals(c.psmID) &&
                (((c.scoreRatio == this.scoreRatio || (Double.isNaN(c.scoreRatio) && Double.isNaN(this.scoreRatio))) && c.peptide1.equals(this.peptide1) && c.peptide2.equals(this.peptide2) && c.pepsite1 == pepsite1 && c.pepsite2 == pepsite2) 
                || /* to be safe from binary inaccuracy we make an integer-comparison*/
                ((Math.round(100000*c.scoreRatio) == Math.round(100000-100000*this.scoreRatio) || (Double.isNaN(c.scoreRatio) && Double.isNaN(this.scoreRatio))) && c.peptide2.equals(this.peptide1) && c.peptide1.equals(this.peptide2) && c.pepsite2 == pepsite1 && c.pepsite1 == pepsite2));
    }

    /**
     * get a peptide pair represented by this PSM
     * @return 
     */
    public PeptidePair getPeptidePair() {
                             //psmID, id1, id2, site1,    site2,    isDecoy1, isDecoy2, peptide1Positions, peptide2Positions, charge, double score, double scoreRatio
        return new PeptidePair(this);
    }
    
    /**
     * get all peptidepairs represented by this PSM - it's exactly one :)
     * @return 
     */    
    public ArrayList<PeptidePair> getPeptidePairs() {
        ArrayList<PeptidePair> ret =  new ArrayList<PeptidePair>(1);
        ret.add(getPeptidePair());
        return ret;
    }


    /**
     * get all proteins involved with this PSM 
     * @return 
     */
    public ArrayList<Protein> getProteins() {
        HashSet<Protein> prots = new HashSet<Protein>();
        prots.addAll(peptide1.getProteins());
        prots.addAll(peptide2.getProteins());
        ArrayList<Protein> ret = new ArrayList<Protein>(prots);
        return ret;
    }
    
    /**
     * get the proteingroups for both peptides
     * @return 
     */
    public ArrayList<ProteinGroup> getProteinGroupss() {
        HashSet<ProteinGroup> prots = new HashSet<ProteinGroup>();
        for (Peptide p : getPeptides()) {
            if (p != Peptide.NOPEPTIDE) {
                prots.add(p.getProteinGroup());
            }
        }
        ArrayList<ProteinGroup> ret = new ArrayList<ProteinGroup>(prots);
        return ret;
    }
    
    /**
     * Return the supporting proteins, but filter them first through the 
     * given {@link SelfAddHashSet} to make sure that all PSMs return the same
     * instance of the proteins
     * @param allProteins
     * @return 
     */
    public ArrayList<Protein> getProteins(SelfAddHashSet<Protein> allProteins) {
        return allProteins.registerAll(getProteins());
    }


//    public int compareTo(PSM o) {
//        return Double.compare(o.getScore(), this.getScore());
//    }

    /**
     * adds the information of another PSM to this one.
     * @param p 
     */
    @Override
    public void add(PSM p) {
        boolean invertScoreRatio = false;
        if (p == this)
            return;
        if (!equals(p))
            throw new UnsupportedOperationException("Can only add the same PSM together (e.g. for assigning other protein positions)");
        
        
        if (p.isInternal && !isInternal) {
            isInternal = true;
          //  setFDRGroup();
        }
    }

    /**
     * @return the psmID
     */
    public String getPsmID() {
        return psmID;
    }

    /**
     * @return the id of the first peptide
     */
    public Peptide getPeptide1() {
        return peptide1;
    }

    /**
     * @return the id of the second peptide
     */
    public Peptide getPeptide2() {
        return peptide2;
    }

    /**
     * @return the link site of the first peptide
     */
    public int getPeptideLinkSite1() {
        return pepsite1;
    }

    /**
     * @return the link site of the second peptide
     */
    public int getPeptideLinkSite2() {
        return pepsite2;
    }

    /**
     * @param peptide
     * @return the link site of the given peptide
     */
    public int getPeptideLinkSite(int peptide) {
        return peptide == 0 ? pepsite1 : pepsite2;
    }
    

    /**
     * @return the score
     */
    @Override
    public double getScore() {
        return score;
    }

    /**
     * @return the precursor charge state of the PSM
     */
    public int getCharge() {
        return charge;
    }

    /**
     * Score ratio is used to split the score for the as support for the matched 
     * protein groups
     * @return the scoreRatio
     */
    public double getScoreRatio() {
        return scoreRatio;
    }

    /**
     * @return is the first peptide a decoy
     */
    public boolean isDecoy1() {
        return isDecoy1;
    }

    /**
     * @return is the second peptide a decoy
     */
    public boolean isDecoy2() {
        return isDecoy2;
    }

    /**
     * This returns whether the PSM could be a protein-internal match. Meaning 
     * is there an overlap between the Protein-groups for each peptide.
     * @return whether the PSM could be a protein-internal match.
     * 
     */
    @Override
    public boolean isInternal() {
        return isInternal;
    }

    /**
     * Are all peptides target peptides
     * @return the isTT
     */
    @Override
    public boolean isTT() {
        return isTT;
    }

    /**
     * is one peptide from the decoy database?
     * @return TD
     */
    @Override
    public boolean isTD() {
        return isTD;
    }

    /**
     * are all peptides decoy peptides
     * @return the isDD
     */
    @Override
    public boolean isDD() {
        return isDD;
    }

    /**
     * numeric id for the FDR-group of this PSM
     * @return 
     */
    @Override
    public int getFDRGroup() {
        return fdrGroup;
    }

    /**
     * a name for the FDR-group
     * @return 
     */
    @Override
    public String getFDRGroupName() {
        return PeptidePair.fdrGroupNames.get(fdrGroup);
    }
    
    /**
     * Returns the name for the given FDR-group-id
     * @param fdrgroup
     * @return 
     */
    public static String getFDRGroupName(int fdrgroup) {
        return PeptidePair.getFDRGroupName(fdrgroup);
    }    
    
    /**
     * Returns a list of links supported by this PSM.
     * As the ambiguity is handled in the Protein groups this will report 
     * exactly one link.
     * @return 
     */
    public Collection<ProteinGroupLink> getLinks() {
        return this.getPeptidePair().getLinks();
    }

    
    /**
     * return all peptides involved in this PSM
     * @return 
     */
    public Collection<Peptide> getPeptides() {
        ArrayList<Peptide> ret = new ArrayList<Peptide>(2);
        if (getPeptide2() != Peptide.NOPEPTIDE) {
            ret.add(getPeptide2());
        }
        ret.add(getPeptide1());
        return ret;
    }
    
    /**
     * A textual representation of this PSM
     * @return 
     */
    @Override
    public String toString() {
        
        if (isTT())
            return "TT - " + peptide1.toString() + " xl " + peptide2.toString();
        if (isTD())
            return "TD - " + peptide1.toString() + " xl " + peptide2.toString();
        
        return "DD - " + peptide1.toString() + " xl " + peptide2.toString();
        
    }
    
    /**
     * is at least one of the peptides a decoy
     * @return 
     */
    @Override
    public boolean isDecoy() {
        return isDecoy1 || isDecoy2;
    }
    
    
    /**
     * set the FDR group according to the information on this PSM
     */
    public void setFDRGroup() {
        
        fdrGroup = PeptidePair.getFDRGroup(peptide1, peptide2, isLinear(), isInternal, specialcase);
       
    }     
    
    /**
     * store the FDR assigned to the score of this PSM
     * @param fdr 
     */
    @Override
    public void setFDR(double fdr) {
        m_fdr = fdr;
    }

    
    /**
     * The FDR assigned to the score of this PSM
     * @return the score FDR
     */
    @Override
    public double getFDR() {
        return m_fdr;
    }


    /**
     * @return the length of peptide 1
     */
    public int getPeptideLength1() {
        return peptide1.length();
    }

    /**
     * @return the length of peptide 2
     */
    public int getPeptideLength2() {
        return peptide2.length();
    }

    /**
     * @param peptide
     * @return the length of the given peptide
     */
    public int getPeptideLength(int peptide) {
        return peptide ==0 ? peptide1.length() : peptide2.length() ;
    }

    

    /**
     * Set the peptide pair that was supported by this PSM and passed the 
     * peptide pair FDR
     * @param pp
     */
    public void setFdrPeptidePair(PeptidePair pp) {
        this.m_peppair = pp;
    }

    /**
     * Get the peptide pair that was supported by this PSM and passed the 
     * peptide pair FDR
     * @return 
     */
    public PeptidePair getFdrPeptidePair() {
        return this.m_peppair;
    }
    
    /**
     * reset the stored information, that might be different for a new FDR-calculation
     */
    public void reset() {
        resetFdrProteinGroup();
        setFdrPeptidePair(null);
        setLinkedSupport(1);
    }
    
    /**
     * unset the fdrproteingroups
     */
    public void resetFdrProteinGroup() {
            this.fdrProteinGroup1 = null;
            this.fdrProteinGroup2 = null;
        
    }
    
    
    /**
     * set the ProteinGroup that passed the fdr and is supported by this PSM
     * @param pg
    */
    public void setFdrProteinGroup(ProteinGroup pg) {
        if (pg.equals(peptide1.getProteinGroup()))
            this.fdrProteinGroup1 = pg;
        if (pg.equals(peptide2.getProteinGroup()))
            this.fdrProteinGroup2 = pg;
    }
    
    /**
     * get the ProteinGroup that passed the fdr and is supported by the first
     * peptide of this PSM 
     * @return 
     */
    public ProteinGroup getFdrProteinGroup1() {
        return this.fdrProteinGroup1;
    }
    

    /**
     * get the ProteinGroup that passed the fdr and is supported by the second
     * peptide of this PSM 
     * @return 
     */
    public ProteinGroup getFdrProteinGroup2() {
        return this.fdrProteinGroup2;
    }

    /** 
     * number of supporting peptidepairs - is for a psm always 1.
     * @return  
     */
    @Override
    public int getPeptidePairCount() {
        return 1;
    }

    /**
     * Is this a special case. E.g. PSMs that where matched with unknown 
     * charge-state have a inherently higher chance to be wrong. Therefore it
     * makes sense to calculate the FDR for these separately
     * @return the specialcase
     */
    public boolean isSpecialcase() {
        return specialcase;
    }

    /**
     * Is this a special case. E.g. PSMs that where matched with unknown 
     * charge-state have a inherently higher chance to be wrong. Therefore it
     * makes sense to calculate the FDR for these separately
     * @param specialcase the specialcase to set
     */
    public void setSpecialcase(boolean specialcase) {
        this.specialcase = specialcase;
    }

    /**
     * Is this a non-cross-linked PSM?
     * @return the isLinear
     */
    @Override
    public boolean isLinear() {
        return isLinear;
    }

    /**
     * Is this a non-cross-linked PSM?
     * @param isLinear 
     */
    public void setLinear(boolean isLinear) {
        this.isLinear = isLinear;
    }

    /**
     * the run-name for the PSM
     * @return the run
     */
    public String getRun() {
        return run;
    }

    /**
     * the run-name for the PSM
     * @param run the run to set
     */
    public void setRun(String run) {
        this.run = run;
    }

    /**
     * the scan-id for the PSM
     * @return the scan
     */
    public String getScan() {
        return scan;
    }

    /**
     * the scan-id for the PSM
     * @param scan the scan to set
     */
    public void setScan(String scan) {
        this.scan = scan;
    }

    /**
     * Get the Peptide for the given site.
     * @param n
     * @return 
     */
    @Override
    public Object getSite(int n) {
        return n==0 ? getPeptide1() : getPeptide2();
    }

    /**
     * how many sites does this PSM have.
     * Is either 2 for cross-linked PSMs or 1 for linear PSMs
     * @return 
     */
    @Override
    public int getSites() {
        if (isLinear)
            return 1;
        else 
            return 2;
    }

    /**
     * is this a cross-linked PSM where both peptides come definitely from two 
     * different Proteins
     * @return 
     */
    @Override
    public boolean isBetween() {
        return !isLinear && !isInternal;
    }

    /**
     * is this a linear loop-linked peptide PSM
     * @return the isLoop
     */
    public boolean isLoop() {
        return isLoop;
    }

    /**
     * is this a linear loop-linked peptide PSM
     * @param isLoop the isLoop to set
     */
    public void setLoop(boolean isLoop) {
        this.isLoop = isLoop;
    }

    /**
     * M/Z value of the spectrum
     * @return the experimentalMZ
     */
    public double getExperimentalMZ() {
        return experimentalMZ;
    }

    /**
     * M/Z value of the spectrum
     * @param experimentalMZ the experimentalMZ to set
     */
    public void setExperimentalMZ(double experimentalMZ) {
        this.experimentalMZ = experimentalMZ;
    }

    /**
     * Mass value of the matched peptide(pair)
     * @return the calcMZ
     */
    public double getCalcMass() {
        return calcMass;
    }

    /**
     * Mass value of the matched peptide(pair)
     * @param calcMass
     */
    public void setCalcMass(double calcMass) {
        this.calcMass = calcMass;
    }

    /**
     * Charge value of the matched peptide(pair)
     * @return the match_charge
     */
    public int getMatchedCharge() {
        return match_charge;
    }

    /**
     * Charge value of the matched peptide(pair)
     * @param match_charge the match_charge to set
     */
    public void setMatchedCharge(int match_charge) {
        this.match_charge = match_charge;
    }
    
}
