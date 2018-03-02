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
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.rappsilber.fdr.groups.ProteinGroup;
import org.rappsilber.fdr.utils.AbstractFDRElement;
import org.rappsilber.fdr.utils.FDRGroupNames;
import org.rappsilber.utils.DoubleArrayList;
import org.rappsilber.utils.RArrayUtils;
import org.rappsilber.utils.SelfAddHashSet;

/**
 * represents a peptide pair and stores the support for the given peptide pair.
 * <br/>
 * Linear peptides will also be presented by PeptidePair (to ease a unified FDR)
 * by setting the second peptide to the {@link Peptide#NOPEPTIDE} place-holder.
 */
public class PeptidePair extends AbstractFDRElement<PeptidePair> {//implements Comparable<PeptidePair>{
    /** counts up each instance of PeptidePairs used to assign unique IDs */
    public static int PEPTIDEPAIRCOUNT = 0;
    /** a link represented by this peptide pair that passed the fdr */
    private ProteinGroupLink m_link;
    /** unique id for this peptide pair */
    protected int peptidePairID = PEPTIDEPAIRCOUNT++;
    /** sorry is used in the rappsilber group for a search type specific FDR */
    public static boolean ISTARGETED = false;
    /** length groups used for doing a length depended FDR-grouping */
    protected static int[] lenghtGroup;
//    /** meaningful names for the FDR-groups */
//    public static HashMap<Integer, String> fdrGroupNames = new HashMap<Integer, String>();
    /** is used in the rappsilber group internally for a a specific search type */ 
    private static Pattern targetMod = Pattern.compile("X([0-9]+(\\.[0-9]+)?)");
    
    /**
     * id of peptide1
     */
    private Peptide peptide1;
    /**
     * id of peptide2
     */
    private Peptide peptide2;
    /**
     * link-site in peptide1
     */
    private byte pepsite1;
    /**
     * link-site in peptide2
     */
    private byte pepsite2;

    /**
     * overall score of this pair
     */
    private double score;
    
    
    /**
     * Cross-linker
     */
    private String crosslinker;
    /** top-score per charge state */
    private DoubleArrayList chargeTopScoresRatios = new DoubleArrayList();
    /**
     * top match per chare state
     */
    //private HashMap<Integer, PSM> chargeTopScoresPSM = new HashMap<Integer, PSM>();
    private PSM[] chargeTopScoresPSM;
    /**
     * PSMs supporting this peptide pair
     */
    private ArrayList<PSM> psms = new ArrayList<PSM>();
    /**
     * hash code for fast access in hashmaps
     */
    int hashcode;
    /**
     * first peptide comes from a decoy sequence
     */
    private boolean isDecoy1;
    /**
     * second peptide comes from a decoy sequence
     */
    private boolean isDecoy2;
    /**
     * it could be a protein internal link
     */
    private boolean isInternal = false;
    /** only one peptide */
    protected boolean isLinear = false;
    /** single loop-linked peptide */
    private boolean isLoop = false;
    /** no decoy */
    private boolean isTT = false;
    /** one decoy */
    private boolean isTD = false;
    /** two decoy */
    private boolean isDD = false;

    /** FDR group assigned to this peptide pair */
    private String fdrGroup;
    /** fdr assigned to the score of this peptidepair */
    public double m_fdr = -1;
    /** the fdr of the link that is supported by this peptide pair*/
    protected double m_LinkFdr = -1;
    /** the fdr of the protein pair that is supported by this peptide pair */
    protected double m_PPIFdr = -1;
    /** proteingroup supported by this peptidepair that passed the fdr */
    public ProteinGroup fdrProteinGroup1 = null;
    /** proteingroup supported by this peptidepair that passed the fdr */
    public ProteinGroup fdrProteinGroup2 = null;
    /** are all supporting PSMs special cases? */
    private ArrayList<String> specialcase = null;

    /** indicates that the two peptides where in the same spectrum but non-covalently linked */
    private boolean isNonCovalent = false;
    private String validated;
    
    
    /**
     * constructor
     * @param psm 
     */
    public PeptidePair(PSM psm) {
        this.peptide1 = psm.getPeptide1();
        this.peptide2 = psm.getPeptide2();
        this.pepsite1 = psm.getPeptideLinkSite1();
        this.pepsite2 = psm.getPeptideLinkSite2();
//        this.score = psm.getScore() * psm.getScore();
        this.isDecoy1 = psm.isDecoy1();
        this.isDecoy2 = psm.isDecoy2();
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
        crosslinker = psm.getCrosslinker();
//        this.chargeTopScores.set(psm.getCharge(), psm.getScore());
//        this.chargeTopScoresPSMId.set(psm.getCharge(), psm.getPsmID());
        this.psms.add(psm);
        this.chargeTopScoresPSM = new PSM[psm.getCharge()+1];
        this.chargeTopScoresPSM[psm.getCharge()] = psm;
        this.chargeTopScoresRatios.set(psm.getCharge(), psm.getScoreRatio());

        this.hashcode = (peptide1.hashCode() + peptide2.hashCode()) % 100000 * 10 + (pepsite1 + pepsite2) % 10;

        this.isLinear = peptide1 == Peptide.NOPEPTIDE || peptide2 == Peptide.NOPEPTIDE;
        this.isLoop = isLinear && pepsite1 >= 0 && pepsite2 >= 0;
        this.isInternal = peptide1.sameProtein(peptide2);
        if (psm.hasNegativeGrouping()) {
            this.specialcase = new ArrayList<>();
            this.specialcase.add(psm.getNegativeGrouping());
        }
        //this.score = psm.getScore();
        isNonCovalent = psm.isNonCovalent();

        setFDRGroup();
        this.validated = psm.getPositiveGrouping();

    }

    @Override
    public int hashCode() {
        return hashcode;
    }

    @Override
    public boolean equals(Object l) {
        PeptidePair c = (PeptidePair) l;
        if (isNonCovalent != c.isNonCovalent)
            return false;
        
        return  crosslinker == c.crosslinker && 
                ((c.peptide1.equals(peptide1) && c.peptide2.equals(peptide2) && c.pepsite1 == pepsite1 && c.pepsite2 == pepsite2)
                || (c.peptide2.equals(peptide1) && c.peptide1.equals(peptide2) && c.pepsite2 == pepsite1 && c.pepsite1 == pepsite2));
    }

    /** 
     * adds the information of another peptide pair to this one.
     * Is used to join the PeptidePairs that are generated from different PSMs 
     * but are actually the same PeptidePairs.
     * @param p 
     */
    public void add(PeptidePair p) {
        if (p == this) {
            return;
        }
        this.psms.addAll(p.psms);
        boolean invertScoreRatio = false;
        // transfer topScores and match-ids
        if (p.chargeTopScoresPSM.length > chargeTopScoresPSM.length) {
            PSM[] dummy = chargeTopScoresPSM;
            chargeTopScoresPSM = new PSM[p.chargeTopScoresPSM.length];
            System.arraycopy(dummy, 0, chargeTopScoresPSM, 0, dummy.length);
        }
        for (int i = 0 ; i< p.chargeTopScoresPSM.length;i++) {
            PSM ppsm = p.chargeTopScoresPSM[i];
            if (ppsm == null)
                continue;

            PSM cpsm =null; 
            if (chargeTopScoresPSM.length>i)
                cpsm= chargeTopScoresPSM[i];
            
            if (cpsm == null || ppsm.getScore() > cpsm.getScore()) {

                chargeTopScoresPSM[i]= ppsm;
                if (invertScoreRatio) {
                    chargeTopScoresRatios.set(i, 1 - p.chargeTopScoresRatios.get(i));
                } else {
                    chargeTopScoresRatios.set(i, p.chargeTopScoresRatios.get(i));
                }
                
            }
            
        }
        boolean setFDR = false;
        if (specialcase!=null && !p.hasNegativeGrouping()) {
            specialcase = null;
            setFDR = true;
        }
        
        if (p.isInternal && !isInternal) {
            isInternal = true;
            setFDR = true;
        } 
        
        if (p.hasPositiveGrouping()) {
            if (this.validated == null) {
                this.validated = p.getPositiveGrouping();
                setFDR = true;
            } else if (!this.validated.contentEquals(p.getPositiveGrouping())) {
                this.validated += " " + p.getPositiveGrouping();
                setFDR = true;
            }
        }

        if (setFDR)
            setFDRGroup();
        
    }

    /**
     * calculate the score of this Peptide pair based on the scores of the 
     * supporting PSMs
     */
    public void setScore() {
        score = 0;
        for (PSM psm : chargeTopScoresPSM) {
            if (psm!=null)
                score += psm.getScore() * psm.getScore();
        }
        score = Math.sqrt(score);
    }

    /**
     * returns the link that is supported by this peptidepair encapsulated in a 
     * ArrayList.
     * @return 
     */
    public ArrayList<ProteinGroupLink> getLinks() {
        ArrayList<ProteinGroupLink> links = new ArrayList<ProteinGroupLink>();
        links.add(new ProteinGroupLink(this));
        return links;
    }

    /**
     * returns the link that is supported by this peptidepair 
     * @return 
     */
    public ProteinGroupLink getLink() {
        return new ProteinGroupLink(this);
    }

    /**
     * returns the link supported by this peptidepair but filtered through the 
     * supplied SelfAddHashSet. This ensures that all peptidepairs that would 
     * support the same link will actually return references to the same 
     * instance of the ProteinGroupLink.
     * @param allLinks
     * @return 
     */
    public ArrayList<ProteinGroupLink> getLinks(SelfAddHashSet<ProteinGroupLink> allLinks) {
        return allLinks.registerAll(getLinks());
    }

    /**
     * returns a list of proteins supported by this peptide pair.
     * @return 
     */
    public ArrayList<Protein> getProteins() {
        HashSet<Protein> prots = new HashSet<Protein>();
        int sc1 = 0;
        int sc2 = 0;
        prots.addAll(peptide1.getProteins());
        prots.addAll(peptide2.getProteins());
        ArrayList<Protein> ret = new ArrayList<Protein>(prots);
        return ret;
    }

    /**
     * returns the proteins supported by this peptidepair but filtered through the 
     * supplied SelfAddHashSet. This ensures that all peptidepairs that would 
     * support the same proteins will actually return references to the same 
     * instances of the Protein.
     * @param allProteins 
     * @return 
     */
    public ArrayList<Protein> getProteins(SelfAddHashSet<Protein> allProteins) {
        return allProteins.registerAll(getProteins());
    }


    /**
     * calculate the score for the given protein based on the score of the 
     * current PeptidePair
     * @param prot
     * @return 
     */
    public double getProteinScore(Protein prot) {
        double score = 0;
        // is the protein among the first site of proteins
        ArrayList<Protein> prot1 = peptide1.getProteins();
        int prot1size = prot1.size();
        for (Protein p : prot1) {
            if (p.equals(prot)) {
                for (int c =0; c<chargeTopScoresPSM.length;c++) {
                    if (chargeTopScoresPSM[c]==null)
                        continue;
                    double cscore = chargeTopScoresPSM[c].getScore() * chargeTopScoresRatios.get(c, 0) / prot1size;
                    score += cscore * cscore;
                }
                break;
            }
        }
        // is the protein among the second site of proteins
        ArrayList<Protein> prot2 = peptide2.getProteins();
        int prot2size = prot2.size();
        for (Protein p : prot2) {
            if (p.equals(prot)) {
                for (int c =0; c<chargeTopScoresPSM.length;c++) {
                    if (chargeTopScoresPSM[c]==null)
                        continue;
                    double cscore = chargeTopScoresPSM[c].getScore() * chargeTopScoresRatios.get(c, 0) / prot2size;
                    score += cscore * cscore;
                }
                break;
            }
        }
        return Math.sqrt(score);
    }

    /**
     * calculate the score of a Peptide based on the current PeptidePair
     * @param pep
     * @return 
     */
    public double getPeptideScore(Peptide pep) {
        double score = 0;
        // is the protein among the first site of proteins
        if (pep.equals(peptide1)) {
                for (int c =0; c<chargeTopScoresPSM.length;c++) {
                    if (chargeTopScoresPSM[c]==null)
                        continue;
                double cscore = chargeTopScoresPSM[c].getScore() * chargeTopScoresRatios.get(c, 0);
                score += cscore * cscore;
            }
        }

        if (pep.equals(peptide2)) {
                for (int c =0; c<chargeTopScoresPSM.length;c++) {
                    if (chargeTopScoresPSM[c]==null)
                        continue;
                double cscore = chargeTopScoresPSM[c].getScore() * (1 - chargeTopScoresRatios.get(c, 0));
                score += cscore * cscore;
            }
        }

        return Math.sqrt(score);
    }

    /**
     * @return the first peptide of the pair
     */
    public Peptide getPeptide1() {
        return peptide1;
    }

    /**
     * the second peptide of the pair
     * @return 
     */
    public Peptide getPeptide2() {
        return peptide2;
    }

    /**
     * @return the link site in the first peptide
     */
    public int getPeptideLinkSite1() {
        return pepsite1;
    }

    /**
     * @return the link site in the second peptide
     */
    public int getPeptideLinkSite2() {
        return pepsite2;
    }

    /**
     * @param peptide
     * @return the link site for the given peptide
     */
    public int getPeptideLinkSite(int peptide) {
        return peptide == 0 ? pepsite1 : pepsite2;
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
     * @return the score
     */
    public double getScore() {
        if (score == 0) {
            setScore();
        }
        return score;
    }

    /**
     * whether the first peptide is decoy
     * @return 
     */
    public boolean isDecoy1() {
        return isDecoy1;
    }

    /**
     * whether the second peptide is decoy
     * @return 
     */
    public boolean isDecoy2() {
        return isDecoy2;
    }

    /**
     * whether this could be an protein internal link.<br/>
     * Meaning whether there is an overlap between the proteins that the first 
     * and the second peptide could originate from.
     * @return the isInternal
     */
    @Override
    public boolean isInternal() {
        return isInternal;
    }

    /**
     * Whether it is definitely not an internal and also not a linear peptide pair
     * @return the isInternal
     */
    @Override
    public boolean isBetween() {
        return (!isInternal) && (!isLinear);
    }
    
    /**
     * all peptides are target
     * @return the isTT
     */
    @Override
    public boolean isTT() {
        return isTT;
    }

    /**
     * one peptide is decoy
     * @return the isTD
     */
    @Override
    public boolean isTD() {
        return isTD;
    }

    /**
     * two peptides are decoy
     * @return the isDD
     */
    @Override
    public boolean isDD() {
        return isDD;
    }
    
    /**
     * the assigned FDR group for this peptide pair
     * @return 
     */
    @Override
    public String getFDRGroup() {
        return fdrGroup;
    }


//    /**
//     * translates FDR-group IDs into FDR-group names
//     * @param fdrgroup
//     * @return 
//     */
//    public static String getFDRGroupName(int fdrgroup) {
//        String name = fdrGroupNames.get(fdrgroup);
//        
//        if (name == null) {
//            return Integer.toString(fdrgroup);
//        }
//        return fdrGroupNames.get(fdrgroup);
//    }

    /**
     * returns all supporting peptide pairs - meaning itself
     * @return 
     */
    public Collection<PeptidePair> getPeptidePairs() {
        ArrayList<PeptidePair> ret = new ArrayList<PeptidePair>(1);
        ret.add(this);
        return ret;
    }

    /**
     * list of the peptides of this peptide pair
     * @return 
     */
    public ArrayList<Peptide> getPeptides() {
        ArrayList<Peptide> ret = new ArrayList<Peptide>(2);
        ret.add(getPeptide1());
        if (getPeptide2() != Peptide.NOPEPTIDE) {
            ret.add(getPeptide2());
        }
        return ret;
    }

    /**
     * returns a list of IDs of the top scoring PSMs for each precursor charge state 
     * of the PSMs that support this peptide pair.
     * @return 
     */
    public String getTopPSMIDs() {
        return RArrayUtils.toStringNoNull(chargeTopScoresPSM, ";");
    }
    
    /**
     * returns a list of the top scoring PSMs for each precursor charge state 
     * of the PSMs that support this peptide pair.
     * @return 
     */
    public Collection<PSM> getTopPSMs() {
        ArrayList<PSM> psms = new ArrayList<>();
        for (PSM p : chargeTopScoresPSM) 
            if (p!=null)
                psms.add(p);
        return psms;
    }

    /**
     * A string representation of this peptide pair
     * @return 
     */
    @Override
    public String toString() {

        if (isTT()) {
            return "TT - " + peptide1.toString() + " xl " + peptide2.toString();
        }
        if (isTD()) {
            return "TD - " + peptide1.toString() + " xl " + peptide2.toString();
        }

        return "DD - " + peptide1.toString() + " xl " + peptide2.toString();

    }

    /**
     * is any peptide a decoy
     * @return 
     */
    public boolean isDecoy() {
        return isDecoy1 || isDecoy2;
    }

    /**
     * @return length based grouping for fdr
     */
    public static int[] getLenghtGroup() {
        return lenghtGroup;
    }

    /**
     * set maximum length groups for fdr
     * @param aLenghtGroup 
     */
    public static void setLenghtGroup(int[] aLenghtGroup) {
        // make sure every value exists only ones
        HashSet<Integer> lengthSet = new HashSet<Integer>();
        boolean hasZero = false;
        for (int i = 0; i < aLenghtGroup.length; i++) {
            if (aLenghtGroup[i] == 0) {
                hasZero = true;
            }
            lengthSet.add(aLenghtGroup[i]);
        }

        // and that we have a zero in ther
        if (!hasZero) {
            lengthSet.add(0);
        }

        // turn it into a sorted array
        Integer[] dl = new Integer[lengthSet.size()];
        dl = lengthSet.toArray(dl);
        Arrays.sort(dl);
        lenghtGroup = new int[dl.length];
        for (int i = 0; i < dl.length; i++) {
            lenghtGroup[i] = dl[dl.length - 1 - i];
        }


//        // give groups a name
//        String[] groupNames = new String[]{"Linear", "Within", "Between", "Linear Special", "Within Special", "Between Special"};
//        fdrGroupNames.put(-1, "all combined");
//        for (int gn = 0; gn < groupNames.length; gn++) {
//            int g = gn * (lenghtGroup.length);
//            for (int i = 0; i < lenghtGroup.length; i++) {
//                fdrGroupNames.put(g + i, groupNames[gn] + "  >" + lenghtGroup[i]);
//            }
//        }

    }

    /** 
     * define a fdr-group id based on the  inputs
     * @param pep1
     * @param pep2
     * @param isLinear
     * @param isInternal
     * @param specialCase
     * @return 
     */
    public static String getFDRGroup(Peptide pep1, Peptide pep2, boolean isLinear, boolean isInternal, String specialCase, String groupExt) {
        String group = (isLinear ? "linear" : (isInternal ? "internal" :"between"));
        groupExt=" " + groupExt;
        //int metaGroup = (isLinear ? 0 : 1);//(isInternal ? 1 :2));
        if (specialCase != null) {
            groupExt+=" " + specialCase;
        }
        

        
       // int fdrGroup = metaGroup * lenghtGroup.length;
        

        if (ISTARGETED) { // for targetd modifications add the mass to the fdr group
            // sorry is used rappsilber internally for targeted modification search 
            Matcher m = targetMod.matcher(pep1.getSequence());
            if (m.matches()) {
                int mass = (int)(Math.round(Double.parseDouble(m.group(1))*10));
//                fdrGroup +=100 * mass;
                groupExt = Math.round(mass*100)/100 + groupExt;
                
            } else {
                m = targetMod.matcher(pep2.getSequence());
                if (m.matches()) {
                    int mass = (int)(Math.round(Double.parseDouble(m.group(1))*10));
                    groupExt = Math.round(mass*100)/100 + groupExt;
                }
            }
        }  
        
        if (isLinear) {
            // and the next decision is based on minimum peptide length
            for (int lg = 0; lg < lenghtGroup.length; lg++) {
                if (pep1.length() > lenghtGroup[lg]) {
                    group += " > " + lenghtGroup[lg];
                    break;
                }
            }

        } else {
            // and the next decision is based on minimum peptide length
            for (int lg = 0; lg < lenghtGroup.length; lg++) {
                if (pep1.length() > lenghtGroup[lg] && pep2.length() > lenghtGroup[lg]) {
                    group += " > " +  lenghtGroup[lg];
                    break;
                }
            }
        }
        
        String g = FDRGroupNames.get(group+groupExt);
        return g;
        
    }
    
    /**
     * define the fdr-group for this match
     */
    public void setFDRGroup() {
        String sc = null;
        if (specialcase != null)
            sc = RArrayUtils.toString(specialcase, " ");
        fdrGroup = getFDRGroup(peptide1, peptide2, isLinear, isInternal, sc, (isNonCovalent ? "NonCovalent":"") + (validated == null ? "" : validated));
        if (hasPositiveGrouping()) {
            fdrGroup += " " + validated;
            fdrGroup = FDRGroupNames.get(fdrGroup);
        }
    }
    /**
     * define the fdr-group for this match
     */
    public void setFDRGroup(String fdrGroup) {
        this.fdrGroup = FDRGroupNames.get(fdrGroup);
    }

    /**
     * is this actually a single peptide
     * @return the isLinear
     */
    @Override
    public boolean isLinear() {
        return isLinear;
    }

    /**
     * is this actually a single peptide
     * @param isLinear the isLinear to set
     */
    public void setLinear(boolean isLinear) {
        this.isLinear = isLinear;
    }

    /**
     * The fdr assigned to the score of this peptide pair
     * @param fdr 
     */
    public void setFDR(double fdr) {
        m_fdr = fdr;
        for (PSM psm : chargeTopScoresPSM) {
            if (psm != null)
                psm.setFdrPeptidePair(this);
        }

    }
    /**
     * The fdr assigned to the score of this peptide pair
     * @return  
     */
    public double getFDR() {
        return m_fdr;
    }

    /**
     * a link supported by this peptidepair that passed the fdr
     * @param l
     */
    public void setFdrLink(ProteinGroupLink l) {
        this.m_link = l;
        if (l != null) {
            for (PSM psm : chargeTopScoresPSM) {
                if (psm != null)
                    psm.setFdrPeptidePair(this);
            }
        }
    }


    /**
     * a ProteinGroup supported by this peptidepair that passed the fdr
     * @param pg 
     */
    public void setFdrProteinGroup(ProteinGroup pg) {

        if (pg == null) {
            this.fdrProteinGroup1 = null;
            this.fdrProteinGroup2 = null;
            for (PSM psm : chargeTopScoresPSM) {
                if (psm != null)
                    psm.setFdrProteinGroup(null);
            }
        } else {
            if (pg.equals(peptide1.getProteinGroup())) {
                this.fdrProteinGroup1 = pg;
            }

            if (pg.equals(peptide2.getProteinGroup())) {
                this.fdrProteinGroup2 = pg;
            }

            for (PSM psm : chargeTopScoresPSM) {
                if (psm != null)
                    psm.setFdrProteinGroup(pg);
            }
        }
            
    }

    /**
     * the first ProteinGroup supported by this peptidepair that passed the fdr
     * @return 
     */
    public ProteinGroup getFdrProteinGroup1() {
        return this.fdrProteinGroup1;
    }

    /**
     * the second ProteinGroup supported by this peptidepair that passed the fdr
     * @return 
     */
    public ProteinGroup getFdrProteinGroup2() {
        return this.fdrProteinGroup2;
    }

    /**
     * PSM id for the supporting PSM (only top scoring PSM per Charge state)
     * @return 
     */
    public String[] getPSMids() {
        //String[] ret = new String[chargeTopScoresPSM.size()];
        ArrayList<String> ret = new ArrayList<>(chargeTopScoresPSM.length);
        
        int c = 0;
        for (PSM psm : chargeTopScoresPSM) {
            if (psm!=null) {
                ret.add(psm.getPsmID());
            }
        }

        return ret.toArray(new String[ret.size()]);
    }

    /**
     * Unique id for this peptidepair
     * @return the peptidePairID
     */
    public int getPeptidePairID() {
        return peptidePairID;
    }

    /**
     * a link supported by this peptidepair that passed the FDR
     * @return 
     */
    public ProteinGroupLink getFdrLink() {
        return this.m_link;
    }

    /**
     * number of supporting peptides. Meaning 1 (itself)
     * @return 
     */
    public int getPeptidePairCount() {
        return 1;
    }

    /**
     * are all supporting PSMs "special" cases?
     * @return the specialcase
     */
    public boolean hasNegativeGrouping() {
        return specialcase!=null;
    }

    /**
     * are all supporting PSMs "special" cases?
     * @param specialcase 
     */
    public void setNegativeGrouping(boolean specialcase) {
        if (specialcase) {
            this.specialcase = new ArrayList<>();
            this.specialcase.add("Special");
                    
        } else {
            this.specialcase = null;
        }
    }

    /**
     * are all supporting PSMs "special" cases?
     * @param specialcase 
     */
    @Override
    public void setNegativeGrouping(String cause) {
        if (cause == null) {
            this.specialcase = null;
        } else {
            this.specialcase = new ArrayList<>();
            this.specialcase.add(cause);
        }
    }
    
    @Override
    public String getNegativeGrouping() {
        if (this.specialcase == null)
            return null;
        
        return RArrayUtils.toString(this.specialcase, " ");
    }
    
    
    /**
     * get the peptide at the given site
     * @param n
     * @return 
     */
    public Object getSite(int n) {
        return n==0 ? getPeptide1() : getPeptide2();
    }

    /** 
     * how many site do we have
     * <p> <ul><li>2 - cross-linked peptides</li>
     * <li> 1 - linear peptide </li></ul></p>
     * @return 
     */
    public int getSites() {
        if (isLinear)
            return 1;
        else 
            return 2;
    }    

    /**
     * is this a loop linked peptide
     * @return the isLoop
     */
    public boolean isLoop() {
        return isLoop;
    }

    /**
     * is this a loop linked peptide
     * @param isLoop the isLoop to set
     */
    public void setLoop(boolean isLoop) {
        this.isLoop = isLoop;
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
     * @return the psms
     */
    public ArrayList<PSM> getAllPSMs() {
        return psms;
    }

    /**
     * @return the crosslinker
     */
    public String getCrosslinker() {
        return crosslinker;
    }

    @Override
    public ProteinGroup getProteinGroup1() {
        return peptide1.getProteinGroup();
    }

    @Override
    public ProteinGroup getProteinGroup2() {
        return peptide2.getProteinGroup();
    }

    /**
     * indicates that the two peptides where in the same spectrum but non-covalently linked
     * @return the isNonCovalent
     */
    public boolean isNonCovalent() {
        return isNonCovalent;
    }

    /**
     * indicates that the two peptides where in the same spectrum but non-covalently linked
     * @param isNonCovalent the isNonCovalent to set
     */
    public void setNonCovalent(boolean isNonCovalent) {
        this.isNonCovalent = isNonCovalent;
    }
    
    @Override
    public boolean hasPositiveGrouping() {
        return this.validated != null;
    }
    
    @Override
    public void setPositiveGrouping(String av) {
        this.validated = av;
    }

    @Override
    public String getPositiveGrouping() {
        return this.validated;
    }
    
}
