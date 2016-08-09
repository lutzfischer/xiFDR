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
import org.rappsilber.utils.SelfAdd;
import java.util.HashSet;
import org.rappsilber.fdr.groups.ProteinGroup;
import org.rappsilber.fdr.utils.AbstractFDRElement;
import org.rappsilber.fdr.utils.FDRSelfAdd;
import org.rappsilber.utils.IntArrayList;

/**
 *
 * @author lfischer
 */
public class ProteinGroupPair extends AbstractFDRElement<ProteinGroupPair> { //implements Comparable<ProteinGroupPair> {

    /**
     * counts all instances of this class.
     * used for assigning unique IDs
     */
    public static int PROTEINGROUPPAIRCOUNT = 0;
    
    
    protected int proteinGroupPairID = PROTEINGROUPPAIRCOUNT++;
    private HashSet<ProteinGroupLink> links = new HashSet<ProteinGroupLink>();
    private ProteinGroup protein1;
    private ProteinGroup protein2;
    private double score;
    int hashcode;
    /** first peptide comes from a decoy sequence */
    private boolean isDecoy1;
    /** second peptide comes from a decoy sequence */
    private boolean isDecoy2;
    /** it could be a protein internal link */
    private boolean isInternal = false;
    private boolean isTT = false;
    private boolean isTD = false;
    private boolean isDD = false;
    private int fdrGroup = 0;
    public double m_fdr = -1;
    private boolean m_specialOnly = true;

    
    public ProteinGroupPair(ProteinGroup Prot1, ProteinGroup Prot2, double score, boolean isSpecialOnly) {
        this.protein1 = Prot1;
        this.protein2 = Prot2;
        
        hashcode = (Prot1.hashCode() + Prot2.hashCode()) % 10000;
        
        this.isDecoy1 = !Prot1.isTT();;
        this.isDecoy2 = !Prot2.isTT();;
        this.isInternal = Prot1.hasOverlap(Prot2);
        this.score = score;
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
        
        fdrGroup = 0;
        if (Prot1 == Prot2 || Prot1.hasOverlap(Prot2))
            fdrGroup = 1;
        m_specialOnly = isSpecialOnly;
        if (m_specialOnly)
            fdrGroup +=2;
                
    }
    
    
    public ProteinGroupPair(PeptidePair pp) {
        this(pp.getPeptide1().getProteinGroup(), pp.getPeptide2().getProteinGroup(), pp.getScore(), pp.isSpecialcase());
    }

    public ProteinGroupPair(ProteinGroupLink l) {
        this(l.getProteinGroup1(), l.getProteinGroup2(),  l.getScore(), l.isSpecialOnly());
        this.links.add(l);
    }

    @Override
    public Site getLinkSite1() {
        return new ProteinGroupSite(protein1);
    }

    @Override
    public Site getLinkSite2() {
        return new ProteinGroupSite(protein2);
    }

    @Override
    public int hashCode() {
        return hashcode;
    }

    @Override
    public boolean equals(Object l) {
        if (l ==  this)
            return true;
        ProteinGroupPair c = (ProteinGroupPair) l;
        return (c.protein1.equals(protein1) && c.protein2.equals(protein2)) || (c.protein1.equals(protein2) && c.protein2.equals(protein1));
    }

    public int compareTo(ProteinGroupPair o) {
        return Double.compare(o.getScore(), this.getScore());
    }

    public void add(ProteinGroupPair pp) {
        if (pp == this)
            return;
        this.links.addAll(pp.links);
        this.score = Math.sqrt(this.score * this.score + pp.score * pp.score);
        m_specialOnly &= pp.m_specialOnly;
        if (!m_specialOnly) {
            fdrGroup = fdrGroup & 1;
        }
        
    }

    /**
     * @return the links
     */
    public HashSet<ProteinGroupLink> getLinks() {
        return links;
    }

    /**
     * @return the protein1
     */
    public ProteinGroup getProtein1() {
        return protein1;
    }

    /**
     * @return the protein2
     */
    public ProteinGroup getProtein2() {
        return protein2;
    }

    /**
     * @return the score
     */
    public double getScore() {
        return score;
    }

//    /**
//     * @return the fdrGroup
//     */
//    public int getFdrGroup() {
//        return fdrGroup;
//    }


    /**
     * @return the isDecoy1
     */
    public boolean isDecoy1() {
        return isDecoy1;
    }

    /**
     * @return the isDecoy2
     */
    public boolean isDecoy2() {
        return isDecoy2;
    }

    /**
     * @return the isInternal
     */
    public boolean isInternal() {
        return isInternal;
    }

    /**
     * @return the isTT
     */
    public boolean isTT() {
        return isTT;
    }

    /**
     * @return the isTD
     */
    public boolean isTD() {
        return isTD;
    }

    /**
     * @return the isDD
     */
    public boolean isDD() {
        return isDD;
    }

    public int getFDRGroup() {
        return fdrGroup;
    }

    public String getFDRGroupName() {
        return getFDRGroupName(fdrGroup);

    }

    public static String getFDRGroupName(int fdrGroup) {
        switch (fdrGroup) {
            case -1 : return "all combined";
            case 0 : return "Between";
            case 1 : return "Within";
            case 2 : return "Special Between";
            case 3 : return "Special Within";
            default : return "Unknown";
        }
    }
    
    
    public Collection<PeptidePair> getPeptidePairs() {
        HashSet<PeptidePair> ret = new HashSet<PeptidePair>();
        for (ProteinGroupLink l : links)
            ret.addAll(l.getPeptidePairs());
        return ret;
    }

    public int[] getPeptidePairID() {
        IntArrayList retal = new IntArrayList();
        for (ProteinGroupLink l : links)
            retal.addAll(l.getPeptidePairIDs());
        return retal.toIntArray();
    }
    
    public Collection<ProteinGroupPair> getProteinPairs() {
        ArrayList<ProteinGroupPair> ret = new ArrayList<ProteinGroupPair>();
        ret.add(this);
        return ret;
        
    }

    public Collection<Protein> getProteins() {
        ArrayList<Protein> ret = new ArrayList<Protein>(2);
        ret.addAll(getProtein1().getProteins());
        ret.addAll(getProtein2().getProteins());
        return ret;
    }

    
    
    public String toString() {
        if (isTT()) 
            return "TT - " + getProtein1().toString() + " xl " + getProtein2().toString();
        if (isTD()) 
            return "TD - " + getProtein1().toString() + " xl " + getProtein2().toString();
        return "DD - " + getProtein1().toString() + " xl " + getProtein2().toString();
    }
    
    public boolean isDecoy() {
        return isDecoy1 || isDecoy2;
    }
 
    
    public void setFDR(double fdr) {
        m_fdr = fdr;
        // flag up Links
        for (ProteinGroupLink pgl : links) {
            pgl.setFdrPPI(this);
        }
    }

    public double getFDR() {
        return m_fdr;
    }

    
    /**
     * @return the peptidePairID
     */
    public int getProteinGroupPairID() {
        return proteinGroupPairID;
    }
    
    /**
     * @return the peptidePairID
     */
    public int[] getLinkIDs() {
        int[] ret = new int[links.size()];
        int c =0;
        for (ProteinGroupLink l : links) {
            ret[c++] = l.getLinkID();
        }
        return ret;
    }

    /**
     * @return the peptidePairID
     */
    public int[] getPeptidePairIDs() {
        IntArrayList ret = new IntArrayList();
        for (ProteinGroupLink l : links) {
            ret.addAll(l.getPeptidePairIDs());
        }
        return ret.toIntArray();
    }
    
    public String[] getPSMIDs() {

        ArrayList<String> ret = new ArrayList<String>();
        for (ProteinGroupLink l : links) {
            for (String id : l.getPSMIDs())
                ret.add(id);
        }
        String[] reta = new String[ret.size()];
        return ret.toArray(reta);
        
    }

    public int getPeptidePairCount() {
        int ret =0;
        for (ProteinGroupLink l : links)
            ret+=l.getPeptidePairCount();
        return ret;
    }

    public Object getSite(int n) {
        return n==0 ? getProtein1() : getProtein2();
    }

    public int getSites() {
        return 2;
    }    

    public boolean isLinear() {
        return false;
    }

    public boolean isBetween() {
        return !isInternal;
    }
    
}
