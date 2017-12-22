/*
 * Copyright 2017 Lutz Fischer <lfischer@staffmail.ed.ac.uk>.
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
import java.util.Iterator;
import org.rappsilber.fdr.groups.ProteinGroup;

/**
 *
 * @author Lutz Fischer <lfischer@staffmail.ed.ac.uk>
 */
public class UniquePSM extends PSM implements Iterable<PSM>{
    
    private ArrayList<PSM> support;
    private PSM topPSM;
    
    public UniquePSM(PSM psm) {
        super(psm.getPsmID(), psm.getRun(), psm.getScan(), psm.getPeptide1(), psm.getPeptide2(), psm.getPeptideLinkSite1(), psm.getPeptideLinkSite2(), psm.isDecoy1(), psm.isDecoy2(), psm.getCharge(), psm.getScore(), psm.getScoreRatio());
        support= new ArrayList<>();
        support.add(psm);
        topPSM=psm;
        this.crosslinker = psm.getCrosslinker();
        psm.setPartOfUniquePSM(this);
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
        return this.getCrosslinker() == c.getCrosslinker() && this.getCharge() == c.getCharge() &&
//                (c.getPeptide1().equals(this.getPeptide1()) && c.getPeptide2().equals(this.getPeptide2()) && ((PeptideSite)c.getLinkSite1() == this.getLinkSite1() && c.getLinkSite2() == this.getLinkSite2()) 
                (c.getPeptide1().equals(this.getPeptide1()) && c.getPeptide2().equals(this.getPeptide2()) && c.getPeptideLinkSite1()== this.getPeptideLinkSite1() && c.getPeptideLinkSite2() == this.getPeptideLinkSite2()) 
                || /* to be safe from binary inaccuracy we make an integer-comparison*/
                (c.getPeptide2().equals(this.getPeptide1()) && c.getPeptide1().equals(this.getPeptide2()) && c.getPeptideLinkSite2()== this.getPeptideLinkSite1() && c.getPeptideLinkSite1() == this.getPeptideLinkSite2());
    }
    
    /**
     * adds the information of another PSM to this one.
     * @param p 
     */
    @Override
    public void add(PSM p) {
        if (p == this) {
            return;
        }
        if (p instanceof UniquePSM) {
            for (PSM psm : ((UniquePSM) p).support) {
                support.add(psm);
                psm.setPartOfUniquePSM(this);
            }
            //support.addAll(((UniquePSM)p).support);
            
            if (p.getScore() > getScore()) {
                setScore(p.getScore());
                topPSM=((UniquePSM) p).topPSM;
            }
        } else {
            support.add(p);
            p.setPartOfUniquePSM(this);
            if (p.getScore() > getScore()) {
                setScore(p.getScore());
                topPSM=p;
            }
        }
    }    

    @Override
    public boolean isSpecialcase() {
        return getTopPSM().isSpecialcase(); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public int getFDRGroup() {
        return getTopPSM().getFDRGroup(); //To change body of generated methods, choose Tools | Templates.
    }
    
    public int hashCode() {
        return topPSM.hashCode();
    }    

    /**
     * @return the topPSM
     */
    public PSM getTopPSM() {
        return topPSM;
    }

    /**
     * @param topPSM the topPSM to set
     */
    public void setTopPSM(PSM topPSM) {
        this.topPSM = topPSM;
    }

    @Override
    public Iterator<PSM> iterator() {
        return support.iterator();
    }

    @Override
    public void setFDR(double fdr) {
        super.setFDR(fdr); //To change body of generated methods, choose Tools | Templates.
        topPSM.setFDR(fdr);
    }

    @Override
    public void setFDRGroup() {
        super.setFDRGroup(); //To change body of generated methods, choose Tools | Templates.
        topPSM.setFDRGroup();
    }

    @Override
    public void setFdrPeptidePair(PeptidePair pp) {
        super.setFdrPeptidePair(pp); //To change body of generated methods, choose Tools | Templates.
        topPSM.setFdrPeptidePair(pp);
    }

    @Override
    public void setFdrProteinGroup(ProteinGroup pg) {
        super.setFdrProteinGroup(pg); //To change body of generated methods, choose Tools | Templates.
        topPSM.setFdrProteinGroup(pg);
    }

    
}
