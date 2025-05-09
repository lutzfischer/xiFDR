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

import java.util.Collection;
import java.util.HashSet;
import org.rappsilber.fdr.entities.PeptidePair;
import org.rappsilber.fdr.entities.Site;
import org.rappsilber.fdr.entities.ProteinGroup;
import org.rappsilber.utils.SelfAdd;

/**
 *
 * @author lfischer
 */
public abstract class AbstractFDRElement<T extends SelfAdd<T>> implements FDRSelfAdd<T>{


    protected double m_higherFDR;
    protected double m_lowerFDR;
    /**
     * The next Target-Decoy match with a lower FDR
     */
    private AbstractFDRElement<T> m_lowerFDRTD;
    /**
     * The next Target-Decoy match with a higher FDR
     */
    private AbstractFDRElement<T> m_HigherFDRTD;
    
    protected double m_linkedSupport = 1;
    protected HashSet<String> m_negativeGroups;
    protected HashSet<String> m_positiveGroups;
    private HashSet<String> m_additionalGroups = new HashSet<>();
    private double m_pep;
    
    public abstract Site getLinkSite1();
    public abstract Site getLinkSite2();
    
    /**
     * @return the higherFDR
     */
    public double getHigherFDR() {
        return m_higherFDR;
    }

    /**
     * @param higherFDR the higherFDR to set
     */
    public void setHigherFDR(double higherFDR) {
        this.m_higherFDR = higherFDR;
    }

    /**
     * @return the m_linkedSupport
     */
    public double getLinkedSupport() {
        return m_linkedSupport;
    }

    /**
     * @param linkedSupport the linkedSupport to set
     */
    public void setLinkedSupport(double linkedSupport) {
        this.m_linkedSupport = linkedSupport;
    }

    /**
     * @return the lowerFDR
     */
    public double getLowerFDR() {
        return m_lowerFDR;
    }

    /**
     * @param lowerFDR the lowerFDR to set
     */
    public void setLowerFDR(double lowerFDR) {
        this.m_lowerFDR = lowerFDR;
    }
    
    public abstract ProteinGroup getProteinGroup1();
    public abstract ProteinGroup getProteinGroup2();

    public boolean hasPositiveGrouping() {
        return this.m_positiveGroups!=null;
    }
    
    public void setPositiveGrouping(String av) {
        if (av == null) {
            this.m_positiveGroups = null;
        }else {
            this.m_positiveGroups = new HashSet<String>(1);
            this.m_positiveGroups.add(av);
        }
    }

    public HashSet<String> getPositiveGrouping() {
        return m_positiveGroups;
    }

    public void addPositiveGrouping(String av) {
        if (this.m_positiveGroups == null)
            this.m_positiveGroups = new HashSet<String>(1);
        this.m_positiveGroups.add(av);
    }


    /** 
     * indicates this passed some form of Negative value that makes this 
     * inherently less likely to be true
     */
    public boolean hasNegativeGrouping() {
        return m_negativeGroups != null;
    }
    
    /**
     * are all supporting PSMs "special" cases?
     * @param specialcase 
     */
    public void setNegativeGrouping(String cause) {
        if (cause == null) {
            this.m_negativeGroups = null;
        } else {
            this.m_negativeGroups = new HashSet<>(1);
            this.m_negativeGroups.add(cause);
        }
    }
    
    public HashSet<String> getNegativeGrouping() {
        return this.m_negativeGroups;
    }
 
    public void addNegativeGrouping(String cause) {
        if (this.m_negativeGroups == null) {
            this.m_negativeGroups = new HashSet<>(1);
        }
        this.m_negativeGroups.add(cause);
    }
    
    
    public boolean addFDRGroups(AbstractFDRElement e) {
        boolean changed = false;
        if (this.m_negativeGroups != null) {
            if (e.m_negativeGroups == null) {
                this.m_negativeGroups = null;
            } else {
                m_negativeGroups.retainAll(e.m_negativeGroups);
                changed = true;
            }
        }
        if (e.hasPositiveGrouping()) {
            if (this.m_positiveGroups == null) {
                this.m_positiveGroups = e.getPositiveGrouping();
                changed = true;
            } else if (!this.m_positiveGroups.containsAll(e.getPositiveGrouping())) {
                this.m_positiveGroups.addAll(e.getPositiveGrouping());
                changed = true;
            }
        }
        return changed;
    }

    public abstract Collection<PeptidePair> getPeptidePairs();

    /**
     * @return the m_additionalGroups
     */
    public HashSet<String> getAdditionalFDRGroups() {
        return m_additionalGroups;
    }
    
    
    /**
     * The next Target-Decoy match with a lower FDR
     * @return the m_lowerTD
     */
    public AbstractFDRElement<T> getLowerFDRTD() {
        return m_lowerFDRTD;
    }

    /**
     * The next Target-Decoy match with a lower FDR
     * @param m_lowerTD the m_lowerTD to set
     */
    public void setLowerFDRTD(AbstractFDRElement<T> TD) {
        this.m_lowerFDRTD = TD;
    }

    /**
     * The next Target-Decoy match with a higher FDR
     * @return the m_higherTD
     */
    public AbstractFDRElement<T> getHigherFDRTD() {
        return m_HigherFDRTD;
    }

    /**
     * The next Target-Decoy match with a higher FDR
     * @param m_higherTD the m_higherTD to set
     */
    public void setHigherFDRTD(AbstractFDRElement<T> TD) {
        this.m_HigherFDRTD = TD;
    }    
    
    
    public double getEstimatedFDR() {
        if (this.getLowerFDRTD() == null || this.getHigherFDRTD() == null) {
            return Double.NaN;
        }
        double lowerFDRScore = this.getLowerFDRTD().getScore();
        double higherFDRScore = this.getHigherFDRTD().getScore();
        double lowerFDR = this.getLowerFDRTD().getFDR();
        double higherFDR = Math.max(this.getHigherFDRTD().getFDR(), lowerFDR);
        double stepScore = lowerFDRScore - higherFDRScore;
        double myStepScore = lowerFDRScore - getScore();
        if (stepScore == 0) {
            return lowerFDR + (higherFDR-lowerFDR)/2;
        } else {
            double stepScoreRatio = myStepScore/stepScore;
            double stepFDR = higherFDR - lowerFDR;
            double myStepFDR = stepFDR*stepScoreRatio;
            return lowerFDR + myStepFDR;
        }
        
    }

    public double getOriginalScore() {
        return getScore();
    }

    @Override
    public void setPEP(double pep) {
        this.m_pep = pep;
    }

    @Override
    public Double getPEP() {
        return this.m_pep;
    }
    
    
    
}
