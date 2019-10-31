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

import java.awt.event.ActionListener;

/**
 *
 * @author lfischer
 */
public interface FDRSettings {
    public static int DEFAULT_MIN_TD_COUNT=2;
    
    void addCalcListener(ActionListener al);

    OfflineFDR.FDRLevel doOptimize();
    void doOptimize(OfflineFDR.FDRLevel level);

    int getMaxLinkAmbiguity();


    int getMaxProteinAmbiguity();

    int getMinLinkPepCount();

    int getMinPPIPepCount();

    int getMinPeptideLength();

    int getMinProteinPepCount();

    double getPSMFDR();

    double getPeptidePairFDR();

    double getProteinGroupFDR();

    double getProteinGroupLinkFDR();

    double getProteinGroupPairFDR();

    void setPSMFDR(Double fdr);
    
    void setPeptidePairFDR(Double fdr);

    void setProteinGroupFDR(Double fdr);

    void setProteinGroupLinkFDR(Double fdr);
    
    
    int getBoostingSteps();
    public void setBoostingSteps(int steps);

    boolean getBoostBetween();
    void setBoostBetween(boolean between);
    
    boolean isLinkDirectional();

    boolean isPPIDirectional();

    boolean isPSMDirectional();

    boolean isPeptidePairDirectional();

    void removeCalcListener(ActionListener al);

    void setLinkDirectional(boolean directional);

    void setMaxLinkAmbiguity(Integer maxAmbiguity);


    void setMaxProteinAmbiguity(Integer maxAmbiguity);

    void setMinLinkPepCount(Integer minPep);

    void setMinPPIPepCount(Integer minPep);

    void setMinPeptideLength(Integer minLength);

    void setMinProteinPepCount(Integer minPep);

    void setPPIDirectional(boolean directional);

    void setPSMDirectional(boolean directional);

    void setPeptidePairDirectional(boolean directional);


    void setProteinGroupPairFDR(Double fdr);
    
    double getReportFactor();
    
    public void setReportFactor(double factor);       
    
    boolean filterToUniquePSM();
    
    void setFilterToUniquePSM(boolean filterToUnique);
    
    boolean filterConsecutivePeptides();
    
    void setFilterConsecutivePeptides(boolean filterConsecutive);

    
    public void setAll(FDRSettings settings);
    
    public boolean boostProteins();
    public void boostProteins(boolean boost);

    public boolean boostPepCoverage();    
    public void boostPepCoverage(boolean boost);    

    public boolean boostMinFragments();
    public void boostMinFragments(boolean boost);
    
    public boolean boostDeltaScore();
    public void boostDeltaScore(boolean boost);

    public boolean boostPSMs();
    public void boostPSMs(boolean boost);
//    public boolean boostSubScores();
//    public void boostSubScores(boolean boost);
    public boolean boostPeptidePairs();
    public void boostPeptidePairs(boolean boost);
    public boolean boostLinks();

    public void boostLinks(boolean boost);
    public boolean isGroupByPSMCount();
    public void setGroupByPSMCount(boolean groupByPPI);
    
    public double getMinPeptideCoverageFilter();
    public void setMinPeptideCoverageFilter(double d);
    
    public int getMinPeptideFragmentsFilter();
    public void setMinPeptideFragmentsFilter(int frags);

    public double getMinDeltaScoreFilter();
    public void setMinDeltaScoreFilter(double d);
    
    public void setMinTD(Integer c) ;

    public int getMinTD();

    public Boolean psmLocalFDR();
    public Boolean peppairLocalFDR();
    public Boolean protLocalFDR();
    public Boolean linkLocalFDR();
    public Boolean ppiLocalFDR();
    public void psmLocalFDR(Boolean local);
    public void peppairLocalFDR(Boolean local);
    public void protLocalFDR(Boolean local);
    public void linkLocalFDR(Boolean local);
    public void ppiLocalFDR(Boolean local);
    public boolean combineScoreAndDelta();
    public void combineScoreAndDelta(boolean c);
    public boolean ignoreValidityChecks();
    public void ignoreValidityChecks(boolean ignore);
//    public double getSubScoreCutOff();
//    public void setSubScoreCutOff(double localfdr);
    public void setGroupByCrosslinkerStubs(boolean group);
    public boolean getGroupByCrosslinkerStubs();
    
}
