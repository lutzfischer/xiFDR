/*
 * Copyright 2016 lfischer.
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


public class FDRSettingsImpl implements FDRSettings {

    OfflineFDR.FDRLevel doOptimize;
    int MaxLinkAmbiguity;
    int MaxProteinAmbiguity;
    int MinLinkPepCount;
    int MinPPIPepCount;
    int MinPeptideLength;
    int MinProteinPepCount;
    double PSMFDR;
    double PeptidePairFDR;
    double ProteinGroupFDR;
    double ProteinGroupLinkFDR;
    double ProteinGroupPairFDR;
    int BoostingSteps = 4;
    boolean BoostBetween;
    boolean LinkDirectional;
    boolean PPIDirectional;
    boolean PSMDirectional;
    boolean PeptidePairDirectional;
    double ReportFactor = 1000000;
    boolean filterToUniquePSM;
    protected boolean boostPSM = true;
    protected boolean boostPeptidePairs = true;
    protected boolean boostProteins = true;
    protected boolean boostLinks = true;
//    private boolean boostSubScores = false;
    protected boolean groupByPPI = false;
    int minTD = DEFAULT_MIN_TD_COUNT;
    private boolean filterConsectutivePeptides;
    private double subScoreCutOff = 1;
    private boolean boostPepCoverage = true;
    private Boolean psmLocalFDR;
    private Boolean peppairLocalFDR;
    private Boolean protLocalFDR;
    private Boolean linkLocalFDR;
    private Boolean ppiLocalFDR;
    private double minPeptideCoverageFilter;
    private double minDeltaScoreFilter;
    private boolean boostPepDeltaScore =true;
    private boolean combineScoreAndDelta;
    private int minFragments;
    private boolean boostMinFragments = false;
    private boolean ignoreValidityChecks = true;
    private boolean groubByCrosslinkerStubs;
    
    @Override
    public void addCalcListener(ActionListener al) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public OfflineFDR.FDRLevel doOptimize() {
        return doOptimize;
    }

    @Override
    public void doOptimize(OfflineFDR.FDRLevel level) {
        doOptimize = level;
    }

    @Override
    public int getMaxLinkAmbiguity() {
        return MaxLinkAmbiguity;
    }

    @Override
    public int getMaxProteinAmbiguity() {
        return MaxProteinAmbiguity;
    }

    @Override
    public int getMinLinkPepCount() {
        return MinLinkPepCount;
    }

    @Override
    public int getMinPPIPepCount() {
        return MinPPIPepCount;
    }

    @Override
    public int getMinPeptideLength() {
        return MinPeptideLength;
    }

    @Override
    public int getMinProteinPepCount() {
        return MinProteinPepCount;
    }

    @Override
    public double getPSMFDR() {
        return PSMFDR; 
    }

    @Override
    public double getPeptidePairFDR() {
        return PeptidePairFDR;
    }

    @Override
    public double getProteinGroupFDR() {
        return ProteinGroupFDR;
    }

    @Override
    public double getProteinGroupLinkFDR() {
        return ProteinGroupLinkFDR;
    }

    @Override
    public double getProteinGroupPairFDR() {
        return ProteinGroupPairFDR;
    }

    @Override
    public void setPSMFDR(Double fdr) {
        PSMFDR = fdr;
    }

    @Override
    public void setPeptidePairFDR(Double fdr) {
        PeptidePairFDR=fdr;
    }

    @Override
    public void setProteinGroupFDR(Double fdr) {
        ProteinGroupFDR=fdr;
    }

    @Override
    public void setProteinGroupLinkFDR(Double fdr) {
        ProteinGroupLinkFDR=fdr;
    }

    @Override
    public int getBoostingSteps() {
        return BoostingSteps;
    }

    @Override
    public void setBoostingSteps(int steps) {
        BoostingSteps=steps;
    }

    @Override
    public boolean getBoostBetween() {
        return BoostBetween;
    }

    @Override
    public void setBoostBetween(boolean between) {
        BoostBetween=between;
    }

    @Override
    public boolean isLinkDirectional() {
        return LinkDirectional;
    }

    @Override
    public boolean isPPIDirectional() {
        return PPIDirectional;
    }

    @Override
    public boolean isPSMDirectional() {
        return PSMDirectional;
    }

    @Override
    public boolean isPeptidePairDirectional() {
        return PeptidePairDirectional;
    }

    @Override
    public void removeCalcListener(ActionListener al) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void setLinkDirectional(boolean directional) {
        LinkDirectional = directional;
    }

    @Override
    public void setMaxLinkAmbiguity(Integer maxAmbiguity) {
        MaxLinkAmbiguity = maxAmbiguity;
    }

    @Override
    public void setMaxProteinAmbiguity(Integer maxAmbiguity) {
        MaxProteinAmbiguity = maxAmbiguity;
    }

    @Override
    public void setMinLinkPepCount(Integer minPep) {
        MinLinkPepCount = minPep;
    }

    @Override
    public void setMinPPIPepCount(Integer minPep) {
        MinPPIPepCount = minPep;
    }

    @Override
    public void setMinPeptideLength(Integer minLength) {
        MinPeptideLength = minLength;
    }

    @Override
    public void setMinProteinPepCount(Integer minPep) {
        MinProteinPepCount = minPep;
    }

    @Override
    public void setPPIDirectional(boolean directional) {
        PPIDirectional=directional;
    }

    @Override
    public void setPSMDirectional(boolean directional) {
        PSMDirectional = directional;
    }

    @Override
    public void setPeptidePairDirectional(boolean directional) {
        PeptidePairDirectional = directional;
    }

    @Override
    public void setProteinGroupPairFDR(Double fdr) {
        ProteinGroupPairFDR = fdr;
    }

    @Override
    public double getReportFactor() {
        return ReportFactor;
    }

    @Override
    public void setReportFactor(double factor) {
        ReportFactor =factor;
    }

    @Override
    public boolean filterToUniquePSM() {
        return filterToUniquePSM;
    }

    @Override
    public void setFilterToUniquePSM(boolean filterToUnique) {
        filterToUniquePSM=filterToUnique;
    }

    @Override
    public void setAll(FDRSettings settings) {
        transferSettings(settings, this);
    }
    
    public static void transferSettings(FDRSettings from, FDRSettings to) {
        to.setBoostingSteps(from.getBoostingSteps());
        to.setMaxLinkAmbiguity(from.getMaxLinkAmbiguity());
        to.setMaxProteinAmbiguity(from.getMaxProteinAmbiguity());
        to.setMinLinkPepCount(from.getMinLinkPepCount());
        to.setMinPPIPepCount(from.getMinPPIPepCount());
        to.setMinPeptideLength(from.getMinPeptideLength());
        to.setMinProteinPepCount(from.getMinProteinPepCount());
        to.setMinTD(from.getMinTD());
        to.setPSMFDR(from.getPSMFDR());
        to.setPSMDirectional(from.isPSMDirectional());
        to.setPPIDirectional(from.isPPIDirectional());
        to.setLinkDirectional(from.isLinkDirectional());
        to.setPeptidePairDirectional(from.isPeptidePairDirectional());
        to.setPeptidePairFDR(from.getPeptidePairFDR());
        to.setProteinGroupFDR(from.getProteinGroupFDR());
        to.setProteinGroupLinkFDR(from.getProteinGroupLinkFDR());
        to.setProteinGroupPairFDR(from.getProteinGroupPairFDR());
        to.setReportFactor(from.getReportFactor());
        to.doOptimize(from.doOptimize());
        to.setFilterToUniquePSM(from.filterToUniquePSM());
        to.setBoostBetween(from.getBoostBetween());
        to.boostLinks(from.boostLinks());
        to.boostPSMs(from.boostPSMs());
        to.boostPeptidePairs(from.boostPeptidePairs());
        to.boostProteins(from.boostProteins());
        to.boostMinFragments(from.boostMinFragments());
//        to.boostSubScores(from.boostSubScores());
//        to.setSubScoreCutOff(from.getSubScoreCutOff());
        to.setFilterConsecutivePeptides(from.filterConsecutivePeptides());
        to.psmLocalFDR(from.psmLocalFDR());    
        to.peppairLocalFDR(from.peppairLocalFDR());    
        to.protLocalFDR(from.protLocalFDR());    
        to.linkLocalFDR(from.linkLocalFDR());    
        to.ppiLocalFDR(from.ppiLocalFDR());    
        to.setMinDeltaScoreFilter(from.getMinDeltaScoreFilter());    
        to.setMinPeptideFragmentsFilter(from.getMinPeptideFragmentsFilter()); 
        to.setMinPeptideCoverageFilter(from.getMinPeptideCoverageFilter()); 
        to.boostPepCoverage(from.boostPepCoverage());
        to.boostDeltaScore(from.boostDeltaScore());
        to.combineScoreAndDelta(from.combineScoreAndDelta());
        to.ignoreValidityChecks(from.ignoreValidityChecks());
        to.setGroupByCrosslinkerStubs(from.getGroupByCrosslinkerStubs());
    }

    public boolean combineScoreAndDelta() {
        return this.combineScoreAndDelta;
    }
    public void combineScoreAndDelta(boolean c) {
        this.combineScoreAndDelta = c;
    }


    public boolean boostProteins() {
        return boostProteins;
    }

    public void boostProteins(boolean boost) {
        boostProteins = boost;

    }

    
//    public boolean boostSubScores() {
//        return boostSubScores;
//    }
//
//    public void boostSubScores(boolean boost) {
//        boostSubScores = boost;
//    }
    
    public boolean boostPSMs() {
        return boostPSM;
    }

    public void boostPSMs(boolean boost) {
        boostPSM = boost;

    }

    public boolean boostPeptidePairs() {
        return boostPeptidePairs;
    }

    public void boostPeptidePairs(boolean boost) {
        boostPeptidePairs = boost;
    }

    public boolean boostLinks() {
        return boostLinks;
    }

    public void boostLinks(boolean boost) {
        boostLinks = boost;
    }    

    @Override
    public boolean isGroupByPSMCount() {
        return groupByPPI;
    }

    @Override
    public void setGroupByPSMCount(boolean groupByPPI) {
        this.groupByPPI=groupByPPI;
    }

    @Override
    public void setMinTD(Integer c) {
        minTD =c;
    }

    @Override
    public int getMinTD() {
        return minTD;
    }

    @Override
    public boolean filterConsecutivePeptides() {
        return filterConsectutivePeptides;
    }

    @Override
    public void setFilterConsecutivePeptides(boolean filterConsecutive) {
        this.filterConsectutivePeptides=filterConsecutive;
    }
    
    public double getSubScoreCutOff() {
        return this.subScoreCutOff;
    }
    public void setSubScoreCutOff(double localfdr) {
        this.subScoreCutOff = localfdr;
    }
    
    @Override
    public boolean boostDeltaScore(){
        return this.boostPepDeltaScore;
    }

    @Override
    public void boostDeltaScore(boolean boost){
        this.boostPepDeltaScore = boost;
    }
    
    @Override
    public boolean boostPepCoverage(){
        return this.boostPepCoverage;
    }

    @Override
    public void boostPepCoverage(boolean boost){
        this.boostPepCoverage = boost;
    }

    @Override
    public boolean boostMinFragments(){
        return this.boostMinFragments;
    }

    @Override
    public void boostMinFragments(boolean boost){
        this.boostMinFragments = boost;
    }

    @Override
    public Boolean psmLocalFDR() {
        return this.psmLocalFDR;
    }

    @Override
    public Boolean peppairLocalFDR() {
        return this.peppairLocalFDR;
    }

    @Override
    public Boolean protLocalFDR() {
        return this.protLocalFDR;
    }

    @Override
    public Boolean linkLocalFDR() {
        return this.linkLocalFDR;
    }

    @Override
    public Boolean ppiLocalFDR() {
        return this.ppiLocalFDR;
    }

    @Override
    public void psmLocalFDR(Boolean local) {
       this.psmLocalFDR = local;
    }

    @Override
    public void peppairLocalFDR(Boolean local) {
        this.peppairLocalFDR = local;
    }

    @Override
    public void protLocalFDR(Boolean local) {
        this.protLocalFDR = local;
    }

    @Override
    public void linkLocalFDR(Boolean local) {
        this.linkLocalFDR = local;
    }

    @Override
    public void ppiLocalFDR(Boolean local) {
        this.ppiLocalFDR = local;
    }
    
    public double getMinPeptideCoverageFilter() {
        return this.minPeptideCoverageFilter;
    }
    public void setMinPeptideCoverageFilter(double d) {
        this.minPeptideCoverageFilter = d;
    }

    public double getMinDeltaScoreFilter() {
        return this.minDeltaScoreFilter;
    }
    
    public void setMinDeltaScoreFilter(double d) {
        this.minDeltaScoreFilter = d;
    }

    void setMinAbsolutePepCoverFilter(double d) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public int getMinPeptideFragmentsFilter() {
        return this.minFragments;
    }

    @Override
    public void setMinPeptideFragmentsFilter(int frags) {
        this.minFragments = frags;
    }

    @Override
    public boolean ignoreValidityChecks() {
        return this.ignoreValidityChecks;
    }

    @Override
    public void ignoreValidityChecks(boolean ignore) {
        this.ignoreValidityChecks = ignore;
    }

    @Override
    public void setGroupByCrosslinkerStubs(boolean group) {
        this.groubByCrosslinkerStubs = group;
    }

    @Override
    public boolean getGroupByCrosslinkerStubs() {
        return this.groubByCrosslinkerStubs;
    }
    
    
}
