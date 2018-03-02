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
    int BoostingSteps;
    boolean BoostBetween;
    boolean LinkDirectional;
    boolean PPIDirectional;
    boolean PSMDirectional;
    boolean PeptidePairDirectional;
    double ReportFactor;
    boolean filterToUniquePSM;
    protected boolean boostPSM = true;
    protected boolean boostPeptidePairs = true;
    protected boolean boostProteins = true;
    protected boolean boostLinks = true;
    protected boolean groupByPPI = false;
    int minTD = 0;
    
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
        this.setBoostingSteps(settings.getBoostingSteps());
        this.setMaxLinkAmbiguity(settings.getMaxLinkAmbiguity());
        this.setMaxProteinAmbiguity(settings.getMaxProteinAmbiguity());
        this.setMinLinkPepCount(settings.getMinLinkPepCount());
        this.setMinPPIPepCount(settings.getMinPPIPepCount());
        this.setMinPeptideLength(settings.getMinPeptideLength());
        this.setMinProteinPepCount(settings.getMinProteinPepCount());
        this.setPSMFDR(settings.getPSMFDR());
        this.setPSMDirectional(settings.isPSMDirectional());
        this.setPPIDirectional(settings.isPPIDirectional());
        this.setLinkDirectional(settings.isLinkDirectional());
        this.setPeptidePairDirectional(settings.isPeptidePairDirectional());
        this.setPeptidePairFDR(settings.getPeptidePairFDR());
        this.setProteinGroupFDR(settings.getProteinGroupFDR());
        this.setProteinGroupLinkFDR(settings.getProteinGroupLinkFDR());
        this.setProteinGroupPairFDR(settings.getProteinGroupPairFDR());
        this.setReportFactor(settings.getReportFactor());
        this.doOptimize(settings.doOptimize());
        this.setFilterToUniquePSM(settings.filterToUniquePSM());
        this.setBoostBetween(settings.getBoostBetween());
        this.boostLinks(settings.boostLinks());
        this.boostPSMs(settings.boostPSMs());
        this.boostPeptidePairs(settings.boostPeptidePairs());
        this.boostProteins(settings.boostProteins());
    }


    public boolean boostProteins() {
        return boostProteins;
    }

    public void boostProteins(boolean boost) {
        boostProteins = boost;

    }

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
}
