/*
 * Copyright 2020 Lutz Fischer <lfischer@staffmail.ed.ac.uk>.
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
package org.rappsilber.fdr.calculation;

import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.rappsilber.fdr.FDRSettings;
import org.rappsilber.fdr.FDRSettingsImpl;
import org.rappsilber.fdr.OfflineFDR;
import org.rappsilber.fdr.filter.PSMFilter;
import org.rappsilber.fdr.filter.SingleSubScoreFilter;
import org.rappsilber.fdr.result.FDRResult;
import org.rappsilber.fdr.utils.MaximisingStatus;
import org.rappsilber.fdr.utils.MaximizeLevelInfo;
import org.rappsilber.fdr.utils.MaximizeLevelInfoInteger;
import org.rappsilber.fdr.utils.MaximizingUpdate;

/**
 *
 * @author Lutz Fischer <lfischer@staffmail.ed.ac.uk>
 */
public class Boost {

    private OfflineFDR fdr;
    private boolean stopMaximizing;
    
    public MaximisingStatus maximise(FDRSettings fdrSettings, OfflineFDR.FDRLevel level, final boolean between, final MaximizingUpdate stateUpdate) {
        return maximise(fdrSettings, level, between, stateUpdate, true);
    }

    public MaximisingStatus maximise(FDRSettings fdrSettings, OfflineFDR.FDRLevel level, final boolean between, final MaximizingUpdate stateUpdate, boolean twoStep) {
        if (twoStep && 
                (fdrSettings.boostDeltaScore() || fdrSettings.boostMinFragments()|| fdrSettings.boostPepCoverage()) &&
                (fdrSettings.boostPSMs() || fdrSettings.boostPeptidePairs() || fdrSettings.boostProteins() || fdrSettings.boostLinks())) {
            FDRSettingsImpl step = new FDRSettingsImpl();
            step.setAll(fdrSettings);
            step.boostDeltaScore(false);
            step.boostMinFragments(false);
            step.boostPepCoverage(false);
            
            MaximisingStatus res = maximiseInner(step, level, between, stateUpdate);
            step.setPSMFDR(res.showPSMFDR);
            step.boostPSMs(false);
            step.setPeptidePairFDR(res.showPepFDR);
            step.boostPeptidePairs(false);
            step.setProteinGroupFDR(res.showProtFDR);
            step.boostProteins(false);
            step.setProteinGroupLinkFDR(res.showLinkFDR);
            step.boostLinks(false);
            step.boostDeltaScore(fdrSettings.boostDeltaScore());
            step.boostMinFragments(fdrSettings.boostMinFragments());
            step.boostPepCoverage(fdrSettings.boostPepCoverage());
            
            return maximiseInner(step, level, between, stateUpdate);
        
        } else {
            return maximiseInner(fdrSettings, level, between, stateUpdate);
        }
    }
    
    
    public MaximisingStatus maximiseInner(FDRSettings fdrSettings, OfflineFDR.FDRLevel level, final boolean between, final MaximizingUpdate stateUpdate) {

        final FDRSettingsImpl settings = new FDRSettingsImpl();
        settings.setAll(fdrSettings);
        int maxCountDown=4;

        try {
            int steps = settings.getBoostingSteps();
            StringBuffer sb = new StringBuffer();
            double maxDelta = 1;
            double maxPeptideCoverage = 1;

            // get some settings, that are constant for all calculations
            boolean ignoreGroups = fdr.isIgnoreGroupsSetting();

            // if minpepcoverage is 0 and we boost on it we can probably savly start ad 0.1
            if (settings.getMinPeptideCoverageFilter() == 0 && settings.boostPepCoverage()) {
//                settings.setMinPeptideCoverageFilter(0.01);
                maxPeptideCoverage = 0.8;
            }

            // also for delta score we can probably savly start as 0.05
            if (settings.boostDeltaScore()) {
//                if (settings.getMinDeltaScoreFilter() == 0) {
//                    settings.setMinDeltaScoreFilter(0.01);
//                }
                // there is probably no point in being more restrictive then 0.7 (even that is overkill)
                if (settings.getMinDeltaScoreFilter() < 0.6) {
                    maxDelta = 0.7;
                }
            }

            final MaximizeLevelInfo absPepCover = new MaximizeLevelInfoInteger(10, settings.boostMinFragments(), steps);
            final MaximizeLevelInfo deltaScore = new MaximizeLevelInfo(1 - settings.getMinDeltaScoreFilter(), 1 - maxDelta, settings.boostDeltaScore(), steps);
            final MaximizeLevelInfo pepCoverage = new MaximizeLevelInfo(1 - settings.getMinPeptideCoverageFilter(), 1 - maxPeptideCoverage, settings.boostPepCoverage(), steps);
            final MaximizeLevelInfo psmFDRInfo = new MaximizeLevelInfo(settings.getPSMFDR(), settings.boostPSMs(), steps);
            final MaximizeLevelInfo pepFDRInfo = new MaximizeLevelInfo(settings.getPeptidePairFDR(), settings.boostPeptidePairs() && level.compareTo(OfflineFDR.FDRLevel.PEPTIDE_PAIR) > 0, steps);
            final MaximizeLevelInfo protFDRInfo = new MaximizeLevelInfo(settings.getProteinGroupFDR(), settings.boostProteins() && level.compareTo(OfflineFDR.FDRLevel.PROTEINGROUP) > 0, steps);
            final MaximizeLevelInfo linkFDRInfo = new MaximizeLevelInfo(settings.getProteinGroupLinkFDR(), settings.boostLinks() && level.compareTo(OfflineFDR.FDRLevel.PROTEINGROUPLINK) > 0, steps);
            final MaximizeLevelInfo ppiFDRInfo = new MaximizeLevelInfo(settings.getProteinGroupPairFDR(), false, steps);
            final ArrayList<MaximizeLevelInfo> subScoreFilter = new ArrayList<>();
            MaximizeLevelInfo targetInfoFirst;

            switch (level) {
                case PEPTIDE_PAIR:
                    targetInfoFirst = pepFDRInfo;
                    break;
                case PROTEINGROUP:
                    targetInfoFirst = protFDRInfo;
                    break;
                case PROTEINGROUPLINK:
                    targetInfoFirst = linkFDRInfo;
                    break;
                case PROTEINGROUPPAIR:
                    targetInfoFirst = ppiFDRInfo;
                    break;
                default:
                    targetInfoFirst = null;
                    break;
            }
            final MaximizeLevelInfo targetInfo = targetInfoFirst;

            boolean filterToUniquePSM = settings.filterToUniquePSM();

            double linkfdr;

            boolean optimizing = true;

            int maxCount = 0;
            int maxCountBetween = 0;

            int countDown = maxCountDown;

            int optimizingRound = 1;

            while (optimizing) {
                int lastMaxCount = maxCount;
                Logger.getLogger(this.getClass().getName()).log(Level.INFO, "Round " + optimizingRound++);

                Logger.getLogger(this.getClass().getName()).log(Level.INFO, "{0}{1}{2}{3}{4}{5}{6}Steps : {7}{8}", new Object[]{
                    deltaScore.boost ? "deltaFilter from  :        " + (1 - deltaScore.fromFDR) + " to " + (1 - deltaScore.toFDR) + "\n" : "", 
                    pepCoverage.boost ? "PepCoverage from  :        " + (1 - pepCoverage.fromFDR) + " to " + (1 - pepCoverage.toFDR) + "\n" : "", 
                    absPepCover.boost ? "Min Pep Frags from  :        " + (10 - absPepCover.fromFDR) + " to " + (10 - absPepCover.toFDR) + "\n" : "", 
                    psmFDRInfo.boost ? "PSM fdr from  :        " + psmFDRInfo.fromFDR + " to " + psmFDRInfo.toFDR + "\n" : "", 
                    pepFDRInfo.boost ? "Peptide pair fdr from  " + pepFDRInfo.fromFDR + " to " + pepFDRInfo.toFDR + "\n" : "", 
                    protFDRInfo.boost ? "Protein-groupfdr from  " + protFDRInfo.fromFDR + " to " + protFDRInfo.toFDR + "\n" : "", 
                    linkFDRInfo.boost ? "linkfdr from  " + linkFDRInfo.fromFDR + " to " + linkFDRInfo.toFDR + "\n" : "", 
                    psmFDRInfo.steps + psmFDRInfo.stepChange, ""});
                FDRResult result = new FDRResult();

                // initialise subscore filter
                if (subScoreFilter.size() > 0) {
                    for (int c = 0; c < subScoreFilter.size(); c++) {
                        subScoreFilter.get(c).firstStep();
                    }
                }
                boolean optimzeSubScores = true;

                // find the combinations with the maximum number of ppis
                deltaScoreloop:
                for (deltaScore.firstStep(); deltaScore.doThisStep(); deltaScore.nextStep()) {
                    settings.setMinDeltaScoreFilter(1 - deltaScore.getCurrentFDR());
                    pepFragsloop:
                    for (absPepCover.firstStep(); absPepCover.doThisStep(); absPepCover.nextStep()) {
                        settings.setMinPeptideFragmentsFilter(10 - (int) Math.round(absPepCover.getCurrentFDR()));
                        pepCoverageloop:
                        for (pepCoverage.firstStep(); pepCoverage.doThisStep(); pepCoverage.nextStep()) {
                            settings.setMinPeptideCoverageFilter(1 - pepCoverage.getCurrentFDR());

//                        double highestPSMFDR = Double.POSITIVE_INFINITY;
//                        boolean firstPSMFDR = true;
                            psmloop:
                            for (psmFDRInfo.firstStep(); psmFDRInfo.doThisStep(); psmFDRInfo.nextStep()) {
//                            if (highestPSMFDR < psmFDRInfo.currentFDR)
//                                continue;

                                settings.setPSMFDR(psmFDRInfo.getCurrentFDR());
                                fdr.calculatePSMFDR(true, fdr.isIgnoreGroupsSetting(), result, settings);
                                psmFDRInfo.setCountsPrefilter(result.psmFDR);
//                            if ( firstPSMFDR) {
//                                firstPSMFDR  =false;
//                                double maxInFDR = 0;
//                                highestPSMFDR = getHighestSubGroupInputFDR(result.psmFDR.getGroups());
//                            }

                                if (stopMaximizing()) {
                                    break deltaScoreloop;
                                }
                                // if we don't get PSM - stop looking at later stages
                                if (result.psmFDR.getResultCount() == 0) {
                                    break psmloop;
                                }

//                            double highestPepFDR = Double.POSITIVE_INFINITY;
//                            boolean firstPepPairFDR = true;
                                peploop:
                                for (pepFDRInfo.firstStep(); pepFDRInfo.doThisStep(); pepFDRInfo.nextStep()) {
                                    // we can skip any steps that would not reduce the peptide FDR below what is already in the input
                                    // calculate peptide level fdr
//                                if (highestPepFDR < pepFDRInfo.currentFDR)
//                                    continue;

                                    //double highestProtFDR = Double.POSITIVE_INFINITY;
                                    //boolean firstProtFDR = true;
                                    protloop:
                                    for (protFDRInfo.firstStep(); protFDRInfo.doThisStep(); protFDRInfo.nextStep()) {
                                        // we can skip any steps that would not reduce the peptide FDR below what is already in the input
                                        // calculate peptide level fdr
//                                    if (highestProtFDR < protFDRInfo.currentFDR)
//                                        continue;

                                        settings.setPeptidePairFDR(pepFDRInfo.getCurrentFDR());
                                        fdr.calculatePeptidePairFDR(true, result, settings, ignoreGroups);

//                                    if ( firstPepPairFDR) {
//                                        firstPepPairFDR  =false;
//                                        double maxInFDR = 0;
//                                        highestPepFDR = getHighestSubGroupInputFDR(result.peptidePairFDR.getGroups());
//                                    }
//                                    
                                        // if we don't get peptide pairs - stop looking at later stages
                                        if (result.peptidePairFDR.getResultCount() == 0) {
                                            break peploop;
                                        }
                                        pepFDRInfo.setCountsPrefilter(result.peptidePairFDR);

                                        if (stopMaximizing()) {
                                            break deltaScoreloop;
                                        }

                                        // calculate protein level fdr
                                        settings.setProteinGroupFDR(protFDRInfo.getCurrentFDR());
                                        fdr.calculateProteinGroupFDR(fdr.isIgnoreGroupsSetting(), true, settings, result);

//                                    if ( firstProtFDR) {
//                                        firstProtFDR  =false;
//                                        double maxInFDR = 0;
//                                        highestProtFDR = getHighestSubGroupInputFDR(result.proteinGroupFDR.getGroups());
//                                    }
                                        if (result.proteinGroupFDR.getResultCount() == 0) {
                                            break protloop;
                                        }
                                        protFDRInfo.setCountsPrefilter(result.proteinGroupFDR);

                                        // cut down the peptides by proteins                           
                                        fdr.filterFDRPeptidePairsByFDRProteinGroups(result);

                                        if (stopMaximizing()) {
                                            break deltaScoreloop;
                                        }

//                                    double highestLinkFDR = Double.POSITIVE_INFINITY;
//                                    boolean firstLinkFDR = true;
                                        linkloop:
                                        for (linkFDRInfo.firstStep(); linkFDRInfo.doThisStep(); linkFDRInfo.nextStep()) {
                                            // we can skip any steps that would not reduce the peptide FDR below what is already in the input
                                            // calculate peptide level fdr
//                                        if (highestLinkFDR < linkFDRInfo.currentFDR)
//                                            continue;

                                            // calculate links
                                            settings.setProteinGroupLinkFDR(linkFDRInfo.getCurrentFDR());
                                            fdr.calculateLinkFDR(ignoreGroups, true, settings, result);
                                            linkFDRInfo.setCountsPrefilter(result.proteinGroupLinkFDR);

//                                        if ( firstLinkFDR && linkFDRInfo.boost) {
//                                            Collection<SubGroupFdrInfo<ProteinGroupLink>> c = result.proteinGroupLinkFDR.getGroups();
//                                            firstLinkFDR  =false;
//                                            // detect the highest fdr that we can see - filtering for anything higher is pointless
//                                            highestLinkFDR = getHighestSubGroupInputFDR(c);
//                                        }
                                            if (result.proteinGroupLinkFDR.getResultCount() == 0) {
                                                break linkloop;
                                            }

                                            if (stopMaximizing()) {
                                                break deltaScoreloop;
                                            }

                                            settings.setProteinGroupPairFDR(ppiFDRInfo.getCurrentFDR());
                                            fdr.calculateProteinGroupPairFDR(ignoreGroups, true, settings, result);

                                            if (result.proteinGroupPairFDR.getResultCount() == 0) {
                                                break linkloop;
                                            }
                                            // now we need to filter down to the required level
                                            //                                if (level.compareTo(level.PROTEINGROUPPAIR)!=0) {
                                            fdr.filterFDRLinksByFDRProteinGroupPairs(result);
                                            //                                }
                                            //                                if (level.compareTo(level.PROTEINGROUPLINK)!=0){
                                            fdr.filterFDRPeptidePairsByFDRProteinGroupLinks(result);
                                            //                                }

                                            //                                if (level.compareTo(level.PROTEINGROUP)==0){
                                            fdr.filterFDRProteinGroupsByFDRPeptidePairs(result);
                                            //                                }

                                            // how many links do we now have?
                                            pepCoverage.setCounts(result.psmFDR);
                                            absPepCover.setCounts(result.psmFDR);
                                            deltaScore.setCounts(result.psmFDR);
                                            psmFDRInfo.setCounts(result.psmFDR);
                                            pepFDRInfo.setCounts(result.peptidePairFDR);
                                            protFDRInfo.setCounts(result.proteinGroupFDR);
                                            linkFDRInfo.setCounts(result.proteinGroupLinkFDR);
                                            ppiFDRInfo.setCounts(result.proteinGroupPairFDR);

                                            int count = targetInfo.count;
                                            int countBetween = targetInfo.countBetween;
                                            if (count == 0 && maxCount <= targetInfo.countPreFilter) {
                                                count = targetInfo.countPreFilter;
                                            }
                                            if (countBetween == 0 && maxCountBetween <= targetInfo.countBetweenPreFilter) {
                                                countBetween = targetInfo.countBetweenPreFilter;
                                            }

                                            // is it a new best?
                                            if (((between && (countBetween > maxCountBetween || (countBetween == maxCountBetween && count > maxCount)))
                                                    || (!between && count > maxCount))) {
                                                maxCount = count;
                                                maxCountBetween = countBetween;

                                                deltaScore.setNewMaxFDR();
                                                absPepCover.setNewMaxFDR();
                                                pepCoverage.setNewMaxFDR();
                                                psmFDRInfo.setNewMaxFDR();
                                                pepFDRInfo.setNewMaxFDR();
                                                protFDRInfo.setNewMaxFDR();
                                                linkFDRInfo.setNewMaxFDR();
                                                protFDRInfo.setNewMaxFDR();

                                                // record that we found a new top
                                                String message = "link count, " + linkFDRInfo.count + "(" + linkFDRInfo.countBetween + " between), Protein Pairs, " + ppiFDRInfo.count + "(" + ppiFDRInfo.countBetween + " between)";
                                                if (linkFDRInfo.boost) {
                                                    message = "link fdr, " + (linkFDRInfo.getCurrentFDR()) + ", " + message;
                                                }
                                                if (protFDRInfo.boost) {
                                                    message = "prot fdr, " + (protFDRInfo.getCurrentFDR()) + ", " + message;
                                                }
                                                if (pepFDRInfo.boost) {
                                                    message = "pep fdr, " + (pepFDRInfo.getCurrentFDR()) + ", " + message;
                                                }
                                                if (psmFDRInfo.boost) {
                                                    message = "psm fdr, " + (psmFDRInfo.getCurrentFDR()) + ", " + message;
                                                }
                                                if (pepCoverage.boost) {
                                                    message = "PepCoverage, " + (1 - pepCoverage.getCurrentFDR()) + ", " + message;
                                                }
                                                if (absPepCover.boost) {
                                                    message = "min Fragmenst, " + (10 - absPepCover.getCurrentFDR()) + ", " + message;
                                                }
                                                if (deltaScore.boost) {
                                                    message = "deltascore, " + (1 - deltaScore.getCurrentFDR()) + ", " + message;
                                                }

                                                // String message = "psmfdr, " + psmFDRInfo.currentFDR + " , pepfdr, " + pepFDRInfo.currentFDR + " ,protfdr, " + protFDRInfo.currentFDR + ", link count, " + linkFDRInfo.count + "(" + linkFDRInfo.countBetween + " between), Protein Pairs, " + ppiFDRInfo.count + "(" + ppiFDRInfo.countBetween + " between)";
                                                sb.append(message + "\n");
                                                Logger.getLogger(this.getClass().getName()).log(Level.INFO, sb.toString());

                                                // forward the values to the gui
                                                final double showDelta = (1 - deltaScore.getCurrentFDR());
                                                final double showPepCoverage = (1 - pepCoverage.getCurrentFDR());
                                                final int showMinFrags = (int) Math.round((10 - absPepCover.getCurrentFDR()));
                                                final double showPSMFDR = psmFDRInfo.getCurrentFDR();
                                                final double showPepFDR = pepFDRInfo.getCurrentFDR();
                                                final double showProtFDR = protFDRInfo.getCurrentFDR();
                                                final double showLinkFDR = linkFDRInfo.getCurrentFDR();
                                                final double showPSMCount = psmFDRInfo.count;
                                                final double showPepCount = pepFDRInfo.count;
                                                final double showProtCount = protFDRInfo.countLinear;
                                                final double showLinkCount = linkFDRInfo.count;
                                                final double showPPICount = ppiFDRInfo.count;
                                                final double showLinkCountBetween = linkFDRInfo.countBetween;
                                                final double showPPICountBetween = ppiFDRInfo.countBetween;
                                                final double showLinkCountBetweenPreFilter = linkFDRInfo.countBetweenPreFilter;
                                                final double showLinkCountPreFilter = linkFDRInfo.countPreFilter;
                                                final MaximisingStatus state = new MaximisingStatus();
                                                state.showDelta = showDelta;
                                                state.showMinFrags = showMinFrags;
                                                state.showPepCoverage = showPepCoverage;
                                                state.showPSMFDR = showPSMFDR;
                                                state.showPepFDR = showPepFDR;
                                                state.showProtFDR = showProtFDR;
                                                state.showLinkFDR = showLinkFDR;

                                                state.showPSMCount = showPSMCount + "";
                                                state.showPepCount = showPepCount + "";
                                                state.showProtCount = showProtCount + "";
                                                state.showLinkCount = showLinkCount + (showLinkCount == 0 ? "(before PPI:" + showLinkCountPreFilter + ")" : "");
                                                state.showLinkCountBetween = showLinkCountBetween + (showLinkCountBetween == 0 ? "(before PPI:" + showLinkCountBetweenPreFilter + ")" : "");

                                                state.showPPICount = showPPICount + "";
                                                state.showPPICountBetween = showPPICountBetween + "";
                                                state.resultCount = maxCount;
                                                state.resultCountBetween = maxCountBetween;
                                                stateUpdate.setStatus(state);
                                            } else if (count == maxCount || (between && countBetween == maxCountBetween)) {
                                                deltaScore.setEqualFDR();
                                                absPepCover.setEqualFDR();
                                                pepCoverage.setEqualFDR();
                                                psmFDRInfo.setEqualFDR();
                                                pepFDRInfo.setEqualFDR();
                                                protFDRInfo.setEqualFDR();
                                                linkFDRInfo.setEqualFDR();
                                                ppiFDRInfo.setEqualFDR();
                                            }

                                            if (stopMaximizing()) {
                                                break deltaScoreloop;
                                            }
                                        }

                                    }

                                }

                            }
                        }
                    }
                }

                stateUpdate.setStatusText("Max Round: " + optimizingRound + " - " + maxCount + " matches");

                // no improvement for the last few rounds?
                if ((maxCount == lastMaxCount && --countDown == 0) || stopMaximizing()) {
                    optimizing = false;
                    Logger.getLogger(this.getClass().getName()).log(Level.INFO, sb.toString());
//                    FDRResult ret = this.calculateFDR(psmFDRInfo.maximumFDR, pepFDRInfo.maximumFDR, protFDRInfo.maximumFDR, linkFDRInfo.maximumFDR, ppiFDRInfo.maximumFDR, 10000, ignoreGroups, true, filterToUniquePSM);
                    settings.setMinDeltaScoreFilter(1 - deltaScore.maximumFDR);
                    settings.setMinPeptideCoverageFilter(1 - pepCoverage.maximumFDR);
                    settings.setMinPeptideFragmentsFilter(10 - (int) Math.round(absPepCover.maximumFDR));
                    settings.setPSMFDR(psmFDRInfo.maximumFDR);
                    settings.setPeptidePairFDR(pepFDRInfo.maximumFDR);
                    settings.setProteinGroupFDR(protFDRInfo.maximumFDR);
                    settings.setProteinGroupLinkFDR(linkFDRInfo.maximumFDR);
                    settings.setProteinGroupPairFDR(ppiFDRInfo.maximumFDR);
                    FDRResult ret = fdr.calculateFDR(settings, true);

                    final int foundCount = maxCount;
                    MaximisingStatus res = new MaximisingStatus();
                    res.showDelta = (1 - deltaScore.maximumFDR);
                    res.showMinFrags = (10 - (int) Math.round(absPepCover.maximumFDR));
                    res.showPepCoverage = (1 - pepCoverage.maximumFDR);
                    res.showPSMFDR = psmFDRInfo.maximumFDR;
                    res.showPepFDR = pepFDRInfo.maximumFDR;
                    res.showProtFDR = protFDRInfo.maximumFDR;
                    res.showLinkFDR = linkFDRInfo.maximumFDR;
                    res.result = ret;
                    res.resultCount = maxCount;
                    res.resultCountBetween = maxCountBetween;

                    if (pepCoverage.maximumFDR < 1) {
                        PSMFilter f = new SingleSubScoreFilter("minPepCoverage", 1 - pepCoverage.maximumFDR, true);
                        fdr.setPrefilteredPSMs(f.filter(fdr.getAllPSMs()));
                        if (fdr.getPrefilteredPSMs().size() < 50) {
                            break;
                        }
                    } else {
                        fdr.setPrefilteredPSMs(null);
                    }

                    stopMaximizing(false);

                    stateUpdate.setStatus(res);

                    return res;

                } else {
                    if (maxCount > lastMaxCount) {
                        // yes we improved
                        countDown = maxCountDown;
                    } else {
                        stateUpdate.setStatusText("Max Link Round: " + optimizingRound + " - " + maxCount + " links - count down: " + countDown);

                    }
                    // so see if we make the resoltuion finer
                    // can we get a better result?
                    boolean startLow = false;
                    int stepChange = 0;
                    if (countDown == 2) {
                        startLow = true;
                        stepChange = 1;
                    }
                    if (countDown == 1) {
                        stepChange = 1;
                    }
                    deltaScore.calcNextFDRRange(startLow, stepChange);
                    absPepCover.calcNextFDRRange(startLow, stepChange);
                    pepCoverage.calcNextFDRRange(startLow, stepChange);
                    psmFDRInfo.calcNextFDRRange(startLow, stepChange);
                    pepFDRInfo.calcNextFDRRange(startLow, stepChange);
                    linkFDRInfo.calcNextFDRRange(startLow, stepChange);
                    protFDRInfo.calcNextFDRRange(startLow, stepChange);
                    ppiFDRInfo.calcNextFDRRange(startLow, stepChange);

                }

            }
        } catch (Exception ex) {
            Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, "Error maximizing links", ex);
            stateUpdate.reportError("Error maximizing", ex);
        }
        return null;

    }

    public MaximisingStatus maximiseSequential(FDRSettings fdrSettings, OfflineFDR.FDRLevel level, final boolean between, final MaximizingUpdate stateUpdate) {

        final FDRSettingsImpl settings = new FDRSettingsImpl();
        settings.setAll(fdrSettings);

        try {
            int steps = settings.getBoostingSteps();
            StringBuffer sb = new StringBuffer();
            double maxDelta = 1;
            double maxPeptideCoverage = 1;

            // get some settings, that are constant for all calculations
            boolean ignoreGroups = fdr.isIgnoreGroupsSetting();

            // if minpepcoverage is 0 and we boost on it we can probably savly start ad 0.1
            if (settings.getMinPeptideCoverageFilter() == 0 && settings.boostPepCoverage()) {
//                settings.setMinPeptideCoverageFilter(0.01);
                maxPeptideCoverage = 0.8;
            }

            // also for delta score we can probably savly start as 0.05
            if (settings.boostDeltaScore()) {
//                if (settings.getMinDeltaScoreFilter() == 0) {
//                    settings.setMinDeltaScoreFilter(0.01);
//                }
                // there is probably no point in being more restrictive then 0.7 (even that is overkill)
                if (settings.getMinDeltaScoreFilter() < 0.6) {
                    maxDelta = 0.7;
                }
            }

            final MaximizeLevelInfo absPepCover = new MaximizeLevelInfoInteger(10, true, steps);
            final MaximizeLevelInfo deltaScore = new MaximizeLevelInfo(1 - settings.getMinDeltaScoreFilter(), 1 - maxDelta, settings.boostDeltaScore(), steps);
            final MaximizeLevelInfo pepCoverage = new MaximizeLevelInfo(1 - settings.getMinPeptideCoverageFilter(), 1 - maxPeptideCoverage, settings.boostPepCoverage(), steps);
            final MaximizeLevelInfo psmFDRInfo = new MaximizeLevelInfo(settings.getPSMFDR(), settings.boostPSMs(), steps);
            final MaximizeLevelInfo pepFDRInfo = new MaximizeLevelInfo(settings.getPeptidePairFDR(), settings.boostPeptidePairs() && level.compareTo(OfflineFDR.FDRLevel.PEPTIDE_PAIR) > 0, steps);
            final MaximizeLevelInfo protFDRInfo = new MaximizeLevelInfo(settings.getProteinGroupFDR(), settings.boostProteins() && level.compareTo(OfflineFDR.FDRLevel.PROTEINGROUP) > 0, steps);
            final MaximizeLevelInfo linkFDRInfo = new MaximizeLevelInfo(settings.getProteinGroupLinkFDR(), settings.boostLinks() && level.compareTo(OfflineFDR.FDRLevel.PROTEINGROUPLINK) > 0, steps);
            final MaximizeLevelInfo ppiFDRInfo = new MaximizeLevelInfo(settings.getProteinGroupPairFDR(), false, steps);
            absPepCover.name = "min frags";
            deltaScore.name = "delta";
            pepCoverage.name = "pep coverage";
            psmFDRInfo.name = "psm";
            pepFDRInfo.name = "pep";
            protFDRInfo.name = "prot";
            linkFDRInfo.name = "link";
            ppiFDRInfo.name = "ppi";

            final ArrayList<MaximizeLevelInfo> subScoreFilter = new ArrayList<>();
            MaximizeLevelInfo targetInfoFirst;

            ArrayList<MaximizeLevelInfo> boostLevels = new ArrayList<>();
            if (settings.boostMinFragments()) {
                boostLevels.add(absPepCover);
            }
            if (settings.boostPepCoverage()) {
                boostLevels.add(pepCoverage);
            }
            if (settings.boostDeltaScore()) {
                boostLevels.add(deltaScore);
            }
            switch (level) {
                case PROTEINGROUPPAIR:
                    if (settings.boostLinks()) {
                        boostLevels.add(linkFDRInfo);
                    }
                case PROTEINGROUPLINK:
                    if (settings.boostProteins()) {
                        boostLevels.add(protFDRInfo);
                    }
                    if (settings.boostPeptidePairs()) {
                        boostLevels.add(pepFDRInfo);
                    }
                case PEPTIDE_PAIR:
                    if (settings.boostPSMs()) {
                        boostLevels.add(psmFDRInfo);
                    }
            }

            switch (level) {
                case PSM:
                    targetInfoFirst = psmFDRInfo;
                    break;
                case PEPTIDE_PAIR:
                    targetInfoFirst = pepFDRInfo;
                    break;
                case PROTEINGROUP:
                    targetInfoFirst = protFDRInfo;
                    break;
                case PROTEINGROUPLINK:
                    targetInfoFirst = linkFDRInfo;
                    break;
                case PROTEINGROUPPAIR:
                    targetInfoFirst = ppiFDRInfo;
                    break;
                default:
                    targetInfoFirst = null;
                    break;
            }
            final MaximizeLevelInfo targetInfo = targetInfoFirst;

            boolean filterToUniquePSM = settings.filterToUniquePSM();

            MaximisingStatus res = new MaximisingStatus();

            int maxCount = 0;
            int maxCountBetween = 0;


            for (int i = 2; i>0;i--) {
                for (MaximizeLevelInfo ml : boostLevels) {
                    ml.fromFDR = ml.minimumFDR;
                    boolean optimizing = true;

                    int countDown = 5;

                    int optimizingRound = 0;
                    if (ml == absPepCover) {
                        Logger.getLogger(this.getClass().getName()).log(Level.INFO,
                                "Min Pep Frags from  :        " + (10 - ml.fromFDR) + " to " + (10 - ml.toFDR));
                    } else if (ml == pepCoverage) {
                        Logger.getLogger(this.getClass().getName()).log(Level.INFO,
                                "Min Pep Coverage from  :        " + (1 - ml.fromFDR) + " to " + (1 - ml.toFDR));
                    } else if (ml == deltaScore) {
                        Logger.getLogger(this.getClass().getName()).log(Level.INFO,
                                "delta score from  :        " + (1 - ml.fromFDR) + "*score to " + (1 - ml.toFDR) + "*score");
                    } else if (ml == psmFDRInfo) {
                        Logger.getLogger(this.getClass().getName()).log(Level.INFO,
                                "psmFDR from  :        " + ml.fromFDR + " to " + ml.toFDR);
                    } else if (ml == pepFDRInfo) {
                        Logger.getLogger(this.getClass().getName()).log(Level.INFO,
                                "pep FDR from  :        " + ml.fromFDR + " to " + ml.toFDR);
                    } else if (ml == protFDRInfo) {
                        Logger.getLogger(this.getClass().getName()).log(Level.INFO,
                                "prot FDR from  :        " + ml.fromFDR + " to " + ml.toFDR);
                    } else if (ml == linkFDRInfo) {
                        Logger.getLogger(this.getClass().getName()).log(Level.INFO,
                                "link FDR from  :        " + ml.fromFDR + " to " + ml.toFDR);
                    }
                    while (optimizing && !stopMaximizing()) {
                        optimizingRound++;
                        int lastMaxCount = maxCount;
                        for (ml.firstStep(); ml.doThisStep(); ml.nextStep()) {
                            FDRResult result = new FDRResult();
                            settings.setMinDeltaScoreFilter(1 - deltaScore.getCurrentFDR());
                            settings.setMinPeptideFragmentsFilter(10 - (int) Math.round(absPepCover.getCurrentFDR()));
                            settings.setMinPeptideCoverageFilter(1 - pepCoverage.getCurrentFDR());
                            settings.setPeptidePairFDR(pepFDRInfo.getCurrentFDR());
                            settings.setProteinGroupFDR(protFDRInfo.getCurrentFDR());
                            settings.setProteinGroupLinkFDR(linkFDRInfo.getCurrentFDR());
                            fdr.calculatePSMFDR(true, ignoreGroups, result, settings);
                            if (result.psmFDR.getResultCount() == 0 || stopMaximizing()) {
                                break;
                            }
                            psmFDRInfo.setCountsPrefilter(result.psmFDR);

                            fdr.calculatePeptidePairFDR(true, result, settings, ignoreGroups);
                            if (result.peptidePairFDR.getResultCount() == 0 || stopMaximizing()) {
                                break;
                            }
                            pepFDRInfo.setCountsPrefilter(result.peptidePairFDR);
                            fdr.calculateProteinGroupFDR(ignoreGroups, true, settings, result);
                            if (result.proteinGroupFDR.getResultCount() == 0 || stopMaximizing()) {
                                break;
                            }
                            protFDRInfo.setCountsPrefilter(result.proteinGroupFDR);
                            // cut down the peptides by proteins                           
                            fdr.filterFDRPeptidePairsByFDRProteinGroups(result);
                            fdr.calculateLinkFDR(ignoreGroups, true, settings, result);
                            linkFDRInfo.setCountsPrefilter(result.proteinGroupLinkFDR);
                            if (result.proteinGroupLinkFDR.getResultCount() == 0 || stopMaximizing()) {
                                break;
                            }
                            settings.setProteinGroupPairFDR(ppiFDRInfo.getCurrentFDR());
                            fdr.calculateProteinGroupPairFDR(ignoreGroups, true, settings, result);

                            if (result.proteinGroupPairFDR.getResultCount() == 0 || stopMaximizing()) {
                                break;
                            }
                            // now we need to filter down 
                            fdr.filterFDRLinksByFDRProteinGroupPairs(result);
                            if (stopMaximizing()) {
                                break;
                            }
                            fdr.filterFDRPeptidePairsByFDRProteinGroupLinks(result);
                            if (stopMaximizing()) {
                                break;
                            }
                            fdr.filterFDRProteinGroupsByFDRPeptidePairs(result);
                            if (stopMaximizing()) {
                                break;
                            }
                            fdr.filterFDRPSMByFDRPeptidePairs(result);

                            // how many links do we now have?
                            pepCoverage.setCounts(result.psmFDR);
                            absPepCover.setCounts(result.psmFDR);
                            deltaScore.setCounts(result.psmFDR);
                            psmFDRInfo.setCounts(result.psmFDR);
                            pepFDRInfo.setCounts(result.peptidePairFDR);
                            protFDRInfo.setCounts(result.proteinGroupFDR);
                            linkFDRInfo.setCounts(result.proteinGroupLinkFDR);
                            ppiFDRInfo.setCounts(result.proteinGroupPairFDR);

                            int count = targetInfo.count;
                            int countBetween = targetInfo.countBetween;
                            if (count == 0 && maxCount <= targetInfo.countPreFilter) {
                                count = targetInfo.countPreFilter;
                            }
                            if (countBetween == 0 && maxCountBetween <= targetInfo.countBetweenPreFilter) {
                                countBetween = targetInfo.countBetweenPreFilter;
                            }

                            // is it a new best?
                            if (((between && (countBetween > maxCountBetween || (countBetween == maxCountBetween && count > maxCount)))
                                    || (!between && count > maxCount))) {
                                maxCount = count;
                                maxCountBetween = countBetween;
                                ml.setNewMaxFDR();
                                String message = " , counts: " + count + ", of wich between : " + countBetween;
                                if (ml == absPepCover) {
                                    message = "Min Pep Frags : " + (10 - ml.getCurrentFDR()) + " , counts: " + count + ", of wich between : " + countBetween;
                                } else if (ml == pepCoverage) {
                                    message = "Min Pep Coverage : " + (1 - ml.getCurrentFDR()) + " , counts: " + count + ", of wich between : " + countBetween;
                                } else if (ml == deltaScore) {
                                    message = "min delta score : " + (1 - ml.getCurrentFDR()) + "*score , counts: " + count + ", of wich between : " + countBetween;
                                } else if (ml == psmFDRInfo) {
                                    message = "psmFDR : " + ml.getCurrentFDR() + "*score , counts: " + count + ", of wich between : " + countBetween;
                                } else if (ml == pepFDRInfo) {
                                    message = "peptide pair FDR : " + ml.getCurrentFDR() + "*score , counts: " + count + ", of wich between : " + countBetween;
                                    Logger.getLogger(this.getClass().getName()).log(Level.INFO,
                                            "pep FDR from  :        " + ml.fromFDR + " to " + ml.toFDR);
                                } else if (ml == protFDRInfo) {
                                    message = "protein FDR : " + ml.getCurrentFDR() + "*score , counts: " + count + ", of wich between : " + countBetween;
                                } else if (ml == linkFDRInfo) {
                                    message = "link FDR : " + ml.getCurrentFDR() + "*score , counts: " + count + ", of wich between : " + countBetween;
                                }
                                Logger.getLogger(this.getClass().getName()).log(Level.INFO, sb.toString());

                                // forward the values to the gui
                                final double showDelta = (1 - deltaScore.getCurrentFDR());
                                final double showPepCoverage = (1 - pepCoverage.getCurrentFDR());
                                final int showMinFrags = (int) Math.round((10 - absPepCover.getCurrentFDR()));
                                final double showPSMFDR = psmFDRInfo.getCurrentFDR();
                                final double showPepFDR = pepFDRInfo.getCurrentFDR();
                                final double showProtFDR = protFDRInfo.getCurrentFDR();
                                final double showLinkFDR = linkFDRInfo.getCurrentFDR();
                                final double showPSMCount = psmFDRInfo.count;
                                final double showPepCount = pepFDRInfo.count;
                                final double showProtCount = protFDRInfo.countLinear;
                                final double showLinkCount = linkFDRInfo.count;
                                final double showPPICount = ppiFDRInfo.count;
                                final double showLinkCountBetween = linkFDRInfo.countBetween;
                                final double showPPICountBetween = ppiFDRInfo.countBetween;
                                final double showLinkCountBetweenPreFilter = linkFDRInfo.countBetweenPreFilter;
                                final double showLinkCountPreFilter = linkFDRInfo.countPreFilter;
                                final MaximisingStatus state = new MaximisingStatus();
                                state.showDelta = showDelta;
                                state.showMinFrags = showMinFrags;
                                state.showPepCoverage = showPepCoverage;
                                state.showPSMFDR = showPSMFDR;
                                state.showPepFDR = showPepFDR;
                                state.showProtFDR = showProtFDR;
                                state.showLinkFDR = showLinkFDR;

                                state.showPSMCount = showPSMCount + "";
                                state.showPepCount = showPepCount + "";
                                state.showProtCount = showProtCount + "";
                                state.showLinkCount = showLinkCount + (showLinkCount == 0 ? "(before PPI:" + showLinkCountPreFilter + ")" : "");
                                state.showLinkCountBetween = showLinkCountBetween + (showLinkCountBetween == 0 ? "(before PPI:" + showLinkCountBetweenPreFilter + ")" : "");

                                state.showPPICount = showPPICount + "";
                                state.showPPICountBetween = showPPICountBetween + "";
                                state.resultCount = maxCount;
                                state.resultCountBetween = maxCountBetween;
                                stateUpdate.setStatus(state);
                            } else if (count == maxCount || (between && countBetween == maxCountBetween)) {
                                ml.setEqualFDR();
                            }

                        }

                        stateUpdate.setStatusText("Max Round: " + optimizingRound + " ("+ ml.name+") - " + maxCount + " matches");

                        // no improvement for the last few rounds?
                        if ((maxCount == lastMaxCount && --countDown == 0) || stopMaximizing()) {
                            optimizing = false;
                            Logger.getLogger(this.getClass().getName()).log(Level.INFO, sb.toString());
                            //                    FDRResult ret = this.calculateFDR(psmFDRInfo.maximumFDR, pepFDRInfo.maximumFDR, protFDRInfo.maximumFDR, linkFDRInfo.maximumFDR, ppiFDRInfo.maximumFDR, 10000, ignoreGroups, true, filterToUniquePSM);
                            settings.setMinDeltaScoreFilter(1 - deltaScore.maximumFDR);
                            settings.setMinPeptideCoverageFilter(1 - pepCoverage.maximumFDR);
                            settings.setMinPeptideFragmentsFilter(10 - (int) Math.round(absPepCover.maximumFDR));
                            settings.setPSMFDR(psmFDRInfo.maximumFDR);
                            settings.setPeptidePairFDR(pepFDRInfo.maximumFDR);
                            settings.setProteinGroupFDR(protFDRInfo.maximumFDR);
                            settings.setProteinGroupLinkFDR(linkFDRInfo.maximumFDR);
                            //settings.setProteinGroupPairFDR(settings.getProteinGroupPairFDR());
                            res.showDelta = (1 - deltaScore.maximumFDR);
                            res.showMinFrags = (10 - (int) Math.round(absPepCover.maximumFDR));
                            res.showPepCoverage = (1 - pepCoverage.maximumFDR);
                            res.showPSMFDR = psmFDRInfo.maximumFDR;
                            res.showPepFDR = pepFDRInfo.maximumFDR;
                            res.showProtFDR = protFDRInfo.maximumFDR;
                            res.showLinkFDR = linkFDRInfo.maximumFDR;
                            res.resultCount = maxCount;
                            res.resultCountBetween = maxCountBetween;

                            if (pepCoverage.maximumFDR < 1) {
                                PSMFilter f = new SingleSubScoreFilter("minPepCoverage", 1 - pepCoverage.maximumFDR, true);
                                fdr.setPrefilteredPSMs(f.filter(fdr.getAllPSMs()));
                                if (fdr.getPrefilteredPSMs().size() < 50) {
                                    break;
                                }
                            } else {
                                fdr.setPrefilteredPSMs(null);
                            }

                            stateUpdate.setStatus(res);
                            optimizing = false;
                        } else {
                            if (maxCount > lastMaxCount) {
                                // yes we improved
                                countDown = 5;
                            } else {
                                stateUpdate.setStatusText("Max Link Round: " + optimizingRound + " - " + maxCount + " links - count down: " + countDown);

                            }
                            // so see if we make the resoltuion finer
                            // can we get a better result?
                            boolean startLow = false;
                            int stepChange = 0;
                            if (countDown == 2) {
                                startLow = true;
                                stepChange = 1;
                            }
                            if (countDown == 1) {
                                stepChange = 1;
                            }
                            ml.calcNextFDRRange(startLow, stepChange);

                        }

                    }
                    ml.setCurrentFDR(ml.maximumFDR);

                }
            }
            FDRResult ret = fdr.calculateFDR(settings, true);
            res.result = ret;

            final int foundCount = maxCount;
            
            stopMaximizing(false);
            return res;
        } catch (Exception ex) {
            Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, "Error maximizing links", ex);
            stateUpdate.reportError("Error maximizing", ex);
        }
        return null;

    }

    /**
     * @return the stopMaximizing
     */
    public boolean stopMaximizing() {
        return stopMaximizing;
    }

    /**
     * @param stopMaximizing the stopMaximizing to set
     */
    public void stopMaximizing(boolean stopMaximizing) {
        this.stopMaximizing = stopMaximizing;
    }
    
}
