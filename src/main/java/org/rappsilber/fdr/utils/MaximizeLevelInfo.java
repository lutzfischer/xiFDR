/*
 * Copyright 2018 Lutz Fischer <lfischer@staffmail.ed.ac.uk>.
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
package org.rappsilber.fdr.utils;

import org.rappsilber.fdr.result.FDRResultLevel;

/**
 * used during maximising results
 * @author Lutz Fischer <lfischer@staffmail.ed.ac.uk>
 */
public class MaximizeLevelInfo {
    public String name;
    public double firstFDR;
    public double maximumFDR = 0;
    public double minimumFDR = 0;
    public double fromFDR;
    public double toFDR;
    public double steps;
    public double stepWidth;
    public double smallestEqualFDR = 0;
    public double largestEqualFDR = 0;
    public boolean boost;
    private double currentFDR;

    public int count = 0;
    public int countBetween = 0;
    public int countLinear = 0;

    public int countPreFilter = 0;
    public int countBetweenPreFilter = 0;
    public int countLinearPreFilter = 0;

    public MaximizeLevelInfo(double maximumFDR, boolean boost, int steps) {
        this(maximumFDR, 0, boost, steps);
    }

    public MaximizeLevelInfo(double maximumFDR, double minimumFDR, boolean boost, int steps) {
        this.maximumFDR = maximumFDR;
        this.firstFDR = maximumFDR;
        this.currentFDR = maximumFDR;
        this.steps = steps;
        this.toFDR = maximumFDR;
        this.boost = boost;
        this.minimumFDR = minimumFDR;
        if (boost) {
            double min = Math.max(minimumFDR, 0.0005);
            fromFDR = Math.min(min, toFDR / steps);
            stepWidth = (toFDR - fromFDR) / steps;
        } else {
            fromFDR = toFDR;
            stepWidth = 100;
        }
    }

    public void firstStep() {
        currentFDR = toFDR;
    }

    public boolean doThisStep() {
        return currentFDR >= fromFDR;
    }

    public void nextStep() {
        currentFDR -= stepWidth;
    }

    public void setCounts(FDRResultLevel l) {
        count = l.getBetween() + l.getWithin();
        countLinear = l.getLinear();
        countBetween = l.getBetween();
    }

    public void setCountsPrefilter(FDRResultLevel l) {
        countPreFilter = l.getBetween() + l.getWithin();
        countBetweenPreFilter = l.getBetween();
        countLinearPreFilter = l.getLinear();
    }

    public void setNewMaxFDR() {
        maximumFDR = currentFDR;
        smallestEqualFDR = currentFDR;
        largestEqualFDR = currentFDR;
    }

    public void setEqualFDR() {
        smallestEqualFDR = Math.min(smallestEqualFDR, currentFDR);
        largestEqualFDR = Math.max(largestEqualFDR, currentFDR);
    }

    public void calcNextFDRRange(boolean startLow, int stepChange) {
        if (boost) {
            toFDR = Math.min(largestEqualFDR + (stepWidth * 3 / 4), firstFDR);
            
            double newFrom = Math.max(smallestEqualFDR - (stepWidth * 3 / 4), 0);
            
            if (newFrom == fromFDR) {
                if (smallestEqualFDR > 0) {
                    newFrom=smallestEqualFDR-Math.random()*(smallestEqualFDR);
                } else {
                    newFrom=Math.random() * stepWidth;
                }
            }
            if (startLow)  {
                newFrom = Math.min(newFrom, Math.max(minimumFDR,0.001));
            }
            fromFDR = newFrom;
            steps += stepChange;
            stepWidth = (toFDR - fromFDR) / steps;
        }
    }

    /**
     * @return the currentFDR
     */
    public double getCurrentFDR() {
        return currentFDR;
    }

    /**
     * @return the currentFDR
     */
    public void setCurrentFDR(double fdr) {
        currentFDR = fdr;
    }

}
