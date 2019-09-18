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
public class MaximizeLevelInfoInteger extends MaximizeLevelInfo{

    
    public MaximizeLevelInfoInteger(double maximumValue, boolean boost, int steps) {
        super(maximumValue, boost, steps);
        correctStepWidth();
    }

    protected void correctStepWidth() {
        stepWidth = Math.round(stepWidth);
        if (stepWidth < 1) {
            stepWidth = 1;
        }
    }
 
    public void calcNextFDRRange(boolean startLow, int stepChange) {
        super.calcNextFDRRange(startLow, stepChange);
        correctStepWidth();
    }

    @Override
    public double getCurrentFDR() {
        return Math.round(super.getCurrentFDR()); //To change body of generated methods, choose Tools | Templates.
    }
    
    public void setNewMaxFDR() {
        maximumFDR = getCurrentFDR();
        smallestEqualFDR = maximumFDR;
        largestEqualFDR = maximumFDR;
    }
    

}
