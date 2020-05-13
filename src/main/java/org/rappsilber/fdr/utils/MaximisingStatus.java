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

import org.rappsilber.fdr.result.FDRResult;

/**
 *
 * maximising 
 * @author Lutz Fischer <lfischer@staffmail.ed.ac.uk>
 */
public class MaximisingStatus {
        public double showDelta;
        public double showPepCoverage;
        public int    showMinFrags;
        public int    showMinStubs;
        public int    showMinDoublets;
        public double showPSMFDR;
        public double showPepFDR;
        public double showProtFDR;
        public double showLinkFDR;
        public double showProtPairFDR;

        public String showPSMCount;
        public String showPepCount;
        public String showProtCount;
        public String showLinkCount;
        public String showPPICount;

        public String showLinkCountBetween;
        public String showPPICountBetween;
        
        public long resultCount;
        public long resultCountBetween;
        
        public FDRResult result;
    
}
