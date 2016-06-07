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
package org.rappsilber.fdr.result;

import java.util.Collection;
import org.rappsilber.fdr.entities.PSM;
import org.rappsilber.fdr.entities.PeptidePair;
import org.rappsilber.fdr.entities.ProteinGroupLink;
import org.rappsilber.fdr.entities.ProteinGroupPair;
import org.rappsilber.fdr.groups.ProteinGroup;

/**
 *
 * @author lfischer
 */
public  class FDRResult {
    public Collection<PSM> input;
    public FDRResultLevel<PSM> psmFDR;   
    public FDRResultLevel<PeptidePair> peptidePairFDR;   
    public FDRResultLevel<ProteinGroup> proteinGroupFDR;   
    public FDRResultLevel<ProteinGroupLink> proteinGroupLinkFDR;   
    public FDRResultLevel<ProteinGroupPair> proteinGroupPairFDR;   

    public int minPeptideLength = 0;
    public int maximumProteinAmbiguity = 0;
    public int maximumLinkAmbiguity = 0;
    public boolean uniquePSMs = false;
    


}
