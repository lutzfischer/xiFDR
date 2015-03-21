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

import org.rappsilber.fdr.groups.ProteinGroup;

/**
 *
 * @author lfischer
 */
public class ProteinGroupDirectionalPair extends ProteinGroupPair {
    
    public ProteinGroupDirectionalPair(ProteinGroup Prot1, ProteinGroup Prot2, double score, boolean isSpecialOnly) {
        super(Prot1, Prot2, score, isSpecialOnly);
    }

    public ProteinGroupDirectionalPair(PeptidePair dpp) {
        super(dpp);
    }

    public ProteinGroupDirectionalPair(ProteinGroupLink l) {
        super(l);
    }
    
    
}
