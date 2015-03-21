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

/**
 * represents a directional peptide pair and stores the support for the given 
 * peptide pair.
 * <p>The only difference to the {@link PeptidePair} is that the equal function 
 * only returns true if the peptide in the peptide pair are in the same order.</p>
 * <p> by declaring only the same order as equal different orders will be summed 
 * up to two separate entries.</p>
 */
public class DirectionalPeptidePair extends PeptidePair{

    /**
     * Constructor - just forwards to {@link PeptidePair}
     * @param psm 
     */
    public DirectionalPeptidePair(PSM psm) {
        super(psm);
    }
    
    /**
     * Returns true if both peptides are the same - including the same order - 
     * and the link sites are the same.
     * @param l
     * @return 
     */
    @Override
    public boolean equals(Object l) {
        DirectionalPeptidePair c = (DirectionalPeptidePair) l;
        return (c.getPeptide1().equals(getPeptide1()) && c.getPeptide2().equals(getPeptide2()) && c.getPeptideLinkSite1() == c.getPeptideLinkSite1() && c.getPeptideLinkSite2() == c.getPeptideLinkSite2());
    }    

}
