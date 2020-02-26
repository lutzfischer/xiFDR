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

import org.rappsilber.utils.MapUtils;

/**
 *
 * @author lfischer
 */
public class ProteinGroupDirectionalLink extends ProteinGroupLink {

    public ProteinGroupDirectionalLink(PeptidePair pp) {
        super(pp);
        
    }

    /**
     * a link group is considered equal, if it contains all the same links
     * @param lg
     * @return
     */
    @Override
    public boolean equals(Object o) {
        // if it is actaually the same object
        if (this == o)
            return true;
        
        if (o instanceof ProteinGroupDirectionalLink) {
            
            ProteinGroupDirectionalLink pgl = (ProteinGroupDirectionalLink) o;

            // is it a complete internal link? - meaning a link within the same protein-group (as opposed links between groups that contain a comon protein)
            if (getProteinGroup1().equals(pgl.getProteinGroup1()) && getProteinGroup2().equals(pgl.getProteinGroup2()))
                return (MapUtils.sameMapCollection(getPosition1(), pgl.getPosition1()) && MapUtils.sameMapCollection(getPosition2(), pgl.getPosition2()));
            return false;
        }
        return false;

    }
    
    
}
