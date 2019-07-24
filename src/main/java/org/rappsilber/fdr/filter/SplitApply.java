/*
 * Copyright 2019 Lutz Fischer <lfischer@staffmail.ed.ac.uk>.
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
package org.rappsilber.fdr.filter;

import java.util.ArrayList;
import org.rappsilber.fdr.entities.PSM;

/**
 *
 * @author Lutz Fischer <lfischer@staffmail.ed.ac.uk>
 */
public class SplitApply implements PSMFilter{
    PSMFilter inner;

    public SplitApply(PSMFilter inner) {
        this.inner = inner;
    }
    
    
    @Override
    public ArrayList<PSM> filter(ArrayList<PSM> psm) {
        ArrayList<PSM>  linear = new ArrayList<PSM>(psm.size()/3);
        ArrayList<PSM>  within = new ArrayList<PSM>(psm.size()/3);
        ArrayList<PSM>  between = new ArrayList<PSM>(psm.size()/3);
        
        for (PSM p : psm) {
            if (p.isBetween())
                between.add(p);
            else if (p.isInternal())
                within.add(p);
            else
                linear.add(p);
        }
        linear = inner.filter(linear);
        linear.addAll(inner.filter(within));
        linear.addAll(inner.filter(between));
        return linear;
    }
    
}
