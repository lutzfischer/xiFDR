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
import java.util.Collection;
import org.rappsilber.fdr.entities.PSM;

/**
 *
 * @author Lutz Fischer <lfischer@staffmail.ed.ac.uk>
 */
public class SingleSubScoreFilter  implements PSMFilter{
    /**
     * what subscore to filter on
     */
    String subscore;
    /**
     * what is the value above or below we discard PSMs
     */
    double splitValue;
    /**
     * keep PSMs with value greater (true) or smaller (false) then splitValue
     */
    boolean keepAbove;

    public SingleSubScoreFilter(String subscore, double splitValue, boolean keepAbove) {
        this.subscore = subscore;
        this.splitValue = splitValue;
        this.keepAbove = keepAbove;
    }

    @Override
    public ArrayList<PSM> filter(Collection<PSM> psm) {
        ArrayList<PSM> ret = new ArrayList<>();
        if (keepAbove)
            for (PSM p : psm) {
                if (p.getDoubleInfo(subscore)>= splitValue) {
                    ret.add(p);
                }
            }
        else
            for (PSM p : psm) {
                if ((Double)p.getDoubleInfo(subscore)<= splitValue) {
                    ret.add(p);
                }
            }
        return ret;
    }
    
    
    
    
}
