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
package org.rappsilber.fdr.utils;

import org.rappsilber.utils.SelfAdd;

/**
 *
 * @author lfischer
 */
public abstract class AbstractFDRElement<T extends SelfAdd<T>> implements FDRSelfAdd<T>{
    protected double m_higherFDR;
    protected double m_lowerFDR;

    /**
     * @return the higherFDR
     */
    public double getHigherFDR() {
        return m_higherFDR;
    }

    /**
     * @param higherFDR the higherFDR to set
     */
    public void setHigherFDR(double higherFDR) {
        this.m_higherFDR = higherFDR;
    }

    /**
     * @return the lowerFDR
     */
    public double getLowerFDR() {
        return m_lowerFDR;
    }

    /**
     * @param lowerFDR the lowerFDR to set
     */
    public void setLowerFDR(double lowerFDR) {
        this.m_lowerFDR = lowerFDR;
    }
    

    
}
