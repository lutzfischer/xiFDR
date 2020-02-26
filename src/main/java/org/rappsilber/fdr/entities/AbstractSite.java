/*
 * Copyright 2015 lfischer.
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

import org.rappsilber.utils.SelfAdd;

/**
 *
 * @author lfischer
 */
public abstract class AbstractSite implements Site {
    protected double m_connetcedness = 1;

    @Override
    public void add(Site o) {
        setConnectedness(getConnectedness() + o.getConnectedness());
    }

    @Override
    public double getConnectedness() {
        return m_connetcedness;
    }

    @Override
    public void setConnectedness(double connectedness) {
        m_connetcedness = connectedness;
    }
    
}
