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

import org.rappsilber.fdr.groups.ProteinGroup;

/**
 *
 * @author lfischer
 */
public class ProteinSite extends AbstractSite {
    private Protein m_protein;

    public ProteinSite(Protein p) {
        m_protein = p;
    }


    @Override
    public int hashCode() {
        return m_protein.hashCode(); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public boolean equals(Object obj) {
        if (obj instanceof ProteinSite)
            return m_protein.equals(((ProteinSite)obj).m_protein); //To change body of generated methods, choose Tools | Templates.
        return false;
    }
    
}
