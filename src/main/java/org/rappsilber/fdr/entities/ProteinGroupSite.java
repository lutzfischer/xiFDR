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

/**
 *
 * @author lfischer
 */
public class ProteinGroupSite extends AbstractSite {
    private ProteinGroup m_proteinGroup;

    public ProteinGroupSite(ProteinGroup p) {
        m_proteinGroup = p;
    }


    @Override
    public int hashCode() {
        return m_proteinGroup.hashCode(); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public boolean equals(Object obj) {
        if (obj instanceof ProteinGroupSite)
            return m_proteinGroup.equals(((ProteinGroupSite)obj).m_proteinGroup); //To change body of generated methods, choose Tools | Templates.
        return false;
    }
    
}
