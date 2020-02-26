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
public class PeptideSite extends AbstractSite{
    private Peptide m_peptide;
    private int m_link;

    public PeptideSite(Peptide p, int link) {
        m_peptide = p;
        m_link = link;
    }


    @Override
    public int hashCode() {
        return m_peptide.hashCode(); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public boolean equals(Object obj) {
        if (obj instanceof PeptideSite) {
            PeptideSite e = (PeptideSite) obj;
            return m_peptide.equals(e.m_peptide) && m_link == e.m_link ; //To change body of generated methods, choose Tools | Templates.
        }
        return false;
    }
    
    
}
