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

import java.util.HashMap;
import org.rappsilber.utils.IntArrayList;
import org.rappsilber.utils.SelfAdd;

/**
 *
 * @author lfischer
 */
public class LinkSite extends AbstractSite {
    private HashMap<Protein, IntArrayList> m_position;
    int hashCode;

    public LinkSite(HashMap<Protein, IntArrayList> positions) {
        m_position = positions;
        setHashCode();
    }


    public LinkSite(ProteinGroupLink pgl, int site) {
        m_position = site==0?pgl.getPosition1() : pgl.getPosition2();
        setHashCode();
    }
    
    
    private void setHashCode() {
        for (Protein p : m_position.keySet()) {
            hashCode+=p.hashCode();
            for (int i : m_position.get(p)) {
                hashCode+=i;
            }
        }
    }
    protected HashMap<Protein, IntArrayList> getPosition() {
        return m_position;
    }

    @Override
    public boolean equals(Object obj) {
        if (!(obj instanceof LinkSite)) {
            return false;
        }
        LinkSite c = (LinkSite) obj;
        HashMap<Protein, IntArrayList> csites = c.getPosition();
        for (Protein cp : csites.keySet()) {
            IntArrayList s = m_position.get(cp);
            if (s == null) {
                return false;
            }
            IntArrayList cs = csites.get(cp);
            if (cs.size() != s.size() || !cs.containsAll(s)) 
                return false;
        }
        return true;
    }

    @Override
    public int hashCode() {
        return hashCode; //To change body of generated methods, choose Tools | Templates.
    }
    
    
    
}
