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
package org.rappsilber.data.csv;

import java.util.ArrayList;
import java.util.HashMap;


/**
 * Provides a hashed access for a joined condition.
 * <br/>Usage:<pre>{@code
 * CSVConditionLookUp csvl = new CSVConditionLookUp(csv,new String{"Run","Scan"});
 * ArrayList<Integer> rows = getRows(new String{"E123456.raw","123"});
 * }</pre>
 * 
 * It works, by setting up a nested HashMap for the values in the defined columns
 * @author lfischer
 */
public class CSVConditionLookUp {
    private int[] m_columns;
    HashMap<String, Object> m_toplevel;
    CSVRandomAccess m_csv;


    /**
     * Sets up a new hashed lookup for the given CSV
     * @param csv 
     * @param columns 
     */
    public CSVConditionLookUp(CSVRandomAccess csv, int ... columns) {
        m_csv=csv;
        m_columns = columns;
        build();
    }    
    
    /**
     * Returned if nothing is matched
     */
    private final static ArrayList<Integer> NORESULT = new ArrayList<Integer>(){
                @Override
                public boolean add(Integer e) {
                    throw new UnsupportedOperationException("This is an always empty ArrayList");
                }
                
                @Override
                public boolean isEmpty() {
                    return true;
                }
                
            };

    
    /**
     * Sub-routine used to buildup the hashmaps
     * @param parent
     * @param value
     * @return 
     */
    protected HashMap<String, Object> getAddSubMap(HashMap<String, Object> parent, String value) {
        HashMap<String, Object> m = (HashMap<String, Object>) parent.get(value);
        if (m==null) {
            m = new HashMap<String, Object>();
            parent.put(value, m);
        }
        return m;
    }

    /**
     * sets up a nested hash-map for the defined columns
     */
    public void build() {
        m_toplevel = new HashMap<String, Object>();
        for (int r = m_csv.getRowCount()-1; r>= 0 ; r--) {
            HashMap<String, Object> level = m_toplevel;
            int c = 0;
            for (; c< m_columns.length-1; c ++) {
                level = getAddSubMap(level, m_csv.getValue(m_columns[c], r));
            }
            ArrayList<Integer> lines = (ArrayList<Integer>) level.get(m_csv.getValue(m_columns[c], r));
            if ( lines == null) {
                lines = new ArrayList<Integer>();
                level.put(m_csv.getValue(m_columns[c], r),lines);
            }
            lines.add(r);
        }
    }
    
    /**
     * get a list of rows, that have the given values
     * @param values
     * @return 
     */
    public ArrayList<Integer> getRows(String ... values) {
        HashMap<String, Object> level = m_toplevel;
        int c = 0;
        for (; c< values.length-1; c ++) {
            level = getAddSubMap(level, values[c]);
        }
        ArrayList<Integer> lines = (ArrayList<Integer>) level.get(values[c]);
        if (lines == null)
            return NORESULT;
        return lines;
    }
}
