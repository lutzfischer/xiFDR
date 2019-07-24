/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.rappsilber.fdr;

import java.util.ArrayList;
import java.util.Collection;
import org.rappsilber.utils.IntArrayList;
import rappsilber.config.DBRunConfig;
import rappsilber.config.RunConfig;

/**
 *
 * @author Lutz Fischer <lfischer@staffmail.ed.ac.uk>
 */
public interface XiInFDR {

    /**
     * @return the m_conf
     */
    RunConfig getConfig();

    RunConfig getConfig(int searchid);

    public IntArrayList getSearchIDs();


    int getFastas(int searchID, ArrayList<Integer> dbIDs, ArrayList<String> names);
}
