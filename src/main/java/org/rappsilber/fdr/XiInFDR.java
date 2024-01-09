/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.rappsilber.fdr;

import java.util.ArrayList;
import java.util.Collection;
import org.rappsilber.fdr.dataimport.Xi2Config;
import org.rappsilber.utils.IntArrayList;
import rappsilber.config.DBRunConfig;
import rappsilber.config.RunConfig;

/**
 *
 * @author Lutz Fischer <lfischer@staffmail.ed.ac.uk>
 */
public interface XiInFDR {
    int getxiMajorVersion();
    
    /**
     * @return the m_conf
     */
    RunConfig getConfig();

    RunConfig getConfig(String searchid);
    
    Xi2Config getXi2Config();
    
    Xi2Config getXi2Config(String searchid);

    public ArrayList<String> getSearchIDs();


    int getFastas(String searchID, ArrayList<String> dbIDs, ArrayList<String> names);
}
