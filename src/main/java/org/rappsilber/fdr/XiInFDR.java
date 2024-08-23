/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.rappsilber.fdr;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import org.rappsilber.fdr.dataimport.Xi2Xi1Config;
import org.rappsilber.utils.IntArrayList;
import org.rappsilber.utils.Version;
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

    HashMap<String,? extends RunConfig> getConfigs();
    
    RunConfig getConfig(String searchid);
    
    Version getXiVersion();
    
    Version getXiVersion(String searchid);
    HashMap<String,Version> getXiVersions();

    public ArrayList<String> getSearchIDs();


    int getFastas(String searchID, ArrayList<String> dbIDs, ArrayList<String> names);
}
