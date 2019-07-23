/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.rappsilber.fdr.utils;

import org.rappsilber.fdr.entities.DBPSM;
import org.rappsilber.fdr.entities.PSM;

/**
 *
 * @author Lutz Fischer <lfischer@staffmail.ed.ac.uk>
 */
public interface PSMToMzIdentScanId {
    String getID(DBPSM psm);
}
