/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.rappsilber.fdr.entities;

import org.rappsilber.fdr.utils.FDRGroupNames;

/**
 *
 * @author lfischer
 */
public class DBPSM extends PSM{

    protected int searchID;
    protected int scanID;
    private boolean isAutoValidated;
    private boolean hasVarMods;
    private boolean hasFixedMods;

    /**
     * creates a new instance of a PSM.
     * @param psmID a unique ID for the PSM (e.g. used when a CSV-file contains 
     * several rows for a single PSM that represent each possible origin of that 
     * PSM).
     * @param run run-name
     * @param scan scan-id
     * @param peptide1 first peptide
     * @param peptide2 second peptide
     * @param site1 link-site in the first peptide
     * @param site2 link site in the second peptide
     * @param isDecoy1 is the first peptide a decoy?
     * @param isDecoy2 is the second peptide a decoy?
     * @param charge precursor charge state
     * @param score score for this PSM
     * @param scoreRatio how to split the score between the peptides - this is 
     * only used, if a protein FDR (each single protein) is to be calculated.
     */
    public DBPSM(String psmID, String run, String scan, Peptide peptide1, Peptide peptide2, byte site1, byte site2, boolean isDecoy1, boolean isDecoy2, byte charge, double score, double peptide1score, double peptide2score) {
        super(psmID, run, scan, peptide1, peptide2, site1, site2, isDecoy1, isDecoy2, charge, score, peptide1score,peptide2score);
    }
    /**
     * creates a new instance of a PSM.
     * @param psmID a unique ID for the PSM (e.g. used when a CSV-file contains 
     * several rows for a single PSM that represent each possible origin of that 
     * PSM).
     * @param run run-name
     * @param scan scan-id
     * @param peptide1 first peptide
     * @param peptide2 second peptide
     * @param site1 link-site in the first peptide
     * @param site2 link site in the second peptide
     * @param isDecoy1 is the first peptide a decoy?
     * @param isDecoy2 is the second peptide a decoy?
     * @param charge precursor charge state
     * @param score score for this PSM
     * @param scoreRatio how to split the score between the peptides - this is 
     * only used, if a protein FDR (each single protein) is to be calculated.
     */
    public DBPSM(String psmID, String run, String scan, Peptide peptide1, Peptide peptide2, byte site1, byte site2, boolean isDecoy1, boolean isDecoy2, byte charge, double score, double peptide1score, double peptide2score, int searchID) {
        this(psmID, run, scan, peptide1, peptide2, site1, site2, isDecoy1, isDecoy2, charge, score, peptide1score, peptide1score);
        this.searchID = searchID;
    }
    
    /**
     * creates a new instance of a PSM.
     * @param psmID a unique ID for the PSM (e.g. used when a CSV-file contains 
     * several rows for a single PSM that represent each possible origin of that 
     * PSM).
     * @param peptide1 first peptide
     * @param peptide2 second peptide
     * @param site1 link-site in the first peptide
     * @param site2 link site in the second peptide
     * @param isDecoy1 is the first peptide a decoy?
     * @param isDecoy2 is the second peptide a decoy?
     * @param charge precursor charge state
     * @param score score for this PSM
     * @param scoreRatio how to split the score between the peptides - this is 
     * only used, if a protein FDR (each single protein) is to be calculated.
     */
    public DBPSM(String psmID, Peptide peptide1, Peptide peptide2, byte site1, byte site2, boolean isDecoy1, boolean isDecoy2, byte charge, double score, double peptide1score, double peptide2score) {
        super(psmID, peptide1, peptide2, site1, site2, isDecoy1, isDecoy2, charge, score, peptide1score, peptide2score);
    }
    
    /**
     * creates a new instance of a PSM.
     * @param psmID a unique ID for the PSM (e.g. used when a CSV-file contains 
     * several rows for a single PSM that represent each possible origin of that 
     * PSM).
     * @param peptide1 first peptide
     * @param peptide2 second peptide
     * @param site1 link-site in the first peptide
     * @param site2 link site in the second peptide
     * @param isDecoy1 is the first peptide a decoy?
     * @param isDecoy2 is the second peptide a decoy?
     * @param charge precursor charge state
     * @param score score for this PSM
     * @param scoreRatio how to split the score between the peptides - this is 
     * only used, if a protein FDR (each single protein) is to be calculated.
     */
    public DBPSM(String psmID, Peptide peptide1, Peptide peptide2, byte site1, byte site2, boolean isDecoy1, boolean isDecoy2, byte charge, double score, double peptide1score, double peptide2score, int searchID) {
        this(psmID, peptide1, peptide2, site1, site2, isDecoy1, isDecoy2, charge, score, peptide1score, peptide2score);
        this.searchID = searchID;
    }

    
    /**
     * @return the searchID
     */
    public int getSearchID() {
        return searchID;
    }

    /**
     * @param searchID the searchID to set
     */
    public void setSearchID(int searchID) {
        this.searchID = searchID;
    }

    /**
     * @return the scanID
     */
    public int getScanID() {
        return scanID;
    }

    /**
     * @param scanID the scanID to set
     */
    public void setScanID(int scanID) {
        this.scanID = scanID;
    }

    
    /**
     * @return the hasVarMods
     */
    public boolean hasVarMods() {
        return hasVarMods;
    }

    /**
     * @param hasVarMods the hasVarMods to set
     */
    public void setHasVarMods(boolean hasVarMods) {
        this.hasVarMods = hasVarMods;
    }

    /**
     * @return the hasFixedMods
     */
    public boolean hasFixedMods() {
        return hasFixedMods;
    }

    /**
     * @param hasFixedMods the hasFixedMods to set
     */
    public void setHasFixedMods(boolean hasFixedMods) {
        this.hasFixedMods = hasFixedMods;
    }

    /**
     * @return the isAutoValidated
     */
    public boolean isAutoValidated() {
        return isAutoValidated;
    }

    /**
     * @param isAutoValidated the isAutoValidated to set
     */
    public void setAutoValidated(boolean isAutoValidated) {
        this.isAutoValidated = isAutoValidated;
    }
    
}
