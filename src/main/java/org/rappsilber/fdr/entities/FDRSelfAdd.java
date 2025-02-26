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
package org.rappsilber.fdr.entities;

import org.rappsilber.utils.SelfAdd;

/**
 *
 * @author lfischer
 */
public interface FDRSelfAdd<T> extends SelfAdd<T> { //, Comparable<T> {
    /**
     * The score of the match, but only consider the topN scores.
     * @return 
     */
    double getScore(int topN);
    /**
     * The score of the match
     * @return 
     */
    double getScore();
    /**
     * an ID for the FDRgroup the match belongs to.
     * @return 
     */
    String getFDRGroup();
    
    /**
     * defines the fdr-group for a single entry
     */
    void setFDRGroup(String fdrGroup);
//    /**
//     * Name of the FDR-group this match belongs to
//     * @return 
//     */
//    String getFDRGroupName();
    /**
     * Non of the sites are decoy
     * @return 
     */
    boolean isTT();
    /**
     * one of the sites is decoy
     * @return 
     */
    boolean isTD();
    /**
     * all sites are decoy
     * @return 
     */
    boolean isDD();
    /**
     * is any site decoy
     * @return 
     */
    boolean isDecoy();
    /**
     * is this a linear match
     * @return 
     */
    boolean isLinear();
    /**
     * could this be a protein internal match
     * <p>This does refer to protein type not necessarily a within one molecule
     * </p> 
     * @return 
     */
    boolean isInternal();
    /**
     * is this definitely not a protein internal match
     * @return 
     */
    boolean isBetween();
    /**
     * set the FDR, that is assigned to the score of the match
     * @param fdr 
     */
    void setFDR(double fdr);
    /**
     * get the FDR that is assigned to the score of this match
     * @return 
     */
    double getFDR();

    /**
     * set the local FDR, of a match
     * @param fdr 
     */
    void setPEP(double pep);
    /**
     * get the local FDR that is assigned to the score of this match
     * @return 
     */
    Double getPEP();
    
    
    /**
     * how many unique peptides support this match
     * @return 
     */
    int getPeptidePairCount();
    /**
     * What is the next higher FDR then the score assigned FDR
     * @return 
     */
    double getHigherFDR();
    /**
     * What is the next higher FDR then the score assigned FDR
     * @param fdr
     */
    void setHigherFDR(double fdr);
    /**
     * What is the next lower FDR then the score assigned FDR     
     * @return 
     */
    double getLowerFDR();
    /**
     * What is the next lower FDR then the score assigned FDR
     * @param fdr 
     */
    void setLowerFDR(double fdr);
    /**
     * get the peptide/protein group at the given site of the match
     * @param n
     * @return 
     */
    Object getSite(int n);
    /**
     * how many sites does this match have.
     * Currently 1 for linear and 2 for cross-linked
     * @return 
     */
    int getSites();
    
    /**
     * returns if this is derived from a non-covalent Peptide-Pair
     */
    public boolean isNonCovalent();
    
    
}
