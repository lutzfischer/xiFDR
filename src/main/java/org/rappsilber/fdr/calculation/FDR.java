/*
 * Copyright 2020 Lutz Fischer <lfischer@staffmail.ed.ac.uk>.
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
package org.rappsilber.fdr.calculation;

import java.util.Collection;
import org.rappsilber.fdr.FDRSettings;
import org.rappsilber.fdr.entities.AbstractFDRElement;
import org.rappsilber.fdr.result.FDRResultLevel;

/**
 *
 * @author Lutz Fischer <lfischer@staffmail.ed.ac.uk>
 */
public interface FDR {

    //    private <T extends FDRSelfAdd<T>> HashedArrayList<T> fdr(double fdr, double safetyfactor, Collection<T> c, HashMap<Integer, SubGroupFdrInfo<T>> groupInfo, double tCount, double dCount, int minPepCount, boolean ignoreGroups,boolean setElementFDR) {
    /**
     * Splits up data into groups and forward each group to
     * {@link #subFDR(java.util.ArrayList, java.util.ArrayList, boolean, rappsilber.fdr.result.SubGroupFdrInfo) subFDR}.
     * results will then be collected
     *
     * @param <T> What type of Information is the FDRImplement to be applied
     * @param fdr the target FDRImplement that is supossed to be returned
     * @param safetyfactor don't report anything if the next higher calculatable
     * @param fdrInput the input data that should be filtered by FDRImplement
     * @param groupInfo will hold all the information for the currently
     * calculated FDRs
     * @param tCount number of entries in the target database. This is used to
    normalise the FDRImplement calculation for imbalenced databases
     * @param dCount number of entries in the decoy database. This is used to
    normalise the FDRImplement calculation for imbalenced databases
     * @param minPepCount if an entry has less then these number of peptides
    supporting it - it will not be considered for the FDRImplement-calculation
     * @param ignoreGroups should we ignore groups cmpletely and just calculate
    a joined FDRImplement?
     * @param setElementFDR should each element be flaged up with the FDRImplement that
    it group has at the given score? - is used to speed up the maximation (by
    not setting these)
     */
    <T extends AbstractFDRElement<T>> void fdr(double fdr, FDRSettings settings, Collection<T> fdrInput, FDRResultLevel<T> groupInfo, double tCount, double dCount, int minPepCount, boolean ignoreGroups, boolean setElementFDR, Boolean localFDR, boolean groupByProteinPair);
    
    void setValidityCheck(CheckValid c);
}
