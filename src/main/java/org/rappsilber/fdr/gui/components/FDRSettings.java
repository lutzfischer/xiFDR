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
package org.rappsilber.fdr.gui.components;

import java.awt.event.ActionListener;
import org.rappsilber.fdr.OfflineFDR;

/**
 *
 * @author lfischer
 */
public interface FDRSettings {

    void addCalcListener(ActionListener al);

    OfflineFDR.FDRLevel doOptimize();

    int getMaxLinkAmbiguity();


    int getMaxProteinAmbiguity();

    int getMinLinkPepCount();

    int getMinPPIPepCount();

    int getMinPeptideLength();

    int getMinProteinPepCount();

    double getPSMFDR();

    double getPeptidePairFDR();

    double getProteinGroupFDR();

    double getProteinGroupLinkFDR();

    double getProteinGroupPairFDR();

    void setPSMFDR(Double fdr);
    
    void setPeptidePairFDR(Double fdr);

    void setProteinGroupFDR(Double fdr);

    void setProteinGroupLinkFDR(Double fdr);
    
    
    int getBoostingSteps();

    boolean getBoostBetween();
    void setBoostBetween(boolean between);
    
    boolean isLinkDirectional();

    boolean isPPIDirectional();

    boolean isPSMDirectional();

    boolean isPeptidePairDirectional();

    void removeCalcListener(ActionListener al);

    void setLinkDirectional(boolean directional);

    void setMaxLinkAmbiguity(Integer maxAmbiguity);


    void setMaxProteinAmbiguity(Integer maxAmbiguity);

    void setMinLinkPepCount(Integer minPep);

    void setMinPPIPepCount(Integer minPep);

    void setMinPeptideLength(Integer minLength);

    void setMinProteinPepCount(Integer minPep);

    void setPPIDirectional(boolean directional);

    void setPSMDirectional(boolean directional);

    void setPeptidePairDirectional(boolean directional);


    void setProteinGroupPairFDR(Double fdr);
    
    double getReportFactor();
    
    public void setReportFactor(double factor);       
    
    boolean filterToUniquePSM();
    
    void setFilterToUniquePSM(boolean filterToUnique);
    
    
    public void setAll(FDRSettings settings);
    
    
}
