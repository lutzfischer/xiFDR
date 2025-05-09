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

import org.rappsilber.fdr.entities.AbstractFDRElement;
import org.rappsilber.fdr.result.SubGroupFdrInfo;

/**
 *
 * @author Lutz Fischer <lfischer@staffmail.ed.ac.uk>
 */
public interface CheckValid {

    /**
     * test whether a result for a subgroup should be considered valid
     *
     * @param <T>
     * @param info
     * @return null pass; otherwise reason
     */
    <T extends AbstractFDRElement<T>> String checkValid(SubGroupFdrInfo<T> info);

    <T extends AbstractFDRElement<T>> String checkValid(SubGroupFdrInfo<T> info, int minTD, double factor);
    
}
