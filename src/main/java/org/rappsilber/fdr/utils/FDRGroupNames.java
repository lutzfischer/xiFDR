/*
 * Copyright 2018 Lutz Fischer <lfischer@staffmail.ed.ac.uk>.
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
package org.rappsilber.fdr.utils;

import java.util.HashMap;

/**
 *
 * @author Lutz Fischer <lfischer@staffmail.ed.ac.uk>
 */
public abstract class FDRGroupNames {
    private static HashMap<String,String> ALL_FDR_GROUP_NAMES = new HashMap<>();
    
    public static String get(String g) {
        String ret = ALL_FDR_GROUP_NAMES.get(g);
        if (ret == null) {
            ret = g;
            ALL_FDR_GROUP_NAMES.put(ret, ret);
        }
        return ret;
    }

}
