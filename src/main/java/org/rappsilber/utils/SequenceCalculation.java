/*
 * Copyright 2017 lfischer.
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
package org.rappsilber.utils;

import java.util.ArrayList;

/**
 *
 * @author lfischer
 */
public abstract class SequenceCalculation {
    public static interface MassCalculator {
        double mass(String sequence);
    }
    
    static MassCalculator m_masscalc;
    
    public static void setCalculator(MassCalculator cal) {
        m_masscalc = cal;
    }
    
    public static MassCalculator getMassCalculator() {
        return m_masscalc;
    }
}
