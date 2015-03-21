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
package org.rappsilber.data.csv;

/**
 * interface that is used to calculate default-values for a column based on other values in that row
 * @author lfischer
 */
public interface CSVValueCalc {
    /**
     * calculate a value for the current row of the CSVParser
     * @param csv
     * @return 
     */
    public String getValue(CsvParser csv);
    /**
     * calculate a value for the given row of the CSVParser
     * @param csv
     * @param row
     * @return 
     */
    public String getValue(CSVRandomAccess csv, int row);
}
