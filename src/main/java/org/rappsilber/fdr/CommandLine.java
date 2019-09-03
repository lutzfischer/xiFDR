/*
 * Copyright 2019 Lutz Fischer <lfischer@staffmail.ed.ac.uk>.
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
package org.rappsilber.fdr;

import java.io.FileNotFoundException;
import java.sql.SQLException;
import org.rappsilber.fdr.gui.FDRGUI;

/**
 * Proxy class that looks at the given arguments and then forwards the whole 
 * command-line to the appropriate class.
 * @author Lutz Fischer <lfischer@staffmail.ed.ac.uk>
 */
public class CommandLine {
    
    public static void main(String args[]) throws SQLException, FileNotFoundException {
        boolean xicsv  = false;
        boolean csv = false;
        boolean gui = true;
        boolean db = false;
        for (String a : args) {
            gui = false;
            if (a.startsWith("--xiconfig=")) {
                xicsv=true;
            } else if (a.startsWith("--dbids=")) {
                db = true;
            } else if (!a.startsWith("-")) {
                csv = true;
            }
        }
        if (db && !csv) {
            DBinFDR.main(args);
            return;
        }
        if (xicsv && !db) {
            XiCSVinFDR.main(args);
            return;
        }
        if (csv && (!db) && (!xicsv)) {
            CSVinFDR.main(args);
            return;
        }
        if (gui) {
            FDRGUI.main(args);
            return;
        }
        new XiCSVinFDR().printUsage();
    }
}
