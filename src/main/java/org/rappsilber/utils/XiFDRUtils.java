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
package org.rappsilber.utils;

import java.io.File;

/**
 *
 * @author Lutz Fischer <lfischer@staffmail.ed.ac.uk>
 */
public class XiFDRUtils {
    public static class FDRCSVOUT {
        public String basename;
        public String folder;
        public String extension;
        public boolean tsv;
        public boolean csv;
    }
    
    public static FDRCSVOUT splitFilename(String filename) {
        FDRCSVOUT out = new FDRCSVOUT();
        File file = new File(filename);
        String name = file.getName();
        out.extension = "";
        if (name.contains(".")) {
            if (name.toLowerCase().endsWith(".tsv")
                    || name.toLowerCase().endsWith(".txt")) {
                out.tsv = true;
            } else if (name.toLowerCase().endsWith(".csv")) {
                out.csv = true;
            }
            out.extension = name.substring(name.lastIndexOf("."));
        }
        out.folder = file.getParent();
        out.basename = name.replaceAll("\\.[^\\.]*$", "");
        if (out.basename.matches("^.*_[^_]*_xiFDR.*$")) {
            out.basename = out.basename.substring(0, out.basename.indexOf("_xiFDR"));
            if (out.basename.contains("_Linear_")) {
                out.basename = out.basename.substring(0, out.basename.lastIndexOf("_Linear_"));
            } else {
                out.basename = out.basename.substring(0, out.basename.lastIndexOf("_"));
            }
        }
        return out;
    }
}
