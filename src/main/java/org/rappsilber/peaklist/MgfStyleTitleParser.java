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
package org.rappsilber.peaklist;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.regex.Pattern;

/**
 *
 * @author Lutz Fischer <lfischer@staffmail.ed.ac.uk>
 */
public class MgfStyleTitleParser {
    public class ParseEnrtry{
        public int index;
        public double expMZ;
        public byte expCharge;
        public String peakFileName;
    }
    String header="TITLE=";
    String pepMZ="PEPMASS=";
    String pepCharge="CHARGE=";
    String scanEnd="END IONS";
    
    private File parsedFile;
    public String fileName;
    HashMap<String,ArrayList<String>> parts2titles = new HashMap<String, ArrayList<String>>();
    HashMap<Integer,ArrayList<String>> int2titles = new HashMap<Integer, ArrayList<String>>();
    Pattern isNummeric = Pattern.compile("^[\\+\\-]?[0-9]+(:?\\.[0-9]+)?$");
    
    HashMap<String,ParseEnrtry> titlesToID = new HashMap<String, ParseEnrtry>();

    public int parseFile(File f) throws IOException {
        if (f.getName().endsWith(".apl"))
            return parseFile(f,"header=", "mz=", "charge=","peaklist end", null);
        return parseFile(f,header, pepMZ, pepCharge, scanEnd, null);
    }

    private void registerTitleParts(String key, String title, ParseEnrtry entry) {
        titlesToID.put(title, entry);
        ArrayList<String> titles = parts2titles.get(key);
        if (titles == null) {
            titles = new ArrayList<>();
            parts2titles.put(key, titles);
        }
        if (isNummeric.matcher(key).matches()) {
            Integer ikey =Integer.parseInt(key);
            titles = int2titles.get(ikey);
            if (titles == null) {
                titles = new ArrayList<>();
                int2titles.put(ikey, titles);
            }
        }
            
        titles.add(title);
    }
    
    /**
     * reads in the file
     * @param f file to parse
     * @param header how does the header line starts
     * @param ignoreTitle a pattern that indicates what scans to ignore (i.g. apl-files wit
     * @return
     * @throws IOException 
     */
    public int parseFile(File f, String header,  String pepMZ, String pepCharge, String scanEnd, Pattern ignoreTitle) throws IOException {
        parsedFile = f;
        fileName = f.getName();
        this.header=header;
        this.pepMZ = pepMZ;
        this.pepCharge = pepCharge;
        this.scanEnd = scanEnd;
        BufferedReader in = new BufferedReader(new FileReader(f));
        String line;
        int index=-1;
        byte charge= 0;
        double mz = Double.NaN;
        String title=null;
        ArrayList<String> titelParts = new ArrayList<>();
        while ((line =in.readLine())!= null) {
            if (line.startsWith(scanEnd)) {
                index++;
                ParseEnrtry e = new ParseEnrtry();
                e.index = index;
                e.expMZ = mz;
                e.expCharge = charge;
                for (String p : titelParts)
                    registerTitleParts(p, title, e);
                mz = Double.NaN;
                charge = -1;
                titelParts.clear();
                title=null;
            } else if (line.startsWith(pepMZ)) {
                String smz = line.substring(pepMZ.length()).trim();
                smz = smz.split("\\s")[0];
                if (isNummeric.matcher(smz).matches()) {
                    mz = Double.parseDouble(smz);
                }
            } else if (line.startsWith(pepCharge)) {
                String scharge = line.substring(pepCharge.length()).trim();
                scharge = scharge.replaceAll("[\\+\\s]", "");
                if (isNummeric.matcher(scharge).matches()) {
                    charge = Byte.parseByte(scharge);
                }
                
            } else if (line.startsWith(header)) {
                title = line;
                if (ignoreTitle!= null && ignoreTitle.matcher(line).matches())
                    continue;
                String[] parts=line.split("\\b");
                for (String p : parts) {
                    if (p.trim().length()>0) {
                        titelParts.add(p);
                        String pt  = p.trim();
                        if (pt.length()!=p.length()) 
                            titelParts.add(pt);
                    }
                }
                parts=line.split("\\.");
                for (String p : parts) {
                    if (p.trim().length()>0) {
                        titelParts.add(p);
                        String pt  = p.trim();
                        if (pt.length()!=p.length()) 
                            titelParts.add(pt);
                    }
                }
                parts=line.split("\\s+");
                for (String p : parts) {
                    if (p.trim().length()>0)
                        titelParts.add(p);
                }
            }
        }
        return index+1;
    }
    
    /**
     * tries to guess the right index of a scan based on run and scan
     * @param run
     * @param scan
     * @return null if nothing could be found; -1 if no unambiguous scan was found; index of matching scan otherwise
     */
    public ParseEnrtry findScanIndex(String run, String scan) {
        HashSet<String> runTitles = new HashSet<String>();
        HashSet<String> scanTitles = new HashSet<String>();
        ArrayList<String> rt = parts2titles.get(run);
        if (rt != null) {
            runTitles.addAll(rt);
        } else if (run.contains(".")) {
            run = run.substring(0,run.indexOf("."));
            rt = parts2titles.get(run);
            if (rt!=null)
                runTitles.addAll(rt);
            else
                return null;
        }

        if (runTitles.isEmpty())
            return null;
        
        ArrayList<String> st = null;
        if (!isNummeric.matcher(scan).matches()) {
            st = parts2titles.get(scan);
        } else {
            st = int2titles.get(Integer.parseInt(scan));
        }
        if (st != null) {
            HashSet<String> runTitlesScan=new HashSet<>(runTitles);
            scanTitles.addAll(st);

            runTitlesScan.retainAll(scanTitles);

            if (runTitlesScan.size()==1) {
                return titlesToID.get(runTitlesScan.iterator().next());
            }

        }

        Pattern sp = Pattern.compile("^.*[^0-9]"+scan+"([^0-9].*)?$");
        for (String t : runTitles) {
            String tr = t.replace(run, "");
            if (sp.matcher(tr).matches())
                return titlesToID.get(t);
        }
        
        return null;
    }

    /**
     * @return the parsedFile
     */
    public File getParsedFile() {
        return parsedFile;
    }
}
