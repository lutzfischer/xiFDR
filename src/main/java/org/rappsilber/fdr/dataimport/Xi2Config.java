/*
 * Copyright 2021 Lutz Fischer <lfischer@staffmail.ed.ac.uk>.
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
package org.rappsilber.fdr.dataimport;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.json.simple.parser.ParseException;
import org.json.simple.parser.JSONParser;


/**
 * Parses the minimal needed information from a xi2config json
 * @author Lutz Fischer <lfischer@staffmail.ed.ac.uk>
 */
public class Xi2Config {
    HashMap<String, Double> default_xl_masses = new HashMap<>();
    {
        default_xl_masses.put("BS3", 138.06807961);
        default_xl_masses.put("EDC", -18.01056027);
        default_xl_masses.put("SDA", 82.04186484);
        default_xl_masses.put("DSSO", 158.0038);
    }
    
    public class Crosslinker {
        public String name;
        public Double mass ;

        public Crosslinker(String name) {
            this.name = name;
            if (default_xl_masses.containsKey(name.toUpperCase())) {
                this.mass = default_xl_masses.get(name.toUpperCase());
            } else {
                this.mass = Double.NaN;
            }
        }
        
        public Crosslinker(String name, double mass) {
            this.name = name;
            this.mass = mass;
        }
        
    }
    public class Modification {
        String symbol;
        String type;

        public Modification(String symbol, String modtype) {
            this.symbol = symbol;
            this.type = modtype;
        }
        
    }
    public ArrayList<Crosslinker> crosslinker;
    public ArrayList<Modification> modifications;
    
    public String textConfig;

    public Xi2Config(String config) {
        this.crosslinker = new ArrayList<>(1);
        this.modifications = new ArrayList<>();
        this.textConfig = config;
                
        JSONParser prsr = new JSONParser();
        try {
            Map json = (Map)prsr.parse(config);
            List crosslinker = parseListSetting(json.get("crosslinker"));
            for (Object xlo : (List)crosslinker) {
                String name;
                if (xlo instanceof Map) {
                    name = ((Map)xlo).get("name").toString();
                    double mass = Double.parseDouble(((Map)xlo).get("mass").toString());
                    this.crosslinker.add(new Crosslinker(name, mass));
                } else {
                    name = xlo.toString();
                    this.crosslinker.add(new Crosslinker(name));
                }
            }
            List modifications = parseListSetting(((Map)json.get("modification")).get("modifications"));
            for (Object mo : modifications) {
                Map m = (Map)mo;
                this.modifications.add(new Modification(m.get("name").toString(), m.get("type").toString()));
            }
        } catch (ParseException ex) {
            Logger.getLogger(Xi2Config.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public Xi2Config(Xi2Config config) {
        this.crosslinker = new ArrayList<>(config.crosslinker);
        this.modifications = new ArrayList<>(config.modifications);
    }

    
    public Xi2Config add(Xi2Config o) {
        Xi2Config ret = new Xi2Config(this);
        ArrayList<Modification> addNewMod = new ArrayList<>();
        oloop: for (Modification mo : o.modifications) {
            for (Modification mt : ret.modifications) {
                if (mo.symbol.contentEquals(mt.symbol) && mo.type.contentEquals(mt.type))
                    continue oloop;
            }
            addNewMod.add(mo);
        }
        ret.modifications.addAll(addNewMod);
        
        ArrayList<Crosslinker> addNewXL = new ArrayList<>();
        oloop: for (Crosslinker co : o.crosslinker) {
            for (Crosslinker ct : ret.crosslinker) {
                if (co.name.contentEquals(ct.name) && Math.round(co.mass * 1000) == Math.round(ct.mass * 1000))
                    continue oloop;
            }
            addNewXL.add(co);
        }
        ret.crosslinker.addAll(addNewXL);
        return ret;
    }
    
    public List parseListSetting(Object json) {
        if (json instanceof List)
            return (List)json;
        List ret = new ArrayList();
        ret.add(json);
        return ret;
    }
    
    public String peptPeptide(String basepeptide, int[] mods) {
        StringBuffer ret = new StringBuffer();
        if (mods[0] != 0)
            ret.append(this.modifications.get(mods[0]).symbol);
        for (int i = 0; i<basepeptide.length(); i++) {
            ret.append(basepeptide.substring(i,i+1));
            if (mods[2+i] != 0)
                ret.append(this.modifications.get(mods[0]).symbol);
        }
        if (mods[1] != 0)
            ret.append(this.modifications.get(mods[1]).symbol);
        return ret.toString();
        
    }
    

    public static void main(String[] args) {
        String conf = "{\n" +
"  \"reporting_requirements\": {\n" +
"    \"report_top_ranking_only\": true\n" +
"  },\n" +
"  \"digestion\": {\n" +
"    \"enzymes\": [\"trypsin\"],\n" +
"    \"missed_cleavages\": 2,\n" +
"    \"min_peptide_length\": 2\n" +
"  },\n" +
"  \"isotope_error_ximpa\": 2,\n" +
"  \"ms1_tol\": \"3ppm\",\n" +
"  \"ms2_tol\": \"5ppm\",\n" +
"  \"top_n_alpha_scores\": 10,\n" +
"  \"top_n_alpha_beta_scores\": 10,\n" +
"  \"crosslinker\": \"BS3\",\n" +
"  \"conservative_n_multi_loss\": 3,\n" +
"  \"denoise_alpha\": {\n" +
"    \"top_n\": 10,\n" +
"    \"bin_size\": 100\n" +
"  },\n" +
"  \"denoise_alpha_beta\": {\n" +
"    \"top_n\": 20,\n" +
"    \"bin_size\": 100\n" +
"  },\n" +
"  \"fragmentation\": {\n" +
"    \"nterm_ions\": [\"b\"],\n" +
"    \"cterm_ions\": [\"y\"],\n" +
"    \"add_precursor\": true,\n" +
"    \"max_nloss\": 4,\n" +
"    \"match_missing_monoisotopic\": true,\n" +
"    \"losses\": [\n" +
"      {\n" +
"        \"name\": \"H2O\",\n" +
"        \"specificity\": [\"S\", \"T\", \"D\", \"E\", \"cterm\"],\n" +
"        \"composition\": \"\",\n" +
"        \"mass\": 18.01056027\n" +
"      },\n" +
"      {\n" +
"        \"name\": \"NH3\",\n" +
"        \"specificity\": [\"R\", \"K\", \"N\", \"Q\", \"nterm\"],\n" +
"        \"composition\": \"\",\n" +
"        \"mass\": 17.02654493\n" +
"      },\n" +
"      {\n" +
"        \"name\": \"CH3SOH\",\n" +
"        \"specificity\": [\"oxM\"],\n" +
"        \"composition\": \"\",\n" +
"        \"mass\": 63.99828547\n" +
"      }\n" +
"    ]\n" +
"  },\n" +
"  \"modification\": {\n" +
"    \"max_var_protein_mods\": 3,\n" +
"    \"max_modified_peps\": 20,\n" +
"    \"modifications\": [\n" +
"      {\n" +
"        \"name\": \"cm\",\n" +
"        \"specificity\": [\"C\"],\n" +
"        \"type\": \"fixed\",\n" +
"        \"composition\": \"C2H3N1O1\"\n" +
"      },\n" +
"      {\n" +
"        \"name\": \"ox\",\n" +
"        \"specificity\": [\"M\"],\n" +
"        \"type\": \"variable\",\n" +
"        \"composition\": \"O1\"\n" +
"      }\n" +
"    ]\n" +
"  }\n" +
"}";
        Xi2Config xiconf = new Xi2Config(conf);
        
    }
}
