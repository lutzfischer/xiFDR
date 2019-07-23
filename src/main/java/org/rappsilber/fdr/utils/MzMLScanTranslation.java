/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.rappsilber.fdr.utils;

import java.util.regex.Matcher;
import java.util.regex.Pattern;
import org.rappsilber.fdr.entities.DBPSM;

/**
 *
 * @author Lutz Fischer <lfischer@staffmail.ed.ac.uk>
 */
public class MzMLScanTranslation implements PSMToMzIdentScanId{
    private String template;
    private String toReplace="%s%";
    private int offset = 0;
    private static Pattern hasOffset = Pattern.compile("(.*%s)([+\\-][0-9]+)(%.*)");

    public MzMLScanTranslation(String template) {
        this.template = template;
        Matcher m = hasOffset.matcher(template);
        if (m.matches()) {
            offset = Integer.parseInt(m.group(2));
            this.toReplace = "%s"+m.group(2)+"%";
        }
    }
    
    String getTemplate() {
        return template;
    }
    
    @Override
    public String getID(DBPSM psm) {
        return template.replace(toReplace, ""+(Integer.parseInt(psm.getScan())+offset));
    }
    
}
