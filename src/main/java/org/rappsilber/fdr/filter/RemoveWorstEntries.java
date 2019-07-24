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
package org.rappsilber.fdr.filter;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Map;
import org.rappsilber.fdr.entities.PSM;

/**
 *
 * @author Lutz Fischer <lfischer@staffmail.ed.ac.uk>
 */
public class RemoveWorstEntries implements PSMFilter{
    /**
     * try to remove the worst part that on itself has a worse FDR then
     */
    private double targetPercent = 0.95;
    
    private boolean keep = true;

    private double maxReducedFDR = 0.50;

    private int windowSize = 1000;
    
    /**
     * only remove if we don't remove more then this relative amount of PSM
     */
    private double maxReduction = 0.25;

    public RemoveWorstEntries() {
    }
    
    
    public RemoveWorstEntries(double targetPercent) {
        this.targetPercent = targetPercent;
    }
    
    @Override
    public ArrayList<PSM> filter(ArrayList<PSM> psms) {
        if (targetPercent == 0 ||  targetPercent == 1)
            return psms;
        if (psms.size() == 0)
            return psms;
        PSM first =  psms.get(0);
        ArrayList<PSM> all = new ArrayList<>(psms.size());
        HashMap<String,Object> infos = first.getOtherInfo();
        ArrayList<String> comparableInfos = new ArrayList<>(infos.size());
        HashMap<String,Object> lowcutoffs = new HashMap<>();
        HashMap<String,Object> highcutoffs = new HashMap<>();
        
        int allTT = 0;
        int allTD = 0;
        int allDD = 0;
        for (PSM p : psms) {
            if (p.isTT())
                allTT++;
            else if (p.isTD())
                allTD++;
            else
                allDD++;
            all.add(p);
        }
        
        double allFDR=(allTD - allDD)/(double)allTT;
        
        for (Map.Entry<String,Object> e : infos.entrySet()) {
            Object o = e.getValue();
            if (o instanceof Comparable) {
                // sort by this value
                final String name = e.getKey();
                java.util.Collections.sort(all, new Comparator<PSM>() {
                    @Override
                    public int compare(PSM arg0, PSM arg1) {
                        return ((Comparable) arg0.getOtherInfo().get(name)).compareTo(arg1.getOtherInfo().get(name));
                    }
                });
                // take the higest portion and the lowest portion of the data
                int lowLimit = (int)(psms.size()*maxReduction);


                defineLimits(all, name, false, keep?1-targetPercent:targetPercent, lowcutoffs, allFDR);
                
                defineLimits(all, name, true, keep?1-targetPercent:targetPercent, highcutoffs, allFDR);
                
            }
        }
        
        if (lowcutoffs.size() == 0 && highcutoffs.size() == 0)
            return psms;
        ArrayList<PSM> ret = new ArrayList<>();

        // delete everything in that is below the lowcutoffs
        for (PSM psm : psms) {
            boolean keep = true;
            for (Map.Entry<String,Object> e : lowcutoffs.entrySet()) {
                if (((Comparable) psm.getOtherInfo().get(e.getKey())).compareTo(e.getValue())<=0) {
                    keep = false;
                    break;
                }
            }

            if (keep) {
                for (Map.Entry<String,Object> e : highcutoffs.entrySet()) {
                    if (((Comparable) psm.getOtherInfo().get(e.getKey())).compareTo(e.getValue())>=0) {
                        keep = false;
                        break;
                    }
                }
            }
            if (keep) {
                ret.add(psm);
            }
        }
        return ret;
        
    }

    protected void defineLimits(ArrayList<PSM> all, final String name, boolean fromHigh, double percentDelete, HashMap<String, Object> cutoffs, double allfdr) {
        // count TT TD and DD in both sets
        double TT = 0;
        double TD = 0;
        double DD = 0;
        Boolean cutSensible  = null;
        double cutSensibleFDR = 0;
        
        LinkedList<PSM> window = new LinkedList<>();
        int increment = -1;
        int cutpoint =(int) (all.size()*percentDelete);;
        if (fromHigh) {
            if (cutpoint == 0)
                return;
            cutpoint = all.size()-cutpoint;
            increment=1;
        }
        
        Comparable o = (Comparable) all.get(cutpoint).getOtherInfo().get(name);
        Comparable c = null;
        
        while (o.equals(c = (Comparable) all.get(cutpoint).getOtherInfo().get(name))) {
            cutpoint += increment;
            if (cutpoint<0 || cutpoint>=all.size())
                return;
        }
        
        for(;cutpoint<0 || cutpoint>=all.size();cutpoint += increment) {
            PSM p = all.get(cutpoint);
            if (p.isTT())
                TT++;
            else if (p.isTD())
                TT++;
            else if (p.isDD())
                TT++;
        }
        // only accept the cutoff if we actually remove proportionally more decoys then targets
        if ((TD-DD)/TT > allfdr)
            cutoffs.put(name, c);
    }

    /**
     * try to remove the worst part that on itself has a worse FDR then
     * @return the targetPercent
     */
    public double getTargetPercent() {
        return targetPercent;
    }

    /**
     * try to remove the worst part that on itself has a worse FDR then
     * @param targetPercent the targetPercent to set
     */
    public void setTargetPercent(double targetPercent) {
        this.targetPercent = targetPercent;
    }

    /**
     * only remove if we don't remove more then this relative amount of PSM
     * @return the maxReduction
     */
    public double getMaxReduction() {
        return maxReduction;
    }

    /**
     * only remove if we don't remove more then this relative amount of PSM
     * @param maxReduction the maxReduction to set
     */
    public void setMaxReduction(double maxReduction) {
        this.maxReduction = maxReduction;
    }
}
