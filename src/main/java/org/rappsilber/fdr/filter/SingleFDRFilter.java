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
import java.util.LinkedList;
import java.util.Map;
import org.rappsilber.fdr.entities.PSM;

/**
 *
 * @author Lutz Fischer <lfischer@staffmail.ed.ac.uk>
 */
public class SingleFDRFilter implements PSMFilter{
    /**
     * try to remove the worst part that on itself has a worse FDR then
     */
    private double targetPercent = 0.95;

    private double maxReducedFDR = 0.50;

    private int windowSize = 1000;
    
    private String name;
    
    /**
     * only remove if we don't remove more then this relative amount of PSM
     */
    private double maxReduction = 0.25;

     public SingleFDRFilter(String name, double targetPercent) {
        this.targetPercent = targetPercent;
    }
    
    @Override
    public ArrayList<PSM> filter(ArrayList<PSM> psms) {
        PSM first =  psms.get(0);
        ArrayList<PSM> all = new ArrayList<>(psms.size());
        HashMap<String,Object> infos = first.getOtherInfo();
        ArrayList<String> comparableInfos = new ArrayList<>(infos.size());
        HashMap<String,Object[]> lowcutoffs = new HashMap<>();
        HashMap<String,Object[]> highcutoffs = new HashMap<>();
        
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


                defineLimits(all, name, 0, lowLimit, 1, lowcutoffs);

                defineLimits(all, name, all.size()-1, all.size()-lowLimit, -1, highcutoffs);
                if (lowcutoffs.containsKey(name) && highcutoffs.containsKey(name)) {
                    if (lowcutoffs.get(name).length > highcutoffs.get(name).length) {
                        highcutoffs.remove(name);
                    } else {
                        lowcutoffs.remove(name);
                    }
                }
                
            }
        }
        
        if (lowcutoffs.size() == 0 && highcutoffs.size() == 0)
            return psms;
        ArrayList<PSM> ret = new ArrayList<>();
        // delete everything in that is below the lowcutoffs
        for (PSM psm : psms) {
            boolean keep = true;
            for (Map.Entry<String,Object[]> e : lowcutoffs.entrySet()) {
                if (((Comparable) psm.getOtherInfo().get(e.getKey())).compareTo(e.getValue()[0])<=0) {
                    keep = false;
                    break;
                }
            }
            if (keep) {
                for (Map.Entry<String,Object[]> e : highcutoffs.entrySet()) {
                    if (((Comparable) psm.getOtherInfo().get(e.getKey())).compareTo(e.getValue()[0])>=0) {
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

    protected void defineLimits(ArrayList<PSM> all, final String name, int start, int lowLimit, int increment, HashMap<String, Object[]> lowcutoffs) {
        // count TT TD and DD in both sets
        double windowTT = 0;
        double windowTD = 0;
        double windowDD = 0;
        Boolean cutSensible  = null;
        double cutSensibleFDR = 0;
        
        LinkedList<PSM> window = new LinkedList<>();
        
        Object last = all.get(0).getOtherInfo().get(name);
        int c = 0;
        for (int i = start; i!=lowLimit; i+=increment) {
            PSM l = all.get(i);
            c++;
            // current score
            Object s = l.getOtherInfo().get(name);
            if (l.isTT())
                windowTT++;
            if (l.isTD())
                windowTD++;
            if (l.isDD())
                windowDD++;
            window.add(l);
            if (window.size()>windowSize) {
                if (cutSensible == null) {
                    cutSensibleFDR = (windowTD-windowDD)/windowTT;
                    cutSensible = cutSensibleFDR > targetPercent;
                    if (!cutSensible) {
                        break;
                    }
                }
                double windowFDR = (windowTD-windowDD)/windowTT;
                if (cutSensible && windowFDR<= targetPercent && windowFDR >maxReducedFDR  && !last.equals(s)) {
                    lowcutoffs.put(name,new Object[]{last,c});
                    break;
                }
                PSM r = window.removeFirst();
                if (r.isTT())
                    windowTT--;
                if (r.isTD())
                    windowTD--;
                if (r.isDD())
                    windowDD--;
            }
            last = l.getOtherInfo().get(name);
        }
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
