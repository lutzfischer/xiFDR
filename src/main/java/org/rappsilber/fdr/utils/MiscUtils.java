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
package org.rappsilber.fdr.utils;

import java.util.ArrayList;

/**
 *
 * @author lfischer
 */
public abstract class MiscUtils {
    
    public static double minDiff(ArrayList<Double> list, boolean gz) {
        int l = list.size();
        double minDiff = Double.POSITIVE_INFINITY;
        if (gz) {
            for (int f=0; f<l; f++) {
                for (int s = f+1; s<l;s++) {
                    double diff = Math.abs(list.get(f) - list.get(s));
                    if ( diff >0)
                        minDiff=Math.min(minDiff, diff);
                }
            }
            if (minDiff==Double.POSITIVE_INFINITY)
                return 0;
            
            return minDiff;
        }
        // we don't care about zeros
        for (int f=0; f<l; f++) {
            for (int s = f+1; s<l;s++) {
                double diff = Math.abs(list.get(f) - list.get(s));
                minDiff=Math.min(minDiff, diff);
            }
        }

        return minDiff;
    }

    /**
     * minimum difference between any of the values in the array
     * if 
     * @param list
     * @param greaterZero
     * @return 
     */
    public static double minDiff(double[] list, boolean greaterZero) {
        int l = list.length;
        double minDiff = Double.POSITIVE_INFINITY;
        if (greaterZero) {
            for (int f=0; f<l; f++) {
                for (int s = f+1; s<l;s++) {
                    double diff = Math.abs(list[f] - list[s]);
                    if ( diff >0)
                        minDiff=Math.min(minDiff, diff);
                }
            }
            if (minDiff==Double.POSITIVE_INFINITY)
                return 0;
            
            return minDiff;
        }
        // we don't care about zeros
        for (int f=0; f<l; f++) {
            for (int s = f+1; s<l;s++) {
                double diff = Math.abs(list[f] - list[s]);
                minDiff=Math.min(minDiff, diff);
            }
        }

        return minDiff;
    }

    public static String formatStringForPrettyPrintingRelatedValues(double[] values, int minDigits) {
        
        double diff = minDiff(values, true);

        // turn that into the decimal by using log10 and round that down
        double roundingFactor = Math.round(Math.log(diff)/Math.log(10)-0.5);

        // we don't want to round away any digits before the dot
        if (roundingFactor>-minDigits)
            roundingFactor=-minDigits;
        if (roundingFactor<-10)
            roundingFactor=-10;
        if (roundingFactor>10)
            roundingFactor=10;

        double normFactor = Math.pow(10, roundingFactor);
        //if (roundingFactor<0)
            roundingFactor=-roundingFactor;

        return "%."+ String.format("%.0f",roundingFactor) + "f";
    }
    
    public static String formatStringForPrettyPrintingRelatedValues(double[] values, int minDigits, int maxDigits) {

        double[] roundedvalues = new double[values.length];
        double factor = Math.pow(10, maxDigits);

        // make sure that no difference smaller then the number of maximum 
        // digits is still there, by rounding all values to these digits.
        for (int i = 0 ; i< values.length; i++) {
            roundedvalues[i] = Math.round(values[i]*factor)/factor;
        }
        
        // get the smallest difference greater than 0
        double diff = minDiff(roundedvalues, true);
        

        // turn that into the decimal by using log10 and round that down
        double roundingFactor = Math.round(Math.log(diff)/Math.log(10)-0.5);

        // we don't want to round away any digits before the dot
        if (roundingFactor>-minDigits)
            roundingFactor=-minDigits;

        double normFactor = Math.pow(10, roundingFactor);
        //if (roundingFactor<0)
            roundingFactor=-roundingFactor;

        return "%."+ String.format("%.0f",roundingFactor) + "f";
    }
    


    
    
    public static String[] arrayToStringWithDifferenceOrientedFormat(double[] values, int minDigits) {
        
        String formatString = formatStringForPrettyPrintingRelatedValues(values, minDigits);

        String[] ret = new String[values.length];
        for (int i = 0; i< values.length;i++) {
            ret[i] = String.format(formatString, values[i]);
        }
        return ret;
        
    }
    
}
