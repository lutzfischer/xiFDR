/*
 * Copyright 2020 lfischer.
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
package org.rappsilber.fdr.gui.components;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import javax.swing.AbstractSpinnerModel;

//    private ArrayList<java.awt.event.ActionListener> m_calc_listener = new ArrayList<ActionListener>();

//
//    private OfflineFDR.FDRLevel m_optimizeWhat;
//
public class FDRSpinnerModel extends AbstractSpinnerModel {

    private double value = 100.0;
    private double increment = 1;
    private double max = 100;
    private double min = 0;
    String sep;

    public FDRSpinnerModel() {
        //this.value=value;
        DecimalFormat format = (DecimalFormat) DecimalFormat.getInstance();
        DecimalFormatSymbols symbols = format.getDecimalFormatSymbols();
        sep = symbols.getDecimalSeparator() + "";
    } //this.value=value;

    public FDRSpinnerModel(Double value) {
        this();
        this.value = value;
    }

    @Override
    public Object getNextValue() {
        setValue(Math.round((value + increment) / increment) * increment);
        return value + ""; //return as string to avoid round
    }

    @Override
    public Object getPreviousValue() {
        setValue(Math.round((value - increment) / increment) * increment);
        return value + ""; //return as string to avoid round
    }

    @Override
    public Object getValue() {
        if (value == (int) value) {
            return ((int) value) + "";
        }
        return value + ""; //return as string to avoid round
    }

    public double getFDR() {
        return value/100; //return as string to avoid round
    }

    public void setFDR(double fdr) {
        value = fdr*100; //return as string to avoid round
        fireStateChanged();
    }

    @Override
    public void setValue(Object o) {
        double newValue = Double.NaN;
        if (o instanceof Number) {
            newValue = ((Number) o).doubleValue();
        } else {
            try {
                newValue = Double.parseDouble(o.toString());
            } catch (Exception e) {
                newValue = value;
            }
        }
        if (newValue < min) {
            newValue = min;
        }
        if (newValue > max) {
            newValue = max;
        }
        String svalue = value + "";
        // remove the odd rounding errors of somethng .000000000001
        if (svalue.contains("00000")) {
            newValue = Double.parseDouble(svalue.substring(0, svalue.indexOf("00000")));
        }
        if (newValue - value != 0) {
            value = newValue;
            fireStateChanged();
        } else if (o.toString().contains(sep) && !(newValue + "").contentEquals(o.toString())) {
            fireStateChanged();
        }
    }

    /**
     * @return the increment
     */
    public double getIncrement() {
        return increment;
    }

    /**
     * @param increment the increment to set
     */
    public void setIncrement(double increment) {
        this.increment = increment;
    }
    
}
