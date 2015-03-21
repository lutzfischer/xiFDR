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
package org.rappsilber.fdr.gui.components;

import java.awt.event.FocusAdapter;
import java.awt.event.FocusEvent;
import java.text.Format;
import java.text.NumberFormat;
import java.text.ParseException;
import javax.swing.JComponent;
import javax.swing.JSpinner;
import javax.swing.JTextField;
import javax.swing.SpinnerModel;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.text.DefaultFormatterFactory;
import javax.swing.text.NumberFormatter;

/**
 * Just an adaption of the JSpinner that replaces 0 with "unrestricted"
 * @author lfischer
 */
public class SingleTextValueNumericSpinner extends javax.swing.JSpinner {
    /** what to display for the special value */
    private String specialValueText = "Unrestricted";

    /** what value to replace with a text */
    private Number specialValue = 0.0;
    
    
    

    /** 
     * The model for the spinner 
     * The Spinner NumberModel does something fancy, that breaks if I don't overwrite some functions
     * 
     */
    public class AmbiguitySpinnerModel extends javax.swing.SpinnerNumberModel {
        Number value = 0;
        Comparable minimum = 0 ;
        Comparable maximum = null;
        Number stepSize = 1;
        
        public AmbiguitySpinnerModel(Number value, Comparable minimum, Comparable maximum, Number stepSize) {
            super(value, minimum, maximum, stepSize);
            this.minimum = minimum;
            this.maximum = maximum;
            this.value = value;
            this.stepSize = stepSize;
        }
        public <T>AmbiguitySpinnerModel() {
            this(Integer.valueOf(0), new Integer(0), null, Integer.valueOf(1));
        }
        
        public Object getValue() {
            return value;
        }
        
        public void setValue(Object value) {
            Number newValue = this.value;
            if (value.toString().contains(getSpecialValueText()))
                setValue(specialValue);
            try {
                newValue = (Number)value;
            } catch (Exception ne) {
//                value = minimum;
            }
            
            if ((minimum != null &&  ((Number)minimum).doubleValue() > newValue.doubleValue()))
                newValue = (Number)minimum;

            if ((maximum != null &&  ((Number)maximum).doubleValue() < newValue.doubleValue()))
                newValue = (Number)maximum;
            
//            if (!newValue.equals(this.value)) {
                
                if (this.value instanceof Double) {
                    newValue = Double.valueOf(newValue.doubleValue());
                } 
                else if (this.value instanceof Float) {
                    newValue = Float.valueOf(newValue.floatValue());
                } 
                else if (this.value instanceof Long) {
                    newValue = Long.valueOf(newValue.longValue());
                }
                else if (this.value instanceof Integer) {
                    newValue = Integer.valueOf(newValue.intValue());
                }
                else if (this.value instanceof Short) {
                    newValue = Short.valueOf(newValue.shortValue());
                }
                else {
                    newValue = Byte.valueOf(newValue.byteValue());
                }
                
                this.value = newValue;
                fireStateChanged();
                
  //          }
        }        

        public Object getNextValue() {
            Double ret = value.doubleValue()+stepSize.doubleValue();
            if (maximum == null || ((Number)maximum).doubleValue() > ret)
                return ret;
            else
                return maximum;
        }

        public Object getPreviousValue() {
            Double ret = value.doubleValue()-stepSize.doubleValue();
            if (minimum == null || ((Number)minimum).doubleValue() < ret)
                return ret;
            else
                return minimum;
        }
    }

    public class AmbiguitySpinnerEditor extends NumberEditor {

        
        
        public AmbiguitySpinnerEditor(JSpinner spinner) {
            super(spinner);
        // You have to do in your contructor
            AmbiguitySpinnerFormatter formatter = new AmbiguitySpinnerFormatter (NumberFormat.getInstance());
            getTextField().setFormatterFactory(new DefaultFormatterFactory(formatter));
        }


        
    }    
    
    public class AmbiguitySpinnerFormatter extends NumberFormatter {

        public AmbiguitySpinnerFormatter(NumberFormat format) {
            super(format);
        }
        
        public Object stringToValue(String text) {
            if (text.contains(getSpecialValueText()))
                return specialValue;
            Object ret;
            try {
                ret = super.stringToValue(text);
            } catch (ParseException e) {
                return specialValue;
            }
            return ret;
        }
        
    
        /**
         * Returns a String representation of the Object <code>value</code>.
         * This invokes <code>format</code> on the current <code>Format</code>.
         *
         * @throws ParseException if there is an error in the conversion
         * @param value Value to convert
         * @return String representation of value
         */
        public String valueToString(Object value) throws ParseException {
            
            if (value == null) {
                return "";
            }
            
            if (value.toString().contains(getSpecialValueText()))
                return getSpecialValueText();
            if (((Number)value).doubleValue() == getSpecialValue().doubleValue())
                return getSpecialValueText();
            return super.valueToString(value);
        }
        
    
    
        
        // Just override the methos StringToValue and ValueToString.
        // You can check here if the value is special
        // i.e you must display its special text instead. e.g. : "Auto" instead of -1
    }    

    public SingleTextValueNumericSpinner() {
        setModel(new AmbiguitySpinnerModel());
        AmbiguitySpinnerEditor editor = new AmbiguitySpinnerEditor(this);
        setEditor(editor);
        final JTextField tf = editor.getTextField();
        tf.addFocusListener(new FocusAdapter() {
            @Override
            public void focusGained(FocusEvent e) {
                if (((Number)getModel().getValue()).doubleValue() == specialValue.doubleValue()) {
                    JTextField tf = (JTextField) e.getSource();
                    if (tf.getText().contentEquals(specialValueText)) {
                        tf.setText(getModel().getValue().toString());
                        tf.selectAll();
                    }
                }
            }
        });
    }
    
    
    
    protected JComponent createEditor(SpinnerModel model) {
        return new AmbiguitySpinnerEditor(this);
    }    

    

    /**
     * what to display for the special value
     * @return the specialValueText
     */
    public String getSpecialValueText() {
        if (specialValueText == null)
            return "";
        return specialValueText;
    }

    /**
     * what to display for the special value
     * @param specialValueText the specialValueText to set
     */
    public void setSpecialValueText(String specialValueText) {
        this.specialValueText = specialValueText;
    }

    /**
     * what value to replace with a text
     * @return the specialValue
     */
    public Number getSpecialValue() {
        return specialValue;
    }

    /**
     * what value to replace with a text
     * @param specialValue the specialValue to set
     */
    public void setSpecialValue(Number specialValue) {
        this.specialValue = specialValue;
        this.fireStateChanged();
    }
        
    
    public void setRanges(Number value, Comparable minimum, Comparable maximum, Number stepSize) {
        AmbiguitySpinnerEditor ed = (AmbiguitySpinnerEditor) getEditor();
        AmbiguitySpinnerModel model = new AmbiguitySpinnerModel(value, minimum, maximum, stepSize);
        super.setModel(model );
        setEditor(ed);
        setValue(model.getNextValue());
        setValue(model.getPreviousValue());
        setValue(value);
    }
    
    @Override
    public void setModel(SpinnerModel model) {
        if (model instanceof javax.swing.SpinnerNumberModel) {
            javax.swing.SpinnerNumberModel snm = (javax.swing.SpinnerNumberModel) model;
            setRanges(snm.getNumber(), snm.getMinimum(), snm.getMaximum(), snm.getStepSize());
        }
    }
    
}
