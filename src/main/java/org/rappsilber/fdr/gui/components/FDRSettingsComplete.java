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

import javax.swing.AbstractSpinnerModel;
import javax.swing.JCheckBox;
import javax.swing.JFormattedTextField;
import javax.swing.JSpinner;
import javax.swing.SpinnerModel;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import org.rappsilber.fdr.OfflineFDR;

/**
 *
 * @author lfischer
 */
public class FDRSettingsComplete extends FDRSettingsPanel {

//    private ArrayList<java.awt.event.ActionListener> m_calc_listener = new ArrayList<ActionListener>();
//    
//    private OfflineFDR.FDRLevel m_optimizeWhat; 
//    

    public static class CustomSpinnerModel extends AbstractSpinnerModel {


        private double value = 100.0;
        private double increment = 1;
        private double max = 100;
        private double min  = 0;
                

        public CustomSpinnerModel(){
             this.value=value;
        }

        public CustomSpinnerModel(Double value){
             this.value=value;
        }

        @Override
        public Object getNextValue() {
            setValue(Math.round((value+increment)/increment)*increment);
            return value+"";//return as string to avoid round
        }

        @Override
        public Object getPreviousValue() {
            setValue(Math.round((value-increment)/increment)*increment);
            return value+"";//return as string to avoid round
        }

        @Override
        public Object getValue() {
            if (value == (int)value) {
                return ((int) value) + "";
            }
            return value+"";//return as string to avoid round
        }

        @Override
        public void setValue(Object o) {
            double newValue = Double.NaN;
            if (o instanceof Number) {
                newValue = ((Number)o).doubleValue();
            } else {
                try{
                   newValue=Double.parseDouble(o.toString());
                } catch(Exception e){
                   newValue = value;
                }
            }
            if (newValue < min)
                newValue = min;
            if (newValue>max) 
                newValue = max;
            String svalue = value+"";
            if (svalue.contains("00000")) {
                newValue = Double.parseDouble(svalue.substring(0,svalue.indexOf("0000")));
            }
            if (newValue != value) {
                value = newValue;
                fireStateChanged();
            } else if (newValue+"" != o.toString()) {
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
    
    /**
     * Creates new form FDRSettingsComplete
     */
    public FDRSettingsComplete() {
        initComponents();
        ChangeListener max100Listener = new ChangeListener() {
            public void stateChanged(ChangeEvent e) {
                SpinnerModel sp = (SpinnerModel) e.getSource();
                if (((Double)sp.getValue()) >100)
                    sp.setValue(100d);
            }
        };    
        setSpinnerModel(spPsmFDR, 100);
        setSpinnerModel(spPepFDR, 100);
        setSpinnerModel(spProteinFDR, 100);
        setSpinnerModel(spLinkFDR, 5);
        setSpinnerModel(spPPIFdr, 100);
        
        
        spMinProteinPepCount.setSpecialValue(1);
        spMinLinkPepCount.setSpecialValue(1);
        spMinPPIPepCount.setSpecialValue(1);
        spReportFactor.setVisible(false);
        lblReportFactor.setVisible(false);
        setMinTD(DEFAULT_MIN_TD_COUNT);
        spOtherFilter.setVisible(ckMoreOptions.isSelected());
    }

    public static void setSpinnerModel(JSpinner sp, double intialValue) {
        sp.setModel(new CustomSpinnerModel(intialValue));
        JFormattedTextField tf = ((JSpinner.DefaultEditor) sp.getEditor()).getTextField();
        tf.setEditable(true);
        tf.setHorizontalAlignment(SwingConstants.RIGHT);
    }
    
    
    @Override
    public boolean getBoostBetween() {
        return ckBoostBetween.isSelected();
    }

    @Override
    public void setBoostBetween(final boolean between) {
        SwingUtilities.invokeLater(new Runnable() {

            public void run() {
                ckBoostBetween.setSelected(between);
            }
        });
    }


    
    private void setValueLater(final JSpinner sp, final Object value) {
        SwingUtilities.invokeLater(new Runnable() {

            public void run() {
                sp.setValue(value);
            }
        });
            
    }
    
    private void setValueLater(final JCheckBox ck, final boolean value) {
        SwingUtilities.invokeLater(new Runnable() {

            public void run() {
                ck.setSelected(value);
            }
        });
    }      
        
    
    private void doCalc() {
        btnStopBoost.setEnabled(ckMaximize.isSelected());        
        raiseStartCalc(ckMaximize.isSelected());
        
//        if (ckMaximize.isSelected()) {
//            for (ActionListener al : m_calc_listener) {
//                al.actionPerformed(new ActionEvent(this, 0, "boost"));
//            }
//        } else {
//            for (ActionListener al : m_calc_listener) {
//                al.actionPerformed(new ActionEvent(this, -1, "calc"));
//            }
//        }
    }    
    
    @Override
    public double getPSMFDR() {
        return Double.parseDouble(spPsmFDR.getValue().toString())/100.0;
    }

    
    @Override
    public double getPeptidePairFDR() {
        return Double.parseDouble(spPepFDR.getValue().toString())/100.0;
    }

    @Override
    public double getProteinGroupFDR() {
        return Double.parseDouble(spProteinFDR.getValue().toString())/100.0;
    }

    @Override
    public double getProteinGroupLinkFDR() {
        return Double.parseDouble(spLinkFDR.getValue().toString())/100.0;
    }

    @Override
    public double getProteinGroupPairFDR() {
        return Double.parseDouble(spPPIFdr.getValue().toString())/100.0;
    }

    @Override
    public void setPSMFDR( Double fdr) {
        setValueLater(spPsmFDR,fdr*100);
    }

    
    @Override
    public void setPeptidePairFDR( Double fdr) {
        setValueLater(spPepFDR,fdr*100);
    }

    @Override
    public void setProteinGroupFDR( Double fdr) {
        setValueLater(spProteinFDR,fdr*100);
    }

    @Override
    public void setProteinGroupLinkFDR( Double fdr) {
        setValueLater(spLinkFDR,fdr*100);
    }

    @Override
    public void setProteinGroupPairFDR( Double fdr) {
        setValueLater(spPPIFdr,fdr*100);
    }

    @Override
    public int getMinProteinPepCount() {
        return (Integer) spMinProteinPepCount.getValue();
    }

    @Override
    public int getMinLinkPepCount() {
        return (Integer) spMinLinkPepCount.getValue();
    }

    @Override
    public int getMinPPIPepCount() {
        return (Integer) spMinPPIPepCount.getValue();
    }

    @Override
    public void setMinProteinPepCount(Integer minPep) {
        setValueLater(spMinProteinPepCount,minPep);
    }

    @Override
    public void setMinLinkPepCount(Integer minPep) {
        setValueLater(spMinLinkPepCount,minPep);
    }

    @Override
    public void setMinPPIPepCount(Integer minPep) {
        setValueLater(spMinPPIPepCount,minPep);
    }

    
    

    @Override
    public int getMinPeptideLength() {
        return (Integer) spMinPeptideLength.getValue();
    }

    @Override
    public void setMinPeptideLength(Integer minLength) {
        setValueLater(spMinPeptideLength,minLength);
    }

    
    @Override
    public int getMaxLinkAmbiguity() {
        return (Integer) spMaxLinkAmbiguity.getValue();
    }

    @Override
    public void setMaxLinkAmbiguity(Integer maxAmbiguity) {
        setValueLater(spMaxLinkAmbiguity,maxAmbiguity);
    }

    @Override
    public int getMaxProteinAmbiguity() {
        return (Integer) spMaxProteinAmbiguity.getValue();
    }

    @Override
    public void setMaxProteinAmbiguity(Integer maxAmbiguity) {
        setValueLater(spMaxProteinAmbiguity,maxAmbiguity);
    }
    
    @Override
    public boolean isPSMDirectional() {
        return false;
    }
    
    @Override
    public boolean isPeptidePairDirectional() {
        return false;
    }

    @Override
    public boolean isLinkDirectional() {
        return false;
    }

    @Override
    public boolean isPPIDirectional() {
        return false;
    }

    @Override
    public void setPSMDirectional(boolean directional) {
        if (directional)
            throw new UnsupportedOperationException("Can't switch to directional");
    }
    
    @Override
    public void setPeptidePairDirectional(boolean directional) {
        if (directional)
            throw new UnsupportedOperationException("Can't switch to directional");
    }

    @Override
    public void setLinkDirectional(boolean directional) {
        if (directional)
            throw new UnsupportedOperationException("Can't switch to directional");
    }

    @Override
    public void setPPIDirectional(boolean directional) {
        if (directional)
            throw new UnsupportedOperationException("Can't switch to directional");
    }
    
    
    @Override
    public OfflineFDR.FDRLevel doOptimize() {
        if (!ckMaximize.isSelected())
            return null;
       
        return (OfflineFDR.FDRLevel) cbBoostWhat.getSelectedItem();
    }
    
    @Override
    public void doOptimize(OfflineFDR.FDRLevel level) {
        if (level == null)
            ckMaximize.setSelected(false);
        else { 
            ckMaximize.setSelected(true);
            cbBoostWhat.getModel().setSelectedItem(level);
        }
    }

    public int getBoostingSteps() {
        return (Integer) spMaximizeSteps.getValue();
    }

    public void setBoostingSteps(int steps) {
        spMaximizeSteps.setValue(steps);
    }
    
    public double getReportFactor() {
        return (Double) spReportFactor.getValue();
    }       
    
    public void setReportFactor(double factor) {
        spReportFactor.setValue(factor);
    }       
    
    @Override
    public void setEnabled(boolean e) {
        super.setEnabled(e);
        btnCalc.setEnabled(e);
        if (e)
            btnStopBoost.setEnabled(false);        
    }    

    public boolean filterToUniquePSM() {
        return ckUniquePSM.isSelected();
    }

    public void setFilterToUniquePSM(boolean filterToUniquePSM) {
        ckUniquePSM.setSelected(filterToUniquePSM);
    }
    
    public void setMinTD(Integer c) {
        spMinTDChance.setValue(c);
    }

    @Override
    public int getMinTD() {
        return (Integer) spMinTDChance.getValue();
    }

    @Override
    public Boolean psmLocalFDR() {
        return this.tsLocalPSMFDR.getSelectionState();
    }

    @Override
    public Boolean peppairLocalFDR() {
        return this.tsLocalPepPairFDR.getSelectionState();
    }

    @Override
    public Boolean protLocalFDR() {
        return this.tsLocalProtFDR.getSelectionState();
    }

    @Override
    public Boolean linkLocalFDR() {
        return this.tsLocalLinkFDR.getSelectionState();
    }

    @Override
    public Boolean ppiLocalFDR() {
        return this.tsLocalPPIFDR.getSelectionState();
    }

    @Override
    public void psmLocalFDR(Boolean local) {
        this.tsLocalPSMFDR.setSelected(local);
    }

    @Override
    public void peppairLocalFDR(Boolean local) {
        this.tsLocalPepPairFDR.setSelected(local);
    }

    @Override
    public void protLocalFDR(Boolean local) {
        this.tsLocalProtFDR.setSelected(local);
    }

    @Override
    public void linkLocalFDR(Boolean local) {
        this.tsLocalLinkFDR.setSelected(local);
    }

    @Override
    public void ppiLocalFDR(Boolean local) {
        this.tsLocalPPIFDR.setSelected(local);
    }    

    public double getMinPeptideCoverageFilter() {
        return this.otherFilter.getMinPepCoverage();
    }
    public void setMinPeptideCoverageFilter(double d) {
        this.otherFilter.setMinPepCoverage(d);
    }

    public double getMinDeltaScoreFilter() {
        return this.otherFilter.getMinDeltaScore();
    }
    
    public void setMinDeltaScoreFilter(double d) {
        this.otherFilter.setMinDeltaScore(d);
    }

    
    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
        /* Set the Nimbus look and feel */
        //<editor-fold defaultstate="collapsed" desc=" Look and feel setting code (optional) ">
        /* If Nimbus (introduced in Java SE 6) is not available, stay with the default look and feel.
         * For details see http://download.oracle.com/javase/tutorial/uiswing/lookandfeel/plaf.html 
         */
        try {
            for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels()) {
                if ("Nimbus".equals(info.getName())) {
                    javax.swing.UIManager.setLookAndFeel(info.getClassName());
                    break;
                }
            }
        } catch (ClassNotFoundException ex) {
            java.util.logging.Logger.getLogger(FDRSettingsComplete.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (InstantiationException ex) {
            java.util.logging.Logger.getLogger(FDRSettingsComplete.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (IllegalAccessException ex) {
            java.util.logging.Logger.getLogger(FDRSettingsComplete.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (javax.swing.UnsupportedLookAndFeelException ex) {
            java.util.logging.Logger.getLogger(FDRSettingsComplete.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>

        /* Create and display the form */
        java.awt.EventQueue.invokeLater(new Runnable() {
            public void run() {
                new FDRSettingsComplete().setVisible(true);
            }
        });
    }    
    
    public boolean combineScoreAndDelta() {
        return this.otherFilter.combineScoreAndDelta();
    }
    public void combineScoreAndDelta(boolean c) {
        this.otherFilter.combineScoreAndDelta(c);
    }

    @Override
    public int getMinPeptideFragmentsFilter() {
        return Integer.parseInt(this.otherFilter.txtMinPepFrags.getText());
    }

    @Override
    public void setMinPeptideFragmentsFilter(int frags) {
        this.otherFilter.txtMinPepFrags.setText(Integer.toString(frags));
    }    

    @Override
    public boolean ignoreValidityChecks() {
        return this.ckIgnoreValidity.isSelected();
    }

    @Override
    public void ignoreValidityChecks(boolean ignore) {
        ckIgnoreValidity.setSelected(ignore);
    }


    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jLabel5 = new javax.swing.JLabel();
        spPsmFDR = new javax.swing.JSpinner();
        jLabel2 = new javax.swing.JLabel();
        spPepFDR = new javax.swing.JSpinner();
        jLabel6 = new javax.swing.JLabel();
        spProteinFDR = new javax.swing.JSpinner();
        jLabel3 = new javax.swing.JLabel();
        spLinkFDR = new javax.swing.JSpinner();
        jLabel4 = new javax.swing.JLabel();
        spPPIFdr = new javax.swing.JSpinner();
        lblReportFactor = new javax.swing.JLabel();
        spReportFactor = new javax.swing.JSpinner();
        jLabel8 = new javax.swing.JLabel();
        jLabel28 = new javax.swing.JLabel();
        spMaxLinkAmbiguity = new org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner();
        spMaxProteinAmbiguity = new org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner();
        jLabel1 = new javax.swing.JLabel();
        spMinPeptideLength = new org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner();
        jLabel18 = new javax.swing.JLabel();
        jLabel9 = new javax.swing.JLabel();
        btnCalc = new javax.swing.JButton();
        ckMaximize = new javax.swing.JCheckBox();
        spMaximizeSteps = new javax.swing.JSpinner();
        jLabel23 = new javax.swing.JLabel();
        cbBoostWhat = new org.rappsilber.fdr.gui.components.FDRLevelComboBox();
        spMinProteinPepCount = new org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner();
        spMinLinkPepCount = new org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner();
        spMinPPIPepCount = new org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner();
        btnStopBoost = new javax.swing.JButton();
        ckUniquePSM = new javax.swing.JCheckBox();
        ckBoostBetween = new javax.swing.JCheckBox();
        btnBoostIgnores = new javax.swing.JButton();
        ckGroupByPSMs = new javax.swing.JCheckBox();
        spMinTDChance = new org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner();
        jLabel7 = new javax.swing.JLabel();
        ckFilterConsecutives = new javax.swing.JCheckBox();
        tsLocalPSMFDR = new org.rappsilber.gui.components.TriStateCheckBox();
        tsLocalPepPairFDR = new org.rappsilber.gui.components.TriStateCheckBox();
        tsLocalProtFDR = new org.rappsilber.gui.components.TriStateCheckBox();
        tsLocalLinkFDR = new org.rappsilber.gui.components.TriStateCheckBox();
        tsLocalPPIFDR = new org.rappsilber.gui.components.TriStateCheckBox();
        ckMoreOptions = new javax.swing.JCheckBox();
        spOtherFilter = new javax.swing.JScrollPane();
        otherFilter = new org.rappsilber.fdr.gui.components.OtherFilter();
        ckIgnoreValidity = new javax.swing.JCheckBox();

        jLabel5.setText("PSM");

        spPsmFDR.setToolTipText("FDR value accepted for PSMs");

        jLabel2.setText("Peptide Pair");

        spPepFDR.setToolTipText("FDR value accepted for Peptide Pairs (including linksite within the peptide)");

        jLabel6.setText("Protein Group");

        spProteinFDR.setToolTipText("FDR value accepted for protein-groups");

        jLabel3.setText("Residue Pairs");

        spLinkFDR.setToolTipText("FDR value accepted for Links (Protein-group-links)");

        jLabel4.setText("Protein Pairs");

        spPPIFdr.setToolTipText("FDR value accepted for protein-pairs");

        lblReportFactor.setText("Report Factor");
        lblReportFactor.setToolTipText("maximum factor the next step in fdr is permited to exced the target fdr");

        spReportFactor.setModel(new javax.swing.SpinnerNumberModel(100000.0d, 1.0d, null, 0.1d));
        spReportFactor.setToolTipText("maximum factor the next step in fdr is permited to exced the target fdr");

        jLabel8.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
        jLabel8.setText("Ambiguity");

        jLabel28.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
        jLabel28.setText("Min supporting peptide-pairs");

        spMaxLinkAmbiguity.setToolTipText("maximum ambiguity for a peptide-pair giving raise to a residue pair");

        spMaxProteinAmbiguity.setToolTipText("maximum number of proteins a peptide could originate from (this ignores if a peptide could come from several places within a protein)");

        jLabel1.setText("Local FDR");
        jLabel1.setToolTipText("ticked calculate and filtert by PEP; unticked: do not calculate calculate PEP; square: calculate PEP but do not filter (default)");

        spMinPeptideLength.setToolTipText("matches to involving\u0000 peptides shorter then this will be ignored");

        jLabel18.setHorizontalAlignment(javax.swing.SwingConstants.RIGHT);
        jLabel18.setText("Min Pep. Length:");

        jLabel9.setText("Max FDRs");

        btnCalc.setText("Calculate");
        btnCalc.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnCalcActionPerformed(evt);
            }
        });

        ckMaximize.setSelected(true);
        ckMaximize.setText("boost");
        ckMaximize.setToolTipText("should we try to boost the results reported for a specified maximum FDR");
        ckMaximize.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                ckMaximizeActionPerformed(evt);
            }
        });

        spMaximizeSteps.setModel(new javax.swing.SpinnerNumberModel(4, 2, 20, 1));
        spMaximizeSteps.setToolTipText("in how many steps is each FDR-level checked for each round of optimisations");

        jLabel23.setHorizontalAlignment(javax.swing.SwingConstants.RIGHT);
        jLabel23.setText("Steps:");

        cbBoostWhat.setSelectedIndex(3);
        cbBoostWhat.setToolTipText("Boost the results for the given level by modifying the lower-level FDRs");

        spMinProteinPepCount.setModel(new javax.swing.SpinnerNumberModel(1, 1, null, 1));
        spMinProteinPepCount.setSpecialValue(1);

        spMinLinkPepCount.setModel(new javax.swing.SpinnerNumberModel(1, 1, null, 1));
        spMinLinkPepCount.setSpecialValue(1);

        spMinPPIPepCount.setModel(new javax.swing.SpinnerNumberModel(1, 1, null, 1));
        spMinPPIPepCount.setSpecialValue(1);

        btnStopBoost.setText("stop");
        btnStopBoost.setEnabled(false);
        btnStopBoost.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnStopBoostActionPerformed(evt);
            }
        });

        ckUniquePSM.setText("Unique PSMs");
        ckUniquePSM.setToolTipText("only use the top match for each combination of charge peptide1, peptide2 and charge state as a PSM");

        ckBoostBetween.setText("Between");
        ckBoostBetween.setToolTipText("optimize settings for maximizing between matches");

        btnBoostIgnores.setText("Boost Includes");
        btnBoostIgnores.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnBoostIgnoresActionPerformed(evt);
            }
        });

        ckGroupByPSMs.setText("Group by PSMs");
        ckGroupByPSMs.setToolTipText("Groups  matches baessed on how many psms support a given protein pair.");
        ckGroupByPSMs.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                ckGroupByPSMsActionPerformed(evt);
            }
        });

        spMinTDChance.setModel(new javax.swing.SpinnerNumberModel(0, 0, null, 1));
        spMinTDChance.setToolTipText("The total number of matches times the fdr should be larger then the given number - otherwise the subgroup will be considered unreliable");
        spMinTDChance.setSpecialValue(0);

        jLabel7.setText("Min TD Chance");

        ckFilterConsecutives.setText("no consecutive");
        ckFilterConsecutives.setToolTipText("Filter out cross-linked PSMs to consecutive peptides");
        ckFilterConsecutives.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                ckFilterConsecutivesActionPerformed(evt);
            }
        });

        tsLocalPSMFDR.setToolTipText("calculate (PEP) (square); calculate and filter by PEP (tick); or do not calculate PEP");

        tsLocalPepPairFDR.setToolTipText("calculate (PEP) (square); calculate and filter by PEP (tick); or do not calculate PEP");

        tsLocalProtFDR.setToolTipText("calculate (PEP) (square); calculate and filter by PEP (tick); or do not calculate PEP");

        tsLocalLinkFDR.setToolTipText("calculate (PEP) (square); calculate and filter by PEP (tick); or do not calculate PEP");

        tsLocalPPIFDR.setToolTipText("calculate (PEP) (square); calculate and filter by PEP (tick); or do not calculate PEP");

        ckMoreOptions.setText("More");
        ckMoreOptions.setToolTipText("some addtional optional prefilter");
        ckMoreOptions.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                ckMoreOptionsActionPerformed(evt);
            }
        });

        spOtherFilter.setViewportView(otherFilter);

        ckIgnoreValidity.setText("Ignore Validity Checks");
        ckIgnoreValidity.setToolTipText("Sub-groups with e.g. negative calculated FDR or to few target matches for the requested FDR)");
        ckIgnoreValidity.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                ckIgnoreValidityActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(spOtherFilter)
                        .addContainerGap())
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(layout.createSequentialGroup()
                                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(jLabel5)
                                    .addComponent(jLabel2)
                                    .addComponent(jLabel6)
                                    .addComponent(jLabel3)
                                    .addComponent(jLabel4)
                                    .addComponent(lblReportFactor))
                                .addGap(12, 12, 12)
                                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                                    .addComponent(spReportFactor, javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(spPPIFdr, javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(spLinkFDR, javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(spProteinFDR, javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(spPepFDR, javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(spPsmFDR, javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(jLabel9, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, 92, Short.MAX_VALUE))
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(jLabel1)
                                    .addComponent(jLabel7)
                                    .addGroup(layout.createSequentialGroup()
                                        .addGap(25, 25, 25)
                                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                            .addComponent(tsLocalPPIFDR, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                            .addComponent(tsLocalLinkFDR, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                            .addComponent(tsLocalPepPairFDR, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                            .addComponent(tsLocalPSMFDR, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                            .addComponent(tsLocalProtFDR, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))))
                            .addComponent(ckMoreOptions))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(ckIgnoreValidity)
                                .addGap(0, 0, Short.MAX_VALUE))
                            .addComponent(spMinPPIPepCount, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(spMinLinkPepCount, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(spMinProteinPepCount, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(jLabel28, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(jLabel18, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(spMinTDChance, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(spMaxLinkAmbiguity, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(jLabel8, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(spMinPeptideLength, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(spMaxProteinAmbiguity, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addGroup(layout.createSequentialGroup()
                                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(ckGroupByPSMs)
                                    .addComponent(ckUniquePSM)
                                    .addComponent(ckFilterConsecutives))
                                .addGap(0, 0, Short.MAX_VALUE))))
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(ckMaximize)
                                .addGap(12, 12, 12)
                                .addComponent(cbBoostWhat, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                            .addComponent(btnStopBoost))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(jLabel23)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(layout.createSequentialGroup()
                                .addGap(70, 70, 70)
                                .addComponent(btnBoostIgnores)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                .addComponent(btnCalc))
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(spMaximizeSteps)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                                .addComponent(ckBoostBetween, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))))))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGap(19, 19, 19)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addGroup(layout.createSequentialGroup()
                                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                                    .addGroup(layout.createSequentialGroup()
                                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                            .addComponent(jLabel9)
                                            .addComponent(jLabel1))
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                                                .addComponent(spPsmFDR, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                                .addComponent(jLabel5))
                                            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                                                .addComponent(spMinPeptideLength, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                                .addComponent(jLabel18)
                                                .addComponent(tsLocalPSMFDR, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                                        .addGap(6, 6, 6)
                                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                                                .addComponent(spPepFDR, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                                .addComponent(jLabel2))
                                            .addComponent(jLabel28, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.PREFERRED_SIZE, 19, javax.swing.GroupLayout.PREFERRED_SIZE)
                                            .addComponent(jLabel8, javax.swing.GroupLayout.Alignment.TRAILING)))
                                    .addComponent(tsLocalPepPairFDR, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                                .addGap(6, 6, 6)
                                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                                    .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                                        .addComponent(spProteinFDR, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                        .addComponent(jLabel6))
                                    .addComponent(spMinProteinPepCount, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                    .addComponent(spMaxProteinAmbiguity, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                            .addComponent(tsLocalProtFDR, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addGap(6, 6, 6)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(spMaxLinkAmbiguity, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(spMinLinkPepCount, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                                .addComponent(spLinkFDR, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addComponent(jLabel3))))
                    .addComponent(tsLocalLinkFDR, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        .addComponent(spMinPPIPepCount, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addComponent(ckUniquePSM))
                    .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        .addComponent(spPPIFdr, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addComponent(jLabel4))
                    .addComponent(tsLocalPPIFDR, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(spReportFactor, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(lblReportFactor)
                    .addComponent(ckGroupByPSMs)
                    .addComponent(spMinTDChance, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel7))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(ckFilterConsecutives)
                    .addComponent(ckMoreOptions)
                    .addComponent(ckIgnoreValidity))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(spOtherFilter, javax.swing.GroupLayout.DEFAULT_SIZE, 96, Short.MAX_VALUE)
                .addGap(18, 18, 18)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(ckMaximize)
                    .addComponent(cbBoostWhat, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel23)
                    .addComponent(spMaximizeSteps, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(ckBoostBetween))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(btnStopBoost)
                    .addComponent(btnCalc)
                    .addComponent(btnBoostIgnores))
                .addContainerGap())
        );
    }// </editor-fold>//GEN-END:initComponents

    private void btnCalcActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnCalcActionPerformed
        doCalc();
    }//GEN-LAST:event_btnCalcActionPerformed

    private void btnStopBoostActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnStopBoostActionPerformed
        raiseStopMaximizing();
    }//GEN-LAST:event_btnStopBoostActionPerformed

    private void ckMaximizeActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_ckMaximizeActionPerformed
                ckBoostBetween.setEnabled(ckMaximize.isSelected());
                cbBoostWhat.setEnabled(ckMaximize.isSelected());
                spMaximizeSteps.setEnabled(ckMaximize.isSelected());
    }//GEN-LAST:event_ckMaximizeActionPerformed

    private void btnBoostIgnoresActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnBoostIgnoresActionPerformed
        setBoostIgnores();
    }//GEN-LAST:event_btnBoostIgnoresActionPerformed

    private void ckGroupByPSMsActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_ckGroupByPSMsActionPerformed
        setGroupByPSMCount(ckGroupByPSMs.isSelected());
        if (ckGroupByPSMs.isSelected() && ((Integer)spMinTDChance.getValue())==0 ) {
            spMinTDChance.setValue(1);
        }
    }//GEN-LAST:event_ckGroupByPSMsActionPerformed

    private void ckFilterConsecutivesActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_ckFilterConsecutivesActionPerformed
        setFilterConsecutivePeptides(ckFilterConsecutives.isSelected());
    }//GEN-LAST:event_ckFilterConsecutivesActionPerformed

    private void ckMoreOptionsActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_ckMoreOptionsActionPerformed
        spOtherFilter.setVisible(ckMoreOptions.isSelected());
        this.repaint();
        this.revalidate();
    }//GEN-LAST:event_ckMoreOptionsActionPerformed

    private void ckIgnoreValidityActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_ckIgnoreValidityActionPerformed
        //spMinTDChance.setEnabled(!ckIgnoreValidity.isSelected());
    }//GEN-LAST:event_ckIgnoreValidityActionPerformed

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton btnBoostIgnores;
    private javax.swing.JButton btnCalc;
    public javax.swing.JButton btnStopBoost;
    private org.rappsilber.fdr.gui.components.FDRLevelComboBox cbBoostWhat;
    private javax.swing.JCheckBox ckBoostBetween;
    private javax.swing.JCheckBox ckFilterConsecutives;
    private javax.swing.JCheckBox ckGroupByPSMs;
    private javax.swing.JCheckBox ckIgnoreValidity;
    private javax.swing.JCheckBox ckMaximize;
    private javax.swing.JCheckBox ckMoreOptions;
    private javax.swing.JCheckBox ckUniquePSM;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel18;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel23;
    private javax.swing.JLabel jLabel28;
    private javax.swing.JLabel jLabel3;
    private javax.swing.JLabel jLabel4;
    private javax.swing.JLabel jLabel5;
    private javax.swing.JLabel jLabel6;
    private javax.swing.JLabel jLabel7;
    private javax.swing.JLabel jLabel8;
    private javax.swing.JLabel jLabel9;
    public javax.swing.JLabel lblReportFactor;
    private org.rappsilber.fdr.gui.components.OtherFilter otherFilter;
    private javax.swing.JSpinner spLinkFDR;
    private org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner spMaxLinkAmbiguity;
    private org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner spMaxProteinAmbiguity;
    private javax.swing.JSpinner spMaximizeSteps;
    private org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner spMinLinkPepCount;
    private org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner spMinPPIPepCount;
    private org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner spMinPeptideLength;
    private org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner spMinProteinPepCount;
    private org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner spMinTDChance;
    private javax.swing.JScrollPane spOtherFilter;
    private javax.swing.JSpinner spPPIFdr;
    private javax.swing.JSpinner spPepFDR;
    private javax.swing.JSpinner spProteinFDR;
    private javax.swing.JSpinner spPsmFDR;
    public javax.swing.JSpinner spReportFactor;
    private org.rappsilber.gui.components.TriStateCheckBox tsLocalLinkFDR;
    private org.rappsilber.gui.components.TriStateCheckBox tsLocalPPIFDR;
    private org.rappsilber.gui.components.TriStateCheckBox tsLocalPSMFDR;
    private org.rappsilber.gui.components.TriStateCheckBox tsLocalPepPairFDR;
    private org.rappsilber.gui.components.TriStateCheckBox tsLocalProtFDR;
    // End of variables declaration//GEN-END:variables




}
