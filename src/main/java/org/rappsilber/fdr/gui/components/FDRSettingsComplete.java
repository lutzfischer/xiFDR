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

import javax.swing.JCheckBox;
import javax.swing.JSpinner;
import javax.swing.SpinnerModel;
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
        spPepFDR.getModel().addChangeListener(max100Listener);
        spProteinFDR.getModel().addChangeListener(max100Listener);
        spLinkFDR.getModel().addChangeListener(max100Listener);
        spPsmFDR.getModel().addChangeListener(max100Listener);
        spPPIFdr.getModel().addChangeListener(max100Listener);
        
        spMinProteinPepCount.setSpecialValue(1);
        spMinLinkPepCount.setSpecialValue(1);
        spMinPPIPepCount.setSpecialValue(1);
        
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
        return ((Double) spPsmFDR.getValue())/100.0;
    }

    
    @Override
    public double getPeptidePairFDR() {
        return ((Double) spPepFDR.getValue())/100.0;
    }

    @Override
    public double getProteinGroupFDR() {
        return ((Double) spProteinFDR.getValue())/100.0;
    }

    @Override
    public double getProteinGroupLinkFDR() {
        return ((Double) spLinkFDR.getValue())/100.0;
    }

    @Override
    public double getProteinGroupPairFDR() {
        return ((Double) spPPIFdr.getValue())/100.0;
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
        return ckDirectionalPSM.isSelected();
    }
    
    @Override
    public boolean isPeptidePairDirectional() {
        return ckDirectionalPeptidePair.isSelected();
    }

    @Override
    public boolean isLinkDirectional() {
        return ckDirectionalLink.isSelected();
    }

    @Override
    public boolean isPPIDirectional() {
        return ckDirectionalPPI.isSelected();
    }

    @Override
    public void setPSMDirectional(boolean directional) {
        setValueLater(ckDirectionalPSM,directional);
    }
    
    @Override
    public void setPeptidePairDirectional(boolean directional) {
        setValueLater(ckDirectionalPeptidePair,directional);
    }

    @Override
    public void setLinkDirectional(boolean directional) {
        setValueLater(ckDirectionalLink,directional);
    }

    @Override
    public void setPPIDirectional(boolean directional) {
        setValueLater(ckDirectionalPPI,directional);
    }
    
    
    @Override
    public OfflineFDR.FDRLevel doOptimize() {
        if (!ckMaximize.isSelected())
            return null;
       
        return (OfflineFDR.FDRLevel) cbBoostWhat.getSelectedItem();
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
    
    public void setAll(FDRSettings settings) {
        this.setBoostingSteps(settings.getBoostingSteps());
        this.setMaxLinkAmbiguity(settings.getMaxLinkAmbiguity());
        this.setMaxProteinAmbiguity(settings.getMaxProteinAmbiguity());
        this.setMinLinkPepCount(settings.getMinLinkPepCount());
        this.setMinPPIPepCount(settings.getMinPPIPepCount());
        this.setMinPeptideLength(settings.getMinPeptideLength());
        this.setMinProteinPepCount(settings.getMinProteinPepCount());
        this.setPSMFDR(settings.getPSMFDR());
        this.setPSMDirectional(settings.isPSMDirectional());
        this.setPPIDirectional(settings.isPPIDirectional());
        this.setLinkDirectional(settings.isLinkDirectional());
        this.setPeptidePairDirectional(settings.isPeptidePairDirectional());
        this.setPeptidePairFDR(settings.getPeptidePairFDR());
        this.setProteinGroupFDR(settings.getProteinGroupFDR());
        this.setProteinGroupLinkFDR(settings.getProteinGroupLinkFDR());
        this.setProteinGroupPairFDR(settings.getProteinGroupPairFDR());
        this.setReportFactor(settings.getReportFactor());
        cbBoostWhat.getModel().setSelectedItem(settings.doOptimize());
        this.setFilterToUniquePSM(settings.filterToUniquePSM());
        this.setBoostBetween(settings.getBoostBetween());
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
        jLabel7 = new javax.swing.JLabel();
        spReportFactor = new javax.swing.JSpinner();
        jLabel8 = new javax.swing.JLabel();
        jLabel28 = new javax.swing.JLabel();
        spMaxLinkAmbiguity = new org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner();
        spMaxProteinAmbiguity = new org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner();
        ckDirectionalPeptidePair = new javax.swing.JCheckBox();
        ckDirectionalLink = new javax.swing.JCheckBox();
        ckDirectionalPPI = new javax.swing.JCheckBox();
        jLabel1 = new javax.swing.JLabel();
        ckDirectionalPSM = new javax.swing.JCheckBox();
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

        jLabel5.setText("PSM");

        spPsmFDR.setModel(new javax.swing.SpinnerNumberModel(Double.valueOf(100.0d), Double.valueOf(0.0d), null, Double.valueOf(1.0d)));
        spPsmFDR.setToolTipText("FDR value accepted for PSMs");

        jLabel2.setText("Peptide Pair");

        spPepFDR.setModel(new javax.swing.SpinnerNumberModel(Double.valueOf(100.0d), Double.valueOf(0.0d), null, Double.valueOf(1.0d)));
        spPepFDR.setToolTipText("FDR value accepted for Peptide Pairs (including linksite within the peptide)");

        jLabel6.setText("Protein Group");

        spProteinFDR.setModel(new javax.swing.SpinnerNumberModel(Double.valueOf(100.0d), Double.valueOf(0.0d), null, Double.valueOf(1.0d)));
        spProteinFDR.setToolTipText("FDR value accepted for protein-groups");

        jLabel3.setText("Link");

        spLinkFDR.setModel(new javax.swing.SpinnerNumberModel(Double.valueOf(5.0d), Double.valueOf(0.0d), null, Double.valueOf(1.0d)));
        spLinkFDR.setToolTipText("FDR value accepted for Links (Protein-group-links)");

        jLabel4.setText("Protein Pairs");

        spPPIFdr.setModel(new javax.swing.SpinnerNumberModel(Double.valueOf(100.0d), Double.valueOf(0.0d), null, Double.valueOf(1.0d)));
        spPPIFdr.setToolTipText("FDR value accepted for protein-pairs");

        jLabel7.setText("Report Factor");
        jLabel7.setToolTipText("maximum factor the next step in fdr is permited to exced the target fdr");

        spReportFactor.setModel(new javax.swing.SpinnerNumberModel(Double.valueOf(1.1d), Double.valueOf(1.0d), null, Double.valueOf(0.1d)));
        spReportFactor.setToolTipText("maximum factor the next step in fdr is permited to exced the target fdr");

        jLabel8.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
        jLabel8.setText("Ambiguity");

        jLabel28.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
        jLabel28.setText("Min supporting peptide-pairs");

        ckDirectionalPeptidePair.setToolTipText("Are peptide-pairs directional (E.g. is peptide1 linked to peptide2 distinct from peptide 2 linked to peptide1)");

        ckDirectionalLink.setToolTipText("Are residue-pairs directional (E.g. is residue 1 linked to residue 2 distinct from residue 2 linked to residue 1)");

        ckDirectionalPPI.setToolTipText("Are protein pairs directional (E.g. is protein 1 linked to protein 2 distinct from protein 2 linked to protein 1)");

        jLabel1.setText("Directional");

        ckDirectionalPSM.setToolTipText("Are PSM-matches directional (E.g. is peptide1 linked to peptide2 distinct from peptide 2 linked to peptide1)");

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

        spMinProteinPepCount.setModel(new javax.swing.SpinnerNumberModel(Integer.valueOf(1), Integer.valueOf(1), null, Integer.valueOf(1)));
        spMinProteinPepCount.setSpecialValue(1);

        spMinLinkPepCount.setModel(new javax.swing.SpinnerNumberModel(Integer.valueOf(1), Integer.valueOf(1), null, Integer.valueOf(1)));
        spMinLinkPepCount.setSpecialValue(1);

        spMinPPIPepCount.setModel(new javax.swing.SpinnerNumberModel(Integer.valueOf(1), Integer.valueOf(1), null, Integer.valueOf(1)));
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

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jLabel5)
                            .addComponent(jLabel2)
                            .addComponent(jLabel6)
                            .addComponent(jLabel3)
                            .addComponent(jLabel4)
                            .addComponent(jLabel7)
                            .addComponent(ckMaximize))
                        .addGap(12, 12, 12)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                            .addComponent(spReportFactor, javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(spPPIFdr, javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(spLinkFDR, javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(spProteinFDR, javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(spPepFDR, javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(spPsmFDR, javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jLabel9, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, 92, Short.MAX_VALUE)
                            .addComponent(cbBoostWhat, javax.swing.GroupLayout.PREFERRED_SIZE, 1, Short.MAX_VALUE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jLabel23)
                            .addComponent(ckDirectionalPPI)
                            .addComponent(ckDirectionalLink)
                            .addComponent(ckDirectionalPeptidePair)
                            .addComponent(ckDirectionalPSM)
                            .addComponent(jLabel1))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(spMaximizeSteps)
                            .addComponent(spMinPPIPepCount, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(spMinLinkPepCount, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(spMinProteinPepCount, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(jLabel28, javax.swing.GroupLayout.PREFERRED_SIZE, 0, Short.MAX_VALUE)
                            .addComponent(jLabel18, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(spMaxLinkAmbiguity, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(jLabel8, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(spMinPeptideLength, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(spMaxProteinAmbiguity, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(ckUniquePSM)
                                .addGap(0, 55, Short.MAX_VALUE))
                            .addComponent(ckBoostBetween, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(btnStopBoost)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(btnCalc))))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGap(19, 19, 19)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jLabel9)
                    .addComponent(jLabel1))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        .addComponent(spPsmFDR, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addComponent(jLabel5))
                    .addComponent(ckDirectionalPSM, javax.swing.GroupLayout.Alignment.TRAILING)
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        .addComponent(spMinPeptideLength, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addComponent(jLabel18)))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        .addComponent(spPepFDR, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addComponent(jLabel2))
                    .addComponent(ckDirectionalPeptidePair, javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(jLabel28, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.PREFERRED_SIZE, 19, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel8, javax.swing.GroupLayout.Alignment.TRAILING))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        .addComponent(spProteinFDR, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addComponent(jLabel6))
                    .addComponent(spMinProteinPepCount, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(spMaxProteinAmbiguity, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(spMaxLinkAmbiguity, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(spMinLinkPepCount, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(ckDirectionalLink, javax.swing.GroupLayout.Alignment.TRAILING)
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        .addComponent(spLinkFDR, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addComponent(jLabel3)))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(spMinPPIPepCount, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(ckDirectionalPPI, javax.swing.GroupLayout.Alignment.TRAILING)
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        .addComponent(spPPIFdr, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addComponent(jLabel4)))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(spReportFactor, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel7)
                    .addComponent(ckUniquePSM))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 6, Short.MAX_VALUE)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(ckMaximize)
                    .addComponent(cbBoostWhat, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel23)
                    .addComponent(spMaximizeSteps, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(ckBoostBetween))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(btnStopBoost)
                    .addComponent(btnCalc))
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

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton btnCalc;
    private javax.swing.JButton btnStopBoost;
    private org.rappsilber.fdr.gui.components.FDRLevelComboBox cbBoostWhat;
    private javax.swing.JCheckBox ckBoostBetween;
    private javax.swing.JCheckBox ckDirectionalLink;
    private javax.swing.JCheckBox ckDirectionalPPI;
    private javax.swing.JCheckBox ckDirectionalPSM;
    private javax.swing.JCheckBox ckDirectionalPeptidePair;
    private javax.swing.JCheckBox ckMaximize;
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
    private javax.swing.JSpinner spLinkFDR;
    private org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner spMaxLinkAmbiguity;
    private org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner spMaxProteinAmbiguity;
    private javax.swing.JSpinner spMaximizeSteps;
    private org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner spMinLinkPepCount;
    private org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner spMinPPIPepCount;
    private org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner spMinPeptideLength;
    private org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner spMinProteinPepCount;
    private javax.swing.JSpinner spPPIFdr;
    private javax.swing.JSpinner spPepFDR;
    private javax.swing.JSpinner spProteinFDR;
    private javax.swing.JSpinner spPsmFDR;
    private javax.swing.JSpinner spReportFactor;
    // End of variables declaration//GEN-END:variables




}
