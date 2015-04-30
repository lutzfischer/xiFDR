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
public class FDRSettingsSimple extends FDRSettingsPanel  {

//    private ArrayList<java.awt.event.ActionListener> m_calc_listener = new ArrayList<ActionListener>();
//    
//    private OfflineFDR.FDRLevel m_optimizeWhat; 
    
    private int m_minPepLength = 0;
    private int m_boostingSteps = 4;
    private double m_psmfdr = 100;
    private double m_pepfdr = 100;
    private double m_protfdr = 100;
    private double m_linkfdr = 100;
    private double m_ppifdr = 100;
    
    private boolean m_filterToUniquePSM = true;

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
    
    /**
     * Creates new form FDRSettingsComplete
     */
    public FDRSettingsSimple() {
        initComponents();
        
        ChangeListener max100Listener = new ChangeListener() {
            public void stateChanged(ChangeEvent e) {
                SpinnerModel sp = (SpinnerModel) e.getSource();
                if (((Double)sp.getValue()) >100)
                    sp.setValue(100d);
            }
        };        
        spFDR.getModel().addChangeListener(max100Listener);
                
    }
    
    private void doCalc() {
 //       if (ckMaximize.isSelected()) {
            switch ((OfflineFDR.FDRLevel) cbFDRLevel.getSelectedItem()) {
                case PSM: 
                    m_pepfdr = m_linkfdr = m_protfdr = m_ppifdr = 1;
                    break;
                case PEPTIDE_PAIR: 
                    m_psmfdr = m_linkfdr = m_protfdr = m_ppifdr = 1;
                    break;
                case PROTEINGROUP: 
                    m_pepfdr = m_linkfdr = m_psmfdr = m_ppifdr = 1;
                    break;
                case PROTEINGROUPLINK: 
                    m_pepfdr = m_psmfdr = m_protfdr = m_ppifdr = 1;
                    break;
                case PROTEINGROUPPAIR: 
                    m_pepfdr = m_linkfdr = m_protfdr = m_psmfdr = 1;
                    break;
            }
            btnStopBoost.setEnabled(ckMaximize.isSelected());            
            raiseStartCalc(ckMaximize.isSelected());
///        } else {
//            for (ActionListener al : m_calc_listener) {
//                al.actionPerformed(new ActionEvent(this, -1, "calc"));
//            }
//        }

    }
    
 
    
    @Override
    public double getPSMFDR() {
        if (cbFDRLevel.getSelectedItem() == OfflineFDR.FDRLevel.PSM)
            return ((Double) spFDR.getValue())/100.0;
        else 
            return m_psmfdr;
    }

    
    @Override
    public double getPeptidePairFDR() {
        if (cbFDRLevel.getSelectedItem() == OfflineFDR.FDRLevel.PEPTIDE_PAIR)
            return ((Double) spFDR.getValue())/100.0;
        else 
            return m_pepfdr;
    }

    @Override
    public double getProteinGroupFDR() {
        if (cbFDRLevel.getSelectedItem() == OfflineFDR.FDRLevel.PROTEINGROUP)
            return ((Double) spFDR.getValue()) /100;
        else 
            return m_protfdr;
    }

    @Override
    public double getProteinGroupLinkFDR() { 
        if (cbFDRLevel.getSelectedItem() == OfflineFDR.FDRLevel.PROTEINGROUPLINK)
            return ((Double) spFDR.getValue())/100.0;
        else 
            return m_ppifdr;
    }

    @Override
    public double getProteinGroupPairFDR() {
        if (cbFDRLevel.getSelectedItem() == OfflineFDR.FDRLevel.PROTEINGROUPPAIR)
            return ((Double) spFDR.getValue())/100.0;
        else 
            return 1;
    }

    @Override
    public void setPSMFDR( Double fdr) {
        if (cbFDRLevel.getSelectedItem() == OfflineFDR.FDRLevel.PSM)
            setValueLater(spFDR, fdr*100.0);
        m_psmfdr = fdr;
    }

    
    @Override
    public void setPeptidePairFDR( Double fdr) {
        if (cbFDRLevel.getSelectedItem() == OfflineFDR.FDRLevel.PEPTIDE_PAIR)
            setValueLater(spFDR, fdr*100.0);
        m_pepfdr = fdr;
    }

    @Override
    public void setProteinGroupFDR( Double fdr) {
        if (cbFDRLevel.getSelectedItem() == OfflineFDR.FDRLevel.PROTEINGROUP)
            setValueLater(spFDR, fdr*100.0);
        m_protfdr = fdr;
    }

    @Override
    public void setProteinGroupLinkFDR( Double fdr) {
        if (cbFDRLevel.getSelectedItem() == OfflineFDR.FDRLevel.PROTEINGROUPLINK)
            setValueLater(spFDR, fdr*100.0);
        m_linkfdr = fdr;
    }

    @Override
    public void setProteinGroupPairFDR( Double fdr) {
        if (cbFDRLevel.getSelectedItem() == OfflineFDR.FDRLevel.PROTEINGROUPPAIR)
            setValueLater(spFDR, fdr*100.0);
        m_ppifdr = fdr;
    }

    @Override
    public int getMinProteinPepCount() {
        return 1;
    }

    @Override
    public int getMinLinkPepCount() {
        return 1;
    }

    @Override
    public int getMinPPIPepCount() {
        return 1;
    }

    @Override
    public void setMinProteinPepCount(Integer minPep) {
    }

    @Override
    public void setMinLinkPepCount(Integer minPep) {
    }

    @Override
    public void setMinPPIPepCount(Integer minPep) {
    }

    
    

    @Override
    public int getMinPeptideLength() {
        return m_minPepLength;
    }

    @Override
    public void setMinPeptideLength(Integer minLength) {
        m_minPepLength = minLength;
    }

    
    @Override
    public int getMaxLinkAmbiguity() {
        return 0;
    }

    @Override
    public void setMaxLinkAmbiguity(Integer maxAmbiguity) {
        
    }

    @Override
    public int getMaxProteinAmbiguity() {
        return 0;
    }

    @Override
    public void setMaxProteinAmbiguity(Integer maxAmbiguity) {
    }
    
    @Override
    public boolean isPSMDirectional() {
        return ckDirectional.isSelected();
    }
    
    @Override
    public boolean isPeptidePairDirectional() {
        return ckDirectional.isSelected();
    }

    @Override
    public boolean isLinkDirectional() {
        return ckDirectional.isSelected();
    }

    @Override
    public boolean isPPIDirectional() {
        return ckDirectional.isSelected();
    }

    @Override
    public void setPSMDirectional(boolean directional) {
        setValueLater(ckDirectional,directional);
    }
    
    @Override
    public void setPeptidePairDirectional(boolean directional) {
        setValueLater(ckDirectional,directional);
    }

    @Override
    public void setLinkDirectional(boolean directional) {
        setValueLater(ckDirectional,directional);
    }

    @Override
    public void setPPIDirectional(boolean directional) {
        setValueLater(ckDirectional,directional);
    }
    
    
    @Override
    public OfflineFDR.FDRLevel doOptimize() {
        if (!ckMaximize.isSelected())
            return null;
        else return (OfflineFDR.FDRLevel) cbFDRLevel.getSelectedItem();
    }
    
    
    public int getBoostingSteps() {
        return m_boostingSteps;
    }

    public void setBoostingSteps(int steps) {
        m_boostingSteps = steps;
    }
    
    public double getReportFactor() {
        return (Double) spReportFactor.getValue();
    }    
    
    public void setReportFactor(double factor) {
        spReportFactor.setValue(factor);
    }       

    
    public boolean filterToUniquePSM() {
        return m_filterToUniquePSM;
    }

    public void setFilterToUniquePSM(boolean filterToUniquePSM) {
        m_filterToUniquePSM = filterToUniquePSM;
    }
    
    public void setAll(FDRSettings settings) {
        cbFDRLevel.getModel().setSelectedItem(settings.doOptimize());
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
        this.setFilterToUniquePSM(settings.filterToUniquePSM());
        this.setBoostBetween(settings.getBoostBetween());
    }
    
    
    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jLabel7 = new javax.swing.JLabel();
        spReportFactor = new javax.swing.JSpinner();
        jLabel1 = new javax.swing.JLabel();
        ckDirectional = new javax.swing.JCheckBox();
        btnCalc = new javax.swing.JButton();
        ckMaximize = new javax.swing.JCheckBox();
        lblDirectional = new javax.swing.JLabel();
        lblBoost = new javax.swing.JLabel();
        spFDR = new javax.swing.JSpinner();
        cbFDRLevel = new org.rappsilber.fdr.gui.components.FDRLevelComboBox();
        btnStopBoost = new javax.swing.JButton();
        ckBoostBetween = new javax.swing.JCheckBox();

        jLabel7.setText("Report Factor");
        jLabel7.setToolTipText("maximum factor the next step in fdr is permited to exced the target fdr");

        spReportFactor.setModel(new javax.swing.SpinnerNumberModel(Double.valueOf(1.1d), Double.valueOf(1.0d), null, Double.valueOf(0.1d)));
        spReportFactor.setToolTipText("maximum factor the next step in fdr is permited to exced the target fdr");

        jLabel1.setText("FDR");

        ckDirectional.setToolTipText("Is the Cross-linking considered directional (A links to B is different to B links to A)");

        btnCalc.setText("Calculate");
        btnCalc.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnCalcActionPerformed(evt);
            }
        });

        ckMaximize.setSelected(true);
        ckMaximize.setToolTipText("Use prefiltering on lower level to boost the results on the chossen level of information");
        ckMaximize.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                ckMaximizeActionPerformed(evt);
            }
        });

        lblDirectional.setText("Is directional");
        lblDirectional.setToolTipText("Is the Cross-linking considered directional (A links to B is different to B links to A)");
        lblDirectional.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                lblDirectionalMouseClicked(evt);
            }
        });

        lblBoost.setText("Boost result");
        lblBoost.setToolTipText("Use prefiltering on lower level to boost the results on the chossen level of information");
        lblBoost.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                lblBoostMouseClicked(evt);
            }
        });

        spFDR.setModel(new javax.swing.SpinnerNumberModel(Double.valueOf(100.0d), Double.valueOf(0.0d), null, Double.valueOf(1.0d)));
        spFDR.setToolTipText("What is the maximal acceptable FDR");

        cbFDRLevel.setSelectedIndex(3);
        cbFDRLevel.setToolTipText("To what to apply the FDR");

        btnStopBoost.setText("stop boost");
        btnStopBoost.setEnabled(false);
        btnStopBoost.setMargin(new java.awt.Insets(2, 7, 2, 7));
        btnStopBoost.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnStopBoostActionPerformed(evt);
            }
        });

        ckBoostBetween.setText("Between");

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(btnStopBoost, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(jLabel1)
                    .addComponent(lblDirectional)
                    .addComponent(lblBoost)
                    .addComponent(jLabel7, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addComponent(spReportFactor, javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(spFDR)
                            .addGroup(layout.createSequentialGroup()
                                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(ckDirectional)
                                    .addComponent(ckMaximize))
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 58, Short.MAX_VALUE)
                                .addComponent(ckBoostBetween)))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(cbFDRLevel, javax.swing.GroupLayout.DEFAULT_SIZE, 140, Short.MAX_VALUE))
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(btnCalc)))
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel1)
                    .addComponent(spFDR, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(cbFDRLevel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(ckDirectional)
                            .addComponent(lblDirectional))
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(ckMaximize)
                            .addComponent(lblBoost)))
                    .addComponent(ckBoostBetween, javax.swing.GroupLayout.Alignment.TRAILING))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel7)
                    .addComponent(spReportFactor, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 32, Short.MAX_VALUE)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(btnCalc)
                    .addComponent(btnStopBoost))
                .addContainerGap())
        );
    }// </editor-fold>//GEN-END:initComponents

    private void btnCalcActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnCalcActionPerformed
        // TODO add your handling code here:
        doCalc();
        
    }//GEN-LAST:event_btnCalcActionPerformed

    private void lblDirectionalMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_lblDirectionalMouseClicked
        ckDirectional.setSelected(!ckDirectional.isSelected());
    }//GEN-LAST:event_lblDirectionalMouseClicked

    private void lblBoostMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_lblBoostMouseClicked
        ckMaximize.setSelected(!ckMaximize.isSelected());
        
    }//GEN-LAST:event_lblBoostMouseClicked

    private void btnStopBoostActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnStopBoostActionPerformed
        raiseStopMaximizing();
    }//GEN-LAST:event_btnStopBoostActionPerformed

    private void ckMaximizeActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_ckMaximizeActionPerformed
        ckBoostBetween.setEnabled(ckMaximize.isSelected());
    }//GEN-LAST:event_ckMaximizeActionPerformed

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton btnCalc;
    private javax.swing.JButton btnStopBoost;
    private org.rappsilber.fdr.gui.components.FDRLevelComboBox cbFDRLevel;
    private javax.swing.JCheckBox ckBoostBetween;
    private javax.swing.JCheckBox ckDirectional;
    private javax.swing.JCheckBox ckMaximize;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel7;
    private javax.swing.JLabel lblBoost;
    private javax.swing.JLabel lblDirectional;
    private javax.swing.JSpinner spFDR;
    private javax.swing.JSpinner spReportFactor;
    // End of variables declaration//GEN-END:variables

    @Override
    public void setEnabled(boolean e) {
        super.setEnabled(e);
        btnCalc.setEnabled(e);
        if (e) 
            btnStopBoost.setEnabled(false);
    }

//    @Override
//    public boolean isEnabeled(boolean e) {
//    }


}
