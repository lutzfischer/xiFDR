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

import org.rappsilber.fdr.FDRSettings;
import java.awt.Component;
import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import javax.swing.JDialog;
import javax.swing.JFrame;

/**
 *
 * @author lfischer
 */
public abstract class FDRSettingsPanel extends javax.swing.JPanel implements FDRSettings {

    protected ArrayList<ActionListener> m_StopMaximizingListener = new ArrayList<ActionListener>();

    protected ArrayList<java.awt.event.ActionListener> m_calc_listener = new ArrayList<ActionListener>();

//    protected OfflineFDR.FDRLevel m_optimizeWhat; 
    protected boolean boostPSM = true;
    protected boolean boostPeptidePairs = true;
    protected boolean boostProteins = true;
    protected boolean boostLinks = true;

    /**
     * Creates new form FDRSettingsPanel
     */
    public FDRSettingsPanel() {
        initComponents();
    }

    public void raiseStopMaximizing() {
        ActionEvent e = new ActionEvent(this, 0, "STOP");
        for (ActionListener al : m_StopMaximizingListener) {
            al.actionPerformed(e);
        }
    }

    public void addStopMaximizingListener(java.awt.event.ActionListener al) {
        m_StopMaximizingListener.add(al);
    }

    public void removeStopMaximizingListener(java.awt.event.ActionListener al) {
        m_StopMaximizingListener.remove(al);
    }

    @Override
    public void addCalcListener(java.awt.event.ActionListener al) {
        this.m_calc_listener.add(al);
    }

    protected void raiseStartCalc(boolean boost) {
        ActionEvent e;
        if (boost) {
            e = new ActionEvent(this, 0, "boost");
        } else {
            e = new ActionEvent(this, -1, "calc");
        }

        for (ActionListener al : m_calc_listener) {
            al.actionPerformed(e);
        }
    }

    @Override
    public void removeCalcListener(java.awt.event.ActionListener al) {
        this.m_calc_listener.remove(al);
    }


    public boolean boostProteins() {
        return boostProteins;
    }

    public void boostProteins(boolean boost) {
        boostProteins = boost;

    }

    public boolean boostPSMs() {
        return boostPSM;
    }

    public void boostPSMs(boolean boost) {
        boostPSM = boost;

    }

    public boolean boostPeptidePairs() {
        return boostPeptidePairs;
    }

    public void boostPeptidePairs(boolean boost) {
        boostPeptidePairs = boost;
    }

    public boolean boostLinks() {
        return boostLinks;
    }

    public void boostLinks(boolean boost) {
        boostLinks = boost;
    }

    public void setBoostIgnores() {
        Component c = this.getParent();
        while (!(c==null || c instanceof Frame) ) {
            c=c.getParent();
        }
        if (c!= null) {
            final JDialog fi = new JDialog((Frame) c,"Boost Ignores",true);
            BoostIgnores bi = new BoostIgnores(this);
            fi.getContentPane().add(bi);
            fi.pack();
            bi.addActionListener(new ActionListener() {

                @Override
                public void actionPerformed(ActionEvent e) {
                    fi.setVisible(false);
                    fi.dispose();
                }
            });
            fi.setVisible(true);
        } else {
            final JFrame fi = new JFrame("Boost Ignores");
            BoostIgnores bi = new BoostIgnores(this);
            fi.getContentPane().add(bi);
            fi.pack();
            bi.addActionListener(new ActionListener() {

                @Override
                public void actionPerformed(ActionEvent e) {
                    fi.setVisible(false);
                    fi.dispose();
                }
            });
            fi.setVisible(true);
        }
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
        this.doOptimize(settings.doOptimize());
        this.setFilterToUniquePSM(settings.filterToUniquePSM());
        this.setBoostBetween(settings.getBoostBetween());
        this.boostLinks(settings.boostLinks());
        this.boostPSMs(settings.boostPSMs());
        this.boostPeptidePairs(settings.boostPeptidePairs());
        this.boostProteins(settings.boostProteins());
    }
    
//    public abstract void setEnabled(boolean e);
//
//    public abstract boolean isEnabeled(boolean e);
    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 400, Short.MAX_VALUE)
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGap(0, 300, Short.MAX_VALUE)
        );
    }// </editor-fold>//GEN-END:initComponents
    // Variables declaration - do not modify//GEN-BEGIN:variables
    // End of variables declaration//GEN-END:variables
}
