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
package org.rappsilber.fdr.gui.components.settings;

import org.rappsilber.fdr.FDRSettings;
import java.awt.Component;
import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JSpinner;
import javax.swing.SwingUtilities;
import org.rappsilber.fdr.FDRSettingsImpl;
import org.rappsilber.fdr.gui.components.BoostIncludes;
import org.rappsilber.fdr.gui.components.FDRSpinnerModel;

/**
 *
 * @author lfischer
 */
public abstract class FDRSettingsPanel extends javax.swing.JPanel implements FDRSettings {



    protected ArrayList<ActionListener> m_StopMaximizingListener = new ArrayList<ActionListener>();

    protected ArrayList<java.awt.event.ActionListener> m_calc_listener = new ArrayList<ActionListener>();

//    protected OfflineFDR.FDRLevel m_optimizeWhat; 
//    protected boolean boostSubScores = false;
    protected boolean boostPSM = true;
    protected boolean boostPeptidePairs = true;
    protected boolean boostProteins = true;
    protected boolean boostLinks = true;
    protected boolean groupByPSMCount = false;
    private boolean filterConsectutivePeptides = false;
//    private double subScoreCutOff = 1;
    private boolean boostPepCoverage = true;
    private boolean boostPepDeltaScore = true;
    private boolean boostMinFragments = false;
    private boolean groubByCrosslinkerStubs;
    private boolean twoStepBoost;
    private boolean boostPeptideDoublets;
    private boolean boostPeptideStubs;
    private double minMinPeptideDoubpletFilter;
    private double minPeptideStubFilter;
    private boolean filterBySelfAndMono = false;
    private boolean boostMinScore;
    private Double minScore = 0d;
    private Integer scoreTopNAggregate;

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

//    public boolean boostSubScores() {
//        return boostSubScores;
//    }
//
//    public void boostSubScores(boolean boost) {
//        boostSubScores = boost;
//    }

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
            BoostIncludes bi = new BoostIncludes(this);
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
            BoostIncludes bi = new BoostIncludes(this);
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
    /**
     * @return the groupByPPI
     */
    public boolean isGroupByPSMCount() {
        return groupByPSMCount;
    }

    /**
     * @param groupByPPI the groupByPPI to set
     */
    public void setGroupByPSMCount(boolean groupBySMCount) {
        this.groupByPSMCount = groupBySMCount;
    }    
    
    public void setAll(FDRSettings settings) {
        FDRSettingsImpl.transferSettings(settings, this);
    }


    protected static void setFDRLater(final JSpinner sp, final double value) {
        SwingUtilities.invokeLater(new Runnable() {

            public void run() {
                ((FDRSpinnerModel)sp.getModel()).setFDR(value);
            }
        });
            
    }

    protected void setValueLater(final JSpinner sp, final Object value) {
        SwingUtilities.invokeLater(new Runnable() {

            public void run() {
                sp.setValue(value);
            }
        });
            
    }
    
    protected void setValueLater(final JCheckBox ck, final boolean value) {
        SwingUtilities.invokeLater(new Runnable() {

            public void run() {
                ck.setSelected(value);
            }
        });
        
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

    @Override
    public boolean filterConsecutivePeptides() {
        return filterConsectutivePeptides;
    }

    @Override
    public void setFilterConsecutivePeptides(boolean filterConsecutive) {
        this.filterConsectutivePeptides=filterConsecutive;
    }
    
  
    @Override
    public boolean boostPepCoverage(){
        return this.boostPepCoverage;
    }
    
    @Override
    public void boostPepCoverage(boolean boost){
        this.boostPepCoverage = boost;
    }

    @Override
    public boolean boostDeltaScore(){
        return this.boostPepDeltaScore;
    }

    @Override
    public void boostDeltaScore(boolean boost){
        this.boostPepDeltaScore = boost;
    }
    
    
    @Override
    public boolean boostMinFragments(){
        return this.boostMinFragments;
    }

    @Override
    public void boostMinFragments(boolean boost){
        this.boostMinFragments = boost;
    }
    
//    public double getSubScoreCutOff() {
//        return this.subScoreCutOff;
//    }
//    public void setSubScoreCutOff(double localfdr) {
//        this.subScoreCutOff = localfdr;
//    }    
    
    @Override
    public void setGroupByCrosslinkerStubs(boolean group) {
        this.groubByCrosslinkerStubs = group;
    }

    @Override
    public boolean getGroupByCrosslinkerStubs() {
        return this.groubByCrosslinkerStubs;
    }

    @Override
    public boolean twoStepOptimization() {
        return this.twoStepBoost;
    }

    @Override
    public void twoStepOptimization(boolean stepped) {
        this.twoStepBoost = stepped;
    }

    @Override
    public boolean boostPeptideStubs(){
        return this.boostPeptideStubs;
    }

    @Override
    public void boostPeptideStubs(boolean boost){
        this.boostPeptideStubs = boost;
    }
    
    @Override
    public boolean boostPeptideDoublets(){
        return this.boostPeptideDoublets;
    }

    @Override
    public void boostPeptideDoublets(boolean boost){
        this.boostPeptideDoublets = boost;
    }


    @Override
    public double getMinPeptideStubFilter() {
        return this.minPeptideStubFilter;
    }
    
    @Override
    public void setMinPeptideStubFilter(double d) {
        this.minPeptideStubFilter = d;
    }
    
    @Override
    public double getMinPeptideDoubletFilter() {
        return this.minMinPeptideDoubpletFilter;
    }
    
    @Override
    public void setMinPeptideDoubletFilter(double d) {
        this.minMinPeptideDoubpletFilter = d;
    }

    @Override
    public boolean filterBySelfAndMono() {
        return filterBySelfAndMono;
    }

    @Override
    public void setfilterBySelfAndMono(boolean filter) {
        filterBySelfAndMono = filter;
    }

    
    
    @Override
    public boolean boostMinScore(){
        return this.boostMinScore;
    }

    @Override
    public void boostMinScore(boolean boost){
        this.boostMinScore = boost;
    }

    
    @Override
    public Double minScore(){
        return this.minScore;
    }

    @Override
    public void minScore(Double minScore){
        this.minScore = minScore;
    }

    @Override
    public Integer getScoreTopNAggregate() {
        return this.scoreTopNAggregate;
    }

    @Override
    public void setScoreTopNAggregate(Integer n) {
        this.scoreTopNAggregate = n;
    }
    
}
