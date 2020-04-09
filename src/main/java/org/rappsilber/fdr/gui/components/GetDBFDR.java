/*
 * Copyright 2020 Lutz Fischer <lfischer@staffmail.ed.ac.uk>.
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

import java.awt.EventQueue;
import java.sql.Connection;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JFrame;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JOptionPane;
import org.rappsilber.fdr.DBinFDR;
import org.rappsilber.fdr.OfflineFDR;
import org.rappsilber.fdr.entities.PSM;
import org.rappsilber.fdr.gui.FDRGUI;
import org.rappsilber.utils.RArrayUtils;
import rappsilber.config.RunConfig;
import rappsilber.ms.sequence.AminoModification;

/**
 *
 * @author Lutz Fischer <lfischer@staffmail.ed.ac.uk>
 */
public class GetDBFDR extends javax.swing.JPanel {
    OfflineFDR m_fdr;
    
    private FDRGUI fdrgui;
    
    
    

    /**
     * Creates new form GetDBFDR2
     */
    public GetDBFDR() {
        initComponents();
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
                JFrame getDBFDR_test = new JFrame("DBFDR TEST");
                getDBFDR_test.getContentPane().add(new GetDBFDR());
                getDBFDR_test.pack();
                getDBFDR_test.setVisible(true);
            }
        });
    }    


    protected void readNewFromDB() {
        try {

            setStatus("Start");
            Logger.getLogger(this.getClass().getName()).log(Level.INFO, "Start");
//            String searchId;
//            int[] selectedlines = lstSearches.getSelectedIndices();
            final int[] searchIds = getSearch.getSelectedSearchIds();
            if (searchIds.length == 0) {
                setStatus("Nothing Selected");
                return;
            }

            String titel = RArrayUtils.toString(searchIds,",") + " - " + getSearch.getSelectedNames()[0];
            if (searchIds.length > 1) {
                titel += " ...";
            }

            this.setTitle(titel);

            

            String txtconnection = null;
//            txtconnection = cmbConnection.getModel().getSelectedItem().toString();
//
//            Object o = Class.forName("org.postgresql.Driver");
            // Establish network connection to database
            Connection connection = getSearch.getConnection();
            final DBinFDR ofdr = new DBinFDR();
            if (!txtSubScorePattern.getText().trim().isEmpty()) {
                ofdr.setSubScoresToForward(txtSubScorePattern.getText());
            }
            ofdr.setDatabaseProvider(getSearch);
            m_fdr = ofdr;

            readData(ofdr, searchIds, (OfflineFDR.Normalisation) cbNormalize.getSelectedItem());
        } catch (Exception ex) {
            Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, null, ex);
            setStatus("error:" + ex);
        }
    }

    protected void addNormalizedFromDB() {
        try {
            final OfflineFDR.Normalisation normalize = (OfflineFDR.Normalisation) cbNormalize.getSelectedItem();
            setStatus("Start");
            Logger.getLogger(this.getClass().getName()).log(Level.INFO, "Start");
//            String searchId;
//            int[] selectedlines = lstSearches.getSelectedIndices();
            final int[] searchIds = getSearch.getSelectedSearchIds();
            if (searchIds.length == 0) {
                setStatus("Nothing Selected");
                return;
            }

            String titel = RArrayUtils.toString(searchIds,",") + " - " + getSearch.getSelectedNames()[0];
            if (searchIds.length > 1) {
                titel += " ...";
            }

            this.setTitle(titel);

            String txtconnection = null;
//            txtconnection = cmbConnection.getModel().getSelectedItem().toString();
//
//            Object o = Class.forName("org.postgresql.Driver");
            // Establish network connection to database
            Connection connection = getSearch.getConnection();
            final DBinFDR ofdr = new DBinFDR();
            ofdr.setDatabaseProvider(getSearch);

            readData(ofdr, searchIds, m_fdr,normalize);
            
        } catch (Exception ex) {
            Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, null, ex);
            setStatus("error:" + ex);
        }
    }
    

    protected void readData(final DBinFDR ofdr, final int[] searchIds, OfflineFDR.Normalisation normalize) {
        this.readData(ofdr, searchIds, null, normalize);
    }
    
    protected void readData(final DBinFDR ofdr, final int[] searchIds, final OfflineFDR addTo, final OfflineFDR.Normalisation normalize) {
        fdrgui.setEnableCalc(false);
        fdrgui.setEnableRead(false);
        if (additionalFDRGroups.getFlagAutoValidated())
            ofdr.setFlagAutoValidated(true);
        if (additionalFDRGroups.getFlagModified())
            ofdr.setMarkModifications(true);
        if (additionalFDRGroups.getFlagSearchID())
            ofdr.setMarkSearchID(true);
        if (additionalFDRGroups.getFlagRun())
            ofdr.setMarkRun(true);
        if (addTo == null) {
            PSM.resetAdditionalColumnNames();
        }
        
        Runnable runnable = new Runnable() {
            public void run() {
                try {
                    setStatus("Read from db");
                    Logger.getLogger(this.getClass().getName()).log(Level.INFO, "read from db");
                    ofdr.readDB(searchIds, txtReadFilter.getText(), ckToprankingOnly.isSelected());
                    
                    EventQueue.invokeLater(new Runnable() {
                        
                        public void run() {
                            fdrgui.ckPrePostAA.setEnabled(true);
                        }
                    });
                    if (fdrgui.getFdr() != null)
                        fdrgui.getFdr().cleanup();
                    if (m_fdr != null) {
                        fdrgui.setFdr(m_fdr);
                    }
                    fdrgui.setEnableRead(true);
                    fdrgui.setEnableCalc(true);
                    if (addTo != null) {
                        if (normalize != OfflineFDR.Normalisation.None) {
                            ofdr.coNormalizePSMs(addTo, normalize);
                        }
                        addTo.getAllPSMs().addAll(ofdr.getAllPSMs());
                        if (addTo instanceof DBinFDR) {
                            RunConfig prev_conf = ((DBinFDR)addTo).getConfig();
                            for (rappsilber.ms.sequence.AminoAcid aa:  ofdr.getConfig().getAllAminoAcids()) {
                                boolean missing = true;
                                for (rappsilber.ms.sequence.AminoAcid aaold:  prev_conf.getAllAminoAcids()) {
                                    if (aaold.SequenceID.equals(aa.SequenceID)) {
                                        missing = false;
                                        break;
                                    }
                                }
                                if (missing) {
                                    if (aa instanceof rappsilber.ms.sequence.AminoModification) {
                                        prev_conf.addKnownModification((AminoModification) aa);
                                    }
                                }
                            }
                        }
                        
                    }  else {
                        if (normalize != OfflineFDR.Normalisation.None) {
                            ofdr.normalizePSMs(normalize);
                        }
                    }
                    
                    
                    setStatus("finished reading");
//                    setEnableAdd(true);
                } catch (Exception ex) {
                    fdrgui.setEnableRead(true);
                    Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, null, ex);
                    setStatus("error:" + ex);
                    fdrgui.setEnableRead(true);
                }
                //  setEnableWrite(true);
            }
        };
        Thread t = new Thread(runnable);
        t.setName("Reading From DB");
        t.start();
    }

    
    private void setStatus(String message) {
        fdrgui.setStatus(message);
    }

    private void setTitle(String titel) {
        fdrgui.setTitle(titel);
    }

    public void setEnableRead(final boolean enable) {
        Runnable setModel = new Runnable() {
            public void run() {
                btnRead.setEnabled(enable);
            }
        };
        javax.swing.SwingUtilities.invokeLater(setModel);
    }    

    public void setEnableAdd(final boolean enable) {
        Runnable setModel = new Runnable() {
            public void run() {
                btnAdd.setEnabled(enable);
            }
        };
        javax.swing.SwingUtilities.invokeLater(setModel);
    }    
    
    /**
     * @return the fdrgui
     */
    public FDRGUI getFdrgui() {
        return fdrgui;
    }

    /**
     * @param fdrgui the fdrgui to set
     */
    public void setFdrgui(FDRGUI fdrgui) {
        this.fdrgui = fdrgui;
    }
    

    
    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        additionalFDRGroups = new org.rappsilber.fdr.gui.components.AdditionalFDRGroups();
        getSearch = new rappsilber.gui.components.db.GetSearch();
        txtReadFilter = new javax.swing.JTextField();
        jButton1 = new javax.swing.JButton();
        ckToprankingOnly = new javax.swing.JCheckBox();
        btnReadFilter = new javax.swing.JButton();
        btnSelectIDS = new javax.swing.JButton();
        lblSubScorePattern = new javax.swing.JLabel();
        txtSubScorePattern = new javax.swing.JTextField();
        btnRead = new javax.swing.JButton();
        btnAdd = new javax.swing.JButton();
        cbNormalize = new javax.swing.JComboBox<>();
        jLabel2 = new javax.swing.JLabel();

        txtReadFilter.setText("(site1 > 0 OR pepSeq2 isnull)");
        txtReadFilter.setToolTipText("Texttual representation of the filter");

        jButton1.setText("Additional Groups");
        jButton1.setToolTipText("define some additional grouping for PSMs");
        jButton1.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jButton1ActionPerformed(evt);
            }
        });

        ckToprankingOnly.setSelected(true);
        ckToprankingOnly.setText("Top Only");

        btnReadFilter.setText("Filter");
        btnReadFilter.setToolTipText("GUI for selecting filter");
        btnReadFilter.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnReadFilterActionPerformed(evt);
            }
        });

        btnSelectIDS.setText("Select By IDs");
        btnSelectIDS.setToolTipText("Seelct searches by ID");
        btnSelectIDS.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnSelectIDSActionPerformed(evt);
            }
        });

        lblSubScorePattern.setText("Forward");

        txtSubScorePattern.setToolTipText("regular expression for what columns should be forwarded to the PSM-files");

        btnRead.setText("read");
        btnRead.setToolTipText("Read in data from DB");
        btnRead.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnReadActionPerformed(evt);
            }
        });

        btnAdd.setText("add");
        btnAdd.setToolTipText("add the selected searches to an already read in search");
        btnAdd.setEnabled(false);
        btnAdd.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnAddActionPerformed(evt);
            }
        });

        cbNormalize.setModel(new DefaultComboBoxModel<OfflineFDR.Normalisation>(OfflineFDR.Normalisation.values()));
        cbNormalize.setToolTipText("How to normalize the input data");

        jLabel2.setText("Normalize");

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(txtReadFilter)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(ckToprankingOnly)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jButton1))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(btnReadFilter)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(btnSelectIDS)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(lblSubScorePattern)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(txtSubScorePattern, javax.swing.GroupLayout.DEFAULT_SIZE, 104, Short.MAX_VALUE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jLabel2)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(cbNormalize, javax.swing.GroupLayout.PREFERRED_SIZE, 121, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(btnAdd)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(btnRead)))
                .addContainerGap())
            .addComponent(getSearch, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addComponent(getSearch, javax.swing.GroupLayout.PREFERRED_SIZE, 266, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(txtReadFilter, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jButton1)
                    .addComponent(ckToprankingOnly))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(btnReadFilter)
                    .addComponent(btnSelectIDS)
                    .addComponent(lblSubScorePattern)
                    .addComponent(txtSubScorePattern, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(btnRead)
                    .addComponent(btnAdd)
                    .addComponent(cbNormalize, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel2))
                .addContainerGap(17, Short.MAX_VALUE))
        );
    }// </editor-fold>//GEN-END:initComponents

    private void btnReadActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnReadActionPerformed
        readNewFromDB();
    }//GEN-LAST:event_btnReadActionPerformed

    private void btnAddActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnAddActionPerformed
        m_fdr = fdrgui.getFdr();
        addNormalizedFromDB();
        
    }//GEN-LAST:event_btnAddActionPerformed

    private void btnReadFilterActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnReadFilterActionPerformed
        new Thread() {
            public void run () {
                txtReadFilter.setText(DBFIlters.showAndGetFilter(txtReadFilter.getText()));
            }
        }.start();            

    }//GEN-LAST:event_btnReadFilterActionPerformed

    private void btnSelectIDSActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnSelectIDSActionPerformed
        String idsString = JOptionPane.showInputDialog(this, "Database IDs to select (comma-separated list)");
        String[] ids = idsString.trim().split("\\s*,\\s*");
        int[] iids = new int[ids.length];
        for (int i = 0; i<ids.length;i++) {
            iids[i]=Integer.parseInt(ids[i]);
        }
        getSearch.setSelectedSearchIds(iids);
    }//GEN-LAST:event_btnSelectIDSActionPerformed

    private void jButton1ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jButton1ActionPerformed
        additionalFDRGroups.showWindow();
    }//GEN-LAST:event_jButton1ActionPerformed


    // Variables declaration - do not modify//GEN-BEGIN:variables
    private org.rappsilber.fdr.gui.components.AdditionalFDRGroups additionalFDRGroups;
    private javax.swing.JButton btnAdd;
    private javax.swing.JButton btnRead;
    private javax.swing.JButton btnReadFilter;
    private javax.swing.JButton btnSelectIDS;
    private javax.swing.JComboBox<OfflineFDR.Normalisation> cbNormalize;
    private javax.swing.JCheckBox ckToprankingOnly;
    public rappsilber.gui.components.db.GetSearch getSearch;
    private javax.swing.JButton jButton1;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel lblSubScorePattern;
    private javax.swing.JTextField txtReadFilter;
    private javax.swing.JTextField txtSubScorePattern;
    // End of variables declaration//GEN-END:variables
}
