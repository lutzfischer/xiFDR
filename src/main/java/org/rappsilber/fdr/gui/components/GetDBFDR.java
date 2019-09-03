/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.rappsilber.fdr.gui.components;

import java.awt.EventQueue;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.sql.Connection;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JOptionPane;
import org.rappsilber.fdr.DBinFDR;
import org.rappsilber.fdr.gui.FDRGUI;
import org.rappsilber.utils.IntArrayList;
import org.rappsilber.utils.RArrayUtils;

/**
 *
 * @author lfischer
 */
public class GetDBFDR extends javax.swing.JPanel {
    javax.swing.JButton btnReadFilter = new javax.swing.JButton();
    javax.swing.JButton btnSelectIDS = new javax.swing.JButton();
    public javax.swing.JButton btnRead = new javax.swing.JButton();
    public javax.swing.JCheckBox ckAdd = new javax.swing.JCheckBox();
    public javax.swing.JCheckBox ckNormalize = new javax.swing.JCheckBox();
    javax.swing.JButton btnGetDBSize = new javax.swing.JButton();
    javax.swing.JTextField txtReadFilter = new javax.swing.JTextField();
    AdditionalFDRGroups additionalFDRGroups= new AdditionalFDRGroups();
    javax.swing.JButton btnGroups = new javax.swing.JButton("Additional Groups");
    javax.swing.JCheckBox ckToprankingOnly = new javax.swing.JCheckBox();
    public rappsilber.gui.components.db.GetSearch getSearch = new rappsilber.gui.components.db.GetSearch();
    DBinFDR m_fdr;
    
    private FDRGUI fdrgui;
    
    private double m_DBSizePSMTarget = 999999999;
    private double m_DBSizePSMDecoy = 999999999;
    private double m_DBSizeLinkTarget = 999999999;
    private double m_DBSizeLinkDecoy = 999999999;
    private double m_DBSizeProteinTarget = 999999999;
    private double m_DBSizeProteinDecoy = 999999999;
    
    /**
     * Creates new form GetDBFDR
     */
    public GetDBFDR() {
        initComponents();

        btnRead.setText("Read");
        btnRead.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                if (ckAdd.isSelected()) {
                    
                    addNormalizedFromDB();
                } else {
                    readNewFromDB();
                }
                
            }
        });

        ckAdd.setText("Add");
        ckNormalize.setText("Normalize");

        ckAdd.setEnabled(false);

        btnGetDBSize.setText("Estimate DB-Size");
        btnGetDBSize.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                readSizeFromDB();
            }
        });

        txtReadFilter.setText("");

        btnReadFilter.setText("Filter");
        btnReadFilter.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                new Thread() {
                    public void run () {
//                        txtReadFilter.setText(DBFIlters.showAndGetFilter());
                        txtReadFilter.setText(DBFIlters.showAndGetFilter(txtReadFilter.getText()));
                    }
                }.start();            
            }
        });

        btnSelectIDS.setText("Select By IDs");
        btnSelectIDS.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                String idsString = JOptionPane.showInputDialog(this, "Database IDs to select (comma-separated list)");
                String[] ids = idsString.trim().split("\\s*,\\s*");
                int[] iids = new int[ids.length];
                for (int i = 0; i<ids.length;i++) {
                    iids[i]=Integer.parseInt(ids[i]);
                }
                getSearch.setSelectedSearchIds(iids);
            }
        });

        ckToprankingOnly.setEnabled(true);        
        ckToprankingOnly.setSelected(true);        
        ckToprankingOnly.setText("Top Ranking only");
        btnGroups.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                additionalFDRGroups.showWindow();
            }
        });
        //ckToprankingOnly.setEnabled(false);        
        
        javax.swing.GroupLayout pInputDBLayout = new javax.swing.GroupLayout(this);
        this.setLayout(pInputDBLayout);
        pInputDBLayout.setHorizontalGroup(pInputDBLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(pInputDBLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(pInputDBLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(pInputDBLayout.createSequentialGroup()
                        .addComponent(btnReadFilter)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(btnSelectIDS)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(ckNormalize)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(ckAdd)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(btnRead)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(btnGetDBSize))
                    .addGroup(pInputDBLayout.createSequentialGroup()
                        .addComponent(txtReadFilter)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(ckToprankingOnly)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(btnGroups))
                    .addComponent(getSearch, javax.swing.GroupLayout.DEFAULT_SIZE, 690, Short.MAX_VALUE))
                .addContainerGap())
        );
        pInputDBLayout.setVerticalGroup(pInputDBLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(pInputDBLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(getSearch, javax.swing.GroupLayout.DEFAULT_SIZE, 335, Short.MAX_VALUE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(pInputDBLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(txtReadFilter, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(ckToprankingOnly)
                    .addComponent(btnGroups))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(pInputDBLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(btnRead)
                    .addComponent(ckAdd)
                    .addComponent(ckNormalize)
                    .addComponent(btnGetDBSize)
                    .addComponent(btnReadFilter)
                    .addComponent(btnSelectIDS))
                .addContainerGap())
        );
        
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
            ofdr.setDatabaseProvider(getSearch);
            m_fdr = ofdr;
            btnGetDBSize.setEnabled(true);

            readData(ofdr, searchIds, ckNormalize.isSelected());
        } catch (Exception ex) {
            Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, null, ex);
            setStatus("error:" + ex);
        }
    }

    protected void addNormalizedFromDB() {
        try {
            final boolean normalize = ckNormalize.isSelected();
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
            btnGetDBSize.setEnabled(true);

            readData(ofdr, searchIds, m_fdr, ckNormalize.isSelected());
            
        } catch (Exception ex) {
            Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, null, ex);
            setStatus("error:" + ex);
        }
    }
    

    protected void addFromDB() {
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

            IntArrayList prevIds = ((DBinFDR) m_fdr).getSearchIDs();
            String titel = RArrayUtils.toString(prevIds,",") +"," + RArrayUtils.toString(searchIds,",") + " - " + getSearch.getSelectedNames()[0];
            if (searchIds.length + prevIds.size() > 1) {
                titel += " ...";
            }

            this.setTitle(titel);

            

            String txtconnection = null;
//            txtconnection = cmbConnection.getModel().getSelectedItem().toString();
//
//            Object o = Class.forName("org.postgresql.Driver");
            // Establish network connection to database
            btnGetDBSize.setEnabled(false);

            readData(m_fdr, searchIds, ckNormalize.isSelected());
        } catch (Exception ex) {
            Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, null, ex);
            setStatus("error:" + ex);
        }
    }

    protected void readData(final DBinFDR ofdr, final int[] searchIds, boolean normalize) {
        this.readData(ofdr, searchIds, null, normalize);
    }
    
    protected void readData(final DBinFDR ofdr, final int[] searchIds, final DBinFDR addTo, final boolean normalize) {
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
                    fdrgui.setFdr(m_fdr);
                    fdrgui.setEnableRead(true);
                    fdrgui.setEnableCalc(true);
                    if (normalize) {
                        ofdr.normalizePSMsToFDR();
                    }
                    
                    if (addTo != null) {
                        if (normalize && !addTo.isNormalized()) {
                            addTo.normalizePSMsToFDR();
                        }
                        addTo.getAllPSMs().addAll(ofdr.getAllPSMs());
//                        if (normalize) {
//                            
//                            addTo.addNormalisedPsmList(ofdr.getAllPSMs(), ofdr.getPsmNormalizationOffset());
//                        } else {
//                            addTo.getAllPSMs().addAll(ofdr.getAllPSMs());
//                        }
//                        ofdr.normalizePSMs();
//                        if (!addTo.isNormalized())
//                            addTo.normalizePSMs();
//                        addTo.addNormalisedPsmList(ofdr.getAllPSMs(), ofdr.getPsmNormalizationOffset());
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
    
    
    protected void readSizeFromDB() {
        try {

            if (!(m_fdr instanceof DBinFDR)) {
                JOptionPane.showMessageDialog(this, "Data are not read from the db!");
                return;
            }
            Logger.getLogger(this.getClass().getName()).log(Level.INFO, "Start");

            btnGetDBSize.setEnabled(false);

            final DBinFDR ofdr = (DBinFDR) m_fdr;

            Runnable runnable = new Runnable() {
                public void run() {
                    try {
                        setStatus("Read db sizes from db");
                        Logger.getLogger(this.getClass().getName()).log(Level.INFO, "Read db sizes from db");

                        ofdr.getDBSizes();

                        m_DBSizePSMTarget =  ofdr.getDecoyDBSize();
                        m_DBSizePSMDecoy =  ofdr.getTargetDBSize();

                        m_DBSizeLinkTarget =  ofdr.getDecoyLinkDBSize();
                        m_DBSizeLinkDecoy =  ofdr.getTargetLinkDBSize();
                        
                        m_DBSizeProteinDecoy =  ofdr.getDecoyProtDBSize();
                        m_DBSizeProteinDecoy =  ofdr.getTargetProtDBSize();

                        setStatus("finished reading db sizes");
                    } catch (Exception ex) {
                        Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, null, ex);
                        setStatus("error:" + ex);
                    }

                    btnGetDBSize.setEnabled(m_fdr instanceof DBinFDR);

                    //  setEnableWrite(true);
                }
            };
            Thread t = new Thread(runnable);
            t.setName("Reading sizes From DB");
            t.start();
        } catch (Exception ex) {
            Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, null, ex);
            setStatus("error:" + ex);
        }
    }
    
    
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
                ckAdd.setEnabled(enable);
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
}
