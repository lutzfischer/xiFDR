/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.rappsilber.fdr.gui.components;

import java.sql.SQLException;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JOptionPane;
import org.rappsilber.fdr.DBinFDR;
import org.rappsilber.fdr.gui.FDRGUI;

/**
 *
 * @author lfischer
 */
public class WriteToDB extends javax.swing.JPanel {

    FDRGUI m_gui;
    /**
     * Creates new form WriteToDB
     */
    public WriteToDB() {
        initComponents();
    }

    
    private void writeDB() {
        final DBinFDR ofdr = (DBinFDR) m_gui.getFdr();
        
        m_gui.setEnableCalc(false);
        this.setEnableWrite(false);
        m_gui.setEnableRead(false);
        final boolean overwrite = ckResValidateOverwrite.isSelected();
        final String validate = ckResValidate.isSelected() ? txtResValidate.getText() : null;
        final boolean writefdr = ckResWriteFDR.isSelected();
        final boolean within = rbWithin.isSelected() || rbAll.isSelected();
        final boolean between = rbBetween.isSelected()  || rbAll.isSelected();
        

        Runnable runnable = new Runnable() {
            public void run() {
                try {
                    m_gui.setStatus("Writing to db: " + ofdr.summaryString(m_gui.getResult()));
                    ofdr.writeDB(validate, overwrite, writefdr,m_gui.getResult(),within, between);
                    m_gui.setStatus("finished writing: " + ofdr.summaryString(m_gui.getResult()));

                } catch (SQLException sex) {
                    Logger.getLogger(FDRGUI.class.getName()).log(Level.SEVERE, null, sex);
                    m_gui.setStatus(sex.toString());
                    SQLException ex = sex;
                    int chain = 1;
                    while ((ex = ex.getNextException()) != null && chain++ <10) {
                        Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, "next exception", sex);
                    }

                }
                m_gui.setEnableCalc(true);
                WriteToDB.this.setEnableWrite(true);
                m_gui.setEnableRead(true);
            }
        };

        new Thread(runnable).start();


    }
 
    
    private void clearValidation() {
        clearValidation(ckResValidate.isSelected() ? txtResValidate.getText() : "");
        
    }
    
    public void setEnableWrite(final boolean enable) {
        Runnable setModel = new Runnable() {
            public void run() {
                btnWriteDB.setEnabled(enable);
                ckResValidateActionPerformed(null);
            }
        };
        javax.swing.SwingUtilities.invokeLater(setModel);
    }

    /**
     *  CLEAR A SPECIFIC VALIDATION FROM THE DATABSE
     * @param validate what validation to clear
     */
    private void clearValidation(final String validate) {

        final DBinFDR ofdr = (DBinFDR)m_gui.getFdr();
        m_gui.setEnableCalc(false);
        this.setEnableWrite(false);
        m_gui.setEnableRead(false);

        Runnable runnable = new Runnable() {
            public void run() {
                try {
                    if (validate.isEmpty()) {
                        if (JOptionPane.showConfirmDialog(m_gui, "Erase all validations ?", "Clear Validation", JOptionPane.OK_CANCEL_OPTION) == JOptionPane.OK_OPTION) {
                            ofdr.clearDBValidation("");
                            m_gui.setStatus("cleared   ");
                        }
                    } else {
                        if (JOptionPane.showConfirmDialog(m_gui, "Erase all validations of type " + validate + "?", "Clear Validation", JOptionPane.OK_CANCEL_OPTION) == JOptionPane.OK_OPTION) {
                            ofdr.clearDBValidation(validate);
                            m_gui.setStatus("cleared   ");
                        }
                    }
                } catch (SQLException ex) {
                    Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, null, ex);
                    m_gui.setStatus(ex.toString());
                    while ((ex = ex.getNextException()) != null) {
                        Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, "Next exception ", ex);
                    }

                }
                m_gui.setEnableCalc(true);
                WriteToDB.this.setEnableWrite(true);
                m_gui.setEnableRead(true);
            }
        };

        new Thread(runnable).start();

    }
    /** 
     * clear the validation from the database  
     */
    private void clearDBFDR() {

        final DBinFDR ofdr = (DBinFDR) m_gui.getFdr();
        m_gui.setEnableCalc(false);
        this.setEnableWrite(false);
        m_gui.setEnableRead(false);
        final String validate = ckResValidate.isSelected() ? txtResValidate.getText() : null;


        Runnable runnable = new Runnable() {
            public void run() {
                try {
                    if (JOptionPane.showConfirmDialog(m_gui, "delete previously created fdr values in the DB?", "Clear FDR", JOptionPane.OK_CANCEL_OPTION) == JOptionPane.OK_OPTION) {
                        ofdr.clearDBFDR();
                        m_gui.setStatus("finished CLEANING ");
                    }
                } catch (SQLException ex) {
                    Logger.getLogger(FDRGUI.class.getName()).log(Level.SEVERE, null, ex);
                    m_gui.setStatus(ex.toString());

                }
                m_gui.setEnableCalc(true);
                WriteToDB.this.setEnableWrite(true);
                m_gui.setEnableRead(true);
            }
        };

        new Thread(runnable).start();

    }    
    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        buttonGroup1 = new javax.swing.ButtonGroup();
        btnResFDRClear = new javax.swing.JButton();
        btnResValidateClear = new javax.swing.JButton();
        btnWriteDB = new javax.swing.JButton();
        ckResWriteFDR = new javax.swing.JCheckBox();
        ckResValidateOverwrite = new javax.swing.JCheckBox();
        txtResValidate = new javax.swing.JTextField();
        ckResValidate = new javax.swing.JCheckBox();
        btnResValidateClearAll = new javax.swing.JButton();
        rbAll = new javax.swing.JRadioButton();
        rbWithin = new javax.swing.JRadioButton();
        rbBetween = new javax.swing.JRadioButton();

        btnResFDRClear.setText("clear");
        btnResFDRClear.setEnabled(false);
        btnResFDRClear.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnResFDRClearActionPerformed(evt);
            }
        });

        btnResValidateClear.setText("clear");
        btnResValidateClear.setEnabled(false);
        btnResValidateClear.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnResValidateClearActionPerformed(evt);
            }
        });

        btnWriteDB.setText("Write");
        btnWriteDB.setEnabled(false);
        btnWriteDB.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnWriteDBActionPerformed(evt);
            }
        });

        ckResWriteFDR.setText("Write fdrs");
        ckResWriteFDR.setEnabled(false);
        ckResWriteFDR.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                ckResWriteFDRActionPerformed(evt);
            }
        });

        ckResValidateOverwrite.setText("Overwrite");
        ckResValidateOverwrite.setEnabled(false);

        txtResValidate.setText("A");
        txtResValidate.setEnabled(false);

        ckResValidate.setText("Validate");
        ckResValidate.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                ckResValidateActionPerformed(evt);
            }
        });

        btnResValidateClearAll.setText("clear all");
        btnResValidateClearAll.setEnabled(false);
        btnResValidateClearAll.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnResValidateClearAllActionPerformed(evt);
            }
        });

        buttonGroup1.add(rbAll);
        rbAll.setSelected(true);
        rbAll.setText("All");

        buttonGroup1.add(rbWithin);
        rbWithin.setText("Within");

        buttonGroup1.add(rbBetween);
        rbBetween.setText("Between");
        rbBetween.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                rbBetweenActionPerformed(evt);
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
                        .addComponent(ckResWriteFDR)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(btnResFDRClear)
                        .addGap(0, 0, Short.MAX_VALUE))
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                        .addComponent(rbAll)
                        .addGap(18, 18, 18)
                        .addComponent(rbWithin)
                        .addGap(18, 18, 18)
                        .addComponent(rbBetween)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(btnWriteDB))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(ckResValidate)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(txtResValidate, javax.swing.GroupLayout.PREFERRED_SIZE, 25, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(ckResValidateOverwrite)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 75, Short.MAX_VALUE)
                        .addComponent(btnResValidateClear)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(btnResValidateClearAll))))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(ckResValidate)
                    .addComponent(txtResValidate, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(ckResValidateOverwrite)
                    .addComponent(btnResValidateClear)
                    .addComponent(btnResValidateClearAll))
                .addGap(18, 18, 18)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(btnWriteDB)
                        .addGap(1, 1, 1))
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(rbAll)
                            .addComponent(rbWithin)
                            .addComponent(rbBetween))
                        .addGap(18, 18, 18)))
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(ckResWriteFDR)
                    .addComponent(btnResFDRClear))
                .addGap(0, 45, Short.MAX_VALUE))
        );
    }// </editor-fold>//GEN-END:initComponents

    private void btnResValidateClearActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnResValidateClearActionPerformed

        clearValidation();
    }//GEN-LAST:event_btnResValidateClearActionPerformed

    private void btnWriteDBActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnWriteDBActionPerformed

        writeDB();

    }//GEN-LAST:event_btnWriteDBActionPerformed

    private void ckResValidateActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_ckResValidateActionPerformed

        ckResValidateOverwrite.setEnabled(ckResValidate.isSelected());
        txtResValidate.setEnabled(ckResValidate.isSelected());
        btnResValidateClear.setEnabled(ckResValidate.isSelected() && btnWriteDB.isEnabled());
        btnResValidateClearAll.setEnabled(ckResValidate.isSelected() && btnWriteDB.isEnabled());
    }//GEN-LAST:event_ckResValidateActionPerformed

    private void btnResValidateClearAllActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnResValidateClearAllActionPerformed
        clearValidation("");
    }//GEN-LAST:event_btnResValidateClearAllActionPerformed

    private void btnResFDRClearActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnResFDRClearActionPerformed

        clearDBFDR();
    }//GEN-LAST:event_btnResFDRClearActionPerformed

    private void ckResWriteFDRActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_ckResWriteFDRActionPerformed

        btnResFDRClear.setEnabled(ckResWriteFDR.isSelected() && m_gui.getFdr() instanceof DBinFDR);
    }//GEN-LAST:event_ckResWriteFDRActionPerformed

    private void rbBetweenActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_rbBetweenActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_rbBetweenActionPerformed

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JButton btnResFDRClear;
    private javax.swing.JButton btnResValidateClear;
    private javax.swing.JButton btnResValidateClearAll;
    private javax.swing.JButton btnWriteDB;
    private javax.swing.ButtonGroup buttonGroup1;
    private javax.swing.JCheckBox ckResValidate;
    private javax.swing.JCheckBox ckResValidateOverwrite;
    private javax.swing.JCheckBox ckResWriteFDR;
    private javax.swing.JRadioButton rbAll;
    private javax.swing.JRadioButton rbBetween;
    private javax.swing.JRadioButton rbWithin;
    private javax.swing.JTextField txtResValidate;
    // End of variables declaration//GEN-END:variables

    public void setFdrgui(FDRGUI gui) {
        this.m_gui= gui;
    }
}
