/*
 * Copyright 2019 Lutz Fischer <lfischer@staffmail.ed.ac.uk>.
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

import java.awt.BorderLayout;
import java.awt.GraphicsEnvironment;
import java.awt.HeadlessException;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.lang.invoke.MethodHandles;
import java.util.Locale;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;
import org.rappsilber.config.LocalProperties;
import org.rappsilber.fdr.utils.MZIdentMLOwner;

/**
 *
 * @author Lutz Fischer <lfischer@staffmail.ed.ac.uk>
 */
public class MZIdentMLOwnerGUI extends javax.swing.JPanel {

    /**
     * Creates new form MzIdenMLOwner
     */
    public MZIdentMLOwnerGUI() {
        initComponents();
        setMZIdentMLOwner();
    }

    
    public void setMZIdentMLOwner() {
        txtmzIdentOwnerFirst.setText(LocalProperties.getProperty("mzIdenMLOwnerFirst", txtmzIdentOwnerFirst.getText()));
        txtmzIdentOwnerLast.setText(LocalProperties.getProperty("mzIdenMLOwnerLast", txtmzIdentOwnerLast.getText()));
        txtmzIdentOwnerEmail.setText(LocalProperties.getProperty("mzIdenMLOwnerEmail", txtmzIdentOwnerEmail.getText()));
        txtmzIdentAddress.setText(LocalProperties.getProperty("mzIdenMLOwnerAddress", txtmzIdentAddress.getText()));
        txtmzIdentOwnerOrg.setText(LocalProperties.getProperty("mzIdenMLOwnerOrg", txtmzIdentOwnerOrg.getText()));
    }

    public void setMZIdentMLOwner(MZIdentMLOwner owner) {
        txtmzIdentOwnerFirst.setText(owner.first);
        txtmzIdentOwnerLast.setText(owner.last);
        txtmzIdentOwnerEmail.setText(owner.email);
        txtmzIdentAddress.setText(owner.address);
        txtmzIdentOwnerOrg.setText(owner.org);
    }
    
    public static MZIdentMLOwner getLastOwner() {
        return new MZIdentMLOwner(
                LocalProperties.getProperty("mzIdenMLOwnerFirst", ""),
                LocalProperties.getProperty("mzIdenMLOwnerLast", ""),
                LocalProperties.getProperty("mzIdenMLOwnerEmail", ""),
                LocalProperties.getProperty("mzIdenMLOwnerOrg", ""),
                LocalProperties.getProperty("mzIdenMLOwnerAddress", ""));
        
    }

    public static void saveDefaultOwner(MZIdentMLOwner owner) {
        LocalProperties.setProperty("mzIdenMLOwnerFirst", owner.first);
        LocalProperties.setProperty("mzIdenMLOwnerLast", owner.last);
        LocalProperties.setProperty("mzIdenMLOwnerEmail", owner.email);
        LocalProperties.setProperty("mzIdenMLOwnerAddress", owner.address);
        LocalProperties.setProperty("mzIdenMLOwnerOrg", owner.org);
    }


    public static void resetDefaultOwner() {
        LocalProperties.setProperty("mzIdenMLOwnerFirst", null);
        LocalProperties.setProperty("mzIdenMLOwnerLast", null);
        LocalProperties.setProperty("mzIdenMLOwnerEmail", null);
        LocalProperties.setProperty("mzIdenMLOwnerAddress", null);
        LocalProperties.setProperty("mzIdenMLOwnerOrg", null);
    }

    
    public MZIdentMLOwner getOwnerInfo() {
        return new MZIdentMLOwner(txtmzIdentOwnerFirst.getText(), txtmzIdentOwnerLast.getText(), txtmzIdentOwnerEmail.getText(), txtmzIdentOwnerOrg.getText(), txtmzIdentAddress.getText());
    }
    
    public static MZIdentMLOwner getDocumentOwner() {
        MZIdentMLOwner o = getLastOwner();
        if (o.isEmpty()) {
            return getNewOwner(o);
        }
        return o;
    }

    public static MZIdentMLOwner getSetOwner(MZIdentMLOwner defaultOwner) throws HeadlessException, UnsupportedOperationException {
        MZIdentMLOwner o = getNewOwner(defaultOwner);
        saveDefaultOwner(o);
        return o;
    }
    
    public static MZIdentMLOwner getNewOwner(MZIdentMLOwner defaultOwner) throws HeadlessException, UnsupportedOperationException {
        if (GraphicsEnvironment.isHeadless()) {
            Logger.getLogger(MethodHandles.lookup().lookupClass().getName()).log(Level.WARNING, "NO OWNER PREDEFINED AND CAN'T ASK FOR ONE!");
            return defaultOwner;
        }
        if (SwingUtilities.isEventDispatchThread())  {
            Logger.getLogger(MethodHandles.lookup().lookupClass().getName()).log(Level.WARNING, "NO OWNER PREDEFINED AND CAN'T ASK FOR ONE!");
            throw new UnsupportedOperationException("This Method should not be run from the EventDispachThread");
        }
        final Object lock = new Object();
        final JFrame window = new JFrame("mzIdentML Document Owner");
        JLabel desc = new JLabel("mzIdentML-files need to provide information about the owner of the document");
        MZIdentMLOwnerGUI og = new MZIdentMLOwnerGUI();
        window.getContentPane().setLayout(new BorderLayout());
        window.add(desc,BorderLayout.NORTH);
        window.add(og,BorderLayout.CENTER);
        JPanel pButtons = new JPanel(new BorderLayout());
        window.getContentPane().add(pButtons,BorderLayout.SOUTH);
        JButton close = new JButton("OK");
        pButtons.add(close,BorderLayout.EAST);
        window.pack();
        close.addActionListener(new ActionListener(){
            @Override
            public void actionPerformed(ActionEvent e) {
                window.setVisible(false);
                synchronized (lock) {
                    lock.notify();
                }
            }
        });
        window.addWindowListener(new WindowAdapter() {
            
            @Override
            public void windowClosing(WindowEvent arg0) {
                synchronized (lock) {
                    //window.setVisible(false);
                    lock.notify();
                }
            }
            
        });
        window.setVisible(true);
        Thread t = new Thread(new Runnable() {
            public void run() {
                while (window.isVisible()) {
                    try {
                        synchronized(lock) {
                            lock.wait();
                        }
                    } catch (InterruptedException e) {
                        e.printStackTrace();
                    }
                    Logger.getLogger(MZIdentMLOwnerGUI.class.getName()).log(Level.INFO,"Window Closed");
                }
            }
        });
        t.start();
        try {
            t.join();
        } catch (InterruptedException ex) {
            Logger.getLogger(MZIdentMLOwnerGUI.class.getName()).log(Level.SEVERE, null, ex);
        }
        MZIdentMLOwner ret = og.getCurrentDocumentOwner();
        return ret;
    }
    
    public MZIdentMLOwner getCurrentDocumentOwner() {
        return new MZIdentMLOwner(
                txtmzIdentOwnerFirst.getText(),
                txtmzIdentOwnerLast.getText(),
                txtmzIdentOwnerEmail.getText(),
                txtmzIdentOwnerOrg.getText(),
                txtmzIdentAddress.getText());
    }
    
    public static void main(String[] args) {
        getNewOwner(getLastOwner());
    }
    
    
    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jLabel2 = new javax.swing.JLabel();
        txtmzIdentOwnerFirst = new javax.swing.JTextField();
        txtmzIdentOwnerLast = new javax.swing.JTextField();
        txtmzIdentOwnerEmail = new javax.swing.JTextField();
        jLabel3 = new javax.swing.JLabel();
        jLabel5 = new javax.swing.JLabel();
        txtmzIdentOwnerOrg = new javax.swing.JTextField();
        jLabel4 = new javax.swing.JLabel();
        jScrollPane1 = new javax.swing.JScrollPane();
        txtmzIdentAddress = new javax.swing.JTextArea();

        jLabel2.setText("Document Owner");

        txtmzIdentOwnerFirst.setText("First");
        txtmzIdentOwnerFirst.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                txtmzIdentOwnerFirstActionPerformed(evt);
            }
        });

        txtmzIdentOwnerLast.setText("Last");
        txtmzIdentOwnerLast.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                txtmzIdentOwnerLastActionPerformed(evt);
            }
        });

        txtmzIdentOwnerEmail.setText("reseacher@organisation.org");
        txtmzIdentOwnerEmail.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                txtmzIdentOwnerEmailActionPerformed(evt);
            }
        });

        jLabel3.setText("E-Mail");

        jLabel5.setText("Organisation");

        txtmzIdentOwnerOrg.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                txtmzIdentOwnerOrgActionPerformed(evt);
            }
        });

        jLabel4.setText("Address");

        txtmzIdentAddress.setColumns(20);
        txtmzIdentAddress.setRows(5);
        jScrollPane1.setViewportView(txtmzIdentAddress);

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jLabel2)
                    .addComponent(jLabel3)
                    .addComponent(jLabel4)
                    .addComponent(jLabel5))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jScrollPane1)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(txtmzIdentOwnerFirst)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(txtmzIdentOwnerLast)
                        .addGap(3, 3, 3))
                    .addComponent(txtmzIdentOwnerEmail)
                    .addComponent(txtmzIdentOwnerOrg)))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel2)
                    .addComponent(txtmzIdentOwnerFirst, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(txtmzIdentOwnerLast, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel3)
                    .addComponent(txtmzIdentOwnerEmail, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(txtmzIdentOwnerOrg, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel5))
                .addGap(9, 9, 9)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jScrollPane1)
                    .addComponent(jLabel4)))
        );
    }// </editor-fold>//GEN-END:initComponents

    private void txtmzIdentOwnerFirstActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_txtmzIdentOwnerFirstActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_txtmzIdentOwnerFirstActionPerformed

    private void txtmzIdentOwnerLastActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_txtmzIdentOwnerLastActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_txtmzIdentOwnerLastActionPerformed

    private void txtmzIdentOwnerEmailActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_txtmzIdentOwnerEmailActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_txtmzIdentOwnerEmailActionPerformed

    private void txtmzIdentOwnerOrgActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_txtmzIdentOwnerOrgActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_txtmzIdentOwnerOrgActionPerformed


    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel3;
    private javax.swing.JLabel jLabel4;
    private javax.swing.JLabel jLabel5;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JTextArea txtmzIdentAddress;
    private javax.swing.JTextField txtmzIdentOwnerEmail;
    private javax.swing.JTextField txtmzIdentOwnerFirst;
    private javax.swing.JTextField txtmzIdentOwnerLast;
    private javax.swing.JTextField txtmzIdentOwnerOrg;
    // End of variables declaration//GEN-END:variables
}
