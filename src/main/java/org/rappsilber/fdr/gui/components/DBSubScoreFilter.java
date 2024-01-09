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

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Timer;
import java.util.TimerTask;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JFrame;
import org.rappsilber.config.LocalProperties;
import org.rappsilber.fdr.DBinFDR;
import org.rappsilber.utils.RArrayUtils;

/**
 *
 * @author Lutz Fischer <lfischer@staffmail.ed.ac.uk>
 */
public class DBSubScoreFilter extends javax.swing.JPanel {
    
    ArrayList<String> subscores;
    
    public static Pattern filterPattern = Pattern.compile("\\(\\s*(@?)\\[%((?:[^%])*)%\\]\\s*([<=>]+|is(?=\\s+NaN)|is not(?=\\s+NaN))\\s*((?<=is|is not)\\s*NaN|[+-]?[0-9.]*\\s*)\\)\\s*");

    /**
     * Creates new form DBSubScoreFilters
     */
    public DBSubScoreFilter() {
        initComponents();
    }

    /**
     * Creates new form DBSubScoreFilters
     */
    public DBSubScoreFilter(String filter) {
        initComponents();
        setFilter(filter);
               
    }

    public void setFilter(String filter) throws NumberFormatException {
        Matcher m = filterPattern.matcher(filter);
        if (m.matches()) {
            m.group(1);
            this.cmbScore.setSelectedItem(m.group(2));
            String op = m.group(3);
            if (op.contentEquals("is")) {
                this.cmbOperator.setSelectedItem("is NaN");
                this.spValue.setEnabled(false);
            } else if (op.contentEquals("is not")) {
                this.cmbOperator.setSelectedItem("is not NaN");
                this.spValue.setEnabled(false);
            } else {
                this.cmbOperator.setSelectedItem(op);
                this.spValue.setValue(Double.parseDouble(m.group(4)));
            }
            ckAbsolute.setSelected(m.group(1).length()>0);
        }
    }
    
  
    public void setSubScores(ArrayList<String> subscores) {
        Object cs =  cmbScore.getSelectedItem();
        String[] prevSel = LocalProperties.getProperty(DBinFDR.subscoreroperty,"").split(";");
        final HashMap<String,Integer>  prevSelected = new HashMap<>();
        for (int i = 0; i<prevSel.length;i++) {
            prevSelected.put(prevSel[i], i);
        }
        
        if (subscores != null) {
            subscores = new ArrayList<>(subscores);
            Collections.sort(subscores, new Comparator<String>(){
                @Override
                public int compare(String o1, String o2) {
                    
                    int s1 = 0;
                    int s2 = 0;
                    if (o1.toLowerCase().contains("error")) {
                        s1-=1;
                    }
                    if (o2.toLowerCase().contains("error")) {
                        s2-=1;
                    }
                    if (o1.toLowerCase().matches(".*preco?uo?rsor")) {
                        s1-=1;
                    }
                    if (o2.toLowerCase().matches(".*preco?uo?rsor")) {
                        s2-=1;
                    }
                    
                    // if something was previously selected
                    Integer p1 = prevSelected.get(o1);
                    if (p1 != null)
                        s1 -= 100 + prevSelected.size() -  p1;
                    Integer p2 = prevSelected.get(o2);
                    if (p2 != null)
                        s2 -= 100  + prevSelected.size() -  p2;
                    
                    
                    
                    if (s1==s2)
                        return o1.toLowerCase().compareTo(o2.toLowerCase());
                    else
                        return s1 - s2;
                    
                }
                
            });
            cmbScore.setModel(new DefaultComboBoxModel<String>(subscores.toArray(new String[subscores.size()])));
        } else {
            cmbScore.setModel(new DefaultComboBoxModel<String>(new String[0]));
        }
        if (cs != null && ((String)cs).trim().length() > 0) {
            cmbScore.setSelectedItem(cs);
        }
    }
    
    public String getFilter() {
        StringBuilder filter = new StringBuilder("(");
        if (ckAbsolute.isSelected()) {
            filter.append("@");
        }
        filter.append("[%").append(cmbScore.getSelectedItem().toString()).append("%] ");
        String op = cmbOperator.getSelectedItem().toString();
        if (op.contentEquals("is NaN")) {
            filter.append("= 'NaN'");
        } else if (op.contentEquals("is not NaN")) {
            filter.append("<> 'NaN'");
        } else {
            filter.append(op).append(" ");
            filter.append(spValue.getValue().toString());
        }
        filter.append(")");
        
        return filter.toString();
    }
    

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        cmbScore = new javax.swing.JComboBox<>();
        cmbOperator = new javax.swing.JComboBox<>();
        ckAbsolute = new javax.swing.JCheckBox();
        spValue = new javax.swing.JSpinner();

        cmbScore.setEditable(true);

        cmbOperator.setModel(new javax.swing.DefaultComboBoxModel<>(new String[] { "=", ">", ">=", "<", "<=", "<>", "is NaN", "is not NaN" }));
        cmbOperator.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                cmbOperatorActionPerformed(evt);
            }
        });

        ckAbsolute.setText("absolute");

        spValue.setModel(new javax.swing.SpinnerNumberModel(0.0f, null, null, 0.5f));

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGap(0, 0, 0)
                .addComponent(ckAbsolute)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(cmbScore, 0, 168, Short.MAX_VALUE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(cmbOperator, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(spValue, javax.swing.GroupLayout.DEFAULT_SIZE, 53, Short.MAX_VALUE)
                .addGap(2, 2, 2))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGap(0, 0, 0)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(cmbScore, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(cmbOperator, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(ckAbsolute)
                    .addComponent(spValue, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
        );
    }// </editor-fold>//GEN-END:initComponents

    private void cmbOperatorActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_cmbOperatorActionPerformed
        spValue.setEnabled(!cmbOperator.getSelectedItem().toString().contains("NaN"));
    }//GEN-LAST:event_cmbOperatorActionPerformed


    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JCheckBox ckAbsolute;
    private javax.swing.JComboBox<String> cmbOperator;
    private javax.swing.JComboBox<String> cmbScore;
    private javax.swing.JSpinner spValue;
    // End of variables declaration//GEN-END:variables

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
            java.util.logging.Logger.getLogger(DBFIlters.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (InstantiationException ex) {
            java.util.logging.Logger.getLogger(DBFIlters.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (IllegalAccessException ex) {
            java.util.logging.Logger.getLogger(DBFIlters.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (javax.swing.UnsupportedLookAndFeelException ex) {
            java.util.logging.Logger.getLogger(DBFIlters.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>
        //</editor-fold>
        System.out.println("\\(\\s*\\[%((?:[^%][^\\]])*)%\\]\\s*([<=>]+|is|is not)\\s*((?<=is|is not)\\s*NaN|[+-][0-9.]*\\s*\\))\\)");
        final DBSubScoreFilter dbsf = new DBSubScoreFilter();

        /* Create and display the form */
        java.awt.EventQueue.invokeLater(new Runnable() {
            public void run() {
                JFrame testframe = new JFrame("Test");
                dbsf.setSubScores(new ArrayList<String>(RArrayUtils.toCollection(new String[]{"A", "B"})));
                dbsf.setFilter("(@[%C%] > -1)");
                testframe.getContentPane().add(dbsf);
                testframe.pack();
                testframe.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
                testframe.setVisible(true);
            }
        });
        
        Timer t = new Timer(true);
        t.schedule(new TimerTask() {
            @Override
            public void run() {
                System.out.println(dbsf.getFilter());
            }
        }, 1000, 1000);
    }

    
}
