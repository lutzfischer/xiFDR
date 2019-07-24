/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.rappsilber.fdr.gui.components;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.concurrent.Semaphore;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import rappsilber.config.LocalProperties;
import rappsilber.utils.MyArrayUtils;

/**
 *
 * @author lfischer
 */
public class DBFIlters extends javax.swing.JFrame {

    public Semaphore finishedsp = new Semaphore(0);
    private ArrayList<String> runHistory = new ArrayList<String>();
    
    private String withinFilter = "(pepSeq2 is null OR ("
                    + "(accession1 = accession2) "
                    + "OR (lower(accession1) = 'rev_' || lower(accession2) OR lower(accession2) = 'rev_' || lower(accession1) )"
                    + "OR (lower(accession1) = 'ran_' || lower(accession2) OR lower(accession2) = 'ran_' || lower(accession1) )"
                    + "))";
    private String betweenFilter = "(pepSeq2 is null OR NOT ("
                    + "(accession1 = accession2) "
                    + "OR (lower(accession1) = 'rev_' || lower(accession2) OR lower(accession2) = 'rev_' || lower(accession1) )"
                    + "OR (lower(accession1) = 'ran_' || lower(accession2) OR lower(accession2) = 'ran_' || lower(accession1) )"
                    + "))";
    
    /**
     * Creates new form DBFIlters
     */
    public DBFIlters() {
        initComponents();
        stspMinDelta.setRanges(new Double(0), new Double(0), (Double)null, new Double(0.1));
        stspMinLinkSiteDelta.setRanges(new Double(0.0), new Double(0), (Double)null, new Double(0.01));
        stspMinFrags.setRanges(new Integer(0), new Integer(0), (Integer)null, new Integer(1));
        String prevruns =  LocalProperties.getProperty(this.getClass().getName()+".PreviousRuns", "");
        HashSet<String> ur = new HashSet<String>();
        String[] runs =prevruns.split("\\|\\|\\|");
        for (String r : runs ) {
            r=r.trim();
            if (!r.isEmpty() && !ur.contains(r)) {
                cbRun.addItem(r);
                runHistory.add(r);
                ur.add(r);
            }
        }
        cbRun.setSelectedIndex(-1);
        this.pWithinBetween.setVisible(false);
    }

    /**
     * Creates new form DBFIlters
     */
    public DBFIlters(String prevFilter) {
        this();
        Pattern p = Pattern.compile("\\(\\s*(?:run_name|lower\\s*\\(\\s*run_name\\s*\\))\\s+like\\s+ '%([^']*)%'\\s*\\)",Pattern.CASE_INSENSITIVE);
        Matcher m = p.matcher(prevFilter);
        if (m.find()) {
            cbRun.setSelectedItem(m.group(1));
        }
        
        String pep1 = "\\(\\s*pepSeq1\\s+like\\s+'%([^']*)%'\\s*\\)";
        String pep2 = "\\(\\s*pepSeq2\\slike\\s+'%([^']*)%'\\s*\\)";
        
        p = Pattern.compile(pep1 +"\\s+OR\\s+" + pep2,Pattern.CASE_INSENSITIVE);
        m = p.matcher(prevFilter);
        if (m.find()) {
            txtPep12.setText(m.group(1));
        }

        p = Pattern.compile(pep1 +"(?:\\s+AND\\s+|\\s*$)",Pattern.CASE_INSENSITIVE);
        m = p.matcher(prevFilter);
        if (m.find()) {
            txtPep1.setText(m.group(1));
        }

        p = Pattern.compile("(?:AND\\s|^)\\s*" + pep2,Pattern.CASE_INSENSITIVE);
        m = p.matcher(prevFilter);
        if (m.find()) {
            txtPep2.setText(m.group(1));
        }

        String prot1 = protinfilter(1, "{[^%]*}").replaceAll("\\(", "\\\\(").replaceAll("\\)", "\\\\)").replaceAll("\\{", "(").replaceAll("\\}", ")");
        String prot2 = protinfilter(2, "{[^%]*}").replaceAll("\\(", "\\\\(").replaceAll("\\)", "\\\\)").replaceAll("\\{", "(").replaceAll("\\}", ")");
        
        p = Pattern.compile(prot1 +"\\s+OR\\s+" + prot2,Pattern.CASE_INSENSITIVE);
        m = p.matcher(prevFilter);
        if (m.find()) {
            txtProt12.setText(m.group(1));
        }

        p = Pattern.compile(prot1 +"(?:\\s+AND\\s+|\\s*$)",Pattern.CASE_INSENSITIVE);
        m = p.matcher(prevFilter);
        if (m.find()) {
            txtProt1.setText(m.group(1));
        }

        p = Pattern.compile("(?:AND\\s|^)\\s*" + prot2,Pattern.CASE_INSENSITIVE);
        m = p.matcher(prevFilter);
        if (m.find()) {
            txtProt2.setText(m.group(1));
        }

        p = Pattern.compile("\\(\\s*scoreP1Coverage\\s+>=\\s+([0-9]+)\\s+AND\\s+",Pattern.CASE_INSENSITIVE);
        m = p.matcher(prevFilter);
        if (m.find()) {
            stspMinFrags.setValue(Integer.parseInt(m.group(1)));
        }
        p = Pattern.compile("\\(\\s*deltaScore\\s+>=\\s+([0-9]+(?:.?[0-9]*)?)(\\s*\\*\\s*score\\s*)?\\s*\\)",Pattern.CASE_INSENSITIVE);
        m = p.matcher(prevFilter);
        if (m.find()) {
            stspMinDelta.setValue(Double.parseDouble(m.group(1)));
            if (m.group(2) != null) {
                ckDeltaTimesScore.setSelected(true);
            }
        }

        p = Pattern.compile("\\(\\s*deltaScore\\s+>=\\s+score\\s*/\\s*([0-9]+(?:.?[0-9]*)?)\\s*\\)",Pattern.CASE_INSENSITIVE);
        m = p.matcher(prevFilter);
        if (m.find()) {
            stspMinDelta.setValue(1/Double.parseDouble(m.group(1)));
            ckDeltaTimesScore.setSelected(true);
        }

        p = Pattern.compile("\\(\\s*LinkSiteDelta\\s+>=\\s+([0-9]+(.?[0-9]*))\\s*\\)",Pattern.CASE_INSENSITIVE);
        m = p.matcher(prevFilter);
        if (m.find()) {
            stspMinLinkSiteDelta.setValue(Double.parseDouble(m.group(1)));
        }
        
        p = Pattern.compile("\\(\\s*precursor_charge\\s+>\\s+0\\s*\\)",Pattern.CASE_INSENSITIVE);
        m = p.matcher(prevFilter);
        if (m.find()) {
            ckSpecial.setSelected(true);
        }
        
        p = Pattern.compile("\\(\\s*autovalidated\\s*\\)",Pattern.CASE_INSENSITIVE);
        m = p.matcher(prevFilter);
        if (m.find()) {
            ckAutovalidated.setSelected(true);
        }
        
        p = Pattern.compile("\\(\\s*pepSeq2\\s+is\\s+not\\s+nul\\s*l\\)",Pattern.CASE_INSENSITIVE);
        m = p.matcher(prevFilter);
        if (m.find()) {
            ckXLOnly.setSelected(true);
        }
        
        p = Pattern.compile("\\(\\s*calc_charge\\s+>=\\s+([0-9]+)\\s*\\) ",Pattern.CASE_INSENSITIVE);
        m = p.matcher(prevFilter);
        if (m.find()) {
            stspMinPrecCharge.setValue(Integer.parseInt(m.group(1)));
        }

        p = Pattern.compile("\\(\\s*calc_charge\\s+<=\\s+([0-9]+)\\s*\\) ",Pattern.CASE_INSENSITIVE);
        m = p.matcher(prevFilter);
        if (m.find()) {
            stspMaxPrecCharge.setValue(Integer.parseInt(m.group(1)));
        }
        
        if (prevFilter.contains(withinFilter)) {
            rbWithin.setSelected(true);
        }
        
        if (prevFilter.contains(betweenFilter)) {
            rbBetween.setSelected(true);
        }

        p = Pattern.compile("\\((?:crosslinker|lower\\s*\\(\\s*crosslinker\\s*\\))\\s+like\\s+'%([^%]+)%'\\)",Pattern.CASE_INSENSITIVE);
        m = p.matcher(prevFilter);
        if (m.find()) {
            txtCrosslinker.setText(m.group(1));
        }

        p = Pattern.compile("\\(cleavclpep1fragmatched::int\\s+\\+\\s+cleavclpep1fragmatched::int\\s+>=\\s+([0-9]+)\\s*\\)",Pattern.CASE_INSENSITIVE);
        m = p.matcher(prevFilter);
        if (m.find()) {
            stspMinDSSOPep.setValue(Integer.parseInt(m.group(1)));
        }

        p = Pattern.compile("\\(mgxrank\\s+<=\\s+([0-9]+)\\s*\\)",Pattern.CASE_INSENSITIVE);
        m = p.matcher(prevFilter);
        if (m.find()) {
            stspMaxMGXRank.setValue(Integer.parseInt(m.group(1)));
        }

        p = Pattern.compile("\\(\\s*site1\\s*>\\s*0\\s+OR\\s+pepSeq2\\s+isnull\\s*\\)",Pattern.CASE_INSENSITIVE);
        m = p.matcher(prevFilter);
        if (m.find()) {
            cbNonCovalent.setSelected(true);
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

        bgCrosslinkTypes = new javax.swing.ButtonGroup();
        ckXLOnly = new javax.swing.JCheckBox();
        txtPep1 = new javax.swing.JTextField();
        txtPep2 = new javax.swing.JTextField();
        lblPep1 = new javax.swing.JLabel();
        lblPep2 = new javax.swing.JLabel();
        lblPep12 = new javax.swing.JLabel();
        txtPep12 = new javax.swing.JTextField();
        lblProt1 = new javax.swing.JLabel();
        lblProt2 = new javax.swing.JLabel();
        lblProt12 = new javax.swing.JLabel();
        txtProt12 = new javax.swing.JTextField();
        txtProt2 = new javax.swing.JTextField();
        txtProt1 = new javax.swing.JTextField();
        jButton1 = new javax.swing.JButton();
        stspMinFrags = new org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner();
        jLabel1 = new javax.swing.JLabel();
        jLabel2 = new javax.swing.JLabel();
        stspMinDelta = new org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner();
        jLabel3 = new javax.swing.JLabel();
        stspMinLinkSiteDelta = new org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner();
        ckSpecial = new javax.swing.JCheckBox();
        ckAutovalidated = new javax.swing.JCheckBox();
        jLabel4 = new javax.swing.JLabel();
        cbRun = new javax.swing.JComboBox();
        stspMinPrecCharge = new org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner();
        stspMaxPrecCharge = new org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner();
        lblPrecCharge = new javax.swing.JLabel();
        lblMinPrecCharge = new javax.swing.JLabel();
        lblMaxPrecCharge = new javax.swing.JLabel();
        pWithinBetween = new javax.swing.JPanel();
        rbAll = new javax.swing.JRadioButton();
        rbWithin = new javax.swing.JRadioButton();
        rbBetween = new javax.swing.JRadioButton();
        txtCrosslinker = new javax.swing.JTextField();
        lblPep3 = new javax.swing.JLabel();
        ckDeltaTimesScore = new javax.swing.JCheckBox();
        lblPrecCharge1 = new javax.swing.JLabel();
        stspMinDSSOPep = new org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner();
        lblPrecCharge2 = new javax.swing.JLabel();
        stspMaxMGXRank = new org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner();
        cbNonCovalent = new javax.swing.JCheckBox();

        addWindowListener(new java.awt.event.WindowAdapter() {
            public void windowClosed(java.awt.event.WindowEvent evt) {
                formWindowClosed(evt);
            }
            public void windowClosing(java.awt.event.WindowEvent evt) {
                formWindowClosing(evt);
            }
        });

        ckXLOnly.setText("Crosslinks Only");
        ckXLOnly.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                ckXLOnlyActionPerformed(evt);
            }
        });

        lblPep1.setText("Peptide1");

        lblPep2.setText("Peptide2");

        lblPep12.setText("Either");

        lblProt1.setText("Protein1");

        lblProt2.setText("Protein2");

        lblProt12.setText("Either");

        jButton1.setText("Ok");
        jButton1.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jButton1ActionPerformed(evt);
            }
        });

        jLabel1.setText("Min frags per peptide");

        jLabel2.setText("Minimum DeltaScore");

        stspMinDelta.setToolTipText("minimum delta score");

        jLabel3.setText("Minimum LinkSiteDelta");

        ckSpecial.setText("Exclude unknown charge state");

        ckAutovalidated.setText("autovalidated only");

        jLabel4.setText("Run name");

        cbRun.setEditable(true);

        lblPrecCharge.setText("Precursor Charge");

        lblMinPrecCharge.setText("min");

        lblMaxPrecCharge.setText("max");

        pWithinBetween.setBorder(javax.swing.BorderFactory.createTitledBorder("Cross-link types"));

        bgCrosslinkTypes.add(rbAll);
        rbAll.setSelected(true);
        rbAll.setText("All");

        bgCrosslinkTypes.add(rbWithin);
        rbWithin.setText("Within Only");

        bgCrosslinkTypes.add(rbBetween);
        rbBetween.setText("Between only");

        javax.swing.GroupLayout pWithinBetweenLayout = new javax.swing.GroupLayout(pWithinBetween);
        pWithinBetween.setLayout(pWithinBetweenLayout);
        pWithinBetweenLayout.setHorizontalGroup(
            pWithinBetweenLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(pWithinBetweenLayout.createSequentialGroup()
                .addComponent(rbAll)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(rbWithin)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(rbBetween))
        );
        pWithinBetweenLayout.setVerticalGroup(
            pWithinBetweenLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(pWithinBetweenLayout.createSequentialGroup()
                .addGroup(pWithinBetweenLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(rbAll)
                    .addComponent(rbWithin)
                    .addComponent(rbBetween))
                .addContainerGap(10, Short.MAX_VALUE))
        );

        lblPep3.setText("Cross-linker");

        ckDeltaTimesScore.setText("x score");
        ckDeltaTimesScore.setToolTipText("dektascore >= x * score");

        lblPrecCharge1.setText("DSSO Peps Found");

        stspMinDSSOPep.setModel(new javax.swing.SpinnerNumberModel(0, null, 3, 1));

        lblPrecCharge2.setText("max Candidate Pair Rank");

        stspMaxMGXRank.setModel(new javax.swing.SpinnerNumberModel(0, null, 3, 1));

        cbNonCovalent.setText("exclude non-covalent");

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                        .addGap(0, 0, Short.MAX_VALUE)
                        .addComponent(jButton1))
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(lblPep1)
                            .addComponent(lblProt1))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(txtPep1, javax.swing.GroupLayout.DEFAULT_SIZE, 84, Short.MAX_VALUE)
                            .addComponent(txtProt1))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(lblPep2)
                                .addGap(0, 0, Short.MAX_VALUE))
                            .addGroup(layout.createSequentialGroup()
                                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(jLabel1)
                                    .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                                        .addGroup(layout.createSequentialGroup()
                                            .addComponent(lblPrecCharge)
                                            .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                            .addComponent(lblMinPrecCharge))
                                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                            .addComponent(jLabel2, javax.swing.GroupLayout.PREFERRED_SIZE, 160, javax.swing.GroupLayout.PREFERRED_SIZE)
                                            .addComponent(jLabel3, javax.swing.GroupLayout.PREFERRED_SIZE, 169, javax.swing.GroupLayout.PREFERRED_SIZE)
                                            .addComponent(jLabel4)
                                            .addComponent(lblProt2)
                                            .addComponent(lblPep3)))
                                    .addComponent(lblPrecCharge1)
                                    .addComponent(lblPrecCharge2))
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                                    .addComponent(stspMaxMGXRank, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                    .addComponent(stspMinDSSOPep, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                    .addComponent(stspMinLinkSiteDelta, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, 110, Short.MAX_VALUE)
                                    .addComponent(stspMinDelta, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                    .addComponent(stspMinFrags, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                    .addComponent(txtProt2, javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(txtPep2, javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(txtCrosslinker)
                                    .addComponent(cbRun, 0, 1, Short.MAX_VALUE)
                                    .addComponent(stspMinPrecCharge, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addGroup(layout.createSequentialGroup()
                                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                            .addComponent(lblMaxPrecCharge, javax.swing.GroupLayout.Alignment.TRAILING)
                                            .addComponent(lblProt12, javax.swing.GroupLayout.Alignment.TRAILING)
                                            .addComponent(lblPep12, javax.swing.GroupLayout.Alignment.TRAILING))
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                            .addComponent(stspMaxPrecCharge, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, 109, Short.MAX_VALUE)
                                            .addComponent(txtProt12, javax.swing.GroupLayout.DEFAULT_SIZE, 109, Short.MAX_VALUE)
                                            .addComponent(txtPep12, javax.swing.GroupLayout.DEFAULT_SIZE, 109, Short.MAX_VALUE)))
                                    .addGroup(layout.createSequentialGroup()
                                        .addComponent(ckDeltaTimesScore)
                                        .addGap(0, 0, Short.MAX_VALUE))))))
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(ckSpecial)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                                .addComponent(ckAutovalidated, javax.swing.GroupLayout.PREFERRED_SIZE, 169, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(cbNonCovalent))
                            .addComponent(ckXLOnly)
                            .addComponent(pWithinBetween, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addGap(0, 0, Short.MAX_VALUE)))
                .addContainerGap())
        );

        layout.linkSize(javax.swing.SwingConstants.HORIZONTAL, new java.awt.Component[] {jLabel1, jLabel2, jLabel3, jLabel4, lblPep2, lblPep3, lblProt2});

        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(ckXLOnly)
                    .addComponent(lblPep3)
                    .addComponent(txtCrosslinker, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(18, 18, 18)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(txtPep12, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(lblPep12))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(txtProt12, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(lblProt12)))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(txtPep2, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(txtProt2, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(lblProt2)
                            .addComponent(txtProt1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(lblProt1))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(stspMinFrags, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(jLabel1))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(stspMinDelta, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(jLabel2)
                            .addComponent(ckDeltaTimesScore))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(stspMinLinkSiteDelta, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(jLabel3))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(cbRun, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(jLabel4))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(stspMinPrecCharge, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(lblPrecCharge)
                            .addComponent(lblMinPrecCharge)
                            .addComponent(stspMaxPrecCharge, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(lblMaxPrecCharge)))
                    .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        .addComponent(txtPep1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addComponent(lblPep1)
                        .addComponent(lblPep2)))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(stspMinDSSOPep, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(lblPrecCharge1))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(stspMaxMGXRank, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(lblPrecCharge2))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 50, Short.MAX_VALUE)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(ckSpecial)
                    .addComponent(ckAutovalidated)
                    .addComponent(cbNonCovalent))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(pWithinBetween, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jButton1))
                .addContainerGap())
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void ckXLOnlyActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_ckXLOnlyActionPerformed
        
    }//GEN-LAST:event_ckXLOnlyActionPerformed

    private void jButton1ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jButton1ActionPerformed
        this.setVisible(false);
        finishedsp.release();
    }//GEN-LAST:event_jButton1ActionPerformed

    private void formWindowClosed(java.awt.event.WindowEvent evt) {//GEN-FIRST:event_formWindowClosed
        finishedsp.release();
    }//GEN-LAST:event_formWindowClosed

    private void formWindowClosing(java.awt.event.WindowEvent evt) {//GEN-FIRST:event_formWindowClosing
        finishedsp.release();
    }//GEN-LAST:event_formWindowClosing

    
    public String getWhere() {
        ArrayList<String> filters = new ArrayList<String>();
        if (ckXLOnly.isSelected()) {
            filters.add("(pepSeq2 is not null)");
        }
        if (!txtPep1.getText().isEmpty()) {
            filters.add("(pepSeq1 LIKE '%" + txtPep1.getText() + "%')");
        }
        if (!txtPep2.getText().isEmpty()) {
            filters.add("(pepSeq2 LIKE '%" + txtPep2.getText() + "%')");
        }
        if (!txtPep12.getText().isEmpty()) {
            filters.add("((pepSeq1 LIKE '%" + txtPep12.getText() + "%') OR (pepSeq2 LIKE '%" + txtPep12.getText() + "%'))");
        }
        
        if (!txtProt1.getText().isEmpty()) {
            filters.add(protinfilter(1,txtProt1.getText()));
        }
        if (!txtProt2.getText().isEmpty()) {
            filters.add(protinfilter(2,txtProt2.getText()));
        }
        if (!txtProt12.getText().isEmpty()) {
            filters.add("(" + protinfilter(1,txtProt12.getText()) +" OR " + protinfilter(2,txtProt12.getText()) + ")");
        }

        if ((Integer)stspMinFrags.getValue() > 0) {
            filters.add("(scoreP1Coverage >= " + ((Integer)stspMinFrags.getValue()) + " AND (pepSeq2 is null OR  scoreP2Coverage >= " + ((Integer)stspMinFrags.getValue())+ " ))");
        }

        if ((Double)stspMinDelta.getValue() > 0) {
            filters.add("(deltaScore >= " + ((Double)stspMinDelta.getValue()) + (ckDeltaTimesScore.isSelected()? "* score":"") + ")");
        }
        if ((Double)stspMinLinkSiteDelta.getValue() > 0) {
            filters.add("(LinkSiteDelta >= " + ((Double)stspMinLinkSiteDelta.getValue()) +")");
        }
        
        if(ckSpecial.isSelected()) {
            filters.add("(precursor_charge > 0)");
        }
        
        if(ckAutovalidated.isSelected()) {
            filters.add("(autovalidated)");
        }
        if(cbRun.getSelectedItem() != null && !cbRun.getSelectedItem().toString().isEmpty()) {
            String txt = cbRun.getSelectedItem().toString();
            if (txt.toLowerCase().contentEquals(txt)) {
                filters.add("(LOWER(run_name) LIKE '%"+txt+"%')");
            }else 
                filters.add("(run_name LIKE '%"+txt+"%')");
        }

        if ((Integer)stspMinPrecCharge.getValue() > 0) {
            filters.add("(calc_charge >= " + ((Integer)stspMinPrecCharge.getValue()) +")");
        }

        if ((Integer)stspMaxPrecCharge.getValue() > 0) {
            filters.add("(calc_charge <= " + ((Integer)stspMaxPrecCharge.getValue()) +")");
        }

        if (txtCrosslinker.getText() != null && !txtCrosslinker.getText().isEmpty()) {
            String txt = txtCrosslinker.getText();
            if (txt.toLowerCase().contentEquals(txt)) {
                filters.add("(LOWER(crosslinker) LIKE '%"+txt +"%')");        
            } else {
                filters.add("(crosslinker LIKE '%"+txt +"%')");        
            }
        }
        
        if ((Integer)stspMinDSSOPep.getValue() > 0) {
            filters.add("(cleavclpep1fragmatched::int + cleavclpep1fragmatched::int >= " + ((Integer)stspMinDSSOPep.getValue()) +")");
        }
        
        
        if ((Integer)stspMaxMGXRank.getValue() > 0) {
            filters.add("(mgxrank <= " + ((Integer)stspMaxMGXRank.getValue()) +")");
        }

        if (rbWithin.isSelected()) {
            filters.add(withinFilter);
        }
        
        if (rbBetween.isSelected()) {
            filters.add(betweenFilter);
        }
        
        if (cbNonCovalent.isSelected()) {
            filters.add("(site1>0 OR pepSeq2 isnull)");
        }
        
        return MyArrayUtils.toString(filters, " AND ");

    }
    
    
    public static String showAndGetFilter() {
        DBFIlters dbf = new DBFIlters();
        dbf.setVisible(true);
        boolean ok = false;
        while (!ok)
            try {
                dbf.finishedsp.acquire();
                ok = true;
            } catch (InterruptedException ex) {
                Logger.getLogger(DBFIlters.class.getName()).log(Level.SEVERE, null, ex);
            }
        
        String ret = dbf.getWhere();
        HashSet<String> runset = new HashSet<String>();
        ArrayList<String> rh = new ArrayList<String>();
        if (dbf.cbRun.getSelectedItem() != null) {
            dbf.runHistory.add(dbf.cbRun.getSelectedItem().toString());
            while (dbf.runHistory.size() > 10) {
                dbf.runHistory.remove(0);
            }
            LocalProperties.setProperty(dbf.getClass().getName()+".PreviousRuns", MyArrayUtils.toString(dbf.runHistory, "|||"));
        }
        dbf.dispose();
        return ret;
    }

    public static String showAndGetFilter(String prevFilter) {
        DBFIlters dbf = new DBFIlters(prevFilter);
        dbf.setVisible(true);
        boolean ok = false;
        while (!ok)
            try {
                dbf.finishedsp.acquire();
                ok = true;
            } catch (InterruptedException ex) {
                Logger.getLogger(DBFIlters.class.getName()).log(Level.SEVERE, null, ex);
            }

        if (dbf.cbRun.getSelectedItem() != null) {
            dbf.runHistory.add(dbf.cbRun.getSelectedItem().toString());
            while (dbf.runHistory.size() > 10) {
                dbf.runHistory.remove(0);
            }
            LocalProperties.setProperty(dbf.getClass().getName()+".PreviousRuns", MyArrayUtils.toString(dbf.runHistory, "|||"));
        }        
        String ret = dbf.getWhere();
        dbf.dispose();
        return ret;
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

        /* Create and display the form */
        java.awt.EventQueue.invokeLater(new Runnable() {
            public void run() {
                new DBFIlters().setVisible(true);
            }
        });
    }
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.ButtonGroup bgCrosslinkTypes;
    private javax.swing.JCheckBox cbNonCovalent;
    private javax.swing.JComboBox cbRun;
    private javax.swing.JCheckBox ckAutovalidated;
    private javax.swing.JCheckBox ckDeltaTimesScore;
    private javax.swing.JCheckBox ckSpecial;
    private javax.swing.JCheckBox ckXLOnly;
    private javax.swing.JButton jButton1;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel3;
    private javax.swing.JLabel jLabel4;
    private javax.swing.JLabel lblMaxPrecCharge;
    private javax.swing.JLabel lblMinPrecCharge;
    private javax.swing.JLabel lblPep1;
    private javax.swing.JLabel lblPep12;
    private javax.swing.JLabel lblPep2;
    private javax.swing.JLabel lblPep3;
    private javax.swing.JLabel lblPrecCharge;
    private javax.swing.JLabel lblPrecCharge1;
    private javax.swing.JLabel lblPrecCharge2;
    private javax.swing.JLabel lblProt1;
    private javax.swing.JLabel lblProt12;
    private javax.swing.JLabel lblProt2;
    private javax.swing.JPanel pWithinBetween;
    private javax.swing.JRadioButton rbAll;
    private javax.swing.JRadioButton rbBetween;
    private javax.swing.JRadioButton rbWithin;
    private org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner stspMaxMGXRank;
    private org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner stspMaxPrecCharge;
    private org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner stspMinDSSOPep;
    private org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner stspMinDelta;
    private org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner stspMinFrags;
    private org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner stspMinLinkSiteDelta;
    private org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner stspMinPrecCharge;
    private javax.swing.JTextField txtCrosslinker;
    private javax.swing.JTextField txtPep1;
    private javax.swing.JTextField txtPep12;
    private javax.swing.JTextField txtPep2;
    private javax.swing.JTextField txtProt1;
    private javax.swing.JTextField txtProt12;
    private javax.swing.JTextField txtProt2;
    // End of variables declaration//GEN-END:variables

    protected String protinfilter(int Protein,String text) {
        String filter = "(accession"+ Protein + " like '%" + text + "%' OR protein"+ Protein + "name like '%" + text + "%' OR description"+ Protein + " like '%" + text + "%')";
        return filter;
    }
}
