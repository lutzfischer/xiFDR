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
package org.rappsilber.fdr.gui;

import java.awt.Component;
import java.awt.Container;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.io.File;
import java.io.FileNotFoundException;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.EventObject;
import java.util.HashSet;
import java.util.logging.Filter;
import java.util.logging.Level;
import java.util.logging.LogRecord;
import java.util.logging.Logger;
import javax.swing.DefaultCellEditor;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.table.TableModel;
import org.rappsilber.config.LocalProperties;
import org.rappsilber.data.csv.CsvParser;
import org.rappsilber.fdr.CSVinFDR;
import org.rappsilber.fdr.result.FDRResult;
import org.rappsilber.fdr.result.FDRResultLevel;
import org.rappsilber.fdr.MZIdentXLFDR;
import org.rappsilber.fdr.OfflineFDR;
import org.rappsilber.fdr.result.SubGroupFdrInfo;
import org.rappsilber.fdr.entities.PSM;
import org.rappsilber.fdr.entities.PeptidePair;
import org.rappsilber.fdr.entities.ProteinGroupLink;
import org.rappsilber.fdr.entities.ProteinGroupPair;
import org.rappsilber.fdr.gui.components.FDRSettingsPanel;
import org.rappsilber.fdr.utils.MiscUtils;
import org.rappsilber.gui.components.AutoAddTableModelListener;
import org.rappsilber.gui.components.GenericTextPopUpMenu;
import org.rappsilber.gui.components.JoinedThreadedTextOuput;
import org.rappsilber.utils.MyArrayUtils;
import org.rappsilber.gui.logging.JTextAreaHandle;
        

/**
 *
 * @author lfischer
 */
public class FDRGUI extends javax.swing.JFrame {

    private boolean stopMaximizing = false;
    
    private JoinedThreadedTextOuput m_status = new JoinedThreadedTextOuput();
    
    private JTextAreaHandle loggingOutput;
    
    
    private FDRResult m_result;
    
    private FDRSettingsPanel fdrSettings;
    
    private OfflineFDR m_fdr;
//    /** denotes required fields in the CSV */
//    public String missingColumn = "!! !MISSING! !!";
//    public String optionalColumn = "  OPTIONAL   ";
//    public String[] csvColumns = new String[]{missingColumn};
//    public String[] csvColumnsOptional = new String[]{optionalColumn};


    /**
     * Creates new form
     */
    public FDRGUI() {
        initComponents();

        // setup the logging to the text-field
        loggingOutput = new JTextAreaHandle(txtLog);
        loggingOutput.setFilter(new Filter() {

            public boolean isLoggable(LogRecord record) {
                return true;
            }
        });
        loggingOutput.setLevel(Level.ALL);

        // make sure we display one of the possible fdrsettings-panel
        changeFDRSettings(null);

        
        //Logger.getLogger("rappsilber").addHandler(loggingOutput);
        Logger.getLogger("org.rappsilber").setLevel(Level.ALL);
        Logger.getLogger("org.rappsilber").addHandler(loggingOutput);
        
        
        this.setTitle("");   
        txtXiFDRVersion.setText(OfflineFDR.version.toString());
        //tpInput.removeTabAt(0);
        tblPepLength.getModel().addTableModelListener(new AutoAddTableModelListener());


        // setup a context-menu with copy past-functions for all test-fields
        GenericTextPopUpMenu gtpm = new GenericTextPopUpMenu();
        gtpm.installContextMenu(this);
        gtpm.installContextMenu(fdrSettingsComplete);

        this.toFront();
        
        // register, where we want to have status messages printed
        m_status.addLoggerOutput(Logger.getLogger(this.getClass().getName()));
        m_status.addTextOutput(txtStatus);
        

        
//        fdrSettings=fdrSettingsSimple;
        
        
        
        // make spinners select all test on gaining the focus
        // it's a bit of a dirty hack - but the best I could find
        ArrayList<Container> con = new ArrayList<Container>();
        HashSet<Container> spinSet = new HashSet<Container>();
        ArrayList<Container> spins = new ArrayList<Container>();
        con.add(this);
        con.add(fdrSettingsComplete);
        con.add(fdrSettingsSimple);
        // get all the spinners in all containers recursivly
        for (int i =0 ; i< con.size(); i++) {
            for (Component c : con.get(i).getComponents()) {
                if (c instanceof JSpinner){
                    if (!spinSet.contains(c)) {
                        spins.add((Container) c);
                        spinSet.add((Container) c);
                    }
                } else if (c instanceof Container) {
                    con.add((Container) c);              
                }
            }
        }
        // now install the a handler on all text-fields within all spinners that selects all 
        for (int i = 0; i< spins.size(); i++) {
            Container spin = spins.get(i);
            for (Component c : spin.getComponents()) {
                if (c instanceof JTextField){
                    ((JTextField) c).addFocusListener(new java.awt.event.FocusAdapter() {
                        public void focusGained(final java.awt.event.FocusEvent evt) {
                            javax.swing.SwingUtilities.invokeLater(new Runnable() {
                                public void run() {
                                    // on focus gain select all text
                                    ((JTextField)evt.getSource()).selectAll();
                                }
                            });
                        }
                    });
                } else if (c instanceof Container) {
                    if (!spinSet.contains(c)) {
                        spins.add((Container) c);
                        spinSet.add((Container) c);
                    }
                }
            }
        }
        
        // hide the DB-sizes panel
        ckDBSizeActionPerformed(null);
        // hide the FDR-groups panel
        ckDefineGroupsActionPerformed(null);
        
        // preset minimum peptide-length
        fdrSettingsComplete.setMinPeptideLength(6);
        fdrSettingsSimple.setMinPeptideLength(6);
        
        ActionListener stopMaximizingAl = new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                stopMaximizing = true;
            }
        };
          
        fdrSettingsComplete.addStopMaximizingListener(stopMaximizingAl);
        fdrSettingsSimple.addStopMaximizingListener(stopMaximizingAl);
        
        fdrSettingsComplete.setMinPeptideLength(6);
        fdrSettingsSimple.setMinPeptideLength(6);
        fdrSettingsComplete.setReportFactor(100000);
        fdrSettingsSimple.setReportFactor(100000);
        
        
        ActionListener calcListener = new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                basicCalc();
            }
        };
        
        fdrSettingsSimple.addCalcListener(calcListener);

        fdrSettingsComplete.addCalcListener(calcListener);
        
        
        spDistanceGroup.setSpecialValueText("No distance group");
        spDistanceGroup.setValue((Integer)0);

        
        txtmzIdentOwnerFirst.setText(LocalProperties.getProperty("mzIdenMLOwnerFirst",txtmzIdentOwnerFirst.getText()));
        txtmzIdentOwnerLast.setText(LocalProperties.getProperty("mzIdenMLOwnerLast",txtmzIdentOwnerLast.getText()));
        txtmzIdentOwnerEmail.setText(LocalProperties.getProperty("mzIdenMLOwnerEmail",txtmzIdentOwnerEmail.getText()));
        txtmzIdentAdress.setText(LocalProperties.getProperty("mzIdenMLOwnerAddress",txtmzIdentAdress.getText()));
        txtmzIdentOwnerOrg.setText(LocalProperties.getProperty("mzIdenMLOwnerAddress",txtmzIdentOwnerOrg.getText()));
        
    }
    
    
    public void setTitle(String title) {
        super.setTitle("xiFDR (" + OfflineFDR.version + ")" + title);
        
    }
    
    public void prepareFDRCalculation() throws NumberFormatException {
        ProteinGroupLink.MIN_DISTANCE_FOR_LONG = ((Number)spDistanceGroup.getValue()).intValue();
        getFdr().setMaximumLinkAmbiguity(fdrSettings.getMaxLinkAmbiguity());
        getFdr().setMaximumProteinAmbiguity(fdrSettings.getMaxProteinAmbiguity());
        getFdr().setMinimumPeptideLength(fdrSettings.getMinPeptideLength());
        if (ckDBSize.isSelected()) {
            getFdr().setTargetDBSize((Double) spTargetDB.getValue());
            getFdr().setDecoyDBSize((Double) spDecoyDB.getValue());
            getFdr().setTargetProtDBSize((Double) spTargetDBProt.getValue());
            getFdr().setDecoyProtDBSize((Double) spDecoyDBProt.getValue());
            getFdr().setTargetProtDBSize((Double) spTargetDBProt.getValue());
            getFdr().setDecoyProtDBSize((Double) spDecoyDBProt.getValue());
        } else {
            getFdr().setTargetDBSize(999999999);
            getFdr().setDecoyDBSize(999999999);
            getFdr().setTargetProtDBSize(999999999);
            getFdr().setDecoyProtDBSize(999999999);
        }
        
        getFdr().setPsm_directional(fdrSettings.isPSMDirectional());
        getFdr().setPeptides_directional(fdrSettings.isPeptidePairDirectional());
        getFdr().setLinks_directional(fdrSettings.isLinkDirectional());
        getFdr().setPpi_directional(fdrSettings.isPPIDirectional());
        
        

        getFdr().setMinPepPerProteinGroup(fdrSettings.getMinProteinPepCount());
        getFdr().setMinPepPerProteinGroupLink(fdrSettings.getMinLinkPepCount());
        getFdr().setMinPepPerProteinGroupPair(fdrSettings.getMinPPIPepCount());

        int lengthSteps = tblPepLength.getRowCount();
        int[] pepLength = new int[lengthSteps];
        for (int i = 0; i < lengthSteps - 1; i++) {
            pepLength[i] = Integer.parseInt(tblPepLength.getValueAt(i, 0).toString());
        }

        getFdr().setLengthGroups(pepLength);
    }

    protected String fdrLevelSummary(FDRResultLevel fdrl) {
        StringBuilder sbSummary = new StringBuilder();
        for (Object fg: fdrl.keySet()) {
            SubGroupFdrInfo sg = (SubGroupFdrInfo) fdrl.get(fg);
            String gn = sg.fdrGroupName;
            sbSummary.append(gn);
            sbSummary.append(" : ");
            sbSummary.append(String.format("%6d",sg.results.size()));
            sbSummary.append("[");
            sbSummary.append(String.format("%6.2f",sg.higherFDR*100));
            sbSummary.append(",");
            sbSummary.append(String.format("%6.2f",sg.lowerFDR*100));
            sbSummary.append("]");
            sbSummary.append(sg.filteredResult.size());
            sbSummary.append(" \n");
        }
        return sbSummary.toString();
    }

    /**
     * @return the m_fdr
     */
    public OfflineFDR getFdr() {
        return m_fdr;
    }

    /**
     * @param m_fdr the m_fdr to set
     */
    public void setFdr(OfflineFDR m_fdr) {
        this.m_fdr = m_fdr;
    }

    private class NeededOptionalComboBoxCellEditor extends DefaultCellEditor {

        EditorDelegate neededDelegate;
        EditorDelegate optionalDelegate;
        JComponent neededEditorComponent;
        JComponent optionalEditorComponent;

        public NeededOptionalComboBoxCellEditor(final JComboBox needed, final JComboBox optional) {
            super(needed);
            this.neededDelegate = delegate;
            this.neededEditorComponent = needed;

            optionalDelegate = new EditorDelegate() {
                public void setValue(Object value) {
                    optional.setSelectedItem(value);
                }

                public Object getCellEditorValue() {
                    return optional.getSelectedItem();
                }

                public boolean shouldSelectCell(EventObject anEvent) {
                    if (anEvent instanceof MouseEvent) {
                        MouseEvent e = (MouseEvent) anEvent;
                        return e.getID() != MouseEvent.MOUSE_DRAGGED;
                    }
                    return true;
                }

                public boolean stopCellEditing() {
                    if (optional.isEditable()) {
                        // Commit edited value.
                        optional.actionPerformed(new ActionEvent(NeededOptionalComboBoxCellEditor.this, 0, ""));
                    }
                    return super.stopCellEditing();
                }
            };
            optional.addActionListener(optionalDelegate);
            optionalEditorComponent = optional;

        }

        @Override
        public Component getTableCellEditorComponent(JTable table, Object value,
                boolean isSelected,
                int row, int column) {
            TableModel tm = table.getModel();
            if (Boolean.TRUE.equals(tm.getValueAt(row, 1))) {
                delegate = optionalDelegate;
                editorComponent = optionalEditorComponent;
            } else {
                delegate = neededDelegate;
                editorComponent = neededEditorComponent;
            }
            return super.getTableCellEditorComponent(table, value, isSelected, row, column);
        }
    }

    protected void writeMZIdentML() {
        try {
            File f = fbMzIdentMLOut.getFile();
            // if (f.canWrite())
            if ( getFdr() instanceof MZIdentXLFDR) {
                ((MZIdentXLFDR) getFdr()).writeMZIdentML(f.getAbsolutePath(), m_result);
            } else {
                JOptionPane.showMessageDialog(this, "mzIdentML export currently only supported with mzIdentML source files");
            }
        } catch (Exception ex) {
            JOptionPane.showMessageDialog(this, "Error while writing the xml-file:" + ex + "\n" + MyArrayUtils.toString(ex.getStackTrace(), "\n"));
            Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, "Error while writing the xml-file", ex);
            setStatus("error:" + ex);
            
            
        }
    }

    private boolean mzScoreDirectionChanged = false;

    protected void writeResults() {
        final OfflineFDR ofdr = getFdr();

        setEnableRead(false);
        setEnableCalc(false);
        setEnableWrite(false);
        final boolean prepostaa = ckPrePostAA.isSelected();
        final String sep = rbTSV.isSelected() ? "\t" : ",";

        Runnable runnable = new Runnable() {
            public void run() {
                try {
                    Logger.getLogger(this.getClass().getName()).log(Level.INFO, "start writing");
                    //                    ofdr.writeFiles(txtFolder.getText(), txtBaseName.getText(), 0.05, 0.05, 0.05, 0.05, new int[]{4, 9});
                    ofdr.writeFiles(fbFolder.getText(), txtBaseName.getText(), sep, m_result);
                    setStatus("finished: " + ofdr.summaryString(m_result));

                } catch (FileNotFoundException ex) {
                    Logger.getLogger(FDRGUI.class.getName()).log(Level.SEVERE, null, ex);
                    setStatus("error: " + ex);
                } catch (Exception e) {
                    Logger.getLogger(FDRGUI.class.getName()).log(Level.SEVERE, null, e);
                    setStatus("error: " + e);
                }
                setEnableRead(true);
                setEnableCalc(true);
                setEnableWrite(true);
            }
        };
        new Thread(runnable).start();
    }


    public void readCSV() {
        try {
            setEnableRead(false);
            setEnableCalc(false);
            setEnableWrite(false);

            setStatus("Start");
            final CSVinFDR ofdr = new CSVinFDR();



            setFdr(ofdr);

            getFdr().setPSMScoreHighBetter(csvSelect.higherIsBetter());


            Runnable runnable = new Runnable() {
                public void run() {
                    try {
                        setStatus("Read from " + csvSelect.getFile().getName());
                        CsvParser csv = csvSelect.getCsv();
                        Logger.getLogger(this.getClass().getName()).log(Level.INFO, "read from " + csvSelect.getFile().getAbsolutePath());
                        ofdr.readCSV(csv);
                        setEnableCalc(true);
                        setStatus("finished reading");
                    } catch (Exception ex) {
                        Logger.getLogger(FDRGUI.class.getName()).log(Level.SEVERE, null, ex);
                        setStatus("error:" + ex);
                    }
                    setEnableRead(true);
                    //  setEnableWrite(true);
                }


            };
            new Thread(runnable).start();
        } catch (Exception ex) {
            Logger.getLogger(FDRGUI.class.getName()).log(Level.SEVERE, null, ex);
            setStatus("error:" + ex);
            setEnableRead(true);
        }

    }

    public void readMZIdentML() {
        try {
            setEnableRead(false);
            setEnableCalc(false);
            setEnableWrite(false);

            final File mzIdentIn = fbMZIdentMLIn.getFile();
            if (mzIdentIn == null) {
                setStatus("no file selected");
                setEnableRead(true);
                //btnRead.setEnabled(true);
                return;
            } else if (!mzIdentIn.canRead()) {
                setStatus("Can't read selected file");
                setEnableRead(true);
                return;
            }

            setStatus("Start");
            final MZIdentXLFDR ofdr = new MZIdentXLFDR();


            setFdr(ofdr);

            getFdr().setPSMScoreHighBetter(rbMZHighBetter.isSelected());
//            ofdr.setCrosslinkedDonorModAcc(txtCrossLinkedDonorModCvParam.getText());
//            ofdr.setCrosslinkedReceptorModAcc(txtCrossLinkedReceptorModCvParam.getText());
//            ofdr.setCrosslinkedSIIAcc(txtCrossLinkedPepCvParam.getText());
            ofdr.setPSMScore(cbMZMatchScoreName.getSelectedItem().toString());


            Runnable runnable = new Runnable() {
                public void run() {
                    try {
                        setStatus("Read from " + mzIdentIn.getName());
                        Logger.getLogger(this.getClass().getName()).log(Level.INFO, "read from " + mzIdentIn.getAbsolutePath());
                        ofdr.readMzIdentML(mzIdentIn);
                        setEnableCalc(true);
                        setStatus("finished reading");
                    } catch (Exception ex) {
                        Logger.getLogger(FDRGUI.class.getName()).log(Level.SEVERE, null, ex);
                        setStatus("error:" + ex);
                    }
                    setEnableRead(true);
                    //  setEnableWrite(true);
                }
            };
            new Thread(runnable).start();
        } catch (Exception ex) {
            Logger.getLogger(FDRGUI.class.getName()).log(Level.SEVERE, null, ex);
            setStatus("error:" + ex);
            setEnableRead(true);
        }

    }

    protected void calculateFDR() {
        final OfflineFDR ofdr = getFdr();
        setEnableCalc(false);
        setEnableWrite(false);
        setEnableRead(false);
        final Double saftyfactor = (Double) fdrSettings.getReportFactor();

        prepareFDRCalculation();


        Runnable runnable = new Runnable() {
            final Double psmfdr = (Double) fdrSettings.getPSMFDR() ;
            final Double pepfdr = (Double) fdrSettings.getPeptidePairFDR() ;
            final Double protfdr = (Double) fdrSettings.getProteinGroupFDR() ;
            final Double linkfdr = (Double) fdrSettings.getProteinGroupLinkFDR();
            final Double ppifdr = (Double) fdrSettings.getProteinGroupPairFDR();
            final boolean filterToUniquePSM = fdrSettings.filterToUniquePSM();

            public void run() {
                try {
                    innerFDRCalculation(psmfdr, pepfdr, protfdr, linkfdr, ppifdr, ofdr, saftyfactor,filterToUniquePSM);

                } catch (Exception e) {
                    Logger.getLogger(FDRGUI.class.getName()).log(Level.SEVERE, null, e);
                    setStatus("error:" + e);
                }
                setEnableCalc(true);
                setEnableRead(true);
                //setEnableWrite(true);
            }
        };
        new Thread(runnable).start();
    }

    protected void innerFDRCalculation(Double psmfdr, Double pepfdr, Double protfdr, Double linkfdr, Double ppifdr, OfflineFDR ofdr, Double saftyfactor, boolean filterToUniquePSM) {
        setStatus("Start");
        setStatus("Calculating fdr");
        Logger.getLogger(this.getClass().getName()).log(Level.INFO, "Calculating fdr");
        FDRResult result = ofdr.calculateFDR(psmfdr, pepfdr, protfdr, linkfdr, ppifdr, saftyfactor, ckIgnoreGroups1.isSelected(), true, filterToUniquePSM);
        m_result = result;
        setStatus("finished");
        setStatus(ofdr.summaryString(m_result));
        setEnableWrite(true);
//        HashMap<String, Double> summary = ofdr.summaryString(m_result);
        txtSumInput.setText(Integer.toString(ofdr.getInputPSMs().size()));


        int sumXLPSM = 0;
        int sumLinearPSM = 0;

        int targetXLPSM = 0;
        int targetLinearPSM = 0;
        int decoyXLPSM = 0;
        int decoyLinearPSM = 0;
        for (PSM psm : m_result.psmFDR) {
            if (psm.isLinear()) {
                sumLinearPSM++;
                if (!psm.isDecoy()) {
                    targetLinearPSM++;
                }

            } else {
                sumXLPSM++;
                if (!psm.isDecoy()) {
                    targetXLPSM++;
                }
            }
        }
        int sumPSM = sumXLPSM + sumLinearPSM;


        int sumXLPepPairs = 0;
        int sumLinearPepPairs = 0;

        int targetXLPepPairs = 0;
        int targetLinearPepPairs = 0;
        int decoyXLPepPairs = 0;
        int decoyLinearPepPairs = 0;
        for (PeptidePair pp : result.peptidePairFDR) {
            if (pp.isLinear()) {
                sumLinearPepPairs++;
                if (!pp.isDecoy()) {
                    targetLinearPepPairs++;
                }

            } else {
                sumXLPepPairs++;
                if (!pp.isDecoy()) {
                    targetXLPepPairs++;
                }
            }
        }
        int sumPepPairs = sumXLPepPairs + sumLinearPepPairs;



        int sumLinksBetweenDecoy = 0;
        int sumLinksInternalDecoy = 0;
        int sumLinksInternalTarget = 0;
        int sumLinksBetweenTarget = 0;
        int sumLinks = m_result.proteinGroupLinkFDR.getResultCount();
        for (ProteinGroupLink l : m_result.proteinGroupLinkFDR) {
            if (l.isDecoy()) {
                if (l.isInternal) {
                    sumLinksInternalDecoy++;
                } else {
                    sumLinksBetweenDecoy++;
                }
            } else if (l.isInternal) {
                sumLinksInternalTarget++;
            } else {
                sumLinksBetweenTarget++;
            }

        }


        int sumProteinGroupPairs = m_result.proteinGroupPairFDR.getResultCount();
        int sumProteinGroupPairsBetweenDecoy = 0;
        int sumProteinGroupPairsInternalDecoy = 0;
        int sumProteinGroupPairsInternalTarget = 0;
        int sumProteinGroupPairsBetweenTarget = 0;
        for (ProteinGroupPair pgl : m_result.proteinGroupPairFDR) {
            if (pgl.isDecoy()) {
                if (pgl.isInternal()) {
                    sumProteinGroupPairsInternalDecoy++;
                } else {
                    sumProteinGroupPairsBetweenDecoy++;
                }
            } else if (pgl.isInternal()) {
                sumProteinGroupPairsInternalTarget++;
            } else {
                sumProteinGroupPairsBetweenTarget++;
            }

        }
//        Integer sumLinksProtGroups =  ofdr.getFDRProteinGroups().size();
        
        String[] nice =  MiscUtils.arrayToStringWithDifferenceOrientedFormat(new double[] {result.psmFDR.getHigherFDR()*100, result.psmFDR.getLowerFDR()*100},1);
        txtSumPSM.setText(sumPSM + " [" + nice[0] + "%," + nice[1] + "%]");
        txtSumPSM.setToolTipText(fdrLevelSummary(result.psmFDR));
        txtSumPSMXL.setText(sumXLPSM + " (" + (targetXLPSM) + " Target)");
        txtSumPSMLinear.setText(sumLinearPSM + " (" + (targetLinearPSM) + " Target)");
        
        nice =  MiscUtils.arrayToStringWithDifferenceOrientedFormat(new double[] {result.peptidePairFDR.getHigherFDR()*100, result.peptidePairFDR.getLowerFDR()*100},1);
        txtSumPepPairs.setText(sumPepPairs + " [" + nice[0] + "%," + nice[1] + "%]");
        txtSumPepPairs.setToolTipText(fdrLevelSummary(result.peptidePairFDR));
        txtSumPepPairsXL.setText(sumXLPepPairs + " (" + targetXLPepPairs + " Target)");
        txtSumPepPairsLinear.setText(sumLinearPepPairs + " (" + targetLinearPepPairs + " Target)");
        
        nice =  MiscUtils.arrayToStringWithDifferenceOrientedFormat(new double[] {result.proteinGroupLinkFDR.getHigherFDR()*100, result.proteinGroupLinkFDR.getLowerFDR()*100},1);
        txtSumLinks.setText(sumLinks + " [" +  nice[0] + "%," + nice[1] + "%]");
        txtSumLinks.setToolTipText(fdrLevelSummary(result.proteinGroupLinkFDR));
        txtSumLinksBetween.setText((sumLinksBetweenTarget + sumLinksBetweenDecoy) + " (" + sumLinksBetweenTarget + " Target)");
        txtSumLinksInternal.setText((sumLinksInternalDecoy + sumLinksInternalTarget) + " (" + sumLinksInternalTarget + ")");
        
        nice =  MiscUtils.arrayToStringWithDifferenceOrientedFormat(new double[] {result.proteinGroupFDR.getHigherFDR()*100, result.proteinGroupFDR.getLowerFDR()*100},1);
        txtSumProtGroups.setText(Integer.toString(m_result.proteinGroupFDR.getResultCount()) + " [" +  nice[0] + "%," + nice[1] + "%]");
        txtSumProtGroups.setToolTipText(fdrLevelSummary(result.proteinGroupFDR));
        
        nice =  MiscUtils.arrayToStringWithDifferenceOrientedFormat(new double[] {result.proteinGroupPairFDR.getHigherFDR()*100, result.proteinGroupPairFDR.getLowerFDR()*100},1);
        txtSumProtGroupPairs.setText(sumProteinGroupPairs + " [" +  nice[0] + "%," + nice[1] + "%]");
        txtSumProtGroupPairs.setToolTipText(fdrLevelSummary(result.proteinGroupPairFDR));
        txtSumProtGroupPairsBetween.setText((sumProteinGroupPairsBetweenTarget + sumProteinGroupPairsBetweenDecoy) + " (" + sumProteinGroupPairsBetweenTarget + " Target)");
        txtSumProtGroupPairsInternal.setText((sumProteinGroupPairsInternalDecoy + sumProteinGroupPairsInternalTarget) + " (" + sumProteinGroupPairsInternalTarget + " Target)");

        // Logger.getLogger(this.getClass().getName()).log(Level.INFO, "finished writing");
    }


    

    
    private void basicCalc() {
        if (getFdr() !=null) {
            fdrSettingsComplete.setPSMDirectional(fdrSettings.isPSMDirectional());
            fdrSettingsComplete.setPeptidePairDirectional(fdrSettings.isPeptidePairDirectional());
            fdrSettingsComplete.setLinkDirectional(fdrSettings.isLinkDirectional());
            fdrSettingsComplete.setPPIDirectional(fdrSettings.isPPIDirectional());

            fdrSettingsComplete.setPSMFDR(fdrSettings.getPSMFDR());
            fdrSettingsComplete.setPeptidePairFDR(fdrSettings.getPeptidePairFDR());
            fdrSettingsComplete.setProteinGroupFDR(fdrSettings.getProteinGroupFDR());
            fdrSettingsComplete.setProteinGroupLinkFDR(fdrSettings.getProteinGroupLinkFDR());
            fdrSettingsComplete.setProteinGroupPairFDR(fdrSettings.getProteinGroupPairFDR());

            fdrSettingsComplete.setMaxLinkAmbiguity(fdrSettings.getMaxLinkAmbiguity());
            fdrSettingsComplete.setMaxProteinAmbiguity(fdrSettings.getMaxProteinAmbiguity());
            fdrSettingsComplete.setMinLinkPepCount(fdrSettings.getMinLinkPepCount());
            fdrSettingsComplete.setMinPPIPepCount(fdrSettings.getMinPPIPepCount());
            fdrSettingsComplete.setMinPeptideLength(fdrSettings.getMinPeptideLength());
            fdrSettingsComplete.setMinProteinPepCount(fdrSettings.getMinProteinPepCount());

            fdrSettingsSimple.setPSMDirectional(fdrSettings.isPSMDirectional());
            fdrSettingsSimple.setPeptidePairDirectional(fdrSettings.isPeptidePairDirectional());
            fdrSettingsSimple.setLinkDirectional(fdrSettings.isLinkDirectional());
            fdrSettingsSimple.setPPIDirectional(fdrSettings.isPPIDirectional());

            fdrSettingsSimple.setPSMFDR(fdrSettings.getPSMFDR());
            fdrSettingsSimple.setPeptidePairFDR(fdrSettings.getPeptidePairFDR());
            fdrSettingsSimple.setProteinGroupFDR(fdrSettings.getProteinGroupFDR());
            fdrSettingsSimple.setProteinGroupLinkFDR(fdrSettings.getProteinGroupLinkFDR());
            fdrSettingsSimple.setProteinGroupPairFDR(fdrSettings.getProteinGroupPairFDR());

            fdrSettingsSimple.setMaxLinkAmbiguity(fdrSettings.getMaxLinkAmbiguity());
            fdrSettingsSimple.setMaxProteinAmbiguity(fdrSettings.getMaxProteinAmbiguity());
            fdrSettingsSimple.setMinLinkPepCount(fdrSettings.getMinLinkPepCount());
            fdrSettingsSimple.setMinPPIPepCount(fdrSettings.getMinPPIPepCount());
            fdrSettingsSimple.setMinPeptideLength(fdrSettings.getMinPeptideLength());
            fdrSettingsSimple.setMinProteinPepCount(fdrSettings.getMinProteinPepCount());


            if (fdrSettings.doOptimize() != null) {
                switch (fdrSettings.doOptimize()) {
                    case PROTEINGROUPLINK:
                        Thread ml = new Thread() {
                            public void run() {
                                maximiseLink();
                            }
                        };
                        ml.start();
                        return;

                    case PROTEINGROUPPAIR:
                        Thread mp = new Thread() {
                            public void run() {
                                maximisePPI();
                            }
                        };
                        mp.start();
                        return;
                }
                JOptionPane.showMessageDialog(this, "Boosting of that Level currently not supported!");
            } else {
                calculateFDR();
            }
        }
    }


    public void setStatus(final String status) {
        m_status.write(status);
//        m_statusmessages.put(Thread.currentThread(),status);
//        Runnable setModel = new Runnable() {
//            public void run() {
//                String message = "";
//                for (Thread t : m_statusmessages.keySet())  {
//                    message+=" - " + m_statusmessages.get(t);
//                    if (!t.isAlive()) {
//                        m_statusmessages.remove(t);
//                    }
//                }
//                txtStatus.setText(message.substring(3));
//                set
//                
//            }
//        };
//        javax.swing.SwingUtilities.invokeLater(setModel);
    }

    public void setEnableRead(final boolean enable) {
        Runnable setModel = new Runnable() {
            public void run() {
                csvSelect.setEnabled(enable);
                btnReadMZIdent.setEnabled(enable && ((JPanel) fdrSettings).isEnabled());

            }
        };
        javax.swing.SwingUtilities.invokeLater(setModel);
    }
    
    public void setEnableCalc(final boolean enable) {
        Runnable setModel = new Runnable() {
            public void run() {
                ((JPanel) fdrSettings).setEnabled(enable);
//                btnMaxLink.setEnabled(enable);
//                btnMaxPPI.setEnabled(enable);
                fdrSettingsComplete.setEnabled(enable);
                fdrSettingsSimple.setEnabled(enable);
            }
        };
        javax.swing.SwingUtilities.invokeLater(setModel);
    }

    protected void setEnableWrite(final boolean enable) {
        final boolean isMzIdent = getFdr() instanceof MZIdentXLFDR;
        Runnable setModel = new Runnable() {
            public void run() {
                btnWrite.setEnabled(enable);

                btnWriteMzIdentML.setEnabled(enable && (isMzIdent));
            }
        };
        javax.swing.SwingUtilities.invokeLater(setModel);
    }

    public void clearResults() {
        txtSumLinks.setText("");
        txtSumLinksBetween.setText("");
        txtSumLinksInternal.setText("");
        txtSumPSM.setText("");
        txtSumPSMLinear.setText("");
        txtSumPSMXL.setText("");
        txtSumPepPairs.setText("");
        txtSumPepPairsLinear.setText("");
        txtSumPepPairsXL.setText("");
        txtSumProtGroupPairs.setText("");
        txtSumProtGroupPairsBetween.setText("");
        txtSumProtGroupPairsInternal.setText("");
        txtSumProtGroups.setText("");
    }

//    public void maximisePPIold() {
//        setEnableRead(false);
//        setEnableCalc(false);
//        setEnableWrite(false);
//
//        double maxpsmfdr = (Double) spPsmFDR.getValue() / 100.0;
//        double maxpepfdr = (Double) spPepFDR.getValue() / 100.0;
//        double maxlinkfdr = (Double) spLinkFDR.getValue() / 100.0;
//        double maxprotfdr = (Double) spProteinFDR.getValue() / 100.0;
//        double maxppifdr = (Double) spPPIFdr.getValue() / 100.0;
//        double reportFactor = (Double) spReportFactor.getValue();
//
//        double minpsmfdr = Math.min(0.05, maxpsmfdr / 4);
//        double minpepfdr = Math.min(0.05, maxpepfdr / 4);;
//        double minlinkfdr = Math.min(0.05, maxlinkfdr / 4);;
//        double minprotfdr = Math.min(0.05, maxprotfdr / 4);;
//
//        double psmfdr = 1;
//        double pepfdr = 1;
//        double linkfdr = 1;
//        double protfdr = 1;
//
//        double psmfdrStep = (maxpsmfdr - minpsmfdr) / 4;
//        double pepfdrStep = (maxpepfdr - minpepfdr) / 4;
//        double linkfdrStep = (maxlinkfdr - minlinkfdr) / 4;
//        double protfdrStep = (maxprotfdr - minprotfdr) / 4;
//
//        boolean optimizing = true;
//
//        int maxPPICount = 0;
//
//
//        double maxPPIPsmfdr = 0;
//        double maxPPIPepfdr = 0;
//        double maxPPILinkfdr = 0;
//        double maxPPIProtfdr = 0;
//
//        double maxPPIminPsmfdr = 0;
//        double maxPPIminPepfdr = 0;
//        double maxPPIminLinkfdr = 0;
//        double maxPPIminProtfdr = 0;
//        double maxPPImaxPsmfdr = 0;
//        double maxPPImaxPepfdr = 0;
//        double maxPPImaxLinkfdr = 0;
//        double maxPPImaxProtfdr = 0;
//        StringBuffer sb = new StringBuffer();
//
//        int countDown = 5;
//
//        prepareFDRCalculation();
//
//        while (optimizing) {
//            int lastMaxPPICount = maxPPICount;
//
//            // find the combinations with the maximum number of ppis
//            psmloop:
//            for (psmfdr = maxpsmfdr; psmfdr > 0 && psmfdr >= minpsmfdr - psmfdrStep / 2; psmfdr -= psmfdrStep) {
//
//                for (pepfdr = maxpepfdr; pepfdr > 0 && pepfdr >= minpepfdr - pepfdrStep / 2; pepfdr -= pepfdrStep) {
//
//                    for (linkfdr = maxlinkfdr; linkfdr > 0 && linkfdr >= minlinkfdr - linkfdrStep / 2; linkfdr -= linkfdrStep) {
//
//                        for (protfdr = maxprotfdr; protfdr > 0 & protfdr >= minprotfdr - protfdrStep / 2; protfdr -= protfdrStep) {
//
//                            innerFDRCalculation(psmfdr, pepfdr, protfdr, linkfdr, maxppifdr, m_fdr, reportFactor);
//
//                            int ppiCount = m_result.proteinGroupPairFDR.getResultCount();
//
//                            if (ppiCount > maxPPICount) {
//
//                                maxPPIPsmfdr = psmfdr;
//                                maxPPIPepfdr = pepfdr;
//                                maxPPILinkfdr = linkfdr;
//                                maxPPIProtfdr = protfdr;
//
//                                maxPPIminPsmfdr = psmfdr;
//                                maxPPIminPepfdr = pepfdr;
//                                maxPPIminLinkfdr = linkfdr;
//                                maxPPIminProtfdr = protfdr;
//                                maxPPImaxPsmfdr = psmfdr;
//                                maxPPImaxPepfdr = pepfdr;
//                                maxPPImaxLinkfdr = linkfdr;
//                                maxPPImaxProtfdr = protfdr;
//
//                                maxPPICount = ppiCount;
//                                String message = "psmfdr, " + psmfdr + " ,pepfdr, " + pepfdr + " ,protfdr, " + protfdr + ", linkfdr, " + linkfdr + ", ppi count, " + ppiCount;
//                                sb.append(message + "\n");
//                                Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, message);
//
//                            } else if (ppiCount == maxPPICount) {
//                                maxPPIminPsmfdr = Math.min(psmfdr, maxPPIminPsmfdr);
//                                maxPPIminPepfdr = Math.min(pepfdr, maxPPIminPepfdr);
//                                maxPPIminLinkfdr = Math.min(linkfdr, maxPPIminLinkfdr);
//                                maxPPIminProtfdr = Math.min(protfdr, maxPPIminProtfdr);
//                                maxPPImaxPsmfdr = Math.max(psmfdr, maxPPIminPsmfdr);
//                                maxPPImaxPepfdr = Math.max(pepfdr, maxPPIminPepfdr);;
//                                maxPPImaxLinkfdr = Math.max(linkfdr, maxPPIminLinkfdr);
//                                maxPPImaxProtfdr = Math.max(protfdr, maxPPIminProtfdr);;
//                            }
//
//                            if (stopMaximizing) {
//                                break psmloop;
//                            }
//
//                            if (m_result.proteinGroupFDR.getResultCount() == 0) {
//                                break;
//                            }
//
//                        }
//
//                    }
//
//                }
//
//            }
//
//            // no improvement for the last few rounds?
//            if ((maxPPICount == lastMaxPPICount && --countDown == 0) || stopMaximizing) {
//                optimizing = false;
//                Logger.getLogger(this.getClass().getName()).log(Level.INFO, sb.toString());
//                final int foundPPICount = maxPPICount;
//                final Double setPsmfdr = maxPPIPsmfdr;
//                final Double setPepfdr = maxPPIPepfdr;
//                final Double setLinkfdr = maxPPILinkfdr;
//                final Double setProtfdr = maxPPIProtfdr;
//                boolean isSet = false;
//                int trySet = 10;
//                while (!isSet && trySet > 0) {
//                    try {
//                        javax.swing.SwingUtilities.invokeAndWait(new Runnable() {
//                            public void run() {
//                                spPsmFDR.setValue(setPsmfdr);
//                                spPepFDR.setValue(setPepfdr);
//                                spLinkFDR.setValue(setLinkfdr);
//                                spProteinFDR.setValue(setProtfdr);
//                            }
//                        });
//                        isSet = true;
//                    } catch (InterruptedException ex) {
//                        trySet--;
//                        Logger.getLogger(this.getClass().getName()).log(Level.SEVERE, null, ex);
//                    } catch (InvocationTargetException ex) {
//                        trySet--;
//                        Logger.getLogger(FDRGUI.class.getName()).log(Level.SEVERE, null, ex);
//                    }
//                }
//                if (trySet < 0) {
//                    JOptionPane.showMessageDialog(this, "Could not set the values for fdr-calculation\n"
//                            + "\nPSM fdr:    " + setPsmfdr
//                            + "\nPeptide fdr:" + setPepfdr
//                            + "\nLink fdr:   " + setLinkfdr
//                            + "\nProtein FDR:" + setProtfdr, "Failed to set Parameters", JOptionPane.ERROR_MESSAGE);
//                    setEnableRead(true);
//                    setEnableCalc(true);
//                    setEnableWrite(true);
//                    stopMaximizing = false;
//                    return;
//                }
//
//                javax.swing.SwingUtilities.invokeLater(new Runnable() {
//                    public void run() {
//                        JOptionPane.showMessageDialog(rootPane, "found " + foundPPICount + " protein-group-pairs for following settings \n"
//                                + "\nPSM fdr:    " + setPsmfdr
//                                + "\nPeptide fdr:" + setPepfdr
//                                + "\nLink fdr:   " + setLinkfdr
//                                + "\nProtein FDR:" + setProtfdr, "best parameters found for max protein group pairs", JOptionPane.INFORMATION_MESSAGE);
//                    }
//                });
//
//                innerFDRCalculation(setPsmfdr, setPepfdr, setProtfdr, setLinkfdr, maxppifdr, m_fdr, reportFactor);
//
//                stopMaximizing = false;
//
//                setEnableRead(true);
//                setEnableCalc(true);
//                setEnableWrite(true);
//
//            } else { // yes we improved
//                countDown = 5;
//
//                // so see if we make the resoltuion finer
//                // can we get a better result?
//                maxpsmfdr = Math.min(maxPPImaxPsmfdr + psmfdrStep / 2.0, 1);
//                maxpepfdr = Math.min(maxPPImaxPepfdr + pepfdrStep / 2.0, 1);
//                maxlinkfdr = Math.min(maxPPImaxLinkfdr + linkfdrStep / 2.0, 1);
//                maxprotfdr = Math.min(maxPPImaxProtfdr + protfdrStep / 2.0, 1);
//
//                minpsmfdr = Math.max(maxPPIminPsmfdr - psmfdrStep / 2.0, 0);
//                minpepfdr = Math.max(maxPPIminPepfdr - pepfdrStep / 2.0, 0);
//                minlinkfdr = Math.max(maxPPIminLinkfdr - linkfdrStep / 2.0, 0);
//                minprotfdr = Math.max(maxPPIminProtfdr - protfdrStep / 2.0, 0);
//
//
//                psmfdrStep = (maxpsmfdr - minpsmfdr) / 4;
//                pepfdrStep = (maxpepfdr - minpepfdr) / 4;
//                linkfdrStep = (maxlinkfdr - minlinkfdr) / 4;
//                protfdrStep = (maxprotfdr - minprotfdr) / 4;
//
//            }
//
//
//        }
//
//
//    }

    public void maximisePPI() {
        clearResults();
        double steps = fdrSettings.getBoostingSteps();
        setEnableRead(false);
        setEnableCalc(false);
        setEnableWrite(false);
        StringBuilder sb = new StringBuilder();

        // get some settings, that are constant for all calculations
        boolean ignoreGroups = ckIgnoreGroups1.isSelected();

        int minLinkPeptide = fdrSettings.getMinLinkPepCount();
        int minProteinPeptide = fdrSettings.getMinProteinPepCount();
        int minProteinPairPeptide = fdrSettings.getMinPPIPepCount();

        int maxLinkAmbiguity = fdrSettings.getMaxLinkAmbiguity();
        int maxProteinAmbiguity = fdrSettings.getMaxProteinAmbiguity();

        double reportFactor = fdrSettings.getReportFactor();
        
        boolean filterToUniquePSM = fdrSettings.filterToUniquePSM();


        // get the maximum fdr values, that we can play with
        double maxpsmfdr = fdrSettings.getPSMFDR() / 100.0;
        double maxpepfdr = fdrSettings.getPeptidePairFDR() / 100.0;
        double maxprotfdr = fdrSettings.getProteinGroupFDR() / 100.0;
        double maxlinkfdr = fdrSettings.getProteinGroupLinkFDR() / 100.0;

        // the fixed fdr values
        double maxppifdr = fdrSettings.getProteinGroupPairFDR() / 100.0;

        // in the first round we start looking on these values to the maximum values
        double minpsmfdr = Math.min(0.05, maxpsmfdr / steps);
        double minpepfdr = Math.min(0.05, maxpepfdr / steps);
        double minprotfdr = Math.min(0.05, maxprotfdr / steps);
        double minlinkfdr = Math.min(0.05, maxlinkfdr / steps);

        // and we run through with these steps
        double psmfdrStep = (maxpsmfdr - minpsmfdr) / steps;
        double pepfdrStep = (maxpepfdr - minpepfdr) / steps;
        double protfdrStep = (maxprotfdr - minprotfdr) / steps;
        double linkfdrStep = (maxlinkfdr - minlinkfdr) / steps;

        double psmfdr = 1;
        double pepfdr = 1;
        double protfdr = 1;
        double linkfdr = 1;


        boolean optimizing = true;

        int maxPPICount = 0;


        double maxPPIPsmfdr = 0;
        double maxPPIPepfdr = 0;
        double maxPPIProtfdr = 0;
        double maxPPILinkfdr = 0;

        double maxPPIminPsmfdr = 0;
        double maxPPIminPepfdr = 0;
        double maxPPIminProtfdr = 0;
        double maxPPIminLinkfdr = 0;
        double maxPPImaxPsmfdr = 0;
        double maxPPImaxPepfdr = 0;
        double maxPPImaxProtfdr = 0;
        double maxPPImaxLinkfdr = 0;

        int countDown = 5;

        prepareFDRCalculation();
        

        int optimizingRound=1;
        while (optimizing) {
            FDRResult result = new FDRResult();
            int lastMaxPPICount = maxPPICount;
            Logger.getLogger(this.getClass().getName()).log(Level.INFO, "Round {0}", optimizingRound++);
            
            Logger.getLogger(this.getClass().getName()).log(Level.INFO, "PSM fdr from  :        {0} to {1}\nPeptide pair fdr from  {2} to {3}\nProtein-groupfdr from  {4} to {5}\nLink fdr from          {6} to {7}\nSteps : {8}", new Object[]{minpsmfdr, maxpsmfdr, minpepfdr, maxpepfdr, minprotfdr, maxprotfdr, minlinkfdr, maxlinkfdr, steps});
            // find the combinations with the maximum number of ppis
            psmloop:
            for (psmfdr = maxpsmfdr; psmfdr > minpsmfdr - psmfdrStep / 2; psmfdr -= psmfdrStep) {
                // we only need to calculate the link fdr ones for each
                getFdr().calculatePSMFDR(psmfdr, reportFactor, ignoreGroups, false, result, getFdr().isPsm_directional(),filterToUniquePSM);
                Logger.getLogger(this.getClass().getName()).log(Level.INFO, "Round {0}", optimizingRound++);
                if (stopMaximizing) {
                    break psmloop;
                }      
                // if we don't get PSM - stop looking at later stages
                if (result.psmFDR.getResultCount() == 0) {
                    break psmloop;
                }

                peploop:
                for (pepfdr = maxpepfdr; pepfdr > minpepfdr - pepfdrStep / 2; pepfdr -= pepfdrStep) {

                    protloop:
                    for (protfdr = maxprotfdr; protfdr > minprotfdr - protfdrStep / 2; protfdr -= protfdrStep) {

                        // calculate peptide level fdr
                        getFdr().calculatePeptidePairFDR(pepfdr, reportFactor, ignoreGroups, false,result,getFdr().isPeptides_directional());
                        if (stopMaximizing) {
                            break psmloop;
                        }      
                        // if we don't get peptide pairs - stop looking at later stages
                        if (result.peptidePairFDR.getResultCount() == 0) {
                            break peploop;
                        }

                        
                        // calculate protein level fdr
                        if (protfdr < 1) {
                            getFdr().calculateProteinGroupFDR(protfdr, reportFactor, ignoreGroups, minProteinPeptide, maxProteinAmbiguity, false,result);
                            
                            if (result.proteinGroupFDR.getResultCount() == 0) {
                                break protloop;
                            }
                            // cut down the peptides by proteins
                            getFdr().filterFDRPeptidePairsByFDRProteinGroups(result);
                        }

                        if (stopMaximizing) {
                            break psmloop;
                        }                              
                        
                        linkloop:
                        for (linkfdr = maxlinkfdr; linkfdr > minlinkfdr - linkfdrStep / 2; linkfdr -= linkfdrStep) {

                            if (stopMaximizing) {
                                break psmloop;
                            }      
                            
                            // calculate links
                            getFdr().calculateLinkFDR(linkfdr, reportFactor, ignoreGroups, minLinkPeptide, maxLinkAmbiguity, false, result, getFdr().isLinks_directional());

                            if (result.proteinGroupLinkFDR.getResultCount() == 0) {
                                break linkloop;
                            }

                            if (stopMaximizing) {
                                break psmloop;
                            }      
                            getFdr().calculateProteinGroupPairFDR(maxppifdr, reportFactor, ignoreGroups, minProteinPairPeptide, 0, false, result, getFdr().isPpi_directional());
                            getFdr().filterFDRLinksByFDRProteinGroupPairs(result);

                            // how many links do we now have?
                            final int ppiCount = result.proteinGroupPairFDR.getResultCount();

                            // is it a new best?
                            if (ppiCount > maxPPICount) {

                                // store the values for this fdr
                                maxPPIPsmfdr = psmfdr;
                                maxPPIPepfdr = pepfdr;
                                maxPPIProtfdr = protfdr;
                                maxPPILinkfdr = linkfdr;

                                // and we found it just here
                                maxPPIminPsmfdr = psmfdr;
                                maxPPIminPepfdr = pepfdr;
                                maxPPIminProtfdr = protfdr;
                                maxPPIminLinkfdr = linkfdr;
                                maxPPImaxPsmfdr = psmfdr;
                                maxPPImaxPepfdr = pepfdr;
                                maxPPImaxProtfdr = protfdr;
                                maxPPImaxLinkfdr = linkfdr;

                                maxPPICount = ppiCount;
                                // record that we found a new top
                                String message = "psmfdr, " + psmfdr + " ,pepfdr, " + pepfdr + " ,protfdr, " + protfdr + " ,linkfdr, " + linkfdr + ", ppi count, " + ppiCount;
                                sb.append(message + "\n");
                                Logger.getLogger(this.getClass().getName()).log(Level.INFO, sb.toString());

                                // forward the values to the gui
                                final double schowPSMFDR = psmfdr * 100;
                                final double schowPepFDR = pepfdr * 100;
                                final double schowProtFDR = protfdr * 100;
                                final double schowLinkFDR = linkfdr * 100;

                                fdrSettingsComplete.setPSMFDR(schowPSMFDR);
                                fdrSettingsComplete.setPeptidePairFDR(schowPepFDR);
                                fdrSettingsComplete.setProteinGroupFDR(schowProtFDR);
                                fdrSettingsComplete.setProteinGroupLinkFDR(schowLinkFDR);
                                fdrSettingsSimple.setPSMFDR(schowPSMFDR);
                                fdrSettingsSimple.setPeptidePairFDR(schowPepFDR);
                                fdrSettingsSimple.setProteinGroupFDR(schowProtFDR);
                                fdrSettingsSimple.setProteinGroupLinkFDR(schowLinkFDR);
                                
                                javax.swing.SwingUtilities.invokeLater(new Runnable() {
                                    public void run() {
                                        txtSumProtGroupPairs.setText(ppiCount + "");
                                    }
                                });

                            } else if (ppiCount == maxPPICount) {
                                maxPPIminPsmfdr = Math.min(psmfdr, maxPPIminPsmfdr);
                                maxPPIminPepfdr = Math.min(pepfdr, maxPPIminPepfdr);
                                maxPPIminProtfdr = Math.min(protfdr, maxPPIminProtfdr);
                                maxPPImaxPsmfdr = Math.max(psmfdr, maxPPIminPsmfdr);
                                maxPPImaxPepfdr = Math.max(pepfdr, maxPPIminPepfdr);;
                                maxPPImaxProtfdr = Math.max(protfdr, maxPPIminProtfdr);;
                            }

                            if (stopMaximizing) {
                                break psmloop;
                            }

                        }
                    }


                }



            }
            
            setStatus("Max Protein-Pairs Round: " + optimizingRound + " - " + maxPPICount + " pairs");
            // no improvement for the last few rounds?
            if ((maxPPICount == lastMaxPPICount && --countDown == 0) || stopMaximizing) {
                optimizing = false;
                Logger.getLogger(this.getClass().getName()).log(Level.INFO, sb.toString());

                final int foundPPICount = maxPPICount;
                final Double setPsmfdr = maxPPIPsmfdr;
                final Double setPepfdr = maxPPIPepfdr;
                final Double setProtfdr = maxPPIProtfdr;
                final Double setLinkfdr = maxPPILinkfdr;
                boolean isSet = false;
                int trySet = 10;
                while (!isSet && trySet > 0) {
                    try {
                        javax.swing.SwingUtilities.invokeAndWait(new Runnable() {
                            public void run() {
                                fdrSettingsComplete.setPSMFDR(setPsmfdr);
                                fdrSettingsComplete.setPeptidePairFDR(setPepfdr);
                                fdrSettingsComplete.setProteinGroupFDR(setProtfdr);
                                fdrSettingsComplete.setProteinGroupLinkFDR(setLinkfdr);

                                fdrSettingsSimple.setPSMFDR(setPsmfdr);
                                fdrSettingsSimple.setPeptidePairFDR(setPepfdr);
                                fdrSettingsSimple.setProteinGroupFDR(setProtfdr);
                                fdrSettingsSimple.setProteinGroupLinkFDR(setLinkfdr);
                                
//                                spPsmFDR.setValue(setPsmfdr*100);
//                                spPepFDR.setValue(setPepfdr*100);
//                                spProteinFDR.setValue(setProtfdr*100);
//                                spLinkFDR.setValue(setLinkfdr*100);
                                txtSumProtGroupPairs.setText(foundPPICount + "");
                            }
                        });
                        isSet = true;
                    } catch (InterruptedException ex) {
                        trySet--;
                        Logger.getLogger(FDRGUI.class.getName()).log(Level.SEVERE, null, ex);
                    } catch (InvocationTargetException ex) {
                        trySet--;
                        Logger.getLogger(FDRGUI.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }
                if (trySet < 0) {
                    JOptionPane.showMessageDialog(this, "Could not set the values for fdr-calculation\n"
                            + "\nPSM fdr:    " + setPsmfdr*100 + "%"
                            + "\nPeptide fdr:" + setPepfdr*100 + "%"
                            + "\nProtein FDR:" + setProtfdr*100 + "%"
                            + "\nLink FDR:" + setLinkfdr *100 +"%", "Failed to set Parameters", JOptionPane.ERROR_MESSAGE);
                    setEnableRead(true);
                    setEnableCalc(true);
                    setEnableWrite(true);
                    stopMaximizing = false;
                    return;
                }

                javax.swing.SwingUtilities.invokeLater(new Runnable() {
                    public void run() {
                        JOptionPane.showMessageDialog(rootPane, "found " + foundPPICount + " protein group pairs for following settings \n"
                            + "\nPSM fdr:    " + setPsmfdr*100 + "%"
                            + "\nPeptide fdr:" + setPepfdr*100 + "%"
                            + "\nProtein FDR:" + setProtfdr*100 + "%"
                                + "\nLink FDR:" + setLinkfdr*100 +"%", "best parameters found for max protein group pairs", JOptionPane.INFORMATION_MESSAGE);
                    }
                });

                innerFDRCalculation(setPsmfdr, setPepfdr, setProtfdr, maxlinkfdr, maxppifdr, getFdr(), reportFactor, filterToUniquePSM);

                setEnableRead(true);
                setEnableCalc(true);
                setEnableWrite(true);

                stopMaximizing = false;

            } else {
                if (maxPPICount > lastMaxPPICount) {
                    // yes we improved
                    countDown = 5;
                } else {
                    setStatus("Max Protein-Pairs Round: " + optimizingRound + " - " + maxPPICount + " pairs - count down:" + countDown );                    
                }

                // so see if we make the resoltuion finer
                // can we get a better result?
                maxpsmfdr = Math.min(maxPPImaxPsmfdr + psmfdrStep / 2.0, 1);
                maxpepfdr = Math.min(maxPPImaxPepfdr + pepfdrStep / 2.0, 1);
                maxprotfdr = Math.min(maxPPImaxProtfdr + protfdrStep / 2.0, 1);

                minpsmfdr = Math.max(maxPPIminPsmfdr - psmfdrStep / 2.0, 0);
                minpepfdr = Math.max(maxPPIminPepfdr - pepfdrStep / 2.0, 0);
                minprotfdr = Math.max(maxPPIminProtfdr - protfdrStep / 2.0, 0);


                psmfdrStep = (maxpsmfdr - minpsmfdr) / steps;
                pepfdrStep = (maxpepfdr - minpepfdr) / steps;
                protfdrStep = (maxprotfdr - minprotfdr) / steps;

            }


        }


    }
//
//    public void maximise(int maxWhat, boolean fixedPSM, boolean fixedPep, boolean fixedProt, boolean fixedLink, boolean fixedPPI) {
//        clearResults();
//        double steps = (Integer)spMaximizeSteps.getValue();
//        setEnableRead(false);
//        setEnableCalc(false);
//        setEnableWrite(false);
//        StringBuilder sb = new StringBuilder();
//
//        // get some settings, that are constant for all calculations
//        boolean ignoreGroups = ckIgnoreGroups.isSelected();
//
//        int minLinkPeptide = (Integer) spMinLinkPepCount.getValue();
//        int minProteinPeptide = (Integer) spMinProteinPepCount.getValue();
//        int minProteinPairPeptide = (Integer) spMinPPIPepCount.getValue();
//
//        int maxLinkAmbiguity = (Integer) spMaxLinkAmbiguity.getValue();
//        int maxProteinAmbiguity = (Integer) spMaxProteinAmbiguity.getValue();
//
//        double reportFactor = (Double) spReportFactor.getValue();
//        
//        // define the maximum fdr to start from
//        double[] maxfdr = new double[] {
//            (Double) spPsmFDR.getValue() / 100.0,
//            (Double) spPepFDR.getValue() / 100.0,
//            (Double) spProteinFDR.getValue() / 100.0,
//            (Double) spLinkFDR.getValue() / 100.0,
//            (Double) spPPIFdr.getValue() / 100.0
//        };
//        // define the minimum fdr to start from
//        double[] minfdr = new double[maxfdr.length];
//        for (int i = 0; i< minfdr.length; i++) {
//            minfdr[i] = Math.min(0.05, maxfdr[i] / steps);
//        };
//        double[] fdrStep  = new double[maxfdr.length];
//        for (int i = 0; i< minfdr.length; i++) {
//            fdrStep[i] = (maxfdr[i]-minfdr[i]) / steps;
//        }
//
//        boolean optimizing = true;
//
//        int maxLinearCount = 0;
//        int maxWithinCount = 0;
//        int maxBetweenCount = 0;
//
//        for (int i = 0; i< minfdr.length; i++) {
//            fdrStep[i] = (maxfdr[i]-minfdr[i]) / steps;
//        }
//
//        for (int i = 0; i< minfdr.length; i++) {
//            fdrStep[i] = (maxfdr[i]-minfdr[i]) / steps;
//        }
//        
//        for (int i = 0; i< minfdr.length; i++) {
//            fdrStep[i] = (maxfdr[i]-minfdr[i]) / steps;
//        }
//        
//        
//        double maxPPIPsmfdr = 0;
//        double maxPPIPepfdr = 0;
//        double maxPPIProtfdr = 0;
//        double maxPPILinkfdr = 0;
//
//        double maxPPIminPsmfdr = 0;
//        double maxPPIminPepfdr = 0;
//        double maxPPIminProtfdr = 0;
//        double maxPPIminLinkfdr = 0;
//        double maxPPImaxPsmfdr = 0;
//        double maxPPImaxPepfdr = 0;
//        double maxPPImaxProtfdr = 0;
//        double maxPPImaxLinkfdr = 0;
//        
//
//        int countDown = 5;
//
//        prepareFDRCalculation();
//        
//
//        int optimizingRound=1;
//        while (optimizing) {
//            FDRResult result = new FDRResult();
//            int lastMaxPPICount = maxPPICount;
//            Logger.getLogger(this.getClass().getName()).log(Level.INFO, "Round {0}", optimizingRound++);
//            
//            Logger.getLogger(this.getClass().getName()).log(Level.INFO, "PSM fdr from  :        {0} to {1}\nPeptide pair fdr from  {2} to {3}\nProtein-groupfdr from  {4} to {5}\nLink fdr from          {6} to {7}\nSteps : {8}", new Object[]{minpsmfdr, maxpsmfdr, minpepfdr, maxpepfdr, minprotfdr, maxprotfdr, minlinkfdr, maxlinkfdr, steps});
//            // find the combinations with the maximum number of ppis
//            psmloop:
//            for (psmfdr = maxpsmfdr; psmfdr > minpsmfdr - psmfdrStep / 2; psmfdr -= psmfdrStep) {
//                // we only need to calculate the link fdr ones for each
//                m_fdr.calculatePSMFDR(psmfdr, reportFactor, ignoreGroups, false, result, m_fdr.isPsm_directional());
//                Logger.getLogger(this.getClass().getName()).log(Level.INFO, "Round {0}", optimizingRound++);
//                if (stopMaximizing) {
//                    break psmloop;
//                }      
//                // if we don't get PSM - stop looking at later stages
//                if (result.psmFDR.getResultCount() == 0) {
//                    break psmloop;
//                }
//                
//                if (maxWhat == 1) {
//                    
//                    m_fdr.calculatePeptidePairFDR(maxpepfdr, reportFactor, ignoreGroups, false,result,m_fdr.isPeptides_directional());
//                    if (protfdr < 1) {
//                        m_fdr.calculateProteinGroupFDR(protfdr, reportFactor, ignoreGroups, minProteinPeptide, maxProteinAmbiguity, false,result);
//
//                        if (result.proteinGroupFDR.getResultCount() == 0) {
//                            continue psmloop;
//                        }
//                        // cut down the peptides by proteins
//                        m_fdr.filterFDRPeptidePairsByFDRProteinGroups(result);
//                    }
//                    if (maxlinkfdr < 1 || maxppifdr <1) {
//                        m_fdr.calculateLinkFDR(maxlinkfdr, reportFactor, ignoreGroups, minLinkPeptide, maxLinkAmbiguity, false,result,m_fdr.isPeptides_directional());
//                        if (maxppifdr < 1) {
//                            m_fdr.calculateProteinGroupPairFDR(maxppifdr, reportFactor, ignoreGroups, minProteinPairPeptide, 0, false,result,m_fdr.isPeptides_directional());
//                            m_fdr.filterFDRLinksByFDRProteinGroupPairs(result);
//                        }
//                        m_fdr.filterFDRPeptidePairsByFDRProteinGroupLinks(result);
//                    }
//                    if (result.peptidePairFDR.getBetween() + result.peptidePairFDR.getWithin()){
//                        
//                    }
//                    
//                    
//                } else {
//                    peploop:
//                    for (pepfdr = maxpepfdr; pepfdr > minpepfdr - pepfdrStep / 2; pepfdr -= pepfdrStep) {
//
//                        protloop:
//                        for (protfdr = maxprotfdr; protfdr > minprotfdr - protfdrStep / 2; protfdr -= protfdrStep) {
//
//                            // calculate peptide level fdr
//                            m_fdr.calculatePeptidePairFDR(pepfdr, reportFactor, ignoreGroups, false,result,m_fdr.isPeptides_directional());
//                            if (stopMaximizing) {
//                                break psmloop;
//                            }      
//                            // if we don't get peptide pairs - stop looking at later stages
//                            if (result.peptidePairFDR.getResultCount() == 0) {
//                                break peploop;
//                            }
//
//
//                            // calculate protein level fdr
//                            if (protfdr < 1) {
//                                m_fdr.calculateProteinGroupFDR(protfdr, reportFactor, ignoreGroups, minProteinPeptide, maxProteinAmbiguity, false,result);
//
//                                if (result.proteinGroupFDR.getResultCount() == 0) {
//                                    break protloop;
//                                }
//                                // cut down the peptides by proteins
//                                m_fdr.filterFDRPeptidePairsByFDRProteinGroups(result);
//                            }
//
//                            if (stopMaximizing) {
//                                break psmloop;
//                            }                              
//
//                            linkloop:
//                            for (linkfdr = maxlinkfdr; linkfdr > minlinkfdr - linkfdrStep / 2; linkfdr -= linkfdrStep) {
//
//                                if (stopMaximizing) {
//                                    break psmloop;
//                                }      
//
//                                // calculate links
//                                m_fdr.calculateLinkFDR(linkfdr, reportFactor, ignoreGroups, minLinkPeptide, maxLinkAmbiguity, false, result, m_fdr.isLinks_directional());
//
//                                if (result.proteinGroupLinkFDR.getResultCount() == 0) {
//                                    break linkloop;
//                                }
//
//                                if (stopMaximizing) {
//                                    break psmloop;
//                                }      
//                                m_fdr.calculateProteinGroupPairFDR(maxppifdr, reportFactor, ignoreGroups, minProteinPairPeptide, 0, false, result, m_fdr.isPpi_directional());
//                                m_fdr.filterFDRLinksByFDRProteinGroupPairs(result);
//
//                                // how many links do we now have?
//                                final int ppiCount = result.proteinGroupPairFDR.getResultCount();
//
//                                // is it a new best?
//                                if (ppiCount > maxPPICount) {
//
//                                    // store the values for this fdr
//                                    maxPPIPsmfdr = psmfdr;
//                                    maxPPIPepfdr = pepfdr;
//                                    maxPPIProtfdr = protfdr;
//                                    maxPPILinkfdr = linkfdr;
//
//                                    // and we found it just here
//                                    maxPPIminPsmfdr = psmfdr;
//                                    maxPPIminPepfdr = pepfdr;
//                                    maxPPIminProtfdr = protfdr;
//                                    maxPPIminLinkfdr = linkfdr;
//                                    maxPPImaxPsmfdr = psmfdr;
//                                    maxPPImaxPepfdr = pepfdr;
//                                    maxPPImaxProtfdr = protfdr;
//                                    maxPPImaxLinkfdr = linkfdr;
//
//                                    maxPPICount = ppiCount;
//                                    // record that we found a new top
//                                    String message = "psmfdr, " + psmfdr + " ,pepfdr, " + pepfdr + " ,protfdr, " + protfdr + " ,linkfdr, " + linkfdr + ", ppi count, " + ppiCount;
//                                    sb.append(message + "\n");
//                                    Logger.getLogger(this.getClass().getName()).log(Level.INFO, sb.toString());
//
//                                    // forward the values to the gui
//                                    final double schowPSMFDR = psmfdr * 100;
//                                    final double schowPepFDR = pepfdr * 100;
//                                    final double schowProtFDR = protfdr * 100;
//                                    final double schowLinkFDR = linkfdr * 100;
//
//                                    javax.swing.SwingUtilities.invokeLater(new Runnable() {
//                                        public void run() {
//                                            spPsmFDR.setValue(schowPSMFDR);
//                                            spPepFDR.setValue(schowPepFDR);
//                                            spProteinFDR.setValue(schowProtFDR);
//                                            spLinkFDR.setValue(schowLinkFDR);
//                                            txtSumProtGroupPairs.setText(ppiCount + "");
//                                        }
//                                    });
//
//                                } else if (ppiCount == maxPPICount) {
//                                    maxPPIminPsmfdr = Math.min(psmfdr, maxPPIminPsmfdr);
//                                    maxPPIminPepfdr = Math.min(pepfdr, maxPPIminPepfdr);
//                                    maxPPIminProtfdr = Math.min(protfdr, maxPPIminProtfdr);
//                                    maxPPImaxPsmfdr = Math.max(psmfdr, maxPPIminPsmfdr);
//                                    maxPPImaxPepfdr = Math.max(pepfdr, maxPPIminPepfdr);;
//                                    maxPPImaxProtfdr = Math.max(protfdr, maxPPIminProtfdr);;
//                                }
//
//                                if (stopMaximizing) {
//                                    break psmloop;
//                                }
//
//                            }
//                        }
//
//
//                    }
//
//                }
//
//            }
//            
//            setStatus("Max Protein-Pairs Round: " + optimizingRound + " - " + maxPPICount + " pairs");
//            // no improvement for the last few rounds?
//            if ((maxPPICount == lastMaxPPICount && --countDown == 0) || stopMaximizing) {
//                optimizing = false;
//                Logger.getLogger(this.getClass().getName()).log(Level.INFO, sb.toString());
//
//                final int foundPPICount = maxPPICount;
//                final Double setPsmfdr = maxPPIPsmfdr;
//                final Double setPepfdr = maxPPIPepfdr;
//                final Double setProtfdr = maxPPIProtfdr;
//                final Double setLinkfdr = maxPPILinkfdr;
//                boolean isSet = false;
//                int trySet = 10;
//                while (!isSet && trySet > 0) {
//                    try {
//                        javax.swing.SwingUtilities.invokeAndWait(new Runnable() {
//                            public void run() {
//                                spPsmFDR.setValue(setPsmfdr*100);
//                                spPepFDR.setValue(setPepfdr*100);
//                                spProteinFDR.setValue(setProtfdr*100);
//                                spLinkFDR.setValue(setLinkfdr*100);
//                                txtSumProtGroupPairs.setText(foundPPICount + "");
//                            }
//                        });
//                        isSet = true;
//                    } catch (InterruptedException ex) {
//                        trySet--;
//                        Logger.getLogger(FDRGUI.class.getName()).log(Level.SEVERE, null, ex);
//                    } catch (InvocationTargetException ex) {
//                        trySet--;
//                        Logger.getLogger(FDRGUI.class.getName()).log(Level.SEVERE, null, ex);
//                    }
//                }
//                if (trySet < 0) {
//                    JOptionPane.showMessageDialog(this, "Could not set the values for fdr-calculation\n"
//                            + "\nPSM fdr:    " + setPsmfdr*100 + "%"
//                            + "\nPeptide fdr:" + setPepfdr*100 + "%"
//                            + "\nProtein FDR:" + setProtfdr*100 + "%"
//                            + "\nLink FDR:" + setLinkfdr *100 +"%", "Failed to set Parameters", JOptionPane.ERROR_MESSAGE);
//                    setEnableRead(true);
//                    setEnableCalc(true);
//                    setEnableWrite(true);
//                    stopMaximizing = false;
//                    return;
//                }
//
//                javax.swing.SwingUtilities.invokeLater(new Runnable() {
//                    public void run() {
//                        JOptionPane.showMessageDialog(rootPane, "found " + foundPPICount + " protein group pairs for following settings \n"
//                            + "\nPSM fdr:    " + setPsmfdr*100 + "%"
//                            + "\nPeptide fdr:" + setPepfdr*100 + "%"
//                            + "\nProtein FDR:" + setProtfdr*100 + "%"
//                                + "\nLink FDR:" + setLinkfdr*100 +"%", "best parameters found for max protein group pairs", JOptionPane.INFORMATION_MESSAGE);
//                    }
//                });
//
//                innerFDRCalculation(setPsmfdr, setPepfdr, setProtfdr, maxlinkfdr, maxppifdr, m_fdr, reportFactor);
//
//                setEnableRead(true);
//                setEnableCalc(true);
//                setEnableWrite(true);
//
//                stopMaximizing = false;
//
//            } else {
//                if (maxPPICount > lastMaxPPICount) {
//                    // yes we improved
//                    countDown = 5;
//                } else {
//                    setStatus("Max Protein-Pairs Round: " + optimizingRound + " - " + maxPPICount + " pairs - count down:" + countDown );                    
//                }
//
//                // so see if we make the resoltuion finer
//                // can we get a better result?
//                maxpsmfdr = Math.min(maxPPImaxPsmfdr + psmfdrStep / 2.0, 1);
//                maxpepfdr = Math.min(maxPPImaxPepfdr + pepfdrStep / 2.0, 1);
//                maxprotfdr = Math.min(maxPPImaxProtfdr + protfdrStep / 2.0, 1);
//
//                minpsmfdr = Math.max(maxPPIminPsmfdr - psmfdrStep / 2.0, 0);
//                minpepfdr = Math.max(maxPPIminPepfdr - pepfdrStep / 2.0, 0);
//                minprotfdr = Math.max(maxPPIminProtfdr - protfdrStep / 2.0, 0);
//
//
//                psmfdrStep = (maxpsmfdr - minpsmfdr) / steps;
//                pepfdrStep = (maxpepfdr - minpepfdr) / steps;
//                protfdrStep = (maxprotfdr - minprotfdr) / steps;
//
//            }
//
//
//        }
//
//
//    }
//    
    
    public void maximiseLink() {
        clearResults();
        
        double steps = fdrSettings.getBoostingSteps();
        setEnableRead(false);
        setEnableCalc(false);
        setEnableWrite(false);
        StringBuffer sb = new StringBuffer();

        // get some settings, that are constant for all calculations
        boolean ignoreGroups = ckIgnoreGroups1.isSelected();

        int minLinkPeptide = fdrSettings.getMinLinkPepCount();
        int minProteinPeptide = fdrSettings.getMinProteinPepCount();
        int minProteinPairPeptide = fdrSettings.getMinPPIPepCount();

        int maxLinkAmbiguity = fdrSettings.getMaxLinkAmbiguity();
        int maxProteinAmbiguity = fdrSettings.getMaxProteinAmbiguity();

        double reportFactor = fdrSettings.getReportFactor();


        // get the maximum fdr values, that we can play with
        double maxpsmfdr = fdrSettings.getPSMFDR();
        double maxpepfdr = fdrSettings.getPeptidePairFDR();
        double maxprotfdr = fdrSettings.getProteinGroupFDR();
        double maxlinkfdr = fdrSettings.getProteinGroupLinkFDR();

        // the fixed fdr values
        double maxppifdr = fdrSettings.getProteinGroupPairFDR();
        boolean filterToUniquePSM = fdrSettings.filterToUniquePSM();
//        int minLinkPeptide = (Integer) spMinLinkPepCount.getValue();
//        int minProteinPeptide = (Integer) spMinProteinPepCount.getValue();
//        int minProteinPairPeptide = (Integer) spMinPPIPepCount.getValue();
//
//        int maxLinkAmbiguity = (Integer) spMaxLinkAmbiguity.getValue();
//        int maxProteinAmbiguity = (Integer) spMaxProteinAmbiguity.getValue();
//
//        double reportFactor = (Double) spReportFactor.getValue();
//
//
//        // get the maximum fdr values, that we can play with
//        double maxpsmfdr = (Double) spPsmFDR.getValue() / 100.0;
//        double maxpepfdr = (Double) spPepFDR.getValue() / 100.0;
//        double maxprotfdr = (Double) spProteinFDR.getValue() / 100.0;

        double maxpsmfdrabs = maxpsmfdr;
        double maxpepfdrabs = maxpepfdr;
        double maxprotfdrabs = maxprotfdr;

        
//        // the fixed fdr values
//        double maxppifdr = (Double) spPPIFdr.getValue() / 100.0;
//        double maxlinkfdr = (Double) spLinkFDR.getValue() / 100.0;

        // in the first round we start looking on these values to the maximum values
        double minpsmfdr = Math.min(0.005, maxpsmfdr / steps);
        double minpepfdr = Math.min(0.005, maxpepfdr / steps);
        double minprotfdr = Math.min(0.005, maxprotfdr / steps);

        // and we run through with these steps
        double psmfdrStep = (maxpsmfdr - minpsmfdr) / steps;
        double pepfdrStep = (maxpepfdr - minpepfdr) / steps;
        double protfdrStep = (maxprotfdr - minprotfdr) / steps;

        double psmfdr = 1;
        double pepfdr = 1;
        double protfdr = 1;


        boolean optimizing = true;

        int maxLinkCount = 0;
        int maxPPICount = 0;


        double maxLinkPsmfdr = 0;
        double maxLinkPepfdr = 0;
        double maxLinkProtfdr = 0;

        double maxLinkminPsmfdr = 0;
        double maxLinkminPepfdr = 0;
        double maxLinkminProtfdr = 0;
        double maxLinkmaxPsmfdr = 0;
        double maxLinkmaxPepfdr = 0;
        double maxLinkmaxProtfdr = 0;

        int countDown = 5;

        prepareFDRCalculation();

        int optimizingRound=1;
        while (optimizing) {
            int lastMaxLinkCount = maxLinkCount;
            int lastMaxPPICount = maxPPICount;
            Logger.getLogger(this.getClass().getName()).log(Level.INFO,"Round " +optimizingRound++);
            
            Logger.getLogger(this.getClass().getName()).log(Level.INFO,"PSM fdr from  :        " + minpsmfdr +" to " + maxpsmfdr +
                                                                     "\nPeptide pair fdr from  " + minpepfdr +" to " + maxpepfdr +
                                                                     "\nProtein-groupfdr from  " + minprotfdr +" to " + maxprotfdr +            
                                                                     "\nSteps : " + steps);
            FDRResult result = new FDRResult();
            // find the combinations with the maximum number of ppis
            psmloop:
            for (psmfdr = maxpsmfdr; psmfdr > minpsmfdr - psmfdrStep / 2; psmfdr -= psmfdrStep) {
                // we only need to calculate the link fdr ones for each
                getFdr().calculatePSMFDR(psmfdr, reportFactor, ignoreGroups, false,result, getFdr().isPsm_directional(),filterToUniquePSM);
                if (stopMaximizing) {
                    break psmloop;
                }      
                // if we don't get PSM - stop looking at later stages
                if (result.psmFDR.getResultCount() == 0) {
                    break psmloop;
                }

                peploop:
                for (pepfdr = maxpepfdr; pepfdr > minpepfdr - pepfdrStep; pepfdr -= pepfdrStep) {

                    protloop:
                    for (protfdr = maxprotfdr; protfdr > minprotfdr - protfdrStep; protfdr -= protfdrStep) {
   
                        
                        // calculate peptide level fdr
                        getFdr().calculatePeptidePairFDR(pepfdr, reportFactor, ignoreGroups, false,result, getFdr().isPeptides_directional() );
                        // if we don't get peptide pairs - stop looking at later stages
                        if (result.peptidePairFDR.getResultCount() == 0) {
                            break peploop;
                        }

                        if (stopMaximizing) {
                            break psmloop;
                        }
                        
                        // calculate protein level fdr
                        if (protfdr < 1) {
                            getFdr().calculateProteinGroupFDR(protfdr, reportFactor, ignoreGroups, minProteinPeptide, maxProteinAmbiguity, false, result);

                            if (result.proteinGroupFDR.getResultCount() == 0) {
                                break protloop;
                            }
                            
                            // cut down the peptides by proteins                           
                            getFdr().filterFDRPeptidePairsByFDRProteinGroups(result);

                            if (stopMaximizing) {
                                break psmloop;
                            }                            
                        }

                        // calculate links
                        getFdr().calculateLinkFDR(maxlinkfdr, reportFactor, ignoreGroups, minLinkPeptide, maxLinkAmbiguity, false,result, getFdr().isLinks_directional());


                        
                        if (result.proteinGroupLinkFDR.getResultCount() == 0) {
                            break protloop;
                        }

                        // do we need to calculate protein pairs?
                        if (maxppifdr < 1) {
                            if (stopMaximizing) {
                                break psmloop;
                            }                                  
                            getFdr().calculateProteinGroupPairFDR(maxppifdr, reportFactor, ignoreGroups, minProteinPairPeptide, 0, false,result, getFdr().isPpi_directional());

                            if (result.proteinGroupPairFDR.getResultCount() == 0) {
                                break protloop;
                            }
                            
                            getFdr().filterFDRLinksByFDRProteinGroupPairs(result);

                  
                        }

                        
                        // how many links do we now have?
                        final int linkCount = result.proteinGroupLinkFDR.getResultCount();
                        final int ppiCount =  result.proteinGroupPairFDR ==null ? 0 : result.proteinGroupPairFDR.getResultCount();

                        // is it a new best?
                        if (linkCount > maxLinkCount || (linkCount == maxLinkCount && ppiCount > maxPPICount)) {
                            
                            // store the values for this fdr
                            maxLinkPsmfdr = psmfdr;
                            maxLinkPepfdr = pepfdr;
                            maxLinkProtfdr = protfdr;

                            // and we found it just here
                            maxLinkminPsmfdr = psmfdr;
                            maxLinkminPepfdr = pepfdr;
                            maxLinkminProtfdr = protfdr;
                            maxLinkmaxPsmfdr = psmfdr;
                            maxLinkmaxPepfdr = pepfdr;
                            maxLinkmaxProtfdr = protfdr;

                            maxLinkCount = linkCount;
                            maxPPICount = ppiCount;
                            // record that we found a new top
                            String message = "psmfdr, " + psmfdr + " ,pepfdr, " + pepfdr + " ,protfdr, " + protfdr + ", link count, " + linkCount;
                            if (ppiCount > 0 )
                                message += ", Protein Pairs, " +  ppiCount;
                            sb.append(message + "\n");
                            Logger.getLogger(this.getClass().getName()).log(Level.INFO, sb.toString());

                            // forward the values to the gui
                            final double schowPSMFDR = psmfdr;
                            final double schowPepFDR = pepfdr;
                            final double schowProtFDR = protfdr;

                            javax.swing.SwingUtilities.invokeLater(new Runnable() {
                                public void run() {
                                    fdrSettingsComplete.setPSMFDR(schowPSMFDR);
                                    fdrSettingsComplete.setPeptidePairFDR(schowPepFDR);
                                    fdrSettingsComplete.setProteinGroupFDR(schowProtFDR);

                                    fdrSettingsSimple.setPSMFDR(schowPSMFDR);
                                    fdrSettingsSimple.setPeptidePairFDR(schowPepFDR);
                                    fdrSettingsSimple.setProteinGroupFDR(schowProtFDR);
                                    
                                    txtSumLinks.setText(linkCount + "");
                                    
                                    if (ppiCount > 0)
                                        txtSumProtGroupPairs.setText(ppiCount + "");
                                    
                                }
                            });

                        } else if (linkCount == maxLinkCount) {
                            maxLinkminPsmfdr = Math.min(psmfdr, maxLinkminPsmfdr);
                            maxLinkminPepfdr = Math.min(pepfdr, maxLinkminPepfdr);
                            maxLinkminProtfdr = Math.min(protfdr, maxLinkminProtfdr);
                            maxLinkmaxPsmfdr = Math.max(psmfdr, maxLinkminPsmfdr);
                            maxLinkmaxPepfdr = Math.max(pepfdr, maxLinkminPepfdr);;
                            maxLinkmaxProtfdr = Math.max(protfdr, maxLinkminProtfdr);;
                        }

                        if (stopMaximizing) {
                            break psmloop;
                        }

                    }


                }



            }

            setStatus("Max Link Round: " + optimizingRound + " - " + maxLinkCount + " links");
            
            // no improvement for the last few rounds?
            if ((maxLinkCount == lastMaxLinkCount && --countDown == 0) || stopMaximizing) {
                optimizing = false;
                Logger.getLogger(this.getClass().getName()).log(Level.INFO, sb.toString());

                final int foundPPICount = maxLinkCount;
                final Double setPsmfdr = maxLinkPsmfdr;
                final Double setPepfdr = maxLinkPepfdr;
                final Double setProtfdr = maxLinkProtfdr;
                boolean isSet = false;
                int trySet = 10;
                while (!isSet && trySet > 0) {
                    try {
                        javax.swing.SwingUtilities.invokeAndWait(new Runnable() {
                            public void run() {
                                fdrSettingsComplete.setPSMFDR(setPsmfdr);
                                fdrSettingsComplete.setPeptidePairFDR(setPepfdr);
                                fdrSettingsComplete.setProteinGroupFDR(setProtfdr);

                                fdrSettingsSimple.setPSMFDR(setPsmfdr);
                                fdrSettingsSimple.setPeptidePairFDR(setPepfdr);
                                fdrSettingsSimple.setProteinGroupFDR(setProtfdr);
                            }
                        });
                        isSet = true;
                    } catch (InterruptedException ex) {
                        trySet--;
                        Logger.getLogger(FDRGUI.class.getName()).log(Level.SEVERE, null, ex);
                    } catch (InvocationTargetException ex) {
                        trySet--;
                        Logger.getLogger(FDRGUI.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }
                if (trySet < 0) {
                    JOptionPane.showMessageDialog(this, "Could not set the values for fdr-calculation\n"
                            + "\nPSM fdr:    " + setPsmfdr
                            + "\nPeptide fdr:" + setPepfdr
                            + "\nProtein FDR:" + setProtfdr, "Failed to set Parameters", JOptionPane.ERROR_MESSAGE);
                    setEnableRead(true);
                    setEnableCalc(true);
                    setEnableWrite(true);
                    stopMaximizing = false;
                    return;
                }

                javax.swing.SwingUtilities.invokeLater(new Runnable() {
                    public void run() {
                        JOptionPane.showMessageDialog(rootPane, "found " + foundPPICount + " links for following settings \n"
                                + "\nPSM fdr:    " + setPsmfdr
                                + "\nPeptide fdr:" + setPepfdr
                                + "\nProtein FDR:" + setProtfdr, "best parameters found for max protein group pairs", JOptionPane.INFORMATION_MESSAGE);
                    }
                });

                innerFDRCalculation(setPsmfdr, setPepfdr, setProtfdr, maxlinkfdr, maxppifdr, getFdr(), reportFactor, filterToUniquePSM);

                setEnableRead(true);
                setEnableCalc(true);
                setEnableWrite(true);

                stopMaximizing = false;

            } else {
                if (maxLinkCount > lastMaxLinkCount) {
                    // yes we improved
                    countDown = 5;
                } else {
                    setStatus("Max Link Round: " + optimizingRound + " - " + maxLinkCount + " links - count down: " + countDown);
                    
                }

                // so see if we make the resoltuion finer
                // can we get a better result?
                maxpsmfdr = Math.min(maxLinkmaxPsmfdr + psmfdrStep, maxpsmfdrabs);
                maxpepfdr = Math.min(maxLinkmaxPepfdr + pepfdrStep, maxpepfdrabs);
                maxprotfdr = Math.min(maxLinkmaxProtfdr + protfdrStep, maxprotfdrabs);

                minpsmfdr = Math.max(maxLinkminPsmfdr - psmfdrStep, 0);
                minpepfdr = Math.max(maxLinkminPepfdr - pepfdrStep, 0);
                minprotfdr = Math.max(maxLinkminProtfdr - protfdrStep, 0);


                psmfdrStep = (maxpsmfdr - minpsmfdr) / steps;
                pepfdrStep = (maxpepfdr - minpepfdr) / steps;
                protfdrStep = (maxprotfdr - minprotfdr) / steps;
                
                

            }


        }


    }
    
    private void changeFDRSettings(java.awt.event.ActionEvent evt) {                                            

        if (rbFDRSimple.isSelected()) {
            if (fdrSettings != null)
                fdrSettingsSimple.setAll(fdrSettings);
            if (spFDRSettingsWrapper.getComponentCount() > 3) {
                spFDRSettingsWrapper.remove(spFDRSettingsWrapper.getComponent(3));
            }
            spFDRSettingsWrapper.setViewportView(fdrSettingsSimple);
            fdrSettings = fdrSettingsSimple;
        } else if (rbFDRComplete.isSelected()) {
            if (fdrSettings != null)
                fdrSettingsComplete.setAll(fdrSettings);
            if (spFDRSettingsWrapper.getComponentCount() > 3) {
                spFDRSettingsWrapper.remove(spFDRSettingsWrapper.getComponent(3));
            }
            spFDRSettingsWrapper.setViewportView(fdrSettingsComplete);
            fdrSettings = fdrSettingsComplete;
        }
        // TODO add your handling code here:
    }      
    

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        bgSeparator = new javax.swing.ButtonGroup();
        bgScoreDirectionMzIdentML = new javax.swing.ButtonGroup();
        cbCSVHeaders = new javax.swing.JComboBox();
        cbCSVHeaderOptional = new javax.swing.JComboBox();
        bgFDRSettingType = new javax.swing.ButtonGroup();
        fdrSettingsComplete = new org.rappsilber.fdr.gui.components.FDRSettingsComplete();
        fdrSettingsSimple = new org.rappsilber.fdr.gui.components.FDRSettingsSimple();
        txtStatus = new javax.swing.JTextField();
        jScrollPane6 = new javax.swing.JScrollPane();
        jTabbedPane1 = new javax.swing.JTabbedPane();
        jPanel2 = new javax.swing.JPanel();
        tpInput = new javax.swing.JTabbedPane();
        jPanel6 = new javax.swing.JPanel();
        csvSelect = new org.rappsilber.fdr.gui.components.CSVSelection();
        jPanel11 = new javax.swing.JPanel();
        fbMZIdentMLIn = new org.rappsilber.gui.components.FileBrowser();
        jLabel19 = new javax.swing.JLabel();
        btnReadMZIdent = new javax.swing.JButton();
        jLabel22 = new javax.swing.JLabel();
        rbMZHighBetter = new javax.swing.JRadioButton();
        rbMZLowBetter = new javax.swing.JRadioButton();
        cbMZMatchScoreName = new javax.swing.JComboBox();
        spExtraLong = new javax.swing.JPanel();
        pDatabseSize = new javax.swing.JPanel();
        spDecoyDBProt = new javax.swing.JSpinner();
        spTargetDBProt = new javax.swing.JSpinner();
        lblProtein = new javax.swing.JLabel();
        lblPeptide = new javax.swing.JLabel();
        spTargetDB = new javax.swing.JSpinner();
        spDecoyDB = new javax.swing.JSpinner();
        lblDecoyDB = new javax.swing.JLabel();
        lblTargetDB = new javax.swing.JLabel();
        spDecoyDBLinks = new javax.swing.JSpinner();
        spTargetDBLinks = new javax.swing.JSpinner();
        lblLinkDB = new javax.swing.JLabel();
        ckDBSize = new javax.swing.JCheckBox();
        rbFDRSimple = new javax.swing.JRadioButton();
        rbFDRComplete = new javax.swing.JRadioButton();
        spFDRSettingsWrapper = new javax.swing.JScrollPane();
        pFDRGroups = new javax.swing.JPanel();
        spPepLength = new javax.swing.JScrollPane();
        tblPepLength = new javax.swing.JTable();
        ckIgnoreGroups1 = new javax.swing.JCheckBox();
        spDistanceGroup = new org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner();
        jLabel1 = new javax.swing.JLabel();
        ckDefineGroups = new javax.swing.JCheckBox();
        pResult = new javax.swing.JPanel();
        jTabbedPane3 = new javax.swing.JTabbedPane();
        jPanel8 = new javax.swing.JPanel();
        jLabel12 = new javax.swing.JLabel();
        txtSumInput = new javax.swing.JTextField();
        jPanel9 = new javax.swing.JPanel();
        btnPSMInfo = new javax.swing.JButton();
        jPanel18 = new javax.swing.JPanel();
        jLabel34 = new javax.swing.JLabel();
        jLabel35 = new javax.swing.JLabel();
        lblSumBetween1 = new javax.swing.JLabel();
        lblSumInternal1 = new javax.swing.JLabel();
        jPanel15 = new javax.swing.JPanel();
        jLabel13 = new javax.swing.JLabel();
        txtSumPSM = new javax.swing.JTextField();
        txtSumPSMXL = new javax.swing.JTextField();
        txtSumPSMLinear = new javax.swing.JTextField();
        jPanel16 = new javax.swing.JPanel();
        jLabel14 = new javax.swing.JLabel();
        txtSumPepPairs = new javax.swing.JTextField();
        txtSumPepPairsXL = new javax.swing.JTextField();
        txtSumPepPairsLinear = new javax.swing.JTextField();
        btnPepInfo = new javax.swing.JButton();
        jPanel20 = new javax.swing.JPanel();
        jLabel15 = new javax.swing.JLabel();
        txtSumProtGroups = new javax.swing.JTextField();
        jLabel30 = new javax.swing.JLabel();
        jLabel31 = new javax.swing.JLabel();
        btnProtInfo = new javax.swing.JButton();
        jPanel12 = new javax.swing.JPanel();
        jLabel33 = new javax.swing.JLabel();
        jLabel32 = new javax.swing.JLabel();
        lblSumBetween = new javax.swing.JLabel();
        lblSumInternal = new javax.swing.JLabel();
        jPanel21 = new javax.swing.JPanel();
        jLabel16 = new javax.swing.JLabel();
        txtSumLinks = new javax.swing.JTextField();
        txtSumLinksBetween = new javax.swing.JTextField();
        txtSumLinksInternal = new javax.swing.JTextField();
        btnLinkInfo = new javax.swing.JButton();
        jPanel22 = new javax.swing.JPanel();
        jLabel17 = new javax.swing.JLabel();
        txtSumProtGroupPairs = new javax.swing.JTextField();
        txtSumProtGroupPairsBetween = new javax.swing.JTextField();
        txtSumProtGroupPairsInternal = new javax.swing.JTextField();
        btnPPIInfo = new javax.swing.JButton();
        jPanel10 = new javax.swing.JPanel();
        jLabel10 = new javax.swing.JLabel();
        rbTSV = new javax.swing.JRadioButton();
        rbCSV = new javax.swing.JRadioButton();
        btnWrite = new javax.swing.JButton();
        fbFolder = new org.rappsilber.gui.components.FileBrowser();
        jLabel9 = new javax.swing.JLabel();
        txtBaseName = new javax.swing.JTextField();
        ckPrePostAA = new javax.swing.JCheckBox();
        jPanel14 = new javax.swing.JPanel();
        fbMzIdentMLOut = new org.rappsilber.gui.components.FileBrowser();
        btnWriteMzIdentML = new javax.swing.JButton();
        txtmzIdentOwnerFirst = new javax.swing.JTextField();
        jLabel2 = new javax.swing.JLabel();
        txtmzIdentOwnerLast = new javax.swing.JTextField();
        jLabel3 = new javax.swing.JLabel();
        txtmzIdentOwnerEmail = new javax.swing.JTextField();
        jScrollPane1 = new javax.swing.JScrollPane();
        txtmzIdentAdress = new javax.swing.JTextArea();
        jLabel4 = new javax.swing.JLabel();
        jLabel5 = new javax.swing.JLabel();
        txtmzIdentOwnerOrg = new javax.swing.JTextField();
        pLog = new javax.swing.JPanel();
        spLog = new javax.swing.JScrollPane();
        txtLog = new javax.swing.JTextArea();
        memory2 = new org.rappsilber.gui.components.memory.Memory();
        jButton3 = new javax.swing.JButton();
        pAbout = new javax.swing.JPanel();
        jTabbedPane4 = new javax.swing.JTabbedPane();
        pVersion = new javax.swing.JPanel();
        jLabel29 = new javax.swing.JLabel();
        txtXiFDRVersion = new javax.swing.JTextField();
        jScrollPane3 = new javax.swing.JScrollPane();
        editAbout = new javax.swing.JEditorPane();
        jScrollPane5 = new javax.swing.JScrollPane();
        editAboutCSV = new javax.swing.JEditorPane();
        jScrollPane4 = new javax.swing.JScrollPane();
        editAboutMzIdentML = new javax.swing.JEditorPane();

        cbCSVHeaders.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Item 1", "Item 2", "Item 3", "Item 4" }));

        cbCSVHeaderOptional.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Item 1", "Item 2", "Item 3", "Item 4" }));

        fdrSettingsSimple.setBorder(javax.swing.BorderFactory.createEtchedBorder());

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);
        setTitle("Cross-Link FDR");
        addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                formMouseClicked(evt);
            }
        });
        addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusGained(java.awt.event.FocusEvent evt) {
                formFocusGained(evt);
            }
        });

        txtStatus.setEditable(false);
        txtStatus.setText("status");

        jTabbedPane1.setPreferredSize(new java.awt.Dimension(700, 400));

        csvSelect.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                csvSelectActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout jPanel6Layout = new javax.swing.GroupLayout(jPanel6);
        jPanel6.setLayout(jPanel6Layout);
        jPanel6Layout.setHorizontalGroup(
            jPanel6Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel6Layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(csvSelect, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
        );
        jPanel6Layout.setVerticalGroup(
            jPanel6Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel6Layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(csvSelect, javax.swing.GroupLayout.DEFAULT_SIZE, 391, Short.MAX_VALUE)
                .addContainerGap())
        );

        tpInput.addTab("CSV", jPanel6);

        fbMZIdentMLIn.setDescription("MZIdentML-Files");
        fbMZIdentMLIn.setExtensions(new String[] {"mzid"});

        jLabel19.setText("mzid-File");

        btnReadMZIdent.setText("Read");
        btnReadMZIdent.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnReadMZIdentActionPerformed(evt);
            }
        });

        jLabel22.setText("PSM-Score");

        bgScoreDirectionMzIdentML.add(rbMZHighBetter);
        rbMZHighBetter.setSelected(true);
        rbMZHighBetter.setText("High better");
        rbMZHighBetter.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                rbMZHighBetterActionPerformed(evt);
            }
        });

        bgScoreDirectionMzIdentML.add(rbMZLowBetter);
        rbMZLowBetter.setText("Low better");
        rbMZLowBetter.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                rbMZLowBetterActionPerformed(evt);
            }
        });

        cbMZMatchScoreName.setEditable(true);
        cbMZMatchScoreName.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "Score", "pvalue", "FDRScore", "xi:score", "MS:1001143", "search engine specific score for PSMs" }));
        cbMZMatchScoreName.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                cbMZMatchScoreNameActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout jPanel11Layout = new javax.swing.GroupLayout(jPanel11);
        jPanel11.setLayout(jPanel11Layout);
        jPanel11Layout.setHorizontalGroup(
            jPanel11Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel11Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel11Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(jPanel11Layout.createSequentialGroup()
                        .addComponent(jLabel19)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(fbMZIdentMLIn, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(btnReadMZIdent))
                    .addGroup(jPanel11Layout.createSequentialGroup()
                        .addComponent(jLabel22)
                        .addGap(239, 239, 239)
                        .addGroup(jPanel11Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(jPanel11Layout.createSequentialGroup()
                                .addGroup(jPanel11Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addComponent(rbMZLowBetter)
                                    .addComponent(rbMZHighBetter))
                                .addGap(88, 270, Short.MAX_VALUE))
                            .addComponent(cbMZMatchScoreName, javax.swing.GroupLayout.Alignment.TRAILING, 0, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))))
                .addContainerGap())
        );
        jPanel11Layout.setVerticalGroup(
            jPanel11Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel11Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel11Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel22)
                    .addComponent(cbMZMatchScoreName, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addGap(1, 1, 1)
                .addComponent(rbMZHighBetter)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(rbMZLowBetter)
                .addGap(18, 18, 18)
                .addGroup(jPanel11Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addGroup(jPanel11Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        .addComponent(btnReadMZIdent)
                        .addComponent(jLabel19))
                    .addComponent(fbMZIdentMLIn, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(289, Short.MAX_VALUE))
        );

        tpInput.addTab("mzIdentML", jPanel11);

        javax.swing.GroupLayout jPanel2Layout = new javax.swing.GroupLayout(jPanel2);
        jPanel2.setLayout(jPanel2Layout);
        jPanel2Layout.setHorizontalGroup(
            jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(tpInput)
        );
        jPanel2Layout.setVerticalGroup(
            jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel2Layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(tpInput)
                .addContainerGap())
        );

        jTabbedPane1.addTab("Input", jPanel2);

        pDatabseSize.setBorder(javax.swing.BorderFactory.createTitledBorder("Database Sizes"));
        pDatabseSize.setEnabled(false);

        spDecoyDBProt.setModel(new javax.swing.SpinnerNumberModel(Double.valueOf(9.99999999E8d), Double.valueOf(0.001d), null, Double.valueOf(1.0d)));
        spDecoyDBProt.setEnabled(false);

        spTargetDBProt.setModel(new javax.swing.SpinnerNumberModel(Double.valueOf(9.99999999E8d), Double.valueOf(0.001d), null, Double.valueOf(1.0d)));
        spTargetDBProt.setEnabled(false);

        lblProtein.setText("Protein");
        lblProtein.setEnabled(false);

        lblPeptide.setText("Peptide");
        lblPeptide.setEnabled(false);

        spTargetDB.setModel(new javax.swing.SpinnerNumberModel(Double.valueOf(9.99999999E8d), Double.valueOf(0.001d), null, Double.valueOf(1.0d)));
        spTargetDB.setEnabled(false);

        spDecoyDB.setModel(new javax.swing.SpinnerNumberModel(Double.valueOf(9.99999999E8d), Double.valueOf(0.001d), null, Double.valueOf(1.0d)));
        spDecoyDB.setEnabled(false);

        lblDecoyDB.setText("Decoy");
        lblDecoyDB.setEnabled(false);

        lblTargetDB.setText("Target");
        lblTargetDB.setEnabled(false);

        spDecoyDBLinks.setModel(new javax.swing.SpinnerNumberModel(Double.valueOf(9.99999999E8d), Double.valueOf(0.001d), null, Double.valueOf(1.0d)));
        spDecoyDBLinks.setEnabled(false);

        spTargetDBLinks.setModel(new javax.swing.SpinnerNumberModel(Double.valueOf(9.99999999E8d), Double.valueOf(0.001d), null, Double.valueOf(1.0d)));
        spTargetDBLinks.setEnabled(false);

        lblLinkDB.setText("Links");
        lblLinkDB.setEnabled(false);

        javax.swing.GroupLayout pDatabseSizeLayout = new javax.swing.GroupLayout(pDatabseSize);
        pDatabseSize.setLayout(pDatabseSizeLayout);
        pDatabseSizeLayout.setHorizontalGroup(
            pDatabseSizeLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(pDatabseSizeLayout.createSequentialGroup()
                .addGap(78, 78, 78)
                .addGroup(pDatabseSizeLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(lblDecoyDB, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(lblTargetDB, javax.swing.GroupLayout.PREFERRED_SIZE, 56, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(pDatabseSizeLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(spDecoyDB, javax.swing.GroupLayout.PREFERRED_SIZE, 100, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(lblPeptide, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(spTargetDB, javax.swing.GroupLayout.PREFERRED_SIZE, 100, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(pDatabseSizeLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(spDecoyDBLinks, javax.swing.GroupLayout.PREFERRED_SIZE, 100, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(spTargetDBLinks, javax.swing.GroupLayout.PREFERRED_SIZE, 100, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(lblLinkDB, javax.swing.GroupLayout.PREFERRED_SIZE, 100, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(pDatabseSizeLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(spDecoyDBProt, javax.swing.GroupLayout.PREFERRED_SIZE, 100, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(spTargetDBProt, javax.swing.GroupLayout.PREFERRED_SIZE, 100, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(lblProtein, javax.swing.GroupLayout.PREFERRED_SIZE, 100, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(209, Short.MAX_VALUE))
        );
        pDatabseSizeLayout.setVerticalGroup(
            pDatabseSizeLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(pDatabseSizeLayout.createSequentialGroup()
                .addGroup(pDatabseSizeLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(pDatabseSizeLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                        .addComponent(lblPeptide)
                        .addComponent(lblProtein))
                    .addComponent(lblLinkDB))
                .addGroup(pDatabseSizeLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(spTargetDBLinks, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addGroup(pDatabseSizeLayout.createSequentialGroup()
                        .addGroup(pDatabseSizeLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(spTargetDB, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(lblTargetDB)
                            .addComponent(spTargetDBProt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(pDatabseSizeLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(lblDecoyDB)
                            .addComponent(spDecoyDB, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                    .addGroup(pDatabseSizeLayout.createSequentialGroup()
                        .addGap(26, 26, 26)
                        .addGroup(pDatabseSizeLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(spDecoyDBLinks, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(spDecoyDBProt, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))))
        );

        ckDBSize.setText("Define databse size");
        ckDBSize.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                ckDBSizeActionPerformed(evt);
            }
        });

        bgFDRSettingType.add(rbFDRSimple);
        rbFDRSimple.setSelected(true);
        rbFDRSimple.setText("Simple FDR");
        rbFDRSimple.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                changFDRSettingsInterface(evt);
            }
        });

        bgFDRSettingType.add(rbFDRComplete);
        rbFDRComplete.setText("complete FDR");
        rbFDRComplete.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                rbFDRCompleteActionPerformed(evt);
            }
        });

        spPepLength.setMinimumSize(new java.awt.Dimension(50, 50));
        spPepLength.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                spPepLengthFocusLost(evt);
            }
        });

        tblPepLength.setModel(new javax.swing.table.DefaultTableModel(
            new Object [][] {
                { new Integer(4)},
                {null}
            },
            new String [] {
                "Peptide Length"
            }
        ) {
            Class[] types = new Class [] {
                java.lang.Integer.class
            };

            public Class getColumnClass(int columnIndex) {
                return types [columnIndex];
            }
        });
        spPepLength.setViewportView(tblPepLength);

        ckIgnoreGroups1.setText("Ignore Groups");
        ckIgnoreGroups1.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                ckIgnoreGroups1ActionPerformed(evt);
            }
        });

        spDistanceGroup.setSpecialValueText("No distance group");

        jLabel1.setText("Sequnce Distance");

        javax.swing.GroupLayout pFDRGroupsLayout = new javax.swing.GroupLayout(pFDRGroups);
        pFDRGroups.setLayout(pFDRGroupsLayout);
        pFDRGroupsLayout.setHorizontalGroup(
            pFDRGroupsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(spPepLength, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.PREFERRED_SIZE, 0, Short.MAX_VALUE)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, pFDRGroupsLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(ckIgnoreGroups1, javax.swing.GroupLayout.DEFAULT_SIZE, 134, Short.MAX_VALUE)
                .addContainerGap())
            .addComponent(spDistanceGroup, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
            .addComponent(jLabel1, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
        );
        pFDRGroupsLayout.setVerticalGroup(
            pFDRGroupsLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(pFDRGroupsLayout.createSequentialGroup()
                .addComponent(ckIgnoreGroups1)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(spPepLength, javax.swing.GroupLayout.DEFAULT_SIZE, 226, Short.MAX_VALUE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jLabel1)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(spDistanceGroup, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
        );

        ckDefineGroups.setText("Define Groups");
        ckDefineGroups.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                ckDefineGroupsActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout spExtraLongLayout = new javax.swing.GroupLayout(spExtraLong);
        spExtraLong.setLayout(spExtraLongLayout);
        spExtraLongLayout.setHorizontalGroup(
            spExtraLongLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(spExtraLongLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(spExtraLongLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(pDatabseSize, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addGroup(spExtraLongLayout.createSequentialGroup()
                        .addComponent(ckDBSize)
                        .addGap(0, 0, Short.MAX_VALUE))
                    .addGroup(spExtraLongLayout.createSequentialGroup()
                        .addComponent(rbFDRSimple)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(rbFDRComplete)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                        .addComponent(ckDefineGroups, javax.swing.GroupLayout.PREFERRED_SIZE, 134, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, spExtraLongLayout.createSequentialGroup()
                        .addComponent(spFDRSettingsWrapper)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(pFDRGroups, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addContainerGap())
        );
        spExtraLongLayout.setVerticalGroup(
            spExtraLongLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(spExtraLongLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(spExtraLongLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(rbFDRSimple)
                    .addComponent(rbFDRComplete)
                    .addComponent(ckDefineGroups))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(spExtraLongLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(spFDRSettingsWrapper)
                    .addComponent(pFDRGroups, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(ckDBSize)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(pDatabseSize, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(20, 20, 20))
        );

        jTabbedPane1.addTab("FDR Settings", spExtraLong);

        jLabel12.setText("PSM input");

        jPanel9.setBorder(javax.swing.BorderFactory.createTitledBorder("AfterFDR"));

        btnPSMInfo.setText("+");
        btnPSMInfo.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnPSMInfoActionPerformed(evt);
            }
        });

        jPanel18.setLayout(new java.awt.GridLayout(1, 4, 15, 0));
        jPanel18.add(jLabel34);
        jPanel18.add(jLabel35);

        lblSumBetween1.setText("Cross-Linked");
        lblSumBetween1.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jPanel18.add(lblSumBetween1);

        lblSumInternal1.setText("Linear");
        jPanel18.add(lblSumInternal1);

        jPanel15.setLayout(new java.awt.GridLayout(1, 0, 10, 0));

        jLabel13.setText("PSMs");
        jPanel15.add(jLabel13);

        txtSumPSM.setMaximumSize(null);
        txtSumPSM.setMinimumSize(null);
        txtSumPSM.setPreferredSize(null);
        txtSumPSM.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                txtSumPSMActionPerformed(evt);
            }
        });
        jPanel15.add(txtSumPSM);

        txtSumPSMXL.setMaximumSize(null);
        txtSumPSMXL.setMinimumSize(null);
        txtSumPSMXL.setPreferredSize(null);
        txtSumPSMXL.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                txtSumPSMXLActionPerformed(evt);
            }
        });
        jPanel15.add(txtSumPSMXL);

        txtSumPSMLinear.setMaximumSize(null);
        txtSumPSMLinear.setMinimumSize(null);
        txtSumPSMLinear.setPreferredSize(null);
        txtSumPSMLinear.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                txtSumPSMLinearActionPerformed(evt);
            }
        });
        jPanel15.add(txtSumPSMLinear);

        jPanel16.setLayout(new java.awt.GridLayout(1, 0, 10, 0));

        jLabel14.setText("Peptides");
        jPanel16.add(jLabel14);

        txtSumPepPairs.setMaximumSize(null);
        txtSumPepPairs.setMinimumSize(null);
        txtSumPepPairs.setPreferredSize(null);
        jPanel16.add(txtSumPepPairs);

        txtSumPepPairsXL.setMaximumSize(null);
        txtSumPepPairsXL.setMinimumSize(null);
        txtSumPepPairsXL.setPreferredSize(null);
        jPanel16.add(txtSumPepPairsXL);

        txtSumPepPairsLinear.setMaximumSize(null);
        txtSumPepPairsLinear.setMinimumSize(null);
        txtSumPepPairsLinear.setPreferredSize(null);
        jPanel16.add(txtSumPepPairsLinear);

        btnPepInfo.setText("+");
        btnPepInfo.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnPepInfoActionPerformed(evt);
            }
        });

        jPanel20.setLayout(new java.awt.GridLayout(1, 0, 10, 0));

        jLabel15.setText("Protein Groups");
        jPanel20.add(jLabel15);

        txtSumProtGroups.setMaximumSize(null);
        txtSumProtGroups.setMinimumSize(null);
        txtSumProtGroups.setPreferredSize(null);
        jPanel20.add(txtSumProtGroups);
        jPanel20.add(jLabel30);
        jPanel20.add(jLabel31);

        btnProtInfo.setText("+");
        btnProtInfo.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnProtInfoActionPerformed(evt);
            }
        });

        jPanel12.setLayout(new java.awt.GridLayout(1, 4, 15, 0));
        jPanel12.add(jLabel33);
        jPanel12.add(jLabel32);

        lblSumBetween.setText("Between");
        lblSumBetween.setHorizontalTextPosition(javax.swing.SwingConstants.CENTER);
        jPanel12.add(lblSumBetween);

        lblSumInternal.setText("Within");
        jPanel12.add(lblSumInternal);

        jPanel21.setLayout(new java.awt.GridLayout(1, 0, 10, 0));

        jLabel16.setText("Links");
        jPanel21.add(jLabel16);

        txtSumLinks.setMaximumSize(null);
        txtSumLinks.setMinimumSize(null);
        txtSumLinks.setPreferredSize(null);
        jPanel21.add(txtSumLinks);

        txtSumLinksBetween.setMaximumSize(null);
        txtSumLinksBetween.setMinimumSize(null);
        txtSumLinksBetween.setPreferredSize(null);
        jPanel21.add(txtSumLinksBetween);

        txtSumLinksInternal.setMaximumSize(null);
        txtSumLinksInternal.setMinimumSize(null);
        txtSumLinksInternal.setPreferredSize(null);
        jPanel21.add(txtSumLinksInternal);

        btnLinkInfo.setText("+");
        btnLinkInfo.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnLinkInfoActionPerformed(evt);
            }
        });

        jPanel22.setLayout(new java.awt.GridLayout(1, 0, 10, 0));

        jLabel17.setText("Protein Group Pairs");
        jPanel22.add(jLabel17);

        txtSumProtGroupPairs.setMaximumSize(null);
        txtSumProtGroupPairs.setMinimumSize(null);
        txtSumProtGroupPairs.setPreferredSize(null);
        jPanel22.add(txtSumProtGroupPairs);

        txtSumProtGroupPairsBetween.setMaximumSize(null);
        txtSumProtGroupPairsBetween.setMinimumSize(null);
        txtSumProtGroupPairsBetween.setPreferredSize(null);
        jPanel22.add(txtSumProtGroupPairsBetween);

        txtSumProtGroupPairsInternal.setMaximumSize(null);
        txtSumProtGroupPairsInternal.setMinimumSize(null);
        txtSumProtGroupPairsInternal.setPreferredSize(null);
        jPanel22.add(txtSumProtGroupPairsInternal);

        btnPPIInfo.setText("+");
        btnPPIInfo.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnPPIInfoActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout jPanel9Layout = new javax.swing.GroupLayout(jPanel9);
        jPanel9.setLayout(jPanel9Layout);
        jPanel9Layout.setHorizontalGroup(
            jPanel9Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel9Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel9Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(jPanel9Layout.createSequentialGroup()
                        .addGroup(jPanel9Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addComponent(jPanel22, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.PREFERRED_SIZE, 0, Short.MAX_VALUE)
                            .addComponent(jPanel21, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(jPanel12, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(jPanel9Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(btnLinkInfo)
                            .addComponent(btnPPIInfo)))
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, jPanel9Layout.createSequentialGroup()
                        .addGroup(jPanel9Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addComponent(jPanel18, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(jPanel20, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, 523, Short.MAX_VALUE)
                            .addComponent(jPanel16, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(jPanel15, javax.swing.GroupLayout.Alignment.LEADING, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(jPanel9Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(btnPSMInfo)
                            .addComponent(btnPepInfo)
                            .addComponent(btnProtInfo))))
                .addContainerGap())
        );
        jPanel9Layout.setVerticalGroup(
            jPanel9Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel9Layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jPanel18, javax.swing.GroupLayout.PREFERRED_SIZE, 20, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel9Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(btnProtInfo)
                    .addGroup(jPanel9Layout.createSequentialGroup()
                        .addGroup(jPanel9Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                            .addComponent(btnPSMInfo, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(jPanel15, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addGroup(jPanel9Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                            .addComponent(jPanel16, javax.swing.GroupLayout.PREFERRED_SIZE, 25, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(btnPepInfo))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(jPanel20, javax.swing.GroupLayout.PREFERRED_SIZE, 25, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jPanel12, javax.swing.GroupLayout.PREFERRED_SIZE, 18, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel9Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                    .addComponent(btnLinkInfo, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(jPanel21, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED, 32, Short.MAX_VALUE)
                .addGroup(jPanel9Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(btnPPIInfo, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(jPanel22, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addContainerGap())
        );

        javax.swing.GroupLayout jPanel8Layout = new javax.swing.GroupLayout(jPanel8);
        jPanel8.setLayout(jPanel8Layout);
        jPanel8Layout.setHorizontalGroup(
            jPanel8Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel8Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel8Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(jLabel12, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(txtSumInput))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jPanel9, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addContainerGap())
        );
        jPanel8Layout.setVerticalGroup(
            jPanel8Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel8Layout.createSequentialGroup()
                .addGroup(jPanel8Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(jPanel8Layout.createSequentialGroup()
                        .addComponent(jLabel12)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(txtSumInput, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addComponent(jPanel9, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(133, Short.MAX_VALUE))
        );

        jTabbedPane3.addTab("Summary", jPanel8);

        jLabel10.setText("Base Name");

        bgSeparator.add(rbTSV);
        rbTSV.setText("Tab Separated");

        bgSeparator.add(rbCSV);
        rbCSV.setSelected(true);
        rbCSV.setText("Comma Separated");

        btnWrite.setText("Write");
        btnWrite.setEnabled(false);
        btnWrite.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnWriteActionPerformed(evt);
            }
        });

        fbFolder.setDirectoryOnly(true);
        fbFolder.setLoad(false);

        jLabel9.setText("Folder");

        txtBaseName.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                txtBaseNameActionPerformed(evt);
            }
        });

        ckPrePostAA.setText("Pre and post amino-acids");

        javax.swing.GroupLayout jPanel10Layout = new javax.swing.GroupLayout(jPanel10);
        jPanel10.setLayout(jPanel10Layout);
        jPanel10Layout.setHorizontalGroup(
            jPanel10Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel10Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel10Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jLabel9)
                    .addComponent(jLabel10))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(jPanel10Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(fbFolder, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addComponent(txtBaseName)
                    .addGroup(jPanel10Layout.createSequentialGroup()
                        .addComponent(rbTSV)
                        .addGap(95, 95, 95)
                        .addComponent(ckPrePostAA)
                        .addGap(0, 0, Short.MAX_VALUE))
                    .addGroup(jPanel10Layout.createSequentialGroup()
                        .addComponent(rbCSV)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 354, Short.MAX_VALUE)
                        .addComponent(btnWrite)))
                .addGap(23, 23, 23))
        );
        jPanel10Layout.setVerticalGroup(
            jPanel10Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel10Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel10Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(fbFolder, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel9))
                .addGap(18, 18, 18)
                .addGroup(jPanel10Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(txtBaseName, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel10))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addGroup(jPanel10Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(rbTSV)
                    .addComponent(ckPrePostAA))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel10Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(rbCSV)
                    .addComponent(btnWrite))
                .addContainerGap(283, Short.MAX_VALUE))
        );

        jTabbedPane3.addTab("CSV/TSV", jPanel10);

        fbMzIdentMLOut.setDescription("mzIdentML files");
        fbMzIdentMLOut.setExtensions(new String[] {"mzid"});
        fbMzIdentMLOut.setLoad(false);
        fbMzIdentMLOut.setLocalPropertyKey("LastAccessedmzIdentMLFolder");

        btnWriteMzIdentML.setText("Write");
        btnWriteMzIdentML.setEnabled(false);
        btnWriteMzIdentML.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                btnWriteMzIdentMLActionPerformed(evt);
            }
        });

        txtmzIdentOwnerFirst.setText("First");
        txtmzIdentOwnerFirst.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                txtmzIdentOwnerFirstActionPerformed(evt);
            }
        });

        jLabel2.setText("Document Owner");

        txtmzIdentOwnerLast.setText("Last");
        txtmzIdentOwnerLast.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                txtmzIdentOwnerLastActionPerformed(evt);
            }
        });

        jLabel3.setText("E-Mail");

        txtmzIdentOwnerEmail.setText("reseacher@organisation.org");
        txtmzIdentOwnerEmail.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                txtmzIdentOwnerEmailActionPerformed(evt);
            }
        });

        txtmzIdentAdress.setColumns(20);
        txtmzIdentAdress.setRows(5);
        jScrollPane1.setViewportView(txtmzIdentAdress);

        jLabel4.setText("Address");

        jLabel5.setText("Organisation");

        txtmzIdentOwnerOrg.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                txtmzIdentOwnerOrgActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout jPanel14Layout = new javax.swing.GroupLayout(jPanel14);
        jPanel14.setLayout(jPanel14Layout);
        jPanel14Layout.setHorizontalGroup(
            jPanel14Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel14Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel14Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(fbMzIdentMLOut, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .addGroup(jPanel14Layout.createSequentialGroup()
                        .addGroup(jPanel14Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jLabel2)
                            .addComponent(jLabel3)
                            .addComponent(jLabel4)
                            .addComponent(jLabel5))
                        .addGap(18, 18, 18)
                        .addGroup(jPanel14Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                            .addComponent(jScrollPane1)
                            .addGroup(jPanel14Layout.createSequentialGroup()
                                .addComponent(txtmzIdentOwnerFirst, javax.swing.GroupLayout.PREFERRED_SIZE, 80, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addGap(18, 18, 18)
                                .addComponent(txtmzIdentOwnerLast, javax.swing.GroupLayout.PREFERRED_SIZE, 80, javax.swing.GroupLayout.PREFERRED_SIZE))
                            .addComponent(txtmzIdentOwnerEmail)
                            .addComponent(txtmzIdentOwnerOrg))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 254, Short.MAX_VALUE)
                        .addComponent(btnWriteMzIdentML)))
                .addContainerGap())
        );
        jPanel14Layout.setVerticalGroup(
            jPanel14Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel14Layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(fbMzIdentMLOut, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel14Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(btnWriteMzIdentML)
                    .addComponent(jLabel2)
                    .addComponent(txtmzIdentOwnerFirst, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(txtmzIdentOwnerLast, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel14Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel3)
                    .addComponent(txtmzIdentOwnerEmail, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(jPanel14Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(txtmzIdentOwnerOrg, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel5))
                .addGap(9, 9, 9)
                .addGroup(jPanel14Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(jScrollPane1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel4))
                .addContainerGap(210, Short.MAX_VALUE))
        );

        jTabbedPane3.addTab("mzIdentML", jPanel14);

        javax.swing.GroupLayout pResultLayout = new javax.swing.GroupLayout(pResult);
        pResult.setLayout(pResultLayout);
        pResultLayout.setHorizontalGroup(
            pResultLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jTabbedPane3, javax.swing.GroupLayout.Alignment.TRAILING, javax.swing.GroupLayout.DEFAULT_SIZE, 719, Short.MAX_VALUE)
        );
        pResultLayout.setVerticalGroup(
            pResultLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(pResultLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jTabbedPane3)
                .addContainerGap())
        );

        jTabbedPane1.addTab("Results", pResult);

        txtLog.setColumns(20);
        txtLog.setRows(5);
        spLog.setViewportView(txtLog);

        jButton3.setText("Paper-Interface");
        jButton3.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jButton3ActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout pLogLayout = new javax.swing.GroupLayout(pLog);
        pLog.setLayout(pLogLayout);
        pLogLayout.setHorizontalGroup(
            pLogLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(pLogLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(pLogLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(spLog)
                    .addComponent(memory2, javax.swing.GroupLayout.DEFAULT_SIZE, 695, Short.MAX_VALUE)
                    .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, pLogLayout.createSequentialGroup()
                        .addGap(0, 0, Short.MAX_VALUE)
                        .addComponent(jButton3)))
                .addContainerGap())
        );
        pLogLayout.setVerticalGroup(
            pLogLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, pLogLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(memory2, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(spLog, javax.swing.GroupLayout.DEFAULT_SIZE, 386, Short.MAX_VALUE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(jButton3)
                .addContainerGap())
        );

        jTabbedPane1.addTab("Log", pLog);

        jLabel29.setText("Version");

        txtXiFDRVersion.setEditable(false);

        javax.swing.GroupLayout pVersionLayout = new javax.swing.GroupLayout(pVersion);
        pVersion.setLayout(pVersionLayout);
        pVersionLayout.setHorizontalGroup(
            pVersionLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(pVersionLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jLabel29)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(txtXiFDRVersion, javax.swing.GroupLayout.PREFERRED_SIZE, 119, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(517, Short.MAX_VALUE))
        );
        pVersionLayout.setVerticalGroup(
            pVersionLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(pVersionLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(pVersionLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel29)
                    .addComponent(txtXiFDRVersion, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(384, Short.MAX_VALUE))
        );

        jTabbedPane4.addTab("Version", pVersion);

        editAbout.setEditable(false);
        jScrollPane3.setViewportView(editAbout);

        jTabbedPane4.addTab("General", jScrollPane3);

        jScrollPane5.setViewportView(editAboutCSV);

        jTabbedPane4.addTab("CSV-Import", jScrollPane5);

        jScrollPane4.setViewportView(editAboutMzIdentML);

        jTabbedPane4.addTab("mzIdentML-import", jScrollPane4);

        javax.swing.GroupLayout pAboutLayout = new javax.swing.GroupLayout(pAbout);
        pAbout.setLayout(pAboutLayout);
        pAboutLayout.setHorizontalGroup(
            pAboutLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jTabbedPane4)
        );
        pAboutLayout.setVerticalGroup(
            pAboutLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(pAboutLayout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jTabbedPane4)
                .addContainerGap())
        );

        jTabbedPane1.addTab("About", pAbout);

        jScrollPane6.setViewportView(jTabbedPane1);

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jScrollPane6, javax.swing.GroupLayout.DEFAULT_SIZE, 727, Short.MAX_VALUE)
            .addComponent(txtStatus)
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addComponent(jScrollPane6, javax.swing.GroupLayout.DEFAULT_SIZE, 496, Short.MAX_VALUE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(txtStatus, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void btnWriteActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnWriteActionPerformed
        writeResults();
    }//GEN-LAST:event_btnWriteActionPerformed

    private void txtBaseNameActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_txtBaseNameActionPerformed
    }//GEN-LAST:event_txtBaseNameActionPerformed

    private void txtSumPSMActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_txtSumPSMActionPerformed
    }//GEN-LAST:event_txtSumPSMActionPerformed

    private void btnReadMZIdentActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnReadMZIdentActionPerformed
        readMZIdentML();
        ckPrePostAA.setVisible(false);
        
    }//GEN-LAST:event_btnReadMZIdentActionPerformed

    private void cbMZMatchScoreNameActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_cbMZMatchScoreNameActionPerformed

        if (!mzScoreDirectionChanged) {
            String name = cbMZMatchScoreName.getModel().getSelectedItem().toString();
            if (name.toLowerCase().contains("pvalue") || name.toLowerCase().contains("fdr")) {
                rbMZLowBetter.setSelected(true);
            }
        }

    }//GEN-LAST:event_cbMZMatchScoreNameActionPerformed

    private void rbMZHighBetterActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_rbMZHighBetterActionPerformed
        mzScoreDirectionChanged = true;
    }//GEN-LAST:event_rbMZHighBetterActionPerformed

    private void rbMZLowBetterActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_rbMZLowBetterActionPerformed
        mzScoreDirectionChanged = true;
    }//GEN-LAST:event_rbMZLowBetterActionPerformed

    private void btnWriteMzIdentMLActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnWriteMzIdentMLActionPerformed
        writeMZIdentML();
        LocalProperties.setProperty("mzIdenMLOwnerFirst",txtmzIdentOwnerFirst.getText());
        LocalProperties.setProperty("mzIdenMLOwnerLast",txtmzIdentOwnerLast.getText());
        LocalProperties.setProperty("mzIdenMLOwnerEmail",txtmzIdentOwnerEmail.getText());
        LocalProperties.setProperty("mzIdenMLOwnerAddress",txtmzIdentAdress.getText());
        LocalProperties.setProperty("mzIdenMLOwnerAddress",txtmzIdentOwnerOrg.getText());

    }//GEN-LAST:event_btnWriteMzIdentMLActionPerformed

    private void txtSumPSMXLActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_txtSumPSMXLActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_txtSumPSMXLActionPerformed

    private void txtSumPSMLinearActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_txtSumPSMLinearActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_txtSumPSMLinearActionPerformed

    private void spPepLengthFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_spPepLengthFocusLost
    }//GEN-LAST:event_spPepLengthFocusLost

    private void ckDefineGroupsActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_ckDefineGroupsActionPerformed

        pFDRGroups.setVisible(ckDefineGroups.isSelected());
        tblPepLength.setEnabled(!ckIgnoreGroups1.isSelected());
        spPepLength.setEnabled(!ckIgnoreGroups1.isSelected());

    }//GEN-LAST:event_ckDefineGroupsActionPerformed

    private void formFocusGained(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_formFocusGained
        
    
        
    }//GEN-LAST:event_formFocusGained

    private void formMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_formMouseClicked
        // TODO add your handling code here:
        this.toFront();
    }//GEN-LAST:event_formMouseClicked

    private void btnPSMInfoActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnPSMInfoActionPerformed
        if (m_result!= null)
            new FDRLevelInformations(m_result.psmFDR, "PSM FDR").setVisible(true);
    }//GEN-LAST:event_btnPSMInfoActionPerformed

    private void btnPepInfoActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnPepInfoActionPerformed
        if (m_result!= null)
            new FDRLevelInformations(m_result.peptidePairFDR, "Peptide-Pair FDR").setVisible(true);
    }//GEN-LAST:event_btnPepInfoActionPerformed

    private void btnProtInfoActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnProtInfoActionPerformed
        if (m_result!= null)
            new FDRLevelInformations(m_result.proteinGroupFDR, "Protein Group FDR").setVisible(true);
    }//GEN-LAST:event_btnProtInfoActionPerformed

    private void btnLinkInfoActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnLinkInfoActionPerformed
        if (m_result!= null)
            new FDRLevelInformations(m_result.proteinGroupLinkFDR, "Protein Group Link FDR").setVisible(true);
    }//GEN-LAST:event_btnLinkInfoActionPerformed

    private void btnPPIInfoActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_btnPPIInfoActionPerformed
        if (m_result!= null)
            new FDRLevelInformations(m_result.proteinGroupPairFDR, "Protein Group Pairs FDR").setVisible(true);
    }//GEN-LAST:event_btnPPIInfoActionPerformed

    private void ckDBSizeActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_ckDBSizeActionPerformed
        
        pDatabseSize.setVisible(ckDBSize.isSelected());
        
        spTargetDB.setEnabled(ckDBSize.isSelected());
        spDecoyDB.setEnabled(ckDBSize.isSelected());
        spTargetDBLinks.setEnabled(ckDBSize.isSelected());
        spDecoyDBLinks.setEnabled(ckDBSize.isSelected());
        spTargetDBProt.setEnabled(ckDBSize.isSelected());
        spDecoyDBProt.setEnabled(ckDBSize.isSelected());
        lblDecoyDB.setEnabled(ckDBSize.isSelected());
        lblTargetDB.setEnabled(ckDBSize.isSelected());
        lblLinkDB.setEnabled(ckDBSize.isSelected());
        lblPeptide.setEnabled(ckDBSize.isSelected());
        lblProtein.setEnabled(ckDBSize.isSelected());

        spTargetDB.setVisible(ckDBSize.isSelected());
        spDecoyDB.setVisible(ckDBSize.isSelected());
        spTargetDBLinks.setVisible(ckDBSize.isSelected());
        spDecoyDBLinks.setVisible(ckDBSize.isSelected());
        spTargetDBProt.setVisible(ckDBSize.isSelected());
        spDecoyDBProt.setVisible(ckDBSize.isSelected());
        lblDecoyDB.setVisible(ckDBSize.isSelected());
        lblTargetDB.setVisible(ckDBSize.isSelected());
        lblLinkDB.setVisible(ckDBSize.isSelected());
        lblPeptide.setVisible(ckDBSize.isSelected());
        lblProtein.setVisible(ckDBSize.isSelected());
        
    }//GEN-LAST:event_ckDBSizeActionPerformed

    private void jButton3ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jButton3ActionPerformed
        
        ckPrePostAA.setVisible(false);
        ckPrePostAA.getParent().remove(ckPrePostAA);

    }//GEN-LAST:event_jButton3ActionPerformed

    private void changFDRSettingsInterface(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_changFDRSettingsInterface

        changeFDRSettings(evt);
        
    }//GEN-LAST:event_changFDRSettingsInterface

    private void rbFDRCompleteActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_rbFDRCompleteActionPerformed
        changeFDRSettings(evt);
    }//GEN-LAST:event_rbFDRCompleteActionPerformed

    private void ckIgnoreGroups1ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_ckIgnoreGroups1ActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_ckIgnoreGroups1ActionPerformed

    private void csvSelectActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_csvSelectActionPerformed
        readCSV();
    }//GEN-LAST:event_csvSelectActionPerformed

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

    
    
//    private void fdrSpinnerMaximumCheck(JSpinner sp, double max) {                                   
//        SpinnerModel sm = sp.getModel();
//        Double fdr = (Double) sm.getValue();
//        if (fdr > max)
//            sm.setValue(max);
//    }                                  

    
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
            java.util.logging.Logger.getLogger(FDRGUI.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (InstantiationException ex) {
            java.util.logging.Logger.getLogger(FDRGUI.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (IllegalAccessException ex) {
            java.util.logging.Logger.getLogger(FDRGUI.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (javax.swing.UnsupportedLookAndFeelException ex) {
            java.util.logging.Logger.getLogger(FDRGUI.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>
        //</editor-fold>

        /* Create and display the form */
        java.awt.EventQueue.invokeLater(new Runnable() {
            public void run() {
                new FDRGUI().setVisible(true);
            }
        });
    }
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.ButtonGroup bgFDRSettingType;
    private javax.swing.ButtonGroup bgScoreDirectionMzIdentML;
    private javax.swing.ButtonGroup bgSeparator;
    private javax.swing.JButton btnLinkInfo;
    private javax.swing.JButton btnPPIInfo;
    private javax.swing.JButton btnPSMInfo;
    private javax.swing.JButton btnPepInfo;
    private javax.swing.JButton btnProtInfo;
    private javax.swing.JButton btnReadMZIdent;
    private javax.swing.JButton btnWrite;
    private javax.swing.JButton btnWriteMzIdentML;
    private javax.swing.JComboBox cbCSVHeaderOptional;
    private javax.swing.JComboBox cbCSVHeaders;
    private javax.swing.JComboBox cbMZMatchScoreName;
    private javax.swing.JCheckBox ckDBSize;
    private javax.swing.JCheckBox ckDefineGroups;
    private javax.swing.JCheckBox ckIgnoreGroups1;
    private javax.swing.JCheckBox ckPrePostAA;
    private org.rappsilber.fdr.gui.components.CSVSelection csvSelect;
    private javax.swing.JEditorPane editAbout;
    private javax.swing.JEditorPane editAboutCSV;
    private javax.swing.JEditorPane editAboutMzIdentML;
    private org.rappsilber.gui.components.FileBrowser fbFolder;
    private org.rappsilber.gui.components.FileBrowser fbMZIdentMLIn;
    private org.rappsilber.gui.components.FileBrowser fbMzIdentMLOut;
    private org.rappsilber.fdr.gui.components.FDRSettingsComplete fdrSettingsComplete;
    private org.rappsilber.fdr.gui.components.FDRSettingsSimple fdrSettingsSimple;
    private javax.swing.JButton jButton3;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel10;
    private javax.swing.JLabel jLabel12;
    private javax.swing.JLabel jLabel13;
    private javax.swing.JLabel jLabel14;
    private javax.swing.JLabel jLabel15;
    private javax.swing.JLabel jLabel16;
    private javax.swing.JLabel jLabel17;
    private javax.swing.JLabel jLabel19;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel22;
    private javax.swing.JLabel jLabel29;
    private javax.swing.JLabel jLabel3;
    private javax.swing.JLabel jLabel30;
    private javax.swing.JLabel jLabel31;
    private javax.swing.JLabel jLabel32;
    private javax.swing.JLabel jLabel33;
    private javax.swing.JLabel jLabel34;
    private javax.swing.JLabel jLabel35;
    private javax.swing.JLabel jLabel4;
    private javax.swing.JLabel jLabel5;
    private javax.swing.JLabel jLabel9;
    private javax.swing.JPanel jPanel10;
    private javax.swing.JPanel jPanel11;
    private javax.swing.JPanel jPanel12;
    private javax.swing.JPanel jPanel14;
    private javax.swing.JPanel jPanel15;
    private javax.swing.JPanel jPanel16;
    private javax.swing.JPanel jPanel18;
    private javax.swing.JPanel jPanel2;
    private javax.swing.JPanel jPanel20;
    private javax.swing.JPanel jPanel21;
    private javax.swing.JPanel jPanel22;
    private javax.swing.JPanel jPanel6;
    private javax.swing.JPanel jPanel8;
    private javax.swing.JPanel jPanel9;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JScrollPane jScrollPane3;
    private javax.swing.JScrollPane jScrollPane4;
    private javax.swing.JScrollPane jScrollPane5;
    private javax.swing.JScrollPane jScrollPane6;
    private javax.swing.JTabbedPane jTabbedPane1;
    private javax.swing.JTabbedPane jTabbedPane3;
    private javax.swing.JTabbedPane jTabbedPane4;
    private javax.swing.JLabel lblDecoyDB;
    private javax.swing.JLabel lblLinkDB;
    private javax.swing.JLabel lblPeptide;
    private javax.swing.JLabel lblProtein;
    private javax.swing.JLabel lblSumBetween;
    private javax.swing.JLabel lblSumBetween1;
    private javax.swing.JLabel lblSumInternal;
    private javax.swing.JLabel lblSumInternal1;
    private javax.swing.JLabel lblTargetDB;
    private org.rappsilber.gui.components.memory.Memory memory2;
    private javax.swing.JPanel pAbout;
    private javax.swing.JPanel pDatabseSize;
    private javax.swing.JPanel pFDRGroups;
    private javax.swing.JPanel pLog;
    private javax.swing.JPanel pResult;
    private javax.swing.JPanel pVersion;
    private javax.swing.JRadioButton rbCSV;
    private javax.swing.JRadioButton rbFDRComplete;
    private javax.swing.JRadioButton rbFDRSimple;
    private javax.swing.JRadioButton rbMZHighBetter;
    private javax.swing.JRadioButton rbMZLowBetter;
    private javax.swing.JRadioButton rbTSV;
    private javax.swing.JSpinner spDecoyDB;
    private javax.swing.JSpinner spDecoyDBLinks;
    private javax.swing.JSpinner spDecoyDBProt;
    private org.rappsilber.fdr.gui.components.SingleTextValueNumericSpinner spDistanceGroup;
    private javax.swing.JPanel spExtraLong;
    private javax.swing.JScrollPane spFDRSettingsWrapper;
    private javax.swing.JScrollPane spLog;
    private javax.swing.JScrollPane spPepLength;
    private javax.swing.JSpinner spTargetDB;
    private javax.swing.JSpinner spTargetDBLinks;
    private javax.swing.JSpinner spTargetDBProt;
    private javax.swing.JTable tblPepLength;
    protected javax.swing.JTabbedPane tpInput;
    private javax.swing.JTextField txtBaseName;
    private javax.swing.JTextArea txtLog;
    private javax.swing.JTextField txtStatus;
    private javax.swing.JTextField txtSumInput;
    private javax.swing.JTextField txtSumLinks;
    private javax.swing.JTextField txtSumLinksBetween;
    private javax.swing.JTextField txtSumLinksInternal;
    private javax.swing.JTextField txtSumPSM;
    private javax.swing.JTextField txtSumPSMLinear;
    private javax.swing.JTextField txtSumPSMXL;
    private javax.swing.JTextField txtSumPepPairs;
    private javax.swing.JTextField txtSumPepPairsLinear;
    private javax.swing.JTextField txtSumPepPairsXL;
    private javax.swing.JTextField txtSumProtGroupPairs;
    private javax.swing.JTextField txtSumProtGroupPairsBetween;
    private javax.swing.JTextField txtSumProtGroupPairsInternal;
    private javax.swing.JTextField txtSumProtGroups;
    private javax.swing.JTextField txtXiFDRVersion;
    private javax.swing.JTextArea txtmzIdentAdress;
    private javax.swing.JTextField txtmzIdentOwnerEmail;
    private javax.swing.JTextField txtmzIdentOwnerFirst;
    private javax.swing.JTextField txtmzIdentOwnerLast;
    private javax.swing.JTextField txtmzIdentOwnerOrg;
    // End of variables declaration//GEN-END:variables
}
