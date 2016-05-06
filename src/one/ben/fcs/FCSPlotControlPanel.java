package one.ben.fcs;

import java.awt.Font;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.io.FilenameFilter;
import java.text.Format;
import java.text.NumberFormat;
import java.util.TreeSet;

import javax.swing.BorderFactory;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFormattedTextField;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.border.BevelBorder;
import javax.swing.table.DefaultTableModel;

import net.miginfocom.swing.MigLayout;
import one.ben.fcs.AbstractPanel2.AXIS_SCALE;
import one.ben.fcs.AbstractPanel2.PLOT_TYPE;
import cnv.gui.JAccordionPanel;

import common.Files;
import common.ext;

public class FCSPlotControlPanel extends JPanel {

    private FCSPlot plot;
    
    private JComboBox<PLOT_TYPE> cbType;
    private JComboBox<String> cbYData;
    private JComboBox<String> cbXData;
    private JComboBox<AXIS_SCALE> cbYScale;
    private JComboBox<AXIS_SCALE> cbXScale;
    
    private Font lblFont = new Font("Arial", 0, 12);
    private JFormattedTextField yBndsMin;
    private JFormattedTextField yBndsMax;
    private JFormattedTextField xBndsMin;
    private JFormattedTextField xBndsMax;
    
    volatile boolean progSet = false;

    private JCheckBox chckbxShowMedianX;
    private JCheckBox chckbxShowSdX;
    private JCheckBox chckbxShowMedianY;
    private JCheckBox chckbxShowSdY;

    private JProgressBar progressBar;

    private JAccordionPanel plotControlPanel;
    private JTextField fileDirField;
    private JTable dataFileTable;
    private JButton btnLoadFile;
    private JButton btnMoveUp;
    private JButton btnMoveDown;
    private JButton btnRemove;

    /**
     * Create the panel.
     */
    public FCSPlotControlPanel(final FCSPlot plot) {
        this.plot = plot;
        
        setLayout(new MigLayout("ins 0", "[grow]", "[grow][]"));
        
        panel_1 = new JPanel();
        panel_1.setLayout(new MigLayout("ins 0", "[grow]", "[]0[]0[]0[grow]"));
        
        plotControlPanel = new JAccordionPanel();
        panel_1.add(plotControlPanel, "cell 0 0,grow");
        JLabel ctrlLabel = new JLabel("<html><u>Plot Controls</u></html>");
        ctrlLabel.setFont(ctrlLabel.getFont().deriveFont(Font.PLAIN, 14));
        plotControlPanel.topPanel.add(ctrlLabel, "pad 0 10 0 0, cell 0 0, grow");
        plotControlPanel.topPanel.setBorder(BorderFactory.createSoftBevelBorder(BevelBorder.RAISED));
        plotControlPanel.contentPanel.setBorder(BorderFactory.createSoftBevelBorder(BevelBorder.LOWERED));
        JPanel panel = plotControlPanel.contentPanel;
        
        panel.setLayout(new MigLayout("", "[][][grow][]", "[][][][][][][][][][][]"));
        
        JLabel lblPlotType = new JLabel("Plot Type:");
        panel.add(lblPlotType, "cell 0 0,alignx trailing");
        lblPlotType.setFont(lblFont);
        
        cbType = new JComboBox<PLOT_TYPE>(PLOT_TYPE.values());
        panel.add(cbType, "cell 1 0 3 1,growx");
        
        JLabel lblYaxisData = new JLabel("Y-Axis Data:");
        panel.add(lblYaxisData, "cell 0 2,alignx trailing");
        lblYaxisData.setFont(lblFont);
        
        cbYData = new JComboBox<String>();
        panel.add(cbYData, "cell 1 2 3 1,growx");
        cbYData.addItemListener(new ItemListener() {
            @Override
            public void itemStateChanged(ItemEvent arg0) {
                if (arg0.getStateChange() == ItemEvent.SELECTED) {
                    plot.setYDataName(arg0.getItem().toString());
                    plot.updateGUI();
                }
            }
        });
        cbYData.setMaximumRowCount(15);
        
        JLabel lblScale = new JLabel("Scale:");
        panel.add(lblScale, "cell 0 3 2 1,alignx trailing");
        lblScale.setFont(lblFont);
        
        cbYScale = new JComboBox<AbstractPanel2.AXIS_SCALE>(AXIS_SCALE.values());
        panel.add(cbYScale, "cell 2 3 2 1,growx");
        cbYScale.addItemListener(new ItemListener() {
            @Override
            public void itemStateChanged(ItemEvent arg0) {
                if (arg0.getStateChange() == ItemEvent.SELECTED) {
                    plot.setYScale((AXIS_SCALE) arg0.getItem());
                    plot.updateGUI();
                }
            }
        });
        
        chckbxShowMedianY = new JCheckBox("Show Median", plot.showMedian(true));
        panel.add(chckbxShowMedianY, "cell 0 4 3 1,alignx trailing");
        chckbxShowMedianY.setHorizontalAlignment(SwingConstants.TRAILING);
        
        chckbxShowSdY = new JCheckBox("Show SD", plot.showSD(true));
        panel.add(chckbxShowSdY, "cell 3 4,alignx trailing");
        chckbxShowSdY.setHorizontalAlignment(SwingConstants.TRAILING);
        
        JLabel lblYbounds = new JLabel("Y-Bounds:");
        panel.add(lblYbounds, "cell 0 5 2 1,alignx trailing");
        lblYbounds.setFont(lblFont);

        Format numberFormat = NumberFormat.getNumberInstance();
        
        yBndsMin = new JFormattedTextField(numberFormat);
        panel.add(yBndsMin, "cell 2 5 2 1");
        yBndsMin.addPropertyChangeListener("value", pcl);
        yBndsMin.setColumns(10);
        yBndsMin.setValue(0);
        yBndsMin.setEditable(false);

        yBndsMax = new JFormattedTextField(numberFormat);
        panel.add(yBndsMax, "cell 2 5 2 1");
        yBndsMax.addPropertyChangeListener("value", pcl);
        yBndsMax.setColumns(10);
        yBndsMax.setEditable(false);
        
        JLabel lblXaxisData = new JLabel("X-Axis Data:");
        panel.add(lblXaxisData, "cell 0 7,alignx trailing");
        lblXaxisData.setFont(lblFont);
        
        cbXData = new JComboBox<String>();
        panel.add(cbXData, "cell 1 7 3 1,growx");
        cbXData.addItemListener(new ItemListener() {
            @Override
            public void itemStateChanged(ItemEvent arg0) {
                if (arg0.getStateChange() == ItemEvent.SELECTED) {
                    plot.setXDataName(arg0.getItem().toString());
                    plot.updateGUI();
                }
            }
        });
        cbXData.setMaximumRowCount(15);
        
        JLabel lblScale_1 = new JLabel("Scale:");
        panel.add(lblScale_1, "cell 0 8 2 1,alignx trailing");
        lblScale_1.setFont(lblFont);
        
        cbXScale = new JComboBox<AbstractPanel2.AXIS_SCALE>(AXIS_SCALE.values());
        panel.add(cbXScale, "cell 2 8 2 1,growx");
        
        chckbxShowMedianX = new JCheckBox("Show Median", plot.showMedian(false));
        panel.add(chckbxShowMedianX, "cell 0 9 3 1,alignx right");
        chckbxShowMedianX.setHorizontalAlignment(SwingConstants.TRAILING);
        
        chckbxShowSdX = new JCheckBox("Show SD", plot.showSD(false));
        panel.add(chckbxShowSdX, "cell 3 9,alignx right");
        chckbxShowSdX.setHorizontalAlignment(SwingConstants.TRAILING);
        
        JLabel lblXbounds = new JLabel("X-Bounds:");
        panel.add(lblXbounds, "cell 0 10 2 1,alignx trailing");
        lblXbounds.setFont(lblFont);
        
        xBndsMin = new JFormattedTextField(numberFormat);
        panel.add(xBndsMin, "cell 2 10 2 1");
        xBndsMin.addPropertyChangeListener("value", pcl);
        xBndsMin.setColumns(10);
        xBndsMin.setValue(0);
        xBndsMin.setEditable(false);

        xBndsMax = new JFormattedTextField(numberFormat);
        panel.add(xBndsMax, "cell 2 10 2 1");
        xBndsMax.addPropertyChangeListener("value", pcl);
        xBndsMax.setColumns(10);
        xBndsMax.setEditable(false);
        
        chckbxShowSdX.addItemListener(new ItemListener() {
            @Override
            public void itemStateChanged(ItemEvent arg0) {
                boolean show = arg0.getStateChange() == ItemEvent.SELECTED;
                plot.setSDVisible(show, false);
                plot.updateGUI();
            }
        });
        chckbxShowMedianX.addItemListener(new ItemListener() {
            @Override
            public void itemStateChanged(ItemEvent arg0) {
                boolean show = arg0.getStateChange() == ItemEvent.SELECTED;
                plot.setMedianVisible(show, false);
                plot.updateGUI();
            }
        });
        cbXScale.addItemListener(new ItemListener() {
            @Override
            public void itemStateChanged(ItemEvent arg0) {
                if (arg0.getStateChange() == ItemEvent.SELECTED) {
                    plot.setXScale((AXIS_SCALE) arg0.getItem());
                    plot.updateGUI();
                }
            }
        });
        chckbxShowSdY.addItemListener(new ItemListener() {
            @Override
            public void itemStateChanged(ItemEvent arg0) {
                boolean show = arg0.getStateChange() == ItemEvent.SELECTED;
                plot.setSDVisible(show, true);
                plot.updateGUI();
            }
        });
        chckbxShowMedianY.addItemListener(new ItemListener() {
            @Override
            public void itemStateChanged(ItemEvent arg0) {
                boolean show = arg0.getStateChange() == ItemEvent.SELECTED;
                plot.setMedianVisible(show, true);
                plot.updateGUI();
            }
        });
        cbType.addItemListener(new ItemListener() {
            @Override
            public void itemStateChanged(ItemEvent arg0) {
                if (arg0.getStateChange() == ItemEvent.SELECTED) {
                    plot.setPlotType((PLOT_TYPE) arg0.getItem());
                    plot.updateGUI();
                }
            }
        });
        
        gateControlPanel = new JAccordionPanel();
        panel_1.add(gateControlPanel, "cell 0 1,grow");
        gateControlPanel.shrink();
        gateControlPanel.contentPanel.setBorder(BorderFactory.createSoftBevelBorder(BevelBorder.LOWERED));
        gateControlPanel.topPanel.setBorder(BorderFactory.createSoftBevelBorder(BevelBorder.RAISED));
        JLabel gCtrlLabel = new JLabel("<html><u>Gate Controls</u></html>");
        gCtrlLabel.setFont(gCtrlLabel.getFont().deriveFont(Font.PLAIN, 14));
        gateControlPanel.topPanel.add(gCtrlLabel, "pad 0 10 0 0, cell 0 0, grow");
        
        dataControlsPanel = new JAccordionPanel();
        panel_1.add(dataControlsPanel, "cell 0 2,grow");
        dataControlsPanel.shrink();
        dataControlsPanel.contentPanel.setBorder(BorderFactory.createSoftBevelBorder(BevelBorder.LOWERED));
        dataControlsPanel.topPanel.setBorder(BorderFactory.createSoftBevelBorder(BevelBorder.RAISED));
        JLabel dCtrlLabel = new JLabel("<html><u>Data Controls</u></html>");
        dCtrlLabel.setFont(dCtrlLabel.getFont().deriveFont(Font.PLAIN, 14));
        dataControlsPanel.topPanel.add(dCtrlLabel, "pad 0 10 0 0, cell 0 0, grow");
        
        JPanel dataPanel = dataControlsPanel.contentPanel;
        dataPanel.setLayout(new MigLayout("", "[][grow][]", "[][grow][][]"));
        JLabel lblFileDirectory = new JLabel("FCS Dir:");
        dataPanel.add(lblFileDirectory, "cell 0 0,alignx trailing");
        
        fileDirField = new JTextField();
        dataPanel.add(fileDirField, "cell 1 0,growx");
        fileDirField.setColumns(10);
        
        dirSelectBtn = new JButton(">");
        dirSelectBtn.setMargin(new Insets(0, 0, 0, 0));
        dirSelectBtn.addActionListener(dirSelectListener);
        dataPanel.add(dirSelectBtn, "cell 2 0");
        
        dataFileTable = new JTable();
        dataPanel.add(dataFileTable, "cell 0 1 3 1,grow");
        
        JPanel dataBtnPanel = new JPanel();
        dataBtnPanel.setBorder(null);
        dataPanel.add(dataBtnPanel, "cell 0 2 3 1,grow");
        dataBtnPanel.setLayout(new MigLayout("insets 0", "[60px,grow][60px,grow]", "[][][]"));
        
        btnLoadFile = new JButton("Load");
        dataBtnPanel.add(btnLoadFile, "cell 0 0 2 1,alignx center,growx");
        
        btnMoveUp = new JButton("Move Up");
        dataBtnPanel.add(btnMoveUp, "cell 0 1,growx");
        
        btnMoveDown = new JButton("Move Down");
        dataBtnPanel.add(btnMoveDown, "cell 1 1,growx");
        
        btnRemove = new JButton("Remove");
        dataBtnPanel.add(btnRemove, "cell 0 2 2 1,alignx center,growx");
        
        add(panel_1, "cell 0 0,grow");
        
        progressBar = new JProgressBar();
        add(progressBar, "cell 0 1,growx, pad -3 3 -3 -3");
        
    }
    
    private void listFiles() {
        String fileDir = ext.verifyDirFormat(fileDirField.getText());
        String[] files = new File(fileDir).list(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                return name.endsWith(".fcs");
            }
        });
        TreeSet<String> fileSet = new TreeSet<String>();
        for (String f : files) {
            fileSet.add(f);
        }
        DefaultTableModel dtm = new DefaultTableModel(new String[]{"File", "Size", "Loaded?"}, 0);
        for (String f : fileSet) {
            dtm.addRow(new Object[]{f, Files.getSizeScaledString(fileDir + f, false), false});
        }
        dataFileTable.setModel(dtm);
        dataFileTable.invalidate();
        repaint();
    }
    
    ActionListener dirSelectListener = new ActionListener() {
        @Override
        public void actionPerformed(ActionEvent e) {
            String curr = fileDirField.getText();
            JFileChooser jfc = new JFileChooser(curr);
            jfc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
            jfc.setDialogTitle("Select FCS File Directory");
            jfc.setMultiSelectionEnabled(false);
            int resp = jfc.showOpenDialog(FCSPlotControlPanel.this);
            if (resp == JFileChooser.APPROVE_OPTION) {
                String newPath = jfc.getSelectedFile().getAbsolutePath();
                fileDirField.setText(newPath);
                listFiles();
            }
        }
    };
    
    PropertyChangeListener pcl = new PropertyChangeListener() {
        @Override
        public void propertyChange(PropertyChangeEvent evt) {
            if (evt == null || plot == null || progSet) return;
            JFormattedTextField jftf = (JFormattedTextField) evt.getSource();
            String prop = evt.getPropertyName();
            if (jftf == xBndsMin) {
                prop = AbstractPanel2.X_MIN;
            } else if (jftf == xBndsMax) {
                prop = AbstractPanel2.X_MAX;
            } else if (jftf == yBndsMin) {
                prop = AbstractPanel2.Y_MIN;
            } else if (jftf == yBndsMax) {
                prop = AbstractPanel2.Y_MAX;
            }
            Object oldV = evt.getOldValue();
            Object newV = evt.getNewValue();
            Double oldV2 = oldV == null ? null : ((Number) oldV).doubleValue(); 
            Double newV2 = newV == null ? null : ((Number) newV).doubleValue();
            /*FCSPlotControlPanel.this.*/firePropertyChange(prop, oldV2, newV2);
//            plot.firePropertyChange(prop, oldV2, newV2);
//            PropertyChangeEvent pce = new PropertyChangeEvent(FCSPlotControlPanel.this, prop, oldV2, newV2);
//            plot.propertyChange(pce);
            
        }
    };
    private JAccordionPanel dataControlsPanel;
    private JAccordionPanel gateControlPanel;
    private JPanel panel_1;

    private JButton dirSelectBtn;
    
    public void setPlotType(PLOT_TYPE typ) {
        cbType.setSelectedItem(typ);
    }
    
    public void setScale(AXIS_SCALE scl, boolean x) {
        (x ? cbXScale : cbYScale).setSelectedItem(scl);
    }
    
    public void setColumns(String[] dataNames, boolean x, int selected) {
        (x ? cbXData : cbYData).setModel(new DefaultComboBoxModel<String>(dataNames));
        (x ? cbXData : cbYData).setSelectedIndex(selected);
        (x ? cbXData : cbYData).repaint();
    }
    
    private void resetProgSet() {
        SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
                progSet = false;
            }
        });
    }
    
    public void setXMin(double xMin) {
        progSet = true;
        xBndsMin.setValue(Math.floor(xMin));
        resetProgSet();
    }
    public void setXMax(double xMax) {
        progSet = true;
        xBndsMax.setValue(Math.ceil(xMax));
        resetProgSet();
    }
    public void setYMin(double yMin) {
        progSet = true;
        yBndsMin.setValue(Math.floor(yMin));
        resetProgSet();
    }
    public void setYMax(double yMax) {
        progSet = true;
        yBndsMax.setValue(Math.ceil(yMax));
        resetProgSet();
    }

    // TODO doesn't work very well - no updates until data gets displayed
//    public void startFileLoading(FCSDataLoader newDataLoader) {
//        new Thread(new Runnable() {
//            @Override
//            public void run() {
//                LOAD_STATE state = null;
//                while((state = newDataLoader.getLoadState()) != LOAD_STATE.LOADED) {
//                    final LOAD_STATE finalState = state;
//                    SwingUtilities.invokeLater(new Runnable() {
//                        @Override
//                        public void run() {
//                            switch(finalState) {
//                            case LOADED:
//                                progressBar.setStringPainted(false);
//                                progressBar.setString(null);
//                                progressBar.setIndeterminate(false);
//                                // hide or set to complete or reset
//                                break;
//                            case LOADING:
//                                progressBar.setIndeterminate(true);
//                                progressBar.setStringPainted(true);
//                                progressBar.setString("Starting File Load...");
//                                progressBar.setVisible(true);
//                                // set to indeterminate
//                                break;
//                            case PARTIALLY_LOADED:
//                            case LOADING_REMAINDER:
//                                progressBar.setIndeterminate(false);
//                                progressBar.setStringPainted(true);
//                                int[] stat = newDataLoader.getLoadedStatus();
//                                progressBar.setMinimum(0);
//                                progressBar.setMaximum(stat[1]);
//                                progressBar.setString(null);
////                            progressBar.setString("Loading File: " + stat[0] + "/" + stat[1]);
//                                progressBar.setVisible(true);
//                                // set to determinate, wait for updates
//                                break;
//                            case UNLOADED:
//                                // what??
//                                break;
//                            default:
//                                // what??
//                                break;
//                            }
//                        }
//                    });
//                    Thread.yield();
//                }
//                progressBar.setStringPainted(false);
//                progressBar.setString(null);
//                progressBar.setIndeterminate(false);
//            }
//        }).start();
//    }
    
}
