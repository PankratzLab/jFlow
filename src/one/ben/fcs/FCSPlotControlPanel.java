package one.ben.fcs;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Insets;
import java.awt.LayoutManager;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.io.FilenameFilter;
import java.text.Format;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.TreeSet;

import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFormattedTextField;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.ScrollPaneConstants;
import javax.swing.Scrollable;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.border.BevelBorder;
import javax.swing.table.DefaultTableModel;

import net.miginfocom.swing.MigLayout;
import one.ben.fcs.AbstractPanel2.AXIS_SCALE;
import one.ben.fcs.AbstractPanel2.PLOT_TYPE;
import one.ben.fcs.FCSDataLoader.LOAD_STATE;
import scala.collection.mutable.HashMap;
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

    public JProgressBar progressBar;

    private JAccordionPanel plotControlPanel;
    private JAccordionPanel dataControlsPanel;
    private JAccordionPanel gateControlPanel;
    private JPanel panel_1;

    private JButton dirSelectBtn;
    private JTextField fileDirField;

    /**
     * Create the panel.
     */
    public FCSPlotControlPanel(final FCSPlot plot) {
        this.plot = plot;
        
        setLayout(new MigLayout("ins 0", "[grow]", "[grow][]"));
        
        panel_1 = new JPanel();
        panel_1.setLayout(new MigLayout("ins 0", "[grow]", "[]0[]0[]0[grow]"));
        
        ArrayList<JAccordionPanel> bg = new ArrayList<JAccordionPanel>();
         
        plotControlPanel = new JAccordionPanel();
        panel_1.add(plotControlPanel, "cell 0 0,grow");
        JLabel ctrlLabel = new JLabel("<html><u>Plot Controls</u></html>");
        ctrlLabel.setFont(ctrlLabel.getFont().deriveFont(Font.PLAIN, 14));
        plotControlPanel.topPanel.add(ctrlLabel, "pad 0 10 0 0, cell 0 0, grow");
        plotControlPanel.topPanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED));
        plotControlPanel.contentPanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.LOWERED));
        plotControlPanel.addToGroup(bg);
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
        gateControlPanel.contentPanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.LOWERED));
        gateControlPanel.topPanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED));
        gateControlPanel.addToGroup(bg);
        JLabel gCtrlLabel = new JLabel("<html><u>Gate Controls</u></html>");
        gCtrlLabel.setFont(gCtrlLabel.getFont().deriveFont(Font.PLAIN, 14));
        gateControlPanel.topPanel.add(gCtrlLabel, "pad 0 10 0 0, cell 0 0, grow");
        
        dataControlsPanel = new JAccordionPanel();
        panel_1.add(dataControlsPanel, "cell 0 2,grow");
        dataControlsPanel.shrink();
        dataControlsPanel.contentPanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.LOWERED));
        dataControlsPanel.topPanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED));
        dataControlsPanel.addToGroup(bg);
        JLabel dCtrlLabel = new JLabel("<html><u>Data Controls</u></html>");
        dCtrlLabel.setFont(dCtrlLabel.getFont().deriveFont(Font.PLAIN, 14));
        dataControlsPanel.topPanel.add(dCtrlLabel, "pad 0 10 0 0, cell 0 0, grow");
        
        JPanel dataPanel = dataControlsPanel.contentPanel;
        dataPanel.setLayout(new MigLayout("hidemode 3,ins 0", "[grow]", "[grow]0px[]"));
        
        fileDirField = new JTextField();
        dataPanel.add(fileDirField, "cell 0 0, growx, split 2");
        dirSelectBtn = new JButton(">");
        dirSelectBtn.addActionListener(dirSelectListener);
        dirSelectBtn.setMargin(new Insets(0, 2, 0, 2));
        dataPanel.add(dirSelectBtn, "cell 0 0");
        
        final JScrollPane scrollPane = new JScrollPane();
        actualDataPanel = new ScrollablePanel(new MigLayout("", "", ""));
        scrollPane.setViewportView(actualDataPanel);
        scrollPane.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
        scrollPane.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS);
        dataPanel.add(scrollPane, "cell 0 1, grow");
        
        // TODO mouse keys (up/down for changing selection, alt/ctrl-up/down for displaying prev/next file data, etc)
        // TODO listener for move up
        // TODO listener for move down
        // TODO listener for information
        // TODO listener for delete
        // TODO listener for load
        // TODO memory management and warnings
        
//        {
//            int ind = 0;
//            int test = 19;
//            HashMap<String, DataControlPanel> fileCtrl = new HashMap<String, DataControlPanel>();
//            final ArrayList<DataControlPanel> ctrlList = new ArrayList<DataControlPanel>();
//            for (int i = 0; i < test; i++) {
//                final int index = i;
//                DataControlPanel dcp = new DataControlPanel();
//                dcp.addMouseListener(new MouseAdapter() {
//                    @Override
//                    public void mouseClicked(MouseEvent e) {
//                        super.mouseClicked(e);
//                        for (int j = 0; j < ctrlList.size(); j++) {
//                            ctrlList.get(j).setSelected(j == index);
//                        }
//                    }
//                });
//                fileCtrl.put(i + "", dcp);
//                ctrlList.add(dcp);
//                actualDataPanel.add(dcp, "cell 0 " + (ind++));
//                actualDataPanel.add(new JSeparator(JSeparator.HORIZONTAL), "grow, cell 0 " + (ind++));
//            }
//        }
        
        dataMsgPanel = new JPanel();
        dataPanel.add(dataMsgPanel, "cell 0 1, grow");
        dataMsgPanel.setVisible(false);
        
        add(panel_1, "cell 0 0,grow");
        
        progressBar = new JProgressBar();
        add(progressBar, "cell 0 1,growx, pad -3 3 -3 -3");
        
    }
    
    public class ScrollablePanel extends JPanel implements Scrollable {
        
        public ScrollablePanel(LayoutManager layout) {
            super(layout);
        }
        
        public Dimension getPreferredScrollableViewportSize() {
            return getPreferredSize();
        }

        public int getScrollableUnitIncrement(Rectangle visibleRect, int orientation, int direction) {
           return 10;
        }

        public int getScrollableBlockIncrement(Rectangle visibleRect, int orientation, int direction) {
            return ((orientation == SwingConstants.VERTICAL) ? visibleRect.height : visibleRect.width) - 10;
        }

        public boolean getScrollableTracksViewportWidth() {
            return true;
        }

        public boolean getScrollableTracksViewportHeight() {
            return false;
        }
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
        
        int index = 0;
        actualDataPanel.removeAll();
        actualDataPanel.add(new JSeparator(JSeparator.HORIZONTAL), "grow, cell 0 " + (index++));
        for (String f : fileSet) {
            String sz = Files.getSizeScaledString(fileDir + f, false);
            String dt = "";
            DataControlPanel dcp = new DataControlPanel(f, sz, dt);
            actualDataPanel.add(dcp, "cell 0 " + (index++));
            actualDataPanel.add(new JSeparator(JSeparator.HORIZONTAL), "grow, cell 0 " + (index++));
        }
        
        revalidate();
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

    private JPanel dataMsgPanel;

    private JPanel actualDataPanel;
    
    public void setPlotType(PLOT_TYPE typ) {
        cbType.setSelectedItem(typ);
    }
    
    public PLOT_TYPE getPlotType() {
        return (PLOT_TYPE) cbType.getSelectedItem();
    }
    
    public void setScale(AXIS_SCALE scl, boolean x) {
        (x ? cbXScale : cbYScale).setSelectedItem(scl);
    }
    
    public void setColumns(String[] dataNames, boolean x, int selected) {
        (x ? cbXData : cbYData).setModel(new DefaultComboBoxModel<String>(dataNames));
        (x ? cbXData : cbYData).setSelectedIndex(selected);
        (x ? cbXData : cbYData).repaint();
    }
    
    public String getSelectedX() { return (String) cbXData.getSelectedItem(); }
    public String getSelectedY() { return (String) cbYData.getSelectedItem(); }
    
    public void setYData(int index) {
        cbYData.setSelectedIndex(index);
        cbYData.repaint();
    }
    
    public void setYData(String name) {
        cbYData.setSelectedItem(name);
        cbYData.repaint();
    }
    
    public void setXData(String name) {
        cbXData.setSelectedItem(name);
        cbXData.repaint();
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

//    public void startFileLoading(FCSDataLoader newDataLoader) {
//        // TODO Auto-generated method stub
//        
//    }

    // TODO doesn't work very well - no updates until data gets displayed
    public void startFileLoading(FCSDataLoader newDataLoader) {
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
    }
    
}
