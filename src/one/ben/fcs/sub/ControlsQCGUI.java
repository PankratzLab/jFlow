package one.ben.fcs.sub;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.EventQueue;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.File;
import java.io.FileFilter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FilenameFilter;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.Properties;
import java.util.TreeMap;
import java.util.TreeSet;

import javax.swing.DefaultComboBoxModel;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JSpinner;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.ScrollPaneConstants;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingConstants;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.table.DefaultTableColumnModel;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumnModel;

import net.miginfocom.swing.MigLayout;
import one.ben.fcs.AbstractPanel2;
import one.ben.fcs.FCSPlot;
import one.ben.fcs.sub.MeanCtrlPanel.LabelPresenter;
import one.ben.fcs.sub.OneDPanel.PLOT_TYPE;
import cnv.gui.IncludeExcludeGUI;

import common.Array;
import common.HashVec;
import common.Matrix;
import common.ext;


public class ControlsQCGUI extends JFrame {

    private JPanel contentPane;
    private JTextField txtFldBaseFile;
    private JButton btnBaseFileSelect;
    private JTable table;
    private JScrollPane scrollPane;

    HashMap<String, ArrayList<Float>> paramMeanLists = new HashMap<String, ArrayList<Float>>();
//    HashMap<String, HashMap<String, Float>> fileParamMeanMap = new HashMap<String, HashMap<String,Float>>();
    HashSet<Integer> boldRows = new HashSet<Integer>();
    HashSet<Integer> statRows = new HashSet<Integer>();
    HashMap<String, Float> paramMeans = new HashMap<String, Float>();
    HashMap<String, Float> paramSDs = new HashMap<String, Float>();
    HashMap<String, Float> paramCVs = new HashMap<String, Float>();

    HashSet<String> hiddenCols = new HashSet<String>();
    
    JFrame meanFrame = new JFrame("Genvisis - Controls QC - Overall Mean/SD");
    OneDPanel meanPanel = new OneDPanel();
    MeanCtrlPanel meanCtrlPanel = new MeanCtrlPanel();
    private JLabel lblCompareFiles;
    private JTextField txtFldCompDir;
    private JButton btnCompDirSelect;
    private DefaultTableModel dtmMean;
    private JButton btnHideshowColumns;
    Color SD1_COLOR = Color.YELLOW;
    Color SD2_COLOR = Color.RED;
    Color ABOVE_3_CV_COLOR = Color.CYAN;
    Color ABOVE_5_CV_COLOR = Color.CYAN.darker();

    CompDir compDirData;
    DataFile baseData;
    private final static String[] EXCLUDED_ROW_HEADERS = {
            "mean",
            "average",
            "sd",
            "cv",
            "cv (%)"
    };
    private JPanel panel_1;
    private JComboBox<String> comboControl;
    private JSeparator separator;
    private JLabel lblControlGroup;
    private JSeparator separator_1;
    private JLabel lblPanel;
    private JComboBox<String> comboPanel;
    private JSeparator separator_2;

    protected void loadFileOrDir(String fileOrDir, boolean baseline) {
        if (baseline) {
            if (baseData != null && baseData.fileName.equals(fileOrDir)) {
                return; // same as already loaded
            }
            baseData = DataFile.construct(fileOrDir);
            String[] pnls = Array.addStrToArray("All", baseData.getPanels(), 0);
            String[] ctrls = Array.addStrToArray("All", baseData.getControlGroups(), 0);
            
            String currPanel, currCtrl;
            currCtrl = (String) comboControl.getSelectedItem();
            DefaultComboBoxModel<String> ctrlModel = new DefaultComboBoxModel<String>(ctrls);
            int ind = ext.indexOfStr(currCtrl, ctrls);
            comboControl.setModel(ctrlModel);
            comboControl.setSelectedIndex(ind == -1 ? 0 : ind);
            
            currPanel = (String) comboPanel.getSelectedItem();
            DefaultComboBoxModel<String> pnlModel = new DefaultComboBoxModel<String>(pnls);
            ind = ext.indexOfStr(currPanel, pnls);
            comboPanel.setModel(pnlModel);
            comboPanel.setSelectedIndex(ind == -1 ? 0 : ind);
        } else {
            if (compDirData != null && compDirData.dir.equals(ext.verifyDirFormat(fileOrDir))) {
                return; // same as already loaded
            }
            compDirData = CompDir.construct(fileOrDir);
        }
        reCalcTableData();
    }



    /**
     * Create the frame.
     */
    public ControlsQCGUI() {
        super("Controls QC");
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setBounds(50, 50, 850, 400);
        contentPane = new JPanel();
        contentPane.setBorder(null);
        setContentPane(contentPane);
        contentPane.setLayout(new MigLayout("ins 7 7 3 7,hidemode 3", "[][grow][]", "[][][grow][]"));
        
        JLabel lblFileDir = new JLabel("<html><u>B</u>aseline File:</html>");
        contentPane.add(lblFileDir, "cell 0 0,alignx trailing");
        
        txtFldBaseFile = new JTextField();
        
        txtFldBaseFile.setEditable(false);
        contentPane.add(txtFldBaseFile, "cell 1 0,growx");
        txtFldBaseFile.setColumns(10);
        
        Insets btnInsets = new Insets(0, 14, 0, 14);
        
        btnBaseFileSelect = new JButton(">");
        btnBaseFileSelect.setMargin(btnInsets);
        btnBaseFileSelect.setMnemonic(KeyEvent.VK_B);
        btnBaseFileSelect.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent arg0) {
                String curr = txtFldBaseFile.getText();
                if (curr.equals("")) {
                    curr = "./";
                }
                JFileChooser jfc = new JFileChooser(curr);
                jfc.setFileSelectionMode(JFileChooser.FILES_ONLY);
                jfc.setDialogTitle("Select Baseline File");
                jfc.setMultiSelectionEnabled(false);
                jfc.setFileFilter(new FileNameExtensionFilter("CSV Files", "csv"));
                int resp = jfc.showOpenDialog(ControlsQCGUI.this);
                if (resp == JFileChooser.APPROVE_OPTION) {
                    String newPath = jfc.getSelectedFile().getAbsolutePath();
                    txtFldBaseFile.setText(newPath);
                    loadFileOrDir(newPath, true);
                    saveProps();
                }
            }
        });
        contentPane.add(btnBaseFileSelect, "cell 2 0");
        
        lblCompareFiles = new JLabel("<html>Compare <u>D</u>ir/Files:</html>");
        contentPane.add(lblCompareFiles, "cell 0 1,alignx trailing");
        
        txtFldCompDir = new JTextField();
        txtFldCompDir.setEditable(false);
        contentPane.add(txtFldCompDir, "cell 1 1,growx");
        txtFldCompDir.setColumns(10);
        
        btnCompDirSelect = new JButton(">");
        btnCompDirSelect.setMargin(new Insets(0, 14, 0, 14));
        btnCompDirSelect.setMnemonic(KeyEvent.VK_D);
        btnCompDirSelect.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent arg0) {
                String curr = txtFldCompDir.getText();
                if (curr.equals("")) {
                    curr = "./";
                }
                JFileChooser jfc = new JFileChooser(curr);
                jfc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
                jfc.setDialogTitle("Select Compare File Directory");
                jfc.setMultiSelectionEnabled(true);
                int resp = jfc.showOpenDialog(ControlsQCGUI.this);
                if (resp == JFileChooser.APPROVE_OPTION) {
                    String newPath = ext.verifyDirFormat(jfc.getSelectedFile().getAbsolutePath());
                    txtFldCompDir.setText(newPath);
                    loadFileOrDir(newPath, false);
                    saveProps();
                }
            }
        });
        contentPane.add(btnCompDirSelect, "cell 2 1");
        
        scrollPane = new JScrollPane();
        scrollPane.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_ALWAYS);
        scrollPane.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS);
        contentPane.add(scrollPane, "cell 0 2 3 1,grow");
                
        table = new JTable() {
            
            @Override
            public Component prepareRenderer(TableCellRenderer renderer, int row, int column) {
                Component c = super.prepareRenderer(renderer, row, column);
                
                if (column == 0 && boldRows.contains(row)) {
                    c.setFont(c.getFont().deriveFont(Font.BOLD));
                }
                if (isCellSelected(row, column)) return c;
                
                Color col = Color.WHITE;
                if (column > 0 && !statRows.contains(row)) {
                    Object val = table.getModel().getValueAt(table.convertRowIndexToModel(row), table.convertColumnIndexToModel(column));
                    if (val != null) {
                        if (val instanceof Float) {
                            String colNm = table.getModel().getColumnName(table.convertColumnIndexToModel(column));
                            Float value = (Float) val;
                            if (paramMeans.containsKey(colNm)) {
                                if (value > (paramMeans.get(colNm) + paramSDs.get(colNm)) || value < (paramMeans.get(colNm) - paramSDs.get(colNm))) {
                                    col = SD1_COLOR;
                                }  
                                if (value > (paramMeans.get(colNm) + 2*paramSDs.get(colNm)) || value < (paramMeans.get(colNm) - 2*paramSDs.get(colNm))) {
                                    col = SD2_COLOR;
                                }
                            }
                        }
                    } 
                }
                c.setBackground(col);
                
                return c;
            }
        };

        table.addMouseListener(new MouseAdapter() {
            public void mouseClicked(MouseEvent e) {
                if (e.getClickCount() == 2) {
                    JTable target = (JTable) e.getSource();
                    int row = target.getSelectedRow();
                    int column = target.getSelectedColumn();
                    if (column == 0) return;
                    if (row != 0 && boldRows.contains(row)) return;
                    Object o1 = target.getValueAt(row, column);
                    if (o1 == null) return;
                    
                    final String col = (String) target.getValueAt(0, column);
                    if (row == 0) {
                        ControlsQCGUI.this.showMeanPanel(col);
                    }
                }
            }
        });
        table.getTableHeader().addMouseListener(new MouseAdapter() {
            @Override
            public void mouseClicked(MouseEvent e) {
                if (e.getClickCount() == 2) {
                    int colInd = table.columnAtPoint(e.getPoint());
                    String col = (String) table.getValueAt(0, colInd);
                    ControlsQCGUI.this.showMeanPanel(col);
                }
            }
        });
        table.setShowVerticalLines(false);
        table.setShowHorizontalLines(false);
        scrollPane.setViewportView(table);
//        table.setCellSelectionEnabled(false);
        table.setColumnModel(new DefaultTableColumnModel() {
            @Override
            public void moveColumn(int columnIndex, int newIndex) {
                if (columnIndex == 0 || newIndex == 0) return; 
                super.moveColumn(columnIndex, newIndex);
            }
        });
        table.setRowSelectionAllowed(true);
        table.setColumnSelectionAllowed(true);
        
        table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        
        panel_1 = new JPanel();
        panel_1.setBorder(null);
        contentPane.add(panel_1, "cell 0 3 3 1,grow");
        panel_1.setLayout(new MigLayout("ins 0", "[][][][][][][][][][][][][][grow]", "[]"));
        
        btnHideshowColumns = new JButton("Hide/Show Columns");
        panel_1.add(btnHideshowColumns, "cell 0 0");
        
        separator = new JSeparator();
        separator.setOrientation(SwingConstants.VERTICAL);
        panel_1.add(separator, "cell 1 0,growy");
        
        lblControlGroup = new JLabel("Control:");
        panel_1.add(lblControlGroup, "cell 2 0,alignx trailing");
        
        ActionListener reCalcListener = new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent arg0) {
                reCalcTableData();
            }
        };
        
        comboControl = new JComboBox<String>();
        comboControl.addActionListener(reCalcListener);
        panel_1.add(comboControl, "cell 3 0,growx");
        
        separator_1 = new JSeparator();
        separator_1.setOrientation(SwingConstants.VERTICAL);
        panel_1.add(separator_1, "cell 4 0,growy");
        
        lblPanel = new JLabel("Panel:");
        panel_1.add(lblPanel, "cell 5 0,alignx trailing");
        
        comboPanel = new JComboBox<String>();
        comboPanel.addActionListener(reCalcListener);
        panel_1.add(comboPanel, "cell 6 0,growx");
        
        separator_2 = new JSeparator();
        separator_2.setOrientation(SwingConstants.VERTICAL);
        panel_1.add(separator_2, "cell 7 0,growy");
        
        chckbxShowInternalControl = new JCheckBox();
        chckbxShowInternalControl.addActionListener(reCalcListener);
        chckbxShowInternalControl.setText("Show Internal Control Data: ");
        chckbxShowInternalControl.setHorizontalTextPosition(SwingConstants.LEFT);
        chckbxShowInternalControl.setSelected(true);
        panel_1.add(chckbxShowInternalControl, "cell 8 0");
        
        separator_3 = new JSeparator();
        separator_3.setOrientation(SwingConstants.VERTICAL);
        panel_1.add(separator_3, "cell 9 0,growy");
        
        lblOfControls = new JLabel("# of Controls (0 = all):");
        panel_1.add(lblOfControls, "cell 10 0");
        
        spinner = new JSpinner();
        spinner.setModel(new SpinnerNumberModel(new Integer(0), new Integer(0), null, new Integer(1)));
        panel_1.add(spinner, "cell 11 0");
        
        separator_4 = new JSeparator();
        separator_4.setOrientation(SwingConstants.VERTICAL);
        panel_1.add(separator_4, "cell 12 0,growy");
        btnHideshowColumns.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent arg0) {
                if (baseData == null) return;
                String[] opts = new String[baseData.getAllParams().size()];
                boolean[] incl = new boolean[baseData.getAllParams().size()];
                int ind = 0;
                for (String p : baseData.getAllParams()) {
                    opts[ind] = p;
                    incl[ind] = !hiddenCols.contains(p);
                    ind++;
                }
                IncludeExcludeGUI dialog = new IncludeExcludeGUI(ControlsQCGUI.this, opts, incl);
                dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
                dialog.pack();
                dialog.setVisible(true);
                int code = dialog.getCloseCode();
                if (code == JOptionPane.OK_OPTION) {
                    boolean[] inc = dialog.getSelected();
                    hiddenCols.clear();
                    for (int i = 0; i < opts.length; i++) {
                        if (!inc[i]) {
                            hiddenCols.add(opts[i]);
                        }
                    }
                    saveProps();
                    reCalcTableData();
                }
            }
        });

        meanPanel.setOpaque(true);
        meanFrame.getContentPane().add(meanPanel, BorderLayout.CENTER);
        meanFrame.getContentPane().add(meanCtrlPanel, BorderLayout.SOUTH);
        ActionListener prevLst = new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent arg0) {
                String newCol = arg0.getActionCommand();
                showMeanPanel(newCol);
            }
        };
        ActionListener plotLst = new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent arg0) {
                String newCol = arg0.getActionCommand();
                if ("DOT".equalsIgnoreCase(newCol)) {
                    meanPanel.setOpaque(true);
                    meanPanel.setPlotType(OneDPanel.PLOT_TYPE.DOT_LINE_PLOT);
                    meanPanel.setAxisXHeight(AbstractPanel2.HEIGHT_X_AXIS - AbstractPanel2.HEIGHT_X_AXIS/5);
                    meanPanel.setXAxisLabel("File by Date");
                } else if ("BOX".equalsIgnoreCase(newCol)) {
                    meanPanel.setPlotType(PLOT_TYPE.BOX_PLOT);
                    meanPanel.setAxisXHeight(AbstractPanel2.HEIGHT_X_AXIS - AbstractPanel2.HEIGHT_X_AXIS / 2);
                    meanPanel.setAxisYWidth(AbstractPanel2.WIDTH_Y_AXIS - AbstractPanel2.WIDTH_Y_AXIS / 3);
                    meanPanel.setXAxisLabel("");// pts[0].trim().replaceAll("/", " /\n");
                    meanPanel.setYAxisLabel(meanPanel.plotLabel.split("\\|")[1].trim());
                }
                meanPanel.invalidate();
                meanPanel.repaint();
//                showMeanPanel(newCol);
            }
        };
        meanCtrlPanel.setChangeListener(prevLst);
        meanCtrlPanel.setPlotChangeListener(plotLst);
        meanPanel.setPlotType(PLOT_TYPE.BOX_PLOT);
        meanPanel.setAxisXHeight(AbstractPanel2.HEIGHT_X_AXIS - AbstractPanel2.HEIGHT_X_AXIS / 2);
        meanPanel.setAxisYWidth(AbstractPanel2.WIDTH_Y_AXIS - AbstractPanel2.WIDTH_Y_AXIS / 3);
        meanPanel.setXAxisLabel("");// pts[0].trim().replaceAll("/", " /\n");
//        meanPanel.setYAxisLabel(meanPanel.plotLabel.split("\\|")[1].trim());
        meanCtrlPanel.setLabelPresenter(new LabelPresenter() {
            @Override
            public String getPresentationView(String label) {
                String[] pts = label.split("\\|")[0].trim().split("/");
                String col = pts[pts.length - 1];
                return col;
            }
        });
        meanFrame.setBounds(FCSPlot.START_X, FCSPlot.START_Y, FCSPlot.START_WIDTH, FCSPlot.START_HEIGHT);
        
        loadProps();
    }
    
    private void showMeanPanel(String col) {
        // Get all means (yData), get Files&Dates (xData), get Mean/SD (meanSD)
        double[] yDataBase, yDataComp;
        String[] baseLbls, compLbls;
        
        String currPanel = (String) comboPanel.getSelectedItem();
        if (currPanel.equals("All")) {
            currPanel = null;
        }
        String currCtrl = (String) comboControl.getSelectedItem();
        if (currCtrl.equals("All")) {
            currCtrl = null;
        }
        
        ArrayList<CtrlFileMetaData> metaData = baseData.getMetaDataFor(currCtrl, currPanel);

        metaData.sort(CtrlFileMetaData.COMPARATOR);
        
        yDataBase = new double[metaData.size()];
        baseLbls = new String[metaData.size()];
        int ind = 0;    
        float[] data = baseData.getParamData(col, currCtrl, currPanel);
        ArrayList<String> files = baseData.getAllFiles();
        for (CtrlFileMetaData cfm : metaData) {
            yDataBase[ind] = data[files.indexOf(cfm.file)];
            baseLbls[ind] = cfm.file;
            ind++;
        }
        
        int cnt = (Integer) spinner.getValue();
        ArrayList<DataFile> dataFiles = compDirData == null ? new ArrayList<ControlsQCGUI.DataFile>() : compDirData.getRecentFiles(currCtrl, currPanel, cnt);
        yDataComp = new double[dataFiles.size()];
        compLbls = new String[dataFiles.size()];
        ind = 0;
        for (int i = 0; i < dataFiles.size(); i++) {
            yDataComp[i] = Array.mean(dataFiles.get(i).getParamData(col, currCtrl, currPanel), true);
            compLbls[i] = dataFiles.get(i).fileName;
            
        }

        ArrayList<String> cols = new ArrayList<String>();
        ind = 0;
        for (int i = 1; i < table.getColumnModel().getColumnCount(); i++) {
            String hdr = (String) table.getValueAt(0, i);
            cols.add(hdr);
            if (hdr.equals(col)) {
                ind = i;
            }
        }
        
        meanCtrlPanel.setColumns(cols, ind-1);

        meanPanel.setXAxisLabel("");// pts[0].trim().replaceAll("/", " /\n");
        meanPanel.setYAxisLabel(col.split("\\|")[1].trim());
        meanPanel.setData(col, new String[][]{baseLbls, compLbls}, new double[][]{yDataBase, yDataComp});
        meanPanel.paintAgain();
        meanFrame.setTitle("Genvisis - Control QC - Overall Mean/SD - " + col);
        meanFrame.setVisible(true);
    }
    
    private static final String PROP_FILE = "controlsQC.properties";
    private static final String PROPKEY_COMPAREDIR = "COMPARE_DIR";
    private static final String PROPKEY_RANGEFILE = "RANGE_FILE";
    private static final String PROPKEY_COLS = "HIDDEN_COLUMNS";
    private JCheckBox chckbxShowInternalControl;
    private JSeparator separator_3;
    private JSpinner spinner;
    private JLabel lblOfControls;
    private JSeparator separator_4;

    private void saveProps() {
        try {
            Properties props = new Properties();
            props.setProperty(PROPKEY_RANGEFILE, txtFldBaseFile.getText());
            props.setProperty(PROPKEY_COMPAREDIR, ext.verifyDirFormat(txtFldCompDir.getText()));
            StringBuilder cols = new StringBuilder();
            int ind = 0;
            for (String c : hiddenCols) {
                cols.append(c);
                if (ind < hiddenCols.size() - 1) {
                    cols.append(";");
                }
                ind++;
            }
            props.setProperty(PROPKEY_COLS, cols.toString());
            File f = new File(PROP_FILE);
            OutputStream out = new FileOutputStream( f );
            props.store(out, "");
        } catch (Exception e ) {
            e.printStackTrace();
        }
    }
    
    private void loadProps() {
        Properties props = new Properties();
        InputStream is = null;
     
        try {
            File f = new File(PROP_FILE);
            if (!f.exists()) return;
            is = new FileInputStream(f);
            props.load(is);
            String base = props.getProperty(PROPKEY_RANGEFILE, "");
            String comp = props.getProperty(PROPKEY_COMPAREDIR, "");
            String colsTemp = props.getProperty(PROPKEY_COLS, "");
            String[] cols = colsTemp.split(";");
            
            if (!base.equals("")) {
                txtFldBaseFile.setText(base);
                loadFileOrDir(base, true);
            }
            if (!comp.equals("")) {
                txtFldCompDir.setText(comp);
                loadFileOrDir(comp, false);
            }
            this.hiddenCols.clear();
            for (String c : cols) {
                this.hiddenCols.add(c);
            }
            reCalcTableData();
        }
        catch ( Exception e ) {
            e.printStackTrace();
            is = null; 
        }
    }
    
    private void resizeColumnWidth() {
        final TableColumnModel columnModel = table.getColumnModel();
        FontMetrics fm = table.getFontMetrics(table.getFont());
        for (int column = 0; column < table.getColumnCount(); column++) {
            int width = fm.stringWidth(columnModel.getColumn(column).getHeaderValue().toString()) + 20; // Min width
            for (int row = 1; row < table.getRowCount(); row++) {
                TableCellRenderer renderer = table.getCellRenderer(row, column);
                Component comp = table.prepareRenderer(renderer, row, column);
                width = Math.max(comp.getPreferredSize().width +1 , width);
            }
            columnModel.getColumn(column).setPreferredWidth(width);
        }
    }
    
    static class CtrlFileMetaData {
        String file;
        String fileDateStr;
        Date fileDate;
        private String panel;
        private String runGroup;
        private String ctrlGroup;
        private String initials;
        private String number;

        public static final Comparator<CtrlFileMetaData> COMPARATOR = new Comparator<CtrlFileMetaData>() {
            @Override
            public int compare(CtrlFileMetaData o1, CtrlFileMetaData o2) {
                int comp;
                comp = o1.fileDate.compareTo(o2.fileDate);
                if (comp == 0) {
                    comp = o1.number.compareTo(o2.number);
                }
                return comp;
            }
        };
        
        private CtrlFileMetaData() {}
        
        public static CtrlFileMetaData parse(String filename) {
            CtrlFileMetaData cfmd = new CtrlFileMetaData();
            cfmd.file = filename;
            String[] pts = ext.rootOf(filename, true).split("_");
            if (pts.length != 6) {
                System.err.println("Error - file " + filename + " is an unexpected format!");
            }
            if (pts.length > 0) {
                cfmd.fileDateStr = pts[0];
                String[] dtPts = pts[0].split("-");
                try {
                    cfmd.fileDate = new GregorianCalendar(Integer.parseInt(dtPts[0]), Integer.parseInt(dtPts[1]), Integer.parseInt(dtPts[2])).getTime();
                } catch (NumberFormatException e) {
                    System.err.println("Error - filename " + filename + " does not contain a date string as the first token!");
                    cfmd.fileDate = null;
                }
            }
            if (pts.length > 1)
                cfmd.panel = pts[1].toUpperCase().replaceAll("[ \\-_/]", "").replaceAll("PANEL", "");
            if (pts.length > 2)
                cfmd.initials = pts[2].toUpperCase().replaceAll("[ \\-_/]", "");
            if (pts.length > 3)
                cfmd.runGroup = pts[3].toUpperCase().replaceAll("[ \\-_/]", "");
            if (pts.length > 4)
                cfmd.ctrlGroup = pts[4].toUpperCase().replaceAll("[ \\-_/]", "").replaceAll("CTL", "");
            if (pts.length > 5)
                cfmd.number = pts[5];
            return cfmd;
        }

        public String getPanel() {
            return panel;
        }

        public boolean isControlGroup(String ctrl) {
            if (ctrl == null && this.ctrlGroup != null) return false;
            if (this.ctrlGroup.equals(ctrl.replaceAll("[ \\-_/]", "")))
                return true;
            return false;
        }

        public boolean isPanel(String panel2) {
            if (panel2 == null && this.panel != null) return false;
            if (this.panel.equals(panel2.replaceAll("[ \\-_/]", "")))
                return true;
            return false;
        }
    }
    
    private void resetShownColumns() {
        TableColumnModel tcm = table.getColumnModel();
        ArrayList<Integer> toRemove = new ArrayList<Integer>();
        for (int i = 0; i < tcm.getColumnCount(); i++) {
            String hdr = (String) tcm.getColumn(i).getHeaderValue();
            if (!hdr.equals("") && hiddenCols.contains(hdr)) {
                toRemove.add(i);
            }
        }
        for (int i = toRemove.size() - 1; i >= 0; i--) {
            tcm.removeColumn(tcm.getColumn(toRemove.get(i)));
        }
    }
    
    static class DataFile {
        
        private String fileName;
        private ArrayList<String> allParams = new ArrayList<String>();
        private ArrayList<String> internalFiles = new ArrayList<String>();
        private float[][] fileData;
        private HashMap<String, float[]> paramData = new HashMap<String, float[]>();
        private HashMap<String, CtrlFileMetaData> internalMetaData = new HashMap<String, ControlsQCGUI.CtrlFileMetaData>();
        
        private DataFile() {}
        
        public static DataFile construct(String dataFile) {
            String[][] data = HashVec.loadFileToStringMatrix(dataFile, false, null, false);
            
            DataFile bd = new DataFile();
            bd.fileName = dataFile;
            
            for (int i = 1; i < data[0].length; i++) {
                if (!"".equals(data[0][i])) {
                    bd.allParams.add(data[0][i]);
                }
            }
            HashMap<String, Integer> fileInd = new HashMap<String, Integer>();
            for (int i = 1; i < data.length; i++) {
                int ind = ext.indexOfStr(data[i][0], EXCLUDED_ROW_HEADERS, false, false);
                if ("".equals(data[i][0]) || ind >= 0) {
                    continue;
                }
                bd.internalFiles.add(data[i][0]);
                bd.internalMetaData.put(data[i][0], CtrlFileMetaData.parse(data[i][0]));
                fileInd.put(data[i][0], i);
            }
            
            bd.fileData = new float[bd.internalFiles.size()][bd.allParams.size()];
            
            for (int f = 0; f < bd.internalFiles.size(); f++) {
                int rowInd = fileInd.get(bd.internalFiles.get(f));
                float[] arr = new float[bd.allParams.size()];
                for (int i = 0; i < bd.allParams.size(); i++) {
                    arr[i] = "".equals(data[rowInd][i+1]) ? Float.NaN : Float.parseFloat(data[rowInd][i+1]);
                }
                bd.fileData[f] = arr;
            }
            
            return bd;
        }
        
        public String[] getPanels() {
            TreeSet<String> panelSet = new TreeSet<String>();
            for (CtrlFileMetaData cfmd : this.internalMetaData.values()) {
                panelSet.add(cfmd.getPanel());
            }
            return panelSet.toArray(new String[panelSet.size()]);
        }

        public String[] getControlGroups() {
            TreeSet<String> ctrlSet = new TreeSet<String>();
            for (CtrlFileMetaData cfmd : this.internalMetaData.values()) {
                ctrlSet.add(cfmd.ctrlGroup);
            }
            return ctrlSet.toArray(new String[ctrlSet.size()]);
        }
        
        public float[] getFileData(String internalFile) {
            float[] retArr = null;
            int ind = this.internalFiles.indexOf(internalFile);
            if (ind > -1) {
                return this.fileData[ind];
            }
            return retArr;
        }
        
        public float[] getParamData(String paramName) {
            float[] retArr = null;
            if (this.paramData.containsKey(paramName)) {
                retArr = this.paramData.get(paramName);
            } else {
                int ind = this.allParams.indexOf(paramName);
                if (ind > -1) {
                    retArr = Matrix.extractColumn(fileData, ind);
                    this.paramData.put(paramName, retArr);
                }
            }
            return retArr;
        }
        
        public float[] getParamData(String param, String ctrl) {
            float[] retArr = this.paramData.get(param + "\t" + ctrl);
            if (retArr == null) {
                retArr = getParamData(param);
                if (retArr != null) { 
                    boolean[] incl = new boolean[this.internalFiles.size()];
                    for (int i = 0; i < incl.length; i++) {
                        incl[i] = internalMetaData.get(this.internalFiles.get(i)).ctrlGroup.equals(ctrl);
                    }
                    retArr = Array.subArray(retArr, incl);
                    paramData.put(param + "\t" + ctrl, retArr);
                }
            }
            return retArr;
        }

        public float[] getParamData(String param, String ctrl, String panel) {
            float[] retArr = this.paramData.get(param + "\t" + ctrl + "\t" + panel);
            if (retArr == null) {
                retArr = getParamData(param);
                if (retArr != null) { 
                    boolean[] incl = new boolean[this.internalFiles.size()];
                    for (int i = 0; i < incl.length; i++) {
                        CtrlFileMetaData cfmd = internalMetaData.get(this.internalFiles.get(i));
                        incl[i] = (ctrl == null || cfmd.isControlGroup(ctrl)) && (panel == null || cfmd.isPanel(panel));
                    }
                    retArr = Array.subArray(retArr, incl);
                    paramData.put(param + "\t" + ctrl + "\t" + panel, retArr);
                }
            }
            return retArr;
        }

        @SuppressWarnings("unchecked")
        public ArrayList<String> getAllParams() {
            return (ArrayList<String>) allParams.clone();
        }

        public ArrayList<CtrlFileMetaData> getMetaDataFor(String currCtrl, String currPanel) {
            ArrayList<CtrlFileMetaData> retList = new ArrayList<ControlsQCGUI.CtrlFileMetaData>();
            for (CtrlFileMetaData cfmd : this.internalMetaData.values()) {
                if ((currCtrl == null || cfmd.isControlGroup(currCtrl)) && (currPanel == null || cfmd.isPanel(currPanel))) { 
                    retList.add(cfmd);
                }
            }
            return retList;
        }
        
        @SuppressWarnings("unchecked")
        public ArrayList<String> getAllFiles() {
            return (ArrayList<String>) internalFiles.clone();
        }
    
    }
    
    
    static class CompDir {
        
        private static final FilenameFilter csvFileFilter = new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                return name.endsWith(".csv");
            }
        };
        private static final FileFilter dirFilter = new FileFilter() {
            @Override
            public boolean accept(File pathname) {
                return pathname.isDirectory();
            }
        };
        
        private CompDir() {}
        
        private String dir;
        private HashMap<String, DataFile> dataFiles = new HashMap<String, ControlsQCGUI.DataFile>();
        private HashMap<String, CtrlFileMetaData> fileMetaData = new HashMap<String, ControlsQCGUI.CtrlFileMetaData>();
        private ArrayList<CompDir> subDirs = new ArrayList<ControlsQCGUI.CompDir>();
        
        public static CompDir construct(String dir) {
            File dirFile = new File(dir);
            String[] rootCSV = dirFile.list(csvFileFilter);
            File[] subDirs = dirFile.listFiles(dirFilter);
            
            CompDir cd = new CompDir();
            cd.dir = ext.verifyDirFormat(dir);
            
            if (rootCSV != null) {
                for (String f : rootCSV) {
                    cd.dataFiles.put(f, DataFile.construct(dir + f));
                    cd.fileMetaData.put(f, CtrlFileMetaData.parse(f));
                }
            }
            if (subDirs != null) {
                for (File f : subDirs) {
                    cd.subDirs.add(CompDir.construct(ext.verifyDirFormat(f.getPath())));
                }
            }
            
            return cd;
        }
        
        private TreeMap<CtrlFileMetaData, DataFile> getFilesMap(String currCtrl, String currPanel) {
            TreeMap<CtrlFileMetaData, DataFile> map = new TreeMap<ControlsQCGUI.CtrlFileMetaData, ControlsQCGUI.DataFile>(CtrlFileMetaData.COMPARATOR);
            for (Entry<String, DataFile> ent : dataFiles.entrySet()) {
                CtrlFileMetaData cfmd = fileMetaData.get(ent.getKey()); 
                if ((currCtrl == null || cfmd.isControlGroup(currCtrl)) && (currPanel == null || cfmd.isPanel(currPanel))) {
                    map.put(cfmd, ent.getValue());
                }
            }
            for (CompDir cd : subDirs) {
                map.putAll(cd.getFilesMap(currCtrl, currPanel));
            }
            return map;
        }
        
        public ArrayList<DataFile> getRecentFiles(String currCtrl, String currPanel, int numOfMostRecentToAdd) {
            ArrayList<DataFile> recentFiles = new ArrayList<ControlsQCGUI.DataFile>();
            
            TreeMap<CtrlFileMetaData, DataFile> map = getFilesMap(currCtrl, currPanel);
            int ind = 0;
            for (CtrlFileMetaData cfmd : map.keySet()) {
                recentFiles.add(map.get(cfmd));
                ind++;
                if (numOfMostRecentToAdd > 0 && numOfMostRecentToAdd <= ind) {
                    break;
                }
            }
            
            return recentFiles;
        }
        
    }
    
    private void reCalcTableData() {
        if (baseData == null /*|| compDirData == null*/) return;
        
        TreeSet<String> paramSet = new TreeSet<String>();
        
        paramSet.addAll(baseData.getAllParams()); // don't have to worry about panel/ctrl for params
        
        String[] paramNames = paramSet.toArray(new String[paramSet.size()]);
        HashMap<String, String> colLookup = new HashMap<String, String>();
        String[] colNames = new String[paramNames.length + 1];
        colNames[0] = "";
        for (int i = 0; i < paramNames.length; i++) {
            String col = paramNames[i];
            String[] pts = col.split("\\|")[0].trim().split("/");
            col = pts[pts.length - 1];
            colNames[i+1] = col;
            colLookup.put(col, paramNames[i]);
            colLookup.put(paramNames[i], col);
        }
        for (String p : paramNames) {
            paramMeanLists.put(p, new ArrayList<Float>());
        }
        
        dtmMean = new DefaultTableModel(colNames, 0) {
            @Override
            public boolean isCellEditable(int row, int column) {
                return false;
            }
        };
        
        boldRows.clear();
        
        String[] firstRow = Array.addStrToArray("", paramNames.clone(), 0);
        firstRow[0] = "Source";
        dtmMean.addRow(firstRow);
        boldRows.add(0);
        
        addDataToModel(baseData);
        
        dtmMean.addRow(new Object[colNames.length]);
        statRows.add(dtmMean.getRowCount());
        
        Object[] meanRow = new Object[colNames.length];
        meanRow[0] = "Mean of means";
        for (int i = 1; i < colNames.length; i++) {
            String colNm = colLookup.get(colNames[i]);
            if (paramMeanLists.containsKey(colNm)) {
                Float mn = Array.mean(Array.toFloatArray(paramMeanLists.get(colNm)), true);
                paramMeans.put(colNm, mn);
                meanRow[i] = mn;
            }
        }
        dtmMean.addRow(meanRow);
        statRows.add(dtmMean.getRowCount());
        
        Object[] sdRow = new Object[colNames.length];
        sdRow[0] = "StdDev of means";
        for (int i = 1; i < colNames.length; i++) {
            String colNm = colLookup.get(colNames[i]);
            if (paramMeanLists.containsKey(colNm)) {
                Float sd = Array.stdev(paramMeanLists.get(colNm).toArray(new Float[0]), true);
                paramSDs.put(colNm, sd);
                sdRow[i] = sd;
            }
        }
        dtmMean.addRow(sdRow);
        statRows.add(dtmMean.getRowCount());

        Object[] cvRow = new Object[colNames.length];
        cvRow[0] = "cV ( = 100 * SD / Mean)";
        for (int i = 1; i < colNames.length; i++) {
            String colNm = colLookup.get(colNames[i]);
            if (paramMeanLists.containsKey(colNm)) {
                Float cv = 100 * (paramSDs.get(colNm) / paramMeans.get(colNm));
                paramCVs.put(colNm, cv);
                cvRow[i] = cv;
            }
        }
        dtmMean.addRow(cvRow);
        statRows.add(dtmMean.getRowCount());
        
        Object[] rgRow = new Object[colNames.length];
        rgRow[0] = "Mean - 1SD";
        for (int i = 1; i < colNames.length; i++) {
            String colNm = colLookup.get(colNames[i]);
            if (paramMeans.containsKey(colNm)) {
                rgRow[i] = paramMeans.get(colNm) - paramSDs.get(colNm);
            }
        }
        dtmMean.addRow(rgRow);
        statRows.add(dtmMean.getRowCount());
        
        rgRow = new Object[colNames.length];
        rgRow[0] = "Mean + 1SD";
        for (int i = 1; i < colNames.length; i++) {
            String colNm = colLookup.get(colNames[i]);
            if (paramMeans.containsKey(colNm)) {
                rgRow[i] = paramMeans.get(colNm) + paramSDs.get(colNm);
            }
        }
        dtmMean.addRow(rgRow);
        statRows.add(dtmMean.getRowCount());
        
        dtmMean.addRow(new Object[colNames.length]);

        int numRecent = (Integer) spinner.getValue();
//        addFilesToModel(compDirData, compDirData.dir, numRecent);
        
        table.setModel(dtmMean);
        table.revalidate();
        table.repaint();
        resizeColumnWidth();
        resetShownColumns();
        
    }
    
    private void addDataToModel(DataFile df) {
        String currPanel = (String) comboPanel.getSelectedItem();
        if ("All".equals(currPanel)) {
            currPanel = null;
        }
        String currCtrl = (String) comboControl.getSelectedItem();
        if ("All".equals(currCtrl)) {
            currCtrl = null;
        }
        ArrayList<CtrlFileMetaData> metaData = df.getMetaDataFor(currCtrl, currPanel);
        metaData.sort(CtrlFileMetaData.COMPARATOR);
        for (CtrlFileMetaData cfmd : metaData) {
            String f = cfmd.file;
            
            ArrayList<String> params = df.getAllParams();
            float[] fileData = df.getFileData(f);
            Object[] rowDataM = new Object[fileData.length + 1];
            rowDataM[0] = ext.rootOf(f);
            for (int i = 0; i < fileData.length; i++) {
                paramMeanLists.get(params.get(i)).add(fileData[i]);
                rowDataM[i + 1] = fileData[i];
            }
            dtmMean.addRow(rowDataM);
        }
    }
    
    private void addFilesToModel(CompDir compData, String removePrep, int numOfMostRecentToAdd) {
        String currPanel = (String) comboPanel.getSelectedItem();
        if (currPanel.equals("All")) {
            currPanel = null;
        }
        String currCtrl = (String) comboControl.getSelectedItem();
        if (currCtrl.equals("All")) {
            currCtrl = null;
        }
        ArrayList<String> params = baseData.getAllParams();
        boolean showInternalData = chckbxShowInternalControl.isSelected();
        
        ArrayList<DataFile> files = compData.getRecentFiles(currCtrl, currPanel, numOfMostRecentToAdd);
        
        for (int d = 0; d < (numOfMostRecentToAdd == 0 ? files.size() : Math.min(files.size(), numOfMostRecentToAdd)); d++) {
            DataFile df = files.get(d);
            
            if (showInternalData) {
                addDataToModel(df);
            }
            
            Object[] rowDataM = new Object[params.size() + 1];
            rowDataM[0] = ext.rootOf(df.fileName, true);
            
            for (int i = 0; i < params.size(); i++) {
                float[] pData = df.getParamData(params.get(i), currCtrl, currPanel);
                Float mn = pData == null ? Float.NaN : Array.mean(pData);
                rowDataM[i + 1] = mn;
            }
            dtmMean.addRow(rowDataM);
            
        }
    }
    
    /**
     * Launch the application.
     */
    public static void main(String[] args) {
        EventQueue.invokeLater(new Runnable() {
            public void run() {
                try {
                    ControlsQCGUI frame = new ControlsQCGUI();
                    frame.setVisible(true);
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        });
    }
    
    
}
