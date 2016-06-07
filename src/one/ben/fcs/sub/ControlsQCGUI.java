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
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.Properties;
import java.util.TreeMap;
import java.util.TreeSet;

import javax.swing.AbstractAction;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.ScrollPaneConstants;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.table.DefaultTableColumnModel;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumnModel;
import javax.xml.parsers.ParserConfigurationException;

import net.miginfocom.swing.MigLayout;
import one.ben.fcs.AbstractPanel2.PLOT_TYPE;
import one.ben.fcs.FCSDataLoader;
import one.ben.fcs.FCSDataLoader.DATA_SET;
import one.ben.fcs.FCSPlot;
import one.ben.fcs.gating.Gate;
import one.ben.fcs.gating.GateFileReader;
import one.ben.fcs.gating.GatingStrategy;

import org.xml.sax.SAXException;

import cnv.gui.IncludeExcludeGUI;
import common.Array;
import common.HashVec;
import common.ext;


public class ControlsQCGUI extends JFrame {

    private JPanel contentPane;
    private JTextField txtFldBaseFile;
    private JTextField txtFldGatingFile;
    private JButton btnBaseDirSelect;
    private JButton btnGatingFileSelect;
    private JTable table;
    private JScrollPane scrollPane;

//    private HashMap<String, FCSDataLoader> baseFiles = new HashMap<String, FCSDataLoader>();
//    private HashMap<String, FCSDataLoader> compFiles = new HashMap<String, FCSDataLoader>();
    
    HashMap<String, ArrayList<Float>> paramMeanLists = new HashMap<String, ArrayList<Float>>();
    HashMap<String, ArrayList<Float>> paramSDLists = new HashMap<String, ArrayList<Float>>();
    HashMap<String, ArrayList<Float>> paramCVLists = new HashMap<String, ArrayList<Float>>();
    HashMap<String, HashMap<String, Float>> fileParamMeanMap = new HashMap<String, HashMap<String,Float>>();
    HashSet<Integer> boldRows = new HashSet<Integer>();
    HashSet<Integer> statRows = new HashSet<Integer>();
    HashMap<String, Float> paramMeans = new HashMap<String, Float>();
    HashMap<String, Float> paramSDs = new HashMap<String, Float>();
    HashMap<String, Float> paramCVs = new HashMap<String, Float>();
    
    TreeSet<String> paramNames = new TreeSet<String>();
    HashSet<String> hiddenCols = new HashSet<String>();
    
    JFrame meanFrame = new JFrame("Genvisis - FCS Overall Mean/SD");
    MeanPanel meanPanel = new MeanPanel();
    MeanCtrlPanel meanCtrlPanel = new MeanCtrlPanel();
    GatingStrategy gateStrat;
//    String baseDir;
//    String[] baseFCSFiles;
    String baseFile;
    DirFile dirStruct;
    private JLabel lblCompareFiles;
    private JTextField txtFldCompDir;
    private JButton btnCompDirSelect;
    private JRadioButton rdbtnMean;
    private JRadioButton rdbtnSd;
    private JRadioButton rdbtnCv;
    private final ButtonGroup buttonGroup = new ButtonGroup();

    protected void loadFileOrDir(String fileOrDir, boolean baseline) {
        if (baseline) {
            if (baseFile != null && baseFile.equals(fileOrDir)) {
                return; // same as already loaded
            }
            baseFile = fileOrDir;
        } else {
            if (dirStruct != null && dirStruct.dir.equals(ext.verifyDirFormat(fileOrDir))) {
                return; // same as already loaded
            }
            dirStruct = loadCompFiles(fileOrDir);
        }
        reCalcTableData();
    }


    private DefaultTableModel dtmMean;
    private DefaultTableModel dtmSD;
    private DefaultTableModel dtmCV;
    private JSeparator separator_1;
    private JButton btnHideshowColumns;
    Color SD1_COLOR = Color.YELLOW;
    Color SD2_COLOR = Color.RED;
    Color ABOVE_3_CV_COLOR = Color.CYAN;
    Color ABOVE_5_CV_COLOR = Color.CYAN.darker();

    final HashSet<String> lessThan4CV = new HashSet<String>();
    final HashSet<String> lessThan6CV = new HashSet<String>();
    
    {
        String[] lessThan4 = {"BB515", "PE", "PE-CF594", "PE-Cy7", "BV 421", "BV 510", "BV 605", "BV 711"};
        String[] lessThan6 = {"BUV 395", "BUV 737", "APC", "APC-Cy7"};
        for (String s : lessThan4) {
            lessThan4CV.add(s);
            lessThan4CV.add(s + "-A");
            lessThan4CV.add("Comp-" + s);
            lessThan4CV.add("Comp-" + s + "-A");
        }
        for (String s : lessThan6) {
            lessThan6CV.add(s);
            lessThan6CV.add("Comp-" + s);
        }
    }

    
    {
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
        meanCtrlPanel.setChangeListener(prevLst);
        meanFrame.setBounds(FCSPlot.START_X, FCSPlot.START_Y, FCSPlot.START_WIDTH, FCSPlot.START_HEIGHT);
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
        contentPane.setLayout(new MigLayout("ins 7 7 3 7,hidemode 3", "[][grow][]", "[][][][grow][][]"));
        
        JLabel lblFileDir = new JLabel("<html><u>B</u>aseline File:</html>");
        contentPane.add(lblFileDir, "cell 0 0,alignx trailing");
        
        txtFldBaseFile = new JTextField();
        txtFldBaseFile.setEditable(false);
        contentPane.add(txtFldBaseFile, "cell 1 0,growx");
        txtFldBaseFile.setColumns(10);
        
        Insets btnInsets = new Insets(0, 14, 0, 14);
        
        btnBaseDirSelect = new JButton(">");
        btnBaseDirSelect.setMargin(btnInsets);
        btnBaseDirSelect.setMnemonic(KeyEvent.VK_B);
        btnBaseDirSelect.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent arg0) {
                String curr = txtFldBaseFile.getText();
                if (curr.equals("")) {
                    curr = "./";
                }
                JFileChooser jfc = new JFileChooser(curr);
                jfc.setFileSelectionMode(JFileChooser.FILES_ONLY);
                jfc.setDialogTitle("Select Baseline File");
                jfc.setMultiSelectionEnabled(false);
                jfc.setFileFilter(new FileNameExtensionFilter("CSV Files", ".csv"));
                int resp = jfc.showOpenDialog(ControlsQCGUI.this);
                if (resp == JFileChooser.APPROVE_OPTION) {
                    String newPath = jfc.getSelectedFile().getAbsolutePath();
                    txtFldBaseFile.setText(newPath);
                    loadFileOrDir(newPath, true);
                    saveProps();
                }
            }
        });
        contentPane.add(btnBaseDirSelect, "cell 2 0");
        
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
                jfc.setFileSelectionMode(JFileChooser.FILES_ONLY);
                jfc.setDialogTitle("Select Compare Files");
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
        contentPane.add(scrollPane, "cell 0 3 3 1,grow");
                
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
                            if (rdbtnMean.isSelected()) {
                                if (paramMeans.containsKey(colNm)) {
                                    if (value > (paramMeans.get(colNm) + paramSDs.get(colNm)) || value < (paramMeans.get(colNm) - paramSDs.get(colNm))) {
                                        col = SD1_COLOR;
                                    }  
                                    if (value > (paramMeans.get(colNm) + 2*paramSDs.get(colNm)) || value < (paramMeans.get(colNm) - 2*paramSDs.get(colNm))) {
                                        col = SD2_COLOR;
                                    }
                                }
                            } else if (rdbtnSd.isSelected()) {
                                
                            } else if (rdbtnCv.isSelected()) {
//                                    <4% :
//                                    BLue laser: BB515, PE, PE-CF594 and PE-Cy7
//                                    Violet laser: BV421, BV510, BV605 and BV711
//
//                                    <6%
//                                    UV laser: BUV395 and BUV737
//                                    Red laser: APC and APC-Cy7
                                if (paramCVs.containsKey(colNm)) {
                                    if (value > 3 && lessThan4CV.contains(colNm)) {
                                        col = ABOVE_3_CV_COLOR;
                                    }
                                    if (value > 5 && lessThan6CV.contains(colNm)) {
                                        col = ABOVE_5_CV_COLOR;
                                    }
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
        
        btnHideshowColumns = new JButton("Hide/Show Columns");
        btnHideshowColumns.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent arg0) {
                String[] opts = new String[paramNames.size()];
                boolean[] incl = new boolean[paramNames.size()];
                int ind = 0;
                for (String p : paramNames) {
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
        contentPane.add(btnHideshowColumns, "flowx,cell 0 4 2 1");
        
        separator_1 = new JSeparator();
        separator_1.setOrientation(SwingConstants.VERTICAL);
        contentPane.add(separator_1, "cell 0 4,growy");

        rdbtnMean = new JRadioButton();
        rdbtnMean.setAction(new AbstractAction() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if (dtmMean != null) {
                    table.setModel(dtmMean);
                    table.revalidate();
                    table.repaint();
                    resizeColumnWidth();
                }
            }
        });
        rdbtnMean.setText("Mean");
        rdbtnMean.setMnemonic(KeyEvent.VK_M);
        rdbtnMean.setSelected(true);
        buttonGroup.add(rdbtnMean);
        contentPane.add(rdbtnMean, "cell 0 4");
        
        rdbtnSd = new JRadioButton();
        rdbtnSd.setAction(new AbstractAction() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if (dtmSD != null) {
                    table.setModel(dtmSD);
                    table.revalidate();
                    table.repaint();
                    resizeColumnWidth();
                } else {
                    rdbtnMean.setSelected(true);
                }
            }
        });
        rdbtnSd.setText("SD");
        rdbtnSd.setMnemonic(KeyEvent.VK_S);
        buttonGroup.add(rdbtnSd);
        contentPane.add(rdbtnSd, "cell 0 4");
        
        rdbtnCv = new JRadioButton();
        rdbtnCv.setAction(new AbstractAction() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if (dtmCV != null) {
                    table.setModel(dtmCV);
                    table.revalidate();
                    table.repaint();
                    resizeColumnWidth();
                } else {
                    rdbtnMean.setSelected(true);
                }
            }
        });
        rdbtnCv.setText("cV");
        rdbtnCv.setMnemonic(KeyEvent.VK_C);
        buttonGroup.add(rdbtnCv);
        contentPane.add(rdbtnCv, "cell 0 4");
        
        separator_2 = new JSeparator();
        separator_2.setOrientation(SwingConstants.VERTICAL);
        contentPane.add(separator_2, "cell 0 4,growy");
        
        loadProps();
    }
    
    private void showMeanPanel(String col) {
        // Get all means (yData), get Files&Dates (xData), get Mean/SD (meanSD)
        float[] xDataBase, xDataComp, yDataBase, yDataComp;
        String[] baseLbls, compLbls;
        
        TreeMap<Date, Float> meanMap = new TreeMap<Date, Float>();
        TreeMap<Date, String> fileMap = new TreeMap<Date, String>();
        for(Entry<String, FCSDataLoader> l : baseFiles.entrySet()) {
            meanMap.put(l.getValue().getRunDate(), fileParamMeanMap.get(l.getKey()).get(col));
            fileMap.put(l.getValue().getRunDate(), l.getKey());
        }
        
        xDataBase = new float[meanMap.size()];
        yDataBase = new float[meanMap.size()];
        baseLbls = new String[meanMap.size()];
        int ind = 0;
        for (Entry<Date, Float> etr : meanMap.entrySet()) {
            xDataBase[ind] = ind+1; // TODO should be time?
            yDataBase[ind] = etr.getValue();
            baseLbls[ind] = fileMap.get(etr.getKey());
            ind++;
        }
        
        meanMap = new TreeMap<Date, Float>();
        fileMap = new TreeMap<Date, String>();
        for(Entry<String, FCSDataLoader> l : compFiles.entrySet()) {
            meanMap.put(l.getValue().getRunDate(), fileParamMeanMap.get(ext.removeDirectoryInfo(l.getKey())).get(col));
            fileMap.put(l.getValue().getRunDate(), l.getKey());
        }
        xDataComp = new float[meanMap.size()];
        yDataComp = new float[meanMap.size()];
        compLbls = new String[meanMap.size()];
        for (Entry<Date, Float> etr : meanMap.entrySet()) {
            xDataComp[ind - xDataBase.length] = ind+1; // TODO should be time?
            yDataComp[ind - xDataBase.length] = etr.getValue();
            compLbls[ind - xDataBase.length] = fileMap.get(etr.getKey());
            ind++;
        }
        
        ArrayList<String> cols = new ArrayList<String>();
        ind = 0;
        for (int i = 1; i < table.getColumnModel().getColumnCount(); i++) {
            String hdr = (String) table.getColumnModel().getColumn(i).getHeaderValue();
            cols.add(hdr);
            if (hdr.equals(col)) {
                ind = i;
            }
        }
        
        meanCtrlPanel.setColumns(cols, ind-1);
        
        meanPanel.setData(col, baseLbls, xDataBase, xDataComp, compLbls, yDataBase, yDataComp);
        meanPanel.paintAgain();
        meanFrame.setTitle("Genvisis - FCS Overall Mean/SD - " + col);
        meanFrame.setVisible(true);
    }
    
    private static final String PROP_FILE = "controlsQC.properties";
    private static final String PROPKEY_COMPAREDIR = "COMPARE_DIR";
    private static final String PROPKEY_RANGEFILE = "RANGE_FILE";
    private static final String PROPKEY_COLS = "HIDDEN_COLUMNS";

    private void saveProps() {
        try {
            Properties props = new Properties();
            props.setProperty(PROPKEY_RANGEFILE, ext.verifyDirFormat(txtFldBaseFile.getText()));
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
        catch ( Exception e ) { is = null; }
    }
    
    private void resizeColumnWidth() {
        final TableColumnModel columnModel = table.getColumnModel();
        FontMetrics fm = table.getFontMetrics(table.getFont());
        for (int column = 0; column < table.getColumnCount(); column++) {
            int width = fm.stringWidth(columnModel.getColumn(column).getHeaderValue().toString()) + 20; // Min width
            for (int row = 0; row < table.getRowCount(); row++) {
                TableCellRenderer renderer = table.getCellRenderer(row, column);
                Component comp = table.prepareRenderer(renderer, row, column);
                width = Math.max(comp.getPreferredSize().width +1 , width);
            }
            columnModel.getColumn(column).setPreferredWidth(width);
        }
    }
    
    static class DirFile {
        String dir;
        String[] files;
        DirFile[] subDirs;
        
        public String[] getAllFiles() {
            ArrayList<String> files = new ArrayList<String>();
            for (String f : this.files) {
                files.add(dir + f);
            }
            for (DirFile df : subDirs) {
                for (String f : df.getAllFiles()) {
                    files.add(f);
                }
            }
            return files.toArray(new String[files.size()]);
        }
        
        public String[] getAllFiles(String removePrep) {
            ArrayList<String> files = new ArrayList<String>();
            for (String f : this.files) {
                files.add((dir + f).replaceFirst(removePrep, ""));
            }
            for (DirFile df : subDirs) {
                for (String f : df.getAllFiles(removePrep)) {
                    files.add(f);
                }
            }
            return files.toArray(new String[files.size()]);
        }
        
    }
    

    private final FilenameFilter csvFileFilter = new FilenameFilter() {
        @Override
        public boolean accept(File dir, String name) {
            return name.endsWith(".csv");
        }
    };
    private final FileFilter dirFilter = new FileFilter() {
        @Override
        public boolean accept(File pathname) {
            return pathname.isDirectory();
        }
    };
    private JSeparator separator_2;
    
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
    
    HashMap<String, String[]> compRowIDs = new HashMap<String, String[]>();
    HashMap<String, float[][]> compData = new HashMap<String, float[][]>();
    String[] baseRowIDs;
    HashMap<String, float[]> baseData = new HashMap<String, float[]>(); 
    private final String[] EXCLUDED_ROW_HEADERS = {
            "mean",
            "average",
            "sd",
            "cv",
            "cv (%)"
    };
    
    private void reCalcTableData() {
        
        if (baseFile == null || "".equals(baseFile)) return;
        if (dirStruct == null) return;
        
        paramNames.clear();
        baseData.clear();
        TreeSet<String> paramSet = new TreeSet<String>();
        String[][] data = HashVec.loadEntireFileToStringMatrixSplittingCommasIntelligently(baseFile);
        
        ArrayList<String> cols = new ArrayList<String>();
        ArrayList<String> rows = new ArrayList<String>();
        HashMap<String, Integer> rowInds = new HashMap<String, Integer>();
        
        for (int i = 1; i < data[0].length; i++) {
            cols.add(data[0][i]);
        }
        for (int i = 1; i < data.length; i++) {
            int ind = ext.indexOfStr(data[i][0], EXCLUDED_ROW_HEADERS, false, false);
            if ("".equals(data[i][0]) || ind >= 0) {
                continue;
            }
            rows.add(data[i][0]);
            rowInds.put(data[i][0], i);
        }
        for (int i = 0; i < cols.size(); i++) {
            float[] arr = new float[rows.size()];
            int ind = 0;
            for (String r : rows) {
                arr[ind] = "".equals(data[rowInds.get(r)][i+1]) ? Float.NaN : Float.parseFloat(data[rowInds.get(r)][i+1]);
            }
            baseData.put(cols.get(i), arr);
        }
        
        baseRowIDs = rows.toArray(new String[rows.size()]);
        
        String[] allFilesFullPaths = dirStruct.getAllFiles();
        for (String f : allFilesFullPaths) {
            if (!this.compData.containsKey(f)) {
                data = HashVec.loadEntireFileToStringMatrixSplittingCommasIntelligently(f);
                baseRow = new ArrayList<String>();
                baseTemp = new ArrayList<float[]>();
                for (String s : data[0]) {
                    if (!"".equals(s)) {
                        paramNames.add(s);
                        paramSet.add(s);
                    }
                }
            }
            for (int i = 1; i < data.length; i++) {
                int ind = ext.indexOfStr(data[i][0], EXCLUDED_ROW_HEADERS, false, false);
                if ("".equals(data[i][0]) || ind >= 0) {
                    continue;
                }
                baseRow.add(data[i][0]);
                float[] tempData = new float[data[i].length - 1];
                for (int j = 1; j < data[i].length; j++) {
                    tempData[j-1] = Float.parseFloat(data[i][j]);
                }
                baseTemp.add(tempData);
            }
            compRowIDs.put(f, baseRow.toArray(new String[baseRow.size()]));
            compData.put(f, baseTemp.toArray(new float[baseTemp.size()][]));
        }
        
        String[] paramNames = paramSet.toArray(new String[paramSet.size()]);
        String[] colNames = Array.addStrToArray("", paramNames, 0);
        for (String p : paramNames) {
            paramMeanLists.put(p, new ArrayList<Float>());
            paramSDLists.put(p, new ArrayList<Float>());
            paramCVLists.put(p, new ArrayList<Float>());
        }
        
        dtmMean = new DefaultTableModel(colNames, 0) {
            @Override
            public boolean isCellEditable(int row, int column) {
                return false;
            }
        };
        dtmSD = new DefaultTableModel(colNames, 0) {
            @Override
            public boolean isCellEditable(int row, int column) {
                return false;
            }
        };
        dtmCV = new DefaultTableModel(colNames, 0) {
            @Override
            public boolean isCellEditable(int row, int column) {
                return false;
            }
        };
        
        boldRows.clear();
        
        String[] firstRow = colNames.clone();
        firstRow[0] = "Source";
        dtmMean.addRow(firstRow);
        dtmSD.addRow(firstRow);
        dtmCV.addRow(firstRow);
        boldRows.add(0);
        
        addBaseToModel(paramNames);
        
        dtmMean.addRow(new Object[colNames.length]);
        dtmSD.addRow(new Object[colNames.length]);
        dtmCV.addRow(new Object[colNames.length]);
        statRows.add(dtmMean.getRowCount());
        
        Object[] meanRow = new Object[colNames.length];
        meanRow[0] = "Mean of means";
        for (int i = 1; i < colNames.length; i++) {
            String colNm = colNames[i];
            if (paramMeanLists.containsKey(colNm)) {
                Float mn = Array.mean(Array.toFloatArray(paramMeanLists.get(colNm)), true);
                paramMeans.put(colNm, mn);
                meanRow[i] = mn;
            }
        }
        dtmMean.addRow(meanRow);
        dtmSD.addRow(meanRow);
        dtmCV.addRow(meanRow);
        statRows.add(dtmMean.getRowCount());
        
        Object[] sdRow = new Object[colNames.length];
        sdRow[0] = "StdDev of means";
        for (int i = 1; i < colNames.length; i++) {
            String colNm = colNames[i];
            if (paramMeanLists.containsKey(colNm)) {
                Float sd = Array.stdev(paramMeanLists.get(colNm).toArray(new Float[0]), true);
                paramSDs.put(colNm, sd);
                sdRow[i] = sd;
            }
        }
        dtmMean.addRow(sdRow);
        dtmSD.addRow(sdRow);
        dtmCV.addRow(sdRow);
        statRows.add(dtmMean.getRowCount());

        Object[] cvRow = new Object[colNames.length];
        cvRow[0] = "cV ( = 100 * SD / Mean)";
        for (int i = 1; i < colNames.length; i++) {
            String colNm = colNames[i];
            if (paramMeanLists.containsKey(colNm)) {
                Float cv = 100 * (paramSDs.get(colNm) / paramMeans.get(colNm));
                paramCVs.put(colNm, cv);
                cvRow[i] = cv;
            }
        }
        dtmMean.addRow(cvRow);
        dtmSD.addRow(cvRow);
        dtmCV.addRow(cvRow);
        statRows.add(dtmMean.getRowCount());
        
        Object[] rgRow = new Object[colNames.length];
        rgRow[0] = "Mean - 1SD";
        for (int i = 1; i < colNames.length; i++) {
            String colNm = colNames[i];
            if (paramMeans.containsKey(colNm)) {
                rgRow[i] = paramMeans.get(colNm) - paramSDs.get(colNm);
            }
        }
        dtmMean.addRow(rgRow);
        dtmSD.addRow(rgRow);
        dtmCV.addRow(rgRow);
        statRows.add(dtmMean.getRowCount());
        
        rgRow = new Object[colNames.length];
        rgRow[0] = "Mean + 1SD";
        for (int i = 1; i < colNames.length; i++) {
            String colNm = colNames[i];
            if (paramMeans.containsKey(colNm)) {
                rgRow[i] = paramMeans.get(colNm) + paramSDs.get(colNm);
            }
        }
        dtmMean.addRow(rgRow);
        dtmSD.addRow(rgRow);
        dtmCV.addRow(rgRow);
        statRows.add(dtmMean.getRowCount());
        
        dtmMean.addRow(new Object[colNames.length]);
        dtmSD.addRow(new Object[colNames.length]);
        dtmCV.addRow(new Object[colNames.length]);
        addFilesToModel(dirStruct, paramNames, dirStruct.dir);
        
        if (rdbtnMean.isSelected()) {
            table.setModel(dtmMean);
        } else if (rdbtnSd.isSelected()) {
            table.setModel(dtmSD);
        } else if (rdbtnCv.isSelected()) {
            table.setModel(dtmCV);
        }
        resizeColumnWidth();
        resetShownColumns();
    }
    
    private void addBaseToModel(String[] paramNames) {
//        TreeMap<Date, String> map = new TreeMap<Date, String>();
//        for (String f : baseFCSFiles) {
//            map.put(baseFiles.get(f).getRunDate(), f);
//        }
        for (String f : map.values()) {
            fileParamMeanMap.put(f, new HashMap<String, Float>());
            Object[] rowDataM = new Object[paramNames.length + 1];
            Object[] rowDataS = new Object[paramNames.length + 1];
            Object[] rowDataC = new Object[paramNames.length + 1];
            rowDataM[0] = ext.rootOf(f);
            rowDataS[0] = ext.rootOf(f);
            rowDataC[0] = ext.rootOf(f);
            for (int i = 0; i < paramNames.length; i++) {
                FCSDataLoader loader = baseFiles.get(f);
                
                float[] data = loader.getData(paramNames[i], true);
                Float mn = Array.mean(data, true);
                Float sd = Array.stdev(data, true);
                Float cv = 100 * (sd / mn);
                paramMeanLists.get(paramNames[i]).add(mn);
                paramSDLists.get(paramNames[i]).add(sd);
                paramCVLists.get(paramNames[i]).add(cv);
                fileParamMeanMap.get(f).put(paramNames[i], mn);
                rowDataM[i + 1] = mn;
                rowDataS[i + 1] = sd;
                rowDataC[i + 1] = cv;
            }
            dtmMean.addRow(rowDataM);
            dtmSD.addRow(rowDataS);
            dtmCV.addRow(rowDataC);
        }
    }
    
    private void addFilesToModel(DirFile df, String[] paramNames, String removePrep) {
        int colCnt = dtmMean.getColumnCount();
        int rows = dtmMean.getRowCount();
        
        Object[] newRow = new Object[colCnt];
        newRow[0] = df.dir.replaceFirst(removePrep, "");
        dtmMean.addRow(newRow);
        dtmSD.addRow(newRow);
        dtmCV.addRow(newRow);
        boldRows.add(rows);
        rows++;
        
        for (String f : df.files) {
            fileParamMeanMap.put(f, new HashMap<String, Float>());
            Object[] rowDataM = new Object[paramNames.length + 1];
            Object[] rowDataS = new Object[paramNames.length + 1];
            Object[] rowDataC = new Object[paramNames.length + 1];
            rowDataM[0] = ext.rootOf(f);
            rowDataS[0] = ext.rootOf(f);
            rowDataC[0] = ext.rootOf(f);
            for (int i = 0; i < paramNames.length; i++) {
                FCSDataLoader loader = compFiles.get(df.dir + f);

                float[] data = loader.getData(paramNames[i], true);
                Float mn = Array.mean(data, true);
                fileParamMeanMap.get(f).put(paramNames[i], mn);
                Float sd = Array.stdev(data, true);
                Float cv = 100 * (sd / mn);
                rowDataM[i + 1] = mn;
                rowDataS[i + 1] = sd;
                rowDataC[i + 1] = cv;
            }
            dtmMean.addRow(rowDataM);
            dtmSD.addRow(rowDataS);
            dtmCV.addRow(rowDataC);
        }
        if (df.files.length > 0 && df.getAllFiles().length > 0) {
            dtmMean.addRow(new Object[colCnt]);
            dtmSD.addRow(new Object[colCnt]);
            dtmCV.addRow(new Object[colCnt]);
        }
        
        for (DirFile sub : df.subDirs) {
            addFilesToModel(sub, paramNames, removePrep);
        }
        
//        dtmMean.addRow(new Object[dtmMean.getColumnCount()]);
//        dtmSD.addRow(new Object[dtmMean.getColumnCount()]);
//        dtmCV.addRow(new Object[dtmMean.getColumnCount()]);
        
    }
    
    private DirFile loadCompFiles(String dir) {
        File dirFile = new File(dir);
        String[] rootCSV = dirFile.list(csvFileFilter);
        File[] subDirs = dirFile.listFiles(dirFilter);
        
        DirFile df = new DirFile();
        df.dir = ext.verifyDirFormat(dir);
        df.files = rootCSV;
        df.subDirs = new DirFile[subDirs.length];
        for (int i = 0; i < subDirs.length; i++) {
            df.subDirs[i] = loadCompFiles(ext.verifyDirFormat(subDirs[i].getPath()));
        }
        
        return df;
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
