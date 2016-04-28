package one.ben.fcs;

import java.awt.Color;
import java.awt.Component;
import java.awt.EventQueue;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FileFilter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeMap;
import java.util.TreeSet;

import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.ScrollPaneConstants;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.table.DefaultTableColumnModel;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumnModel;
import javax.xml.parsers.ParserConfigurationException;

import net.miginfocom.swing.MigLayout;
import one.ben.fcs.FCSDataLoader.DATA_SET;
import one.ben.fcs.gating.GateFileReader;
import one.ben.fcs.gating.GatingStrategy;

import org.xml.sax.SAXException;

import common.Array;
import common.ext;


public class RainbowTestGUI extends JFrame {

    private JPanel contentPane;
    private JTextField txtFldFileDir;
    private JTextField txtFldGatingFile;
    private JButton btnFileDirSelect;
    private JButton btnGatingFileSelect;
    private JTable table;
    private JScrollPane scrollPane;

    private String[] gateFileExts = {"wsp", "wspt"};
    
    private HashMap<String, FCSDataLoader> files = new HashMap<String, FCSDataLoader>();
    
    private JButton btnEditColumns;
    
    HashMap<String, ArrayList<Float>> paramMeanLists = new HashMap<String, ArrayList<Float>>();
    HashSet<Integer> boldRows = new HashSet<Integer>();
    HashMap<String, Float> paramMeans = new HashMap<String, Float>();
    HashMap<String, Float> paramSDs = new HashMap<String, Float>();
    
    Color ABOVE_1SD_COLOR = Color.RED;
    Color BELOW_1SD_COLOR = Color.RED;
    

    /**
     * Create the frame.
     */
    public RainbowTestGUI() {
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setBounds(50, 50, 850, 400);
        contentPane = new JPanel();
        contentPane.setBorder(null);
        setContentPane(contentPane);
        contentPane.setLayout(new MigLayout("ins 7 7 3 7", "[][grow][]", "[][][grow][][]"));
        
        JLabel lblFileDir = new JLabel("File Dir:");
        contentPane.add(lblFileDir, "cell 0 0,alignx trailing");
        
        txtFldFileDir = new JTextField();
        txtFldFileDir.setEditable(false);
        contentPane.add(txtFldFileDir, "cell 1 0,growx");
        txtFldFileDir.setColumns(10);
        
        Insets btnInsets = new Insets(0, 14, 0, 14);
        
        btnFileDirSelect = new JButton(">");
        btnFileDirSelect.setMargin(btnInsets);
        btnFileDirSelect.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent arg0) {
                String curr = txtFldFileDir.getText();
                if (curr.equals("")) {
                    curr = "./";
                }
                JFileChooser jfc = new JFileChooser(curr);
                jfc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
                jfc.setDialogTitle("Select FCS File Directory");
                jfc.setMultiSelectionEnabled(false);
                int resp = jfc.showOpenDialog(RainbowTestGUI.this);
                if (resp == JFileChooser.APPROVE_OPTION) {
                    String newPath = ext.verifyDirFormat(jfc.getSelectedFile().getAbsolutePath());
                    txtFldFileDir.setText(newPath);
                    loadFCSDir(newPath);
                }
            }
        });
        contentPane.add(btnFileDirSelect, "cell 2 0");
        
        JLabel lblGatingFile = new JLabel("Gating File:");
        contentPane.add(lblGatingFile, "cell 0 1,alignx trailing");
        
        txtFldGatingFile = new JTextField();
        txtFldGatingFile.setEditable(false);
        contentPane.add(txtFldGatingFile, "cell 1 1,growx");
        txtFldGatingFile.setColumns(10);
        
        btnGatingFileSelect = new JButton(">");
        btnGatingFileSelect.setMargin(btnInsets);
        btnGatingFileSelect.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                String curr = txtFldFileDir.getText();
                if (curr.equals("")) {
                    curr = "./";
                }
                JFileChooser jfc = new JFileChooser(curr);
                jfc.setFileSelectionMode(JFileChooser.FILES_ONLY);
                jfc.setFileFilter(new FileNameExtensionFilter("Gating-ML Files", gateFileExts));
                jfc.setDialogTitle("Select Gating-ML File");
                jfc.setMultiSelectionEnabled(false);
                int resp = jfc.showOpenDialog(RainbowTestGUI.this);
                if (resp == JFileChooser.APPROVE_OPTION) {
                    String newPath = jfc.getSelectedFile().getAbsolutePath();
                    txtFldGatingFile.setText(newPath);
                    setGateFile(newPath);
                }
            }
        });
        contentPane.add(btnGatingFileSelect, "cell 2 1");
        
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
                if (column > 0 && row < getRowCount() - 4) {
                    Object val = table.getModel().getValueAt(table.convertRowIndexToModel(row), table.convertColumnIndexToModel(column));
                    if (val != null) {
                        if (val instanceof Float) {
                            String colNm = table.getModel().getColumnName(table.convertColumnIndexToModel(column));
                            Float value = (Float) val;
                            if (paramMeans.containsKey(colNm)) {
                                if (value > (paramMeans.get(colNm) + paramSDs.get(colNm))) {
                                    col = ABOVE_1SD_COLOR;
                                } else if (value < (paramMeans.get(colNm) - paramSDs.get(colNm))) {
                                    col = BELOW_1SD_COLOR;
                                }
                            }
                        }
                    } 
                }
                c.setBackground(col);
                
                return c;
            }
        };
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
        
        btnEditColumns = new JButton("Edit Columns");
        btnEditColumns.setVisible(false);
        contentPane.add(btnEditColumns, "cell 0 3");
        
        setGateFile("F:\\Flow\\rainbow\\URB.wspt");
    }
    
    private void resizeColumnWidth(JTable table) {
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
    
    private void setGateFile(String filePath) {
        try {
            GatingStrategy gateStrat = GateFileReader.readGateFile(filePath);
            
            
            
        } catch (ParserConfigurationException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (SAXException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
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
    

    private final FilenameFilter fcsFileFilter = new FilenameFilter() {
        @Override
        public boolean accept(File dir, String name) {
            return name.endsWith(".fcs");
        }
    };
    private final FileFilter dirFilter = new FileFilter() {
        @Override
        public boolean accept(File pathname) {
            return pathname.isDirectory();
        }
    };
    
    DirFile dirStruct;
    protected void loadFCSDir(String dir) {
        if (dirStruct != null && dirStruct.dir.equals(ext.verifyDirFormat(dir))) {
            return; // exact same as already loaded
        }
        dirStruct = load(dir);
        this.files.clear();
        reCalcTableData();
    }
    
    private void reCalcTableData() {
        String[] allFilesFullPaths = dirStruct.getAllFiles();
        
        TreeMap<Date, String> dateMap = new TreeMap<Date, String>();
        TreeSet<String> paramSet = new TreeSet<String>();
        
        for (String f : allFilesFullPaths) {
            if (!this.files.containsKey(f)) {
                this.files.put(f, loadFCSFile(f));
            }
            dateMap.put(files.get(f).lastModified, f);
            paramSet.addAll(files.get(f).getAllDisplayableNames(DATA_SET.COMPENSATED));
        }
        String[] paramNames = paramSet.toArray(new String[paramSet.size()]);
        String[] colNames = Array.addStrToArray("", paramNames, 0);
        for (String p : paramNames) {
            paramMeanLists.put(p, new ArrayList<Float>());
        }
        
        DefaultTableModel dtm = new DefaultTableModel(colNames, 0) {
            @Override
            public boolean isCellEditable(int row, int column) {
                return false;
            }
        };
        
        boldRows.clear();
        
        String[] firstRow = colNames.clone();
        firstRow[0] = "Source";
        dtm.addRow(firstRow);
        addFilesToModel(dirStruct, paramNames, dtm, dirStruct.dir);
        
        Object[] meanRow = new Object[dtm.getColumnCount()];
        meanRow[0] = "Mean";
        for (int i = 1; i < dtm.getColumnCount(); i++) {
            String colNm = dtm.getColumnName(i);
            if (paramMeanLists.containsKey(colNm)) {
                Float mn = Array.mean(paramMeanLists.get(colNm).toArray(new Float[0]));
                paramMeans.put(colNm, mn);
                meanRow[i] = mn;
            }
        }
        dtm.addRow(meanRow);
        
        Object[] sdRow = new Object[dtm.getColumnCount()];
        sdRow[0] = "StdDev";
        for (int i = 1; i < dtm.getColumnCount(); i++) {
            String colNm = dtm.getColumnName(i);
            if (paramMeanLists.containsKey(colNm)) {
                Float sd = Array.stdev(paramMeanLists.get(colNm).toArray(new Float[0]), false);
                paramSDs.put(colNm, sd);
                sdRow[i] = sd;
            }
        }
        dtm.addRow(sdRow);

        Object[] cvRow = new Object[dtm.getColumnCount()];
        cvRow[0] = "cV ( = 100 * SD / Mean)";
        for (int i = 1; i < dtm.getColumnCount(); i++) {
            String colNm = dtm.getColumnName(i);
            if (paramMeanLists.containsKey(colNm)) {
                cvRow[i] = 100 * (paramSDs.get(colNm) / paramMeans.get(colNm));
            }
        }
        dtm.addRow(cvRow);
        
        Object[] rgRow = new Object[dtm.getColumnCount()];
        rgRow[0] = "Mean - 1SD";
        for (int i = 1; i < dtm.getColumnCount(); i++) {
            String colNm = dtm.getColumnName(i);
            if (paramMeans.containsKey(colNm)) {
                rgRow[i] = paramMeans.get(colNm) - paramSDs.get(colNm);
            }
        }
        dtm.addRow(rgRow);
        
        rgRow = new Object[dtm.getColumnCount()];
        rgRow[0] = "Mean + 1SD";
        for (int i = 1; i < dtm.getColumnCount(); i++) {
            String colNm = dtm.getColumnName(i);
            if (paramMeans.containsKey(colNm)) {
                rgRow[i] = paramMeans.get(colNm) + paramSDs.get(colNm);
            }
        }
        dtm.addRow(rgRow);
        
        table.setModel(dtm);
        resizeColumnWidth(table);
    }
    
    private void addFilesToModel(DirFile df, String[] paramNames, DefaultTableModel dtm, String removePrep) {
        int colCnt = dtm.getColumnCount();
        int rows = dtm.getRowCount();
        
        Object[] newRow = new Object[colCnt];
        newRow[0] = df.dir.replaceFirst(removePrep, "");
        dtm.addRow(newRow);
        boldRows.add(rows);
        rows++;
        
        for (String f : df.files) {
            Object[] rowData = new Object[paramNames.length + 1];
            rowData[0] = ext.rootOf(f);
            for (int i = 0; i < paramNames.length; i++) {
                // TODO apply Gating here
                Float mn = Array.mean(files.get(df.dir + f).getData(paramNames[i], true));
                paramMeanLists.get(paramNames[i]).add(mn);
                rowData[i + 1] = mn;
            }
            dtm.addRow(rowData);
        }
        if (df.files.length > 0 && df.getAllFiles().length > 0) {
            dtm.addRow(new Object[colCnt]);
        }
        
        for (DirFile sub : df.subDirs) {
            addFilesToModel(sub, paramNames, dtm, removePrep);
        }
        
        dtm.addRow(new Object[dtm.getColumnCount()]);
        
    }
    
    private DirFile load(String dir) {
        File dirFile = new File(dir);
        String[] rootFCS = dirFile.list(fcsFileFilter);
        File[] subDirs = dirFile.listFiles(dirFilter);
        
        DirFile df = new DirFile();
        df.dir = ext.verifyDirFormat(dir);
        df.files = rootFCS;
        df.subDirs = new DirFile[subDirs.length];
        for (int i = 0; i < subDirs.length; i++) {
            df.subDirs[i] = load(ext.verifyDirFormat(subDirs[i].getPath()));
        }
        
        return df;
    }
    
    private FCSDataLoader loadFCSFile(String file) {
        FCSDataLoader fdl = new FCSDataLoader();
        try {
            fdl.loadData(file);
        } catch (IOException e) {
            e.printStackTrace();
            // TODO handle this
            return null;
        }
        return fdl;
    }

    /**
     * Launch the application.
     */
    public static void main(String[] args) {
        EventQueue.invokeLater(new Runnable() {
            public void run() {
                try {
                    RainbowTestGUI frame = new RainbowTestGUI();
                    frame.setVisible(true);
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        });
    }
    
    
}
