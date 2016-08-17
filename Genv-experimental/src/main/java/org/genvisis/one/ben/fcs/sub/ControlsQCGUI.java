package org.genvisis.one.ben.fcs.sub;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dialog;
import java.awt.EventQueue;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Insets;
import java.awt.Window;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.HierarchyEvent;
import java.awt.event.HierarchyListener;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Properties;
import java.util.TreeSet;

import javax.swing.AbstractAction;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JButton;
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
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.table.DefaultTableColumnModel;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumnModel;

import org.genvisis.cnv.gui.IncludeExcludeGUI;
import org.genvisis.common.Array;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.ext;
import org.genvisis.one.ben.fcs.AbstractPanel2;
import org.genvisis.one.ben.fcs.FCSPlot;
import org.genvisis.one.ben.fcs.sub.MeanCtrlPanel.LabelPresenter;
import org.genvisis.one.ben.fcs.sub.OneDPanel.PLOT_TYPE;

import com.google.common.primitives.Doubles;
import com.google.common.primitives.Floats;

import net.miginfocom.swing.MigLayout;


public class ControlsQCGUI extends JFrame {

  static class CompData {

    // compData.fetchRecentData(params, numOfMostRecentToAdd);
    // String[][] fData = compData.getRecentSources();
    // float[][] pData = compData.getRecentData();

    // private static final FilenameFilter csvFileFilter = new FilenameFilter() {
    // @Override
    // public boolean accept(File dir, String name) {
    // return name.endsWith(".csv");
    // }
    // };
    // private static final FileFilter dirFilter = new FileFilter() {
    // @Override
    // public boolean accept(File pathname) {
    // return pathname.isDirectory();
    // }
    // };

    // public static CompData load(String[] files) {
    // files may not be named properly, but may contain proper data
    // check for duplicated data in files [[allow??]]
    // -- check marks in table to include/exclude data?
    // problems with:
    // - mixed dates between files [split up view? i.e. two lines from
    // one file, one line from another file, another line from first
    // file]
    // --- Option for "Most recent regardless of source file" /
    // "Most recent within source files"

    // - files named properly may not contain the data for which they're
    // named [i.e. named 2016-07-26 but contain data from different
    // date, or contain range of data...]
    // --- should file names be used at all??? (other than as bolded
    // table headers) [No reason to use filenames, as proper data has
    // FCS source file listed]

    // CompData cd = new CompData(files);
    //
    // }


    String[][] source;

    float[][][] data;

    HashMap<String, DataFile> files;

    public CompData(String[] filesToLoad, Logger log) {
      files = new HashMap<String, ControlsQCGUI.DataFile>();
      for (String s : filesToLoad) {
        files.put(s, DataFile.construct(s, log));
      }
    }

    public void fetchRecentData(ArrayList<String> params, String currCtrl, String currPanel,
        int numOfMostRecentToAdd, Logger log) {
      HashMap<CtrlFileMetaData, DataFile> metaMap =
          new HashMap<ControlsQCGUI.CtrlFileMetaData, ControlsQCGUI.DataFile>();
      HashMap<CtrlFileMetaData, float[]> allData =
          new HashMap<ControlsQCGUI.CtrlFileMetaData, float[]>();
      for (DataFile df : files.values()) {
        ArrayList<CtrlFileMetaData> metaData = df.getMetaDataFor(currCtrl, currPanel);
        for (int i = 0; i < metaData.size(); i++) {
          if (allData.containsKey(metaData.get(i))) {
            log.reportError("Error - duplicate data entry found: " + metaData.toString()
                + " in files: " + df.fileName + " || " + metaMap.get(metaData.get(i)).fileName);
          }
          allData.put(metaData.get(i), df.getFileData(metaData.get(i).file));
          metaMap.put(metaData.get(i), df);
        }
      }
      ArrayList<CtrlFileMetaData> sorted =
          new ArrayList<ControlsQCGUI.CtrlFileMetaData>(allData.keySet());
      sorted.sort(CtrlFileMetaData.COMPARATOR);

      ArrayList<String[]> sourceLines = new ArrayList<String[]>();
      ArrayList<float[][]> sourceData = new ArrayList<float[][]>();
      String prevSrc = null;
      ArrayList<String> sourceLine = null;
      ArrayList<float[]> dataLines = null;
      for (int c = 0, count = (numOfMostRecentToAdd == 0 ? sorted.size()
          : Math.min(sorted.size(), numOfMostRecentToAdd)); c < count; c++) {
        CtrlFileMetaData cfmd = sorted.get(c);
        String src = metaMap.get(cfmd).fileName;
        if (sourceLine == null) {
          sourceLine = new ArrayList<String>();
          sourceLine.add(src);
          prevSrc = src;
          dataLines = new ArrayList<float[]>();
          dataLines.add(null);
        } else {
          if (src.equals(prevSrc)) {
            // nuffin
          } else {
            sourceLines.add(sourceLine.toArray(new String[sourceLine.size()]));
            sourceData.add(dataLines.toArray(new float[dataLines.size()][]));
            sourceLine = new ArrayList<String>();
            sourceLine.add(src);
            prevSrc = src;
            dataLines = new ArrayList<float[]>();
            dataLines.add(null);
          }
        }
        sourceLine.add(cfmd.file);
        ArrayList<String> dataParams = metaMap.get(cfmd).allParams;
        float[] fileData = allData.get(cfmd);
        float[] lineData = new float[params.size()];
        for (int i = 0; i < lineData.length; i++) {
          int paramInd = dataParams
              .indexOf(params.get(i).trim().toUpperCase().replaceAll("\\(%\\)", "").trim());
          lineData[i] = paramInd == -1 ? Float.NaN : fileData[paramInd];
        }
        dataLines.add(lineData);
      }
      if (sourceLine != null) {
        sourceLines.add(sourceLine.toArray(new String[sourceLine.size()]));
        sourceData.add(dataLines.toArray(new float[dataLines.size()][]));
      } else {
        log.reportError(
            "Error - no data found for CtrlGroup " + currCtrl + " and Panel " + currPanel);
      }

      source = sourceLines.toArray(new String[sourceLines.size()][]);
      data = sourceData.toArray(new float[sourceData.size()][][]);
    }

    public boolean fileSetEquals(String[] fls) {
      if (fls.length != files.size()) {
        return false;
      }
      for (String s : fls) {
        if (!files.containsKey(s)) {
          return false;
        }
      }
      return true;
    }

    public float[][][] getRecentData() {
      return data;
    }

    public String[][] getRecentSources() {
      return source;
    }

  }
  static class CtrlFileMetaData {
    public static final Comparator<CtrlFileMetaData> COMPARATOR =
        new Comparator<CtrlFileMetaData>() {
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

    public static CtrlFileMetaData parse(String filename, Logger log) {
      CtrlFileMetaData cfmd = new CtrlFileMetaData();
      cfmd.file = filename;
      String[] pts = ext.rootOf(filename, true).split("_");
      if (pts.length != 6) {
        log.reportError("Error - file " + filename + " is an unexpected format!");
        return null;
      }
      if (pts.length > 0) {
        cfmd.fileDateStr = pts[0];
        String[] dtPts = pts[0].split("-");
        try {
          cfmd.fileDate = new GregorianCalendar(Integer.parseInt(dtPts[0]),
              Integer.parseInt(dtPts[1]) - 1, Integer.parseInt(dtPts[2])).getTime();
        } catch (NumberFormatException e) {
          log.reportError("Error - filename " + filename
              + " does not contain a date (YYYY-MM-DD format) as the first token!");
          return null;
        }
      }
      if (pts.length > 1) {
        cfmd.panel = pts[1].toUpperCase();
      }
      cfmd.panel = cfmd.panel.replaceAll("[ \\-_/]", "");
      cfmd.panel = cfmd.panel.replaceAll("PANEL", "");
      cfmd.panel = cfmd.panel.replaceAll("P", "");
      cfmd.panel = cfmd.panel.trim();
      if (pts.length > 2) {
        cfmd.initials = pts[2].toUpperCase().replaceAll("[ \\-_/]", "");
      }
      if (pts.length > 3) {
        cfmd.runGroup = pts[3].toUpperCase().replaceAll("[ \\-_/]", "");
      }
      if (pts.length > 4) {
        cfmd.ctrlGroup = pts[4].toUpperCase().replaceAll("[ \\-_/]", "").replaceAll("CTL", "");
      }
      if (pts.length > 5) {
        cfmd.number = pts[5];
      }
      return cfmd;
    }

    String file;
    String fileDateStr;
    Date fileDate;
    private String panel;
    private String runGroup;
    private String ctrlGroup;

    private String initials;

    private String number;

    private CtrlFileMetaData() {}

    public String getPanel() {
      return panel;
    }

    public boolean isControlGroup(String ctrl) {
      if (ctrl == null && ctrlGroup != null) {
        return false;
      }
      if (ctrlGroup.equals(ctrl.replaceAll("[ \\-_/]", ""))) {
        return true;
      }
      return false;
    }

    public boolean isPanel(String panel2) {
      if (panel2 == null && panel != null) {
        return false;
      }
      if (panel.equals(panel2.replaceAll("[ \\-_/]", ""))) {
        return true;
      }
      return false;
    }
  }
  static class DataFile {

    public static DataFile construct(String dataFile, Logger log) {
      String[][] data = HashVec.loadFileToStringMatrix(dataFile, false, null, false);

      DataFile df = new DataFile();
      df.fileName = dataFile;

      for (int i = 1; i < data[0].length; i++) {
        if (!"".equals(data[0][i])) {
          df.allParams.add(data[0][i].trim().toUpperCase().replaceAll("\\(%\\)", "").trim());
        }
      }
      HashMap<String, Integer> fileInd = new HashMap<String, Integer>();
      for (int i = 1; i < data.length; i++) {
        int ind = ext.indexOfStr(data[i][0], EXCLUDED_ROW_HEADERS, false, false);
        if ("".equals(data[i][0]) || ind >= 0) {
          continue;
        }
        CtrlFileMetaData cfmd = CtrlFileMetaData.parse(data[i][0], log);
        if (cfmd == null) {
          log.reportError("Data on line " + i + " in file " + dataFile + " will be excluded.");
          continue; // couldn't parse filename - skip line;
        }
        df.internalFiles.add(data[i][0]);
        df.internalMetaData.put(data[i][0], cfmd);
        fileInd.put(data[i][0], i);
      }

      df.fileData = new float[df.internalFiles.size()][df.allParams.size()];

      for (int f = 0; f < df.internalFiles.size(); f++) {
        int rowInd = fileInd.get(df.internalFiles.get(f));
        float[] arr = new float[df.allParams.size()];
        for (int i = 0; i < df.allParams.size(); i++) {
          arr[i] = "".equals(data[rowInd][i + 1]) ? Float.NaN
              : Float.parseFloat(data[rowInd][i + 1].replace("%", "").trim());
        }
        df.fileData[f] = arr;
      }

      return df;
    }

    private String fileName;
    private final ArrayList<String> allParams = new ArrayList<String>();
    private final ArrayList<String> internalFiles = new ArrayList<String>();
    private float[][] fileData;
    private final HashMap<String, float[]> paramDataCache = new HashMap<String, float[]>();

    private final HashMap<String, CtrlFileMetaData> internalMetaData =
        new HashMap<String, ControlsQCGUI.CtrlFileMetaData>();

    private DataFile() {}

    @SuppressWarnings("unchecked")
    public ArrayList<String> getAllFiles() {
      return (ArrayList<String>) internalFiles.clone();
    }

    @SuppressWarnings("unchecked")
    public ArrayList<String> getAllParams() {
      return (ArrayList<String>) allParams.clone();
    }

    public String[] getControlGroups() {
      TreeSet<String> ctrlSet = new TreeSet<String>();
      for (CtrlFileMetaData cfmd : internalMetaData.values()) {
        ctrlSet.add(cfmd.ctrlGroup);
      }
      return ctrlSet.toArray(new String[ctrlSet.size()]);
    }

    public float[] getFileData(String internalFile) {
      float[] retArr = null;
      int ind = internalFiles.indexOf(internalFile);
      if (ind > -1) {
        return fileData[ind];
      }
      return retArr;
    }

    public ArrayList<CtrlFileMetaData> getMetaDataFor(String currCtrl, String currPanel) {
      ArrayList<CtrlFileMetaData> retList = new ArrayList<ControlsQCGUI.CtrlFileMetaData>();
      for (CtrlFileMetaData cfmd : internalMetaData.values()) {
        if ((currCtrl == null || cfmd.isControlGroup(currCtrl))
            && (currPanel == null || cfmd.isPanel(currPanel))) {
          retList.add(cfmd);
        }
      }
      return retList;
    }

    public String[] getPanels() {
      TreeSet<String> panelSet = new TreeSet<String>();
      for (CtrlFileMetaData cfmd : internalMetaData.values()) {
        panelSet.add(cfmd.getPanel());
      }
      return panelSet.toArray(new String[panelSet.size()]);
    }

    public float[] getParamData(String param, String ctrl, String panel) {
      float[] retArr = paramDataCache.get(param + "\t" + ctrl + "\t" + panel);
      if (retArr == null) {
        int ind = allParams.indexOf(param.trim().toUpperCase().replaceAll("\\(%\\)", "").trim());
        if (ind > -1) {
          retArr = Matrix.extractColumn(fileData, ind);
          boolean[] incl = new boolean[internalFiles.size()];
          for (int i = 0; i < incl.length; i++) {
            CtrlFileMetaData cfmd = internalMetaData.get(internalFiles.get(i));
            incl[i] = (ctrl == null || cfmd.isControlGroup(ctrl))
                && (panel == null || cfmd.isPanel(panel));

          }
          retArr = Array.subArray(retArr, incl);
          paramDataCache.put(param + "\t" + ctrl + "\t" + panel, retArr);
        }
      }
      return retArr;
    }

  }

  /**
  * 
  */
  private static final long serialVersionUID = 1L;
  private final static String[] EXCLUDED_ROW_HEADERS = {"mean", "average", "sd", "cv", "cv (%)"};
  private static final String PROP_FILE = "controlsQC.properties";

  private static final String PROPKEY_COMPARE_FILES = "COMPARE_FILES";
  private static final String PROPKEY_RANGEFILE = "RANGE_FILE";
  private static final String PROPKEY_COLS = "HIDDEN_COLUMNS";
  private static final String PROPKEY_CTRL = "CONTROL_GROUP";
  private static final String PROPKEY_PANEL = "PANEL_GROUP";
  private static final int TREND_ABOVE_1SD_THRESH = 5;

  private static final int TREND_ABOVE_2SD_THRESH = 2;

  private static final double PCT_OF_EVENTS_DEV_TREND = 0.25; // 1-quarter of events outside of 1SD
                                                              // will result in file being reported

  /**
   * Launch the application.
   */
  public static void main(String[] args) {
    EventQueue.invokeLater(new Runnable() {
      @Override
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

  private JPanel contentPane;
  private JTextField txtFldBaseFile;
  private JButton btnBaseFileSelect;
  private JTable table;
  private JScrollPane scrollPane;
  HashMap<String, ArrayList<Float>> paramMeanLists = new HashMap<String, ArrayList<Float>>();
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
  private DefaultTableModel dtmMean;
  private JButton btnHideshowColumns;

  Color SD1_COLOR = Color.YELLOW;

  Color SD2_COLOR = Color.RED;


  Color ABOVE_3_CV_COLOR = Color.CYAN;

  Color ABOVE_5_CV_COLOR = Color.CYAN.darker();

  CompData compData;
  DataFile baseData;
  private JPanel panel_1;
  private JComboBox<String> comboControl;
  private JSeparator separator;
  private JLabel lblControlGroup;
  private JSeparator separator_1;
  private JLabel lblPanel;
  private JComboBox<String> comboPanel;
  private JSeparator separator_3;
  private JSpinner spinner;
  private JLabel lblOfControls;
  private JSeparator separator_4;
  private JButton btnAddComp;

  private JButton btnRemoveComp;

  private Logger log;

  private JLabel lblWarning;

  private JButton btnMore;

  /**
   * Create the frame.
   */
  public ControlsQCGUI() {
    super("Controls QC");
    setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    setBounds(50, 50, 850, 400);
    log = new Logger();
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
      @Override
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
          loadBaseline(newPath);
          saveProps();
        }
      }
    });
    contentPane.add(btnBaseFileSelect, "cell 2 0,growx");

    lblCompareFiles = new JLabel("<html>Compare <u>D</u>ir/Files:</html>");
    contentPane.add(lblCompareFiles, "cell 0 1,alignx trailing");

    txtFldCompDir = new JTextField();
    txtFldCompDir.setEditable(false);
    contentPane.add(txtFldCompDir, "cell 1 1,growx");
    txtFldCompDir.setColumns(10);

    btnAddComp = new JButton("+");
    btnAddComp.setMnemonic(KeyEvent.VK_O);
    btnAddComp.setMargin(new Insets(0, 6, 0, 6));
    btnAddComp.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent arg0) {
        String curr = txtFldCompDir.getText();
        if (curr.equals("")) {
          curr = txtFldBaseFile.getText();
          if (!"".equals(curr)) {
            curr = ext.parseDirectoryOfFile(curr);
          } else {
            curr = "./";
          }
        } else {
          String[] pts = curr.split(";");
          curr = ext.parseDirectoryOfFile(pts[pts.length - 1]);
        }
        JFileChooser jfc = new JFileChooser(curr);
        jfc.setFileSelectionMode(JFileChooser.FILES_ONLY);
        jfc.setDialogTitle("Select Additional Compare Files");
        jfc.setMultiSelectionEnabled(true);
        int resp = jfc.showOpenDialog(ControlsQCGUI.this);
        if (resp == JFileChooser.APPROVE_OPTION) {
          File[] addlFiles = jfc.getSelectedFiles();
          curr = txtFldCompDir.getText();
          if (!"".equals(curr) && addlFiles.length > 0) {
            curr += ";";
          } else if (addlFiles.length == 0) {
            return;
          }
          StringBuilder path = new StringBuilder(curr);
          for (int i = 0; i < addlFiles.length; i++) {
            String pt = addlFiles[i].getAbsolutePath();
            pt = ext.verifyDirFormat(ext.parseDirectoryOfFile(pt)) + ext.removeDirectoryInfo(pt);
            path.append(pt);
            if (i < addlFiles.length - 1) {
              path.append(";");
            }
          }
          txtFldCompDir.setText(path.toString());
          loadCompare(path.toString().split(";")); // TODO replace split call
          saveProps();
        }
      }
    });

    btnRemoveComp = new JButton("-");
    btnRemoveComp.setMnemonic(KeyEvent.VK_O);
    btnRemoveComp.setMargin(new Insets(0, 6, 0, 6));
    btnRemoveComp.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent arg0) {
        String curr = txtFldCompDir.getText();
        String[] pts;
        if (curr.equals("")) {
          return;
        } else {
          pts = curr.split(";");
        }
        IncludeExcludeGUI dialog =
            new IncludeExcludeGUI(ControlsQCGUI.this, pts, Array.booleanArray(pts.length, true));
        dialog.setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        dialog.pack();
        dialog.setVisible(true);
        int code = dialog.getCloseCode();
        if (code == JOptionPane.OK_OPTION) {
          StringBuilder sb = new StringBuilder();
          boolean[] inc = dialog.getSelected();
          for (int i = 0; i < pts.length; i++) {
            if (inc[i]) {
              sb.append(pts[i]).append(";");
            }
          }
          if (sb.length() > 0) {
            sb.deleteCharAt(sb.length() - 1); // remove trailing semicolon
          }
          txtFldCompDir.setText(sb.toString());
          loadCompare(sb.toString().split(";")); // TODO replace split call
          saveProps();
        }

      }
    });
    contentPane.add(btnRemoveComp, "cell 2 1");

    contentPane.add(btnAddComp, "flowx,cell 2 1");
    contentPane.add(btnRemoveComp, "cell 2 1");

    scrollPane = new JScrollPane();
    scrollPane.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_ALWAYS);
    scrollPane.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS);
    contentPane.add(scrollPane, "cell 0 2 3 1,grow");

    table = new JTable() {

      /**
      * 
      */
      private static final long serialVersionUID = 1L;

      @Override
      public Component prepareRenderer(TableCellRenderer renderer, int row, int column) {
        Component c = super.prepareRenderer(renderer, row, column);

        if (/* column == 0 && */ boldRows.contains(row)) {
          c.setFont(c.getFont().deriveFont(Font.BOLD));
        }
        if (isCellSelected(row, column)) {
          return c;
        }

        Color col = Color.WHITE;
        if (column > 0 && !statRows.contains(row)) {
          Object val = table.getModel().getValueAt(table.convertRowIndexToModel(row),
              table.convertColumnIndexToModel(column));
          if (val != null) {
            if (val instanceof Float) {
              String colNm =
                  (String) table.getModel().getValueAt(0, table.convertColumnIndexToModel(column));
              Float value = (Float) val;
              if (paramMeans.containsKey(colNm)) {
                if (value > (paramMeans.get(colNm) + paramSDs.get(colNm))
                    || value < (paramMeans.get(colNm) - paramSDs.get(colNm))) {
                  col = SD1_COLOR;
                }
                if (value > (paramMeans.get(colNm) + 2 * paramSDs.get(colNm))
                    || value < (paramMeans.get(colNm) - 2 * paramSDs.get(colNm))) {
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
      @Override
      public void mouseClicked(MouseEvent e) {
        if (e.getClickCount() == 2) {
          JTable target = (JTable) e.getSource();
          int row = target.getSelectedRow();
          int column = target.getSelectedColumn();
          if (column == 0) {
            return;
          }
          if (row != 0 && boldRows.contains(row)) {
            return;
          }
          Object o1 = target.getValueAt(row, column);
          if (o1 == null) {
            return;
          }

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
    // table.setCellSelectionEnabled(false);
    table.setColumnModel(new DefaultTableColumnModel() {
      /**
      * 
      */
      private static final long serialVersionUID = 1L;

      @Override
      public void moveColumn(int columnIndex, int newIndex) {
        if (columnIndex == 0 || newIndex == 0) {
          return;
        }
        super.moveColumn(columnIndex, newIndex);
      }
    });
    table.setRowSelectionAllowed(true);
    table.setColumnSelectionAllowed(true);

    table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);

    panel_1 = new JPanel();
    panel_1.setBorder(null);
    contentPane.add(panel_1, "cell 0 3 3 1,grow");
    panel_1.setLayout(new MigLayout("hidemode 3, ins 0", "[][][][][][][][][][][][grow]", "[]"));

    btnHideshowColumns = new JButton("Hide/Show Columns");
    panel_1.add(btnHideshowColumns, "cell 0 0");

    separator = new JSeparator();
    separator.setOrientation(SwingConstants.VERTICAL);
    panel_1.add(separator, "cell 1 0,growy");

    lblControlGroup = new JLabel("Control:");
    panel_1.add(lblControlGroup, "cell 2 0,alignx trailing");

    ActionListener reCalcListener = new ActionListener() {
      @SuppressWarnings("unchecked")
      @Override
      public void actionPerformed(ActionEvent arg0) {
        if (((JComboBox<String>) arg0.getSource()).isPopupVisible()) {
          saveProps();
        }
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

    separator_3 = new JSeparator();
    separator_3.setOrientation(SwingConstants.VERTICAL);
    panel_1.add(separator_3, "cell 7 0,growy");

    lblOfControls = new JLabel("# of Controls (0 = all):");
    panel_1.add(lblOfControls, "cell 8 0");

    spinner = new JSpinner();
    spinner.setModel(new SpinnerNumberModel(new Integer(0), new Integer(0), null, new Integer(1)));
    spinner.addChangeListener(new ChangeListener() {
      @Override
      public void stateChanged(ChangeEvent arg0) {
        reCalcTableData();
      }
    });
    panel_1.add(spinner, "cell 9 0");

    separator_4 = new JSeparator();
    separator_4.setOrientation(SwingConstants.VERTICAL);
    panel_1.add(separator_4, "cell 10 0,growy");

    lblWarning = new JLabel("");
    panel_1.add(lblWarning, "flowx,cell 11 0,alignx right");

    btnMore = new JButton("More");
    btnMore.setVisible(false);
    panel_1.add(btnMore, "cell 11 0,alignx right");

    btnHideshowColumns.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent arg0) {
        if (baseData == null) {
          return;
        }
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
          meanPanel.setAxisXHeight(AbstractPanel2.HEIGHT_X_AXIS - AbstractPanel2.HEIGHT_X_AXIS / 5);
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
        // showMeanPanel(newCol);
      }
    };
    meanCtrlPanel.setChangeListener(prevLst);
    meanCtrlPanel.setPlotChangeListener(plotLst);
    meanPanel.setPlotType(PLOT_TYPE.BOX_PLOT);
    meanPanel.setAxisXHeight(AbstractPanel2.HEIGHT_X_AXIS - AbstractPanel2.HEIGHT_X_AXIS / 2);
    meanPanel.setAxisYWidth(AbstractPanel2.WIDTH_Y_AXIS - AbstractPanel2.WIDTH_Y_AXIS / 3);
    meanPanel.setXAxisLabel("");// pts[0].trim().replaceAll("/", " /\n");
    // meanPanel.setYAxisLabel(meanPanel.plotLabel.split("\\|")[1].trim());
    meanCtrlPanel.setLabelPresenter(new LabelPresenter() {
      @Override
      public String getPresentationView(String label) {
        String[] pts = label.split("\\|")[0].trim().split("/");
        String col = pts[pts.length - 1];
        return col;
      }
    });
    meanFrame.setBounds(FCSPlot.START_X, FCSPlot.START_Y, FCSPlot.START_WIDTH,
        FCSPlot.START_HEIGHT);

    loadProps();
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
      Object[] rowDataM = new Object[fileData.length + 2];
      rowDataM[0] = null;
      rowDataM[1] = ext.rootOf(f);
      for (int i = 0; i < fileData.length; i++) {
        ArrayList<Float> paramMeans = paramMeanLists.get(params.get(i));
        paramMeans.add(fileData[i]);
        rowDataM[i + 2] = fileData[i];
      }
      dtmMean.addRow(rowDataM);
    }
  }

  private void addFilesToModel(CompData compData, int numOfMostRecentToAdd) {
    String currPanel = (String) comboPanel.getSelectedItem();
    if (currPanel.equals("All")) {
      currPanel = null;
    }
    String currCtrl = (String) comboControl.getSelectedItem();
    if (currCtrl.equals("All")) {
      currCtrl = null;
    }
    ArrayList<String> params = baseData.getAllParams();

    compData.fetchRecentData(params, currCtrl, currPanel, numOfMostRecentToAdd, log);
    String[][] fData = compData.getRecentSources();
    float[][][] pData = compData.getRecentData();

    for (int i = 0; i < fData.length; i++) {
      String[] fileSourceAndInternals = fData[i];
      Object[] rowDataM = new Object[params.size() + 2];
      rowDataM[0] = null;
      rowDataM[1] = fileSourceAndInternals[0];
      dtmMean.addRow(rowDataM.clone());
      boldRows.add(dtmMean.getRowCount() - 1);

      for (int f = 1; f < fileSourceAndInternals.length; f++) {
        rowDataM = new Object[params.size() + 2];
        rowDataM[0] = null; // TODO check-boxes for inclusion (should be in model/view, not
                            // underlying data structure [i.e., this class, not the data classes])
        rowDataM[1] = fileSourceAndInternals[f];
        for (int p = 0; p < pData[i][f].length; p++) {
          rowDataM[p + 2] = pData[i][f][p];
        }
        dtmMean.addRow(rowDataM.clone());
      }
    }
  }



  private void checkWarnings() {
    if (baseData != null && compData != null) {
      final ArrayList<String> warnings = new ArrayList<String>();

      ArrayList<String> params = baseData.getAllParams();
      float[][][] currData = compData.getRecentData();

      HashMap<String, Integer> paramTrendsAbove1 = new HashMap<String, Integer>();
      HashMap<String, Integer> paramTrendsAbove2 = new HashMap<String, Integer>();

      for (int p = 0; p < params.size(); p++) {
        int trendAbove1 = 0;
        int trendAbove2 = 0;
        float mean = paramMeans.get(params.get(p));
        float sd = paramSDs.get(params.get(p));

        for (float[][] element : currData) {
          for (int s = 1; s < element.length; s++) {
            float param = element[s][p];
            if (param < mean - sd || param > mean + sd) {
              trendAbove1++;
              if (param < mean - 2 * sd || param > mean + 2 * sd) {
                trendAbove2++;
                if (trendAbove2 >= TREND_ABOVE_2SD_THRESH) {
                  paramTrendsAbove2.put(params.get(p), trendAbove2);
                }
              } else {
                trendAbove2 = 0;
                if (trendAbove1 >= TREND_ABOVE_1SD_THRESH) {
                  paramTrendsAbove1.put(params.get(p), trendAbove1);
                }
              }
            } else {
              trendAbove1 = 0;
              trendAbove2 = 0;
            }

          }
        }
      }
      int max = 0;
      for (Integer val : paramTrendsAbove1.values()) {
        max = Math.max(max, val);
      }
      for (Integer val : paramTrendsAbove2.values()) {
        max = Math.max(max, val);
      }

      String maxWarn = null;
      String maxWarnParam = null;
      if (!paramTrendsAbove1.isEmpty()) {
        for (String s : paramTrendsAbove1.keySet()) {
          String[] p = s.split("\\|")[0].trim().split("/");
          String warn = paramTrendsAbove1.get(s) + "-count trend in parameter " + p[p.length - 1]
              + " greater than 1SD from the mean.";
          warnings.add(warn);
          if (paramTrendsAbove1.get(s) == max) {
            maxWarn = warn;
            maxWarnParam = s;
          }
        }
      }
      if (!paramTrendsAbove2.isEmpty()) {
        for (String s : paramTrendsAbove2.keySet()) {
          String[] p = s.split("\\|")[0].trim().split("/");
          String warn = paramTrendsAbove2.get(s) + "-count trend in parameter " + p[p.length - 1]
              + " greater than 2SD from the mean.";
          warnings.add(warn);
          if (paramTrendsAbove2.get(s) == max) {
            maxWarn = warn;
            maxWarnParam = s;
          }
        }
      }

      String[][] srcs = compData.getRecentSources();

      for (int i = 0; i < currData.length; i++) {
        for (int s = 1; s < currData[i].length; s++) {
          int cnt = 0;
          for (int p = 0; p < params.size(); p++) {
            float mean = paramMeans.get(params.get(p));
            float sd = paramSDs.get(params.get(p));
            float param = currData[i][s][p];
            if (param < mean - sd || param > mean + sd) {
              cnt++;
            }
          }

          if (cnt >= params.size() * PCT_OF_EVENTS_DEV_TREND) {
            warnings.add("Source " + srcs[i][s] + " has " + cnt
                + " events greater than 1SD from parameter means.");
          }

        }
      }

      if (warnings.size() > 0) {
        lblWarning.setIcon(UIManager.getIcon("OptionPane.warningIcon"));
        lblWarning.setText(maxWarn);
        final String finalMaxWarnParam = maxWarnParam;
        lblWarning.addMouseListener(new MouseAdapter() {
          @Override
          public void mouseClicked(MouseEvent e) {
            if (e.getClickCount() == 2) {
              showMeanPanel(finalMaxWarnParam);
            }
          }
        });

        if (warnings.size() > 1) {
          btnMore.setVisible(true);
          btnMore.setAction(new AbstractAction() {
            /**
            * 
            */
            private static final long serialVersionUID = 1L;

            @Override
            public void actionPerformed(ActionEvent e) {
              final JPanel msgPane = new JPanel(new MigLayout("", "", ""));
              JScrollPane scroll = new JScrollPane(msgPane);
              for (int i = 0; i < warnings.size(); i++) {
                msgPane.add(new JLabel("" + warnings.get(i)), "cell 0 " + i);
              }
              msgPane.addHierarchyListener(new HierarchyListener() {
                @Override
                public void hierarchyChanged(HierarchyEvent e) {
                  Window window = SwingUtilities.getWindowAncestor(msgPane);
                  if (window instanceof Dialog) {
                    Dialog dialog = (Dialog) window;
                    if (!dialog.isResizable()) {
                      dialog.setResizable(true);
                    }
                  }
                }
              });
              JOptionPane.showMessageDialog(ControlsQCGUI.this, scroll, "Trend Warnings",
                  JOptionPane.WARNING_MESSAGE);
            }
          });
          btnMore.setText("All");
        }
      }
    }
  }

  protected void loadBaseline(String fileOrDir) {
    if (baseData != null && baseData.fileName.equals(fileOrDir)) {
      return; // same as already loaded
    }
    baseData = DataFile.construct(fileOrDir, log);
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
    reCalcTableData();
  }

  protected void loadCompare(String[] files) {
    if (files == null || files.length == 0 || (files.length == 1 && "".equals(files[0]))) {
      if (compData != null) {
        compData = null;
        reCalcTableData();
      }
      return;
    }
    if (compData != null && compData.fileSetEquals(files)) {
      return; // same as already loaded
    }
    compData = new CompData(files, log);
    reCalcTableData();
  }

  private void loadProps() {
    Properties props = new Properties();
    InputStream is = null;

    try {
      File f = new File(PROP_FILE);
      if (!f.exists()) {
        return;
      }
      is = new FileInputStream(f);
      props.load(is);
      String base = props.getProperty(PROPKEY_RANGEFILE, "");
      String comp = props.getProperty(PROPKEY_COMPARE_FILES, "");
      String colsTemp = props.getProperty(PROPKEY_COLS, "");
      String selCtrl = props.getProperty(PROPKEY_CTRL, "");
      String selPnl = props.getProperty(PROPKEY_PANEL, "");
      String[] cols = colsTemp.split(";");
      String[] compAll = comp.split(";");

      if (!base.equals("")) {
        txtFldBaseFile.setText(base);
        loadBaseline(base);
      }
      if (!comp.equals("")) {
        txtFldCompDir.setText(comp);
        loadCompare(compAll);
      }
      hiddenCols.clear();
      for (String c : cols) {
        hiddenCols.add(c);
      }
      if (selCtrl != null && !"".equals(selCtrl)) {
        comboControl.setSelectedItem(selCtrl);
      }
      if (selPnl != null && !"".equals(selPnl)) {
        comboPanel.setSelectedItem(selPnl);
      }
      reCalcTableData();
    } catch (Exception e) {
      e.printStackTrace();
      is = null;
    }
  }

  private void reCalcTableData() {
    if (baseData == null) {
      return;
    }

    TreeSet<String> paramSet = new TreeSet<String>();

    paramSet.addAll(baseData.getAllParams()); // don't have to worry about panel/ctrl for params

    String[] paramNames = paramSet.toArray(new String[paramSet.size()]);
    String[] colNames = new String[paramNames.length + 2];
    colNames[0] = "";
    colNames[1] = "";
    for (int i = 0; i < paramNames.length; i++) {
      String col = paramNames[i];
      String[] pts = col.split("\\|")[0].trim().split("/");
      col = pts[pts.length - 1];
      colNames[i + 2] = col;
    }
    for (String p : paramNames) {
      paramMeanLists.put(p, new ArrayList<Float>());
    }

    dtmMean = new DefaultTableModel(colNames, 0) {
      /**
      * 
      */
      private static final long serialVersionUID = 1L;

      @Override
      public boolean isCellEditable(int row, int column) {
        return false;
      }
    };

    boldRows.clear();

    String[] firstRow = Array.addStrToArray("", paramNames.clone(), 0);
    firstRow = Array.addStrToArray("", firstRow, 0);
    firstRow[0] = "Include?";
    firstRow[1] = "Source";
    dtmMean.addRow(firstRow);
    boldRows.add(0);

    addDataToModel(baseData);

    dtmMean.addRow(new Object[colNames.length]);
    statRows.add(dtmMean.getRowCount());

    Object[] meanRow = new Object[colNames.length];
    meanRow[0] = null;
    meanRow[1] = "Mean of means";
    for (int i = 2; i < colNames.length; i++) {
      String colNm = paramNames[i - 2];
      if (paramMeanLists.containsKey(colNm)) {
        Float mn = Array.mean(Floats.toArray(paramMeanLists.get(colNm)), true);
        paramMeans.put(colNm, mn);
        meanRow[i] = mn;
      }
    }
    dtmMean.addRow(meanRow);
    statRows.add(dtmMean.getRowCount());

    Object[] sdRow = new Object[colNames.length];
    sdRow[0] = null;
    sdRow[1] = "StdDev of means";
    for (int i = 2; i < colNames.length; i++) {
      String colNm = paramNames[i - 2];
      if (paramMeanLists.containsKey(colNm)) {
        Float sd = Array.stdev(paramMeanLists.get(colNm).toArray(new Float[0]), true);
        paramSDs.put(colNm, sd);
        sdRow[i] = sd;
      }
    }
    dtmMean.addRow(sdRow);
    statRows.add(dtmMean.getRowCount());

    Object[] cvRow = new Object[colNames.length];
    cvRow[0] = null;
    cvRow[1] = "cV ( = 100 * SD / Mean)";
    for (int i = 2; i < colNames.length; i++) {
      String colNm = paramNames[i - 2];
      if (paramMeanLists.containsKey(colNm)) {
        Float cv = 100 * (paramSDs.get(colNm) / paramMeans.get(colNm));
        paramCVs.put(colNm, cv);
        cvRow[i] = cv;
      }
    }
    dtmMean.addRow(cvRow);
    statRows.add(dtmMean.getRowCount());

    Object[] rgRow = new Object[colNames.length];
    rgRow[0] = null;
    rgRow[1] = "Mean - 1SD";
    for (int i = 2; i < colNames.length; i++) {
      String colNm = paramNames[i - 2];
      if (paramMeans.containsKey(colNm)) {
        rgRow[i] = paramMeans.get(colNm) - paramSDs.get(colNm);
      }
    }
    dtmMean.addRow(rgRow);
    statRows.add(dtmMean.getRowCount());

    rgRow = new Object[colNames.length];
    rgRow[0] = null;
    rgRow[1] = "Mean + 1SD";
    for (int i = 2; i < colNames.length; i++) {
      String colNm = paramNames[i - 2];
      if (paramMeans.containsKey(colNm)) {
        rgRow[i] = paramMeans.get(colNm) + paramSDs.get(colNm);
      }
    }
    dtmMean.addRow(rgRow);
    statRows.add(dtmMean.getRowCount());

    dtmMean.addRow(new Object[colNames.length]);

    if (compData != null) {
      int numRecent = (Integer) spinner.getValue();
      addFilesToModel(compData, numRecent);
    }

    checkWarnings();

    table.setModel(dtmMean);
    table.revalidate();
    table.repaint();
    resizeColumnWidth();
    resetShownColumns();
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

  private void resizeColumnWidth() {
    final TableColumnModel columnModel = table.getColumnModel();
    FontMetrics fm = table.getFontMetrics(table.getFont());
    for (int column = 0; column < table.getColumnCount(); column++) {
      int width = fm.stringWidth(columnModel.getColumn(column).getHeaderValue().toString()) + 20; // Min
                                                                                                  // width
      for (int row = 1; row < table.getRowCount(); row++) {
        TableCellRenderer renderer = table.getCellRenderer(row, column);
        Component comp = table.prepareRenderer(renderer, row, column);
        width = Math.max(comp.getPreferredSize().width + 1, width);
      }
      columnModel.getColumn(column).setPreferredWidth(width);
    }
  }

  private void saveProps() {
    try {
      Properties props = new Properties();
      props.setProperty(PROPKEY_RANGEFILE, txtFldBaseFile.getText());
      props.setProperty(PROPKEY_COMPARE_FILES, txtFldCompDir.getText());
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
      props.setProperty(PROPKEY_CTRL, (String) comboControl.getSelectedItem());
      props.setProperty(PROPKEY_PANEL, (String) comboPanel.getSelectedItem());
      File f = new File(PROP_FILE);
      OutputStream out = new FileOutputStream(f);
      props.store(out, "");
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  private void showMeanPanel(String col) {
    // Get all means (yData), get Files&Dates (xData), get Mean/SD (meanSD)
    double[] yDataBase;
    String[] baseLbls;

    String currPanel = (String) comboPanel.getSelectedItem();
    if (currPanel.equals("All")) {
      currPanel = null;
    }
    String currCtrl = (String) comboControl.getSelectedItem();
    if (currCtrl.equals("All")) {
      currCtrl = null;
    }

    ArrayList<CtrlFileMetaData> metaData = baseData.getMetaDataFor(currCtrl, currPanel);
    HashMap<CtrlFileMetaData, Integer> metaIndMap =
        new HashMap<ControlsQCGUI.CtrlFileMetaData, Integer>();
    for (int i = 0; i < metaData.size(); i++) {
      metaIndMap.put(metaData.get(i), i);
    }
    metaData.sort(CtrlFileMetaData.COMPARATOR);

    yDataBase = new double[metaData.size()];
    baseLbls = new String[metaData.size()];
    int ind = 0;
    float[] data = baseData.getParamData(col, currCtrl, currPanel);
    for (CtrlFileMetaData cfm : metaData) {
      yDataBase[ind] = data[metaIndMap.get(cfm)];
      baseLbls[ind] = cfm.file;
      ind++;
    }

    ArrayList<String> compLbls = new ArrayList<String>();
    ArrayList<Double> compDbls = new ArrayList<Double>();
    if (compData != null) {
      String[][] sources = compData.getRecentSources();
      float[][][] compDataArr = compData.getRecentData();
      ind = baseData.allParams.indexOf(col);

      for (int i = 0; i < sources.length; i++) {
        for (int s = 1; s < sources[i].length; s++) {
          compLbls.add(sources[i][s]);
          compDbls.add(Double.valueOf(compDataArr[i][s][ind]));
        }
      }
    }

    ArrayList<String> cols = new ArrayList<String>();
    ind = 0;
    for (int i = 2; i < table.getColumnModel().getColumnCount(); i++) {
      String hdr = (String) table.getValueAt(0, i);
      cols.add(hdr);
      if (hdr.equals(col)) {
        ind = i;
      }
    }

    meanCtrlPanel.setColumns(cols, ind - 2);

    meanPanel.setXAxisLabel("");// pts[0].trim().replaceAll("/", " /\n");
    meanPanel.setYAxisLabel(col.split("\\|")[1].trim());
    meanPanel.setData(col, new String[][] {baseLbls, compLbls.toArray(new String[compLbls.size()])},
        new double[][] {yDataBase, Doubles.toArray(compDbls)});
    meanPanel.paintAgain();
    meanFrame.setTitle("jFlow - Control QC - Overall Mean/SD - " + col);
    meanFrame.setVisible(true);
  }


}


