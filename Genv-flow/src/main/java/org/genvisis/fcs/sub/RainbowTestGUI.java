package org.genvisis.fcs.sub;

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
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.io.FileFilter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
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
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
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
import javax.swing.UIManager;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.table.DefaultTableColumnModel;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumnModel;
import javax.xml.parsers.ParserConfigurationException;

import org.genvisis.fcs.AbstractPanel2;
import org.genvisis.fcs.FCSDataLoader;
import org.genvisis.fcs.JFlow;
import org.genvisis.fcs.AbstractPanel2.PLOT_TYPE;
import org.genvisis.fcs.FCSDataLoader.DATA_SET;
import org.genvisis.fcs.gating.Gate;
import org.genvisis.fcs.gating.GateFileReader;
import org.genvisis.fcs.gating.Gating;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.SerializedFiles;
import org.pankratzlab.common.ext;
import org.pankratzlab.utils.gui.IncludeExcludeGUI;
import org.xml.sax.SAXException;

import com.google.common.primitives.Doubles;

import net.miginfocom.swing.MigLayout;

public class RainbowTestGUI extends JFrame {

  /** */
  private static final long serialVersionUID = 1L;

  private JPanel contentPane;
  private JTextField txtFldBaseDir;
  private JTextField txtFldGatingFile;
  private JButton btnBaseDirSelect;
  private JButton btnGatingFileSelect;
  private JTable table;
  private JScrollPane scrollPane;

  private final String[] gateFileExts = {"wsp", "wspt"};

  private final HashMap<String, FCSDataLoader> baseFiles = new HashMap<>();
  private final HashMap<String, FCSDataLoader> compFiles = new HashMap<>();

  HashMap<String, ArrayList<Double>> paramMeanLists = new HashMap<>();
  HashMap<String, ArrayList<Double>> paramSDLists = new HashMap<>();
  HashMap<String, ArrayList<Double>> paramCVLists = new HashMap<>();
  HashMap<String, HashMap<String, Double>> fileParamMeanMap = new HashMap<>();
  HashSet<Integer> boldRows = new HashSet<>();
  HashSet<Integer> statRows = new HashSet<>();
  HashMap<String, Double> paramMeans = new HashMap<>();
  HashMap<String, Double> paramSDs = new HashMap<>();
  HashMap<String, Double> paramCVs = new HashMap<>();

  TreeSet<String> paramNames = new TreeSet<>();
  HashSet<String> hiddenCols = new HashSet<>();

  JFrame meanFrame = new JFrame("Genvisis - FCS Overall Mean/SD");
  OneDPanel meanPanel = new OneDPanel();
  MeanCtrlPanel meanCtrlPanel = new MeanCtrlPanel();
  Gating gateStrat;
  String baseDir;
  String[] baseFCSFiles;
  DirFile dirStruct;
  private JLabel lblCompareFcsDir;
  private JTextField txtFldCompDir;
  private JButton btnCompDirSelect;
  private JRadioButton rdbtnMean;
  private JRadioButton rdbtnSd;
  private JRadioButton rdbtnCv;
  private final ButtonGroup buttonGroup = new ButtonGroup();

  protected void loadFCSDir(String dir, boolean baseline) {
    if (baseline) {
      if (baseDir != null && baseDir.equals(dir)) {
        return; // same as already loaded
      }
      baseDir = dir;
      loadBase(dir);
      baseFiles.clear();
    } else {
      if (dirStruct != null && dirStruct.dir.equals(ext.verifyDirFormat(dir))) {
        return; // same as already loaded
      }
      dirStruct = load(dir);
      compFiles.clear();
    }
    reCalcTableData();
  }

  private DefaultTableModel dtmMean;
  private DefaultTableModel dtmSD;
  private DefaultTableModel dtmCV;
  private JSeparator separator;
  private JRadioButton rdbtnUngated;
  private JRadioButton rdbtnGated;
  private final ButtonGroup buttonGroup_1 = new ButtonGroup();
  private JButton button;
  private JSeparator separator_1;
  private JButton btnHideshowColumns;
  Color SD1_COLOR = Color.YELLOW;
  Color SD2_COLOR = Color.RED;
  Color ABOVE_3_CV_COLOR = Color.CYAN;
  Color ABOVE_5_CV_COLOR = Color.CYAN.darker();

  /** Create the frame. */
  public RainbowTestGUI() {
    super("Rainbow Bead Testing");
    setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
    setBounds(50, 50, 850, 400);

    addWindowListener(
        new WindowAdapter() {

          @Override
          public void windowClosing(WindowEvent e) {
            RainbowTestGUI.this.setVisible(false);
            RainbowTestGUI.this.saveCache();
            RainbowTestGUI.this.dispose();
            System.exit(0);
          }
        });

    contentPane = new JPanel();
    contentPane.setBorder(null);
    setContentPane(contentPane);
    contentPane.setLayout(
        new MigLayout("ins 7 7 3 7,hidemode 3", "[][grow][]", "[][][][grow][][]"));

    JLabel lblFileDir = new JLabel("<html><u>B</u>aseline FCS Dir:</html>");
    contentPane.add(lblFileDir, "cell 0 0,alignx trailing");

    txtFldBaseDir = new JTextField();
    txtFldBaseDir.setEditable(false);
    contentPane.add(txtFldBaseDir, "cell 1 0,growx");
    txtFldBaseDir.setColumns(10);

    Insets btnInsets = new Insets(0, 14, 0, 14);

    btnBaseDirSelect = new JButton(">");
    btnBaseDirSelect.setMargin(btnInsets);
    btnBaseDirSelect.setMnemonic(KeyEvent.VK_B);
    btnBaseDirSelect.addActionListener(
        new ActionListener() {

          @Override
          public void actionPerformed(ActionEvent arg0) {
            String curr = txtFldBaseDir.getText();
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
              txtFldBaseDir.setText(newPath);
              loadFCSDir(newPath, true);
              saveProps();
            }
          }
        });
    contentPane.add(btnBaseDirSelect, "cell 2 0");

    lblCompareFcsDir = new JLabel("<html>Compare FCS <u>D</u>ir:</html>");
    contentPane.add(lblCompareFcsDir, "cell 0 1,alignx trailing");

    txtFldCompDir = new JTextField();
    txtFldCompDir.setEditable(false);
    contentPane.add(txtFldCompDir, "cell 1 1,growx");
    txtFldCompDir.setColumns(10);

    btnCompDirSelect = new JButton(">");
    btnCompDirSelect.setMargin(new Insets(0, 14, 0, 14));
    btnCompDirSelect.setMnemonic(KeyEvent.VK_D);
    btnCompDirSelect.addActionListener(
        new ActionListener() {

          @Override
          public void actionPerformed(ActionEvent arg0) {
            String curr = txtFldCompDir.getText();
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
              txtFldCompDir.setText(newPath);
              loadFCSDir(newPath, false);
              saveProps();
            }
          }
        });
    contentPane.add(btnCompDirSelect, "cell 2 1");

    JLabel lblGatingFile = new JLabel("<html><u>G</u>ating File:</html>");
    contentPane.add(lblGatingFile, "cell 0 2,alignx trailing");

    txtFldGatingFile = new JTextField();
    txtFldGatingFile.setEditable(false);
    contentPane.add(txtFldGatingFile, "cell 1 2,growx");
    txtFldGatingFile.setColumns(10);

    btnGatingFileSelect = new JButton(">");
    btnGatingFileSelect.setMargin(new Insets(0, 3, 0, 3));
    btnGatingFileSelect.setMnemonic(KeyEvent.VK_G);
    btnGatingFileSelect.addActionListener(
        new ActionListener() {

          @Override
          public void actionPerformed(ActionEvent e) {
            String curr = txtFldBaseDir.getText();
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
              saveProps();
            }
          }
        });
    contentPane.add(btnGatingFileSelect, "flowx,cell 2 2");

    scrollPane = new JScrollPane();
    scrollPane.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_ALWAYS);
    scrollPane.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS);
    contentPane.add(scrollPane, "cell 0 3 3 1,grow");

    final HashSet<String> lessThan4CV = new HashSet<>();
    final HashSet<String> lessThan6CV = new HashSet<>();

    String[] lessThan4 = {
      "BB515", "PE", "PE-CF594", "PE-Cy7", "BV 421", "BV 510", "BV 605", "BV 711"
    };
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

    table =
        new JTable() {

          /** */
          private static final long serialVersionUID = 1L;

          @Override
          public Component prepareRenderer(TableCellRenderer renderer, int row, int column) {
            Component c = super.prepareRenderer(renderer, row, column);

            if (column == 0 && boldRows.contains(row)) {
              c.setFont(c.getFont().deriveFont(Font.BOLD));
            }
            if (isCellSelected(row, column)) {
              return c;
            }

            Color col = Color.WHITE;
            if (column > 0 && !statRows.contains(row)) {
              Object val =
                  table
                      .getModel()
                      .getValueAt(
                          table.convertRowIndexToModel(row),
                          table.convertColumnIndexToModel(column));
              if (val != null) {
                if (val instanceof Number) {
                  String colNm =
                      table.getModel().getColumnName(table.convertColumnIndexToModel(column));
                  double value = ((Number) val).doubleValue();
                  if (rdbtnMean.isSelected()) {
                    if (paramMeans.containsKey(colNm)) {
                      double mn = paramMeans.get(colNm);
                      if (value > (mn + .15 * mn) || value < (mn - .15 * mn)) {
                        col = SD2_COLOR;
                      }
                      // if (value > (paramMeans.get(colNm) + paramSDs.get(colNm)) || value <
                      // (paramMeans.get(colNm) - paramSDs.get(colNm))) {
                      // col = SD1_COLOR;
                      // }
                      // if (value > (paramMeans.get(colNm) + 2*paramSDs.get(colNm)) || value <
                      // (paramMeans.get(colNm) - 2*paramSDs.get(colNm))) {
                      // col = SD2_COLOR;
                      // }
                    }
                  } else if (rdbtnSd.isSelected()) {

                  } else if (rdbtnCv.isSelected()) {
                    // <4% :
                    // BLue laser: BB515, PE, PE-CF594 and PE-Cy7
                    // Violet laser: BV421, BV510, BV605 and BV711
                    //
                    // <6%
                    // UV laser: BUV395 and BUV737
                    // Red laser: APC and APC-Cy7
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

    final JFlow fcp = JFlow.createGUI(false);

    table.addMouseListener(
        new MouseAdapter() {

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
                RainbowTestGUI.this.showMeanPanel(col);
              }
              String file = (String) target.getValueAt(row, 0) + ".fcs";
              FCSDataLoader loader = null;
              if (baseFiles.containsKey(file)) {
                loader = baseFiles.get(file);
              } else {
                for (String f : compFiles.keySet()) {
                  if (f.endsWith(file)) { // TODO warning, will break if duplicate files exist
                    loader = compFiles.get(f);
                    break;
                  }
                }
              }
              if (loader == null) {
                return;
              }
              fcp.setData(loader);
              if (gateStrat != null) {
                fcp.setGating(gateStrat);
              }
              SwingUtilities.invokeLater(
                  new Runnable() {

                    @Override
                    public void run() {
                      fcp.setXDataName(col);
                      fcp.setPlotType(PLOT_TYPE.HISTOGRAM);
                      SwingUtilities.getWindowAncestor(fcp).setVisible(true);
                    }
                  });
            }
          }
        });
    table
        .getTableHeader()
        .addMouseListener(
            new MouseAdapter() {

              @Override
              public void mouseClicked(MouseEvent e) {
                if (e.getClickCount() == 2) {
                  int colInd = table.columnAtPoint(e.getPoint());
                  String col = (String) table.getValueAt(0, colInd);
                  RainbowTestGUI.this.showMeanPanel(col);
                }
              }
            });
    table.setShowVerticalLines(false);
    table.setShowHorizontalLines(false);
    scrollPane.setViewportView(table);
    // table.setCellSelectionEnabled(false);
    table.setColumnModel(
        new DefaultTableColumnModel() {

          /** */
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

    btnHideshowColumns = new JButton("Hide/Show Columns");
    btnHideshowColumns.addActionListener(
        new ActionListener() {

          @Override
          public void actionPerformed(ActionEvent arg0) {
            String[] opts = new String[paramNames.size()];
            boolean[] incl = new boolean[paramNames.size()];
            int ind = 0;
            for (String p : paramNames) {
              opts[ind] = p;
              incl[ind] = !hiddenCols.contains(p);
              ind++;
            }
            IncludeExcludeGUI dialog = new IncludeExcludeGUI(RainbowTestGUI.this, opts, incl);
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
    contentPane.add(separator_1, "cell 0 4 2 1,growy");

    AbstractAction gateAction =
        new AbstractAction() {

          /** */
          private static final long serialVersionUID = 1L;

          @Override
          public void actionPerformed(ActionEvent e) {
            reCalcTableData();
          }
        };

    rdbtnUngated = new JRadioButton();
    rdbtnUngated.setAction(gateAction);
    rdbtnUngated.setText("Ungated");
    rdbtnUngated.setSelected(true);
    rdbtnUngated.setMnemonic(KeyEvent.VK_U);
    buttonGroup_1.add(rdbtnUngated);
    contentPane.add(rdbtnUngated, "cell 0 4 2 1");

    rdbtnGated = new JRadioButton();
    rdbtnGated.setAction(gateAction);
    rdbtnGated.setText("Gated");
    rdbtnGated.setEnabled(false);
    rdbtnGated.setMnemonic(KeyEvent.VK_T);
    buttonGroup_1.add(rdbtnGated);
    contentPane.add(rdbtnGated, "cell 0 4 2 1");

    separator = new JSeparator();
    separator.setOrientation(SwingConstants.VERTICAL);
    contentPane.add(separator, "cell 0 4 2 1,growy");

    rdbtnMean = new JRadioButton();
    rdbtnMean.setAction(
        new AbstractAction() {

          /** */
          private static final long serialVersionUID = 1L;

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
    contentPane.add(rdbtnMean, "cell 0 4 2 1");

    rdbtnSd = new JRadioButton();
    rdbtnSd.setAction(
        new AbstractAction() {

          /** */
          private static final long serialVersionUID = 1L;

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
    contentPane.add(rdbtnSd, "cell 0 4 2 1");

    rdbtnCv = new JRadioButton();
    rdbtnCv.setAction(
        new AbstractAction() {

          /** */
          private static final long serialVersionUID = 1L;

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
    contentPane.add(rdbtnCv, "cell 0 4 2 1");

    button = new JButton("X");
    button.addActionListener(
        new ActionListener() {

          @Override
          public void actionPerformed(ActionEvent e) {
            gateStrat = null;
            txtFldGatingFile.setText("");
            reCalcTableData();
          }
        });
    button.setMargin(new Insets(0, 2, 0, 2));
    contentPane.add(button, "cell 2 2");

    separator_2 = new JSeparator();
    separator_2.setOrientation(SwingConstants.VERTICAL);
    contentPane.add(separator_2, "cell 0 4 2 1,growy");

    rdbtnCompensated =
        new JRadioButton(
            new AbstractAction() {

              /** */
              private static final long serialVersionUID = 1L;

              @Override
              public void actionPerformed(ActionEvent arg0) {
                currData = DATA_SET.COMPENSATED;
                reCalcTableData();
              }
            });
    rdbtnCompensated.setText("Compensated");
    buttonGroup_2.add(rdbtnCompensated);
    contentPane.add(rdbtnCompensated, "cell 0 4 2 1");

    rdbtnUncompensated =
        new JRadioButton(
            new AbstractAction() {

              /** */
              private static final long serialVersionUID = 1L;

              @Override
              public void actionPerformed(ActionEvent arg0) {
                currData = DATA_SET.UNCOMPENSATED;
                reCalcTableData();
              }
            });
    rdbtnUncompensated.setText("Uncompensated");
    rdbtnUncompensated.setSelected(true);
    buttonGroup_2.add(rdbtnUncompensated);
    contentPane.add(rdbtnUncompensated, "cell 0 4 2 1");

    btnWarning = new JButton("");
    btnWarning.setIcon(UIManager.getIcon("OptionPane.warningIcon"));
    btnWarning.setMargin(new Insets(0, 2, 0, 2));
    btnWarning.setVisible(false);
    contentPane.add(btnWarning, "cell 2 4");

    meanPanel.setXAxisLabel("File by Date");
    meanPanel.setOpaque(true);
    meanFrame.getContentPane().add(meanPanel, BorderLayout.CENTER);
    meanFrame.getContentPane().add(meanCtrlPanel, BorderLayout.SOUTH);
    ActionListener prevLst =
        new ActionListener() {

          @Override
          public void actionPerformed(ActionEvent arg0) {
            String newCol = arg0.getActionCommand();
            showMeanPanel(newCol);
          }
        };
    meanCtrlPanel.setChangeListener(prevLst);
    meanFrame.setBounds(JFlow.START_X, JFlow.START_Y, JFlow.START_WIDTH, JFlow.START_HEIGHT);
    // meanPanel.setPlotType(OneDPanel.PLOT_TYPE.BOX_PLOT);
    meanPanel.setPlotType(OneDPanel.PLOT_TYPE.DOT_LINE_PLOT);
    meanPanel.setAxisXHeight(AbstractPanel2.HEIGHT_X_AXIS - AbstractPanel2.HEIGHT_X_AXIS / 5);

    meanPanel.setShow1SDLines(false);
    meanPanel.setShow2SDLines(true);
    meanPanel.setShowMean15Line(true);
    meanPanel.setShowRegressionLine(true);

    meanPanel.addPropertyChangeListener(
        "REGRESSION",
        new PropertyChangeListener() {

          @Override
          public void propertyChange(PropertyChangeEvent evt) {
            saveProps();
          }
        });

    meanFrame.setJMenuBar(createMenuBar());
    meanCtrlPanel.link(meanPanel);

    loadCache();
    loadProps();
  }

  private JMenuBar createMenuBar() {
    JMenuBar menuBar = new JMenuBar();

    JMenu actMenu = new JMenu("Actions");
    JMenuItem screenshotItem = new JMenuItem("Screenshot");
    screenshotItem.addActionListener(
        new ActionListener() {

          @Override
          public void actionPerformed(ActionEvent e) {
            JFileChooser fileChooser = new JFileChooser(".");
            int fileOpenActionSelected = fileChooser.showSaveDialog(RainbowTestGUI.this);
            if (fileOpenActionSelected == JFileChooser.APPROVE_OPTION) {
              File fileToOpen = fileChooser.getSelectedFile();
              meanPanel.screenCapture(fileToOpen.toString() + ".png");
            }
          }
        });

    actMenu.add(screenshotItem);
    menuBar.add(actMenu);

    return menuBar;
  }

  private void showMeanPanel(String col) {
    // Get all means (yData), get Files&Dates (xData), get Mean/SD (meanSD)
    double[] yDataBase, yDataComp;
    String[] baseLbls, compLbls;

    TreeMap<Date, ArrayList<Double>> meanMap = new TreeMap<>();
    TreeMap<Date, ArrayList<String>> fileMap = new TreeMap<>();
    int count = 0;
    for (Entry<String, FCSDataLoader> l : baseFiles.entrySet()) {
      Date key;
      if (cache.containsKey(l.getKey())) {
        key = cache.get(l.getKey()).runDate;
      } else {
        key = l.getValue().getRunDate();
      }
      ArrayList<Double> means = meanMap.get(key);
      ArrayList<String> files = fileMap.get(key);
      if (means == null) {
        means = new ArrayList<>();
        meanMap.put(key, means);
      }
      if (files == null) {
        files = new ArrayList<>();
        fileMap.put(key, files);
      }
      means.add(fileParamMeanMap.get(l.getKey()).get(col));
      files.add(l.getKey());
      count++;
    }

    yDataBase = new double[count];
    baseLbls = new String[count];
    int ind = 0;
    for (Entry<Date, ArrayList<Double>> etr : meanMap.entrySet()) {
      for (int i = 0; i < etr.getValue().size(); i++) {
        yDataBase[ind] = etr.getValue().get(i);
        baseLbls[ind] = fileMap.get(etr.getKey()).get(i);
        ind++;
      }
    }

    meanMap = new TreeMap<>();
    fileMap = new TreeMap<>();
    count = 0;
    for (Entry<String, FCSDataLoader> l : compFiles.entrySet()) {
      Date key;
      if (cache.containsKey(l.getKey())) {
        key = cache.get(l.getKey()).runDate;
      } else {
        key = l.getValue().getRunDate();
      }
      ArrayList<Double> means = meanMap.get(key);
      ArrayList<String> files = fileMap.get(key);
      if (means == null) {
        means = new ArrayList<>();
        meanMap.put(key, means);
      }
      if (files == null) {
        files = new ArrayList<>();
        fileMap.put(key, files);
      }
      means.add(fileParamMeanMap.get(l.getKey()).get(col));
      files.add(l.getKey());
      count++;
    }
    yDataComp = new double[count];
    compLbls = new String[count];
    for (Entry<Date, ArrayList<Double>> etr : meanMap.entrySet()) {
      for (int i = 0; i < etr.getValue().size(); i++) {
        yDataComp[ind - yDataBase.length] = etr.getValue().get(i);
        compLbls[ind - yDataBase.length] = fileMap.get(etr.getKey()).get(i);
        ind++;
      }
    }

    ArrayList<String> cols = new ArrayList<>();
    ind = 0;
    for (int i = 1; i < table.getColumnModel().getColumnCount(); i++) {
      String hdr = (String) table.getColumnModel().getColumn(i).getHeaderValue();
      cols.add(hdr);
      if (hdr.equals(col)) {
        ind = i;
      }
    }

    meanCtrlPanel.setColumns(cols, ind - 1);

    meanPanel.setData(
        col, new String[][] {baseLbls, compLbls}, new double[][] {yDataBase, yDataComp});
    meanPanel.setYAxisLabel("Mean - " + col);
    meanFrame.setTitle("Genvisis - FCS Overall Mean/SD - " + col);
    meanFrame.setVisible(true);
    meanPanel.paintAgain();
  }

  private static final String CACHE_FILE = ".URB_DATA_CACHE.ser.gz";

  @SuppressWarnings("unchecked")
  private synchronized void loadCache() {
    if (Files.exists(CACHE_FILE) && (new File(CACHE_FILE).length() > 0)) {
      cache = (HashMap<String, CacheObject>) SerializedFiles.readSerial(CACHE_FILE, false);
    }
    if (cache == null) {
      cache = new HashMap<>();
    }
  }

  private synchronized void saveCache() {
    SerializedFiles.writeSerial(cache, CACHE_FILE, true);
  }

  private static final String PROP_FILE = "rainbow.properties";
  private static final String PROPKEY_COMPAREDIR = "COMPARE_DIR";
  private static final String PROPKEY_BASEDIR = "BASE_DIR";
  private static final String PROPKEY_GATEFILE = "GATING_FILE";
  private static final String PROPKEY_COLS = "HIDDEN_COLUMNS";
  private static final String PROPKEY_REGRESSION_L = "REGRESSION_DROPS";
  private static final String PROPKEY_REGRESSION_G = "REGRESSION_DROPS_GLOBAL";

  private void saveProps() {
    try {
      Properties props = new Properties();
      props.setProperty(PROPKEY_BASEDIR, ext.verifyDirFormat(txtFldBaseDir.getText()));
      props.setProperty(PROPKEY_COMPAREDIR, ext.verifyDirFormat(txtFldCompDir.getText()));
      props.setProperty(PROPKEY_GATEFILE, txtFldGatingFile.getText());
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

      HashMap<String, ArrayList<String>> localDrops = meanPanel.locallyDroppedPoints;
      cols = new StringBuilder();
      int cnt = 0;
      for (Entry<String, ArrayList<String>> drops : localDrops.entrySet()) {
        if (!drops.getValue().isEmpty()) {
          cols.append(drops.getKey()).append("|");
          for (int i = 0; i < drops.getValue().size(); i++) {
            if (i > 0) {
              cols.append(",");
            }
            cols.append(drops.getValue().get(i));
          }
          cnt++;
          if (cnt < localDrops.size()) {
            cols.append(";");
          }
        }
      }
      props.setProperty(PROPKEY_REGRESSION_L, cols.toString());

      cols = new StringBuilder();
      ind = 0;
      for (String s : meanPanel.globallyDroppedPoints) {
        cols.append(s);
        ind++;
        if (ind < meanPanel.globallyDroppedPoints.size()) {
          cols.append(";");
        }
      }
      props.setProperty(PROPKEY_REGRESSION_G, cols.toString());

      File f = new File(PROP_FILE);
      OutputStream out = new FileOutputStream(f);
      props.store(out, "");
    } catch (Exception e) {
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
      String base = props.getProperty(PROPKEY_BASEDIR, "");
      String comp = props.getProperty(PROPKEY_COMPAREDIR, "");
      String gate = props.getProperty(PROPKEY_GATEFILE, "");
      String colsTemp = props.getProperty(PROPKEY_COLS, "");
      String regDrops = props.getProperty(PROPKEY_REGRESSION_L, "");
      String regDropsG = props.getProperty(PROPKEY_REGRESSION_G, "");
      String[] cols = colsTemp.split(";");

      if (!gate.equals("")) {
        txtFldGatingFile.setText(gate);
        setGateFile(gate);
      }
      if (!base.equals("")) {
        txtFldBaseDir.setText(base);
        loadFCSDir(base, true);
      }
      if (!comp.equals("")) {
        txtFldCompDir.setText(comp);
        loadFCSDir(comp, false);
      }
      hiddenCols.clear();
      for (String c : cols) {
        hiddenCols.add(c);
      }
      if (!"".equals(regDrops)) {
        String[] drops = regDrops.split(";");
        for (String drop : drops) {
          String[] pts = drop.split("\\|");
          String lbl = pts[0];
          String[] ptIds = pts[1].split(",");
          ArrayList<String> ids = new ArrayList<>();
          for (String s : ptIds) {
            ids.add(s);
          }
          meanPanel.locallyDroppedPoints.put(lbl, ids);
        }
      }
      if (!"".equals(regDropsG)) {
        String[] drops = regDropsG.split(";");
        for (String s : drops) {
          meanPanel.globallyDroppedPoints.add(s);
        }
      }
      reCalcTableData();
    } catch (Exception e) {
      is = null;
    }
  }

  private void resizeColumnWidth() {
    final TableColumnModel columnModel = table.getColumnModel();
    FontMetrics fm = table.getFontMetrics(table.getFont());
    for (int column = 0; column < table.getColumnCount(); column++) {
      int width =
          fm.stringWidth(columnModel.getColumn(column).getHeaderValue().toString()) + 20; // Min
      // width
      for (int row = 0; row < table.getRowCount(); row++) {
        TableCellRenderer renderer = table.getCellRenderer(row, column);
        Component comp = table.prepareRenderer(renderer, row, column);
        width = Math.max(comp.getPreferredSize().width + 1, width);
      }
      columnModel.getColumn(column).setPreferredWidth(width);
    }
  }

  private void setGateFile(String filePath) {
    try {
      gateStrat = GateFileReader.readGateFile(filePath);
      rdbtnGated.setEnabled(true);
      rdbtnGated.setSelected(true);
      reCalcTableData();
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
      ArrayList<String> files = new ArrayList<>();
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
  }

  private final FilenameFilter fcsFileFilter =
      new FilenameFilter() {

        @Override
        public boolean accept(File dir, String name) {
          return name.endsWith(".fcs");
        }
      };
  private final FileFilter dirFilter =
      new FileFilter() {

        @Override
        public boolean accept(File pathname) {
          return pathname.isDirectory();
        }
      };
  private JSeparator separator_2;
  private JRadioButton rdbtnCompensated;
  private JRadioButton rdbtnUncompensated;
  private DATA_SET currData = DATA_SET.UNCOMPENSATED;
  private final ButtonGroup buttonGroup_2 = new ButtonGroup();
  private JButton btnWarning;

  private void resetShownColumns() {
    TableColumnModel tcm = table.getColumnModel();
    ArrayList<Integer> toRemove = new ArrayList<>();
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

  private void reCalcTableData() {
    if (baseDir == null || dirStruct == null) {
      return;
    }
    TreeSet<String> paramSet = new TreeSet<>();

    for (String f : baseFCSFiles) {
      if (!baseFiles.containsKey(f) && !cache.containsKey(f)) {
        baseFiles.put(f, loadFCSFile(baseDir + f));
      } else if (cache.containsKey(f)) {
        baseFiles.put(f, null);
      }
      List<String> p =
          (cache.containsKey(f)
              ? Arrays.asList(cache.get(f).params)
              : baseFiles.get(f).getAllDisplayableNames(currData));
      paramNames.addAll(p);
      paramSet.addAll(p);
    }

    String[] allFilesFullPaths = dirStruct.getAllFiles();
    for (String f : allFilesFullPaths) {
      if (!compFiles.containsKey(f) && !cache.containsKey(f)) {
        compFiles.put(f, loadFCSFile(f));
      } else if (cache.containsKey(f)) {
        compFiles.put(f, null);
      }
      List<String> p =
          (cache.containsKey(f)
              ? Arrays.asList(cache.get(f).params)
              : compFiles.get(f).getAllDisplayableNames(currData));
      paramNames.addAll(p);
      paramSet.addAll(p);
    }
    String[] paramNames = paramSet.toArray(new String[paramSet.size()]);
    String[] colNames = ArrayUtils.addStrToArray("", paramNames, 0);
    for (String p : paramNames) {
      paramMeanLists.put(p, new ArrayList<Double>());
      paramSDLists.put(p, new ArrayList<Double>());
      paramCVLists.put(p, new ArrayList<Double>());
    }

    dtmMean =
        new DefaultTableModel(colNames, 0) {

          /** */
          private static final long serialVersionUID = 1L;

          @Override
          public boolean isCellEditable(int row, int column) {
            return false;
          }
        };
    dtmSD =
        new DefaultTableModel(colNames, 0) {

          /** */
          private static final long serialVersionUID = 1L;

          @Override
          public boolean isCellEditable(int row, int column) {
            return false;
          }
        };
    dtmCV =
        new DefaultTableModel(colNames, 0) {

          /** */
          private static final long serialVersionUID = 1L;

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
        Double mn = ArrayUtils.mean(Doubles.toArray(paramMeanLists.get(colNm)), true);
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
        Double sd = ArrayUtils.stdev(paramMeanLists.get(colNm).toArray(new Double[0]), true);
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
        Double cv = 100 * (paramSDs.get(colNm) / paramMeans.get(colNm));
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

    TreeMap<Date, ArrayList<String>> fileMap = new TreeMap<>();
    for (Entry<String, FCSDataLoader> l : compFiles.entrySet()) {
      Date key =
          cache.containsKey(l.getKey()) ? cache.get(l.getKey()).runDate : l.getValue().getRunDate();
      ArrayList<String> files = fileMap.get(key);
      if (files == null) {
        files = new ArrayList<>();
        fileMap.put(key, files);
      }
      files.add(l.getKey());
    }
    ArrayList<String> fileList = new ArrayList<>();
    for (Entry<Date, ArrayList<String>> entry : fileMap.entrySet()) {
      fileList.addAll(entry.getValue());
    }

    addFilesToModel(fileList, paramNames, dirStruct.dir);

    checkWarnings(paramNames);

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

  private static final int TREND_ABOVE_1SD_THRESH = 5;
  private static final int TREND_ABOVE_2SD_THRESH = 2;

  private static final double PCT_OF_EVENTS_DEV_TREND = 0.25; // 1-quarter of events outside of 1SD
  // will result in file being reported

  private ArrayList<String> getFileTrendWarnings(String[] params) {
    ArrayList<String> warnings = new ArrayList<>();

    for (String f : compFileList) {
      int cnt = 0;
      HashMap<String, Double> paramVals = fileParamMeanMap.get(f);
      for (String s : params) {
        double param = paramVals.get(s);
        double mean = paramMeans.get(s);
        double sd = paramSDs.get(s);

        if (param < mean - sd || param > mean + sd) {
          cnt++;
        }
      }

      if (cnt >= params.length * PCT_OF_EVENTS_DEV_TREND) {
        warnings.add(
            "Source " + f + " has " + cnt + " events greater than 1SD from parameter means.");
      }
    }

    return warnings;
  }

  private ArrayList<String> check(String[] params) {
    HashMap<String, ArrayList<ArrayList<String>>[]> trendMap = new HashMap<>();
    for (String param : params) {
      ArrayList<ArrayList<String>> all1TrendsForParam = new ArrayList<>();
      ArrayList<ArrayList<String>> all2TrendsForParam = new ArrayList<>();
      ArrayList<ArrayList<String>> all15TrendsForParam = new ArrayList<>();

      if (!paramMeans.containsKey(param) || hiddenCols.contains(param)) {
        continue;
      }

      double mean = paramMeans.get(param);
      double mean15 = (float) (mean * .15);
      double sd = paramSDs.get(param);

      ArrayList<String> trend1 = new ArrayList<>();
      ArrayList<String> trend2 = new ArrayList<>();
      ArrayList<String> any15 = new ArrayList<>();
      for (String f : compFileList) {
        double paramValue = fileParamMeanMap.get(f).get(param);

        if (paramValue < mean - sd || paramValue > mean + sd) {
          trend1.add(f);
          if (trend1.size() >= TREND_ABOVE_1SD_THRESH && !all1TrendsForParam.contains(trend1)) {
            all1TrendsForParam.add(trend1);
          }
          if (paramValue < mean - 2 * sd || paramValue > mean + 2 * sd) {
            trend2.add(f);
            if (trend2.size() >= TREND_ABOVE_2SD_THRESH && !all2TrendsForParam.contains(trend2)) {
              all2TrendsForParam.add(trend2);
            }
          } else {
            trend2 = new ArrayList<>();
          }
        } else {
          trend1 = new ArrayList<>();
          trend2 = new ArrayList<>();
        }
        if (paramValue < mean - mean15 || paramValue > mean + mean15) {
          any15.add(f);
        }
      }
      for (ArrayList<String> s : all1TrendsForParam) {
        for (String ss : s) {
          any15.remove(ss);
        }
      }
      for (ArrayList<String> s : all2TrendsForParam) {
        for (String ss : s) {
          any15.remove(ss);
        }
      }
      all15TrendsForParam.add(any15);
      trendMap.put(
          param, new ArrayList[] {all1TrendsForParam, all2TrendsForParam, all15TrendsForParam});
    }

    ArrayList<String> trendWarnings = new ArrayList<>();

    for (Entry<String, ArrayList<ArrayList<String>>[]> entry : trendMap.entrySet()) {

      if (!entry.getValue()[0].isEmpty()) {
        for (ArrayList<String> trend : entry.getValue()[0]) {
          if (trend.size() == 0) {
            continue;
          }
          String[] p = entry.getKey().split("\\|")[0].trim().split("/");
          String warn =
              trend.size()
                  + "-count trend in parameter "
                  + p[p.length - 1]
                  + " greater than 1SD from the mean.";
          trendWarnings.add(warn);
        }
        for (ArrayList<String> trend : entry.getValue()[1]) {
          if (trend.size() == 0) {
            continue;
          }
          String[] p = entry.getKey().split("\\|")[0].trim().split("/");
          String warn =
              trend.size()
                  + "-count trend in parameter "
                  + p[p.length - 1]
                  + " greater than 2SD from the mean.";
          trendWarnings.add(warn);
        }
        for (ArrayList<String> trend : entry.getValue()[2]) {
          if (trend.size() == 0) {
            continue;
          }
          String[] p = entry.getKey().split("\\|")[0].trim().split("/");
          for (String s : trend) {
            String warn =
                "Source "
                    + s
                    + " has a deviant (>15% of mean) value for parameter "
                    + p[p.length - 1];
            trendWarnings.add(warn);
          }
        }
      }
    }

    return trendWarnings;
  }

  // private ArrayList<String> getParameterTrendWarnings(String[] params) {
  // HashMap<String, Integer> paramTrendsAbove1 = new HashMap<String, Integer>();
  // HashMap<String, Integer> paramTrendsAbove2 = new HashMap<String, Integer>();
  //
  // for (String colNm : params) {
  // if (!paramMeans.containsKey(colNm)) continue;
  //
  // int trendAbove1 = 0;
  // int trendAbove2 = 0;
  // float mean = paramMeans.get(colNm);
  // float sd = paramSDs.get(colNm);
  //
  // for (String f : compFileList) {
  // float param = fileParamMeanMap.get(f).get(colNm);
  // if (param < mean - sd || param > mean + sd) {
  // trendAbove1++;
  // if (param < mean - 2 * sd || param > mean + 2 * sd) {
  // trendAbove2++;
  // if (trendAbove2 >= TREND_ABOVE_2SD_THRESH) {
  // paramTrendsAbove2.put(colNm, trendAbove2);
  // }
  // } else {
  // trendAbove2 = 0;
  // if (trendAbove1 >= TREND_ABOVE_1SD_THRESH) {
  // paramTrendsAbove1.put(colNm, trendAbove1);
  // }
  // }
  // } else {
  // trendAbove1 = 0;
  // trendAbove2 = 0;
  // }
  // }
  // }
  //
  // final ArrayList<String> warnings = new ArrayList<String>();
  // if (!paramTrendsAbove1.isEmpty()) {
  // for (String s : paramTrendsAbove1.keySet()) {
  // String[] p = s.split("\\|")[0].trim().split("/");
  // String warn = paramTrendsAbove1.get(s) + "-count trend in parameter " + p[p.length - 1] +
  // " greater than 1SD from the mean.";
  // warnings.add(warn);
  // }
  // }
  // if (!paramTrendsAbove2.isEmpty()) {
  // for (String s : paramTrendsAbove2.keySet()) {
  // String[] p = s.split("\\|")[0].trim().split("/");
  // String warn = paramTrendsAbove2.get(s) + "-count trend in parameter " + p[p.length - 1] +
  // " greater than 2SD from the mean.";
  // warnings.add(warn);
  // }
  // }
  //
  // return warnings;
  // }

  private void checkWarnings(String[] params) {
    final ArrayList<String> warnings = new ArrayList<>();

    warnings.addAll(check(params));
    warnings.addAll(getFileTrendWarnings(params));

    if (warnings.size() > 0) {
      btnWarning.setVisible(true);
      btnWarning.setAction(
          new AbstractAction() {

            /** */
            private static final long serialVersionUID = 1L;

            @Override
            public void actionPerformed(ActionEvent e) {
              final JPanel msgPane = new JPanel(new MigLayout("", "", ""));
              for (int i = 0; i < warnings.size(); i++) {
                msgPane.add(new JLabel("" + warnings.get(i)), "cell 0 " + i);
              }
              JScrollPane scroll = new JScrollPane();
              scroll.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
              scroll.setViewportView(msgPane);
              JPanel chk = new JPanel();
              chk.add(scroll);
              msgPane.addHierarchyListener(
                  new HierarchyListener() {

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
              JOptionPane.showMessageDialog(
                  RainbowTestGUI.this, scroll, "Trend Warnings", JOptionPane.WARNING_MESSAGE);
            }
          });
      btnWarning.setIcon(UIManager.getIcon("OptionPane.warningIcon"));
    }
  }

  private static final class CacheObject implements Serializable {

    /** */
    private static final long serialVersionUID = 1L;

    String[] params;
    Date runDate;
    String gateFileApplied;
    HashMap<String, Double> paramMeans = new HashMap<>();
    HashMap<String, Double> paramSDs = new HashMap<>();
    HashMap<String, Double> paramCVs = new HashMap<>();
  }

  private HashMap<String, CacheObject> cache = new HashMap<>();

  private void addBaseToModel(String[] paramNames) {
    TreeMap<Date, ArrayList<String>> map = new TreeMap<>();

    for (String f : baseFCSFiles) {
      Date runDate;
      if (cache.containsKey(f)) {
        runDate = cache.get(f).runDate;
      } else {
        runDate = baseFiles.get(f).getRunDate();
      }
      ArrayList<String> fls = map.get(runDate);
      if (fls == null) {
        fls = new ArrayList<>();
        map.put(runDate, fls);
      }
      fls.add(f);
    }
    for (ArrayList<String> fls : map.values()) {
      for (String f : fls) {

        fileParamMeanMap.put(f, new HashMap<String, Double>());
        Object[] rowDataM = new Object[paramNames.length + 1];
        Object[] rowDataS = new Object[paramNames.length + 1];
        Object[] rowDataC = new Object[paramNames.length + 1];
        rowDataM[0] = ext.rootOf(f);
        rowDataS[0] = ext.rootOf(f);
        rowDataC[0] = ext.rootOf(f);

        CacheObject co = cache.get(f);
        if (co != null
            && ((gateStrat == null && "".equals(co.gateFileApplied))
                || gateStrat.getFile().equals(co.gateFileApplied))
            && ArrayUtils.equals(paramNames, co.params, false)) {
          System.out.println("Loading data for " + f + " from cache");
          for (int i = 0; i < paramNames.length; i++) {
            Double mn = co.paramMeans.get(paramNames[i]);
            Double sd = co.paramSDs.get(paramNames[i]);
            Double cv = co.paramCVs.get(paramNames[i]);
            paramMeanLists.get(paramNames[i]).add(mn);
            paramSDLists.get(paramNames[i]).add(sd);
            paramCVLists.get(paramNames[i]).add(cv);
            fileParamMeanMap.get(f).put(paramNames[i], mn);

            rowDataM[i + 1] = mn;
            rowDataS[i + 1] = sd;
            rowDataC[i + 1] = cv;
          }
        } else {
          FCSDataLoader loader = baseFiles.get(f);

          CacheObject newCO = new CacheObject();
          newCO.runDate = loader.getRunDate();
          newCO.gateFileApplied = gateStrat != null ? gateStrat.getFile() : "";
          newCO.params = paramNames;

          for (int i = 0; i < paramNames.length; i++) {

            boolean[] gating = null;
            if (gateStrat != null && rdbtnGated.isSelected()) {
              ArrayList<Gate> gates = gateStrat.getGatesForParamOnly(paramNames[i]);
              System.out.println(
                  "Applying "
                      + gates.size()
                      + " gates (not including parent-gates) to parameter "
                      + paramNames[i]);
              for (Gate g : gates) {
                if (gating == null) {
                  gating = g.gate(loader);
                } else {
                  ArrayUtils.booleanArrayAndInPlace(gating, g.gate(loader));
                }
              }
            }

            double[] data = loader.getData(paramNames[i], true);
            if (gating != null) {
              data = ArrayUtils.subArray(data, gating);
            }
            Double mn = ArrayUtils.mean(data, true);
            Double sd = ArrayUtils.stdev(data, true);
            Double cv = 100 * (sd / mn);
            paramMeanLists.get(paramNames[i]).add(mn);
            paramSDLists.get(paramNames[i]).add(sd);
            paramCVLists.get(paramNames[i]).add(cv);
            fileParamMeanMap.get(f).put(paramNames[i], mn);

            newCO.paramMeans.put(paramNames[i], mn);
            newCO.paramSDs.put(paramNames[i], sd);
            newCO.paramCVs.put(paramNames[i], cv);

            rowDataM[i + 1] = mn;
            rowDataS[i + 1] = sd;
            rowDataC[i + 1] = cv;
          }

          cache.put(f, newCO);
        }

        dtmMean.addRow(rowDataM);
        dtmSD.addRow(rowDataS);
        dtmCV.addRow(rowDataC);
      }
    }
  }

  private ArrayList<String> compFileList = new ArrayList<>();

  private void addFilesToModel(List<String> files, String[] paramNames, String removePrep) {
    // int colCnt = dtmMean.getColumnCount();
    // int rows = dtmMean.getRowCount();

    // Object[] newRow = new Object[colCnt];
    // newRow[0] = df.dir.replaceFirst(removePrep, "");
    // dtmMean.addRow(newRow);
    // dtmSD.addRow(newRow);
    // dtmCV.addRow(newRow);
    // boldRows.add(rows);
    // rows++;

    compFileList = new ArrayList<>();
    for (String f : files) {
      compFileList.add(f);
      fileParamMeanMap.put(f, new HashMap<String, Double>());
      Object[] rowDataM = new Object[paramNames.length + 1];
      Object[] rowDataS = new Object[paramNames.length + 1];
      Object[] rowDataC = new Object[paramNames.length + 1];
      rowDataM[0] = ext.rootOf(f);
      rowDataS[0] = ext.rootOf(f);
      rowDataC[0] = ext.rootOf(f);

      // load from cache
      CacheObject co = cache.get(f);
      if (co != null
          && ((gateStrat == null && "".equals(co.gateFileApplied))
              || gateStrat.getFile().equals(co.gateFileApplied))
          && ArrayUtils.equals(paramNames, co.params, false)) {
        System.out.println("Loading data for " + f + " from cache");
        for (int i = 0; i < paramNames.length; i++) {
          Double mn = co.paramMeans.get(paramNames[i]);
          Double sd = co.paramSDs.get(paramNames[i]);
          Double cv = co.paramCVs.get(paramNames[i]);
          paramMeanLists.get(paramNames[i]).add(mn);
          paramSDLists.get(paramNames[i]).add(sd);
          paramCVLists.get(paramNames[i]).add(cv);
          fileParamMeanMap.get(f).put(paramNames[i], mn);

          rowDataM[i + 1] = mn;
          rowDataS[i + 1] = sd;
          rowDataC[i + 1] = cv;
        }
      } else {
        FCSDataLoader loader = compFiles.get(/* df.dir + */ f);

        CacheObject newCO = new CacheObject();
        newCO.runDate = loader.getRunDate();
        newCO.gateFileApplied = gateStrat != null ? gateStrat.getFile() : "";
        newCO.params = paramNames;

        for (int i = 0; i < paramNames.length; i++) {

          boolean[] gating = null;
          if (gateStrat != null && rdbtnGated.isSelected()) {
            ArrayList<Gate> gates = gateStrat.getGatesForParamOnly(paramNames[i]);
            for (Gate g : gates) {
              if (gating == null) {
                gating = g.gate(loader);
              } else {
              }
              ArrayUtils.booleanArrayAndInPlace(gating, g.gate(loader));
            }
          }

          double[] data = loader.getData(paramNames[i], true);
          if (gating != null) {
            data = ArrayUtils.subArray(data, gating);
          }
          Double mn = ArrayUtils.mean(data, true);
          fileParamMeanMap.get(f).put(paramNames[i], mn);
          Double sd = ArrayUtils.stdev(data, true);
          Double cv = 100 * (sd / mn);
          rowDataM[i + 1] = mn;
          rowDataS[i + 1] = sd;
          rowDataC[i + 1] = cv;

          newCO.paramMeans.put(paramNames[i], mn);
          newCO.paramSDs.put(paramNames[i], sd);
          newCO.paramCVs.put(paramNames[i], cv);

          cache.put(f, newCO);
        }
      }

      dtmMean.addRow(rowDataM);
      dtmSD.addRow(rowDataS);
      dtmCV.addRow(rowDataC);
    }
    // if (df.files.length > 0 && df.getAllFiles().length > 0) {
    // dtmMean.addRow(new Object[colCnt]);
    // dtmSD.addRow(new Object[colCnt]);
    // dtmCV.addRow(new Object[colCnt]);
    // }
    //
    // for (DirFile sub : df.subDirs) {
    // addFilesToModel(sub, paramNames, removePrep);
    // }

  }

  private void loadBase(String dir) {
    File dirFile = new File(dir);
    baseFCSFiles = dirFile.list(fcsFileFilter);
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

  /** Launch the application. */
  public static void main(String[] args) {
    EventQueue.invokeLater(
        new Runnable() {

          @Override
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
