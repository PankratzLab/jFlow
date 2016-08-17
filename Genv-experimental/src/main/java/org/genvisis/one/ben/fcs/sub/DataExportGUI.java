package org.genvisis.one.ben.fcs.sub;

import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.HashMap;

import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.JTree;
import javax.swing.ScrollPaneConstants;
import javax.swing.SwingUtilities;
import javax.swing.border.BevelBorder;
import javax.swing.border.EmptyBorder;
import javax.swing.event.CaretEvent;
import javax.swing.event.CaretListener;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.TreePath;

import org.genvisis.cnv.gui.JAccordionPanel;
import org.genvisis.common.Files;
import org.genvisis.one.ben.fcs.FCSPlot;
import org.genvisis.one.ben.fcs.gating.Gate;
import org.genvisis.one.ben.fcs.gating.GatingStrategy;

import net.miginfocom.swing.MigLayout;

public class DataExportGUI extends JDialog {

  /**
  * 
  */
  private static final long serialVersionUID = 1L;
  private final JPanel contentPane;
  private final JTree tree;
  private final HashMap<DefaultMutableTreeNode, Gate> gateMap =
      new HashMap<DefaultMutableTreeNode, Gate>();
  private GatingStrategy gating;
  private final FCSPlot plot;
  private ArrayList<JCheckBox> boxes;
  private volatile boolean cancelled = true;
  private final JTextField outputField;
  private final JRadioButton btnCounts;
  private final JRadioButton btnPcts;
  private final JButton btnExport;

  public DataExportGUI(FCSPlot fcsPlot) {
    super(SwingUtilities.getWindowAncestor(fcsPlot), "Select a set of files and gates to export",
        ModalityType.APPLICATION_MODAL);
    plot = fcsPlot;
    setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
    setBounds(100, 100, 650, 600);
    contentPane = new JPanel();
    contentPane.setBorder(null);
    contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
    setContentPane(contentPane);
    contentPane.setLayout(new MigLayout("ins 0", "[grow]", "[grow]"));

    JScrollPane scrollPane = new JScrollPane();
    scrollPane.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS);
    contentPane.add(scrollPane, "cell 0 0,grow");

    JPanel panel_1 = new JPanel();
    scrollPane.setViewportView(panel_1);
    panel_1.setLayout(new MigLayout("", "[grow]", "[]"));
    JAccordionPanel fileAP = new JAccordionPanel();
    JAccordionPanel gateAP = new JAccordionPanel();
    JAccordionPanel optAP = new JAccordionPanel();

    ArrayList<JAccordionPanel> bg = new ArrayList<JAccordionPanel>();

    JLabel fileLabel = new JLabel("<html><u>Select Files...</u></html>");
    fileLabel.setFont(fileLabel.getFont().deriveFont(Font.PLAIN, 14));
    fileAP.topPanel.add(fileLabel, "pad 0 10 0 0, cell 0 0, grow");
    fileAP.topPanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED));
    fileAP.contentPanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.LOWERED));
    fileAP.addToGroup(bg);

    JLabel gateLabel = new JLabel("<html><u>Select Gating...</u></html>");
    gateLabel.setFont(gateLabel.getFont().deriveFont(Font.PLAIN, 14));
    gateAP.topPanel.add(gateLabel, "pad 0 10 0 0, cell 0 0, grow");
    gateAP.topPanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED));
    gateAP.contentPanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.LOWERED));
    gateAP.addToGroup(bg);

    JLabel optLabel = new JLabel("<html><u>Export Options...</u></html>");
    optLabel.setFont(optLabel.getFont().deriveFont(Font.PLAIN, 14));
    optAP.topPanel.add(optLabel, "pad 0 10 0 0, cell 0 0, grow");
    optAP.topPanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.RAISED));
    optAP.contentPanel.setBorder(BorderFactory.createBevelBorder(BevelBorder.LOWERED));
    optAP.addToGroup(bg);

    gateAP.shrink();
    optAP.shrink();

    panel_1.add(fileAP, "cell 0 0, growx");
    panel_1.add(gateAP, "cell 0 1, growx");
    panel_1.add(optAP, "cell 0 2, growx");

    Font lblFont = new Font("Arial", 0, 9);

    JPanel pnl = fileAP.contentPanel;
    pnl.setLayout(new MigLayout("", "[grow]", ""));
    JLabel lbl = new JLabel("(Currently-loaded file marked in bold)");
    lbl.setFont(lblFont);
    pnl.add(lbl, "cell 0 0");

    if (plot != null) {
      ArrayList<String> files = plot.getAddedFiles();
      String curr = plot.getCurrentFile();
      boxes = new ArrayList<JCheckBox>();
      for (int i = 0; i < files.size(); i++) {
        JCheckBox box = new JCheckBox(files.get(i));
        box.addActionListener(new ActionListener() {
          @Override
          public void actionPerformed(ActionEvent e) {
            updateValidity();
          }
        });
        if (curr == null || !files.get(i).equals(curr)) {
          box.setFont(box.getFont().deriveFont(Font.PLAIN));
        }
        pnl.add(box, "growx, cell 0 " + (i + 1));
        boxes.add(box);
      }
    }

    pnl = gateAP.contentPanel;
    pnl.setLayout(new MigLayout("", "[grow]", ""));

    tree = new JTree(new DefaultTreeModel(null));
    tree.addTreeSelectionListener(new TreeSelectionListener() {
      @Override
      public void valueChanged(TreeSelectionEvent arg0) {
        updateValidity();
      }
    });
    tree.setBorder(new BevelBorder(BevelBorder.LOWERED, null, null, null, null));
    tree.setRootVisible(false);
    tree.setShowsRootHandles(true);

    lbl =
        new JLabel("(Hold Ctrl to select multiple gates | Hold Shift to select a range of gates)");
    lbl.setFont(lblFont);
    pnl.add(lbl, "cell 0 0");
    pnl.add(tree, "cell 0 1, grow");

    if (plot != null) {
      gating = plot.getGatingStrategy();
      resetGating();
    }

    pnl = optAP.contentPanel;
    pnl.setLayout(new MigLayout("", "[][][grow][]", "[][][][]"));

    JLabel lblOutputType = new JLabel("<html><u>Output Type:</u></html>");
    pnl.add(lblOutputType, "cell 0 0 2 1");

    ButtonGroup buttonGroup = new ButtonGroup();
    btnPcts = new JRadioButton("Percents");
    btnPcts.setSelected(true);
    buttonGroup.add(btnPcts);
    pnl.add(btnPcts, "cell 0 1");

    btnCounts = new JRadioButton("Counts");
    buttonGroup.add(btnCounts);
    pnl.add(btnCounts, "cell 1 1");

    JLabel lblOutputFilename = new JLabel("Output Filename:");
    pnl.add(lblOutputFilename, "cell 0 2");

    outputField = new JTextField();
    outputField.addCaretListener(new CaretListener() {
      @Override
      public void caretUpdate(CaretEvent arg0) {
        updateValidity();
      }
    });
    pnl.add(outputField, "cell 0 3 3 1,growx");

    JButton btnFileSelect = new JButton(">");
    btnFileSelect.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent arg0) {
        JFileChooser jfc = new JFileChooser();
        jfc.setMultiSelectionEnabled(false);
        int code = jfc.showSaveDialog(DataExportGUI.this);
        if (code == JFileChooser.APPROVE_OPTION) {
          outputField.setText(jfc.getSelectedFile().getAbsolutePath());
        }
        updateValidity();
      }
    });
    pnl.add(btnFileSelect, "cell 3 3");

    JPanel panel = new JPanel();
    contentPane.add(panel, "south,grow");
    panel.setLayout(new MigLayout("ins 5 0 0 0", "[grow][][]", "[]"));

    btnExport = new JButton("Export");
    btnExport.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        cancelled = false;
        setVisible(false);
      }
    });
    btnExport.setEnabled(false);
    panel.add(btnExport, "cell 1 0");

    JButton btnCancel = new JButton("Cancel");
    btnCancel.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent arg0) {
        setVisible(false);
        dispose();
      }
    });
    panel.add(btnCancel, "cell 2 0");
  }

  private void addGatesToTree(DefaultMutableTreeNode root, Gate g) {
    DefaultMutableTreeNode child = new DefaultMutableTreeNode(
        g.getName() == null || "".equals(g.getName()) ? g.getID() : g.getName(), true);
    gateMap.put(child, g);
    root.add(child);
    for (Gate childGate : g.getChildGates()) {
      addGatesToTree(child, childGate);
    }
  }

  private void expandAllNodes() {
    int j = tree.getRowCount();
    int i = 0;
    while (i < j) {
      tree.expandRow(i);
      i += 1;
      j = tree.getRowCount();
    }
  }

  public boolean getCountsExportSelected() {
    return btnCounts.isSelected();
  }

  public String getOutputFile() {
    return outputField.getText();
  }

  public ArrayList<String> getSelectedFiles() {
    ArrayList<String> files = new ArrayList<String>();

    for (JCheckBox box : boxes) {
      if (box.isSelected()) {
        files.add(box.getText());
      }
    }

    return files;
  }

  public ArrayList<Gate> getSelectedGates() {
    ArrayList<Gate> gates = new ArrayList<Gate>();

    TreePath[] selected = tree.getSelectionPaths();
    for (TreePath path : selected) {
      DefaultMutableTreeNode node = (DefaultMutableTreeNode) path.getLastPathComponent();
      gates.add(gateMap.get(node));
    }

    return gates;
  }

  private void resetGating() {
    if (gating == null || tree == null) {
      return;
    }
    ArrayList<Gate> roots = gating.getRootGates();

    DefaultMutableTreeNode rootNode = new DefaultMutableTreeNode();
    for (Gate g : roots) {
      addGatesToTree(rootNode, g);
    }

    tree.setModel(new DefaultTreeModel(rootNode, true));
    expandAllNodes();
  }

  private void updateValidity() {
    boolean anySel = false;
    for (JCheckBox chk : boxes) {
      anySel = chk.isSelected();
      if (anySel) {
        break;
      }
    }
    boolean anyGate = tree.getSelectionCount() > 0;
    boolean file = false;
    String fl = outputField.getText();
    if (!"".equals(fl) && !Files.exists(fl)) {
      file = true;
    }
    btnExport.setEnabled(anySel && anyGate && file);
  }

  public boolean wasCancelled() {
    return cancelled;
  }
}
