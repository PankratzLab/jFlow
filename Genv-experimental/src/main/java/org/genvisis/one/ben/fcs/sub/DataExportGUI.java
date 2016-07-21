package org.genvisis.one.ben.fcs.sub;

import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.HashMap;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTree;
import javax.swing.ScrollPaneConstants;
import javax.swing.SwingUtilities;
import javax.swing.border.BevelBorder;
import javax.swing.border.EmptyBorder;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeModel;

import net.miginfocom.swing.MigLayout;

import org.genvisis.cnv.gui.JAccordionPanel;
import org.genvisis.one.ben.fcs.FCSPlot;
import org.genvisis.one.ben.fcs.gating.Gate;
import org.genvisis.one.ben.fcs.gating.GatingStrategy;

public class DataExportGUI extends JDialog {

    private JPanel contentPane;
    private JTree tree;
    private HashMap<DefaultMutableTreeNode, Gate> gateMap = new HashMap<DefaultMutableTreeNode, Gate>();
    private GatingStrategy gating;

    /**
     * Launch the application.
     */
//    public static void main(String[] args) {
//        EventQueue.invokeLater(new Runnable() {
//            public void run() {
//                try {
//                    DataExportGUI frame = new DataExportGUI(null);
//                    frame.setVisible(true);
//                } catch (Exception e) {
//                    e.printStackTrace();
//                }
//            }
//        });
//    }

    FCSPlot plot;
    ArrayList<JCheckBox> boxes;
    
    public DataExportGUI(FCSPlot fcsPlot) {
        super(SwingUtilities.getWindowAncestor(fcsPlot), "Select a set of files and gates to export", ModalityType.APPLICATION_MODAL);
        this.plot = fcsPlot;
        setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
        setBounds(100, 100, 450, 300);
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
        
        JPanel pnl = fileAP.contentPanel;
        pnl.setLayout(new MigLayout("", "[grow]", ""));
        
        if (plot != null) {
            ArrayList<String> files = plot.getAddedFiles();
            boxes = new ArrayList<JCheckBox>();
            for (int i = 0; i < files.size(); i++) {
                JCheckBox box = new JCheckBox(files.get(i));
                pnl.add(box, "growx, cell 0 " + i);
                boxes.add(box);
            }
        }
        
        pnl = gateAP.contentPanel;
        pnl.setLayout(new MigLayout("", "[grow]", ""));
        
        tree = new JTree(new DefaultTreeModel(null));
        tree.setBorder(new BevelBorder(BevelBorder.LOWERED, null, null, null, null));
        tree.setRootVisible(false);
        tree.setShowsRootHandles(true);
        
        pnl.add(tree, "cell 0 0, grow");
        
        if (plot != null) {
            this.gating = plot.getGatingStrategy();
            resetGating();
        }
        
        
        JPanel panel = new JPanel();
        contentPane.add(panel, "south,grow");
        panel.setLayout(new MigLayout("ins 5 0 0 0", "[grow][][]", "[]"));
        
        JButton btnExport = new JButton("Export");
        panel.add(btnExport, "cell 1 0");
        
        JButton btnCancel = new JButton("Cancel");
        btnCancel.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent arg0) {
                setVisible(false);
                dispose();
            }
        });
        panel.add(btnCancel, "cell 2 0");
    }

    private void resetGating() {
        if (this.gating == null || this.tree == null) return;
        ArrayList<Gate> roots = gating.getRootGates();
        
        DefaultMutableTreeNode rootNode = new DefaultMutableTreeNode();
        for (Gate g : roots) {
            addGatesToTree(rootNode, g);
        }
        
        tree.setModel(new DefaultTreeModel(rootNode));        
        expandAllNodes();
    }
    
    private void expandAllNodes() {
        int j = tree.getRowCount();
        int i = 0;
        while(i < j) {
            tree.expandRow(i);
            i += 1;
            j = tree.getRowCount();
        }
    }
    
    private void addGatesToTree(DefaultMutableTreeNode root, Gate g) {
        DefaultMutableTreeNode child = new DefaultMutableTreeNode(g.getName());
        gateMap.put(child, g);
        root.add(child);
        for (Gate childGate : g.getChildGates()) {
            addGatesToTree(child, childGate);
        }
    }
}
