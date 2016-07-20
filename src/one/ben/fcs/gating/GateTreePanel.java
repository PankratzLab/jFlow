package one.ben.fcs.gating;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Font;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JScrollPane;
import javax.swing.JTree;
import javax.swing.ScrollPaneConstants;
import javax.swing.SwingConstants;
import javax.swing.WindowConstants;
import javax.swing.border.BevelBorder;
import javax.swing.border.Border;
import javax.swing.border.EtchedBorder;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.TreePath;
import javax.xml.parsers.ParserConfigurationException;

import net.miginfocom.swing.MigLayout;
import one.ben.fcs.FCSPlot;

import org.xml.sax.SAXException;

public class GateTreePanel extends JPanel {

    GatingStrategy gating;
    private JTree tree;
    private JPanel topPanel;
    private JScrollPane scrollPane;
    private FCSPlot plot;
    private HashMap<DefaultMutableTreeNode, Gate> gateMap = new HashMap<DefaultMutableTreeNode, Gate>();
    private volatile boolean showing = false;
    /**
     * Create the panel.
     */
    public GateTreePanel(FCSPlot plot) {
        setLayout(new BorderLayout(0, 0));
        
        this.plot = plot;
        
        topPanel = new JPanel();
        topPanel.setBorder(new EtchedBorder(EtchedBorder.LOWERED, null, null));
        topPanel.setBackground(Color.WHITE);
        add(topPanel, BorderLayout.NORTH);
        topPanel.setLayout(new MigLayout("ins 0, hidemode 3", "[grow]", "[grow,center]"));

        treePopup = new JPopupMenu();
        treePopup.setLightWeightPopupEnabled(true);
        treePopup.setOpaque(true);
        final JButton button = new JButton("\\/");
        button.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent arg0) {
                if (showing) {
                    treePopup.setVisible(false);
                    showing = false;
                } else {
                    treePopup.setPopupSize(topPanel.getWidth(), topPanel.getHeight() * 8);
                    Point p = topPanel.getLocationOnScreen();
                    treePopup.show(button, 0, 0);
                    treePopup.setLocation((int) p.getX(), (int) p.getY() + topPanel.getHeight());
                    showing = true;
                }
            }
        });
        topPanel.add(button, "east");
        
        scrollPane_1 = new JScrollPane();
        scrollPane_1.setBorder(null);
        scrollPane_1.setViewportBorder(null);
        scrollPane_1.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_NEVER);
        scrollPane_1.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
        topPanel.add(scrollPane_1, "cell 0 0,grow");
        
        breadcrumbPanel = new JPanel();
        scrollPane_1.setViewportView(breadcrumbPanel);
        breadcrumbPanel.setBorder(null);
        breadcrumbPanel.setBackground(Color.WHITE);
        breadcrumbPanel.setLayout(new MigLayout("ins 0", "[]", "[grow]"));
        addComponentListener(new ComponentAdapter() {
            @Override
            public void componentResized(ComponentEvent e) {
                super.componentResized(e);
                Rectangle rect = new Rectangle(breadcrumbPanel.getWidth() - 10, 10, breadcrumbPanel.getWidth(), 10);
                ((JComponent) breadcrumbPanel.getParent()).scrollRectToVisible(rect);
                breadcrumbPanel.repaint();
                
                if (treePopup.isVisible()) {
                    treePopup.setPopupSize(topPanel.getWidth(), topPanel.getHeight() * 8);
                }
            }
        });
        
        scrollPane = new JScrollPane();
        
        treePopup.add(scrollPane, BorderLayout.CENTER);
        tree = new JTree(new DefaultTreeModel(null));
        scrollPane.setViewportView(tree);
        tree.setBorder(new BevelBorder(BevelBorder.LOWERED, null, null, null, null));
        tree.setRootVisible(false);
        tree.setShowsRootHandles(true);
        tree.addTreeSelectionListener(treeListener);
    }
    
    public void setGating(GatingStrategy gating) {
        this.gating = gating;
        resetGating();
    }
    
    private void resetGating() {
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
    
    Font breadcrumbFont = new Font("Arial", 0, 10);
    
    TreeSelectionListener treeListener = new TreeSelectionListener() {
        @Override
        public void valueChanged(TreeSelectionEvent e) {
            if (treePopup.isVisible()) {
                treePopup.setVisible(false);
            }
            if (!e.isAddedPath()) {
                return;
            }
            
            TreePath newPath = e.getNewLeadSelectionPath();
            DefaultTreeModel tm = (DefaultTreeModel) tree.getModel();
            
            Object[] path = newPath.getPath();
            
            breadcrumbPanel.removeAll();
            int index = 0;
            for (int i = 0; i < path.length; i++) {
                DefaultMutableTreeNode obj = (DefaultMutableTreeNode) path[i];
                String nm = (String) obj.getUserObject();
                String chNm = (String) (i == path.length - 1 ? "" : ((DefaultMutableTreeNode) path[i+1]).getUserObject());
                final JLabel lbl = new JLabel(nm);
                lbl.setBorder(breadcrumbBorder);
                lbl.setFont(breadcrumbFont);
                lbl.setVerticalAlignment(SwingConstants.CENTER);
                breadcrumbPanel.add(lbl, "cell " + index + " 0");
                index++;
                
                ArrayList<DefaultMutableTreeNode> sibs = new ArrayList<DefaultMutableTreeNode>();
                int cnt = tm.getChildCount(obj);
                for (int c = 0; c < cnt; c++) {
                    sibs.add(((DefaultMutableTreeNode) tm.getChild(obj, c)));
                }
                
                breadcrumbPanel.add(getBreadcrumbLabel(">", chNm, sibs, lbl), "cell " + index + " 0");
                index++;
                
            }
            breadcrumbPanel.revalidate();
            breadcrumbPanel.repaint();
            Rectangle rect = new Rectangle(breadcrumbPanel.getWidth() - 10, 10, breadcrumbPanel.getWidth(), 10);
            ((JComponent) breadcrumbPanel.getParent()).scrollRectToVisible(rect);
            plot.gateSelected(gateMap.get(path[path.length - 1]));
        }
    };
    private JPanel breadcrumbPanel;
    private Border mouseoverBorder = BorderFactory.createLineBorder(Color.black, 1);
    private Border breadcrumbBorder = BorderFactory.createEmptyBorder(0, 1, 0, 1);
    
    private JLabel getBreadcrumbLabel(String txt, String me, final ArrayList<DefaultMutableTreeNode> sibs, JLabel prevLbl) {
        final JLabel lbl = new JLabel(txt);
        lbl.setBorder(breadcrumbBorder);
        lbl.setFont(breadcrumbFont);
        lbl.setVerticalAlignment(SwingConstants.CENTER);
        
        JPopupMenu popup = new JPopupMenu();
        StringBuilder sb = new StringBuilder("<html>");
        sb.append("<b>").append(me).append("</b><br />");
        for (int c = 0; c < sibs.size(); c++) {
            final int ind = c; 
            sb.append("--").append((String) sibs.get(c).getUserObject());
            Action act = new AbstractAction((String) sibs.get(c).getUserObject()) {
                @Override
                public void actionPerformed(ActionEvent e) {
                    tree.setSelectionPath(new TreePath(((DefaultTreeModel)tree.getModel()).getPathToRoot(sibs.get(ind))));
                }
            };
            JMenuItem jmi = new JMenuItem(act);
            popup.add(jmi);
            if (c < sibs.size() - 1) sb.append("<br />");
        }
        sb.append("</html>");
        
        lbl.setToolTipText(sb.toString());
        
        lbl.setComponentPopupMenu(popup);
        prevLbl.setComponentPopupMenu(popup);
        lbl.addMouseListener(mouseOver);
        prevLbl.addMouseListener(mouseOver);
        return lbl;
    }
    
    MouseAdapter mouseOver = new MouseAdapter() {
        @Override
        public void mouseEntered(MouseEvent e) {
            super.mouseEntered(e);
            ((JComponent) e.getSource()).setBorder(mouseoverBorder);
        }
        @Override
        public void mouseExited(MouseEvent e) {
            super.mouseExited(e);
            ((JComponent) e.getSource()).setBorder(breadcrumbBorder);
        }
        public void mouseClicked(MouseEvent e) {
            super.mouseClicked(e);
            ((JComponent) e.getSource()).getComponentPopupMenu().show(((JComponent) e.getSource()), e.getX(), e.getY());
        };
    };
    private JScrollPane scrollPane_1;
    private JPopupMenu treePopup;
    
    public static void main(String[] args) throws ParserConfigurationException, SAXException, IOException {
        JFrame frame = new JFrame();
        frame.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        
        GateTreePanel gtp = new GateTreePanel(null);
        frame.getContentPane().add(gtp);
        frame.setVisible(true);
        
        GatingStrategy gs = GateFileReader.readGateFile("F:/Flow/PBMC A&C comparison HB ZF DHS 20-Mar-2016.wspt");
        gtp.setGating(gs);
        frame.repaint();
    }

}
