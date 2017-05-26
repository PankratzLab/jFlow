package org.genvisis.one.ben.fcs.gating;

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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JScrollPane;
import javax.swing.JTree;
import javax.swing.ScrollPaneConstants;
import javax.swing.SwingConstants;
import javax.swing.border.BevelBorder;
import javax.swing.border.Border;
import javax.swing.border.EtchedBorder;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.TreeNode;
import javax.swing.tree.TreePath;

import org.genvisis.one.ben.fcs.FCSPlot;

import net.miginfocom.swing.MigLayout;

public class GateTreePanel extends JPanel {

	/**
	* 
	*/
	private static final long serialVersionUID = 1L;
	Gating gating;
	private final JTree tree;
	private final JPanel topPanel;
	private final JScrollPane scrollPane;
	private final FCSPlot plot;
	private final HashMap<DefaultMutableTreeNode, Gate> gateMap = new HashMap<DefaultMutableTreeNode, Gate>();
	private final HashMap<Gate, DefaultMutableTreeNode> nodeMap = new HashMap<Gate, DefaultMutableTreeNode>();
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
			@Override
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
				Rectangle rect = new Rectangle(breadcrumbPanel.getWidth() - 10, 10,
																			 breadcrumbPanel.getWidth(), 10);
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
		tree.setRootVisible(true);
		tree.setShowsRootHandles(true);
		tree.addTreeSelectionListener(treeListener);
	}

	volatile boolean selecting = false;

	public void selectGate(Gate selected) {
		if (selected != null) {
			selecting = true;
			expandAllNodes();
			DefaultMutableTreeNode selNode = nodeMap.get(selected);
			TreeNode[] path = ((DefaultTreeModel) tree.getModel()).getPathToRoot(selNode);
			tree.setSelectionPath(new TreePath(path));
			selecting = false;
		}
	}

	public void resetGating(Gating gating, Gate selectedGate) {
		this.gating = gating;
		ArrayList<Gate> roots = gating.getRootGates();

		DefaultMutableTreeNode rootNode = new DefaultMutableTreeNode();
		for (Gate g : roots) {
			addGatesToTree(rootNode, g);
		}

		breadcrumbPanel.removeAll();
		DefaultTreeModel newModel = new DefaultTreeModel(rootNode, true);
		tree.setModel(newModel);
		expandAllNodes();

		if (selectedGate != null) {
			selectGate(selectedGate);
		}

		repaint();
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

	private void addGatesToTree(DefaultMutableTreeNode root, Gate g) {
		StringBuilder ident = new StringBuilder();
		ident.append(g.getName() == null || "".equals(g.getName()) ? g.getID() : g.getName());

		ident.append(" (");
		ident.append(g.getXDimension().paramName);
		if (g.getYDimension() != null) {
			ident.append(" v ");
			ident.append(g.getYDimension().paramName);
		}
		ident.append(")");
		DefaultMutableTreeNode child = new DefaultMutableTreeNode(ident.toString(), true);
		gateMap.put(child, g);
		nodeMap.put(g, child);
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
				showing = false;
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
				String chNm = (String) (i == path.length
																		 - 1 ? ""
																				 : ((DefaultMutableTreeNode) path[i + 1]).getUserObject());
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
			Rectangle rect = new Rectangle(breadcrumbPanel.getWidth() - 10, 10,
																		 breadcrumbPanel.getWidth(), 10);
			((JComponent) breadcrumbPanel.getParent()).scrollRectToVisible(rect);
			if (!selecting) {
				plot.gateSelected(gateMap.get(path[path.length - 1]), false);
			}
		}
	};
	private final JPanel breadcrumbPanel;
	private final Border mouseoverBorder = BorderFactory.createLineBorder(Color.black, 1);
	private final Border breadcrumbBorder = BorderFactory.createEmptyBorder(0, 1, 0, 1);

	private JLabel getBreadcrumbLabel(String txt, String me,
																		final List<DefaultMutableTreeNode> sibs, JLabel prevLbl) {
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
				/**
				* 
				*/
				private static final long serialVersionUID = 1L;

				@Override
				public void actionPerformed(ActionEvent e) {
					tree.setSelectionPath(new TreePath(((DefaultTreeModel) tree.getModel()).getPathToRoot(sibs.get(ind))));
				}
			};
			JMenuItem jmi = new JMenuItem(act);
			popup.add(jmi);
			if (c < sibs.size() - 1) {
				sb.append("<br />");
			}
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

		@Override
		public void mouseClicked(MouseEvent e) {
			super.mouseClicked(e);
			((JComponent) e.getSource()).getComponentPopupMenu().show(((JComponent) e.getSource()),
																																e.getX(), e.getY());
		};
	};
	private final JScrollPane scrollPane_1;
	private final JPopupMenu treePopup;

}
