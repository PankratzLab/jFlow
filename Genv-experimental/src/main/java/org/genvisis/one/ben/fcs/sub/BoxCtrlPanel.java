package org.genvisis.one.ben.fcs.sub;

import java.awt.Component;
import java.awt.Font;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTree;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeCellRenderer;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.DefaultTreeSelectionModel;
import javax.swing.tree.TreeNode;
import javax.swing.tree.TreePath;
import javax.swing.tree.TreeSelectionModel;

import org.genvisis.common.ArrayUtils;

import net.miginfocom.swing.MigLayout;

public class BoxCtrlPanel extends JPanel {

	/**
	* 
	*/
	private static final long serialVersionUID = 1L;
	public JTree tree;
	HashSet<String> actualData = new HashSet<String>();

	public void setData(String[] hdrs) {
		nodes.clear();
		ArrayList<String[]> headers = new ArrayList<String[]>();
		for (String s : hdrs) {
			if (s.startsWith("\"")) {
				s = s.substring(1);
			}
			if (s.endsWith("\"")) {
				s = s.substring(0, s.length() - 1);
			}
			s = s.split("\\|")[0].trim();
			headers.add(s.split("/"));
		}
		ArrayList<DefaultMutableTreeNode> rootNodes = new ArrayList<DefaultMutableTreeNode>();
		for (String[] hdr : headers) {
			String parent = hdr.length > 1 ? hdr[hdr.length - 2] : null;
			if (parent == null) {
				DefaultMutableTreeNode dmtn = nodes.get(hdr[0]);
				if (dmtn != null) {
					System.err.println("Error - duplicate root node: " + hdr[0]);
				} else {
					dmtn = new DefaultMutableTreeNode(hdr[0]);
					nodes.put(hdr[0], dmtn);
					rootNodes.add(dmtn);
					actualData.add(hdr[0]);
				}
			} else {
				String parentKey = ArrayUtils.toStr(ArrayUtils.subArray(hdr, 0, hdr.length - 1), "\t");
				DefaultMutableTreeNode dmtnParent = nodes.get(parentKey);
				if (dmtnParent == null) {
					dmtnParent = new DefaultMutableTreeNode(hdr[hdr.length - 2]);
					nodes.put(parentKey, dmtnParent);
					for (int i = hdr.length - 3; i >= 0; i--) {
						String key = ArrayUtils.toStr(ArrayUtils.subArray(hdr, 0, i + 1), "\t");
						DefaultMutableTreeNode nd = nodes.get(key);
						boolean found = true;
						if (nd == null) {
							nd = new DefaultMutableTreeNode(hdr[i]);
							nodes.put(key, nd);
							found = false;
						}
						key = ArrayUtils.toStr(ArrayUtils.subArray(hdr, 0, i + 2), "\t");
						nd.add(nodes.get(key));
						if (found) {
							break;
						}
					}
				}
				DefaultMutableTreeNode dmtn = new DefaultMutableTreeNode(hdr[hdr.length - 1]);
				dmtnParent.add(dmtn);
				nodes.put(ArrayUtils.toStr(hdr, "\t"), dmtn);
				actualData.add(ArrayUtils.toStr(hdr, "\t"));
			}
		}
		if (rootNodes.size() == 0) {
			System.err.println("Error - at least one root node must be specified!");
			return;
		}
		DefaultMutableTreeNode root = rootNodes.get(0);
		if (rootNodes.size() > 1) {
			DefaultMutableTreeNode dmtn = new DefaultMutableTreeNode();
			for (DefaultMutableTreeNode dmtnCh : rootNodes) {
				dmtn.add(dmtnCh);
			}
			root = dmtn;
		}
		DefaultTreeModel dtm = new DefaultTreeModel(root);
		DefaultTreeSelectionModel dtsm = new DefaultTreeSelectionModel() {
			/**
			* 
			*/
			private static final long serialVersionUID = 1L;

			@Override
			public void setSelectionPaths(TreePath[] pPaths) {
				ArrayList<TreePath> validPaths = new ArrayList<TreePath>();
				for (TreePath path : pPaths) {
					Object[] pathObjs = path.getPath();
					if (actualData.contains(ArrayUtils.toStr(pathObjs, "\t"))) { // only allow selections that have
																																	// data
						validPaths.add(path);
					}
				}
				super.setSelectionPaths(validPaths.toArray(new TreePath[validPaths.size()]));
			}

			@Override
			public void setSelectionPath(TreePath path) {
				Object[] pathObjs = path.getPath();
				if (actualData.contains(ArrayUtils.toStr(pathObjs, "\t"))) { // only allow selections that have
																																// data
					super.setSelectionPath(path);
				}
			}
		};
		dtsm.setSelectionMode(TreeSelectionModel.DISCONTIGUOUS_TREE_SELECTION);
		dtsm.addTreeSelectionListener(new TreeSelectionListener() {
			@Override
			public void valueChanged(TreeSelectionEvent e) {
				if (e.getNewLeadSelectionPath() == null) {
					return;
				}
				String newKey = ArrayUtils.toStr(e.getNewLeadSelectionPath().getPath(), "\t");
				if (actualData.contains(newKey)) {
					return;
				} else {
					int oldRow = tree.getRowForPath(e.getOldLeadSelectionPath());
					int newRow = tree.getRowForPath(e.getNewLeadSelectionPath());
					if (oldRow == -1) {
						boolean foundValid = false;
						do {
							newRow--;
							foundValid = actualData.contains(ArrayUtils.toStr(	tree.getPathForRow(newRow).getPath(),
																														"\t"));
						} while (tree.getPathForRow(newRow).getPathCount() >= e	.getNewLeadSelectionPath()
																																		.getPathCount()
											|| (!foundValid && newRow >= 0));
					} else if (oldRow > newRow) { // moving up
						boolean foundValid = false;
						do {
							newRow--;
							foundValid = actualData.contains(ArrayUtils.toStr(	tree.getPathForRow(newRow).getPath(),
																														"\t"));
						} while (!foundValid && newRow >= 0);
					} else { // moving down
						boolean foundValid = false;
						do {
							newRow++;
							foundValid = actualData.contains(ArrayUtils.toStr(	tree.getPathForRow(newRow).getPath(),
																														"\t"));
						} while (!foundValid && newRow < tree.getRowCount());
					}
					if (newRow < 0 || newRow == tree.getRowCount()) {
						// error!
						return;
					} else {
						tree.setSelectionRow(newRow);
					}
				}
			}
		});
		tree.setToggleClickCount(0);
		tree.setSelectionModel(dtsm);
		tree.setModel(dtm);
		tree.setCellRenderer(new CustomDefaultTreeCellRenderer());
		for (int i = 0; i < tree.getRowCount(); i++) {
			tree.expandRow(i);
		}
		tree.invalidate();
		tree.repaint();
	}

	public TreePath[] getSelectedPaths() {
		return tree.getSelectionPaths();
	}

	public void addTreeSelectionListener(TreeSelectionListener tsl) {
		tree.addTreeSelectionListener(tsl);
	}

	class CustomDefaultTreeCellRenderer extends DefaultTreeCellRenderer {
		/**
		* 
		*/
		private static final long serialVersionUID = 1L;

		@Override
		public Component getTreeCellRendererComponent(JTree tree, Object value, boolean sel,
																									boolean expanded, boolean leaf, int row,
																									boolean hasFocus) {
			boolean enabled = actualData == null
												|| tree.getPathForRow(row) == null	? true
																														: actualData.contains(ArrayUtils.toStr(tree.getPathForRow(row)
																																																	.getPath(),
																																															"\t")); // <--
																																																			// here
																																																			// is
																																																			// your
																																																			// logic
																																																			// for
																																																			// enable/disable
																																																			// cell

			Component treeCellRendererComponent = super.getTreeCellRendererComponent(	tree, value, sel,
																																								expanded, leaf, row,
																																								hasFocus);
			treeCellRendererComponent.setEnabled(enabled);

			return treeCellRendererComponent;
		}
	}

	/**
	 * Create the panel.
	 */
	public BoxCtrlPanel() {
		setLayout(new MigLayout("", "[grow]", "[grow]"));

		JScrollPane scrollPane = new JScrollPane();
		add(scrollPane, "cell 0 0,grow");

		tree = new JTree(new DefaultTreeModel(null));
		tree.setExpandsSelectedPaths(true);
		scrollPane.setViewportView(tree);
		tree.setFont(new Font("Arial", Font.PLAIN, 9));

	}

	HashMap<String, DefaultMutableTreeNode> nodes = new HashMap<String, DefaultMutableTreeNode>();

	public TreeNode getNodeForKey(String key) {
		return nodes.get(key);
	}

}
