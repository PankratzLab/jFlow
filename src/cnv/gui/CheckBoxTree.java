package cnv.gui;

import java.awt.*;
import java.awt.event.*;
import java.util.*;
import javax.swing.*;
import javax.swing.JTree.DynamicUtilTreeNode;
import javax.swing.event.*;
import javax.swing.tree.*;

public class CheckBoxTree extends JTree implements ItemListener {
	public static final long serialVersionUID = 1L;
	
	private JCheckBox[] selections;
//	private Branch root;
	private int count = 0;
	private boolean dynamic;
	
	public CheckBoxTree(String[][] names, int maxSelectable) {
		super(createTreeStructure(names));
		
		selections = new JCheckBox[maxSelectable];
		dynamic = false;

		setCellRenderer(new CheckBoxNodeRenderer());
		setCellEditor(new CheckBoxCellEditor(this));
		setEditable(true);
	}
	
	public CheckBoxTree(String[] namesOfBranches, String[] branchHandles, String[][] namesOfNodes, boolean[] active, int maxSelectable) {
		super(createTreeStructure(namesOfBranches, branchHandles, namesOfNodes, active));
		
		selections = new JCheckBox[maxSelectable];
		dynamic = false;

		setCellRenderer(new CheckBoxNodeRenderer());
		setCellEditor(new CheckBoxCellEditor(this));
		setEditable(true);
	}
	
	public void itemStateChanged(ItemEvent itemEvent) {
		JCheckBox checkbox, deselect;
		boolean found;
		int index;
		
		deselect = null;
		checkbox = (JCheckBox)itemEvent.getSource();
		System.out.println("State change ("+(++count)+"): "+checkbox.getName()+" is "+checkbox.isSelected());
		System.out.println("There are "+selections.length+" elements in selections");
		if (checkbox.isSelected()) {
			found = false;
			for (index = 0; index<selections.length && !found; index++) {
				if (selections[index] == null) {
					found = true;
				}
	        }
			if (!found) {
				deselect = selections[0];
				selections[0] = null;
				collapseSelections(selections);
			}
			selections[index-1] = checkbox;
		} else {
			for (int i = 0; i<selections.length; i++) {
				if (selections[i] == checkbox) {
					selections[i] = null;
				}
            }
			collapseSelections(selections);
		}
		if (deselect != null) {
			deselect.setSelected(false);
			repaint();
		}
		
		fireValueChanged(new TreeSelectionEvent(this, getPathForRow(0), true, getPathForRow(0), getPathForRow(0)));
	}
	
	public void addNode(String nameOfBranch, String branchHandle, String[] namesOfNodes, boolean[] active) {
		TreeModel model;
		Font font;
		Boolean booleanValue;
		boolean focusPainted;
		JCheckBox selection1, selection2;
		int row1, row2;
		
		selection1 = selections[0];
		selection2 = selections[1];
//		getExpandedState();
		row1=-1;
		row2=-1;
		for (int i=0; i<getRowCount(); i++) {
			if (isExpanded(i)) {
				if (row1==-1) {
					row1=i;
				} else {
					row2=i;
				}
			}
		}
		setSelectionPath(null);

		font = UIManager.getFont("Tree.font");
		booleanValue = (Boolean)UIManager.get("Tree.drawsFocusBorderAroundIcon");
		focusPainted = (booleanValue!=null)&&(booleanValue.booleanValue());
		
		model = getModel();
		
		Object root = model.getRoot();
		
		JCheckBox[] boxes = new JCheckBox[namesOfNodes.length];
		for (int j = 0; j<boxes.length; j++) {
			boxes[j] = new JCheckBox(namesOfNodes[j], false);
			boxes[j].setFont(font);
			boxes[j].setName(branchHandle+" "+j);
			boxes[j].setFocusPainted(focusPainted);
			boxes[j].setEnabled(active[j]);
        }
//        DynamicUtilTreeNode.createChildren(root, new Branch(nameOfBranch, boxes));
//		Object ob = new Branch(nameOfBranch, boxes);
        DynamicUtilTreeNode.createChildren((DefaultMutableTreeNode) root, new Branch[] {new Branch(nameOfBranch, boxes)});

        dynamic = true;
        ((DefaultTreeModel)model).reload();
        
		if (row1>=0) {
			expandRow(row1);
		}
		if (row2>=0) {
			expandRow(row2);
		}
		if (selection1 != null) {
			selection1.setSelected(true);
			// reexpand tree
//			setExpandedState(anchorPath, focusPainted);
//			for (int i=0; i<getRowCount(); i++) {
//				expandRow(i);
//			}
//		    DefaultMutableTreeNode  rooot;
//		    rooot = (DefaultMutableTreeNode) getModel().getRoot();
//		    scrollPathToVisible(new TreePath(rooot.getLastLeaf().getPath()));
		}
		if (selection2 != null) {
			selection2.setSelected(true);
		}
	}
	
//	public void deleteNode(String nameOfBranch, String branchHandle) {
//		int index=-1;
//		DefaultMutableTreeNode root, branch;
//		DefaultTreeModel model;
//		
//		model = (DefaultTreeModel)getModel();		
//		root = (DefaultMutableTreeNode)model.getRoot();
//		for (int i=0; i<root.getChildCount(); i++ ) {
//			branch = (DefaultMutableTreeNode)getModel().getChild(getModel().getRoot(), i);
//			if (branch.toString().equals(nameOfBranch)) {
//				Branch br = (Branch)branch.getUserObject();
//				Object ob = br.firstElement();
//				if (((JCheckBox)ob).getName().split("[\\s]+")[0].startsWith(branchHandle)) {
//					System.out.println("Branch '"+branch+"' index="+i);
//					index = i;
//					break;
//				}
//			}
//		}
//		if (index != -1) {
//			root.remove(index);
//		} else {
//			System.err.println("Branch "+nameOfBranch+" not found.");
//		}
//		
//        model.reload();
//	}
	
	public void deleteSelectedNode() {
		DefaultMutableTreeNode selectedNode;
		JCheckBox selection1, selection2;
		int row1, row2, index;
		DefaultMutableTreeNode root;	//, branch;
		DefaultTreeModel model;
		
		selectedNode = (DefaultMutableTreeNode)getLastSelectedPathComponent();
		if (selectedNode == null) {
			System.err.println("Error - nothing selected");
			return;
		}

		selection1 = selections[0];
		selection2 = selections[1];
		setSelectionPath(null);
		row1=-1;
		row2=-1;
		for (int i=0; i<getRowCount(); i++) {
			if (isExpanded(i)) {
				if (row1==-1) {
					row1=i;
				} else {
					row2=i;
				}
			}
		}


//		index = -1;
		index = getSelectedPathComponent();
		model = (DefaultTreeModel)getModel();		
		root = (DefaultMutableTreeNode)model.getRoot();
//		for (int i=0; i<root.getChildCount(); i++ ) {
////			branch = (DefaultMutableTreeNode)getModel().getChild(getModel().getRoot(), i);
//			branch = (DefaultMutableTreeNode)model.getChild(root, i);
//			if (branch == selectedNode) {
//				index = i;
//			}
//			
//			if (selectedNode.getUserObject() instanceof JCheckBox) {
//				Branch br = (Branch)branch.getUserObject();
//				for (int j=0; j<branch.getChildCount(); j++ ) {
//					JCheckBox box = (JCheckBox)br.elementAt(j);
//					if (box == (JCheckBox)(selectedNode.getUserObject())){
//						index = i;
//					}
//				}
//			}
//		}
		if (index != -1) {
			root.remove(index);
		} else {
			System.err.println("Branch "+selectedNode+" not found.");
		}
		
        model.reload();
        
        if (row1>=-1) {
        	if (row1<index) {
        		expandRow(row1);
        	} else if (row1>index) {
        		expandRow(row1-1);
        	}
        }
        
        if (row2>=-1) {
        	if (row2<index) {
        		expandRow(row2);
        	} else if (row2>index) {
        		expandRow(row2-1);
        	}
        }
        
        if (selection1 != null) {
			selection1.setSelected(true);
		}
		if (selection2 != null) {
			selection2.setSelected(true);
		}
	}
	
	public static void collapseSelections(JCheckBox[] selections) {
		JCheckBox trav;
		int index;
		
		index = 0;
		for (int i = 0; i<selections.length; i++) {
			if (selections[i] != null) {
				trav = selections[i];
				selections[i] = null;
				selections[index++] = trav;
			}
        }
	}
	
	public int[][] getSelectionIndices() {
		int[][] indices;
		String[] line;
		
		if (dynamic) {
			System.err.println("Error - addNode was used; cannot getSelectionIndices");
			return null;
		}
		
		indices = new int[selections.length][2];
		for (int i = 0; i<selections.length; i++) {
			if (selections[i] == null) {
				indices[i][0] = indices[i][1] = -1;
			} else {
				line = selections[i].getName().split("[\\s]+");
				indices[i][0] = Integer.parseInt(line[0]);
				indices[i][1] = Integer.parseInt(line[1]);
			}
        }
		
		return indices;
	}

	public String[][] getSelectionValues() {
		String[][] selectedValues;
		String[] line;
		
		selectedValues = new String[selections.length][2];
		for (int i = 0; i<selections.length; i++) {
			if (selections[i] == null) {
				selectedValues[i][0] = selectedValues[i][1] = null;
			} else {
				line = selections[i].getName().split("[\\s]+");
				selectedValues[i][0] = line[0];
				selectedValues[i][1] = line[1];
			}
        }
		
		return selectedValues;
	}

	public int getSelectedPathComponent() {
		DefaultMutableTreeNode selectedNode;
		int index;
		DefaultMutableTreeNode root, branch;
		DefaultTreeModel model;
		
		selectedNode = (DefaultMutableTreeNode)getLastSelectedPathComponent();
		if (selectedNode == null) {
			System.err.println("Error - nothing selected");
			return -1;
		}

		index=-1;
		model = (DefaultTreeModel)getModel();		
		root = (DefaultMutableTreeNode)model.getRoot();
		for (int i=0; i<root.getChildCount(); i++ ) {
			branch = (DefaultMutableTreeNode)model.getChild(root, i);
			if (branch == selectedNode) {
				index = i;
			}
			
			if (selectedNode.getUserObject() instanceof JCheckBox) {
				Branch br = (Branch)branch.getUserObject();
				for (int j=0; j<branch.getChildCount(); j++ ) {
					JCheckBox box = (JCheckBox)br.elementAt(j);
					if (box == (JCheckBox)(selectedNode.getUserObject())){
						index = i;
					}
				}
			}
		}
		return index;
	}
	
	public String getSelectedPathComponentName() {
		int index = getSelectedPathComponent();
		if (index>=0) {
			return selections[index].getName().split("[\\s]+")[0];
		} else {
			return null;
		}
	}
	
//	public static void main(String[] args) {
//		JFrame frame = new JFrame("CheckBox Tree");
//
//		String[][] options = {{"Accessibility", "Move system caret with focus/selection changes", "Always expand alt text for images"},
//		{"Browsing", "Notify when downloads complete", "Disable script debugging", "Use AutoComplete", "Browse in a new process"}};
//
//		JScrollPane scrollPane = new JScrollPane(new CheckBoxTree(options, 2));
//		frame.getContentPane().add(scrollPane, BorderLayout.CENTER);
//		frame.setSize(300, 150);
//		frame.setVisible(true);
//		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
//	}
	
	public static Branch createTreeStructure(String[][] names) {
		Branch[] branches = new Branch[names.length];
		JCheckBox[] boxes;
		Font font;
		Boolean booleanValue;
		boolean focusPainted;

		font = UIManager.getFont("Tree.font");
		booleanValue = (Boolean)UIManager.get("Tree.drawsFocusBorderAroundIcon");
		focusPainted = (booleanValue!=null)&&(booleanValue.booleanValue());
		for (int i = 0; i<names.length; i++) {
			boxes = new JCheckBox[names[i].length-1];
			for (int j = 0; j<boxes.length; j++) {
				boxes[j] = new JCheckBox(names[i][j+1], false);
				boxes[j].setFont(font);
				boxes[j].setName(i+" "+j);
				boxes[j].setFocusPainted(focusPainted);
            }
			branches[i] = new Branch(names[i][0], boxes);
        }

		return new Branch("Root", branches);
	}
	
	public static Branch createTreeStructure(String[] namesOfBranches, String[] branchHandles, String[][] namesOfNodes, boolean[] active) {
		if (namesOfBranches==null || namesOfBranches.length==0) {
			return null;
		}
		Branch[] branches = new Branch[namesOfNodes.length];
		JCheckBox[] boxes;
		Font font;
		Boolean booleanValue;
		boolean focusPainted;

		font = UIManager.getFont("Tree.font");
		booleanValue = (Boolean)UIManager.get("Tree.drawsFocusBorderAroundIcon");
		focusPainted = (booleanValue!=null)&&(booleanValue.booleanValue());
		for (int i = 0; i<namesOfBranches.length; i++) {
			boxes = new JCheckBox[namesOfNodes[i].length];
			for (int j = 0; j<boxes.length; j++) {
				boxes[j] = new JCheckBox(namesOfNodes[i][j], false);
				boxes[j].setFont(font);
				boxes[j].setName(branchHandles[i]+" "+j);
				boxes[j].setFocusPainted(focusPainted);
				boxes[j].setEnabled(active[j]);
            }
			branches[i] = new Branch(namesOfBranches[i], boxes);
        }

		return new Branch("Root", branches);
	}
}


class CheckBoxNodeRenderer implements TreeCellRenderer {
//	private JCheckBox leafRenderer = new JCheckBox();
	private JCheckBox latestLeaf;
	private DefaultTreeCellRenderer branchRenderer = new DefaultTreeCellRenderer();

	private Color selectionForeground, selectionBackground, textForeground, textBackground;

	protected JCheckBox getLeafRenderer() {
		return latestLeaf;
	}

	public CheckBoxNodeRenderer() {
		selectionForeground = UIManager.getColor("Tree.selectionForeground");
		selectionBackground = UIManager.getColor("Tree.selectionBackground");
		textForeground = UIManager.getColor("Tree.textForeground");
		textBackground = UIManager.getColor("Tree.textBackground");

		branchRenderer.setOpenIcon(null);
		branchRenderer.setClosedIcon(null);
	}

	public Component getTreeCellRendererComponent(JTree tree, Object value, boolean selected, boolean expanded, boolean leaf, int row, boolean hasFocus) {
		if (leaf) {
//			String stringValue = tree.convertValueToText(value, selected, expanded, leaf, row, false);
			
//			leafRenderer.setText(stringValue);
//			leafRenderer.setSelected(false);
//
//			leafRenderer.setEnabled(tree.isEnabled());
//

			if ((value!=null)&&(value instanceof DefaultMutableTreeNode)) {
				Object userObject = ((DefaultMutableTreeNode)value).getUserObject();
				if (userObject instanceof JCheckBox) {
					JCheckBox node = (JCheckBox)userObject;
//					leafRenderer.setText(node.getText());
//					leafRenderer.setSelected(node.isSelected());
					
					if (selected) {
						node.setForeground(selectionForeground);
						node.setBackground(selectionBackground);
					} else {
						node.setForeground(textForeground);
						node.setBackground(textBackground);
					}
					latestLeaf = node;
					return node;
				}
			}
			return new JCheckBox("messed", true);
		} else {
			return branchRenderer.getTreeCellRendererComponent(tree, value, selected, expanded, leaf, row, hasFocus);
		}
	}
}

class CheckBoxCellEditor extends AbstractCellEditor implements TreeCellEditor {
	public static final long serialVersionUID = 1L;

	CheckBoxNodeRenderer renderer = new CheckBoxNodeRenderer();
	ChangeEvent changeEvent = null;
	CheckBoxTree checkBoxTree;
	int count;

	public CheckBoxCellEditor(CheckBoxTree checkBoxTree) {
		this.checkBoxTree = checkBoxTree;
	}

	public Object getCellEditorValue() {
//		JCheckBox checkbox = renderer.getLeafRenderer();
//		return new JCheckBox(checkbox.getText(), checkbox.isSelected());
//		return new JCheckBox("damnit", false);
		return renderer.getLeafRenderer();
	}

	@Override
	public boolean isCellEditable(EventObject event) {
		if (event instanceof MouseEvent) {
			MouseEvent mouseEvent = (MouseEvent)event;
			TreePath path = checkBoxTree.getPathForLocation(mouseEvent.getX(), mouseEvent.getY());
			return ((DefaultMutableTreeNode)path.getLastPathComponent()).isLeaf();
		}
		return false;
	}

	public Component getTreeCellEditorComponent(JTree tree, Object value, boolean selected, boolean expanded, boolean leaf, int row) {
		Component editor = renderer.getTreeCellRendererComponent(tree, value, true, expanded, leaf, row, true);

//		ItemListener itemListener = new ItemListener() {
//			public void itemStateChanged(ItemEvent itemEvent) {
//				if (stopCellEditing()) {
//					fireEditingStopped();
//					System.out.println(++count);
//				} else {
//					System.out.println("didn't stop");
//				}
//			}
//		};
		if (editor instanceof JCheckBox) {
			JCheckBox checkbox = ((JCheckBox)editor);
			if (checkbox.getItemListeners().length == 0) {
//				checkbox.addItemListener(itemListener);
				checkbox.addItemListener(checkBoxTree);
			}
		}

		return editor;
	}
}

class Branch extends Vector<Object> {
	public static final long serialVersionUID = 1L;
	String name;

	public Branch(String name) {
		this.name = name;
	}

	public Branch(String name, Object[] elements) {
		this.name = name;
		for (int i = 0, n = elements.length; i<n; i++) {
			add(elements[i]);
		}
	}

	@Override
	public String toString() {
		return name;
	}
}
