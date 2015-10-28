package cnv.gui;

import java.awt.*;
import java.awt.event.*;
import java.util.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.tree.*;

import common.IntVector;

public class CheckBoxTree extends JTree implements ItemListener {
	public static final long serialVersionUID = 1L;
	
	private JCheckBox[] selections;
//	private Branch root;
//	private int count = 0;
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

	/**
	 * Function to update the maximum number of checkbox which can be selected after creating the tree
	 * 
	 * @param maxSelectable
	 *            the new maximum number of selectable checkboxes
	 */
	public void setMaxSelections(int maxSelectable) {
		// get new selection with specified size
		selections = Arrays.copyOf(selections, maxSelectable);
	}

	public void itemStateChanged(ItemEvent itemEvent) {
		JCheckBox checkbox, deselect;
		boolean found;
		int index;
		
		deselect = null;
		checkbox = (JCheckBox)itemEvent.getSource();
//		System.out.println("State change ("+(++count)+"): "+checkbox.getName()+" is "+checkbox.isSelected());
//		System.out.println("There are "+selections.length+" elements in selections");
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


	@SuppressWarnings("unchecked")
	public DefaultMutableTreeNode searchNode(String fileName, String nodeStr)
	{
		DefaultTreeModel model = (DefaultTreeModel)getModel();
		DefaultMutableTreeNode node = (DefaultMutableTreeNode) model.getRoot();
		Enumeration<DefaultMutableTreeNode> topLevel = node.children();
		while(topLevel.hasMoreElements()) {
			DefaultMutableTreeNode tpNd = topLevel.nextElement();
			if(((Branch)tpNd.getUserObject()).toString().equals(fileName)) {
				Enumeration<DefaultMutableTreeNode> subElements = tpNd.breadthFirstEnumeration();
				while(subElements.hasMoreElements()) {
					node = subElements.nextElement();
					if(node.getUserObject() instanceof JCheckBox){
						JCheckBox thisNode= (JCheckBox) node.getUserObject();
						//System.out.println("this:" + thisNode.getText());
						if(nodeStr.equals(thisNode.getText())) {
							return node;
						}
					}
				}
			}
		}
		
		//tree node with string node found return null
		return null;
	}
	
	@SuppressWarnings("unchecked")
	public DefaultMutableTreeNode searchNode(String nodeStr)
	{
		DefaultTreeModel model = (DefaultTreeModel)getModel();
		DefaultMutableTreeNode node = (DefaultMutableTreeNode) model.getRoot();		
		Enumeration<DefaultMutableTreeNode> enumeration = node.breadthFirstEnumeration();
		while(enumeration.hasMoreElements()) {
			node = enumeration.nextElement();
			if(node.getUserObject() instanceof JCheckBox){
				JCheckBox thisNode= (JCheckBox) node.getUserObject();
				//System.out.println("this:" + thisNode.getText());
				if(nodeStr.equals(thisNode.getText())) {
					return node;
				}
			}
		}
		//tree node with string node found return null
		return null;
	}

	/**
	 * Function to perform a given action on a given checkbox
	 *
	 * @param checkboxName
	 *            the name of the checkbox on which action has to be performed
	 * @param action
	 *            the action to be performed which is either SELECTED or DESELECTED from {@link ItemEvent}
	 */
	public void performCheckBoxAction(String checkboxName, int action) {
		// try to get the checkbox from the CheckBoxTree
		DefaultMutableTreeNode searchNode = searchNode(checkboxName);
		checkBoxAction(searchNode, action);
	}
	
	public void performCheckBoxAction(String fileName, String checkboxName, int action) {
		DefaultMutableTreeNode searchNode = searchNode(fileName, checkboxName);
		checkBoxAction(searchNode, action);
	}
	
	private void checkBoxAction(DefaultMutableTreeNode searchNode, int action) {
		if (searchNode != null) {
			JCheckBox thisCheckBox = (JCheckBox) searchNode.getUserObject();
			if (action == ItemEvent.SELECTED) {
				if (thisCheckBox.isSelected()) { // if action is selected and the checkbox is already selected then
					// deselect first
					thisCheckBox.setSelected(false);
					this.itemStateChanged(new ItemEvent(thisCheckBox, ItemEvent.ITEM_LAST, thisCheckBox, ItemEvent.DESELECTED));
				}
				thisCheckBox.setSelected(true); // then select the checkbox again
			} else if (action == ItemEvent.DESELECTED) {
				thisCheckBox.setSelected(false);
			}
			itemStateChanged(new ItemEvent(thisCheckBox, ItemEvent.ITEM_LAST, thisCheckBox, action));
			repaint();
		}
	}

	/**
	 *
 	 * @param nameOfBranch
	 * @param branchHandle
	 * @param namesOfNodes
	 * @param active
	 * @param mouseListener: the mouse listener for the checkboxes. Should be null if we don't want any action to be performed on mouse click. itemStateChange will work irrespective of mouseListener
	 */
	public void addNode(String nameOfBranch, String branchHandle, String[] namesOfNodes, boolean[] active, MouseListener mouseListener) {
		TreeModel model;
		Font font;
		Boolean booleanValue;
		boolean focusPainted;
		JCheckBox[] travSelections;
		IntVector expansions;
		
		travSelections = new JCheckBox[selections.length];
		for (int i = 0; i < selections.length; i++) {
			travSelections[i] = selections[i];
		}

		expansions = new IntVector();
		for (int i=0; i<getRowCount(); i++) {
			if (isExpanded(i)) {
				expansions.add(i);
				collapseRow(i);
			}
		}

		setSelectionPath(null);

		font = UIManager.getFont("Tree.font");
		booleanValue = (Boolean)UIManager.get("Tree.drawsFocusBorderAroundIcon");
		focusPainted = (booleanValue!=null)&&(booleanValue.booleanValue());
		
		model = getModel();
		
		Object root = model.getRoot();
//		DynamicUtilTreeNode root = (DynamicUtilTreeNode)model.getRoot();
//		root.setAllowsChildren(true);
		
		JCheckBox[] boxes = new JCheckBox[namesOfNodes.length];
		for (int j = 0; j<boxes.length; j++) {
			boxes[j] = new JCheckBox(namesOfNodes[j], false);
			boxes[j].setFont(font);
			boxes[j].setName(branchHandle+" "+j);
			boxes[j].setFocusPainted(focusPainted);
			boxes[j].setEnabled(active[j]);
			boxes[j].addMouseListener(mouseListener);
        }
//        DynamicUtilTreeNode.createChildren(root, new Branch(nameOfBranch, boxes));
//		Object ob = new Branch(nameOfBranch, boxes);
		((DefaultMutableTreeNode)root).setAllowsChildren(true);
        DynamicUtilTreeNode.createChildren((DefaultMutableTreeNode) root, new Branch[] {new Branch(nameOfBranch, boxes)});

        dynamic = true;
        ((DefaultTreeModel)model).reload();
        
        for (int i = expansions.size()-1; i >= 0 ; i--) {
    		expandRow(expansions.elementAt(i));
		}

		for (int i = 0; i < selections.length; i++) {
			if (travSelections[i] != null) {
				travSelections[i].setSelected(true);
			}
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
	
//	public void deleteSelectedNode() {
//		DefaultMutableTreeNode selectedNode;
//		int[][] selectionTmp;
//		String[] selects;
//		int row1, row2, index;
//		DefaultMutableTreeNode root;	//, branch;
//		DefaultTreeModel model;
//		
//		selectedNode = (DefaultMutableTreeNode)getLastSelectedPathComponent();
//		if (selectedNode == null) {
//			System.err.println("Error - No root of a branch is selected.");
//			return;
//		}
//		index = getSelectedPathComponent();
//		model = (DefaultTreeModel)getModel();		
//		root = (DefaultMutableTreeNode)model.getRoot();
//
//		selectionTmp = new int[selections.length][2];
//		selects = new String[2];
//		for (int i=0; i<selectionTmp.length; i++) {
//			selects[i] = selections[i].getName();
//			
//			selectionTmp[i][0] = selections[i].getX();
//			selectionTmp[i][1] = selections[i].getY();
//			//TODO What's the relationship between selections[] and getSelectedPathComponent()???
//			if ( getSelectionValues()[i][0].equals(getSelectedPathComponent()) ) {
//				selectionTmp[i] = null;
//			}
//		}
//		setSelectionPath(null);
//		selections[0]=null;
//		selections[1]=null;
//		row1=-1;
//		row2=-1;
//		for (int i=0; i<getRowCount(); i++) {
//			if (isExpanded(i)) {
//				if (row1==-1) {
//					row1=i;
//				} else {
//					row2=i;
//				}
//			}
//		}
//
////		for (int i=0; i<root.getChildCount(); i++ ) {
//////			branch = (DefaultMutableTreeNode)getModel().getChild(getModel().getRoot(), i);
////			branch = (DefaultMutableTreeNode)model.getChild(root, i);
////			if (branch == selectedNode) {
////				index = i;
////			}
////			
////			if (selectedNode.getUserObject() instanceof JCheckBox) {
////				Branch br = (Branch)branch.getUserObject();selections
////				for (int j=0; j<branch.getChildCount(); j++ ) {
////					JCheckBox box = (JCheckBox)br.elementAt(j);
////					if (box == (JCheckBox)(selectedNode.getUserObject())){
////						index = i;
////					}
////		index = -1;
////				}
////			}
////		}
//		if (index != -1) {
//			root.remove(index);
//		} else {
//			System.err.println("Branch "+selectedNode+" not foungetSelectionIndicesd.");
//		}
//		
//        model.reload();
//        
//        if (row1>=-1) {
//        	if (row1<index) {
//        		expandRow(row1);
//        	} else if (row1>index) {
//        		expandRow(row1-1);
//        	}
//        }
//        
//        if (row2>=-1) {
//        	if (row2<index) {
//        		expandRow(row2);
//        	} else if (row2>index) {
//        		expandRow(row2-1);
//        	}
//        }
//        
//		for (int i=0; i<selectionTmp.length; i++) {
//	        if (selectionTmp[i] != null) {
//	        	setSelectionPath(getPathForLocation(selectionTmp[i][0],selectionTmp[i][1]));
//				selections[i].setSelected(true);//TODO
//			}
//	        if (selects[i] != null) {
//	        	
//	        	setSelectionPath(getPathForLocation(selectionTmp[i][0],selectionTmp[i][1]));
//				selections[i].setSelected(true);//TODO
//			}
//		}
//	}
//	

	// TODO if node to be deleted has selected values, then it may fail down stream in the parent appication (e.g., TwoDPlot), tried to make a work around below (currently commented out), but it does not work 
	public void deleteSelectedNode() {
		DefaultMutableTreeNode selectedNode;
		int index;
		DefaultMutableTreeNode root;	//, branch;
		DefaultTreeModel model;
		IntVector expansions;
//		JCheckBox[] travSelections;
		
		selectedNode = (DefaultMutableTreeNode)getLastSelectedPathComponent();
		if (selectedNode == null) {
			System.err.println("Error - No root of a branch is selected.");
			return;
		}
		index = getSelectedPathComponent();
		model = (DefaultTreeModel)getModel();		
		root = (DefaultMutableTreeNode)model.getRoot();
		
//		travSelections = new JCheckBox[selections.length];
//		for (int i = selections.length-1; i >= 0 ; i--) {
//			travSelections[i] = selections[i];
//			if (selections[i] != null) {
//				selections[i].setSelected(false);
////				selections[i] = null;
//			}
//		}

		setSelectionPath(null);

		// Reserve the Tree Expansion/Collapse status
		expansions = new IntVector();
		for (int i=0; i<getRowCount(); i++) {
			if (isExpanded(i)) {
				if (i != index) {
					expansions.add(i);
				}
				collapseRow(i);
			}
		}

		root.remove(index);
        model.reload();
        
//		for (int i = 0; i < selections.length; i++) {
//			if (travSelections[i] != null) {
//				selections[i] = travSelections[i];
//				selections[i].setSelected(true);
//			} else {
//				System.out.println("nothing to select for index "+i);
//			}
//		}
        
        // Restore the tree expansion/collapse status
        for (int i = expansions.size()-1; i >= 0 ; i--) {
    		expandRow(expansions.elementAt(i)-(expansions.elementAt(i)<index?0:1));
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
			System.err.println("Error - No root of a branch is selected");
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
	private DefaultTreeCellRenderer branchRenderer;

	private Color selectionForeground, selectionBackground, textForeground, textBackground;

	protected JCheckBox getLeafRenderer() {
		return latestLeaf;
	}

	public CheckBoxNodeRenderer() {
		branchRenderer = new DefaultTreeCellRenderer();
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

	private CheckBoxNodeRenderer renderer;
//	private ChangeEvent changeEvent;
	private CheckBoxTree checkBoxTree;
//	private int count;

	public CheckBoxCellEditor(CheckBoxTree newCheckBoxTree) {
		renderer = new CheckBoxNodeRenderer();
		changeEvent = null;
		checkBoxTree = newCheckBoxTree;
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
	private String name;

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
