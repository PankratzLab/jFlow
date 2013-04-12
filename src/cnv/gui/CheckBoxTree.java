package cnv.gui;

import java.awt.*;
import java.awt.event.*;
import java.util.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.tree.*;

public class CheckBoxTree extends JTree implements ItemListener {
	public static final long serialVersionUID = 1L;
	
	private JCheckBox[] selections;
	
	public CheckBoxTree(String[][] names, int maxSelectable) {
		super(createTreeStructure(names));
		
		selections = new JCheckBox[maxSelectable];

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
//		System.out.println("State change ("+(++count)+"): "+checkbox.getName()+" is "+checkbox.isSelected());
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

	public static void main(String[] args) {
		JFrame frame = new JFrame("CheckBox Tree");

		String[][] options = {{"Accessibility", "Move system caret with focus/selection changes", "Always expand alt text for images"},
		{"Browsing", "Notify when downloads complete", "Disable script debugging", "Use AutoComplete", "Browse in a new process"}};

		JScrollPane scrollPane = new JScrollPane(new CheckBoxTree(options, 2));
		frame.getContentPane().add(scrollPane, BorderLayout.CENTER);
		frame.setSize(300, 150);
		frame.setVisible(true);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}
	
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
