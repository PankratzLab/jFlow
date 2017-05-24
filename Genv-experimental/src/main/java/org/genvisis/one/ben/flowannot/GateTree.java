package org.genvisis.one.ben.flowannot;

import java.util.HashMap;

import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.MutableTreeNode;

public final class GateTree {

	private GateTree() {}

	public static MutableTreeNode constructTree() {
		DefaultMutableTreeNode rootNode = null;
		HashMap<String, DefaultMutableTreeNode> map = new HashMap<>();
		for (String[] node : GATE_TREE) {
			if (node.length == 1) {
				AnnotatedImage ai = new AnnotatedImage(node[0], true);
				rootNode = new DefaultMutableTreeNode(ai);
				map.put(node[0], rootNode);
			} else {
				AnnotatedImage ai = new AnnotatedImage(node[1], false);
				DefaultMutableTreeNode tn = new DefaultMutableTreeNode(ai);
				map.get(node[0]).add(tn);
				map.put(node[1], tn);
			}
		}
		return rootNode;
	}

	static final String[][] GATE_DIMS = {
																			 {"SSC-A", "FSC-A"},
																			 {"FSC-H", "FSC-W"},
																			 {"PE-A"},
																			 {"CD3", "CD19"},
																			 {"CD27", "IgD"},
																			 {"CD8", "CD4"},
	};

	static final String[][] GATE_TREE = { // mapping[i][0] = parent, mapping[i][1] = child,
																				// if length == 1, @ root
			{"SSC-A v FSC-A"},
			{"SSC-A v FSC-A", "FSC-H v FSC-W"},
			{"FSC-H v FSC-W", "PE-A"},
			{"PE-A", "CD3 v CD19"},
			{"CD3 v CD19", "CD27 v IgD"},
			{"CD3 v CD19", "CD8 v CD4"},

	};

}
