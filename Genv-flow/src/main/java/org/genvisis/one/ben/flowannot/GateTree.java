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
				AnnotatedImage ai = new AnnotatedImage(node[0], false);
				DefaultMutableTreeNode tn = new DefaultMutableTreeNode(ai);
				map.get(node[1]).add(tn);
				map.put(node[0], tn);
			}
		}
		return rootNode;
	}


	static final String[][] GATE_TREE = { // mapping[i][0] = parent, mapping[i][1] = child,
																				// if length == 1, @ root
			{"boundary"},
			{"nonDebris", "boundary"},
			{"lymph", "nonDebris"},
			{"FSC-W+", "lymph"},
			{"Singlets", "lymph"},
			{"PE-A", "Singlets"},
			{"CD3+", "PE-A"},
			{"CD4", "CD3+"},
			{"ActivatedCD4", "CD4"},
			{"CCR7+CD45RA+", "CD4"},
			{"CCR7+CD45RA-", "CD4"},
			{"CCR7-CD45RA+", "CD4"},
			{"CCR7-CD45RA-", "CD4"},
			{"CD45RA+", "CD4"},
			{"CCR7gate", "CD45RA+"},
			{"CD45RA-", "CD4"},
			{"HLA-DR+", "CD4"},
			{"CD4+", "CD3+"},
			{"CD4-", "CD3+"},
			{"CD8+", "CD4-"},
			{"CD8", "CD4"},
			{"ActivatedCD8", "CD8"},
			{"CCR7+CD45RA+", "CD8"},
			{"CCR7+CD45RA-", "CD8"},
			{"CCR7-CD45RA+", "CD8"},
			{"CD8EFFCD27Gate", "CCR7-CD45RA+"},
			{"CD8EFFCD28Gate", "CCR7-CD45RA+"},
			{"CD28+CD27+", "CCR7-CD45RA+"},
			{"CD28+CD27-", "CCR7-CD45RA+"},
			{"CD28-CD27+", "CCR7-CD45RA+"},
			{"CD28-CD27-", "CCR7-CD45RA+"},
			{"CCR7-CD45RA-", "CD8"},
			{"CD8MEMCD27Gate", "CCR7-CD45RA-"},
			{"CD8MEMCD28Gate", "CCR7-CD45RA-"},
			{"CD28+CD27+", "CCR7-CD45RA-"},
			{"CD28+CD27-", "CCR7-CD45RA-"},
			{"CD28-CD27+", "CCR7-CD45RA-"},
			{"CD28-CD27-", "CCR7-CD45RA-"},
			{"CD45RA+", "CD8"},
			{"CCR7gate", "CD45RA+"},
			{"CD45RA-", "CD8"},
			{"CCR7gate", "CD45RA-"},
			{"HLA-DR+", "CD8"},
			{"CD8-", "CD3+"},
			{"CD3-", "PE-A"},
			{"CD19gate", "CD3-"},
			{"CD27gate", "CD19gate"},
			{"IgD+CD27+", "CD19gate"},
			{"IgD+CD27-", "CD19gate"},
			{"IgD-CD27+", "CD19gate"},
			{"IgD-CD27-", "CD19gate"},
			{"IgDgate", "CD19gate"},
			{"PE-A+", "PE-A"},
	};

}
