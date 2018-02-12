package org.genvisis.one.ben.flowannot;

import java.util.HashMap;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.MutableTreeNode;
import org.genvisis.one.ben.flowannot.IAnnotator.PANEL;

public final class GateTree {

  private GateTree() {}

  public static MutableTreeNode constructTree(PANEL panel) {
    String[][] tree = panel == PANEL.PANEL_1 ? GATE_TREE_PANEL_1 : GATE_TREE_PANEL_2;
    DefaultMutableTreeNode rootNode = null;
    HashMap<String, DefaultMutableTreeNode> map = new HashMap<>();
    for (String[] node : tree) {
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

  static final String[][] GATE_TREE_PANEL_2 = {{"boundary"}, {"CD45-", "boundary"},
                                               {"CD45+", "boundary"}, {"PE-A+", "CD45+"},
                                               {"Live immune cells (CD45+ PE-)", "CD45+"},
                                               {"FSC-H+", "Live immune cells (CD45+ PE-)"},
                                               {"SingletsH", "Live immune cells (CD45+ PE-)"},
                                               {"FSC-W+", "Live immune cells (CD45+ PE-)"},
                                               {"SingletsW", "Live immune cells (CD45+ PE-)"},
                                               {"Live Single immune cells(FSC-H_FSC-W)",
                                                "Live immune cells (CD45+ PE-)"},
                                               {"nonDebris",
                                                "Live Single immune cells(FSC-H_FSC-W)"},
                                               {"Live Single PBMCs (SSC-A_FSC-A)", "nonDebris"},
                                               {"CD19+", "Live Single PBMCs (SSC-A_FSC-A)"},
                                               {"CD19-", "Live Single PBMCs (SSC-A_FSC-A)"},
                                               {"CD3+", "Live Single PBMCs (SSC-A_FSC-A)"},
                                               {"CD3-", "Live Single PBMCs (SSC-A_FSC-A)"},
                                               {"DC NK MONOCYTES (CD3- CD19-)",
                                                "Live Single PBMCs (SSC-A_FSC-A)"},
                                               {"MONOCYTES (CD14+)",
                                                "DC NK MONOCYTES (CD3- CD19-)"},
                                               {"Non classical monocytes (CD16+ CD14+)",
                                                "MONOCYTES (CD14+)"},
                                               {"Classical monocytes (CD16- CD14+)",
                                                "MONOCYTES (CD14+)"},
                                               {"CD14-", "DC NK MONOCYTES (CD3- CD19-)"},
                                               {"DC NK (CD20- CD14-)", "CD14-"},
                                               {"HLA-DR+", "DC NK (CD20- CD14-)"},
                                               {"DC (HLA-DR+)", "DC NK (CD20- CD14-)"},
                                               {"BB515-A-BV 711-A+", "DC (HLA-DR+)"},
                                               {"BB515-A+BV 711-A+", "DC (HLA-DR+)"},
                                               {"Myeloid DC (CD11c+ CD123-)", "DC (HLA-DR+)"},
                                               {"Plasmacytoid DC (CD11c- CD123+)", "DC (HLA-DR+)"},
                                               {"CD56+", "DC NK (CD20- CD14-)"},
                                               {"NK (CD16+)", "DC NK (CD20- CD14-)"},
                                               {"CD16-CD56-", "DC NK (CD20- CD14-)"},
                                               {"CD16+CD56-", "DC NK (CD20- CD14-)"},
                                               {"CD16-CD56+", "DC NK (CD20- CD14-)"},
                                               {"NK CD56LO", "DC NK (CD20- CD14-)"},
                                               {"NK CD56HI", "NK CD56LO"},};

  static final String[][] GATE_TREE_PANEL_1 = { // mapping[i][0] = child, mapping[i][1] = parent,
                                               // if length == 1, @ root
                                               {"boundary"}, {"nonDebris", "boundary"},
                                               {"Lymphocytes (SSC-A v FSC-A)", "nonDebris"},
                                               {"FSC-W+", "Lymphocytes (SSC-A v FSC-A)"},
                                               {"Single Cells (FSC-H v FSC-W)",
                                                "Lymphocytes (SSC-A v FSC-A)"},
                                               {"PE-A+", "Single Cells (FSC-H v FSC-W)"},
                                               {"Live cells (PE-)", "Single Cells (FSC-H v FSC-W)"},
                                               {"Tcells (CD3+ CD19-)", "Live cells (PE-)"},
                                               {"CD8-", "Tcells (CD3+ CD19-)"},
                                               {"CD4+", "Tcells (CD3+ CD19-)"},
                                               {"Helper Tcells-CD4+", "Tcells (CD3+ CD19-)"},
                                               {"CD45RA+", "Helper Tcells-CD4+"},
                                               {"CCR7gate", "CD45RA+"},
                                               {"effector helper Tcells (CCR7- CD45RA+)",
                                                "Helper Tcells-CD4+"},
                                               {"effector helper Tcells (CD95/CD28)",
                                                "effector helper Tcells (CCR7- CD45RA+)"},
                                               {"naive helper Tcells (CCR7+ CD45RA+)",
                                                "Helper Tcells-CD4+"},
                                               {"naive helper Tcells (CD95/CD28)",
                                                "naive helper Tcells (CCR7+ CD45RA+)"},
                                               {"CD45RA-", "Helper Tcells-CD4+"},
                                               {"CCR7gate", "CD45RA-"},
                                               {"effector memory helper Tcells (CCR7- CD45RA-)",
                                                "Helper Tcells-CD4+"},
                                               {"effector memory helper Tcells (CD95/CD28)",
                                                "effector memory helper Tcells (CCR7- CD45RA-)"},
                                               {"central memory helper Tcells (CCR7+ CD45RA-)",
                                                "Helper Tcells-CD4+"},
                                               {"central memory helper Tcells (CD95/CD28)",
                                                "central memory helper Tcells (CCR7+ CD45RA-)"},
                                               {"HLA-DR+", "Helper Tcells-CD4+"},
                                               {"activated helper Tcells (CD4+ HLA-DR+)",
                                                "Helper Tcells-CD4+"},
                                               {"naive helper Tcells (CD95- CD28+)",
                                                "Helper Tcells-CD4+"},
                                               {"effector memory helper Tcells (CD95- CD28-)",
                                                "Helper Tcells-CD4+"},
                                               {"central memory helper Tcells (CD95+ CD28+)",
                                                "Helper Tcells-CD4+"},
                                               {"CD4-", "Tcells (CD3+ CD19-)"}, {"CD8+", "CD4-"},
                                               {"cytotoxic Tcells-CD8+", "Tcells (CD3+ CD19-)"},
                                               {"CD45RA+", "cytotoxic Tcells-CD8+"},
                                               {"CCR7gate", "CD45RA+"},
                                               {"effector cytotoxic Tcells  (CCR7-  CD45RA+)",
                                                "cytotoxic Tcells-CD8+"},
                                               {"effector cytotoxic Tcells (CD95/CD28)",
                                                "effector cytotoxic Tcells  (CCR7-  CD45RA+)"},
                                               {"CD8EFFCD27Gate",
                                                "effector cytotoxic Tcells  (CCR7-  CD45RA+)"},
                                               {"CD8EFFCD28Gate",
                                                "effector cytotoxic Tcells  (CCR7-  CD45RA+)"},
                                               {"pE cytotoxic Tcells (CD27-  CD28-)",
                                                "effector cytotoxic Tcells  (CCR7-  CD45RA+)"},
                                               {"CD28+CD27-",
                                                "effector cytotoxic Tcells  (CCR7-  CD45RA+)"},
                                               {"pE2 cytotoxic Tcells (CD27+ , CD28-)",
                                                "effector cytotoxic Tcells  (CCR7-  CD45RA+)"},
                                               {"pE1 cytotoxic Tcells (CD27+  CD28+)",
                                                "effector cytotoxic Tcells  (CCR7-  CD45RA+)"},
                                               {"naive cytotoxic Tcells (CCR7+ , CD45RA+)",
                                                "cytotoxic Tcells-CD8+"},
                                               {"naive cytotoxic Tcells (CD95/CD28)",
                                                "naive cytotoxic Tcells (CCR7+ , CD45RA+)"},
                                               {"CD45RA-", "cytotoxic Tcells-CD8+"},
                                               {"CCR7gate", "CD45RA-"},
                                               {"effector memory cytotoxic Tcells (CCR7- , CD45RA-)",
                                                "cytotoxic Tcells-CD8+"},
                                               {"effector memory cytotoxic Tcells (CD95/CD28)",
                                                "effector memory cytotoxic Tcells (CCR7- , CD45RA-)"},
                                               {"CD8MEMCD27Gate",
                                                "effector memory cytotoxic Tcells (CCR7- , CD45RA-)"},
                                               {"CD8MEMCD28Gate",
                                                "effector memory cytotoxic Tcells (CCR7- , CD45RA-)"},
                                               {"EM3 cytotoxic Tcells (CD27-  CD28-)",
                                                "effector memory cytotoxic Tcells (CCR7- , CD45RA-)"},
                                               {"EM4 cytotoxic Tcells (CD27-  CD28+)",
                                                "effector memory cytotoxic Tcells (CCR7- , CD45RA-)"},
                                               {"EM2 cytotoxic Tcells (CD27+  CD28-)",
                                                "effector memory cytotoxic Tcells (CCR7- , CD45RA-)"},
                                               {"EM1 cytotoxic Tcells (CD27+  CD28+)",
                                                "effector memory cytotoxic Tcells (CCR7- , CD45RA-)"},
                                               {"central memory cytotoxic Tcells (CCR7+ , CD45RA-)",
                                                "cytotoxic Tcells-CD8+"},
                                               {"central memory cytotoxic Tcells (CD95/CD28)",
                                                "central memory cytotoxic Tcells (CCR7+ , CD45RA-)"},
                                               {"HLA-DR+", "cytotoxic Tcells-CD8+"},
                                               {"activated cytotoxic Tcells (CD8+ HLA-DR+)",
                                                "cytotoxic Tcells-CD8+"},
                                               {"naive cytotoxic Tcells (CD95- , CD28+)",
                                                "cytotoxic Tcells-CD8+"},
                                               {"central memory cytotoxic Tcells (CD95+ , CD28+)",
                                                "cytotoxic Tcells-CD8+"},
                                               {"effector memory cytotoxic Tcells (CD95- , CD28-)",
                                                "cytotoxic Tcells-CD8+"},
                                               {"pE cytotoxic Tcells (CD95-CD28-,CD27-  CD28-)",
                                                "effector memory cytotoxic Tcells (CD95- , CD28-)"},
                                               {"pE2 cytotoxic Tcells (CD95-CD28-,CD27+ , CD28-)",
                                                "effector memory cytotoxic Tcells (CD95- , CD28-)"},
                                               {"pE1 cytotoxic Tcells (CD95-CD28-,CD27+  CD28+)",
                                                "effector memory cytotoxic Tcells (CD95- , CD28-)"},
                                               {"CD3-", "Live cells (PE-)"},
                                               {"B cells (CD3- CD19+)", "CD3-"},
                                               {"IgDgate", "B cells (CD3- CD19+)"},
                                               {"CD27gate", "B cells (CD3- CD19+)"},
                                               {"IgD-CD27-", "B cells (CD3- CD19+)"},
                                               {"naive Bcells (CD27- IgD+)",
                                                "B cells (CD3- CD19+)"},
                                               {"IgD- memory Bcells (CD27+)",
                                                "B cells (CD3- CD19+)"},
                                               {"IgD+ memory Bcells (CD27+)",
                                                "B cells (CD3- CD19+)"},

  };

}
