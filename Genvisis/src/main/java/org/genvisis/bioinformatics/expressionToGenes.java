package org.genvisis.bioinformatics;

import java.io.*;
import java.util.*;

import org.genvisis.common.*;

public class expressionToGenes {
	public static final int[] GENE_COLS = {0, 1, 2, 3, 4, 5, 6, 7};
	public static final String[] HS_HEADER = {"NM_ID", "ID", "GENE ", "TITLE", "GENE_ID", "CHROMOSOME", "CYTOBAND", "EXPRESS"};
	public static final String[] GENES_HEADER = {"Name", "reference_GeneID", "reference_chr", "reference_start", "reference_stop", "reference_sense", "reference_placed", "reference_imputed", "Celera_GeneID", "Celera_chr", "Celera_start", "Celera_stop", "Celera_sense", "Celera_placed", "Celera_imputed"};
	public static final String[] TISSUE_TERMS_CHOSEN = {"brain"};
	public static final String[] TISSUE_TYPES = {"fetal brain", "whole brain", "temporal lobe", "parietal lobe", "occipital lobe", "prefrontal cortex", "cingulate cortex", "cerebellum", "cerebellum peduncles", "amygdala", "hypothalamus", "thalamus", "subthalamic nucleus", "caudate nucleus", "globus pallidus", "olfactory bulb", "pons", "medulla oblongata", "spinal cord", "ciliary ganglion", "trigeminal ganglion", "superior cervical ganglion", "dorsal root ganglion", "thymus", "tonsil", "lymph node", "bone marrow", "BM-CD71+ early erythroid", "BM-CD33+ myeloid", "BM-CD105+ endothelial", "BM-CD34+", "whole blood", "PB-BDCA4+ dentritic cells", "PB-CD14+ monocytes", "PB-CD56+ NKCells", "PB-CD4+ Tcells", "PB-CD8+ Tcells", "PB-CD19+ Bcells", "leukemia lymphoblastic(molt4)", "721 B lymphoblasts", "lymphoma Burkitts Raji", "leukemia promyelocytic(hl60)", "lymphoma Burkitts Daudi", "leukemia chronic myelogenous(k562)", "colorectal adenocarcinoma", "appendix", "skin", "adipocyte", "fetal thyroid", "thyroid", "pituitary gland", "adrenal gland", "adrenal cortex", "prostate", "salivary gland", "pancreas", "pancreatic islets", "atrioventricular node", "heart", "cardiac myocytes", "skeletal muscle", "tongue", "smooth muscle", "uterus", "uterus corpus", "trachea", "bronchial epithelial cells", "fetal lung", "lung", "kidney", "fetal liver", "liver", "placenta", "testis", "testis Leydig cell", "testis germ cell", "testis interstitial", "testis seminiferous tubule", "ovary"};
	public static final String[] TISSUES_CHOSEN = {"whole brain", "prefrontal cortex", "cingulate cortex", "cerebellum", "amygdala", "hypothalamus", "thalamus", "subthalamic nucleus", "caudate nucleus", "globus pallidus", "olfactory bulb"};
	public static final String DELIMITER = "\t";

	public static class HS_DataPoint {
		public String[] NM_IDs;
		public String HS_ID;
		public String GENE;
		public String title;
		public String geneID;
		public String chromosome;
		public String cytoband;
		public String[] express;

		public HS_DataPoint(String[] line) {
			this.NM_IDs = line[0].split(";");
			for (int i = 0; i<this.NM_IDs.length; i++) {
				this.NM_IDs[i] = this.NM_IDs[i].indexOf(".")>0?this.NM_IDs[i].substring(0, this.NM_IDs[i].indexOf(".")):this.NM_IDs[i];
			}
			this.HS_ID = line[1];
			this.GENE = line[2];
			this.title = line[3];
			this.geneID = line[4];
			this.chromosome = line[5];
			this.cytoband = line[6];
			line = line[7].split("\\|");
			this.express = new String[line.length];
			for (int i = 0; i<line.length; i++) {
				this.express[i] = line[i].trim();
			}
		}
	}

	public static void runExpressionToGenes(String genes, String express, String tag_prop, String prop_ref, String hs_data) {
		BufferedReader reader = null;
		PrintWriter writer = null;
		PrintWriter geneids = null;
		String[] line, header;
		String temp, trav;
		Hashtable<String,Vector<String>> tagPropLookup = new Hashtable<String,Vector<String>>();
		Hashtable<String,String> propRefLookup = new Hashtable<String,String>();
		Hashtable<String,String> geneData = new Hashtable<String,String>();
		Hashtable<String,HS_DataPoint> hsData = new Hashtable<String,HS_DataPoint>();
		Vector<String> v = new Vector<String>();
		HS_DataPoint hsdp;
		Vector<String> nms;

		System.out.println("Populating tagPropLookup...");
		try {
			reader = new BufferedReader(new FileReader(tag_prop));
			reader.readLine();
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				HashVec.addToHashVec(tagPropLookup, line[1], line[0], false);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+tag_prop+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+tag_prop+"\"");
			System.exit(2);
		}

		System.out.println("Populating propRefLookup...");
		try {
			reader = new BufferedReader(new FileReader(prop_ref));
			reader.readLine();
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				propRefLookup.put(line[0], line[1]);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+prop_ref+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+prop_ref+"\"");
			System.exit(2);
		}

		System.out.println("Populating Hs.data...");
		try {
			reader = new BufferedReader(new FileReader(hs_data));
			header = reader.readLine().split("\t");
			ext.checkHeader(header, HS_HEADER, true);
			while (reader.ready()) {
				line = reader.readLine().split("\t");
				if (!line[0].equals(".")) {
					hsdp = new HS_DataPoint(line);
					for (int i = 0; i<hsdp.NM_IDs.length; i++) {
						hsData.put(hsdp.NM_IDs[i], hsdp);
					}
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+hs_data+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+hs_data+"\"");
			System.exit(2);
		}

		System.out.println("Populating genes.xls data...");
		try {
			reader = new BufferedReader(new FileReader(genes));
			ext.checkHeader(line = reader.readLine().split("[\\s]+"), GENES_HEADER, true);
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				geneData.put(line[0], line[1]+"\t"+line[2]);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+genes+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+genes+"\"");
			System.exit(2);
		}

		System.out.println("Merging with expression data...");
		try {
			reader = new BufferedReader(new FileReader(express));
			reader.readLine();
			writer = new PrintWriter(new FileWriter("reverse.xls"));
			geneids = new PrintWriter(new FileWriter("lookup.xls"));
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				trav = line[0];
				writer.print(trav);
				nms = new Vector<String>();
				if (tagPropLookup.containsKey(trav)) {
					v = tagPropLookup.get(trav);
					for (int i = 0; i<v.size(); i++) {
						if (propRefLookup.containsKey(v.elementAt(i))) {
							temp = propRefLookup.get(v.elementAt(i));
							if (hsData.containsKey(temp)) {
								HashVec.addIfAbsent(temp, nms);
							} else {
								System.out.println(temp);
							}
						}
					}
					for (int i = 0; i<nms.size(); i++) {
						hsdp = (HS_DataPoint)hsData.get(nms.elementAt(i));
						writer.print("\t"+nms.elementAt(i)+"\t"+hsdp.title+"\t"+hsdp.chromosome+"\t"+hsdp.geneID);
						if (geneData.containsKey(hsdp.geneID)) {
							writer.println("\t"+geneData.get(hsdp.geneID));
							// geneids.println(trav+"\t"+hsdp.geneID);
						} else {
							writer.println("\tno gene.xls data");
						}
					}
					for (int i = 0; i<nms.size(); i++) {
						hsdp = (HS_DataPoint)hsData.get(nms.elementAt(i));
						writer.print("\t"+nms.elementAt(i)+"\t"+hsdp.title+"\t"+hsdp.chromosome+"\t"+hsdp.geneID);
						if (geneData.containsKey(hsdp.geneID)) {
							writer.println("\t"+geneData.get(hsdp.geneID));
							geneids.println(trav+"\t"+hsdp.geneID);
						} else {
							writer.println("\tno gene.xls data");
						}
					}
					if (nms.size()==0) {
						writer.println("\tno refSeq matchup");
					}
				} else {
					writer.println("\tno prop gene name");
				}
			}
			writer.close();
			geneids.close();
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+express+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+express+"\"");
			System.exit(2);
		}

	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		String genes = "genes.xls";
		String express = "expression_levels.prn";
		// String tissueTypes = "79_tissue_names.prn";
		String tag_prop = "knownToAll.prn";
		String prop_ref = "propgenenames_refseqgenenames.prn";
		String hs_data = "Hs.data.xls";

		String usage = "\n"+"park.genesByTissue requires 0-5 arguments\n"+"   (1) filename of genes file (i.e. genes="+genes+" (default))\n"+"   (2) filename of expression data (i.e. express="+express+" (default))\n"+
		// " (3) tissue_types lookup file (i.e. tissueTypes="+tag_prop+"
		// (default))\n"+
		"   (4) tag_prop lookup file (i.e. tagProp="+tag_prop+" (default))\n"+"   (5) prop_ref lookup file (i.e. propRef="+prop_ref+" (default))\n"+"   (6) Hs.data file (i.e. hsData="+hs_data+" (default))\n"+"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("genes=")) {
				genes = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("express=")) {
				express = args[i].split("=")[1];
				numArgs--;
				// } else if (args[i].startsWith("tissueTypes=")) {
				// tissueTypes = args[i].split("=")[1];
				// numArgs--;
			} else if (args[i].startsWith("tagProp=")) {
				tag_prop = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("propRef=")) {
				prop_ref = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("hsData=")) {
				hs_data = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			runExpressionToGenes(genes, express, tag_prop, prop_ref, hs_data);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
