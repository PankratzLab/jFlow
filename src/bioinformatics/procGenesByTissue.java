// -Xms1024M -Xmx1024M
package bioinformatics;

import java.io.*;
import java.util.*;

import common.*;

public class procGenesByTissue {
	public static final int[] GENE_COLS = {0, 1, 2, 3, 4, 5, 6, 7};

	public static final String[] HS_HEADER = {"NM_ID", "ID", "GENE ", "TITLE", "GENE_ID", "CHROMOSOME", "CYTOBAND", "EXPRESS"};

	public static final String[] GENES_HEADER = {"GeneID", "reference_name", "reference_chr", "reference_start", "reference_stop", "reference_sense", "reference_placed", "reference_imputed", "Celera_name", "Celera_chr", "Celera_start", "Celera_stop", "Celera_sense", "Celera_placed", "Celera_imputed"};

	public static final String[] TISSUE_TERMS_CHOSEN = {"brain"};

	public static final String[] TISSUE_TYPES = {"fetal brain", "whole brain", "temporal lobe", "parietal lobe", "occipital lobe", "prefrontal cortex", "cingulate cortex", "cerebellum", "cerebellum peduncles", "amygdala", "hypothalamus", "thalamus", "subthalamic nucleus", "caudate nucleus", "globus pallidus", "olfactory bulb", "pons", "medulla oblongata", "spinal cord", "ciliary ganglion", "trigeminal ganglion", "superior cervical ganglion", "dorsal root ganglion", "thymus", "tonsil", "lymph node", "bone marrow", "BM-CD71+ early erythroid", "BM-CD33+ myeloid", "BM-CD105+ endothelial", "BM-CD34+", "whole blood", "PB-BDCA4+ dentritic cells", "PB-CD14+ monocytes", "PB-CD56+ NKCells", "PB-CD4+ Tcells", "PB-CD8+ Tcells", "PB-CD19+ Bcells", "leukemia lymphoblastic(molt4)", "721 B lymphoblasts", "lymphoma Burkitts Raji", "leukemia promyelocytic(hl60)", "lymphoma Burkitts Daudi", "leukemia chronic myelogenous(k562)", "colorectal adenocarcinoma", "appendix", "skin", "adipocyte", "fetal thyroid", "thyroid", "pituitary gland", "adrenal gland", "adrenal cortex", "prostate", "salivary gland", "pancreas", "pancreatic islets", "atrioventricular node", "heart", "cardiac myocytes", "skeletal muscle", "tongue", "smooth muscle", "uterus", "uterus corpus", "trachea", "bronchial epithelial cells", "fetal lung", "lung", "kidney", "fetal liver", "liver", "placenta", "testis", "testis Leydig cell", "testis germ cell", "testis interstitial", "testis seminiferous tubule", "ovary"};

	// public static final String[] TISSUES_CHOSEN = {"whole brain", "prefrontal
	// cortex", "cingulate cortex", "cerebellum", "amygdala", "hypothalamus",
	// "thalamus", "subthalamic nucleus", "caudate nucleus", "globus pallidus",
	// "olfactory bulb"};
	public static final String[] TISSUES_CHOSEN = {"whole brain"};

	public static final String DELIMITER = "\t";

	// public static final String CHROMOSOME = "13";
	// public static final String START = "62710392";
	// public static final String STOP = "101179302";

	// public static final String CHROMOSOME = "10";
	// public static final String START = "3050550";
	// public static final String STOP = "7202570";

	// public static final String CHROMOSOME = "12";
	// public static final String START = "10891047";
	// public static final String STOP = "42011271";

	// public static final String CHROMOSOME = "1";
	// public static final String START = "57178309";
	// public static final String STOP = "62293893";

	// public static final String CHROMOSOME = "1";
	// public static final String START = "200700000";
	// public static final String STOP = "200800000";

	// public static final String CHROMOSOME = "5";
	// public static final String START = "3900000";
	// public static final String STOP = "4050000";

	// public static final String CHROMOSOME = "13";
	// public static final String START = "18000000";
	// public static final String STOP = "20000000";

	// public static final String CHROMOSOME = "2"; // 2q region
	// public static final String START = "227241128";
	// public static final String STOP = "235934686";

	public static final String CHROMOSOME = "20";

	public static final String START = "14875484";

	public static final String STOP = "15275484";

	// public static final String OMIMLINK =
	// "http://www.ncbi.nlm.nih.gov/entrez/dispomim.cgi?id=";
	public static final String OMIMLINK = "http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&db=OMIM&dopt=Detailed&tmpl=dispomimTemplate&list_uids=";

	public static class HS_DataPoint {
		public String GENE;

		public String title;

		public String geneID;

		public String chromosome;

		public String cytoband;

		public String[] express;

		public HS_DataPoint(String[] line) {
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

	public static void runProcGenesByTissue(String genes, String express, String tagGeneID, String hs_data, String omimGeneID, String omimNames) {
		BufferedReader reader = null;
		PrintWriter writer = null;
		String[] line, header;
		String trav;
		Hashtable<String,HS_DataPoint> hash = new Hashtable<String,HS_DataPoint>();
		Hashtable<String,Vector<String>> tagGeneIDLookup, geneidOmimLookup;
		Hashtable<String,double[]> expressionData = new Hashtable<String,double[]>();
		Hashtable<String,String> omimNameLookup = new Hashtable<String,String>();
		Vector<String> v;
		int count;
		HS_DataPoint hsdp;
		int[] tissueIndices;
		double[] levels;

		tissueIndices = new int[TISSUES_CHOSEN.length];
		for (int i = 0; i<tissueIndices.length; i++) {
			tissueIndices[i] = ext.indexOfStr(TISSUES_CHOSEN[i], TISSUE_TYPES);
			if (tissueIndices[i]==-1) {
				System.err.println("Error - '"+TISSUES_CHOSEN[i]+"' is not a valid tissue type");
				System.exit(1);
			}
		}

		System.out.println("Populating tagToGeneID lookup table...");
		tagGeneIDLookup = HashVec.loadFileToHashVec(tagGeneID, 1, new int[] {0}, null, true, true);

		System.out.println("Populating expressionData...");
		try {
			reader = new BufferedReader(new FileReader(express));
			reader.readLine();
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				trav = line[0];
				levels = new double[tissueIndices.length];
				line = line[2].split(",");
				for (int i = 0; i<tissueIndices.length; i++) {
					levels[i] = Double.parseDouble(line[tissueIndices[i]]);
				}
				expressionData.put(trav, levels);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+express+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+express+"\"");
			System.exit(2);
		}

		System.out.println("Populating Hs.data...");
		try {
			reader = new BufferedReader(new FileReader(hs_data));
			header = reader.readLine().split("\t");
			ext.checkHeader(header, HS_HEADER, true);
			while (reader.ready()) {
				line = reader.readLine().split("\t");
				if (!line[4].equals(".")) {
					hash.put(line[4], new HS_DataPoint(line));
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

		System.out.println("Populating OMIM GeneID Lookup...");
		geneidOmimLookup = HashVec.loadFileToHashVec(omimGeneID, 0, new int[] {1}, null, false, true);

		System.out.println("Populating OMIM descriptions Lookup...");
		omimNameLookup = HashVec.loadFileToHashString(omimNames, false);

		System.out.println("Merging data with gene positions...");
		try {
			reader = new BufferedReader(new FileReader(genes));
			ext.checkHeader(line = reader.readLine().split("[\\s]+"), GENES_HEADER, true);
			writer = new PrintWriter(new FileWriter("expression.xls"));
			writer.print(line[0]+DELIMITER+line[1]+DELIMITER+line[2]+DELIMITER+line[3]+DELIMITER+line[4]+DELIMITER+"CytoBand"+DELIMITER+"Title");
			for (int i = 0; i<TISSUE_TERMS_CHOSEN.length; i++) {
				writer.print(DELIMITER+TISSUE_TERMS_CHOSEN[i]);
			}
			writer.print(DELIMITER+"ExpressionTagName");
			for (int i = 0; i<TISSUES_CHOSEN.length; i++) {
				writer.print(DELIMITER+TISSUES_CHOSEN[i]);
			}
			writer.print(DELIMITER+"OMIM entries");
			writer.println();
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				if (line[2].equals(CHROMOSOME)&&(Integer.parseInt(line[3])>Integer.parseInt(START)||Integer.parseInt(line[4])>Integer.parseInt(START))&&(Integer.parseInt(line[3])<Integer.parseInt(STOP)||Integer.parseInt(line[4])<Integer.parseInt(STOP))) {
					writer.print(line[0]+DELIMITER+line[1]+DELIMITER+line[2]+DELIMITER+line[3]+DELIMITER+line[4]);
					if (hash.containsKey(line[0])) {
						hsdp = (HS_DataPoint)hash.get(line[0]);
						writer.print(DELIMITER+hsdp.cytoband+DELIMITER+hsdp.title);
						for (int i = 0; i<TISSUE_TERMS_CHOSEN.length; i++) {
							writer.print(DELIMITER+(ext.indexOfStr(TISSUE_TERMS_CHOSEN[i], hsdp.express, false, true)>=0?1:0));
						}
					} else {
						writer.print(DELIMITER+Array.toStr(Array.stringArray(2+TISSUE_TERMS_CHOSEN.length, ".."), DELIMITER));
					}
					count = 0;
					if (tagGeneIDLookup.containsKey(line[0])) {
						v = tagGeneIDLookup.get(line[0]);
						for (int i = 0; i<v.size(); i++) {
							trav = v.elementAt(i);
							if (expressionData.containsKey(trav)) {
								levels = expressionData.get(trav);
								if (count>0) {
									writer.print("\n"+Array.toStr(Array.stringArray(5+2+TISSUE_TERMS_CHOSEN.length), DELIMITER));
								}
								writer.print(DELIMITER+trav);
								for (int j = 0; j<levels.length; j++) {
									writer.print(DELIMITER+levels[j]);
								}
								count++;
							} else {
								System.err.println("My don't we feel crunchy");
							}
						}
					} else {
						writer.print(DELIMITER+Array.toStr(Array.stringArray(TISSUES_CHOSEN.length+1, "."), DELIMITER));
					}
					if (geneidOmimLookup.containsKey(line[0])) {
						v = geneidOmimLookup.get(line[0]);
						for (int i = 0; i<v.size(); i++) {
							if (omimNameLookup.containsKey(v.elementAt(i))) {
								writer.print(DELIMITER+"=HYPERLINK(\""+OMIMLINK+v.elementAt(i)+"\", \""+omimNameLookup.get(v.elementAt(i))+"\")");
							} else {
								writer.print(DELIMITER+"=HYPERLINK(\""+OMIMLINK+v.elementAt(i)+"\", \""+v.elementAt(i)+"\")");
							}
						}
					} else {
						// writer.print(DELIMITER+"no OMIM entry");
					}

					writer.println(DELIMITER+" ");

				}
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+genes+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+genes+"\"");
			System.exit(2);
		}

	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		String genes = "genes.xls";
		String express = "expression_levels.prn";
		String tagGeneID = "tagsToGeneIDs.prn";
		String hs_data = "Hs.data.xls";
		String omimNames = "omim_names.dat";
		String omimGeneID = "geneid_omim.dat";

		String usage = "\n"+"park.genesByTissue requires 0-4 arguments\n"+"   (1) filename of genes file (i.e. genes="+genes+" (default))\n"+"   (2) filename of expression data (i.e. express="+express+" (default))\n"+"   (3) Tag_to_GeneID_table lookup file (i.e. tagGeneID="+tagGeneID+" (default))\n"+"   (4) Hs.data file (i.e. hsData="+hs_data+" (default))\n"+"";

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
			} else if (args[i].startsWith("tagGeneID=")) {
				tagGeneID = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("hsData=")) {
				hs_data = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("geneidOmim=")) {
				omimGeneID = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("omimNames=")) {
				omimNames = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			runProcGenesByTissue(genes, express, tagGeneID, hs_data, omimGeneID, omimNames);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
