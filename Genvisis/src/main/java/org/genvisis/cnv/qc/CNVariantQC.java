package org.genvisis.cnv.qc;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;

import org.genvisis.cnv.filesys.CNVQC;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Logger;
import org.genvisis.filesys.CNVariant;

// TODO extends?
public class CNVariantQC implements Serializable {
	public static final long serialVersionUID = 1L;
	public static final String[] PLINK_CNV_QC_HEADER = {"FID", "IID", "CHR", "BP1", "BP2", "TYPE",
																											"SCORE", "SITES", "HEIGHT", "BAFQC",
																											"NUMHETSBAF", "NUMHETSGENOTYPE", "NUMPOLYMAF",
																											"SUM2pq", "SOURCEFILE"};
	private CNVariant cnVariant;
	private double height;
	private double bafQC;
	private int numHetsBAF;
	private int numHetsGenotype;
	private int numPolyMAF;
	private Hashtable<String, Integer> markersIncnVariant;
	private String[] markerList;
	private double[] mafs;
	private double[] bafs;
	private double[] lrrs;
	private byte[] genotypes;
	private boolean[] polyMarkers;
	private double sum2pq;
	private String sourceFile;
	private double cnvLRRStdev;
	private double sampleCallRate;

	public CNVariantQC(CNVariant cnVariant) {
		this.cnVariant = cnVariant;
		height = Double.NaN;
		bafQC = Double.NaN;
		cnvLRRStdev = Double.NaN;
		sum2pq = Double.NaN;
		sampleCallRate = Double.NaN;
		numHetsBAF = 0;
		numHetsGenotype = 0;
		numPolyMAF = 0;
		markersIncnVariant = new Hashtable<String, Integer>();
		markerList = new String[cnVariant.getNumMarkers()];
		mafs = new double[cnVariant.getNumMarkers()];
		bafs = new double[cnVariant.getNumMarkers()];
		lrrs = new double[cnVariant.getNumMarkers()];
		genotypes = new byte[cnVariant.getNumMarkers()];
		polyMarkers = new boolean[cnVariant.getNumMarkers()];
		sourceFile = "NA";
	}

	public CNVariantQC(String[] plinkPlusQCFormatLine) {
		cnVariant = new CNVariant(ArrayUtils.subArray(plinkPlusQCFormatLine, 0, 7));
		height = Double.parseDouble(plinkPlusQCFormatLine[8]);
		bafQC = Double.parseDouble(plinkPlusQCFormatLine[9]);
		numHetsBAF = Integer.parseInt(plinkPlusQCFormatLine[10]);
		numHetsGenotype = Integer.parseInt(plinkPlusQCFormatLine[11]);
		numPolyMAF = Integer.parseInt(plinkPlusQCFormatLine[12]);
		sum2pq = Double.parseDouble(plinkPlusQCFormatLine[13]);
	}


	public void findMarkerNamesinCNV(Project proj, int[][] indices, int[] positions,
																	 String[] markerNames) {
		int numMarkers = 0;
		Logger log = proj.getLog();

		for (int i = 0; i < indices[cnVariant.getChr()].length; i++) {
			int position = positions[indices[cnVariant.getChr()][i]];
			if (inCNV(position, cnVariant)) {
				// will exit if found too many markers
				checkNumMarkers(numMarkers, false, log);
				markersIncnVariant.put(markerNames[indices[cnVariant.getChr()][i]],
															 indices[cnVariant.getChr()][i]);
				markerList[numMarkers] = markerNames[indices[cnVariant.getChr()][i]];
				numMarkers++;
			}
		}
		// will exit if markers found do not match number in cnv file
		checkNumMarkers(numMarkers, true, log);
	}

	public void assignMAFs(Hashtable<String, Double> markerBAFs, Logger log) {
		for (int i = 0; i < markerList.length; i++) {
			if (!markerBAFs.containsKey(markerList[i])) {
				log.reportError("Error - did not find marker " + markerList[i] + " in markerFreq file");
				System.exit(1);
			} else {
				mafs[i] = markerBAFs.get(markerList[i]);
			}
		}
	}

	public static String[] getIDList(CNVariantQC[] cnVariantQCs,
																	 Hashtable<String, Hashtable<String, Integer>> defineCompHash) {
		ArrayList<String> al = new ArrayList<String>();
		Hashtable<String, Integer> tracker = new Hashtable<String, Integer>();
		for (CNVariantQC cnVariantQC : cnVariantQCs) {
			String FID = cnVariantQC.getCnVariant().getFamilyID();
			// TODO
			// change to ind id when i update input file
			if (!tracker.containsKey(FID)
					&& (defineCompHash == null || defineCompHash.containsKey(FID))) {
				tracker.put(FID, 1);
				al.add(FID);
			}
		}
		return toStringArray(al);
	}

	public static Hashtable<String, CNVariantQC[]> getIndCNVQCs(String[] inds,
																															CNVariantQC[] cnVariantQCs) {
		Hashtable<String, ArrayList<CNVariantQC>> IndCNVQCs = new Hashtable<String, ArrayList<CNVariantQC>>();
		for (int i = 0; i < cnVariantQCs.length; i++) {
			if (!IndCNVQCs.containsKey(cnVariantQCs[i].getCnVariant().getFamilyID())) {
				// TODO
				// change to ind id when i update input file
				IndCNVQCs.put(cnVariantQCs[i].getCnVariant().getFamilyID(), new ArrayList<CNVariantQC>());
			}
			IndCNVQCs.get(cnVariantQCs[i].getCnVariant().getFamilyID()).add(cnVariantQCs[i]);
		}
		return getIndividualCNVQCArrays(inds, IndCNVQCs);
	}

	public static CNVariantQC[][][] prepCNVQCsForComparison(Project proj, String plinkCnvQCs,
																													Hashtable<String, Hashtable<String, Integer>> defineCompHash) {
		CNVariantQC[] cnVariantQCs = CNVQC.load(proj.PROJECT_DIRECTORY.getValue() + plinkCnvQCs)
																			.getCnVariantQCs();
		String[] inds = getIDList(cnVariantQCs, defineCompHash);
		if (inds.length < 2) {
			proj.getLog()
					.reportError("Error - the cnvQC file " + proj.PROJECT_DIRECTORY.getValue() + plinkCnvQCs
											 + " does not contain any matched IDs found in duplicates file");
			System.exit(1);
		}
		Hashtable<String, CNVariantQC[]> indCNVQCssArrays = getIndCNVQCs(inds, cnVariantQCs);
		return assignCNVComparisions(inds, defineCompHash, indCNVQCssArrays, proj.getLog());

	}

	public static void filterCNVQCsByComparison(Project proj, String plinkCnvQCs,
																							Hashtable<String, Hashtable<String, Integer>> defineCompHash) {
		CNVariantQC[] cnVariantQCs = CNVQC.load(proj.PROJECT_DIRECTORY.getValue() + plinkCnvQCs)
																			.getCnVariantQCs();
		String[] inds = getIDList(cnVariantQCs, defineCompHash);
		if (inds.length < 2) {
			proj.getLog()
					.reportError("Error - the cnvQC file " + proj.PROJECT_DIRECTORY.getValue() + plinkCnvQCs
											 + " does not contain any matched IDs found in duplicates file");
			System.exit(1);
		}
		Hashtable<String, CNVariantQC[]> indCNVQCssArrays = getIndCNVQCs(inds, cnVariantQCs);
		ArrayList<CNVariantQC> filteredByComparision = new ArrayList<CNVariantQC>();
		for (String ind : inds) {
			if (defineCompHash.containsKey(ind)) {
				CNVariantQC[] indcnVariantQCs = indCNVQCssArrays.get(ind);
				for (CNVariantQC indcnVariantQC : indcnVariantQCs) {
					filteredByComparision.add(indcnVariantQC);
				}
			}
		}
		new CNVQC(filteredByComparision.toArray(new CNVariantQC[filteredByComparision.size()])).serialize(proj.PROJECT_DIRECTORY.getValue()
																																																			+ plinkCnvQCs.replaceAll(".ser",
																																																															 ".comp.ser"));
	}

	public static CNVariantQC[] getCNVariantQCFromPlinkFile(Project proj, String plinkCnvs) {
		CNVariant[] fileCNVs = CNVariant.loadPlinkFile(proj.PROJECT_DIRECTORY.getValue() + plinkCnvs);
		CNVariantQC[] cnVariantQCs = new CNVariantQC[fileCNVs.length];
		for (int i = 0; i < fileCNVs.length; i++) {
			cnVariantQCs[i] = new CNVariantQC(fileCNVs[i]);
		}
		return cnVariantQCs;
	}

	public void assigncnvLRRStdev() {
		cnvLRRStdev = ArrayUtils.median(lrrs);
	}

	public double getSampleCallRate() {
		return sampleCallRate;
	}

	public void setSampleCallRate(double sampleCallRate) {
		this.sampleCallRate = sampleCallRate;
	}

	public double getCnvLRRStdev() {
		return cnvLRRStdev;
	}

	public void setCnvLRRStdev(double cnvLRRStdev) {
		this.cnvLRRStdev = cnvLRRStdev;
	}

	public double[] getLrrs() {
		return lrrs;
	}

	public void setLrrs(double[] lrrs) {
		this.lrrs = lrrs;
	}

	public byte[] getGenotypes() {
		return genotypes;
	}

	public void setGenotypes(byte[] genotypes) {
		this.genotypes = genotypes;
	}

	public double[] getMafs() {
		return mafs;
	}

	public double[] getBafs() {
		return bafs;
	}

	public void setBafs(double[] bafs) {
		this.bafs = bafs;
	}

	public boolean[] getPolyMarkers() {
		return polyMarkers;
	}

	public String toPlinkPlusQCFormat() {
		return cnVariant.toPlinkFormat() + "\t" + height + "\t" + bafQC + "\t" + numHetsBAF + "\t"
					 + numHetsGenotype + "\t" + numPolyMAF + "\t" + sum2pq + "\t" + sourceFile;
	}

	public int getNumPolyMAF() {
		return numPolyMAF;
	}

	public void setNumPolyMAF(int numPolyMAF) {
		this.numPolyMAF = numPolyMAF;
	}

	public void setSourceFile(String sourceFile) {
		this.sourceFile = sourceFile;
	}

	public double getMaf2pq() {
		return sum2pq;
	}

	public void setMaf2pq(double maf2pq) {
		sum2pq = maf2pq;
	}

	public Hashtable<String, Integer> getMarkersIncnVariant() {
		return markersIncnVariant;
	}

	public String[] getMarkerList() {
		return markerList;
	}

	public CNVariant getCnVariant() {
		return cnVariant;
	}

	public int getNumHetsBAF() {
		return numHetsBAF;
	}

	public void setNumHetsBAF(int numHetsBAF) {
		this.numHetsBAF = numHetsBAF;
	}

	public int getNumHetsGenotype() {
		return numHetsGenotype;
	}

	public void setNumHetsGenotype(int numHetsGenotype) {
		this.numHetsGenotype = numHetsGenotype;
	}

	public void setCnVariant(CNVariant cnVariant) {
		this.cnVariant = cnVariant;
	}

	public double getHeight() {
		return height;
	}

	public void setHeight(double height) {
		this.height = height;
	}

	public double getBafQC() {
		return bafQC;
	}

	public void setBafQC(double bafQC) {
		this.bafQC = bafQC;
	}

	private void checkNumMarkers(int numMarkers, boolean last, Logger log) {
		if (numMarkers < cnVariant.getNumMarkers() && last) {
			log.reportError(numMarkers + "Error - mismatched number of markers contained in cnv "
											+ cnVariant.toPlinkFormat() + " and markers in the position file ");
			System.exit(1);
		} else if (numMarkers == cnVariant.getNumMarkers() && !last) {
			log.reportError("Error - mismatched number of markers contained in cnv "
											+ cnVariant.toPlinkFormat() + " and markers in the position file ");
			System.exit(1);
		}
	}

	private boolean inCNV(int position, CNVariant cnv) {
		return position >= cnv.getStart() && position <= cnv.getStop();
	}

	private static String[] toStringArray(List<String> stringList) {
		return stringList.toArray(new String[stringList.size()]);
	}

	private static Hashtable<String, CNVariantQC[]> getIndividualCNVQCArrays(String[] inds,
																																					 Hashtable<String, ArrayList<CNVariantQC>> allIndCNVQCs) {
		Hashtable<String, CNVariantQC[]> allIndCNVQCsArray = new Hashtable<String, CNVariantQC[]>();
		for (String ind : inds) {
			allIndCNVQCsArray.put(ind, toCNVQCArray(allIndCNVQCs.get(ind)));
		}
		return allIndCNVQCsArray;
	}

	public static CNVariantQC[] toCNVQCArray(List<CNVariantQC> cnvQCs) {
		return cnvQCs.toArray(new CNVariantQC[cnvQCs.size()]);
	}

	private static CNVariantQC[][][] assignCNVComparisions(String[] inds,
																												 Hashtable<String, Hashtable<String, Integer>> defineCompHash,
																												 Hashtable<String, CNVariantQC[]> indCNVQCssArrays,
																												 Logger log) {
		Hashtable<String, Boolean> compared = new Hashtable<String, Boolean>();
		ArrayList<CNVariantQC[]> toCompare1 = new ArrayList<CNVariantQC[]>();
		ArrayList<CNVariantQC[]> toCompare2 = new ArrayList<CNVariantQC[]>();
		for (int j = 0; j < inds.length; j++) {
			for (int k = 0; k < inds.length; k++) {
				if (doComparision(defineCompHash, compared, inds[j], inds[k], indCNVQCssArrays) && j != k) {
					if (indCNVQCssArrays.get(inds[j]).length > 0
							&& indCNVQCssArrays.get(inds[k]).length > 0) {
						toCompare1.add(indCNVQCssArrays.get(inds[j]));
						toCompare2.add(indCNVQCssArrays.get(inds[k]));
					} else {
						log.reportError("Error -  CNVs were not found for comparison between " + inds[j]
														+ " and " + inds[k]);
						System.exit(1);
					}
				}
			}
		}
		return assignedToArray(toCompare1, toCompare2);
	}

	private static boolean doComparision(Hashtable<String, Hashtable<String, Integer>> defineCompHash,
																			 Hashtable<String, Boolean> compared, String ind1,
																			 String ind2,
																			 Hashtable<String, CNVariantQC[]> allIndCNVsArrays) {
		boolean doIt = false;
		if (defineCompHash.containsKey(ind1) && defineCompHash.get(ind1).containsKey(ind2)) {
			if (allIndCNVsArrays.containsKey(ind1) && allIndCNVsArrays.containsKey(ind2)) {
				if (!isCompared(compared, ind1, ind2)) {
					doIt = true;
				}
			}
		}
		return doIt;
	}

	private static boolean isCompared(Hashtable<String, Boolean> compared, String ind1, String ind2) {
		if (compared.containsKey(ind1 + "\t" + ind2) || compared.containsKey(ind2 + "\t" + ind1)) {
			return true;
		} else {
			compared.put(ind1 + "\t" + ind2, true);
			compared.put(ind2 + "\t" + ind1, true);
			return false;
		}
	}

	private static CNVariantQC[][][] assignedToArray(List<CNVariantQC[]> toCompare1,
																									 List<CNVariantQC[]> toCompare2) {
		CNVariantQC[][][] cnvQCsAssigned = new CNVariantQC[2][toCompare2.size()][];
		if (toCompare1.size() == toCompare2.size()) {
			for (int i = 0; i < toCompare2.size(); i++) {
				cnvQCsAssigned[0][i] = toCompare1.get(i);
				cnvQCsAssigned[1][i] = toCompare2.get(i);
			}
		} else {
			System.err.println("Warning -  Unmatched number of individuals found for comparisions");
			System.exit(1);
		}
		return cnvQCsAssigned;
	}

	public String getSourceFile() {
		return sourceFile;
	}
}

// private static CNVariantQC[] loadCNVariantQCFile(String plinkCnvQCs, Hashtable<String,
// Hashtable<String, Integer>> defineCompHash) {
// ArrayList<CNVariantQC> v = null;
// String[] line;
//
// v = new ArrayList<CNVariantQC>();
// try {
// BufferedReader reader = new BufferedReader(new FileReader(plinkCnvQCs));
//
// reader.mark(1000);
// line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
// if (!line[2].toLowerCase().equals("chr") && Positions.chromosomeNumber(line[2]) != -1) {
// reader.reset();
// }
// while (reader.ready()) {
// line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
// if (defineCompHash == null || defineCompHash.containsKey(line[0])) {
// v.add(new CNVariantQC(line));
// }
// }
// reader.close();
//
// return toCNVQCArray(v);
// } catch (FileNotFoundException fnfe) {
// System.err.println("Error: file \"" + plinkCnvQCs + "\" not found in current directory");
// fnfe.printStackTrace();
// } catch (IOException ioe) {
// System.err.println("Error reading file \"" + plinkCnvQCs + "\"");
// ioe.printStackTrace();
// }
//
// return null;
// }
