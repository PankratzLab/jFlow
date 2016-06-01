package cnv.qc;

import java.util.HashMap;
import java.util.HashSet;

import cnv.filesys.ClusterFilterCollection;
import cnv.filesys.Project;
import cnv.filesys.Sample;
import common.Array;
import common.Files;
import common.HashVec;
import common.Logger;
import common.ext;

public class DuplicateConcordance {
	
	private double projectConcordance;
	private double[] markerConcordance;
	
	private DuplicateConcordance(double projectConcordance, double[] markerConcordance) {
		this.projectConcordance = projectConcordance;
		this.markerConcordance = markerConcordance;
	}
	
	
	
	public double getProjectConcordance() {
		return projectConcordance;
	}

	public double[] getMarkerConcordance() {
		return markerConcordance;
	}

	public static DuplicateConcordance calculateMarkerConcordances(Project proj, String targetMarkersFile) {
		Logger log = proj.getLog();
		ClusterFilterCollection clusterFilterCollection = proj.getClusterFilterCollection();
		String[] markerNames;
		int[] markerIndices;
		if (targetMarkersFile == null) {
			markerNames = proj.getMarkerNames();
			markerIndices = null;
		} else {
			markerNames = proj.getTargetMarkers(targetMarkersFile);
			markerIndices = ext.indexFactors(markerNames, proj.getMarkerNames(), true, false);
		}
		String sampleData = proj.SAMPLE_DATA_FILENAME.getValue();
		if (sampleData == null) {
			log.reportTimeError("Project Sample Data file is not defined, cannot determine duplicates");
			return null;
		}
		if (!Files.exists(sampleData)) {
			log.reportTimeError("Project Sample Data file, " + sampleData + " could not be found, cannot determine duplicates");
			return null;
		}
		String[] sampleDataHeader = Files.getHeaderOfFile(sampleData, log);
		String[] sampleDataCols = new String[] {"DNA", "CLASS=Exclude", "DuplicateId"};
		int[] sampleDataIndices = ext.indexFactors(sampleDataCols, sampleDataHeader, false, log, false, false);
		for (int i = 0; i < sampleDataIndices.length; i++) {
			if (sampleDataIndices[i] == -1) {
				log.reportTimeError("Could not find " + sampleDataCols[i] + " in Sample Data file, cannot determine duplicates");
				return null;
			}
		}
		String[][] sampleInfo = HashVec.loadFileToStringMatrix(sampleData, true, sampleDataIndices, proj.JAR_STATUS.getValue());
		HashVec.loadFileToStringArray(sampleData, true, sampleDataIndices, false);
		HashMap<String, HashSet<String>> duplicateSets = new HashMap<String, HashSet<String>>();
		for (String[] sampleLine : sampleInfo) {
			String dna = sampleLine[0];
			boolean exclude = sampleLine[1].equals("1");
			String duplicateID = sampleLine[2];
			
			if (!duplicateID.equals(".") && !exclude) {
				HashSet<String> duplicateSet = duplicateSets.get(duplicateID);
				if (duplicateSet == null) {
					duplicateSet = new HashSet<String>();
					duplicateSets.put(duplicateID, duplicateSet);
				}
				duplicateSet.add(dna);
			}
		}
		
		
		int[] markerHits = Array.intArray(markerNames.length, 0);
		int totalHits = 0;
		int pairsChecked = 0;
		
		for (HashSet<String> duplicateSet : duplicateSets.values()) {
			HashSet<String> loopDuplicateSet = new HashSet<String>(duplicateSet);
			for (String dna1 : loopDuplicateSet) {
				duplicateSet.remove(dna1);
				if (!duplicateSet.isEmpty()) {
					Sample sample1 = proj.getFullSampleFromRandomAccessFile(dna1);
					if (sample1 == null) {
						log.reportTimeError("Could not find data for Sample " + dna1 + ", will not be used to calculate concordance");
						continue;
					}
					for (String dna2 : duplicateSet) {
						Sample sample2 = proj.getFullSampleFromRandomAccessFile(dna2);
						if (sample2 == null) {
							log.reportTimeError("Could not find data for Sample " + dna2 + ", will not be used to calculate concordance");
							continue;
						}
						pairsChecked++;
						byte[] s1Genotypes, s2Genotypes;
						if (clusterFilterCollection == null) {
							s1Genotypes = sample1.getAB_Genotypes();
							s2Genotypes = sample2.getAB_Genotypes();
						} else {
							s1Genotypes = sample1.getAB_GenotypesAfterFilters(markerNames, markerIndices, clusterFilterCollection, 0.0f);
							s2Genotypes = sample2.getAB_GenotypesAfterFilters(markerNames, markerIndices, clusterFilterCollection, 0.0f);
						}
						for (int j = 0; j < s1Genotypes.length; j++) {
							if (s1Genotypes[j] == s2Genotypes[j]) {
								markerHits[j]++;
								totalHits++;
							}
						}
					}
				}
			}
		}
		
		if (pairsChecked == 0) {
			log.reportTimeError("No duplicates could be compared, duplicate concordance cannot be calculated");
			return null;
		}
		
		double[] markerConcordance = new double[markerHits.length];
		for (int i = 0; i < markerHits.length; i++) {
			markerConcordance[i] = (double)markerHits[i] / pairsChecked;
		}
		
		double projectConcordance = (double)totalHits / (markerHits.length * pairsChecked);
		
		return new DuplicateConcordance(projectConcordance, markerConcordance);
			
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		Project proj = null;
		String targetMarkersFile = null;

		String usage = "\n" +
		"cnv.qc.DuplicateConcordance requires 1-2 arguments\n" +
		"   (1) Project properties filename (i.e. proj=" + cnv.Launch.getDefaultDebugProjectFile(false) + " (not the default))\n"+
		"   (2) Target markers file (i.e. targetMarkersFile=" + targetMarkersFile + " (default))\n" +
		"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				proj = new Project(args[i].split("=")[1], false);
                numArgs--;
			} else if (args[i].startsWith("targetMarkersFile=")) {
				targetMarkersFile = args[i].split("=")[1];
                numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0 || args.length < 1) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			DuplicateConcordance duplicateConcordance = calculateMarkerConcordances(proj, targetMarkersFile);
			if (duplicateConcordance != null) {
				proj.getLog().report("Project duplicate concordance is " + duplicateConcordance.getProjectConcordance());
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
}
