package cnv.qc;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

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
	private int markersChecked;
	private int duplicatePairsChecked;
	
	
	/**
	 * @param discordantCalls number of marker calls that did not match across
	 * @param markersChecked number of markers checked
	 * @param duplicatePairsChecked number of duplicate pairs checked
	 */
	private DuplicateConcordance(int discordantCalls, int markersChecked, int duplicatePairsChecked) {
		super();
		int totalChecks = markersChecked * duplicatePairsChecked;
		this.projectConcordance = (double)(totalChecks - discordantCalls) / totalChecks;
		this.markersChecked = markersChecked;
		this.duplicatePairsChecked = duplicatePairsChecked;
	}

	public double getProjectConcordance() {
		return projectConcordance;
	}

	public int getMarkersChecked() {
		return markersChecked;
	}

	public int getDuplicatePairsChecked() {
		return duplicatePairsChecked;
	}
	
	public String getConcordanceString() {
		return "Duplicate Concordance was calculated to be " + projectConcordance + " using " + duplicatePairsChecked + " pairs of duplicates at " + markersChecked + " markers.";
	}

	/**
	 * 
	 * @param proj Project to calculate duplicate concordance for
	 * @param targetMarkers Markers to use in concordance checks or null to check all markers
	 * @return
	 */
	public static DuplicateConcordance calculateDuplicateConcordances(Project proj, String[] targetMarkers) {
		Logger log = proj.getLog();
		ClusterFilterCollection clusterFilterCollection = proj.getClusterFilterCollection();
		String[] markerNames;
		int[] markerIndices;
		if (targetMarkers == null) {
			markerNames = proj.getMarkerNames();
			markerIndices = null;
		} else {
			markerNames = targetMarkers;
			markerIndices = ext.indexLargeFactors(markerNames, proj.getMarkerNames(), true, log, false, false);
			for (int i = 0; i < markerIndices.length; i++) {
				if (markerIndices[i] == -1) {
					log.reportTimeError("Marker " + markerNames[i] + " could not be found in project");
					return null;
				}
			}
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
		int discordantCalls = 0;
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
							s1Genotypes = sample1.getAB_Genotypes(markerIndices);
							s2Genotypes = sample2.getAB_Genotypes(markerIndices);
						} else {
							s1Genotypes = sample1.getAB_GenotypesAfterFilters(markerNames, markerIndices, clusterFilterCollection, 0.0f);
							s2Genotypes = sample2.getAB_GenotypesAfterFilters(markerNames, markerIndices, clusterFilterCollection, 0.0f);
						}
						for (int i = 0; i < s1Genotypes.length; i++) {
							if (s1Genotypes[i] != s2Genotypes[i]) {
								discordantCalls++;
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
		
		return new DuplicateConcordance(discordantCalls, markerNames.length, pairsChecked);
			
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		Project proj = null;
		String markerKeeps = null;
		String markerDrops = null;

		String usage = "\n" +
		"cnv.qc.DuplicateConcordance requires 1-2 arguments\n" +
		"   (1) Project properties filename (i.e. proj=" + cnv.Launch.getDefaultDebugProjectFile(false) + " (not the default))\n"+
		"AND\n" + 
		"   (2) File of markers to use (i.e. markerKeeps=keeps.txt (not the default))\n" +
		"OR\n" +
		"   (2) File of markers to not use (i.e. markerDrops=drops.txt (not the default))\n" +
		"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				proj = new Project(args[i].split("=")[1], false);
                numArgs--;
			} else if (args[i].startsWith("markerKeeps=")) {
				markerKeeps = args[i].split("=")[1];
                numArgs--;
			} else if (args[i].startsWith("markerDrops=")) {
				markerDrops = args[i].split("=")[1];
                numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			if (proj == null) {
				System.err.println("Project must be defined");
				System.err.println(usage);
				System.exit(1);
			}
			if (markerKeeps != null && markerDrops != null) {
				System.err.println("Include a marker keeps or drops file but not both");
				System.err.println(usage);
				System.exit(1);
			}
			String[] targetMarkers;
			if (markerKeeps != null) {
				targetMarkers = proj.getTargetMarkers(markerKeeps);
			} else if (markerDrops != null) {
				Set<String> excludes = HashVec.loadFileToHashSet(markerDrops, false);
				ArrayList<String> markers = new ArrayList<String>();
				for (String marker : proj.getMarkerNames()) {
					if (!excludes.contains(marker)) {
						markers.add(marker);
					}
				}
				targetMarkers = Array.toStringArray(markers);
			} else {
				targetMarkers = null;
			}
			
			DuplicateConcordance duplicateConcordance = calculateDuplicateConcordances(proj, targetMarkers);
			if (duplicateConcordance != null) {
				proj.getLog().report(duplicateConcordance.getConcordanceString());
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
}
