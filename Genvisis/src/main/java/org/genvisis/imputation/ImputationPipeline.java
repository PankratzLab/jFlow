package org.genvisis.imputation;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.PlinkData;
import org.genvisis.cnv.manage.VCFData;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Sort;
import org.genvisis.common.ext;
import org.genvisis.gwas.FurtherAnalysisQc;
import org.genvisis.gwas.Qc;

import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Sets;

/**
 * Export genotypes to either PLINK or VCF format and, depending on selected options, generate
 * scripts to run ShapeIt and Minimac programs externally. <br />
 * <br />
 * Imputation TODO's:<br />
 * <ul>
 * <li>Make content aware (i.e., if PLINK files already exist, skip chrs that exist (or overwrite,
 * or fail))<br />
 * </li>
 * <li>Multithread data file export<br />
 * </li>
 * <li>Use ScriptExecutor for ShapeIt/Minimac scripts</li>
 * </ul>
 */
public class ImputationPipeline {

	Project proj;
	Set<String> dropMarkers = new HashSet<String>();
	Set<String> dropSamples = new HashSet<String>();
	Set<String> keepMarkers = new HashSet<String>();
	Set<String> keepSamples = new HashSet<String>();
	Map<String, Marker> prepMarkers = new HashMap<String, Marker>();

	Set<String> markersToExport;

	public ImputationPipeline(Project proj, String referenceFile) {
		this.proj = proj;
		ImputationPrep prep = new ImputationPrep(proj, referenceFile);
		Set<Marker> markers = prep.getMatchingMarkers();
		// PrintWriter writer =
		// Files.getAppropriateWriter("/home/pankrat2/cole0482/imputationMarkers.txt");
		for (Marker m : markers) {
			prepMarkers.put(m.getName(), m);
			// writer.println(m.getName());
		}
		// writer.flush();
		// writer.close();
	}

	public void loadDefaultDropFiles(String plinkDir) {
		String dir = plinkDir + Qc.QC_SUBDIR + FurtherAnalysisQc.FURTHER_ANALYSIS_DIR;
		String mark = dir + FurtherAnalysisQc.MARKER_QC_DROPS;
		String samp = dir + FurtherAnalysisQc.SAMPLE_QC_DROPS;
		setMarkersToDropFile(mark);
		setSamplesToDropFile(samp);
	}

	public void setSamplesToDropFile(String samplesToDropFile) {
		if (!Files.exists(samplesToDropFile)) {
			proj.getLog().reportTimeWarning("Sample drop file doesn't exist: " + samplesToDropFile);
			return;
		}
		dropSamples = HashVec.loadFileToHashSet(samplesToDropFile, new int[] {0, 1}, "\t", false);
		Set<String> missingIDs = Sets.difference(dropSamples, proj.getSampleData(false)
																															.getSampleIDLookup().keySet());
		if (!missingIDs.isEmpty()) {
			proj.getLog().reportError("Not all Samples to drop could be found in SampleData: "
																+ ext.listWithCommas(missingIDs));
		}
	}

	public void setSamplesToKeepFile(String samplesToKeepFile) {
		if (!Files.exists(samplesToKeepFile)) {
			proj.getLog().reportTimeWarning("Sample keep file doesn't exist: " + samplesToKeepFile);
			return;
		}
		keepSamples = HashVec.loadFileToHashSet(samplesToKeepFile, new int[] {0, 1}, "\t", false);
		Set<String> missingIDs = Sets.difference(keepSamples, proj.getSampleData(false)
																															.getSampleIDLookup().keySet());
		if (!missingIDs.isEmpty()) {
			proj.getLog().reportError("Not all Samples to keep could be found in SampleData: "
																+ ext.listWithCommas(missingIDs));
		}
	}

	public void setMarkersToDropFile(String markersToDropFile) {
		if (!Files.exists(markersToDropFile)) {
			proj.getLog().reportTimeWarning("Marker drop file doesn't exist: " + markersToDropFile);
			return;
		}
		dropMarkers = HashVec.loadFileToHashSet(markersToDropFile, false);
		Set<String> missingIDs = Sets.difference(dropMarkers,
																						 proj.getMarkerSet().getMarkerNameMap().keySet());
		if (!missingIDs.isEmpty()) {
			proj.getLog().reportError("Not all Markers to drop could be found in Project: "
																+ ext.listWithCommas(missingIDs));
		}
	}

	public void setMarkersToKeepFile(String markersToKeepFile) {
		if (!Files.exists(markersToKeepFile)) {
			proj.getLog().reportTimeWarning("Marker keep file doesn't exist: " + markersToKeepFile);
			return;
		}
		keepMarkers = HashVec.loadFileToHashSet(markersToKeepFile, false);
		Set<String> missingIDs = Sets.difference(keepMarkers,
																						 proj.getMarkerSet().getMarkerNameMap().keySet());
		if (!missingIDs.isEmpty()) {
			proj.getLog().reportError("Not all Markers to keep could be found in Project: "
																+ ext.listWithCommas(missingIDs));
		}
	}

	private Set<String> getSamplesToExport() {
		HashSet<String> samplesToExport = new HashSet<String>();
		if (keepSamples == null) {
			for (String s : proj.getSamples()) {
				if (dropSamples == null || !dropSamples.contains(s)) {
					samplesToExport.add(s);
				}
			}
		} else {
			samplesToExport.addAll(keepSamples);
		}
		return samplesToExport;
	}

	private Set<String> getMarkersToExport() {
		if (markersToExport == null) {
			markersToExport = new HashSet<String>(Arrays.asList(proj.getMarkerNames()));
			if (keepMarkers != null && keepMarkers.size() > 0) {
				markersToExport.retainAll(keepMarkers);
			}
			markersToExport.removeAll(dropMarkers);
			markersToExport.retainAll(prepMarkers.keySet());
		}
		return markersToExport;
	}

	private Set<String> getChrMarkers(int chr) {
		Set<String> markersToExport = getMarkersToExport();
		Set<String> chrMarkers = new HashSet<String>();
		for (String m : markersToExport) {
			if (prepMarkers.get(m).getChr() == (byte) chr) {
				chrMarkers.add(m);
			}
		}
		return chrMarkers;
	}

	private ArrayList<String> getMarkersSortedNoDupes(int chr) {
		Set<String> chrMarkers = getChrMarkers(chr);

		int[] pos = new int[chrMarkers.size()];
		String[] mkr = chrMarkers.toArray(new String[chrMarkers.size()]);
		for (int i = 0; i < mkr.length; i++) {
			pos[i] = prepMarkers.get(mkr[i]).getPosition();
		}

		int[] indices = Sort.getSortedIndices(pos);

		ArrayList<String> mkrs = new ArrayList<String>();
		for (int i = 0; i < indices.length; i++) {
			if (i == 0 || pos[indices[i]] != pos[indices[i - 1]]) { // skip if prev (in sorted array) was
																															// same position
				mkrs.add(mkr[indices[i]]);
			}
		}
		return mkrs;
	}

	public void exportToPlink(String plinkDirAndRoot, int[] chrs) {
		// TODO (??) Only alphanumeric characters in FID/IID
		(new File(ext.parseDirectoryOfFile(plinkDirAndRoot + ".bim"))).mkdirs();
		HashSet<String> toDrop = new HashSet<String>();
		if (keepSamples != null && keepSamples.size() > 0) {
			for (String s : proj.getSamples()) {
				if (keepSamples.contains(s) && !dropSamples.contains(s)) {
					continue;
				}
				toDrop.add(s);
			}
		} else {
			toDrop.addAll(dropSamples);
		}
		String[] writtenDNAs = PlinkData.createFamFile(proj, plinkDirAndRoot, toDrop, true);
		if (writtenDNAs == null) {
			// TODO error
			return;
		}
		int[] indicesOfTargetSamplesInProj = PlinkData.getIndicesOfTargetSamplesInProj(proj,
																																									 writtenDNAs,
																																									 proj.getLog());

		// UNUSED - could potentially apply
		String clusterFilterFileName = proj.CLUSTER_FILTER_COLLECTION_FILENAME.getValue();

		// TODO multi-thread
		for (int chr : chrs) {
			proj.getLog().report("Exporting chr" + chr);
			ArrayList<String> mkrs = getMarkersSortedNoDupes(chr);

			String[] targetMarkers = mkrs.toArray(new String[mkrs.size()]);
			int[] indicesOfTargetMarkersInProj = new int[targetMarkers.length];
			HashMap<String, Byte> chrsOfTargetMarkers = new HashMap<String, Byte>();
			HashMap<String, Integer> posOfTargetMarkers = new HashMap<String, Integer>();
			PlinkData.getIndicesOfTargetMarkers(proj, targetMarkers, indicesOfTargetMarkersInProj,
																					chrsOfTargetMarkers, posOfTargetMarkers, proj.getLog());

			String dirAndRoot = plinkDirAndRoot + "_chr" + chr;
			boolean success = PlinkData.createBedFileSnpMajor10KperCycle(proj,
																																	 ImmutableSet.copyOf(targetMarkers),
																																	 chrsOfTargetMarkers,
																																	 posOfTargetMarkers,
																																	 indicesOfTargetSamplesInProj,
																																	 clusterFilterFileName,
																																	 dirAndRoot, proj.getLog());

			if (success) {
				PrintWriter refWriter = Files.getAppropriateWriter(dirAndRoot + "_alleles.ref");
				for (String s : targetMarkers) {
					refWriter.println(s + "\t" + prepMarkers.get(s).getRef());
				}
				refWriter.flush();
				refWriter.close();
				Files.copyFile(plinkDirAndRoot + ".fam", dirAndRoot + ".fam");
			}
		}
	}

	public void exportToVCF(String vcfDirAndRoot, int[] chrs, boolean useGRC) {
		String[] samplesToExport = getSamplesToExport().toArray(new String[0]);
		String[] markersToExport = getMarkersToExport().toArray(new String[0]);

		VCFData.exportGenvisisToVCF(proj, samplesToExport, markersToExport, true, useGRC, vcfDirAndRoot);
	}

	private static class ImputationPipeRunner {

		public static void runVCF(String projPropFile, int[] chrs, String refFile, String plinkSubdir,
															String vcfDirAndRoot, boolean useGRC) {
			ImputationPipeline ip = setupPipe(projPropFile, refFile, plinkSubdir);
			ip.exportToVCF(vcfDirAndRoot, chrs, useGRC);
		}

		public static void runPlink(String projPropFile, int[] chrs, String refFile,
																String plinkSubdir, String outputDirAndRoot) {
			ImputationPipeline ip = setupPipe(projPropFile, refFile, plinkSubdir);
			ip.exportToPlink(outputDirAndRoot, chrs);
		}

		public static void runShapeIt(String projPropFile, int[] chrs, String plinkFileDir,
																	String plinkPrefix, String outDir) {
			Project proj = new Project(projPropFile, false);
			new ImputationImpl.ShapeIt(proj, chrs, plinkFileDir, plinkPrefix, outDir).createScripts();
		}

		public static void runMinimac(String projPropFile, int[] chrs, String hapsDir, String outDir) {
			Project proj = new Project(projPropFile, false);
			new ImputationImpl.MiniMac(proj, chrs, hapsDir, outDir).createScripts();
		}

		public static void runPlinkAndShapeIt(String projPropFile, int[] chrs, String refFile,
																					String plinkSubdir, String outputDir) {
			runPlink(projPropFile, chrs, refFile, plinkSubdir, outputDir + "/plink/plink");
			runShapeIt(projPropFile, chrs, outputDir + "/plink/", "plink_chr", outputDir + "/haps/");
		}

		public static void runPlinkShapeItAndMinimac(String projPropFile, int[] chrs, String refFile,
																								 String plinkSubdir, String outputDir) {
			runPlinkAndShapeIt(projPropFile, chrs, refFile, plinkSubdir, outputDir);
			runMinimac(projPropFile, chrs, outputDir + "/haps/", outputDir);
		}

		private static ImputationPipeline setupPipe(String projPropFile, String refFile,
																								String plinkSubdir) {
			Project proj = new Project(projPropFile, false);
			ImputationPipeline ip = new ImputationPipeline(proj, refFile);
			ip.loadDefaultDropFiles(proj.PROJECT_DIRECTORY.getValue() + plinkSubdir);
			return ip;
		}

	}

	private enum IMPUTATION_PIPELINE_PATH {
		VCF_ONLY,
		PLINK_ONLY,
		PLINK_SHAPE,
		PLINK_SHAPE_MINI,
		SHAPE,
		MINI;
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String projFile = null;
		String refFile = null;
		String plinkSubdir = null;
		String outDirAndRoot = null;
		String hapsDir = null;
		String outDir = null;
		String plinkPrefix = null;
		int[] chrs = null;
		boolean useGRC = true;
		IMPUTATION_PIPELINE_PATH path = null;

		String usage = "\n" +
									 "org.genvisis.imputation.ImputationPipeline requires 0-1 arguments\n" +
									 "   (1) Project properties filename (i.e. proj=" + projFile + " (default))\n" +
									 "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h")
					|| args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				projFile = ext.parseStringArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("ref=")) {
				refFile = ext.parseStringArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("plinkDir=")) {
				plinkSubdir = ext.parseStringArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("plinkPrefix=")) {
				plinkPrefix = ext.parseStringArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("outDirAndRoot=")) {
				outDirAndRoot = ext.parseStringArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("outDir=")) {
				outDir = ext.parseStringArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("hapsDir=")) {
				hapsDir = ext.parseStringArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("chrs=")) {
				chrs = ext.parseIntArrayArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("type=")) {
				path = IMPUTATION_PIPELINE_PATH.valueOf(ext.parseStringArg(args[i]));
				numArgs--;
			} else if (args[i].startsWith("useGRC=")) {
				useGRC = ext.parseBooleanArg(args[i]);
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
			switch (path) {
				case VCF_ONLY:
					ImputationPipeRunner.runVCF(projFile, chrs, refFile, plinkSubdir, outDirAndRoot, useGRC);
					break;
				case PLINK_ONLY:
					ImputationPipeRunner.runPlink(projFile, chrs, refFile, plinkSubdir, outDirAndRoot);
					break;
				case PLINK_SHAPE:
					ImputationPipeRunner.runPlinkAndShapeIt(projFile, chrs, refFile, plinkSubdir, outDir);
					break;
				case PLINK_SHAPE_MINI:
					ImputationPipeRunner.runPlinkShapeItAndMinimac(projFile, chrs, refFile, plinkSubdir,
																												 outDir);
					break;
				case MINI:
					ImputationPipeRunner.runMinimac(projFile, chrs, hapsDir, outDir);
					break;
				case SHAPE:
					ImputationPipeRunner.runShapeIt(projFile, chrs, plinkSubdir, plinkPrefix, outDir);
					break;
				default:
					System.err.println("Error - unrecognized imputation path: " + path);
					break;
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
