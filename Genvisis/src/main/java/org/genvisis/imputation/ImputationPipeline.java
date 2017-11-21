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
import org.genvisis.common.ArrayUtils;
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

	public static final String PROJ_ARG = "proj=";
	public static final String REF_ARG = "ref=";
	public static final String PLINK_DIR_ARG = "plinkDir=";
	public static final String PLINK_PREFIX_ARG = "plinkPrefix=";
	public static final String OUT_DIR_AND_ROOT_ARG = "outDirAndRoot=";
	public static final String OUT_DIR_ARG = "outDir=";
	public static final String HAPS_DIR_ARG = "hapsDir=";
	public static final String CHRS_ARG = "chrs=";
	public static final String RUN_TYPE_ARG = "type=";
	public static final String USE_GRC_ARG = "useGRC=";

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
		for (Marker m : markers) {
			prepMarkers.put(m.getName(), m);
		}
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
																																	 proj.GC_THRESHOLD.getValue()
																																										.floatValue(),
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

		VCFData.exportGenvisisToVCF(proj, samplesToExport, markersToExport, true, useGRC, chrs,
																vcfDirAndRoot);
	}

	protected static class ImputationPipeRunner {

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
			Project proj = new Project(projPropFile);
			new ImputationImpl.ShapeIt(proj, chrs, plinkFileDir, plinkPrefix, outDir).createScripts();
		}

		public static void runMinimac(String projPropFile, int[] chrs, String hapsDir, String outDir) {
			Project proj = new Project(projPropFile);
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
			Project proj = new Project(projPropFile);
			ImputationPipeline ip = new ImputationPipeline(proj, refFile);
			if (Files.exists(plinkSubdir)) {
				ip.loadDefaultDropFiles(proj.PROJECT_DIRECTORY.getValue() + plinkSubdir);
			}
			return ip;
		}

	}

	public enum IMPUTATION_PIPELINE_PATH {
		VCF_ONLY,
		PLINK_ONLY,
		PLINK_SHAPEIT,
		PLINK_SHAPEIT_MINIMAC,
		SHAPEIT,
		MINIMAC;
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

		String usage = "\n"
									 +
									 "org.genvisis.imputation.ImputationPipeline requires 3+ arguments\n"
									 +
									 "   (1) Project properties filename (i.e. " + PROJ_ARG + projFile
									 + " (default))\n"
									 +
									 "   (2) A comma-separated list of chromosomes to export, or null for all (may be omitted) (i.e. "
									 + CHRS_ARG + chrs + " (default))\n"
									 +
									 "   (3) Imputation Pipeline path (i.e. " + RUN_TYPE_ARG + " one of "
									 + ArrayUtils.toStr(IMPUTATION_PIPELINE_PATH.values(), ", ") + "))\n"
									 +
									 "   --------------------- \n"
									 +
									 "   The following arguments may be necessary depending on chosen pipeline:\n"
									 +
									 "   (a) Reference Panel / Site List file, with mkr, chr, pos, ref, and alt columns (i.e. "
									 + REF_ARG + refFile + " (default))\n" +
									 "   (b) Subdirectory in which to create PLINK files (i.e. " + PLINK_DIR_ARG
									 + plinkSubdir + " (default))\n" +
									 "   (c) PLINK output prefix (i.e. " + PLINK_PREFIX_ARG + plinkPrefix
									 + " (default))\n" +
									 "   (d) Output directory and fileroot (i.e " + OUT_DIR_AND_ROOT_ARG
									 + outDirAndRoot
									 + " (default))\n" +
									 "   (e) Output directory (i.e " + OUT_DIR_ARG + outDir + " (default))\n" +
									 "   (f) Export contigs as 'chr1' instead of '1' (i.e. " + USE_GRC_ARG + useGRC
									 + " (default))\n" +
									 "   (g) Directory with output from ShapeIt (i.e. " + HAPS_DIR_ARG + hapsDir
									 + " (default))\n" +
									 "   --------------------- \n" +
									 "   Additional pipeline argument requirements are as follows:\n" +
									 "\tVCF_ONLY:\n" +
									 "\t\trefFile, plinkSubdir, outDirAndRoot, useGRC\n" +
									 "\tPLINK_ONLY:\n" +
									 "\t\trefFile, plinkSubdir, outDirAndRoot\n" +
									 "\tPLINK_SHAPEIT:\n" +
									 "\t\trefFile, plinkSubdir, outDir\n" +
									 "\tPLINK_SHAPEIT_MINIMAC:\n" +
									 "\t\trefFile, plinkSubdir\n" +
									 "\tMINIMAC:\n" +
									 "\t\thapsDir, outDir\n" +
									 "\tSHAPEIT:\n" +
									 "\t\tplinkSubdir, plinkPrefix, outDir\n" +
									 "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h")
					|| args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith(PROJ_ARG)) {
				projFile = ext.parseStringArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith(REF_ARG)) {
				refFile = ext.parseStringArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith(PLINK_DIR_ARG)) {
				plinkSubdir = ext.parseStringArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith(PLINK_PREFIX_ARG)) {
				plinkPrefix = ext.parseStringArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith(OUT_DIR_AND_ROOT_ARG)) {
				outDirAndRoot = ext.parseStringArg(args[i]);
				new File(ext.parseDirectoryOfFile(outDirAndRoot)).mkdirs();
				numArgs--;
			} else if (args[i].startsWith(OUT_DIR_ARG)) {
				outDir = ext.verifyDirFormat(ext.parseStringArg(args[i]));
				new File(outDir).mkdirs();
				numArgs--;
			} else if (args[i].startsWith(HAPS_DIR_ARG)) {
				hapsDir = ext.parseStringArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith(CHRS_ARG)) {
				chrs = ext.parseIntArrayArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith(RUN_TYPE_ARG)) {
				path = IMPUTATION_PIPELINE_PATH.valueOf(ext.parseStringArg(args[i]));
				numArgs--;
			} else if (args[i].startsWith(USE_GRC_ARG)) {
				useGRC = ext.parseBooleanArg(args[i]);
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}

		if (numArgs != 0 || path == null) {
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
				case PLINK_SHAPEIT:
					ImputationPipeRunner.runPlinkAndShapeIt(projFile, chrs, refFile, plinkSubdir, outDir);
					break;
				case PLINK_SHAPEIT_MINIMAC:
					ImputationPipeRunner.runPlinkShapeItAndMinimac(projFile, chrs, refFile, plinkSubdir,
																												 outDir);
					break;
				case MINIMAC:
					ImputationPipeRunner.runMinimac(projFile, chrs, hapsDir, outDir);
					break;
				case SHAPEIT:
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
