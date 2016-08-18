package org.genvisis.seq.analysis;

import java.io.FileWriter;
import java.io.PrintWriter;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.ext;
import org.genvisis.seq.manage.VCFOps;

public class PlinkSeqIteration {

	public static void runPlinkSeq(String dir, String phenoFile, String resourceDirectory, String outputRoot, String[] locGroups, boolean overwriteExisting, int macStart, int macStop, int numThreads, Logger log) {
		String[] vcfs = Files.listFullPaths(dir, ".vcf.gz", false);
		String finalSummary = dir + "finalBurden.summary";
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(finalSummary));
			writer.println("NUMVAR\tMAC_GENE_CUT\tANALYSIS\tPASSING_UNITS\tBonferroni\tNumSig");
			for (int i = 0; i < vcfs.length; i++) {
				int numVar = VCFOps.getNumberOfVariants(vcfs[i]);
				int max = 500;
				for (int j = macStart; j < macStop; j++) {
					String projName = ext.rootOf(vcfs[i]) + "_mac_" + j;
					PlinkSeq.runPlinkSeq(projName, vcfs[i], phenoFile, resourceDirectory, outputRoot, locGroups, overwriteExisting, true, true, j == 0 ? "0" : j + ":" + max, numThreads, log);
					String summaryFile = dir + "pseqProj_" + projName + "/assoc/burden.summary";
					String[] start = Files.getHeaderOfFile(summaryFile, new Logger());
					int[] indices = new int[start.length];
					for (int k = 0; k < indices.length; k++) {
						indices[k] = k;
					}
					String[] summary = HashVec.loadFileToStringArray(summaryFile, false, null, false);
					for (int k = 0; k < summary.length; k++) {
						writer.println(numVar + "\t" + j + "\t" + summary[k]);
					}
					writer.flush();
				}
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + finalSummary);
			log.reportException(e);
		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String dir = null;
		String phenoFile = null;
		String resourceDirectory = null;
		String projName = "pseq_project";
		String outputRoot = "pseq_analysis";
		String[] locGroups = new String[] { "refseq" };
		int numThreads = 4;
		int macStart = 0;
		int macStop = 5;

		boolean overwriteExisting = false;

		// TODO, use just a project file Name to load everything, loc.group args (refseq)
		String usage = "\n" + "seq.analysis.PlinkSeq requires 3 arguments\n";
		usage += "   (1) full path to a dir with vcfs (i.e. vcf= (no default))\n" + "";
		usage += "   (2) full path to a plinkseq phenotype file (i.e. phenoFile= (no default))\n" + "";
		usage += "   (3) full path to a plinkseq resource directory  (i.e. resourceDirectory= (no default))\n" + "";
		usage += "   (4) overwrite existing files  (i.e. -overwriteExisting (default))\n" + "";
		usage += "   (5) analysis output root  (i.e. outputRoot=" + outputRoot + " (default))\n" + "";
		usage += "   (6) project name  (i.e. projName=" + projName + " (default))\n" + "";
		usage += "   (7) overwrite existing files  (i.e. -overwriteExisting (default))\n" + "";
		usage += "   (8) comma -delimted loc groups to test in the association  (i.e. locGroups=" + Array.toStr(locGroups, ",") + " (default))\n" + "";
		usage += "   (9) mac start  (i.e. macStart=" + macStart + " (default))\n" + "";
		usage += "   (9) mac stop  (i.e. macStop=" + macStop + " (default))\n" + "";

		usage += "   (10) " + PSF.Ext.getNumThreadsCommand(10, numThreads);

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("dir=")) {
				dir = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("phenoFile=")) {
				phenoFile = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("outputRoot=")) {
				outputRoot = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("overwriteExisting=")) {
				overwriteExisting = true;
				numArgs--;
			} else if (args[i].startsWith("resourceDirectory=")) {
				resourceDirectory = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("locGroups=")) {
				locGroups = ext.parseStringArg(args[i], "").split(",");
				numArgs--;
			} else if (args[i].startsWith("projName=")) {
				projName = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
				numThreads = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("macStart=")) {
				macStart = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("macStop=")) {
				macStop = ext.parseIntArg(args[i]);
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
			runPlinkSeq(dir, phenoFile, resourceDirectory, outputRoot, locGroups, overwriteExisting, macStart, macStop, numThreads, new Logger());
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
