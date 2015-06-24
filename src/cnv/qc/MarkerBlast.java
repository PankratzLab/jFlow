package cnv.qc;

import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;

import common.Array;
import common.Logger;
import common.PSF;
import common.Positions;
import common.WorkerHive;
import common.ext;
import cnv.filesys.Project;
import cnv.manage.ExtProjectDataParser;
import filesys.Segment;
import seq.analysis.Blast;
import seq.analysis.Blast.BlastWorker;
import seq.analysis.Blast.FastaEntry;

public class MarkerBlast {
	private enum FILE_SEQUENCE_TYPE {
		/**
		 * like HumanOmni2-5-8-v1-2-A-Strand-Report-FDT.txt
		 */
		STRAND_REPORT, /**
		 * Like HumanExome-12-v1-0-B.csv
		 */
		MANIFEST_FILE;
	}

	private static void blastEm(Project proj, String fileSeq, FILE_SEQUENCE_TYPE type, String fastaDb, int blastWordSize, int reportWordSize, int numThreads, boolean reportToTmp) {
		Blast blast = new Blast(fastaDb, blastWordSize, reportWordSize, proj.getLog(), true, true);
		WorkerHive<Blast.BlastResultsSummary[]> hive = new WorkerHive<Blast.BlastResultsSummary[]>(numThreads, 10, proj.getLog());
		FastaEntry[] fastaEntries = getMarkerFastaEntries(proj, fileSeq, type);
		ArrayList<FastaEntry[]> splits = Array.splitUpArray(fastaEntries, numThreads, proj.getLog());
		BlastWorker[] workers = new BlastWorker[splits.size()];
		String[] tmps = new String[workers.length];
		if (fastaEntries != null && fastaEntries.length > 0) {
			for (int i = 0; i < splits.size(); i++) {
				tmps[i] = ext.rootOf(fileSeq, false) + ".tmp" + i;
				workers[i] = new BlastWorker(blast, splits.get(i), reportToTmp ? tmps[i] : null);
			}
		}
		hive.addCallables(workers);
		hive.setReportEvery(1);
		hive.execute(true);
		ArrayList<Blast.BlastResultsSummary[]> bSummaries = hive.getResults();
		String output = ext.rootOf(fileSeq, false) + ".blasted";

		summarize(proj, bSummaries, output, proj.getLog());
		String outputOneHitWonders = ext.rootOf(fileSeq, false) + ".oneHitWonders";
		summarizeOneHitWonders(proj, bSummaries, outputOneHitWonders, proj.getLog());
		// Files.cat(tmps, ext.rootOf(fileSeq, false) + ".blasted", new int[] {}, proj.getLog());
	}

	private static void summarizeOneHitWonders(Project proj, ArrayList<Blast.BlastResultsSummary[]> bSummaries, String outputFile, Logger log) {
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(outputFile));
			Hashtable<String, Integer> indices = proj.getMarkerIndices();
			byte[] chrs = proj.getMarkerSet().getChrs();
			int[] pos = proj.getMarkerSet().getPositions();
			for (int i = 0; i < bSummaries.size(); i++) {
				for (int j = 0; j < bSummaries.get(i).length; j++) {
					if (bSummaries.get(i)[j].getNumPerfectMatches() == 1) {
						int currentMarkerIndex = indices.get(bSummaries.get(i)[j].getName());
						Segment markerSegment = new Segment(chrs[currentMarkerIndex], pos[currentMarkerIndex]-1, pos[currentMarkerIndex] + 1);
						Segment blastPerfectSegment = bSummaries.get(i)[j].getPerfectMatchSegment();
						if (blastPerfectSegment != null && blastPerfectSegment.overlaps(markerSegment)) {
							writer.println(bSummaries.get(i)[j].getName() + "\t" + chrs[currentMarkerIndex] + "\t" + pos[currentMarkerIndex] + "\t" + blastPerfectSegment.getChr() + "\t" + blastPerfectSegment.getStart() + "\t" + blastPerfectSegment.getStop());
						}
					}
				}
			}

			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + outputFile);
			log.reportException(e);
		}
	}

	private static void summarize(Project proj, ArrayList<Blast.BlastResultsSummary[]> bSummaries, String outputFile, Logger log) {
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(outputFile));
			Hashtable<String, Integer> indices = proj.getMarkerIndices();
			byte[] chrs = proj.getMarkerSet().getChrs();
			writer.print("Marker\tDesignContig\tNumTotalContigs");
			double[] bins = bSummaries.get(0)[0].getDynamicHistogram().getBins();
			HashSet<String> allContigs = new HashSet<String>();
			ArrayList<String> uniqContigs = new ArrayList<String>();
			for (int i = 0; i < bSummaries.size(); i++) {
				for (int j = 0; j < bSummaries.get(i).length; j++) {
					allContigs.addAll(bSummaries.get(i)[j].getHitCounts().keySet());
				}
			}
			uniqContigs.addAll(allContigs);
			for (int i = 0; i < bins.length; i++) {
				writer.print("\t" + ext.formDeci(bins[i], 0) + "_PercentIdentity");
			}
			for (String contig : uniqContigs) {
				writer.print("\t" + contig);
			}
			writer.println();
			for (int i = 0; i < bSummaries.size(); i++) {
				for (int j = 0; j < bSummaries.get(i).length; j++) {
					int[] counts = bSummaries.get(i)[j].getDynamicHistogram().getCounts();
					int contigHits = bSummaries.get(i)[j].getHitCounts().size();
					String chr = Positions.getChromosomeUCSC(chrs[indices.get(bSummaries.get(i)[j].getName())], true);
					writer.print(bSummaries.get(i)[j].getName() + "\t" + chr + "\t" + contigHits + "\t" + Array.toStr(counts));
					for (String contig : uniqContigs) {
						if (bSummaries.get(i)[j].getHitCounts().containsKey(contig)) {
							writer.print("\t" + bSummaries.get(i)[j].getHitCounts().get(contig));
						} else {
							writer.print("\t" + "0");
						}
					}
					writer.println();
				}

			}

			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + outputFile);
			log.reportException(e);
		}
	}

	private static boolean isPerfectlyGoodMatch(Blast.BlastResultsSummary bSummary, byte chr) {
		if (bSummary.getNumPerfectMatches() != 1) {
			return false;
		}

		return true;

	}

	private static FastaEntry[] getMarkerFastaEntries(Project proj, String strandReportFile, FILE_SEQUENCE_TYPE type) {
		FastaEntry[] fastaEntries = new FastaEntry[0];
		ExtProjectDataParser.Builder builder = new ExtProjectDataParser.Builder();

		switch (type) {
		case MANIFEST_FILE:
			builder.separator(",");
			builder.dataKeyColumnName("Name");
			builder.stringDataTitles(new String[] { "AlleleA_ProbeSeq" });
			builder.headerFlags(new String[] { "Name", "AlleleA_ProbeSeq" });
			break;
		case STRAND_REPORT:
			builder.dataKeyColumnName("SNP_Name");
			builder.commentString("##");
			builder.stringDataTitles(new String[] { "Design_Seq" });
			break;
		default:
			break;

		}
		builder.sampleBased(false);
		builder.treatAllNumeric(false);
		builder.requireAll(false);
		builder.verbose(true);
		try {

			ExtProjectDataParser parser = builder.build(proj, strandReportFile);
			parser.determineIndicesFromTitles();
			parser.loadData();
			ArrayList<FastaEntry> entries = new ArrayList<FastaEntry>(Array.booleanArraySum(parser.getDataPresent()));
			for (int i = 0; i < parser.getDataPresent().length; i++) {
				if (parser.getDataPresent()[i]) {

					String seq = parser.getStringData()[0][i];
					// int leftBracket = seq.indexOf("[");
					// int rightBrackt = seq.indexOf("]");
					// int slash = seq.indexOf("/");
					// String A = seq.substring(leftBracket + 1, slash);
					// String B = seq.substring(slash + 1, rightBrackt);
					// String left = seq.substring(0, leftBracket);
					// String right = seq.substring(rightBrackt + 1);
					// String seqA = left + A + right;
					// String seqB = left + B + right;
					// if (!A.equals("-")) {
					// entries.add(new FastaEntry(parser.getDataToLoad()[i] + "_A", seqA));
					// }
					// if (!B.equals("-")) {
					// entries.add(new FastaEntry(parser.getDataToLoad()[i] + "_B", seqB));
					// }
					entries.add(new FastaEntry(parser.getDataToLoad()[i], seq));
				}
			}
			fastaEntries = entries.toArray(new FastaEntry[entries.size()]);
		} catch (FileNotFoundException e) {
			proj.getLog().reportFileNotFound(strandReportFile);
			e.printStackTrace();
		}
		proj.getLog().reportTimeInfo("Found " + fastaEntries.length + " marker sequences to blast");
		return fastaEntries;
	}

	// public static void test() {
	// String file = "/home/pankrat2/shared/aric_exome_chip/Blast/HumanExome-12-v1-0-B.csv";
	// String fastaDb = "/home/pankrat2/public/bin/ref/hg19_canonical.fa";
	// blastEm(new Project("/home/pankrat2/lanej/projects/aric_exome.properties", false), file, FILE_SEQUENCE_TYPE.MANIFEST_FILE, fastaDb, 2);
	// }

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		String fastaDb = "hg19_canonical.fa";
		String fileSeq = "HumanExome-12-v1-0-B.csv";
		int numThreads = 4;
		int blastWordSize = 30;
		int reportWordSize = 40;
		boolean report = false;
		FILE_SEQUENCE_TYPE fSequence_TYPE = FILE_SEQUENCE_TYPE.MANIFEST_FILE;
		String usage = "\n" + "cnv.qc.MarkerBlast requires 3 arguments\n";
		usage += "   (1) Project file name (i.e. proj=" + filename + " (default))\n" + "";
		usage += "   (2) full path to a fasta file to search (i.e. fastaDb=" + fastaDb + " (default))\n" + "";
		usage += "   (3) full path to an Illumina manifest file  (i.e. fileSeq=" + fileSeq + " (default))\n" + "";
		usage += "   (4) number of threads to use  (i.e. " + PSF.Ext.NUM_THREADS_COMMAND + numThreads + " (default))\n" + "";
		usage += "   (5) word size for initial match in blast db  (i.e. blastWordSize=" + blastWordSize + " (default))\n" + "";
		usage += "   (6) number of base pairs with an exact match to report (i.e. reportWordSize=" + reportWordSize + " (default))\n" + "";
		usage += "   (7) report results to temporary files (i.e. -report (not the default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("fastaDb=")) {
				fastaDb = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("fileSeq=")) {
				fileSeq = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
				numThreads = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("blastWordSize=")) {
				blastWordSize = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("reportWordSize=")) {
				reportWordSize = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("-report")) {
				report = true;
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
			Project proj = new Project(filename, false);
			blastEm(proj, fileSeq, fSequence_TYPE, fastaDb, blastWordSize, reportWordSize, numThreads, report);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
