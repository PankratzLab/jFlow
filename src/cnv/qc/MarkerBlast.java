package cnv.qc;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;

import common.Array;
import common.Files;
import common.Logger;
import common.PSF;
import common.Positions;
import common.WorkerHive;
import common.ext;
import cnv.filesys.Project;
import cnv.filesys.Project.ARRAY;
import cnv.manage.ExtProjectDataParser;
import filesys.Segment;
import seq.analysis.Blast;
import seq.analysis.Blast.BlastWorker;
import seq.analysis.Blast.FastaEntry;

public class MarkerBlast {
	public enum FILE_SEQUENCE_TYPE {
		/**
		 * like HumanOmni2-5-8-v1-2-A-Strand-Report-FDT.txt
		 */
		STRAND_REPORT, /**
		 * Like HumanExome-12-v1-0-B.csv
		 */
		MANIFEST_FILE, /**
		 * Like GenomeWideSNP_6.probe_tab
		 */
		AFFY_ANNOT;
	}

	public static MarkerBlastResult blastEm(Project proj, String fileSeq, FILE_SEQUENCE_TYPE type, int blastWordSize, int reportWordSize, int numThreads, boolean reportToTmp) {
		String fastaDb = proj.REFERENCE_GENOME_FASTA_FILENAME.getValue();
		if (!Files.exists(fastaDb)) {
			proj.getLog().reportTimeError("Was not able to find reference genome defined by " + proj.REFERENCE_GENOME_FASTA_FILENAME.getName());
			return null;
		} else {

			Blast blast = new Blast(fastaDb, blastWordSize, reportWordSize, proj.getLog(), true, true);
			String dir = proj.PROJECT_DIRECTORY.getValue() + "Blasts/";
			new File(dir).mkdirs();

			String root = dir + ext.rootOf(fileSeq, true) + ".blasted.ws." + blastWordSize + ".rep." + reportWordSize;
			String output = root + ".blasted";
			String outputOneHitWonders = root + ".oneHitWonders";

			String[] tmps = new String[numThreads];
			for (int i = 0; i < numThreads; i++) {
				tmps[i] = root + ".tmp" + i;
			}

			if (!Files.exists(output) || !Files.exists(outputOneHitWonders)) {
				WorkerHive<Blast.BlastResultsSummary[]> hive = new WorkerHive<Blast.BlastResultsSummary[]>(numThreads, 10, proj.getLog());
				FastaEntry[] fastaEntries = getMarkerFastaEntries(proj, fileSeq, type);
				ArrayList<FastaEntry[]> splits = Array.splitUpArray(fastaEntries, numThreads, proj.getLog());

				BlastWorker[] workers = new BlastWorker[splits.size()];
				if (fastaEntries != null && fastaEntries.length > 0) {
					for (int i = 0; i < splits.size(); i++) {
						workers[i] = new BlastWorker(blast, splits.get(i), reportToTmp ? tmps[i] : null);
					}
				}

				hive.addCallables(workers);
				hive.setReportEvery(1);
				hive.execute(true);
				ArrayList<Blast.BlastResultsSummary[]> bSummaries = hive.getResults();

				summarize(proj, bSummaries, output, proj.getLog());
				summarizeOneHitWonders(proj, bSummaries, outputOneHitWonders, proj.getLog());
			} else {
				proj.getLog().reportFileExists(output);
				proj.getLog().reportFileExists(outputOneHitWonders);
			}
			proj.getLog().reportTimeWarning("Assuming that all sequences have length " + proj.getArrayType().getProbeLength() + " based on array type " + proj.getArrayType());
			MarkerBlastResult result = new MarkerBlastResult(output, outputOneHitWonders, blastWordSize, reportWordSize, proj.getArrayType().getProbeLength());
			result.setTmpFiles(tmps);
			return result;
		}
	}

	public static class MarkerBlastResult {
		private String output;
		private String outputOneHitWonders;
		private String[] tmpFiles;
		private int blastWordSize;
		private int reportWordSize;
		private int sequenceSize;

		public MarkerBlastResult(String output, String outputOneHitWonders, int blastWordSize, int reportWordSize, int sequenceSize) {
			super();
			this.output = output;
			this.outputOneHitWonders = outputOneHitWonders;
			this.blastWordSize = blastWordSize;
			this.reportWordSize = reportWordSize;
			this.sequenceSize = sequenceSize;
		}

		public int getSequenceSize() {
			return sequenceSize;
		}

		public void setTmpFiles(String[] tmpFiles) {
			this.tmpFiles = tmpFiles;
		}

		public String[] getTmpFiles() {
			return tmpFiles;
		}

		public int getBlastWordSize() {
			return blastWordSize;
		}

		public int getReportWordSize() {
			return reportWordSize;
		}

		public String getOutput() {
			return output;
		}

		public String getOutputOneHitWonders() {
			return outputOneHitWonders;
		}

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
						if (Array.countIf(bSummaries.get(i)[j].getPercentIdentityHistogram().getCounts(), 0) == bSummaries.get(i)[j].getPercentIdentityHistogram().getCounts().length - 1) {
							int currentMarkerIndex = indices.get(bSummaries.get(i)[j].getName());
							Segment markerSegment = new Segment(chrs[currentMarkerIndex], pos[currentMarkerIndex] - 1, pos[currentMarkerIndex] + 1);
							Segment blastPerfectSegment = bSummaries.get(i)[j].getPerfectMatchSegment();
							if (blastPerfectSegment != null && blastPerfectSegment.overlaps(markerSegment)) {
								writer.println(bSummaries.get(i)[j].getName() + "\t" + chrs[currentMarkerIndex] + "\t" + pos[currentMarkerIndex] + "\t" + blastPerfectSegment.getChr() + "\t" + blastPerfectSegment.getStart() + "\t" + blastPerfectSegment.getStop());
							}
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
			double[] bins = bSummaries.get(0)[0].getPercentIdentityHistogram().getBins();
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
					int[] counts = bSummaries.get(i)[j].getPercentIdentityHistogram().getCounts();
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

	private static FastaEntry[] getMarkerFastaEntries(Project proj, String strandReportFile, FILE_SEQUENCE_TYPE type) {
		FastaEntry[] fastaEntries = new FastaEntry[0];
		ExtProjectDataParser.Builder builder = new ExtProjectDataParser.Builder();
		switch (type) {
		case MANIFEST_FILE:
			if (proj.getArrayType() != ARRAY.ILLUMINA) {
				proj.getLog().reportTimeError("Array type was set to " + proj.getArrayType() + " and this file is for " + ARRAY.ILLUMINA);
				return null;
			}
			builder.separator(",");
			builder.dataKeyColumnName("Name");
			builder.stringDataTitles(new String[] { "AlleleA_ProbeSeq" });
			builder.headerFlags(new String[] { "Name", "AlleleA_ProbeSeq" });
			break;
		case STRAND_REPORT:
			if (proj.getArrayType() != ARRAY.ILLUMINA) {
				proj.getLog().reportTimeError("Array type was set to " + proj.getArrayType() + " and this file is for " + ARRAY.ILLUMINA);
				return null;

			}
			builder.dataKeyColumnName("SNP_Name");
			builder.commentString("##");
			builder.stringDataTitles(new String[] { "Design_Seq" });
			break;
		case AFFY_ANNOT:
			if (proj.getArrayType() != ARRAY.AFFY_GW6 && proj.getArrayType() != ARRAY.AFFY_GW6_CN) {
				proj.getLog().reportTimeError("Array type was set to " + proj.getArrayType() + " and this file is for " + ARRAY.AFFY_GW6 + " or " + ARRAY.AFFY_GW6_CN);
				return null;
			}
			builder.separator("\t");
			builder.dataKeyColumnName("PROBESET_ID");
			builder.stringDataTitles(new String[] { "PROBE_SEQUENCE" });
			builder.concatMultipleStringEntries(true);

		default:
			break;

		}
		builder.sampleBased(false);
		builder.treatAllNumeric(false);
		builder.requireAll(false);
		builder.verbose(false);
		try {
			int seqLength = proj.ARRAY_TYPE.getValue().getProbeLength();
			ExtProjectDataParser parser = builder.build(proj, strandReportFile);
			parser.determineIndicesFromTitles();
			parser.loadData();
			ArrayList<FastaEntry> entries = new ArrayList<FastaEntry>(Array.booleanArraySum(parser.getDataPresent()));
			// SequenceLookup sequenceLookup = new SequenceLookup(proj.getLog());
			for (int i = 0; i < parser.getDataPresent().length; i++) {
				if (parser.getDataPresent()[i]) {

					String seq = null;

					if (type != FILE_SEQUENCE_TYPE.AFFY_ANNOT) {
						seq = parser.getStringData()[0][i];
						if (seq.length() != seqLength) {
							proj.getLog().reportTimeError("Sequence " + seq + " did not have length " + proj.ARRAY_TYPE.getValue().getProbeLength());
							return null;
						}
						entries.add(new FastaEntry(parser.getDataToLoad()[i], seq));
					} else {
						String[] tmpSeq = Array.unique(parser.getStringData()[0][i].split("\t"));
						if (tmpSeq.length != 2) {
							proj.getLog().reportTimeError("Marker " + parser.getDataToLoad()[i] + " did not have 2 unique probe designs");
							proj.getLog().reportTimeError("found the following " + parser.getDataToLoad()[i] + "\t" + Array.toStr(tmpSeq));
							return null;
						} else {
							if (tmpSeq[0].length() != seqLength || tmpSeq[1].length() != seqLength) {
								proj.getLog().reportTimeError("Sequence " + tmpSeq[0] + " or " + tmpSeq[1] + "  did not have length " + proj.ARRAY_TYPE.getValue().getProbeLength());
								return null;
							}
							entries.add(new FastaEntry(parser.getDataToLoad()[i]+"_A", tmpSeq[0]));
							entries.add(new FastaEntry(parser.getDataToLoad()[i]+"_B", tmpSeq[1]));
						}
					}
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
		// String fastaDb = "hg19_canonical.fa";
		String fileSeq = "HumanExome-12-v1-0-B.csv";
		int numThreads = 4;
		int blastWordSize = 30;
		int reportWordSize = 40;
		boolean report = false;
		FILE_SEQUENCE_TYPE fSequence_TYPE = FILE_SEQUENCE_TYPE.MANIFEST_FILE;
		String usage = "\n" + "cnv.qc.MarkerBlast requires 3 arguments\n";
		usage += "   (1) Project file name (i.e. proj=" + filename + " (default))\n" + "";
		usage += "   (2) full path to an Illumina manifest file  (i.e. fileSeq=" + fileSeq + " (default))\n" + "";
		usage += "   (3) number of threads to use  (i.e. " + PSF.Ext.NUM_THREADS_COMMAND + numThreads + " (default))\n" + "";
		usage += "   (4) word size for initial match in blast db  (i.e. blastWordSize=" + blastWordSize + " (default))\n" + "";
		usage += "   (5) number of base pairs with an exact match to report (i.e. reportWordSize=" + reportWordSize + " (default))\n" + "";
		usage += "   (6) report results to temporary files (i.e. -report (not the default))\n" + "";
		usage += "   (7) sequence file type  (i.e. seqType=" + fSequence_TYPE + " (default))\n" + ", Options are:\n ";
		for (int i = 0; i < FILE_SEQUENCE_TYPE.values().length; i++) {
			usage += FILE_SEQUENCE_TYPE.values()[i] + "\n";
		}

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = ext.parseStringArg(args[i], "");
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
			} else if (args[i].startsWith("seqType=")) {
				fSequence_TYPE = FILE_SEQUENCE_TYPE.valueOf(ext.parseStringArg(args[i], FILE_SEQUENCE_TYPE.MANIFEST_FILE.toString()));
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
			blastEm(proj, fileSeq, fSequence_TYPE, blastWordSize, reportWordSize, numThreads, report);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}

// String alts = null;
// for (int j = 0; j < tmpSeq[0].length(); j++) {
// if (tmpSeq[0].charAt(j) != tmpSeq[1].charAt(j)) {
// if (alts != null) {
// proj.getLog().reportTimeError("found multiple probe differences for " + parser.getDataToLoad()[i] + "\t" + Array.toStr(tmpSeq));
// } else {
// alts = sequenceLookup.getNucForCombo("" + tmpSeq[0].charAt(j) + tmpSeq[1].charAt(j));
// seq = tmpSeq[0].substring(0, j) + alts + tmpSeq[0].substring(j + 1, tmpSeq[0].length());
// System.out.println(alts + "\t" + parser.getDataToLoad()[i] + "\t" + tmpSeq[0]);
// System.out.println(alts + "\t" + parser.getDataToLoad()[i] + "\t" + tmpSeq[1]);
// System.out.println(alts + "\t" + parser.getDataToLoad()[i] + "\t" + seq + "\t" + tmpSeq[0].length() + "\t" + seq.length());
// // }
// // }
// // }
// private static boolean isPerfectlyGoodMatch(Blast.BlastResultsSummary bSummary, byte chr) {
// if (bSummary.getNumPerfectMatches() != 1) {
// return false;
// }
//
// return true;
//
// }
