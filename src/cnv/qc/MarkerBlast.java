package cnv.qc;

import java.io.FileNotFoundException;
import java.util.ArrayList;

import common.Array;
import common.Files;
import common.PSF;
import common.WorkerHive;
import common.ext;
import cnv.filesys.Project;
import cnv.manage.ExtProjectDataParser;
import seq.analysis.Blast;
import seq.analysis.Blast.BlastResults;
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

	private static void blastEm(Project proj, String fileSeq, FILE_SEQUENCE_TYPE type, String fastaDb, int numThreads) {
		Blast blast = new Blast(fastaDb, 50, 110, proj.getLog(), true, true);
		WorkerHive<BlastResults[]> hive = new WorkerHive<Blast.BlastResults[]>(numThreads, 10, proj.getLog());
		FastaEntry[] fastaEntries = getMarkerFastaEntries(proj, fileSeq, type);
		ArrayList<FastaEntry[]> splits = Array.splitUpArray(fastaEntries, numThreads, proj.getLog());
		BlastWorker[] workers = new BlastWorker[splits.size()];
		String[] tmps = new String[workers.length];
		if (fastaEntries != null && fastaEntries.length > 0) {
			for (int i = 0; i < splits.size(); i++) {
				tmps[i] = ext.rootOf(fileSeq, false) + ".tmp" + i;
				workers[i] = new BlastWorker(blast, splits.get(i), tmps[i]);
			}
		}
		hive.addCallables(workers);
		hive.setReportEvery(1);
		hive.execute(true);
		Files.cat(tmps, ext.rootOf(fileSeq, false) + ".blasted", new int[] {}, proj.getLog());
	}

	private static FastaEntry[] getMarkerFastaEntries(Project proj, String strandReportFile, FILE_SEQUENCE_TYPE type) {
		FastaEntry[] fastaEntries = new FastaEntry[0];
		ExtProjectDataParser.Builder builder = new ExtProjectDataParser.Builder();

		switch (type) {
		case MANIFEST_FILE:
			builder.separator(",");
			builder.dataKeyColumnName("Name");
			builder.stringDataTitles(new String[] { "SourceSeq" });
			builder.headerFlags(new String[] { "Name", "SourceSeq" });
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
					int leftBracket = seq.indexOf("[");
					int rightBrackt = seq.indexOf("]");
					int slash = seq.indexOf("/");
					String A = seq.substring(leftBracket + 1, slash);
					String B = seq.substring(slash + 1, rightBrackt);
					String left = seq.substring(0, leftBracket);
					String right = seq.substring(rightBrackt + 1);
					String seqA = left + A + right;
					String seqB = left + B + right;
					if (!A.equals("-")) {
						entries.add(new FastaEntry(parser.getDataToLoad()[i] + "_A", seqA));
					}
					if (!B.equals("-")) {
						entries.add(new FastaEntry(parser.getDataToLoad()[i] + "_B", seqB));
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
		String fastaDb = "hg19_canonical.fa";
		String fileSeq = "HumanExome-12-v1-0-B.csv";
		int numThreads = 4;
		String usage = "\n" + "cnv.qc.MarkerBlast requires 3 arguments\n";
		usage += "   (1) Project file name (i.e. proj=" + filename + " (default))\n" + "";
		usage += "   (2) full path to a fasta file to search (i.e. fastaDb=" + fastaDb + " (default))\n" + "";
		usage += "   (3) full path to an Illumina manifest file  (i.e. fileSeq=" + fileSeq + " (default))\n" + "";
		usage += "   (4) number of threads to use  (i.e. " + PSF.Ext.NUM_THREADS_COMMAND + numThreads + " (default))\n" + "";

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
			blastEm(proj, fileSeq, FILE_SEQUENCE_TYPE.MANIFEST_FILE, fastaDb, numThreads);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
