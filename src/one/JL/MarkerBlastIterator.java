package one.JL;

import cnv.filesys.Project;
import cnv.qc.MarkerBlast;
import cnv.qc.MarkerBlast.FILE_SEQUENCE_TYPE;
import cnv.qc.MarkerBlast.MarkerBlastResult;
import common.PSF;
import common.ext;

public class MarkerBlastIterator {

	
	
	public static void blastIter(Project proj, String fileSeq, FILE_SEQUENCE_TYPE type, int numThreads, boolean reportToTmp) {
		
		
		
		int[] blastWordSizes = new int[] { 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25 };
		int[] reportWordSizes = new int[] { 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25 };
		MarkerBlastResult[] results = new MarkerBlastResult[blastWordSizes.length * reportWordSizes.length];
		int index = 0;
		for (int i = 0; i < blastWordSizes.length; i++) {
			for (int j = 0; j < reportWordSizes.length; j++) {
				results[index] = MarkerBlast.blastEm(proj, fileSeq, type, blastWordSizes[i], reportWordSizes[j], numThreads, reportToTmp);
			}
		}
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		//String fastaDb = "hg19_canonical.fa";
		String fileSeq = "HumanExome-12-v1-0-B.csv";
		int numThreads = 4;
		int blastWordSize = 30;
		int reportWordSize = 40;
		boolean report = false;
		FILE_SEQUENCE_TYPE fSequence_TYPE = FILE_SEQUENCE_TYPE.AFFY_ANNOT;
		String usage = "\n" + "cnv.qc.MarkerBlast requires 3 arguments\n";
		usage += "   (1) Project file name (i.e. proj=" + filename + " (default))\n" + "";
		usage += "   (2) full path to an Illumina manifest file  (i.e. fileSeq=" + fileSeq + " (default))\n" + "";
		usage += "   (3) number of threads to use  (i.e. " + PSF.Ext.NUM_THREADS_COMMAND + numThreads + " (default))\n" + "";
		usage += "   (4) word size for initial match in blast db  (i.e. blastWordSize=" + blastWordSize + " (default))\n" + "";
		usage += "   (5) number of base pairs with an exact match to report (i.e. reportWordSize=" + reportWordSize + " (default))\n" + "";
		usage += "   (6) report results to temporary files (i.e. -report (not the default))\n" + "";
		usage += "   (7) sequence file type  (i.e. seqType=" + fSequence_TYPE + " (default))\n" + ", Options are:\n ";
		

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
			blastIter(proj, fileSeq, fSequence_TYPE, numThreads, report);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
