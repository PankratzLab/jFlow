package seq.telomere;

import java.io.File;
import java.util.ArrayList;

import seq.telomere.SRAUtils.SRAConversionResult;
import common.Logger;
import common.PSF;
import common.ext;

public class TelSeq {

	private static void telSeqIt(String inputBam, String outputDir, String optionalBed, Logger log) {
		String[] input = new String[] { inputBam };
		ArrayList<String> command = new ArrayList<String>();
		command.add("telseq");

	}

	public static void run(String sraDir, String outDir, String optionalBed, int threads) {
		ArrayList<SRAConversionResult> conv = SRAUtils.run(sraDir, outDir, threads);
		String telseqDir = outDir + "telseq/";
		new File(telseqDir).mkdirs();
		Logger log = new Logger(telseqDir + ".telseq.log");
		log.reportTimeInfo("Assuming telseq is on system path");

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String sraDir = "/scratch.global/lanej/aric_raw/sra/";
		String outDir = "/scratch.global/lanej/aric_raw/";
		String captureBed = "/home/pankrat2/public/bin/ref/VCRome_2_1_hg19_capture_targets.bed";
		int threads = 24;

		String usage = "\n" +
				"telomere.SRAUtils requires 0-1 arguments\n" +
				"   (1) SRA directory (i.e. sraDir=" + sraDir + " (default))\n" +
				"   (2) out directory (i.e. outDir=" + outDir + " (default))\n" +
				"   (3) capture bed (i.e. bed=" + outDir + " (default))\n" +

				PSF.Ext.getNumThreadsCommand(3, threads) +
				"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("sraDir")) {
				sraDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("bed=")) {
				captureBed = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("outDir")) {
				outDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
				threads = ext.parseIntArg(args[i]);
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
			run(sraDir, outDir, captureBed, threads);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
