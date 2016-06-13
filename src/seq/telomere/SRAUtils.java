package telomere;

import java.util.ArrayList;

import common.Logger;
import common.PSF;
import common.ext;

public class SRAUtils {

	private static final String SAM_DUMP = "sam-dump.2.6.3";

	private static boolean dumpSra(String inputSra, String outputBam, Logger log) {
		ArrayList<String> command = new ArrayList<String>();

		command.add(SAM_DUMP);
		command.add(ext.rootOf(inputSra));
		command.add(ext.rootOf(inputSra));

		return false;
	}

	// sam-dump.2.6.3 SRR1737697 |samtools view -bS -

	private static void run(String sraDir, String outDir, int threads) {

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "SRAUtils.dat";
		String sraDir = "/scratch.global/lanej/aric_raw/sra/";
		String outDir = "/scratch.global/lanej/aric_raw/";
		int threads = 24;

		String usage = "\n" +
				"telomere.SRAUtils requires 0-1 arguments\n" +
				"   (1) SRA directory (i.e. sraDir=" + sraDir + " (default))\n" +
				"   (2) out directory (i.e. outDir=" + outDir + " (default))\n" +
				PSF.Ext.getNumThreadsCommand(3, threads) +
				"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("sraDir")) {
				sraDir = args[i].split("=")[1];
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
			run(sraDir, outDir, threads);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
