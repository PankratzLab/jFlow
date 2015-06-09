package seq.manage;

import seq.analysis.GATK;
import common.Array;
import common.Files;
import common.Logger;
import common.ext;

/**
 * @author lane0212 Merge a list of vcfs using GATK;
 */
public class VCFMerge {

	public static void merge(String[] vcfs, String mergeOut, String gatkLoc, String referenceGenomeFasta, Logger log) {
		GATK gatk = new GATK(gatkLoc, referenceGenomeFasta, true, false, log);
		if (vcfs != null) {
			if (mergeOut != null) {
				if (!Files.exists("", vcfs)) {
					boolean merged = gatk.mergeVCFs(vcfs, mergeOut, log);
					if (merged) {
						log.reportTimeInfo("Merged " + Array.toStr(vcfs, "\n") + " to " + mergeOut);
						VCFOps.extractSamps(mergeOut, log);
					} else {
						log.reportTimeError("Could not merge " + Array.toStr(vcfs, "\n") + " to " + mergeOut);
					}
				} else {
					log.reportFilesNotFound(vcfs);
				}
			}else{
				log.reportTimeError("Must provide an output file for the merged results");
			}
		} else {
			log.reportFilesNotFound(vcfs);
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String[] vcfs = null;
		String mergeOut = null;
		String referenceGenomeFasta = ReferenceGenome.DEFAULT_REFERENCE;
		String gatk = GATK.DEFAULT_GATK;

		String usage = "\n" + "seq.manage.VCFMerge requires 0-1 arguments\n";
		usage += "   (1) comma delimited full paths to vcf files (i.e. vcfs= (no default))\n" + "";
		usage += "   (2) full path to the merged output (i.e. mergeOut= (no default))\n" + "";
		usage += "   (3) full path to the gatk directory (i.e. gatk=" + gatk + " (default))\n" + "";
		usage += "   (4) full path to the reference genome  (i.e. ref=" + gatk + " (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("vcfs=")) {
				vcfs = args[i].split("=")[1].split(",");
				numArgs--;
			} else if (args[i].startsWith("mergeOut=")) {
				mergeOut = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("gatk=")) {
				gatk = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("ref=")) {
				referenceGenomeFasta = args[i].split("=")[1];
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
			merge(vcfs, mergeOut, gatk, referenceGenomeFasta, new Logger(vcfs == null ? null : ext.addToRoot(vcfs[0], "mergeLog.log")));
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
