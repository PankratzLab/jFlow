package one.JL;

import java.io.IOException;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import seq.manage.BamOps;
import seq.manage.VCFOps;
import common.Files;
import common.Logger;

/**
 * @author lane0212 just read bams and vcfs to use the htsjdk error checking to validate. Can be used after transferring files
 *
 */
public class VerifyNGS {

	public static void verify(String bamDir, String vcfDir, Logger log) {
		String[] bams = Files.list(bamDir, null, ".bam", true, false, true);
		log.reportTimeInfo("found " + bams.length + " bams");
		for (int i = 0; i < bams.length; i++) {
			log.reportTimeInfo("Validating " + bams[i]);
			SamReader reader = BamOps.getDefaultReader(bams[i], ValidationStringency.STRICT);
			for (SAMRecord samRecord : reader) {
				samRecord.getCigar();//
			}
			try {
				reader.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		String[] vcfs = Files.list(vcfDir, null, ".vcfs", true, false, true);
		for (int i = 0; i < vcfs.length; i++) {
			log.reportTimeInfo("Validating " + vcfs[i]);
			VCFFileReader reader = new VCFFileReader(vcfs[i], true);
			for (VariantContext vc : reader) {
				vc.fullyDecode(reader.getFileHeader(), false);
				vc.getGenotypes();
			}
			reader.close();
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String bamDir = "/bams/";
		String vcfDir = "/vcf/";
		String logfile = null;
		Logger log;

		String usage = "\n" + "one.JL.verifyBam requires 0-1 arguments\n";
		usage += "   (1) directory of bams (i.e. bams=" + bamDir + " (default))\n" + "";
		usage += "   (2) directory of vcfs to validate (i.e. vcfs=" + bamDir + " (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("bams=")) {
				bamDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("vcfs=")) {
				vcfDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = args[i].split("=")[1];
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
			log = new Logger(logfile);
			verify(bamDir, vcfDir, log);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
