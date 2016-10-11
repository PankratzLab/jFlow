package org.genvisis.seq.manage;

import java.io.File;
import java.io.IOException;

import org.genvisis.common.Files;
import org.genvisis.common.Logger;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

/**
 * @author lane0212 just read bams and vcfs to use the htsjdk error checking to validate. Can be
 *         used after transferring files
 *
 */
public class VerifyNGS {

	public static void verify(String bamDir, String vcfDir, Logger log) {
		String[] bams = Files.list(bamDir, null, ".bam", true, false, true);
		log.reportTimeInfo("found " + bams.length + " bams");
		for (String bam : bams) {
			log.reportTimeInfo("Validating " + bam);
			SamReader reader = BamOps.getDefaultReader(bam, ValidationStringency.STRICT);
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
		String[] vcfs = Files.list(vcfDir, null, ".vcf", true, false, true);
		for (String vcf : vcfs) {
			log.reportTimeInfo("Validating " + vcf);
			VCFFileReader reader = new VCFFileReader(new File(vcf), true);
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

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("bams=")) {
				bamDir = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("vcfs=")) {
				vcfDir = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("log=")) {
				logfile = arg.split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + arg);
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
