package org.genvisis.one.JL;

import java.io.File;
import java.util.HashSet;
import java.util.concurrent.Callable;

import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.ext;
import org.genvisis.common.WorkerTrain.AbstractProducer;
import org.genvisis.seq.manage.BamOps;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.manage.VCOps;
import org.genvisis.seq.manage.VCFOps.HEADER_COPY_TYPE;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;

public class VerifyBamID {

	private static final String VerifyBamIDLocation = "/panfs/roc/groups/5/pankrat2/public/bin/verifyBamID_1.1.0/verifyBamID/bin/verifyBamID";

	public static void verifyBamID(String bamFile, String vcf) {
		String dir = ext.parseDirectoryOfFile(vcf) + "verifyBamID/";
		new File(dir).mkdirs();
		Logger log = new Logger();
		String sample = BamOps.getSampleName(bamFile);
		HashSet<String> sub = new HashSet<String>();
		String verfifyVCF = ext.rootOf(ext.addToRoot(vcf, ".verify")) + ".gz";

		String output = ext.rootOf(bamFile) + "_verify";
		verfifyVCF = dir + verfifyVCF;
		sub.add(sample);
		// output+=
		if (!Files.exists(verfifyVCF)) {
			VCFFileReader reader = new VCFFileReader(new File(vcf), true);
			VariantContextWriter writer = VCFOps.initWriter(verfifyVCF, null, reader.getFileHeader().getSequenceDictionary());
			VCFOps.copyHeader(reader, writer, null, HEADER_COPY_TYPE.FULL_COPY, log);
			int index = 0;
			for (VariantContext vc : reader) {
				index++;
				if (index % 100000 == 0) {
					log.reportTimeInfo(index + " variants");
				}

				if (vc.isBiallelic() && !vc.isIndel() && VCOps.getSegment(vc).getChr() < 23) {
					writer.add(vc);
				}
			}
			reader.close();
			writer.close();
		}

		String[] command = new String[] { VerifyBamIDLocation, "--vcf", verfifyVCF, "--bam", bamFile, "--out", output, "--verbose", "--ignoreRG", "--chip-none" };
		CmdLine.runCommandWithFileChecks(command, "", null, null, false, false, false, log);
	}

	public static class VerifyProducer extends AbstractProducer<Boolean> {
		private String[] bamFiles;
		private String vcf;
		private int index;

		public VerifyProducer(String[] bamFiles, String vcf) {
			super();
			this.bamFiles = bamFiles;
			this.vcf = vcf;
		}

		@Override
		public boolean hasNext() {

			return index < bamFiles.length;
		}

		@Override
		public Callable<Boolean> next() {
			final String bamFile = bamFiles[index];
			Callable<Boolean> callable = new Callable<Boolean>() {

				@Override
				public Boolean call() throws Exception {
					verifyBamID(bamFile, vcf);
					// TODO Auto-generated method stub
					return true;
				}
			};
			index++;
			// TODO Auto-generated method stub
			return callable;
		}
	}

	public static void runVerifies(String bamDir, String vcf, int nThreads) {
		String[] bams = Files.listFullPaths(bamDir, ".bam", false);

		//VerifyProducer producer = new VerifyProducer(bams, vcf);
		//WorkerTrain<Boolean> train = new WorkerTrain<Boolean>(producer, nThreads, 10, new Logger());
		for (int i = 0; i < bams.length; i++) {
			verifyBamID(bams[i], vcf);
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "VerifyBamID.dat";
		String bamDir = "/bams/";
		int numThreads = 8;
		//String logfile = null;
		//Logger log;

		String usage = "\n" + "one.JL.VerifyBamID requires 0-1 arguments\n";
		usage += "   (1) vcf file name (i.e. file=" + filename + " (default))\n" + "";
		usage += "   (2) bam Directory (i.e. bams=" + filename + " (default))\n" + "";
		usage += PSF.Ext.getNumThreadsCommand(3, numThreads);

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("bams=")) {
				bamDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
				numThreads = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				//logfile = args[i].split("=")[1];
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
			//log = new Logger(logfile);
			runVerifies(bamDir, filename, numThreads);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
