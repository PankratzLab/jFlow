package seq.analysis;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Set;

import filesys.Segment;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import seq.manage.BamExtractor;
import seq.manage.BamExtractor.WorkerExtractor;
import seq.manage.BamOps;
import seq.manage.VCFOps;
import seq.manage.VCOps;
import common.Files;
import common.HashVec;
import common.Logger;
import common.PSF;
import common.WorkerHive;
import common.ext;

/**
 * Indel qc
 *
 */
public class Indelathon {

	@SuppressWarnings("unchecked")
	private static void extractIndels(String vcf, String bamFilesFile, String outDir, Set<String> variantSets, int buffer, int numThreads) {
		new File(outDir).mkdirs();
		Logger log = new Logger(outDir + "indel.log");
		String[] bams = HashVec.loadFileToStringArray(bamFilesFile, false, new int[] { 0 }, true);
		Hashtable<String, String> matchedSamps = BamOps.matchToVcfSamplesToBamFiles(VCFOps.getSamplesInFile(vcf), variantSets, bams, numThreads, log);

		String outIndelVCF = outDir + VCFOps.getAppropriateRoot(vcf, true) + ".indels.vcf.gz";
		String outSegSer = outDir + VCFOps.getAppropriateRoot(vcf, true) + ".indels.seg.ser";
		Hashtable<String, ArrayList<Segment>> sampSegs = new Hashtable<String, ArrayList<Segment>>();

		if (!Files.exists(outIndelVCF) || !Files.exists(outSegSer)) {
			VCFFileReader reader = new VCFFileReader(vcf, true);
			VariantContextWriter writer = VCFOps.initWriterWithHeader(reader, outIndelVCF, VCFOps.DEFUALT_WRITER_OPTIONS, log);
			int numTotal = 0;
			int numIndels = 0;
			for (VariantContext vc : reader) {
				numTotal++;
				if (vc.isIndel()) {
					GenotypesContext gc = vc.getGenotypes();
					for (Genotype g : gc) {
						if (!g.isNoCall() && !g.isHomRef()) {
							if (!sampSegs.containsKey(g.getSampleName())) {
								sampSegs.put(g.getSampleName(), new ArrayList<Segment>());
							}
							sampSegs.get(g.getSampleName()).add(VCOps.getSegment(vc));
						}
					}
					numIndels++;
					writer.add(vc);
				}
			}
			log.reportTimeInfo("Found " + numIndels + " indels for " + numTotal + " variants total");
			writer.close();
			reader.close();
			Files.writeSerial(sampSegs, outSegSer, true);
		}
		sampSegs = (Hashtable<String, ArrayList<Segment>>) Files.readSerial(outSegSer, false, log, false, true);
		WorkerHive<BamExtractor> hive = new WorkerHive<BamExtractor>(numThreads, 10, log);
		ArrayList<String> indelBams = new ArrayList<String>();
		String indelBamDir = outDir + "indel_bams/";
		new File(indelBamDir).mkdirs();
		for (String samp : matchedSamps.keySet()) {
			String outbam = indelBamDir + samp + "_indels_" + buffer + "bp.bam";
			indelBams.add(outbam);
			Segment[] segmentsToExtract = new Segment[] {};
			if (sampSegs.containsKey(samp)) {
				segmentsToExtract = sampSegs.get(samp).toArray(new Segment[sampSegs.get(samp).size()]);
			}
			WorkerExtractor extractor = new WorkerExtractor(segmentsToExtract, matchedSamps.get(samp), false, false, buffer, outbam, log);
			hive.addCallable(extractor);
		}
		hive.execute(true);
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String vcf = "Indelathon.vcf";
		String outDir = "Indelathon/";
		String bams = null;
		HashSet<String> sets = new HashSet<String>();
		int numThreads = 24;
		int buffer = 100;

		String usage = "\n" +
				"seq.analysis.Indelathon requires 0-1 arguments\n" +
				"   (1) filename (i.e. vcf=" + vcf + " (default))\n" +
				"   (2) outDir (i.e. out=" + outDir + " (default))\n" +
				"   (3) bams (i.e. bams=" + bams + " (default))\n" +
				"   (4) comma delimited variant sets (i.e. sets= (default))\n" +

				PSF.Ext.getNumThreadsCommand(4, numThreads) +

				"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("vcf=")) {
				vcf = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("out=")) {
				outDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("bams=")) {
				bams = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("sets=")) {
				String[] tmp = ext.parseStringArg(args[i], "").split(",");
				for (int j = 0; j < tmp.length; j++) {
					sets.add(tmp[j]);
				}
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

			extractIndels(vcf, bams, outDir, sets, buffer, numThreads);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
