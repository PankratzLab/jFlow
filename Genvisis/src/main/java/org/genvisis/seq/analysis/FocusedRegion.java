/**
 * 
 */
package org.genvisis.seq.analysis;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;

import org.genvisis.CLI;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.WorkerTrain.AbstractProducer;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.BamOps;
import org.genvisis.seq.manage.BamPile;
import org.genvisis.seq.manage.BamPileUp;
import org.genvisis.seq.manage.BamPileUp.PILE_TYPE;
import org.genvisis.seq.manage.ReferenceGenome;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.qc.FilterNGS;
import org.genvisis.seq.qc.FilterNGS.SAM_FILTER_TYPE;

/**
 * @author Kitty Widget to perform a focused bp analysis of a region
 */
public class FocusedRegion {
	private FocusedRegion() {

	}


	private static void focus(String bamFile, String ref, String outDir, String vcf, Segment seg,
														String variantSet, int numthreads) {
		new File(outDir).mkdirs();

		String output = outDir	+ ext.replaceWithLinuxSafeCharacters(seg.getUCSClocation(), true)
										+ "_focus.txt";

		Logger log = new Logger(outDir + "focus.log");
		String[] bams = HashVec.loadFileToStringArray(bamFile, false, new int[] {0}, true);
		log.reportTimeInfo(bams.length + " bam files found");
		PrintWriter writer = Files.getAppropriateWriter(output);
		FocusProducer producer = new FocusProducer(seg, new ReferenceGenome(ref, log), bams, log);
		WorkerTrain<FocusResults> train = new WorkerTrain<FocusResults>(producer, numthreads, 10, log);

		writer.println("SAMPLE\t" + BamPile.getBampPileHeader());
		while (train.hasNext()) {

			FocusResults results = train.next();
			for (BamPile bamPile : results.piles) {
				writer.println(results.sample + "\t" + bamPile.getOuput(log));
			}
			log.reportTimeInfo("Finished " + results.sample);
		}
		writer.close();
		String segFile = ext.addToRoot(output, ".segs");

		Files.write(seg.getChromosomeUCSC() + "\t" + seg.getStart() + "\t" + seg.getStop(), segFile);
		String subOut = outDir + "extract/";
		new File(subOut).mkdirs();
		VCFOps.extractSegments(	vcf, segFile, 300, bamFile, subOut, false, true, true, true,
														new String[] {variantSet}, numthreads, log);
	}

	private static class FocusResults {
		private List<BamPile> piles;
		private String sample;

		private FocusResults(List<BamPile> piles, String sample) {
			super();
			this.piles = piles;
			this.sample = sample;
		}


	}

	/**
	 * 
	 *
	 */
	private static class FocusProducer extends AbstractProducer<FocusResults> {
		private final Segment seg;
		private final ReferenceGenome referenceGenome;
		private final String[] bamFiles;
		private int index;
		private final Logger log;

		private FocusProducer(Segment seg, ReferenceGenome referenceGenome, String[] bamFiles,
													Logger log) {
			super();
			this.seg = seg;
			this.referenceGenome = referenceGenome;
			this.bamFiles = bamFiles;
			this.log = log;
			index = 0;
		}

		@Override
		public boolean hasNext() {
			return index < bamFiles.length;
		}

		@Override
		public Callable<FocusResults> next() {
			final int current = index;
			Callable<FocusResults> worker = new Callable<FocusResults>() {

				@Override
				public FocusResults call() throws Exception {
					BamPileUp pileUp = new BamPileUp(	bamFiles[current], referenceGenome, 1, new FilterNGS(),
																						new Segment[] {seg}, PILE_TYPE.REGULAR,
																						SAM_FILTER_TYPE.COPY_NUMBER, true, log);
					ArrayList<BamPile> pile = new ArrayList<BamPile>();
					while (pileUp.hasNext()) {
						pile.add(pileUp.next());
					}
					return new FocusResults(pile, BamOps.getSampleName(bamFiles[current], log));
				}
			};
			index++;

			return worker;
		}
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		String bams = "bams.txt";
		String outDir = "out/";
		String segment = "chr18:20715728-20716571";
		String ref = "hg19.fa";
		String vcf = "a.vcf";
		String variantSet = ".variant";

		int numthreads = 24;
		CLI c = new CLI(FocusedRegion.class);
		c.addArgWithDefault(CLI.ARG_OUTDIR, CLI.DESC_OUTDIR, outDir);
		c.addArgWithDefault("vcf", "vcf for extraction", vcf);

		c.addArgWithDefault("seg", "segment to focus", segment);
		c.addArgWithDefault("bams", "file of bams", bams);
		c.addArgWithDefault("variantSet",
												"variant set mapping bam samples to vcf samples, set to \"\" if none",
												variantSet);

		c.addArgWithDefault("ref", "reference genome", ref);
		c.addArgWithDefault("threads", "number of threads", Integer.toString(numthreads));
		c.parseWithExit(args);
		focus(c.get("bams"), c.get("ref"), c.get(CLI.ARG_OUTDIR), c.get("vcf"),
					new Segment(c.get("seg")), c.get("variantSet"), c.getI("threads"));

	}
}
