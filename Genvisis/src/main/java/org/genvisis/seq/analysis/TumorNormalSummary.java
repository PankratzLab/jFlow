package org.genvisis.seq.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.Set;

import javax.jms.IllegalStateException;

import org.genvisis.CLI;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.manage.VCFOps.HEADER_COPY_TYPE;
import org.genvisis.seq.manage.VCFOps.VcfPopulation;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.RETRIEVE_TYPE;
import org.genvisis.seq.manage.VCOps;
import org.genvisis.seq.manage.VCOps.VC_SUBSET_TYPE;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;

/**
 * @author lane0212 Quick summary of a region geared towards tumor normal
 */
public class TumorNormalSummary {


	private static final String[] BASE_OUT = new String[] {	"CHROM", "POS", "ID", "REF", "FULL_ALT",
																													"ALT", "HIGH||MODERATE||LOW", "TN_PAIR",
																													"NORMAL_SAMPLE", "TUMOR_SAMPLE",
																													"NORMAL_GENOTYPE", "TUMOR_GENOTYPE",
																													"NORMAL_ALLELES", "TUMOR_ALLELES",
																													"NORMAL_AD", "TUMOR_AD", "NORMAL_GQ",
																													"TUMOR_GQ", "NORMAL_HAS_ALT",
																													"TUMOR_HAS_ALT", "MIN_GQ", "TN_MATCH"};

	private TumorNormalSummary() {

	}

	private static void run(String vcf, String vpopFile, String outputDir, Segment seg, String name,
													int buffer, Logger log) {
		VcfPopulation vpop = VcfPopulation.load(vpopFile, POPULATION_TYPE.TUMOR_NORMAL, log);

		vpop.report();
		Set<String> all = new HashSet<String>();
		all.addAll(vpop.getSuperPop().get(VcfPopulation.TUMOR));
		all.addAll(vpop.getSuperPop().get(VcfPopulation.NORMAL));

		log.reportTimeInfo("Detected " + vpop.getSubPop().keySet().size() + " Tumor normal pairs");
		String outRoot = outputDir + ext.rootOf(vpop.getFileName()) + "_" + name;
		String outVCF = outRoot + ".vcf.gz";
		String outSummary = outRoot + ".summary.txt";

		VCFFileReader reader = new VCFFileReader(new File(vcf), true);
		VariantContextWriter writer = VCFOps.initWriter(outVCF, VCFOps.DEFUALT_WRITER_OPTIONS,
																										reader.getFileHeader().getSequenceDictionary());
		VCFOps.copyHeader(reader, writer, null, HEADER_COPY_TYPE.FULL_COPY, log);
		Segment segBuffer = seg.getBufferedSegment(buffer);
		CloseableIterator<VariantContext> iter = reader.query(Positions.getChromosomeUCSC(segBuffer.getChr(),
																																											true),
																													segBuffer.getStart(),
																													segBuffer.getStop());
		try {
			PrintWriter writerSummary = new PrintWriter(new FileWriter(outSummary));
			String[][] annos = VCFOps.getAnnotationKeys(vcf, log);
			writerSummary.println(ArrayUtils.toStr(BASE_OUT) + "\t" + ArrayUtils.toStr(annos[0]));

			while (iter.hasNext()) {
				VariantContext vcFull = iter.next();
				VariantContext vc = VCOps.getSubset(vcFull, all, VC_SUBSET_TYPE.SUBSET_STRICT);
				if (!vc.isMonomorphicInSamples()) {
					writer.add(vc);
					for (String tnPair : vpop.getSubPop().keySet()) {
						Set<String> samps = vpop.getSubPop().get(tnPair);

							String tumor = null;
							String normal = null;
							for (String samp : samps) {
								if (vpop.getPopulationForInd(	samp,
																							RETRIEVE_TYPE.SUPER)[0].equals(VcfPopulation.TUMOR)) {
									tumor = samp;
								} else if (vpop.getPopulationForInd(samp,
																										RETRIEVE_TYPE.SUPER)[0].equals(VcfPopulation.NORMAL)) {
									normal = samp;
								} else {
									writerSummary.close();
									throw new IllegalStateException("Unknown types");
								}
							}

							Genotype gTumor = vc.getGenotype(tumor);
							Genotype gNormal = vc.getGenotype(normal);
							StringBuilder builder = new StringBuilder();
							builder.append(vc.getContig());
							builder.append("\t" + vc.getStart());
							builder.append("\t" + vc.getID());
							builder.append("\t" + vc.getReference().getDisplayString());
							builder.append("\t" + vcFull.getAlternateAlleles());
							builder.append("\t" + vc.getAlternateAlleles());
							String impact = VCOps.getSNP_EFFImpact(vc);
							builder.append("\t" + (impact.equals("HIGH")|| impact.equals("MODERATE")
																			|| impact.equals("LOW")));
							builder.append("\t" + tnPair);
							builder.append("\t" + normal);
							builder.append("\t" + tumor);
							builder.append("\t" + gNormal.toString());
							builder.append("\t" + gTumor.toString());
							builder.append("\t" + gNormal.getGenotypeString());
							builder.append("\t" + gTumor.getGenotypeString());
							builder.append("\t"
															+ ArrayUtils.toStr(ArrayUtils.toStringArray(gNormal.getAD()), ","));
							builder.append("\t"
															+ ArrayUtils.toStr(ArrayUtils.toStringArray(gTumor.getAD()), ","));
							builder.append("\t" + gNormal.getGQ());
							builder.append("\t" + gTumor.getGQ());
							builder.append("\t" + (gNormal.isCalled() && !gNormal.isHomRef()));
							builder.append("\t" + (gTumor.isCalled() && !gTumor.isHomRef()));

							builder.append("\t" + Math.min(gNormal.getGQ(), gTumor.getGQ()));
							builder.append("\t" + gTumor.sameGenotype(gNormal));
							builder.append("\t" + ArrayUtils.toStr(VCOps.getAnnotationsFor(annos[0], vc, ".")));
							writerSummary.println(builder.toString());
					}
				}
			}

			writerSummary.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + outSummary);
			log.reportException(e);
		}
		writer.close();
		reader.close();

	}

	public static void main(String[] args) {
		CLI c = new CLI(TumorNormalSummary.class);

		c.addArgWithDefault("vcf", "vcf to tally", "a.vcf");
		c.addArgWithDefault("vpop", "vpop to use", "a.vpop");
		c.addArgWithDefault("segment", "UCSC segment", "chr18:20714428-20840534");
		c.addArgWithDefault("name", "typically gene name", "CABLES1");
		c.addArgWithDefault("outDir", "output directory", "out");


		c.parseWithExit(args);

		String vcf = c.get("vcf");
		String vpop = c.get("vpop");
		Segment seg = new Segment(c.get("segment"));
		String name = c.get("name");
		String outDir = c.get("outDir");

		int buffer = 300;
		run(vcf, vpop, outDir, seg, name, buffer, new Logger());

	}
}

