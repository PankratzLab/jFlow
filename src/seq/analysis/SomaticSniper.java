package seq.analysis;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Set;
import java.util.concurrent.Callable;

import common.Array;
import common.CmdLine;
import common.Files;
import common.HashVec;
import common.Logger;
import common.WorkerTrain;
import common.ext;
import common.WorkerTrain.Producer;
import seq.manage.BamOps;
import seq.manage.GenotypeOps;
import seq.manage.VCFOps;
import seq.manage.VCFOps.VcfPopulation;
import seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import seq.manage.VCFOps.VcfPopulation.RETRIEVE_TYPE;

/**
 * Wrapper and parser/assembler for Somatic Sniper (https://github.com/genome/somatic-sniper)
 *
 */
public class SomaticSniper {

	private static void run(GATK gatk, SomaticParams params, String bams, String vcfPop, int numThreads, Logger log) {
		VcfPopulation vpop = VcfPopulation.load(vcfPop, POPULATION_TYPE.TUMOR_NORMAL, log);
		String[] bamFiles = null;
		if (Files.isDirectory(bams)) {
			bamFiles = Files.list(bams, ".bam", false);
		} else {
			bamFiles = HashVec.loadFileToStringArray(bams, false, new int[] { 0 }, true);
		}
		log.reportTimeInfo("Detected " + bamFiles.length + " total bam files");

		TNSample[] tnSamples = matchSamples(bamFiles, params.getOutputDir(), vpop, params, log);
		String tnMatch = params.getOutputDir() + "tnMatch.txt";
		ArrayList<String> summaryMatch = new ArrayList<String>();
		summaryMatch.add("NORMAL_BAM\tTUMOR_BAM");
		for (int i = 0; i < tnSamples.length; i++) {
			summaryMatch.add(tnSamples[i].getNormalBam() + "\t" + tnSamples[i].getTumorBam());
		}
		Files.writeList(Array.toStringArray(summaryMatch), tnMatch);
		TNProducer producer = new TNProducer(tnSamples);
		WorkerTrain<TNSample> train = new WorkerTrain<TNSample>(producer, numThreads, 2, log);
		ArrayList<String> finalOuts = new ArrayList<String>();
		while (train.hasNext()) {
			TNSample tmp = train.next();
			finalOuts.add(tmp.getOutputGz());
		}
		String out = params.getOutputDir() + "tn.out.vcf.gz";
		String outFiltRename = params.getOutputDir() + "tn.out.filt.vcf.gz";

		if (gatk.mergeVCFs(Array.toStringArray(finalOuts), out, numThreads, false, log)) {
			VCFFileReader reader = new VCFFileReader(out, true);
			Set<String> samps = new HashSet<String>();
			String[] sampIn = VCFOps.getSamplesInFile(out);
			for (int i = 0; i < sampIn.length; i++) {
				String fix = sampIn[i].replaceAll(".variant.*", "");
				samps.add(fix);
			}

			final VCFHeader outHeader = new VCFHeader(reader.getFileHeader().getMetaDataInInputOrder(), samps);

			VariantContextWriter writer = VCFOps.initWriter(outFiltRename, VCFOps.DEFUALT_WRITER_OPTIONS, reader.getFileHeader().getSequenceDictionary());
			writer.writeHeader(outHeader);
			for (VariantContext vc : reader) {
				VariantContextBuilder builder = new VariantContextBuilder(vc);
				ArrayList<Genotype> renameGeno = new ArrayList<Genotype>();
				for (Genotype g : vc.getGenotypes()) {
					GenotypeBuilder gBuilder = new GenotypeBuilder(g);
					gBuilder.name(g.getSampleName().replaceAll(".variant.*", ""));
					double mapQ = 0;
					double ssc = 0;
					try {
						if (g.hasAnyAttribute("MQ")) {
							mapQ = Double.parseDouble(g.getAnyAttribute("MQ").toString());
						}
						if (g.hasAnyAttribute("SSC")) {
							ssc = Double.parseDouble(g.getAnyAttribute("SSC").toString());
						}
					} catch (NumberFormatException nfe) {

					}
					if (mapQ < 40 || ssc < 40) {
						gBuilder.alleles(GenotypeOps.getNoCall());
					}
					renameGeno.add(gBuilder.make());
				}
				builder.genotypes(renameGeno);
				VariantContext vcFilt = builder.make();
				if (!vcFilt.isMonomorphicInSamples()) {
					writer.add(vcFilt);
				}
			}
			reader.close();
			writer.close();
		}

	}

	private static class TNProducer implements Producer<TNSample> {
		private TNSample[] tnSamples;
		private int index;

		public TNProducer(TNSample[] tnSamples) {
			super();
			this.tnSamples = tnSamples;
			this.index = 0;
		}

		@Override
		public boolean hasNext() {
			return index < tnSamples.length;
		}

		@Override
		public Callable<TNSample> next() {
			TNSample current = tnSamples[index];
			index++;

			return current;
		}

		@Override
		public void shutdown() {
			// TODO Auto-generated method stub

		}

	}

	private static TNSample[] matchSamples(String[] bamFiles, String outputDir, VcfPopulation vpop, SomaticParams somaticParams, Logger log) {
		Hashtable<String, String> all = new Hashtable<String, String>();
		for (int i = 0; i < bamFiles.length; i++) {
			all.put(BamOps.getSampleName(bamFiles[i]), bamFiles[i]);
		}
		ArrayList<String> analysisBams = new ArrayList<String>();
		ArrayList<TNSample> tnSamples = new ArrayList<TNSample>();
		for (String tnPair : vpop.getSubPop().keySet()) {
			Set<String> samps = vpop.getSubPop().get(tnPair);
			String tumor = null;
			String normal = null;
			for (String samp : samps) {
				if (vpop.getPopulationForInd(samp, RETRIEVE_TYPE.SUPER)[0].equals(VcfPopulation.TUMOR)) {
					tumor = samp;
				} else if (vpop.getPopulationForInd(samp, RETRIEVE_TYPE.SUPER)[0].equals(VcfPopulation.NORMAL)) {
					normal = samp;
				} else {
					throw new IllegalArgumentException("Unknown types");
				}
			}

			if (!all.containsKey(tumor) || !all.containsKey(normal)) {
				throw new IllegalArgumentException("Could not find bam file for Tumor " + tumor + " or for Normal " + normal);
			} else {
				TNSample tSample = new TNSample(all.get(normal), normal, all.get(tumor), tumor, somaticParams, log);
				tnSamples.add(tSample);
				analysisBams.add(all.get(normal));
				analysisBams.add(all.get(tumor));

			}
		}
		if (analysisBams.size() < bamFiles.length) {
			Files.writeList(Array.toStringArray(analysisBams), outputDir + ext.rootOf(vpop.getFileName() + ".analysis.bams.txt"));
		}

		log.reportTimeInfo("Matched " + tnSamples.size() + " tumor normals with bams ");
		return tnSamples.toArray(new TNSample[tnSamples.size()]);
	}

	private static class SomaticParams {
		private String somaticSnipeLoc;
		private String refGenome;
		private String outputDir;
		private String samtoolsLoc;
		private boolean overwriteExisting;

		public SomaticParams(String somaticSnipeLoc, String refGenome, String outputDir, boolean overwriteExisting) {
			super();
			this.somaticSnipeLoc = somaticSnipeLoc;
			this.refGenome = refGenome;
			this.outputDir = outputDir;
			this.overwriteExisting = overwriteExisting;
			this.samtoolsLoc = somaticSnipeLoc + "build/vendor/samtools/";

		}

		public String getSamtoolsLoc() {
			return samtoolsLoc;
		}

		public String getOutputDir() {
			return outputDir;
		}

		public boolean isOverwriteExisting() {
			return overwriteExisting;
		}

		public String getRefGenome() {
			return refGenome;
		}

		public String getSnipe() {
			return somaticSnipeLoc + "build/bin/" + "bam-somaticsniper";
		}

		private String[] generateCommand(String normalBam, String tumorBam, String output) {
			ArrayList<String> command = new ArrayList<String>();
			command.add(getSnipe());
			command.add("-F");
			command.add("vcf");

			command.add("-Q");
			command.add("40");
			command.add("-G");
			command.add("-L");
			command.add("-f");
			command.add(refGenome);
			command.add(tumorBam);
			command.add(normalBam);
			command.add(output);
			return Array.toStringArray(command);
		}
	}

	private static class TNSample implements Callable<TNSample> {
		private String normalBam;
		private String normalSample;
		private String tumorBam;
		private String tumorSample;
		private SomaticParams somaticParams;
		private String indelPile;
		private String indelPileFilt;
		private String output;
		private String outputGz;
		private Logger log;

		public TNSample(String normalBam, String normalSample, String tumorBam, String tumorSample, SomaticParams somaticParams, Logger log) {
			super();
			this.normalBam = normalBam;
			this.normalSample = normalSample;
			this.tumorBam = tumorBam;
			this.tumorSample = tumorSample;
			this.somaticParams = somaticParams;
			this.output = somaticParams.getOutputDir() + "T_" + tumorSample + "_N_" + normalSample + ".vcf";
			this.outputGz = VCFOps.getAppropriateRoot(output, false) + "_rename.vcf.gz";
			this.indelPile = ext.rootOf(output, false) + ".indel.pileup";
			this.indelPileFilt = ext.rootOf(output, false) + ".indel.pileup";

			this.log = log;
		}

		public String getOutputGz() {
			return outputGz;
		}

		public String getNormalBam() {
			return normalBam;
		}

		public String getTumorBam() {
			return tumorBam;
		}

		private boolean runSamToolsIndel() {
			String[] inputs = new String[] { tumorBam, somaticParams.getRefGenome() };
			String[] output = new String[] { indelPile };
			ArrayList<String> commandCall = new ArrayList<String>();
			commandCall.add(somaticParams.getSamtoolsLoc() + "samtools");
			commandCall.add("pileup");
			commandCall.add("–cvi");
			commandCall.add("–f");
			commandCall.add(somaticParams.getRefGenome());
			commandCall.add(tumorBam);
			commandCall.add(indelPile);

			boolean progress = CmdLine.runCommandWithFileChecks(Array.toStringArray(commandCall), "", inputs, output, true, somaticParams.isOverwriteExisting(), false, log);
			if (progress) {
				ArrayList<String> commandFilt = new ArrayList<String>();
				commandFilt.add("perl");
				commandFilt.add(somaticParams.getSamtoolsLoc() + "misc/samtools.pl");
				commandFilt.add("varFilter");
				commandFilt.add(indelPile);
				commandFilt.add(">");
				commandFilt.add(indelPileFilt);
				String[] runBat = CmdLine.prepareBatchForCommandLine(Array.toStringArray(commandFilt), indelPileFilt + ".sh", true, log);
				progress = CmdLine.runCommandWithFileChecks(runBat, "", new String[] { indelPile }, new String[] { indelPileFilt }, true, somaticParams.isOverwriteExisting(), false, log);
			}
			return progress;

		}

		@Override
		public TNSample call() throws Exception {
			String[] inputs = new String[] { somaticParams.getRefGenome(), somaticParams.getSnipe(), tumorBam, normalBam };
			String[] outputs = new String[] { output };
			String[] command = somaticParams.generateCommand(normalBam, tumorBam, output);

			boolean progress = CmdLine.runCommandWithFileChecks(command, "", inputs, outputs, true, somaticParams.isOverwriteExisting(), false, log);

			if (progress) {

				progress = runSamToolsIndel();
				if (!Files.exists(outputGz)) {
					VCFFileReader reader = new VCFFileReader(output, false);
					VariantContextWriter writer = VCFOps.initWriter(outputGz, VCFOps.DEFUALT_WRITER_OPTIONS, reader.getFileHeader().getSequenceDictionary());
					Set<String> samps = new HashSet<String>();
					samps.add(normalSample);
					samps.add(tumorSample);
					final VCFHeader outHeader = new VCFHeader(reader.getFileHeader().getMetaDataInInputOrder(), samps);
					writer.writeHeader(outHeader);
					for (VariantContext vc : reader) {
						VariantContextBuilder builder = new VariantContextBuilder(vc);
						ArrayList<Genotype> renamed = new ArrayList<Genotype>();
						Genotype normal = rename(vc.getGenotype("NORMAL"), normalSample);
						Genotype tumor = rename(vc.getGenotype("TUMOR"), tumorSample);
						renamed.add(normal);
						renamed.add(tumor);
						builder.genotypes(renamed);
						if (!renamed.get(0).sameGenotype(vc.getGenotype("NORMAL"))) {
							reader.close();
							writer.close();
							throw new IllegalStateException("Improprer rename");
						}
						builder.genotypes(renamed);
						if (!renamed.get(1).sameGenotype(vc.getGenotype("TUMOR"))) {
							reader.close();
							writer.close();
							throw new IllegalStateException("Improprer rename");
						}

						writer.add(builder.make());
					}
					log.reportTimeInfo("Re-named and indexed " + output + " to " + outputGz);
					reader.close();
					writer.close();
				}

			}
			if (progress) {
				return this;
			} else {
				return null;
			}
		}
	}

	private static Genotype rename(Genotype g, String newName) {
		GenotypeBuilder builder = new GenotypeBuilder(g);
		builder.name(newName);
		return builder.make();
	}

	public static void main(String[] args) {
		String bams = "/scratch.global/lane0212/Project_Tsai_Project_028/151120_D00635_0090_BC81JNANXX/somaticSniper/results/TN.vpop.analysis.bams";
		String vpop = "/scratch.global/lane0212/Project_Tsai_Project_028/151120_D00635_0090_BC81JNANXX/somaticSniper/TN.vpop";
		String outputDir = "/scratch.global/lane0212/Project_Tsai_Project_028/151120_D00635_0090_BC81JNANXX/somaticSniper/results/";
		String refGenome = "/home/pankrat2/public/bin/ref/hg19_canonical.fa";
		String somaticSnipLoc = "/home/pankrat2/public/bin/somaticSniper/somatic-sniper/";
		String gatkLoc = "/home/pankrat2/public/bin/GATK/";

		int numThreads = 1;
		new File(outputDir).mkdirs();
		Logger log = new Logger(outputDir + "somaticSniper.log");
		GATK gatk = new GATK(gatkLoc, refGenome, true, false, log);
		SomaticParams params = new SomaticParams(somaticSnipLoc, gatk.getReferenceGenomeFasta(), outputDir, false);
		run(gatk, params, bams, vpop, numThreads, log);

	}

}
