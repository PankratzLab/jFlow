package seq.manage;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Set;
import java.util.Vector;
import java.util.concurrent.Callable;

import cnv.var.LocusSet;
import cnv.var.LocusSet.TO_STRING_TYPE;
import seq.analysis.GATK;
import seq.analysis.PlinkSeq;
import seq.analysis.PlinkSeq.ANALYSIS_TYPES;
import seq.analysis.PlinkSeq.LOAD_TYPES;
import seq.analysis.PlinkSeq.PlinkSeqWorker;
import seq.analysis.PlinkSeqUtils.PseqProject;
import seq.manage.BEDFileReader.BEDFeatureSeg;
import seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import seq.manage.VCOps.VC_SUBSET_TYPE;
import seq.qc.FilterNGS.VARIANT_FILTER_BOOLEAN;
import seq.qc.FilterNGS.VARIANT_FILTER_DOUBLE;
import seq.qc.FilterNGS.VariantContextFilter;
import stats.Histogram.DynamicAveragingHistogram;
import stats.Histogram.DynamicHistogram;
import filesys.Segment;
import gwas.MatchSamples;
import gwas.MatchesVisualized;
import gwas.MergeDatasets;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndex;
import htsjdk.tribble.util.LittleEndianOutputStream;
import htsjdk.tribble.util.TabixUtils;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import common.Array;
import common.CmdLine;
import common.Files;
import common.HashVec;
import common.Logger;
import common.PSF;
import common.Positions;
import common.WorkerHive;
import common.WorkerTrain;
import common.ext;
import common.WorkerTrain.Producer;

/**
 * Class for common actions on a VCF
 *
 */
public class VCFOps {
	public static final Set<String> BLANK_SAMPLE = new HashSet<String>();
	public static final Options[] DEFUALT_WRITER_OPTIONS = new Options[] { Options.INDEX_ON_THE_FLY };

	private static final String[] ANNO_BASE = new String[] { "CHROM", "POS", "ID", "REF", "ALT", "NUM_HOM_VAR", "NUM_HET", "NUM_HOM_ALT" };

	public enum VCF_EXTENSIONS {
		GZIP_VCF(".vcf.gz"), REG_VCF(".vcf"), BCF(".bcf");

		private String literal;

		private VCF_EXTENSIONS(String literal) {
			this.literal = literal;
		}

		public String getLiteral() {
			return literal;
		}
	}

	public enum UTILITY_TYPE {
		/**
		 * Convert a vcf to plink and run gwas QC
		 */
		CONVERT_PLINK,

		/**
		 * Subset a vcf by a super population
		 */
		SUBSET_SUPER, /**
		 * Extracts var
		 */
		EXTRACT_SEGMENTS,

		/**
		 * Extract segments and make annotation file
		 */
		EXTRACT_SEGMENTS_ANNOTATION,
		/**
		 * Extract rsids from a vcf
		 */
		EXTRACT_IDS,
		/**
		 * Removed variants flagged as filtered
		 */
		REMOVE_FILTERED, /**
		 * gzip and index a vcf file
		 */
		GZIP,

		/**
		 * Creates a vpop file with all cases
		 */
		DUMP_SAMPLES,
		/**
		 * Use plinkSeq to qc a vcf
		 */
		QC,

		/**
		 * Annotate a vcf using a bed file
		 */
		BED_ANNOTATE,
		/**
		 * Determine Homogeneity for populations in a vcf
		 */
		HOMOGENEITY;
	}

	/**
	 * Initialize a writer where all the samples will be present in the new vcf
	 */
	public static VariantContextWriter initWriterWithHeader(final VCFFileReader in, final String output, final Options[] options, Logger log) {
		VariantContextWriter writer = initWriter(output, options, in.getFileHeader().getSequenceDictionary());
		copyHeader(in, writer, null, HEADER_COPY_TYPE.FULL_COPY, log);
		return writer;
	}

	/**
	 * @param output
	 *            the output file to write to
	 * @param options
	 *            Options to be passed to the builder
	 * @param sequenceDictionary
	 *            an existing sequence dictionary- for example:<br>
	 *            ({@link VCFFileReader#getFileHeader() } and then {@link VCFHeader#getSequenceDictionary()}
	 * @return an initialized writer
	 */
	public static VariantContextWriter initWriter(final String output, final Options[] options, final SAMSequenceDictionary sequenceDictionary) {
		VariantContextWriterBuilder builder = new VariantContextWriterBuilder().setOutputFile(output);
		if (options != null) {
			for (int i = 0; i < options.length; i++) {
				builder.setOption(options[i]);

			}
		}
		if (sequenceDictionary != null) {
			builder.setReferenceDictionary(sequenceDictionary);
		}

		VariantContextWriter writer = builder.build();
		return writer;
	}

	/**
	 * @param vcf
	 * @return appropriate index
	 */
	public static String getIndex(String vcf) {
		return vcf.endsWith("gz") ? vcf + GATK.VCF_GZ_INDEX : GATK.VCF_INDEX;

	}

	/**
	 * @param vcf
	 * @return whether the vcf and appropriate index exists
	 */
	public static boolean existsWithIndex(String vcf) {
		return Files.exists(vcf) && Files.exists(getIndex(vcf));
	}

	public enum HEADER_COPY_TYPE {
		/**
		 * Site only header stripping sample info
		 */
		SITE_ONLY, /**
		 * Copy header from input
		 */
		FULL_COPY, /**
		 * Samples not contained in the vcf will be given missing genotypes
		 */
		SUBSET_LOOSE, /**
		 * A check will be performed and only samples present in the input set and the vcf file will be exported
		 */
		SUBSET_STRICT;
	}

	/**
	 * @param vcfFileReader
	 *            taker header from
	 * @param writer
	 *            write header to
	 * @param samples
	 *            subset to these samples only. To obtain a site only output, use {@link VCFOps#BLANK_SAMPLE}
	 * @return
	 */
	public static VariantContextWriter copyHeader(final VCFFileReader vcfFileReader, final VariantContextWriter writer, final Set<String> samples, HEADER_COPY_TYPE type, Logger log) {

		switch (type) {
		case FULL_COPY:
			writer.writeHeader(vcfFileReader.getFileHeader());
			break;
		case SITE_ONLY:
			writer.writeHeader(new VCFHeader(vcfFileReader.getFileHeader().getMetaDataInInputOrder(), BLANK_SAMPLE));
			break;
		case SUBSET_LOOSE:
			writer.writeHeader(new VCFHeader(vcfFileReader.getFileHeader().getMetaDataInInputOrder(), samples));
			break;

		case SUBSET_STRICT:
			ArrayList<String> samplesHave = vcfFileReader.getFileHeader().getSampleNamesInOrder();
			ArrayList<String> newSampleSubset = new ArrayList<String>();
			for (String samp : samplesHave) {
				if (samples.contains(samp)) {
					newSampleSubset.add(samp);
				}
			}
			writer.writeHeader(new VCFHeader(vcfFileReader.getFileHeader().getMetaDataInInputOrder(), newSampleSubset));

			break;
		default:
			break;

		}

		return writer;
	}

	public static String[] getSamplesInFile(String vcf) {
		VCFFileReader reader = new VCFFileReader(vcf, false);
		String[] samples = getSamplesInFile(reader);
		reader.close();
		return samples;
	}

	public static String[] getSamplesInFile(VCFFileReader reader) {
		List<String> samples = reader.getFileHeader().getGenotypeSamples();
		return samples.toArray(new String[samples.size()]);
	}

	public static boolean hasInfoLine(VCFFileReader reader, String anno) {
		return reader.getFileHeader().hasInfoLine(anno);
	}

	/**
	 * Retrieves the sequence dictionary from a reader
	 * 
	 * @param vcfFileReader
	 * @return
	 */
	public static SAMSequenceDictionary getSequenceDictionary(final VCFFileReader vcfFileReader) {
		return vcfFileReader.getFileHeader().getSequenceDictionary();
	}

	/**
	 * Retrieves the sequence dictionary from a reader
	 * 
	 * @param vcfFileReader
	 * @return
	 */
	public static SAMSequenceDictionary getSequenceDictionary(final String vcf) {
		VCFFileReader reader = new VCFFileReader(vcf, false);
		SAMSequenceDictionary samSequenceDictionary = getSequenceDictionary(reader);
		reader.close();
		return samSequenceDictionary;
	}

	public enum PLINK_SET_MODE {
		GWAS_QC, HOMOGENEITY;
	}

	public static boolean generateStripVCF(String inputVCF, String outputVCF, int maxAlleles, Logger log) {
		System.out.println("Why are you using this method");
		System.exit(1);
		if (Files.exists(outputVCF)) {
			log.reportFileExists(outputVCF);
			return true;
		}
		HashSet<String> samps = new HashSet<String>();
		for (int i = 0; i < maxAlleles; i++) {
			samps.add("Allele" + i);
		}
		VCFFileReader reader = new VCFFileReader(inputVCF, false);
		VariantContextWriter writer = initWriter(outputVCF, DEFUALT_WRITER_OPTIONS, getSequenceDictionary(reader));
		copyHeader(reader, writer, samps, HEADER_COPY_TYPE.SUBSET_LOOSE, log);
		int num = 0;
		for (VariantContext vc : reader) {
			num++;
			if (num % 100000 == 0) {
				log.reportTimeInfo(num + " variants stripped");
			}
			// int totalDiff = vc.getAlleles().size();

			GenotypesContext genotypesContext = vc.getGenotypes();
			ArrayList<Genotype> genos = new ArrayList<Genotype>();
			Hashtable<String, String> added = new Hashtable<String, String>();
			int currentIndex = 0;
			for (Genotype g : genotypesContext) {
				List<Allele> alleles = g.getAlleles();
				String aKey = alleles.get(0).getDisplayString();
				for (int i = 1; i < alleles.size(); i++) {
					aKey += "/" + alleles.get(i).getDisplayString();
				}
				if (!added.containsKey(aKey)) {
					GenotypeBuilder builder = new GenotypeBuilder("Allele" + currentIndex, alleles);
					currentIndex++;
					genos.add(builder.make());
					added.put(aKey, aKey);
				}
			}
			if (currentIndex > maxAlleles) {
				throw new IllegalArgumentException("max allele size was too small, found " + currentIndex + " alts for " + vc.toStringWithoutGenotypes());
			}
			while (currentIndex < samps.size()) {
				GenotypeBuilder builder = new GenotypeBuilder("Allele" + currentIndex, GenotypeOps.getNoCall());
				currentIndex++;
				genos.add(builder.make());
			}

			VariantContextBuilder vcBuilder = new VariantContextBuilder(vc);
			vcBuilder.genotypes(GenotypesContext.create(genos));
			VariantContext st = vcBuilder.make();
			if (!vc.isMonomorphicInSamples() && st.isMonomorphicInSamples()) {
				throw new IllegalArgumentException("non to mono");
			}
			if (vc.isMonomorphicInSamples() && !st.isMonomorphicInSamples()) {
				throw new IllegalArgumentException("mono to non");
			}
			if (!vc.hasSameAllelesAs(st)) {
				throw new IllegalArgumentException("diff alleles");
			}

			writer.add(st);

		}
		return Files.exists(outputVCF);
	}

	public static String[] reportCallRateHWEFiltered(String vcf, String outputFile, double callRate, double hweP, Logger log) {
		VARIANT_FILTER_DOUBLE callRateFilter = VARIANT_FILTER_DOUBLE.CALL_RATE_LOOSE;
		VARIANT_FILTER_DOUBLE hwe = VARIANT_FILTER_DOUBLE.HWE;
		hwe.setDFilter(hweP);
		callRateFilter.setDFilter(callRate);
		VariantContextFilter filter = new VariantContextFilter(new VARIANT_FILTER_DOUBLE[] { callRateFilter, hwe }, new VARIANT_FILTER_BOOLEAN[] {}, null, null, log);
		DynamicAveragingHistogram dynamicAveragingHistogram = new DynamicAveragingHistogram(0, 1.1, 2);

		try {
			PrintWriter writer = new PrintWriter(new FileWriter(outputFile));
			PrintWriter writerAll = new PrintWriter(new FileWriter(ext.addToRoot(outputFile, ".all")));
			String segFile = ext.addToRoot(outputFile, ".segment");
			PrintWriter writerFilteredSeg = new PrintWriter(new FileWriter(segFile));

			VCFFileReader reader = new VCFFileReader(vcf, true);
			int count = 0;
			for (VariantContext vc : reader) {
				count++;
				if (count % 100000 == 0) {
					log.reportTimeInfo(count + " variants tested at callrate " + callRate);
				}
				if (!filter.filter(vc).passed()) {
					writer.println((vc.getID().equals(".") ? new VCOps.LocusID(vc).getId() : vc.getID()) + "\t" + vc.getNoCallCount() + "\t" + vc.getNSamples());
					Segment vcSeg = VCOps.getSegment(vc);
					writerFilteredSeg.println(vcSeg.getChr() + "\t" + vcSeg.getStart() + "\t" + vcSeg.getStop());
				}
				writerAll.println((vc.getID().equals(".") ? new VCOps.LocusID(vc).getId() : vc.getID()) + "\t" + vc.getNoCallCount() + "\t" + vc.getNSamples());
				double g1000 = -1;
				try {
					if (vc.getCommonInfo().hasAttribute("g10002014oct_all")) {

						String g = VCOps.getAnnotationsFor(new String[] { "g10002014oct_all" }, vc, ".")[0];
						if (g.equals(".")) {
							g1000 = 0;
						} else {
							g1000 = Double.parseDouble(g);
						}

					} else {
						// notAnnotated.addDataPair(VCOps.getCallRate(vc, null), Double.parseDouble(vc.getCommonInfo().getAttributeAsDouble("AF", "0")));
					}
				} catch (NumberFormatException nfe) {

				}
				dynamicAveragingHistogram.addDataPair(g1000, VCOps.getCallRate(vc, null));
			}
			reader.close();
			writer.close();
			writerAll.close();
			writerFilteredSeg.close();
			String hist = ext.addToRoot(outputFile, ".callrateHist");
			PrintWriter writerHist = new PrintWriter(new FileWriter(hist));
			writerHist.println("g1000Bin\tCount\tAvgCallRate");
			dynamicAveragingHistogram.average();
			for (int i = 0; i < dynamicAveragingHistogram.getBins().length; i++) {
				writerHist.println(dynamicAveragingHistogram.getBins()[i] + "\t" + dynamicAveragingHistogram.getCounts()[i] + "\t" + dynamicAveragingHistogram.getAverages()[i]);
			}
			writerHist.close();
			return new String[] { hist, segFile };
		} catch (Exception e) {
			log.reportError("Error writing to " + outputFile);
			log.reportException(e);
		}
		return null;
	}

	/**
	 * @param vcf
	 *            a vcf file to convert to plink
	 * @param rootOut
	 *            the root output for the plink files
	 * @param log
	 */

	public static String[] convertToPlinkSet(String outputDir, String vcf, String rootOut, PLINK_SET_MODE mode, Logger log) {
		String[] plinkCommand = null;
		String dir = outputDir == null ? ext.parseDirectoryOfFile(vcf) + "plink" + ext.rootOf(vcf) + "/" : outputDir;

		new File(dir).mkdirs();

		rootOut = dir + rootOut;
		String[] outFiles = PSF.Plink.getPlinkBedBimFam(rootOut);
		if (!Files.exists("", outFiles)) {
			plinkCommand = PSF.Plink.getPlinkVCFCommand(vcf, rootOut);
			if (CmdLine.runCommandWithFileChecks(plinkCommand, "", new String[] { vcf }, outFiles, true, true, false, log)) {

			}
		} else {
			log.reportTimeWarning("Detected that the following files already exist " + Array.toStr(outFiles));
		}
		if (Files.exists("", outFiles)) {
			// Hashtable<String, String> newIDS = new Hashtable<String, String>();
			fixFamFile(log, outFiles[2]);
			log.reportTimeInfo("MODE=" + mode);
			if (mode == PLINK_SET_MODE.GWAS_QC) {

				gwas.Qc.fullGamut(dir, rootOut, false, new Logger(dir + "fullGamutOfMarkerAndSampleQC.log"));
				String mdsFile = dir + "genome/mds20.mds";
				if (Files.exists(mdsFile)) {
					// fixMdsFile(log, dir, newIDS, mdsFile);
					// CmdLine.run("runEigenstratWoHapMap", dir + "ancestry/");
					// CmdLine.run("runEigenstrat2", dir + "ancestry/");
					// fixMdsFile(log, dir + "ancestry/", newIDS, combo_fancy_postnormed_eigens.xln);
				}
			} else if (mode == PLINK_SET_MODE.HOMOGENEITY) {

				if (!Files.exists(ext.parseDirectoryOfFile(outFiles[2]) + "hardy.hwe")) {
					CmdLine.run("plink --bfile " + rootOut + " --maf 0 --geno 1 --mind 1 --hardy --out hardy --noweb", ext.parseDirectoryOfFile(outFiles[2]));
				} else {
					log.reportTimeInfo("Found file " + ext.parseDirectoryOfFile(outFiles[2]) + "hardy.hwe");
				}
			}
		}
		if (!Files.exists(dir + ".qc.pbs")) {
			String gwasQC = Array.toStr(PSF.Load.getAllModules(), "\n") + "\njcp gwas.Qc dir=" + dir;
			Files.qsub(dir + "qc.pbs", gwasQC, 62000, 24, 16);
		}
		return outFiles;
	}

	// private static void fixMdsFile(Logger log, String dir, Hashtable<String, String> newIDS, String mdsFile) {
	// if (newIDS.size() > 0 && Files.exists(mdsFile)) {
	// String[] header = Files.getHeaderOfFile(mdsFile, log);
	// int[] indices = new int[header.length];
	// for (int i = 0; i < indices.length; i++) {
	// indices[i] = i;
	// }
	// String[][] mds = HashVec.loadFileToStringMatrix(mdsFile, false, new int[] { 0, 1, 2, 3, 4, 5 }, false);
	// String[][] newMds = new String[mds.length][mds[0].length];
	// for (int i = 0; i < mds.length; i++) {
	// for (int j = 0; j < mds[i].length; j++) {
	// if (j == 0 && j != 1 && newIDS.containsKey(mds[i][0])) {
	// newMds[i][0] = newIDS.get(mds[i][0]);
	// newMds[i][1] = newIDS.get(mds[i][0]);
	// } else {
	// newMds[i][j] = mds[i][j];
	// }
	// }
	// }
	// Files.writeMatrix(newMds, ext.addToRoot(mdsFile, ".fixed"), "\t");
	//
	// }
	// }

	private static Hashtable<String, String> fixFamFile(Logger log, String famFile) {
		Hashtable<String, String> changedIds = new Hashtable<String, String>();
		Files.copyFile(famFile, famFile + ".bak");
		String[][] fam = HashVec.loadFileToStringMatrix(famFile, false, new int[] { 0, 1, 2, 3, 4, 5 }, false);
		String[][] newfam = new String[fam.length][fam[0].length];
		boolean newSex = false;
		// String[][] fam = HashVec.loadFileToStringMatrix(, false, new int[]{1,2,3,4,5,6}, "[\\s]+", false, 1000, false);
		int noSexcount = 0;

		for (int i = 0; i < fam.length; i++) {
			if (fam[i][4].equals("0")) {
				noSexcount++;
			}
		}
		if (noSexcount == fam.length) {
			newSex = true;
			log.reportTimeWarning("Assigning alternating sex specifications");
		}
		Hashtable<String, String> uniqIds = new Hashtable<String, String>();
		boolean swit = true;
		for (int i = 0; i < fam.length; i++) {
			String FidIid = fam[i][0];
			if (FidIid.length() >= 37) {
				String newID = FidIid.substring(0, 38);
				log.reportTimeWarning("Changing " + FidIid + " to " + newID);
				changedIds.put(newID, FidIid);
				// uniqIds.put(FidIid, newID);
				FidIid = newID;
			}
			uniqIds.put(FidIid, FidIid);

			newfam[i][0] = FidIid;
			newfam[i][1] = FidIid;
			if (newSex && swit) {
				newfam[i][4] = "1";
				swit = false;
			} else if (newSex && !swit) {
				newfam[i][4] = "2";
				swit = true;
			} else {
				newfam[i][4] = fam[i][4];
			}
			for (int j = 0; j < newfam[i].length; j++) {
				if (j != 4 && j != 0 && j != 1) {
					newfam[i][j] = fam[i][j];
				}
			}

		}
		if (uniqIds.size() != fam.length) {
			log.reportTimeError("Could not remedy fam file");
		} else {
			log.reportTimeInfo("fixed fam file");
			Files.writeMatrix(newfam, famFile, "\t");
			if (changedIds.size() > 0) {
				try {
					PrintWriter writer = new PrintWriter(new FileWriter(famFile + ".changedIds"));
					for (String newID : changedIds.keySet()) {
						writer.print(newID + "\t" + changedIds.get(newID));
					}
					writer.close();
				} catch (Exception e) {
					log.reportError("Error writing to " + famFile + ".changedIds");
					log.reportException(e);
				}
			}
		}
		return changedIds;
	}

	/**
	 * @param vcf
	 *            runs {@link gwas.Qc#fullGamut(String, boolean)} after converting to plink* files if neccesary
	 * @param log
	 */
	public static void vcfGwasQC(String vcf, Logger log) {
		if (Files.exists(vcf)) {
			String dir = ext.parseDirectoryOfFile(vcf);
			String[] plinkFiles = PSF.Plink.getPlinkBedBimFam("plink");
			if (!Files.exists(dir, plinkFiles)) {
				log.reportTimeInfo("Generating plink files for " + vcf + " in " + dir);
				convertToPlinkSet(null, vcf, "plink", PLINK_SET_MODE.GWAS_QC, log);
			}
			log.reportTimeInfo("Running gwas.qc on the following files in " + dir + ":");
			log.reportTimeInfo("\t" + Array.toStr(plinkFiles, "\n"));
			// gwas.Qc.fullGamut(dir, false, new Logger(dir + "fullGamutOfMarkerAndSampleQC.log"));
		} else {
			log.reportFileNotFound(vcf);
		}
	}

	/**
	 * Class that manages the population structure represented in a vcf<br>
	 * Can be used for HWE tests on sub and super populations etc...
	 *
	 */
	public static class VcfPopulation implements Serializable {
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;
		public static final String ANCHOR = "ANCHOR";
		public static final String BARNACLE = "BARNACLE";
		public static final String CASE = "CASE";
		public static final String CONTROL = "CONTROL";
		public static final String EXCLUDE = "EXCLUDE";
		public static final String TUMOR = "TUMOR";
		public static final String NORMAL = "NORMAL";
		public static final String OFFSPRING = "OFFSPRING";
		public static final String PARENTS = "PARENTS";

		public static final String DETERMINE_ANCESTRY = "DETERMINE_ANCESTRY";
		public static final String[] HEADER = new String[] { "IID", "Population", "SuperPopulation" };
		private static final String SKIP = "#N/A";
		private Hashtable<String, Set<String>> subPop;
		private Hashtable<String, Set<String>> superPop;
		private ArrayList<String> uniqSubPop;
		private ArrayList<String> uniqSuperPop;
		private POPULATION_TYPE type;
		private String fileName;
		private Logger log;

		public enum POPULATION_TYPE {
			CASE_CONTROL, ANY, STRATIFICATION, EXOME_DEPTH, PC_ANCESTRY, ANCHOR_BARNACLE, TUMOR_NORMAL, DENOVO;
		}

		public enum RETRIEVE_TYPE {
			SUPER, SUB;
		}

		public VcfPopulation(Hashtable<String, Set<String>> subPop, Hashtable<String, Set<String>> superPop, POPULATION_TYPE type, Logger log) {
			super();
			this.subPop = subPop;
			this.superPop = superPop;
			this.type = type;
			this.uniqSubPop = new ArrayList<String>();
			uniqSubPop.addAll(subPop.keySet());
			this.uniqSuperPop = new ArrayList<String>();
			uniqSuperPop.addAll(superPop.keySet());
			this.fileName = "inMemoryPop";
			this.log = log;
		}

		public VcfPopulation(POPULATION_TYPE type, Logger log) {
			this.subPop = new Hashtable<String, Set<String>>();
			this.superPop = new Hashtable<String, Set<String>>();
			superPop.put(EXCLUDE, new HashSet<String>());
			this.uniqSubPop = new ArrayList<String>();
			this.uniqSuperPop = new ArrayList<String>();
			this.type = type;
			this.log = log;
		}

		public Logger getLog() {
			return log;
		}

		public String[] getOffP1P2ForFam(String fam) {
			if (type != POPULATION_TYPE.DENOVO) {
				throw new IllegalArgumentException("Method can only be used for denovo");
			}
			if (!subPop.containsKey(fam) || subPop.get(fam).size() != 3) {
				throw new IllegalArgumentException("Invalid fam");

			} else {
				String off = null;
				String p1 = null;
				String p2 = null;
				for (String ind : subPop.get(fam)) {
					if (getPopulationForInd(ind, RETRIEVE_TYPE.SUPER)[0].equals(OFFSPRING)) {
						off = ind;
					} else if (getPopulationForInd(ind, RETRIEVE_TYPE.SUPER)[0].equals(PARENTS)) {
						if (p1 == null) {
							p1 = ind;
						} else {
							p2 = ind;
						}
					}
				}
				return new String[] { off, p1, p2 };
			}
		}

		public String[] getPopulationForInd(String ind, RETRIEVE_TYPE type) {
			ArrayList<String> tmp = new ArrayList<String>();
			Set<String> avail;
			switch (type) {
			case SUB:
				avail = subPop.keySet();
				for (String key : avail) {
					if (subPop.get(key).contains(ind)) {
						tmp.add(key);
					}
				}
				break;
			case SUPER:
				avail = superPop.keySet();
				for (String key : avail) {
					if (superPop.get(key).contains(ind)) {
						tmp.add(key);
					}
				}
				break;
			default:
				log.reportTimeError("Invalid type " + type);
				break;

			}
			return tmp.toArray(new String[tmp.size()]);
		}

		private HashSet<String> getAllInds() {
			HashSet<String> allInds = new HashSet<String>();
			for (String key : superPop.keySet()) {
				allInds.addAll(superPop.get(key));
			}
			for (String key : subPop.keySet()) {
				allInds.addAll(subPop.get(key));
			}
			return allInds;
		}

		public void dump(String file) {
			try {
				PrintWriter writer = new PrintWriter(new FileWriter(file));
				writer.println(Array.toStr(HEADER));
				HashSet<String> inds = getAllInds();
				for (String ind : inds) {
					writer.println(ind + "\t" + getPopulationForInd(ind, RETRIEVE_TYPE.SUB)[0] + "\t" + getPopulationForInd(ind, RETRIEVE_TYPE.SUPER)[0]);
				}
				writer.close();
			} catch (Exception e) {
				log.reportError("Error writing to " + file);
				log.reportException(e);
			}
		}

		public boolean generatePlinkSeqPheno(String output) {
			boolean generated = false;
			if (type != POPULATION_TYPE.CASE_CONTROL) {
				log.reportTimeError("Population type must be set to " + POPULATION_TYPE.CASE_CONTROL);
				return generated;
			} else if (!valid()) {
				return generated;
			} else {

				try {
					PrintWriter writer = new PrintWriter(new FileWriter(output));
					writer.println("##" + CASE + ",Integer,0,Primary disease phenotype");
					writer.println("#ID\t" + CASE);
					Set<String> cases = subPop.get(CASE);
					Set<String> controls = subPop.get(CONTROL);
					for (String aCase : cases) {
						writer.println(aCase + "\t2");
					}
					for (String control : controls) {
						writer.println(control + "\t1");
					}
					writer.close();
					generated = true;
				} catch (Exception e) {
					log.reportError("Error writing to " + output);
					log.reportException(e);
				}
			}
			return generated;
		}

		public boolean valid() {
			boolean valid = true;
			switch (type) {
			case ANY:
				break;
			case CASE_CONTROL:
				if (!subPop.containsKey(CASE) || !subPop.containsKey(CONTROL)) {
					log.reportTimeError("Population type was set to " + type + ", but did not contain " + CASE + " and  " + CONTROL);
				}
				break;
			case STRATIFICATION:
				break;
			case PC_ANCESTRY:
				if (!superPop.containsKey(DETERMINE_ANCESTRY)) {
					log.reportTimeError("Population type was set to " + type + ", but did not contain " + DETERMINE_ANCESTRY);
					log.reportTimeError(DETERMINE_ANCESTRY + " must be present in the " + HEADER[2] + " column as a flag to determine ancestry, all other categories will be used as cluster generators");

					valid = false;
				}
				break;
			case ANCHOR_BARNACLE:
				if (!superPop.containsKey(ANCHOR) || !superPop.containsKey(ANCHOR)) {
					log.reportTimeError("Population type was set to " + type + ", but did not contain " + ANCHOR + " AND " + BARNACLE);
					valid = false;
				}
				break;
			case TUMOR_NORMAL:
				if (superPop.size() != 3) {
					throw new IllegalArgumentException(POPULATION_TYPE.TUMOR_NORMAL + " must have two and only three types and found " + superPop.keySet());

				} else if (!superPop.containsKey(TUMOR) || !superPop.containsKey(NORMAL)) {
					throw new IllegalArgumentException(POPULATION_TYPE.TUMOR_NORMAL + " must only contain " + TUMOR + " and " + NORMAL);
				}
				break;

			case DENOVO:
				if (superPop.size() != 3) {
					throw new IllegalArgumentException(POPULATION_TYPE.DENOVO + " must have two and found " + superPop.keySet());

				} else if (!superPop.containsKey(OFFSPRING) || !superPop.containsKey(PARENTS)) {
					throw new IllegalArgumentException(POPULATION_TYPE.DENOVO + " must only contain " + OFFSPRING + " and " + PARENTS);
				}
				break;
			default:
				break;

			}
			return valid;
		}

		/**
		 * @param IID
		 *            sample name
		 * @param sub
		 *            samples sub population
		 * @param sup
		 *            samples super population
		 */
		public void add(String IID, String sub, String sup) {
			if (!sub.equals(SKIP)) {
				if (!subPop.containsKey(sub)) {
					subPop.put(sub, new HashSet<String>());
					uniqSubPop.add(sub);
				}
				subPop.get(sub).add(IID);
			}
			if (!sup.equals(SKIP)) {
				if (!superPop.containsKey(sup)) {
					superPop.put(sup, new HashSet<String>());
					uniqSuperPop.add(sup);
				}
				superPop.get(sup).add(IID);
			}
		}

		public Hashtable<String, Set<String>> getSubPop() {
			return subPop;
		}

		public Hashtable<String, Set<String>> getSuperPop() {
			return superPop;
		}

		public ArrayList<String> getUniqSuperPop() {
			return uniqSuperPop;
		}

		/**
		 * Breakdown of sample sizes per population
		 */
		public void report() {
			log.reportTimeInfo("Population type :" + type);
			for (String key : subPop.keySet()) {
				log.reportTimeInfo("Sub - population " + key + " had " + subPop.get(key).size() + " individuals");
			}
			for (String key : superPop.keySet()) {
				log.reportTimeInfo("Super - population " + key + " had " + superPop.get(key).size() + " individuals");
			}
		}

		public String getFileName() {
			return fileName;
		}

		public void setFileName(String fileName) {
			this.fileName = fileName;
		}

		private VariantContextWriter[] getWritersForPop(String outputbase, VCFFileReader reader, Logger log) {
			String[] filenames = getFileNamesForPop(outputbase, log);
			if (Files.exists("", filenames)) {
				log.reportTimeWarning("All split vcfs exist, skipping extraction");
				return null;
			} else {

				VariantContextWriter[] writers = new VariantContextWriter[filenames.length];

				for (int i = 0; i < filenames.length; i++) {
					if (Files.exists(filenames[i])) {
						writers[i] = null;
						log.reportTimeWarning(filenames[i] + " exists, skipping");
					} else {
						String sPop = uniqSuperPop.get(i);
						writers[i] = initWriter(filenames[i], DEFUALT_WRITER_OPTIONS, getSequenceDictionary(reader));
						copyHeader(reader, writers[i], superPop.get(sPop), HEADER_COPY_TYPE.SUBSET_STRICT, log);
					}
				}
				return writers;
			}

		}

		public String[] getFileNamesForPop(String outputbase, Logger log) {
			String[] filenames = new String[uniqSuperPop.size()];
			for (int i = 0; i < uniqSuperPop.size(); i++) {
				filenames[i] = outputbase + "." + uniqSuperPop.get(i) + ".vcf.gz";
			}
			return filenames;

		}

		public static String[] splitVcfByPopulation(String vcf, String fullPathToPopFile, boolean removeMonoMorphic, boolean keepIds, Logger log) {
			if (vcf != null && !Files.exists(vcf)) {
				log.reportFileNotFound(vcf);
				return null;
			}
			if (fullPathToPopFile != null && !Files.exists(fullPathToPopFile)) {
				log.reportFileNotFound(fullPathToPopFile);
				return null;
			}
			if (!keepIds) {
				log.reportTimeWarning("Variant context ids that are set to \".\" will be updated");
			}
			VcfPopulation vpop = VcfPopulation.load(fullPathToPopFile, POPULATION_TYPE.ANY, log);
			vpop.report();
			VCFFileReader reader = new VCFFileReader(vcf, true);

			String dir = ext.parseDirectoryOfFile(fullPathToPopFile);
			String root = getAppropriateRoot(vcf, true);
			log.reportTimeInfo("Writing to root split " + dir + root);
			VariantContextWriter[] writers = vpop.getWritersForPop(dir + root, reader, log);
			if (writers != null) {

				int progress = 0;
				int monoMorphic = 0;
				for (VariantContext vc : reader) {
					progress++;
					if (progress % 100000 == 0) {
						log.reportTimeInfo(progress + " variants read...");
					}
					if (!keepIds && vc.getID().equals(".")) {
						VariantContextBuilder builder = new VariantContextBuilder(vc);
						builder.id(new VCOps.LocusID(vc).getId());
						vc = builder.make();
					}
					for (int i = 0; i < writers.length; i++) {
						if (writers[i] != null) {
							VariantContext vcSub = VCOps.getSubset(vc, vpop.getSuperPop().get(vpop.getUniqSuperPop().get(i)));
							if (!removeMonoMorphic || !vcSub.isMonomorphicInSamples()) {
								writers[i].add(vcSub);
							} else {
								monoMorphic++;
							}
						}
					}
				}
				for (int i = 0; i < writers.length; i++) {
					if (writers[i] != null) {
						writers[i].close();
					}
				}
				if (monoMorphic > 0) {
					log.reportTimeInfo("Removed " + monoMorphic + " monomorphic variants ");
				}
			}
			return vpop.getFileNamesForPop(dir + root, log);
		}

		/**
		 * @param fullPathToPopFile
		 *            load this file to an {@link VcfPopulation}
		 * @param log
		 * @return
		 */
		public static VcfPopulation load(String fullPathToPopFile, POPULATION_TYPE type, Logger log) {
			VcfPopulation vcfPopulation = new VcfPopulation(type, log);
			try {
				BufferedReader reader = Files.getAppropriateReader(fullPathToPopFile);
				String[] header = Files.getHeaderOfFile(fullPathToPopFile, log);
				int[] indices = ext.indexFactors(HEADER, header, true, false);
				if (Array.countIf(indices, -1) > 0) {
					log.reportTimeError("Could not find required headers " + Array.toStr(HEADER) + " in " + fullPathToPopFile);
					return null;
				}
				reader.readLine();
				while (reader.ready()) {
					String[] line = reader.readLine().trim().split("[\\s]+");
					vcfPopulation.add(line[indices[0]], line[indices[1]], line[indices[2]]);
				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				log.reportError("Error: file \"" + fullPathToPopFile + "\" not found in current directory");
				return null;
			} catch (IOException ioe) {
				log.reportError("Error reading file \"" + fullPathToPopFile + "\"");
				return null;
			}
			vcfPopulation.setFileName(fullPathToPopFile);
			if (!vcfPopulation.valid()) {
				return null;
			} else {
				return vcfPopulation;
			}
		}

	}

	/**
	 * @param vcf
	 * @param fullPathToPopFiles
	 *            split the vcf by the definitions in this file and run homogeneity tests
	 * @param log
	 */
	public static void runHomoGeneity(String vcf, String[] fullPathToPopFiles, Logger log) {
		HashSet<String> toRemoveHash = new HashSet<String>();
		String[] toRemove = new String[] {};
		Segment[] toRemoveSeg = new Segment[] {};
		double callRate = 0.80;
		double hwe = .00001;
		int numBarnsPerSample = 5;
		String finalSamples = vcf + ".finalSamples";
		String[] matchUpVpops = Files.listFullPaths(ext.parseDirectoryOfFile(vcf), ".homogeneity.vpop", false);
		if (matchUpVpops.length < 1) {
			log.reportTimeError("Required file(s) ending with .homogeneity.vpop in directory " + ext.parseDirectoryOfFile(vcf) + " were not found");
			return;
		}
		if (!Files.exists(finalSamples)) {
			log.reportTimeError("Required file " + finalSamples + " is missing");
			return;
		}
		String finalDir = ext.parseDirectoryOfFile(vcf) + "homogeneity/";
		new File(finalDir).mkdirs();
		String idFile = finalDir + "variants_Removed.txt";
		String mdsFile = finalDir + "/genome/mds20.mds";

		if (!Files.exists(idFile) || !Files.exists(mdsFile)) {

			String[] samples = HashVec.loadFileToStringArray(finalSamples, false, new int[] { 0 }, false);
			log.reportTimeInfo("Found " + samples.length + " samples for the final analysis");
			for (int i = 0; i < fullPathToPopFiles.length; i++) {
				String fullPathToPopFile = fullPathToPopFiles[i];
				String[] splits = VcfPopulation.splitVcfByPopulation(vcf, fullPathToPopFile, false, false, log);
				String[] dirs = new String[splits.length];
				String dir = ext.parseDirectoryOfFile(vcf) + ext.rootOf(fullPathToPopFile) + "/";
				for (int j = 0; j < splits.length; j++) {

					String export = dir + "plink_" + ext.rootOf(splits[j]) + "/";
					dirs[j] = ext.parseDirectoryOfFile(convertToPlinkSet(export, splits[j], "plink", PLINK_SET_MODE.HOMOGENEITY, log)[0]);
					if (VCFOps.getSamplesInFile(splits[j]).length > 50) {
						String callRateFiltered = dir + ext.rootOf(splits[j]) + ".CR." + callRate + ".hwe." + hwe + ".txt";
						if (!Files.exists(callRateFiltered)) {
							reportCallRateHWEFiltered(splits[j], callRateFiltered, callRate, hwe, log);
						}
						// LocusSet<Segment> segs = LocusSet.loadSegmentSetFromFile(ext.addToRoot(callRateFiltered, ".segment"), 0, 1, 2, 0, true, true, 0, log);
						// toRemoveSeg = Array.concatAll(toRemoveSeg, segs.getLoci());
						String[] callRateRemove = HashVec.loadFileToStringArray(callRateFiltered, false, new int[] { 0 }, true);
						for (int k = 0; k < callRateRemove.length; k++) {
							toRemoveHash.add(callRateRemove[k]);
						}
						log.reportTimeInfo(callRateRemove.length + " variants removed from " + splits[j] + " at callrate " + callRateFiltered);
						toRemove = Array.unique(Array.concatAll(toRemove, callRateRemove));
					}
				}
				String lackOfHomoGeneity = dir + "lackOfHomogeneity.dat";
				String problems = dir + "problematic.dat";

				if (!Files.exists(lackOfHomoGeneity)) {
					MergeDatasets.checkForHomogeneity(null, dirs, dir, "ALL", log);

				} else {
					log.reportTimeInfo("Found " + lackOfHomoGeneity + ", assuming this has run to completion");
				}

				String[] lackOfHomoGeneityIDs = HashVec.loadFileToStringArray(lackOfHomoGeneity, false, new int[] { 0 }, true);
				log.reportTimeInfo(lackOfHomoGeneityIDs.length + " markers lacking homogeneity from " + fullPathToPopFiles[i]);
				toRemove = Array.unique(Array.concatAll(toRemove, lackOfHomoGeneityIDs));
				for (int j = 0; j < lackOfHomoGeneityIDs.length; j++) {
					toRemoveHash.add(lackOfHomoGeneityIDs[j]);
				}
				if (Files.exists(problems)) {
					String[] problematic = HashVec.loadFileToStringArray(problems, false, new int[] { 0 }, true);
					log.reportTimeInfo(problematic.length + " markers with problems from " + fullPathToPopFiles[i]);
					// toRemove = Array.unique(Array.concatAll(toRemove, problematic));
				}
			}

			Files.writeList(toRemoveHash.toArray(new String[toRemoveHash.size()]), idFile);
			HashSet<String> sampleHash = new HashSet<String>();
			for (int i = 0; i < samples.length; i++) {
				sampleHash.add(samples[i]);
			}
			LocusSet<Segment> sort = new LocusSet<Segment>(toRemoveSeg, true, log) {

				/**
			 * 
			 */
				private static final long serialVersionUID = 1L;

			};

			log.reportTimeInfo("Removing " + toRemove.length + " variants from " + vcf);
			log.reportTimeInfo("Removing " + toRemoveSeg.length + " segments as well");
			log.reportTimeInfo("Subsetting to " + sampleHash.size() + " samples");
			String extractVCF = extractIDs(vcf, idFile, finalDir, true, true, sampleHash, sort.getLoci(), true, false, log);

			if (!Files.exists(mdsFile)) {
				convertToPlinkSet(finalDir, extractVCF, "plink", PLINK_SET_MODE.GWAS_QC, log);
			}
		} else {
			log.reportTimeWarning("found file " + idFile + " and " + mdsFile + " , assuming processing up to this point has been completed");
		}
		widgets.TabVersion.make(mdsFile);
		mdsFile = mdsFile + ".xln";

		for (int i = 1; i < matchUpVpops.length; i++) {
			VcfPopulation vpop = VcfPopulation.load(matchUpVpops[i], POPULATION_TYPE.ANCHOR_BARNACLE, log);
			vpop.report();
			String matchDir = finalDir + "match_" + ext.rootOf(matchUpVpops[i]) + "/";
			new File(matchDir).mkdirs();

			String factorFile = matchDir + "plink.mds";
			Files.copyFileUsingFileChannels(new File(mdsFile), new File(factorFile), log);
			String[] barnacleIdsPresent = HashVec.loadFileToStringArray(factorFile, true, new int[] { 0 }, true);

			Set<String> anchors = vpop.getSuperPop().get(VcfPopulation.ANCHOR);
			Set<String> barnacles = vpop.getSuperPop().get(VcfPopulation.BARNACLE);

			ArrayList<String> currentBarns = new ArrayList<String>();

			for (int j = 0; j < numBarnsPerSample; j++) {
				HashSet<String> barnaclesPresent = new HashSet<String>();
				String[] currentBarnsA = currentBarns.toArray(new String[currentBarns.size()]);
				System.out.println("SIZE" + Array.unique(currentBarnsA).length);
				for (String barn : barnacles) {
					if (ext.indexOfStr(barn, barnacleIdsPresent) >= 0 && ext.indexOfStr(barn, currentBarnsA) < 0) {
						barnaclesPresent.add(barn);
					} else {
						// log.reportTimeWarning("Missing sample " + barn + " in file " + factorFile);
					}
				}

				System.out.println(barnaclesPresent.size());
				try {
					Thread.sleep(100);
				} catch (InterruptedException ie) {
				}
				String anchorList = matchDir + "anchors.txt";
				String barnacleList = matchDir + j + "barnacles.txt";

				Files.writeList(anchors.toArray(new String[anchors.size()]), anchorList);
				Files.writeList(barnaclesPresent.toArray(new String[barnaclesPresent.size()]), barnacleList);
				String[] run = new String[] { "C1", "C3" };
				System.out.println("RUNNING match1");
				String matchFile = MatchSamples.matchMaker(matchDir, ext.removeDirectoryInfo(anchorList), ext.removeDirectoryInfo(barnacleList), ext.removeDirectoryInfo(factorFile), run, new double[] { 1, 1 }, false);
				System.out.println("RUNNING match2");

				matchFile = MatchSamples.normalizeDistances(matchDir, matchFile, 0, 100);
				System.out.println("RUNNING match3");

				String pairs = matchDir + MatchSamples.matchPairs(matchDir, matchFile, true);
				System.out.println("RUNNING match4");

				Files.copyFileUsingFileChannels(new File(pairs), new File(pairs + j + ".selection"), log);

				String[] barnesPicked = HashVec.loadFileToStringArray(pairs, true, new int[] { 1 }, true);
				String[] deletes = Files.listFullPaths(matchDir, ".xln", false);
				new MatchesVisualized(matchDir, ext.removeDirectoryInfo(anchorList), ext.removeDirectoryInfo(barnacleList), ext.removeDirectoryInfo(factorFile), ext.indexFactors(run, Files.getHeaderOfFile(factorFile, log), true, false), ext.removeDirectoryInfo(pairs));

				for (int k = 0; k < deletes.length; k++) {
					new File(deletes[k]).delete();
				}
				System.out.println(Array.toStr(barnesPicked));
				for (int k = 0; k < barnesPicked.length; k++) {
					currentBarns.add(barnesPicked[k]);
				}

			}
			String finalVpop = matchDir + "barnacle.vpop";
			try {
				PrintWriter writer = new PrintWriter(new FileWriter(finalVpop));
				writer.println(Array.toStr(VcfPopulation.HEADER));
				for (int j = 0; j < currentBarns.size(); j++) {
					writer.println(currentBarns.get(j) + "\t" + VcfPopulation.CONTROL + "\t" + VcfPopulation.CONTROL);
				}
				for (String anchor : anchors) {
					writer.println(anchor + "\t" + VcfPopulation.CASE + "\t" + VcfPopulation.CASE);
				}
				writer.close();
			} catch (Exception e) {
				log.reportError("Error writing to " + finalVpop);
				log.reportException(e);
			}
		}

	}

	public static String getAppropriateRoot(String vcf, boolean removeDirectoryInfo) {
		String root = "";
		if (vcf.endsWith(VCF_EXTENSIONS.GZIP_VCF.getLiteral())) {
			StringBuilder b = new StringBuilder(vcf);
			b.replace(vcf.lastIndexOf(VCF_EXTENSIONS.GZIP_VCF.getLiteral()), vcf.lastIndexOf(VCF_EXTENSIONS.GZIP_VCF.getLiteral()) + VCF_EXTENSIONS.GZIP_VCF.getLiteral().length(), "");
			root = b.toString();

		} else {
			root = ext.rootOf(vcf, false);
		}
		if (removeDirectoryInfo) {
			root = ext.removeDirectoryInfo(root);
		}
		return root;
	}

	/**
	 * @param vcf
	 * @param idFile
	 * @param outputDir
	 * @param skipFiltered
	 * @param gzipOutput
	 * @param samples
	 *            optional if you would like to only keep a subset of samples
	 * 
	 * @param locusID
	 *            used the locus ID
	 * @param keepIDs
	 *            keep the ids in the file, otherwise remove them
	 * @param log
	 * @return
	 */
	public static String extractIDs(String vcf, String idFile, String outputDir, boolean skipFiltered, boolean gzipOutput, HashSet<String> samples, Segment[] segsToExclude, boolean locusID, boolean keepIDs, Logger log) {
		String outputVCF = null;
		if (idFile == null || !Files.exists(idFile)) {
			log.reportFileNotFound(idFile);
			return null;
		}
		if (vcf == null || !Files.exists(vcf)) {
			log.reportFileNotFound(vcf);
			return null;
		} else {
			String[] ids = HashVec.loadFileToStringArray(idFile, false, new int[] { 0 }, true);
			HashSet<String> tmp = new HashSet<String>();
			for (int i = 0; i < ids.length; i++) {
				tmp.add(ids[i]);
			}
			String dir = outputDir == null ? ext.parseDirectoryOfFile(vcf) : outputDir;
			new File(dir).mkdirs();
			String root = getAppropriateRoot(vcf, true);
			outputVCF = outputDir + root + "." + ext.rootOf(idFile) + ".vcf" + (gzipOutput ? ".gz" : "");
			if (!Files.exists(outputVCF)) {
				VCFFileReader reader = new VCFFileReader(vcf, true);
				VariantContextWriter writer = initWriter(outputVCF, DEFUALT_WRITER_OPTIONS, getSequenceDictionary(reader));
				copyHeader(reader, writer, samples == null ? BLANK_SAMPLE : samples, HEADER_COPY_TYPE.SUBSET_STRICT, log);
				int progress = 0;
				int found = 0;

				if (hasInfoLine(reader, "snp138") || locusID) {
					log.reportTimeWarning("If a variant has an ID of \".\", the" + (locusID ? " locusID" : " snp138 ") + "annotation will be added");

					for (VariantContext vc : reader) {
						progress++;
						if (progress % 100000 == 0) {
							log.reportTimeInfo(progress + " variants read...");
							log.reportTimeInfo(found + " variants found...");

						}
						// if (found == 10000) {
						// writer.close();
						// reader.close();
						// return outputVCF;
						// }
						String anno = locusID ? new VCOps.LocusID(vc).getId() : VCOps.getAnnotationsFor(new String[] { "snp138" }, vc, ".")[0];

						if ((!skipFiltered || !vc.isFiltered()) && (keepIDs && tmp.contains(anno)) || (!keepIDs && !tmp.contains(anno))) {
							if (segsToExclude == null || (keepIDs && Segment.binarySearchForAllOverLappingIndices(VCOps.getSegment(vc), segsToExclude) != null) || (!keepIDs && Segment.binarySearchForAllOverLappingIndices(VCOps.getSegment(vc), segsToExclude) == null)) {
								VariantContextBuilder builder = new VariantContextBuilder(vc);
								if (vc.getID().equals(".")) {
									builder.id(anno);
								}
								VariantContext vcAdd = builder.make();

								if (samples != null) {
									vcAdd = VCOps.getSubset(vcAdd, samples, VC_SUBSET_TYPE.SUBSET_STRICT, false);
								}
								writer.add(vcAdd);
								found++;
							}
						}
					}
				} else {
					log.reportTimeError("This method relies on the  \"snp138\" annotation, and none was detected, sorry");
				}
				log.reportTimeInfo(progress + " total variants read...");
				log.reportTimeInfo(found + " variants found...");
				reader.close();
				writer.close();

			} else {
				log.reportFileExists(outputVCF);
			}
		}
		return outputVCF;
	}

	public static String[][] getAnnotationKeys(String vcf, Logger log) {
		VCFFileReader reader = new VCFFileReader(vcf, false);
		String[] annotationKeys = new String[reader.getFileHeader().getInfoHeaderLines().size()];
		String[] descriptions = new String[reader.getFileHeader().getInfoHeaderLines().size()];
		int index = 0;
		for (VCFInfoHeaderLine vcfHeaderLine : reader.getFileHeader().getInfoHeaderLines()) {
			annotationKeys[index] = vcfHeaderLine.getID();
			descriptions[index] = vcfHeaderLine.getDescription();
			index++;
		}
		reader.close();
		return new String[][] { annotationKeys, descriptions };

	}

	/**
	 * just converts the contigs from number (b37 ex) to chr (hg10 ex)
	 */
	public static void b_ToHg(String input, String output, String refGenome, Logger log) {
		ReferenceGenome ref = new ReferenceGenome(refGenome, log);
		VCFFileReader reader = new VCFFileReader(input, false);
		VariantContextWriterBuilder builderWriter = new VariantContextWriterBuilder().setOutputFile(output);
		builderWriter.setReferenceDictionary(ref.getIndexedFastaSequenceFile().getSequenceDictionary());
		builderWriter.setOption(Options.DO_NOT_WRITE_GENOTYPES);
		VariantContextWriter writer = builderWriter.build();
		VCFHeader header = new VCFHeader(reader.getFileHeader());
		header.setSequenceDictionary(ref.getIndexedFastaSequenceFile().getSequenceDictionary());

		writer.writeHeader(header);
		for (VariantContext vc : reader) {
			VariantContextBuilder builder = new VariantContextBuilder(vc);
			builder.chr(Positions.getChromosomeUCSC(Integer.parseInt(vc.getContig()), true));
			builder.attributes(null);
			builder.noGenotypes();
			VariantContext vcNew = builder.make();
			writer.add(vcNew);

		}
		reader.close();
		writer.close();

	}

	public static void createSiteOnlyVcf(String inputVCF, String outputVCF, boolean overWrite, Logger log) {
		VCFFileReader reader = new VCFFileReader(inputVCF, false);
		if (Files.exists(outputVCF) && !overWrite) {
			log.reportTimeWarning(outputVCF + " exists, skipping site only");
		} else {
			VariantContextWriter writer = initWriter(outputVCF, new Options[] { Options.DO_NOT_WRITE_GENOTYPES }, reader.getFileHeader().getSequenceDictionary());
			copyHeader(reader, writer, null, HEADER_COPY_TYPE.SITE_ONLY, log);
			for (VariantContext vc : reader) {
				writer.add(vc);
			}
			reader.close();
			writer.close();
		}
	}

	public static String extractSegments(String vcf, String segmentFile, int bpBuffer, String bams, String outputDir, boolean skipFiltered, boolean gzipOutput, int numThreads, Logger log) {
		return extractSegments(vcf, segmentFile, bpBuffer, bams, outputDir, skipFiltered, gzipOutput, false, numThreads, log);
	}

	public static String extractSegments(String vcf, String segmentFile, int bpBuffer, String bams, String outputDir, boolean skipFiltered, boolean gzipOutput, boolean createAnnotationFile, int numThreads, Logger log) {
		return extractSegments(vcf, segmentFile, bpBuffer, bams, outputDir, skipFiltered, gzipOutput, createAnnotationFile, false, null, numThreads, log);
	}

	public static String extractSegments(String vcf, String segmentFile, int bpBuffer, String bams, String outputDir, boolean skipFiltered, boolean gzipOutput, boolean createAnnotationFile, boolean subToBam, String[] varSets, int numThreads, Logger log) {
		BamExtractor.BamSample bamSample = null;

		if (vcf == null || !Files.exists(vcf)) {
			log.reportFileNotFound(vcf);
			return null;
		}
		if (segmentFile != null && !Files.exists(segmentFile)) {
			log.reportFileNotFound(segmentFile);
			return null;
		}
		HashSet<String> bamSamples = new HashSet<String>();

		String dir = outputDir == null ? ext.parseDirectoryOfFile(vcf) : outputDir;
		new File(dir).mkdirs();
		String root = getAppropriateRoot(vcf, true);

		VCFFileReader reader = new VCFFileReader(vcf, Files.exists(vcf + ".idx"));
		if (bams == null) {
			log.reportTimeInfo("A bam directory was not provided, skipping bam extraction");
		} else {
			log.reportTimeInfo("A bam directory was provided, extracting bams to " + dir);
			if (Files.isDirectory(bams)) {
				bamSample = new BamExtractor.BamSample(Files.listFullPaths(bams, ".bam", false), log, true);
			} else {
				String[] tmp = null;
				if (bams.split(",").length > 1) {
					ArrayList<String> baTmp = new ArrayList<String>();
					String[] split = bams.split(",");
					for (int i = 0; i < split.length; i++) {
						if (Files.isDirectory(split[i])) {
							String[] dirBams = Files.listFullPaths(split[i], ".bam", false);
							for (int j = 0; j < dirBams.length; j++) {
								baTmp.add(dirBams[i]);
							}
						} else {
							String[] fileBams = HashVec.loadFileToStringArray(split[i], false, new int[] { 0 }, false);
							for (int j = 0; j < fileBams.length; j++) {
								baTmp.add(fileBams[i]);
							}
						}
					}
					tmp = Array.toStringArray(baTmp);
				} else {
					tmp = HashVec.loadFileToStringArray(bams, false, new int[] { 0 }, false);
				}
				log.reportTimeInfo("found " + tmp.length + " bam files in " + bams);
				bamSample = new BamExtractor.BamSample(tmp, log, true);
			}
			bamSample.generateMap();
			bamSample.getBamSampleMap();
			bamSample.verify(getSamplesInFile(reader), varSets);
			String[] vcfSamples = getSamplesInFile(reader);
			if (subToBam) {
				for (String abamSample : bamSample.getBamSampleMap().keySet()) {
					if (varSets == null) {
						bamSamples.add(abamSample);
					} else {
						for (int i = 0; i < varSets.length; i++) {
							String tmp = abamSample + varSets[i];
							if (ext.indexOfStr(tmp, vcfSamples) >= 0) {
								bamSamples.add(tmp);
							}
						}
					}
				}
			}
		}
		String output = dir + root + "." + ext.rootOf(segmentFile) + ".vcf" + (gzipOutput ? ".gz" : "");
		String annoFile = dir + root + "." + ext.rootOf(segmentFile) + ".anno.txt";
		if (!Files.exists(output)) {
			Segment[] segsToSearch = null;
			if (segmentFile.endsWith(".in") || segmentFile.endsWith(".bim")) {
				segsToSearch = Segment.loadRegions(segmentFile, 0, 3, 3, false);
			} else {
				segsToSearch = Segment.loadRegions(segmentFile, 0, 1, 2, 0, true, true, true, bpBuffer);
				log.reportTimeInfo("Loaded " + segsToSearch.length + " segments to search");
				segsToSearch = Segment.unique(segsToSearch);
				log.reportTimeInfo(segsToSearch.length + " were unique");
				segsToSearch = Segment.mergeOverlapsAndSortAllChromosomes(segsToSearch, 1);
				log.reportTimeInfo(segsToSearch.length + " after merging");

			}

			VariantContextWriter writer = initWriter(output, DEFUALT_WRITER_OPTIONS, getSequenceDictionary(reader));
			copyHeader(reader, writer, bamSamples, subToBam ? HEADER_COPY_TYPE.SUBSET_STRICT : HEADER_COPY_TYPE.FULL_COPY, log);
			int progress = 0;
			int found = 0;

			PrintWriter annoWriter = null;
			String[][] annotations = getAnnotationKeys(vcf, log);
			if (createAnnotationFile) {
				annoWriter = Files.getAppropriateWriter(annoFile);
				annoWriter.println("##" + Array.toStr(ANNO_BASE) + "\t" + Array.toStr(annotations[1]) + "\t" + Array.toStr(VCFOps.getSamplesInFile(vcf)));
				annoWriter.println(Array.toStr(ANNO_BASE) + "\t" + Array.toStr(annotations[0]) + "\t" + Array.toStr(VCFOps.getSamplesInFile(vcf)));
			}
			for (VariantContext vc : reader) {
				progress++;
				if (progress % 100000 == 0) {
					log.reportTimeInfo(progress + " variants read...");
					log.reportTimeInfo(found + " variants found...");

				}
				if ((!skipFiltered || !vc.isFiltered()) && VCOps.isInTheseSegments(vc, segsToSearch)) {
					if (subToBam) {
						writer.add(VCOps.getSubset(vc, bamSamples));
					} else {
						writer.add(vc);
					}
					if (bamSample != null) {
						bamSample.addSegmentToExtract(new Segment(Positions.chromosomeNumber(vc.getContig()), vc.getStart(), vc.getEnd()));
					}
					if (createAnnotationFile) {
						annoWriter.print(vc.getContig() + "\t" + vc.getStart() + "\t" + vc.getID() + "\t" + vc.getReference().getBaseString() + "\t" + vc.getAlternateAlleles().toString() + "\t" + vc.getHomRefCount() + "\t" + vc.getHetCount() + "\t" + vc.getHomVarCount() + "\t" + Array.toStr(VCOps.getAnnotationsFor(annotations[0], vc, ".")));
						GenotypesContext gc = vc.getGenotypes();

						for (Genotype g : gc) {
							annoWriter.print("\t" + g.getGenotypeString());
						}
						annoWriter.println();
					}
					found++;
				}
			}
			writer.close();

			if (bamSample != null) {
				BamExtractor.extractAll(bamSample, dir, bpBuffer, true, true, numThreads, log);
				bamSample = new BamExtractor.BamSample(Files.listFullPaths(dir, ".bam", false), log, true);
				bamSample.generateMap();
				bamSample.dumpToIGVMap(output, varSets);
			}
			if (createAnnotationFile) {
				annoWriter.close();
			}
		} else {
			log.reportTimeWarning("The file " + output + " already exists, skipping extraction step");
		}
		reader.close();
		return output;
	}

	public static int getNumberOfVariants(String vcf) {
		VCFFileReader vcfFileReader = new VCFFileReader(vcf, true);
		int numVar = 0;
		for (VariantContext vc : vcfFileReader) {
			vc.getID();
			numVar++;
		}
		vcfFileReader.close();
		return numVar;
	}

	/**
	 * @return true if the index exists or was created for all vcfs
	 */
	public static boolean verifyIndices(String[] vcfFiles, Logger log) {
		boolean verify = true;
		for (int i = 0; i < vcfFiles.length; i++) {
			if (!verifyIndex(vcfFiles[i], log)) {
				verify = false;
			}
		}
		return verify;
	}

	public static ChrSplitResults[] splitByChrs(String vcf, int numthreads, boolean onlyWithVariants, Logger log) {
		return splitByChrs(vcf, null, numthreads, onlyWithVariants, log);
	}

	public static ChrSplitResults[] splitByChrs(String vcf, String newDir, int numthreads, boolean onlyWithVariants, Logger log) {
		String[] toSplit = getAllContigs(vcf, log);
		log.reportTimeInfo("Detected " + toSplit.length + " chrs to split");
		VCFSplitProducer producer = new VCFSplitProducer(vcf, newDir, toSplit, log);
		WorkerTrain<ChrSplitResults> train = new WorkerTrain<VCFOps.ChrSplitResults>(producer, numthreads, numthreads, log);
		ArrayList<ChrSplitResults> chrSplitResults = new ArrayList<ChrSplitResults>();
		while (train.hasNext()) {
			ChrSplitResults tmp = train.next();
			if (onlyWithVariants) {
				if (tmp.hasVariants()) {
					chrSplitResults.add(tmp);
				} else {
					log.reportTimeWarning("skipping " + tmp.getChr() + " , no variants were found");
				}
			} else {
				chrSplitResults.add(tmp);
			}
		}
		return chrSplitResults.toArray(new ChrSplitResults[chrSplitResults.size()]);
	}

	private static class VCFSplitProducer implements Producer<ChrSplitResults> {
		private String vcfFile;
		private String[] toSplit;
		private Logger log;
		private int index;
		private String newDir;

		private VCFSplitProducer(String vcfFile, String newDir, String[] toSplit, Logger log) {
			super();
			this.vcfFile = vcfFile;
			this.toSplit = toSplit;
			this.log = log;
			this.newDir = newDir;
			this.index = 0;
		}

		@Override
		public boolean hasNext() {
			return index < toSplit.length;
		}

		@Override
		public Callable<ChrSplitResults> next() {
			final String chr = toSplit[index];
			String dir = newDir == null ? ext.parseDirectoryOfFile(vcfFile) : newDir + chr + "/";
			new File(dir).mkdirs();
			final String output = dir + getAppropriateRoot(vcfFile, true) + "." + chr + VCF_EXTENSIONS.GZIP_VCF.getLiteral();
			Callable<ChrSplitResults> callable = new Callable<VCFOps.ChrSplitResults>() {

				@Override
				public ChrSplitResults call() throws Exception {

					return splitByChr(vcfFile, output, chr, log);
				}

			};
			index++;
			return callable;
		}

		@Override
		public void remove() {

		}

		@Override
		public void shutdown() {

		}

	}

	public static String[] getAllContigs(String vcfFile, Logger log) {
		List<SAMSequenceRecord> sList = getAllSAMSequenceRecords(vcfFile, log);
		String[] all = new String[sList.size()];
		for (int i = 0; i < all.length; i++) {
			all[i] = sList.get(i).getSequenceName();
		}
		return all;
	}

	public static List<SAMSequenceRecord> getAllSAMSequenceRecords(String vcfFile, Logger log) {
		SAMSequenceDictionary samSequenceDictionary = getSequenceDictionary(vcfFile);
		List<SAMSequenceRecord> sList = samSequenceDictionary.getSequences();
		return sList;
	}

	private static boolean hasSomeVariants(String vcfFile, Logger log) {
		VCFFileReader reader = new VCFFileReader(vcfFile, false);

		int count = 0;
		for (VariantContext vc : reader) {
			vc.getContig();
			count++;
			if (count > 0) {
				reader.close();
				return true;
			}
		}
		reader.close();
		return false;
	}

	public static void dumpSnpEffGenes(String vcfFile, String geneFile, Logger log) {
		Hashtable<String, String> track = new Hashtable<String, String>();
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(geneFile));
			VCFFileReader reader = new VCFFileReader(vcfFile, true);
			for (VariantContext vc : reader) {
				String name = VCOps.getSNP_EFFGeneName(vc);
				if (!track.containsKey(name)) {
					writer.println(name);
					track.put(name, name);
				}
			}
			reader.close();
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + geneFile);
			log.reportException(e);
		}

	}

	public static ChrSplitResults splitByChr(String vcfFile, String outputVCF, String chr, Logger log) {
		int numChr = 0;
		int numTotal = 0;
		ChrSplitResults chrSplitResults = null;
		if (!Files.exists(outputVCF) || !verifyIndex(outputVCF, log)) {
			log.reportTimeInfo("Extracting chr " + chr + " from " + vcfFile + " to " + outputVCF);
			SAMSequenceDictionary samSequenceDictionary = getSequenceDictionary(vcfFile);
			VCFFileReader reader = new VCFFileReader(vcfFile, true);
			VariantContextWriter writer = initWriter(outputVCF, null, samSequenceDictionary);
			copyHeader(reader, writer, null, HEADER_COPY_TYPE.FULL_COPY, log);
			// samSequenceDictionary.getSequence(chr).g

			CloseableIterator<VariantContext> iterator = reader.query(chr, 1, samSequenceDictionary.getSequence(chr).getSequenceLength() + 100);

			while (iterator.hasNext()) {
				VariantContext vc = iterator.next();
				if (numTotal % 100000 == 0) {
					log.reportTimeInfo("Scanned " + numTotal + " variants from " + vcfFile + " (" + numChr + " matching " + chr + ")");
				}
				if (vc.getContig().equals(chr)) {
					writer.add(vc);
					numChr++;
				} else {
					throw new IllegalStateException("Iterator returned invalid sequence...");
				}
				numTotal++;
			}

			writer.close();
			reader.close();
			chrSplitResults = new ChrSplitResults(chr, vcfFile, outputVCF, numChr);
		} else {
			log.reportFileExists(outputVCF);
			log.reportFileExists(chr);

			chrSplitResults = new ChrSplitResults(chr, vcfFile, outputVCF, numChr);
			chrSplitResults.setHasVariants(hasSomeVariants(outputVCF, log));
			chrSplitResults.setSkippedExists(true);
		}
		return chrSplitResults;
	}

	public static class ChrSplitResults {
		private String inputVCF;
		private String outputVCF;
		private String chr;
		private int numChr;
		private boolean skippedExists;
		private boolean hasVariants;

		public ChrSplitResults(String chr, String inputVCF, String outputVCF, int numChr) {
			super();
			this.chr = chr;
			this.inputVCF = inputVCF;
			this.outputVCF = outputVCF;
			this.numChr = numChr;
			this.skippedExists = false;
			this.hasVariants = numChr > 0;
		}

		public String getChr() {
			return chr;
		}

		public String getOutputVCF() {
			return outputVCF;
		}

		public boolean isSkippedExists() {
			return skippedExists;
		}

		public boolean hasVariants() {
			return hasVariants;
		}

		public void setHasVariants(boolean hasVariants) {
			this.hasVariants = hasVariants;
		}

		public void setSkippedExists(boolean skippedExists) {
			this.skippedExists = skippedExists;
		}

		public String getSummary() {
			String summary = "Extracted chr " + chr + " from " + inputVCF + " to " + outputVCF;
			summary += "Found " + numChr + " total variants on " + chr;
			return summary;
		}
	}

	/**
	 * @return true if the index exists and was valid, or was created
	 */
	public static boolean verifyIndex(String vcfFile, Logger log) {
		boolean created = false;
		if (!vcfFile.endsWith(VCF_EXTENSIONS.REG_VCF.getLiteral()) && !vcfFile.endsWith(VCF_EXTENSIONS.GZIP_VCF.getLiteral())) {
			log.reportTimeError("We currently can only index the following extensions");
			log.reportTimeError(VCF_EXTENSIONS.REG_VCF.getLiteral());
			log.reportTimeError(VCF_EXTENSIONS.GZIP_VCF.getLiteral());

		} else {
			if (vcfFile.endsWith(VCF_EXTENSIONS.REG_VCF.getLiteral())) {
				File indexFile = Tribble.indexFile(new File(vcfFile));
				if (indexFile.canRead()) {
					log.report("Info - Loading index file " + indexFile);
					IndexFactory.loadIndex(indexFile.getAbsolutePath());
					created = true;
				} else {
					log.report("Info - creating index file " + indexFile);
					try {
						Index index = IndexFactory.createLinearIndex(new File(vcfFile), new VCFCodec());
						LittleEndianOutputStream stream = new LittleEndianOutputStream(new FileOutputStream(indexFile));
						index.write(stream);
						stream.close();
						created = true;
					} catch (IOException e) {
						log.reportError("Error - could not create index file " + indexFile);
						created = false;
					}
				}
			} else if (vcfFile.endsWith(VCF_EXTENSIONS.GZIP_VCF.getLiteral())) {
				// return false;
				String indexFile = vcfFile + TabixUtils.STANDARD_INDEX_EXTENSION;
				if (Files.exists(indexFile)) {
					created = true;
				} else {
					log.reportTimeWarning("Indexing not quite implemented yet for " + VCF_EXTENSIONS.GZIP_VCF.getLiteral() + ", exiting");
					// System.exit(1);
					VCFFileReader readerVcfGz = new VCFFileReader(vcfFile, false);
					TabixIndex index = IndexFactory.createTabixIndex(new File(vcfFile), new VCFCodec(), TabixFormat.VCF, readerVcfGz.getFileHeader().getSequenceDictionary());
					try {
						index.writeBasedOnFeatureFile(new File(vcfFile));
					} catch (IOException e) {
						e.printStackTrace();
					}
					created = false;
					readerVcfGz.close();

				}
				// // TabixIndexCreator idxCreator = new TabixIndexCreator( getSequenceDictionary(vcfFile), TabixFormat.VCF);
				// // VCFFileReader reader = new VCFFileReader(vcfFile, false);
				// // for(VariantContext vc:reader){
				// // idxCreator.addFeature(vc, 1);
				// // }
				//
				// // initWriter(output, options, sequenceDictionary)
				// TabixIndex tbi = IndexFactory.createTabixIndex(new File(vcfFile), new VCFCodec(), TabixFormat.VCF, getSequenceDictionary(vcfFile));
				//
				//
				// VCFFileReader reader = new VCFFileReader(vcfFile, false);
				//

				// }
			}
		}
		return created;
	}

	/**
	 * @return true if the vcf file was gzipped
	 */
	public static String gzipAndIndex(String vcfFile, Logger log) {
		String vcfFileGz = vcfFile + ".gz";

		if (Files.exists(vcfFileGz)) {
			log.reportTimeWarning("Gzipped vcf " + vcfFileGz + " already exists, skipping");
		} else {
			if (verifyIndex(vcfFile, log)) {
				VCFFileReader reader = new VCFFileReader(vcfFile, true);
				VariantContextWriter writer = initWriter(vcfFileGz, null, reader.getFileHeader().getSequenceDictionary());
				copyHeader(reader, writer, null, HEADER_COPY_TYPE.FULL_COPY, log);
				for (VariantContext vc : reader) {
					writer.add(vc);
				}
				reader.close();
				writer.close();

			}
		}
		return vcfFileGz;
	}

	public static VariantContext lookupExactVariant(String vcf, VariantContext vc, Logger log) {
		VCFFileReader reader = new VCFFileReader(vcf, true);
		Segment vcSeg = VCOps.getSegment(vc);
		CloseableIterator<VariantContext> cIterator = reader.query(vcSeg.getChromosomeUCSC(), vcSeg.getStart(), vcSeg.getStop());
		VariantContext vcMatch = null;
		while (cIterator.hasNext()) {
			VariantContext vcTmp = cIterator.next();
			if (vcSeg.equals(VCOps.getSegment(vcTmp))) {
				if (vcTmp.getReference().equals(vc.getReference(), false) && vcTmp.hasSameAlternateAllelesAs(vc)) {
					vcMatch = vcTmp;
					break;
				}
			}
		}
		reader.close();
		return vcMatch;
	}

	/**
	 * Creates a new vcf with filtered variants removed
	 */
	public static String removeFilteredVariants(String vcf, boolean gzipOutput, boolean standardFilters, Logger log) {
		VCFFileReader reader = new VCFFileReader(vcf, true);
		String output = ext.addToRoot(vcf.endsWith(VCF_EXTENSIONS.GZIP_VCF.getLiteral()) ? vcf.replaceAll(".gz", "") : vcf, ".filtered") + (gzipOutput ? ".gz" : "");

		output = getAppropriateRoot(vcf, false) + ".filtered" + (gzipOutput ? VCF_EXTENSIONS.GZIP_VCF.getLiteral() : VCF_EXTENSIONS.REG_VCF.getLiteral());
		log.reportTimeInfo("Will write filtered variants to " + output);
		VariantContextWriter writer = initWriter(output, DEFUALT_WRITER_OPTIONS, getSequenceDictionary(reader));
		VCFOps.copyHeader(reader, writer, null, HEADER_COPY_TYPE.FULL_COPY, log);
		DynamicHistogram dyHistogramVQSLOD = null;
		if (reader.getFileHeader().hasInfoLine("VQSLOD")) {
			log.reportTimeInfo("Detected info header line for VQSLOD, creating a histogram of scores for passing variants");
			dyHistogramVQSLOD = new DynamicHistogram(0, 100, 0);
		}
		VARIANT_FILTER_DOUBLE[] vDoubles = new VARIANT_FILTER_DOUBLE[0];
		if (standardFilters) {
			VARIANT_FILTER_DOUBLE gq = VARIANT_FILTER_DOUBLE.GQ_STRICT;
			VARIANT_FILTER_DOUBLE dp = VARIANT_FILTER_DOUBLE.DP;
			VARIANT_FILTER_DOUBLE maf = VARIANT_FILTER_DOUBLE.MAF;
			VARIANT_FILTER_DOUBLE cr = VARIANT_FILTER_DOUBLE.CALL_RATE;

			vDoubles = new VARIANT_FILTER_DOUBLE[] { cr, maf, gq, dp };
		}
		VariantContextFilter variantContextFilter = new VariantContextFilter(vDoubles, new VARIANT_FILTER_BOOLEAN[] { VARIANT_FILTER_BOOLEAN.FAILURE_FILTER }, null, null, log);
		int count = 0;
		int countFiltered = 0;
		int countPassed = 0;
		for (VariantContext vc : reader) {
			count++;
			if (variantContextFilter.filter(vc).passed()) {
				writer.add(vc);
				if (dyHistogramVQSLOD != null) {
					dyHistogramVQSLOD.addDataPointToHistogram(vc.getCommonInfo().getAttributeAsDouble("VQSLOD", 0.0));
				}
				countPassed++;
			} else {
				countFiltered++;
			}
			if (count % 100000 == 0) {
				log.reportTimeInfo(count + " total variants, " + countPassed + " passed the filters, " + countFiltered + " were filtered");
			}
		}
		reader.close();
		writer.close();
		log.reportTimeInfo(count + " total variants read...");
		log.reportTimeInfo(countPassed + " variants passed the filters...");
		if (dyHistogramVQSLOD != null) {
			String outputHist = ext.addToRoot(output, ".hist.VQSLOD");

			try {
				PrintWriter writerHist = new PrintWriter(new FileWriter(outputHist));
				double[] bins = dyHistogramVQSLOD.getBins();
				int[] counts = dyHistogramVQSLOD.getCounts();
				writerHist.println("VQSLOD_BIN\tVQSLOD_COUNT");
				for (int i = 0; i < counts.length; i++) {
					writerHist.println(bins[i] + "\t" + counts[i]);
				}
				writerHist.close();
			} catch (Exception e) {
				log.reportError("Error writing to " + outputHist);
				log.reportException(e);
			}
		}
		return output;
	}

	public static void qcVCF(String vcf, Logger log) {
		// verifyIndex(vcf, log);
		PlinkSeq plinkSeq = new PlinkSeq(false, true, log);
		PseqProject pseqProject = PlinkSeq.initialize(plinkSeq, ext.rootOf(vcf), null, vcf, null, ext.parseDirectoryOfFile(vcf), true, true, log);
		// plinkSeq.createNewProject(pseqProject);
		plinkSeq.loadData(pseqProject, LOAD_TYPES.VCF, null);
		WorkerHive<PlinkSeqWorker> assocHive = new WorkerHive<PlinkSeq.PlinkSeqWorker>(1, 10, log);
		assocHive.addCallable(PlinkSeq.generateAWorker(pseqProject, ANALYSIS_TYPES.I_SUMMARY, null, null, null, null, -1, "0", pseqProject.getProjectName(), true, log));
		assocHive.execute(true);
	}

	/**
	 * @param vcf
	 *            Samples from this vcf will be dumped to a {@link VcfPopulation} formatted file
	 * @param log
	 */
	public static void extractSamps(String vcf, Logger log) {
		VCFFileReader reader = new VCFFileReader(vcf, false);
		String[] samples = getSamplesInFile(reader);
		reader.close();
		String output = getAppropriateRoot(vcf, false) + ".vpop";
		if (!Files.exists(output)) {
			log.reportTimeInfo("Detected " + samples.length + " samples in vcf " + vcf);
			log.reportTimeInfo("exporting to " + output);
			try {
				PrintWriter writer = new PrintWriter(new FileWriter(output));
				writer.println(Array.toStr(VcfPopulation.HEADER));
				for (int i = 0; i < samples.length; i++) {
					writer.println(samples[i] + "\t" + VcfPopulation.CONTROL + "\t" + VcfPopulation.CONTROL);
				}
				writer.close();
			} catch (Exception e) {
				log.reportError("Error writing to " + output);
				log.reportException(e);
			}

		} else {
			log.reportFileExists(output);
		}
	}

	/**
	 * @param segs
	 *            Genvisis type segments to search
	 * @param vcfHeader
	 *            an {@link VCFHeader}
	 * @param bpBuffer
	 *            bp buffer to be added to the segments
	 * @param optimize
	 *            perform an extra step to make a more efficient query
	 * @param log
	 * @return array of {@link QueryInterval} that can be queried by a bamfile reader
	 */
	public static QueryInterval[] convertSegsToQI(Segment[] segs, VCFHeader vcfHeader, int bpBuffer, boolean optimize, Logger log) {
		QueryInterval[] qIntervals = new QueryInterval[segs.length];
		segs = Segment.sortSegments(segs);
		for (int i = 0; i < qIntervals.length; i++) {
			String sequenceName = Positions.getChromosomeUCSC(segs[i].getChr(), true);

			int referenceIndex = vcfHeader.getSequenceDictionary().getSequenceIndex(sequenceName);
			if (referenceIndex < 0) {
				log.reportError("Error - could not find " + sequenceName + " in the sequence dictionary, halting");
				return null;
			}
			int start = segs[i].getStart();
			int stop = segs[i].getStop();
			if (segs[i].getChr() == 0 && segs[i].getStart() <= 0) {
				start = 1;
				if (segs[i].getStop() <= 0) {
					stop = 1;
				}
			}
			qIntervals[i] = new QueryInterval(referenceIndex, start - bpBuffer, stop + bpBuffer);
		}
		if (optimize) {
			qIntervals = QueryInterval.optimizeIntervals(qIntervals);
		}
		return qIntervals;
	}

	public static void annoWithBed(String vcfFile, String bedFile) {
		System.out.println(bedFile);

		BEDFileReader bedReader = new BEDFileReader(bedFile, false);
		// LocusSet<BEDFeatureSeg> bLocusSet = bedReader.loadAll(new Logger());
		System.out.println(bedFile);
		String key = ext.rootOf(bedFile);
		String description = "Custom annotation provide by genvisis from " + ext.removeDirectoryInfo(bedFile) + " on " + ext.getTimestampForFilename();
		VCFFileReader reader = new VCFFileReader(vcfFile, false);
		String outputVcf = VCFOps.getAppropriateRoot(vcfFile, false) + "_" + ext.rootOf(bedFile) + VCFOps.VCF_EXTENSIONS.GZIP_VCF.getLiteral();
		VariantContextWriter writer = VCFOps.initWriter(outputVcf, DEFUALT_WRITER_OPTIONS, null);
		VCFHeader vcfHeader = reader.getFileHeader();
		VCFInfoHeaderLine vHeaderLine = new VCFInfoHeaderLine(key, 1, VCFHeaderLineType.String, description);
		vcfHeader.addMetaDataLine(vHeaderLine);
		writer.writeHeader(vcfHeader);
		ArrayList<Segment> segs = new ArrayList<Segment>();

		for (VariantContext vc : reader) {
			LocusSet<BEDFeatureSeg> bLocusSet = bedReader.loadSegsFor(VCOps.getSegment(vc), new Logger());
			BEDFeatureSeg[] olaps = bLocusSet.getOverLappingLoci(VCOps.getSegment(vc));
			segs.add(VCOps.getSegment(vc));
			VariantContextBuilder builder = new VariantContextBuilder(vc);
			builder.attribute(key, ".");
			if (olaps != null && olaps.length > 0) {
				if (olaps.length > 1) {
					reader.close();
					bedReader.close();
					throw new IllegalStateException("Bed file contained multiple overlapping segments for " + vc.toStringWithoutGenotypes());
				} else {
					builder.attribute(key, olaps[0].getBedFeature().getName());
				}
			}
			writer.add(builder.make());
		}
		bedReader.close();
		reader.close();
		writer.close();
		LocusSet<Segment> segSet = new LocusSet<Segment>(segs.toArray(new Segment[segs.size()]), true, new Logger()) {

			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

		};
		String outSegs = VCFOps.getAppropriateRoot(vcfFile, false) + "_" + ext.rootOf(bedFile) + ".segs";
		segSet.writeRegions(outSegs, TO_STRING_TYPE.REGULAR, false, new Logger());
		extractSegments(outputVcf, outSegs, 0, null, ext.parseDirectoryOfFile(outSegs), false, true, true, 1, new Logger());
	}

	private static String[] getExtractAnnotationCommand(UTILITY_TYPE type, Logger log) {
		ArrayList<String> params = new ArrayList<String>();
		params.add("#Full path to a vcf (preferably indexed");
		params.add(VCF_COMMAND);
		params.add("#Full path to a .bed file to extract");
		params.add(SEGMENT_FILE_COMMAND);
		params.add("#Full path to a either a file listing .bams, or a directory of .bam files");
		params.add("#" + BAM_COMMAND);
		System.out.println(Array.toStr(Array.toStringArray(params)));
		return Array.toStringArray(params);
	}

	public static void fromParameters(String filename, UTILITY_TYPE type, Logger log) {
		Vector<String> params = new Vector<String>();
		switch (type) {
		case CONVERT_PLINK:
			log.reportTimeInfo("Invalid op type " + type);
			break;
		case DUMP_SAMPLES:
			log.reportTimeInfo("Invalid op type " + type);

			break;
		case EXTRACT_IDS:
			log.reportTimeInfo("Invalid op type " + type);

			break;
		case EXTRACT_SEGMENTS:
			log.reportTimeInfo("Invalid op type " + type);

			break;
		case EXTRACT_SEGMENTS_ANNOTATION:
			params = Files.parseControlFile(filename, COMMAND_VCF_OPS_EXTRACT, getExtractAnnotationCommand(type, log), log);

			break;
		case GZIP:
			log.reportTimeInfo("Invalid op type " + type);

			break;
		case HOMOGENEITY:
			log.reportTimeInfo("Invalid op type " + type);

			break;
		case QC:
			break;
		case REMOVE_FILTERED:
			log.reportTimeInfo("Invalid op type " + type);

			break;
		case SUBSET_SUPER:
			log.reportTimeInfo("Invalid op type " + type);

			break;
		default:
			break;
		}
		if (params != null) {
			params.add(UTILITY_COMMAND + type);

			main(Array.toStringArray(params));
		}
	}

	private static final String VCF_COMMAND = "vcf=";
	private static final String UTILITY_COMMAND = "utility=";
	private static final String SEGMENT_FILE_COMMAND = "segs=";
	private static final String BAM_COMMAND = "bams=";
	public static final String COMMAND_VCF_OPS_EXTRACT = "vcfExtract";
	public static final String COMMAND_VCF_EXTRACT_DESCRIPTION = "extract a file of segments from a vcf and bam files";

	public static void main(String[] args) {
		int numArgs = args.length;
		String vcf = "Avcf.vcf";
		String populationFile = null;
		String logfile = null;
		int bpBuffer = 0;
		UTILITY_TYPE type = UTILITY_TYPE.CONVERT_PLINK;
		String segmentFile = null;
		String idFile = null;
		String bams = null;
		String outDir = null;
		String bedFile = null;
		boolean skipFiltered = false;
		boolean standardFilters = false;
		boolean gzip = true;
		boolean removeMonoMorphic = false;
		boolean keepIds = false;
		int numThreads = 1;
		Logger log;
		boolean subToBam = false;
		String[] varSet = null;

		String usage = "\n" + "seq.analysis.VCFUtils requires 0-1 arguments\n";
		usage += "   (1) full path to a vcf file (i.e. " + VCF_COMMAND + vcf + " (default))\n" + "";
		usage += "   (2) utility type (i.e. " + UTILITY_COMMAND + type + " (default))\n" + "";
		usage += "   (3) full path to a file (can be comma delimited for homogeneity utility) defining a population for the vcf (i.e. pop= (no default))\n" + "";
		usage += "   (4) the type of vcf extension (i.e. pop= (no default))\n" + "";
		usage += "   (5) full path to a file name with chr,start,stop or *.bim to extract (i.e. segs= (no default))\n" + "";
		usage += "   (6) bp buffer for segments to extract (i.e. bp=" + bpBuffer + "(default))\n" + "";
		usage += "   (7) a bam directory to extract associtated reads (i.e. " + BAM_COMMAND + bams + "( no default))\n" + "";
		usage += "   (8) an output directory for extracted vcfs/minibams (i.e. outDir=" + outDir + "( no default))\n" + "";
		usage += "   (9) skip filtered variants when extracting (i.e. -skipFiltered (not the default))\n" + "";
		usage += "   (10) gzip the output when extracting (i.e. -gzip ( the default))\n" + "";
		usage += "   (11) full path to a file of ids (i.e. idFile= (no default))\n" + "";
		usage += "   (12) when removing filtered variants, apply our standard filters as well (i.e. -standardFilters (not the default, GQ >=" + VARIANT_FILTER_DOUBLE.GQ.getDFilter() + " and DP >=" + VARIANT_FILTER_DOUBLE.DP.getDFilter() + "))\n" + "";
		usage += "   (13) when subsetting by samples, remove monomorphic variants (i.e. -removeMonoMorphic (not the default))\n" + "";
		usage += "   (14) when subsetting, keep variant ids if set to \".\" (i.e. -keepIds (not the default))\n" + "";
		usage += "   (15) full path to a segment file  (i.e. " + SEGMENT_FILE_COMMAND + "( no default)\n" + "";
		usage += "   (16) full path to a bed file for annotating  (i.e. bed= (no default)\n" + "";
		usage += "   (17) subset samples to bam files in vcf (i.e. subToBam=false (no default)\n" + "";
		usage += "   (18) varset for sample matching in vcf (i.e. varSet= (no default)\n" + "";

		usage += PSF.Ext.getNumThreadsCommand(14, numThreads);
		usage += "   NOTE: available utilities are:\n";

		for (int i = 0; i < UTILITY_TYPE.values().length; i++) {
			usage += UTILITY_TYPE.values()[i] + "\n";
		}
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith(VCF_COMMAND)) {
				vcf = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith(UTILITY_COMMAND)) {
				type = UTILITY_TYPE.valueOf(ext.parseStringArg(args[i], ""));
				numArgs--;
			} else if (args[i].startsWith("vpopFile=")) {
				populationFile = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("segs=")) {
				segmentFile = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("varSet=")) {
				varSet = ext.parseStringArg(args[i], "").split(",");
				numArgs--;
			} else if (args[i].startsWith("subToBam=")) {
				subToBam = ext.parseBooleanArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("idFile=")) {
				idFile = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith(BAM_COMMAND)) {
				bams = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("outDir=")) {
				outDir = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("bp=")) {
				bpBuffer = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("bed=")) {
				bedFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("-skipFiltered")) {
				skipFiltered = true;
				numArgs--;
			} else if (args[i].startsWith("-standardFilters")) {
				standardFilters = true;
				numArgs--;
			} else if (args[i].startsWith("-removeMonoMorphic")) {
				removeMonoMorphic = true;
				numArgs--;
			} else if (args[i].startsWith("-keepIds")) {
				keepIds = true;
				numArgs--;
			} else if (args[i].startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
				numThreads = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("-gzip")) {
				gzip = true;
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
			log.reportTimeInfo("Running utiltity type: " + type);

			switch (type) {
			case CONVERT_PLINK:
				convertToPlinkSet(null, vcf, "plink", PLINK_SET_MODE.GWAS_QC, log);
				break;
			case SUBSET_SUPER:
				VcfPopulation.splitVcfByPopulation(vcf, populationFile, removeMonoMorphic, keepIds, log);
				break;
			case EXTRACT_SEGMENTS:
				extractSegments(vcf, segmentFile, bpBuffer, bams, outDir, skipFiltered, gzip, false, subToBam, varSet, numThreads, log);
				break;
			case EXTRACT_SEGMENTS_ANNOTATION:
				extractSegments(vcf, segmentFile, bpBuffer, bams, outDir, skipFiltered, gzip, true, numThreads, log);
				break;
			case REMOVE_FILTERED:
				removeFilteredVariants(vcf, gzip, standardFilters, log);
				break;
			case GZIP:
				gzipAndIndex(vcf, log);
				break;
			case QC:
				qcVCF(vcf, log);
				break;
			case EXTRACT_IDS:
				extractIDs(vcf, idFile, outDir, skipFiltered, gzip, null, null, false, true, log);
				break;
			case DUMP_SAMPLES:
				extractSamps(vcf, log);
				break;
			case HOMOGENEITY:
				runHomoGeneity(vcf, populationFile.split(","), log);
				break;
			case BED_ANNOTATE:
				annoWithBed(vcf, bedFile);
				break;
			default:
				System.err.println("Invalid utility type: Available are ->");
				for (int i = 0; i < UTILITY_TYPE.values().length; i++) {
					usage += UTILITY_TYPE.values()[i] + "\n";
				}
				break;
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
