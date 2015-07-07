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
import java.util.TreeSet;
import java.util.concurrent.Callable;

import seq.analysis.PlinkSeq;
import seq.analysis.PlinkSeq.ANALYSIS_TYPES;
import seq.analysis.PlinkSeq.LOAD_TYPES;
import seq.analysis.PlinkSeq.PlinkSeqWorker;
import seq.analysis.PlinkSeqUtils.PseqProject;
import seq.manage.VCOps.VC_SUBSET_TYPE;
import seq.qc.FilterNGS.VARIANT_FILTER_BOOLEAN;
import seq.qc.FilterNGS.VARIANT_FILTER_DOUBLE;
import seq.qc.FilterNGS.VariantContextFilter;
import stats.Histogram.DynamicAveragingHistogram;
import stats.Histogram.DynamicHistogram;
import filesys.Segment;
import gwas.MergeDatasets;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.util.LittleEndianOutputStream;
import htsjdk.tribble.util.TabixUtils;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
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
	public static final Set<String> BLANK_SAMPLE = new TreeSet<String>();
	public static final Options[] DEFUALT_WRITER_OPTIONS = new Options[] { Options.INDEX_ON_THE_FLY };

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

	private enum UTILITY_TYPE {
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
		 * Determine Homogeneity for populations in a vcf
		 */
		HOMOGENEITY;
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

	public static String reportCallRateHWEFiltered(String vcf, String outputFile, double callRate, double hweP, Logger log) {
		VARIANT_FILTER_DOUBLE callRateFilter = VARIANT_FILTER_DOUBLE.CALL_RATE_LOOSE;
		VARIANT_FILTER_DOUBLE hwe = VARIANT_FILTER_DOUBLE.HWE;
		hwe.setDFilter(hweP);
		callRateFilter.setDFilter(callRate);
		VariantContextFilter filter = new VariantContextFilter(new VARIANT_FILTER_DOUBLE[] { callRateFilter, hwe }, new VARIANT_FILTER_BOOLEAN[] {}, null, null, log);
		DynamicAveragingHistogram dynamicAveragingHistogram = new DynamicAveragingHistogram(0, 1.1, 2);

		try {
			PrintWriter writer = new PrintWriter(new FileWriter(outputFile));
			VCFFileReader reader = new VCFFileReader(vcf, true);
			int count = 0;
			for (VariantContext vc : reader) {
				count++;
				if (count % 100000 == 0) {
					log.reportTimeInfo(count + " variants tested at callrate " + callRate);
				}
				if (!filter.filter(vc).passed()) {
					writer.println((vc.getID().equals(".") ? new VCOps.LocusID(vc).getId() : vc.getID()) + "\t" + vc.getNoCallCount() + "\t" + vc.getNSamples());
				}
				double g1000 = -1;
				try {
					g1000 =Double.parseDouble(VCOps.getAnnotationsFor(new String[] { "g10002014oct_all" }, vc, ".")[0]);
				} catch (NumberFormatException nfe) {

				}

				dynamicAveragingHistogram.addDataPair(VCOps.getCallRate(vc, null), g1000);
			}
			reader.close();
			writer.close();
			String hist = ext.addToRoot(outputFile, ".callrateHist");
			PrintWriter writerHist = new PrintWriter(new FileWriter(hist));
			writerHist.println("CallRateBin\tCount\t1000GFreq");
			dynamicAveragingHistogram.average();
			for (int i = 0; i < dynamicAveragingHistogram.getBins().length; i++) {
				writerHist.println(dynamicAveragingHistogram.getBins()[i] + "\t" + dynamicAveragingHistogram.getCounts()[i] + "\t" + dynamicAveragingHistogram.getAverages()[i]);
			}
			writerHist.close();
			return hist;
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
			Hashtable<String, String> newIDS = new Hashtable<String, String>();
			newIDS = fixFamFile(log, outFiles[2]);
			log.reportTimeInfo("MODE=" + mode);
			if (mode == PLINK_SET_MODE.GWAS_QC) {

				gwas.Qc.fullGamut(dir, false, new Logger(dir + "fullGamutOfMarkerAndSampleQC.log"));
				String mdsFile = dir + "genome/mds20.mds";
				if (Files.exists(mdsFile)) {
					fixMdsFile(log, dir, newIDS, mdsFile);
					CmdLine.run("runEigenstratWoHapMap", dir);
					CmdLine.run("runEigenstrat2", dir + "ancestry/");
					// fixMdsFile(log, dir + "ancestry/", newIDS, combo_fancy_postnormed_eigens.xln);
				}
			} else if (mode == PLINK_SET_MODE.HOMOGENEITY) {

				if (!Files.exists(ext.parseDirectoryOfFile(outFiles[2]) + "hardy.hwe")) {
					CmdLine.run("plink --bfile plink --maf 0 --geno 1 --mind 1 --hardy --out hardy --noweb", ext.parseDirectoryOfFile(outFiles[2]));
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

	private static void fixMdsFile(Logger log, String dir, Hashtable<String, String> newIDS, String mdsFile) {
		if (newIDS.size() > 0 && Files.exists(mdsFile)) {
			String[] header = Files.getHeaderOfFile(mdsFile, log);
			int[] indices = new int[header.length];
			for (int i = 0; i < indices.length; i++) {
				indices[i] = i;
			}
			String[][] mds = HashVec.loadFileToStringMatrix(mdsFile, false, new int[] { 0, 1, 2, 3, 4, 5 }, false);
			String[][] newMds = new String[mds.length][mds[0].length];
			for (int i = 0; i < mds.length; i++) {
				for (int j = 0; j < mds[i].length; j++) {
					if (j == 0 && j != 1 && newIDS.containsKey(mds[i][0])) {
						newMds[i][0] = newIDS.get(mds[i][0]);
						newMds[i][1] = newIDS.get(mds[i][0]);
					} else {
						newMds[i][j] = mds[i][j];
					}
				}
			}
			Files.writeMatrix(newMds, ext.addToRoot(mdsFile, ".fixed"), "\t");

		}
	}

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
		public static final String CASE = "CASE";
		public static final String CONTROL = "CONTROL";
		public static final String EXCLUDE = "EXCLUDE";
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
			CASE_CONTROL, ANY, STRATIFICATION, EXOME_DEPTH, PC_ANCESTRY;
		}

		public enum RETRIEVE_TYPE {
			SUPER, SUB;
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

		public static String[] splitVcfByPopulation(String vcf, String fullPathToPopFile, Logger log) {
			if (vcf != null && !Files.exists(vcf)) {
				log.reportFileNotFound(vcf);
				return null;
			}
			if (fullPathToPopFile != null && !Files.exists(fullPathToPopFile)) {
				log.reportFileNotFound(fullPathToPopFile);
				return null;
			}

			log.reportTimeWarning("Variant context ids that are set to \".\" will be updated");
			VcfPopulation vpop = VcfPopulation.load(fullPathToPopFile, POPULATION_TYPE.ANY, log);
			vpop.report();
			VCFFileReader reader = new VCFFileReader(vcf, true);

			String dir = ext.parseDirectoryOfFile(vcf);
			String root = getAppropriateRoot(vcf, true);

			VariantContextWriter[] writers = vpop.getWritersForPop(dir + root, reader, log);
			if (writers != null) {

				int progress = 0;
				for (VariantContext vc : reader) {
					progress++;
					if (progress % 100000 == 0) {
						log.reportTimeInfo(progress + " variants read...");
					}
					if (vc.getID().equals(".")) {
						VariantContextBuilder builder = new VariantContextBuilder(vc);
						builder.id(new VCOps.LocusID(vc).getId());
						vc = builder.make();
					}
					for (int i = 0; i < writers.length; i++) {
						if (writers[i] != null) {
							VariantContext vcSub = VCOps.getSubset(vc, vpop.getSuperPop().get(vpop.getUniqSuperPop().get(i)));

							// if (vcSub.getHomVarCount() > 0 || vcSub.getHetCount() > 0) {
							writers[i].add(vcSub);
						}
						// System.out.println(vpop.getFileNamesForPop(dir + root, log)[i]);
						// }
					}
				}
				for (int i = 0; i < writers.length; i++) {
					if (writers[i] != null) {
						writers[i].close();
					}
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
		String[] toRemove = new String[] {};
		double callRate = 0.80;
		double hwe = .00001;
		String finalSamples = vcf + ".finalSamples";
		if (!Files.exists(finalSamples)) {
			log.reportTimeError("Required file " + finalSamples + " is missing");
			return;
		}
		String[] samples = HashVec.loadFileToStringArray(finalSamples, false, new int[] { 0 }, false);
		log.reportTimeInfo("Found " + samples.length + " samples for the final analysis");
		for (int i = 0; i < fullPathToPopFiles.length; i++) {
			String fullPathToPopFile = fullPathToPopFiles[i];
			String[] splits = VcfPopulation.splitVcfByPopulation(vcf, fullPathToPopFile, log);
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
					String[] callRateRemove = HashVec.loadFileToStringArray(callRateFiltered, false, new int[] { 0 }, true);
					log.reportTimeInfo(callRateRemove.length + " variants removed from " + splits[j] + " at callrate " + callRateFiltered);
					// toRemove = Array.unique(Array.concatAll(toRemove, callRateRemove));
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

			if (Files.exists(problems)) {
				String[] problematic = HashVec.loadFileToStringArray(problems, false, new int[] { 0 }, true);
				log.reportTimeInfo(problematic.length + " markers with problems from " + fullPathToPopFiles[i]);
				toRemove = Array.unique(Array.concatAll(toRemove, problematic));
			}

		}

		String finalDir = ext.parseDirectoryOfFile(vcf) + "homogeneity/";
		new File(finalDir).mkdirs();
		String idFile = finalDir + "variants_Removed.txt";
		Files.writeList(toRemove, idFile);
		HashSet<String> sampleHash = new HashSet<String>();
		for (int i = 0; i < samples.length; i++) {
			sampleHash.add(samples[i]);
		}
		log.reportTimeInfo("Removing " + toRemove.length + " variants from " + vcf);
		String extractVCF = extractIDs(vcf, idFile, finalDir, true, true, sampleHash, true, false, log);
		convertToPlinkSet(finalDir, extractVCF, "plink", PLINK_SET_MODE.GWAS_QC, log);

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
	public static String extractIDs(String vcf, String idFile, String outputDir, boolean skipFiltered, boolean gzipOutput, HashSet<String> samples, boolean locusID, boolean keepIDs, Logger log) {
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
				copyHeader(reader, writer, BLANK_SAMPLE, HEADER_COPY_TYPE.FULL_COPY, log);
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
						String anno = locusID ? new VCOps.LocusID(vc).getId() : VCOps.getAnnotationsFor(new String[] { "snp138" }, vc, ".")[0];
						if ((!skipFiltered || !vc.isFiltered()) && (keepIDs && tmp.contains(anno)) || (!keepIDs && !tmp.contains(anno))) {
							VariantContextBuilder builder = new VariantContextBuilder(vc);
							if (vc.getID().equals(".")) {
								builder.id(anno);
							} else {
								builder.id(new VCOps.LocusID(vc).getId());
							}
							writer.add(samples == null ? builder.make() : VCOps.getSubset(builder.make(), samples, VC_SUBSET_TYPE.SUBSET_STRICT, false));
							found++;
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

	public static String extractSegments(String vcf, String segmentFile, int bpBuffer, String bams, String outputDir, boolean skipFiltered, boolean gzipOutput, int numThreads, Logger log) {
		BamExtractor.BamSample bamSample = null;

		if (vcf == null || !Files.exists(vcf)) {
			log.reportFileNotFound(vcf);
			return null;
		}
		if (segmentFile == null || !Files.exists(segmentFile)) {
			log.reportFileNotFound(segmentFile);
			return null;
		}

		String dir = outputDir == null ? ext.parseDirectoryOfFile(vcf) : outputDir;
		new File(dir).mkdirs();
		String root = getAppropriateRoot(vcf, true);

		VCFFileReader reader = new VCFFileReader(vcf, true);
		if (bams == null) {
			log.reportTimeInfo("A bam directory was not provided, skipping bam extraction");
		} else {
			log.reportTimeInfo("A bam directory was provided, extracting bams to " + dir);
			if (Files.isDirectory(bams)) {
				bamSample = new BamExtractor.BamSample(Files.listFullPaths(bams, ".bam", false), log, true);
			} else {
				bamSample = new BamExtractor.BamSample(HashVec.loadFileToStringArray(bams, false, new int[] { 0 }, false), log, true);
			}
			bamSample.generateMap();
			bamSample.getBamSampleMap();
			bamSample.verify(getSamplesInFile(reader));
		}
		String output = dir + root + "." + ext.rootOf(segmentFile) + ".vcf" + (gzipOutput ? ".gz" : "");
		if (!Files.exists(output)) {
			Segment[] segsToSearch = null;
			if (segmentFile.endsWith(".in") || segmentFile.endsWith(".bim")) {
				segsToSearch = Segment.loadRegions(segmentFile, 0, 3, 3, false);
			} else {
				segsToSearch = Segment.loadRegions(segmentFile, 0, 1, 2, 0, true, true, true, bpBuffer);
			}
			log.reportTimeInfo("Loaded " + segsToSearch.length + " segments to search");
			VariantContextWriter writer = initWriter(output, DEFUALT_WRITER_OPTIONS, getSequenceDictionary(reader));
			copyHeader(reader, writer, BLANK_SAMPLE, HEADER_COPY_TYPE.FULL_COPY, log);
			int progress = 0;
			int found = 0;

			for (VariantContext vc : reader) {
				progress++;
				if (progress % 100000 == 0) {
					log.reportTimeInfo(progress + " variants read...");
					log.reportTimeInfo(found + " variants found...");

				}
				if ((!skipFiltered || !vc.isFiltered()) && VCOps.isInTheseSegments(vc, segsToSearch)) {
					writer.add(vc);
					if (bamSample != null) {
						bamSample.addSegmentToExtract(new Segment(Positions.chromosomeNumber(vc.getChr()), vc.getStart(), vc.getEnd()));
					}
					found++;
				}
			}
			if (bamSample != null) {
				BamExtractor.extractAll(bamSample, dir, bpBuffer, true, true, numThreads, log);
				bamSample = new BamExtractor.BamSample(Files.listFullPaths(dir, ".bam", false), log, true);
				bamSample.generateMap();
				bamSample.dumpToIGVMap(output);
			}
			writer.close();
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
		String[] toSplit = getAllContigs(vcf, log);
		log.reportTimeInfo("Detected " + toSplit.length + " chrs to split");
		log.reportTimeInfo(Array.toStr(toSplit, "\n"));
		VCFSplitProducer producer = new VCFSplitProducer(vcf, toSplit, log);
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

		private VCFSplitProducer(String vcfFile, String[] toSplit, Logger log) {
			super();
			this.vcfFile = vcfFile;
			this.toSplit = toSplit;
			this.log = log;
			this.index = 0;
		}

		@Override
		public boolean hasNext() {
			return index < toSplit.length;
		}

		@Override
		public Callable<ChrSplitResults> next() {
			final String chr = toSplit[index];
			String dir = ext.parseDirectoryOfFile(vcfFile) + chr + "/";
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
			vc.getChr();
			count++;
			if (count > 0) {
				reader.close();
				return true;
			}
		}
		reader.close();
		return false;
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
				if (vc.getChr().equals(chr)) {
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
					created = false;
					log.reportTimeWarning("Indexing not quite implemented yet for " + VCF_EXTENSIONS.GZIP_VCF.getLiteral());

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
	public static boolean gzipAndIndex(String vcfFile, Logger log) {
		boolean created = false;
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
				created = true;

			}
		}
		return created;
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
		boolean skipFiltered = false;
		boolean standardFilters = false;
		boolean gzip = true;
		Logger log;

		String usage = "\n" + "seq.analysis.VCFUtils requires 0-1 arguments\n";
		usage += "   (1) full path to a vcf file (i.e. vcf=" + vcf + " (default))\n" + "";
		usage += "   (2) utility type (i.e. utility=" + type + " (default))\n" + "";
		usage += "   (3) full path to a file (can be comma delimited for homogeneity utility) defining a population for the vcf (i.e. pop= (no default))\n" + "";
		usage += "   (4) the type of vcf extension (i.e. pop= (no default))\n" + "";
		usage += "   (5) full path to a file name with chr,start,stop or *.bim to extract (i.e. segs= (no default))\n" + "";
		usage += "   (6) bp buffer for segments to extract (i.e. bp=" + bpBuffer + "(default))\n" + "";
		usage += "   (7) a bam directory to extract associtated reads (i.e. bams=" + bams + "( no default))\n" + "";
		usage += "   (8) an output directory for extracted vcfs/minibams (i.e. outDir=" + outDir + "( no default))\n" + "";
		usage += "   (9) skip filtered variants when extracting (i.e. -skipFiltered (not the default))\n" + "";
		usage += "   (10) gzip the output when extracting (i.e. -gzip ( the default))\n" + "";
		usage += "   (11) full path to a file of ids (i.e. idFile= (no default))\n" + "";
		usage += "   (12) when removing filtered variants, apply our standard filters as well (i.e. -standardFilters (not the default, GQ >=" + VARIANT_FILTER_DOUBLE.GQ_LOOSE.getDFilter() + " and DP >=" + VARIANT_FILTER_DOUBLE.DP.getDFilter() + "))\n" + "";

		usage += "   NOTE: available utilities are:\n";

		for (int i = 0; i < UTILITY_TYPE.values().length; i++) {
			usage += UTILITY_TYPE.values()[i] + "\n";
		}
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("vcf=")) {
				vcf = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("utility=")) {
				type = UTILITY_TYPE.valueOf(ext.parseStringArg(args[i], ""));
				numArgs--;
			} else if (args[i].startsWith("vpopFile=")) {
				populationFile = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("segs=")) {
				segmentFile = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("idFile=")) {
				idFile = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("bams=")) {
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
			} else if (args[i].startsWith("-skipFiltered")) {
				skipFiltered = true;
				numArgs--;
			} else if (args[i].startsWith("-standardFilters")) {
				standardFilters = true;
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
				VcfPopulation.splitVcfByPopulation(vcf, populationFile, log);
				break;
			case EXTRACT_SEGMENTS:
				extractSegments(vcf, segmentFile, bpBuffer, bams, outDir, skipFiltered, gzip, 1, log);
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
				extractIDs(vcf, idFile, outDir, skipFiltered, gzip, null, false, true, log);
				break;
			case DUMP_SAMPLES:
				extractSamps(vcf, log);
				break;
			case HOMOGENEITY:
				runHomoGeneity(vcf, populationFile.split(","), log);
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
