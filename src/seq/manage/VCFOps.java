package seq.manage;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Set;
import java.util.TreeSet;

import filesys.Segment;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import common.Array;
import common.CmdLine;
import common.Files;
import common.Logger;
import common.PSF;
import common.ext;

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
		 * Run a full gwas qc on a vcf file
		 */
		GWAS_QC, /**
		 * Conver a vcf to plink
		 */
		CONVERT_PLINK,

		/**
		 * Subset a vcf by a super population
		 */
		SUBSET_SUPER, /**
		 * 
		 */
		EXTRACT_SEGMENTS;
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
	 * @param vcf
	 *            a vcf file to convert to plink
	 * @param rootOut
	 *            the root output for the plink files
	 * @param log
	 */
	public static void convertToPlinkSet(String vcf, String rootOut, Logger log) {
		String[] plinkCommand = null;
		String dir = ext.parseDirectoryOfFile(vcf) + "plink" + ext.rootOf(vcf) + "/";
		new File(dir).mkdirs();
		rootOut = dir + rootOut;
		String[] outFiles = PSF.Plink.getPlinkBedBimFam(rootOut);
		if (!Files.exists(rootOut, outFiles)) {
			plinkCommand = PSF.Plink.getPlinkVCFCommand(vcf, rootOut);
			CmdLine.runCommandWithFileChecks(plinkCommand, "", new String[] { vcf }, outFiles, true, true, false, log);
		} else {
			log.reportTimeWarning("Detected that the following files already exist " + Array.toStr(outFiles));
		}
		if (!Files.exists(dir + ".qc.pbs")) {
			String gwasQC = Array.toStr(PSF.Load.getAllModules(), "\n") + "\njcp gwas.Qc dir=" + dir;
			Files.qsub(dir + "qc.pbs", gwasQC, 62000, 24, 16);
		}
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
				convertToPlinkSet(vcf, dir + "plink", log);
			}
			log.reportTimeInfo("Running gwas.qc on the following files in " + dir + ":");
			log.reportTimeInfo("\t" + Array.toStr(plinkFiles, "\n"));
			gwas.Qc.fullGamut(dir, false);
		} else {
			log.reportFileNotFound(vcf);
		}
	}

	/**
	 * Class that manages the population structure represented in a vcf<br>
	 * Can be used for HWE tests on sub and super populations etc...
	 *
	 */
	public static class VcfPopulation {
		private static final String[] HEADER = new String[] { "IID", "Population", "SuperPopulation" };
		private static final String SKIP = "#N/A";
		private Hashtable<String, Set<String>> subPop;
		private Hashtable<String, Set<String>> superPop;
		private ArrayList<String> uniqSubPop;
		private ArrayList<String> uniqSuperPop;

		public VcfPopulation() {
			this.subPop = new Hashtable<String, Set<String>>();
			this.superPop = new Hashtable<String, Set<String>>();
			this.uniqSubPop = new ArrayList<String>();
			this.uniqSuperPop = new ArrayList<String>();
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
		public void report(Logger log) {
			for (String key : subPop.keySet()) {
				log.reportTimeInfo("Sub - population " + key + " had " + subPop.get(key).size() + " individuals");
			}
			for (String key : superPop.keySet()) {
				log.reportTimeInfo("Super - population " + key + " had " + superPop.get(key).size() + " individuals");
			}
		}

		private VariantContextWriter[] getWritersForPop(String outputbase, VCFFileReader reader, Logger log) {
			VariantContextWriter[] writers = new VariantContextWriter[uniqSuperPop.size()];
			for (int i = 0; i < uniqSuperPop.size(); i++) {
				String sPop = uniqSuperPop.get(i);
				String out = outputbase + "." + uniqSuperPop.get(i) + ".vcf.gz";
				writers[i] = initWriter(out, DEFUALT_WRITER_OPTIONS, getSequenceDictionary(reader));
				copyHeader(reader, writers[i], superPop.get(sPop), HEADER_COPY_TYPE.SUBSET_STRICT, log);
			}
			return writers;
		}

		public static void splitVcfByPopulation(String vcf, String fullPathToPopFile, Logger log) {
			if (vcf != null && !Files.exists(vcf)) {
				log.reportFileNotFound(vcf);
				return;
			}
			if (fullPathToPopFile != null && !Files.exists(fullPathToPopFile)) {
				log.reportFileNotFound(fullPathToPopFile);
				return;
			}
			VcfPopulation vpop = VcfPopulation.load(fullPathToPopFile, log);
			vpop.report(log);
			VCFFileReader reader = new VCFFileReader(vcf, true);

			String dir = ext.parseDirectoryOfFile(vcf);
			String root = ext.rootOf(vcf).replaceFirst(VCF_EXTENSIONS.REG_VCF.getLiteral(), "");

			VariantContextWriter[] writers = vpop.getWritersForPop(dir + root, reader, log);
			int progress = 0;
			for (VariantContext vc : reader) {
				progress++;
				if (progress % 100000 == 0) {
					log.reportTimeInfo(progress + " variants read...");
				}
				for (int i = 0; i < writers.length; i++) {
					writers[i].add(VCOps.getSubset(vc, vpop.getSuperPop().get(vpop.getUniqSuperPop().get(i))));
				}
			}
			for (int i = 0; i < writers.length; i++) {
				writers[i].close();
			}
		}

		/**
		 * @param fullPathToPopFile
		 *            load this file to an {@link VcfPopulation}
		 * @param log
		 * @return
		 */
		public static VcfPopulation load(String fullPathToPopFile, Logger log) {
			VcfPopulation vcfPopulation = new VcfPopulation();
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
			return vcfPopulation;
		}

	}

	public static void extractSegments(String vcf, String segmentFile, Logger log) {
		if (vcf != null && !Files.exists(vcf)) {
			log.reportFileNotFound(vcf);
			return;
		}
		if (segmentFile != null && !Files.exists(segmentFile)) {
			log.reportFileNotFound(segmentFile);
			return;
		}
		Segment[] segsToSearch = null;
		if (segmentFile.endsWith(".in") || segmentFile.endsWith(".bim")) {
			segsToSearch = Segment.loadRegions(segmentFile, 0, 3, 3, false);
		} else {
			segsToSearch = Segment.loadRegions(segmentFile, 0, 1, 2, false);
		}
		String dir = ext.parseDirectoryOfFile(vcf);
		String root = ext.rootOf(vcf).replaceFirst(VCF_EXTENSIONS.REG_VCF.getLiteral(), "");

		VCFFileReader reader = new VCFFileReader(vcf, true);
		String output = dir + root + "." + ext.rootOf(segmentFile) + ".vcf.gz";
		VariantContextWriter writer = initWriter(output, DEFUALT_WRITER_OPTIONS, getSequenceDictionary(reader));
		copyHeader(reader, writer, BLANK_SAMPLE, HEADER_COPY_TYPE.FULL_COPY, log);
		int progress =0;
		int found =0;

		for (VariantContext vc : reader) {
			progress++;
			if (progress % 100000 == 0) {
				log.reportTimeInfo(progress + " variants read...");
				log.reportTimeInfo(found + " variants found...");

			}
			if (VCOps.isInTheseSegments(vc, segsToSearch)) {
				writer.add(vc);
				found++;
			}
		}
		reader.close();
		writer.close();
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String vcf = "Avcf.vcf";
		String populationFile = null;
		String logfile = null;
		UTILITY_TYPE type = UTILITY_TYPE.GWAS_QC;
		String segmentFile = null;
		Logger log;

		String usage = "\n" + "seq.analysis.VCFUtils requires 0-1 arguments\n";
		usage += "   (1) full path to a vcf file (i.e. vcf=" + vcf + " (default))\n" + "";
		usage += "   (2) utility type (i.e. utility=" + type + " (default))\n" + "";
		usage += "   (3) full path to a file defining a population for the vcf (i.e. pop= (no default))\n" + "";
		usage += "   (4) the type of vcf extension (i.e. pop= (no default))\n" + "";
		usage += "   (5) full path to a file name with chr,start,stop or *.bim to extract (i.e. segs= (no default))\n" + "";

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
			log.reportTimeInfo("Running utiltity type: " + type);
			switch (type) {
			case GWAS_QC:
				vcfGwasQC(vcf, log);
				break;
			case CONVERT_PLINK:
				convertToPlinkSet(vcf, "plink", log);
				break;
			case SUBSET_SUPER:
				VcfPopulation.splitVcfByPopulation(vcf, populationFile, log);
				break;
			case EXTRACT_SEGMENTS:
				extractSegments(vcf, segmentFile, log);
				break;
			default:
				break;
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
