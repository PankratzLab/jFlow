package seq.manage;

import java.io.File;

import common.Files;
import common.Logger;
import common.PSF;
import common.ext;
import seq.analysis.VcfQuery.Location;
import seq.manage.VCFOps.HEADER_COPY_TYPE;
import seq.manage.VCFOps.VcfPopulation;
import seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import seq.qc.FilterNGS.VARIANT_FILTER_BOOLEAN;
import seq.qc.FilterNGS.VARIANT_FILTER_DOUBLE;
import seq.qc.FilterNGS.VariantContextFilter;
import seq.qc.FilterNGS.VariantContextFilterPass;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;

/**
 * https://www.cog-genomics.org/plink2/dev#pigz
 *
 */
public class PopGen {
	public static final String DEFAULT_FTP_DIR = "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/";
	public static final String DEFAULT_FTP_EXT = ".genotypes.vcf.gz";

	public static void filterForPopGen(String directory, String outputVCF, String vcfSuffix, Location location, final VCFOps.VcfPopulation hwePopTests, int numThreads, Logger log) {
		String[] vcfs = null;
		switch (location) {
		case LOCAL:
			vcfs = Files.listFullPaths(directory, vcfSuffix, false);
			break;
		case REMOTE:
			vcfs = Files.parseRemoteFTPFiles(directory, vcfSuffix, log);
			break;
		default:
			log.reportTimeError("Invalid file location " + location);
			break;
		}
		VARIANT_FILTER_DOUBLE[] fullPopFilters = hwePopTests == null ? VARIANT_FILTER_DOUBLE.values() : VARIANT_FILTER_DOUBLE.getFiltersExcluding(new VARIANT_FILTER_DOUBLE[] { VARIANT_FILTER_DOUBLE.HWE });
		final VariantContextFilter vcfFilter = new VariantContextFilter(fullPopFilters, VARIANT_FILTER_BOOLEAN.values(), null, null, log);
		final VariantContextFilter superPopFilter = new VariantContextFilter(new VARIANT_FILTER_DOUBLE[] { VARIANT_FILTER_DOUBLE.HWE, VARIANT_FILTER_DOUBLE.CALL_RATE_LOOSE }, new VARIANT_FILTER_BOOLEAN[] {}, null, null, log);

		VCFFileReader tmp = new VCFFileReader(vcfs[0], true);
		VariantContextWriter writer = VCFOps.initWriter(outputVCF, VCFOps.DEFUALT_WRITER_OPTIONS, VCFOps.getSequenceDictionary(tmp));
		VCFOps.copyHeader(tmp, writer, VCFOps.BLANK_SAMPLE, HEADER_COPY_TYPE.FULL_COPY, log);

		tmp.close();
		int pass = 0;
		int fail = 0;
		int hweSuperPopFail = 0;
		for (int i = 0; i < vcfs.length; i++) {
			log.reportTimeInfo("Initializing reader for " + vcfs[i]);
			tmp = new VCFFileReader(vcfs[i], true);
			log.reportTimeInfo("Finished initializing reader for " + vcfs[i]);
			long time = System.currentTimeMillis();
			for (final VariantContext variantContext : tmp) {
				VariantContextFilterPass vcfp = vcfFilter.filter(variantContext);
				VariantContextFilterPass superPopPass = new VariantContextFilterPass(true, "Did not need to test super populations");
				if (vcfp.passed() && hwePopTests != null) {
					superPopPass = testSuperPopulations(hwePopTests, superPopFilter, variantContext, log);
				}
				if (vcfp.passed() && superPopPass.passed() && !variantContext.isFiltered()) {
					try {
						writer.add(variantContext);
					} catch (IllegalStateException is) {
						log.reportTimeError("A variant violated vcf convention, skipping file " + vcfs[i]);
						tmp.close();
						break;
					} catch (IllegalArgumentException ia) {
						log.reportTimeError("A variant violated vcf convention, skipping file " + vcfs[i]);
						tmp.close();
						break;
					}
					pass++;
				} else {
					fail++;
					if (!superPopPass.passed()) {
						hweSuperPopFail++;
					}
				}
				if ((pass + fail) % 10000 == 0) {
					log.report("");
					log.reportTimeInfo(ext.getTimeElapsed(time) + " since last update");
					log.reportTimeInfo("Currently reading file (" + (i + 1) + " of " + vcfs.length + ") : " + vcfs[i]);
					log.reportTimeInfo(pass + " variants passed");
					log.reportTimeInfo(fail + " variants failed");
					log.reportTimeInfo(hweSuperPopFail + " variants failed for super population HWE tests");
					time = System.currentTimeMillis();
				}
			}
			tmp.close();
		}
		writer.close();
	}

	private static VariantContextFilterPass testSuperPopulations(final VCFOps.VcfPopulation hwePopTests, final VariantContextFilter hweFilter, final VariantContext variantContext, Logger log) {
		VariantContextFilterPass hweContextFilterPass = new VariantContextFilterPass(true, "Passed");
		if (hwePopTests != null) {
			for (String superPopulation : hwePopTests.getSuperPop().keySet()) {
				VariantContext vcSub = VCOps.getSubset(variantContext, hwePopTests.getSuperPop().get(superPopulation));
				VariantContextFilterPass vcfpHWE = hweFilter.filter(vcSub);
				if (!vcfpHWE.passed()) {
					hweContextFilterPass = vcfpHWE;
					break;
				}
			}
		} else {
			hweContextFilterPass = new VariantContextFilterPass(true, "Super population not provided");
		}
		return hweContextFilterPass;
	}

	public static void popGen(String inputDirectory, String output, String populationFile, Location location, String extension, int numThreads, Logger log) {
		VcfPopulation vpop = VCFOps.VcfPopulation.load(populationFile,POPULATION_TYPE.ANY, log);
		if (vpop != null) {
			vpop.report();
		}
		filterForPopGen(inputDirectory, output, extension, location, vpop, numThreads, log);
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String inputDirectory = "/Pop/";
		String logfile = null;
		Logger log;
		String output = "Pop.vcf.gz";
		String populationFile = "vpop.txt";
		String extension = ".vcf.gz";
		Location location = Location.LOCAL;
		int numThreads = 3;
		String usage = "\n" + "seq.manage.PopGen requires 0-1 arguments\n";
		usage += "   (1) input directory containing vcf files (i.e. inputDirectory=" + inputDirectory + " (default))\n" + "";
		usage += "   (2) full path to an output file (i.e. output=" + output + " (default))\n" + "";
		usage += "   (3) full path to file defining population sub and super structure (i.e. vpopFile=" + populationFile + " (default))\n" + "";
		usage += "   (4) the directory is a remote directory (i.e. -remote (default))\n" + "";
		usage += "   (5) extension of vcf files to use(i.e. ext=" + extension + " (default))\n" + "";
		usage += PSF.Ext.getNumThreadsCommand(6, numThreads);
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("inputDirectory=")) {
				inputDirectory = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("output=")) {
				output = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("vpopFile=")) {
				populationFile = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("ext=")) {
				extension = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
				numThreads = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("-remote")) {
				location = Location.REMOTE;
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
			new File(ext.parseDirectoryOfFile(output)).mkdirs();
			log = new Logger(logfile == null ? ext.rootOf(output, false) + ".log" : logfile);
			popGen(inputDirectory, output, populationFile, location, extension, numThreads, log);
			// parse(filename, log);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}

// private static int ReportPassingVariants(Logger log, String vcfs, VCFFileReader tmp, VariantContextWriter writer, int pass, CallableBufferSet<VariantContext, VariantContext[]> bufferSet) {
// ArrayList<VariantContext[]> passingVariants = bufferSet.execute();
// for (int i = 0; i < passingVariants.size(); i++) {
// for (int j = 0; j < passingVariants.get(i).length; j++) {
// try {
// writer.add(passingVariants.get(i)[j]);
// } catch (IllegalStateException is) {
// log.reportTimeError("A variant violated vcf convention, skipping file " + vcfs);
// tmp.close();
//
// return -1;
// // log.reportException(is);
// } catch (IllegalArgumentException ia) {
// log.reportTimeError("A variant violated vcf convention, skipping file " + vcfs);
// tmp.close();
//
// return -1;
//
// // log.reportException(ia);
// }
// pass++;
// }
// }
// return pass;
// }
// private static PopGenBuffer[] initBuffers(int numThreads, int bufferSize, VariantContextFilter vcfFilter, VariantContextFilter hweFilter, VCFOps.VcfPopulation hwePopTests, VCFHeader header, Logger log) {
// PopGenBuffer[] popGenBuffers = new PopGenBuffer[numThreads];
// for (int i = 0; i < numThreads; i++) {
// popGenBuffers[i] = new PopGenBuffer(bufferSize, vcfFilter, hweFilter, hwePopTests, new ArrayList<VariantContext>(), header, log);
// }
// return popGenBuffers;
// }
//
// private static class PopGenBuffer extends CallableBufferSet.CallableBuffer<VariantContext, VariantContext[]> {
// private VariantContextFilter vcfFilter;
// private VariantContextFilter hweFilter;
// private VCFOps.VcfPopulation hwePopTests;
// private ArrayList<VariantContext> otherBuffer;
// private VCFHeader header;
// private Logger log;
//
// public PopGenBuffer(int bufferSize, VariantContextFilter vcfFilter, VariantContextFilter hweFilter, VCFOps.VcfPopulation hwePopTests, ArrayList<VariantContext> vcBuffer, VCFHeader header, Logger log) {
// super(bufferSize, vcBuffer);
// this.otherBuffer = vcBuffer;
// this.header = header;
// this.vcfFilter = vcfFilter;
// this.hweFilter = hweFilter;
// this.hwePopTests = hwePopTests;
// this.log = log;
// }
//
// public void addToOther(VariantContext vc) {
// otherBuffer.add(vc);
// }
//
// @Override
// public VariantContext[] call() throws Exception {
// ArrayList<VariantContext> variantsPassing = new ArrayList<VariantContext>(getBuffer().size());
// VariantContextFilterPass[] passes = new VariantContextFilterPass[getBuffer().size()];
//
// for (int i = 0; i < otherBuffer.size(); i++) {
// VariantContext vc = otherBuffer.get(i).fullyDecode(header, false);
// // System.out.println(Array.toStr(Array.toStringArray(header.getSampleNamesInOrder())));
// System.out.println(vc.isFullyDecoded());
// VariantContextFilterPass vcfp = vcfFilter.filter(vc);
// VariantContextFilterPass hwefp = new VariantContextFilterPass(true, "Did not need to test super populations");
// if (vcfp.passed() && hwePopTests != null) {
// hwefp = testSuperPopulations(hwePopTests, hweFilter, vc, log);
// }
// if (vcfp.passed() && hwefp.passed()) {
// variantsPassing.add(vc);
// System.out.println("HDSFSDF");
// }
// // System.out.println(vcfp.getTestPerformed());
// passes[i] = new VariantContextFilterPass(vcfp.passed() && hwefp.passed(), "All Sample Filter: " + vcfp.getTestPerformed() + " SuperPopFilter: " + hwefp.getTestPerformed());
// }
// // TODO Auto-generated method stub
// return variantsPassing.toArray(new VariantContext[variantsPassing.size()]);
// }
// }
// for (int i = 0; i < vcfs.length; i++) {
// // if (i >= 22) {
// log.reportTimeInfo("Initializing reader for " + vcfs[i]);
// tmp = new VCFFileReader(vcfs[i], true);
// log.reportTimeInfo("Finished initializing reader for " + vcfs[i]);
// long time = System.currentTimeMillis();
// for (final VariantContext variantContext : tmp) {
// VariantContextFilterPass vcfp = vcfFilter.filter(variantContext);
// VariantContextFilterPass hwefp = new VariantContextFilterPass(true, "Did not need to test super populations");
// if (vcfp.passed() && hwePopTests != null) {
// hwefp = testSuperPopulations(hwePopTests, hweFilter, variantContext, log);
// }
//
// if (vcfp.passed() && hwefp.passed()) {
// try {
// writer.add(variantContext);
// } catch (IllegalStateException is) {
// log.reportTimeError("A variant violated vcf convention, skipping file " + vcfs[i]);
// tmp.close();
//
// break;
// // log.reportException(is);
// } catch (IllegalArgumentException ia) {
// log.reportTimeError("A variant violated vcf convention, skipping file " + vcfs[i]);
// tmp.close();
//
// break;
// // log.reportException(ia);
// }
// pass++;
// } else {
// fail++;
// if (!hwefp.passed()) {
// hweSuperPopFail++;
// }
// // log.reportTimeInfo(vcfp.getTestPerformed());
// // log.reportTimeInfo(hwefp.getTestPerformed());
// }
// if ((pass + fail) % 50000 == 0) {
// log.report("");
// log.reportTimeInfo(ext.getTimeElapsed(time) + " since last update");
// log.reportTimeInfo("Currently reading file (" + (i + 1) + " of " + vcfs.length + ") : " + vcfs[i]);
// log.reportTimeInfo(pass + " variants passed");
// log.reportTimeInfo(fail + " variants failed");
// log.reportTimeInfo(hweSuperPopFail + " variants failed for super population HWE tests");
// time = System.currentTimeMillis();
// }
// }
// tmp.close();
// // }
// }
// writer.close();
// }

// for (int i = 0; i < vcfs.length; i++) {
//
// // if (i >= 22) {
// log.reportTimeInfo("Initializing reader for " + vcfs[i]);
// tmp = new VCFFileReader(vcfs[i], true);
//
// VCFHeader header = tmp.getFileHeader();
// PopGenBuffer[] buffers = initBuffers(numThreads, 2000, vcfFilter, hweFilter, hwePopTests, header, log);
// CallableBufferSet<VariantContext, VariantContext[]> bufferSet = new CallableBufferSet<VariantContext, VariantContext[]>(numThreads, 10, buffers, log);
//
// log.reportTimeInfo("Finished initializing reader for " + vcfs[i]);
// long time = System.currentTimeMillis();
// PopGenBuffer popGenBuffer = new PopGenBuffer(50000, vcfFilter, hweFilter, hwePopTests, new ArrayList<VariantContext>(), tmp.getFileHeader(), log);
// PopGenBuffer popGenBuffer2 = new PopGenBuffer(50000, vcfFilter, hweFilter, hwePopTests, new ArrayList<VariantContext>(), tmp.getFileHeader(), log);
// boolean swiadf = false;
// for (final VariantContext variantContext : tmp) {
// // VariantContext variantContext = variantContext2.fullyDecode(header, false);
// if (swiadf) {
// swiadf = false;
// popGenBuffer.addToOther(variantContext);
// } else {
// swiadf = true;
// popGenBuffer2.addToOther(variantContext);
// }
// totalVariants++;
// if ((totalVariants) % 50000 == 0) {
// log.report("");
// log.reportTimeInfo(ext.getTimeElapsed(time) + " since last update");
// log.reportTimeInfo("Currently reading file (" + (i + 1) + " of " + vcfs.length + ") : " + vcfs[i]);
// log.reportTimeInfo(pass + " variants passed");
// log.reportTimeInfo(totalVariants + " total Variants");
//
// time = System.currentTimeMillis();
// }
// if (popGenBuffer.isFull()) {
//
// WorkerHive<VariantContext[]> hive = new WorkerHive<VariantContext[]>(numThreads, 10, log);
// popGenBuffer = new PopGenBuffer(50000, vcfFilter, hweFilter, hwePopTests, popGenBuffer.getBuffer(), tmp.getFileHeader(), log);
// popGenBuffer2 = new PopGenBuffer(50000, vcfFilter, hweFilter, hwePopTests, popGenBuffer2.getBuffer(), tmp.getFileHeader(), log);
//
// hive.addCallable(popGenBuffer);
// hive.addCallable(popGenBuffer2);
//
// hive.execute(true);
// // ArrayList<VariantContext[]> results = hive.getResults();
//
// // bufferSet.execute();
// System.exit(1);
// // pass = ReportPassingVariants(log, vcfs[i], tmp, writer, pass, bufferSet);
// bufferSet.clearBuffers();
// }
// if (pass < 0) {
// break;
// }
// // bufferSet.addToBuffers(variantContext);
// }
// tmp.close();
// // }
// }
