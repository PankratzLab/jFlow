package org.genvisis.seq.cnv;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.concurrent.Callable;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.manage.Markers;
import org.genvisis.cnv.manage.TransposeData;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.WorkerTrain.AbstractProducer;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.cnv.ExomeDepth.ExomeDepthAnalysis;
import org.genvisis.seq.manage.BamOps;
import org.genvisis.seq.manage.VCFOps.VcfPopulation;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import org.genvisis.stats.Maths;

/**
 * @author lane0212 Handles sample and global exclusions for the reference set, need in particular
 *         when related samples are run together (i.e trios)
 */
public class ExomeDepthRun {

	public static void runExomeDepth(	String bams, String vpopFile, String outputDir,
																		String outputRoot, String rLoc, boolean somaticMode,
																		int numthreads, Logger log) {
		VcfPopulation vpop = null;
		String[] allReferenceBamFiles = Files.isDirectory(bams)
																															? Files.listFullPaths(bams,
																																									BamOps.BAM_EXT,
																																									false)
																														: HashVec.loadFileToStringArray(bams,
																																														false,
																																														new int[] {0},
																																														true);
		String outputResultsDir = (outputDir == null ? ext.parseDirectoryOfFile(bams) : outputDir)
															+ "results/";
		new File(outputResultsDir).mkdirs();
		ExomeDepth exomeDepth = new ExomeDepth(	allReferenceBamFiles, allReferenceBamFiles,
																						outputResultsDir, outputRoot, rLoc, log);
		if (!Files.exists(exomeDepth.getCountFile())) {
			log.reportTimeWarning("Did not find "	+ exomeDepth.getCountFile()
														+ ", generating it now (takes a long time)");
			exomeDepth.generateCountFile();
		} else {
			log.reportTimeWarning("Using existing count file " + exomeDepth.getCountFile());
		}
		if (vpopFile == null) {
			log.reportTimeWarning("A vpopulation file was not provided, sample specific and global exclusions will not be applied");
		} else {
			if (somaticMode) {
				vpop = VcfPopulation.load(vpopFile, POPULATION_TYPE.TUMOR_NORMAL, log);
			} else {
				vpop = VcfPopulation.load(vpopFile, POPULATION_TYPE.EXOME_DEPTH, log);
			}
			vpop.report();
			exomeDepth.parseVpop(vpop);
		}
		ExomeDepthAnalysis[] eDepthAnalysis = ExomeDepth.callCNVs(exomeDepth, outputResultsDir,
																															outputRoot, numthreads, log);
		log.reportTimeInfo("Finished running exome depth for " + eDepthAnalysis.length + " .bam files");

		log.reportTimeInfo("Generating project in " + outputResultsDir);
		String projectDir = (outputDir == null ? ext.parseDirectoryOfFile(bams) : outputDir)
												+ "project_" + outputRoot + "/";
		new File(projectDir).mkdirs();
		String projectFile = projectDir + outputRoot + ".properties";
		Files.write("", projectFile);
		Project proj = new Project(projectFile, false);
		proj.PROJECT_DIRECTORY.setValue(projectDir);
		proj.saveProperties(projectFile);
		generateMarkerPositions(proj, eDepthAnalysis[0], log);
		String currentCnvFile = outputResultsDir + outputRoot + ".all.cnvs";
		proj.CNV_FILENAMES.setValue(new String[] {currentCnvFile});
		GenvisisSampleProducer producer = new GenvisisSampleProducer(	proj, eDepthAnalysis, outputRoot,
																																	log);
		WorkerTrain<ExomeSample> train = new WorkerTrain<ExomeDepthRun.ExomeSample>(producer,
																																								numthreads, 10,
																																								log);
		Hashtable<String, Float> allOutliers = new Hashtable<String, Float>();
		ExomeSample[] exomeSamples = new ExomeSample[eDepthAnalysis.length];
		int index = 0;
		while (train.hasNext()) {
			ExomeSample eSample = train.next();
			exomeSamples[index] = eSample;
			allOutliers.putAll(eSample.getOutliers());
			index++;
			if (index % 100 == 0) {
				log.reportTimeInfo("Imported " + index + " samples");
			}
		}
		String outliersSer = proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser";
		SerializedFiles.writeSerial(allOutliers, outliersSer);
		TransposeData.transposeData(proj, 2000000000, false);
		proj.saveProperties();
		SampleData.createMinimalSampleData(proj);

	}

	private static void generateMarkerPositions(Project proj, ExomeDepthAnalysis first, Logger log) {
		String[] header = Files.getHeaderOfFile(first.getExomeDepthRawDataOutput(), log);
		int chrCol = ext.indexOfStr("anno.chromosome", header);
		int startCol = ext.indexOfStr("anno.start", header);
		int stopCol = ext.indexOfStr("anno.end", header);

		Segment[] segments = Segment.loadRegions(	first.getExomeDepthRawDataOutput(), chrCol, startCol,
																							stopCol, true);
		String positions = proj.MARKER_POSITION_FILENAME.getValue();
		proj.getLog().reportTimeInfo("Postions will be set to the midpoint of each segment");
		String[] markerNames = new String[segments.length];
		if (!Files.exists(positions) || !Files.exists(proj.MARKERSET_FILENAME.getValue(true, true))) {
			try {
				PrintWriter writer = new PrintWriter(new FileWriter(positions));
				int markerIndex = 0;
				writer.println("BinName\tChr\tPosition");
				for (Segment bFeatureSeg : segments) {
					String markerName = bFeatureSeg.getUCSClocation();
					markerNames[markerIndex] = markerName;
					int diff = bFeatureSeg.getStop() - bFeatureSeg.getStart();
					int mid = Math.round((float) diff / 2);
					int pos = bFeatureSeg.getStart() + mid;
					writer.println(markerName + "\t" + bFeatureSeg.getChr() + "\t" + pos);
					markerIndex++;
				}

				writer.close();
			} catch (Exception e) {
				proj.getLog().reportError("Error writing to " + positions);
				proj.getLog().reportException(e);
			}

			Markers.orderMarkers(	markerNames, proj.MARKER_POSITION_FILENAME.getValue(),
														proj.MARKERSET_FILENAME.getValue(true, true), proj.getLog());
		}
	}

	private static class ExomeSample implements Callable<ExomeSample> {
		private final Project proj;
		private final ExomeDepthAnalysis eDepthAnalysis;
		// private String outputSampleRoot;
		private Hashtable<String, Float> outliers;
		private final Logger log;

		public ExomeSample(	Project proj, ExomeDepthAnalysis eDepthAnalysis, String outputSampleRoot,
												Logger log) {
			super();
			this.proj = proj;
			this.eDepthAnalysis = eDepthAnalysis;
			// this.outputSampleRoot = outputSampleRoot;
			this.log = log;
			outliers = new Hashtable<String, Float>();
		}

		public Hashtable<String, Float> getOutliers() {
			return outliers;
		}

		@Override
		public ExomeSample call() throws Exception {
			String sample = BamOps.getSampleName(eDepthAnalysis.getInputBam());
			String sampFile = proj.SAMPLE_DIRECTORY.getValue(true, false)	+ sample
												+ Sample.SAMPLE_FILE_EXTENSION;

			if (!Files.exists(sampFile)) {
				String input = eDepthAnalysis.getExomeDepthRawDataOutput();
				int ratioIndex = ext.indexOfStr("ratio", Files.getHeaderOfFile(input, log));
				String[] ratios = HashVec.loadFileToStringArray(input, true, new int[] {ratioIndex}, false);
				float[] ratioLrr = new float[ratios.length];
				for (int i = 0; i < ratios.length; i++) {
					try {
						float tmp = Float.parseFloat(ratios[i]);
						ratioLrr[i] = (float) Maths.log2(tmp);
					} catch (NumberFormatException nfe) {
						ratioLrr[i] = Float.NaN;
					}
				}
				byte[] genos = ArrayUtils.byteArray(ratioLrr.length, (byte) 1);
				float[] zeroArray = ArrayUtils.floatArray(ratioLrr.length, 0);
				Sample samp = new Sample(	sample, proj.getMarkerSet().getFingerprint(), zeroArray, zeroArray,
																	zeroArray, zeroArray, ratioLrr, genos, genos, false);
				samp.saveToRandomAccessFile(sampFile, outliers, sample);
			} else {
				outliers = Sample.loadOutOfRangeValuesFromRandomAccessFile(sampFile);
			}
			return this;
		}
	}

	private static class GenvisisSampleProducer extends AbstractProducer<ExomeSample> {
		private final Project proj;
		private final String outputSampleRoot;
		private final ExomeDepthAnalysis[] eDepthAnalysis;
		private final Logger log;
		private int index;

		public GenvisisSampleProducer(Project proj, ExomeDepthAnalysis[] eDepthAnalysis,
																	String outputSampleRoot, Logger log) {
			super();
			this.proj = proj;
			this.eDepthAnalysis = eDepthAnalysis;
			this.outputSampleRoot = outputSampleRoot;
			this.log = log;
			index = 0;
		}

		@Override
		public boolean hasNext() {
			return index < eDepthAnalysis.length;
		}

		@Override
		public Callable<ExomeSample> next() {
			ExomeSample exomeSample = new ExomeSample(proj, eDepthAnalysis[index], outputSampleRoot, log);
			index++;
			return exomeSample;
		}
	}

	// private Sample parseToSample(ExomeDepthAnalysis exomeDepthAnalysis, Logger log) {
	// String input = exomeDepthAnalysis.getExomeDepthRawDataOutput();
	// int ratioIndex = ext.indexOfStr("ratio", Files.getHeaderOfFile(input, log));
	// String[] ratios = HashVec.loadFileToStringArray(input, true, new int[] { ratioIndex }, false);
	// float[] ratioLrr = new float[ratios.length];
	// for (int i = 0; i < ratios.length; i++) {
	// try {
	// float tmp = Float.parseFloat(ratios[i]);
	// ratioLrr[i] = tmp;
	// } catch (NumberFormatException nfe) {
	// ratioLrr[i] = Float.NaN;
	// }
	// }
	//
	// return null;
	// }

	public static void main(String[] args) {
		int numArgs = args.length;
		String bams = "bams/";
		String outputDir = null;
		String outputRoot = "ExomeDepth";
		int numthreads = 1;
		String vpopFile = null;
		String logfile = null;
		String Rloc = null;
		boolean somaticMode = false;
		Logger log;

		String usage = "\n" + "seq.analysis.ExomeDepth requires 0-1 arguments\n";
		usage += "   (1) full path to a directory of or file of bams (i.e. bams="	+ bams
							+ " (default))\n" + "";
		usage += PSF.Ext.getOutputDirCommand(2, outputDir);
		usage += "   (3) output root command (i.e. root=" + outputRoot + " (default))\n" + "";
		usage += PSF.Ext.getNumThreadsCommand(4, numthreads);
		usage +=
					"   (5) full path to a v population file, individuals with the same population will not be used as ref(i.e. vpop= (no default))\n"
							+ "";
		usage += "   (6) alternative R location (i.e. rDir= (no default))\n" + "";
		usage += "   (7) somatic mode (i.e. somaticMode=" + somaticMode + " (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("bams=")) {
				bams = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("vpop=")) {
				vpopFile = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
				numthreads = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith(PSF.Ext.OUTPUT_DIR_COMMAND)) {
				outputDir = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("rDir=")) {
				Rloc = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("root=")) {
				outputRoot = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("somaticMode=")) {
				somaticMode = ext.parseBooleanArg(arg);
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
			runExomeDepth(bams, vpopFile, outputDir, outputRoot, Rloc, somaticMode, numthreads, log);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
