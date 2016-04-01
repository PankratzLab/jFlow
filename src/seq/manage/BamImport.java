package seq.manage;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.concurrent.Callable;

import seq.manage.BEDFileReader.BEDFeatureSeg;
import seq.manage.BamOps.BamIndexStats;
import seq.manage.BamSegPileUp.BamPileResult;
import seq.manage.BamSegPileUp.PileupProducer;
import seq.qc.FilterNGS;
import common.Array;
import common.Files;
import common.HashVec;
import common.Logger;
import common.Positions;
import common.WorkerTrain;
import common.WorkerTrain.Producer;
import common.ext;
import common.PSF.Ext;
import cnv.analysis.CentroidCompute;
import cnv.analysis.CentroidCompute.CentroidBuilder;
import cnv.analysis.PennCNVPrep;
import cnv.analysis.pca.PrincipalComponentsCrossTabs;
import cnv.analysis.Mosaicism;
import cnv.filesys.Centroids;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.filesys.Sample;
import cnv.filesys.Project.ARRAY;
import cnv.manage.Markers;
import cnv.manage.MitoPipeline;
import cnv.manage.TransposeData;
import cnv.qc.LrrSd;
import cnv.var.LocusSet;
import cnv.var.SampleData;
import filesys.Segment;

public class BamImport {
	// public static final String OFF_TARGET_FLAG = "OFF_TARGET";
	// public static final String VARIANT_SITE_FLAG = "VARIANT_SITE";

	public enum NGS_MARKER_TYPE {
		/**
		 * Typically representing exons
		 */
		ON_TARGET("ON_TARGET"),
		/**
		 * Off target bins
		 */
		OFF_TARGET("OFF_TARGET"),
		/**
		 * Corresponding to actual variant sites
		 */
		VARIANT_SITE("VARIANT_SITE");

		private String flag;

		private NGS_MARKER_TYPE(String flag) {
			this.flag = flag;
		}

		public String getFlag() {
			return flag;
		}

		public static NGS_MARKER_TYPE getType(String markerName) {
			NGS_MARKER_TYPE type = null;
			for (int i = 0; i < NGS_MARKER_TYPE.values().length; i++) {
				if (markerName.contains(NGS_MARKER_TYPE.values()[i].getFlag())) {
					if (type != null) {
						throw new IllegalArgumentException("Multiple types for " + markerName);

					}
					type = NGS_MARKER_TYPE.values()[i];
				}
			}

			if (type == null) {
				throw new IllegalArgumentException("Could not determine type for " + markerName);
			}
			return type;
		}

	}

	private static class BamPileConversionResults implements Callable<BamPileConversionResults> {
		private Project proj;
		private BamPileResult result;
		private String sample;
		private BamIndexStats bamIndexStats;
		private Hashtable<String, Float> outliers;
		private Logger log;
		private long fingerPrint;

		public BamPileConversionResults(Project proj, BamPileResult result, long fingerPrint, Logger log) {
			super();
			this.proj = proj;
			this.result = result;
			this.outliers = new Hashtable<String, Float>();
			this.fingerPrint = fingerPrint;
			this.log = log;
		}

		@Override
		public BamPileConversionResults call() throws Exception {
			String sampleFile = proj.SAMPLE_DIRECTORY.getValue() + BamOps.getSampleName(result.getBam()) + Sample.SAMPLE_DATA_FILE_EXTENSION;
			if (!Files.exists(sampleFile)) {
				BamSample bamSample = new BamSample(proj, result.getBam(), result.loadResults(log));
				sample = bamSample.getSampleName();
				bamIndexStats = BamOps.getBamIndexStats(result.getBam());
				outliers = bamSample.writeSample(fingerPrint);
			} else {
				log.reportFileExists(sampleFile);
				sample = BamOps.getSampleName(result.getBam());
				bamIndexStats = BamOps.getBamIndexStats(result.getBam());
				outliers = null;
				// outliers = Sample.loadOutOfRangeValuesFromRandomAccessFile(sampleFile);
			}
			return this;
		}

		public String getSample() {
			return sample;
		}

		public BamIndexStats getBamIndexStats() {
			return bamIndexStats;
		}

		public Hashtable<String, Float> getOutliers() {
			return outliers;
		}

	}

	private static class BamPileConverterProducer implements Producer<BamPileConversionResults> {
		private Project proj;
		private BamPileResult[] pileResults;
		private long fingerPrint;
		private Logger log;
		private int index;

		public BamPileConverterProducer(Project proj, BamPileResult[] pileResults, long fingerPrint, Logger log) {
			super();
			this.proj = proj;
			this.pileResults = pileResults;
			this.fingerPrint = fingerPrint;
			this.log = log;
		}

		@Override
		public boolean hasNext() {
			return index < pileResults.length;
		}

		@Override
		public Callable<BamPileConversionResults> next() {
			BamPileResult current = pileResults[index];
			BamPileConversionResults conv = new BamPileConversionResults(proj, current, fingerPrint, log);
			index++;
			return conv;
		}

		@Override
		public void shutdown() {

		}

	}

	private static class VariantSeg extends Segment {

		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;
		private String tag;

		public VariantSeg(byte chr, int start, int stop, String tag) {
			super(chr, start, stop);
			this.tag = tag;
		}

		public String getTag() {
			return tag;
		}

	}

	private static LocusSet<VariantSeg> extractVCF(Project proj, String vcf) {
		LocusSet<VariantSeg> varLocusSet = new LocusSet<VariantSeg>(new VariantSeg[] {}, true, proj.getLog()) {

			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

		};

		if (vcf != null) {

			String outDir = proj.PROJECT_DIRECTORY.getValue() + "vcf/";
			new File(outDir).mkdirs();
			String out = outDir + VCFOps.getAppropriateRoot(vcf, true) + ".regions.txt";
			String outVCF = outDir + VCFOps.getAppropriateRoot(vcf, true) + ".site.only.vcf.gz";
			VCFOps.createSiteOnlyVcf(vcf, outVCF, false, proj.getLog());
			if (!Files.exists(out)) {
				try {
					PrintWriter writer = new PrintWriter(new FileWriter(out));
					VCFFileReader reader = new VCFFileReader(outVCF, true);

					int num = 0;
					for (VariantContext vc : reader) {
						num++;
						if (num % 100000 == 0) {
							proj.getLog().reportTimeInfo(num + " variant sites read, writing to " + out);
						}
						String name = "REF_" + vc.getReference().getDisplayString();
						int i = 0;
						for (Allele allele : vc.getAlternateAlleles()) {
							i++;
							name += "_ALT" + i + "_" + allele.getDisplayString();
						}
						Segment vcSeg = VCOps.getSegment(vc);
						writer.println(vcSeg.getChromosomeUCSC() + "\t" + vcSeg.getStart() + "\t" + vcSeg.getStop() + "\t" + name);
					}
					reader.close();

					writer.close();
				} catch (Exception e) {
					proj.getLog().reportError("Error writing to " + out);
					proj.getLog().reportException(e);
				}

			}

			ArrayList<VariantSeg> segs = new ArrayList<BamImport.VariantSeg>();

			try {

				BufferedReader reader = Files.getAppropriateReader(out);
				while (reader.ready()) {
					String[] line = reader.readLine().trim().split("\t");
					byte chr = Positions.chromosomeNumber(line[0]);
					int start = Integer.parseInt(line[1]);
					int stop = Integer.parseInt(line[2]);
					String tag = line[3];
					segs.add(new VariantSeg(chr, start, stop, tag));
				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				proj.getLog().reportError("Error: file \"" + out + "\" not found in current directory");
				return null;
			} catch (IOException ioe) {
				proj.getLog().reportError("Error reading file \"" + out + "\"");
				return null;
			}
			varLocusSet = new LocusSet<VariantSeg>(segs.toArray(new VariantSeg[segs.size()]), true, proj.getLog()) {

				/**
				 * 
				 */
				private static final long serialVersionUID = 1L;

			};
			// BEDFileReader readerVarBed = new BEDFileReader(out, false);
			// varLocusSet = readerVarBed.loadAll(proj.getLog());
			// readerVarBed.close();

		}

		return varLocusSet;

	}

	public static void importTheWholeBamProject(Project proj, String binBed, String captureBed, String optionalVCF, int captureBuffer, int correctionPCs, int numthreads) {

		if (proj.getArrayType() == ARRAY.NGS) {
			Logger log = proj.getLog();

			String serDir = proj.PROJECT_DIRECTORY.getValue() + "tmpBamSer/";
			String[] bamsToImport = null;
			if (Files.isDirectory(proj.SOURCE_DIRECTORY.getValue())) {
				bamsToImport = Files.listFullPaths(proj.SOURCE_DIRECTORY.getValue(), proj.SOURCE_FILENAME_EXTENSION.getValue(), false);
			} else {
				bamsToImport = HashVec.loadFileToStringArray(proj.SOURCE_DIRECTORY.getValue(), false, new int[] { 0 }, true);
			}

			log.reportTimeInfo("Found " + bamsToImport.length + " bam files to import");
			if (BedOps.verifyBedIndex(binBed, log)) {
				BEDFileReader readerBin = new BEDFileReader(binBed, false);
				LocusSet<BEDFeatureSeg> bLocusSet = readerBin.loadAll(log);
				readerBin.close();
				BEDFileReader readerCapture = new BEDFileReader(captureBed, false);
				readerCapture.close();
				if (!bLocusSet.hasNoOverlap()) {
					ReferenceGenome referenceGenome = new ReferenceGenome(proj.REFERENCE_GENOME_FASTA_FILENAME.getValue(), log);
					log.memoryFree();
					// TODO, skip centromeres
					LocusSet<Segment> genomeBinsMinusBinsCaputure = referenceGenome.getBins(20000).removeThese(LocusSet.combine(bLocusSet, readerCapture.loadAll(log), true, log).mergeOverlapping(true), 4000);//
					log.reportTimeInfo(genomeBinsMinusBinsCaputure.getBpCovered() + " bp covered by reference bins int the anti-on-target regions");
					log.memoryFree();

					LocusSet<VariantSeg> varFeatures = extractVCF(proj, optionalVCF);// get variant sites, piled up according to actual location, not bins
					if (varFeatures.getLoci().length > 0) {
						log.reportTimeInfo(varFeatures.getBpCovered() + " bp covered by known variant sites");
					}

					ArrayList<MarkerFileType> markerTypes = generateMarkerPositions(proj, bLocusSet, genomeBinsMinusBinsCaputure, varFeatures);
					log.memoryFree();
					LocusSet<Segment> analysisSet = LocusSet.combine(bLocusSet.getStrictSegmentSet(), genomeBinsMinusBinsCaputure, true, log);
					analysisSet = LocusSet.combine(analysisSet, varFeatures.getStrictSegmentSet(), true, log);
					dumpLikelyOffTargetProblems(proj);
					log.memoryFree();
					if (!analysisSet.verifyPositiveLength()) {
						throw new IllegalArgumentException("all import segments must be gte length 1");
					}

					FilterNGS filterNGS = new FilterNGS(20, 20, null);// TODO, args for mapQ/phred
					PileupProducer pileProducer = new PileupProducer(bamsToImport, serDir, referenceGenome.getReferenceFasta(), filterNGS, analysisSet.getStrictSegments(), log);
					WorkerTrain<BamPileResult> pileTrain = new WorkerTrain<BamPileResult>(pileProducer, numthreads, 2, log);
					int index = 0;
					proj.SAMPLE_DIRECTORY.getValue(true, false);
					proj.XY_SCALE_FACTOR.setValue((double) 10);

					BamPileResult[] results = new BamPileResult[bamsToImport.length];
					while (pileTrain.hasNext()) {// creating temporary bam pileup of read counts for positions/segments of interest
						results[index] = pileTrain.next();
						index++;

					}
					String[] mappedReadCounts = new String[bamsToImport.length + 1];
					mappedReadCounts[0] = "Sample\tAlignedReadCount\tUnalignedReadCount";
					long fingerPrint = proj.getMarkerSet().getFingerprint();
					BamPileConverterProducer conversionProducer = new BamPileConverterProducer(proj, results, fingerPrint, log);
					WorkerTrain<BamPileConversionResults> conversionTrain = new WorkerTrain<BamImport.BamPileConversionResults>(conversionProducer, numthreads, 10, log);

					Hashtable<String, Float> allOutliers = new Hashtable<String, Float>();

					int convIndex = 0;
					while (conversionTrain.hasNext()) {// normalize read counts and dump to sampRAF, special care for variant sites
						BamPileConversionResults conversionResult = conversionTrain.next();
						BamIndexStats bamIndexStats = conversionResult.getBamIndexStats();
						int numAligned = bamIndexStats.getAlignedRecordCount();
						int numNotAligned = bamIndexStats.getUnalignedRecordCount();
						mappedReadCounts[convIndex + 1] = conversionResult.getSample() + "\t" + numAligned + "\t" + numNotAligned;
						convIndex++;
						if (conversionResult.getOutliers() != null && conversionResult.getOutliers().size() > 0) {
							allOutliers.putAll(conversionResult.getOutliers());
						}
					}
					String readCountFile = proj.PROJECT_DIRECTORY.getValue() + "sample.readCounts.txt";

					Files.writeList(mappedReadCounts, readCountFile);

					if (allOutliers.size() > 0 || !Files.exists(proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser")) {// currently do to all the skipping
						Files.writeSerial(allOutliers, proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser");
					}

					if (!Files.exists(proj.CUSTOM_CENTROIDS_FILENAME.getValue())) {// compute Log R ratio, since its not immediately available
						TransposeData.transposeData(proj, 2000000000, false);
						CentroidCompute.computeAndDumpCentroids(proj, proj.CUSTOM_CENTROIDS_FILENAME.getValue(), new CentroidBuilder(), numthreads, 2);
						Centroids.recompute(proj, proj.CUSTOM_CENTROIDS_FILENAME.getValue(), true, numthreads);
						TransposeData.transposeData(proj, 2000000000, false);
					} else {
						proj.getLog().reportTimeWarning(proj.CUSTOM_CENTROIDS_FILENAME.getValue() + " exists, and currently is the proxy for LRR computation being completed");
					}
					proj.saveProperties();
					// All below stuff is just for fun...

					SampleData.createMinimalSampleData(proj);
					if (!Files.exists(proj.SAMPLE_QC_FILENAME.getValue())) {
						LrrSd.init(proj, null, null, numthreads);
					} else {
						log.reportFileExists(proj.SAMPLE_QC_FILENAME.getValue());
					}
					if (!Files.exists(proj.MOSAIC_RESULTS_FILENAME.getValue())) {
						Mosaicism.findOutliers(proj, numthreads);
					}
					ArrayList<ProjectCorrected> correcteds = correctifyProject(proj, markerTypes, correctionPCs, numthreads);// Generates and corrects the project for each marker type

					String newSampleDir = proj.PROJECT_DIRECTORY.getValue() + "samplesCorrected/";
					String newtransposedDir = proj.PROJECT_DIRECTORY.getValue() + "transposedCorrected/";

					RecompileProducer producer = new RecompileProducer(proj, proj.getSamples(), newSampleDir, proj.getMarkerSet(), correcteds);
					WorkerTrain<Hashtable<String, Float>> train = new WorkerTrain<Hashtable<String, Float>>(producer, numthreads, 10, proj.getLog());
					Hashtable<String, Float> recompallOutliers = new Hashtable<String, Float>();

					while (train.hasNext()) {// consolidate the pc corrected projects back into a single sample
						Hashtable<String, Float> tmp = train.next();
						if (tmp != null) {
							recompallOutliers.putAll(tmp);
						}
					}

					proj.SAMPLE_DIRECTORY.setValue(newSampleDir);// Final resting place
					if (recompallOutliers.size() > 0 || !Files.exists(proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser")) {// currently do to all the skipping
						Files.writeSerial(allOutliers, proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser");
					}
					proj.MARKER_DATA_DIRECTORY.setValue(newtransposedDir);
					TransposeData.transposeData(proj, 2000000000, false);// already recomputed with the correction

					proj.MOSAIC_RESULTS_FILENAME.setValue(ext.addToRoot(proj.MOSAIC_RESULTS_FILENAME.getValue(), ".pcCorrected"));
					if (!Files.exists(proj.MOSAIC_RESULTS_FILENAME.getValue())) {
						Mosaicism.findOutliers(proj, numthreads);
					}

					// generateGCModel(proj, analysisSet, referenceGenome, 50);
					// generateGCModel(proj, analysisSet, referenceGenome, 100);
					// generateGCModel(proj, analysisSet, referenceGenome, 200);
					// generateGCModel(proj, analysisSet, referenceGenome, 500000);
					// generateGCModel(proj, analysisSet, referenceGenome, 1000000);

					//
					// GCAdjustorBuilder gAdjustorBuilder = new GCAdjustorBuilder();
					// GcAdjustorParameter.generate(proj, "GC_ADJUSTMENT/", proj.REFERENCE_GENOME_FASTA_FILENAME.getValue(), gAdjustorBuilder, false, GcAdjustor.GcModel.DEFAULT_GC_MODEL_BIN_FASTA, numthreads);
				} else {
					log.reportTimeError("The bed file " + binBed + " had overlapping segments, currently non -overlapping segments are required");
				}
			}
		} else {
			proj.getLog().reportTimeError(proj.ARRAY_TYPE.getName() + " must be set to " + ARRAY.NGS);
		}
	}

	private static String generateGCModel(Project proj, LocusSet<Segment> analysisSet, ReferenceGenome referenceGenome, int buffer) {
		String gcFile = ext.addToRoot(proj.GC_MODEL_FILENAME.getValue(), ".buffer_" + buffer);
		if (!Files.exists(gcFile)) {
			MarkerSet markerSet = proj.getMarkerSet();
			String[] markerNames = markerSet.getMarkerNames();

			try {
				PrintWriter writer = new PrintWriter(new FileWriter(gcFile));
				String[] header = new String[] { "Name", "Chr", "Position", "GC" };
				writer.println(Array.toStr(header));
				for (int i = 0; i < markerNames.length; i++) {
					if (i % 1000 == 0) {
						proj.getLog().reportTimeInfo("Loaded gc content for " + (i + 1) + " bins");
					}
					writer.println(markerNames[i] + "\t" + markerSet.getChrs()[i] + "\t" + markerSet.getPositions()[i] + "\t" + ReferenceGenome.getPercentGC(referenceGenome.getSequenceFor(analysisSet.getLoci()[i].getBufferedSegment(buffer))));
				}
				writer.close();
			} catch (Exception e) {
				proj.getLog().reportError("Error writing to " + gcFile);
				proj.getLog().reportException(e);
			}
		} else {
			proj.getLog().reportFileExists(gcFile);
		}
		return gcFile;
	}

	private static class ProjectCorrected {
		private Project proj;
		private NGS_MARKER_TYPE type;

		public ProjectCorrected(Project proj, NGS_MARKER_TYPE type) {
			super();
			this.proj = proj;
			this.type = type;
		}

		private Project getProj() {
			return proj;
		}

		private NGS_MARKER_TYPE getType() {
			return type;
		}

	}

	private static ArrayList<ProjectCorrected> correctifyProject(Project proj, ArrayList<MarkerFileType> types, int correctionPCs, int numthreads) {
		proj.SAMPLE_CALLRATE_THRESHOLD.setValue(0.0);
		proj.LRRSD_CUTOFF.setValue(.40);
		proj.INTENSITY_PC_NUM_COMPONENTS.setValue(20);
		String mediaMarks = ext.addToRoot(proj.INTENSITY_PC_MARKERS_FILENAME.getValue(), ".median");
		ArrayList<ProjectCorrected> correctedProjects = new ArrayList<ProjectCorrected>();
		Files.writeList(Array.subArray(proj.getMarkerNames(), 0, 1000), mediaMarks);
		for (MarkerFileType type : types) {
			String base = "";
			if (type.getType() == null) {
				base = base + "BAM_PCS_ALL_MARKERS";
			} else {
				base = base + "BAM_PCS_" + type.getType().getFlag();
			}
			String markerfile = proj.PROJECT_DIRECTORY.getValue() + base + "_inputMarkers.txt";
			Files.writeList(HashVec.loadFileToStringArray(type.getFile(), true, new int[] { 0 }, true), markerfile);
			proj.INTENSITY_PC_MARKERS_FILENAME.setValue(markerfile);
			MitoPipeline.catAndCaboodle(proj, numthreads, mediaMarks, 20, base, false, true, 0, null, null, null, false, false, false, true, false, null, -1, -1);
			// PrincipalComponentsCrossTabs.crossTabulate(proj, proj.INTENSITY_PC_NUM_COMPONENTS.getValue(), null, true);

			String PCCorrected = ext.addToRoot(proj.getPropertyFilename(), "." + proj.INTENSITY_PC_NUM_COMPONENTS.getValue() + "_pc_corrected_" + base);
			String newName = proj.PROJECT_NAME.getValue() + "_" + proj.INTENSITY_PC_NUM_COMPONENTS.getValue() + "_pc_corrected_" + base;
			Files.copyFileUsingFileChannels(proj.getPropertyFilename(), PCCorrected, proj.getLog());
			Project pcCorrected = new Project(PCCorrected, false);
			pcCorrected.PROJECT_DIRECTORY.setValue(proj.PROJECT_DIRECTORY.getValue() + newName + "/");
			pcCorrected.PROJECT_NAME.setValue(newName);
			proj.copyBasicFiles(pcCorrected, true);
			pcCorrected.SAMPLE_DIRECTORY.setValue(pcCorrected.PROJECT_DIRECTORY.getValue() + "shadowSamples/");

			String[] correctedSamps = Array.tagOn(proj.getSamples(), pcCorrected.SAMPLE_DIRECTORY.getValue(), Sample.SAMPLE_DATA_FILE_EXTENSION);
			if (!Files.exists("", correctedSamps)) {
				proj.getLog().reportTimeInfo("PC correcting project using " + proj.INTENSITY_PC_NUM_COMPONENTS.getValue() + " components ");

				PennCNVPrep.exportSpecialPennCNV(proj, "correction/", pcCorrected.PROJECT_DIRECTORY.getValue() + "tmpPCCorrection/", correctionPCs, null, numthreads, 1, false, false, false, -1, true);
				// TODO, auto adjust batch size by memory
				PennCNVPrep.exportSpecialPennCNV(pcCorrected, "correction/", pcCorrected.PROJECT_DIRECTORY.getValue() + "tmpPCCorrection/", correctionPCs, null, 1, 1, false, true, false, 2, true);
			}
			pcCorrected.saveProperties();
			if (type.getType() != null) {
				correctedProjects.add(new ProjectCorrected(pcCorrected, type.getType()));
			}
		}
		return correctedProjects;
	}

	private static class RecompileProducer implements Producer<Hashtable<String, Float>> {
		private Project proj;
		private String[] samples;
		private String newSampleDirectory;
		private MarkerSet markerSet;
		private ArrayList<ProjectCorrected> correctedProjects;
		private int index;

		public RecompileProducer(Project proj, String[] samples, String newSampleDirectory, MarkerSet markerSet, ArrayList<ProjectCorrected> correctedProjects) {
			super();
			this.proj = proj;
			this.samples = samples;
			this.newSampleDirectory = newSampleDirectory;
			this.markerSet = markerSet;
			this.correctedProjects = correctedProjects;
			this.index = 0;
		}

		@Override
		public boolean hasNext() {
			return index < samples.length;

		}

		@Override
		public Callable<Hashtable<String, Float>> next() {
			RecompileWorker current = new RecompileWorker(proj, samples[index], newSampleDirectory, markerSet, correctedProjects);
			index++;
			return current;
		}

		@Override
		public void shutdown() {
			// TODO Auto-generated method stub

		}

	}

	private static class RecompileWorker implements Callable<Hashtable<String, Float>> {
		private Project proj;
		private String sampleName;
		private String newSampleDirectory;
		private MarkerSet markerSet;
		private ArrayList<ProjectCorrected> correctedProjects;

		public RecompileWorker(Project proj, String sampleName, String newSampleDirectory, MarkerSet markerSet, ArrayList<ProjectCorrected> correctedProjects) {
			super();
			this.proj = proj;
			this.sampleName = sampleName;
			this.newSampleDirectory = newSampleDirectory;
			this.markerSet = markerSet;
			this.correctedProjects = correctedProjects;
		}

		@Override
		public Hashtable<String, Float> call() throws Exception {
			return recompileSample(proj, sampleName, newSampleDirectory, markerSet, correctedProjects);
		}

	}

	private static Hashtable<String, Float> recompileSample(Project proj, String sampleName, String newSampleDirectory, MarkerSet markerSet, ArrayList<ProjectCorrected> correctedProjects) {
		String sampleFile = newSampleDirectory + sampleName + Sample.SAMPLE_DATA_FILE_EXTENSION;
		proj.getLog().reportTimeInfo("Sample file = " + sampleFile);
		Hashtable<String, Float> outliers = new Hashtable<String, Float>();
		if (!Files.exists(sampleName)) {
			Sample sampleOriginal = proj.getFullSampleFromRandomAccessFile(sampleName);
			int numAccountedFor = 0;
			float[] gcs = sampleOriginal.getGCs();
			float[] intensity = Array.floatArray(markerSet.getMarkerNames().length, Float.NaN);

			float[] bafs = sampleOriginal.getBAFs();// preserve these;
			float[] lrrs = Array.floatArray(markerSet.getMarkerNames().length, Float.NaN);

			String[] markerNames = markerSet.getMarkerNames();
			for (ProjectCorrected corrected : correctedProjects) {
				Sample typeCorrected = corrected.getProj().getFullSampleFromRandomAccessFile(sampleName);
				for (int i = 0; i < markerNames.length; i++) {
					NGS_MARKER_TYPE marker_TYPE = NGS_MARKER_TYPE.getType(markerNames[i]);
					if (marker_TYPE == corrected.getType()) {
						intensity[i] = typeCorrected.getXs()[i];
						lrrs[i] = typeCorrected.getLRRs()[i];
						numAccountedFor++;
					}
				}
			}
			Sample sampleCorrected = new Sample(sampleName, markerSet.getFingerprint(), gcs, intensity, intensity, bafs, lrrs, sampleOriginal.getForwardGenotypes(), sampleOriginal.getAB_Genotypes(), sampleOriginal.getCanXYBeNegative());
			sampleCorrected.saveToRandomAccessFile(sampleFile, outliers, sampleName);
			if (numAccountedFor != markerNames.length) {
				throw new IllegalArgumentException("Not all markers accounted for in corrections");
			}
		} else {
			proj.getLog().reportFileExists(sampleFile);
		}

		return outliers;
	}

	private static void dumpLikelyOffTargetProblems(Project proj) {
		MarkerSet markerSet = proj.getMarkerSet();
		String problemFile = ext.addToRoot(proj.MARKER_POSITION_FILENAME.getValue(), ".likelyOffTargetProblems");
		String noproblemFile = ext.addToRoot(proj.MARKER_POSITION_FILENAME.getValue(), ".withoutlikelyOffTargetProblems");
		String allFile = ext.addToRoot(proj.MARKER_POSITION_FILENAME.getValue(), ".OffTargetProblemsFlagged");

		ArrayList<String> problems = new ArrayList<String>();
		ArrayList<String> noProblems = new ArrayList<String>();
		ArrayList<String> all = new ArrayList<String>();

		problems.add("BinName\tCLASS=MARKER_COLOR;OFF_TARGET_OK=Blue;LIKELY_OFF_TARGET_PROBLEM=RED;OTHER_TYPE=Green");
		noProblems.add("BinName\tCLASS=MARKER_COLOR;OFF_TARGET_OK=Blue;LIKELY_OFF_TARGET_PROBLEM=RED;OTHER_TYPE=Green");
		all.add("BinName\tCLASS=MARKER_COLOR;OFF_TARGET_OK=Blue;LIKELY_OFF_TARGET_PROBLEM=RED;OTHER_TYPE=Green");
		int[][] indices = markerSet.getIndicesByChr();
		String[] names = markerSet.getMarkerNames();
		for (int i = 0; i < indices.length; i++) {
			for (int j = 0; j < indices[i].length; j++) {
				int compLeft = Math.max(j - 1, 0);
				int compRight = Math.min(j + 1, indices[i].length - 1);
				NGS_MARKER_TYPE current = NGS_MARKER_TYPE.getType(names[indices[i][j]]);
				if (current == NGS_MARKER_TYPE.OFF_TARGET) {
					NGS_MARKER_TYPE left = NGS_MARKER_TYPE.getType(names[indices[i][compLeft]]);
					NGS_MARKER_TYPE right = NGS_MARKER_TYPE.getType(names[indices[i][compRight]]);
					if ((compLeft != j && left != NGS_MARKER_TYPE.OFF_TARGET) || (compRight != j && right != NGS_MARKER_TYPE.OFF_TARGET)) {
						problems.add(names[indices[i][j]] + "\tLIKELY_OFF_TARGET_PROBLEM");
						all.add(names[indices[i][j]] + "\tLIKELY_OFF_TARGET_PROBLEM");
					} else {
						noProblems.add(names[indices[i][j]] + "\tOFF_TARGET_OK");
						all.add(names[indices[i][j]] + "\tOFF_TARGET_OK");
					}
				} else {
					noProblems.add(names[indices[i][j]] + "\tOTHER_TYPE");
					all.add(names[indices[i][j]] + "\tOTHER_TYPE");
				}
			}
		}
		proj.getLog().reportTimeInfo("Dumping " + problems.size() + " off target markers that will likely be biased to " + problemFile);
		Files.writeArrayList(problems, problemFile);
		Files.writeArrayList(noProblems, noproblemFile);
		Files.writeArrayList(all, allFile);

	}

	private static ArrayList<MarkerFileType> generateMarkerPositions(Project proj, LocusSet<BEDFeatureSeg> bLocusSet, LocusSet<Segment> genomeBinsMinusBinsCaputure, LocusSet<VariantSeg> varFeatures) {
		String positions = proj.MARKER_POSITION_FILENAME.getValue();
		proj.getLog().reportTimeInfo("Postions will be set to the midpoint of each segment");
		String[] markerNames = new String[bLocusSet.getLoci().length + genomeBinsMinusBinsCaputure.getLoci().length + varFeatures.getLoci().length];
		String header = "BinName\tChr\tPosition\tCLASS=MARKER_COLOR;OFF_TARGET=Blue;VARIANT_SITE=RED;ON_TARGET=Green";
		ArrayList<String> onTMarkers = new ArrayList<String>();
		onTMarkers.add(header);

		ArrayList<String> offTMarkers = new ArrayList<String>();
		offTMarkers.add(header);

		ArrayList<String> variantSiteMarkers = new ArrayList<String>();
		variantSiteMarkers.add(header);

		try {
			PrintWriter writer = new PrintWriter(new FileWriter(positions));
			int markerIndex = 0;
			writer.println(header);

			for (int i = 0; i < bLocusSet.getLoci().length; i++) {
				BEDFeatureSeg bFeatureSeg = bLocusSet.getLoci()[i];
				String markerName = bFeatureSeg.getUCSClocation();
				String name = bFeatureSeg.getBedFeature().getName();
				if (name != null) {
					markerName += "|" + name;
				}
				markerName += "|" + NGS_MARKER_TYPE.ON_TARGET.getFlag();
				markerNames[markerIndex] = markerName;
				int diff = bFeatureSeg.getStop() - bFeatureSeg.getStart();
				int mid = Math.round((float) diff / 2);
				int pos = bFeatureSeg.getStart() + mid;
				String out = markerName + "\t" + bFeatureSeg.getChr() + "\t" + pos + "\t" + NGS_MARKER_TYPE.ON_TARGET.getFlag();
				onTMarkers.add(out);
				writer.println(out);
				markerIndex++;
			}
			for (int i = 0; i < genomeBinsMinusBinsCaputure.getLoci().length; i++) {
				Segment binnedSeg = genomeBinsMinusBinsCaputure.getLoci()[i];
				String markerName = binnedSeg.getUCSClocation() + "|" + NGS_MARKER_TYPE.OFF_TARGET.getFlag();

				markerNames[markerIndex] = markerName;
				int diff = binnedSeg.getStop() - binnedSeg.getStart();
				int mid = Math.round((float) diff / 2);
				int pos = binnedSeg.getStart() + mid;
				String out = markerName + "\t" + binnedSeg.getChr() + "\t" + pos + "\t" + NGS_MARKER_TYPE.OFF_TARGET.getFlag();
				offTMarkers.add(out);
				writer.println(out);

				markerIndex++;
			}

			for (int i = 0; i < varFeatures.getLoci().length; i++) {
				VariantSeg variantFeatureSeg = varFeatures.getLoci()[i];
				String markerName = variantFeatureSeg.getUCSClocation();
				String name = variantFeatureSeg.getTag();
				if (name != null) {
					markerName += "|" + name;
				}
				markerName = markerName + "|" + NGS_MARKER_TYPE.VARIANT_SITE.getFlag();
				markerNames[markerIndex] = markerName;
				int pos = variantFeatureSeg.getStart();
				String out = markerName + "\t" + variantFeatureSeg.getChr() + "\t" + pos + "\t" + NGS_MARKER_TYPE.VARIANT_SITE.getFlag();
				variantSiteMarkers.add(out);

				writer.println(out);
				markerIndex++;
			}

			writer.close();
		} catch (Exception e) {
			proj.getLog().reportError("Error writing to " + positions);
			proj.getLog().reportException(e);
		}
		String onTargetFile = ext.addToRoot(proj.MARKER_POSITION_FILENAME.getValue(), "." + NGS_MARKER_TYPE.ON_TARGET.getFlag());
		String offTargetFile = ext.addToRoot(proj.MARKER_POSITION_FILENAME.getValue(), "." + NGS_MARKER_TYPE.OFF_TARGET.getFlag());
		String variantSiteTargetFile = ext.addToRoot(proj.MARKER_POSITION_FILENAME.getValue(), "." + NGS_MARKER_TYPE.VARIANT_SITE.getFlag());
		String allMarkerFile = ext.addToRoot(proj.MARKER_POSITION_FILENAME.getValue(), ".allMarkers");
		ArrayList<MarkerFileType> markerTypes = new ArrayList<MarkerFileType>();

		if (onTMarkers.size() > 0) {
			Files.writeList(Array.toStringArray(onTMarkers), onTargetFile);
			markerTypes.add(new MarkerFileType(NGS_MARKER_TYPE.ON_TARGET, onTargetFile));
		} else {
			proj.getLog().reportTimeWarning("No " + NGS_MARKER_TYPE.ON_TARGET.getFlag() + " markers detected");
		}
		if (offTMarkers.size() > 0) {
			Files.writeList(Array.toStringArray(offTMarkers), offTargetFile);
			markerTypes.add(new MarkerFileType(NGS_MARKER_TYPE.OFF_TARGET, offTargetFile));
		} else {
			proj.getLog().reportTimeWarning("No " + NGS_MARKER_TYPE.OFF_TARGET.getFlag() + " markers detected");
		}

		if (variantSiteMarkers.size() > 0) {
			Files.writeList(Array.toStringArray(variantSiteMarkers), variantSiteTargetFile);
			markerTypes.add(new MarkerFileType(NGS_MARKER_TYPE.VARIANT_SITE, variantSiteTargetFile));
		} else {
			proj.getLog().reportTimeWarning("No " + NGS_MARKER_TYPE.VARIANT_SITE.getFlag() + " markers detected");
		}

		Files.writeList(HashVec.loadFileToStringArray(proj.MARKER_POSITION_FILENAME.getValue(), true, new int[] { 0 }, true), allMarkerFile);
		markerTypes.add(new MarkerFileType(null, allMarkerFile));

		Markers.orderMarkers(markerNames, proj.MARKER_POSITION_FILENAME.getValue(), proj.MARKERSET_FILENAME.getValue(true, true), proj.getLog());
		return markerTypes;
	}

	private static class MarkerFileType {
		private NGS_MARKER_TYPE type;
		private String file;

		public MarkerFileType(NGS_MARKER_TYPE type, String file) {
			super();
			this.type = type;
			this.file = file;
		}

		private NGS_MARKER_TYPE getType() {
			return type;
		}

		private String getFile() {
			return file;
		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		String binBed = "binsToImport.bed";
		String captureBed = "AgilentCaptureRegions.txt";
		int numthreads = 24;
		int captureBuffer = 400;
		String vcf = null;
		int correctionPCs = 4;
		// String referenceGenomeFasta = "hg19_canonical.fa";
		// String logfile = null;
		// Logger log;

		String usage = "\n" + "seq.manage.BamImport requires 0-1 arguments\n";
		usage += "(1) filename (i.e. proj= ( nodefault))\n" + "";
		usage += "(2) bed file to import  (i.e. importBed=" + binBed + " ( no default))\n" + "";
		usage += Ext.getNumThreadsCommand(3, numthreads);
		usage += "(4) bed file to import  (i.e. captureBed=" + captureBed + " ( no default))\n" + "";
		usage += "(5) a vcf, if provided the variants will be imported with bp resolution  (i.e. vcf= ( no default))\n" + "";
		usage += "(6) number of PCs to correct with  (i.e. correctionPCs= ( no default))\n" + "";

		// usage += "(3) reference genome  (i.e. ref=" + referenceGenomeFasta + " ( default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("importBed=")) {
				binBed = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("captureBed=")) {
				captureBed = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("vcf=")) {
				vcf = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith(Ext.NUM_THREADS_COMMAND)) {
				numthreads = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("correctionPCs")) {
				correctionPCs = ext.parseIntArg(args[i]);
				numArgs--;
			}
			// else if (args[i].startsWith("log=")) {
			// logfile = args[i].split("=")[1];
			// numArgs--;
			// }
			else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			// log = new Logger(logfile);
			Project proj = new Project(filename, false);
			importTheWholeBamProject(proj, binBed, captureBed, vcf, captureBuffer, correctionPCs, numthreads);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	//
	// private BamPileUp bamPileUp;
	// private String vcf;
	// private String bam;
	// private FilterNGS readDepthFilter;
	// private LocusSet<Segment> intervals;
	// private ReferenceGenome referenceGenome;
	// private Logger log;
	//
	// public BamImport(String vcf, String bam, FilterNGS readDepthFilter, LocusSet<Segment> intervals, ReferenceGenome referenceGenome, Logger log) {
	// super();
	// this.vcf = vcf;
	// this.bam = bam;
	// this.readDepthFilter = readDepthFilter;
	// this.referenceGenome = referenceGenome;
	// this.bamPileUp = new BamPileUp(bam, referenceGenome, 1, readDepthFilter, intervals.getLoci(), PILE_TYPE.REGULAR, SAM_FILTER_TYPE.COPY_NUMBER, true, log);
	// this.intervals = new LocusSet<Segment>(BamOps.converQItoSegs(bamPileUp.getQueryIntervals(), BamOps.getHeader(bam), log), true, log) {
	//
	// /**
	// *
	// */
	// private static final long serialVersionUID = 1L;
	// };
	//
	// this.log = log;
	// }
	//
	// public void importBam() {
	// VCFFileReader reader = new VCFFileReader(vcf, true);
	// String bamSample = BamOps.getSampleName(bam);
	//
	// if (ext.indexOfStr(bamSample, VCFOps.getSamplesInFile(reader)) < 0) {
	// log.reportTimeError("Could not find sample " + bamSample + " in the vcf " + vcf);
	// return;
	// } else {
	// log.reportTimeInfo("Detected sample " + bamSample + " in vcf " + vcf);
	// }
	// log.reportTimeWarning("Only un-ambigous and biallelic variants will be imported from " + vcf);
	// FilterNGS.VariantContextFilter niceAllele = new FilterNGS.VariantContextFilter(new VARIANT_FILTER_DOUBLE[] {}, new VARIANT_FILTER_BOOLEAN[] { VARIANT_FILTER_BOOLEAN.BIALLELIC_FILTER, VARIANT_FILTER_BOOLEAN.AMBIGUOUS_FILTER }, null, null, log);
	// SampleNGS ngsSample = new SampleNGS(bamSample);
	// TmpBin tmpBin = null;
	// int scanIndex = -1;
	// while (bamPileUp.hasNext()) {
	// BamPile bamPile = bamPileUp.next();
	// Segment curSeg = bamPile.getBin();
	// int[] indices = intervals.getOverlappingIndices(curSeg);
	// System.out.println("HDI\t" + curSeg.getUCSClocation() + "\t" + intervals.getLoci()[0].getUCSClocation());
	// if (indices == null || indices.length == 0) {
	// // log.reportTimeInfo("Un-matched segments");
	// } else {
	// if (indices.length > 1) {
	// log.reportTimeInfo("Non-unique (overlapping) segments were supplied, halting");
	// return;
	// }
	// Segment[] overlaps = Array.subArray(intervals.getLoci(), indices);
	// if (tmpBin == null) {
	// tmpBin = new TmpBin(overlaps[0]);
	// }
	// if (!tmpBin.getCurrentPile().equals(overlaps[0])) {
	// ngsSample.addGeno(null, tmpBin.developFakeGenotype(), log);
	// tmpBin = new TmpBin(overlaps[0]);
	// CloseableIterator<VariantContext> reg = reader.query(Positions.getChromosomeUCSC(overlaps[0].getChr(), true), overlaps[0].getStart(), overlaps[0].getStop());
	// while (reg.hasNext()) {
	// VariantContext vc = reg.next();
	// if (niceAllele.filter(vc).passed()) {
	// ngsSample.addGeno(vc, VCOps.getGenotypeFor(vc, bamSample, VC_SUBSET_TYPE.SUBSET_STRICT), log);
	// }
	// }
	// }
	// tmpBin.setNumRef(tmpBin.getNumRef() + bamPile.getNumRef(log));
	// tmpBin.setNumAlt(tmpBin.getNumAlt() + bamPile.getNumAlt(log));
	// }
	// }
	// }
	//
	// // int index = 0;
	// // boolean newIndex = true;
	// // while (bamPileUp.hasNext() && index < intervals.getLoci().length) {
	// // BamPile bamPile = bamPileUp.next();
	// // Segment curSeg = bamPile.getBin();
	// // Segment[] overlaps = intervals.getOverLappingLoci(curSeg);
	// // if (newIndex) {
	// // CloseableIterator<VariantContext> reg = reader.query(Positions.getChromosomeUCSC(intervals.getLoci()[index].getChr(), true), intervals.getLoci()[index].getStart(), intervals.getLoci()[index].getStop());
	// // while (reg.hasNext()) {
	// // VariantContext vc = reg.next();
	// // if (niceAllele.filter(vc).passed()) {
	// // ngsSample.addGeno(VCOps.getGenotypeFor(vc, bamSample, VC_SUBSET_TYPE.SUBSET_STRICT), log);
	// // }
	// // }
	// // newIndex = false;
	// // }
	// // if (overlaps == null || overlaps.length == 0) {
	// // if (curSeg.getChr() == intervals.getLoci()[index].getChr() && curSeg.getStop() < intervals.getLoci()[index].getStart()) {// before current index
	// //
	// // } else {
	// // while (curSeg.getChr() == intervals.getLoci()[index].getChr() && curSeg.getStart() > intervals.getLoci()[index].getStop()) {
	// // index++;
	// // }
	// // }
	// // }
	// //
	// // }
	// //
	// //
	// private static class TmpBin {
	// private Segment currentPile;
	// private int numRef;
	// private int numAlt;
	//
	// public Segment getCurrentPile() {
	// return currentPile;
	// }
	//
	// public TmpBin(Segment currentPile) {
	// super();
	// this.currentPile = currentPile;
	// this.numRef = 0;
	// this.numAlt = 0;
	// }
	//
	// public int getNumRef() {
	// return numRef;
	// }
	//
	// public void setNumRef(int numRef) {
	// this.numRef = numRef;
	// }
	//
	// public int getNumAlt() {
	// return numAlt;
	// }
	//
	// public void setNumAlt(int numAlt) {
	// this.numAlt = numAlt;
	// }
	//
	// public void setCurrentPile(Segment currentPile) {
	// this.currentPile = currentPile;
	// }
	//
	// public Genotype developFakeGenotype() {
	// GenotypeBuilder builder = new GenotypeBuilder();
	// builder.AD(new int[] { numRef, numAlt });
	// builder.GQ(100);
	// ArrayList<Allele> alleles = new ArrayList<Allele>();
	// alleles.add(Allele.create("N", true));
	// builder.alleles(alleles);
	// return builder.make();
	// }
	//
	// }
	//
	// public static void test() {
	// String bam = "D:/data/Project_Tsai_Project_021/testPileUp/rrd_lane_CONTROL_4_CTCTCTAC-CTCTCTAT.merge.sorted.dedup.realigned.bam";
	// String ref = "C:/bin/ref/hg19_canonical.fa";
	// String segs = "C:/bin/Agilent/captureLibraries/SureSelectHumanAllExonV5UTRs/AgilentCaptureRegions_chr1.txt";
	// String vcf = "D:/data/Project_Tsai_Spector_Joint/joint_genotypes_tsai_21_25_spector_mt.AgilentCaptureRegions.SNP.recal.INDEL.recal.hg19_multianno.eff.gatk.vcf.gz";
	// Logger log = new Logger();
	// Segment[] q = segs == null ? null : Segment.loadRegions(segs, 0, 1, 2, 0, true, true, true, 100);
	// FilterNGS filterNGS = new FilterNGS(0, 0, new int[] { 0, 0 });
	// ReferenceGenome referenceGenome = ref == null ? null : new ReferenceGenome(ref, log);
	// LocusSet<Segment> locusSet = new LocusSet<Segment>(q, true, log) {
	//
	// /**
	// *
	// */
	// private static final long serialVersionUID = 1L;
	//
	// };
	// BamImport bamImport = new BamImport(vcf, bam, filterNGS, locusSet, referenceGenome, log);
	// bamImport.importBam();
	//
	// }
	//
	// public static void main(String[] args) {
	// String bam = "D:/data/Project_Tsai_Project_021/testPileUp/rrd_lane_CONTROL_4_CTCTCTAC-CTCTCTAT.merge.sorted.dedup.realigned.bam";
	// String ref = "C:/bin/ref/hg19_canonical.fa";
	// String segs = "C:/bin/Agilent/captureLibraries/SureSelectHumanAllExonV5UTRs/AgilentCaptureRegions_chr1.txt";
	// String vcf = "D:/data/Project_Tsai_Spector_Joint/joint_genotypes_tsai_21_25_spector_mt.AgilentCaptureRegions.SNP.recal.INDEL.recal.hg19_multianno.eff.gatk.vcf.gz";
	// test();
	// }
}
