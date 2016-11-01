package org.genvisis.seq.manage;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;
import java.util.concurrent.Callable;

import org.genvisis.cnv.analysis.CentroidCompute;
import org.genvisis.cnv.analysis.CentroidCompute.CentroidBuilder;
import org.genvisis.cnv.analysis.Mosaicism;
import org.genvisis.cnv.analysis.PennCNVPrep;
import org.genvisis.cnv.filesys.Centroids;
import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.manage.Markers;
import org.genvisis.cnv.manage.MitoPipeline;
import org.genvisis.cnv.manage.Resources.GENOME_BUILD;
import org.genvisis.cnv.manage.TransposeData;
import org.genvisis.cnv.qc.LrrSd;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF.Ext;
import org.genvisis.common.Positions;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.WorkerTrain.AbstractProducer;
import org.genvisis.common.ext;
import org.genvisis.filesys.LocusSet;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.SeqVariables.ASSAY_TYPE;
import org.genvisis.seq.SeqVariables.ASSEMBLY_NAME;
import org.genvisis.seq.manage.BEDFileReader.BEDFeatureSeg;
import org.genvisis.seq.manage.BamOps.BamIndexStats;
import org.genvisis.seq.manage.BamSegPileUp.BamPileResult;
import org.genvisis.seq.manage.BamSegPileUp.PileupProducer;
import org.genvisis.seq.qc.FilterNGS;
import org.genvisis.stats.LeastSquares.LS_TYPE;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

/**
 * Class for converting bam files to Genvisis compatable format
 *
 */
public class BamImport {
	public static final int CAPTURE_BUFFER = 4000;
	public static final int BIN_SIZE_WXS = 20000;
	public static final int BIN_SIZE_WGS = 2000;



	/**
	 * Some enums that define NGS "marker" types
	 *
	 */
	/**
	 * @author Kitty
	 *
	 */
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

		private final String flag;

		private NGS_MARKER_TYPE(String flag) {
			this.flag = flag;
		}

		public String getFlag() {
			return flag;
		}

		/**
		 * @param markerName
		 * @return the {@link NGS_MARKER_TYPE} this marker represents
		 */
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
		private final Project proj;
		private final BamPileResult result;
		private String sample;
		private BamIndexStats bamIndexStats;
		private Hashtable<String, Float> outliers;
		private final Logger log;
		private final long fingerPrint;

		public BamPileConversionResults(Project proj, BamPileResult result, long fingerPrint,
																		Logger log) {
			super();
			this.proj = proj;
			this.result = result;
			outliers = new Hashtable<String, Float>();
			this.fingerPrint = fingerPrint;
			this.log = log;
		}

		@Override
		public BamPileConversionResults call() throws Exception {
			String sampleFile = proj.SAMPLE_DIRECTORY.getValue()
													+ BamOps.getSampleName(result.getBam(), log)
													+ Sample.SAMPLE_FILE_EXTENSION;
			if (!Files.exists(sampleFile)) {
				BamSample bamSample = new BamSample(proj, result.getBam(), result.loadResults(log));
				sample = bamSample.getSampleName();
				bamIndexStats = BamOps.getBamIndexStats(result.getBam());
				outliers = bamSample.writeSample(fingerPrint);
			} else {
				log.reportFileExists(sampleFile);
				sample = BamOps.getSampleName(result.getBam(), log);
				bamIndexStats = BamOps.getBamIndexStats(result.getBam());
				outliers = null;
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

	private static class BamPileConverterProducer extends AbstractProducer<BamPileConversionResults> {
		private final Project proj;
		private final BamPileResult[] pileResults;
		private final long fingerPrint;
		private final Logger log;
		private int index;

		private BamPileConverterProducer(	Project proj, BamPileResult[] pileResults, long fingerPrint,
																			Logger log) {
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
	}

	private static class VariantSeg extends Segment {

		/**
		 *
		 */
		private static final long serialVersionUID = 1L;
		private final String tag;

		private VariantSeg(byte chr, int start, int stop, String tag) {
			super(chr, start, stop);
			this.tag = tag;
		}

		private String getTag() {
			return tag;
		}

	}

	private static LocusSet<VariantSeg> extractVCF(	Project proj, String vcf,
																									boolean removeSingletons) {
		LocusSet<VariantSeg> varLocusSet = new LocusSet<VariantSeg>(new VariantSeg[] {}, true,
																																proj.getLog()) {

			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

		};

		if (vcf != null) {
			proj.getLog().reportTimeInfo("Extracting variant sites from " + vcf);
			String outDir = proj.PROJECT_DIRECTORY.getValue() + "vcf/";
			new File(outDir).mkdirs();
			String out = outDir + VCFOps.getAppropriateRoot(vcf, true) + ".regions.txt";
			String outVCF = outDir + VCFOps.getAppropriateRoot(vcf, true) + ".site.only.vcf.gz";
			VCFOps.createSiteOnlyVcf(vcf, outVCF, removeSingletons, false, proj.getLog());
			if (!Files.exists(out)) {
				try {
					PrintWriter writer = new PrintWriter(new FileWriter(out));
					VCFFileReader reader = new VCFFileReader(new File(outVCF), true);

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
						writer.println(vcSeg.getChromosomeUCSC()+ "\t" + vcSeg.getStart() + "\t"
														+ vcSeg.getStop() + "\t" + name);
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
				proj.getLog().reportException(fnfe);
				return null;
			} catch (IOException ioe) {
				proj.getLog().reportError("Error reading file \"" + out + "\"");
				proj.getLog().reportException(ioe);

				return null;
			}
			varLocusSet = new LocusSet<VariantSeg>(	segs.toArray(new VariantSeg[segs.size()]), true,
																							proj.getLog()) {

				/**
				 * 
				 */
				private static final long serialVersionUID = 1L;

			};

		}

		return varLocusSet;

	}

	public static void importTheWholeBamProject(Project proj, String binBed, String captureBed,
																							String optionalVCF, int captureBuffer,
																							int correctionPCs, boolean compileProject,
																							ASSAY_TYPE atType, ASSEMBLY_NAME aName,
																							String[] bamsToImport, int numthreads) {

		if (proj.getArrayType() == ARRAY.NGS) {
			Logger log = proj.getLog();

			String serDir = proj.PROJECT_DIRECTORY.getValue() + "tmpBamSer/";
			if (bamsToImport == null) {
				if (Files.isDirectory(proj.SOURCE_DIRECTORY.getValue())) {
					bamsToImport = Files.listFullPaths(	proj.SOURCE_DIRECTORY.getValue(),
																							proj.SOURCE_FILENAME_EXTENSION.getValue(), false);
				} else {
					bamsToImport = HashVec.loadFileToStringArray(	proj.SOURCE_DIRECTORY.getValue(), false,
																												new int[] {0}, true);
				}
			}
			ReferenceGenome referenceGenome =
																			new ReferenceGenome(proj.getReferenceGenomeFASTAFilename(),
																													log);
			log.reportTimeInfo("Found " + bamsToImport.length + " bam files to import");
			AnalysisSets analysisSet = generateAnalysisSet(	proj, binBed, captureBed, optionalVCF,
																											captureBuffer, atType, log, referenceGenome);
			FilterNGS filterNGS = new FilterNGS(20, 20, null);
			PileupProducer pileProducer = new PileupProducer(	bamsToImport, serDir,
																												referenceGenome.getReferenceFasta(),
																												filterNGS,
																												analysisSet.analysisSet.getStrictSegments(),
																												aName, log);
			WorkerTrain<BamPileResult> pileTrain = new WorkerTrain<BamPileResult>(pileProducer,
																																						numthreads, 2, log);
			int index = 0;
			proj.SAMPLE_DIRECTORY.getValue(true, false);
			proj.XY_SCALE_FACTOR.setValue((double) 10);

			BamPileResult[] results = new BamPileResult[bamsToImport.length];
			while (pileTrain.hasNext()) {// creating temporary bam
																		// pileup of read counts for
																		// positions/segments of
																		// interest
				results[index] = pileTrain.next();
				index++;

			}
			if (compileProject) {
				compileProject(	proj, correctionPCs, numthreads, log, bamsToImport, referenceGenome,
												analysisSet.markerTypes, analysisSet.analysisSet,
												analysisSet.offTargetsToUse, results, aName);
			}


		} else {
			proj.getLog().reportTimeError(proj.ARRAY_TYPE.getName() + " must be set to " + ARRAY.NGS);
		}
	}

	/**
	 * @param proj determining where everything will live
	 * @param binBed the bins of interest
	 * @param captureBed capture library
	 * @param optionalVCF vcf for variants
	 * @param log
	 * @param referenceGenome
	 * @return
	 */
	public static AnalysisSets generateAnalysisSet(	Project proj, String binBed, String captureBed,
																									String optionalVCF, int captureBuffer,
																									ASSAY_TYPE atype, Logger log,

																									ReferenceGenome referenceGenome) {
		LocusSet<Segment> analysisSet;
		LocusSet<BEDFeatureSeg> bLocusSet = null;
		LocusSet<Segment> genomeBinsMinusBinsCaputure = null;
		if (atype == ASSAY_TYPE.WXS) {
			if (binBed == null || BedOps.verifyBedIndex(binBed, log)) {
				BEDFileReader readerBin = new BEDFileReader(binBed, false);
				bLocusSet = readerBin.loadAll(log);
				readerBin.close();
				BEDFileReader readerCapture = new BEDFileReader(captureBed, false);
				readerCapture.close();
				if (!bLocusSet.hasNoOverlap()) {

					log.memoryFree();
					genomeBinsMinusBinsCaputure = referenceGenome	.getBins(BIN_SIZE_WXS)
																												.removeThese(	LocusSet.combine(bLocusSet,
																																											readerCapture.loadAll(log),
																																											true, log)
																																							.mergeOverlapping(true),
																																			captureBuffer);//
				} else {
					log.reportTimeError("The bed file "+ binBed
															+ " had overlapping segments, currently non -overlapping segments are required");
					throw new IllegalArgumentException();
				}
			}
		} else if (atype == ASSAY_TYPE.WGS) {
			genomeBinsMinusBinsCaputure = referenceGenome.getBins(BIN_SIZE_WGS);
			String tmp = proj.PROJECT_DIRECTORY.getValue() + "tmp.bed";
			Files.write("0\t1\t100000", tmp);
			BEDFileReader readerBin = new BEDFileReader(tmp, false);
			bLocusSet = readerBin.loadAll(log);
			readerBin.close();
		} else {
			throw new IllegalArgumentException("Invalid assay type " + atype);
		}
		log.reportTimeInfo(genomeBinsMinusBinsCaputure.getBpCovered()+ " bp covered by " + atype
												+ " analysis  set");
		log.memoryFree();

		String vcfToCount = optionalVCF;
		if (atype == ASSAY_TYPE.WXS) {
			String subsetVcf = VCFOps.extractSegments(vcfToCount, captureBed, 100, null,
																								proj.PROJECT_DIRECTORY.getValue() + "vcf/", false,
																								false, 1, log);
			vcfToCount = subsetVcf;
		}
		LocusSet<VariantSeg> varFeatures = extractVCF(proj, vcfToCount, true);

		if (varFeatures.getLoci().length > 0) {
			log.reportTimeInfo(varFeatures.getBpCovered() + " bp covered by known variant sites");
		}

		ArrayList<MarkerFileType> markerTypes = generateMarkerPositions(proj, bLocusSet,
																																		genomeBinsMinusBinsCaputure,
																																		varFeatures);
		log.memoryFree();
		analysisSet = LocusSet.combine(	bLocusSet.getStrictSegmentSet(), genomeBinsMinusBinsCaputure,
																		true, log);
		analysisSet = LocusSet.combine(analysisSet, varFeatures.getStrictSegmentSet(), true, log);
		String[] offTargetsToUse = dumpLikelyOffTargetProblems(proj, atype);
		log.memoryFree();
		if (!analysisSet.verifyPositiveLength()) {
			throw new IllegalArgumentException("all import segments must be gte length 1");
		}
		return new AnalysisSets(analysisSet, offTargetsToUse, markerTypes);
	}

	private static class AnalysisSets {
		private final LocusSet<Segment> analysisSet;
		private final String[] offTargetsToUse;
		private final ArrayList<MarkerFileType> markerTypes;

		private AnalysisSets(	LocusSet<Segment> analysisSet, String[] offTargetsToUse,
													ArrayList<MarkerFileType> markerTypes) {
			super();
			this.analysisSet = analysisSet;
			this.offTargetsToUse = offTargetsToUse;
			this.markerTypes = markerTypes;
		}


	}



	private static void compileProject(	Project proj, int correctionPCs, int numthreads, Logger log,
																			String[] bamsToImport, ReferenceGenome referenceGenome,
																			ArrayList<MarkerFileType> markerTypes,
																			LocusSet<Segment> analysisSet, String[] offTargetsToUse,
																			BamPileResult[] results, ASSEMBLY_NAME aName) {
		String[] mappedReadCounts = new String[bamsToImport.length + 1];
		mappedReadCounts[0] = "Sample\tAlignedReadCount\tUnalignedReadCount";
		long fingerPrint = proj.getMarkerSet().getFingerprint();
		BamPileConverterProducer conversionProducer = new BamPileConverterProducer(	proj, results,
																																								fingerPrint, log);
		WorkerTrain<BamPileConversionResults> conversionTrain =
																													new WorkerTrain<BamImport.BamPileConversionResults>(conversionProducer,
																																																							numthreads,
																																																							10,
																																																							log);

		Hashtable<String, Float> allOutliers = new Hashtable<String, Float>();

		int convIndex = 0;
		while (conversionTrain.hasNext()) {// normalize read counts
																				// and dump to sampRAF,
																				// special care for
																				// variant sites
			BamPileConversionResults conversionResult = conversionTrain.next();
			BamIndexStats bamIndexStats = conversionResult.getBamIndexStats();
			int numAligned = bamIndexStats.getAlignedRecordCount();
			int numNotAligned = bamIndexStats.getUnalignedRecordCount();
			mappedReadCounts[convIndex + 1] = conversionResult.getSample()+ "\t" + numAligned + "\t"
																				+ numNotAligned;
			convIndex++;
			if (conversionResult.getOutliers() != null && conversionResult.getOutliers().size() > 0) {
				allOutliers.putAll(conversionResult.getOutliers());
			}
		}
		String readCountFile = proj.PROJECT_DIRECTORY.getValue() + "sample.readCounts.txt";

		Files.writeArray(mappedReadCounts, readCountFile);

		if (allOutliers.size() > 0
				|| !Files.exists(proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser")) {// currently
																																												// do
																																												// to
																																												// all
																																												// the
																																												// skipping
			SerializedFiles.writeSerial(allOutliers,
																	proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser");
		}

		if (!Files.exists(proj.CUSTOM_CENTROIDS_FILENAME.getValue())) {// compute
																																		// Log
																																		// R
																																		// ratio,
																																		// since
																																		// its
																																		// not
																																		// immediately
																																		// available
			TransposeData.transposeData(proj, 2000000000, false);
			CentroidCompute.computeAndDumpCentroids(proj, proj.CUSTOM_CENTROIDS_FILENAME.getValue(),
																							new CentroidBuilder(), numthreads, 2);
			Centroids.recompute(proj, proj.CUSTOM_CENTROIDS_FILENAME.getValue(), true, numthreads);
			TransposeData.transposeData(proj, 2000000000, false);
		} else {
			proj.getLog()
					.reportTimeWarning(proj.CUSTOM_CENTROIDS_FILENAME.getValue()
															+ " exists, and currently is the proxy for LRR computation being completed");
		}
		proj.saveProperties();
		// All below stuff is just for fun...

		SampleData.createMinimalSampleData(proj);
		if (!Files.exists(proj.SAMPLE_QC_FILENAME.getValue())) {
			LrrSd.init(proj, null, null, numthreads);
		} else {
			log.reportFileExists(proj.SAMPLE_QC_FILENAME.getValue());
		}
		ArrayList<ProjectCorrected> correcteds = correctifyProject(	proj, markerTypes, offTargetsToUse,
																																correctionPCs, numthreads);// Generates
																																														// and
																																														// corrects
																																														// the
																																														// project
																																														// for
																																														// each
																																														// marker
																																														// type

		String newSampleDir = proj.PROJECT_DIRECTORY.getValue() + "samplesCorrected/";
		String newtransposedDir = proj.PROJECT_DIRECTORY.getValue() + "transposedCorrected/";

		RecompileProducer producer = new RecompileProducer(	proj, proj.getSamples(), newSampleDir,
																												proj.getMarkerSet(), correcteds);
		WorkerTrain<Hashtable<String, Float>> train =
																								new WorkerTrain<Hashtable<String, Float>>(producer,
																																													numthreads,
																																													10,
																																													proj.getLog());
		Hashtable<String, Float> recompallOutliers = new Hashtable<String, Float>();

		while (train.hasNext()) {// consolidate the pc corrected
															// projects back into a single
															// sample
			Hashtable<String, Float> tmp = train.next();
			if (tmp != null) {
				recompallOutliers.putAll(tmp);
			}
		}

		proj.SAMPLE_DIRECTORY.setValue(newSampleDir);// Final
																									// resting
																									// place
		if (recompallOutliers.size() > 0
				|| !Files.exists(proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser")) {// currently
																																												// do
																																												// to
																																												// all
																																												// the
																																												// skipping
			SerializedFiles.writeSerial(allOutliers,
																	proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser");
		}
		proj.MARKER_DATA_DIRECTORY.setValue(newtransposedDir);
		TransposeData.transposeData(proj, 2000000000, false);// already
																													// recomputed
																													// with
																													// the
																													// correction

		proj.MOSAIC_RESULTS_FILENAME.setValue(ext.addToRoot(proj.MOSAIC_RESULTS_FILENAME.getValue(),
																												".pcCorrected"));
		if (!Files.exists(proj.MOSAIC_RESULTS_FILENAME.getValue())) {
			Mosaicism.findOutliers(proj, numthreads);
		}

		generateGCModel(proj, analysisSet, referenceGenome, aName, 100);
		generateGCModel(proj, analysisSet, referenceGenome, aName, 10000);
		// generateGCModel(proj, analysisSet, referenceGenome, aName,
		// GcModel.DEFAULT_GC_MODEL_BIN_FASTA);
	}

	private static String generateGCModel(Project proj, LocusSet<Segment> analysisSet,
																				ReferenceGenome referenceGenome, ASSEMBLY_NAME aName,
																				int buffer) {
		String gcFile = ext.addToRoot(proj.GC_MODEL_FILENAME.getValue(), ".buffer_" + buffer);
		if (!Files.exists(gcFile)) {
			MarkerSet markerSet = proj.getMarkerSet();
			String[] markerNames = markerSet.getMarkerNames();

			try {
				PrintWriter writer = new PrintWriter(new FileWriter(gcFile));
				String[] header = new String[] {"Name", "Chr", "Position", "GC"};
				writer.println(Array.toStr(header));
				for (int i = 0; i < markerNames.length; i++) {
					if (i % 1000 == 0) {
						proj.getLog().reportTimeInfo("Loaded gc content for " + (i + 1) + " bins");
					}
					writer.println(markerNames[i]+ "\t" + markerSet.getChrs()[i] + "\t"
													+ markerSet.getPositions()[i] + "\t"
													+ ReferenceGenome.getPercentGC(referenceGenome.getSequenceFor(analysisSet.getLoci()[i].getBufferedSegment(buffer),
																																												aName, buffer > 100000)));
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
		private final Project proj;
		private final NGS_MARKER_TYPE type;

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

	private static ArrayList<ProjectCorrected> correctifyProject(	Project proj,
																																ArrayList<MarkerFileType> types,
																																String[] offTargetsToUse,
																																int correctionPCs, int numthreads) {
		proj.SAMPLE_CALLRATE_THRESHOLD.setValue(0.0);
		proj.LRRSD_CUTOFF.setValue(.60);
		proj.INTENSITY_PC_NUM_COMPONENTS.setValue(Math.max(20, correctionPCs));
		String mediaMarks = ext.addToRoot(proj.INTENSITY_PC_MARKERS_FILENAME.getValue(), ".median");
		ArrayList<ProjectCorrected> correctedProjects = new ArrayList<ProjectCorrected>();
		Files.writeArray(Array.subArray(proj.getMarkerNames(), 0, 1000), mediaMarks);
		String[] autoMarks = proj.getAutosomalMarkers();
		for (MarkerFileType type : types) {
			String base = "";
			if (type.getType() == null) {
				base = base + "BAM_PCS_ALL_MARKERS";
			} else {
				base = base + "BAM_PCS_" + type.getType().getFlag();
			}
			String markerfile = proj.PROJECT_DIRECTORY.getValue() + base + "autosomal_inputMarkers.txt";

			String[] tmpList;
			if (type.getType() != null && type.getType() == NGS_MARKER_TYPE.OFF_TARGET) {
				proj.getLog().reportTimeInfo("Detected "+ offTargetsToUse.length
																			+ " off target regions to use for pca of");
				tmpList = offTargetsToUse;
			} else {
				tmpList = HashVec.loadFileToStringArray(type.getFile(), true, new int[] {0}, true);
			}
			ArrayList<String> autosomalToUse = new ArrayList<String>();
			int[] indices = ext.indexLargeFactors(tmpList, autoMarks, true, proj.getLog(), false, false);
			for (int indice : indices) {
				if (indice >= 0) {
					autosomalToUse.add(autoMarks[indice]);
				}
			}
			Files.writeIterable(autosomalToUse, markerfile);
			if (!autosomalToUse.isEmpty()) {
				proj.INTENSITY_PC_MARKERS_FILENAME.setValue(markerfile);
				MitoPipeline.catAndCaboodle(proj, numthreads, mediaMarks,
																		proj.INTENSITY_PC_NUM_COMPONENTS.getValue(), base, false, true,
																		0, null, null, null, null, false, false, true, false, true,
																		false, -1, -1, GENOME_BUILD.HG19, MitoPipeline.DEFAULT_PVAL_OPTS,
																		null, false, true);

				String pcCorrectedFile = ext.addToRoot(	proj.getPropertyFilename(),
																								"."													+ proj.INTENSITY_PC_NUM_COMPONENTS.getValue()
																																						+ "_pc_corrected_"
																																						+ base);
				String newName = proj.PROJECT_NAME.getValue()+ "_"
													+ proj.INTENSITY_PC_NUM_COMPONENTS.getValue() + "_pc_corrected_" + base;
				Files.copyFileUsingFileChannels(proj.getPropertyFilename(), pcCorrectedFile, proj.getLog());
				Project pcCorrected = new Project(pcCorrectedFile, false);
				pcCorrected.PROJECT_DIRECTORY.setValue(proj.PROJECT_DIRECTORY.getValue() + newName + "/");
				pcCorrected.PROJECT_NAME.setValue(newName);
				proj.copyBasicFiles(pcCorrected, true);
				pcCorrected.SAMPLE_DIRECTORY.setValue(pcCorrected.PROJECT_DIRECTORY.getValue()
																							+ "shadowSamples/");

				String[] correctedSamps = Array.tagOn(proj.getSamples(),
																							pcCorrected.SAMPLE_DIRECTORY.getValue(),
																							Sample.SAMPLE_FILE_EXTENSION);
				if (!Files.exists("", correctedSamps) && type.getType() != null) {
					proj.getLog()
							.reportTimeInfo("PC correcting project using " + correctionPCs + " components ");

					PennCNVPrep.exportSpecialPennCNV(	proj, "correction/",
																						pcCorrected.PROJECT_DIRECTORY.getValue()
																																	+ "tmpPCCorrection/",
																						correctionPCs, null, numthreads, 1, false,
																						LS_TYPE.REGULAR, -1, true);
					// Warning currently set up for 24 threads..
					// TODO
					PennCNVPrep.exportSpecialPennCNV(	pcCorrected,
																						"correction/",
																						pcCorrected.PROJECT_DIRECTORY.getValue()
																														+ "tmpPCCorrection/",
																						correctionPCs, null, 1, 24, true, LS_TYPE.REGULAR, 5,
																						true);
				}
				pcCorrected.saveProperties();
				if (type.getType() != null) {
					correctedProjects.add(new ProjectCorrected(pcCorrected, type.getType()));
				}
			} else {
				proj.getLog().reportTimeWarning("Did not find any markers of type " + type);
			}
		}
		return correctedProjects;
	}

	private static class RecompileProducer extends AbstractProducer<Hashtable<String, Float>> {
		private final Project proj;
		private final String[] samples;
		private final String newSampleDirectory;
		private final MarkerSet markerSet;
		private final List<ProjectCorrected> correctedProjects;
		private int index;

		public RecompileProducer(	Project proj, String[] samples, String newSampleDirectory,
															MarkerSet markerSet, List<ProjectCorrected> correctedProjects) {
			super();
			this.proj = proj;
			this.samples = samples;
			this.newSampleDirectory = newSampleDirectory;
			this.markerSet = markerSet;
			this.correctedProjects = correctedProjects;
			index = 0;
		}

		@Override
		public boolean hasNext() {
			return index < samples.length;

		}

		@Override
		public Callable<Hashtable<String, Float>> next() {
			RecompileWorker current = new RecompileWorker(proj, samples[index], newSampleDirectory,
																										markerSet, correctedProjects);
			index++;
			return current;
		}
	}

	private static class RecompileWorker implements Callable<Hashtable<String, Float>> {
		private final Project proj;
		private final String sampleName;
		private final String newSampleDirectory;
		private final MarkerSet markerSet;
		private final List<ProjectCorrected> correctedProjects;

		public RecompileWorker(	Project proj, String sampleName, String newSampleDirectory,
														MarkerSet markerSet, List<ProjectCorrected> correctedProjects) {
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

	private static Hashtable<String, Float> recompileSample(Project proj, String sampleName,
																													String newSampleDirectory,
																													MarkerSet markerSet,
																													List<ProjectCorrected> correctedProjects) {
		String sampleFile = newSampleDirectory + sampleName + Sample.SAMPLE_FILE_EXTENSION;
		proj.getLog().reportTimeInfo("Sample file = " + sampleFile);
		Hashtable<String, Float> outliers = new Hashtable<String, Float>();
		if (!Files.exists(sampleName)) {
			Sample sampleOriginal = proj.getFullSampleFromRandomAccessFile(sampleName);
			int numAccountedFor = 0;
			float[] gcs = sampleOriginal.getGCs();
			float[] intensity = Array.floatArray(markerSet.getMarkerNames().length, Float.NaN);

			float[] bafs = sampleOriginal.getBAFs();// preserve the bafs
			float[] lrrs = Array.floatArray(markerSet.getMarkerNames().length, Float.NaN);

			String[] markerNames = markerSet.getMarkerNames();
			for (ProjectCorrected corrected : correctedProjects) {
				Sample typeCorrected = corrected.getProj().getFullSampleFromRandomAccessFile(sampleName);
				for (int i = 0; i < markerNames.length; i++) {
					NGS_MARKER_TYPE markerType = NGS_MARKER_TYPE.getType(markerNames[i]);
					if (markerType == corrected.getType()) {
						intensity[i] = typeCorrected.getXs()[i];
						lrrs[i] = typeCorrected.getLRRs()[i];
						numAccountedFor++;
					}
				}
			}
			Sample sampleCorrected = new Sample(sampleName, markerSet.getFingerprint(), gcs, intensity,
																					intensity, bafs, lrrs,
																					sampleOriginal.getForwardGenotypes(),
																					sampleOriginal.getAB_Genotypes(),
																					sampleOriginal.getCanXYBeNegative());
			sampleCorrected.saveToRandomAccessFile(sampleFile, outliers, sampleName);
			if (numAccountedFor != markerNames.length) {
				throw new IllegalArgumentException("Not all markers accounted for in corrections");
			}
		} else {
			proj.getLog().reportFileExists(sampleFile);
		}

		return outliers;
	}

	private static String[] dumpLikelyOffTargetProblems(Project proj, ASSAY_TYPE aType) {

		MarkerSet markerSet = proj.getMarkerSet();
		String problemFile = ext.addToRoot(	proj.MARKER_POSITION_FILENAME.getValue(),
																				".likelyOffTargetProblems");
		String noproblemFile = ext.addToRoot(	proj.MARKER_POSITION_FILENAME.getValue(),
																					".withoutlikelyOffTargetProblems");
		String allFile = ext.addToRoot(	proj.MARKER_POSITION_FILENAME.getValue(),
																		".OffTargetProblemsFlagged");

		ArrayList<String> problems = new ArrayList<String>();
		ArrayList<String> noProblems = new ArrayList<String>();
		ArrayList<String> all = new ArrayList<String>();
		ArrayList<String> goodOffTargets = new ArrayList<String>();

		String header =
									"BinName\tCLASS=MARKER_COLOR;OFF_TARGET_OK=Blue;LIKELY_OFF_TARGET_PROBLEM=RED;OTHER_TYPE=Green";
		problems.add(header);
		noProblems.add(header);
		all.add(header);
		int[][] indices = markerSet.getIndicesByChr();
		String[] names = markerSet.getMarkerNames();
		for (int[] indice : indices) {
			for (int j = 0; j < indice.length; j++) {
				int compLeft = Math.max(j - 1, 0);
				int compRight = Math.min(j + 1, indice.length - 1);
				NGS_MARKER_TYPE current = NGS_MARKER_TYPE.getType(names[indice[j]]);
				if (current == NGS_MARKER_TYPE.OFF_TARGET) {
					NGS_MARKER_TYPE left = NGS_MARKER_TYPE.getType(names[indice[compLeft]]);
					NGS_MARKER_TYPE right = NGS_MARKER_TYPE.getType(names[indice[compRight]]);
					if (aType != ASSAY_TYPE.WGS
							&& ((compLeft != j && left != NGS_MARKER_TYPE.OFF_TARGET)
									|| (compRight != j && right != NGS_MARKER_TYPE.OFF_TARGET))) {
						goodOffTargets.add(names[indice[j]]);
						problems.add(names[indice[j]] + "\tLIKELY_OFF_TARGET_PROBLEM");
						all.add(names[indice[j]] + "\tLIKELY_OFF_TARGET_PROBLEM");
					} else {
						noProblems.add(names[indice[j]] + "\tOFF_TARGET_OK");
						all.add(names[indice[j]] + "\tOFF_TARGET_OK");
					}
				} else {
					noProblems.add(names[indice[j]] + "\tOTHER_TYPE");
					all.add(names[indice[j]] + "\tOTHER_TYPE");
				}
			}
		}
		proj.getLog()
				.reportTimeInfo("Dumping "+ (problems.size() - 1)
												+ " off target markers that will likely be biased to " + problemFile);
		Files.writeIterable(problems, problemFile);
		Files.writeIterable(noProblems, noproblemFile);
		Files.writeIterable(all, allFile);
		return Array.toStringArray(goodOffTargets);

	}

	private static ArrayList<MarkerFileType> generateMarkerPositions(	Project proj,
																																		LocusSet<BEDFeatureSeg> bLocusSet,
																																		LocusSet<Segment> genomeBinsMinusBinsCaputure,
																																		LocusSet<VariantSeg> varFeatures) {
		String positions = proj.MARKER_POSITION_FILENAME.getValue();
		proj.getLog().reportTimeInfo("Postions will be set to the midpoint of each segment");
		String[] markerNames = new String[bLocusSet.getLoci().length
																			+ genomeBinsMinusBinsCaputure.getLoci().length
																			+ varFeatures.getLoci().length];
		String header =
									"BinName\tChr\tPosition\tCLASS=MARKER_COLOR;OFF_TARGET=Blue;VARIANT_SITE=RED;ON_TARGET=Green";
		ArrayList<String> onTMarkers = new ArrayList<String>();
		onTMarkers.add(header);

		ArrayList<String> offTMarkers = new ArrayList<String>();
		offTMarkers.add(header);

		ArrayList<String> variantSiteMarkers = new ArrayList<String>();
		variantSiteMarkers.add(header);

		ArrayList<String> allMarkerColors = new ArrayList<String>();
		allMarkerColors.add(header);

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
				String out = markerName+ "\t" + bFeatureSeg.getChr() + "\t" + pos + "\t"
											+ NGS_MARKER_TYPE.ON_TARGET.getFlag();
				onTMarkers.add(out);
				allMarkerColors.add(out);
				writer.println(out);
				markerIndex++;
			}
			for (int i = 0; i < genomeBinsMinusBinsCaputure.getLoci().length; i++) {
				Segment binnedSeg = genomeBinsMinusBinsCaputure.getLoci()[i];
				String markerName =
													binnedSeg.getUCSClocation() + "|" + NGS_MARKER_TYPE.OFF_TARGET.getFlag();

				markerNames[markerIndex] = markerName;
				int diff = binnedSeg.getStop() - binnedSeg.getStart();
				int mid = Math.round((float) diff / 2);
				int pos = binnedSeg.getStart() + mid;
				String out = markerName+ "\t" + binnedSeg.getChr() + "\t" + pos + "\t"
											+ NGS_MARKER_TYPE.OFF_TARGET.getFlag();
				offTMarkers.add(out);
				allMarkerColors.add(out);
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
				String out = markerName+ "\t" + variantFeatureSeg.getChr() + "\t" + pos + "\t"
											+ NGS_MARKER_TYPE.VARIANT_SITE.getFlag();
				variantSiteMarkers.add(out);
				allMarkerColors.add(out);
				writer.println(out);
				markerIndex++;
			}

			writer.close();
		} catch (Exception e) {
			proj.getLog().reportError("Error writing to " + positions);
			proj.getLog().reportException(e);
		}
		String onTargetFile = ext.addToRoot(proj.MARKER_POSITION_FILENAME.getValue(),
																				"." + NGS_MARKER_TYPE.ON_TARGET.getFlag());
		String offTargetFile = ext.addToRoot(	proj.MARKER_POSITION_FILENAME.getValue(),
																					"." + NGS_MARKER_TYPE.OFF_TARGET.getFlag());
		String variantSiteTargetFile = ext.addToRoot(	proj.MARKER_POSITION_FILENAME.getValue(),
																									"." + NGS_MARKER_TYPE.VARIANT_SITE.getFlag());

		String allMarkerColorFile = ext.addToRoot(proj.MARKER_POSITION_FILENAME.getValue(),
																							".allMarkers.color");
		proj.MARKER_COLOR_KEY_FILENAMES.setValue(new String[] {	onTargetFile, offTargetFile,
																														variantSiteTargetFile,
																														allMarkerColorFile});
		String allMarkerFile = ext.addToRoot(proj.MARKER_POSITION_FILENAME.getValue(), ".allMarkers");


		ArrayList<MarkerFileType> markerTypes = new ArrayList<MarkerFileType>();

		if (!onTMarkers.isEmpty()) {
			Files.writeArray(Array.toStringArray(onTMarkers), onTargetFile);
			markerTypes.add(new MarkerFileType(NGS_MARKER_TYPE.ON_TARGET, onTargetFile));
		} else {
			proj.getLog()
					.reportTimeWarning("No " + NGS_MARKER_TYPE.ON_TARGET.getFlag() + " markers detected");
		}
		if (!offTMarkers.isEmpty()) {
			Files.writeArray(Array.toStringArray(offTMarkers), offTargetFile);
			markerTypes.add(new MarkerFileType(NGS_MARKER_TYPE.OFF_TARGET, offTargetFile));
		} else {
			proj.getLog()
					.reportTimeWarning("No " + NGS_MARKER_TYPE.OFF_TARGET.getFlag() + " markers detected");
		}

		if (!variantSiteMarkers.isEmpty()) {
			Files.writeArray(Array.toStringArray(variantSiteMarkers), variantSiteTargetFile);
			markerTypes.add(new MarkerFileType(NGS_MARKER_TYPE.VARIANT_SITE, variantSiteTargetFile));
		} else {
			proj.getLog()
					.reportTimeWarning("No " + NGS_MARKER_TYPE.VARIANT_SITE.getFlag() + " markers detected");
		}

		Files.writeArray(Array.toStringArray(allMarkerColors), allMarkerColorFile);

		Files.writeArray(	HashVec.loadFileToStringArray(proj.MARKER_POSITION_FILENAME.getValue(), true,
																										new int[] {0}, true),
											allMarkerFile);
		markerTypes.add(new MarkerFileType(null, allMarkerFile));

		if (!Files.exists(proj.MARKERSET_FILENAME.getValue(true, true))) {
			Markers.orderMarkers(	markerNames, proj.MARKER_POSITION_FILENAME.getValue(),
														proj.MARKERSET_FILENAME.getValue(true, true), proj.getLog());
		} else {
			proj.getLog().reportFileExists(proj.MARKERSET_FILENAME.getValue(true, true));
		}
		return markerTypes;
	}

	private static class MarkerFileType {
		private final NGS_MARKER_TYPE type;
		private final String file;

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
		int captureBuffer = CAPTURE_BUFFER;
		String vcf = null;
		int correctionPCs = 4;


		String usage = "\n" + "seq.manage.BamImport requires 0-1 arguments\n";
		usage += "(1) filename (i.e. proj= ( nodefault))\n" + "";
		usage += "(2) bed file to import  (i.e. importBed=" + binBed + " ( no default))\n" + "";
		usage += Ext.getNumThreadsCommand(3, numthreads);
		usage += "(4) bed file to import  (i.e. captureBed=" + captureBed + " ( no default))\n" + "";
		usage +=
					"(5) a vcf, if provided the variants will be imported with bp resolution  (i.e. vcf= ( no default))\n"
							+ "";
		usage += "(6) number of PCs to correct with  (i.e. correctionPCs= ( no default))\n" + "";


		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("proj=")) {
				filename = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("importBed=")) {
				binBed = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("captureBed=")) {
				captureBed = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("vcf=")) {
				vcf = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith(Ext.NUM_THREADS_COMMAND)) {
				numthreads = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("correctionPCs")) {
				correctionPCs = ext.parseIntArg(arg);
				numArgs--;
			}

			else {
				System.err.println("Error - invalid argument: " + arg);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			Project proj = new Project(filename, false);
			importTheWholeBamProject(	proj, binBed, captureBed, vcf, captureBuffer, correctionPCs, true,
																ASSAY_TYPE.WXS, ASSEMBLY_NAME.HG19, null, numthreads);
		} catch (Exception e) {
			new Logger().reportException(e);
		}
	}
}
