package seq.manage;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
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
import common.WorkerTrain;
import common.WorkerTrain.Producer;
import common.ext;
import common.PSF.Ext;
import cnv.analysis.CentroidCompute;
import cnv.analysis.CentroidCompute.CentroidBuilder;
import cnv.filesys.Centroids;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.filesys.Project.ARRAY;
import cnv.manage.Markers;
import cnv.manage.MitoPipeline;
import cnv.manage.TransposeData;
import cnv.qc.GcAdjustor;
import cnv.qc.GcAdjustorParameter;
import cnv.qc.GcAdjustor.GCAdjustorBuilder;
import cnv.var.LocusSet;
import filesys.Segment;

public class BamImport {
	public static final String OFF_TARGET_FLAG = "OFF_TARGET";
	public static final String VARIANT_SITE_FLAG = "VARIANT_SITE";

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
			BamSample bamSample = new BamSample(proj, result.getBam(), result.loadResults(log));
			sample = bamSample.getSampleName();
			bamIndexStats = BamOps.getBamIndexStats(result.getBam());
			outliers = bamSample.writeSample(fingerPrint);
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

	private static LocusSet<BEDFeatureSeg> extractVCF(Project proj, String vcf) {
		LocusSet<BEDFeatureSeg> varLocusSet = new LocusSet<BEDFeatureSeg>(new BEDFeatureSeg[] {}, true, proj.getLog()) {

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
						writer.println(vc.getContig() + "\t" + vc.getStart() + "\t" + vc.getEnd() + "\t" + name);
					}
					reader.close();

					writer.close();
				} catch (Exception e) {
					proj.getLog().reportError("Error writing to " + out);
					proj.getLog().reportException(e);
				}

			}

			BEDFileReader readerVarBed = new BEDFileReader(out, false);
			varLocusSet = readerVarBed.loadAll(proj.getLog());
			readerVarBed.close();
		}
		return varLocusSet;

	}

	public static void importTheWholeBamProject(Project proj, String binBed, String captureBed, String optionalVCF, int captureBuffer, int numthreads) {
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
					LocusSet<Segment> genomeBinsMinusBinsCaputure = referenceGenome.getBins(20000).removeThese(LocusSet.combine(bLocusSet, readerCapture.loadAll(log), true, log).mergeOverlapping(true), 4000);
					log.reportTimeInfo(genomeBinsMinusBinsCaputure.getBpCovered() + " bp covered by reference bins int the anti-on-target regions");
					log.memoryFree();

					LocusSet<BEDFeatureSeg> varFeatures = extractVCF(proj, optionalVCF);
					if(varFeatures.getLoci().length>0){
					log.reportTimeInfo(varFeatures.getBpCovered() + " bp covered by known variant sites");
					}

					generateMarkerPositions(proj, bLocusSet, genomeBinsMinusBinsCaputure, varFeatures);
					log.memoryFree();
					LocusSet<Segment> analysisSet = LocusSet.combine(bLocusSet.getStrictSegmentSet(), genomeBinsMinusBinsCaputure, true, log);
					analysisSet = LocusSet.combine(analysisSet, varFeatures.getStrictSegmentSet(), true, log);

					log.memoryFree();

					log.reportTimeInfo(analysisSet.getLoci().length + " segments to pile");
					System.exit(1);

					generateGCModel(proj, analysisSet, referenceGenome);

					FilterNGS filterNGS = new FilterNGS(20, 20, null);
					PileupProducer pileProducer = new PileupProducer(bamsToImport, serDir, referenceGenome.getReferenceFasta(), filterNGS, analysisSet.getStrictSegments(), log);
					WorkerTrain<BamPileResult> pileTrain = new WorkerTrain<BamPileResult>(pileProducer, numthreads, 2, log);
					int index = 0;
					proj.SAMPLE_DIRECTORY.getValue(true, false);
					proj.XY_SCALE_FACTOR.setValue((double) 10);

					BamPileResult[] results = new BamPileResult[bamsToImport.length];
					while (pileTrain.hasNext()) {
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
					while (conversionTrain.hasNext()) {
						BamPileConversionResults conversionResult = conversionTrain.next();
						BamIndexStats bamIndexStats = conversionResult.getBamIndexStats();
						int numAligned = bamIndexStats.getAlignedRecordCount();
						int numNotAligned = bamIndexStats.getUnalignedRecordCount();
						mappedReadCounts[convIndex + 1] = conversionResult.getSample() + "\t" + numAligned + "\t" + numNotAligned;
						convIndex++;
						if (conversionResult.getOutliers().size() > 0) {
							allOutliers.putAll(conversionResult.getOutliers());
						}

					}
					String readCountFile = proj.PROJECT_DIRECTORY.getValue() + "sample.readCounts.txt";

					Files.writeList(mappedReadCounts, readCountFile);

					Files.writeSerial(allOutliers, proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser");

					TransposeData.transposeData(proj, 2000000000, false);
					CentroidCompute.computeAndDumpCentroids(proj, proj.CUSTOM_CENTROIDS_FILENAME.getValue(), new CentroidBuilder(), numthreads, 2);
					Centroids.recompute(proj, proj.CUSTOM_CENTROIDS_FILENAME.getValue(), true, numthreads);
					TransposeData.transposeData(proj, 2000000000, false);

					GCAdjustorBuilder gAdjustorBuilder = new GCAdjustorBuilder();
					GcAdjustorParameter.generate(proj, "GC_ADJUSTMENT/", proj.REFERENCE_GENOME_FASTA_FILENAME.getValue(), gAdjustorBuilder, false, GcAdjustor.GcModel.DEFAULT_GC_MODEL_BIN_FASTA, numthreads);
					// GcAdjustorParameters params =
					generatePCFile(proj, numthreads);
					proj.INTENSITY_PC_NUM_COMPONENTS.setValue(5);
					proj.saveProperties();

					// String PCCorrected = ext.addToRoot(proj.getPropertyFilename(), "." + proj.INTENSITY_PC_NUM_COMPONENTS.getValue() + "_pc_corrected");
					// String newName = proj.PROJECT_NAME.getValue() + "_" + proj.INTENSITY_PC_NUM_COMPONENTS.getValue() + "_pc_corrected";
					// Files.copyFileUsingFileChannels(proj.getPropertyFilename(), PCCorrected, log);
					// Project pcCorrected = new Project(PCCorrected, false);
					// pcCorrected.PROJECT_DIRECTORY.setValue(proj.PROJECT_DIRECTORY.getValue() + newName + "/");
					// pcCorrected.PROJECT_NAME.setValue(newName);
					// proj.copyBasicFiles(pcCorrected, true);
					//
					// log.reportTimeInfo("PC correcting project using " + proj.INTENSITY_PC_NUM_COMPONENTS.getValue() + " components ");
					// PennCNVPrep.exportSpecialPennCNV(proj, "correction/", pcCorrected.PROJECT_DIRECTORY.getValue() + "tmpPCCorrection/", proj.INTENSITY_PC_NUM_COMPONENTS.getValue(), null, 1, numthreads, false, false, false, -1, true);
					// PennCNVPrep.exportSpecialPennCNV(pcCorrected, "correction/", pcCorrected.PROJECT_DIRECTORY.getValue() + "tmpPCCorrection/", proj.INTENSITY_PC_NUM_COMPONENTS.getValue(), null, 1, numthreads, false, true, false, -1, true);
					// pcCorrected.SAMPLE_DIRECTORY.setValue(pcCorrected.PROJECT_DIRECTORY.getValue() + "shadowSamples3/");
					// pcCorrected.saveProperties();
					// TransposeData.transposeData(pcCorrected, 2000000000, false);

					// MDL

					// Project proj =
				} else {
					log.reportTimeError("The bed file " + binBed + " had overlapping segments, currently non -overlapping segments are required");
				}
			}
		} else {
			proj.getLog().reportTimeError(proj.ARRAY_TYPE.getName() + " must be set to " + ARRAY.NGS);
		}
	}

	private static void generateGCModel(Project proj, LocusSet<Segment> analysisSet, ReferenceGenome referenceGenome) {
		String gcFile = proj.GC_MODEL_FILENAME.getValue();
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
					writer.println(markerNames[i] + "\t" + markerSet.getChrs()[i] + "\t" + markerSet.getPositions()[i] + "\t" + ReferenceGenome.getPercent(referenceGenome.getSequenceFor(analysisSet.getLoci()[i])));
				}
				writer.close();
			} catch (Exception e) {
				proj.getLog().reportError("Error writing to " + gcFile);
				proj.getLog().reportException(e);
			}
		} else {
			proj.getLog().reportFileExists(gcFile);
		}
	}

	private static void generatePCFile(Project proj, int numthreads) {
		Files.writeList(proj.getMarkerNames(), proj.TARGET_MARKERS_FILENAMES.getValue()[0]);
		String mediaMarks = ext.addToRoot(proj.TARGET_MARKERS_FILENAMES.getValue()[0], ".median");
		Files.writeList(Array.subArray(proj.getMarkerNames(), 0, 1000), mediaMarks);
		String base = "BAM_PCS";
		MitoPipeline.catAndCaboodle(proj, numthreads, "0", mediaMarks, proj.getSamples().length - 1, base, false, false, 0, null, null, null, false, false, false, true, false, null, -1, -1);
	}

	private static void generateMarkerPositions(Project proj, LocusSet<BEDFeatureSeg> bLocusSet, LocusSet<Segment> genomeBinsMinusBinsCaputure, LocusSet<BEDFeatureSeg> varFeatures) {
		String positions = proj.MARKER_POSITION_FILENAME.getValue();
		proj.getLog().reportTimeInfo("Postions will be set to the midpoint of each segment");
		String[] markerNames = new String[bLocusSet.getLoci().length + genomeBinsMinusBinsCaputure.getLoci().length + varFeatures.getLoci().length];
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(positions));
			int markerIndex = 0;
			writer.println("BinName\tChr\tPosition");
			for (int i = 0; i < bLocusSet.getLoci().length; i++) {
				BEDFeatureSeg bFeatureSeg = bLocusSet.getLoci()[i];
				String markerName = bFeatureSeg.getUCSClocation();
				String name = bFeatureSeg.getBedFeature().getName();
				if (name != null) {
					markerName += "|" + name;
				}
				markerNames[markerIndex] = markerName;
				int diff = bFeatureSeg.getStop() - bFeatureSeg.getStart();
				int mid = Math.round((float) diff / 2);
				int pos = bFeatureSeg.getStart() + mid;
				writer.println(markerName + "\t" + bFeatureSeg.getChr() + "\t" + pos);
				markerIndex++;
			}
			for (int i = 0; i < genomeBinsMinusBinsCaputure.getLoci().length; i++) {
				Segment binnedSeg = genomeBinsMinusBinsCaputure.getLoci()[i];
				String markerName = binnedSeg.getUCSClocation() + "|" + OFF_TARGET_FLAG;

				markerNames[markerIndex] = markerName;
				int diff = binnedSeg.getStop() - binnedSeg.getStart();
				int mid = Math.round((float) diff / 2);
				int pos = binnedSeg.getStart() + mid;
				writer.println(markerName + "\t" + binnedSeg.getChr() + "\t" + pos);
				markerIndex++;
			}

			for (int i = 0; i < varFeatures.getLoci().length; i++) {
				BEDFeatureSeg variantFeatureSeg = varFeatures.getLoci()[i];
				String markerName = variantFeatureSeg.getUCSClocation();
				String name = variantFeatureSeg.getBedFeature().getName();
				if (name != null) {
					markerName += "|" + name;
				}
				markerName = markerName + "|" + VARIANT_SITE_FLAG;
				markerNames[markerIndex] = markerName;
				int pos = variantFeatureSeg.getStart();
				writer.println(markerName + "\t" + variantFeatureSeg.getChr() + "\t" + pos);
				markerIndex++;
			}

			writer.close();
		} catch (Exception e) {
			proj.getLog().reportError("Error writing to " + positions);
			proj.getLog().reportException(e);
		}

		Markers.orderMarkers(markerNames, proj.MARKER_POSITION_FILENAME.getValue(), proj.MARKERSET_FILENAME.getValue(true, true), proj.getLog());
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		String binBed = "binsToImport.bed";
		String captureBed = "AgilentCaptureRegions.txt";
		int numthreads = 24;
		int captureBuffer = 400;
		String vcf = null;
		// String referenceGenomeFasta = "hg19_canonical.fa";
		// String logfile = null;
		// Logger log;

		String usage = "\n" + "seq.manage.BamImport requires 0-1 arguments\n";
		usage += "(1) filename (i.e. proj= ( nodefault))\n" + "";
		usage += "(2) bed file to import  (i.e. importBed=" + binBed + " ( no default))\n" + "";
		usage += Ext.getNumThreadsCommand(3, numthreads);
		usage += "(4) bed file to import  (i.e. captureBed=" + captureBed + " ( no default))\n" + "";
		usage += "(5) a vcf, if provided the variants will be imported with bp resolution  (i.e. vcf= ( no default))\n" + "";

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
			importTheWholeBamProject(proj, binBed, captureBed, vcf, captureBuffer, numthreads);
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
