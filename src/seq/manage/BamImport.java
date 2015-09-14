package seq.manage;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Hashtable;

import seq.manage.BEDFileReader.BEDFeatureSeg;
import seq.manage.BamSegPileUp.PileupProducer;
import seq.qc.FilterNGS;
import common.Array;
import common.Files;
import common.Logger;
import common.WorkerTrain;
import common.ext;
import common.PSF.Ext;
import cnv.analysis.CentroidCompute;
import cnv.analysis.CentroidCompute.Builder;
import cnv.filesys.Centroids;
import cnv.filesys.Project;
import cnv.filesys.Project.ARRAY;
import cnv.manage.Markers;
import cnv.manage.MitoPipeline;
import cnv.manage.TransposeData;
import cnv.var.LocusSet;

public class BamImport {

	public static void importTheWholeBamProject(Project proj, String binBed, int numthreads) {
		if (proj.getArrayType() == ARRAY.NGS) {
			Logger log = proj.getLog();
			String serDir = proj.PROJECT_DIRECTORY.getValue() + "tmpBamSer/";
			String[] bamsToImport = Files.listFullPaths(proj.SOURCE_DIRECTORY.getValue(), proj.SOURCE_FILENAME_EXTENSION.getValue(), false);
			log.reportTimeInfo("Found " + bamsToImport.length + " bam files to import");
			if (BedOps.verifyBedIndex(binBed, log)) {
				BEDFileReader reader = new BEDFileReader(binBed, true);
				LocusSet<BEDFeatureSeg> bLocusSet = reader.loadAll(log);
				reader.close();
				if (!bLocusSet.hasNoOverlap()) {
					generateMarkerPositions(proj, bLocusSet);
					ReferenceGenome referenceGenome = new ReferenceGenome(proj.REFERENCE_GENOME_FASTA_FILENAME.getValue(), log);
					String projectSeqSer = proj.PROJECT_DIRECTORY.getValue() + ext.removeDirectoryInfo(binBed) + "seq.ser";
					String[][] refForIntervals = null;
					if (!Files.exists(projectSeqSer)) {
						refForIntervals = referenceGenome.getSequencesFor(bLocusSet.getLoci(), 1000);
						Files.writeSerial(refForIntervals, projectSeqSer, true);
						log.reportTimeInfo("Writing serialized bin sequences to " + projectSeqSer);

					} else {
						log.reportTimeInfo("Loading bin sequences from " + projectSeqSer);
						refForIntervals = (String[][]) Files.readSerial(projectSeqSer, false, log, false, true);
					}
					if (refForIntervals.length != bLocusSet.getLoci().length) {
						String error = "Have " + refForIntervals.length + " sequence intervals for " + bLocusSet.getLoci().length + " bins";
						log.reportTimeError(error);
						return;
					}

					long fingerPrint = proj.getMarkerSet().getFingerprint();

					log.reportTimeInfo(bLocusSet.getLoci().length + " segments to pile");
					reader.close();
					FilterNGS filterNGS = new FilterNGS(20, 20, null);

					PileupProducer producer = new PileupProducer(bamsToImport, serDir, refForIntervals, filterNGS, bLocusSet.getStrictSegments(), log);
					WorkerTrain<BamPile[]> train = new WorkerTrain<BamPile[]>(producer, numthreads, 2, log);
					int index = 0;
					Hashtable<String, Float> allOutliers = new Hashtable<String, Float>();
					proj.SAMPLE_DIRECTORY.getValue(true, false);
					while (train.hasNext()) {
						BamSample bamSample = new BamSample(proj, bamsToImport[index], train.next());
						Hashtable<String, Float> outliers = bamSample.writeSample(fingerPrint);
						if (outliers.size() > 0) {
							allOutliers.putAll(allOutliers);
						}
						index++;
					}
					Files.writeSerial(allOutliers, proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser");
					TransposeData.transposeData(proj, 2000000000, false);
					CentroidCompute.computeAndDumpCentroids(proj, null, proj.CUSTOM_CENTROIDS_FILENAME.getValue(), new Builder(), numthreads, 2);
					Centroids.recompute(proj, proj.CUSTOM_CENTROIDS_FILENAME.getValue(), true);
					TransposeData.transposeData(proj, 2000000000, false);

					generatePCFile(proj, numthreads);

				} else {
					log.reportTimeError("The bed file " + binBed + " had overlapping segments, currently non -overlapping segments are required");
				}
			}
		} else {
			proj.getLog().reportTimeError(proj.ARRAY_TYPE.getName() + " must be set to " + ARRAY.NGS);
		}
	}

	private static void generatePCFile(Project proj, int numthreads) {
		Files.writeList(proj.getMarkerNames(), proj.TARGET_MARKERS_FILENAME.getValue());
		String mediaMarks = ext.addToRoot(proj.TARGET_MARKERS_FILENAME.getValue(), ".median");
		Files.writeList(Array.subArray(proj.getMarkerNames(), 0, 10), mediaMarks);
		String base = ext.rootOf(proj.INTENSITY_PC_FILENAME.getValue());
		MitoPipeline.catAndCaboodle(proj, numthreads, "0", mediaMarks, proj.getSamples().length - 1, base, false, false, 0, null, null, null, false, false, false, true);
	}

	private static void generateMarkerPositions(Project proj, LocusSet<BEDFeatureSeg> bLocusSet) {
		String positions = proj.MARKER_POSITION_FILENAME.getValue();
		String[] markerNames = new String[bLocusSet.getLoci().length];
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(positions));
			writer.println("BinName\tChr\tPosition");
			for (int i = 0; i < bLocusSet.getLoci().length; i++) {
				BEDFeatureSeg bFeatureSeg = bLocusSet.getLoci()[i];
				String markerName = bFeatureSeg.getUCSClocation();
				String name = bFeatureSeg.getBedFeature().getName();
				if (name != null) {
					markerName += "|" + name;
				}
				markerNames[i] = markerName;
				writer.println(markerName + "\t" + bFeatureSeg.getChr() + "\t" + bFeatureSeg.getStart());
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
		int numthreads = 24;
		// String referenceGenomeFasta = "hg19_canonical.fa";
		String logfile = null;
		Logger log;

		String usage = "\n" + "seq.manage.BamImport requires 0-1 arguments\n";
		usage += "(1) filename (i.e. proj= ( nodefault))\n" + "";
		usage += "(2) bed file to import  (i.e. bed=" + binBed + " ( no default))\n" + "";
		usage += Ext.getNumThreadsCommand(3, numthreads);
		// usage += "(3) reference genome  (i.e. ref=" + referenceGenomeFasta + " ( default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("bed=")) {
				binBed = args[i].split("=")[1];
				numArgs--;
			}
			// else if (args[i].startsWith("ref=")) {
			// referenceGenomeFasta = args[i].split("=")[1];
			// numArgs--;
			// }
			else if (args[i].startsWith(Ext.NUM_THREADS_COMMAND)) {
				numthreads = ext.parseIntArg(args[i]);
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
			Project proj = new Project(filename, false);
			importTheWholeBamProject(proj, binBed, numthreads);
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
