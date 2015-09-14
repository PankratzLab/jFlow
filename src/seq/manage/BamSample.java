package seq.manage;

import java.util.Hashtable;

import seq.manage.BamOps.BamIndexStats;
import stats.Maths;
import cnv.analysis.BeastScore;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.filesys.Sample;
import common.Array;
import filesys.Segment;

public class BamSample {
	private static final double MAX_MAPQ = 60;
	private static final double SCALE_FACTOR_NUM_READS = 1000000;
	private String bamFile;
	private String sampleName;
	private BamPile[] bamPiles;
	private double[] rawDepth;
	private double[] normDepth;
	private double[] mapQs;
	private double[] percentWithMismatch;
	private Project proj;

	public BamSample(Project proj, String bamFile, BamPile[] bamPiles) {
		super();
		this.proj = proj;
		this.bamFile = bamFile;
		this.sampleName = BamOps.getSampleName(bamFile);
		this.bamPiles = bamPiles;
		init();
	}

	private static double computeRPKM(int numMappedReads, Segment bin, int numTotalMappedReads) {
		double data = (double) numMappedReads / bin.getSize();
		double scale = numTotalMappedReads > 0 ? (double) SCALE_FACTOR_NUM_READS / numTotalMappedReads : 0;
		return data * scale;
	}

	private void init() {
		MarkerSet markerSet = proj.getMarkerSet();
		BamIndexStats bamIndexStats = BamOps.getBamIndexStats(bamFile);
		System.out.println(bamIndexStats.getAlignedRecordCount());
		// double scaleFactor = (double) bamIndexStats.getAlignedRecordCount()
		this.rawDepth = new double[bamPiles.length];
		this.mapQs = new double[bamPiles.length];
		this.percentWithMismatch = new double[bamPiles.length];
		byte currentChr = 0;
		int currentPos = 0;
		for (int i = 0; i < bamPiles.length; i++) {
			if (currentChr > bamPiles[i].getBin().getChr() || (bamPiles[i].getBin().getChr() <= currentChr && currentPos > bamPiles[i].getBin().getStart())) {
				String error = "BUG, segments are unsorted";
				error += " current CHR: " + currentChr + " CurrentPOS: " + currentPos;
				error += "Bin: " + bamPiles[i].getBin().getUCSClocation();
				error += "Index " + i;
				proj.getLog().reportTimeError(error);
				throw new IllegalStateException(error);
			} else {
				currentChr = bamPiles[i].getBin().getChr();
				currentPos = bamPiles[i].getBin().getStart();
			}
			// rawDepth[i] = bamPiles[i].getOverallAvgDepth();
			// System.out.println(bamPiles[i].getNumOverlappingReads());
			// rawDepth[i] = 10 + computeRPKM(bamPiles[i].getNumOverlappingReads(), bamPiles[i].getBin(), bamIndexStats.getAlignedRecordCount());
			rawDepth[i] = computeRPKM(bamPiles[i].getNumOverlappingReads(), bamPiles[i].getBin(), bamIndexStats.getAlignedRecordCount());
			// rawDepth[i] = Maths.log2(rawDepth[i]);
			mapQs[i] = Math.min(bamPiles[i].getOverallAvgMapQ() / MAX_MAPQ, 1);
			int currentSize = bamPiles[i].getBin().getSize();
			double normBasesOverlap = (double) bamPiles[i].getNumBasesOverlap() / currentSize;
			double normBasesMiss = (double) bamPiles[i].getNumBasesWithMismatch() / currentSize;
			double percentMiss = 0;
			if (normBasesOverlap > 0) {
				percentMiss = normBasesMiss / normBasesOverlap;
			}
			percentWithMismatch[i] = percentMiss;
		}

		BeastScore beastScore = new BeastScore(Array.toFloatArray(rawDepth), markerSet.getIndicesByChr(), null, proj.getLog());
		float[] scaleMAD = beastScore.getScaleMadRawData();
		int[][] chrIndices = markerSet.getIndicesByChr();
		for (int i = 0; i < chrIndices.length; i++) {
			boolean error = false;
			for (int j = 0; j < chrIndices[i].length; j++) {
				if (Double.isNaN(scaleMAD[chrIndices[i][j]]) && rawDepth[chrIndices[i][j]] == 0) {// should only happen if the MAD is NaN
					if (!error) {
						String warning = "Found invalid scale MAD depth for " + bamFile + ", bin " + markerSet.getMarkerNames()[chrIndices[i][j]];
						warning += "Setting all of chr" + i + " to 0";
						proj.getLog().reportTimeWarning(warning);
						error = true;
					}
				}
			}
			if (error) {
				for (int j = 0; j < chrIndices[i].length; j++) {
					scaleMAD[chrIndices[i][j]] = 0;
				}

			}
		}

		for (int j = 0; j < scaleMAD.length; j++) {
			if (Double.isNaN(scaleMAD[j])) {
				String error = "Found invalid scale MAD depth for " + bamFile + ", bin " + markerSet.getMarkerNames()[j];
				proj.getLog().reportTimeError(error);
				// throw new IllegalStateException(error);
			}
		}

		this.normDepth = Array.scaleMinTo(Array.toDoubleArray(scaleMAD), 1);
		double scale = proj.XY_SCALE_FACTOR.getValue();
		normDepth = Array.multiply(normDepth, (double) 1 / scale);
		
		
		percentWithMismatch = Array.scale(percentWithMismatch);
		System.out.println(Array.min(percentWithMismatch));

		//
		// // // if(Array.min(normDepth)<0){
		// // // System.err.println("less dan 0");
		// // // System.exit(1);
		// // // }
		// double maxNorm = Array.max(normDepth);
		// double minNorm = Array.min(normDepth);
		// for (int i = 0; i < normDepth.length; i++) {
		// normDepth[i] = (normDepth[i] - minNorm) / (maxNorm - minNorm);
		// normDepth[i] += 1;
		// // normDepth[i] *= 10;
		// }
//		System.out.println(Array.mean(normDepth));
//		System.out.println(Array.min(normDepth));
//		System.out.println(Array.max(normDepth));
		for (int j = 0; j < normDepth.length; j++) {
			if (Double.isNaN(normDepth[j])) {
				String error = "Found invalid normalized depth for " + bamFile + ", bin " + markerSet.getMarkerNames()[j];
				proj.getLog().reportTimeError(error);
				throw new IllegalStateException(error);
			}
		}

	}

	public Hashtable<String, Float> writeSample(long fingerprint) {
		Hashtable<String, Float> outliers = new Hashtable<String, Float>();
		byte[] genos = Array.byteArray(bamPiles.length, (byte) 1);
		float[] blankLRRs = Array.floatArray(bamPiles.length, 0);
		String sampleFile = proj.SAMPLE_DIRECTORY.getValue() + sampleName + Sample.SAMPLE_DATA_FILE_EXTENSION;
		Sample sample = new Sample(sampleFile, fingerprint, Array.toFloatArray(mapQs), Array.toFloatArray(normDepth), Array.toFloatArray(normDepth), Array.toFloatArray(percentWithMismatch), blankLRRs, genos, genos, false);
		sample.saveToRandomAccessFile(sampleFile, outliers, sampleFile);
		return outliers;
	}
}
