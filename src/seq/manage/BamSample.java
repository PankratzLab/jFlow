package seq.manage;

import java.util.Hashtable;

import cnv.analysis.BeastScore;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.filesys.Sample;
import common.Array;

public class BamSample {
	private static final double MAX_MAPQ = 60;
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

	private void init() {
		MarkerSet markerSet = proj.getMarkerSet();
		this.rawDepth = new double[bamPiles.length];
		this.mapQs = new double[bamPiles.length];
		this.percentWithMismatch = new double[bamPiles.length];
		for (int i = 0; i < bamPiles.length; i++) {
			rawDepth[i] = bamPiles[i].getOverallAvgDepth();
			mapQs[i] = Math.min(bamPiles[i].getOverallAvgMapQ() / MAX_MAPQ, 1);
			double percentMiss = bamPiles[i].getNumReadsOverlap() <= 0 ? 0 : (double) bamPiles[i].getNumReadsWithMismatch() / bamPiles[i].getNumReadsOverlap();

			percentWithMismatch[i] = percentMiss;

		}
		BeastScore beastScore = new BeastScore(Array.toFloatArray(rawDepth), markerSet.getIndicesByChr(), null, proj.getLog());
		this.normDepth = Array.toDoubleArray(beastScore.getinverseTransformedDataScaleMAD());

//		if(Array.min(normDepth)<0){
//			System.err.println("less dan 0");
//			System.exit(1);
//		}
		double maxNorm = Array.max(normDepth);
		double minNorm = Array.min(normDepth);
		for (int i = 0; i < normDepth.length; i++) {
			normDepth[i] = (normDepth[i] - minNorm) / (maxNorm - minNorm);
			normDepth[i] += 1;
			normDepth[i] *= 10;
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
