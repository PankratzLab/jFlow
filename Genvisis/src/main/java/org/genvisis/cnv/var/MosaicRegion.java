package org.genvisis.cnv.var;

import java.io.Serializable;
import java.util.ArrayList;

import org.genvisis.common.ArrayUtils;
import org.genvisis.filesys.CNVariant;

public class MosaicRegion extends CNVariant implements Serializable {

	/**
	 *
	 */
	private static final long serialVersionUID = 1L;
	public static final String[] ADD_HEADER =
																					new String[] {"bpWeightedScore", "nearestStatScore",
																												"pdfScore", "delta", "f", "customF",
																												"numFMarkers", "numCustomFMarkers",
																												"beastScore", "beastHeight", "beastLength"};
	private final double nearestStateScore;
	private final double bpWeightedScore;
	private final double pdfScore;
	private final double delta;
	private double f;
	private final double customF;
	private int numFMarkers;
	private int numCustomFMarkers;
	private int beastLength;
	private double beastHeight;
	private double beastScore;

	public MosaicRegion(CNVariant cnv, double bpWeightedScore, double nearestStateScore,
											double pdfScore, double delta, double f, double customF) {
		super(cnv);
		this.bpWeightedScore = bpWeightedScore;
		this.nearestStateScore = nearestStateScore;
		this.pdfScore = pdfScore;
		this.delta = delta;
		this.f = f;
		this.customF = customF;
		numFMarkers = 0;
		numCustomFMarkers = 0;
		beastScore = 0;
	}

	public MosaicRegion(CNVariant cnv, MosaicRegion another) {
		super(cnv);
		bpWeightedScore = another.bpWeightedScore;
		nearestStateScore = another.nearestStateScore;
		pdfScore = another.pdfScore;
		delta = another.delta;
		f = another.f;
		customF = another.customF;
	}

	public void setF(double f) {
		this.f = f;
	}



	public void setBeastLength(int beastLength) {
		this.beastLength = beastLength;
	}

	public void setBeastHeight(double beastHeight) {
		this.beastHeight = beastHeight;
	}

	public void setBeastScore(double beastScore) {
		this.beastScore = beastScore;
	}

	public void setNumFMarkers(int numFMarkers) {
		this.numFMarkers = numFMarkers;
	}

	public void setNumCustomFMarkers(int numCustomFMarkers) {
		this.numCustomFMarkers = numCustomFMarkers;
	}

	@Override
	public String toAnalysisString() {
		ArrayList<String> tmp = new ArrayList<String>();
		tmp.add(bpWeightedScore + "");
		tmp.add(nearestStateScore + "");
		tmp.add(pdfScore + "");
		tmp.add(delta + "");
		tmp.add(f + "");
		tmp.add(customF + "");
		tmp.add(numFMarkers + "");
		tmp.add(numCustomFMarkers + "");
		tmp.add(beastScore + "");
		tmp.add(beastHeight + "");
		tmp.add(beastLength + "");

		String[] s = ArrayUtils.concatAll(toPlinkFormat().split("\t"), ArrayUtils.toStringArray(tmp));
		return ArrayUtils.toStr(s);
	}

	@Override
	public String[] getHeader() {
		return ArrayUtils.concatAll(PLINK_CNV_HEADER, ADD_HEADER);

	}


	public double getCustomF() {
		return customF;
	}

	public double getPdfScore() {
		return pdfScore;
	}

	public double getNearestStateScore() {
		return nearestStateScore;
	}

	public double getBpWeightedScore() {
		return bpWeightedScore;
	}

}
