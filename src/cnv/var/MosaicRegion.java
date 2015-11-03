package cnv.var;

import java.util.ArrayList;

import common.Array;

public class MosaicRegion extends CNVariant {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	public static final String[] ADD_HEADER = new String[] { "bpWeightedScore", "nearestStatScore", "pdfScore", "delta", "f", "customF" };
	private double nearestStateScore;
	private double bpWeightedScore;
	private double pdfScore;
	private double delta;
	private double f;
	private double customF;

	public MosaicRegion(CNVariant cnv, double bpWeightedScore, double nearestStateScore, double pdfScore, double delta, double f, double customF) {
		super(cnv);
		this.bpWeightedScore = bpWeightedScore;
		this.nearestStateScore = nearestStateScore;
		this.pdfScore = pdfScore;
		this.delta = delta;
		this.f = f;
		this.customF = customF;
	}

	public MosaicRegion(CNVariant cnv, MosaicRegion another) {
		super(cnv);
		this.bpWeightedScore = another.bpWeightedScore;
		this.nearestStateScore = another.nearestStateScore;
		this.pdfScore = another.pdfScore;
		this.delta = another.delta;
		this.f = another.f;
		this.customF = another.customF;
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

		String[] s = Array.concatAll(toPlinkFormat().split("\t"), Array.toStringArray(tmp));
		return Array.toStr(s);
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
