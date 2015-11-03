package cnv.var;

public class MosaicRegion extends CNVariant {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	private double nearestStateScore;
	private double bpWeightedScore;
	private double pdfScore;

	public MosaicRegion(CNVariant cnv, double bpWeightedScore, double nearestStateScore, double pdfScore) {
		super(cnv);
		this.bpWeightedScore = bpWeightedScore;
		this.nearestStateScore = nearestStateScore;
		this.pdfScore = pdfScore;
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
