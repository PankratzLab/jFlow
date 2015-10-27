package cnv.var;

public class MosaicRegion extends CNVariant {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	private double nearestStateScore;
	private double bpWeightedScore;

	public MosaicRegion(CNVariant cnv, double bpWeightedScore, double nearestStateScore) {
		super(cnv);
		this.bpWeightedScore = bpWeightedScore;
		this.nearestStateScore = nearestStateScore;
	}

	public double getNearestStateScore() {
		return nearestStateScore;
	}

	public double getBpWeightedScore() {
		return bpWeightedScore;
	}

}
