package cnv.var;

public class MosaicRegion extends CNVariant {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	private double nearestStateScore;

	public MosaicRegion(CNVariant cnv, double nearestStateScore) {
		super(cnv);
		this.nearestStateScore = nearestStateScore;
	}

	public double getNearestStateScore() {
		return nearestStateScore;
	}

}
