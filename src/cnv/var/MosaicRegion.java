package cnv.var;

import filesys.Segment;

public class MosaicRegion extends Segment {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	private String sample;
	private double trisomyDisomyF;
	private double monsomyDisomyF;
	private double shift;

	public MosaicRegion(String sample, double trisomyDisomyF, double monsomyDisomyF, double shift, byte chr, int start, int stop) {
		super(chr, start, stop);
		this.sample = sample;
		this.trisomyDisomyF = trisomyDisomyF;
		this.monsomyDisomyF = monsomyDisomyF;
		this.shift = shift;
	}

}
