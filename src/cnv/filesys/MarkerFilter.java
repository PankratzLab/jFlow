package cnv.filesys;

import java.io.Serializable;

public class MarkerFilter implements Serializable {
	private static final long serialVersionUID = 1L;
	
	private byte coordType; // e.g. X/Y, Theta/R, BAF/LRR
	private float[][] coords; // already translated from pixels into values
	private byte newClass;
	
	public MarkerFilter() {}
	

}
