package cnv.filesys;

import java.io.Serializable;
import common.Files;
import filesys.Segment;

public class MeanLRRset implements Serializable {
	public static final long serialVersionUID = 1L;

	private long sampleFingerprint;
	private Segment[] regions;
	private int[] numberOfMarkers;
	private float[][][] data;
	
	public MeanLRRset(long sampleFingerprint, Segment[] regions, int[] numberOfMarkers, float[][][] data) {
		this.sampleFingerprint = sampleFingerprint;
		this.numberOfMarkers = numberOfMarkers;
		this.regions = regions;
		this.data = data;
	}

	public long getSampleFingerprint() {
		return sampleFingerprint;
	}

	public Segment[] getRegions() {
		return regions;
	}

	public int[] getNumerOfMarkersPerRegion() {
		return numberOfMarkers;
	}

	public float[][][] getData() {
		return data;
	}
	
	public void serialize(String filename) {
		Files.writeSerial(this, filename);
	}

	public static MeanLRRset load(String filename, boolean jar) {
		return (MeanLRRset)Files.readSerial(filename, jar, true);
	}
}
