package cnv.filesys;

import java.io.Serializable;
import common.Files;
import filesys.Segment;

public class MeanLRRset implements Serializable {
	public static final long serialVersionUID = 1L;

	private long sampleFingerprint;
	private Segment[] cnps;
	private int[] numberOfMarkers;
	private float[][] data;
	
	public MeanLRRset(long sampleFingerprint, Segment[] cnps, int[] numberOfMarkers, float[][] data) {
		this.sampleFingerprint = sampleFingerprint;
		this.numberOfMarkers = numberOfMarkers;
		this.cnps = cnps;
		this.data = data;
	}

	public long getSampleFingerprint() {
		return sampleFingerprint;
	}

	public Segment[] getCnps() {
		return cnps;
	}

	public int[] getNumerOfMarkersPerCNV() {
		return numberOfMarkers;
	}

	public float[][] getData() {
		return data;
	}
	
	public void serialize(String filename) {
		Files.writeSerial(this, filename);
	}

	public static MeanLRRset load(String filename, boolean jar) {
		return (MeanLRRset)Files.readSerial(filename, jar, true);
	}
}
