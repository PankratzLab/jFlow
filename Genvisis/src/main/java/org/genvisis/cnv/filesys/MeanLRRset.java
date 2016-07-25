package org.genvisis.cnv.filesys;

import java.io.Serializable;

import org.genvisis.common.SerializedFiles;
import org.genvisis.filesys.Segment;

public class MeanLRRset implements Serializable {
	public static final long serialVersionUID = 1L;

	private long sampleFingerprint;
	private Segment[] regions;
	private int[] numberOfMarkers;
	private float[][][] data;
	private String[] transformations;
	
	public MeanLRRset(long sampleFingerprint, Segment[] regions, int[] numberOfMarkers, float[][][] data, String[] transformations) {
		this.sampleFingerprint = sampleFingerprint;
		this.numberOfMarkers = numberOfMarkers;
		this.regions = regions;
		this.data = data;
		this.transformations = transformations;
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
	
	public String[] getTransformations() {
		return transformations;
	}
	
	public void serialize(String filename) {
		SerializedFiles.writeSerial(this, filename);
	}

	public static MeanLRRset load(String filename, boolean jar) {
		return (MeanLRRset)SerializedFiles.readSerial(filename, jar, true);
	}
}
