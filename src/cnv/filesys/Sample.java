package cnv.filesys;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.Serializable;

import common.Array;
import common.DoubleVector;
import common.Files;

public class Sample implements Serializable {
	public static final long serialVersionUID = 1L;

	private long fingerprint;
	private float[] lrrs;
	private float[] bafs;
	private byte[] genotypes;

	public Sample(long fingerprint, float[] lrrs, float[] bafs, byte[] genotypes) {
		this.fingerprint = fingerprint;
		this.lrrs = lrrs;
		this.bafs = bafs;
		this.genotypes = genotypes;
	}

	public long getFingerprint() {
		return fingerprint;
	}

	public float[] getLRRs() {
		return lrrs;
	}

	public float[] getBAFs() {
		return bafs;
	}

	public byte[] getGenotypes() {
		return genotypes;
	}

	public void writeToFile(String[] markerNames, String filename) {
		PrintWriter writer;

		if (markerNames.length!=lrrs.length) {
			System.err.println("Error - MarkerNames (n="+markerNames.length+") do not match up with the number of LRRs/BAFs (n="+lrrs.length+")");
			System.exit(1);
		}
		try {
			writer = new PrintWriter(new FileWriter(filename));
			writer.println("SNP\tLRR\tBAF\tGenotype");
			for (int i = 0; i<markerNames.length; i++) {
				writer.println(markerNames[i]+"\t"+lrrs[i]+"\t"+bafs[i]+"\t"+genotypes[i]);
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing "+filename);
			e.printStackTrace();
		}
	}

	public void serialize(String filename) {
		Files.writeSerial(this, filename);
	}

	public static Sample load(String filename, boolean jar) {
		return (Sample)Files.readSerial(filename, jar, true);
	}

	public boolean hasBimodalBAF(byte chr, int startPosition, int endPosition) {
		DoubleVector dv = new DoubleVector();
		for (int i=0; i<bafs.length; i++) {
			if (bafs[i]>0.15 && bafs[i]<0.85) {
				dv.add(bafs[i]);
			}
		}
		return Array.isBimodal(dv.toArray());
	}

}
