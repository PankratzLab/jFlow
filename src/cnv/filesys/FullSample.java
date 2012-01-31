package cnv.filesys;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.Serializable;
import common.Files;

public class FullSample implements Serializable {
	public static final long serialVersionUID = 1L;
	public static final String[][] DATA_FIELDS = {{"GC Score", "GCscore"}, {"X Raw"}, {"Y Raw"}, {"X", "Xvalue"}, {"Y", "Yvalue"}, {"Theta"}, {"R"}, {"B Allele Freq"}, {"Log R Ratio"}};
	public static final String[][] GENOTYPE_FIELDS = {{"Allele1 - Forward", "Allele1"}, {"Allele2 - Forward", "Allele2"}, {"Allele1 - AB"}, {"Allele2 - AB"}};
	public static final String[] ALL_STANDARD_GENOTYPE_FIELDS = {"Allele1 - AB", "Allele2 - AB", "Allele1 - Forward", "Allele2 - Forward", "Allele1 - Top", "Allele2 - Top", "Allele1 - Design", "Allele2 - Design"};
	public static final String[] ALLELE_PAIRS = {"--", "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT", "DD", "DI", "II", "ID"};
	public static final String[] ALT_NULL = {"-", "0"};
	public static final String[] ALT_NULLS = {"--", "00"};
	public static final String[] AB_PAIRS = {"AA", "AB", "BB"};

	private String sampleName;
	private long fingerprint;
	private float[] gcs;
	private float[] xs;
	private float[] ys;
	private float[] xRaws;
	private float[] yRaws;
	private float[] thetas;
	private float[] rs;
	private float[] lrrs;
	private float[] bafs;
	private byte[] forwardGenotypes;
	private byte[] abGenotypes;

	public FullSample(String sampleName, long fingerprint, float[] gcs, float[] xRaws, float[] yRaws, float[] xs, float[] ys, float[] rs, float[] thetas, float[] bafs, float[] lrrs, byte[] forwardGenotypes, byte[] abGenotypes) {
		this.sampleName = sampleName;
		this.fingerprint = fingerprint;
		this.gcs = gcs;
		this.xRaws = xRaws;
		this.yRaws = yRaws;
		this.xs = xs;
		this.ys = ys;
		this.rs = rs;
		this.thetas = thetas;
		this.rs = rs;
		this.bafs = bafs;
		this.lrrs = lrrs;
		this.forwardGenotypes = forwardGenotypes;
		this.abGenotypes = abGenotypes;
	}

	public FullSample(String sampleName, long fingerprint, float[][] data, byte[][] genotypes) {
		this.sampleName = sampleName;
		this.fingerprint = fingerprint;
		this.gcs = data[0];
		this.xRaws = data[1];
		this.yRaws = data[2];
		this.xs = data[3];
		this.ys = data[4];
		this.thetas = data[5];
		this.rs = data[6];
		this.bafs = data[7];
		this.lrrs = data[8];
		this.forwardGenotypes = genotypes[0];
		this.abGenotypes = genotypes[1];
	}

	public float[][] getAllData() {
		return new float[][] {gcs, xRaws, yRaws, xs, ys, thetas, rs, bafs, lrrs};
	}

	public byte[][] getAllGenotypes() {
		return new byte[][] {forwardGenotypes, abGenotypes};
	}

	public String getSampleName() {
		return sampleName;
	}

	public long getFingerprint() {
		return fingerprint;
	}

	public float[] getGCs() {
		return gcs;
	}

	public float[] getX_Raws() {
		return xRaws;
	}

	public float[] getY_Raws() {
		return yRaws;
	}

	public float[] getXs() {
		return xs;
	}

	public float[] getYs() {
		return ys;
	}

	public float[] getThetas() {
		if (thetas == null) {
			thetas = new float[xs.length];
			for (int i = 0; i<xs.length; i++) {
				thetas[i] = Centroids.calcTheta(xs[i], ys[i]);
            }
		}
		return thetas;
	}

	public float[] getRs() {
		if (rs == null) {
			rs = new float[xs.length];
			for (int i = 0; i<xs.length; i++) {
				rs[i] = Centroids.calcR(xs[i], ys[i]);
            }
		}
		return rs;
	}

	public float[] getBAFs() {
		return bafs;
	}

	public float[] getBAFs(float[][][] centroids) {
		float[] thetas, bafs;

		thetas = getThetas();
		bafs = new float[xs.length];
		for (int i = 0; i<xs.length; i++) {
			bafs[i] = Centroids.calcBAF(thetas[i], centroids[i]);
        }

		return bafs;
	}

	public float[] getLRRs() {
		return lrrs;
	}

	public float[] getLRRs(float[][][] centroids) {
		float[] thetas, rs, lrrs;

		thetas = getThetas();
		rs = getRs();
		lrrs = new float[xs.length];
		for (int i = 0; i<xs.length; i++) {
			lrrs[i] = Centroids.calcLRR(thetas[i], rs[i], centroids[i]);
        }

		return lrrs;
	}

	public byte[] getForwardGenotypes() {
		return forwardGenotypes;
	}

	public byte[] getAB_Genotypes() {
		return abGenotypes;
	}

	public void writeToFile(String[] markerNames, String filename) {
		PrintWriter writer;

		if (markerNames.length!=lrrs.length) {
			System.err.println("Error - MarkerNames (n="+markerNames.length+") do not match up with the number of LRRs/BAFs/etc (n="+lrrs.length+")");
			System.exit(1);
		}
		try {
			writer = new PrintWriter(new FileWriter(filename));
			writer.println("SNP\tGC Score\tRaw X\tRaw Y\tX\tY\tTheta\tR\tLRR\tBAF\tGenotypes");
			for (int i = 0; i<markerNames.length; i++) {
				writer.println(markerNames[i]+"\t"+gcs[i]+"\t"+xRaws[i]+"\t"+yRaws[i]+"\t"+xs[i]+"\t"+ys[i]+"\t"+thetas[i]+"\t"+rs[i]+"\t"+lrrs[i]+"\t"+bafs[i]+"\t"+ALLELE_PAIRS[forwardGenotypes[i]]);
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing "+filename);
			e.printStackTrace();
		}
	}

	public float[][][] loadCentroids(Project proj) {
		Centroids centroids;
		
		centroids = Centroids.load(proj.getFilename(Project.CUSTOM_CENTROIDS_FILENAME), proj.getJarStatus());
		if (centroids.getFingerprint() != getFingerprint()) {
			System.err.println("Error - mismatched fingerprint for "+sampleName);
		}
		return centroids.getCentroids();
	}
	
	public void compareCalculationsFile(Project proj, String[] markerNames, String filename) {
		PrintWriter writer;
		float[][][] centroids;
		float[] compBAFs, compLRRs;

		centroids = loadCentroids(proj);
		compBAFs = getBAFs(centroids);
		compLRRs = getLRRs(centroids);
		
		if (markerNames.length!=lrrs.length) {
			System.err.println("Error - MarkerNames (n="+markerNames.length+") do not match up with the number of LRRs/BAFs/etc (n="+lrrs.length+")");
			System.exit(1);
		}
		try {
			writer = new PrintWriter(new FileWriter(filename));
			writer.println("SNP\tX\tY\tTheta\tR\tcompTheta\tcompR\tLRR\tBAF\tcompLRR\tcompBAF");
			for (int i = 0; i<markerNames.length; i++) {
//			for (int i = 0; i<1000; i++) {
				writer.println(markerNames[i]+"\t"+xs[i]+"\t"+ys[i]+"\t"+thetas[i]+"\t"+rs[i]+"\t"+bafs[i]+"\t"+lrrs[i]+"\t"+compBAFs[i]+"\t"+compLRRs[i]);
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing "+filename);
			e.printStackTrace();
		}
	}

	public void compLRRs(Project proj) {
		float[] compLRRs, diffs;

		compLRRs = getLRRs(loadCentroids(proj));
		diffs = new float[lrrs.length];
		for (int i = 0; i<lrrs.length; i++) {
			diffs[i] = lrrs[i]-compLRRs[i];
		}
		Files.writeSerial(diffs, proj.getProjectDir()+"comps/"+sampleName+".comp");
	}

	public void serialize(String filename) {
		Files.writeSerial(this, filename);
	}

	public Sample convertToSample() {
		return new Sample(fingerprint, lrrs, bafs, abGenotypes);
	}

	public static FullSample load(String filename, boolean jar) {
		return (FullSample)Files.readSerial(filename, jar, true);
	}
	
	public static void main(String[] args) {
		Project proj = new Project(Project.DEFAULT_CURRENT, false);
		String[] samples = proj.getSamples();
		FullSample fsamp;
		
		for (int i = 0; i<samples.length; i++) {
			fsamp = proj.getFullSample(samples[i]);
			fsamp.compareCalculationsFile(proj, proj.getMarkerNames(), proj.getProjectDir()+samples[i]+"_comp.xln");
        }
    }
}
