package seq.manage;

import htsjdk.variant.variantcontext.Genotype;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Hashtable;

import cnv.filesys.Project;
import cnv.filesys.Sample;
import common.Array;
import common.Logger;

public class SampleNGS {
	private enum DATA_TYPE {
		GENO, X, Y, GC;
	}

	private String sampleName;
	private ArrayList<Byte> geno;
	private ArrayList<Float> xs;
	private ArrayList<Float> ys;
	private ArrayList<Float> gcs;

	public SampleNGS(String sampleName) {
		super();
		this.sampleName = sampleName;
		this.geno = new ArrayList<Byte>(1000000);
		this.xs = new ArrayList<Float>(1000000);
		this.ys = new ArrayList<Float>(1000000);
		this.gcs = new ArrayList<Float>(1000000);
	}

	public String getSampleName() {
		return sampleName;
	}

	public void addGeno(Genotype geno, Logger log) {
		if (!geno.getSampleName().equals(sampleName)) {
			log.reportTimeError("Sample names do not match");
		} else {
			if (!geno.hasAD()) {
				log.reportTimeError("Genotype does not have AD annotation,setting X and Y to 0");
				addFloat(0, DATA_TYPE.X, log);
				addFloat(0, DATA_TYPE.Y, log);
			} else {
				addFloat(geno.getAD()[0], DATA_TYPE.X, log);
				addFloat(geno.getAD()[1], DATA_TYPE.Y, log);
			}
			if (!geno.hasGQ()) {
				addFloat(0, DATA_TYPE.GC, log);
				// log.reportTimeError("Genotype does not have GQ annotation, setting to 0");
			} else {
				addFloat(geno.getGQ(), DATA_TYPE.GC, log);
			}
			if (geno.isHomRef()) {
				addByte((byte) 0, DATA_TYPE.GENO, log);
			} else if (geno.isHet()) {
				addByte((byte) 1, DATA_TYPE.GENO, log);
			} else if (geno.isHomVar()) {
				addByte((byte) 2, DATA_TYPE.GENO, log);
			} else {
				addByte((byte) -1, DATA_TYPE.GENO, log);

			}

		}
	}

	private float convertGQ(float GQ) {
		return (float) GQ / 100;
	}

	private float convertXY(float xy) {
		return (float) xy / 100;
	}

	private void addFloat(float f, SampleNGS.DATA_TYPE type, Logger log) {
		switch (type) {
		case GC:
			gcs.add(convertGQ(f));
			break;
		case GENO:
			log.reportTimeError("Invalid data type for " + type);
			break;
		case X:
			xs.add(convertXY(f));
			break;
		case Y:
			ys.add(convertXY(f));
			break;
		default:
			log.reportTimeError("Invalid data type for " + type);
			break;
		}
	}

	private void addByte(byte b, SampleNGS.DATA_TYPE type, Logger log) {
		switch (type) {
		case GC:
			log.reportTimeError("Invalid data type for " + type);
			break;
		case GENO:
			geno.add(b);
			break;
		case X:
			log.reportTimeError("Invalid data type for " + type);
			break;
		case Y:
			log.reportTimeError("Invalid data type for " + type);
			break;
		default:
			log.reportTimeError("Invalid data type for " + type);
			break;
		}
	}

	public boolean verify(Project proj) {
		int numMarkers = proj.getMarkerNames().length;
		boolean verified = true;
		if (numMarkers != gcs.size() || numMarkers != xs.size() || numMarkers != ys.size() || numMarkers != geno.size()) {
			verified = false;
		}
		return verified;
	}

	public Hashtable<String, Float> dump(Project proj, Hashtable<String, Float> allOutliers, long fingerprint, Logger log) {
		if (!verify(proj)) {
			log.reportTimeError("Could not verify that all data has been added for sample " + sampleName);
		} else {

			String dir = proj.SAMPLE_DIRECTORY.getValue();
			float[] fakeLRRs = new float[gcs.size()];
			float[] fakeBAFS = new float[gcs.size()];
			Arrays.fill(fakeLRRs, -1);
			Arrays.fill(fakeBAFS, 0);
			Sample samp = new Sample(sampleName, fingerprint, Array.toFloatArray(gcs), Array.toFloatArray(xs), Array.toFloatArray(ys), fakeBAFS, fakeLRRs, Array.toByteArray(geno), Array.toByteArray(geno), false);
			samp.saveToRandomAccessFile(dir + sampleName + Sample.SAMPLE_DATA_FILE_EXTENSION, allOutliers, sampleName);
		}
		return allOutliers;
	}

	public static SampleNGS[] getSamples(HashSet<String> samples) {
		SampleNGS[] vcfSamples = new SampleNGS[samples.size()];
		int index = 0;
		for (String sample : samples) {
			vcfSamples[index] = new SampleNGS(sample);
			index++;
		}
		return vcfSamples;
	}

}