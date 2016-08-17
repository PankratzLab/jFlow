package org.genvisis.seq.manage;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Hashtable;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.common.Logger;

import com.google.common.primitives.Bytes;
import com.google.common.primitives.Floats;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public class SampleNGS {
  private enum DATA_TYPE {
                          GENO, X, Y, GC;
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

  private final String sampleName;
  private final ArrayList<Byte> geno;
  private final ArrayList<Float> xs;
  private final ArrayList<Float> ys;

  private final ArrayList<Float> gcs;

  public SampleNGS(String sampleName) {
    super();
    this.sampleName = sampleName;
    geno = new ArrayList<Byte>(2000000);
    xs = new ArrayList<Float>(2000000);
    ys = new ArrayList<Float>(2000000);
    gcs = new ArrayList<Float>(2000000);
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

  public void addGeno(VariantContext vc, Genotype geno, Logger log) {
    if (!geno.getSampleName().equals(sampleName)) {
      log.reportTimeError("Sample names do not match");
    } else {
      try {
        int[] ad =
            vc == null ? geno.getAD() : VCOps.getAppropriateAlleleDepths(vc, geno, true, log);
        addFloat(ad[0], DATA_TYPE.X, log);
        addFloat(ad[1], DATA_TYPE.Y, log);
      } catch (IllegalStateException ils) {
        log.reportTimeWarning("setting invalid allele depths to 0");
        addFloat(0, DATA_TYPE.X, log);
        addFloat(0, DATA_TYPE.Y, log);
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
    return GQ / 100;
  }

  private float convertXY(float xy) {
    return xy / 100;
  }

  public Hashtable<String, Float> dump(Project proj, Hashtable<String, Float> allOutliers,
                                       long fingerprint, Logger log) {
    if (!verify(proj)) {
      log.reportTimeError("Could not verify that all data has been added for sample " + sampleName);
    } else {

      String dir = proj.SAMPLE_DIRECTORY.getValue();
      float[] fakeLRRs = new float[gcs.size()];
      float[] fakeBAFS = new float[gcs.size()];
      Arrays.fill(fakeLRRs, -1);
      Arrays.fill(fakeBAFS, 0);
      Sample samp = new Sample(sampleName, fingerprint, Floats.toArray(gcs), Floats.toArray(xs),
                               Floats.toArray(ys), fakeBAFS, fakeLRRs, Bytes.toArray(geno),
                               Bytes.toArray(geno), false);
      samp.saveToRandomAccessFile(dir + sampleName + Sample.SAMPLE_FILE_EXTENSION, allOutliers,
                                  sampleName);
    }
    return allOutliers;
  }

  public String getSampleName() {
    return sampleName;
  }

  public boolean verify(Project proj) {
    int numMarkers = proj.getMarkerNames().length;
    boolean verified = true;
    if (numMarkers != gcs.size() || numMarkers != xs.size() || numMarkers != ys.size()
        || numMarkers != geno.size()) {
      verified = false;
    }
    return verified;
  }

}
