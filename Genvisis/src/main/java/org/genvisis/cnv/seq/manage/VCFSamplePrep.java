package org.genvisis.cnv.seq.manage;

import java.util.Arrays;
import java.util.Hashtable;
import java.util.Map;
import java.util.concurrent.Callable;
import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.qc.GcAdjustor;
import org.genvisis.cnv.qc.GcAdjustor.GC_CORRECTION_METHOD;
import org.genvisis.cnv.qc.GcAdjustor.GcModel;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.WorkerTrain.AbstractProducer;

public class VCFSamplePrep {

  private final Project proj;
  private final Sample samp;
  private GcModel gcModel;

  public enum PREPPED_SAMPLE_TYPE {
    NORMALIZED_GC_CORRECTED, NORMALIZED;
  }

  public VCFSamplePrep(Project proj, Sample samp, GcModel gcModel) {
    super();
    this.proj = proj;
    this.samp = samp;
    this.gcModel = gcModel;
    if (gcModel == null && Files.exists(proj.GC_MODEL_FILENAME.getValue())) {
      this.gcModel = GcAdjustor.GcModel.populateFromFile(proj.GC_MODEL_FILENAME.getValue(), false,
                                                         proj.getLog());
    }
  }

  public Sample getPreppedSample(PREPPED_SAMPLE_TYPE type) {
    float[] totalDepth = getTotalDepth(samp.getXs(), samp.getYs());
    double[] normDepth = normDepth(totalDepth);
    switch (type) {
      case NORMALIZED:
        break;
      case NORMALIZED_GC_CORRECTED:
        if (gcModel != null) {
          Map<Marker, Double> correctedIntensities = GcAdjustor.getComputedAdjustor(proj,
                                                                                    proj.getMarkerSet()
                                                                                        .mapProjectOrderData(normDepth),
                                                                                    gcModel,
                                                                                    GC_CORRECTION_METHOD.GENVISIS_GC,
                                                                                    true, true,
                                                                                    false)
                                                               .getCorrectedIntensities();
          normDepth = proj.getMarkerSet().markersAsList().stream()
                          .mapToDouble(correctedIntensities::get).toArray();
        } else {
          proj.getLog()
              .reportError("Projects gcmodel file must be valid for this method, skipping gc correction");

        }
        break;
      default:
        proj.getLog().reportError("Invalid sample type");
        break;
    }
    float[] xs = scale(samp.getXs(), totalDepth, normDepth);
    float[] ys = scale(samp.getYs(), totalDepth, normDepth);
    float[] fake = new float[xs.length];
    Arrays.fill(fake, 0);
    Sample preppedSamp = new Sample(samp.getSampleName(), samp.getFingerprint(), samp.getGCs(), xs,
                                    ys, fake, fake, samp.getForwardGenotypes(),
                                    samp.getAB_Genotypes(), samp.getCanXYBeNegative());
    return preppedSamp;
  }

  private static float[] scale(float[] d, float[] totalDepth, double[] normDepth) {
    double minNorm = ArrayUtils.min(normDepth);
    float[] scaleNorm = new float[d.length];
    for (int i = 0; i < scaleNorm.length; i++) {
      if (!Float.isFinite((float) normDepth[i])) {
        System.out.println(normDepth[i]);
        System.exit(1);
      }
      if (totalDepth[i] == 0) {
        scaleNorm[i] = (float) minNorm;
      } else if (d[i] == 0) {
        scaleNorm[i] = (float) minNorm;
      } else {
        scaleNorm[i] = (float) ((d[i] * normDepth[i]) / totalDepth[i]);
      }
    }
    double min = Math.abs(ArrayUtils.min(scaleNorm));
    for (int i = 0; i < scaleNorm.length; i++) {
      scaleNorm[i] += min;
    }

    return scaleNorm;
  }

  private static double[] normDepth(float[] totalDepth) {
    return ArrayUtils.normalize(ArrayUtils.toDoubleArray(totalDepth));
  }

  private static float[] getTotalDepth(float[] xs, float[] ys) {
    float[] totalDepth = new float[xs.length];
    for (int i = 0; i < totalDepth.length; i++) {
      totalDepth[i] = xs[i] + ys[i];
    }
    return totalDepth;
  }

  public static class VCFSamplePrepWorker extends AbstractProducer<Hashtable<String, Float>> {

    private final Project proj;
    private final PREPPED_SAMPLE_TYPE type;
    private final String[] samples;
    private int index;
    private final String sampleDir;
    private final GcModel gcModel;

    public VCFSamplePrepWorker(Project proj, String sampleDir, PREPPED_SAMPLE_TYPE type,
                               GcModel gcModel) {
      super();
      this.proj = proj;
      this.sampleDir = sampleDir;
      this.type = type;
      samples = proj.getSamples();
      this.gcModel = gcModel;
      index = 0;
    }

    @Override
    public boolean hasNext() {
      return index < samples.length;
    }

    @Override
    public Callable<Hashtable<String, Float>> next() {
      final String sample = samples[index];
      Callable<Hashtable<String, Float>> sampPrep = new Callable<Hashtable<String, Float>>() {

        @Override
        public Hashtable<String, Float> call() throws Exception {
          Hashtable<String, Float> outliers = new Hashtable<>();
          VCFSamplePrep prep = new VCFSamplePrep(proj,
                                                 proj.getFullSampleFromRandomAccessFile(sample),
                                                 gcModel);
          Sample prepped = prep.getPreppedSample(type);
          prepped.saveToRandomAccessFile(sampleDir + prepped.getSampleName()
                                         + Sample.SAMPLE_FILE_EXTENSION, outliers,
                                         prepped.getSampleName());
          byte[] genos = prepped.getAB_Genotypes();
          for (byte geno : genos) {
            if (geno == 3) {
              System.out.println("genotype 3 in sample after normalizing "
                                 + prepped.getSampleName());
              System.exit(1);
            }
          }
          return outliers;
        }
      };
      index++;
      return sampPrep;
    }
  }
}
