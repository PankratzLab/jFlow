package org.genvisis.one;

import java.io.File;
import java.io.PrintWriter;
import java.util.concurrent.Callable;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.WorkerTrain.AbstractProducer;

/**
 * @author lane0212 Taking a look at the intensity bias of X and Y probes
 */
public class XYIntensityBias {

  private static final String[] HEADER_BASE = new String[] {"Mean", "Median", "SD"};

  private static void run(Project proj) {

    XYProducer producer = new XYProducer(proj, proj.getSamples());
    try (WorkerTrain<double[][]> train = new WorkerTrain<double[][]>(producer, 6, 2,
                                                                     proj.getLog())) {
      String outDir = proj.PROJECT_DIRECTORY.getValue() + "xyComp/";
      new File(outDir).mkdirs();
      String out = outDir + "xyComp.txt";
      try {
        PrintWriter writer = Files.openAppropriateWriter(out);
        String[] xh = ArrayUtils.tagOn(HEADER_BASE, "X_", null);
        String[] yh = ArrayUtils.tagOn(HEADER_BASE, "Y_", null);

        writer.println("Sample\t" + ArrayUtils.toStr(xh) + "\t" + ArrayUtils.toStr(yh));
        int index = 0;
        while (train.hasNext()) {
          double[][] vals = train.next();
          writer.println(proj.getSamples()[index] + "\t" + ArrayUtils.toStr(vals[0]) + "\t"
                         + ArrayUtils.toStr(vals[1]));
          index++;
          proj.getLog().reportTimeInfo(index + "");
        }
        writer.close();
      } catch (Exception e) {
        proj.getLog().reportError("Error writing to " + out);
        proj.getLog().reportException(e);
      }
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = null;

    String usage = "\n" + "one.XYIntensityBias requires 0-1 arguments\n" + "   (1) proj (i.e. proj="
                   + filename + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("proj=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      Project proj = new Project(filename);
      run(proj);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  private static class XYProducer extends AbstractProducer<double[][]> {

    private final Project proj;
    private final String[] samples;
    private int index;

    public XYProducer(Project proj, String[] samples) {
      super();
      this.proj = proj;
      this.samples = samples;
      index = 0;
    }

    @Override
    public boolean hasNext() {
      return index < samples.length;
    }

    @Override
    public Callable<double[][]> next() {
      final String currentSamp = samples[index];
      Callable<double[][]> callable = new Callable<double[][]>() {

        @Override
        public double[][] call() throws Exception {
          return compute(proj, currentSamp);
        }

      };
      index++;
      // TODO Auto-generated method stub
      return callable;
    }
  }

  private static double[][] compute(Project proj, String sample) {
    double[][] vals = new double[2][3];
    Sample samp = proj.getFullSampleFromRandomAccessFile(sample);
    int[] auto = proj.getAutosomalMarkerIndices();
    float[] xauto = ArrayUtils.subArray(samp.getXs(), auto);
    float[] yauto = ArrayUtils.subArray(samp.getYs(), auto);
    vals[0] = getVals(ArrayUtils.toDoubleArray(xauto));
    vals[1] = getVals(ArrayUtils.toDoubleArray(yauto));
    return vals;
  }

  private static double[] getVals(double[] in) {
    double[] vals = new double[3];
    double[] tmp = ArrayUtils.removeNaN(in);
    vals[0] = ArrayUtils.mean(tmp);
    vals[1] = ArrayUtils.median(tmp);
    vals[2] = ArrayUtils.stdev(tmp);
    return vals;
  }

}
