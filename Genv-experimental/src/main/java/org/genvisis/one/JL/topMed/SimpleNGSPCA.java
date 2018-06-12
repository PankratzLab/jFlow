package org.genvisis.one.JL.topMed;

import java.io.File;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;
import org.apache.commons.math3.linear.DiagonalMatrix;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.ejml.alg.dense.decomposition.svd.SvdImplicitQrDecompose_D64;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.SingularOps;
import org.genvisis.CLI;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.VCFOps.VcfPopulation;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import org.genvisis.stats.Maths;

/**
 * A simplified version of BamImport that uses MosDepth output to generate PCS
 */
public class SimpleNGSPCA implements Serializable {

  /**
   * 
   */
  private static final long serialVersionUID = 5517545544997813383L;
  private final RealMatrix m;
  private final String[] colNames;
  private final String[] rowNames;

  //  With V,W, and original data M we can always compute U
  private RealMatrix v;
  private DiagonalMatrix w;

  /**
   * @param m
   * @param colNames
   * @param rowNames
   */
  private SimpleNGSPCA(RealMatrix m, String[] colNames, String[] rowNames) {
    super();
    this.m = m;
    this.colNames = colNames;
    this.rowNames = rowNames;
    if (m.getColumnDimension() != colNames.length) {
      throw new IllegalArgumentException("Mismatched column lengths");
    }
    if (m.getRowDimension() != rowNames.length) {
      throw new IllegalArgumentException("Mismatched row lengths");
    }
  }

  private void computeSVD(Logger log) {
    log.reportTimeInfo("computing SVD base");

    DenseMatrix64F a = new DenseMatrix64F(1, 1);
    a.reshape(m.getRowDimension(), m.getColumnDimension());
    for (int row = 0; row < m.getRowDimension(); row++) {
      for (int col = 0; col < m.getColumnDimension(); col++) {
        a.add(row, col, m.getEntry(row, col));
      }
    }
    log.reportTimeInfo("Computing EJML PCs");
    SvdImplicitQrDecompose_D64 svd = new SvdImplicitQrDecompose_D64(false, false, true, false);
    svd.decompose(a);
    log.reportTimeInfo("Finished Computing EJML PCs");

    log.reportTimeInfo("finished computing SVD base");

    DenseMatrix64F tv = svd.getV(null, true);

    DenseMatrix64F tmpW = svd.getW(null);
    SingularOps.descendingOrder(null, false, tmpW, tv, true);
    int numSingular = Math.min(tmpW.numRows, tmpW.numCols);
    double[] singularValues = new double[numSingular];
    for (int i = 0; i < numSingular; i++) {
      singularValues[i] = tmpW.get(i, i);
    }
    this.w = new DiagonalMatrix(singularValues);
    this.v = MatrixUtils.createRealMatrix(tv.numRows, tv.numCols);
    for (int row = 0; row < tv.numRows; row++) {
      for (int col = 0; col < tv.numCols; col++) {
        v.addToEntry(row, col, tv.get(row, col));
      }
    }
  }

  private enum SCALE_METHOD {
    FC_MEDIAN, CENTER_SCALE_SAMPLE, CENTER_SCALE_SAMPLE_SCALE_MARKER;
  }

  private void dumpPCsToText(String file) {

    PrintWriter writer = Files.getAppropriateWriter(file);
    StringJoiner joiner = new StringJoiner("\t");
    joiner.add("SAMPLE");
    for (int i = 0; i < v.getColumnDimension(); i++) {
      joiner.add("PC" + (i + 1));
    }
    writer.println(joiner.toString());

    for (int i = 0; i < v.getRowDimension(); i++) {
      StringJoiner sample = new StringJoiner("\t");
      sample.add(colNames[i]);
      for (int j = 0; j < v.getColumnDimension(); j++) {
        sample.add(Double.toString(v.getEntry(j, i)));
      }
      writer.println(sample.toString());

    }
    writer.close();
  }

  /**
   * compute fold change for sample and center to marker median
   * 
   * @param m an {@link RealMatrix}
   */
  private static void foldChangeAndCenter(RealMatrix m) {
    double[] medians = new double[m.getColumnDimension()]; //Samples are columns

    //    convert samples to log2 fold-change from median
    for (int j = 0; j < m.getColumnDimension(); j++) {
      double[] tmp = new double[m.getRowDimension()];//Markers are rows
      for (int i = 0; i < m.getRowDimension(); i++) {
        tmp[i] += m.getEntry(i, j);
      }
      medians[j] = ArrayUtils.median(tmp);
    }
    for (int i = 0; i < m.getRowDimension(); i++) {
      for (int j = 0; j < m.getColumnDimension(); j++) {
        double entry = m.getEntry(i, j);
        if (entry > 0) {
          double standard = Maths.log2(entry / medians[j]);
          m.setEntry(i, j, standard);
        } else {
          m.setEntry(i, j, 0);
        }
      }
    }
    //    convert markers to distance from median
    scaleMarkersToMedian(m);
  }

  private static void scaleMarkersToMedian(RealMatrix m) {
    for (int i = 0; i < m.getRowDimension(); i++) {
      double[] tmp = m.getRow(i);
      double median = ArrayUtils.median(tmp);
      for (int j = 0; j < m.getColumnDimension(); j++) {
        m.setEntry(i, j, tmp[j] - median);
      }
    }
  }

  private static void scale(RealMatrix m) {
    double[] sds = new double[m.getColumnDimension()];
    double[] mean = new double[m.getColumnDimension()];

    for (int column = 0; column < mean.length; column++) {
      double[] tmp = new double[m.getRowDimension()];
      for (int i = 0; i < m.getRowDimension(); i++) {
        tmp[i] += m.getEntry(i, column);
      }
      mean[column] = ArrayUtils.mean(tmp);
      sds[column] = ArrayUtils.stdev(tmp);
    }
    for (int row = 0; row < m.getRowDimension(); row++) {
      for (int column = 0; column < mean.length; column++) {
        double standard = m.getEntry(row, column) - mean[column];
        standard /= sds[column];
        m.setEntry(row, column, standard);
      }
    }
  }

  public static void main(String[] args) {
    CLI c = new CLI(SimpleNGSPCA.class);
    c.addArg(CLI.ARG_OUTDIR, CLI.DESC_OUTDIR, true);
    c.addArg(CLI.ARG_INDIR, CLI.DESC_INDIR, true);
    c.addArg("vpop",
             "population definition, samples with the same super population will be PC'ed together",
             true);

    c.parseWithExit(args);
    String rootoutDir = c.get(CLI.ARG_OUTDIR);
    new File(rootoutDir).mkdirs();
    Logger log = new Logger(rootoutDir + "log.log");

    VcfPopulation vpop = VcfPopulation.load(c.get("vpop"), POPULATION_TYPE.ANY, log);
    vpop.report();

    for (SCALE_METHOD method : SCALE_METHOD.values()) {
      Set<String> allSamps = new HashSet<>();
      for (String superPop : vpop.getSuperPop().keySet()) {
        Set<String> samps = vpop.getSuperPop().get(superPop);
        if (samps.size() > 1) {
          runGroup(c, method, rootoutDir, log, superPop, samps);
        } else {
          log.reportTimeWarning("Skipping super population " + superPop + ", not enough samples");
        }
        allSamps.addAll(samps);
      }
      runGroup(c, method, rootoutDir, log, "ALL_SAMPLES", allSamps);
    }

  }

  static void runGroup(CLI c, SCALE_METHOD method, String rootoutDir, Logger log, String superPop,
                       Set<String> samps) {
    String outDir = rootoutDir + superPop + "/" + method + "/";
    new File(outDir).mkdirs();
    log.reportTimeInfo("Operating in " + outDir);
    String matFile = outDir + "mat.ser.gz";
    String svdFile = ext.addToRoot(matFile, ".svd");
    String pcFile = outDir + "pcs.gz";

    if (!Files.exists(matFile)) {
      List<String> files = new ArrayList<>();

      List<String> tmpFiles = Arrays.asList(Files.listFullPaths(c.get(CLI.ARG_INDIR), "bed.gz"));
      for (String tmpFile : tmpFiles) {
        for (String sample : samps) {
          if (tmpFile.contains(sample)) {
            files.add(tmpFile);
            break;
          }
        }
      }

      log.reportTimeInfo("Found " + files.size() + " files to PCA for " + superPop);
      String[] availableRows = HashVec.loadFileToStringArray(files.get(0), false,
                                                             new int[] {0, 1, 2}, false);
      List<String> rows = new ArrayList<>();
      List<Integer> indices = new ArrayList<>();
      for (int i = 0; i < availableRows.length; i++) {
        String[] split = availableRows[i].trim().split("\t");

        Segment seg = new Segment(split[0], split[1], split[2]);
        //        or other filter logic
        if (seg.getChr() > 0 && seg.getChr() < 23) {
          rows.add(seg.getUCSClocation());
          indices.add(i);
        }
      }
      log.reportTimeInfo("Preparing matrix with " + rows.size() + " rows and " + files.size()
                         + "columns");

      RealMatrix m = MatrixUtils.createRealMatrix(rows.size(), files.size());
      String[] colNames = new String[files.size()];
      String[] rowNames = ArrayUtils.toStringArray(rows);
      for (int i = 0; i < colNames.length; i++) {
        colNames[i] = ext.rootOf(files.get(i));
        log.reportTime("Parsing sample " + (i + 1) + ", " + colNames[i]);
        String[] data = HashVec.loadFileToStringArray(files.get(i), false, new int[] {3}, false);
        for (int j = 0; j < indices.size(); j++) {
          m.addToEntry(j, i, Double.parseDouble(data[indices.get(j)]));
        }
      }
      //      new SingularValueDecomposition(

      log.reportTime("Saving progress");

      SerializedFiles.writeSerial(new SimpleNGSPCA(m, colNames, rowNames), matFile, true);
    }
    if (!Files.exists(svdFile)) {
      log.reportTime("Loading data");
      SimpleNGSPCA tm = (SimpleNGSPCA) SerializedFiles.readSerial(matFile, log, false);
      log.reportTime("Scaling data");
      switch (method) {
        case CENTER_SCALE_SAMPLE:
          scale(tm.m);
          break;
        case CENTER_SCALE_SAMPLE_SCALE_MARKER:
          scale(tm.m);
          scaleMarkersToMedian(tm.m);
          break;
        case FC_MEDIAN:
          foldChangeAndCenter(tm.m);

          break;
        default:
          throw new IllegalArgumentException("Invalid method " + method);

      }

      log.reportTimeInfo("Rows=" + tm.m.getRowDimension());
      log.reportTimeInfo("Cols=" + tm.m.getColumnDimension());
      log.reportTime("Computing SVD");
      tm.computeSVD(log);

      log.reportTime("Finished Computing SVD");

      log.reportTime("Saving progress");

      SerializedFiles.writeSerial(tm, svdFile, true);
    }
    log.reportTimeInfo("Loading " + svdFile);
    SimpleNGSPCA tm = (SimpleNGSPCA) SerializedFiles.readSerial(svdFile, log, false);
    tm.dumpPCsToText(pcFile);
  }
}