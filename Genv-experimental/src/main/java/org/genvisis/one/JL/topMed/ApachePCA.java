package org.genvisis.one.JL.topMed;

import java.io.File;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.StringJoiner;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.QRDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.ejml.alg.dense.decomposition.svd.SvdImplicitQrDecompose_D64;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;

public class ApachePCA {

  private static class TOPMedPCAMatrix implements Serializable {

    /**
     * 
     */
    private static final long serialVersionUID = 5517545544997813383L;
    private final RealMatrix m;
    private final String[] colNames;
    private final String[] rowNames;
    RealMatrix s;
    RealMatrix u;
    RealMatrix v;

    /**
     * @param m
     * @param colNames
     * @param rowNames
     */
    private TOPMedPCAMatrix(RealMatrix m, String[] colNames, String[] rowNames) {
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
      SingularValueDecomposition svd = new SingularValueDecomposition(m);

      log.reportTimeInfo("finished computing SVD base");

      s = svd.getS();
      u = svd.getU();
      v = svd.getV();
    }

    private void dumpPCsToText(String file) {

      PrintWriter writer = Files.getAppropriateWriter(file);
      StringJoiner joiner = new StringJoiner("\t");
      joiner.add("DNA");
      for (int i = 0; i < v.getColumnDimension(); i++) {
        joiner.add("PC" + (i + 1));
      }
      writer.println(joiner.toString());

      for (int i = 0; i < v.getRowDimension(); i++) {
        StringJoiner sample = new StringJoiner("\t");
        sample.add(colNames[i]);
        for (int j = 0; j < v.getColumnDimension(); j++) {
          sample.add(Double.toString(v.getEntry(i, j)));
        }
        writer.println(sample.toString());

      }
      writer.close();
    }
  }

  private static void scale(RealMatrix m) {
    double[] sds = new double[m.getColumnDimension()];
    double[] mean = new double[m.getColumnDimension()];

    for (int j = 0; j < mean.length; j++) {
      double[] tmp = new double[m.getRowDimension()];
      for (int i = 0; i < m.getRowDimension(); i++) {
        tmp[i] += m.getEntry(i, j);
      }
      mean[j] = ArrayUtils.mean(tmp);
      sds[j] = ArrayUtils.stdev(tmp);
    }
    for (int i = 0; i < m.getRowDimension(); i++) {
      for (int j = 0; j < mean.length; j++) {
        double standard = m.getEntry(i, j) - mean[j];
        standard /= sds[j];
        m.setEntry(i, j, standard);
      }
    }
  }

  public static void main(String[] args) {
    String outDir = args[1];
    new File(outDir).mkdirs();
    Logger log = new Logger(outDir + "log.log");
    String matFile = outDir + "mat.ser.gz";
    String svdFile = ext.addToRoot(matFile, ".svd");
    String pcFile = outDir + "pcs.gz";

    log.reportTimeInfo("Operating in " + outDir);

    if (!Files.exists(matFile)) {
      String[] files = Files.listFullPaths(args[0], ".gz");
      log.reportTimeInfo("Found " + files.length + " files to PCA");
      String[] availableRows = HashVec.loadFileToStringArray(files[0], false, new int[] {0, 1, 2},
                                                             false);
      List<String> rows = new ArrayList<>();
      List<Integer> indices = new ArrayList<>();
      for (int i = 0; i < availableRows.length; i++) {
        String[] split = availableRows[i].trim().split("\t");

        Segment seg = new Segment(split[0], split[1], split[2]);
        //        or other filter logic
        if (seg.getChr() > 0 && seg.getChr() < 3) {
          rows.add(seg.getUCSClocation());
          indices.add(i);
        }
      }
      log.reportTimeInfo("Preparing matrix with " + rows.size() + " rows and " + files.length
                         + "columns");

      RealMatrix m = MatrixUtils.createRealMatrix(rows.size(), files.length);
      String[] colNames = new String[files.length];
      String[] rowNames = ArrayUtils.toStringArray(rows);
      for (int i = 0; i < colNames.length; i++) {
        colNames[i] = ext.rootOf(files[i]);
        log.reportTime("Parsing sample " + i + ", " + colNames[i]);
        String[] data = HashVec.loadFileToStringArray(files[i], false, new int[] {3}, false);
        for (int j = 0; j < indices.size(); j++) {
          m.addToEntry(j, i, Double.parseDouble(data[indices.get(j)]));
        }
      }
      //      new SingularValueDecomposition(
      log.reportTime("Scaling data");

      scale(m);
      log.reportTime("Saving progress");

      SerializedFiles.writeSerial(new TOPMedPCAMatrix(m, colNames, rowNames), matFile, true);
    }
    if (!Files.exists(svdFile)) {
      log.reportTime("Loading data");
      TOPMedPCAMatrix tm = (TOPMedPCAMatrix) SerializedFiles.readSerial(matFile, log, false);
      log.reportTimeInfo("Rows=" + tm.m.getRowDimension());
      log.reportTimeInfo("Cols=" + tm.m.getColumnDimension());
      log.reportTime("Computing SVD");
      tm.computeSVD(log);

      log.reportTime("Finished Computing SVD");

      log.reportTime("Saving progress");

      SerializedFiles.writeSerial(tm, svdFile, true);
    }
    log.reportTimeInfo("Loading " + svdFile);
    TOPMedPCAMatrix tm = (TOPMedPCAMatrix) SerializedFiles.readSerial(svdFile, log, false);

    tm.dumpPCsToText(pcFile);

    DenseMatrix64F a = new DenseMatrix64F(1, 1);
    a.reshape(tm.m.getRowDimension(), tm.m.getColumnDimension());
    for (int row = 0; row < tm.m.getRowDimension(); row++) {
      for (int col = 0; col < tm.m.getColumnDimension(); col++) {
        a.add(row, col, tm.m.getEntry(row, col));
      }
    }
    log.reportTimeInfo("Computing EJML PCs");
    SvdImplicitQrDecompose_D64 svd = new SvdImplicitQrDecompose_D64(false, false, true, false);
    svd.decompose(a);
    log.reportTimeInfo("Finished Computing EJML PCs");

  }
}
