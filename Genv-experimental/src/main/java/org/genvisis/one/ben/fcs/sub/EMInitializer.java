package org.genvisis.one.ben.fcs.sub;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.distribution.MixtureMultivariateNormalDistribution;
import org.apache.commons.math3.distribution.MultivariateNormalDistribution;
import org.apache.commons.math3.distribution.fitting.MultivariateNormalMixtureExpectationMaximization;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.Pair;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

import edu.stanford.facs.logicle.Logicle;

public class EMInitializer {
  
  private static double count(double max) {
    double dec = max;
    double decCnt = 0;
    while (dec > 1) {
      dec = dec / 10d;
      decCnt++;
    }
    return decCnt;
  }

  private static Logicle getBiexScale(double max) {
    double w = 2; // linear decades // W=1 is weird, W=3 -> "scale() didn't
                  // converge"
    double a = 0; // negative decades
    double decCnt = count(max);

    return new Logicle(max, Math.min(w, decCnt / 2), decCnt, a);
  }

  private static Logicle getLogicle(double[] data) {
    return getBiexScale(Array.max(data));
  }

  private static double[][] transpose(double[][] matrixData) {
    RealMatrix m = MatrixUtils.createRealMatrix(matrixData);
    return m.transpose().getData();
  }
  
  static final String[] DATA_COLUMNS = {
      "Comp-BB515-A (CD27)",
      "Comp-PE-CF594-A (HLA-DR)",
      "Comp-PE-Cy7-A (CD19)",
      "Comp-BUV 395-A (CD8)",
      "Comp-BUV 737-A (IgD)",
      "Comp-APC-A (CD3)",
      "Comp-BV 421-A (CCR7)",
      "Comp-BV 510-A (CD28)",
      "Comp-BV 605-A (CD95)",
      "Comp-BV 711-A (CD45RA)"
  };
  
  static final int CLUSTERS = 18;
  
  static final String DATA_DIR = "F:/Flow/controlFCS/";
  static final String DIST_DIR = "F:/Flow/controlFCS/dists/";
  static final String DATA_EXT = ".xln";
  static final Logger log = new Logger();
  
  static void run() {
    String[] files = Files.list(DATA_DIR, DATA_EXT, false);
    
    int cnt = 0;
    for (String f : files) {
      cnt += Files.countLines(DATA_DIR + f, 1);
    }
    
    System.out.println("Reading " + cnt + " lines");
    double[][] data = new double[DATA_COLUMNS.length][cnt];
    
    int index = 0;
    for (String f : files) {
      try {
        BufferedReader reader = Files.getAppropriateReader(DATA_DIR + f);
        String[] parts = reader.readLine().split("\t");
        int[] fact = ext.indexFactors(DATA_COLUMNS, parts, false, false);
        String line = null;
        while((line = reader.readLine()) != null) {
          parts = line.split("\t");
          for (int i = 0; i < DATA_COLUMNS.length; i++) {
            data[i][index] = Double.parseDouble(parts[fact[i]]);
          }
          index++;
        }
        reader.close();
      } catch (FileNotFoundException e) {
        // TODO Auto-generated catch block
        e.printStackTrace();
      } catch (IOException e) {
        // TODO Auto-generated catch block
        e.printStackTrace();
      }
    }
    
    Logicle[] scales = new Logicle[DATA_COLUMNS.length];
    for (int i = 0; i < scales.length; i++) {
      scales[i] = getLogicle(data[i]);
    }

    double[][] transformed = new double[data.length][data[0].length];
    for (int i = 0; i < scales.length; i++) {
      for (int j = 0; j < data[0].length; j++) {
        transformed[i][j] = scales[i].scale(data[i][j]);
      }
    }

    data = transpose(transformed);
    
    MixtureMultivariateNormalDistribution initialMix = MultivariateNormalMixtureExpectationMaximization.estimate(data, CLUSTERS);
    
    List<Pair<Double, MultivariateNormalDistribution>> init = initialMix.getComponents();
    int len = ("" + init.size()).length() + 1;
    
    log.report("Writing output");
    
    for (int i = 0, count = init.size(); i < count; i++) {
      String filename = i + "";
      while (filename.length() < len) {
        filename = "0" + filename;
      }
      
      Pair<Double, MultivariateNormalDistribution> p = init.get(i);
      Files.write("" + p.getFirst(), DIST_DIR + filename + ".comp");
      
      MultivariateNormalDistribution dist = p.getSecond();
      double[] means = dist.getMeans();
      double[][] covars = dist.getCovariances().getData();
      
      PrintWriter writer2 = Files.getAppropriateWriter(DIST_DIR + filename + ".comp" + ".means");
      for (double d : means) {
        writer2.println("" + d);
      }
      writer2.flush();
      writer2.close();

      writer2 = Files.getAppropriateWriter(DIST_DIR + filename + ".comp" + ".covars");
      for (double[] d : covars) {
        writer2.println("" + Array.toStr(d, "\t"));
      }
      writer2.flush();
      writer2.close();
      
    }
    
  }
  
  public static MixtureMultivariateNormalDistribution load(String dir) {
    String[] compFiles = Files.list(dir, ".comp", false);
    
    List<Pair<Double, MultivariateNormalDistribution>> ret = new ArrayList<Pair<Double, MultivariateNormalDistribution>>();
    double[] means;
    double[][] covars;
    
    for (int i = 0; i < compFiles.length; i++) {
      String[] cmp = HashVec.loadFileToStringArray(dir + compFiles[i], false, null, false);
      String[] mns = HashVec.loadFileToStringArray(dir + compFiles[i] + ".means", false, null, false);
      String[][] cvrs = HashVec.loadFileToStringMatrix(dir + compFiles[i] + ".covars", false, null, "\t", false, 10, false);
      
      double comps = Double.parseDouble(cmp[0]);
      means = new double[mns.length];
      for (int d = 0; d < mns.length; d++) {
        means[d] = Double.parseDouble(mns[d]);
      }
      covars = new double[cvrs.length][];
      for (int d = 0; d < cvrs.length; d++) {
        covars[d] = new double[cvrs[d].length];
        for (int d2 = 0; d2 < cvrs[d].length; d2++) {
          covars[d][d2] = Double.parseDouble(cvrs[d][d2]);
        }
      }
      
      ret.add(new Pair<Double, MultivariateNormalDistribution>(comps, new MultivariateNormalDistribution(means, covars)));
    }
    
    MixtureMultivariateNormalDistribution mmnd = new MixtureMultivariateNormalDistribution(ret);
    return mmnd;
  }
  
  
  public static void main(String[] args) {
    run();
  }
  
}
