package org.genvisis.expression;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.DoubleVector;
import org.genvisis.common.HashVec;
import org.genvisis.common.ext;

import com.google.common.primitives.Doubles;

public class RelativeQuantification {
  public static final String[] HEADER = {"Plate", "Well ID", "Well", "Sample", "Detector", "Task",
      "Ct", "Ct Std Err", "Avg Ct", "Avg dCt", "dCt Std Err", "ddCt", "RQ", "RQ Min", "RQ Max",
      "Omit", "Filtered", "Threshold", "Auto Ct", "Baseline", "Start", "End"};

  public static void main(String[] args) {
    int numArgs = args.length;
    String dir =
        "C:\\Documents and Settings\\npankrat\\My Documents\\Expression\\RelativeQuantification\\";
    String subdir = "source/";
    String filename = "RelativeQuantification.dat";

    String usage = "\\n" + "park.expression.RelativeQuantification requires 0-1 arguments\n"
        + "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("file=")) {
        filename = arg.split("=")[1];
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      parse(dir, subdir);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  // public static void parse(String dir, String subdir) {
  // BufferedReader reader;
  // PrintWriter writer, summary;
  // String[] line, indKeys, regionKeys, probeKeys;
  // Hashtable<String, Hashtable<String, Hashtable<String, DoubleVector>>>
  // probes = new Hashtable<String,
  // Hashtable<String,Hashtable<String,DoubleVector>>>();
  // Hashtable<String, Hashtable<String, DoubleVector>> regions;
  // Hashtable<String, DoubleVector> inds;
  // DoubleVector values, regionMeans, regionStdevs, allObs, allObsAroundZero;
  // IntVector regionCounts;
  // String ind, region, probeName;
  // Hashtable<String, String> lookup = new Hashtable<String, String>();
  // double mean, stdev, z, t;
  //
  // File[] files = new File(dir+subdir).listFiles(new FilenameFilter() {
  // public boolean accept(File file, String filename) {
  // return filename.endsWith(".csv");
  // }
  // });
  //
  // for (int i = 0; i < files.length; i++) {
  // try {
  // reader = new BufferedReader(new FileReader(files[i]));
  // ext.checkHeader(reader.readLine().trim().split(","), HEADER, true);
  // while (reader.ready()) {
  // line = reader.readLine().trim().split(",");
  // ind = line[3].substring(0, line[3].length()-2);
  // region = line[3].substring(line[3].length()-2, line[3].length()-1);
  // probeName = line[4];
  // if (!lookup.containsKey(probeName)) {
  // lookup.put(probeName, ext.rootOf(files[i].getName()));
  // }
  // if (probes.containsKey(probeName)) {
  // regions = probes.get(probeName);
  // } else {
  // probes.put(probeName, regions = new Hashtable<String,
  // Hashtable<String,DoubleVector>>());
  // }
  // if (regions.containsKey(region)) {
  // inds = regions.get(region);
  // } else {
  // regions.put(region, inds = new Hashtable<String, DoubleVector>());
  // }
  // if (inds.containsKey(ind)) {
  // values = inds.get(ind);
  // } else {
  // inds.put(ind, values = new DoubleVector());
  // }
  // if (line.length > 2 && !line[16].equals("Outlier")){
  // values.add(Double.parseDouble(line[6]));
  // }
  // }
  // reader.close();
  // } catch (FileNotFoundException fnfe) {
  // System.err.println("Error: file \""+files[i].getName()+"\" not found in
  // current directory");
  // System.exit(1);
  // } catch (IOException ioe) {
  // System.err.println("Error reading file \""+files[i].getName()+"\"");
  // System.exit(2);
  // }
  // }
  //
  //
  // try {
  // writer = new PrintWriter(new FileWriter(dir+"RQ.xln"));
  // writer.println("Gene\tSNP\tRegion\tIndividual\tN\tMean\tStdev\tStderr\tRegionMean\tRegionStdev\tZ-score\tpval
  // from Zdist\tT\tpval from Tdist");
  // summary = new PrintWriter(new FileWriter(dir+"allelotype_summary.xln"));
  // summary.println("gene\tsnp\tN\tproportion\tstd dev\tt-value df\tregion
  // A\tregion B\tregion C");
  // probeKeys = HashVec.getKeys(probes);
  // for (int i = 0; i < probeKeys.length; i++) {
  // regions = probes.get(probeKeys[i]);
  // regionKeys = HashVec.getKeys(regions);
  // for (int j = 0; j < regionKeys.length; j++) {
  // inds = regions.get(regionKeys[j]);
  // indKeys = HashVec.getKeys(inds);
  // regionMeans = new DoubleVector();
  // regionStdevs = new DoubleVector();
  // regionCounts = new IntVector();
  // allObs = new DoubleVector();
  // allObsAroundZero = new DoubleVector();
  // for (int k = 0; k < indKeys.length; k++) {
  // values = inds.get(indKeys[k]);
  // regionMeans.add(Array.mean(values.toArray()));
  // regionStdevs.add(Array.stdev(values.toArray()));
  // regionCounts.add(values.size());
  // for (int l = 0; l < values.size(); l++) {
  // allObs.add(values.elementAt(l));
  // allObsAroundZero.add(values.elementAt(l)-0.5);
  // }
  // writer.println((j==0 &&
  // k==0?lookup.get(probeKeys[i])+"\t"+probeKeys[i]:"\t")+"\t"+(k==0?regionKeys[j]:"")+"\t"+indKeys[k]+"\t"+values.size()+"\t"+ext.formDeci(Array.mean(values.toArray()),
  // 3)+"\t"+ext.formDeci(Array.stdev(values.toArray()),
  // 3)+"\t"+ext.formDeci(Array.stdev(values.toArray())/Math.sqrt(values.size()),
  // 3));
  // //
  // writer.println(lookup.get(snpKeys[i])+"\t"+snpKeys[i]+"\t"+regionKeys[j]+"\t"+indKeys[k]+"\t"+values.size()+"\t"+ext.formDeci(Array.mean(values.toArray()),
  // 3)+"\t"+ext.formDeci(Array.stdev(values.toArray()), 3));
  // }
  // mean = Array.mean(regionMeans.toArray());
  // stdev = Array.stdev(regionMeans.toArray());
  // z = Stats.ztest(mean, stdev, 0.5);
  // t = Stats.ttestOneSample(mean, stdev, regionMeans.size(), 0.5);
  // writer.println("\t\t\t\t\t\t\t\t"+mean+"\t"+stdev+"\t"+z+"\t"+ProbDist.NormDist(z)+"\t"+t+"\t"+ProbDist.TDist(t,
  // indKeys.length-1)+"\tOne observation (mean) per sample");
  // if (j==0) {
  // summary.print(lookup.get(probeKeys[i])+"\t"+probeKeys[i]+"\t"+indKeys.length+"\t"+mean+"\t"+stdev+"\t"+t+"\t"+(indKeys.length-1));
  // }
  // summary.print("\t"+ext.prettyP(ProbDist.TDist(t, indKeys.length-1)));
  // stdev = Array.stdev(allObs.toArray());
  // z = Stats.ztest(mean, stdev, 0.5);
  // t = Stats.ttestOneSample(mean, stdev, regionMeans.size(), 0.5);
  // writer.println("\t\t\t\t\t\t\t\t"+mean+"\t"+stdev+"\t"+z+"\t"+ProbDist.NormDist(z)+"\t"+t+"\t"+ProbDist.TDist(t,
  // indKeys.length-1)+"\tAll observations");
  // stdev = pool(regionStdevs, regionCounts);
  // z = Stats.ztest(mean, stdev, 0.5);
  // t = Stats.ttestOneSample(mean, stdev, regionMeans.size(), 0.5);
  // writer.println("\t\t\t\t\t\t\t\t"+mean+"\t"+stdev+"\t"+z+"\t"+ProbDist.NormDist(z)+"\t"+t+"\t"+ProbDist.TDist(t,
  // indKeys.length-1)+"\tWith pooled variance estimator");
  // }
  // summary.println();
  // }
  // writer.close();
  // summary.close();
  // } catch (Exception e) {
  // System.err.println("Error writing to allelptypes.out in "+dir);
  // e.printStackTrace();
  // }
  // }
  //
  // // pooled two sample variance estimator
  // public static double pool(DoubleVector stdevs, IntVector counts) {
  // double sumVars = 0;
  // int totalN = 0;
  //
  // for (int i = 0; i < stdevs.size(); i++) {
  // sumVars += stdevs.elementAt(i)*stdevs.elementAt(i)*(counts.elementAt(i) -
  // 1);
  // totalN += counts.elementAt(i) - 1;
  // }
  //
  // return Math.sqrt(sumVars/totalN);
  // }

  public static void parse(String dir, String subdir) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line, inds, regionNames, probeNames, values;
    Hashtable<String, Hashtable<String, Hashtable<String, String[]>>> individuals =
        new Hashtable<String, Hashtable<String, Hashtable<String, String[]>>>();
    Hashtable<String, Hashtable<String, String[]>> probes;
    Hashtable<String, String[]> regions;
    double[][][] means, stderrs;
    DoubleVector sampleStderrs;
    String ind, region, probeName;
    Hashtable<String, String> probeGeneLookup = new Hashtable<String, String>();
    Vector<String> vProbes = new Vector<String>();
    Vector<String> vRegions = new Vector<String>();

    File[] files = new File(dir + subdir).listFiles(new FilenameFilter() {
      @Override
      public boolean accept(File file, String filename) {
        return filename.endsWith(".csv");
      }
    });

    for (File file : files) {
      try {
        reader = new BufferedReader(new FileReader(file));
        ext.checkHeader(reader.readLine().trim().split(","), HEADER, true);
        while (reader.ready()) {
          line = reader.readLine().trim().split(",");
          ind = line[3].substring(0, line[3].length() - 2);
          region = line[3].substring(line[3].length() - 2, line[3].length() - 1).toUpperCase();
          probeName = line[4];
          HashVec.addIfAbsent(region, vRegions);
          HashVec.addIfAbsent(probeName, vProbes);
          if (!probeGeneLookup.containsKey(probeName)) {
            probeGeneLookup.put(probeName, ext.rootOf(file.getName()));
          }
          if (individuals.containsKey(ind)) {
            probes = individuals.get(ind);
          } else {
            individuals.put(ind, probes = new Hashtable<String, Hashtable<String, String[]>>());
          }
          if (probes.containsKey(probeName)) {
            regions = probes.get(probeName);
          } else {
            probes.put(probeName, regions = new Hashtable<String, String[]>());
          }
          if (regions.containsKey(region)) {
            values = regions.get(region);
            if (!values[0].equals(line[8]) || !values[1].equals(line[7])
                || !values[2].equals(line[12]) || !values[3].equals(line[13])
                || !values[4].equals(line[14])) {
              System.err.println("Error - mismatch with " + ind);
            }
          } else {
            // AvgCt, CtStderr, RQ, RQ Min, RQ Max
            regions.put(region,
                values = new String[] {line[8], line[7], line[12], line[13], line[14]});
          }
        }
        reader.close();
      } catch (FileNotFoundException fnfe) {
        System.err.println("Error: file \"" + file.getName() + "\" not found in current directory");
        System.exit(1);
      } catch (IOException ioe) {
        System.err.println("Error reading file \"" + file.getName() + "\"");
        System.exit(2);
      }
    }

    try {
      inds = HashVec.getKeys(individuals);
      probeNames = Array.toStringArray(vProbes);
      regionNames = Array.toStringArray(vRegions);

      writer = new PrintWriter(new FileWriter(dir + "RQ.xln"));
      writer.print("Sample");
      for (String probeName2 : probeNames) {
        for (String regionName : regionNames) {
          writer.print("\t" + regionName + "_" + probeGeneLookup.get(probeName2) + "-" + probeName2
              + "_AvgCt\t" + regionName + "_" + probeGeneLookup.get(probeName2) + "-" + probeName2
              + "_CtStderr\t" + regionName + "_" + probeGeneLookup.get(probeName2) + "-"
              + probeName2 + "_RQ\t" + regionName + "_" + probeGeneLookup.get(probeName2) + "-"
              + probeName2 + "_RQ Min\t" + regionName + "_" + probeGeneLookup.get(probeName2) + "-"
              + probeName2 + "_RQ Max");
        }
      }
      writer.println("\tMeanSampleStderr");
      means = new double[inds.length][probeNames.length][regionNames.length];
      stderrs = new double[inds.length][probeNames.length][regionNames.length];
      for (int i = 0; i < inds.length; i++) {
        probes = individuals.get(inds[i]);
        sampleStderrs = new DoubleVector();
        writer.print(inds[i]);
        for (int j = 0; j < probeNames.length; j++) {
          regions = probes.get(probeNames[j]);
          for (int k = 0; k < regionNames.length; k++) {
            values = regions.get(regionNames[k]);
            if (values == null) {
              writer.print("\t" + Array.toStr(Array.stringArray(5, ".")));
            } else {
              means[i][j][k] = Double.parseDouble(values[0]);
              stderrs[i][j][k] = Double.parseDouble(values[1]);
              sampleStderrs.add(Double.parseDouble(values[1]));
              writer.print("\t" + Array.toStr(values));
            }
          }
        }
        writer.println("\t" + Array.mean(Doubles.toArray(sampleStderrs)));
      }
      writer.close();

      writer = new PrintWriter(new FileWriter(dir + "RQ_db.xln"));
      writer.print("Sample");
      for (String probeName2 : probeNames) {
        for (String regionName : regionNames) {
          writer.print("\t" + regionName + "_" + probeGeneLookup.get(probeName2) + "-" + probeName2
              + "_RQ_log");
        }
      }
      writer.println("\tMeanSampleStderr");
      means = new double[inds.length][probeNames.length][regionNames.length];
      stderrs = new double[inds.length][probeNames.length][regionNames.length];
      for (String ind2 : inds) {
        probes = individuals.get(ind2);
        sampleStderrs = new DoubleVector();
        writer.print(ind2);
        for (String probeName2 : probeNames) {
          regions = probes.get(probeName2);
          for (String regionName : regionNames) {
            values = regions.get(regionName);
            if (values == null) {
              writer.print("\t.");
            } else {
              writer.print("\t" + Math.log(Double.parseDouble(values[2])) / Math.log(2));
            }
          }
        }
        writer.println("\t" + Array.mean(Doubles.toArray(sampleStderrs)));
      }
      writer.close();

      for (String probeName2 : probeNames) {
        writer = new PrintWriter(
            new FileWriter(dir + probeGeneLookup.get(probeName2) + "_" + probeName2 + ".xln"));
        writer.print("Sample");
        for (String regionName : regionNames) {
          writer.print(
              "\t" + regionName + "_" + "_AvgCt\t" + regionName + "_" + "_CtStderr\t" + regionName
                  + "_" + "_RQ\t" + regionName + "_" + "_RQ Min\t" + regionName + "_" + "_RQ Max");
        }
        writer.println();
        for (String ind2 : inds) {
          regions = individuals.get(ind2).get(probeName2);
          writer.print(ind2);
          for (String regionName : regionNames) {
            values = regions.get(regionName);
            if (values == null) {
              writer.print("\t" + Array.toStr(Array.stringArray(5, ".")));
            } else {
              values[3] = (Double.parseDouble(values[2]) - Double.parseDouble(values[3])) + "";
              values[4] = (Double.parseDouble(values[4]) - Double.parseDouble(values[2])) + "";
              writer.print("\t" + Array.toStr(values));
            }
          }
          writer.println();
        }
        writer.close();
      }

    } catch (Exception e) {
      System.err.println("Error writing to RQ.xln in " + dir);
      e.printStackTrace();
    }
  }

}
