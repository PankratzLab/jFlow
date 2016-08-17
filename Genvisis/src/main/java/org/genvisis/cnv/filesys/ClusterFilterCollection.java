package org.genvisis.cnv.filesys;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Map.Entry;
import java.util.Vector;

import javax.swing.JOptionPane;

import org.genvisis.cnv.manage.TextExport;
import org.genvisis.cnv.plots.GenericRectangle;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.ext;

/**
 * This is data structure to hold a group of filters that can be used to screen data points.
 * 
 * @author npankratz and zxu
 *
 */
public class ClusterFilterCollection implements Serializable, TextExport {
  private static final long serialVersionUID = 1L;

  public static void describe(String filename) {
    PrintWriter writer;
    String[] markerNames;
    ClusterFilterCollection trav;
    int count;

    trav = load(filename, false);

    count = 0;
    try {
      writer = new PrintWriter(new FileWriter(ext.rootOf(filename) + "_described.xln"));
      writer.println("MarkerName\t#filters");
      markerNames = trav.getMarkerNames();
      for (String markerName : markerNames) {
        writer.println(markerName + "\t" + trav.getClusterFilters(markerName).size());
        count += trav.getClusterFilters(markerName).size();
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + ext.rootOf(filename) + "_described.xln");
      e.printStackTrace();
    }

    System.out.println(
        "There were a total of " + count + " clusterFilters across " + trav.getSize() + " markers");
  }

  public static void dump(String filename) {
    ClusterFilterCollection collection;
    String outFile;

    if (!Files.exists(filename)) {
      System.err.println("Error - collection '" + filename + "' does not exist");
      return;
    }

    try {
      collection = load(filename, false);
      outFile = ext.rootOf(filename, false) + "_dump.xln";
      collection.exportToText(null, outFile);
    } catch (Exception e) {
      System.err.println("Error writing to " + ext.rootOf(filename) + "_dump.xln");
      e.printStackTrace();
    }
  }

  public static String getClusterFilterFilenameSelection(Project proj) {
    String result;
    result = (String) JOptionPane.showInputDialog(null, "Please select a cluster filter file:",
        "Apply Cluster Filters", JOptionPane.QUESTION_MESSAGE, null,
        Array.addStrToArray("(--Do not apply any cluster filter--)",
            Files.list(proj.DATA_DIRECTORY.getValue(false, true), null,
                ext.removeDirectoryInfo(proj.getProperty(proj.CLUSTER_FILTER_COLLECTION_FILENAME)),
                false, proj.JAR_STATUS.getValue())),
        proj.getProperty(proj.CLUSTER_FILTER_COLLECTION_FILENAME));
    if (result == null) {
      result = "cancel";
    } else if (result.equals("(--Do not apply any cluster filter--)")) {
      result = null;
    }
    return result;
  }

  public static String getGenotypeLookupTableSelection(Project proj) {
    String result;
    result = (String) JOptionPane.showInputDialog(null, "Please select a AB genotype lookup table:",
        "Select AB Genotype Lookup Table", JOptionPane.QUESTION_MESSAGE, null,
        new String[] {"Lookup Table 1", "Lookup Table 2", "Lookup Table 3"},
        proj.getProperty(proj.CLUSTER_FILTER_COLLECTION_FILENAME));
    return result;
  }

  public static ClusterFilterCollection load(String filename, boolean jar) {
    return (ClusterFilterCollection) SerializedFiles.readSerial(filename, jar, true);
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String[] filenames = null;
    String out = "merged.ser";
    String filename = "data/clusterFilters.ser";
    boolean describeFile = false;
    boolean exportFile = false;
    boolean importFile = false;
    boolean recodeOld = false;

    String usage = "\n" + "cnv.filesys.ClusterFilterCollection requires 0-1 arguments\n"
        + "   (1) names of clusterFilter files to merge (i.e. files=file1.ser,file2.ser,file3.ser (not the default))\n"
        + "   (2) output filename (i.e. out=" + out + " (default))\n" + " OR:\n"
        + "   (1) describe clusterFilter file (i.e. -describe (not the default))\n"
        + "   (2) name of clusterFilter file to describe (i.e. file=" + filename + " (default))\n"
        + " OR:\n" + "   (1) export clusterFilter file (i.e. -export (not the default))\n"
        + "   (2) name of clusterFilter file to import (i.e. file=" + filename + " (default))\n"
        + " OR:\n"
        + "   (1) import text file into a clusterFilter collection (i.e. -import (not the default; and still needs to be implemented))\n"
        + "   (2) name of text file to import (i.e. file=" + filename + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("file=")) {
        filename = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("files=")) {
        filenames = arg.split("=")[1].split(",");
        numArgs--;
      } else if (arg.startsWith("out=")) {
        out = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("-describe")) {
        describeFile = true;
        numArgs--;
      } else if (arg.startsWith("-export")) {
        exportFile = true;
        numArgs--;
      } else if (arg.startsWith("-import")) {
        importFile = true;
        numArgs--;
      } else if (arg.startsWith("-recode")) {
        recodeOld = true;
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
      if (recodeOld) {
        recode(filename, out);
      } else if (filenames != null) {
        merge(filenames, out);
      } else if (describeFile) {
        describe(filename);
      } else if (exportFile) {
        dump(filename);
      } else if (importFile) {
        // TODO still needs to be implemented
        System.err.println("import still needs to be implemented");
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public static void merge(String[] filenames, String outfile) {
    String[] markerNames;
    ClusterFilterCollection master, trav;
    Vector<String> v = new Vector<String>();
    ArrayList<ClusterFilter> masterArray, travArray;

    master = new ClusterFilterCollection();
    for (String filename : filenames) {
      trav = load(filename, false);
      markerNames = trav.getMarkerNames();
      for (String markerName : markerNames) {
        masterArray = master.getClusterFilters(markerName);
        travArray = trav.getClusterFilters(markerName);
        if (masterArray != null) {
          HashVec.addIfAbsent(markerName, v);
        }
        for (int k = 0; k < travArray.size(); k++) {
          master.addClusterFilter(markerName, travArray.get(k));
        }
      }
    }

    master.serialize(outfile);

    if (v.size() > 0) {
      System.out.println("The following markers had cluster filters in multiple files:");
      System.out.println(Array.toStr(Array.toStringArray(v), "\n"));
    }
  }

  // filterIndividual, worry about this later!!

  private static void recode(String filename, String out) {
    ClusterFilterCollection collection;

    if (!Files.exists(filename)) {
      System.err.println("Error - collection '" + filename + "' does not exist");
      return;
    }

    collection = load(filename, false);

    for (Entry<String, ArrayList<ClusterFilter>> entry : collection.hash.entrySet()) {
      for (ClusterFilter cf : entry.getValue()) {
        cf.setPlotType((byte) (cf.getPlotType() - 1));
      }
    }

    collection.serialize(out);
  }

  /**
   * The structure of "hash" is Hashtable<String markerName, ArrayList<ClusterFilter>
   * clusterFilter>.
   */
  private final Hashtable<String, ArrayList<ClusterFilter>> hash;

  // constructor, get(String markerName), load, serialize
  public ClusterFilterCollection() {
    hash = new Hashtable<String, ArrayList<ClusterFilter>>();
  }

  public void addClusterFilter(String markerName, ClusterFilter filter) {
    // ArrayList<ClusterFilter> filters;
    //
    // if (hash.containsKey(markerName)) {
    // filters = hash.get(markerName);
    // } else {
    // hash.put(markerName, filters = new ArrayList<ClusterFilter>());
    // }
    // filters.add(filter);

    if (!hash.containsKey(markerName)) {
      hash.put(markerName, new ArrayList<ClusterFilter>());
    }
    hash.get(markerName).add(filter);
  }

  public void deleteClusterFilter(String markerName, byte index) {
    if (hash.containsKey(markerName) && hash.get(markerName).size() > index) {
      hash.get(markerName).remove(index);
      if (hash.get(markerName).size() == 0) {
        hash.remove(markerName);
      }
    } else {
      System.err.println(
          "Error deleting the cluster filter: either no cluster filter associate with this marker name, or the index for the cluster filter to be deleted does not exist.");
    }
  }

  @Override
  public void exportToText(Project proj, String outputFile) {
    PrintWriter writer;
    String[] markerNames;
    ArrayList<ClusterFilter> list;
    ClusterFilter filter;

    try {
      writer = new PrintWriter(new FileWriter(outputFile));
      writer.println(
          "MarkerIndex\tMarkerName\tFilterIndex\tPlotType\tGenotype\tminX\tminY\tmaxX\tmaxY");
      markerNames = getMarkerNames();
      for (int i = 0; i < markerNames.length; i++) {
        list = getClusterFilters(markerNames[i]);
        for (int j = 0; j < list.size(); j++) {
          filter = list.get(j);
          writer.println(i + "\t" + markerNames[i] + "\t" + j + "\t" + filter.getPlotType() + "\t"
              + filter.getCluterGenotype() + "\t" + filter.getXMin() + "\t" + filter.getYMin()
              + "\t" + filter.getXMax() + "\t" + filter.getYMax());
        }
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + outputFile + "_dump.xln");
      e.printStackTrace();
    }
  }

  public byte[] filterMarker(MarkerData markerData, float gcThreshold) {
    byte[] result, original;
    float[] realX;
    float[] realY;
    ArrayList<ClusterFilter> clusterFilters;

    original = markerData.getAbGenotypes();
    clusterFilters = hash.get(markerData.getMarkerName());
    result = new byte[original.length];
    for (int j = 0; j < original.length; j++) {
      if (markerData.getGCs()[j] < gcThreshold) {
        result[j] = (byte) -1;
      } else {
        result[j] = original[j];
      }
    }

    for (int i = 0; clusterFilters != null && i < clusterFilters.size(); i++) {
      // get appropriate data (X/Y Theta/R LRR/BAF)
      switch (clusterFilters.get(i).getPlotType()) {
        case 0:
          // realX = markerData.getX_Raws();
          // realY = markerData.getY_Raws();
          // break;
          // case 1:
          realX = markerData.getXs();
          realY = markerData.getYs();
          break;
        case 1:
          realX = markerData.getThetas();
          realY = markerData.getRs();
          break;
        case 2:
          realX = markerData.getBAFs();
          realY = markerData.getLRRs();
          break;
        default:
          realX = markerData.getXs();
          realY = markerData.getYs();
      }
      // iterate through all samples
      for (int j = 0; j < markerData.getAbGenotypes().length; j++) {
        if (realX[j] >= clusterFilters.get(i).getXMin()
            && realY[j] >= clusterFilters.get(i).getYMin()
            && realX[j] <= clusterFilters.get(i).getXMax()
            && realY[j] <= clusterFilters.get(i).getYMax()) {
          result[j] = clusterFilters.get(i).getCluterGenotype();
          if (result[j] < -1) {
            System.err.println("Error - result[" + j + "]=" + result[j]);
          }
        }
      }
    }

    return result;
  }

  public ArrayList<ClusterFilter> getClusterFilters(String markerName) {
    return hash.get(markerName);
  }

  // ??? How to select the last filter???
  public byte getGenotype(String markerName, byte index) {
    if (hash.containsKey(markerName) && hash.get(markerName).size() > index) {
      return hash.get(markerName).get(index).getCluterGenotype();
    } else {
      System.err.println("Error - Trying to get a ClusterFilter that does not exist.");
      return (byte) -1;
    }
  }

  public String[] getMarkerNames() {
    return HashVec.getKeys(hash);
  }

  public GenericRectangle[] getRectangles(String markerName, byte plotType, byte thickness,
      boolean fill, boolean roundedCorners, byte color, byte layer) {
    ArrayList<ClusterFilter> clusterFilters;
    GenericRectangle[] result = new GenericRectangle[getSize(markerName)];
    // ArrayList<GenericRectangle> rectangles = new ArrayList<GenericRectangle>();

    clusterFilters = hash.get(markerName);
    for (int i = 0; clusterFilters != null && i < clusterFilters.size(); i++) {
      if (clusterFilters.get(i).getPlotType() == plotType) {
        result[i] = new GenericRectangle(clusterFilters.get(i).getXMin(),
            clusterFilters.get(i).getYMin(), clusterFilters.get(i).getXMax(),
            clusterFilters.get(i).getYMax(), thickness, fill, roundedCorners, color, layer, true);
      } else {
        result[i] =
            new GenericRectangle(clusterFilters.get(i).getXMin(), clusterFilters.get(i).getYMin(),
                clusterFilters.get(i).getXMax(), clusterFilters.get(i).getYMax(), thickness, fill,
                roundedCorners, color, (byte) -1, true);
      }
    }
    return result;
  }

  public int getSize() {
    return (hash == null ? 0 : hash.size());
  }

  public byte getSize(String markerName) {
    return (byte) (hash.get(markerName) == null ? 0 : hash.get(markerName).size());// ???
  }

  public void serialize(String filename) {
    SerializedFiles.writeSerial(this, filename);
  }

  public void updateGenotype(String markerName, byte index, byte newGenotype) {
    ArrayList<ClusterFilter> temp = hash.get(markerName);
    if (temp == null) {
      return;
    }
    temp.get(index).setClusterGenotype(newGenotype);
  }

}
