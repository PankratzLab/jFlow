package org.genvisis.cnv.filesys;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.lang.ref.Reference;
import java.lang.ref.SoftReference;
import java.util.Arrays;
import java.util.Hashtable;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.SciStringComparator;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.ext;
import com.google.common.collect.ImmutableMap;

public class SampleList implements Serializable {

  public static final long serialVersionUID = 1L;

  private final long fingerprint;
  private final String[] samples;
  private transient Reference<ImmutableMap<String, Integer>> sampleIndicesRef = null;

  public SampleList(String[] samples) {
    this.samples = samples;
    fingerprint = MarkerSet.fingerprint(samples);
  }

  public long getFingerprint() {
    return fingerprint;
  }

  public String[] getSamples() {
    return samples;
  }

  public ImmutableMap<String, Integer> getSampleIndices() {
    ImmutableMap<String, Integer> sampleIndices = sampleIndicesRef == null ? null
                                                                           : sampleIndicesRef.get();
    if (sampleIndices == null) {
      ImmutableMap.Builder<String, Integer> sampleIndicesBuilder = ImmutableMap.builder();
      for (int i = 0; i < samples.length; i++) {
        sampleIndicesBuilder.put(samples[i], i);
      }
      sampleIndices = sampleIndicesBuilder.build();
      sampleIndicesRef = new SoftReference<ImmutableMap<String, Integer>>(sampleIndices);
    }
    return sampleIndices;
  }

  public void writeToTextFile(String filename) {
    PrintWriter writer;

    try {
      writer = Files.openAppropriateWriter(filename);
      for (String sample : samples) {
        writer.println(sample);
      }
      writer.close();
    } catch (IOException ioe) {
      System.err.println("Error writing to " + filename);
      ioe.printStackTrace();
    }
  }

  public void serialize(String filename) {
    SerializedFiles.writeSerial(this, filename);
  }

  public static SampleList load(String filename) {
    return (SampleList) SerializedFiles.readSerial(filename, true);
  }

  public static SampleList generateSampleList(Project proj) {
    String[] files, samples;
    SampleList list;
    int[] keys;
    int countAt;
    Logger log;

    log = proj.getLog();

    // if (Files.list(proj.getDir(Project.MARKER_DATA_DIRECTORY, true),
    // MarkerData.MARKER_DATA_FILE_EXTENSION, proj.getJarStatus()).length>0) {
    // System.err.println("Error - Refusing to create new SampleList until the plots directory is
    // either deleted or emptied; altering the SampleList will invalidate those files");
    // System.exit(1);
    // }

    if (Files.exists(proj.SAMPLE_DIRECTORY.getValue(false, true))) {
      files = Files.list(proj.SAMPLE_DIRECTORY.getValue(false, true), Sample.SAMPLE_FILE_EXTENSION);
    } else {
      log.reportError("Error - failed to find the SAMPLE_DIRECTORY ("
                      + proj.SAMPLE_DIRECTORY.getValue(false, true)
                      + "); no SampleList could be generated");
      return null;
    }

    countAt = 0;
    Arrays.sort(files, new SciStringComparator());
    samples = new String[files.length];
    for (int i = 0; i < samples.length; i++) {
      samples[i] = files[i].substring(0, files[i].lastIndexOf("."));
      if (samples[i].contains("@")) {
        countAt++;
      }
    }
    list = new SampleList(samples);
    if (samples.length > 0) {
      list.serialize(proj.SAMPLELIST_FILENAME.getValue(true, true));
    } else {
      log.reportError("Error - there are no samples in the samples directory; parsing must have failed, so cannot create a SampleList");
    }
    if (countAt > 0) {
      proj.getLog()
          .report("Note - " + countAt + " ("
                  + (Double.parseDouble(ext.prettyP((double) countAt / (double) samples.length))
                     * 100)
                  + "%) of your Sample IDs contain the @ symbol, which is often used when concatenating the sample's bar code. If you would like these to be stripped, then set "
                  + proj.PARSE_AT_AT_SYMBOL.getName()
                  + "=TRUE in the properties file, delete the samples/ directory and reparse the data");
    }

    return list;
  }

  public static void serializeOutliers(Project proj) {
    String[] samples = proj.getSamples();
    String sampleDir = proj.SAMPLE_DIRECTORY.getValue(false, true);

    Hashtable<String, Float> allOutliers = new Hashtable<String, Float>();

    proj.getLog().report("Extracting outliers from " + samples.length + " Sample files");
    for (String sample : samples) {
      try {
        String sampleFile = sampleDir + sample + Sample.SAMPLE_FILE_EXTENSION;
        Hashtable<String, Float> outliers = Sample.loadOutOfRangeValuesFromRandomAccessFile(sampleFile);

        // TODO duplicates?
        if (outliers != null) {
          for (java.util.Map.Entry<String, Float> entry : outliers.entrySet()) {
            String[] keyParts = entry.getKey().split("\t");
            allOutliers.put(keyParts[0] + "\t" + sample + "\t" + keyParts[1], entry.getValue());
          }
        }
        // for duplicates, loop through outliers and use allOutliers.putIfAbsent and compare
        // returned value
      } catch (Exception e) {
        proj.getLog().reportException(e);
      }

    }
    SerializedFiles.writeSerial(allOutliers, sampleDir + "outliers.ser");

  }

  public static void main(String[] args) throws IOException {
    int numArgs = args.length;
    String filename = null;
    boolean outliers = false;

    String usage = "\n" + "filesys.SampleList requires 1+ argument\n"
                   + "   (1) project properties filename (i.e. proj="
                   + org.genvisis.cnv.Launch.getDefaultDebugProjectFile(false) + " (default))\n"
                   + "   (2) (Optional) use -outliers argument to generate an 'outliers.ser' file (not the default)\n"
                   + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("proj=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("-outliers")) {
        outliers = true;
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }

    if (outliers) {
      serializeOutliers(new Project(filename));
    } else {
      try {
        generateSampleList(new Project(filename));
      } catch (Exception e) {
        e.printStackTrace();
      }
    }
  }

}
