package org.pankratzlab.shared.filesys;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.Serializable;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Vector;
import javax.swing.JProgressBar;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.Positions;
import org.pankratzlab.common.ProgressMonitor;
import org.pankratzlab.common.SerializedFiles;

public class CNVariantHash implements Serializable {

  public static final long serialVersionUID = 1L;

  public static final int CONSTRUCT_BY_IND = 1;
  public static final int CONSTRUCT_ALL = 2;

  private String filename; // The file this hash came from
  private final Hashtable<String, Hashtable<String, CNVariant[]>> hashes;

  public CNVariantHash(String filename, int structureType, Logger log) {
    Hashtable<String, Hashtable<String, Vector<CNVariant>>> vHashes;
    Hashtable<String, Vector<CNVariant>> vHash;
    Vector<CNVariant> v;
    Hashtable<String, CNVariant[]> finalHash;
    CNVariant cnv;
    String trav, temp;
    String[] inds, chrs;
    // long time;
    BufferedReader reader;
    String[] line;
    ProgressMonitor progMonitor = new ProgressMonitor(new JProgressBar(), log);
    setFilename(filename);
    String taskName = "CNV_CONVERSION";
    progMonitor.beginDeterminateTask(taskName, "Converting CNVs", Files.countLines(filename, 1),
                                     ProgressMonitor.DISPLAY_MODE.GUI_AND_CONSOLE);

    // time = new Date().getTime();
    vHashes = new Hashtable<>();
    try {
      reader = Files.getReader(filename, true, true);

      reader.mark(1000);
      temp = reader.readLine();
      temp.length();
      line = temp.trim().split(PSF.Regex.GREEDY_WHITESPACE);
      if (!line[2].toLowerCase().equals("chr") && Positions.chromosomeNumber(line[2]) != -1) {
        reader.reset();
      }
      while (reader.ready()) {
        temp = reader.readLine();
        temp.length();
        progMonitor.updateTask(taskName);
        cnv = new CNVariant(temp.trim().split(PSF.Regex.GREEDY_WHITESPACE));
        trav = cnv.getFamilyID() + "\t" + cnv.getIndividualID();
        if (structureType == CONSTRUCT_ALL) {
          // cnv.setFamilyID(null); // CompPanel will still need to link to the proper IDs
          // cnv.setIndividualID(null);
          trav = "all";
        }

        if (vHashes.containsKey(trav)) {
          vHash = vHashes.get(trav);
        } else {
          vHashes.put(trav, vHash = new Hashtable<>());
        }
        if (vHash.containsKey(cnv.getChr() + "")) {
          v = vHash.get(cnv.getChr() + "");
        } else {
          vHash.put(cnv.getChr() + "", v = new Vector<>());
        }
        v.add(cnv);

      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + filename + "\" not found in current directory");
      fnfe.printStackTrace();
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + filename + "\"");
      ioe.printStackTrace();
    }
    progMonitor.endTask(taskName);

    hashes = new Hashtable<>();
    // time = new Date().getTime();

    inds = HashVec.getKeys(vHashes);
    taskName = "CNV_SERIALIZATION";
    progMonitor.beginDeterminateTask(taskName, "Serializing CNVs", inds.length,
                                     ProgressMonitor.DISPLAY_MODE.GUI_AND_CONSOLE);
    for (String ind : inds) {
      progMonitor.updateTask(taskName);
      vHash = vHashes.get(ind);
      finalHash = new Hashtable<>();
      chrs = HashVec.getKeys(vHash);
      for (String chr : chrs) {
        finalHash.put(chr, CNVariant.toCNVariantArray(vHash.get(chr)));
      }
      hashes.put(ind, finalHash);
    }
    progMonitor.endTask(taskName);
  }

  /**
   * @return the filename
   */
  public String getFilename() {
    return filename;
  }

  /**
   * @param filename the filename to set
   */
  public void setFilename(String filename) {
    this.filename = filename;
  }

  public Hashtable<String, CNVariant[]> getDataFor(String famID_indID) {
    if (hashes.containsKey(famID_indID)) {
      return hashes.get(famID_indID);
    } else {
      return new Hashtable<>();
    }
  }

  public void serialize(String filename) {
    SerializedFiles.writeSerial(this, filename, filename.endsWith(".gz"));
  }

  public static CNVariantHash load(String filename, int structureType, Logger log) {
    CNVariantHash hashes = null;
    String suffix;
    String suffixGz;

    if (structureType == CONSTRUCT_BY_IND) {
      suffix = ".ind.ser";
    } else if (structureType == CONSTRUCT_ALL) {
      suffix = ".all.ser";
    } else {
      System.err.println("Error - invalid CONSTRUCT type");
      suffix = null;
    }
    suffixGz = suffix + ".gz";

    boolean parse = Files.exists(filename + suffix);
    boolean parseGz = Files.exists(filename + suffixGz);

    if (parseGz && !parse) {
      hashes = (CNVariantHash) SerializedFiles.readSerial(filename + suffixGz, log, false);
    } else if (parse && !parseGz) {
      log.reportTime("Found unzipped CNVariantHash file, converting to zipped file.");
      hashes = (CNVariantHash) SerializedFiles.readSerial(filename + suffix, log, false);
      if (hashes == null) {
        log.reportTime("Couldn't read CNVariantHash file " + filename + suffix
                       + " - recreating...");
        hashes = new CNVariantHash(filename, structureType, log);
      }
      hashes.serialize(filename + suffixGz);
      // check that we've serialized successfully:
      hashes = (CNVariantHash) SerializedFiles.readSerial(filename + suffixGz, log, false);
      if (hashes != null) {
        new File(filename + suffix).delete();
      }
    } else if (parse && parseGz) {
      log.reportTime("Found both zipped and unzipped CNVariantHash file.  Loading zipped file...");
      hashes = (CNVariantHash) SerializedFiles.readSerial(filename + suffixGz, log, false);
    }
    if ((!parse && !parseGz) || hashes == null) {
      if (hashes == null) {
        log.reportTime("Detected that CNVariantHash needs to be updated from cnv.var.CNVariantHash to filesys.CNVariantHash; reparsing...");
      }
      hashes = new CNVariantHash(filename, structureType, log);
      hashes.serialize(filename + suffixGz);
    }

    hashes.setFilename(filename);

    return hashes;
  }

  /**
   * This method will return all of the CNVs within a specified range
   *
   * @param location range to search
   * @param minProbes Minimum number of probes (markers)
   * @param minSizeKb Minimum size in kilobases
   * @param minQualityScore Minimum quality score
   * @return Sorted array of CNVs in the region
   */
  public CNVariant[] getAllInRegion(Segment location, int minProbes, int minSizeKb,
                                    int minQualityScore) {
    Vector<CNVariant> inRegion = new Vector<>();

    Enumeration<String> e = hashes.keys();
    while (e.hasMoreElements()) {
      Hashtable<String, CNVariant[]> h = hashes.get(e.nextElement());
      Enumeration<String> ee = h.keys();
      while (ee.hasMoreElements()) {
        CNVariant[] cnv = h.get(ee.nextElement());
        for (CNVariant element : cnv) {
          byte myChr = element.getChr();

          // Only look at the specified chromosome
          if (myChr == location.getChr()) {
            int myStart = element.getStart();
            int myStop = element.getStop();

            if ((myStop < location.getStart()) || (myStart > location.getStop())) {
              // We stop before the window, or start after it
            } else {
              int minSizeBases = minSizeKb * 1000;
              if ((element.getNumMarkers() >= minProbes) && (element.getSize() >= minSizeBases)
                  && (element.getScore() >= minQualityScore)) {
                inRegion.add(element);
              }
            }
          }
        }
      }
    }
    return CNVariant.sortCNVsInPlace(CNVariant.toCNVariantArray(inRegion));
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = "C:\\Documents and Settings\\npankrat\\My Documents\\CNV_PD\\data\\penncnv_1SNP.cnv";

    String usage = "\n" + "cnv.var.CNVariantHash requires 0-1 arguments\n"
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
      new CNVariantHash(filename, 1, new Logger());
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
