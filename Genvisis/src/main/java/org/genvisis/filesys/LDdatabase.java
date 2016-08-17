// used to filter out markers because of hwe, but since Haploview doesn't limit itself to controls,
// I truned this off (-hwcutoff 0)
// also added -minGeno 0 and -missingCutoff 1, to include shoddy coverage of regions by a mixture of
// chips
// always filter on hwe and call rate (snp and ind) before using as a reference
package org.genvisis.filesys;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.Date;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.Aliases;
import org.genvisis.common.Array;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.IntVector;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.common.Vectors;
import org.genvisis.common.ext;

public class LDdatabase implements Serializable {
  public static final long serialVersionUID = 1L;
  public static final float MISSING_INFO = Float.NEGATIVE_INFINITY;
  public static final float MONOMORPH_FIRST = -1;
  public static final float MONOMORPH_SECOND = -2;
  public static final float MONOMORPH_BOTH = -3;
  public static final float MISSING_FIRST = -11;
  public static final float MISSING_SECOND = -12;
  public static final float MISSING_BOTH = -13;
  public static final int TYPE_STRING = 1;
  public static final int TYPE_LONG = 2;
  public static final int NUM_CHROMOSOMES = 27;
  public static final int BP_LIMIT = 500000;

  public static final int COMP_FIRST = 1;
  public static final int COMP_MIN = 2;
  public static final int COMP_MAX = 3;

  public static final int REPORT_FIRST = 1;
  public static final int REPORT_LAST = 2;
  public static final int REPORT_VERBOSE = 3;
  public static final int REPORT_UNANIMOUS = 4;
  public static final String REPORT_INCONSISTENT = "inconsistent";

  public static final String[] HAPLOVIEW_LD_HEADER =
      {"L1", "L2", "D'", "LOD", "r^2", "CIlow", "CIhi", "Dist", "T-int"};
  public static final String[] FREQ_HEADER = {"CHR", "SNP", "A1", "A2", "MAF", "NCHROBS"};
  public static final String MASTER_HAPMAP_ROOT =
      ext.rootOf(Files.firstPathToFileThatExists(Aliases.REFERENCE_FOLDERS,
          "HapMap/CEU_founders/CEU_founders" + ".bim", true, false, new Logger()), false);
  public static final String HAPLOVIEW_LOC = "/home/npankrat/Haploview.jar";
  public static final String LDDB_ROOT = "lddb";
  public static final String LDDB_TARGETS = "targets";

  public static Hashtable<String, Hashtable<String, String>> getChrHashes(LDdatabase[] lddbs) {
    Hashtable<String, Hashtable<String, String>> chrHashes;

    chrHashes = new Hashtable<String, Hashtable<String, String>>();
    for (int i = 0; i < lddbs.length; i++) {
      chrHashes.put(i + "", lddbs[i].getChrHash());
    }

    return chrHashes;
  }

  public static LongLDdb[] getChrLDdbs(LDdatabase[] lddbs, int chr) {
    LongLDdb[] chrLDdbs;

    chrLDdbs = new LongLDdb[lddbs.length];
    for (int i = 0; i < lddbs.length; i++) {
      chrLDdbs[i] = lddbs[i].getLongChrLDdb(chr);
    }

    return chrLDdbs;
  }

  public static float multiCheck(LongLDdb[] chrLDdbs, String marker1, String marker2, int compType,
      Logger log) {
    float[] r2s;
    float r2;

    r2s = new float[chrLDdbs.length];
    for (int i = 0; i < chrLDdbs.length; i++) {
      r2s[i] = chrLDdbs[i].get(marker1, marker2);
    }

    r2 = MISSING_INFO;
    switch (compType) {
      case COMP_FIRST:
        for (float r22 : r2s) {
          if (r22 > MISSING_INFO) {
            return r22;
          }
        }
        break;
      case COMP_MIN:
        for (float r22 : r2s) {
          if (r2 == MISSING_INFO || (r2 < 0 && r22 >= 0)) {
            r2 = r22;
          } else if (r22 >= 0) {
            r2 = Math.min(r2, r22);
          }
        }
        break;
      case COMP_MAX:
        for (float r22 : r2s) {
          if (r22 > r2) {
            r2 = r22;
          }
        }
        break;
      default:
        log.reportError("Error - invalid LDdatabase array comparison type: " + compType);
        System.exit(1);
    }

    return r2;
  }

  public static String multiHashCheck(Hashtable<String, Hashtable<String, String>> chrHashes,
      String markerName, int reportType, Logger log) {
    String[] poslar;
    String pos;

    poslar = new String[chrHashes.size()];
    for (int i = 0; i < chrHashes.size(); i++) {
      poslar[i] = chrHashes.get(i + "").get(markerName);
    }

    pos = null;
    switch (reportType) {
      case REPORT_FIRST:
        for (String element : poslar) {
          if (element != null) {
            return element;
          }
        }
        break;
      case REPORT_LAST:
        for (int i = poslar.length - 1; i >= 0; i--) {
          if (poslar[i] != null) {
            return poslar[i];
          }
        }
        break;
      case REPORT_VERBOSE:
      case REPORT_UNANIMOUS:
        for (int i = 0; i < poslar.length; i++) {
          if (pos == null) {
            pos = poslar[i];
          } else if (!poslar[i].equals(pos)) {
            if (reportType == REPORT_VERBOSE) {
              log.reportError("Error - inconsistent position for marker " + markerName + ": " + pos
                  + " and " + poslar[i]);
            }
            return REPORT_INCONSISTENT;
          }
        }
        break;
      default:
        log.reportError("Error - invalid LDdatabase array report type: " + reportType);
        System.exit(1);
    }

    return pos;
  }

  public static LDdatabase[] multiLDdbs(String[] roots, Logger log) {
    LDdatabase[] lddbs;

    lddbs = new LDdatabase[roots.length];
    for (int i = 0; i < roots.length; i++) {
      lddbs[i] = new LDdatabase(roots[i], LDdatabase.TYPE_LONG, log);
    }

    return lddbs;
  }

  public static void multiUpdate(LDdatabase[] lddbs, String[] markers, String listName) {
    for (LDdatabase lddb : lddbs) {
      lddb.updateWithTheseMarkers(markers, listName);
    }
  }

  private final String root;

  private final String dir;

  private final int type;

  private boolean unchanged;

  private Hashtable<String, String> chrHash, subHash;

  private final Logger log;

  public LDdatabase(int type, Logger log) {
    this(MASTER_HAPMAP_ROOT, type, log);
  }

  public LDdatabase(String root, int type, Logger log) {
    this.root = root;
    dir = root + "_" + LDDB_ROOT + "/";
    this.type = type;
    unchanged = false;
    chrHash = null;
    subHash = null;
    this.log = log;

    new File(dir).mkdirs();
  }

  public Hashtable<String, String> getChrHash() {
    if (unchanged) {
      return subHash;
    } else {
      if (chrHash == null) {
        chrHash =
            SnpMarkerSet.loadSnpMarkerSetToChrHash(root + ".bim", SnpMarkerSet.PLINK_BIM_FORMAT);
      }
      return chrHash;
    }
  }

  public String getDir() {
    return dir;
  }

  public LongLDdb getLongChrLDdb(int chr) {
    return LongLDdb.load(dir + "chr" + chr, false, true);
  }

  public StringLDdb getStringChrLDdb(int chr) {
    return StringLDdb.load(dir + "chr" + chr, false, true);
  }

  public int getType() {
    return type;
  }

  public boolean isUnchanged() {
    return unchanged;
  }

  public void updateWithTheseMarkers(String[] targets) {
    updateWithTheseMarkers(targets, null);
  }

  public void updateWithTheseMarkers(String[] targets, String listName) {
    BufferedReader reader;
    PrintWriter writer;
    Vector<String> v, monomorphs;
    String[] line, subset, check;
    int[] positions;
    long time;
    IntVector[] ivs;
    String trav;
    // StringLDdb chrLDdb;
    LongLDdb chrLDdb;
    Hashtable<String, Vector<String>> missing;
    SnpMarkerSet missingMarkers;
    byte[] chrs;

    if (listName != null) {
      if (Files.exists(dir + listName + ".list", false)) {
        check =
            HashVec.loadFileToStringArray(dir + listName + ".list", false, new int[] {0}, false);
        if (Array.equals(targets, check, false)) {
          unchanged = true;
          subHash = SerialHash.loadSerializedStringHash(dir + listName + ".chrHash.ser");
          return;
        } else {
          Files.writeList(targets, dir + listName + ".list_");
        }
      } else {
        Files.writeList(targets, dir + listName + ".list_");
      }
    }

    if (chrHash == null) {
      chrHash =
          SnpMarkerSet.loadSnpMarkerSetToChrHash(root + ".bim", SnpMarkerSet.PLINK_BIM_FORMAT);
    }
    subHash = new Hashtable<String, String>();

    time = new Date().getTime();
    log.report("Sorting markers by chromosome");
    ivs = Vectors.initializedArray(IntVector.class, NUM_CHROMOSOMES);
    missing = new Hashtable<String, Vector<String>>();
    for (int i = 0; i < targets.length; i++) {
      trav = chrHash.get(targets[i]);
      if (trav == null) {
        HashVec.addToHashVec(missing, "0", targets[i], true);
      } else {
        ivs[Integer.parseInt(trav.split("[\\s]+")[0])].add(i);
        subHash.put(targets[i], trav);
      }
    }

    if (listName != null) {
      SerialHash.createSerializedStringHash(dir + listName + ".chrHash.ser", subHash);
    }

    if (missing.containsKey("0") && missing.get("0").size() > 0) {
      v = missing.get("0");
      log.reportError("Warning - the following " + v.size() + " marker" + (v.size() > 1 ? "s" : "")
          + " were not found in the " + root + " dataset:");
      for (int i = 0; i < v.size(); i++) {
        log.reportError(v.elementAt(i));
      }
      missingMarkers = new SnpMarkerSet(Array.toStringArray(v));
      missingMarkers.parseSNPlocations();
      chrs = missingMarkers.getChrs();
      for (int i = 0; i < v.size(); i++) {
        HashVec.addToHashVec(missing, chrs[i] + "", targets[i], true);
      }
    }

    time = new Date().getTime();
    for (int chr = 1; chr < NUM_CHROMOSOMES; chr++) {
      if (ivs[chr].size() > 2) {
        log.report(ext.getTime() + "\tchr" + chr);
        chrLDdb = LongLDdb.load(dir + "chr" + chr, false, true);

        subset = new String[ivs[chr].size()];
        positions = new int[ivs[chr].size()];
        for (int j = 0; j < ivs[chr].size(); j++) {
          subset[j] = targets[ivs[chr].elementAt(j)];
          positions[j] = Integer.parseInt(chrHash.get(subset[j]).split("[\\s]+")[1]);
        }

        log.report("Checking how many chr" + chr + " pairs are present already...", false, true);
        v = new Vector<String>();
        try {
          writer = null;
          for (int i = 0; i < subset.length; i++) {
            for (int j = i + 1; j < subset.length; j++) {
              if (Math.abs(positions[i] - positions[j]) < BP_LIMIT
                  && chrLDdb.get(subset[i], subset[j]) == MISSING_INFO) {
                HashVec.addIfAbsent(subset[i], v);
                HashVec.addIfAbsent(subset[j], v);
                if (writer == null) {
                  writer = new PrintWriter(new FileWriter(dir + "pairs.xln"));
                }
                writer.println(
                    subset[i] + "\t" + subset[j] + "\t" + Math.abs(positions[i] - positions[j]));
              }
            }
          }
          if (writer != null) {
            writer.close();
          }
        } catch (Exception e) {
          log.reportError("Error writing to " + dir + "pairs.xln");
          e.printStackTrace();
        }
        log.report("done");

        if (v.size() == 0) {
          log.report("All calculations have been performed previously");
        } else {
          log.report("Extracting chr" + chr + " markers from " + root, false, true);
          Files.writeList(Array.toStringArray(v), dir + LDDB_TARGETS + "_snplist.txt");
          if (!new File(root + ".chr" + chr + ".bed").exists()) {
            CmdLine.run("plink --bfile " + root + " --chr " + chr + " --noweb --make-bed --out "
                + root + ".chr" + chr, "./");
          }
          CmdLine.run("plink --bfile " + root + ".chr" + chr + " --noweb --extract " + dir
              + LDDB_TARGETS + "_snplist.txt --recode --out " + dir + LDDB_TARGETS, "./");
          CmdLine.run(
              "plink --file " + dir + LDDB_TARGETS + " --noweb --freq --out " + dir + "freqCheck",
              "./");
          try {
            reader = new BufferedReader(new FileReader(dir + "freqCheck.frq"));
            ext.checkHeader(reader.readLine().trim().split("[\\s]+"), FREQ_HEADER, true);
            monomorphs = new Vector<String>();
            while (reader.ready()) {
              line = reader.readLine().trim().split("[\\s]+");
              if (line[4].equals("0")) {
                monomorphs.add(line[1]);
              }
            }
            chrLDdb.addMonomorphs(monomorphs);
            reader.close();
          } catch (FileNotFoundException fnfe) {
            log.reportError(
                "Error: file \"" + dir + "freqCheck.frq" + "\" not found in current directory");
            System.exit(1);
          } catch (IOException ioe) {
            log.reportError("Error reading file \"" + dir + "freqCheck.frq" + "\"");
            System.exit(2);
          }
          new SnpMarkerSet(dir + LDDB_TARGETS + ".map").writeToFile(dir + LDDB_TARGETS + ".info",
              SnpMarkerSet.HAPLOVIEW_INFO_FORMAT, log);
          log.report("  ...done");
          log.report("Computing LD for chr" + chr + "...", false, true);
          if (new File(HAPLOVIEW_LOC).exists()) {
            CmdLine.run("java -jar " + HAPLOVIEW_LOC + " -nogui -pedfile " + dir + LDDB_TARGETS
                + ".ped -info " + dir + LDDB_TARGETS + ".info -dprime -chromosome "
                + Positions.chromosomeNumberInverse(chr)
                + " -maxDistance 500 -hwcutoff 0 -minGeno 0 -missingCutoff 1", "./");
          } else {
            log.reportError("\n\nError - could not find haploview JAR file at: " + HAPLOVIEW_LOC
                + "\nAborting...");
            System.exit(1);
          }

          try {
            reader = new BufferedReader(new FileReader(dir + LDDB_TARGETS + ".ped.LD"));
            ext.checkHeader(reader.readLine().trim().split("[\\s]+"), HAPLOVIEW_LD_HEADER, false);
            while (reader.ready()) {
              line = reader.readLine().trim().split("[\\s]+");
              chrLDdb.add(line[0], line[1], Float.parseFloat(line[4]));
            }
            reader.close();
          } catch (FileNotFoundException fnfe) {
            log.reportError("Error: file \"" + dir + LDDB_TARGETS + ".ped.LD"
                + "\" not found in current directory");
            System.exit(1);
          } catch (IOException ioe) {
            log.reportError("Error reading file \"" + dir + LDDB_TARGETS + ".ped.LD" + "\"");
            System.exit(2);
          }

          log.report("done");
        }
        if (chrLDdb.getChanged()) {
          chrLDdb.serialize(dir + "chr" + chr);
        }
      }
    }
    new File(dir + listName + ".list_").renameTo(new File(dir + listName + ".list"));
    log.report("Finished update in " + ext.getTimeElapsed(time));
  }
}
