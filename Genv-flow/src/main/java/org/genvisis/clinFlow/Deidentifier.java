package org.genvisis.clinFlow;

import java.io.File;
import java.io.FileFilter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.genvisis.jfcs.FCSKeywords;
import org.genvisis.jfcs.FCSReader;
import org.genvisis.jfcs.FCS_KEYWORD;
import org.pankratzlab.common.CLI;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;

public class Deidentifier {

  private static final String KEY_PATIENT_ID = "PATIENT ID";
  private static final String KEY_EXPERIMENT_NAME = "EXPERIMENT NAME";
  private static final String KEY_SAMPLE_ID = "SAMPLE ID";
  private static final String KEY_GUID = "GUID";
  private static final String KEY_SRC = "$SRC";
  private static final String KEY_FIL = "$FIL";

  private static final List<String> keys = new ArrayList<>();

  {
    keys.add(KEY_PATIENT_ID);
    keys.add(KEY_EXPERIMENT_NAME);
    keys.add(KEY_SAMPLE_ID);
    keys.add(KEY_SRC);
    keys.add(KEY_FIL);
  }

  private String rootIn;
  private String rootOut;
  private String linkDir;

  public Deidentifier(String in, String out, String link) {
    this.rootIn = ext.verifyDirFormat(in);
    this.rootOut = ext.verifyDirFormat(out);
    this.linkDir = ext.verifyDirFormat(link);
  }

  public List<Conversion> identify() {
    File rootDir = new File(rootIn);
    File outDir = new File(rootOut);

    List<Conversion> allConvs = processDir(rootDir, outDir);
    // removeExisting(allConvs);
    return allConvs;
  }

  // private void removeExisting(List<Conversion> convs) {
  // for (int i = convs.size() - 1; i >= 0; i--) {
  // try {
  // if (exists(convs.get(i))) {
  // convs.remove(i);
  // }
  // } catch (IOException e) {
  // Conversion c = convs.remove(i);
  // cantOpen(e, path(c.dir) + c.fcs);
  // }
  // }
  // }

  public static boolean exists(Conversion c) throws IOException {
    if (!c.out.exists()) {
      return false;
    }

    FCSKeywords keys;
    try {
      keys = FCSReader.readKeywords(path(c.dir) + c.fcs);
    } catch (IOException e) {
      cantOpen(e, path(c.dir) + c.fcs);
      return true;
    }

    if (!keys.hasKeyword(KEY_GUID)) {
      // can't determine, recreate
      return false;
    }
    String newID = generateNewID(keys);
    String newF = path(c.out) + newID + ".fcs";
    if (Files.exists(newF)) {
      keys = FCSReader.readKeywords(newF);
      if (keys.getKeyword(KEY_GUID).equals(keys.getKeyword(KEY_SAMPLE_ID))) {
        // sample is already set to guid
        return true;
      }
    }
    return false;
  }

  private List<Conversion> processDir(File dir, File out) {
    File[] dirs = listDirs(dir);

    String[] fcss = listFCS(dir);

    List<Conversion> convs = new ArrayList<>();
    for (String f : fcss) {
      convs.add(new Conversion(dir, out, linkDir, f));
    }
    for (File d : dirs) {
      convs.addAll(processDir(d, new File(out, generateNewDirName(d))));
    }
    return convs;
  }

  private File[] listDirs(File dir) {
    return dir.listFiles(
        new FileFilter() {

          @Override
          public boolean accept(File arg0) {
            return arg0.isDirectory();
          }
        });
  }

  private String[] listFCS(File dir) {
    return dir.list(
        new FilenameFilter() {

          @Override
          public boolean accept(File dir, String name) {
            return name.endsWith(".fcs");
          }
        });
  }

  private String generateNewDirName(File dir) {
    String name = dir.getName();

    long trav = 1;
    for (int j = 0; j < name.length(); j++) {
      trav += name.charAt(j) * (j + 1);
    }

    return Long.toString(trav);
  }

  public static void processSingleFCS(Conversion conv) {
    FCSReader reader;
    try {
      reader = FCSReader.open(path(conv.dir) + conv.fcs);
    } catch (IOException e) {
      cantOpen(e, conv.fcs);
      return;
    }

    String newID = generateNewID(reader.getKeywords());
    Map<String, String> idents = getIdentifiers(reader.getKeywords());
    fixKeywords(newID, reader.getKeywords());
    String outPath = path(conv.out);
    conv.out.mkdirs();
    try {
      FCSReader.write(reader, outPath + newID + ".fcs");
    } catch (IOException e) {
      writeFail(e, conv.fcs, outPath + newID + ".fcs");
      return;
    }
    idents.put("SOURCE PATH", path(conv.dir));
    idents.put("OUTPUT PATH", outPath);
    writeLinkFile(path(conv.link == null ? conv.out : new File(conv.link)), newID, idents);
    reader.dispose();
  }

  private static String generateNewID(FCSKeywords ks) {
    String id;
    if (ks.hasKeyword(KEY_GUID)) {
      id = ks.getKeyword(KEY_GUID);
    } else {
      id = java.util.UUID.randomUUID().toString();
      ks.setKeyword(KEY_GUID, id);
    }
    return id;
  }

  private static Map<String, String> getIdentifiers(FCSKeywords ks) {
    Map<String, String> id = new HashMap<>();
    for (String k : keys) {
      id.put(k, ks.getKeyword(k));
    }
    return id;
  }

  private static void fixKeywords(String newID, FCSKeywords ks) {
    for (String k : keys) {
      ks.setKeyword(k, newID);
    }
    ks.setKeyword(FCS_KEYWORD.FIL, newID + ".fcs");
  }

  private static void writeLinkFile(String outDir, String newID, Map<String, String> idents) {
    PrintWriter writer = Files.getAppropriateWriter(outDir + newID + ".txt");

    for (Entry<String, String> e : idents.entrySet()) {
      writer.println(e.getKey() + "\t" + e.getValue());
    }

    writer.close();
  }

  private static String path(File f) {
    return ext.verifyDirFormat(f.getAbsolutePath());
  }

  private static void cantOpen(IOException e, String fcsFile) {
    // TODO
  }

  private static void writeFail(IOException e, String fcsFile, String newFile) {
    // TODO
  }

  public static void main(String[] args) {
    CLI cli = new CLI(Deidentifier.class);

    cli.addArg("in", "Input directory", true);
    cli.addArg("out", "Output directory", true);
    cli.addArg("link", "Link file directory", true);

    cli.parseWithExit(args);

    Logger log = new Logger();
    Deidentifier deid = new Deidentifier(cli.get("in"), cli.get("out"), cli.get("link"));
    List<Conversion> toRun = deid.identify();
    log.report("Found " + toRun.size() + " files to convert.");
    for (Conversion c : toRun) {
      try {
        if (Deidentifier.exists(c)) {
          log.report("Found existing conversion results for " + c.fcs);
          continue;
        }
      } catch (IOException e) {
        // just redo
      }
      Deidentifier.processSingleFCS(c);
      log.report("Converted " + c.fcs);
    }
  }
}
