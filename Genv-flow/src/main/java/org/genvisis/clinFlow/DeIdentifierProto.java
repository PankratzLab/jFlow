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
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import org.genvisis.common.Files;
import org.genvisis.common.ext;
import org.genvisis.jfcs.FCSKeywords;
import org.genvisis.jfcs.FCSReader;
import org.genvisis.jfcs.FCS_KEYWORD;

public class DeIdentifierProto {

  // input directory
  // fcs output directory
  // link output directory

  /*-
   * Process per file:
   *  read fcs file
   *  pull keywords
   *  remove/replace keywords
   *  write to link file
   *  write to new fcs file
   *    - rewrite header (text start / end, data start / end if existed prev, analysis start/end if existed prev - if in keywords, replace values)
   */

  String KEY_PATIENT_ID = "PATIENT ID";
  String KEY_EXPERIMENT_NAME = "EXPERIMENT NAME";
  String KEY_SAMPLE_ID = "SAMPLE ID";
  String KEY_GUID = "GUID";
  String KEY_SRC = "$SRC";
  String KEY_FIL = "$FIL";

  List<String> keys = new ArrayList<>();
  {
    keys.add(KEY_PATIENT_ID);
    keys.add(KEY_EXPERIMENT_NAME);
    keys.add(KEY_SAMPLE_ID);
    keys.add(KEY_SRC);
    keys.add(KEY_FIL);
  }

  String rootIn;
  String rootOut;

  class Conversion {

    public Conversion(File dir2, File out2, String f) {
      this.dir = dir2;
      this.out = out2;
      this.fcs = f;
    }

    File dir;
    File out;
    String fcs;
  }

  private void start() {
    File rootDir = new File(rootIn);
    File outDir = new File(rootOut);

    List<Conversion> allConvs = processDir(rootDir, outDir);
    removeExisting(allConvs);

    ExecutorService service = Executors.newFixedThreadPool(Runtime.getRuntime()
                                                                  .availableProcessors()
                                                           - 1);
    for (Conversion c : allConvs) {
      // create new progress bar for each file, update as run

      service.submit(() -> {
        processSingleFCS(c);
      });
    }
    service.shutdown();
    try {
      service.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
    } catch (InterruptedException e) {}
    Thread.yield();
  }

  private void removeExisting(List<Conversion> convs) {
    for (int i = convs.size() - 1; i >= 0; i++) {
      try {
        if (exists(convs.get(i))) {
          convs.remove(i);
        }
      } catch (IOException e) {
        Conversion c = convs.remove(i);
        cantOpen(e, path(c.dir) + c.fcs);
      }
    }
  }

  private boolean exists(Conversion c) throws IOException {
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
      convs.add(new Conversion(dir, out, f));
    }
    for (File d : dirs) {
      convs.addAll(processDir(d, new File(out, generateNewDirName(d))));
    }
    return convs;
  }

  private File[] listDirs(File dir) {
    return dir.listFiles(new FileFilter() {

      @Override
      public boolean accept(File arg0) {
        return arg0.isDirectory();
      }
    });
  }

  private String[] listFCS(File dir) {
    return dir.list(new FilenameFilter() {

      @Override
      public boolean accept(File dir, String name) {
        return name.endsWith(".fcs");
      }
    });
  }

  private String generateNewDirName(File dir) {
    String name = dir.getName();

    int prime = 31;
    long trav = 0;
    for (int j = 0; j < name.length(); j++) {
      trav += name.charAt(j) * prime;
    }

    return Long.toString(trav);
  }

  private String generateNewID(FCSKeywords ks) {
    String id;
    if (ks.hasKeyword(KEY_GUID)) {
      id = ks.getKeyword(KEY_GUID);
    } else {
      id = java.util.UUID.randomUUID().toString();
      ks.setKeyword(KEY_GUID, id);
    }
    return id;
  }

  private Map<String, String> getIdentifiers(FCSKeywords ks) {
    Map<String, String> id = new HashMap<>();
    for (String k : keys) {
      id.put(k, ks.getKeyword(k));
    }
    return id;
  }

  private void writeLinkFile(String outDir, String newID, Map<String, String> idents) {
    PrintWriter writer = Files.getAppropriateWriter(outDir + newID + ".txt");

    writer.close();
  }

  private void fixKeywords(String newID, FCSKeywords ks) {
    for (String k : keys) {
      ks.setKeyword(k, newID);
    }
    ks.setKeyword(FCS_KEYWORD.FIL, newID + ".fcs");
  }

  private void processSingleFCS(Conversion conv) {
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
    try {
      FCSReader.write(reader, outPath + newID + ".fcs");
    } catch (IOException e) {
      writeFail(e, conv.fcs, outPath + newID + ".fcs");
      return;
    }
    writeLinkFile(outPath, newID, idents);
    reader.dispose();
  }

  private String path(File f) {
    return ext.verifyDirFormat(f.getAbsolutePath());
  }

  private void cantOpen(IOException e, String fcsFile) {
    // TODO
  }

  private void writeFail(IOException e, String fcsFile, String newFile) {
    // TODO
  }

  public void run() {
    rootIn = "F:/Flow_stage2/source/";
    rootOut = "F:/Flow_stage2/deident/";
    start();
  }

  public static void main(String[] args) {
    new DeIdentifierProto().run();
  }

}
