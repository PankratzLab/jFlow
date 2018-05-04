package org.genvisis.one.ben.flowannot;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;
import java.util.stream.Collectors;
import org.apache.commons.io.FileUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.one.ben.flowannot.AnnotatedImage.Annotation;

public class AnnotationParser {

  String dir;
  Logger log;

  private List<String> discover() {
    return FileUtils.listFiles(new File(dir), null, true).stream().map(File::getAbsolutePath)
                    .collect(Collectors.toList());
  }

  private Set<String> loadSamples(List<String> filesPaths) {
    Set<String> files = new HashSet<>();
    for (String f : filesPaths) {
      Annotator annotator = new Annotator();
      try {
        annotator.loadAnnotations(f);
      } catch (IOException e) {
        log.reportError("Failed to load file " + f + " -- " + e.getMessage());
        continue;
      }
      for (Entry<String, HashMap<String, AnnotatedImage>> en : annotator.getAnnotationMap()
                                                                        .entrySet()) {
        String fi = en.getKey();
        files.add(fi);
      }
    }
    return files;
  }

  private Set<String> loadAll(List<String> filesPaths, String... caseInsensTags) {
    Set<String> anns = new HashSet<>();
    for (String f : filesPaths) {
      Annotator annotator = new Annotator();
      try {
        annotator.loadAnnotations(f);
      } catch (IOException e) {
        log.reportError("Failed to load file " + f + " -- " + e.getMessage());
        continue;
      }
      for (Entry<String, HashMap<String, AnnotatedImage>> en : annotator.getAnnotationMap()
                                                                        .entrySet()) {
        String fi = en.getKey();
        for (AnnotatedImage ai : en.getValue().values()) {
          for (Annotation a : ai.getAnnotations()) {
            for (String tag : caseInsensTags) {
              if (a.annotation.equalsIgnoreCase(tag)) {
                anns.add(fi);
              }
            }
          }
        }
      }
    }
    return anns;
  }

  private Set<Annotation> loadAllPossible(List<String> filesPaths) {
    Set<Annotation> anns = new HashSet<>();
    for (String f : filesPaths) {
      Annotator annotator = new Annotator();
      try {
        annotator.loadAnnotations(f);
      } catch (IOException e) {
        log.reportError("Failed to load file " + f + " -- " + e.getMessage());
        continue;
      }
      anns.addAll(annotator.getAnnotations());
    }
    return anns;
  }

  private void loadAllMapped(List<String> filesPaths, String out) {
    PrintWriter writer = Files.getAppropriateWriter(out);

    for (String f : filesPaths) {
      String s = ext.rootOf(f, true);
      Annotator annotator = new Annotator();
      try {
        annotator.loadAnnotations(f);
      } catch (IOException e) {
        log.reportError("Failed to load file " + f + " -- " + e.getMessage());
        continue;
      }
      for (Entry<String, HashMap<String, AnnotatedImage>> en : annotator.getAnnotationMap()
                                                                        .entrySet()) {
        String fcsFi = en.getKey();
        for (AnnotatedImage ai : en.getValue().values()) {
          for (Annotation a : ai.getAnnotations()) {
            writer.println(s + "\t" + fcsFi + "\t" + ai.getGateName() + "\t" + a.annotation);
          }
        }
      }
    }

    writer.close();
  }

  public static void main(String[] args) {
    String d = "F:\\Flow\\final results\\annotations\\src\\";
    AnnotationParser ap = new AnnotationParser();
    ap.dir = d;
    ap.log = new Logger();
    //    for (String a : ap.loadAll(ap.discover(), "manual")) {
    //      System.out.println(a);
    //    }
    //    for (String a : ap.loadSamples(ap.discover())) {
    //      System.out.println(a);
    //    }
    ap.loadAllMapped(ap.discover(), d + "parsedAnnots.xln");
  }

}
