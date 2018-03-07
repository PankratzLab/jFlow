package org.genvisis.one.ben.flowannot;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;
import java.util.Set;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.one.ben.flowannot.AnnotatedImage.Annotation;

public class AnnotationParser {

  String dir;
  Logger log;

  private String[] discover() {
    return Files.listFullPaths(dir, "");
  }

  private Set<String> loadAll(String[] filesPaths, String... caseInsensTags) {
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

  private Set<Annotation> loadAllPossible(String[] filesPaths) {
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

  public static void main(String[] args) {
    String d = "F:\\Flow\\boolGating\\annotation\\";
    AnnotationParser ap = new AnnotationParser();
    ap.dir = d;
    ap.log = new Logger();
    for (String a : ap.loadAll(ap.discover(), "manual")) {
      System.out.println(a);
    }
  }

}
