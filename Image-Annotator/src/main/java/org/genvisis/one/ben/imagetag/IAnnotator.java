package org.genvisis.one.ben.imagetag;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.stream.Collectors;

public interface IAnnotator {

  static class LoadedSaved {
    int loaded;
    int missing;
    Map<String, Loaded> unmatched;

    public LoadedSaved() {
      unmatched = new HashMap<>();
    }

    String getString() {
      StringBuilder sb = new StringBuilder("Loaded ").append(loaded + missing)
                                                     .append(" annotated images from file");
      if (missing > 0) {
        sb.append(" (").append(missing).append(" missing)");
      }
      sb.append(".");
      if (unmatched.size() > 0) {
        sb.append("<br />");
        sb.append("<br />");
        sb.append(unmatched.entrySet().stream().map(ent -> {
          return "From " + ent.getKey() + "<br />" + ent.getValue().getString();
        }).collect(Collectors.joining("<br />")));
      }
      return sb.toString();
    }
  }

  static class Loaded {
    Map<String, Integer> unmatched;
    Map<String, Integer> total;

    public Loaded() {
      unmatched = new HashMap<>();
      total = new HashMap<>();
    }

    String getString() {
      if (unmatched.size() == 0) return "";
      StringBuilder sb = new StringBuilder();
      sb.append(unmatched.entrySet().stream().map(ent -> {
        return ent.getValue() + " (of " + total.get(ent.getKey()) + ") in " + ent.getKey()
               + "/ that do not start with " + ent.getKey();
      }).collect(Collectors.joining("<br />")));
      return sb.toString();
    }
  }

  Loaded loadImgDir(String dir);

  LoadedSaved loadAnnotations(String annotFile) throws IOException;

  void saveAnnotations(String annotFile);

  void saveAnnotation(AnnotatedImage.Annotation annotation, String annotFile);

  void addNewAnnotation(AnnotatedImage.Annotation newAnnotation);

  void replaceAnnotation(AnnotatedImage.Annotation prevAnnot, AnnotatedImage.Annotation newAnnot);

  void deleteAnnotation(AnnotatedImage.Annotation annotation);

  HashMap<String, HashMap<String, AnnotatedImage>> getAnnotationMap();

  ArrayList<String> getRoots();

  ArrayList<AnnotatedImage.Annotation> getAnnotations();

  void backupExistingFile(String file);
}
