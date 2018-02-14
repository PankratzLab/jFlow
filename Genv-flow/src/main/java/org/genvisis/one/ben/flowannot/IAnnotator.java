package org.genvisis.one.ben.flowannot;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public interface IAnnotator {

  public static enum PANEL {
    PANEL_1() {

      @Override
      public boolean isPanel(String p) {
        return p.toLowerCase().contains("panel 1") || p.toLowerCase().contains("panel_1");
      }
    },
    PANEL_2() {

      @Override
      public boolean isPanel(String p) {
        return p.toLowerCase().contains("panel 2") || p.toLowerCase().contains("panel_2");
      }
    };

    public abstract boolean isPanel(String p);
  }

  void loadImgDir(String dir);

  void loadAnnotations(String annotFile) throws IOException;

  void saveAnnotations(String annotFile);

  void saveAnnotation(AnnotatedImage.Annotation annotation, String annotFile);

  void addNewAnnotation(AnnotatedImage.Annotation newAnnotation);

  void replaceAnnotation(AnnotatedImage.Annotation prevAnnot, AnnotatedImage.Annotation newAnnot);

  void deleteAnnotation(AnnotatedImage.Annotation annotation);

  HashMap<String, HashMap<String, AnnotatedImage>> getAnnotationMap();

  ArrayList<String> getFCSKeys();

  ArrayList<String> getFCSKeys(PANEL panel);

  ArrayList<AnnotatedImage.Annotation> getAnnotations();

  void backupExistingFile(String file);

}
