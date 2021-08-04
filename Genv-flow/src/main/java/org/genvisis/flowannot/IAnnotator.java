package org.genvisis.flowannot;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.genvisis.fcs.auto.Panel;

public interface IAnnotator {

  void loadImgDir(String dir, List<Panel> panels);

  void loadAnnotations(String annotFile, List<Panel> panels) throws IOException;

  void saveAnnotations(String annotFile);

  void saveAnnotation(AnnotatedImage.Annotation annotation, String annotFile);

  void addNewAnnotation(AnnotatedImage.Annotation newAnnotation);

  void replaceAnnotation(AnnotatedImage.Annotation prevAnnot, AnnotatedImage.Annotation newAnnot);

  void deleteAnnotation(AnnotatedImage.Annotation annotation);

  HashMap<String, HashMap<String, AnnotatedImage>> getAnnotationMap();

  ArrayList<String> getFCSKeys();

  ArrayList<String> getFCSKeys(Panel panel);

  ArrayList<AnnotatedImage.Annotation> getAnnotations();

  void backupExistingFile(String file);
}
