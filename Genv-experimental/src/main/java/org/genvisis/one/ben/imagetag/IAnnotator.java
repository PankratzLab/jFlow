package org.genvisis.one.ben.imagetag;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public interface IAnnotator {

  void loadImgDir(String dir);

  void loadAnnotations(String annotFile) throws IOException;

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
