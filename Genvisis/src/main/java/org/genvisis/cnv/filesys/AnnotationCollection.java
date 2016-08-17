package org.genvisis.cnv.filesys;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Vector;

import javax.swing.JOptionPane;

import org.genvisis.cnv.manage.TextExport;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.Sort;
import org.genvisis.common.ext;

public class AnnotationCollection implements Serializable, TextExport {
  private static final long serialVersionUID = 1L;

  public static void appendFromLists(String filename,
      AnnotationCollection currentAnnotationCollection, String[] markerNamesProj, Logger log) {
    BufferedReader reader;
    String[] line;
    char key;
    char[] keys;
    Hashtable<String, Character> annotationKeys;
    boolean found;

    annotationKeys = new Hashtable<String, Character>();
    keys = currentAnnotationCollection.getKeys();
    for (char key2 : keys) {
      annotationKeys.put(currentAnnotationCollection.getDescriptionForComment(key2, false, false),
          key2);
    }
    try {
      reader = new BufferedReader(new FileReader(filename));
      while (reader.ready()) {
        key = 0;
        line = reader.readLine().split("[\t;]");
        found = false;
        if (markerNamesProj != null) {
          for (String element : markerNamesProj) {
            if (element.equals(line[0])) {
              found = true;
            }
          }
        }
        if (markerNamesProj != null && !found) {
          log.reportError("Skipped importing annotations for marker " + line[0]
              + ", which is found in the annotation file " + filename
              + ", but not in the marker name list of the project.");
        } else {
          for (int i = 1; i < line.length; i++) {
            if (annotationKeys.containsKey(line[i])) {
              key = annotationKeys.get(line[i]);
            } else {
              key = assignKey(line[i], currentAnnotationCollection);
              if (key == 0) {
                if (log == null) {
                  System.out.println(
                      "cannot automatically assign a shortcut to the annotation key '" + line[i]
                          + "'. Skipped importing this annotation for marker " + line[0] + ".");
                } else {
                  log.reportError(
                      "cannot automatically assign a shortcut to the annotation key '" + line[i]
                          + "'. Skipped importing this annotation for marker " + line[0] + ".");
                }
              } else {
                annotationKeys.put(line[i], key);
                currentAnnotationCollection.addAnnotation(key, line[i]);
              }
            }
            currentAnnotationCollection.addAnnotationForMarker(line[0], key);
          }
        }
      }
      reader.close();
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  public static char assignKey(String comment, AnnotationCollection currentAnnotationCollection) {
    char key;
    char[] letters;

    key = 0;
    letters = comment.toLowerCase().toCharArray();
    for (int j = 0; j < letters.length; j++) {
      if (letters[j] >= 97 && letters[j] <= 122
          && !currentAnnotationCollection.containsKey(letters[j])) {
        key = letters[j];
        break;
      }
    }
    if (key == 0) {
      for (int j = 97; j <= 122; j++) {
        if (!currentAnnotationCollection.containsKey((char) j)) {
          key = (char) j;
          break;
        }
      }
    }

    return key;
  }

  public static AnnotationCollection load(String filename, boolean jar) {
    return (AnnotationCollection) SerializedFiles.readSerial(filename, jar, true);
  }

  public static AnnotationCollection loadFromLists(String filename, Logger log) {
    AnnotationCollection annotationCollection;
    annotationCollection = new AnnotationCollection();
    appendFromLists(filename, annotationCollection, null, log);
    return annotationCollection;
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String logfile = null;
    Logger log;
    Project proj;
    String filename = null;
    String exportList = "listOfAnnotations.out";
    String importList = null;
    AnnotationCollection annotationCollection;
    boolean dump = false;
    String recoverDir = null;

    String usage = "\n" + "cnv.filesys.AnnotationCollection requires 0-1 arguments\n"
        + "   (1) project properties filename (i.e. proj="
        + org.genvisis.cnv.Launch.getDefaultDebugProjectFile(false) + " (default))\n" + "  AND\n"
        + "   (2) list annotations (i.e. exportList=" + exportList + " (default))\n" + "  OR\n"
        + "   (2) rename annotations (i.e. importList=importNewList.txt (not the default))\n"
        + "  OR\n"
        + "   (2) dump lists for the project's AnnotationCollection (i.e. -dump (not the default))\n"
        + "  OR\n"
        + "   (1) recover temp annotations from a directory (i.e. recover=C:/data/recover/ (not the default))\n"
        + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("proj=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("exportList=")) {
        exportList = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("importList=")) {
        importList = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("recover=")) {
        recoverDir = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("-dump")) {
        dump = true;
        numArgs--;
      } else if (arg.startsWith("log=")) {
        logfile = arg.split("=")[1];
        numArgs--;
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      if (recoverDir != null) {
        recover(recoverDir);
      } else {
        proj = new Project(filename, logfile, false);
        log = proj.getLog();

        annotationCollection = proj.getAnnotationCollection();
        if (dump) {
          annotationCollection.dumpLists(proj);
        } else if (importList != null) {
          // proj.archiveFile(proj.getFilename(proj.ANNOTATION_FILENAME));
          proj.archiveFile(proj.ANNOTATION_FILENAME.getValue());
          annotationCollection.importList(importList, log);
          // annotationCollection.serialize(proj.getFilename(proj.ANNOTATION_FILENAME, false,
          // false));
          annotationCollection.serialize(proj.ANNOTATION_FILENAME.getValue(false, false));
        } else {
          annotationCollection.exportToText(proj, exportList);
        }
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public static void recover(String dir) {
    AnnotationCollection annotationCollection;
    String[] files;
    String trav;

    files = Files.list(dir, ".tempAnnotation.ser", false);
    for (String file : files) {
      trav = file.substring(0, file.indexOf(".tempAnnotation.ser"));
      new File(dir + trav + "/").mkdirs();
      annotationCollection = load(dir + file, false);
      annotationCollection.dumpLists(dir + trav + "/");
    }
  }

  private final Hashtable<Character, String> commentsHash; // annotation keys

  private final Hashtable<String, Vector<String>> markerAnnotations; // annotations organized by
                                                                     // markers. markerName ->
                                                                     // AnnotationA, AnnotationB,
                                                                     // ...

  private final Hashtable<String, Vector<String>> annotationMarkerLists; // annotations organized by
                                                                         // annotation keys.
                                                                         // Annotation ->
                                                                         // marker1Name,
                                                                         // marker2Name, ...

  public AnnotationCollection() {
    commentsHash = new Hashtable<Character, String>();
    markerAnnotations = new Hashtable<String, Vector<String>>();
    annotationMarkerLists = new Hashtable<String, Vector<String>>();
  }

  public void addAnnotation(char c, String comment) {
    commentsHash.put(c, comment);
    annotationMarkerLists.put(c + "", new Vector<String>());
  }

  public void addAnnotationForMarker(String markerName, char c) {
    HashVec.addToHashVec(markerAnnotations, markerName, c + "", true);
    HashVec.addToHashVec(annotationMarkerLists, c + "", markerName, true);
  }

  public char[] annotationsForMarker(String markerName) {
    if (markerAnnotations.containsKey(markerName)) {
      return Array.toCharArray(Array.toStringArray(markerAnnotations.get(markerName)));
    } else {
      return new char[0];
    }
  }

  public boolean containsKey(char c) {
    return commentsHash.containsKey(c);
  }

  public void dumpLists(Project proj) {
    dumpLists(proj.PROJECT_DIRECTORY.getValue());
  }

  public void dumpLists(String outputDir) {
    String[] list, keysAnnotations, keysMarkers;
    String[][] matrix;
    Vector<String> annotationsVector;
    String annotationsOfTheMarker;

    keysAnnotations = HashVec.getKeys(annotationMarkerLists);
    for (String keysAnnotation : keysAnnotations) {
      Files
          .writeList(Array.toStringArray(annotationMarkerLists.get(keysAnnotation)),
              outputDir + "annotation_" + keysAnnotation + "_"
                  + ext.replaceWithLinuxSafeCharacters(
                      getDescriptionForComment(keysAnnotation.charAt(0), false, false), true)
                  + ".xln");
    }

    list = new String[markerAnnotations.size()];
    matrix = new String[markerAnnotations.size() + 2][];
    matrix[0] = Array.addStrToArray("", keysAnnotations, 0);
    matrix[1] = Array.addStrToArray("MarkerName", keysAnnotations, 0);
    for (int i = 1; i < matrix[1].length; i++) {
      matrix[1][i] = getDescriptionForComment(matrix[1][i].charAt(0), false, false);
    }
    keysMarkers = HashVec.getKeys(markerAnnotations);
    for (int i = 0; i < list.length; i++) {
      annotationsVector = markerAnnotations.get(keysMarkers[i]);
      annotationsOfTheMarker = "";
      matrix[i + 2] = Array.stringArray(keysAnnotations.length + 1, "0");
      matrix[i + 2][0] = keysMarkers[i];
      for (int j = 0; j < annotationsVector.size(); j++) {
        annotationsOfTheMarker +=
            (j == 0 ? "" : ";") + commentsHash.get(annotationsVector.elementAt(j).toCharArray()[0]);
        matrix[i + 2][ext.indexOfStr(annotationsVector.elementAt(j), keysAnnotations) + 1] = "1";
      }
      list[i] = keysMarkers[i] + "\t" + annotationsOfTheMarker;
    }
    Files.writeList(list, outputDir + "annotations_list.xln");
    Files.writeMatrix(matrix, outputDir + "annotations_matrix.xln", "\t");
  }

  @Override
  public void exportToText(Project proj, String outputFile) {
    PrintWriter writer;
    char[] keys;

    keys = getKeys();
    try {
      writer = new PrintWriter(new FileWriter(outputFile));
      for (char key : keys) {
        writer.println(key + "\t" + getDescriptionForComment(key, false, false));
      }
      writer.close();
      proj.getLog().report("Finished exporting annotation list to " + outputFile);
    } catch (Exception e) {
      System.err.println("Error writing to " + outputFile);
      e.printStackTrace();
    }
  }


  public String getDescriptionForComment(char c, boolean includeShortcuts, boolean includeNumbers) {
    return (includeShortcuts ? "'" + c + "' " : "") + commentsHash.get(c)
        + (includeNumbers ? " (n=" + annotationMarkerLists.get(c + "").size() + ")" : "");
  }


  public char[] getKeys() {
    String[] commentDescriptions;
    char[] result;
    Enumeration<Character> keys;
    char key;
    Hashtable<String, Character> temp;

    temp = new Hashtable<String, Character>(commentsHash.size() + 1, (float) 1.00);
    keys = commentsHash.keys();
    commentDescriptions = new String[commentsHash.size()];
    for (int i = 0; i < commentDescriptions.length; i++) {
      key = keys.nextElement();
      commentDescriptions[i] = commentsHash.get(key);
      temp.put(commentDescriptions[i], key);
    }
    commentDescriptions = Sort.putInOrder(commentDescriptions);

    result = new char[commentsHash.size()];
    for (int i = 0; i < commentDescriptions.length; i++) {
      result[i] = temp.get(commentDescriptions[i]);
    }
    return result;
  }

  public String[] getMarkerLists() {
    String[] result;
    Enumeration<String> keys;

    result = new String[markerAnnotations.size()];
    keys = markerAnnotations.keys();
    for (int i = 0; i < result.length; i++) {
      result[i] = keys.nextElement();
    }
    return result;
  }

  public String[] getMarkerLists(char key) {
    return Array.toStringArray(annotationMarkerLists.get(key + ""));
  }

  public void importList(String importList, Logger log) {
    String[][] newMappings;
    char[] keys;

    keys = getKeys();
    newMappings = HashVec.loadFileToStringMatrix(importList, false, new int[] {0, 1},
        Files.determineDelimiter(importList, log), false, keys.length, false);

    for (String[] newMapping : newMappings) {
      if (commentsHash.containsKey(newMapping[0].charAt(0))) {
        renameAnnotation(newMapping[0].charAt(0), newMapping[1]);
      } else {
        log.reportError("Warning - could not find an annotation using the shortcut character '"
            + newMapping[0] + "'; ignoring");
      }
    }
    log.report("Finished importing annotation list from " + importList);
  }

  public boolean markerHasAnnotation(String markerName, char annotation) {
    if (markerAnnotations.containsKey(markerName)) {
      return markerAnnotations.get(markerName).contains(annotation + "");
    } else {
      return false;
    }
  }

  public boolean markerHasAnyAnnotation(String markerName) {
    return markerAnnotations.containsKey(markerName);
  }

  public void removeAnnotation(Project proj, char c) {
    int response;
    Vector<String> markers;

    response = JOptionPane.showConfirmDialog(null,
        "This will remove the annotation '" + commentsHash.get(c) + "' from all markers (n="
            + annotationMarkerLists.get(c + "").size() + ") from the annotation database",
        "Warning", JOptionPane.ERROR_MESSAGE);
    if (response == 0) {
      serialize(proj.BACKUP_DIRECTORY.getValue(true, true) + "annotationsBeforeRemoving_"
          + ext.replaceWithLinuxSafeCharacters(commentsHash.get(c), true) + ".ser."
          + (new SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date())));
      commentsHash.remove(c);
      markers = annotationMarkerLists.get(c + "");
      for (int i = 0; markers != null && i < markers.size(); i++) {
        removeAnnotationForMarker(markers.elementAt(i), c);
      }
      System.out.println(ext.listWithCommas(HashVec.getKeys(annotationMarkerLists, true, false)));
      annotationMarkerLists.remove(c + "");
      System.out.println(ext.listWithCommas(HashVec.getKeys(annotationMarkerLists, true, false)));
      System.out.println();
    }

  }

  public void removeAnnotationForMarker(String markerName, char c) {
    if (markerAnnotations.containsKey(markerName)) {
      if (markerAnnotations.get(markerName).contains(c + "")) {
        markerAnnotations.get(markerName).remove(c + "");
        annotationMarkerLists.get(c + "").remove(markerName);
      } else {
        System.err.println(
            "Error - cannot remove " + c + " from " + markerName + " if it's not already checked");
      }
      if (markerAnnotations.get(markerName).size() == 0) {
        markerAnnotations.remove(markerName);
      }
    } else {
      System.err.println("Error - cannot remove " + c + " from " + markerName
          + " if the marker does not have any annotation");
    }
  }

  public void renameAnnotation(char c, String newComment) {
    commentsHash.put(c, newComment);
  }

  public void renameAnnotationKey(char oldShortcut, char newShortcut) {
    // TODO
  }

  public void serialize(String filename) {
    SerializedFiles.writeSerial(this, filename);
  }
}
