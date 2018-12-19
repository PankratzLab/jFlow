package org.genvisis.one.ben.imagetag;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.Map.Entry;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.filesys.Positions;

public class Annotator implements IAnnotator {

  private ArrayList<String> rootKeys = new ArrayList<>();
  private HashMap<String, HashMap<String, AnnotatedImage>> imageMap = new HashMap<>();
  private ArrayList<AnnotatedImage.Annotation> annotations = new ArrayList<>();

  private static final HashMap<String, Character> mnemonicMap = new HashMap<>();

  public static Character determineMnemonic(String ann) {
    Character mnemonic;
    if (mnemonicMap.containsKey(ann)) {
      mnemonic = mnemonicMap.get(ann);
    } else {
      int ind = 0;
      Character c = ann.toUpperCase().charAt(ind);
      while (mnemonicMap.containsValue(c)) {
        ind++;
        c = ann.toUpperCase().charAt(ind);
      }
      if (ind == ann.length()) {
        System.err.println("Error - all possible mnemonic characters already used for annotation {"
                           + ann + "}.  Using alphanumerics instead.");
        String alphanum = "abcdefghijklmnopqrstuvwxyz0123456789";
        ind = 0;
        while (mnemonicMap.containsValue(alphanum.charAt(ind))) {
          ind++;
          if (ind == alphanum.length()) {
            System.err.println("Error - ran out of alphanumeric mnemonics too!");
            break;
          }
        }
        if (ind < alphanum.length()) {
          mnemonic = alphanum.charAt(ind);
          mnemonicMap.put(ann, mnemonic);
        } else {
          mnemonic = '-';
        }
      } else {
        mnemonicMap.put(ann, c);
        mnemonic = c;
      }
    }
    return mnemonic;
  }

  public ArrayList<String> getRoots() {
    return new ArrayList<>(rootKeys);
  }

  public HashMap<String, HashMap<String, AnnotatedImage>> getAnnotationMap() {
    return imageMap;
  }

  @Override
  public ArrayList<AnnotatedImage.Annotation> getAnnotations() {
    return this.annotations;
  }

  @Override
  public void addNewAnnotation(AnnotatedImage.Annotation newAnnotation) {
    this.annotations.add(newAnnotation);
  }

  @Override
  public void replaceAnnotation(AnnotatedImage.Annotation prevAnnot,
                                AnnotatedImage.Annotation newAnnot) {
    this.annotations.set(this.annotations.indexOf(prevAnnot), newAnnot);
    for (HashMap<String, AnnotatedImage> annMap : imageMap.values()) {
      for (AnnotatedImage ai : annMap.values()) {
        if (ai.getAnnotations().contains(prevAnnot)) {
          ai.getAnnotations().set(ai.getAnnotations().indexOf(prevAnnot), newAnnot);
        }
      }
    }
  }

  @Override
  public void deleteAnnotation(AnnotatedImage.Annotation annotation) {
    this.annotations.remove(annotation);
    for (HashMap<String, AnnotatedImage> annMap : imageMap.values()) {
      for (AnnotatedImage ai : annMap.values()) {
        ai.getAnnotations().remove(annotation);
      }
    }
    mnemonicMap.values().remove(annotation.mnemonic);
  }

  @Override
  public void loadImgDir(String dir) {
    File dFil = new File(dir);
    File[] subDirs = dFil.listFiles(new FileFilter() {

      @Override
      public boolean accept(File pathname) {
        return pathname.isDirectory();
      }
    });
    for (File d : subDirs) {
      String rootIdent = d.getName();
      rootKeys.add(rootIdent);
      HashMap<String, AnnotatedImage> imgs = new HashMap<>();
      imageMap.put(rootIdent, imgs);
      String[] imgFiles = d.list(new FilenameFilter() {

        @Override
        public boolean accept(File dir, String name) {
          return name.startsWith(rootIdent) && (name.endsWith(".png") || name.endsWith(".jpg"));
        }
      });
      for (String img : imgFiles) {
        AnnotatedImage ai = parseAnnotatedImage(ext.verifyDirFormat(d.getAbsolutePath()), rootIdent,
                                                img);
        imgs.put(img, ai);
      }
    }
  }

  private static AnnotatedImage parseAnnotatedImage(String dirPath, String rootID,
                                                    String imgFileName) {
    String f = imgFileName.startsWith("_") ? imgFileName.substring(1) : imgFileName;
    f = ext.rootOf(f);
    String[] parts = (f.startsWith(rootID) ? f.substring(rootID.length() + 1) : f).split("_", -1);
    int chr = Positions.chromosomeNumber(parts[1]);
    int start = Integer.parseInt(parts[2]);
    int stop = Integer.parseInt(parts[3]);
    String name = Positions.getUCSCformat(new int[] {chr, start, stop}) + " ~ " + parts[0];
    AnnotatedImage ai = new AnnotatedImage(name, false);
    ai.setImageFile(dirPath + imgFileName);
    return ai;
  }

  private static final String ANNOT_TOKEN = "@ANNOT";
  private static final String IMAGE_TOKEN = "@IMAGE";

  @Override
  public void loadAnnotations(String annotFile) throws IOException {
    BufferedReader reader = Files.getAppropriateReader(annotFile);
    String line = null;

    while ((line = reader.readLine()) != null) {
      if ("".equals(line)) continue;
      if (line.startsWith(ANNOT_TOKEN)) {
        String[] pts = line.split("\t");
        AnnotatedImage.Annotation a = new AnnotatedImage.Annotation(pts[1], pts[2]);
        this.annotations.add(a);
      } else if (line.startsWith(IMAGE_TOKEN)) {
        String[] pts = line.split("\t")[1].split("\\|", -1);
        AnnotatedImage ai = new AnnotatedImage(pts[0], false);
        String imgFile = pts[1].equals("") ? null : pts[1];
        ai.setImageFile(imgFile);
        ai.setMissing(imgFile == null || (!imgFile.contains(";") && !Files.exists(imgFile)));
        if (pts.length > 2) {
          for (int i = 2; i < pts.length; i++) {
            for (AnnotatedImage.Annotation a : this.annotations) {
              if (a.annotation.equals(pts[i])) {
                ai.getAnnotations().add(a);
              }
            }
          }
        }
        String rootKey = parseRootKey(imgFile);
        if (!rootKeys.contains(rootKey)) {
          rootKeys.add(rootKey);
        }
        HashMap<String, AnnotatedImage> map = imageMap.get(rootKey);
        if (map == null) {
          map = new HashMap<>();
          imageMap.put(rootKey, map);
        }
        map.put(ai.getName(), ai);
      }
    }
    reader.close();
  }

  private String parseRootKey(String fullPathToFile) {
    // top level directory, for this solution
    String path = ext.rootOf(fullPathToFile, false);
    int i = -1;
    if ((i = path.lastIndexOf('/')) >= 0) {
      path = path.substring(0, i);
      i = path.lastIndexOf('/');
      path = path.substring(i + 1);
    } else if ((i = path.lastIndexOf('\\')) >= 0) {
      path = path.substring(0, i);
      i = path.lastIndexOf('\\');
      path = path.substring(i + 1);
    }
    return path;
  }

  public void saveAnnotation(AnnotatedImage.Annotation annot, String file) {
    backupExistingFile(file);
    PrintWriter writer = Files.getAppropriateWriter(file);
    writer.println(ANNOT_TOKEN + "=" + annot.annotation);
    writer.println();
    for (String fcs : rootKeys) {
      for (Entry<String, AnnotatedImage> ent : imageMap.get(fcs).entrySet()) {
        if (ent.getValue().getAnnotations().contains(annot)) {
          writer.println(fcs + "\t" + ent.getKey());
        }
      }
    }
    writer.flush();
    writer.close();
  }

  @Override
  public void saveAnnotations(String annotFile) {
    HashMap<AnnotatedImage.Annotation, ArrayList<String>> map = new HashMap<>();
    for (AnnotatedImage.Annotation a : annotations) {
      map.put(a, new ArrayList<>());
    }
    for (String f : rootKeys) {
      for (Entry<String, AnnotatedImage> ent : imageMap.get(f).entrySet()) {
        for (AnnotatedImage.Annotation a : ent.getValue().getAnnotations()) {
          map.get(a).add(ent.getKey());
        }
      }
    }
    backupExistingFile(annotFile);
    PrintWriter writer = Files.getAppropriateWriter(annotFile);
    for (AnnotatedImage.Annotation a : map.keySet()) {
      StringBuilder sb = new StringBuilder(ANNOT_TOKEN).append("\t").append(a.annotation)
                                                       .append("\t").append(a.mnemonic);
      writer.println(sb.toString());
    }
    writer.println();
    for (HashMap<String, AnnotatedImage> annMap : imageMap.values()) {
      for (Entry<String, AnnotatedImage> ent : annMap.entrySet()) {
        writer.println(IMAGE_TOKEN + "\t" + ent.getValue().exportToString());
      }
    }
    writer.flush();
    writer.close();
  }

  @Override
  public void backupExistingFile(String file) {
    if (file == null || "".equals(file) || !Files.exists(file)) {
      return;
    }
    String root = ext.rootOf(file, false);
    String date = new SimpleDateFormat("yyyyMMdd").format(new Date());
    String exten = file.substring(root.length());

    String newFile = root + "_" + date + exten;
    Files.copyFile(file, newFile);
  }
}
