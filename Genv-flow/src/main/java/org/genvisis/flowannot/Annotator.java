package org.genvisis.flowannot;

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
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.compress.archivers.tar.TarArchiveEntry;
import org.apache.commons.compress.archivers.tar.TarArchiveInputStream;
import org.genvisis.fcs.auto.Panel;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.ImagesPkl;
import org.pankratzlab.common.ext;

public class Annotator implements IAnnotator {

  private ArrayList<String> fcsKeys = new ArrayList<>();
  private HashMap<String, HashMap<String, AnnotatedImage>> imageMap = new HashMap<>();
  private ArrayList<AnnotatedImage.Annotation> annotations = new ArrayList<>();

  public ArrayList<String> getFCSKeys() {
    return getFCSKeys(null);
  }

  public ArrayList<String> getFCSKeys(Panel panel) {
    ArrayList<String> keys = fcsKeys;
    if (panel != null) {
      keys =
          new ArrayList<>(
              fcsKeys.stream().filter(p -> panel.isPanel(p)).collect(Collectors.toList()));
    }
    return keys;
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
  public void replaceAnnotation(
      AnnotatedImage.Annotation prevAnnot, AnnotatedImage.Annotation newAnnot) {
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
  }

  private void loadFromTAR(String tar, List<Panel> panels) throws IOException {
    TarArchiveInputStream tarStream = (TarArchiveInputStream) ImagesPkl.getTarInputStream(tar);
    if (tarStream == null) return;
    TarArchiveEntry entry;
    String dirName = null;
    List<String> dirFiles = new ArrayList<>();
    while ((entry = tarStream.getNextTarEntry()) != null) {
      if (entry.isDirectory()) {
        if (dirName != null && !dirFiles.isEmpty()) {
          checkSpecials(
              panels, tar + File.separatorChar + dirName, dirName, dirFiles.toArray(new String[0]));
        }
        String name = entry.getName();
        if (name.endsWith("/")) {
          name = name.substring(0, name.length() - 1);
        }
        dirName = name;
        dirFiles = new ArrayList<>();
        fcsKeys.add(name);
        HashMap<String, AnnotatedImage> fcsImgs = new HashMap<>();
        imageMap.put(name, fcsImgs);
      } else {
        String fcsName = null;
        String image = entry.getName();
        dirFiles.add(image.substring(dirName.length() + 1));
        int index;
        if ((index = image.lastIndexOf('/')) >= 0) {
          fcsName = image.substring(0, index);
          image = image.substring(index + 1);
        } else if ((index = image.lastIndexOf('\\')) >= 0) {
          fcsName = image.substring(0, index);
          image = image.substring(index + 1);
        } else if (image.startsWith(dirName)) {
          image = image.substring(dirName.length() + 1);
        }
        if (fcsName == null) continue;
        final String img = image;
        Optional<Panel> panel = panels.stream().filter(p -> p.isPanel(img)).findFirst();
        if (!panel.isPresent()) {
          System.err.println("Couldn't identify panel for file: " + img);
          continue;
        }
        String[][] gateTree = panel.get().getGateTree();

        int gateInd = -1;
        String name = img.substring(fcsName.length() + 1, img.length() - 4);
        for (int i = 0; i < gateTree.length; i++) {
          if (gateTree[i][0].equals(name)
              || ext.replaceWithLinuxSafeCharacters(gateTree[i][0]).equals(name)) {
            gateInd = i;
            break;
          }
        }
        if (gateInd >= 0) {
          AnnotatedImage ai = new AnnotatedImage(gateInd + "", gateInd == 0);
          ai.setImageFile(tar + File.separatorChar + entry.getName());
          ai.setGateName(gateTree[gateInd][0]);
          imageMap.get(fcsName).put(gateTree[gateInd][0], ai);
        }
      }
    }
  }

  @Override
  public void loadImgDir(String dir, List<Panel> panels) {
    if (dir.endsWith(".tar") || dir.endsWith(".tar.gz") || dir.endsWith(".tgz")) {
      try {
        loadFromTAR(dir, panels);
      } catch (IOException e) {
        e.printStackTrace();
      }
      return;
    }
    File dFil = new File(dir);
    File[] subDirs =
        dFil.listFiles(
            new FileFilter() {

              @Override
              public boolean accept(File pathname) {
                return pathname.isDirectory();
              }
            });
    for (File d : subDirs) {
      loadSubDir(d, panels);
    }
  }

  private void loadSubDir(File d, List<Panel> panels) {
    String fcsFilename = d.getName();
    fcsKeys.add(fcsFilename);
    HashMap<String, AnnotatedImage> fcsImgs = new HashMap<>();
    imageMap.put(fcsFilename, fcsImgs);
    String[] imgFiles =
        d.list(
            new FilenameFilter() {

              @Override
              public boolean accept(File dir, String name) {
                return name.startsWith(fcsFilename)
                    && (name.endsWith(".png") || name.endsWith(".jpg"));
              }
            });
    for (String img : imgFiles) {
      Optional<Panel> panel = panels.stream().filter(p -> p.isPanel(img)).findFirst();
      if (!panel.isPresent()) {
        System.err.println("Couldn't identify panel for file: " + img);
        continue;
      }
      String[][] gateTree = panel.get().getGateTree();

      int gateInd = -1;
      String name = img.substring(fcsFilename.length() + 1, img.length() - 4);
      for (int i = 0; i < gateTree.length; i++) {
        if (gateTree[i][0].equals(name)
            || ext.replaceWithLinuxSafeCharacters(gateTree[i][0]).equals(name)) {
          gateInd = i;
          break;
        }
      }
      if (gateInd >= 0) {
        AnnotatedImage ai = new AnnotatedImage(gateInd + "", gateTree[gateInd].length == 1);
        ai.setImageFile(ext.verifyDirFormat(d.getAbsolutePath()) + img);
        ai.setGateName(gateTree[gateInd][0]);
        fcsImgs.put(gateTree[gateInd][0], ai);
      }
    }
    checkSpecials(panels, ext.verifyDirFormat(d.getAbsolutePath()), fcsFilename, imgFiles);
  }

  private void checkSpecials(
      List<Panel> panels, String dir, String fcsFilename, String[] imgFiles) {
    for (Panel panel : panels) {
      if (panel.hasSpecials()) {
        String[][] gateTree = panel.getGateTree();
        String[] treeKeys = ArrayUtils.extract(gateTree, 0);
        Map<String, List<String>> specials = panel.getSpecials();
        for (Entry<String, List<String>> special : specials.entrySet()) {
          int gateInd = ext.indexOfStr(special.getKey(), treeKeys);
          if (gateInd >= 0) {
            AnnotatedImage ai = new AnnotatedImage(gateInd + "", gateInd == 0);
            StringBuilder allImgs = new StringBuilder();
            for (String sub : special.getValue()) {
              String subImg = null;
              for (String img : imgFiles) {
                String name = img.substring(fcsFilename.length() + 1, img.length() - 4);
                if (sub.equals(name) || ext.replaceWithLinuxSafeCharacters(sub).equals(name)) {
                  subImg = img;
                  break;
                }
              }
              if (subImg != null) {
                allImgs.append(
                    (allImgs.length() > 0 ? ";" : "") + dir + File.separatorChar + subImg);
              }
            }
            ai.setImageFile(allImgs.toString());
            ai.setGateName(gateTree[gateInd][0]);
            imageMap.get(fcsFilename).put(gateTree[gateInd][0], ai);
          }
        }
      }
    }
  }

  private static final String ANNOT_TOKEN = "@ANNOT";
  private static final String IMAGE_TOKEN = "@IMAGE";

  @Override
  public void loadAnnotations(String annotFile, List<Panel> panels) throws IOException {
    BufferedReader reader = Files.getAppropriateReader(annotFile);
    String line = null;

    Set<String> rootIdents =
        panels.stream().map(Panel::getGateTree).map(g -> g[0][0]).collect(Collectors.toSet());

    while ((line = reader.readLine()) != null) {
      if ("".equals(line)) continue;
      if (line.startsWith(ANNOT_TOKEN)) {
        String[] pts = line.split("\t");
        AnnotatedImage.Annotation a = new AnnotatedImage.Annotation(pts[1], pts[2]);
        this.annotations.add(a);
      } else if (line.startsWith(IMAGE_TOKEN)) {
        String[] pts = line.split("\t")[1].split("\\|", -1);
        AnnotatedImage ai = new AnnotatedImage(pts[0], rootIdents.contains(pts[0]));
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
        String file =
            ext.removeDirectoryInfo(imgFile.contains(";") ? imgFile.split(";")[0] : imgFile);
        String fcsFile = file.substring(0, file.indexOf(".fcs.") + 4);
        if (!fcsKeys.contains(fcsFile)) {
          fcsKeys.add(fcsFile);
        }
        HashMap<String, AnnotatedImage> map = imageMap.get(fcsFile);
        if (map == null) {
          map = new HashMap<>();
          imageMap.put(fcsFile, map);
        }
        map.put(ai.getGateName(), ai);
      }
    }
    reader.close();
  }

  public void saveAnnotation(AnnotatedImage.Annotation annot, String file) {
    backupExistingFile(file);
    PrintWriter writer = Files.getAppropriateWriter(file);
    writer.println(ANNOT_TOKEN + "=" + annot.annotation);
    writer.println();
    for (String fcs : fcsKeys) {
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
    for (String f : fcsKeys) {
      for (Entry<String, AnnotatedImage> ent : imageMap.get(f).entrySet()) {
        for (AnnotatedImage.Annotation a : ent.getValue().getAnnotations()) {
          map.get(a).add(ent.getKey());
        }
      }
    }
    backupExistingFile(annotFile);
    PrintWriter writer = Files.getAppropriateWriter(annotFile);
    for (AnnotatedImage.Annotation a : map.keySet()) {
      StringBuilder sb =
          new StringBuilder(ANNOT_TOKEN)
              .append("\t")
              .append(a.annotation)
              .append("\t")
              .append(a.mnemonic);
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
