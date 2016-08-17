package org.genvisis.cnv.manage;

import java.io.BufferedReader;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Positions;
import org.genvisis.common.ext;

public class PlinkMergePrep {

  private static final int INVALID_ROOT = 0;
  public static final int PEDMAP = 1;
  public static final int BEDBIMFAM = 2;
  public static final String TEMP_MKR_FILE = "./tempMarkersForPlinkMerge.txt";;

  static int isValidRoot(final String plinkRootWithDir) {
    String[] files =
        (new File(ext.parseDirectoryOfFile(plinkRootWithDir + ".bim"))).list(new FilenameFilter() {
          @Override
          public boolean accept(File dir, String name) {
            return name.startsWith(ext.rootOf(plinkRootWithDir));
          }
        });
    boolean foundBed = false;
    boolean foundBim = false;
    boolean foundFam = false;
    boolean foundPed = false;
    boolean foundMap = false;
    for (String f : files) {
      if (f.endsWith(".bim")) {
        foundBim = true;
      } else if (f.endsWith(".bed")) {
        foundBed = true;
      } else if (f.endsWith(".fam")) {
        foundFam = true;
      } else if (f.endsWith(".ped")) {
        foundPed = true;
      } else if (f.endsWith(".map")) {
        foundMap = true;
      }
    }

    if (foundBed && foundBim && foundFam) {
      return BEDBIMFAM;
    } else if (foundPed && foundMap) {
      return PEDMAP;
    }
    return INVALID_ROOT;
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String outRoot = "merge";
    int outType = BEDBIMFAM;
    boolean overwrite = false;
    boolean renameMarkers = false;
    String dir = "./";
    String rgnFile = "";
    String mkrFile = "";
    String[] roots = null;
    String[] prepends = null;

    String usage = "\n" + "cnv.manage.PlinkMergePrep requires 0-1 arguments\n"
        + "   (1) Output Root (i.e. outRoot=" + outRoot + " (default))\n"
        + "   (2) Output Type (i.e. outType=" + outType + " (default, options are (" + PEDMAP
        + ") PED/MAP, or (" + BEDBIMFAM + ") BED/BIM/FAM.))\n"
        + "   (3a) Region file (UCSC-format) (i.e. regions=" + rgnFile + " (default))\n"
        + "   (3b) Markers file (i.e. markers=" + mkrFile + " (default))\n"
        + "   (4) List of Plink fileRoots (either relative to 'dir' argument, or absolute paths) to merge (i.e. roots=plink1,plink2,plink4 (not the default))\n"
        + "   (5) overwrite flag (i.e. -overwrite (not the default))\n"
        + "   (6) rename markers flag (i.e. -renameMarkers (not the default))\n"
        + "   Optional: (7) directory of any relative files (marker list, regions list, or plink files) (i.e. dir="
        + dir + " (default))\n" + "";

    if (numArgs == 0) {
      System.err.println(usage);
      System.exit(1);
    }
    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("outRoot=")) {
        outRoot = ext.parseStringArg(arg, "merge");
        numArgs--;
      } else if (arg.startsWith("outType=")) {
        outType = ext.parseIntArg(arg);
        numArgs--;
      } else if (arg.startsWith("dir=")) {
        dir = ext.parseStringArg(arg, "./");
        numArgs--;
      } else if (arg.startsWith("roots=")) {
        roots = ext.parseStringArg(arg, "").split(",");
        numArgs--;
      } else if (arg.startsWith("regions=")) {
        rgnFile = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith("markers=")) {
        mkrFile = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith("-overwrite")) {
        overwrite = true;
        numArgs--;
      } else if (arg.startsWith("-renameMarkers")) {
        renameMarkers = true;
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
      if (!"./".equals(dir)) {
        dir = ext.verifyDirFormat(dir);
        if (!"".equals(rgnFile) && Files.isRelativePath(rgnFile)) {
          rgnFile = dir + rgnFile;
        }
        if (!"".equals(mkrFile) && Files.isRelativePath(mkrFile)) {
          mkrFile = dir + mkrFile;
        }
        for (int i = 0; i < roots.length; i++) {
          if (Files.isRelativePath(roots[i])) {
            roots[i] = dir + roots[i];
          }
        }
      }
      if (renameMarkers) {
        prepends = new String[roots.length];
        for (int i = 0; i < roots.length; i++) {
          prepends[i] = org.genvisis.common.ext.rootOf(roots[i], true);
        }
      }

      String mergeCommand =
          merge(outType, outRoot, overwrite, renameMarkers, rgnFile, mkrFile, roots, prepends);
      System.out.println(mergeCommand);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  static String merge(int outType, String outDirAndRoot, boolean overwrite, boolean renameMarkers,
      String regionFile, String markersFile, String plinkRootWithDir1, String plinkRootWithDir2,
      String prepend1, String prepend2) {
    int rootType1, rootType2;
    StringBuilder cmds;
    String tempRgnFile;

    if (outType != PEDMAP && outType != BEDBIMFAM) {
      System.err.println("Error - output type must be either (" + PEDMAP + ") PED/MAP, or ("
          + BEDBIMFAM + ") BED/BIM/FAM.");
      return null;
    }
    rootType1 = isValidRoot(plinkRootWithDir1);
    if (rootType1 == INVALID_ROOT) {
      System.err.println(
          "Error - couldn't find a complete set of PED/MAP or BED/BIM/FAM files for plink root '"
              + plinkRootWithDir1 + "'.");
      return null;
    }
    rootType2 = isValidRoot(plinkRootWithDir2);
    if (rootType2 == INVALID_ROOT) {
      System.err.println(
          "Error - couldn't find a complete set of PED/MAP or BED/BIM/FAM files for plink root '"
              + plinkRootWithDir2 + "'.");
      return null;
    }
    if (Files.exists(outDirAndRoot + (outType == PEDMAP ? ".ped" : ".bed"))) {
      if (!overwrite) {
        System.err.println("Error - PLINK files with given output root ('" + outDirAndRoot
            + "') already exist!  Please specify alternate output name or set overwrite flag.");
        return null;
      }
    }
    boolean hasRgns = false;
    if (regionFile != null && !regionFile.equals("")) {
      if (!Files.exists(regionFile)) {
        System.err.println("Error - specified region file ('" + regionFile + "') doesn't exist!");
        return null;
      } else {
        hasRgns = true;
      }
    }
    boolean hasMkrs = false;
    if (markersFile != null && !markersFile.equals("")) {
      if (!Files.exists(markersFile)) {
        System.err.println("Error - specified markers file ('" + markersFile + "') doesn't exist!");
        return null;
      } else {
        hasMkrs = true;
      }
    }
    if (hasRgns && hasMkrs) {
      System.err.println("Error - cannot specify both a regions file and a markers file.");
      return null;
    }
    if (renameMarkers) {
      renameMarkers(plinkRootWithDir1, rootType1, prepend1);
      renameMarkers(plinkRootWithDir2, rootType2, prepend2);
    }

    cmds = new StringBuilder();
    cmds.append("plink2 --noweb --allow-no-sex --").append(rootType1 == PEDMAP ? "file " : "bfile ")
        .append(plinkRootWithDir1);
    if (rootType2 == PEDMAP) {
      cmds.append(" --merge ").append(plinkRootWithDir2).append(".ped ").append(plinkRootWithDir2)
          .append(".map ");
    } else if (rootType2 == BEDBIMFAM) {
      cmds.append(" --bmerge ").append(plinkRootWithDir2).append(".bed ").append(plinkRootWithDir2)
          .append(".bim ").append(plinkRootWithDir2).append(".fam ");
    }
    cmds.append("--").append(outType == PEDMAP ? "recode" : "make-bed");
    if (regionFile != null && !regionFile.equals("") && Files.exists(regionFile)) {
      tempRgnFile = writeTempRegionFileFromUCSCFile(regionFile);
      cmds.append(" --extract range ").append(tempRgnFile);
    }
    if (markersFile != null && !markersFile.equals("") && Files.exists(markersFile)) {
      String mkrFile = markersFile;
      if (renameMarkers) {
        String[] markers =
            HashVec.loadFileToStringArray(mkrFile, false, false, new int[] {0}, true, false, "\t");
        ArrayList<String> newMkrs = new ArrayList<String>();
        String root = ext.rootOf(plinkRootWithDir1, true);
        for (String marker : markers) {
          newMkrs.add(root + "_" + marker);
        }
        root = ext.rootOf(plinkRootWithDir2, true);
        for (String marker : markers) {
          newMkrs.add(root + "_" + marker);
        }
        for (String marker : markers) { // and add all original markers too
          newMkrs.add(marker);
        }
        Files.writeArrayList(newMkrs, TEMP_MKR_FILE);
        mkrFile = (new File(TEMP_MKR_FILE)).getAbsolutePath();
      }
      cmds.append(" --extract ").append(mkrFile);
    }
    cmds.append(" --out ").append(outDirAndRoot);
    return cmds.toString();
  }

  public static String merge(int outType, String outDirAndRoot, boolean overwrite,
      boolean renameMarkers, String regionFile, String markersFile, String[] plinkRootsWithDirs,
      String[] prepends) {
    int[] rootTypes;
    String[] lines;
    String fileListName, tempRgnFile;
    StringBuilder cmds;
    int index;

    if (outType != PEDMAP && outType != BEDBIMFAM) {
      System.err.println("Error - output type must be either (" + PEDMAP + ") PED/MAP, or ("
          + BEDBIMFAM + ") BED/BIM/FAM.");
      return null;
    }
    if (plinkRootsWithDirs.length < 2) {
      System.err
          .println("Error - at least two distinct PLINK file roots must be specified for merging.");
      return null;
    } else if (plinkRootsWithDirs.length == 2) {
      return merge(outType, outDirAndRoot, overwrite, renameMarkers, regionFile, markersFile,
          plinkRootsWithDirs[0], plinkRootsWithDirs[1], prepends[0], prepends[1]);
    }
    if (Files.exists(outDirAndRoot + (outType == PEDMAP ? ".ped" : ".bed"))) {
      if (!overwrite) {
        System.err.println("Error - PLINK files with given output root ('" + outDirAndRoot
            + "') already exist!  Please specify alternate output name or set overwrite flag.");
        return null;
      }
    }
    boolean hasRgns = false;
    if (regionFile != null && !regionFile.equals("")) {
      if (!Files.exists(regionFile)) {
        System.err.println("Error - specified region file ('" + regionFile + "') doesn't exist!");
        return null;
      } else {
        hasRgns = true;
      }
    }
    boolean hasMkrs = false;
    if (markersFile != null && !markersFile.equals("")) {
      if (!Files.exists(markersFile)) {
        System.err.println("Error - specified markers file ('" + markersFile + "') doesn't exist!");
        return null;
      } else {
        hasMkrs = true;
      }
    }
    if (hasRgns && hasMkrs) {
      System.err.println("Error - cannot specify both a regions file and a markers file.");
      return null;
    }

    rootTypes = new int[plinkRootsWithDirs.length];
    for (int i = 0; i < plinkRootsWithDirs.length; i++) {
      rootTypes[i] = isValidRoot(plinkRootsWithDirs[i]);
      if (rootTypes[i] == INVALID_ROOT) {
        System.err.println(
            "Error - couldn't find a complete set of PED/MAP or BED/BIM/FAM files for plink root '"
                + plinkRootsWithDirs[i] + "'.");
        return null;
      }
    }
    lines = new String[plinkRootsWithDirs.length - 1];
    for (int i = 1; i < plinkRootsWithDirs.length; i++) {
      lines[i - 1] =
          rootTypes[i] == PEDMAP ? plinkRootsWithDirs[i] + ".ped " + plinkRootsWithDirs[i] + ".map"
              : plinkRootsWithDirs[i] + ".bed " + plinkRootsWithDirs[i] + ".bim "
                  + plinkRootsWithDirs[i] + ".fam";
    }
    fileListName = "./fileList.txt";
    index = 1;
    while (Files.exists(fileListName)) {
      fileListName = "./fileList_" + index++ + ".txt";
    }
    Files.writeList(lines, fileListName);

    if (renameMarkers) {
      for (int i = 0; i < plinkRootsWithDirs.length; i++) {
        renameMarkers(plinkRootsWithDirs[i], rootTypes[i], prepends[i]);
      }
    }
    fileListName = (new File(fileListName)).getAbsolutePath();
    cmds = new StringBuilder();
    cmds.append("plink2 --noweb --allow-no-sex --")
        .append(rootTypes[0] == PEDMAP ? "file " : "bfile ").append(plinkRootsWithDirs[0])
        .append(" --merge-list ").append(fileListName);
    cmds.append(" --").append(outType == PEDMAP ? "recode" : "make-bed");
    if (regionFile != null && !regionFile.equals("") && Files.exists(regionFile)) {
      tempRgnFile = writeTempRegionFileFromUCSCFile(regionFile);
      cmds.append(" --extract range ").append(tempRgnFile);
    }
    if (markersFile != null && !markersFile.equals("") && Files.exists(markersFile)) {
      String mkrFile = markersFile;
      if (renameMarkers) {
        String[] markers =
            HashVec.loadFileToStringArray(mkrFile, false, false, new int[] {0}, true, false, "\t");
        ArrayList<String> newMkrs = new ArrayList<String>();
        for (String plinkRootsWithDir : plinkRootsWithDirs) { // add all newly-renamed markers
          String root = ext.rootOf(plinkRootsWithDir, true);
          for (String marker : markers) {
            newMkrs.add(root + "_" + marker);
          }
        }
        for (String marker : markers) { // and add all original markers too
          newMkrs.add(marker);
        }
        Files.writeArrayList(newMkrs, TEMP_MKR_FILE);
        mkrFile = (new File(TEMP_MKR_FILE)).getAbsolutePath();
      }
      cmds.append(" --extract ").append(mkrFile);
    }
    cmds.append(" --out ").append(outDirAndRoot);

    return cmds.toString();
  }

  static void renameMarkers(String plinkDirRoot, int rootType, String prepend) {
    BufferedReader reader;
    PrintWriter writer;
    String line, ext, newFile;
    String[] parts;

    ext = rootType == PEDMAP ? ".map" : ".bim";
    newFile = plinkDirRoot + "_orig" + ext;
    if (Files.exists(newFile)) {
      throw new RuntimeException("Error - renamed PLINK .bim files already exist.");
    }
    boolean moved = (new File(plinkDirRoot + ext)).renameTo(new File(newFile));
    if (moved) {
      try {
        reader = Files.getAppropriateReader(plinkDirRoot + "_orig" + ext);
        writer = Files.getAppropriateWriter(plinkDirRoot + ext);
        while ((line = reader.readLine()) != null) {
          parts = line.split("[\\s]+", -1);
          parts[1] = prepend + "_" + parts[1];
          writer.println(Array.toStr(parts, "\t"));
        }
        writer.flush();
        writer.close();
        reader.close();
      } catch (IOException e) {
        e.printStackTrace();
      }
    } else {
      System.err.println("Error - unable to move " + plinkDirRoot + ext + " to " + plinkDirRoot
          + "_orig" + ext + ".");
    }
  }

  static String writeTempRegionFileFromUCSCFile(String regionFile) {
    String outFile = "tempRegions.txt";
    int index = 1;
    while (Files.exists(ext.parseDirectoryOfFile(regionFile) + outFile)) {
      outFile = "tempRegions_" + index++ + ".txt";
    }
    String[] ucscLines =
        HashVec.loadFileToStringArray(regionFile, false, false, new int[] {0}, true, false, "\t");
    String[] plinkLines = new String[ucscLines.length];
    for (int i = 0; i < ucscLines.length; i++) {
      int[] pos = Positions.parseUCSClocation(ucscLines[i]);
      plinkLines[i] = pos[0] + "\t" + pos[1] + "\t" + pos[2] + "\tPLINK_" + i;
    }
    Files.writeList(plinkLines, ext.parseDirectoryOfFile(regionFile) + outFile);
    return ext.verifyDirFormat(new File(ext.parseDirectoryOfFile(regionFile)).getAbsolutePath())
        + outFile;
  }

}
