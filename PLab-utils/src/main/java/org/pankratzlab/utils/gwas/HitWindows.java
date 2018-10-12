package org.pankratzlab.utils.gwas;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Hashtable;
import java.util.List;
import java.util.Vector;
import org.pankratzlab.common.Aliases;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.Matrix;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.Sort;
import org.pankratzlab.common.ext;
import org.pankratzlab.gwas.HitWindowsLD;
import org.pankratzlab.shared.filesys.Positions;
import org.pankratzlab.shared.filesys.Segment;
import org.pankratzlab.utils.filesys.SnpMarkerSet;

public class HitWindows {

  /**
   * @param filename
   * @param indexThreshold
   * @param windowMinSizePerSide
   * @param windowExtensionThreshold
   * @param additionalAnnotationVariableNames
   * @param log
   * @return String[][] of hits table, ready for output including header
   */
  public static String[][] determine(String filename, float indexThreshold,
                                     int windowMinSizePerSide, float windowExtensionThreshold,
                                     String[] additionalAnnotationVariableNames, Logger log) {
    BufferedReader reader;
    String[] line, header;
    int count;
    String[] markerNames;
    byte[] chrs;
    int[] positions;
    double[] pvals;
    int[] indices;
    String temp, delimiter;
    String[][] factors;
    String[][] annotation;

    if (additionalAnnotationVariableNames == null) {
      additionalAnnotationVariableNames = new String[0];
    }
    factors = new String[4 + additionalAnnotationVariableNames.length][];
    factors[0] = Aliases.MARKER_NAMES;
    factors[1] = Aliases.CHRS;
    factors[2] = Aliases.POSITIONS;
    factors[3] = Aliases.PVALUES;
    for (int i = 0; i < additionalAnnotationVariableNames.length; i++) {
      factors[4 + i] = new String[] {additionalAnnotationVariableNames[i]};
    }

    count = Files.countLines(filename, 1);
    // log.report("Parsing "+count+" lines");
    markerNames = new String[count];
    chrs = new byte[count];
    positions = new int[count];
    pvals = new double[count];
    annotation = new String[count][additionalAnnotationVariableNames.length];
    try {
      reader = Files.getAppropriateReader(filename);
      temp = reader.readLine();
      delimiter = ext.determineDelimiter(temp);
      header = temp.trim().split(delimiter);
      indices = ext.indexFactors(factors, header, false, false, true, true, log);
      if (ArrayUtils.min(indices) == -1) {
        log.reportError("Aborting after failing to find the appropriate column headers in file "
                        + filename);
        log.reportError("Missing:");
        for (int i = 0; i < indices.length; i++) {
          if (indices[i] == -1) {
            log.reportError("  " + ArrayUtils.toStr(factors[i], "/"));
          }
        }
        return null;
      }
      count = 0;
      // log.report("Parsing... "+filename);
      while ((temp = reader.readLine()) != null) {
        if (delimiter.equals(",")) {
          line = ext.splitCommasIntelligently(temp, true, log);
        } else {
          line = temp.trim().split(delimiter);
        }
        line = ArrayUtils.subArray(line, indices);
        try {
          markerNames[count] = line[0];
          chrs[count] = Positions.chromosomeNumber(line[1]);
          // positions[count] = Integer.parseInt(line[2]);
          positions[count] = (new java.math.BigDecimal(line[2])).intValueExact();
          pvals[count] = ext.isMissingValue(line[3]) ? 999 : Double.parseDouble(line[3]);
          for (int i = 0; i < additionalAnnotationVariableNames.length; i++) {
            annotation[count][i] = line[4 + i];
          }
        } catch (Exception e) {
          log.reportError("Error - reading file " + filename + " at line " + count + ": " + temp);
          log.reportError("  which was parsed as : " + ArrayUtils.toStr(line));
          log.reportException(e);
          return null;
        }
        count++;
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      log.reportError("Error: file \"" + filename + "\" not found in current directory");
      return null;
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + filename + "\"");
      return null;
    }

    return determineButOrderFirst(markerNames, chrs, positions, pvals, indexThreshold,
                                  windowMinSizePerSide, windowExtensionThreshold,
                                  additionalAnnotationVariableNames, annotation);
  }

  public static String[][] determineButOrderFirst(String[] markerNames, byte[] chrs,
                                                  int[] positions, double[] pvals,
                                                  float indexThreshold, int windowMinSizePerSide,
                                                  float windowExtensionThreshold,
                                                  String[] additionalAnnotationVariableNames,
                                                  String[][] annotation) {
    int[] order;

    order = Sort.getSort2DIndices(chrs, positions);

    markerNames = Sort.getOrdered(markerNames, order);
    chrs = Sort.getOrdered(chrs, order);
    positions = Sort.getOrdered(positions, order);
    pvals = Sort.getOrdered(pvals, order);
    annotation = Sort.getOrdered(annotation, order);

    return determine(markerNames, chrs, positions, pvals, indexThreshold, windowMinSizePerSide,
                     windowExtensionThreshold, additionalAnnotationVariableNames, annotation);
  }

  public static String[][] determine(String[] markerNames, byte[] chrs, int[] positions,
                                     double[] pvals, float indexThreshold, int windowMinSizePerSide,
                                     float windowExtensionThreshold,
                                     String[] additionalAnnotationVariableNames,
                                     String[][] annotation) {
    int startIndex, stopIndex, offset, minIndex;
    double minPval;
    int region;
    int numSig, numSuggestive;
    Vector<String[]> v;
    String[] line;

    v = new Vector<>();
    line = new String[] {"Region", "MarkerName", "Chr", "Position", "p-value", "Region+Window",
                         "RegionStart", "RegionStop", "NumSigMarkers", "NumSuggestiveMarkers",
                         "NumTotalMarkers", "SizeOfRegion"};
    for (String additionalAnnotationVariableName : additionalAnnotationVariableNames) {
      line = ArrayUtils.addStrToArray(additionalAnnotationVariableName, line);
    }
    v.add(line);

    startIndex = -1;
    stopIndex = -1;
    region = 1;
    for (int i = 0; i < markerNames.length; i++) {
      if (pvals[i] < indexThreshold) {
        startIndex = i;
        minIndex = i;
        minPval = pvals[i];
        offset = 0;
        numSig = numSuggestive = 1;
        while (startIndex - offset - 1 >= 0 && chrs[startIndex] == chrs[startIndex - offset - 1]
               && positions[startIndex]
                  // *2 required to ensure that there are no overlapping SNPs 500kb after last hit and 500kb before next hit is technically a 1M region that should be merged:
                  - windowMinSizePerSide * 2 <= positions[startIndex - offset - 1]) {
          offset++;
          if (pvals[startIndex - offset] < windowExtensionThreshold) {
            startIndex -= offset;
            offset = 0;
            numSuggestive++;
          }
        }

        stopIndex = i;
        offset = 0;
        while (stopIndex + offset + 1 < markerNames.length
               && chrs[stopIndex] == chrs[stopIndex + offset + 1]
               && positions[stopIndex]
                  + windowMinSizePerSide >= positions[stopIndex + offset + 1]) { // don't want the 2* here, though
          offset++;
          if (pvals[stopIndex + offset] < indexThreshold) {
            numSig++;
          }
          if (pvals[stopIndex + offset] < windowExtensionThreshold) {
            stopIndex += offset;
            offset = 0;
            numSuggestive++;
          }
          if (pvals[stopIndex] < minPval) {
            minIndex = stopIndex;
            minPval = pvals[stopIndex];
          }
        }

        line = new String[] {region + "", markerNames[minIndex], chrs[minIndex]
                                                                 + "",
                             positions[minIndex] + "", pvals[minIndex]
                                                       + "",
                             "chr" + chrs[startIndex] + ":" + Math.max(positions[startIndex]
                                                                       - windowMinSizePerSide, 1)
                                                             + "-"
                                                             + (positions[stopIndex]
                                                                + windowMinSizePerSide),
                             positions[startIndex] + "", positions[stopIndex] + "", numSig + "",
                             numSuggestive + "", (stopIndex - startIndex + 1) + "",
                             (positions[stopIndex] - positions[startIndex] + 1) + ""};
        for (int j = 0; j < additionalAnnotationVariableNames.length; j++) {
          line = ArrayUtils.addStrToArray(annotation[minIndex][j], line);
        }
        v.add(line);

        i = stopIndex + offset;
        region++;
      }
    }

    return Matrix.toStringArrays(v);
  }

  public static void generateHitsLookup(String inputHits, int window, String outputFile,
                                        String mapfile, Logger log) {
    BufferedReader reader;
    String[] line;
    Hashtable<String, Vector<String>> hash;
    Vector<String> v;
    Segment[][] segs;
    SnpMarkerSet markerSet;
    int[] indices;
    String[] header, traits, markerNames, chrPositions;
    Segment variant;

    v = new Vector<>();
    hash = new Hashtable<>();
    try {
      reader = new BufferedReader(new FileReader(inputHits));
      header = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
      indices = ext.indexFactors(new String[][] {{"Trait"}, Aliases.CHRS, Aliases.POSITIONS},
                                 header, false, true, true);
      if (!ArrayUtils.equals(indices, new int[] {0, 1, 2})) {
        log.reportError("Error - currently expecting format: Trait\tChr\tPosition");
      }
      while (reader.ready()) {
        line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
        HashVec.addIfAbsent(line[0], v);
        HashVec.addToHashVec(hash, line[0], line[1] + "\t" + line[1], false);
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      log.reportError("Error: file \"" + inputHits + "\" not found in current directory");
      return;
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + inputHits + "\"");
      return;
    }

    traits = ArrayUtils.toStringArray(v);
    segs = new Segment[traits.length][];
    for (int i = 0; i < traits.length; i++) {
      v = hash.get(traits[i]);
      segs[i] = new Segment[v.size()];
      for (int j = 0; j < segs[i].length; j++) {
        line = v.elementAt(j).split(PSF.Regex.GREEDY_WHITESPACE);
        segs[i][j] = new Segment(Positions.chromosomeNumber(line[0]),
                                 Integer.parseInt(line[1]) - window,
                                 Integer.parseInt(line[1]) + window);
      }
    }

    markerSet = new SnpMarkerSet(mapfile, false, log);
    markerNames = markerSet.getMarkerNames();
    chrPositions = markerSet.getChrAndPositions();
    for (int m = 0; m < 10; m++) {
      line = chrPositions[m].split(PSF.Regex.GREEDY_WHITESPACE);
      variant = new Segment(Byte.parseByte(line[0]), Integer.parseInt(line[1]),
                            Integer.parseInt(line[1]) + 1);
      v = new Vector<>();
      for (int i = 0; i < segs.length; i++) {
        for (int j = 0; j < segs[i].length; j++) {
          if (variant.overlaps(segs[i][j])) {
            v.add(traits[i]);
          }
        }
      }
      log.report(markerNames[m] + "\t"
                 + (v.size() == 0 ? "." : ArrayUtils.toStr(ArrayUtils.toStringArray(v), "/")));
    }
  }

  public static void fromParameters(String filename, Logger log) {
    List<String> params;

    params = Files.parseControlFile(filename, "hitWindows",
                                    new String[] {"# filename containing the markers/chrs/positions/p-values:",
                                                  "file=plink.assoc.linear",
                                                  "# outfile for the table", "out=table.out",
                                                  "# p-value threshold needed to be an index SNP",
                                                  "indexThresh=0.00000005",
                                                  "# minimum num bp to search per side of window",
                                                  "minWinSize=500000",
                                                  "# p-value threshold needed to extend a window",
                                                  "winThresh=0.000001",
                                                  "# optional additional column headers to include in the table",
                                                  "#annotationCols=betaCol,sdCol,mafCol",},
                                    log);

    if (params != null) {
      params.add("log=" + log.getFilename());
      main(ArrayUtils.toStringArray(params));
    }
  }

  public static void determineLD(String targetFile, String mapFile, String ldFile, int window,
                                 int column, double filter, String output, Logger log) {
    HitWindowsLD.determineLD(targetFile, mapFile, output, ldFile, window, filter, column, log);
  }

  public static void determineLD(String[] targets, String regionName, String mapFile, String output,
                                 String ldFile, int window, double filter, int column,
                                 boolean region, Logger log) {
    HitWindowsLD.determineLD(targets, regionName, mapFile, output, ldFile, window, filter, column,
                             region, log);
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = "input.dat";
    String outfile = "hits.out";
    float indexThreshold = (float) 0.00000005;
    int windowMinSizePerSide = 500000; // 500kb each side is technically a 1M window until the next
                                      // hit region, but we now take this into consideration in
                                      // the
                                      // main algorithm
    int windowMaxSizePerSide = 100000;
    float windowExtensionThreshold = (float) 0.00000005; // (float)0.00001;
    String knownHits = null;
    String ldFile = null;
    int column = 4;
    double ldFilter = 0.5;
    String map = "markers.dat";
    String[][] results;
    String[] annotationCols = null;
    Logger log;
    String logfile = null;

    String usage = "\n" + "gwas.HitWindows requires 0-1 arguments\n"
                   + "   (1) input filename (i.e. file=" + filename + " (default))\n"
                   + "   (2) filename of output (i.e. out=" + outfile
                   + " (default; set to null to print to stdout))\n"
                   + "   (3) p-value threshold for index SNPs (i.e. indexThresh=" + indexThreshold
                   + " (default))\n" + "   (4) minimum num bp per side of window (i.e. minWinSize="
                   + windowMinSizePerSide + " (default))\n"
                   + "   (5) p-value threshold to extend the window (i.e. winThresh="
                   + windowExtensionThreshold + " (default))\n"
                   + "   (6) (optional) additional column headers to include (i.e. annotationCols=betaCol,sdCol,mafCol (not the default))\n"
                   +

                   " OR\n"
                   + "   (1) list of known hits, 3 columns=trait+chr+position (i.e. knownHits=filenameOfKnownHits.dat (not the default))\n"
                   + "   (2) window around hit to extend (i.e. minWinSize=" + windowMinSizePerSide
                   + " (default))\n" + "   (3) map file for lookup (i.e. map=" + map
                   + " (default))\n" +

                   " OR\n"
                   + "   (1) input filename for LD-based hits (plink .flt, plink .ld , or Haploview export (i.e. ldFile="
                   + ldFile + " \n"
                   + "   (2) list of known hits, 1 column=SNP identifier (i.e. knownHits=filenameOfKnownHits.dat (not the default))\n"
                   + "   (3) map file for lookup (i.e. map=" + map + " (default))\n" +

                   "   OPTIONAL: \n"
                   + "   (4) maximum window around hit for LD search (i.e. maxWinSize="
                   + windowMaxSizePerSide + " (default))\n"
                   + "   (5) R2/D'/LOD threshold (.flt files are prefiltered. For .ld files filter will apply to R2 values. For haploview format files filter will apply to column in (6))(i.e. ldFilter="
                   + ldFilter + " (default))\n"
                   + "   (6) if input file is in haploview format, column on which to filter (i.e. column="
                   + column + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        return;
      } else if (arg.startsWith("file=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("ldFile=")) {
        ldFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("out=")) {
        outfile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("indexThresh=")) {
        indexThreshold = ext.parseFloatArg(arg);
        numArgs--;
      } else if (arg.startsWith("minWinSize=")) {
        windowMinSizePerSide = ext.parseIntArg(arg);
        numArgs--;
      } else if (arg.startsWith("maxWinSize=")) {
        windowMaxSizePerSide = ext.parseIntArg(arg);
        numArgs--;
      } else if (arg.startsWith("winThresh=")) {
        windowExtensionThreshold = ext.parseFloatArg(arg);
        numArgs--;
      } else if (arg.startsWith("knownHits=")) {
        knownHits = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("map=")) {
        map = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("log=")) {
        logfile = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("annotationCols=")) {
        annotationCols = ext.parseStringArg(arg, null).split(",");
        numArgs--;
      } else if (arg.startsWith("ldFilter=")) {
        ldFilter = ext.parseDoubleArg(arg);
        numArgs--;
      } else if (arg.startsWith("column=")) {
        column = ext.parseIntArg(arg);
        numArgs--;
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      return;
    }

    // generateHitsLookup("D:/ExomeChip/Hematology/00src/CHARGE-RBC/knownHits.dat", 200000,
    // "D:/ExomeChip/Hematology/00src/CHARGE-RBC/hitLookup.dat",
    // "D:/ExomeChip/Hematology/00src/CHARGE-RBC/ExomeChipV5_wMAF.csv");
    // System.exit(1);

    try {
      log = new Logger(logfile);
      if (ldFile != null) {
        determineLD(knownHits, map, ldFile, windowMaxSizePerSide, column, ldFilter, outfile, log);
      } else if (knownHits != null) {
        generateHitsLookup(knownHits, windowMinSizePerSide, outfile, map, log);
      } else {
        results = determine(filename, indexThreshold, windowMinSizePerSide,
                            windowExtensionThreshold, annotationCols, log);
        if (results != null) {
          if (outfile != null) {
            Files.writeMatrix(results, outfile, "\t");
          } else {
            for (String[] result : results) {
              log.report(ArrayUtils.toStr(result, "\t"));
            }
          }
        }
      }
      // determine("D:/mega/FromMike.11032011/pvals1.txt", "D:/mega/FromMike.11032011/pvals1.out",
      // indexThreshold, windowMinSizePerSide, windowExtensionThreshold);
      // determine("D:/mega/FromMike.11032011/pvals2.txt", "D:/mega/FromMike.11032011/pvals2.out",
      // indexThreshold, windowMinSizePerSide, windowExtensionThreshold);
      // determine("D:/mega/FromMike.11032011/fixed_together.txt",
      // "D:/mega/FromMike.11032011/fixed_together.out", indexThreshold, windowMinSizePerSide,
      // windowExtensionThreshold);
      // determine("D:/mega/FromMike.11032011/all_together.txt",
      // "D:/mega/FromMike.11032011/all_together.out", indexThreshold, windowMinSizePerSide,
      // windowExtensionThreshold);
      // determine("C:/CARe_data/conditionalMeta/uni_input.txt",
      // "C:/CARe_data/conditionalMeta/uniqueRegions.out", (float)0.00001, windowMinSizePerSide,
      // (float)0.0001);
      // determine("D:/tWork/Consortium/Megas/input.txt",
      // "D:/tWork/Consortium/Megas/uniqueRegions.out", indexThreshold, windowMinSizePerSide,
      // windowExtensionThreshold);
      // determine("D:/mega/filtered/OnlyThoseInRef/input.dat",
      // "D:/mega/filtered/OnlyThoseInRef/uniqueRegions.out", indexThreshold, windowMinSizePerSide,
      // windowExtensionThreshold);
      // determine("D:/mega/filtered/OnlyThoseInRef/originalWrongRE2_input.dat",
      // "D:/mega/filtered/OnlyThoseInRef/originalWrongRE2_uniqueRegions.out", indexThreshold,
      // windowMinSizePerSide, windowExtensionThreshold);
      // determine("D:/home/npankrat/jProjects/Ruhi/bioinf/input2.txt",
      // "D:/home/npankrat/jProjects/Ruhi/bioinf/my.out.xln", 0.00001f, 500000, 0.0001f);

      // determine("/home/input.dat", "/home/uniqueRegions.out", indexThreshold,
      // windowMinSizePerSide, windowExtensionThreshold);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
