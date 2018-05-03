package org.genvisis.cnv.filesys;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;
import java.util.Map;
import java.util.Vector;
import org.genvisis.CLI;
import org.genvisis.bioinformatics.Sequence;
import org.genvisis.cnv.annotation.markers.AnnotationFileLoader.QUERY_TYPE;
import org.genvisis.cnv.annotation.markers.AnnotationParser;
import org.genvisis.cnv.annotation.markers.MarkerAnnotationLoader;
import org.genvisis.cnv.annotation.markers.MarkerBlastAnnotation;
import org.genvisis.cnv.annotation.markers.MarkerSeqAnnotation;
import org.genvisis.cnv.manage.MarkerDataLoader;
import org.genvisis.cnv.prop.Property;
import org.genvisis.cnv.qc.IlluminaMarkerBlast;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.ext;
import org.genvisis.seq.manage.StrandOps;
import com.google.common.base.Joiner;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.io.Closeables;
import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.Allele;

public class ABLookup {

  public static final String DEFAULT_AB_FILE = "AB_lookup.dat";
  public static final String MARKER_LABEL = "Marker";
  public static final String A_LABEL = "A";
  public static final String B_LABEL = "B";
  public static final String[] AB_LOOKUP_COLS = new String[] {MARKER_LABEL, A_LABEL, B_LABEL};

  public static final String FLAGS_CLUSTER = "parseFromGenotypeClusterCenters";
  public static final String FLAGS_ORIGIN = "parseFromOriginalGenotypes";
  public static final String FLAGS_VCF = "parseFromAnnotationVCF";
  public static final String ARGS_MANIFEST = "IlluminaManifestFile";
  public static final String ARGS_PARTAB = "incompleteAB";
  public static final String ARGS_MAP = "mapFile";
  public static final String FLAGS_PLINK = "plinkFile";
  public static final String FLAGS_APPLYAB = "applyAB";

  public static final char MISSING_ALLELE = 'N';

  private String[] markerNames;
  private char[][] lookup;

  /**
   * Enum defining potential {@link ABLookup} parse sources.
   */
  public static enum ABSource {
    GENCLUSTER, MANIFEST, VCF, ORIGEN;
  }

  public ABLookup() {
    markerNames = new String[0];
    lookup = new char[0][0];
  }

  public ABLookup(String[] markerNames, char[][] lookup) {
    this.markerNames = markerNames;
    this.lookup = lookup;
  }

  public ABLookup(Project proj, String[] markerNames) {
    this(proj, markerNames, true, false);
  }

  public ABLookup(Project proj, String[] markerNames, boolean verbose, boolean allOrNothing) {
    this.markerNames = markerNames;

    Logger log = proj.getLog();

    if (Files.exists(proj.BLAST_ANNOTATION_FILENAME.getValue())) {
      lookup = generateABLookupFromAnnotationVCF(proj, markerNames, verbose, allOrNothing);
    } else if (Files.exists(proj.AB_LOOKUP_FILENAME.getValue())) {
      lookup = parseLookupFromFile(markerNames, proj.AB_LOOKUP_FILENAME.getValue(), verbose,
                                   allOrNothing, log);
    } else {
      log.reportError("Neither a Marker Blast Annotation nor an AB Lookup file for the project could be found, cannot generate AB Lookup");
      lookup = null;
    }

  }

  public ABLookup(String[] markerNames, String filename, boolean verbose, boolean allOrNothing,
                  Logger log) {
    this.markerNames = markerNames;
    lookup = parseLookupFromFile(markerNames, filename, verbose, allOrNothing, log);
  }

  private static char[][] parseLookupFromFile(String[] markerNames, String abLookupFilename,
                                              boolean verbose, boolean allOrNothing, Logger log) {
    Map<String, char[]> lookupMap;
    if (abLookupFilename.toLowerCase().endsWith(".csv")) {
      lookupMap = generateABLookupHashFromCSV(abLookupFilename, log);
    } else {
      lookupMap = generateABLookupHash(abLookupFilename, log);
    }
    return parseLookupFromMap(markerNames, lookupMap, allOrNothing, log, verbose);

  }

  private static char[][] parseLookupFromMap(String[] markerNames, Map<String, char[]> lookupMap,
                                             boolean allOrNothing, Logger log, boolean verbose) {
    char[][] lookup = new char[markerNames.length][];
    int missing = 0;
    for (int i = 0; i < markerNames.length; i++) {
      lookup[i] = lookupMap.get(markerNames[i]);
      if (lookup[i] == null) {
        if (verbose) {
          log.reportError("Error - no AB value for marker '" + markerNames[i] + "'");
        }
        missing++;
      }
    }
    if (missing > 0) {
      log.reportError("Warning - there " + (missing > 1 ? "were " : "was ") + missing + " marker"
                      + (missing > 1 ? "s" : "") + " without an AB value");
      if (allOrNothing) {
        lookup = null;
      }
    }
    return lookup;
  }

  public static char[] parseABFromMarkerSeqAnnotation(MarkerSeqAnnotation markerSeqAnnotation) {
    Strand strand = markerSeqAnnotation.getStrand();
    Allele A = markerSeqAnnotation.getA();
    Allele B = markerSeqAnnotation.getB();
    char a = MISSING_ALLELE;
    char b = MISSING_ALLELE;
    if (markerSeqAnnotation.isIndel()) {
      if (markerSeqAnnotation.isaDeletion()) {
        a = 'D';
        b = 'I';// actually this is typically just :equals reference;
      } else if (markerSeqAnnotation.isbDeletion()) {
        a = 'I';
        b = 'D';
      } else {
        throw new IllegalStateException("Both A and B should not be deletions");
      }
    } else {
      return parseABFromAlleles(A, B, strand);
    }
    return new char[] {a, b};
  }

  public static char[] parseABFromAlleles(Allele alleleA, Allele alleleB) {
    return parseABFromAlleles(alleleA, alleleB, Strand.POSITIVE);
  }

  public static char[] parseABFromAlleles(Allele alleleA, Allele alleleB, Strand strand) {
    char a = MISSING_ALLELE;
    char b = MISSING_ALLELE;
    if (alleleA.getBaseString().length() != alleleB.getBaseString().length()) {
      if (alleleA.getBaseString().length() < alleleB.getBaseString().length()) {
        a = 'D';
        b = 'I';// actually this is typically just :equals reference;
      } else {
        a = 'I';
        b = 'D';
      }
    } else {
      if (!alleleA.isSymbolic() && alleleA.getDisplayString().length() > 0) {
        a = StrandOps.flipIfNeeded(alleleA.getDisplayString(), strand, false).charAt(0);
      }
      if (!alleleB.isSymbolic() && alleleB.getDisplayString().length() > 0) {
        b = StrandOps.flipIfNeeded(alleleB.getDisplayString(), strand, false).charAt(0);
      }
    }
    return new char[] {a, b};
  }

  public void parseFromGenotypeClusterCenters(Project proj) {
    PrintWriter writer = null;
    PrintWriter writer2 = null;
    String trav;
    String[] samples;
    String[][] genotypesSeen;
    double[][] meanThetasForGenotypes;
    int[][] countsForGenotypes;
    Sample fsamp;
    byte[] genotypes;
    float[] thetas;
    char[] alleles;
    double travD;
    int order;
    MarkerSetInfo markerSet;
    Logger log;

    log = proj.getLog();
    samples = proj.getSamples();
    markerSet = proj.getMarkerSet();
    markerNames = markerSet.getMarkerNames();

    genotypesSeen = new String[markerNames.length][3];
    meanThetasForGenotypes = new double[markerNames.length][3];
    countsForGenotypes = new int[markerNames.length][3];
    log.report("Imputing AB alleles from the genotype cluster centers");
    for (int i = 0; i < samples.length; i++) {
      if (i % 100 == 0) {
        log.report((i + 1) + " of " + samples.length);
      }
      fsamp = proj.getFullSampleFromRandomAccessFile(samples[i]);
      thetas = fsamp.getThetas();
      genotypes = fsamp.getForwardGenotypes();

      for (int j = 0; j < markerNames.length; j++) {
        if (genotypes[j] > 0) {
          trav = Sample.ALLELE_PAIRS[genotypes[j]];
          if (trav.charAt(0) != trav.charAt(1)) {
            if (genotypesSeen[j][1] == null || genotypesSeen[j][1].equals(trav)) {
              genotypesSeen[j][1] = trav;
              meanThetasForGenotypes[j][1] += thetas[j];
              countsForGenotypes[j][1]++;
            } else {
              log.reportError("Error - different heterozygote (" + trav
                              + ") than the on seen previously (" + genotypesSeen[j][1]
                              + ") for marker " + markerNames[j] + " and sample " + samples[i]);
            }
          } else {
            if (genotypesSeen[j][0] == null || genotypesSeen[j][0].equals(trav)) {
              genotypesSeen[j][0] = trav;
              meanThetasForGenotypes[j][0] += thetas[j];
              countsForGenotypes[j][0]++;
            } else if (genotypesSeen[j][2] == null || genotypesSeen[j][2].equals(trav)) {
              genotypesSeen[j][2] = trav;
              meanThetasForGenotypes[j][2] += thetas[j];
              countsForGenotypes[j][2]++;
            } else {
              log.reportError("Error - different genotype (" + trav
                              + ") than anything seen before ("
                              + ArrayUtils.toStr(genotypesSeen[j], "/") + ") for marker "
                              + markerNames[j] + " and sample " + samples[i]);
            }
          }
        }
      }
    }

    try {
      writer = Files.getAppropriateWriter(proj.PROJECT_DIRECTORY.getValue() + "AB_breakdown.xln");
      writer.println("Marker\tG11\tG12\tG22\tG11 counts\tG12 counts\tG22 counts\tMean Theta G11\tMean Theta G12\tMean Theta G22\torder\tA allele\tB allele");
      writer2 = Files.getAppropriateWriter(proj.PROJECT_DIRECTORY.getValue() + "posssible_"
                                           + DEFAULT_AB_FILE);
      writer2.println("Marker\tA\tB");
      // int countMissing;
      lookup = new char[markerNames.length][];
      for (int i = 0; i < markerNames.length; i++) {
        // countMissing = 0;
        for (int j = 0; j < 3; j++) {
          meanThetasForGenotypes[i][j] /= countsForGenotypes[i][j];
          // if (countsForGenotypes[i][j] == 0) {
          // countMissing++;
          // }
        }
        writer.print(markerNames[i] + "\t" + genotypesSeen[i][0] + "\t" + genotypesSeen[i][1] + "\t"
                     + genotypesSeen[i][2] + "\t" + countsForGenotypes[i][0] + "\t"
                     + countsForGenotypes[i][1] + "\t" + countsForGenotypes[i][2] + "\t"
                     + meanThetasForGenotypes[i][0] + "\t" + meanThetasForGenotypes[i][1] + "\t"
                     + meanThetasForGenotypes[i][2]);

        alleles = new char[] {MISSING_ALLELE, MISSING_ALLELE};
        for (int j = 0; j < 3; j++) {
          for (int k = 0; k < 2; k++) {
            if (countsForGenotypes[i][j] > 0) {
              if (alleles[0] == MISSING_ALLELE) {
                alleles[0] = genotypesSeen[i][j].charAt(k);
              } else if (alleles[0] == genotypesSeen[i][j].charAt(k)) {} else if (alleles[1] == MISSING_ALLELE) {
                alleles[1] = genotypesSeen[i][j].charAt(k);
              } else if (alleles[1] == genotypesSeen[i][j].charAt(k)) {} else {
                log.reportError("Error - too many alleles for " + markerNames[i] + " first "
                                + alleles[0] + " and " + alleles[1] + " and now "
                                + genotypesSeen[i][j].charAt(k));
              }
            }
          }
        }
        order = 0;
        travD = -1;
        for (int j = 0; j < 3; j++) {
          if (countsForGenotypes[i][j] > 0) {
            if (travD == -1) {
              travD = meanThetasForGenotypes[i][j];
            } else if (meanThetasForGenotypes[i][j] > travD) {
              if (order == 0) {
                order = 1;
              } else if (order != 1) {
                order = -1;
              }
            } else if (meanThetasForGenotypes[i][j] < travD) {
              if (order == 0) {
                order = 2;
              } else if (order != 2) {
                order = -1;
              }
            } else {
              log.reportError("Error - equal meanThetas?");
            }
          }
        }
        if (order == 0 && meanThetasForGenotypes[i][0] > 0.50) {
          order = 2;
        }
        if (order == -1 && countsForGenotypes[i][0] > 0 && countsForGenotypes[i][2] > 0) {
          order = meanThetasForGenotypes[i][0] < meanThetasForGenotypes[i][2] ? -1 : -2;
        }
        if (order == 2 || order == -2) {
          lookup[i] = new char[] {alleles[1], alleles[0]};
        } else {
          lookup[i] = alleles;
        }
        writer.println("\t" + order + "\t" + lookup[i][0] + "\t" + lookup[i][1]);
        writer2.println(markerNames[i] + "\t" + lookup[i][0] + "\t" + lookup[i][1]);
      }
      writer.close();
      writer2.close();
    } catch (Exception e) {
      log.reportError("Error writing to either " + "AB_breakdown.xln" + " or " + "possible_"
                      + DEFAULT_AB_FILE);
      log.reportException(e);
    } finally {
      if (writer != null) {
        writer.close();
      }
      if (writer2 != null) {
        writer2.close();
      }
    }

  }

  /**
   * {@link ABLookup#parseFromManifest(Project, String)}
   */
  public void parseFromManifest(Project proj, String manifestFile) {
    if (Files.exists(proj.BLAST_ANNOTATION_FILENAME.getValue())) {
      proj.getLog()
          .reportError("Error - annotation file '" + proj.BLAST_ANNOTATION_FILENAME.getValue()
                       + "' already exists; use that instead of the Illumina Manifest file to parse AB values (or rename that file to get around this check)");
      return;
    }
    new IlluminaMarkerBlast(proj, -1, -1, -1, false, false, false, 1, manifestFile, true).blastEm();
    parseFromAnnotationVCF(proj);
  }

  /**
   * @param proj Will try to generate the ab lookup from the annotation vcf. If fails, will use
   *          {@link ABLookup#parseFromOriginalGenotypes(Project)}
   */
  public void parseFromAnnotationVCF(Project proj) {
    if (!Files.exists(proj.BLAST_ANNOTATION_FILENAME.getValue())) {
      proj.getLog().reportTimeWarning("Could not find " + proj.BLAST_ANNOTATION_FILENAME.getValue()
                                      + ", trying to parse from original genotypes instead");
    } else {
      proj.getLog()
          .reportTimeWarning("This method has not been completely tested, you have been warned");
      proj.getLog().reportTimeWarning("This method will convert AB genotypes to positive strand");

      try {
        MarkerDetailSet markerSet = proj.getMarkerSet();
        markerNames = markerSet.getMarkerNames();

        Map<String, MarkerBlastAnnotation> masterMarkerList = MarkerBlastAnnotation.initForMarkers(markerNames);
        MarkerAnnotationLoader annotationLoader = new MarkerAnnotationLoader(null,
                                                                             proj.BLAST_ANNOTATION_FILENAME.getValue(),
                                                                             markerSet, true,
                                                                             proj.getLog());
        annotationLoader.setReportEvery(10000);
        List<Map<String, ? extends AnnotationParser>> parsers = Lists.newArrayList();
        parsers.add(masterMarkerList);
        annotationLoader.fillAnnotations(null, parsers, QUERY_TYPE.ONE_TO_ONE);
        Map<String, Integer> indices = proj.getMarkerIndices();// should be in perfect order
                                                               // but just in case;
        lookup = new char[markerNames.length][2];

        for (String markerName : markerNames) {
          int indx = indices.get(markerName);
          lookup[indx] = parseABFromMarkerSeqAnnotation(masterMarkerList.get(markerName)
                                                                        .getMarkerSeqAnnotation());
        }
        return;

      } catch (Exception e) {
        proj.getLog()
            .reportTimeWarning("Could not generate AB lookup from "
                               + proj.BLAST_ANNOTATION_FILENAME.getValue()
                               + ", trying to parse from original genotypes instead");

      }
    }
    parseFromOriginalGenotypes(proj);
  }

  public void parseFromOriginalGenotypes(Project proj) {
    String fullGenotype;
    String[] samples;
    String[][] genotypesSeen;
    int[][] countsForGenotypes;
    Sample fsamp;
    byte[] abGenotypes, fullGenotypes;
    MarkerSetInfo markerSet;
    Logger log;

    log = proj.getLog();
    samples = proj.getSamples();
    markerSet = proj.getMarkerSet();
    markerNames = markerSet.getMarkerNames();

    genotypesSeen = new String[markerNames.length][3];
    countsForGenotypes = new int[markerNames.length][3];
    log.report("Parsing AB alleles from the original genotypes");
    for (int i = 0; i < samples.length; i++) {
      if (i % 100 == 0) {
        log.report((i + 1) + " of " + samples.length);
      }
      fsamp = proj.getFullSampleFromRandomAccessFile(samples[i]);
      abGenotypes = fsamp.getAB_Genotypes();
      fullGenotypes = fsamp.getForwardGenotypes();

      for (int j = 0; j < markerNames.length; j++) {
        if (abGenotypes[j] >= 0) {
          fullGenotype = Sample.ALLELE_PAIRS[fullGenotypes[j]];
          if (genotypesSeen[j][abGenotypes[j]] == null) {
            genotypesSeen[j][abGenotypes[j]] = fullGenotype;
          } else if (!fullGenotype.equals(genotypesSeen[j][abGenotypes[j]])) {
            log.reportError("Error - different heterozygote (" + fullGenotype
                            + ") than the on seen previously (" + genotypesSeen[j][1]
                            + ") for marker " + markerNames[j] + " and sample " + samples[i]);
          }
          countsForGenotypes[j][abGenotypes[j]]++;
        }
      }
    }

    lookup = new char[markerNames.length][2];
    for (int i = 0; i < markerNames.length; i++) {
      lookup[i] = new char[] {MISSING_ALLELE, MISSING_ALLELE};
      if (countsForGenotypes[i][1] > 0
          && genotypesSeen[i][1].charAt(0) == genotypesSeen[i][1].charAt(1)) {
        log.reportError("Error - impossible heterozygote (" + genotypesSeen[i][1] + ") for marker "
                        + markerNames[i]);
      }
      if (countsForGenotypes[i][0] > 0) {
        if (genotypesSeen[i][0].charAt(0) != genotypesSeen[i][0].charAt(1)) {
          log.reportError("Error - impossible A/A homozygote (" + genotypesSeen[i][0]
                          + ") for marker " + markerNames[i]);
        } else {
          lookup[i][0] = genotypesSeen[i][0].charAt(0);
        }
      }
      if (countsForGenotypes[i][2] > 0) {
        if (genotypesSeen[i][2].charAt(0) != genotypesSeen[i][2].charAt(1)) {
          log.reportError("Error - impossible B/B homozygote (" + genotypesSeen[i][2]
                          + ") for marker " + markerNames[i]);
        } else {
          lookup[i][1] = genotypesSeen[i][2].charAt(0);
        }
      }
      if (countsForGenotypes[i][1] > 0
          && countsForGenotypes[i][0] + countsForGenotypes[i][2] == 0) {
        log.reportError("Error - only heterozygotes present for marker " + markerNames[i]
                        + "; cannot determine which allele is the A allele versus the B allele");
      } else if (countsForGenotypes[i][1] > 0 && countsForGenotypes[i][0] == 0) {
        if (genotypesSeen[i][1].charAt(0) == lookup[i][1]) {
          lookup[i][0] = genotypesSeen[i][1].charAt(1);
        } else if (genotypesSeen[i][1].charAt(1) == lookup[i][1]) {
          lookup[i][0] = genotypesSeen[i][1].charAt(0);
        } else {
          log.reportError("Error - heterozygote for marker " + markerNames[i] + " ("
                          + genotypesSeen[i][1] + ") does not match up with the B allele ("
                          + lookup[i][1] + ")");
        }
      } else if (countsForGenotypes[i][1] > 0 && countsForGenotypes[i][2] == 0) {
        if (genotypesSeen[i][1].charAt(0) == lookup[i][0]) {
          lookup[i][1] = genotypesSeen[i][1].charAt(1);
        } else if (genotypesSeen[i][1].charAt(1) == lookup[i][0]) {
          lookup[i][1] = genotypesSeen[i][1].charAt(0);
        } else {
          log.reportError("Error - heterozygote for marker " + markerNames[i] + " ("
                          + genotypesSeen[i][1] + ") does not match up with the A allele ("
                          + lookup[i][1] + ")");
        }
      }
    }
  }

  /**
   * Static entry point for parsing {@link ABLookup}s from all sources.
   *
   * @param proj Base project for the {@link ABLookup}
   * @param parseSource Where to parse the {@link ABLookup}
   * @param outfile Where to write the parsed {@link ABLookup}
   * @param args Optional parameters if required by parse mechanism (e.g. a manifestFile for
   *          {@link #parseFromManifest(Project, String)}
   */
  public static void parseABLookup(Project proj, ABSource parseSource, String outfile,
                                   Object... args) {
    ABLookup abLookup = new ABLookup();

    if (!Files.exists(outfile, true)) {
      switch (parseSource) {
        case GENCLUSTER:
          abLookup.parseFromGenotypeClusterCenters(proj);
          break;
        case MANIFEST:
          abLookup.parseFromManifest(proj, (String) args[0]);
          break;
        case VCF:
          abLookup.parseFromAnnotationVCF(proj);
          break;
        case ORIGEN:
          abLookup.parseFromOriginalGenotypes(proj);
          break;
      }
      abLookup.writeToFile(outfile, proj.getLog());
    }

  }

  public static Map<String, char[]> generateABLookupHash(String filename, Logger log) {
    BufferedReader reader = null;
    String[] line;
    Map<String, char[]> map;

    map = Maps.newHashMap();
    try {
      reader = new BufferedReader(new FileReader(filename));
      while (reader.ready()) {
        line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
        map.put(line[0], new char[] {line[1].charAt(0), line[2].charAt(0)});
      }
      return map;
    } catch (FileNotFoundException fnfe) {
      return null;
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + filename + "\"");
      return null;
    } finally {
      Closeables.closeQuietly(reader);
    }
  }

  public static Map<String, char[]> generateABLookupHashFromCSV(String filename, Logger log) {
    BufferedReader reader;
    String[] line;
    Map<String, char[]> map;
    int[] indices;
    String temp, prev;
    int count = 0;

    map = Maps.newHashMap();
    try {
      reader = Files.getAppropriateReader(filename);
      do {
        temp = reader.readLine();
        line = temp.trim().split(",");
      } while (reader.ready() && (!temp.contains("Name") || !temp.contains("SNP")));
      indices = ext.indexFactors(new String[][] {{"Name", "MarkerName"}, {"SNP"}}, line, false,
                                 true, true, log);
      prev = temp;
      while (reader.ready()) {
        temp = reader.readLine();
        if (temp.startsWith("[")) {
          count = 0;
          while (reader.ready()) {
            reader.readLine();
            count++;
          }
          log.report("Ended with an ignored tail of " + count + " line(s)");
        } else {
          try {
            line = temp.trim().split(",");
            map.put(line[indices[0]],
                    new char[] {line[indices[1]].charAt(1), line[indices[1]].charAt(3)});
          } catch (ArrayIndexOutOfBoundsException aioobe) {
            log.reportError("Error - could not parse line:");
            log.reportError(temp);
            log.reportError("Where the previous line was:");
            log.reportError(prev);
            log.reportException(aioobe);
          }
          prev = temp;
        }
      }
      reader.close();
      return map;
    } catch (FileNotFoundException fnfe) {
      log.reportError("File not found: \"" + filename + "\"");
      return null;
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + filename + "\"");
      return null;
    }
  }

  public static char[][] generateABLookupFromAnnotationVCF(Project proj, String[] markerNames,
                                                           boolean verbose, boolean allOrNothing) {
    char[][] lookup = new char[markerNames.length][];
    if (!Files.exists(proj.BLAST_ANNOTATION_FILENAME.getValue())) {
      proj.getLog().reportTimeWarning("Could not find " + proj.BLAST_ANNOTATION_FILENAME.getValue()
                                      + ", cannot generate AB Lookup");
      return lookup;
    }

    MarkerDetailSet markerSet = proj.getMarkerSet();
    Map<String, MarkerBlastAnnotation> masterMarkerList = MarkerBlastAnnotation.initForMarkers(markerNames);
    MarkerAnnotationLoader annotationLoader = new MarkerAnnotationLoader(null,
                                                                         proj.BLAST_ANNOTATION_FILENAME.getValue(),
                                                                         markerSet, true,
                                                                         proj.getLog());
    annotationLoader.setReportEvery(10000);
    List<Map<String, ? extends AnnotationParser>> parsers = Lists.newArrayList();
    parsers.add(masterMarkerList);
    annotationLoader.fillAnnotations(null, parsers, QUERY_TYPE.ONE_TO_ONE);

    int missing = 0;
    for (int i = 0; i < markerNames.length; i++) {
      try {
        lookup[i] = parseABFromMarkerSeqAnnotation(masterMarkerList.get(markerNames[i])
                                                                   .getMarkerSeqAnnotation());
      } catch (NullPointerException npe) {
        if (verbose) {
          proj.getLog().reportError("no AB value for marker '" + markerNames[i] + "'");
        }
        missing++;
      }
    }
    if (missing > 0) {
      proj.getLog().reportError("Warning - there " + (missing > 1 ? "were " : "was ") + missing
                                + " marker" + (missing > 1 ? "s" : "") + " without an AB value");
      if (allOrNothing) {
        lookup = null;
      }
    }
    return lookup;

  }

  public char[][] getLookup() {
    return lookup;
  }

  public void writeToFile(String outfile, Logger log) {
    PrintWriter writer = null;
    new File(ext.parseDirectoryOfFile(outfile)).mkdirs();
    try {
      writer = Files.getAppropriateWriter(outfile);
      writer.println(Joiner.on("\t").join(AB_LOOKUP_COLS));
      for (int i = 0; i < markerNames.length; i++) {
        writer.println(markerNames[i] + "\t" + lookup[i][0] + "\t" + lookup[i][1]);
      }
      writer.flush();
      writer.close();
    } catch (Exception e) {
      log.reportError("Error writing to " + outfile);
      log.reportException(e);
    } finally {
      if (writer != null) {
        writer.close();
      }
    }
  }

  public static void applyABLookupToFullSampleFiles(Project proj) {
    applyABLookupToFullSampleFiles(proj, null);
  }

  public static void applyABLookupToFullSampleFiles(Project proj, String abLookupFilename) {
    ABLookup abLookup;
    Sample fsamp;
    String[] samples;
    String genotype;
    byte[] forwardGenotypes, abGenotypes;
    Logger log;
    // Save the given filename if given
    // Otherwise, use the stored value
    if (abLookupFilename != null) {
      proj.AB_LOOKUP_FILENAME.setValue(abLookupFilename);
      proj.saveProperties(new Property[] {proj.AB_LOOKUP_FILENAME});
    } else {
      if (!Files.exists(proj.AB_LOOKUP_FILENAME.getValue())) {
        proj.getLog()
            .reportError("Error - cannot applyABLookupToFullSampleFiles without the AB Lookup file ('"
                         + proj.AB_LOOKUP_FILENAME.getValue() + "').");
        return;
      }
      abLookupFilename = proj.AB_LOOKUP_FILENAME.getValue();
    }
    log = proj.getLog();
    abLookup = new ABLookup(proj.getMarkerNames(), abLookupFilename, false, false, proj.getLog());
    samples = proj.getSamples();
    int alleleMismatches = 0;
    int totalGenotypes = 0;
    for (int i = 0; i < samples.length; i++) {
      if (i % 100 == 0) {
        log.report((i + 1) + " of " + samples.length);
      }
      fsamp = proj.getFullSampleFromRandomAccessFile(samples[i]);
      forwardGenotypes = fsamp.getForwardGenotypes();
      if (forwardGenotypes.length != abLookup.markerNames.length) {
        log.reportError("Error - mismatched array lengths for forwardGenotypes and abLookup markerNames");
      }
      totalGenotypes += abLookup.markerNames.length;
      abGenotypes = new byte[forwardGenotypes.length];
      for (int j = 0; j < abLookup.markerNames.length; j++) {
        if (forwardGenotypes[j] == 0) {
          abGenotypes[j] = -1;
        } else {
          // genotype = FullSample.ALLELE_PAIRS[forwardGenotypes[j]];
          genotype = Sample.ALLELE_PAIRS[forwardGenotypes[j]];
          for (int k = 0; k < 2; k++) {
            if (genotype.charAt(k) == abLookup.lookup[j][1]) {
              abGenotypes[j]++;
            } else if (genotype.charAt(k) != abLookup.lookup[j][0]) {
              alleleMismatches++;
            }
          }
        }
      }
      fsamp.setAB_Genotypes(abGenotypes);
      fsamp.saveToRandomAccessFile(proj.SAMPLE_DIRECTORY.getValue(false, true) + samples[i]
                                   + Sample.SAMPLE_FILE_EXTENSION);
    }
    if (alleleMismatches > 0) {
      log.reportError("Warning - " + alleleMismatches + " total mismatched forward alleles found; ~"
                      + alleleMismatches / samples.length + " per sample; ~"
                      + ext.formPercent(alleleMismatches / (totalGenotypes * 2.0), 4)
                      + " of total alleles");
    }
  }

  /**
   * @return {@code false} if an exception was caught, aborting execution. {@code true} otherwise.
   */
  public static boolean fillInMissingAlleles(Project proj, String incompleteABlookupFilename,
                                             String mapFile, boolean updatingPlinkFile) {
    BufferedReader reader = null;
    PrintWriter writer = null;
    String[] line;
    Map<String, char[]> lookupMap;
    char knownAllele;
    int knownIndex;
    Vector<String> markersWithNoLink;
    MarkerDataLoader markerDataLoader;
    String[] markerNames;
    char[] refAlleles;
    Logger log;
    ClusterFilterCollection clusterFilterCollection;
    String output;
    int markerIndex, first, second;

    log = proj.getLog();

    if (mapFile == null) {
      log.reportError("Error - no Illumina map file provided to ABLookup; aborting...");
      return false;
    }

    if (!Files.exists(mapFile)) {
      log.reportError("Error - could not find Illumina map file '" + mapFile + "'; aborting...");
      return false;
    }

    if (!mapFile.toLowerCase().endsWith(".csv")) {
      log.reportError("Error - expecting an Illumina style format, but this map file '" + mapFile
                      + "' does not end in .csv; aborting...");
      return false;
    }

    if (!Files.exists(incompleteABlookupFilename)) {
      log.reportError("Error - could not find incomplete ABLookup file '"
                      + incompleteABlookupFilename + "'; aborting...");
      return false;
    }

    lookupMap = generateABLookupHashFromCSV(mapFile, proj.getLog());

    markersWithNoLink = new Vector<>();
    try {
      reader = Files.getAppropriateReader(incompleteABlookupFilename);
      writer = Files.getAppropriateWriter(ext.addToRoot(incompleteABlookupFilename, "_filledIn"));
      if (updatingPlinkFile) {
        markerIndex = 1;
        first = 4;
        second = 5;
      } else {

        line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
        int[] colIndices = ext.indexFactors(AB_LOOKUP_COLS, line, false, log, true);

        markerIndex = colIndices[0];
        first = colIndices[1];
        second = colIndices[2];
        if (markerIndex == -1 || first == -1 || second == -1) {
          log.reportError(incompleteABlookupFilename + " does not contain required columns: "
                          + Joiner.on("\t").join(AB_LOOKUP_COLS));
          reader.close();
          writer.close();
          return false;
        }

        writer.println(ArrayUtils.toStr(line));
      }
      while (reader.ready()) {
        line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
        if (line[first].equals("N") && line[second].equals("N")) {
          markersWithNoLink.add(line[0]);
          line[first] = "A";
          line[second] = "B";
        } else if (line[first].equals("N") || line[second].equals("N")) {
          knownIndex = line[first].equals("N") ? 1 : 0;
          knownAllele = line[first + knownIndex].charAt(0);
          refAlleles = lookupMap.get(line[markerIndex]);
          if (refAlleles == null) {
            log.reportError("Error - allele lookup failed for marker " + line[markerIndex]);
          } else {
            for (int i = 0; knownIndex != -9 && i < refAlleles.length; i++) {
              if (knownAllele == refAlleles[i]) {
                line[first + (1 - knownIndex)] = refAlleles[1 - i] + "";
                knownIndex = -9;
              }
            }
            if (knownIndex != -9) {
              for (int i = 0; knownIndex != -9 && i < refAlleles.length; i++) {
                if (knownAllele == Sequence.flip(refAlleles[i])) {
                  line[first + (1 - knownIndex)] = Sequence.flip(refAlleles[1 - i]) + "";
                  knownIndex = -9;
                }
              }
            }
            if (knownIndex != -9) {
              log.reportError("Error - failed to reconcile " + ArrayUtils.toStr(line, "/")
                              + " with reference alleles: " + refAlleles[0] + "/" + refAlleles[1]);
            }
          }

        }
        writer.println(ArrayUtils.toStr(line));
      }

      clusterFilterCollection = proj.getClusterFilterCollection();

      if (!markersWithNoLink.isEmpty()) {
        String markersWithNoLinkFile = Files.getNextAvailableFilename(ext.rootOf(incompleteABlookupFilename,
                                                                                 false)
                                                                      + "_test_markersWithNoLink#.txt");
        log.report("There were " + markersWithNoLink.size()
                   + " markers that were zeroed out in the original export; these are set to alleles 'A' and 'B'; for list check "
                   + markersWithNoLinkFile);
        markerNames = ArrayUtils.toStringArray(markersWithNoLink);
        markerDataLoader = null;
        output = Files.getNextAvailableFilename(markersWithNoLinkFile);
        try {
          if (Files.exists(proj.MARKER_DATA_DIRECTORY.getValue(false, false))) {
            markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj,
                                                                                       markerNames);
            log.reportError("Warning - allele frequencies for any chrX markers will be slightly inaccurate");
          } else {
            log.report("Warning - since " + proj.MARKER_DATA_DIRECTORY.getValue(false, false)
                       + " does not exist, marker data can not be loaded and frequency of B allele will not be reported in "
                       + output
                       + ".\n If you would like to obtain the frequency of B allele for these markers, please transpose the data and then run the following");
            log.report("java -jar /your/path/to/" + org.genvisis.common.PSF.Java.GENVISIS
                       + " cnv.filesys.ABLookup proj=" + proj.getPropertyFilename()
                       + " incompleteAB=" + incompleteABlookupFilename + " mapFile=" + mapFile);
          }
        } catch (NullPointerException nullPointerException) {// MarkerDataLoader will likely throw
          // this if there are other issues
          log.report("Warning - was not able to load marker data, frequency of B allele will not be reported in "
                     + output);
          log.reportException(nullPointerException);
        }
        for (int i = 0; i < markerNames.length; i++) {
          if (markerDataLoader == null) {
            markerNames[i] = markerNames[i];// skip frequency of b allele
          } else {
            MarkerData markerData = markerDataLoader.requestMarkerData(i);
            // markerNames[i] = markerNames[i] + "\t" + markerData.getFrequencyOfB(null, null,
            // clusterFilterCollection, proj.getFloat(proj.GC_THRESHOLD));
            markerNames[i] = markerNames[i] + "\t"
                             + markerData.getFrequencyOfB(null, null, clusterFilterCollection,
                                                          proj.GC_THRESHOLD.getValue().floatValue(),
                                                          log);
            markerDataLoader.releaseIndex(i);
          }
        }
        Files.writeArray(markerNames, output);
      }
    } catch (FileNotFoundException fnfe) {
      log.reportError("Error: file \"" + incompleteABlookupFilename
                      + "\" not found in current directory");
      return false;
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + incompleteABlookupFilename + "\"");
      return false;
    } finally {
      if (writer != null) {
        writer.close();
      }
      Closeables.closeQuietly(reader);
    }
    return true;
  }

  public static void main(String... args) {
    Project proj;
    String outfile = DEFAULT_AB_FILE;
    String mapFile = "SNP_Map.csv";
    String projFile = org.genvisis.cnv.Launch.getDefaultDebugProjectFile(false);
    String illumina = "infiniumomni2-5-8-v1-3-a1-manifest-file-csv.zip";
    String abLookup = "possible_AB_lookup.dat";

    CLI c = new CLI(ABLookup.class);
    c.addArgWithDefault(CLI.ARG_PROJ, CLI.DESC_PROJ, projFile);
    c.addArgWithDefault(CLI.ARG_OUTFILE, CLI.DESC_OUTFILE, outfile);
    c.addFlag(FLAGS_CLUSTER, "parse ABLookup from centroids");
    c.addFlag(FLAGS_ORIGIN, "parse ABLookup from existing original genotypes");
    c.addArg(ARGS_MANIFEST, "parse ABLookup from Illumina manifest file", illumina);
    c.addFlag(FLAGS_VCF, "parse ABLookup from VCF annotation file");
    c.addArg(ARGS_PARTAB, "fill in a partial existing ABLookup file using an Illumina SNP Table",
             abLookup);
    c.addArg(ARGS_MAP, "the filename of the Illumina SNP Tabl", mapFile);
    c.addFlag(FLAGS_PLINK, "use a plink.bim file as input instead of an ABLookup file");
    c.addFlag(FLAGS_APPLYAB, "apply the project's AB lookup to all Sample files in project");

    c.addGroup(FLAGS_APPLYAB, FLAGS_VCF, ARGS_PARTAB, ARGS_MANIFEST, FLAGS_ORIGIN, FLAGS_CLUSTER);

    c.parseWithExit(args);

    proj = new Project(c.get(CLI.ARG_PROJ));
    if (!new File(outfile).isAbsolute()) outfile = proj.PROJECT_DIRECTORY.getValue()
                                                   + c.get(CLI.ARG_OUTFILE);
    if (c.has(FLAGS_APPLYAB)) {
      applyABLookupToFullSampleFiles(proj);
    } else if (c.has(FLAGS_VCF)) {
      parseABLookup(proj, ABSource.VCF, outfile);
    } else if (c.has(ARGS_PARTAB)) {
      fillInMissingAlleles(proj, c.get(ARGS_PARTAB), c.get(ARGS_MAP), c.has(FLAGS_PLINK));
    } else if (c.has(ARGS_MANIFEST)) {
      parseABLookup(proj, ABSource.MANIFEST, outfile, c.get(ARGS_MANIFEST));
    } else if (c.has(FLAGS_ORIGIN)) {
      parseABLookup(proj, ABSource.ORIGEN, outfile);
    } else if (c.has(FLAGS_CLUSTER)) {
      parseABLookup(proj, ABSource.GENCLUSTER, outfile);
    }
  }
}
