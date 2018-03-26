package org.genvisis.cnv.hmm;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;
import java.util.stream.IntStream;
import org.genvisis.CLI;
import org.genvisis.cnv.filesys.MarkerSetInfo;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.manage.ExtProjectDataParser;
import org.genvisis.cnv.manage.ExtProjectDataParser.ProjectDataParserBuilder;
import org.genvisis.cnv.qc.GcAdjustor.GcModel;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.ext;

/**
 * @author lane0212 Handles the pfb data for the hmm in reduced format
 */
public class PFB {

  private final Project proj;
  private final double[] pfbs;

  public PFB(PFB pfb) {
    proj = pfb.proj;
    pfbs = pfb.pfbs;
  }

  private PFB(Project proj, double[] pfbst) {
    super();
    this.proj = proj;
    pfbs = pfbst;
    MarkerSetInfo markerSet = proj.getMarkerSet();
    if (pfbs.length != markerSet.getMarkerNames().length) {
      String error = "Found " + pfbs.length + " pfb entries, but the project has "
                     + markerSet.getMarkerNames().length + " markers";
      proj.getLog().reportError(error);
      throw new IllegalArgumentException(error);
    } else {
      this.proj.getLog().reportTimeInfo("Loaded " + pfbst.length + " pfb entries");
    }
    int problems = 0;
    for (int i = 0; i < pfbs.length; i++) { // what PennCNV does
      if (!Double.isNaN(pfbs[i]) && pfbs[i] >= 0 && pfbs[i] <= 1) {
        if (pfbs[i] < 0.01) {
          pfbs[i] = 0.01;
        }
        if (pfbs[i] > .99) {
          pfbs[i] = .99;
        }
      } else if (!Double.isNaN(pfbs[i]) && pfbs[i] < 1) {
        if (!proj.ARRAY_TYPE.getValue().isCNOnly(markerSet.getMarkerNames()[i])) {
          problems++;
        }
      }
    }
    if (problems > 0) {
      proj.getLog().reportTimeWarning(problems + " markers " + " had a pfb value less than 1 ");
    }
  }

  public double[] getPfbs() {
    return pfbs;
  }

  /**
   * Calculate the population BAF (B Allele Frequency based on all the samples available in the)
   * data. Output is going to be saved on disk. In PennCnv, this file is also called snpFile.
   *
   * @param proj The project you are going to run PennCNV on. The output file looks like the the
   *          following: Name Chr Position PFB rs1000113 5 150220269 0.564615751221256 rs1000115 9
   *          112834321 0.565931333264192 rs10001190 4 6335534 0.5668604380025 rs10002186 4 38517993
   *          0.57141752993563 rs10002743 4 6327482 0.567557695424774
   */
  public static String populationBAF(Project proj) {
    String[] sampleList;
    String output;
  
    Logger log = proj.getLog();
    String filename = proj.SAMPLE_SUBSET_FILENAME.getValue(true, false);
  
    if (ext.rootOf(filename) == null || ext.rootOf(filename).equals("")
        || !Files.exists(filename)) {
      sampleList = proj.getSampleList().getSamples();
      output = proj.CUSTOM_PFB_FILENAME.getValue(true, false);
    } else if (Files.exists(filename)) {
      log.report("filename: " + filename);
      sampleList = HashVec.loadFileToStringArray(filename, false, new int[] {0}, false);
      output = proj.PROJECT_DIRECTORY.getValue() + ext.rootOf(filename) + ".pfb";
    } else {
      proj.message("Failed to load \"" + filename + "\"");
      return null;
    }
  
    MarkerSetInfo markerSet = proj.getMarkerSet();
    String[] markerNames = markerSet.getMarkerNames();
    byte[] chrs = markerSet.getChrs();
    int[] positions = markerSet.getPositions();
    double[] bafSum = new double[chrs.length];
    int[] bafCounts = new int[chrs.length];
    int[] genoCounts = new int[chrs.length];
    int logPer = sampleList.length / 100; // 100 log messages
    for (int i = 0; i < sampleList.length; i++) {
      if (i % logPer == 0) {
        log.reportTime("Loading file " + (i + 1) + " of " + sampleList.length);
      }
      Sample samp = proj.getPartialSampleFromRandomAccessFile(sampleList[i], false, false, true,
                                                              false, true);
      float[] bafs = samp.getBAFs();
      byte[] genotypes = samp.getAB_Genotypes();
      IntStream.range(0, bafSum.length)
               .parallel()
               .forEach((j) -> {
                 if (!Float.isNaN(bafs[j])) {
                   bafSum[j] += bafs[j];
                   bafCounts[j]++;
                   if (genotypes[j] >= 0) {
                     genoCounts[j]++;
                   }
                 }
               });
    }
    double[] bafAverage = new double[chrs.length];
    Set<String> missingGenotypeMarkers = Collections.synchronizedSet(new HashSet<>());
    
    IntStream.range(0, bafSum.length)
             .parallel()
             .forEach((i) -> {
               boolean cnOnly = proj.getArrayType().isCNOnly(markerNames[i]);
               if (genoCounts[i] != 0 && !cnOnly) {// Since mock genotypes can be present, we demand non-CN
                 // only
                 bafAverage[i] = bafSum[i] / bafCounts[i];
               } else if (cnOnly) {
                 bafAverage[i] = 2;
               } else {
                 bafAverage[i] = -1; // This is to more clearly differentiate CN only markers from SNPs
                 // without callrate
                 
                 missingGenotypeMarkers.add(markerNames[i]);
               }
             });
  
    PSF.checkInterrupted();
    try {
  
      PrintWriter writer = Files.openAppropriateWriter(output);
      writer.println("Name\tChr\tPosition\tPFB");
      for (int i = 0; i < markerNames.length; i++) {
        writer.println(markerNames[i] + "\t"
                       + (chrs[i] < 23 ? chrs[i]
                                       : (chrs[i] == 23 ? "X"
                                                        : (chrs[i] == 24 ? "Y"
                                                                         : (chrs[i] == 25 ? "XY"
                                                                                          : (chrs[i] == 26 ? "M"
                                                                                                           : "Un")))))
                       + "\t" + positions[i] + "\t" + bafAverage[i]);
      }
      writer.close();
      log.report("Population BAF file is now ready at: " + output);
      if (!missingGenotypeMarkers.isEmpty()) {
        String missingGenoFile = ext.addToRoot(output, ".missingGenotypes");
        log.reportTimeInfo(missingGenotypeMarkers.size()
                           + " markers had missing genotypes and were set to -1 in " + output
                           + ". These markers can be treated as CN only markers, or removed at your discretion with CNVCaller");
        Files.writeIterable(missingGenotypeMarkers, missingGenoFile);
      }
    } catch (Exception e) {
      log.reportError("Error writing to '" + output + "'");
      log.reportException(e);
    }
    return output;
  }

  public static PFB loadPFB(Project proj) {
    return loadPFB(proj, proj.CUSTOM_PFB_FILENAME.getValue());
  }

  public static PFB loadPFB(Project proj, String fullPathToPfb) {

    ProjectDataParserBuilder builder = new ProjectDataParserBuilder();
    builder.dataKeyColumnName("Name");
    builder.numericDataTitles(new String[] {"PFB"});
    builder.sampleBased(false);
    builder.requireAll(true);
    builder.treatAllNumeric(false);
    builder.verbose(false);

    try {
      ExtProjectDataParser extProjectDataParser = builder.build(proj, fullPathToPfb);
      extProjectDataParser.determineIndicesFromTitles();
      extProjectDataParser.loadData();
      double[] pfbs = extProjectDataParser.getNumericDataForTitle("PFB");
      return new PFB(proj, pfbs);
    } catch (FileNotFoundException e) {
      e.printStackTrace();
      proj.getLog().reportFileNotFound(fullPathToPfb);
      return null;
    }

  }
  
  public static void main(String[] args) {
    CLI cli = new CLI(PFB.class);
    
    cli.addArg("filename", "Project properties file");
    cli.addArg("logfile", "Project log file", false);
    
    cli.parseWithExit(args);
    
    String filename = cli.get("filename");
    String logfile = cli.has("logfile") ? cli.get("logfile") : null;
    
    Project proj = new Project(filename, logfile);

    populationBAF(proj);
  }

}
