package org.genvisis.cnv.annotator;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import org.genvisis.cnv.analysis.BeastScore;
import org.genvisis.cnv.filesys.CNVariant;
import org.genvisis.cnv.filesys.Centroids;
import org.genvisis.cnv.filesys.MarkerDetailSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.var.SampleData;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.ext;

public class CNVannotator {

  /**
   * Create a new CNV file with BEAST height, score, and individual sex annotations for each CNV.
   * 
   * @param proj Project - Genvisis project file (including LRR values, BAFs, genotypes, X/Y, sex,
   *          batch, other annotation, etc.).
   * @param cnvFile CNV File to score, String, .cnv/.csv file (including FID, IID, CHR, BP1, BP2,
   *          TYPE, SCORE, SITES).
   * @param isSexCNVs This flag applies sex-specific centroids to recompute LRRs. If not, original
   *          LRRs will be used for scoring.
   */

  // Define "scoreCNVFile" method:
  public static void scoreCNVFile(Project proj, String cnvFile,
                                  boolean recomputeLRRsFromSexCentroids) {

    // Pull sample data (?) from project file:
    SampleData sd = proj.getSampleData(false);

    // Map with the indices of each marker in the project:
    MarkerDetailSet markerSet = proj.getMarkerSet();
    // Get chromosome number of markers:
    byte[] chr = markerSet.getChrs();
    // Get bp position of markers:
    int[] positions = markerSet.getPositions();
    // Get indices of set, in order of chromosome / position:
    int[][] indicesByChr = markerSet.getIndicesByChr();

    // Array of CNVs (for each individual) corresponding to regions (?) for which to calculate beast
    // scores:
    // Convert from txt to java object:
    // Beast works on single individual, want to break up CNVs
    CNVariant[] cnvs = CNVariant.loadPlinkFile(cnvFile);

    // Individual/sample IDs (?): Goes over every individual from above.
    HashSet<String> inds = new HashSet<>();
    for (CNVariant cnv : cnvs) {
      inds.add(cnv.getIndividualID());
    }
    // CNV IDs (?): Going over , map gets ID, returns list, adds to list. This individual has these
    // CNVs.
    ArrayList<String> ids = new ArrayList<>(inds);
    // 5/15

    HashMap<String, ArrayList<CNVariant>> cnvMap = new HashMap<>();
    // 5/15

    for (String s : ids) {
      cnvMap.put(s, new ArrayList<CNVariant>());
      // 5/15

    }
    // Adding to hashmap (of list of CNV for each INd, return CNV for an individual) existing T CNV,
    // associating w individual
    // Good place to add "new" CNVs.
    for (CNVariant cnv : cnvs) {
      cnvMap.get(cnv.getIndividualID()).add(cnv);
    }
    //
    CNVariant[][] cnvArr = new CNVariant[ids.size()][];

    for (int i = 0; i < cnvArr.length; i++) {
      cnvArr[i] = cnvMap.get(ids.get(i)).toArray(new CNVariant[0]);
    }

    // Initialize 3-dimensional arrays for (Male and Female) centroids:
    float[][][] centFem = null;
    float[][][] centMal = null;

    // If boolean "recomputeLRRsFromSexCentroids" selected true:
    if (recomputeLRRsFromSexCentroids) {
      // If (Female) sex-specific centroids file available:
      if (Files.exists(proj.SEX_CENTROIDS_FEMALE_FILENAME.getValue())) {
        // Sex-specific centroids (Female):
        centFem = Centroids.load(proj.SEX_CENTROIDS_FEMALE_FILENAME.getValue()).getCentroids();
        // No sex-specific (Female) file error:
      } else {
        proj.getLog()
            .reportError("Female-specific centroid file {"
                         + proj.SEX_CENTROIDS_FEMALE_FILENAME.getValue()
                         + "} doesn't exist - LRR correction cannot complete.");
        return;
      }
      // If (Male) sex-specific centroids file available:
      if (Files.exists(proj.SEX_CENTROIDS_MALE_FILENAME.getValue())) {
        // Sex-specific centroids (Male):
        centMal = Centroids.load(proj.SEX_CENTROIDS_MALE_FILENAME.getValue()).getCentroids();
        // No sex-specific (Male) file error:
      } else {
        proj.getLog()
            .reportError("Male-specific centroid file {"
                         + proj.SEX_CENTROIDS_MALE_FILENAME.getValue()
                         + "} doesn't exist - LRR correction cannot complete.");
        return;
      }
    }

    // Initialize array (size is the number of IDs read) to hold beast scores:
    BeastScore[] idScores = new BeastScore[ids.size()];
    // 5/15

    // Iterate over IIDs:
    for (int i = 0; i < ids.size(); i++) {
      // Get (partial) RAF file corresponding to each sample:
      Sample samp = proj.getPartialSampleFromRandomAccessFile(ids.get(i), false, false, false, true,
                                                              false);
      // Calculate LRRs (with sex-specific centroids):
      float[] lrrs = recomputeLRRsFromSexCentroids ? samp.getLRRs((sd.getSexForIndividual(ids.get(i)) == 1) ? centMal
                                                                                                            : centFem)
                                                   : samp.getLRRs();

      // Compute beast scores for CNVariant[] representing a single individual, return BeastScore
      // (for all the individual's cnvs):
      idScores[i] = BeastScore.beastInd(proj, sd, lrrs, cnvArr[i], chr, positions, indicesByChr);
      // 5/15

    }

    // Define output file name/extension:
    String outFile = ext.rootOf(cnvFile, false) + "_beast.cnv";
    // 5/15

    // Write output file:
    PrintWriter writer = Files.getAppropriateWriter(outFile);
    // 5/15

    // Specify/write column headers:
    writer.println(ArrayUtils.toStr(CNVariant.PLINK_CNV_HEADER, "\t") + "\tSEX\tSCORE\tHEIGHT");
    // 5/15

    for (int i = 0; i < ids.size(); i++) {
      // Include IDs (FID and IID (?)):
      CNVariant[] idCnvs = cnvArr[i];
      // 5/15

      // Derive / include beast scores:
      float[] scores = idScores[i].getBeastScores();
      // 5/15

      // Derive / include beast heights:
      float[] heights = idScores[i].getBeastHeights();
      // 5/15

      // Write data above:
      for (int c = 0; c < idCnvs.length; c++) {
        writer.println(idCnvs[c].toPlinkFormat() + "\t" + sd.getSexForIndividual(ids.get(i)) + "\t"
                       + scores[c] + "\t" + heights[c]);

      }

    }
    writer.flush();
    writer.close();
    // 5/15

  }

  public static void main(String[] args) {

    // Give current program directory (reference for self - TK):
    final String dir = System.getProperty("user.dir");
    System.out.println("TK note: current dir = " + dir);

    // CNVHelper to generate a regions file from CNV file (?):
    // CNVHelper.generateRegionsFileFromCNVFile("/Users/taylorkuebler/PennCNV_calls34.cnv");

    // Read in project .properties file:
    Project proj = new Project("C:\\Users\\Mark Hiner\\.genvisis\\projects\\SampleIllumina.properties");

    // Read in .cnv file:
    // String cnvFile = "/Users/taylorkuebler/PennCNV_calls34.cnv";
    // COPY OF ABOVE, FOR TESTING ONLY (SMALL SUBSET OF DATA FOR DEBUGGING)
    String cnvFile = "C:\\Users\\Mark Hiner\\Desktop\\Illumina\\Omni2_5_1000g_Genvisis\\cnvs\\genvisis.cnv";

    long t = System.currentTimeMillis();
    // scoreCNVfile to derive BEAST height, score, and individual sex annotations for each CNV,
    // store in new .cnv file:
    scoreCNVFile(proj, cnvFile, false);
    System.out.println("Took: " + (System.currentTimeMillis() - t) + "ms");

  }
}
