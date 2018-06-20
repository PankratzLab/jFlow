package org.genvisis.cnv.annotator;

import java.io.IOException;
import java.util.List;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.filesys.CNVariant;

public class CNVAnnoTester {

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

    long t = System.currentTimeMillis();
    // Take CNVs in original .cnv file list, convert them to segments:
    // Get a list of CNVariants:
    List<CNVariant> cnvs = CNVariant.loadPlinkFile(cnvFile, null, true);

    // 1) Gather the segments into an AnnotatedCollection:
    AnnotatorConfig config = new AnnotatorConfig().setDoUpstream().setDoDownstream()
                                                  .setUseSegmentSize();
    AnnotatedCollection<CNVariant> annotatedCollection = new AnnotatedCollection<>(cnvs, config);

    // 2) Create a BeastScoreAnnotator and annotate the collection
    BeastScoreAnnotator beastScoreAnnotator = new BeastScoreAnnotator(proj);

    // 3) Use the annotator on the collection:
    beastScoreAnnotator.annotate(annotatedCollection);

    System.out.println("Took: " + (System.currentTimeMillis() - t) + "ms");

    try {
      annotatedCollection.writeResults("C:\\Users\\Mark Hiner\\Desktop\\Illumina\\Omni2_5_1000g_Genvisis\\cnvs\\beastTest.annotated.cnv");
    } catch (IOException e) {
      e.printStackTrace();
    }
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

    // scoreCNVfile to derive BEAST height, score, and individual sex annotations for each CNV,
    // store in new .cnv file:
    scoreCNVFile(proj, cnvFile, false);

  }
}
