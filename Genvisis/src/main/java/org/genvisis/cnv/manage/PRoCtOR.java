package org.genvisis.cnv.manage;

import java.io.File;
import java.io.FilenameFilter;

import org.genvisis.cnv.analysis.PennCNVPrep;
import org.genvisis.cnv.analysis.pca.PCA;
import org.genvisis.cnv.analysis.pca.PCAPrep;
import org.genvisis.cnv.analysis.pca.PrincipalComponentsApply;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.common.ext;
import org.genvisis.stats.LeastSquares.LS_TYPE;

// PRincipal COmponents Residuals - PR [o] C[t]O R
public class PRoCtOR {

  private static final String SHADOW_PREP_DIR = "shadowPrep/";

  private static long getMaxSampleSize(Project proj) {
    File[] sampleFiles =
        (new File(proj.SAMPLE_DIRECTORY.getValue())).listFiles(new FilenameFilter() {
          @Override
          public boolean accept(File dir, String name) {
            return name.endsWith(Sample.SAMPLE_FILE_EXTENSION);
          }
        });
    long max = 0;
    for (File f : sampleFiles) {
      if (f.length() > max) {
        max = f.length();
      }
    }
    return max;
  }

  private static int getSampleChunks(Project proj, int numThreads) {
    /*
     * sampleChunks = (.8 * totalRam * sampleSize) / numThreads;
     */
    long mem = Runtime.getRuntime().maxMemory();
    long samp = getMaxSampleSize(proj);
    double sampleChunks = (0.8 * mem * numThreads) / samp;
    return (int) sampleChunks;
  }

  public static void shadow(Project proj, String tmpDir, String outputBase,
                            double markerCallRateFilter, boolean recomputeLRR_PCs,
                            int numComponents, int totalThreads) {
    int numMarkerThreads = 2;
    int numThreads = (int) Math.ceil((double) totalThreads / (double) numMarkerThreads);
    boolean markerQC = true;
    String useFile = null;
    int sampleChunks = getSampleChunks(proj, numThreads);
    proj.getLog().report("Using " + sampleChunks + " sample chunks");

    int retCode = PCAPrep.prepPCA(proj, numThreads, outputBase, markerQC, markerCallRateFilter,
                                  useFile, proj.getSampleList(), proj.getLog());
    if (retCode != 42) {
      // TODO error
      return;
    }
    PrincipalComponentsApply pcApply =
        PCA.generateFullPCA(proj, numComponents, outputBase, recomputeLRR_PCs, true, null,
                            proj.getLog());
    pcApply.getExtrapolatedPCsFile();
    PennCNVPrep.prepExport(proj, SHADOW_PREP_DIR, tmpDir, numComponents, null, numThreads,
                           numMarkerThreads, LS_TYPE.REGULAR, false);
    PennCNVPrep.exportSpecialPennCNV(proj, SHADOW_PREP_DIR, tmpDir, numComponents, null, numThreads,
                                     numMarkerThreads, true, LS_TYPE.REGULAR, sampleChunks, false);
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = "D:/projects/Poynter.properties";
    String tempDir = null;
    String outputBase = MitoPipeline.FILE_BASE;
    double callrate = MitoPipeline.DEFAULT_MKR_CALLRATE_FILTER;
    boolean recomputeLRR = false;
    int numComponents = MitoPipeline.DEFAULT_NUM_COMPONENTS;
    int numThreads = Runtime.getRuntime().availableProcessors();

    String usage = "\n" + "cnv.manage.PRoCtOR requires 0-1 arguments\n"
                   + "   (1) project properties filename (i.e. proj=" + filename + " (default))\n"
                   + "   (2) Number of principal components for correction (i.e. numComponents="
                   + numComponents + " (default))\n"
                   + "   (3) Output file full path and baseName for principal components correction files (i.e. outputBase="
                   + outputBase + " (default))\n"
                   + "   (4) Call-rate filter for determining high-quality markers (i.e. callrate="
                   + callrate + " (default))\n"
                   + "   (5) Flag specifying whether or not to re-compute Log-R Ratio values (usually false if LRRs already exist) (i.e. recomputeLRR="
                   + recomputeLRR + " (default))\n"
                   + "   (6) Total number of threads to use (i.e. numThreads=" + numThreads
                   + " (default))\n"
                   + "   (7) OPTIONAL: temp directory for intermediate files (which tend to be very large) (i.e. tmp="
                   + tempDir + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("proj=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("numComponents=")) {
        numComponents = ext.parseIntArg(arg);
        numArgs--;
      } else if (arg.startsWith("outputBase=")) {
        outputBase = ext.parseStringArg(arg, outputBase);
        numArgs--;
      } else if (arg.startsWith("callrate=")) {
        callrate = ext.parseDoubleArg(arg);
        numArgs--;
      } else if (arg.startsWith("recomputeLRR=")) {
        recomputeLRR = ext.parseBooleanArg(arg);
        numArgs--;
      } else if (arg.startsWith("numThreads=")) {
        numThreads = ext.parseIntArg(arg);
        numArgs--;
      } else if (arg.startsWith("tmp=")) {
        tempDir = ext.parseStringArg(arg, null);
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
      Project proj = new Project(filename, false);
      shadow(proj, tempDir, outputBase, callrate, recomputeLRR, numComponents, numThreads);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

}
