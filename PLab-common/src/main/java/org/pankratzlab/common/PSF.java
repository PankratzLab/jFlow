package org.pankratzlab.common;

import java.awt.Color;
import java.util.Set;
import com.google.common.collect.Sets;

/**
 * Just some common public static final vars to be used anywhere...used to build external commands,
 * etc.. Might not be worth it but oh well;
 */
public class PSF {

  /**
   * Check if the current thread has been interrupted and, if so, execute the given runnable and
   * throw a RuntimeException (wrapped around an InterruptedException).
   * 
   * @param runIfInterrupt
   */
  public static final void checkInterrupted(Runnable runIfInterrupt) {
    if (Thread.currentThread().isInterrupted()) {
      runIfInterrupt.run();
      throw new RuntimeException(new InterruptedException());
    }
  }

  /**
   * Check if the current thread has been interrupted and, if so, throw a RuntimeException (wrapped
   * around an InterruptedException).
   */
  public static final void checkInterrupted() {
    if (Thread.currentThread().isInterrupted()) {
      throw new RuntimeException(new InterruptedException());
    }
  }

  /**
   * Stores colors rather than recreating them for each use.<br />
   * Names and hue designations are from http://www.color-blindness.com/color-name-hue/<br />
   * Hues are approximate - feel free to move these around if you think they're incorrectly
   * categorized.
   */
  public static final class Colors {

    public static final class REDS {

      public static final Color VENETIAN_RED = new Color(201, 30, 10);
      public static final Color CINNABAR = new Color(230, 73, 39);
      public static final Color FIREBRICK = new Color(178, 34, 34);

    }

    public static final class ORANGES {

      public static final Color SAFETY_ORANGE = new Color(255, 100, 0);
      public static final Color INDOCHINE = new Color(153, 102, 51);
      public static final Color SUNSHADE = new Color(247, 150, 70);
      public static final Color MANGO_TANGO = new Color(227, 108, 9);

    }

    public static final class YELLOWS {

      public static final Color AMBER = new Color(255, 192, 0);

    }

    public static final class GREENS {

      public static final Color GREEN = new Color(33, 87, 0);
      public static final Color GREEN_YELLOW = new Color(189, 243, 61);
      public static final Color FREE_SPEECH_AQUAMARINE = new Color(0, 157, 126);
      public static final Color DARK_CYAN = new Color(0, 150, 150);

    }

    public static final class BLUES {

      public static final Color MIDNIGHT_EXPRESS = new Color(33, 31, 53);
      public static final Color PERSIAN_BLUE = new Color(23, 58, 172);
      public static final Color DODGER_BLUE = new Color(55, 129, 252);
      public static final Color SLATE_BLUE = new Color(94, 88, 214);
      public static final Color NAVY = new Color(0, 0, 128);
      public static final Color CORNFLOWER_BLUE = new Color(100, 149, 237);
      public static final Color DARK_SLATE_BLUE = new Color(72, 61, 139);
      public static final Color MEDIUM_SLATE_BLUE = new Color(123, 104, 238);
      public static final Color LIGHT_SLATE_BLUE = new Color(132, 112, 255);
      public static final Color MEDIUM_BLUE = new Color(0, 0, 205);
      public static final Color ROYAL_BLUE = new Color(65, 105, 225);
      public static final Color DEEP_SKY_BLUE = new Color(0, 191, 255);
      public static final Color LIGHT_SKY_BLUE = new Color(135, 206, 250);
      public static final Color STEEL_BLUE = new Color(70, 130, 180);
      public static final Color LIGHT_STEEL_BLUE = new Color(176, 196, 222);
      public static final Color LIGHT_BLUE = new Color(173, 216, 230);
      public static final Color POWDER_BLUE = new Color(176, 224, 230);
      public static final Color PALE_TURQUOISE = new Color(175, 238, 238);
      public static final Color DARK_TURQUOISE = new Color(0, 206, 209);
      public static final Color MEDIUM_TURQUOISE = new Color(72, 209, 204);
      public static final Color TURQUOISE = new Color(64, 224, 208);
      public static final Color AQUA = Color.CYAN;
      public static final Color LIGHT_CYAN = new Color(224, 255, 255);

    }

    public static final class VIOLETS {

      public static final Color BLUE_VIOLET = new Color(140, 20, 180);
      public static final Color ORCHID = new Color(217, 109, 194);
      public static final Color WINDSOR = new Color(62, 46, 133);
      public static final Color LOLA = new Color(189, 174, 198);
      public static final Color CHRISTALLE = new Color(66, 28, 82);

    }
    public static final class BROWNS {

      public static final Color MAROON = new Color(100, 50, 0);

    }

    public static final class BLACKS {

      /**
       * Slightly more blue than black (rgb 33,31,53)
       */

    }

    public static final class GREYS {

    }

    public static final class WHITES {

    }

  }

  /**
   * related to java commands (-Xmx, -cp etc)
   */
  public static class Java {

    public static final String JAVA = "java";
    public static final String JAR = "-jar";
    public static final String XMX = "-Xmx";
    public static final String CP = "-cp";
    public static final String GENVISIS = "genvisis.jar";

    public static String[] buildJava() {
      return new String[] {JAVA};
    }

    public static String[] buildJavaJar(String jarFile) {
      return new String[] {JAVA, JAR, jarFile};
    }

    public static String[] buildJavaXMXJar(String jarFile, int memoryInMb) {
      return new String[] {JAVA, buildXmxString(memoryInMb), JAR, jarFile};
    }

    public static String[] buildJavaCP(String fullPathTojarFile) {
      return new String[] {JAVA, CP, fullPathTojarFile};
    }

    public static String[] buildJavaCPXMX(String fullPathTojarFile, String cp, int memoryInMb) {
      return new String[] {JAVA, buildXmxString(memoryInMb), CP, cp, JAR, fullPathTojarFile};
    }

    public static String buildXmxString(int memoryInMb) {
      return XMX + memoryInMb + "m";
    }
  }

  /**
   * related to loading modules for qsubs
   */
  public static class Load {

    public static final String MODULE_LOAD_JAVA = "module load java";
    public static final String MODULE_LOAD_RISS_UTIL = "module load riss_util";
    public static final String MODULE_LOAD_R = "module load R";
    public static final String MODULE_LOAD_PERL = "module load perl";
    public static final String MODULE_LOAD_PYTHON = "module load python";

    public static final String MODULE_LOAD_NETCDF = "module load netcdf";
    public static final String MODULE_LOAD_SINGULARITY = "module load singularity";

    public static String[] getAllModules() {
      return new String[] {MODULE_LOAD_JAVA, MODULE_LOAD_PERL, MODULE_LOAD_R, MODULE_LOAD_RISS_UTIL,
                           MODULE_LOAD_PYTHON, MODULE_LOAD_NETCDF, MODULE_LOAD_SINGULARITY};
    }

  }

  /**
   * Common commands to run
   */
  public static class Cmd {

    /**
     * Make sure to include an "&" after any profile.pl usage in a .pbs script It looks like this
     * may actually profile all your jobs on a node
     */
    public static final String PROFILE_PL = "profile.pl";

    public static String getProfilePLWithOutput(String outputFile) {
      return PROFILE_PL + " -o " + outputFile + " " + Ext.AMP;
    }

    public static final String ECHO = "echo";

    /**
     * @param toEcho get command to echo this string
     * @return the echo command
     */
    public static final String echoAString(String toEcho) {
      return ECHO + "\"" + toEcho + "\"";
    }

    public static final String getSedCommand(String fullPathToInput, String fullPathToOutput,
                                             String reg1, String reg2) {
      String sed = "";
      // sed 's/1000g2014oct_all/g10002014oct_all/g'

      return sed;
    }

    public static final String PERL = "perl";

    public static final String SED = "sed";

  }

  /**
   * related to anything
   */
  public static class Ext {

    public static final String CARROT = ">";
    public static final String REV_CARROT = "<";
    public static final String SPACE = " ";

    public static final String AMP = "&";
    public static final String BLANK = "";
    public static final String NUM_THREADS_COMMAND = "threads=";
    public static final String OUTPUT_DIR_COMMAND = "outputdir=";
    public static final String MEMORY_MB = "memoryInMb=";
    public static final int DEFAULT_MEMORY_MB = 61000;
    public static final String WALLTIME_HRS = "wallTimeInHour=";

    public static String getNumThreadsCommand(int argNumber, int numThreads) {
      return "   (" + argNumber + ")" + " number of threads to use (i.e. " + NUM_THREADS_COMMAND
             + numThreads + " (default))\n";
    }

    public static String getWallTimeCommand(int argNumber, int wallTimeInHours) {
      return "   (" + argNumber + ")" + " wall time in hours to use (i.e. " + WALLTIME_HRS
             + wallTimeInHours + " (default))\n";
    }

    public static String getMemoryMbCommand(int argNumber, int memoryInMb) {
      return "   (" + argNumber + ")" + " memory in mb to use (i.e. " + MEMORY_MB + memoryInMb
             + " (default))\n";
    }

    public static String getOutputDirCommand(int argNumber, String defaultDir) {
      return "   (" + argNumber + ")" + " the output directory to use (i.e. " + OUTPUT_DIR_COMMAND
             + (defaultDir == null ? "" : defaultDir) + " (" + (defaultDir == null ? "no" : "")
             + "default))\n";

    }

  }

  public static class Regex {

    public static final String GREEDY_WHITESPACE = "[\\s]+";

    public static String regexSplitPreserveQuoted(String baseRegex) {
      return "[" + baseRegex + "]+(?=([^\"]*\"[^\"]*\")*[^\"]*$)";
    }
  }

  /**
   * For building plink commands
   */
  public static class Plink {

    public static final String PLINK2 = "plink2";
    public static final String PLINK = "plink";

    public static final String BFILE = "--bfile";
    public static final String INDEP_PAIRWISE = "--indep-pairwise";
    public static final String MAKE_BED = "--make-bed";
    public static final String BIALLELIC_ONLY = "--biallelic-only";
    public static final String STRICT = "strict";
    public static final String LIST = "list";

    public static final String DOUBLE_ID = "--double-id";
    public static final String NO_WEB = "--noweb";
    public static final String VCF = "--vcf";
    public static final String OUT = "--out";
    public static final String BIM = ".bim";
    public static final String BED = ".bed";
    public static final String FAM = ".fam";
    public static final String PED = ".ped";
    public static final String MAP = ".map";

    public static final int FAM_FID_INDEX = 0;
    public static final int FAM_IID_INDEX = 1;
    public static final int FAM_FA_INDEX = 2;
    public static final int FAM_MO_INDEX = 3;
    public static final int FAM_SEX_INDEX = 4;
    public static final int FAM_AFF_INDEX = 5;

    public static final int FAM_FIELD_COUNT = 6;

    public static final int AFF_CASE_STATUS = 2;

    public static final int BIM_CHR_INDEX = 0;
    public static final int BIM_NAME_INDEX = 1;
    public static final int BIM_CM_INDEX = 2;
    public static final int BIM_POS_INDEX = 3;
    public static final int BIM_A1_INDEX = 4;
    public static final int BIM_A2_INDEX = 5;

    public static boolean affIsCase(String affField) {
      final int aff;
      try {
        aff = Integer.parseInt(affField);
      } catch (NumberFormatException e) {
        return false;
      }
      return aff == AFF_CASE_STATUS;
    }

    public static boolean bedBimFamExist(String plinkDirAndRoot) {
      return allFilesExist(plinkDirAndRoot, true);
    }

    public static boolean allFilesExist(String plinkDirAndRoot, boolean bed) {
      String[] files;
      if (bed) {
        files = getPlinkBedBimFam(plinkDirAndRoot);
      } else {
        files = new String[] {getPED(plinkDirAndRoot), getMAP(plinkDirAndRoot)};
      }
      return Files.checkAllFiles("", files);
    }

    /**
     * Get a command that will convert a vcf to plink format using plink2
     *
     * @param inputVCF
     * @param outputBase
     * @return
     */
    public static String[] getPlinkVCFCommand(String inputVCF, String outputBase) {
      return new String[] {PLINK2, VCF, inputVCF, DOUBLE_ID, MAKE_BED, OUT, outputBase,
                           BIALLELIC_ONLY, STRICT, LIST, NO_WEB};
    }

    /**
     * @param root likely full path
     * @param outputBase likely full path
     * @return
     */
    // CmdLine.run("plink --bfile plink --indep-pairwise 50 5 0.3", dir+"ldPruning/");

    public static String[] getPlinkLDCommand(String root) {
      return new String[] {PLINK, BFILE, root, INDEP_PAIRWISE, "50", "5", "0.3"};
    }

    /**
     * Get the .bim, .bed, and .fam files associated with this root
     *
     * @param root
     * @return
     */
    public static String[] getPlinkBedBimFam(String root) {
      return new String[] {getBED(root), getBIM(root), getFAM(root)};
    }

    public static Set<String> getPlinkBedBimFamSet(String root) {
      Set<String> bedBimFam = Sets.newHashSetWithExpectedSize(3);
      bedBimFam.add(getBED(root));
      bedBimFam.add(getBIM(root));
      bedBimFam.add(getFAM(root));
      return bedBimFam;
    }

    public static String getPED(String root) {
      return root + PED;
    }

    public static String getMAP(String root) {
      return root + MAP;
    }

    public static String getBED(String root) {
      return root + BED;
    }

    public static String getBIM(String root) {
      return root + BIM;
    }

    public static String getFAM(String root) {
      return root + FAM;
    }

    public static final String[] LOGISTIC_SE_HEADER = {"CHR", "SNP", "BP", "A1", "TEST", "NMISS",
                                                       "OR", "SE", "L95", "U95", "STAT", "P"};
    public static final String[] LINEAR_SE_HEADER = {"CHR", "SNP", "BP", "A1", "TEST", "NMISS",
                                                     "BETA", "SE", "L95", "U95", "STAT", "P"};
    public static final String FAMID = "FID";
    public static final String INDID = "IID";
    public static final String FATHER_ID = "FA_IID";
    public static final String MOTHER_ID = "MO_IID";
    public static final String IND_DNA = "DNA";
    public static final String FATHER_DNA = "FA_DNA";
    public static final String MOTHER_DNA = "MO_DNA";

    public static boolean isValidDNA(String s) {
      return !s.isEmpty() && !s.equals("0") && !s.equals(".");
    }

  }

}
