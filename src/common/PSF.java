package common;

/**
 * Just some common public static final vars to be used anywhere...used to build external commands, etc.. Might not be worth it but oh well;
 *
 */
public class PSF {
	/**
	 * related to java commands (-Xmx, -cp etc)
	 *
	 */
	public static class Java {
		public static final String JAVA = "java";
		public static final String JAR = "-jar";
		public static final String XMX = "-Xmx";
		public static final String CP = "-cp";
		public static final String GENVISIS = "park.jar";

		public static String[] buildJava() {
			return new String[] { JAVA };
		}

		public static String[] buildJavaJar(String jarFile) {
			String[] javaJar = Array.concatAll(buildJava(), new String[] { JAR });
			if (jarFile != null) {
				javaJar = Array.concatAll(javaJar, new String[] { jarFile });
			}
			return javaJar;
		}

		public static String[] buildJavaCP(String fullPathTojarFile) {
			String[] javaJar = Array.concatAll(buildJava(), new String[] { CP });
			if (fullPathTojarFile != null) {
				javaJar = Array.concatAll(javaJar, new String[] { fullPathTojarFile });
			}
			return javaJar;
		}

		public static String[] buildJavaCPXMX(String fullPathTojarFile, String cp, int memoryInMb) {
			String[] javaGenCP = buildJavaCP(fullPathTojarFile);
			String[] add = new String[] { XMX + memoryInMb + "m", cp };
			return Array.concatAll(javaGenCP, add);

		}
	}

	/**
	 * related to loading modules for qsubs
	 *
	 */
	public static class Load {
		public static final String MODULE_LOAD_JAVA = "module load java";
		public static final String MODULE_LOAD_SAMTOOLS = "module load samtools";
		public static final String MODULE_LOAD_RISS_UTIL = "module load riss_util";
		public static final String MODULE_LOAD_R = "module load R";
		public static final String MODULE_LOAD_PERL = "module load perl";

		public static String[] getAllModules() {
			return new String[] { MODULE_LOAD_JAVA, MODULE_LOAD_PERL, MODULE_LOAD_R, MODULE_LOAD_RISS_UTIL, MODULE_LOAD_SAMTOOLS };
		}

	}

	/**
	 * Common commands to run
	 *
	 */
	public static class Cmd {
		/**
		 * Make sure to include an "&" after any profile.pl usage in a .pbs script It looks like this may actually profile all your jobs on a node
		 */
		public static final String PROFILE_PL = "profile.pl";

		public static String getProfilePLWithOutput(String outputFile) {
			return PROFILE_PL + " -o " + outputFile + " " + Ext.AMP;
		}

		public static final String ECHO = "echo";

		/**
		 * @param toEcho
		 *            get command to echo this string
		 * @return the echo command
		 */
		public static final String echoAString(String toEcho) {
			return ECHO + "\"" + toEcho + "\"";
		}

		public static final String getSedCommand(String fullPathToInput, String fullPathToOutput, String reg1, String reg2) {
			String sed = "";
			// sed 's/1000g2014oct_all/g10002014oct_all/g'

			return sed;
		}

		public static final String PERL = "perl";

		public static final String SED = "sed";

	}

	/**
	 * related to anything
	 *
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
		public static final String WALLTIME_HRS = "wallTimeInHour=";

		public static String getNumThreadsCommand(int argNumber, int numThreads) {
			return "   (" + argNumber + ")" + " number of threads to use (i.e. " + NUM_THREADS_COMMAND + numThreads + " (default))\n";
		}

		public static String getWallTimeCommand(int argNumber, int wallTimeInHours) {
			return "   (" + argNumber + ")" + " wall time in hours to use (i.e. " + WALLTIME_HRS + wallTimeInHours + " (default))\n";
		}

		public static String getMemoryMbCommand(int argNumber, int memoryInMb) {
			return "   (" + argNumber + ")" + " memory in mb to use (i.e. " + MEMORY_MB + memoryInMb + " (default))\n";
		}

		public static String getOutputDirCommand(int argNumber, String defaultDir) {
			return "   (" + argNumber + ")" + " the output directory to use (i.e. " + OUTPUT_DIR_COMMAND + (defaultDir == null ? "" : defaultDir) + " (" + (defaultDir == null ? "no" : "") + "default))\n";

		}

	}

	/**
	 * For building plink commands
	 *
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
		
		public static boolean allFilesExist(String plinkDirAndRoot, boolean bed) {
		    String[] files;
		    if (bed) {
		        files = new String[]{getBED(plinkDirAndRoot), getBIM(plinkDirAndRoot), getFAM(plinkDirAndRoot)};
		    } else {
    		    files = new String[]{getPED(plinkDirAndRoot), getMAP(plinkDirAndRoot)};
    		}
		    for (String file : files) {
		        if (!Files.exists(file)) {
		            return false;
		        }
		    }
		    return true;
		}
		
		/**
		 * Get a command that will convert a vcf to plink format using plink2
		 * 
		 * @param inputVCF
		 * @param outputBase
		 * @return
		 */
		public static String[] getPlinkVCFCommand(String inputVCF, String outputBase) {
			return new String[] { PLINK2, VCF, inputVCF, DOUBLE_ID, MAKE_BED, OUT, outputBase, BIALLELIC_ONLY, STRICT, LIST, NO_WEB };
		}

		/**
		 * @param root
		 *            likely full path
		 * @param outputBase
		 *            likely full path
		 * @return
		 */
		// CmdLine.run("plink --bfile plink --indep-pairwise 50 5 0.3", dir+"ldPruning/");

		public static String[] getPlinkLDCommand(String root) {
			return new String[] { PLINK, BFILE, root, INDEP_PAIRWISE, "50", "5", "0.3" };
		}

		/**
		 * Get the .bim, .bed, and .fam files associated with this root
		 * 
		 * @param root
		 * @return
		 */
		public static String[] getPlinkBedBimFam(String root) {
			return new String[] { getBED(root), getBIM(root), getFAM(root) };
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

	}

}
