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
		
		public static final String PERL = "perl";

	}

	/**
	 * related to anything
	 *
	 */
	public static class Ext {
		public static final String CARROT = ">";
		public static final String AMP = "&";
		public static final String BLANK = "";
		public static final String NUM_THREADS_COMMAND = "numthreads=";
		public static final String OUTPUT_DIR_COMMAND = "outputdir=";

	}

	/**
	 * For building plink commands
	 *
	 */
	public static class Plink {
		public static final String PLINK2 = "plink2";
		public static final String BFILE = "--bfile";
		public static final String MAKE_BED = "--make-bed";
		public static final String DOUBLE_ID = "--double-id";
		public static final String NO_WEB = "--noweb";
		public static final String VCF = "--vcf";
		public static final String OUT = "--out";
		public static final String BIM = ".bim";
		public static final String BED = ".bed";
		public static final String FAM = ".fam";

		/**
		 * Get a command that will convert a vcf to plink format using plink2
		 * 
		 * @param inputVCF
		 * @param outputBase
		 * @return
		 */
		public static String[] getPlinkVCFCommand(String inputVCF, String outputBase) {
			return new String[] { PLINK2, VCF, inputVCF, DOUBLE_ID, MAKE_BED, OUT, outputBase, NO_WEB };
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
