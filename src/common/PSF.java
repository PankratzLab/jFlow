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
	}

	/**
	 * Common commands to run
	 *
	 */
	public static class Cmd {
		/**
		 * Make sure to include an "&" after any profile.pl usage in a .pbs script
		 * It looks like this may actually profile all your jobs on a node
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

}
