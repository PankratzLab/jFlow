package stats;

//import java.io.*;
import java.util.*;

import common.*;

public class Rscript {
	public static final HashSet<String> R_INVALID_CHARS = new HashSet<String>(Arrays.asList("-"));
	public static final String R_REPLACEMENT = ".";
	public static final String[] RSCRIPT_EXECS = {
			// "/soft/R/3.0.1/bin/Rscript", // MSI
	"/soft/R/2.15.1/bin/Rscript", // Itasca nodes can only see this
	"/share/apps/R-3.0.1/bin/Rscript", // psych
	"/share/apps/src/R-3.0.1/bin/Rscript", // alcatraz
	};

	public static final String[] R_EXECS = {
			// "/soft/R/3.0.1/bin/R", // MSI
	"/soft/R/2.15.1/bin/R", // Itasca nodes can only see this
	"/share/apps/R-3.0.1/bin/R", // psych
	"/share/apps/src/R-3.0.1/bin/R", // alcatraz
	};

	public static String getRscriptExecutable(Logger log) {
		if (System.getProperty("os.name").startsWith("Windows")) {
			return "Rscript";
		} else {
			for (int i = 0; i < RSCRIPT_EXECS.length; i++) {
				if (Files.exists(RSCRIPT_EXECS[i])) {
					return RSCRIPT_EXECS[i];
				}
			}
		}

		log.report("Warning - could not determine hostname, assuming R executbale is simply 'Rscript'");

		return "Rscript";
	}

	public static String getRExecutable(Logger log) {
		if (System.getProperty("os.name").startsWith("Windows")) {
			return "R";
		} else {
			for (int i = 0; i < R_EXECS.length; i++) {
				if (Files.exists(R_EXECS[i])) {
					return R_EXECS[i];
				}
			}
		}

		log.report("Warning - could not determine hostname, assuming R executbale is simply 'Rscript'");

		return "R";
	}

	private static void batchUp(String dir, Logger log) {
		Vector<String> v;
		String[] files;
		String root;

		dir = ext.verifyDirFormat(dir);

		if (dir.equals("") || dir.equals("./")) {
			dir = ext.pwd();
		}

		files = Files.list(dir, ".R", false);

		v = new Vector<String>();
		for (int i = 0; i < files.length; i++) {
			root = ext.rootOf(dir + files[i], false);
			Files.qsub(root + ".qsub", getRscriptExecutable(log) + " --no-save " + root + ".R", 4000, 2, 1);
			v.add("qsub " + root + ".qsub");
		}

		Files.writeList(Array.toStringArray(v), dir + "master");
		Files.chmod(dir + "master");
	}

	/**
	 * @param toVector
	 *            String[] that will be converted to c("toVector[0]",toVector[1]...);
	 * @return
	 */
	public static String generateRVector(String[] toVector, boolean quotes) {
		String rV = "c(";
		for (int i = 0; i < toVector.length; i++) {
			rV += (i == 0 ? "" : ",");
			if (quotes) {
				rV += "\"" + toVector[i] + "\"";
			} else {
				rV += toVector[i];
			}
		}
		rV += ")";
		return rV;
	}

	/**
	 * @param s
	 * @return s with invalid r charcters replaced
	 */
	public static String makeRSafe(final String s) {
		String rSafe = s;
		for (String toReplace : R_INVALID_CHARS) {
			rSafe = rSafe.replaceAll(toReplace, R_REPLACEMENT);
		}
		return rSafe;
	}

	/**
	 * @param s
	 * @return s with invalid r charcters replaced
	 */
	public static String[] makeRSafe(final String[] s) {
		String[] rSafe = new String[s.length];
		for (int i = 0; i < s.length; i++) {
			rSafe[i] = makeRSafe(s[i]);
		}
		return rSafe;
	}

	/**
	 * Going to base a few R calls off this
	 */
	public interface RCommand {
		/**
		 * Should run an R command, returning true if successful
		 */
		public boolean execute();

		/**
		 * Should develop a String[] R script to be run on the command line
		 */
		public String[] developScript();

		/**
		 * perform any preliminary data checks
		 */
		public boolean validate();
	}

	public enum SCATTER_TYPE{
		LINE("geom_line"),POINT("geom_point");

		private String call;

		private SCATTER_TYPE(String call) {
			this.call = call;
		}

		public String getCall() {
			return call;
		}

		public void setCall(String call) {
			this.call = call;
		}
		

	}
	
	public enum LEGEND_POSITION{
		BOTTOM("legend.position=\"bottom\"");

		private String call;

		private LEGEND_POSITION(String call) {
			this.call = call;
		}

		public String getCall() {
			return call;
		}

		public void setCall(String call) {
			this.call = call;
		}
		

	}
	
	
	/**
	 * @author lane0212 For plotting that mimics Excel's scatter plots using ggplot2
	 */
	public static class RScatter implements RCommand {

		private static final String DATA_TABLE_VAR = "dataTable";
		private static final String DATA_TABLE_MELT_VAR = "dataTableMelt";
		private static final String DATA_TABLE_EXTRACT = "extract";

		private static final String PLOT_VAR = "p";

		private String dataFile, output, rScriptFile;
		private String[] dataYvalueColumns, rSafeYColumns;
		private String dataXvalueColumn, rSafeXColumn;
		private double[] xRange, yRange;
		private String xLabel;
		private String yLabel;
		private Logger log;
		private SCATTER_TYPE sType;
		private LEGEND_POSITION lPosition;
		private boolean valid;
		private int fontsize;
		private String title;

		public RScatter(String dataFile, String rSriptFile, String output, String dataXvalueColumn, String[] dataYvalueColumns, SCATTER_TYPE sType, Logger log) {
			super();
			this.dataFile = dataFile;
			this.rScriptFile = rSriptFile;
			this.dataXvalueColumn = dataXvalueColumn;
			this.dataYvalueColumns = dataYvalueColumns;
			this.log = log;
			this.output = output;
			this.sType = sType;
			this.rSafeXColumn = makeRSafe(dataXvalueColumn);
			this.rSafeYColumns = makeRSafe(dataYvalueColumns);
			this.valid = validate();
			this.fontsize=4;
			
		}
		
		

		public void setTitle(String title) {
			this.title = title;
		}


		public void setlPosition(LEGEND_POSITION lPosition) {
			this.lPosition = lPosition;
		}



		public void setxRange(double[] xRange) {
			this.xRange = xRange;
		}

		public void setyRange(double[] yRange) {
			this.yRange = yRange;
		}

		@Override
		public boolean execute() {
			String[] rScript = developScript();
			log.report(Array.toStr(rScript,"\n"));
			Files.writeList(rScript, rScriptFile);
			boolean ran = CmdLine.runCommandWithFileChecks(new String[]{"Rscript",rScriptFile}, "", new String[]{rScriptFile}, new String[]{output}, true, true, false, log);
			return ran;
		}

		@Override
		public String[] developScript() {

			if (valid) {
				ArrayList<String> rCmd = new ArrayList<String>();

				rCmd.add("library(scales)");
				rCmd.add("library(ggplot2)");
				rCmd.add("require(reshape2)");
				rCmd.add("library(plyr)");
				String dataTable = DATA_TABLE_VAR;
				rCmd.add(dataTable + "=read.table(\"" + dataFile + "\", header=T, sep=\"\\t\")");

				String plot = "";
				if (dataYvalueColumns.length > 1) {
					String[] toExtract = Array.concatAll(new String[] { rSafeXColumn }, rSafeYColumns);
					String extract = DATA_TABLE_EXTRACT + " <- " + DATA_TABLE_VAR + "[," + generateRVector(toExtract, true) + "]";
					rCmd.add(extract);
					String melt = DATA_TABLE_MELT_VAR + "<-  melt(" + DATA_TABLE_EXTRACT + ",id.vars =\"" + rSafeXColumn + "\")";
					rCmd.add(melt);
					dataTable = DATA_TABLE_MELT_VAR;
					plot += PLOT_VAR + " <- ggplot(" + dataTable + ",aes(x=" + rSafeXColumn + ", y=value, group=variable)) ";
					plot += " + " + sType.getCall() + "(aes(colour =variable))";
				
				}else{
					plot += PLOT_VAR + " <- ggplot(" + dataTable + ", aes(x = " + dataTable + "$" + rSafeXColumn + ")) ";
					plot += " + " + sType.getCall() + "(aes(y = " + dataTable + "$" + rSafeYColumns[0] + "))";
				}
				
			
				
//				
//				plot += PLOT_VAR + " <- ggplot(" + dataTable + ", aes(x = " + dataTable + "$" + rSafeXColumn + ")) ";
//				for (int i = 0; i < rSafeYColumns.length; i++) {
//					plot += " + " + sType.getCall() + "(aes(y = " + dataTable + "$" + rSafeYColumns[i] + "))";
//				}
				if (xLabel != null) {
					plot += " + xlab(\"" + xLabel + "\")";
				}
				if (yLabel != null) {
					plot += " + ylab(\"" + yLabel + "\")";
				}
				if (xRange != null) {
					plot += " + xlim(" + xRange[0] + "," + xRange[1] + ")";
				}
				if (yRange != null) {
					plot += " + ylim(" + yRange[0] + "," + yRange[1] + ")";
				}
				plot += " + theme(axis.line = element_line(colour = \"black\"), ";
				plot += "panel.grid.major = element_blank(), ";
				plot += "panel.grid.minor = element_blank(),";
				plot += "panel.border = element_blank(),";
				plot += yLabel == null ? "axis.title.y = element_blank()," : "";
				plot += xLabel == null ? "axis.title.x = element_blank()," : "";
				plot += lPosition == null ? "" : lPosition.getCall() + ",";
				plot += "legend.text=element_text(size=" + fontsize + "),";
				plot += "panel.background = element_blank())";
				plot += title == null ? "" : "+ labs(title = \"" + title + "\")";
				rCmd.add(plot);
				rCmd.add("ggsave(file=\"" + output + "\", width= 8, height = 4)");

				return rCmd.toArray(new String[rCmd.size()]);

			}
			return null;
		}

		@Override
		public boolean validate() {
			boolean valid = true;
			if (!Files.exists(dataFile)) {
				valid = false;
				log.reportFileNotFound(dataFile);
			} else {

				if (!Files.headerOfFileContainsAll(dataFile, dataYvalueColumns, log)) {
					valid = false;
					log.reportTimeError("Could not find all Y value columns in " + dataFile);
				} else {
					if (!Files.headerOfFileContainsAll(dataFile, new String[] { dataXvalueColumn }, log)) {
						valid = false;
						log.reportTimeError("Could not find all X value column in " + dataFile);
					} else {
						if (xRange != null && xRange.length != 2) {
							valid = false;
							log.reportTimeError("x range must be length 2");
						} else {
							if (yRange != null && yRange.length != 2) {
								valid = false;
								log.reportTimeError("y range must be length 2");
							}
						}
					}
				}

			}
			return valid;
		}

		public void setxLabel(String xLabel) {
			this.xLabel = xLabel;
		}

		public void setyLabel(String yLabel) {
			this.yLabel = yLabel;
		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		boolean batchUp = false;
		String dir = "./";
		String logfile = null;
		Logger log;

		String usage = "\n" + "stats.Rscript requires 0-1 arguments\n" + "   (1) create qsub files for all .R files in the directory (i.e. -batchUp (not the default))\n" + "   (2) directory to search for the .R files (i.e. dir=" + dir + " (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("-batchUp")) {
				batchUp = true;
				numArgs--;
			} else if (args[i].startsWith("dir=")) {
				dir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			log = new Logger(logfile);
			if (batchUp) {
				batchUp(dir, log);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
