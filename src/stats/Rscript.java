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

	public enum SCATTER_TYPE {
		LINE("geom_line"), POINT("geom_point");

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

	public enum LEGEND_POSITION {
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

	public enum PLOT_DEVICE {
		PDF("pdf(file="), JPEG("jpeg");

		private String call;

		private PLOT_DEVICE(String call) {
			this.call = call;
		}

		public String getCall() {
			return call;
		}

		public void setCall(String call) {
			this.call = call;
		}
	}

	public enum GEOM_POINT_SIZE {
		GEOM_POINT("size =", 1);

		private String call;
		private double size;

		private GEOM_POINT_SIZE(String call, double size) {
			this.call = call;
			this.size = size;
		}

		public String getCall() {
			return call + size;
		}

		public void setSize(double size) {
			this.size = size;
		}
	}

	public enum COLUMNS_MULTIPLOT {
		COLUMNS_MULTIPLOT_1("cols =", 1), COLUMNS_MULTIPLOT_2("cols =", 2);

		private String call;
		private int cols;

		private COLUMNS_MULTIPLOT(String call, int cols) {
			this.call = call;
			this.cols = cols;
		}

		public String getCall() {
			return call + cols;
		}

		public void setSize(int cols) {
			this.cols = cols;
		}

	}

	public static class RScatters implements RCommand {

		private RScatter[] rScatters;
		private String mergeOutput;
		private String rScriptFile;
		private COLUMNS_MULTIPLOT cMultiplot;
		private PLOT_DEVICE device;
		private Logger log;

		public RScatters(RScatter[] rScatter, String rScriptFile, String mergeOutput, COLUMNS_MULTIPLOT cMultiplot, PLOT_DEVICE device, Logger log) {
			super();
			this.rScatters = rScatter;
			this.mergeOutput = mergeOutput;
			this.cMultiplot = cMultiplot;
			this.device = device;
			this.rScriptFile = rScriptFile;
			this.log = log;
		}

		public COLUMNS_MULTIPLOT getcMultiplot() {
			return cMultiplot;
		}

		@Override
		public boolean execute() {
			String[] rScript = developScript();
			log.report(Array.toStr(rScript, "\n"));
			Files.writeList(rScript, rScriptFile);
			boolean ran = CmdLine.runCommandWithFileChecks(new String[] { "Rscript", rScriptFile }, "", new String[] { rScriptFile }, new String[] { mergeOutput }, true, true, false, log);
			return ran;
		}

		@Override
		public String[] developScript() {
			ArrayList<String> allPlots = new ArrayList<String>();
			allPlots.add(getMultiPlotFunc());

			String plotDevice = device.getCall() + "\"" + mergeOutput + "\" ,onefile = TRUE)";
			allPlots.add(plotDevice);
			for (int i = 0; i < rScatters.length; i++) {
				String[] tmpScript = rScatters[i].developScript();
				for (int j = 0; j < tmpScript.length; j++) {
					allPlots.add(tmpScript[j]);
				}
			}
			for (int i = 0; i < rScatters.length; i++) {
				allPlots.add(rScatters[i].getPlotVar());
			}
			allPlots.add("dev.off()");
			return allPlots.toArray(new String[allPlots.size()]);
		}

		@Override
		public boolean validate() {
			return false;
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
		private double yMin;
		private String xLabel;
		private String yLabel;
		private Logger log;
		private SCATTER_TYPE sType;
		private LEGEND_POSITION lPosition;
		private GEOM_POINT_SIZE gPoint_SIZE;
		private boolean valid;
		private int fontsize;
		private String title;
		private String plotVar;
		private boolean logTransformX;
		private boolean logTransformY;
		private boolean overWriteExisting;

		public RScatter(String dataFile, String rSriptFile, String plotVar, String output, String dataXvalueColumn, String[] dataYvalueColumns, SCATTER_TYPE sType, Logger log) {
			super();
			this.dataFile = dataFile;
			this.rScriptFile = rSriptFile;
			this.dataXvalueColumn = dataXvalueColumn;
			this.dataYvalueColumns = dataYvalueColumns;
			this.log = log;
			this.output = ext.rootOf(output, false) + ".jpeg";
			this.sType = sType;
			this.rSafeXColumn = makeRSafe(dataXvalueColumn);
			this.rSafeYColumns = makeRSafe(dataYvalueColumns);
			this.valid = validate();
			this.fontsize = 4;
			this.plotVar = (plotVar == null ? PLOT_VAR : makeRSafe(plotVar));
			this.logTransformX = false;
			this.logTransformY = false;
			this.overWriteExisting = false;
			this.yMin = Double.NaN;
		}

		public void setOverWriteExisting(boolean overWriteExisting) {
			this.overWriteExisting = overWriteExisting;
		}

		public void setLogTransformX(boolean logTransformX) {
			this.logTransformX = logTransformX;
		}

		public void setLogTransformY(boolean logTransformY) {
			this.logTransformY = logTransformY;
		}

		public void setFontsize(int fontsize) {
			this.fontsize = fontsize;
		}

		public String getPlotVar() {
			return plotVar;
		}

		public void setOutput(String output) {
			this.output = output;
		}

		public String getOutput() {
			return output;
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

		public void setgPoint_SIZE(GEOM_POINT_SIZE gPoint_SIZE) {
			this.gPoint_SIZE = gPoint_SIZE;
		}

		@Override
		public boolean execute() {
			String[] rScript = developScript();
			log.report(Array.toStr(rScript, "\n"));
			Files.writeList(rScript, rScriptFile);
			boolean ran = CmdLine.runCommandWithFileChecks(new String[] { "Rscript", rScriptFile }, "", new String[] { rScriptFile }, new String[] { output }, true, overWriteExisting, false, log);
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
				String dataTable = DATA_TABLE_VAR + plotVar;
				String dataTableExtract = DATA_TABLE_EXTRACT + plotVar;
				String dataTableMelt = DATA_TABLE_MELT_VAR + plotVar;

				rCmd.add(dataTable + "=read.table(\"" + dataFile + "\", header=T, sep=\"\\t\")");

				String plot = "";
				// if (dataYvalueColumns.length > 1) {

				String[] toExtract = Array.concatAll(new String[] { rSafeXColumn }, rSafeYColumns);
				String extract = dataTableExtract + " <- " + dataTable + "[," + generateRVector(toExtract, true) + "]";
				rCmd.add(extract);
				String melt = dataTableMelt + "<-  melt(" + dataTableExtract + ",id.vars =\"" + rSafeXColumn + "\")";
				rCmd.add(melt);
				dataTable = dataTableMelt;
				plot += plotVar + " <- ggplot(" + dataTable + ",aes(x=" + rSafeXColumn + ", y=value, group=variable)) ";
				plot += " + " + sType.getCall() + "(aes(colour =variable)" + (gPoint_SIZE == null ? "" : "," + gPoint_SIZE.getCall()) + ")";

				// } else {
				// plot += plotVar + " <- ggplot(" + dataTable + ", aes(x = " + dataTable + "$" + rSafeXColumn + ")) ";
				// plot += " + " + sType.getCall() + "(aes(y = " + dataTable + "$" + rSafeYColumns[0] + "))";
				// }

				//
				// plot += PLOT_VAR + " <- ggplot(" + dataTable + ", aes(x = " + dataTable + "$" + rSafeXColumn + ")) ";
				// for (int i = 0; i < rSafeYColumns.length; i++) {
				// plot += " + " + sType.getCall() + "(aes(y = " + dataTable + "$" + rSafeYColumns[i] + "))";
				// }
				if (xLabel != null) {
					plot += " + xlab(\"" + xLabel + "\")";
				}
				if (yLabel != null) {
					plot += " + ylab(\"" + yLabel + "\")";
				}
				if (xRange != null) {
					plot += " + xlim(" + xRange[0] + "," + xRange[1] + ")";
				}
				if (yRange != null || !Double.isNaN(yMin)) {
					if (yRange != null) {
						plot += " + ylim(" + yRange[0] + "," + yRange[1] + ")";
					} else {
						plot += " + expand_limits(y=" + yMin + ")";
					}
				}
				if (logTransformX) {
					plot += "+ scale_x_log10()";
				}
				if (logTransformY) {
					plot += "+ scale_y_log10()";
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
				rCmd.add("ggsave(file=\"" + output + "\", width= 11, height = 6)");

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

					log.reportTimeError("Could not find all Y value columns in " + dataFile);
					int[] indices = ext.indexFactors(dataYvalueColumns, Files.getHeaderOfFile(dataFile, log), true, false);
					boolean[] extract = Array.booleanArray(indices.length, false);
					for (int i = 0; i < extract.length; i++) {
						if (indices[i] >= 0) {
							extract[i] = true;
						}
					}
					dataYvalueColumns = Array.subArray(dataYvalueColumns, extract);
					rSafeYColumns = makeRSafe(dataYvalueColumns);
					log.reportTimeWarning("Using available columns\t\t\n" + Array.toStr(dataYvalueColumns, "\n"));
					if (dataYvalueColumns.length < 1) {
						log.reportTimeError("Could not find any Y value columns, invalid plotting");
						valid = false;
					}
				}
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
			return valid;
		}

		public void setxLabel(String xLabel) {
			this.xLabel = xLabel;
		}

		public void setyLabel(String yLabel) {
			this.yLabel = yLabel;
		}

		public void setyMin(double yMin) {
			this.yMin = yMin;
		}

	}

	private static String getMultiPlotFunc() {
		String multiPlot = "";
		multiPlot += "multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {\n";
		multiPlot += "  require(grid)\n";

		multiPlot += " # Make a list from the ... arguments and plotlist\n";
		multiPlot += "  plots <- c(list(...), plotlist)\n";

		multiPlot += "  numPlots = length(plots)\n";

		multiPlot += "# If layout is NULL, then use 'cols' to determine layout\n";
		multiPlot += " if (is.null(layout)) {\n";
		multiPlot += "    # Make the panel\n";
		multiPlot += "   # ncol: Number of columns of plots\n";
		multiPlot += "# nrow: Number of rows needed, calculated from # of cols\n";
		multiPlot += " layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),\n";
		multiPlot += " ncol = cols, nrow = ceiling(numPlots/cols))\n";
		multiPlot += " }\n";

		multiPlot += " if (numPlots==1) {\n";
		multiPlot += "   print(plots[[1]])\n";

		multiPlot += " } else {\n";
		multiPlot += "   # Set up the page\n";
		multiPlot += "  grid.newpage()\n";
		multiPlot += "  pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))\n";

		multiPlot += "   # Make each plot, in the correct location\n";
		multiPlot += "  for (i in 1:numPlots) {\n";
		multiPlot += "    # Get the i,j matrix positions of the regions that contain this subplot\n";
		multiPlot += "    matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))\n";

		multiPlot += "    print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,\n";
		multiPlot += "layout.pos.col = matchidx$col))\n";
		multiPlot += "   }\n";
		multiPlot += " }\n";
		multiPlot += "}\n";

		return multiPlot;
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
