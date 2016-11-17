package org.genvisis.stats;

// import java.io.*;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class Rscript {
	public static final HashSet<String> R_INVALID_CHARS =
																											new HashSet<String>(Arrays.asList("-", "=",
																																												" ", "\\(",
																																												"\\)"));
	public static final String R_REPLACEMENT = ".";
	public static final String[] RSCRIPT_EXECS = {"/panfs/roc/itascasoft/R/3.2.2/bin/Rscript", // good
																																															// for
																																															// Mesabi
																								// "/soft/R/3.0.1/bin/Rscript", // MSI
																								"/soft/R/2.15.1/bin/Rscript", // Itasca nodes can
																																							// only see this
																								"/share/apps/R-3.0.1/bin/Rscript", // psych
																								"/share/apps/src/R-3.0.1/bin/Rscript", // alcatraz
	};

	public static final String[] R_EXECS = {"/panfs/roc/itascasoft/R/3.2.2/bin/R", // good for Mesabi
																					// "/soft/R/3.0.1/bin/R", // MSI
																					"/soft/R/2.15.1/bin/R", // Itasca nodes can only see this
																					"/share/apps/R-3.0.1/bin/R", // psych
																					"/share/apps/src/R-3.0.1/bin/R", // alcatraz
	};

	public static String getRscriptExecutable(Logger log) {
		if (Files.isWindows()) {
			return "Rscript";
		} else {
			for (String element : RSCRIPT_EXECS) {
				if (Files.exists(element)) {
					return element;
				}
			}
		}

		log.report("Warning - could not determine hostname, assuming R executbale is simply 'Rscript'");

		return "Rscript";
	}

	public static String getRExecutable(Logger log) {
		if (Files.isWindows()) {
			return "R";
		} else {
			for (String element : R_EXECS) {
				if (Files.exists(element)) {
					return element;
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
		for (String file : files) {
			root = ext.rootOf(dir + file, false);
			Files.qsub(root	+ ".qsub", getRscriptExecutable(log) + " --no-save " + root + ".R", 4000, 2,
									1);
			v.add("qsub " + root + ".qsub");
		}

		Files.writeArray(Array.toStringArray(v), dir + "master");
		Files.chmod(dir + "master");
	}

	/**
	 * @param toVector String[] that will be converted to c("toVector[0]",toVector[1]...);
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
		try {
			Double.parseDouble(rSafe.charAt(0) + "");
			rSafe = "X" + rSafe;
		} catch (NumberFormatException nfe) {

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
														LINE("geom_line"), POINT("geom_point"), BOX("geom_boxplot"),
														/**
														 *
														 * basic x v y
														 */
														BASIC_POINT("geom_point"), NO_MELT_POINT("geom_point"), BOX_NO_MELT("geom_boxplot"), HIST("geom_histogram");

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
																BOTTOM("legend.position=\"bottom\""), NONE("legend.position=\"none\"");

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

		private final String call;
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
																	COLUMNS_MULTIPLOT_1("cols =", 1), COLUMNS_MULTIPLOT_2("cols =",
																																												2);

		private final String call;
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

		private final RScatter[] rScatters;
		private final String mergeOutput;
		private final String rScriptFile;
		private final COLUMNS_MULTIPLOT cMultiplot;
		private final PLOT_DEVICE device;
		private final Logger log;

		public RScatters(	ArrayList<RScatter> rsScatters, String rScriptFile, String mergeOutput,
											COLUMNS_MULTIPLOT cMultiplot, PLOT_DEVICE device, Logger log) {
			this(	rsScatters.toArray(new RScatter[rsScatters.size()]), rScriptFile, mergeOutput,
						cMultiplot, device, log);
		}

		public RScatters(	RScatter[] rScatter, String rScriptFile, String mergeOutput,
											COLUMNS_MULTIPLOT cMultiplot, PLOT_DEVICE device, Logger log) {
			super();
			rScatters = rScatter;
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
			// log.report(Array.toStr(rScript, "\n"));
			Files.writeArray(rScript, rScriptFile);
			boolean ran = CmdLine.runCommandWithFileChecks(
																											new String[] {rScatters[0].getrScriptLoc(),
																																		rScriptFile},
																											"", new String[] {rScriptFile},
																											new String[] {mergeOutput}, true, true, false,
																											log);
			return ran;
		}

		@Override
		public String[] developScript() {
			if (validate()) {
				ArrayList<String> allPlots = new ArrayList<String>();
				allPlots.add(getMultiPlotFunc());

				String plotDevice = device.getCall() + "\"" + mergeOutput + "\" ,onefile = TRUE)";
				allPlots.add(plotDevice);
				for (RScatter rScatter : rScatters) {
					String[] tmpScript = rScatter.developScript();
					if (tmpScript != null) {
						for (String element : tmpScript) {
							allPlots.add(element);
						}
					}

				}
				for (RScatter rScatter : rScatters) {

					allPlots.add(rScatter.getPlotVar());

				}
				allPlots.add("dev.off()");
				return allPlots.toArray(new String[allPlots.size()]);
			} else {
				return null;
			}
		}

		@Override
		public boolean validate() {
			boolean valid = true;
			Hashtable<String, String> plotVars = new Hashtable<String, String>();
			for (RScatter rScatter : rScatters) {
				if (plotVars.containsKey(rScatter.getPlotVar())) {
					log.reportError("Multiple plot vars " + rScatter.getPlotVar());
					valid = false;
				} else {
					plotVars.put(rScatter.getPlotVar(), rScatter.getPlotVar());
				}
			}
			return valid;
		}
	}

	/**
	 * plot error bars on data points
	 *
	 */
	public static class ErrorBars {

		private String[] rSafeDataColumns;
		private final String[] rSafeErrorColumns;
		private boolean doubleCoded;

		public ErrorBars(String[] dataColumns, String[] stdErrorColumns) {
			super();
			rSafeDataColumns = makeRSafe(dataColumns);
			rSafeErrorColumns = makeRSafe(stdErrorColumns);
			doubleCoded = false;
		}

		public void setDoubleCoded(boolean doubleCoded) {
			this.doubleCoded = doubleCoded;
		}

		public String[] getErrorBarCalls(String dataFrame, String xColumn) {
			String[] errorBars = new String[rSafeDataColumns.length];
			if (doubleCoded) {
				int dataIndex = 0;
				for (int i = 0; i < rSafeErrorColumns.length; i += 2) {
					StringBuilder builder = new StringBuilder();
					builder.append("geom_errorbar(data=" + dataFrame + ",");
					builder.append("aes(x=" + xColumn + ", ymin =" + rSafeErrorColumns[i] + ",");
					builder.append("ymax="	+ rSafeErrorColumns[i + 1] + ",colour=\""
													+ rSafeDataColumns[dataIndex] + "\"))");
					errorBars[dataIndex] = builder.toString();
					dataIndex++;

				}
			} else {

				for (int i = 0; i < errorBars.length; i++) {
					StringBuilder builder = new StringBuilder();
					builder.append("geom_errorbar(data=" + dataFrame + ",");
					builder.append("aes(x="	+ xColumn + ", ymin =" + rSafeDataColumns[i] + "-"
													+ rSafeErrorColumns[i] + ",");
					builder.append("ymax="	+ rSafeDataColumns[i] + "+" + rSafeErrorColumns[i] + ",colour=\""
													+ rSafeDataColumns[i] + "\"))");
					errorBars[i] = builder.toString();
				}
			}
			return errorBars;

		}

		public String[] getrSafeErrorColumns() {
			return rSafeErrorColumns;
		}

		public void setrSafeDataColumns(String[] rSafeDataColumns) {
			this.rSafeDataColumns = rSafeDataColumns;
		}

	}

	public static class GeomText {
		private static final String[] HEADER = new String[] {"X", "Y", "ANGLE", "LABEL", "FONTSIZE"};
		private final double x;
		private double y;
		private final double angle;
		private final String label;
		private String color;

		private final double fontSize;

		public GeomText(double x, double y, double angle, String label, double fontSize) {
			super();
			this.x = x;
			this.y = y;
			this.angle = angle;
			this.label = label;
			this.fontSize = fontSize;
			color = null;
		}

		public double getY() {
			return y;
		}

		public String getColor() {
			return color;
		}

		public void setColor(String color) {
			this.color = color;
		}

		public double getX() {
			return x;
		}

		public String getLabel() {
			return label;
		}

		public void setY(double y) {
			this.y = y;
		}

		public String getCommand() {
			StringBuilder stringBuilder = new StringBuilder();
			stringBuilder.append(" geom_text(data = NULL,aes(x =");
			stringBuilder.append(x + ", y = ");
			stringBuilder.append(y + ", label=");
			stringBuilder.append("\"" + label + "\",angle=" + angle + ",size = " + fontSize + " ))");
			// ,size=" + fontSize + "
			// + (color == null ? "" : ",color=\"" + color + "\"")
			return stringBuilder.toString();
		}

		public static GeomText[] fromFile(String file, Logger log) {
			if (!Files.exists(file)) {
				log.reportFileNotFound(file);
				return null;
			} else if (!Files.headerOfFileContainsAll(file, HEADER, log)) {
				log.reportError("Geom text " + file + "  must have header " + Array.toStr(HEADER));
				return null;
			} else {
				ArrayList<GeomText> geomTexts = new ArrayList<GeomText>();
				try {
					int[] indices = ext.indexFactors(HEADER, Files.getHeaderOfFile(file, log), true, false);
					BufferedReader reader = Files.getAppropriateReader(file);
					reader.readLine();
					while (reader.ready()) {
						String[] line = reader.readLine().trim().split("\t");
						double x = Double.parseDouble(line[indices[0]]);
						double y = Double.parseDouble(line[indices[1]]);
						double angle = Double.parseDouble(line[indices[2]]);
						String label = line[indices[3]];
						int fontSize = Integer.parseInt(line[indices[4]]);
						geomTexts.add(new GeomText(x, y, angle, label, fontSize));

					}
					reader.close();
				} catch (FileNotFoundException fnfe) {
					log.reportError("Error: file \"" + file + "\" not found in current directory");
					return null;
				} catch (IOException ioe) {
					log.reportError("Error reading file \"" + file + "\"");
					return null;
				}
				return geomTexts.toArray(new GeomText[geomTexts.size()]);
			}

		}

	}

	// public static class RBoxPlot implements RCommand{
	//
	// }

	public static class SeriesLabeler {
		private boolean labelMax;
		private boolean labelMin;
		private String indexTag;
		private String valueTag;
		private int textSize;
		private double textAngle;

		public SeriesLabeler(boolean labelMax, boolean labelMin, String indexTag, String valueTag) {
			super();
			this.labelMax = labelMax;
			this.labelMin = labelMin;
			this.indexTag = indexTag;
			this.valueTag = valueTag;
			textSize = 3;
			textAngle = 0;

		}

		public void setLabelMax(boolean labelMax) {
			this.labelMax = labelMax;
		}

		public void setLabelMin(boolean labelMin) {
			this.labelMin = labelMin;
		}

		public void setIndexTag(String indexTag) {
			this.indexTag = indexTag;
		}

		public void setValueTag(String valueTag) {
			this.valueTag = valueTag;
		}

		public void setTextSize(int textSize) {
			this.textSize = textSize;
		}

		public void setTextAngle(double textAngle) {
			this.textAngle = textAngle;
		}

		public int getTextSize() {
			return textSize;
		}

		public double getTextAngle() {
			return textAngle;
		}

		public boolean isLabelMax() {
			return labelMax;
		}

		public boolean isLabelMin() {
			return labelMin;
		}

		public String getIndexTag() {
			return indexTag;
		}

		public String getValueTag() {
			return valueTag;
		}

	}

	public static class SlopeLine {
		private final double xIntercept;
		private final double slope;

		public SlopeLine(double xIntercept, double slope) {
			super();
			this.xIntercept = xIntercept;
			this.slope = slope;
		}

		public String getCall() {
			return "geom_abline(xintercept =" + xIntercept + ",slope=" + slope + ")";
		}

	}

	public static class VertLine {
		private final double xIntercept;

		public VertLine(double xIntercept) {
			super();
			this.xIntercept = xIntercept;
		}

		public String getCall() {
			return "geom_vline(xintercept =" + xIntercept + ")";
		}
	}

	public static class HorizontalLine {
		private final double yIntercept;

		public HorizontalLine(double yIntercept) {
			super();
			this.yIntercept = yIntercept;
		}

		public String getCall() {
			return "geom_hline(yintercept =" + yIntercept + ")";
		}
	}

	public static class Restrictions {
		private final String[] cols;
		private final String[] vals;
		private final String[] ops;
		private final String cond;

		public Restrictions(String cols[], double[] vals, String[] ops, String cond) {
			super();
			this.cols = cols;
			this.vals = Array.toStringArray(vals);
			this.ops = ops;
			this.cond = cond;
		}

		public Restrictions(String cols[], String[] vals, String[] ops, String cond) {
			super();
			this.cols = cols;
			this.vals = vals;
			this.ops = ops;
			this.cond = cond;
		}

		private String getRestriction(String dataFrame) {
			String minStr = "";
			for (int i = 0; i < ops.length; i++) {
				minStr += (i > 0 ? cond : "") + dataFrame + "$" + cols[i] + ops[i] + vals[i];
			}
			String sub = "subset(" + dataFrame + "," + minStr + ")";
			return sub;
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
		private SlopeLine[] slopeLines;
		private VertLine[] vertLines;
		private HorizontalLine[] horizontalLines;

		private final String dataFile;
		private String output;
		private final String rScriptFile;
		private String[] dataYvalueColumns, rSafeYColumns, rSafeAltYColumnNames;
		private final String dataXvalueColumn, rSafeXColumn, colorColumn;
		private String rSafeColorColumn;
		private double[] xRange, yRange;
		private double yMin;
		private String xLabel;
		private String yLabel;
		private final Logger log;
		private final SCATTER_TYPE sType;
		private LEGEND_POSITION lPosition;
		private GEOM_POINT_SIZE gPoint_SIZE;
		private boolean valid;
		private double fontsize;
		private String title;
		private String plotVar;
		private boolean logTransformX;
		private boolean logTransformY;
		private boolean overWriteExisting;
		private int height;
		private final int numHistBins;
		private int width;
		private GeomText[] gTexts;
		private SeriesLabeler seriesLabeler;
		private boolean directLableGtexts, onlyMaxMin;
		private ErrorBars errorBars;
		private boolean factorColor;
		// private String legendName;
		private boolean regLines;
		private boolean scaleDensity;
		private double midScaleDensity;
		private Restrictions[] restrictions;
		private String rScriptLoc;
		private String[][] altLegendTitles;
		private String legendTitle;

		public RScatter(String dataFile, String rSriptFile, String plotVar, String output,
										String dataXvalueColumn, String[] dataYvalueColumns, SCATTER_TYPE sType,
										Logger log) {
			this(	dataFile, rSriptFile, plotVar, output, dataXvalueColumn, dataYvalueColumns, null, sType,
						log);
		}

		public void setSlopeLines(SlopeLine[] slopeLines) {
			this.slopeLines = slopeLines;
		}

		/**
		 * @param dataFile full path to the input data file, where all data to plot will come from
		 * @param rSriptFile full path to the rscript file (which will be generated after calling
		 *        {@link RScatter#execute()}
		 * @param plotVar the variable name that holds the R plot variable, which must be unique per R
		 *        instance
		 * @param output full path to the output plot (extension determines type, i.e .pdf produces
		 *        pdfs, .jpg produces jpegs
		 * @param dataXvalueColumn the column name holding the x-axis data
		 * @param dataYvalueColumns the column names holding the y-axis data
		 * @param colorColumn if not null, points will be colored by this column
		 * @param sType Type of plot eg, {@link SCATTER_TYPE#POINT} or {@link SCATTER_TYPE#BOX}
		 * @param log
		 */
		public RScatter(String dataFile, String rSriptFile, String plotVar, String output,
										String dataXvalueColumn, String[] dataYvalueColumns, String colorColumn,
										SCATTER_TYPE sType, Logger log) {
			super();
			this.dataFile = dataFile;
			rScriptFile = rSriptFile;
			this.dataXvalueColumn = dataXvalueColumn;
			this.dataYvalueColumns = dataYvalueColumns;
			this.log = log;
			this.output = ext.rootOf(output, false) + ".jpeg";
			this.sType = sType;
			rSafeXColumn = makeRSafe(dataXvalueColumn);
			rSafeYColumns = makeRSafe(dataYvalueColumns);
			this.colorColumn = colorColumn;
			if (colorColumn != null) {
				rSafeColorColumn = makeRSafe(colorColumn);
			}
			valid = validate();
			fontsize = 4;
			this.plotVar = (plotVar == null ? PLOT_VAR : makeRSafe(plotVar));
			logTransformX = false;
			logTransformY = false;
			overWriteExisting = false;
			yMin = Double.NaN;
			height = 6;
			width = 11;
			directLableGtexts = false;
			onlyMaxMin = false;
			// this.legendName = "variable";
			regLines = false;
			midScaleDensity = .5;
			rScriptLoc = "Rscript";
			factorColor = true;
			numHistBins = -1;
			altLegendTitles = null;
			legendTitle = null;
		}

		public String getrScriptLoc() {
			return rScriptLoc;
		}

		public void setLegendTitle(String legendTitle) {
			this.legendTitle = legendTitle;
		}

		public void setAltLegendTitles(String[][] altLegendTitles) {
			this.altLegendTitles = altLegendTitles;
		}

		public void setFactorColor(boolean factorColor) {
			this.factorColor = factorColor;
		}

		public void setrScriptLoc(String rScriptLoc) {
			this.rScriptLoc = rScriptLoc;
		}


		// private void setNumHistBins(int numHistBins) {
		// this.numHistBins = numHistBins;
		// }

		public ErrorBars getErrorBars() {
			return errorBars;
		}

		public void setVertLines(VertLine[] vertLines) {
			this.vertLines = vertLines;
		}

		public void setHorizontalLines(HorizontalLine[] horizontalLines) {
			this.horizontalLines = horizontalLines;
		}

		public void setMidScaleDensity(double midScaleDensity) {
			this.midScaleDensity = midScaleDensity;
		}

		public void setScaleDensity(boolean scaleDensity) {
			this.scaleDensity = scaleDensity;
		}

		public void setRestrictions(Restrictions[] restrictions) {
			this.restrictions = restrictions;
		}

		public void setRegLines(boolean regLines) {
			this.regLines = regLines;
		}

		public String[] getrSafeYColumns() {
			return rSafeYColumns;
		}

		public void setErrorBars(ErrorBars errorBars) {
			this.errorBars = errorBars;
		}

		public String[] getrSafeAltYColumnNames() {
			return rSafeAltYColumnNames;
		}

		public void setrSafeAltYColumnNames(String[] rSafeAltYColumnNames) {
			if (rSafeAltYColumnNames.length != rSafeYColumns.length) {
				valid = false;
				log.reportError("alternate y -column names must be the same length as available data");
				rSafeYColumns = null;

			}
			this.rSafeAltYColumnNames = rSafeAltYColumnNames;
		}

		public GeomText[] getgTexts() {
			return gTexts;
		}

		public void setOnlyMaxMin(boolean onlyMaxMin) {
			this.onlyMaxMin = onlyMaxMin;
		}

		public SeriesLabeler getSeriesLabeler() {
			return seriesLabeler;
		}

		public void setDirectLableGtexts(boolean directLableGtexts) {
			this.directLableGtexts = directLableGtexts;
		}

		public void setPlotVar(String plotVar) {
			this.plotVar = plotVar;
		}

		public void setSeriesLabeler(SeriesLabeler seriesLabeler) {
			this.seriesLabeler = seriesLabeler;
		}

		public void setgTexts(GeomText[] gTexts) {
			this.gTexts = gTexts;
		}

		public String getTitle() {
			return title;
		}

		public String getxLabel() {
			return xLabel;
		}

		public String getyLabel() {
			return yLabel;
		}

		public String getDataFile() {
			return dataFile;
		}

		public void setHeight(int height) {
			this.height = height;
		}

		public void setWidth(int width) {
			this.width = width;
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

		public void setFontsize(double fontsize) {
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
			if (rScript != null) {
				// log.report(Array.toStr(rScript, "\n"));
				Files.writeArray(rScript, rScriptFile);
				boolean ran = CmdLine.runCommandWithFileChecks(	new String[] {rScriptLoc, rScriptFile}, "",
																												new String[] {rScriptFile},
																												new String[] {output}, true,
																												overWriteExisting, false, log);
				return ran;
			} else {
				return false;
			}
		}

		private GeomText[] getLabelers() {
			if (seriesLabeler != null) {
				ArrayList<GeomText> geomTexts = new ArrayList<GeomText>();
				double[] xs = Array.toDoubleArray(HashVec.loadFileToStringArray(dataFile, true,
																																				new int[] {ext.indexOfStr(dataXvalueColumn,
																																																	Files.getHeaderOfFile(dataFile,
																																																												log))},
																																				false));
				for (String dataYvalueColumn : dataYvalueColumns) {
					int maxIndex = -1;
					int minIndex = -1;
					double min = Double.MAX_VALUE;
					double max = -1 * Double.MAX_VALUE;

					double[] data = Array.toDoubleArray(HashVec.loadFileToStringArray(dataFile, true,
																																						new int[] {ext.indexOfStr(dataYvalueColumn,
																																																			Files.getHeaderOfFile(dataFile,
																																																														log))},
																																						false));
					for (int j = 0; j < data.length; j++) {
						if (data[j] > max) {
							max = data[j];
							maxIndex = j;
						}
						if (data[j] < min) {
							min = data[j];
							minIndex = j;

						}
					}

					if (seriesLabeler.isLabelMax()) {
						String maxLabel = "MAX_"	+ dataYvalueColumn + "_" + seriesLabeler.getIndexTag() + "-"
															+ Math.round(xs[maxIndex]) + " " + seriesLabeler.getValueTag() + "-"
															+ ext.formDeci(data[maxIndex], 2);
						GeomText maxText = new GeomText(xs[maxIndex], max, seriesLabeler.getTextAngle(),
																						maxLabel, seriesLabeler.getTextSize());
						// maxText.setColor(dataYvalueColumns[i]);
						geomTexts.add(maxText);
					}
					if (seriesLabeler.isLabelMin()) {

						String minLabel = "MIN_"	+ dataYvalueColumn + "_" + seriesLabeler.getIndexTag() + "-"
															+ Math.round(xs[minIndex]) + " " + seriesLabeler.getValueTag() + "-"
															+ ext.formDeci(data[minIndex], 2);
						GeomText minText = new GeomText(xs[minIndex], min, seriesLabeler.getTextAngle(),
																						minLabel, seriesLabeler.getTextSize());
						// minText.setColor(dataYvalueColumns[i]);
						geomTexts.add(minText);
					}
				}

				return geomTexts.toArray(new GeomText[geomTexts.size()]);

			}
			return null;
		}

		private String[] directLabelFrame(String currentPlot, String factorTitel, String xTitle,
																			String meltYTitle) {
			String[] directDataFrame = new String[2];
			String[] types = new String[gTexts.length];
			double[] xs = new double[gTexts.length];
			double[] ys = new double[gTexts.length];
			for (int i = 0; i < gTexts.length; i++) {
				types[i] = gTexts[i].getLabel();
				xs[i] = gTexts[i].getX();
				ys[i] = gTexts[i].getY();
			}
			String typeVector = generateRVector(types, true);
			String xVector = generateRVector(Array.toStringArray(xs), false);
			String yVector = generateRVector(Array.toStringArray(ys), false);

			String df = currentPlot + "_direcLdf";
			directDataFrame[0] = df;
			String dfGenerate = df	+ " <-  data.frame(" + xTitle + "=" + xVector + "," + factorTitel + "="
													+ typeVector + "," + meltYTitle + "=" + yVector + ")";
			directDataFrame[1] = dfGenerate;
			return directDataFrame;
		}

		private static String getRename(String[] original, String[] alts, String dataframe) {
			String[] toVec = new String[original.length];
			for (int i = 0; i < toVec.length; i++) {
				toVec[i] = "\"" + original[i] + "\"=\"" + alts[i] + "\"";
			}

			String rename = dataframe	+ " <- rename(" + dataframe + "," + generateRVector(toVec, false)
											+ ")";
			return rename;
		}

		@Override
		public String[] developScript() {

			if (valid) {
				if (seriesLabeler != null) {
					GeomText[] labels = getLabelers();
					if (gTexts != null) {
						gTexts = Array.concatAll(gTexts, labels);
					} else {
						gTexts = labels;
					}
				}
				ArrayList<String> rCmd = new ArrayList<String>();

				rCmd.add("library(scales)");
				rCmd.add("library(ggplot2)");
				rCmd.add("require(reshape2)");
				rCmd.add("library(RColorBrewer)");
				if (directLableGtexts) {
					rCmd.add(" library(directlabels)");
				}
				rCmd.add("library(plyr)");
				String dataTable = DATA_TABLE_VAR + plotVar;
				String dataTableExtract = DATA_TABLE_EXTRACT + plotVar;
				String dataTableMelt = DATA_TABLE_MELT_VAR + plotVar;

				rCmd.add(dataTable + "=read.table(\"" + dataFile + "\", header=T, sep=\"\\t\")");
				if (restrictions != null) {
					for (Restrictions restriction : restrictions) {
						rCmd.add(dataTable + " <- " + restriction.getRestriction(dataTable));
					}
				}
				String plot = "";
				// if (dataYvalueColumns.length > 1) {
				if (rSafeAltYColumnNames != null
						&& !Array.equals(rSafeYColumns, rSafeAltYColumnNames, true)) {
					String[] currentNames = rSafeYColumns;
					String[] altNames = rSafeAltYColumnNames;
					if (errorBars != null) {
						currentNames = Array.concatAll(currentNames, errorBars.getrSafeErrorColumns());
						altNames = Array.concatAll(altNames, errorBars.getrSafeErrorColumns());
					}
					String rename = getRename(currentNames, altNames, dataTable);
					rCmd.add(rename);
				} else {
					rSafeAltYColumnNames = rSafeYColumns;
				}

				String[] toExtract = Array.concatAll(new String[] {rSafeXColumn}, rSafeAltYColumnNames);
				if (errorBars != null) {
					toExtract = Array.concatAll(toExtract, errorBars.getrSafeErrorColumns());
				}
				if (colorColumn != null) {
					toExtract = Array.concatAll(toExtract, new String[] {rSafeColorColumn});
				}
				String extract = dataTableExtract	+ " <- " + dataTable + "[,"
													+ generateRVector(toExtract, true) + "]";

				rCmd.add(extract);

				String melt = dataTableMelt	+ "<-  melt(" + dataTableExtract + ",id.vars =\"" + rSafeXColumn
											+ "\")";
				rCmd.add(melt);
				int numColors = rSafeYColumns.length;
				if (errorBars != null) {
					numColors += errorBars.getrSafeErrorColumns().length;
				}
				if (colorColumn != null) {
					numColors++;
				}
				StringBuilder pallete = new StringBuilder();
				pallete.append("getPalette = colorRampPalette(brewer.pal(9, \"Set1\"))\n");
				pallete.append("myColors <- getPalette(" + numColors + ")\n");
				if (sType == SCATTER_TYPE.BASIC_POINT) {
					pallete.append("names(myColors)  <- levels(" + dataTableExtract + ")\n");
				} else {
					pallete.append("names(myColors)  <- levels(" + dataTableMelt + "$variable)\n");
				}
				pallete.append("colScale <- scale_colour_manual(name = \"Var\",values = myColors)\n");
				String colScale = pallete.toString();
				if (rSafeColorColumn == null) {
					rCmd.add(colScale);
				}
				// names(myColors) <-
				// levels(dataTableMeltcorrectionEval_WITHOUT_INDEPS_STEPWISE_RANK_R2_WITH_BUILDERS.summary.DK.summary.heritability_summary.parsed$variable)
				// /colScale <- scale_colour_manual(name = "grp",values = myColors)
				dataTable = dataTableMelt;
				if (rSafeColorColumn != null || sType == SCATTER_TYPE.NO_MELT_POINT) {
					plot += plotVar	+ " <- ggplot(" + dataTableExtract + ",aes(x=" + rSafeXColumn
									+ (scaleDensity ? ",color=" + rSafeColorColumn : "") + "))";
					for (String rSafeYColumn : rSafeYColumns) {
						System.out.println(plot);
						if (!scaleDensity) {
							plot += " + "	+ sType.getCall() + "(aes(y=" + rSafeYColumn + ", color="
											+ (factorColor ? "factor(" + rSafeColorColumn + ")" : rSafeColorColumn)
											+ "))";
							if (altLegendTitles != null) {
								plot += " + scale_color_discrete(\""
												+ (legendTitle == null ? "metric" : legendTitle) + "\",";
								plot += "breaks=" + generateRVector(altLegendTitles[0], true) + ",";
								plot += "labels=" + generateRVector(altLegendTitles[1], true) + ")";
								System.out.println(Array.toStr(altLegendTitles[0]));
								System.out.println(Array.toStr(altLegendTitles[1]));
							}

						} else {
							plot += " + " + sType.getCall() + "(aes(y=" + rSafeYColumn + "))";
						}
						if (scaleDensity) {
							System.out.println("DFSSDF");
							System.out.println();
							plot +=
										"+ scale_colour_gradient2(low=\"blue\", high=\"red\",mid=\"yellow\" , midpoint="
											+ midScaleDensity + ", " + "\"" + rSafeColorColumn + "\"" + ")";
						}
					}
				} else if (sType == SCATTER_TYPE.BOX) {

					// String labs = " xlabs <- paste(factor(levels(" + dataTableExtract + "$" + rSafeXColumn
					// + "))," + "\"\\n(N=\",table(" + dataTableExtract + "$" + rSafeXColumn +
					// "),\")\",sep=\"\")";
					// rCmd.add(labs);
					String addN = dataTable	+ " <- ddply(" + dataTable + ",.(" + rSafeXColumn
												+ ",variable),transform,N=length(" + rSafeXColumn + "))";
					rCmd.add(addN);
					String label = dataTable	+ "$label <- paste0(" + dataTable + "$" + rSafeXColumn
													+ ",\"\\n\",\"(n=\"," + dataTable + "$N,\")\")";
					rCmd.add(label);
					String order = dataTable	+ " <- " + dataTable + "[ order(" + dataTable + "$"
													+ rSafeXColumn + "), ]";
					// String sort = dataTableMelt + "$label <- ";
					// System.out.println(sort);
					rCmd.add(order);
					plot += plotVar + " <- ggplot(" + dataTable + ",aes(x=variable, y=value)) ";
					plot += " + " + sType.getCall() + "(aes(fill=factor(" + rSafeXColumn + ")))";

					// plot += " + scale_x_discrete(limits=" + dataTable + "$label)";
					// factor(" + dataTableMelt + "$label, levels = " + dataTableMelt +
					// "$label[order(as.numeric(" + dataTableMelt + "$variable))]))+facet_grid(.~variable)";

				} else if (sType == SCATTER_TYPE.BOX_NO_MELT) {
					plot += plotVar + " <- ggplot(" + dataTable + ",aes(x=variable, y=value)) ";
					plot += " + " + sType.getCall() + "(aes(fill=factor(" + rSafeXColumn + ")))";

				} else if (sType == SCATTER_TYPE.BASIC_POINT) {
					if (rSafeYColumns.length > 1) {
						log.reportError(SCATTER_TYPE.BASIC_POINT
																+ " demands a single x and single y value");
						return null;
					}
					plot += plotVar	+ " <- ggplot(" + dataTableExtract + ",aes(x=" + rSafeXColumn + ", y="
									+ rSafeYColumns[0] + ")) ";
					plot += " + " + sType.getCall() + "(aes(fill=" + rSafeXColumn + "))";
				} else if (sType == SCATTER_TYPE.HIST) {
					plot += plotVar + " <- ggplot(" + dataTableExtract + ",aes(x=" + rSafeYColumns[0] + "))";
					plot += " + "	+ sType.getCall() + "()+stat_bin("
									+ (numHistBins > 0 ? "bins = " + numHistBins + "," : "")
									+ "aes(y=..count.., label=format(round(100*(..count..)/sum(..count..), 0), nsmall = 0)), geom=\"text\", vjust=-.5) ";
				}

				else {
					if (directLableGtexts && gTexts != null) {
						String[] directFrame = directLabelFrame(plotVar, "variable", rSafeXColumn, "value");
						rCmd.add(directFrame[1]);
						String aes = "data="	+ directFrame[0] + " , aes(colour =variable,x=" + rSafeXColumn
													+ ", y=value, group=variable)";
						plot += plotVar	+ " <- direct.label( ggplot(" + aes + ") + " + sType.getCall() + "("
										+ aes + "),  list(\"smart.grid\", cex=" + fontsize + "))";
					} else {
						plot += plotVar + " <- ggplot()";

					}
					// , col=variable
					if (!onlyMaxMin) {
						// data=subset(dataTableMeltcorrectionEval_WITHOUT_INDEPS_STEPWISE_RANK_R2_WITH_BUILDERS.summary.DK.summary.heritability_summary.parsed,variable==c("SOLAR_PERCENT_HERITABLITY","MERLIN_PERCENT_HERITABLITY")
						String data = "";
						if (errorBars != null) {
							data = dataTable	+ "[" + dataTable + "$variable %in% "
											+ generateRVector(rSafeAltYColumnNames, true) + ",]";
							// data = "subset(" + dataTable + ",variable==" + + ")";
						} else {
							data = dataTable;
						}
						plot += " + "	+ sType.getCall() + "(data=" + data + " , aes(colour =variable,x="
										+ rSafeXColumn + ", y=value, group=variable)"
										+ (gPoint_SIZE == null ? "" : "," + gPoint_SIZE.getCall()) + ")";

					}
					System.out.println(plot);
				}

				if (errorBars != null) {
					String[] errors = errorBars.getErrorBarCalls(dataTableExtract, rSafeXColumn);
					for (String error : errors) {
						plot += " + " + error;
					}
				}

				// } else {
				// plot += plotVar + " <- ggplot(" + dataTable + ", aes(x = " + dataTable + "$" +
				// rSafeXColumn + ")) ";
				// plot += " + " + sType.getCall() + "(aes(y = " + dataTable + "$" + rSafeYColumns[0] +
				// "))";
				// }

				//
				// plot += PLOT_VAR + " <- ggplot(" + dataTable + ", aes(x = " + dataTable + "$" +
				// rSafeXColumn + ")) ";
				// for (int i = 0; i < rSafeYColumns.length; i++) {
				// plot += " + " + sType.getCall() + "(aes(y = " + dataTable + "$" + rSafeYColumns[i] +
				// "))";
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
				if (sType != SCATTER_TYPE.BOX && sType != SCATTER_TYPE.HIST && rSafeColorColumn == null) {
					plot += " + colScale";
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
				if (gTexts != null && !directLableGtexts) {
					for (GeomText gText : gTexts) {
						plot += " + " + gText.getCommand();
					}
				}
				if (regLines) {
					plot += " + geom_smooth(method=lm,  se=FALSE)";
				}
				if (vertLines != null) {
					for (VertLine vertLine : vertLines) {
						plot += " + " + vertLine.getCall();
					}
				}
				if (slopeLines != null) {
					for (SlopeLine slopeLine : slopeLines) {
						plot += " + " + slopeLine.getCall();
					}
				}
				if (horizontalLines != null) {
					for (HorizontalLine horizontalLine : horizontalLines) {
						plot += " + " + horizontalLine.getCall();

					}
				}

				if (sType == SCATTER_TYPE.BOX) {

				}
				rCmd.add(plot);

				rCmd.add("ggsave(file=\""	+ output + "\", width= " + width + ", height = " + height
									+ ",limitsize=FALSE)");
				// System.out.println(output);
				// System.out.println(Array.toStr(Array.toStringArray(rCmd),"\n"));
				//
				//
				// System.out.println(plot);
				// System.exit(1);
				// System.exit(1);
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

				if (!Files.headerOfFileContainsAll(dataFile, dataYvalueColumns, false, log)) {

					log.reportError("Could not find all Y value columns in " + dataFile);
					int[] indices = ext.indexFactors(	dataYvalueColumns, Files.getHeaderOfFile(dataFile, log),
																						true, false);
					boolean[] extract = Array.booleanArray(indices.length, false);
					for (int i = 0; i < extract.length; i++) {
						if (indices[i] >= 0) {
							extract[i] = true;
						}
					}
					dataYvalueColumns = Array.subArray(dataYvalueColumns, extract);
					rSafeYColumns = makeRSafe(dataYvalueColumns);

					log.reportTimeWarning("Using available columns\t\t\n"
																+ Array.toStr(dataYvalueColumns, "\n"));
					if (dataYvalueColumns.length < 1) {
						log.reportError("Could not find any Y value columns, invalid plotting");
						valid = false;
					}
				}
				if (!Files.headerOfFileContainsAll(dataFile, new String[] {dataXvalueColumn}, log)) {
					valid = false;
					log.reportError("Could not find all X value column in " + dataFile);
				} else {
					if (xRange != null && xRange.length != 2) {
						valid = false;
						log.reportError("x range must be length 2");
					} else {
						if (yRange != null && yRange.length != 2) {
							valid = false;
							log.reportError("y range must be length 2");
						}
					}
				}
				if (colorColumn != null
						&& !Files.headerOfFileContainsAll(dataFile, new String[] {colorColumn}, log)) {
					log.reportError("Could not find color column in " + dataFile);
					valid = false;
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

		String usage = "\n"	+ "stats.Rscript requires 0-1 arguments\n"
										+ "   (1) create qsub files for all .R files in the directory (i.e. -batchUp (not the default))\n"
										+ "   (2) directory to search for the .R files (i.e. dir=" + dir
										+ " (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("-batchUp")) {
				batchUp = true;
				numArgs--;
			} else if (arg.startsWith("dir=")) {
				dir = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("log=")) {
				logfile = arg.split("=")[1];
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
			log = new Logger(logfile);
			if (batchUp) {
				batchUp(dir, log);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
