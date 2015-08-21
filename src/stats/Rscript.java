package stats;

//import java.io.*;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;

import common.*;

public class Rscript {
	public static final HashSet<String> R_INVALID_CHARS = new HashSet<String>(Arrays.asList("-","="));
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
		LINE("geom_line"), POINT("geom_point"), BOX("geom_boxplot");

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
			// log.report(Array.toStr(rScript, "\n"));
			Files.writeList(rScript, rScriptFile);
			boolean ran = CmdLine.runCommandWithFileChecks(new String[] { "Rscript", rScriptFile }, "", new String[] { rScriptFile }, new String[] { mergeOutput }, true, true, false, log);
			return ran;
		}

		@Override
		public String[] developScript() {
			if (validate()) {
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
			} else {
				return null;
			}
		}

		@Override
		public boolean validate() {
			boolean valid = true;
			Hashtable<String, String> plotVars = new Hashtable<String, String>();
			for (int i = 0; i < rScatters.length; i++) {
				if (plotVars.containsKey(rScatters[i].getPlotVar())) {
					log.reportTimeError("Multiple plot vars " + rScatters[i].getPlotVar());
					valid = false;
				} else {
					plotVars.put(rScatters[i].getPlotVar(), rScatters[i].getPlotVar());
				}
			}
			return valid;
		}
	}

	
	
	public static class GeomText {
		private static final String[] HEADER = new String[] { "X", "Y", "ANGLE", "LABEL","FONTSIZE" };
		private double x;
		private double y;
		private double angle;
		private String label;
		private String color;
		private int fontSize;

		public GeomText(double x, double y,double angle, String label, int fontSize) {
			super();
			this.x = x;
			this.y = y;
			this.angle = angle;
			this.label = label;
			this.fontSize = fontSize;
			this.color =null;
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
			stringBuilder.append(" geom_text(data = NULL,x =");
			stringBuilder.append(x + ", y = ");
			stringBuilder.append(y + ", label=");
			stringBuilder.append("\"" + label + "\",angle=" + angle + ",size=" + fontSize  + " )");
			//+ (color == null ? "" : ",color=\"" + color + "\"")
			return stringBuilder.toString();
		}

		public static GeomText[] fromFile(String file, Logger log) {
			if (!Files.exists(file)) {
				log.reportFileNotFound(file);
				return null;
			} else if (!Files.headerOfFileContainsAll(file, HEADER, log)) {
				log.reportTimeError("Geom text "+file+"  must have header " + Array.toStr(HEADER));
				return null;
			} else {
				ArrayList<GeomText> geomTexts = new ArrayList<GeomText>();
				try {
					int[] indices  = ext.indexFactors(HEADER, Files.getHeaderOfFile(file, log), true, false);
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
			this.textSize = 3;
			this.textAngle=0;
			
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
		private double fontsize;
		private String title;
		private String plotVar;
		private boolean logTransformX;
		private boolean logTransformY;
		private boolean overWriteExisting;
		private int height;
		private int width;
		private GeomText[] gTexts;
		private SeriesLabeler seriesLabeler;
		private boolean directLableGtexts, onlyMaxMin;		

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
			this.height = 6;
			this.width = 11;
			this.directLableGtexts =false;
			this.onlyMaxMin = false;
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
			// log.report(Array.toStr(rScript, "\n"));
			Files.writeList(rScript, rScriptFile);
			boolean ran = CmdLine.runCommandWithFileChecks(new String[] { "Rscript", rScriptFile }, "", new String[] { rScriptFile }, new String[] { output }, true, overWriteExisting, false, log);
			return ran;
		}

		private GeomText[] getLabelers(){
			if (seriesLabeler != null) {
				ArrayList<GeomText> geomTexts = new ArrayList<GeomText>();
				double[] xs = Array.toDoubleArray(HashVec.loadFileToStringArray(dataFile, true, new int[] { ext.indexOfStr(dataXvalueColumn, Files.getHeaderOfFile(dataFile, log)) }, false));
				for (int i = 0; i < dataYvalueColumns.length; i++) {
					int maxIndex = -1;
					int minIndex = -1;
					double min = Double.MAX_VALUE;
					double max = -1*Double.MAX_VALUE;

					double[] data = Array.toDoubleArray(HashVec.loadFileToStringArray(dataFile, true, new int[] { ext.indexOfStr(dataYvalueColumns[i], Files.getHeaderOfFile(dataFile, log)) }, false));
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
						String maxLabel = "MAX_" + dataYvalueColumns[i] + "_" + seriesLabeler.getIndexTag() + "-" + Math.round(xs[maxIndex]) + " " + seriesLabeler.getValueTag() + "-" + ext.formDeci(data[maxIndex], 2);
						GeomText maxText = new GeomText(xs[maxIndex], max, seriesLabeler.getTextAngle(), maxLabel, seriesLabeler.getTextSize());
						//maxText.setColor(dataYvalueColumns[i]);
						geomTexts.add(maxText);
					}
					if (seriesLabeler.isLabelMin()) {

						String minLabel = "MIN_" + dataYvalueColumns[i] + "_" +seriesLabeler.getIndexTag() + "-" + Math.round(xs[minIndex]) + " " + seriesLabeler.getValueTag() + "-" + ext.formDeci(data[minIndex], 2);
						GeomText minText = new GeomText(xs[minIndex], min, seriesLabeler.getTextAngle(), minLabel, seriesLabeler.getTextSize());
					//	minText.setColor(dataYvalueColumns[i]);
						geomTexts.add(minText);
					}
				}

				return geomTexts.toArray(new GeomText[geomTexts.size()]);

			}
			return null;
		}
		
		private String[] directLabelFrame(String currentPlot,String factorTitel,String xTitle,String meltYTitle){
			String[] directDataFrame= new String[2];
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

			String df = currentPlot+"_direcLdf";
			directDataFrame[0]=df;
			String dfGenerate = df + " <-  data.frame(" + xTitle + "=" + xVector + "," + factorTitel + "=" + typeVector + "," + meltYTitle + "=" + yVector + ")";
			directDataFrame[1] = dfGenerate;
			return directDataFrame;
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
				if (directLableGtexts) {
					rCmd.add(" library(directlabels)");
				}
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
				if (sType == SCATTER_TYPE.BOX) {
					plot += plotVar + " <- ggplot(" + dataTable + ",aes(x=variable, y=value)) ";
					plot += " + " + sType.getCall() + "(aes(fill=" + rSafeXColumn + "))";

				} else {
					if (directLableGtexts && gTexts != null) {
						String[] directFrame = directLabelFrame(plotVar, "variable", rSafeXColumn, "value");
						rCmd.add(directFrame[1]);
						String aes = "data=" + directFrame[0] + " , aes(colour =variable,x=" + rSafeXColumn + ", y=value, group=variable)";
						plot += plotVar + " <- direct.label( ggplot(" + aes + ") + " + sType.getCall() + "(" + aes + "),  list(\"smart.grid\", cex=" + fontsize + "))";
					} else {
						plot += plotVar + " <- ggplot()";

					}
					// , col=variable
					if (!onlyMaxMin) {
						plot += " + " + sType.getCall() + "(data=" + dataTable + " , aes(colour =variable,x=" + rSafeXColumn + ", y=value, group=variable)" + (gPoint_SIZE == null ? "" : "," + gPoint_SIZE.getCall()) + ")";
					}
					System.out.println(plot);
				}

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
				if (gTexts != null&&!directLableGtexts) {
					for (int i = 0; i < gTexts.length; i++) {
						plot += " + " + gTexts[i].getCommand();
					}
				}
				rCmd.add(plot);

				rCmd.add("ggsave(file=\"" + output + "\", width= " + width + ", height = " + height + ",limitsize=FALSE)");
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
