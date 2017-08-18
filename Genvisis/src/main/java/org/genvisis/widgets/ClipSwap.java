package org.genvisis.widgets;

import java.awt.Container;
import java.io.File;
import java.util.Arrays;
import java.util.Date;
import java.util.HashSet;

import javax.swing.JFrame;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.plots.TwoDPlot;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.DoubleVector;
import org.genvisis.common.Elision;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.Unique;
import org.genvisis.common.ext;
import org.genvisis.filesys.SerialStringArray;
import org.genvisis.gwas.MetaAnalysis;
import org.genvisis.stats.Histogram;

import com.google.common.primitives.Doubles;


public class ClipSwap {
	public static final String CONTRACT_EXPAND_DELIM = "^";

	public static void contract() {
		String str = ext.getClipboard();

		if (CONTRACT_EXPAND_DELIM.length() == 1) {
			str = str.replace('\t', CONTRACT_EXPAND_DELIM.charAt(0));
		} else {
			str = ext.replaceAllWith(str, "\t", CONTRACT_EXPAND_DELIM);
		}

		ext.setClipboard(str);
	}

	public static void expand() {
		String str = ext.getClipboard();

		if (CONTRACT_EXPAND_DELIM.length() == 1) {
			str = str.replace(CONTRACT_EXPAND_DELIM.charAt(0), '\t');
		} else {
			str = ext.replaceAllWith(str, CONTRACT_EXPAND_DELIM, "\t");
		}

		ext.setClipboard(str);
	}

	public static void slash() {
		String str = ext.getClipboard();
		str = ext.replaceAllWith(str, "\\\\", "\\");
		str = ext.replaceAllWith(str, "\\", "qwerty");
		// str = ext.replaceAllWith(str, "qwerty", "\\\\");
		// ext.setClipboard(str+"\\\\");
		str = ext.replaceAllWith(str, "qwerty", "/");
		ext.setClipboard(str + "/");
	}

	public static void unique() {
		ext.setClipboard(Unique.proc(ext.getClipboard().trim().split("\\n")));
	}

	public static void removeFormatting() {
		ext.setClipboard(ext.getClipboard());
	}

	public static void prettyP() {
		String[] lines, line;
		String result;

		lines = ext.getClipboard().trim().split("\\n");

		result = "";
		for (String line2 : lines) {
			line = line2.split("\t", -1);
			for (int j = 0; j < line.length; j++) {
				result += (j == 0 ? "" : "\t") + "=\"" + ext.prettyP(line[j], 2, 4, 2, true) + "\"";
			}
			result += "\r\n";
		}

		ext.setClipboard(result);
	}

	public static void histogram() {
		Histogram histo = getHistogram();

		String file = "./histograms/clipboard_histogram.png";
		int cnt = 1;
		while (Files.exists(file)) {
			file = "./histograms/clipboard_histogram_" + cnt++ + ".png";
		}
		TwoDPlot tdp = TwoDPlot.createGUI(new Project(), false, false);
		tdp.setHistogram(true);
		tdp.getPanel().overrideAxisLabels("Bins", "");
		tdp.getPanel().setHistogramOverride(true);
		tdp.getPanel().setHistogram(histo);
		tdp.getPanel().createImage();
		tdp.getPanel().screenCapture(file);

		Container p = tdp.getParent();
		while (!(p instanceof JFrame)) {
			p = p.getParent();
		}
		p.setVisible(true);

		System.out.println("Wrote to " + new File(file).getAbsolutePath());
		// TODO set to histogram image location? set to histogram image?
		ext.setClipboard(histo.getSummary());
	}

	/**
	 * @return a {@link Histogram} populated from the current clipboard
	 */
	public static Histogram getHistogram() {
		String[] lines, line;
		DoubleVector dv;
		int countMoreThan, countInvalids;
		Histogram histo;
		int[] sigfigsExtrastep;
		double[] array;

		sigfigsExtrastep = null;
		countMoreThan = countInvalids = 0;
		lines = ext.getClipboard().trim().split("\\n");
		dv = new DoubleVector();
		for (int i = 0; i < lines.length; i++) {
			line = lines[i].split("\t", -1);
			if (line.length > 1) {
				if (countMoreThan < 5) {
					System.out.println("Line # " + (i + 1) + " had more than one column: "
														 + ArrayUtils.toStr(line, " / "));
				}
				countMoreThan++;
			}
			if (line[0].startsWith("sigfigs=")) {
				sigfigsExtrastep = new int[] {ext.parseIntArg(line[0]), 0};
			} else if (line[0].startsWith("step=")) {
				sigfigsExtrastep = Histogram.reverseStep(ext.parseDoubleArg(line[0]));
				System.out.println("Overriding step to be "
													 + ext.formDeci(
																					Histogram.determineStep(sigfigsExtrastep[0],
																																	sigfigsExtrastep[1]),
																					sigfigsExtrastep[0] + 2));
			} else if (ext.isMissingValue(line[0])) {
				if (countInvalids < 5) {
					System.out.println("Line # " + (i + 1) + " had an invalid double: "
														 + ArrayUtils.toStr(line, " / "));
				}
				countInvalids++;
			} else {
				dv.add(Double.parseDouble(line[0]));
			}
		}
		if (countMoreThan >= 5) {
			System.out.println("There were " + countMoreThan + " lines with more than one column");
		}
		if (countInvalids >= 5) {
			System.out.println("There were " + countInvalids + " invalid doubles in the data");
		}

		array = Doubles.toArray(dv);
		if (sigfigsExtrastep == null) {
			histo = new Histogram(array);
		} else {
			histo = new Histogram(array, ArrayUtils.min(array), ArrayUtils.max(array),
														sigfigsExtrastep[0], sigfigsExtrastep[1]);
		}
		return histo;
	}

	public static void inverseVarianceMeta() {
		// System.out.println(Array.toStr(MetaAnalysis.inverseVarianceWeighting(ext.getClipboard().trim().split("\\n"),
		// new Logger()), "/"));
		ext.setClipboard(ArrayUtils.toStr(MetaAnalysis.inverseVarianceWeighting(ext.getClipboard()
																																							 .trim().split("\\n"),
																																						new Logger())));
	}

	public static void convertBetaSE2OR_CI(int sigfigs) {
		String[] lines = ext.getClipboard().trim().split("\\n");
		String str = "";

		for (int i = 0; i < lines.length; i++) {
			String[] line = lines[i].trim().split("\t", -1);
			double beta = 999;
			double se = 999;

			if (line.length != 2) {
				ext.waitForResponse("ClipSwap/Convert beta/SE to OR (95% CI) was called, but did not find two tokens to parse");
				return;
			}
			try {
				beta = Double.parseDouble(line[0]);
			} catch (Exception e) {
				e.printStackTrace();
				ext.waitForResponse("ClipSwap/Convert beta/SE to OR (95% CI) was called, but could not parse beta value from line "
														+ (i + 1) + ": " + line[0]);
				return;
			}

			try {
				se = Double.parseDouble(line[1]);
			} catch (Exception e) {
				e.printStackTrace();
				ext.waitForResponse("ClipSwap/Convert beta/SE to OR (95% CI) was called, but could not parse beta value from line "
														+ (i + 1) + ": " + line[0]);
				return;
			}

			str += ext.formDeci(Math.exp(beta), sigfigs) + " ("
						 + ext.formDeci(Math.exp(beta - 1.96 * se), sigfigs) + ", "
						 + ext.formDeci(Math.exp(beta + 1.96 * se), sigfigs) + ")\n";
		}

		ext.setClipboard(str);
	}

	public static void nominalVariable() {
		String[] lines, values;
		String result;
		HashSet<String> hash;

		lines = ext.getClipboard().trim().split("\\n");

		hash = new HashSet<String>();
		for (int i = 0; i < lines.length; i++) {
			if (!ext.isMissingValue(lines[i])) {
				hash.add(lines[i]);
			}
		}

		values = hash.toArray(new String[hash.size()]);
		Arrays.sort(values);
		result = ArrayUtils.toStr(values) + "\r\n";
		for (String line : lines) {
			for (int j = 0; j < values.length; j++) {
				result += (j == 0 ? "" : "\t")
									+ (ext.isMissingValue(line) ? "." : (line.equals(values[j]) ? "1" : "0"));
			}
			result += "\r\n";
		}

		ext.setClipboard(result);
	}

	public static void saveKeysToFile() {
		long time = new Date().getTime();

		new SerialStringArray(ext.getClipboard().trim()
														 .split("\\n")).serialize("savedKeysForLookup.ser");

		if (new Date().getTime() - time > 1 * 1000) { // report if it took more than a second
			System.out.println("Converted lookup result to text in " + ext.getTimeElapsed(time));
		}
	}

	public static boolean lookupValuesForSavedKeys() {
		String[] keys;
		String[][] result;
		StringBuilder temp;
		long time;
		Logger log;
		String lineEnding;
		int count;

		log = new Logger();
		log.setLevel(11);

		time = new Date().getTime();
		if (!Files.exists("savedKeysForLookup.ser")) {
			log.reportError("Our expected file 'savedKeysForLookup.ser' does not exist; run -saveKeys first");
			return false;
		}

		if (new Date().getTime() - new File("savedKeysForLookup.ser").lastModified() > 10 * 60 * 1000) { // if
																																																		 // file
																																																		 // is
																																																		 // more
																																																		 // than
																																																		 // 10
																																																		 // minutes
																																																		 // old
			log.reportError("Our expected file 'savedKeysForLookup.ser' was created "
											+ ext.getTimeElapsed(new File("savedKeysForLookup.ser").lastModified())
											+ " ago; resave keys using -saveKeys");
			return false;
		}

		keys = SerialStringArray.load("savedKeysForLookup.ser").getArray();
		if (new Date().getTime() - time > 1 * 1000) { // report if it took more than a second
			log.report("Loaded serialized keys in " + ext.getTimeElapsed(time));
		}

		time = new Date().getTime();
		temp = new StringBuilder();
		try {
			result = Files.combineInMemory(keys,
																		 Matrix.toMatrix(ext.getClipboard().trim().split("\\n"), "\t"),
																		 ".", false, false, log);
			if (new Date().getTime() - time > 1 * 1000) { // report if it took more than a second
				log.report("Performed lookup in " + ext.getTimeElapsed(time));
			}

			if (Files.isWindows()) {
				lineEnding = "\r\n";
			} else {
				lineEnding = "\n";
			}

			count = 0;
			time = new Date().getTime();
			for (int i = 0; i < keys.length; i++) {
				if (keys[i].trim().equals("")) {
					temp.append(lineEnding);
					count++;
				} else {
					temp.append(ArrayUtils.toStr(result[i]) + lineEnding);
				}
			}
			if (count > 0) {
				log.report("(There "
									 + (count == 1 ? "was one key that was" : "were " + count + " keys that were")
									 + " blank)" + lineEnding);
			}

			if (new Date().getTime() - time > 1 * 1000) { // report if it took more than a second
				log.report("Converted lookup result to text in " + ext.getTimeElapsed(time));
			}
		} catch (Elision e) {
			log.reportError(e.getMessage());
			return false;
		}

		time = new Date().getTime();
		ext.setClipboard(temp.toString());
		if (new Date().getTime() - time > 1 * 1000) { // report if it took more than a second
			log.report("Posted result to clipboard in " + ext.getTimeElapsed(time));
		}

		return true;
	}

	public static void convertToUCSC() {
		String[][] strs = Matrix.toMatrix(ext.getClipboard().trim().split("\\n"), "\t");
		String[] result = new String[strs.length];
		int numColErrs = 0;
		int numParseErrs = 0;

		for (int i = 0; i < strs.length; i++) {
			if (strs[i].length != 3) {
				if (numColErrs < 10) {
					System.err.println("Error - row #" + (i + 1) + " has " + (strs[i].length)
														 + " columns instead of 3");
					System.err.println("        " + ArrayUtils.toStr(strs[i], "/"));
				}
				numColErrs++;
				result[i] = ".";
			} else {
				if (strs[i][0].toLowerCase().startsWith("chr")) {
					strs[i][0] = strs[i][0].substring(3, strs[i][0].length());
				}
				for (int j = 0; j < 3; j++) {
					try {
						Integer.parseInt(strs[i][j]);
					} catch (NumberFormatException nfe) {
						if (numParseErrs < 10) {
							System.err.println("Error - row #" + (i + 1) + " does not have a valid "
																 + (j == 0 ? "chromosome" : (j == 1 ? "start" : "stop"))
																 + " number: " + strs[i][0]);
						}
						numParseErrs++;
					}
				}
				result[i] = "chr" + Integer.parseInt(strs[i][0]) + ":" + Integer.parseInt(strs[i][1]) + "-"
										+ Integer.parseInt(strs[i][2]);
			}
		}
		if (numColErrs >= 10) {
			System.err.println("There were " + numColErrs + " mismatched column errors");
		}
		if (numParseErrs >= 10) {
			System.err.println("There were " + numColErrs + " invalid number errors");
		}

		ext.setClipboard(ArrayUtils.toStr(result, "\n"));
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		boolean slash = false;
		boolean unique = false;
		boolean contract = false;
		boolean expand = false;
		boolean removeFormatting = false;
		boolean prettyP = false;
		boolean histogram = false;
		boolean inverseVariance = false;
		boolean convertBetaSE2OR_CI = false;
		boolean nominalVariable = false;
		boolean saveKeysToFile = false;
		boolean lookupValuesForSavedKeys = false;
		boolean convertToUCSC = false;
		int sigfigs = 3;

		String usage = "\n"
									 + "widgets.ClipSwap requires 0-1 arguments\n"
									 + "   (1) Fix slashes (i.e. -slash (not the default))\n"
									 + "   (2) Contracts contents of clipboard (i.e. -contract (not the default))\n"
									 + "   (3) Expands contents of clipboard (i.e. -expand (not the default))\n"
									 + "   (4) Find unique set and count counts (i.e. -unique (not the default))\n"
									 + "   (5) Remove formatting, leaving only plain text (i.e. -removeFormatting (not the default))\n"
									 + "   (6) Make p-values pretty (i.e. -prettyP (not the default))\n"
									 + "   (7) Create bins and counts for a histogram (i.e. -histogram (not the default))\n"
									 + "   (8) Perform an inverse-variance weighted meta-analysis on a series of betas/stderrs (i.e. -inverseVariance (not the default))\n"
									 + "   (9) Convert a beta/SE pair to an OR (95% CI) pair (i.e. -convertBetaSE2OR_CI (not the default))\n"
									 + "            using X number of significant digits (i.e. sigfigs="
									 + sigfigs
									 + " (default))\n"
									 + "   (10) Split a nominal variable into binary columns (i.e. -nominalVariable (not the default))\n"
									 + "   (11) Extracts keys from clipboard and saves them to a serialized file, with tabs maintained in lookup values (i.e. -saveKeys (not the default))\n"
									 + "   (12) Lookup the values for the stored keys using the contents of the clipboard (i.e. -lookupValuesForSavedKeys (not the default))\n"
									 + "   (13) convert three columns of chr,start,stop into a single ucsc column chr:start-stop (i.e. -convertToUCSC (not the default))\n"
									 + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("-slash")) {
				slash = true;
				numArgs--;
			} else if (arg.startsWith("-contract")) {
				contract = true;
				numArgs--;
			} else if (arg.startsWith("-expand")) {
				expand = true;
				numArgs--;
			} else if (arg.startsWith("-unique")) {
				unique = true;
				numArgs--;
			} else if (arg.startsWith("-removeFormatting")) {
				removeFormatting = true;
				numArgs--;
			} else if (arg.startsWith("-prettyP")) {
				prettyP = true;
				numArgs--;
			} else if (arg.startsWith("-histogram")) {
				histogram = true;
				numArgs--;
			} else if (arg.startsWith("-inverseVariance")) {
				inverseVariance = true;
				numArgs--;
			} else if (arg.startsWith("-convertBetaSE2OR_CI")) {
				convertBetaSE2OR_CI = true;
				numArgs--;
			} else if (arg.startsWith("sigfigs=")) {
				sigfigs = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("-nominalVariable")) {
				nominalVariable = true;
				numArgs--;
			} else if (arg.startsWith("-saveKeys")) {
				saveKeysToFile = true;
				numArgs--;
			} else if (arg.startsWith("-lookupValuesForSavedKeys")) {
				lookupValuesForSavedKeys = true;
				numArgs--;
			} else if (arg.startsWith("-convertToUCSC")) {
				convertToUCSC = true;
				numArgs--;
			} else {
				System.err.println("Error - don't know what to do with argument '" + arg + "'");
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			if (slash) {
				slash();
			}
			if (contract) {
				contract();
			}
			if (expand) {
				expand();
			}
			if (unique) {
				unique();
			}
			if (removeFormatting) {
				removeFormatting();
			}
			if (prettyP) {
				prettyP();
			}
			if (histogram) {
				histogram();
				ext.waitForResponse("done!");
			}
			if (inverseVariance) {
				inverseVarianceMeta();
			}
			if (convertBetaSE2OR_CI) {
				convertBetaSE2OR_CI(sigfigs);
			}
			if (nominalVariable) {
				nominalVariable();
			}
			if (saveKeysToFile) {
				saveKeysToFile();
				ext.waitForResponse("Keys have been saved; ready to perform lookup\npress any key to continue");
			}
			if (lookupValuesForSavedKeys) {
				if (lookupValuesForSavedKeys()) {
					ext.waitForResponse("Lookup was successful!\npress any key to close");
				} else {
					ext.waitForResponse("Lookup was not successful\npress any key to close");
				}
			}
			if (convertToUCSC) {
				convertToUCSC();
			}

		} catch (Exception e) {
			e.printStackTrace();
			ext.waitForResponse("Some sort of error");
		}
	}
}
