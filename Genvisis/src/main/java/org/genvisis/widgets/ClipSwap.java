package org.genvisis.widgets;

import java.io.File;
import java.util.Date;
import java.util.HashSet;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.plots.TwoDPlot;
import org.genvisis.common.*;
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
//		str = ext.replaceAllWith(str, "qwerty", "\\\\");
//		ext.setClipboard(str+"\\\\");
		str = ext.replaceAllWith(str, "qwerty", "/");
		ext.setClipboard(str+"/");
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
		for (int i = 0; i < lines.length; i++) {
			line = lines[i].split("\t", -1);
			for (int j = 0; j < line.length; j++) {
				result += (j==0?"":"\t") + "=\""+ext.prettyP(line[j], 2, 4, 2, true)+"\"";
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
        tdp.setVisible(true);
//        tdp.windowClosing(null);
//        tdp = null;
		
        // TODO set to histogram image location?  set to histogram image?
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
					System.out.println("Line # " + (i + 1) + " had more than one column: " + Array.toStr(line, " / "));
				}
				countMoreThan++;
			}
			if (line[0].startsWith("sigfigs=")) {
				sigfigsExtrastep = new int[] { ext.parseIntArg(line[0]), 0 };
			} else if (line[0].startsWith("step=")) {
				sigfigsExtrastep = Histogram.reverseStep(ext.parseDoubleArg(line[0]));
				System.out.println("Overriding step to be " + ext.formDeci(Histogram.determineStep(sigfigsExtrastep[0], sigfigsExtrastep[1]), sigfigsExtrastep[0] + 2));
			} else if (ext.isMissingValue(line[0])) {
				if (countInvalids < 5) {
					System.out.println("Line # " + (i + 1) + " had an invalid double: " + Array.toStr(line, " / "));
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
			histo = new Histogram(array, Array.min(array), Array.max(array), sigfigsExtrastep[0], sigfigsExtrastep[1]);
		}
		return histo;
	}

	public static void inverseVarianceMeta() {
//		System.out.println(Array.toStr(MetaAnalysis.inverseVarianceWeighting(ext.getClipboard().trim().split("\\n"), new Logger()), "/"));
		ext.setClipboard(Array.toStr(MetaAnalysis.inverseVarianceWeighting(ext.getClipboard().trim().split("\\n"), new Logger())));
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

		values = HashVec.getKeys(hash);
		result = Array.toStr(values)+"\r\n";
		for (int i = 0; i < lines.length; i++) {
			for (int j = 0; j < values.length; j++) {
				result += (j==0?"":"\t") + (ext.isMissingValue(lines[i])?".":(lines[i].equals(values[j])?"1":"0"));
			}
			result += "\r\n";
		}
		
		ext.setClipboard(result);
	}

	public static void saveKeysToFile() {
		long time = new Date().getTime();
		
		new SerialStringArray(ext.getClipboard().trim().split("\\n")).serialize("savedKeysForLookup.ser");

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

		if (new Date().getTime() - new File("savedKeysForLookup.ser").lastModified() > 10 * 60 * 1000) { // if file is more than 10 minutes old
			log.reportError("Our expected file 'savedKeysForLookup.ser' was created "+ext.getTimeElapsed(new File("savedKeysForLookup.ser").lastModified())+" ago; resave keys using -saveKeys");
			return false;
		}
		
		keys = SerialStringArray.load("savedKeysForLookup.ser", false).getArray();
		if (new Date().getTime() - time > 1 * 1000) { // report if it took more than a second
			log.report("Loaded serialized keys in " + ext.getTimeElapsed(time));
		}

		time = new Date().getTime();
		temp = new StringBuilder();
		try {
			result = Files.combineInMemory(keys, Matrix.toMatrix(ext.getClipboard().trim().split("\\n"), "\t"), ".", false, false, log);
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
					temp.append(Array.toStr(result[i])+lineEnding);
				}
			}
			if (count > 0) {
				log.report("(There "+(count==1?"was one key that was":"were "+count+" keys that were")+" blank)"+lineEnding);
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
	    boolean nominalVariable = false;
    	boolean saveKeysToFile = false;
    	boolean lookupValuesForSavedKeys = false;
    	
	    String usage = "\n"+
	    "widgets.ClipSwap requires 0-1 arguments\n"+
	    "   (1) Fix slashes (i.e. -slash (not the default))\n"+
	    "   (2) Contracts contents of clipboard (i.e. -contract (not the default))\n"+
	    "   (3) Expands contents of clipboard (i.e. -expand (not the default))\n"+
	    "   (4) Find unique set and count counts (i.e. -unique (not the default))\n"+
	    "   (5) Remove formatting, leaving only plain text (i.e. -removeFormatting (not the default))\n"+
	    "   (6) Make p-values pretty (i.e. -prettyP (not the default))\n"+
	    "   (7) Create bins and counts for a histogram (i.e. -histogram (not the default))\n"+
	    "   (8) Perform an inverse-variance weighted meta-analysis on a series of betas/stderrs (i.e. -inverseVariance (not the default))\n"+
	    "   (9) Split a nominal variable into binary columns (i.e. -nominalVariable (not the default))\n"+
	    "   (10) Extracts keys from clipboard and saves them to a serialized file, with tabs maintained in lookup values (i.e. -saveKeys (not the default))\n"+
	    "   (11) Lookup the values for the stored keys using the contents of the clipboard (i.e. -lookupValuesForSavedKeys (not the default))\n"+
	    "";

	    for (int i = 0; i<args.length; i++) {
		    if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
			    System.err.println(usage);
			    System.exit(1);
		    } else if (args[i].startsWith("-slash")) {
			    slash = true;
			    numArgs--;
		    } else if (args[i].startsWith("-contract")) {
			    contract = true;
			    numArgs--;
		    } else if (args[i].startsWith("-expand")) {
			    expand = true;
			    numArgs--;
		    } else if (args[i].startsWith("-unique")) {
			    unique = true;
			    numArgs--;
		    } else if (args[i].startsWith("-removeFormatting")) {
			    removeFormatting = true;
			    numArgs--;
		    } else if (args[i].startsWith("-prettyP")) {
		    	prettyP = true;
			    numArgs--;
		    } else if (args[i].startsWith("-histogram")) {
		    	histogram = true;
			    numArgs--;
		    } else if (args[i].startsWith("-inverseVariance")) {
		    	inverseVariance = true;
			    numArgs--;
		    } else if (args[i].startsWith("-nominalVariable")) {
		    	nominalVariable = true;
			    numArgs--;
		    } else if (args[i].startsWith("-saveKeys")) {
		    	saveKeysToFile = true;
			    numArgs--;
		    } else if (args[i].startsWith("-lookupValuesForSavedKeys")) {
		    	lookupValuesForSavedKeys = true;
			    numArgs--;
		    } else {
			    System.err.println("Error - don't know what to do with argument '"+args[i]+"'");
		    }
	    }
	    if (numArgs!=0) {
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
	    	
	    } catch (Exception e) {
		    e.printStackTrace();
		    ext.waitForResponse("Some sort of error");
	    }
    }
}
