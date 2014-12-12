package common;

import java.io.*;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryUsage;
import java.text.*;
import java.util.*;
import java.awt.Toolkit;
import java.awt.datatransfer.*;

public class ext {
	public static final String[][] VALID_CHARS = { {" ", "_"}, {"?", "_"}, {"'", "_"}, {"/", "_"}, {"\\", "_"}, {"<", "_"}, {">", "_"}, {"|", "_"}, {":", "_"}, {"*", "_"}, {"\"", "_"} };
	public static final String[][] TO_BE_EXTRA_SAFE_CHARS = { {"(", "_"}, {")", "_"}, {"&", "_"} };
	public static final String[] MISSING_VALUES = {"", ".", "NA", "NaN", "x", "#N/A"};
	public static final String[][] META_REPLACEMENTS = { {"{Tab}", "\t"}, {"{Space}", " "}, {"{!}", "!"}, {"{#}", "#"}, {"{+}", "+"}};
	public static final String[] COMMON_IDS = {"id", "IID", "IndID", "gwas_id"};
	
	// safer, reports spurious errors when called a lot 
	public static String replaceAllWithSafer(String str, String from, String to) {
		boolean done;
		String original = str;
		int index;
		
		done = false;
		while (!done) {
			index = str.indexOf(from);
			if (index == -1) {
				done = true;
			} else {
				try {
					str = str.substring(0, str.indexOf(from))+to+str.substring(str.indexOf(from)+from.length());
				} catch (Exception e) {
					System.err.println("Error with "+from+" to "+to+": "+original);
					e.printStackTrace();
				}
			}
		}
		return str;
	}

	public static String replaceAllWith(String str, String from, String to) {
		boolean done;
		int index;

		done = false;
		while (!done) {
			index = str.indexOf(from);
			if (index == -1) {
				done = true;
			} else {
				str = str.substring(0, str.indexOf(from))+to+str.substring(str.indexOf(from)+from.length());
			}
		}
		return str;
	}
	
	public static String replaceQuotesWithSlashQuotes(String original) {
		char[] array;
		String str;
		
		array = original.toCharArray();
		
		str = "";
		for (int i = 0; i < array.length; i++) {
			if (array[i] == '"') {
				str += "\\";
			}
			str += array[i];
		}

		return str;
	}
	
	public static String[][] processMetaCharacters(String[][] replaces) {
		String[][] metaReplaces;
		
		metaReplaces = new String[replaces.length][];
		for (int i = 0; i < replaces.length; i++) {
			metaReplaces[i] = new String[] {replaces[i][0], replaceAllWith(replaces[i][1], META_REPLACEMENTS)};
		}
		
		return metaReplaces;		
	}

	public static String replaceAllWith(String str, String[][] replaces) {
    	for (int i = 0; i<replaces.length; i++) {
    		str = ext.replaceAllWith(str, replaces[i][0], replaces[i][1]);
        }
		return str;
	}

	public static String chrome(int chr) {
		return (chr<10)?"0"+chr:""+chr;
	}

	public static int gender(String gender) {
		if (gender.toLowerCase().equals("m")||gender.equals("1")) {
			return 1;
		} else if (gender.toLowerCase().equals("f")||gender.equals("2")) {
			return 2;
		} else {
			return 0;
		}
	}

	public static String prettyUpDistance(int dist, int numSigFigs) {
		String str = dist+"";
		if (str.length()>6) {
			return formDeci((double)Integer.parseInt(str)/1000000, numSigFigs, true)+" Mb";
		}
		if (str.length()>3) {
			return formDeci((double)Integer.parseInt(str)/1000, numSigFigs, true)+" kb";
		}
		return str;
	}

	public static String prettyUpSize(long size, int numSigFigs) {
		String str = size+"";
		if (str.length()>12) {
			return formDeci((double)size/1048576000000l, numSigFigs, true)+" TB";
		}
		if (str.length()>9) {
			return formDeci((double)size/1048576000, numSigFigs, true)+" GB";
		}
		if (str.length()>6) {
			return formDeci((double)size/1048576, numSigFigs, true)+" MB";
		}
		if (str.length()>3) {
			return formDeci((double)size/1024, numSigFigs, true)+" KB";
		}
		return str;
	}

	public static int indexOfInt(int target, int[] array) {
		int index = -1;

		for (int i = 0; i<array.length; i++) {
			if (array[i]==target) {
				if (index==-1) {
					index = i;
				} else {
					return -2;
				}
			}
		}
		return index;
	}

	public static int indexOfStr(String target, String[][] array, boolean caseSensitive) {
		for (int i = 0; i<array.length; i++) {
			if (caseSensitive?array[i][0].equals(target):array[i][0].toLowerCase().equals(target.toLowerCase())) {
				return i;
			}
		}
		return -1;
	}

	public static int indexOfAnyStr(String target, String[][] array, boolean caseSensitive) {
		for (int i = 0; i<array.length; i++) {
			for (int j = 0; j<array[i].length; j++) {
				if (caseSensitive?array[i][j].equals(target):array[i][j].equalsIgnoreCase(target)) {
					return i;
				}
			}
		}
		return -1;
	}

	public static int indexOfStr(String target, String[] array) {
		return indexOfStr(target, array, true, true, new Logger(), false);
	}

	public static int indexOfStr(String target, String[] array, boolean caseSensitive, boolean exactMatch) {
		return indexOfStr(target, array, caseSensitive, exactMatch, new Logger(), false);
	}

	public static int indexOfStr(String target, String[] array, boolean caseSensitive, boolean exactMatch, Logger log, boolean verbose) {
		int[] indices;
		
		indices = indicesOfStr(target, array, caseSensitive, exactMatch);

		if (indices.length == 0) {
			if (verbose) {
				log.reportError("Error - '"+target+"' was not found in array");
			}
			return -1;
		}
		
		if (indices.length > 1 && verbose) {
			log.reportError("Warning - '"+target+"' was found more than once in the array");
		}
		
		return indices[0];
	}
	
	public static int[] indicesOfStr(String target, String[] array, boolean caseSensitive, boolean exactMatch) {
		IntVector indices;
		
		indices = new IntVector();
		for (int i = 0; i<array.length; i++) {
			if (exactMatch) {
				if (caseSensitive?array[i].equals(target):array[i].toLowerCase().equals(target.toLowerCase())) {
					indices.add(i);
				}
			} else {
				if (caseSensitive?array[i].contains(target)||target.contains(array[i]):array[i].toLowerCase().contains(target.toLowerCase())||target.toLowerCase().contains(array[i].toLowerCase())) {
					indices.add(i);
				}
			}
		}
		
		return indices.toArray();
	}

	public static int[] indicesWithinString(String target, String str) {
		IntVector iv;
		int index, runningIndex;
		
		iv = new IntVector();
		index = -9;
		runningIndex = 0;
		while (index != -1) {
			index = str.indexOf(target);
			if (index != -1) {
				iv.add(runningIndex+index);
				runningIndex += index+target.length();
				str = str.substring(index+target.length());
			}
		}

		return iv.toArray();
	}

	public static int indexOfChar(char target, char[] array) {
		for (int i = 0; i<array.length; i++) {
			if (array[i] == target) {
				return i;
			}
		}

		return -1;
	}

	public static int indexOfStartsWith(String target, String[] array, boolean reverseNotForward) {
		if (reverseNotForward) {
			for (int i = 0; i<array.length; i++) {
				if (target.startsWith(array[i])) {
					return i;
				}
			}
		} else {
			for (int i = 0; i<array.length; i++) {
				if (array[i].startsWith(target)) {
					return i;
				}
			}
		}
		return -1;
	}

	public static int indexOfEndsWith(String target, String[] array, boolean reverseNotForward) {
		if (reverseNotForward) {
			for (int i = 0; i<array.length; i++) {
				if (target.endsWith(array[i])) {
					return i;
				}
			}
		} else {
			for (int i = 0; i<array.length; i++) {
				if (array[i].endsWith(target)) {
					return i;
				}
			}
		}
		return -1;
	}

	/**
	 * Parses the operators, the indices of those operators, and the tokens
	 * 
	 * @param str
	 *            String to be parsed
	 * @param operators
	 *            String of the possible operators
	 * @return a matrix of String[3][] with an array each of the operators, the
	 *         indices of those operators, and the tokens
	 */
	public static String[][] getOperatorsOperatorIndicesAndSplit(String str, String operators) {
		Vector<String> operatorsPresent = new Vector<String>();
		Vector<String> operatorIndices = new Vector<String>();
		Vector<String> tokens = new Vector<String>();
		Vector<String> ops = new Vector<String>();
		int count = 0;
		String hist = "";

		for (int i = 0; i<operators.length(); i++) {
			ops.add(operators.charAt(i)+"");
		}

		while (count<str.length()) {
			if (ops.contains(str.charAt(count)+"")) {
				operatorsPresent.add(str.charAt(count)+"");
				operatorIndices.add(count+"");
				tokens.add(hist);
				hist = "";
			} else {
				hist += str.charAt(count);
			}
			count++;
		}
		tokens.add(hist);

		return new String[][] {Array.toStringArray(operatorsPresent), Array.toStringArray(operatorIndices), Array.toStringArray(tokens)};
	}

	public static boolean eqArrays(Object[] arr1, Object[] arr2) {
		if (arr1.length!=arr2.length) {
			return false;
		}
		for (int i = 0; i<arr1.length; i++) {
			if (!arr1[i].equals(arr2[i])) {
				return false;
			}
		}

		return true;
	}

	public static String formStr(String str, int size) {
		return formStr(str, size, false);
	}

	/**
	 * Adds whitespace to String to fill a specified number of characters
	 * 
	 * @param str
	 *            name
	 * @param size
	 *            number of characters
	 * @param after
	 *            add white space after the string (instead of before)
	 * @return processed string
	 */
	public static String formStr(String str, int size, boolean after) {
		String spaces = "";
		for (int i = 0; i<(size-str.length()); i++) {
			spaces += " ";
		}
		if (after) {
			return str+spaces;
		} else {
			return spaces+str;
		}
	}

	public static String formNum(int num, int size) {
		return formNum(num+"", size);
	}

	public static String formNum(String str, int size) {
		String newNum = "";
		for (int i = 0; i<(size-str.length()); i++) {
			newNum += "0";
		}
		return newNum+str;
	}

	public static String formNum(int num1, int num2, int size1, int size2) {
		return formNum(num1+"", num2+"", size1, size2);
	}

	public static String formNum(String str1, String str2, int size1, int size2) {
		String newNum = "1";
		for (int i = 0; i<(size1-str1.length()); i++) {
			newNum += "0";
		}
		newNum += str1;
		for (int i = 0; i<(size2-str2.length()); i++) {
			newNum += "0";
		}
		return newNum+str2;
	}

	public static String formDeci(double num, int numPlaces) {
		return formDeci(num, numPlaces, false);
	}

	public static String formDeci(double num, int numPlaces, boolean forceDigits) {
		return formDeci(num, forceDigits?numPlaces:0, numPlaces);
	}

	public static String formDeci(double num, int minSigFigs, int maxSigFigs) {
		String theFormat = (maxSigFigs>0?"0.":"0");
		String result;

		if (Double.isNaN(num)) {
			return num+"";
		}

		for (int i = 0; i<minSigFigs; i++) {
			theFormat += "0";
		}

		for (int i = minSigFigs; i<maxSigFigs; i++) {
			theFormat += "#";
		}

		result = new DecimalFormat(theFormat).format(num);
		if (ext.replaceAllWith(result, new String[][] {{"0", ""}, {".", ""}}).equals("-")) {
			result = result.substring(1);
		}
		return result;
	}

	public static String formPercent(double num, int numSigFigs) {
		String theFormat = (numSigFigs>0?"0.":"0");

		if (Double.isNaN(num)) {
			return num+"";
		}

		for (int i = 0; i<numSigFigs; i++) {
			theFormat += "0";
		}

		return new DecimalFormat(theFormat).format(num*100)+"%";
	}

	public static String formSciNot(double num, int numPlaces, boolean forceDigits) {
		DecimalFormat f;
		String theFormat = (numPlaces>0?"0.":"0");

		for (int i = 0; i<numPlaces; i++) {
			theFormat += (forceDigits?"0":"#");
		}

		f = new DecimalFormat(theFormat+"E0");

		return f.format(num);
	}

	public static String prettyP(double d) {
		return prettyP(d+"", 2, 5, 1, false);
	}

	public static String prettyP(String nummer) {
		return prettyP(nummer, 2, 5, 1, false);
	}

	public static String prettyP(double d, int sigfigs, int sigfigWhenConverting, int sigfigWhenConverted, boolean useE) {
		return prettyP(d+"", sigfigs, sigfigWhenConverting, sigfigWhenConverted, useE);
	}

	/**
	 * Converts a String (which is usually a p-value) to a double and then back to a String using scientific notation if beyond a certain magnitude
	 * 
	 * @param nummer
	 *            the original number in String format
	 * @param sigfigs
	 *            the number of significant digits to use if there is no conversion (usually 2)
	 * @param sigfigBeforeConverting
	 *            the number of digits (i.e., the number of zeros after the decimal plus one) required before switching to scientific notation (usually 4 or 5)
	 * @param sigfigWhenConverted
	 *            the number of significant digits to use when using scientific notation (if you use 3, then that means there will be the one digit before the decimal and two after)
	 * @param useE
	 *            use E-0# instead of x 10^-#
	 * @return the pretty version of the number
	 */
	public static String prettyP(String nummer, int sigfigs, int sigfigBeforeConverting, int sigfigWhenConverted, boolean useE) {
		String str = "";
		double d = -999;
		int power;
		String sigdigits, temp;
		
		if (isMissingValue(nummer)) {
			return nummer;
		}

		try {
			d = Double.parseDouble(nummer);
		} catch (Exception e) {
			return "'"+nummer+"' is an ugly number (failed prettyP)";
		}

		if (d>(Math.pow(0.1, sigfigs)-Math.pow(0.1, sigfigs+1)/2)||Double.isNaN(d)) {
			return formDeci(d, sigfigs, true);
		}

		if (d<0) {
			return "0";
		}

		temp = "0";
		for (int i = 1; i<sigfigWhenConverted; i++) {
			if (i==1) {
				temp += ".";
			}
			temp += "0";
		}
		temp += "E0";
		str = new DecimalFormat(temp).format(d);
		sigdigits = str.split("E")[0];
		power = Integer.parseInt(str.split("E")[1])*-1;

		if (power<sigfigBeforeConverting) {
			str = "0.";
			for (int i = 0; i<power-1; i++) {
				str += "0";
			}
			str += ext.replaceAllWith(sigdigits, ".", "");
			return str;
		} else {
			return useE?str:sigdigits+"x10-"+power;
		}
	}

	public static String getDate(Date date, String delimiter) {
		GregorianCalendar gc = new GregorianCalendar();
		gc.setTime(date);
		return gc.get(Calendar.YEAR)+delimiter+chrome((gc.get(Calendar.MONTH)+1))+delimiter+chrome(gc.get(Calendar.DAY_OF_MONTH));
	}

	public static String getDate() {
		return getDate(new Date(), ".");
	}

	public static String getTime() {
		return getTime(new Date().getTime());
	}

	public static String getTime(long time) {
		GregorianCalendar gc = new GregorianCalendar();
		gc.setTime(new Date(time));
		return formNum(gc.get(Calendar.HOUR)==0?12:gc.get(Calendar.HOUR), 2)+":"+formNum(gc.get(Calendar.MINUTE), 2)+":"+formNum(gc.get(Calendar.SECOND), 2)+(gc.get(Calendar.AM_PM)==0?"AM":"PM");
	}
	
	public static void timestamp() {
		System.out.println(getDate()+"\t"+getTime());
	}

	public static String getTimestampForFilename() {
		return replaceAllWith((getDate()+"_"+getTime()), ":", "_");
	}
	
	public static String getTimeElapsed(long startTime) {
		long timeElapsed;
		
		timeElapsed = new Date().getTime()-startTime;
		if (timeElapsed < 10000) {
			return timeElapsed +" ms";
		}
		
		timeElapsed /= 1000;
		
		if (timeElapsed < 60) {
			return timeElapsed +" sec";
		}

		if (timeElapsed < 600) {
			return timeElapsed/60 +" min "+(timeElapsed - (timeElapsed/60)*60)+" sec";
		}

		timeElapsed /= 60;
		
		if (timeElapsed < 60) {
			return timeElapsed +" min";
		}

		if (timeElapsed < 24*60) {
			return timeElapsed/60 +" hrs "+(timeElapsed - (timeElapsed/60)*60)+" min";
		}

		timeElapsed /= 60;
		
		return timeElapsed/24 +" days "+(timeElapsed - (timeElapsed/24)*24)+" hrs";
	}

	public static double getTimeSince(long startTime, char returnType) {
		double timeElapsed;
		
		timeElapsed = new Date().getTime()-startTime;
		timeElapsed /= 1000;

		if (returnType == 'S' || returnType == 's') {
			return timeElapsed;
		}
		
		timeElapsed /= 60;
		
		if (returnType == 'M' || returnType == 'm') {
			return timeElapsed;
		}

		timeElapsed /= 60;
		
		if (returnType == 'H' || returnType == 'h') {
			return timeElapsed;
		}
		
		timeElapsed /= 24;
		
		if (returnType == 'D' || returnType == 'd') {
			return timeElapsed;
		}
		
		System.err.println("Error - '"+returnType+"' is an invalid return type: only S(econds), M(inutes), and H(ours) are allowed");
		
		return -1;
	}

	public static double divide(int num1, int num2) {
		if (num2==0) {
			return 0;
		}
		return (double)num1/(double)num2;
	}

	public static String replaceWithLinuxSafeCharacters(String str, boolean extraSafe) {
		for (int i = 0; i<VALID_CHARS.length; i++) {
			str = replaceAllWith(str, VALID_CHARS[i][0], VALID_CHARS[i][1]);
		}

		if (extraSafe) {
			for (int i = 0; i<TO_BE_EXTRA_SAFE_CHARS.length; i++) {
				str = replaceAllWith(str, TO_BE_EXTRA_SAFE_CHARS[i][0], TO_BE_EXTRA_SAFE_CHARS[i][1]);
			}
		}
		
		return str;
	}

	public static String replaceDirectoryCharsWithUnderscore(String filename, int depth) {
		int index;
		
		filename = replaceAllWith(filename, "\\", "/");
		filename = replaceAllWith(filename, "//", "/");
		
		for (int i = 0; i<depth; i++) {
			index = filename.lastIndexOf("/");
			if (index >= 0) {
				filename = filename.substring(0, index)+"_"+filename.substring(index+1);
			}
		}
		index = filename.lastIndexOf("/");
		if (index >= 0) {
			filename = filename.substring(index+1);
		}
		
		return filename;
	}

//	public static int[] indexFactors(String trait, Vector included, String[] phenoNames) {
//		return indexFactors(trait, included, phenoNames, false);
//	}
//
//	public static int[] indexFactors(String trait, Vector included, String[] phenoNames, boolean possibleEmptySet) {
//		int[] indices;
//		int first, last;
//
//		if (included==null) {
//			included = new Vector();
//		} else if (included.size()==0&&!possibleEmptySet) {
//			throw new RuntimeException("No independent factors listed; flag true if legitimate");
//		}
//
//		if (included.contains(trait)) {
//			System.err.println("Warning - "+trait+" is listed as both a factor and the dependent variable; removed as a factor");
//			included.remove(trait);
//		}
//
//		indices = new int[included.size()+1];
//		for (int i = 0; i<indices.length; i++) {
//			indices[i] = -1;
//		}
//
//		for (int i = 0; i<phenoNames.length; i++) {
//			if (included.contains(phenoNames[i])) {
//				first = included.indexOf(phenoNames[i]);
//				last = included.lastIndexOf(phenoNames[i]);
//				for (int j = first; j<=last; j++) {
//					if (phenoNames[i].equals(included.elementAt(j))) {
//						indices[j+1] = i;
//					}
//				}
//			} else if (phenoNames[i].equals(trait)) {
//				indices[0] = i;
//			}
//		}
//
//		for (int i = 0; i<indices.length; i++) {
//			if (indices[i]==-1) {
//				throw new RuntimeException("Error - variable name '"+(i==0?trait:included.elementAt(i-1))+"' was not found in the database");
//			}
//		}
//
//		return indices;
//	}

	/**
	 * Method very similar to other index factors, but uses a hash based approach. When subset and superset are very large such as extracting marker subsets (i.e. 50K and 2500K), this might be a bit faster
	 * 
	 * @param subset
	 *            : what to extract
	 * @param superset
	 *            : extract indices from
	 * @param casesensitive
	 *            : case sensitive match
	 * @param log
	 * @param verbose
	 *            : report duplicate matches, and no matches
	 * @param kill
	 *            : system.exit(1) if duplicates are found, or subset member is not contained in superset
	 * @return indices of subset in superset
	 * 
	 */
	public static int[] indexLargeFactors(String[] subset, String[] superset, boolean casesensitive, Logger log, boolean verbose, boolean kill) {
		Hashtable<String, Integer> track = new Hashtable<String, Integer>();
		int[] indices = new int[subset.length];
		Arrays.fill(indices, -1);
		boolean err = false;
		for (int i = 0; i < superset.length; i++) {
			track.put(casesensitive ? superset[i] : superset[i].toLowerCase(), i);
		}
		for (int i = 0; i < subset.length; i++) {
			String sub = casesensitive ? subset[i] : subset[i].toLowerCase();
			if (track.containsKey(sub)) {
				int index = track.get(sub);
				if (index >= 0) {
					indices[i] = index;
					track.put(sub, -1);
				} else {
					if (verbose) {
						log.reportError("Error - more than one factor was named '" + subset[i] + "'");
					}
					err = true;
				}
			} else {
				if (verbose) {
					log.reportError("Error - no factor was named '" + subset[i] + "'");
				}
				err = true;
			}
		}
		if (kill && err) {
			System.exit(1);
		}
		return indices;
	}

	public static int[] indexFactors(String[] subset, String[] superset, boolean casesensitive, boolean kill) {
		return indexFactors(subset, superset, casesensitive, new Logger(), true, kill);
	}
	
	public static int[] indexFactors(String[] subset, String[] superset, boolean casesensitive, Logger log, boolean verbose, boolean kill) {
		int[] indices = new int[subset.length];
		boolean err = false;

		for (int i = 0; i<subset.length; i++) {
			indices[i] = -1;
			for (int j = 0; j<superset.length; j++) {
				if (casesensitive?subset[i].equals(superset[j]):subset[i].equalsIgnoreCase(superset[j])) {
					if (indices[i]==-1) {
						indices[i] = j;
					} else {
						log.reportError("Error - more than one factor was named '"+subset[i]+"'");
						err = true;
					}
				}
			}
			if (indices[i]==-1) {
				if (verbose) {
					log.reportError("Error - no factor was named '"+subset[i]+"'");
				}
				err = true;
			}
		}

		if (kill&&err) {
			System.exit(1);
		}

		return indices;
	}

//	public static boolean[] indexPresent(String[] subset, String[] superset, boolean casesensitive) {
//		boolean[] attendance = new boolean[subset.length];
//
//		for (int i = 0; i<subset.length; i++) {
//			attendance[i] = false;
//			for (int j = 0; j<superset.length; j++) {
//				if (casesensitive?subset[i].equals(superset[j]):subset[i].equalsIgnoreCase(superset[j])) {
//					attendance[i] = true;
//				}
//			}
//		}
//
//		return attendance;
//	}
//
	public static int[] indexFactors(String[][] targetsWithAlts, String[] superset, boolean caseSensitive, boolean exactMatch, boolean verbose, boolean kill) {
		return indexFactors(targetsWithAlts, superset, caseSensitive, exactMatch, verbose, new Logger(), kill);
	}
	
	public static int[] indexFactors(String[][] targetsWithAlts, String[] superset, boolean caseSensitive, boolean exactMatch, boolean verbose, Logger log, boolean kill) {
		return indexFactors(targetsWithAlts, superset, false, caseSensitive, exactMatch, verbose, log, kill);
	}
	
	public static int[] indexFactors(String[][] targetsWithAlts, String[] superset, boolean preferFirstInTargetsOverFirstInSuperset, boolean caseSensitive, boolean exactMatch, boolean verbose, Logger log, boolean kill) {
		int[] finalIndices;
		IntVector[] possibleIndices;
		boolean err = false;
		int index, minIndex;

		possibleIndices = IntVector.newIntVectors(targetsWithAlts.length);
		if (preferFirstInTargetsOverFirstInSuperset) {
	 		for (int j = 0; j<targetsWithAlts.length; j++) {
	 			minIndex = Integer.MAX_VALUE;
	 			for (int i = 0; i<superset.length; i++) {
	 				index = indexOfStr(superset[i], targetsWithAlts[j], caseSensitive, exactMatch);
	 				if (index!=-1) {
	 					if (index < minIndex) {
	 						possibleIndices[j].insertElementAt(i, 0);
	 						minIndex = index;
	 					} else {
	 						possibleIndices[j].add(i);
	 					}
	 				}
	 			}
	 		}
		} else {
			for (int i = 0; i<superset.length; i++) {
				for (int j = 0; j<targetsWithAlts.length; j++) {
					if (indexOfStr(superset[i], targetsWithAlts[j], caseSensitive, exactMatch)!=-1) {
						possibleIndices[j].add(i);
					}
				}
			}
		}

		finalIndices = new int[targetsWithAlts.length];
		for (int i = 0; i<targetsWithAlts.length; i++) {
			if (possibleIndices[i].size()==0) {
				if (verbose) {
					log.reportError((kill?"Error":"Warning")+" - no factor named '"+Array.toStr(targetsWithAlts[i], "/")+"'");
				}
				finalIndices[i] = -1;
				err = true;
			} else {
				finalIndices[i] = possibleIndices[i].elementAt(0);
				if (possibleIndices[i].size()>1) {
					if (verbose) {
						log.reportError("Warning - multiple factors matched; using "+superset[possibleIndices[i].elementAt(0)]+" and not ", false, true);
						for (int j = 1; j<possibleIndices[i].size(); j++) {
							log.reportError((j==1?"":" or ")+superset[possibleIndices[i].elementAt(j)], false, true);
						}
						log.reportError("", true, true);
					}
				}
			}
		}
		if (kill&&err) {
			System.exit(1);
		}
		return finalIndices;
	}

	public static boolean checkHeader(String[] observed, String[] expected, boolean kill) {
		return checkHeader(observed, expected, true, kill);
	}
	
	public static boolean checkHeader(String[] observed, String[] expected, boolean caseSensitive, boolean kill) {
		return checkHeader(observed, expected, caseSensitive, new Logger(), kill);
	}
	
	public static boolean checkHeader(String[] observed, String[] expected, boolean caseSensitive, Logger log, boolean kill) {
		boolean kosher = true;

        
        if (observed == null) {
        	System.err.println("Error - null header observed");
        	if (kill) {
        		System.exit(1);
        	}
        	return false;
        }

        if (observed.length!=expected.length) {
			log.reportError("Error - file has an unexpected header; expecting "+expected.length+" columns, found "+observed.length);
			log.reportError("Expected: "+Array.toStr(expected));
			log.reportError("Observed: "+Array.toStr(observed));
			kosher = false;
			if (kill) {
				System.exit(1);
			}
		}
		for (int i = 0; i<(observed.length<expected.length?observed.length:expected.length); i++) {
			if (!observed[i].equalsIgnoreCase(expected[i]) || (caseSensitive && !observed[i].equals(expected[i]))) {
				log.reportError("Error - Expecting "+expected[i]+" in column "+(i+1)+"; got "+observed[i]);
				kosher = false;
			}
		}
		if (kill && !kosher) {
			System.exit(1);
		}

		return kosher;
	}

	public static boolean checkHeader(String[] observed, String[] expected, int[] indicesToCheck, boolean caseSensitive, Logger log, boolean kill) {
        String[] subArray;

        if (indicesToCheck == null) {
        	subArray = observed;
        } else {
        	subArray = new String[indicesToCheck.length];
        	for (int i = 0; i<indicesToCheck.length; i++) {
        		if (indicesToCheck[i] >= observed.length) {
        			System.err.println("Error - indices to checkHeader led to an out of bounds exception");
        			if (kill) {
        				System.exit(1);
        			}
        			return false;
        		}
        		subArray[i] = observed[indicesToCheck[i]];
            }
        }
        
        return checkHeader(subArray, expected, caseSensitive, log, kill);
	}

	public static String getColumn(int col) {
		return (col>=26?(char)(64+(int)(col/26))+"":"")+(char)(64+(col%26)+1);
	}

	public static int getIndexFromColumn(String col) {
		int index = 0;

		if (col.length()>2) {
			System.err.println("Error - column needs to be two characters or less");
			System.exit(1);
		}
		if (col.length()==2) {
			index += ((int)col.charAt(0)-64)*26;
		}
		index += (int)col.charAt(col.length()-1)-64-1;

		return index;
	}

	public static String getClipboard() {
		return getClipboard(true);
	}
	
	public static String getClipboard(boolean verbose) {
		Clipboard systemClipboard;
		Transferable contents;

		try {
			systemClipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
			contents = systemClipboard.getContents(null);
			if (contents==null) {
				return ("Clipboard is empty");
			} else {
				try {
					if (contents.isDataFlavorSupported(DataFlavor.stringFlavor)) {
						return (String)contents.getTransferData(DataFlavor.stringFlavor);
					}
				} catch (UnsupportedFlavorException ufe) {
					ufe.printStackTrace();
				} catch (IOException ioe) {
					ioe.printStackTrace();
				}
			}
		} catch (Exception e) {
			if (verbose) {
				e.printStackTrace();
			}
			return ("Clipboard failure");
		}

		return null;
	}

	public static void setClipboard(String text) {
		setClipboard(text, false);
	}

	public static void setClipboard(String text, boolean html) {
		if (html) {
			Toolkit.getDefaultToolkit().getSystemClipboard().setContents(new HtmlSelection(text), null);
		} else {
			Toolkit.getDefaultToolkit().getSystemClipboard().setContents(new StringSelection(text), null);
		}
	}

	public static String listWithCommas(String[] list) {
		return listWithCommas(list, true);
	}

	public static String listWithCommas(String[] list, boolean finalComma) {
		String str = "";

		for (int i = 0; i<list.length; i++) {
			str += list[i];
			if (i<list.length-1-(finalComma&&list.length>2?0:1)) {
				str += ", ";
			}
			if (list.length>1&&i==list.length-2) {
				str += (list.length==2||!finalComma?" ":"")+"and ";
			}
		}

		return str;
	}

	public static String listRanges(int[] array) {
		String str;
		int[] ordered;
		
		if (array.length == 0) {
			return "";
		}
		
		ordered = Sort.putInOrder(array);
		str = ordered[0]+"";
		for (int i = 1; i<ordered.length; i++) {
			if (ordered[i] - ordered[i-1] > 1) {
				str += ","+ordered[i];
			} else if (i+1 == ordered.length || ordered[i+1] - ordered[i] > 1) {
				str += "-"+ordered[i];
			}
		}

		return str;
	}

	public static String capitalizeFirst(String str) {
		if (str.length()==0) {
			return str;
		} else {
			return str.substring(0, 1).toUpperCase()+str.substring(1);
		}
	}

	public static String uncapitalizeFirst(String str) {
		if (str.length()==0) {
			return str;
		} else {
			return str.substring(0, 1).toLowerCase()+str.substring(1);
		}
	}

	public static double doubleParser(String str) {
		try {
			return Double.parseDouble(str);
		} catch (Exception e) {
			return Double.NEGATIVE_INFINITY;
		}
	}

	public static String getExcelColumn(int index) {
		return ""+(index/26>0?(char)(index/26-1+65):"")+(char)(index%26+65);
	}

	public static String removeQuotes(String str) {
		if (str.startsWith("\"")&&str.endsWith("\"")&&str.length()>1) {
			return str.substring(1, str.length()-1);
		} else {
			return str;
		}
	}

	
	public static String removeQuotesFromExcelToken(String ori) {
		return removeQuotesFromExcelToken(ori, new String[][] {{"\t", "[TAB]"}, {",", "[COMMA]"}});
	}
	
	public static String removeQuotesFromExcelToken(String ori, String[][] theReplacements) {
		String str;
		int start, stop;
		
		str = replaceAllWith(ori, "\"\"", "'");
		while (str.contains("\"")) {
			start = str.indexOf("\"");
			stop = str.substring(start+1).indexOf("\"")+start+1;
			try {
				str = str.substring(0, start)+replaceAllWith(str.substring(start+1, stop), theReplacements)+str.substring(stop+1);
			} catch (Exception e) {
				System.err.println("Error - failed to parse: '"+ori+"'");
				e.printStackTrace();
				ext.waitForResponse();
			}
		}
		return str;
	}

	public static boolean containsAny(String str, String[] array) {
		for (int i = 0; i<array.length; i++) {
			if (str.contains(array[i])) {
				return true;
			}
		}
		return false;
	}

	public static boolean containsAll(String[] targets, String[] array) {
		for (int i = 0; i<targets.length; i++) {
			if (ext.indexOfStr(targets[i], array)==-1) {
				return false;
			}
		}
		return true;
	}

	public static Date parseDate(String str) {
		String[] line;

		line = str.trim().split("/");
		if (line.length!=3) {
			System.err.println("Error parsing date ('"+str+"'); assuming MM/DD/YYYY format");
		}

		return new GregorianCalendar(Integer.parseInt(line[2]), Integer.parseInt(line[0])-1, Integer.parseInt(line[1])).getTime();
	}

	public static int calcDays(Date start, Date end) {
		return (int)((end.getTime()-start.getTime())/1000/60/60/24);
	}

	public static String rootOf(String filename) {
		return rootOf(filename, true);
	}

	public static String rootOf(String filename, boolean trimDirectoryInfo) {
		if (trimDirectoryInfo && filename.indexOf("/")>=0) {
			filename = filename.substring(filename.lastIndexOf("/")+1);
		}
		if (trimDirectoryInfo && filename.indexOf("\\")>=0) {
			filename = filename.substring(filename.lastIndexOf("\\")+1);
		}
		if (filename.lastIndexOf(".") >= 0) {
			filename = filename.substring(0, filename.lastIndexOf(".")); 
		}
		return filename;
	}
	
	public static String rootRootOf(String filename) {
		return filename.substring(0, filename.indexOf("."));
	}
	
	public static String addToRoot(String filename, String addition) {
		return ext.rootOf(filename, false)+addition+(filename.lastIndexOf(".")>0?filename.substring(filename.lastIndexOf(".")):"");
	}
	
	public static String removeDirectoryInfo(String filename) {
		if (filename.indexOf("/")>=0) {
			filename = filename.substring(filename.lastIndexOf("/")+1);
		}
		if (filename.indexOf("\\")>=0) {
			filename = filename.substring(filename.lastIndexOf("\\")+1);
		}
		return filename;
	}
	
	public static String parseDirectoryOfFile(String filename) {
		return parseDirectoryOfFile(filename, false);
	}
	
	public static String parseDirectoryOfFile(String filename, boolean allowBlank) {
		String dir = "";
		
		filename = replaceAllWith(filename, "\\", "/");
		filename = replaceAllWith(filename, "//", "/");
		
		if (filename.indexOf("/")>=0) {
			dir = filename.substring(0, filename.lastIndexOf("/")+1);
		}
		
		if (dir.equals("") && !allowBlank) {
			dir = "./";
		}
		
		return dir;
	}

	public static void writeToAll(String str, PrintWriter[] writers) {
		for (int i = 0; i<writers.length; i++) {
			writers[i].println(str);
		}
	}

	public static String addCommas(int value) {
		String str = value+"";

		for (int i = str.length()-3; i>0; i = i-3) {
			str = str.substring(0, i)+","+str.substring(i);
		}

		return str;
	}
	
	public static String addCommas(long value) {
		String str = value+"";

		for (int i = str.length()-3; i>0; i = i-3) {
			str = str.substring(0, i)+","+str.substring(i);
		}

		return str;
	}
	
	public static void appendToFile(String text, String filename) {
		PrintWriter writer;
		
		try {
	        writer = new PrintWriter(new FileWriter(filename, true));
	        writer.println(text);
	        writer.close();
        } catch (Exception e) {
	        System.err.println("Error writing to "+filename);
	        e.printStackTrace();
        }
	}
	
	public static String verifyDirFormat(String str) {
//		if (str.startsWith("C:")) {
//			str = str.substring(2);
//		}
		str = replaceAllWith(str, "\\", "/");
		if (str.length() > 0 && !str.endsWith("/")) {
			str += "/";
		}
		
		if (str.endsWith("/./")) {
			str = str.substring(0, str.length()-2);
		}
		
		return str;
	}
	
	public static boolean parseBooleanArg(String arg) {
		return parseBooleanArg(arg, new Logger());
	}
	
	public static boolean parseBooleanArg(String arg, Logger log) {
		if (arg.split("=")[1].equalsIgnoreCase("true")) {
			return true;
		} else if (arg.split("=")[1].equalsIgnoreCase("false")) {
			return false;
		} else {
			log.reportError("Error - invalid "+arg.split("=")[0]+"= argument (expecting true/false): "+arg.split("=")[1]);
			System.exit(1);
			return false;
		}
	}

	public static String parseStringArg(String arg, String blankValue) {
		if (arg.split("=").length == 1) {
			return blankValue;
		} else if (arg.split("=")[1].equalsIgnoreCase("null")){
			return null;
		} else {
			return arg.substring(arg.indexOf("=")+1);
		}
	}

	public static int parseIntArg(String arg) {
		try {
			return Integer.parseInt(arg.split("=")[1]);
		} catch (NumberFormatException nfe) {
			System.err.println("Error - invalid "+arg.split("=")[0]+"= argument: "+arg.split("=")[1]);
			System.exit(1);
			return Integer.MIN_VALUE;
		}
	}

	public static byte parseByteArg(String arg) {
		try {
			return Byte.parseByte(arg.split("=")[1]);
		} catch (NumberFormatException nfe) {
			System.err.println("Error - invalid "+arg.split("=")[0]+"= argument: "+arg.split("=")[1]);
			System.exit(1);
			return Byte.MIN_VALUE;
		}
	}

	public static double parseDoubleArg(String arg) {
		try {
			return Double.parseDouble(arg.split("=")[1]);
		} catch (NumberFormatException nfe) {
			System.err.println("Error - invalid "+arg.split("=")[0]+"= argument: "+arg.split("=")[1]);
			System.exit(1);
			return Double.NaN;
		}
	}
	
	public static boolean isMissingValue(String str) {
		for (int i = 0; i < MISSING_VALUES.length; i++) {
			if (str.equalsIgnoreCase(MISSING_VALUES[i])) {
				return true;
			}
		}
		return false;
	}
	
	public static boolean isValidDouble(String str) {
		if (isMissingValue(str)) {
			return false;
		}
		
		try {
			Double.parseDouble(str);
		} catch (NumberFormatException nfe) {
			return false;
		}

		return true;
	}

	public static boolean isValidInteger(String str) {
		if (isMissingValue(str)) {
			return false;
		}
		
		try {
			Integer.parseInt(str);
		} catch (NumberFormatException nfe) {
			return false;
		}

		return true;
	}
	
	public static float parseFloatArg(String arg) {
		try {
			return Float.parseFloat(arg.split("=")[1]);
		} catch (NumberFormatException nfe) {
			System.err.println("Error - invalid "+arg.split("=")[0]+"= argument: "+arg.split("=")[1]);
			System.exit(1);
			return Float.NaN;
		}
	}
	
	public static String insertNumbers(String pattern, int num) {
		return ext.replaceAllWith(ext.replaceAllWith(pattern, "##", chrome(num)), "#", num+"");
	}

	public static String insertNumbers(String pattern, int num, int numDigits) {
		if (numDigits == -1) {
			return ext.replaceAllWith(pattern, "#", num+"");
		} else {
			return ext.replaceAllWith(pattern, "#", ext.formNum(num, numDigits));
		}
	}
	
	public static String right(String str, int numChar) {
		return str.substring(str.length()-numChar, str.length());
	}
	
	public static String link(String linkRoot, String dir) {
//		int count;
		
		linkRoot = verifyDirFormat(linkRoot);

//		count = 0;
		while (dir.startsWith("../")) {
			linkRoot = linkRoot.substring(0, linkRoot.lastIndexOf("/"));
			dir = dir.substring(3);
//			count++;
		}
		
		return linkRoot.substring(0, linkRoot.lastIndexOf("/")+1)+dir;
	}

	public static int countInstancesOf(String str, String pattern) {
		int count;
		String trav;
		
		trav = str;
		count = 0;
		while (trav.contains(pattern)) {
			trav = trav.substring(trav.indexOf(pattern)+pattern.length());
			count++;
		}
		
		return count;		
	}
	
	public static String determineDelimiter(String str) {
		if (str.contains("\t")) {
			return "\t";
		} else if (countInstancesOf(str, ",") > countInstancesOf(str, " ")) {
			return ",";
		} else {
			return "[\\s]+";
		}
	}

	public static String removeAndSimplifyQuotes(String str, Logger log) {
		if (str.startsWith("\"")) {
			if (!str.endsWith("\"")) {
				log.reportError("Error - improperly formed quotes in comma delimited token:");
				log.reportError(str);
			}
			str = str.substring(1, str.length()-1);
		}
		return ext.replaceAllWith(str, "\"\"", "\"");
	}
	
	public static String[] splitCommasIntelligently(String str, boolean removeAndSimplifyQuotes, Logger log) {
		int numCommas;
		int startIndex;
		Vector<String> v;
		boolean insideQuotes;

		insideQuotes = false;
		v = new Vector<String>();
		startIndex = numCommas = 0;
		for (int i = 0; i < str.length(); i++) {
			if (str.charAt(i) == ',') {
				if (numCommas % 2 == 0) {
					v.add(removeAndSimplifyQuotes?removeAndSimplifyQuotes(str.substring(startIndex, i), log):str.substring(startIndex, i));
//					System.out.println(v.elementAt(v.size()-1));
					startIndex = i+1;
					if (insideQuotes && str.charAt(i-1) != '\"') {
						log.reportError("Error - improperly formed quotes in comma delimited string:");
						log.reportError(str);
					}
					insideQuotes = false;
				}
			}
			if (str.charAt(i) == '\"') {
				numCommas++;
				insideQuotes = true;
			}
		}
		v.add(removeAndSimplifyQuotes?removeAndSimplifyQuotes(str.substring(startIndex), log):str.substring(startIndex));
		
		return Array.toStringArray(v);
	}
	
	public static String reportMemoryUsage() {
		MemoryUsage muHeap; // , muNonHeap;
		
		muHeap = ManagementFactory.getMemoryMXBean().getHeapMemoryUsage();
//		muNonHeap = ManagementFactory.getMemoryMXBean().getNonHeapMemoryUsage();
		
//		return "Total: " + ext.prettyUpSize(Runtime.getRuntime().totalMemory(), 1) + 
//				" Free: " + ext.prettyUpSize(Runtime.getRuntime().freeMemory(),1) + 
//				" Max: " + ext.prettyUpSize(Runtime.getRuntime().maxMemory(),1) +
//				" HeapComit: " + ext.prettyUpSize(muHeap.getCommitted(),1) +
//				" HeapInit: " + ext.prettyUpSize(muHeap.getInit(),1) +
//				" HeapMax: " + ext.prettyUpSize(muHeap.getMax(),1) +
//				" HeapUsed: " + ext.prettyUpSize(muHeap.getUsed(),1) +
//				" NonHeapComit: " + ext.prettyUpSize(muNonHeap.getCommitted(),1) +
//				" NonHeapInit: " + ext.prettyUpSize(muNonHeap.getInit(),1) +
//				" NonHeapMax: " + ext.prettyUpSize(muNonHeap.getMax(),1) +
//				" NonHeapUsed: " + ext.prettyUpSize(muNonHeap.getUsed(),1)
//				;

		return "Used: " + ext.prettyUpSize(muHeap.getUsed(),1) +
				"\tInit: " + ext.prettyUpSize(muHeap.getInit(),1) +
				"\tCommit: " + ext.prettyUpSize(muHeap.getCommitted(),1) +
				"\tMax: " + ext.prettyUpSize(Runtime.getRuntime().maxMemory(),1) +
				"\tFree: " + ext.prettyUpSize(Runtime.getRuntime().freeMemory(),1) 
				;
	}
	
	public static String pwd() {
		return ext.verifyDirFormat(new File("./").getAbsolutePath());
	}

	public static String replaceTilde(String path) {
		return path.replace("~", System.getProperty("user.home"));
	}

	public static String[] splitLine(String str, String delimiter, Logger log) {
		if (delimiter.equals("\t")) {
			return str.split(delimiter, -1);
		} else if (delimiter.equals(",")) {
			return ext.splitCommasIntelligently(str, true, log);
		} else {
			return str.trim().split(delimiter);
		}
	}

	public static void waitForResponse() {
		waitForResponse("Press ENTER to continue");
	}
	
	public static void waitForResponse(String message) {
		if (System.getProperty("os.name").startsWith("Windows")) {
			System.out.println(message);
			try {
				new BufferedReader(new InputStreamReader(System.in)).readLine();
			} catch (IOException ioe) {}
		}
	}
	
	public static void main(String[] args) {
		String temp;
		Logger log;
		
		log = new Logger();
		try {
			System.out.println("1\t2\t3\t4\t5\t6\t7");
			temp = "try,this,out";
			System.out.println(temp);
			System.out.println(Array.toStr(splitCommasIntelligently(temp, false, log)));
			System.out.println();
			temp = "\"is this, harder\",o,r,\"not so ,much\"";
			System.out.println(temp);
			System.out.println("1\t2\t3\t4\t5\t6\t7");
			System.out.println(Array.toStr(splitCommasIntelligently(temp, false, log)));
			System.out.println(Array.toStr(splitCommasIntelligently(temp, true, log)));
			System.out.println();
			temp = "\"perhaps this is \"hard,est\" of all\",or,\"\"\"maybe not\"\"\"";
			temp = "\"perhaps this is \"\"hard,est\"\" of all\",or,\"\"\"maybe not\"\"\"";			
			System.out.println(temp);
			System.out.println("1\t2\t3\t4\t5\t6\t7");
			System.out.println(Array.toStr(splitCommasIntelligently(temp, false, log)));
			System.out.println(Array.toStr(splitCommasIntelligently(temp, true, log)));
			System.out.println();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
