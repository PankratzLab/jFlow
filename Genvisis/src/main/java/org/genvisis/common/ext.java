package org.genvisis.common;

import java.awt.Toolkit;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.StringSelection;
import java.awt.datatransfer.Transferable;
import java.awt.datatransfer.UnsupportedFlavorException;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryUsage;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collection;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Vector;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;

import com.google.common.collect.ImmutableSet;
import com.google.common.primitives.Ints;

public class ext {
	public static final Set<Character> UNSAFE_CHARS = ImmutableSet.of(' ', '?', '\'', '/', '\\', '<',
																																		'>',
																																		'|', ':', '*', '\'');
	public static final Set<Character> UNSAFE_CHARS_STRICT;
	static {
		ImmutableSet.Builder<Character> unsafeCharsStrictBuilder = ImmutableSet.builder();
		unsafeCharsStrictBuilder.addAll(UNSAFE_CHARS);
		unsafeCharsStrictBuilder.add('(', ')', '&');
		UNSAFE_CHARS_STRICT = unsafeCharsStrictBuilder.build();
	}

	public static final Character UNSAFE_CHAR_REPLACEMENT = '_';
	public static final String[] MISSING_VALUES = {"", ".", "NA", "NaN", "x", "#N/A", "--", "-"};
	public static final String[][] META_REPLACEMENTS = {{"{Tab}", "\t"}, {"{Space}", " "},
																											{"{!}", "!"}, {"{#}", "#"}, {"{+}", "+"}};
	public static final String[] COMMON_IDS = {"id", "IID", "IndID", "gwas_id"};
	public static final String REGEX_TO_SPLIT_SPACES_NOT_IN_QUOTES = PSF.Regex.regexSplitPreserveQuoted(" ");

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
					str = str.substring(0, str.indexOf(from)) + to
								+ str.substring(str.indexOf(from) + from.length());
				} catch (Exception e) {
					System.err.println("Error with " + from + " to " + to + ": " + original);
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
				str = str.substring(0, str.indexOf(from)) + to
							+ str.substring(str.indexOf(from) + from.length());
			}
		}
		return str;
	}

	public static String replaceQuotesWithSlashQuotes(String original) {
		char[] array;
		String str;

		array = original.toCharArray();

		str = "";
		for (char element : array) {
			if (element == '"') {
				str += "\\";
			}
			str += element;
		}

		return str;
	}

	public static String[][] processMetaCharacters(String[][] replaces) {
		String[][] metaReplaces;

		metaReplaces = new String[replaces.length][];
		for (int i = 0; i < replaces.length; i++) {
			metaReplaces[i] = new String[] {replaces[i][0],
																			replaceAllWith(replaces[i][1], META_REPLACEMENTS)};
		}

		return metaReplaces;
	}

	public static String replaceAllWith(String str, String[][] replaces) {
		for (String[] replace : replaces) {
			str = ext.replaceAllWith(str, replace[0], replace[1]);
		}
		return str;
	}

	public static String chrome(int chr) {
		return (chr < 10) ? "0" + chr : "" + chr;
	}

	public static int gender(String gender) {
		if (gender.toLowerCase().equals("m") || gender.equals("1")) {
			return 1;
		} else if (gender.toLowerCase().equals("f") || gender.equals("2")) {
			return 2;
		} else {
			return 0;
		}
	}

	public static String prettyUpDistance(int dist, int numSigFigs) {
		String str = dist + "";
		if (str.length() > 6) {
			return formDeci((double) Integer.parseInt(str) / 1000000, numSigFigs, true) + " Mb";
		}
		if (str.length() > 3) {
			return formDeci((double) Integer.parseInt(str) / 1000, numSigFigs, true) + " kb";
		}
		return str;
	}

	public static String prettyUpSize(long size, int numSigFigs) {
		String str = size + "";
		if (str.length() > 12) {
			return formDeci((double) size / 1048576000000l, numSigFigs, true) + " TB";
		}
		if (str.length() > 9) {
			return formDeci((double) size / 1048576000, numSigFigs, true) + " GB";
		}
		if (str.length() > 6) {
			return formDeci((double) size / 1048576, numSigFigs, true) + " MB";
		}
		if (str.length() > 3) {
			return formDeci((double) size / 1024, numSigFigs, true) + " KB";
		}
		return str;
	}

	public static int indexOfInt(int target, int[] array) {
		int index = -1;

		for (int i = 0; i < array.length; i++) {
			if (array[i] == target) {
				if (index == -1) {
					index = i;
				} else {
					return -2;
				}
			}
		}
		return index;
	}

	public static int indexOfAnyStr(String target, String[][] array, boolean caseSensitive) {
		for (int i = 0; i < array.length; i++) {
			for (int j = 0; j < array[i].length; j++) {
				if (caseSensitive ? array[i][j].equals(target) : array[i][j].equalsIgnoreCase(target)) {
					return i;
				}
			}
		}
		return -1;
	}

	/**
	 * @see #indexOfStr(String, String[], boolean, boolean, Logger, boolean)
	 */
	public static int indexOfStr(String target, String[] superset) {
		return indexOfStr(target, superset, true, true, new Logger(), false);
	}

	/**
	 * @see #indexOfStr(String, String[], boolean, boolean, Logger, boolean)
	 */
	public static int indexOfStr(String target, String[] superset, boolean caseSensitive,
															 boolean exactMatch) {
		return indexOfStr(target, superset, caseSensitive, exactMatch, new Logger(), false);
	}

	/**
	 * Identify the index of the first string in a superset matching a target string
	 *
	 * @param target String of interest
	 * @param superset Space to search for target string
	 * @param caseSensitive Whether or not to consider case when comparing subset and superset strings
	 * @param exactMatch Whether or not we require a perfect match between a target and superset
	 * @param verbose If true, report any mismatches (default: {@code true})
	 * @param log Logger instance to use when reporting
	 * @return First index of the target string in the superset, or -1 if not found
	 *
	 * @see #indexFactors
	 * @see #indicesOfStr
	 */
	public static int indexOfStr(String target, String[] superset, boolean caseSensitive,
															 boolean exactMatch, Logger log, boolean verbose) {
		// FIXME if not verbose we don't actually use the fact that there could be more than one match,
		// which could potentially shorten search time.
		int[] indices = indicesOfStr(target, superset, caseSensitive, exactMatch);

		if (indices.length == 0) {
			if (verbose) {
				String error = "Error - '" + target + "' was not found in array";
				if (log == null) {
					System.err.println(error);
				} else {
					log.reportError(error);
				}
			}
			return -1;
		}

		if (indices.length > 1 && verbose) {
			String warning = "Warning - '" + target + "' was found more than once in the array";
			if (log == null) {
				System.err.println(warning);
			} else {
				log.reportError(warning);
			}
		}

		return indices[0];
	}

	/**
	 * Identify the indices of all strings in a superset matching a target string
	 *
	 * @param target String of interest
	 * @param superset Space to search for target string
	 * @param caseSensitive Whether or not to consider case when comparing subset and superset strings
	 * @param exactMatch Whether or not we require a perfect match between a target and superset
	 * @return Array of all indices in the superset that match the target
	 *
	 * @see #indexOfStr
	 */
	public static int[] indicesOfStr(String target, String[] superset, boolean caseSensitive,
																	 boolean exactMatch) {
		int[] indices = new int[superset.length];
		int numMatches = 0;

		if (!caseSensitive) {
			target = target.toLowerCase();
		}

		for (int i = 0; i < superset.length; i++) {
			String toCheck = superset[i];
			if (!caseSensitive) {
				toCheck = toCheck.toLowerCase();
			}

			if ((exactMatch && target.equals(toCheck))
					|| (!exactMatch && (target.contains(toCheck) || toCheck.contains(target)))) {
				indices[numMatches++] = i;
			}
		}

		return Arrays.copyOf(indices, numMatches);
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
				iv.add(runningIndex + index);
				runningIndex += index + target.length();
				str = str.substring(index + target.length());
			}
		}

		return Ints.toArray(iv);
	}

	public static int indexOfChar(char target, char[] array) {
		for (int i = 0; i < array.length; i++) {
			if (array[i] == target) {
				return i;
			}
		}

		return -1;
	}

	public static boolean startsWithOneOf(String source, String... prefixes) {
		for (String prefix : prefixes) {
			if (source.startsWith(prefix)) {
				return true;
			}
		}
		return false;
	}

	public static boolean startsWithOneOf(String source, Iterable<String> prefixes) {
		for (String prefix : prefixes) {
			if (source.startsWith(prefix)) {
				return true;
			}
		}
		return false;
	}

	public static int indexOfStartsWith(String target, String[] array, boolean reverseNotForward) {
		return indexOfStartsWith(target, array, reverseNotForward, true);
	}

	/**
	 * 
	 * @param target
	 * @param array
	 * @param reverseNotForward true to match on target starting with array entry, false to match on
	 *        array entry starting with target
	 * @param caseSensitive true for case sensitive match, false to ignore case
	 * @return first index in array that matches criteria or -1 for no match
	 */
	public static int indexOfStartsWith(String target, String[] array, boolean reverseNotForward,
																			boolean caseSensitive) {
		String cleanTarget = maybeLowercase(target, !caseSensitive);
		if (reverseNotForward) {
			for (int i = 0; i < array.length; i++) {
				if (cleanTarget.startsWith(maybeLowercase(array[i], !caseSensitive))) {
					return i;
				}
			}
		} else {
			for (int i = 0; i < array.length; i++) {
				if (maybeLowercase(array[i], !caseSensitive).startsWith(cleanTarget)) {
					return i;
				}
			}
		}
		return -1;
	}

	private static String maybeLowercase(String source, boolean lowercase) {
		return lowercase ? source.toLowerCase() : source;
	}

	public static int indexOfEndsWith(String target, String[] array, boolean reverseNotForward) {
		if (reverseNotForward) {
			for (int i = 0; i < array.length; i++) {
				if (target.endsWith(array[i])) {
					return i;
				}
			}
		} else {
			for (int i = 0; i < array.length; i++) {
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
	 * @param str String to be parsed
	 * @param operators String of the possible operators
	 * @return a matrix of String[3][] with an array each of the operators, the indices of those
	 *         operators, and the tokens
	 */
	public static String[][] getOperatorsOperatorIndicesAndSplit(String str, String operators) {
		Vector<String> operatorsPresent = new Vector<String>();
		Vector<String> operatorIndices = new Vector<String>();
		Vector<String> tokens = new Vector<String>();
		Vector<String> ops = new Vector<String>();
		int count = 0;
		String hist = "";

		for (int i = 0; i < operators.length(); i++) {
			ops.add(operators.charAt(i) + "");
		}

		while (count < str.length()) {
			if (ops.contains(str.charAt(count) + "")) {
				operatorsPresent.add(str.charAt(count) + "");
				operatorIndices.add(count + "");
				tokens.add(hist);
				hist = "";
			} else {
				hist += str.charAt(count);
			}
			count++;
		}
		tokens.add(hist);

		return new String[][] {ArrayUtils.toStringArray(operatorsPresent),
													 ArrayUtils.toStringArray(operatorIndices),
													 ArrayUtils.toStringArray(tokens)};
	}

	public static boolean eqArrays(Object[] arr1, Object[] arr2) {
		if (arr1.length != arr2.length) {
			return false;
		}
		for (int i = 0; i < arr1.length; i++) {
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
	 * @param str name
	 * @param size number of characters
	 * @param after add white space after the string (instead of before)
	 * @return processed string
	 */
	public static String formStr(String str, int size, boolean after) {
		String spaces = "";
		for (int i = 0; i < (size - str.length()); i++) {
			spaces += " ";
		}
		if (after) {
			return str + spaces;
		} else {
			return spaces + str;
		}
	}

	public static String formNum(int num, int size) {
		return formNum(num + "", size);
	}

	public static String formNum(String str, int size) {
		String newNum = "";
		for (int i = 0; i < (size - str.length()); i++) {
			newNum += "0";
		}
		return newNum + str;
	}

	public static String formNum(int num1, int num2, int size1, int size2) {
		return formNum(num1 + "", num2 + "", size1, size2);
	}

	public static String formNum(String str1, String str2, int size1, int size2) {
		String newNum = "1";
		for (int i = 0; i < (size1 - str1.length()); i++) {
			newNum += "0";
		}
		newNum += str1;
		for (int i = 0; i < (size2 - str2.length()); i++) {
			newNum += "0";
		}
		return newNum + str2;
	}

	public static String formDeci(double num, int numPlaces) {
		return formDeci(num, numPlaces, false);
	}

	public static String formDeci(double num, int numPlaces, boolean forceDigits) {
		return formDeci(num, forceDigits ? numPlaces : 0, numPlaces);
	}

	public static String formDeci(double num, int minSigFigs, int maxSigFigs) {
		String theFormat = (maxSigFigs > 0 ? "0." : "0");
		String result;

		if (Double.isNaN(num)) {
			return Double.toString(num);
		}

		for (int i = 0; i < minSigFigs; i++) {
			theFormat += "0";
		}

		for (int i = minSigFigs; i < maxSigFigs; i++) {
			theFormat += "#";
		}

		result = new DecimalFormat(theFormat).format(num);
		if (ext.replaceAllWith(result, new String[][] {{"0", ""}, {".", ""}}).equals("-")) {
			result = result.substring(1);
		}
		return result;
	}

	public static String formDeci(float num, int numPlaces) {
		return formDeci(num, numPlaces, false);
	}

	public static String formDeci(float num, int numPlaces, boolean forceDigits) {
		return formDeci(num, forceDigits ? numPlaces : 0, numPlaces);
	}

	public static String formDeci(float num, int minSigFigs, int maxSigFigs) {
		String theFormat = (maxSigFigs > 0 ? "0." : "0");
		String result;

		if (Double.isNaN(num)) {
			return Float.toString(num);
		}

		for (int i = 0; i < minSigFigs; i++) {
			theFormat += "0";
		}

		for (int i = minSigFigs; i < maxSigFigs; i++) {
			theFormat += "#";
		}

		result = new DecimalFormat(theFormat).format(num);
		if (ext.replaceAllWith(result, new String[][] {{"0", ""}, {".", ""}}).equals("-")) {
			result = result.substring(1);
		}
		return result;
	}

	public static String formPercent(double num, int numSigFigs) {
		String theFormat = (numSigFigs > 0 ? "0." : "0");

		if (Double.isNaN(num)) {
			return num + "";
		}

		for (int i = 0; i < numSigFigs; i++) {
			theFormat += "0";
		}

		return new DecimalFormat(theFormat).format(num * 100) + "%";
	}

	public static String formSciNot(double num, int numPlaces, boolean forceDigits) {
		DecimalFormat f;
		String theFormat = (numPlaces > 0 ? "0." : "0");

		for (int i = 0; i < numPlaces; i++) {
			theFormat += (forceDigits ? "0" : "#");
		}

		f = new DecimalFormat(theFormat + "E0");

		return f.format(num);
	}

	public static String prettyP(double d) {
		return prettyP(d + "", 2, 5, 1, false);
	}

	public static String prettyP(String nummer) {
		return prettyP(nummer, 2, 5, 1, false);
	}

	public static String prettyP(double d, int sigfigs, int sigfigWhenConverting,
															 int sigfigWhenConverted, boolean useE) {
		return prettyP(d + "", sigfigs, sigfigWhenConverting, sigfigWhenConverted, useE);
	}

	/**
	 * Converts a String (which is usually a p-value) to a double and then back to a String using
	 * scientific notation if beyond a certain magnitude
	 *
	 * @param nummer the original number in String format
	 * @param sigfigs the number of significant digits to use if there is no conversion (usually 2)
	 * @param sigfigBeforeConverting the number of digits (i.e., the number of zeros after the decimal
	 *        plus one) required before switching to scientific notation (usually 4 or 5)
	 * @param sigfigWhenConverted the number of significant digits to use when using scientific
	 *        notation (if you use 3, then that means there will be the one digit before the decimal
	 *        and two after)
	 * @param useE use E-0# instead of x 10^-#
	 * @return the pretty version of the number
	 */
	public static String prettyP(String nummer, int sigfigs, int sigfigBeforeConverting,
															 int sigfigWhenConverted, boolean useE) {
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
			return "'" + nummer + "' is an ugly number (failed prettyP)";
		}

		if (d > (Math.pow(0.1, sigfigs) - Math.pow(0.1, sigfigs + 1) / 2) || Double.isNaN(d)) {
			return formDeci(d, sigfigs, true);
		}

		if (d < 0) {
			return "0";
		}

		temp = "0";
		for (int i = 1; i < sigfigWhenConverted; i++) {
			if (i == 1) {
				temp += ".";
			}
			temp += "0";
		}
		temp += "E0";
		str = new DecimalFormat(temp).format(d);
		sigdigits = str.split("E")[0];
		power = Integer.parseInt(str.split("E")[1]) * -1;

		if (power < sigfigBeforeConverting) {
			str = "0.";
			for (int i = 0; i < power - 1; i++) {
				str += "0";
			}
			str += ext.replaceAllWith(sigdigits, ".", "");
			return str;
		} else {
			return useE ? str : sigdigits + "x10-" + power;
		}
	}

	public static String getDate(Date date, String delimiter) {
		GregorianCalendar gc = new GregorianCalendar();
		gc.setTime(date);
		return gc.get(Calendar.YEAR) + delimiter + chrome((gc.get(Calendar.MONTH) + 1)) + delimiter
					 + chrome(gc.get(Calendar.DAY_OF_MONTH));
	}

	public static String getDate(String delimiter) {
		return getDate(new Date(), delimiter);
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
		return formNum(gc.get(Calendar.HOUR) == 0 ? 12 : gc.get(Calendar.HOUR), 2) + ":"
					 + formNum(gc.get(Calendar.MINUTE), 2) + ":" + formNum(gc.get(Calendar.SECOND), 2)
					 + (gc.get(Calendar.AM_PM) == 0 ? "AM" : "PM");
	}

	public static void timestamp() {
		System.out.println(getDate() + "\t" + getTime());
	}

	public static String getTimestampForFilename() {
		return replaceAllWith((getDate() + "_" + getTime()), ":", "_");
	}

	public static String formatTimeElapsed(long timeElapsed, TimeUnit elapsedUnit) {
		String elapsed;
		switch (elapsedUnit) {
			case DAYS:
				elapsed = formatMinutesElapsed(timeElapsed * 60 * 24);
				break;
			case HOURS:
				elapsed = formatMinutesElapsed(timeElapsed * 60);
				break;
			case MICROSECONDS:
				elapsed = TimeUnit.MICROSECONDS.toMillis(timeElapsed) > 9
																																	? formatMillisElapsed(timeElapsed
																																												/ 1000)
																																	: (timeElapsed + " µs");
				break;
			case MILLISECONDS:
				elapsed = formatMillisElapsed(timeElapsed);
				break;
			case MINUTES:
				elapsed = formatMinutesElapsed(timeElapsed);
				break;
			case NANOSECONDS:
				elapsed = TimeUnit.NANOSECONDS.toMillis(timeElapsed) > 9
																																 ? formatMillisElapsed(timeElapsed
																																											 / 1000000)
																																 : (timeElapsed + " ns");
				break;
			case SECONDS:
				elapsed = formatMillisElapsed(timeElapsed * 1000);
				break;
			default:
				elapsed = timeElapsed + " (unknown TimeUnit)";
				break;
		}
		return elapsed;
	}

	/**
	 * 
	 * @param startTimeNanos Start time in nanos (Usually from <code>System.nanoTime()</code>)
	 * @return Formatted string of the most appropriate unit (e.g. "1 hr 22 mins")
	 */
	public static String getTimeElapsedNanos(long startTimeNanos) {
		long timeElapsed = System.nanoTime() - startTimeNanos;
		return formatTimeElapsed(timeElapsed, TimeUnit.NANOSECONDS);
	}

	public static String getTimeElapsed(long startTime) {
		long timeElapsed;
		timeElapsed = new Date().getTime() - startTime;

		return formatMillisElapsed(timeElapsed);
	}

	public static String formatMillisElapsed(long millisElapsed) {
		long timeElapsed = millisElapsed;

		if (timeElapsed < 10000) {
			return timeElapsed + " ms";
		}

		timeElapsed /= 1000;
		if (timeElapsed < 60) {
			return timeElapsed + " sec";
		}

		if (timeElapsed < 600) {
			return timeElapsed / 60 + " min " + (timeElapsed - (timeElapsed / 60) * 60) + " sec";
		}

		timeElapsed /= 60;

		return formatMinutesElapsed(timeElapsed);
	}

	public static String formatMinutesElapsed(long minsElapsed) {
		long timeElapsed = minsElapsed;
		if (timeElapsed < 60) {
			return timeElapsed > 1 ? timeElapsed + " mins" : timeElapsed + " min";
		}

		if (timeElapsed < 24 * 60) {
			int hrs = (int) (timeElapsed / 60);
			return hrs + " hr" + (hrs > 1 ? "s " : " ") + (timeElapsed - (hrs * 60))
						 + " min";
		}

		timeElapsed /= 60;

		int days = (int) (timeElapsed / 24);
		return days + " day" + (days > 1 ? "s " : " ") + (timeElapsed - (days * 24)) + " hrs";
	}

	public static double getTimeSince(long startTime, char returnType) {
		double timeElapsed;

		timeElapsed = new Date().getTime() - startTime;
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

		System.err.println("Error - '"
											 + returnType
											 + "' is an invalid return type: only S(econds), M(inutes), and H(ours) are allowed");

		return -1;
	}

	public static double divide(int num1, int num2) {
		if (num2 == 0) {
			return 0;
		}
		return (double) num1 / (double) num2;
	}

	public static String replaceWithLinuxSafeCharacters(String str) {
		return replaceWithLinuxSafeCharacters(str, true);
	}

	public static String replaceWithLinuxSafeCharacters(String str, boolean extraSafe) {
		String cleanString = str;
		Set<Character> unsafeChars = extraSafe ? UNSAFE_CHARS_STRICT : UNSAFE_CHARS;
		for (char unsafeChar : unsafeChars) {
			cleanString = cleanString.replace(unsafeChar, UNSAFE_CHAR_REPLACEMENT);
		}
		return cleanString;
	}

	public static String replaceDirectoryCharsWithUnderscore(String filename, int depth) {
		int index;

		filename = replaceAllWith(filename, "\\", "/");
		filename = replaceAllWith(filename, "//", "/");

		for (int i = 0; i < depth; i++) {
			index = filename.lastIndexOf("/");
			if (index >= 0) {
				filename = filename.substring(0, index) + "_" + filename.substring(index + 1);
			}
		}
		index = filename.lastIndexOf("/");
		if (index >= 0) {
			filename = filename.substring(index + 1);
		}

		return filename;
	}

	// public static int[] indexFactors(String trait, Vector included, String[] phenoNames) {
	// return indexFactors(trait, included, phenoNames, false);
	// }
	//
	// public static int[] indexFactors(String trait, Vector included, String[] phenoNames, boolean
	// possibleEmptySet) {
	// int[] indices;
	// int first, last;
	//
	// if (included==null) {
	// included = new Vector();
	// } else if (included.size()==0&&!possibleEmptySet) {
	// throw new RuntimeException("No independent factors listed; flag true if legitimate");
	// }
	//
	// if (included.contains(trait)) {
	// System.err.println("Warning - "+trait+" is listed as both a factor and the dependent variable;
	// removed as a factor");
	// included.remove(trait);
	// }
	//
	// indices = new int[included.size()+1];
	// for (int i = 0; i<indices.length; i++) {
	// indices[i] = -1;
	// }
	//
	// for (int i = 0; i<phenoNames.length; i++) {
	// if (included.contains(phenoNames[i])) {
	// first = included.indexOf(phenoNames[i]);
	// last = included.lastIndexOf(phenoNames[i]);
	// for (int j = first; j<=last; j++) {
	// if (phenoNames[i].equals(included.elementAt(j))) {
	// indices[j+1] = i;
	// }
	// }
	// } else if (phenoNames[i].equals(trait)) {
	// indices[0] = i;
	// }
	// }
	//
	// for (int i = 0; i<indices.length; i++) {
	// if (indices[i]==-1) {
	// throw new RuntimeException("Error - variable name '"+(i==0?trait:included.elementAt(i-1))+"'
	// was not found in the database");
	// }
	// }
	//
	// return indices;
	// }

	/**
	 * Method very similar to other index factors, but uses a hash based approach. When subset and
	 * superset are very large such as extracting marker subsets (i.e. 50K and 2500K), this might be a
	 * bit faster
	 *
	 * @param subset : what to extract
	 * @param superset : extract indices from
	 * @param casesensitive : case sensitive match
	 * @param log
	 * @param verbose : report duplicate matches, and no matches
	 * @return indices of subset in superset
	 *
	 */
	public static int[] indexLargeFactors(String[] subset, String[] superset, boolean casesensitive,
																				Logger log, boolean verbose) {
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
		return indices;
	}

	/**
	 * @see #indexFactors(String[], String[], boolean, Logger, boolean)
	 */
	public static int[] indexFactors(String[] subset, String[] superset, boolean caseSensitive) {
		return indexFactors(subset, superset, caseSensitive, new Logger(), true);
	}

	/**
	 * For each string in a subset, identify the index of the first matching string in a superset
	 *
	 * @param subset Query strings to identify indices in superset
	 * @param superset Space to search for subset strings
	 * @param caseSensitive Whether or not to consider case when comparing subset and superset strings
	 * @param log Logger instance to use when reporting
	 * @param verbose If true, report any mismatches (default: {@code true})
	 * @return Parallel array of subset indices in the superset (-1 if not found)
	 *
	 * @see #indexFactors(String[][], String[], boolean, boolean, boolean, boolean, Logger)
	 * @see #indexOfStr
	 */
	public static int[] indexFactors(String[] subset, String[] superset, boolean caseSensitive,
																	 Logger log, boolean verbose) {
		if (subset.length == 1) {
			// If we're looking up a single string, use indexOfStr to avoid computation cost of building a
			// map.
			return new int[] {indexOfStr(subset[0], superset, caseSensitive, true, log, verbose)};
		}

		int[] indices = new int[subset.length];
		boolean err = false;

		// Map the superset strings to their indices in the array
		Map<String, Integer> supersetMap = makeIndexMap(superset, caseSensitive, verbose, log);

		// Loop through our query strings and look up their indices in the map
		for (int i = 0; i < subset.length; i++) {
			String s = subset[i];
			if (!caseSensitive) {
				s = s.toLowerCase();
			}
			if (supersetMap.containsKey(s)) {
				indices[i] = supersetMap.get(s);
			} else {
				err = true;
				indices[i] = -1;
				if (verbose) {
					log.reportError("Error - no factor was named '" + subset[i] + "'");
				}
			}
		}

		return indices;
	}

	// public static boolean[] indexPresent(String[] subset, String[] superset, boolean casesensitive)
	// {
	// boolean[] attendance = new boolean[subset.length];
	//
	// for (int i = 0; i<subset.length; i++) {
	// attendance[i] = false;
	// for (int j = 0; j<superset.length; j++) {
	// if (casesensitive?subset[i].equals(superset[j]):subset[i].equalsIgnoreCase(superset[j])) {
	// attendance[i] = true;
	// }
	// }
	// }
	//
	// return attendance;
	// }
	//
	/**
	 * @see #indexFactors(String[][], String[], boolean, boolean, boolean, boolean, Logger)
	 */
	public static int[] indexFactors(String[][] targetsWithAlts, String[] superset,
																	 boolean caseSensitive, boolean exactMatch, boolean verbose) {
		return indexFactors(targetsWithAlts, superset, caseSensitive, exactMatch, verbose,
												new Logger());
	}

	/**
	 * @see #indexFactors(String[][], String[], boolean, boolean, boolean, boolean, Logger)
	 */
	public static int[] indexFactors(String[][] targetsWithAlts, String[] superset,
																	 boolean caseSensitive, boolean exactMatch, boolean verbose,
																	 Logger log) {
		return indexFactors(targetsWithAlts, superset, false, caseSensitive, exactMatch, verbose, log);
	}

	/**
	 * For each set of aliases in a subset, identify the index of the best matching string in a
	 * superset
	 *
	 * @param targetsWithAlts Array of query strings, each having one or more possible alternatives
	 * @param superset Space to search for subset strings
	 * @param preferFirstInTargetsOverFirstInSuperset For each target, whether we prefer lower indices
	 *        in the alts or the superset array
	 * @param caseSensitive Whether or not to consider case when comparing subset and superset strings
	 * @param exactMatch Whether or not we require a perfect match between a target and superset
	 *        string, or if it is sufficient for one to be a substring of the other
	 * @param verbose If true, report any mismatches (default: {@code true})
	 * @param log Logger instance to use when reporting
	 * @return Parallel array of the best target indices in the superset (-1 indicates target was not
	 *         found)
	 *
	 * @see #indexFactors(String[], String[], boolean, Logger, boolean)
	 * @see #indexOfStr
	 */
	public static int[] indexFactors(String[][] targetsWithAlts, String[] superset,
																	 boolean preferFirstInTargetsOverFirstInSuperset,
																	 boolean caseSensitive, boolean exactMatch, boolean verbose,
																	 Logger log) {
		int[] finalIndices = new int[targetsWithAlts.length];
		boolean err = false;

		// Exact match requires a different algorithm, since there is not an efficient way to create an
		// "inexact" map without generating a complete dictionary.
		if (!exactMatch) {
			if (!caseSensitive) {
				// Since we will be iterating through the superset once for each alt string, we do a single
				// pass here to convert the superset to lowercase. We use a new array to avoid unintended
				// modification of the inputs
				String[] tmp = new String[superset.length];
				for (int i = 0; i < tmp.length; i++) {
					tmp[i] = superset[i].toLowerCase();
				}
				superset = tmp;
			}

			// Each index in the outer array is a target, with a parallel position in finalIndices
			for (int i = 0; i < targetsWithAlts.length; i++) {
				// Track the current superset index for this target
				int index = -1;
				// Matching any alt in the superset is sufficient for this target
				for (int j = 0; j < targetsWithAlts[i].length; j++) {
					String alt = targetsWithAlts[i][j];
					if (!caseSensitive) {
						// Again we do have to compare each alt to each superset string, so we do a single
						// conversion to lowercase to reduce the loop complexity
						alt = alt.toLowerCase();
					}

					// We have to loop through the complete superset to report multiple matches
					for (int k = 0; k < superset.length; k++) {
						if (superset[k].contains(alt) || alt.contains(superset[k])) {
							if (index == -1) {
								// First match so we can just set the index for this target
								index = k;
							} else {
								// We've already found a match for this target
								err = true;
								if (verbose) {
									log.reportError("Error - multiple alts found in superset: " + superset[index]
																	+ ", " + superset[k]);
								}

								// If we prefer the first match in the alt list, then the previously found index
								// value is the preferred value (since we are iterating in order of the alt array)
								// If we prefer the first match in the target list, then we take the smaller of the
								// indices
								if (!preferFirstInTargetsOverFirstInSuperset && k < index) {
									index = k;
								}
							}
						}
					}
				}
				// Set the final value after iterating through all alts
				finalIndices[i] = index;
				if (index == -1) {
					// target not found
					err = true;
					if (verbose) {
						log.reportError("Error - no factor matched alts: "
														+ Arrays.stream(targetsWithAlts[i]).collect(Collectors.joining(", ")));
					}
				}
			}
			return finalIndices;
		}

		// If we are doing an exact match, using a Map over the superset can be much faster than
		// iterating. However, if we are only looking for a single target then creating the map would be
		// a waste
		if (targetsWithAlts.length == 1 && targetsWithAlts[0].length == 1) {
			// If we're looking up a single string, use indexOfStr to avoid computation cost of building a
			// map.
			return new int[] {indexOfStr(targetsWithAlts[0][0], superset, caseSensitive, true, log,
																	 verbose)};
		}

		// Make the map. NB: this method takes care of reporting duplicate strings in the superset
		Map<String, Integer> indexMap = makeIndexMap(superset, caseSensitive, verbose, log);

		List<ImmutableSet<String>> targetsWithAltsList = Arrays.stream(targetsWithAlts)
																													 .map(Arrays::stream)
																													 .map(s -> s.map(s1 -> caseSensitive ? s1
																																															 : s1.toLowerCase())
																																			.collect(ImmutableSet.toImmutableSet()))
																													 .collect(Collectors.toList());
		// For each target, we'll check the superset map for each alt
		for (int i = 0; i < targetsWithAltsList.size(); i++) {
			int index = -1;
			ImmutableSet<String> alts = targetsWithAltsList.get(i);
			for (String alt : alts) {
				if (indexMap.containsKey(alt)) {
					// If the alt is found, we check the index
					int altIndex = indexMap.get(alt);
					if (index == -1) {
						// This is the first alt we've found so use its index
						index = altIndex;
					} else {
						// Multiple matches
						if (verbose) {
							log.reportError("Error - multiple alts found in superset: " + superset[index] + ", "
															+ superset[altIndex]);
						}

						// This logic is the same as for the inexact algorithm
						if (!preferFirstInTargetsOverFirstInSuperset && altIndex < index) {
							index = altIndex;
						}
					}
				}
			}
			finalIndices[i] = index;
			if (index == -1) {
				// Target not found
				err = true;
				if (verbose) {
					log.reportError("Error - no factor matched alts: "
													+ Arrays.stream(targetsWithAlts[i]).collect(Collectors.joining(", ")));
				}
			}
		}
		//
		// if (kill && err) {
		// // System exit if requested
		// System.exit(1);
		// }

		return finalIndices;
	}


	/**
	 * Helper method to turn an array of String to a map of those strings to the first index that
	 * string appears in the array.
	 *
	 * @param superset The array to convert to a map
	 * @param caseSensitive If false, all strings will be {@link String#toLowerCase()}'d
	 * @param verbose Whether or not to print error messages
	 * @param log Log to use if verbose
	 * @param kill If multiple strings are present in the given superset, will exit
	 * @return The constructed map
	 */
	private static Map<String, Integer> makeIndexMap(String[] superset, boolean caseSensitive,
																									 boolean verbose, Logger log) {
		Map<String, Integer> supersetMap = new HashMap<>();
		boolean err = false;
		for (int i = 0; i < superset.length; i++) {
			String s = superset[i];
			if (!caseSensitive) {
				s = s.toLowerCase();
			}
			if (supersetMap.containsKey(s)) {
				if (verbose) {
					log.reportError("Error - more than one factor was named '" + superset[i] + "'");
				}
				err = true;
			} else {
				supersetMap.put(s, i);
			}
		}

		return supersetMap;
	}

	/**
	 * @see #indexMap(String[], String[], boolean, Logger, boolean, boolean)
	 */
	public static Map<String, Integer> indexMap(String[] subset, String[] superset,
																							boolean casesensitive, boolean kill) {
		return indexMap(subset, superset, casesensitive, new Logger(), true, kill);
	}

	/**
	 * @param subset
	 * @param superset
	 * @param casesensitive
	 * @param log
	 * @param verbose
	 * @param kill
	 * @return Map from subset strings to the indices of those strings in the superset
	 */
	public static Map<String, Integer> indexMap(String[] subset, String[] superset,
																							boolean casesensitive, Logger log, boolean verbose,
																							boolean kill) {
		Map<String, Integer> indices = new HashMap<String, Integer>(subset.length);

		if (!casesensitive) {
			// TODO implement CaseInsensitiveMap solution (maybe Apache Commons)
			int[] indexArr = indexFactors(subset, superset, casesensitive, log, verbose);
			for (int i = 0; i < subset.length; i++) {
				indices.put(subset[i], indexArr[i]);
			}
		} else {
			Set<String> targets = new HashSet<String>(Arrays.asList(subset));
			boolean err = false;


			for (int i = 0; i < superset.length; i++) {
				if (targets.contains(superset[i])) {
					if (indices.get(superset[i]) != null) {
						if (verbose) {
							log.reportError("Error - more than one factor was named '" + superset[i] + "'");
						}
						err = true;
					}
					indices.put(superset[i], i);
				}
			}

			if (verbose || kill) {
				for (String key : targets) {
					if (!indices.containsKey(key)) {
						if (verbose) {
							log.reportError("Error - no factor was named '" + key + "'");
						}
						err = true;
					}
				}
			}

			if (kill && err) {
				System.exit(1);
			}
		}

		return indices;
	}

	public static boolean checkHeader(String[] observed, String[] expected, boolean kill) {
		return checkHeader(observed, expected, true, kill);
	}

	public static boolean checkHeader(String[] observed, String[] expected, boolean caseSensitive,
																		boolean kill) {
		return checkHeader(observed, expected, caseSensitive, new Logger(), kill);
	}

	public static boolean checkHeader(String[] observed, String[] expected, boolean caseSensitive,
																		Logger log, boolean kill) {
		boolean kosher = true;


		if (observed == null) {
			System.err.println("Error - null header observed");
			if (kill) {
				System.exit(1);
			}
			return false;
		}

		if (observed.length != expected.length) {
			log.reportError("Error - file has an unexpected header; expecting " + expected.length
											+ " columns, found " + observed.length);
			log.reportError("Expected: " + ArrayUtils.toStr(expected));
			log.reportError("Observed: " + ArrayUtils.toStr(observed));
			kosher = false;
			if (kill) {
				System.exit(1);
			}
		}
		for (int i = 0; i < (observed.length < expected.length ? observed.length
																													 : expected.length); i++) {
			if (!observed[i].equalsIgnoreCase(expected[i])
					|| (caseSensitive && !observed[i].equals(expected[i]))) {
				log.reportError("Error - Expecting " + expected[i] + " in column " + (i + 1) + "; got "
												+ observed[i]);
				kosher = false;
			}
		}
		if (kill && !kosher) {
			System.exit(1);
		}

		return kosher;
	}

	public static boolean checkHeader(String[] observed, String[] expected, int[] indicesToCheck,
																		boolean caseSensitive, Logger log, boolean kill) {
		String[] subArray;

		if (indicesToCheck == null) {
			subArray = observed;
		} else {
			subArray = new String[indicesToCheck.length];
			for (int i = 0; i < indicesToCheck.length; i++) {
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
		return (col >= 26 ? (char) (64 + col / 26) + "" : "") + (char) (64 + (col % 26) + 1);
	}

	public static int getIndexFromColumn(String col) {
		int index = 0;

		if (col.length() > 2) {
			System.err.println("Error - column needs to be two characters or less");
			System.exit(1);
		}
		if (col.length() == 2) {
			index += (col.charAt(0) - 64) * 26;
		}
		index += col.charAt(col.length() - 1) - 64 - 1;

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
			if (contents == null) {
				return ("Clipboard is empty");
			} else {
				try {
					if (contents.isDataFlavorSupported(DataFlavor.stringFlavor)) {
						return (String) contents.getTransferData(DataFlavor.stringFlavor);
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

	public static String listWithCommas(Collection<String> items) {
		return listWithCommas(items.toArray(new String[items.size()]));
	}

	public static String listWithCommas(String[] list) {
		return listWithCommas(list, true);
	}

	public static String listWithCommas(String[] list, boolean finalComma) {
		String str = "";

		for (int i = 0; i < list.length; i++) {
			str += list[i];
			if (i < list.length - 1 - (finalComma && list.length > 2 ? 0 : 1)) {
				str += ", ";
			}
			if (list.length > 1 && i == list.length - 2) {
				str += (list.length == 2 || !finalComma ? " " : "") + "and ";
			}
		}

		return str;
	}

	public static String listRanges(int[] array) {
		if (array.length == 0) {
			return "";
		}

		StringBuilder sb = new StringBuilder();
		int[] ordered = Arrays.copyOf(array, array.length);
		Arrays.sort(ordered);

		sb.append(ordered[0]);
		for (int i = 1; i < ordered.length; i++) {
			if (ordered[i] - ordered[i - 1] > 1) {
				sb.append(",");
				sb.append(ordered[i]);
			} else if (i + 1 == ordered.length || ordered[i + 1] - ordered[i] > 1) {
				sb.append("-");
				sb.append(ordered[i]);
			}
		}

		return sb.toString();
	}

	public static String capitalizeFirst(String str) {
		if (str.length() == 0) {
			return str;
		} else {
			return str.substring(0, 1).toUpperCase() + str.substring(1);
		}
	}

	public static String uncapitalizeFirst(String str) {
		if (str.length() == 0) {
			return str;
		} else {
			return str.substring(0, 1).toLowerCase() + str.substring(1);
		}
	}

	public static String capitalizeWords(String str) {
		String[] arr = str.split(" ");
		StringBuffer sb = new StringBuffer();

		for (String element : arr) {
			sb.append(Character.toUpperCase(element.charAt(0)))
				.append(element.substring(1).toLowerCase())
				.append(" ");
		}
		return sb.toString().trim();
	}

	public static double doubleParser(String str) {
		try {
			return Double.parseDouble(str);
		} catch (Exception e) {
			return Double.NEGATIVE_INFINITY;
		}
	}

	public static String getExcelColumn(int index) {
		return "" + (index / 26 > 0 ? (char) (index / 26 - 1 + 65) : "") + (char) (index % 26 + 65);
	}

	public static String removeQuotes(String str) {
		if (str.startsWith("\"") && str.endsWith("\"") && str.length() > 1) {
			return str.substring(1, str.length() - 1);
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
			stop = str.substring(start + 1).indexOf("\"") + start + 1;
			try {
				str = str.substring(0, start)
							+ replaceAllWith(str.substring(start + 1, stop), theReplacements)
							+ str.substring(stop + 1);
			} catch (Exception e) {
				System.err.println("Error - failed to parse: '" + ori + "'");
				e.printStackTrace();
				ext.waitForResponse();
			}
		}
		return str;
	}

	public static boolean containsAny(String str, String[] targets) {
		for (String element : targets) {
			if (str.contains(element)) {
				return true;
			}
		}
		return false;
	}

	public static boolean containsAny(String str, Iterable<String> targets) {
		for (String element : targets) {
			if (str.contains(element)) {
				return true;
			}
		}
		return false;
	}

	public static boolean containsAnyChar(String str, char[] targets) {
		for (char element : targets) {
			if (str.indexOf(element) != -1) {
				return true;
			}
		}
		return false;
	}

	public static boolean containsAnyChar(String str, Iterable<Character> targets) {
		if (str != null) {
			for (char element : targets) {
				if (str.indexOf(element) != -1) {
					return true;
				}
			}
		}
		return false;
	}

	public static boolean containsAll(String[] targets, String[] array) {
		for (String target : targets) {
			if (ext.indexOfStr(target, array) == -1) {
				return false;
			}
		}
		return true;
	}

	public static Date parseDate(String str) {
		String[] line;

		line = str.trim().split("/");
		if (line.length != 3) {
			System.err.println("Error parsing date ('" + str + "'); assuming MM/DD/YYYY format");
		}

		return new GregorianCalendar(Integer.parseInt(line[2]), Integer.parseInt(line[0]) - 1,
																 Integer.parseInt(line[1])).getTime();
	}

	public static int calcDays(Date start, Date end) {
		return (int) ((end.getTime() - start.getTime()) / 1000 / 60 / 60 / 24);
	}

	public static String rootOf(String filename) {
		return rootOf(filename, true);
	}

	public static String rootOf(String filename, boolean trimDirectoryInfo) {
		if (trimDirectoryInfo && filename.indexOf("/") >= 0) {
			filename = filename.substring(filename.lastIndexOf("/") + 1);
		}
		if (trimDirectoryInfo && filename.indexOf("\\") >= 0) {
			filename = filename.substring(filename.lastIndexOf("\\") + 1);
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
		return ext.rootOf(filename, false) + addition
					 + (filename.lastIndexOf(".") > 0 ? filename.substring(filename.lastIndexOf(".")) : "");
	}

	public static String removeDirectoryInfo(String filename) {
		if (filename.indexOf("/") >= 0) {
			filename = filename.substring(filename.lastIndexOf("/") + 1);
		}
		if (filename.indexOf("\\") >= 0) {
			filename = filename.substring(filename.lastIndexOf("\\") + 1);
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

		if (filename.indexOf("/") >= 0) {
			dir = filename.substring(0, filename.lastIndexOf("/") + 1);
		}

		if (dir.equals("") && !allowBlank) {
			dir = "./";
		}

		return dir;
	}

	public static void writeToAll(String str, PrintWriter[] writers) {
		for (PrintWriter writer : writers) {
			writer.println(str);
		}
	}

	public static String addCommas(int value) {
		String str = value + "";

		for (int i = str.length() - 3; i > 0; i = i - 3) {
			str = str.substring(0, i) + "," + str.substring(i);
		}

		return str;
	}

	public static String addCommas(long value) {
		String str = value + "";

		for (int i = str.length() - 3; i > 0; i = i - 3) {
			str = str.substring(0, i) + "," + str.substring(i);
		}

		return str;
	}

	public static void appendToFile(String text, String filename) {
		PrintWriter writer;

		try {
			writer = Files.openAppropriateWriter(filename, true);
			writer.println(text);
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + filename);
			e.printStackTrace();
		}
	}

	public static String verifyDirFormat(String str) {
		// if (str.startsWith("C:")) {
		// str = str.substring(2);
		// }
		str = replaceAllWith(str, "\\", "/");
		if (str.length() > 0 && !str.endsWith("/")) {
			str += "/";
		}

		if (str.endsWith("/./")) {
			str = str.substring(0, str.length() - 2);
		}

		return str;
	}

	public static boolean parseBooleanArg(String arg) {
		return parseBooleanArg(arg, new Logger());
	}

	public static boolean parseBooleanArg(String arg, Logger log) {
		if (arg.split("=")[1].trim().equalsIgnoreCase("true")) {
			return true;
		} else if (arg.split("=")[1].trim().equalsIgnoreCase("false")) {
			return false;
		} else {
			log.reportError("Error - invalid " + arg.split("=")[0]
											+ "= argument (expecting true/false): "
											+ arg.split("=")[1].trim());
			System.exit(1);
			return false;
		}
	}

	public static String parseStringArg(String arg) {
		return parseStringArg(arg, "");
	}

	public static String parseStringArg(String arg, String blankValue) {
		if (arg.split("=").length == 1) {
			return blankValue;
		} else if (arg.split("=")[1].trim().equalsIgnoreCase("null")) {
			return null;
		} else {
			return arg.substring(arg.indexOf("=") + 1).trim();
		}
	}

	public static int parseIntArg(String arg) {
		try {
			return Integer.parseInt(arg.split("=")[1].trim());
		} catch (NumberFormatException nfe) {
			System.err.println("Error - invalid " + arg.split("=")[0] + "= argument: "
												 + arg.split("=")[1].trim());
			System.exit(1);
			return Integer.MIN_VALUE;
		}
	}

	/**
	 * Parse a string array argument;
	 *
	 * @param arg argument in form of <code>arg=str1,str2,string3</code>, which will return
	 *        <code>String[]{str1,str2,string3}</code>
	 * @return
	 */
	public static String[] parseStringArrayArg(String arg, String[] blankValue) {
		String[] pts = arg.trim().split("=");
		if (pts.length == 1) {
			return blankValue;
		}
		pts = pts[1].split(",");
		return pts;
	}

	/**
	 * Parse an integer array argument;
	 *
	 * @param arg argument in form of <code>arg=2,4-10,20</code>, which will return
	 *        <code>int[]{2,4,5,6,7,8,9,10,20}</code>
	 * @return
	 */
	public static int[] parseIntArrayArg(String arg) {
		String lst = arg.trim().split("=")[1];
		String[] pts = lst.split(",");
		ArrayList<Integer> returnList = new ArrayList<Integer>();
		for (String pt : pts) {
			// check for range
			if (pt.contains("-")) {
				int begin = Integer.parseInt(pt.split("-")[0]);
				int end = Integer.parseInt(pt.split("-")[1]);
				for (int i = begin; i < end + 1; i++) {
					returnList.add(i);
				}
			} else {
				int ent = Integer.parseInt(pt);
				returnList.add(ent);
			}
		}
		return Ints.toArray(returnList);
	}

	public static byte parseByteArg(String arg) {
		try {
			return Byte.parseByte(arg.split("=")[1].trim());
		} catch (NumberFormatException nfe) {
			System.err.println("Error - invalid " + arg.split("=")[0] + "= argument: "
												 + arg.split("=")[1]);
			System.exit(1);
			return Byte.MIN_VALUE;
		}
	}

	public static double parseDoubleArg(String arg) {
		try {
			return Double.parseDouble(arg.split("=")[1].trim());
		} catch (NumberFormatException nfe) {
			System.err.println("Error - invalid " + arg.split("=")[0] + "= argument: "
												 + arg.split("=")[1]);
			System.exit(1);
			return Double.NaN;
		}
	}

	public static boolean isMissingValue(String str) {
		for (String element : MISSING_VALUES) {
			if (str.trim().equalsIgnoreCase(element)) {
				return true;
			}
		}
		return false;
	}

	public static double getValidDouble(String str) {
		if (isMissingValue(str)) {
			return 0.0;
		}
		try {
			return Double.parseDouble(str);
		} catch (NumberFormatException nfe) {
			return 0.0;
		}
	}

	public static int getValidInteger(String str) {
		if (isMissingValue(str)) {
			return 0;
		}
		try {
			return Integer.parseInt(str);
		} catch (NumberFormatException nfe) {
			return 0;
		}
	}

	public static boolean isValidDouble(String str) {
		if (isMissingValue(str)) {
			return false;
		}

		try {
			Double.parseDouble(str.trim());
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
			Integer.parseInt(str.trim());
		} catch (NumberFormatException nfe) {
			return false;
		}

		return true;
	}

	public static boolean isValidByte(String str) {
		return isValidInteger(str) && Integer.parseInt(str) >= Byte.MIN_VALUE
					 && Integer.parseInt(str) <= Byte.MAX_VALUE;
	}

	public static boolean isValidChromosome(String str) {
		return ext.indexOfStr(str, Positions.CHR_CODES) != -1
					 || (isValidInteger(str) && Integer.parseInt(str) >= 1 && Integer.parseInt(str) <= 26);
	}

	public static float parseFloatArg(String arg) {
		try {
			return Float.parseFloat(arg.split("=")[1].trim());
		} catch (NumberFormatException nfe) {
			System.err.println("Error - invalid " + arg.split("=")[0] + "= argument: "
												 + arg.split("=")[1]);
			System.exit(1);
			return Float.NaN;
		}
	}

	public static String insertNumbers(String pattern, int num) {
		return ext.replaceAllWith(ext.replaceAllWith(pattern, "##", chrome(num)), "#", num + "");
	}

	public static String insertNumbers(String pattern, int num, int numDigits) {
		if (numDigits == -1) {
			return ext.replaceAllWith(pattern, "#", num + "");
		} else {
			return ext.replaceAllWith(pattern, "#", ext.formNum(num, numDigits));
		}
	}

	public static String right(String str, int numChar) {
		return str.substring(str.length() - numChar, str.length());
	}

	public static String link(String linkRoot, String dir) {
		// int count;

		linkRoot = verifyDirFormat(linkRoot);

		// count = 0;
		while (dir.startsWith("../")) {
			linkRoot = linkRoot.substring(0, linkRoot.lastIndexOf("/"));
			dir = dir.substring(3);
			// count++;
		}

		return linkRoot.substring(0, linkRoot.lastIndexOf("/") + 1) + dir;
	}

	public static int countInstancesOf(String str, String pattern) {
		int count;
		String trav;

		trav = str;
		count = 0;
		while (trav.contains(pattern)) {
			trav = trav.substring(trav.indexOf(pattern) + pattern.length());
			count++;
		}

		return count;
	}

	public static String determineDelimiter(String str) {
		return determineDelimiter(str, false);
	}

	/**
	 * @param str
	 * @param spaceAsLastResort if true, we return {@link PSF.Regex#GREEDY_WHITESPACE} only if \t and
	 *        , are not found in the header
	 * @return
	 */
	public static String determineDelimiter(String str, boolean spaceAsLastResort) {
		if (str.contains("\t")) {
			return "\t";
		} else if (countInstancesOf(str, ",") > countInstancesOf(str, " ")
							 || (spaceAsLastResort && countInstancesOf(str, ",") > 0)) {
			return ",";
		} else {
			return PSF.Regex.GREEDY_WHITESPACE;
		}
	}

	public static String removeAndSimplifyQuotes(String str, Logger log) {
		if (str.startsWith("\"")) {
			if (!str.endsWith("\"") && log != null) {
				log.reportError("Error - improperly formed quotes in comma delimited token:");
				log.reportError(str);
			}
			str = str.substring(1, str.length() - 1);
		}
		return ext.replaceAllWith(str, "\"\"", "\"");
	}

	public static String enquote(String str) {
		return "\"" + str + "\"";
	}

	public static String[] splitCommasIntelligently(String str, boolean removeAndSimplifyQuotes,
																									Logger log) {
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
					v.add(removeAndSimplifyQuotes
																				? removeAndSimplifyQuotes(str.substring(startIndex, i), log)
																				: str.substring(startIndex, i));
					// System.out.println(v.elementAt(v.size()-1));
					startIndex = i + 1;
					if (insideQuotes && str.charAt(i - 1) != '\"') {
						if (log != null) {
							log.reportError("Error - improperly formed quotes in comma delimited string:");
							log.reportError(str);
						}
					}
					insideQuotes = false;
				}
			}
			if (str.charAt(i) == '\"') {
				numCommas++;
				insideQuotes = true;
			}
		}
		v.add(removeAndSimplifyQuotes ? removeAndSimplifyQuotes(str.substring(startIndex), log)
																	: str.substring(startIndex));

		return ArrayUtils.toStringArray(v);
	}

	public static String reportMemoryUsage() {
		MemoryUsage muHeap; // , muNonHeap;

		muHeap = ManagementFactory.getMemoryMXBean().getHeapMemoryUsage();
		// muNonHeap = ManagementFactory.getMemoryMXBean().getNonHeapMemoryUsage();

		// return "Total: " + ext.prettyUpSize(Runtime.getRuntime().totalMemory(), 1) +
		// " Free: " + ext.prettyUpSize(Runtime.getRuntime().freeMemory(),1) +
		// " Max: " + ext.prettyUpSize(Runtime.getRuntime().maxMemory(),1) +
		// " HeapComit: " + ext.prettyUpSize(muHeap.getCommitted(),1) +
		// " HeapInit: " + ext.prettyUpSize(muHeap.getInit(),1) +
		// " HeapMax: " + ext.prettyUpSize(muHeap.getMax(),1) +
		// " HeapUsed: " + ext.prettyUpSize(muHeap.getUsed(),1) +
		// " NonHeapComit: " + ext.prettyUpSize(muNonHeap.getCommitted(),1) +
		// " NonHeapInit: " + ext.prettyUpSize(muNonHeap.getInit(),1) +
		// " NonHeapMax: " + ext.prettyUpSize(muNonHeap.getMax(),1) +
		// " NonHeapUsed: " + ext.prettyUpSize(muNonHeap.getUsed(),1)
		// ;

		return "Used: " + ext.prettyUpSize(muHeap.getUsed(), 1) + "\tInit: "
					 + ext.prettyUpSize(muHeap.getInit(), 1) + "\tCommit: "
					 + ext.prettyUpSize(muHeap.getCommitted(), 1) + "\tMax: "
					 + ext.prettyUpSize(Runtime.getRuntime().maxMemory(), 1) + "\tFree: "
					 + ext.prettyUpSize(Runtime.getRuntime().freeMemory(), 1);
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
			return str.trim().split(delimiter, -1);
		}
	}

	public static void waitForResponse() {
		waitForResponse("Press ENTER to continue");
	}

	public static void waitForResponse(String message) {
		if (Files.isWindows()) {
			System.out.println(message);
			try {
				new BufferedReader(new InputStreamReader(System.in)).readLine();
			} catch (IOException ioe) {
			}
		}
	}

	public static int getNumSigFig(double num) {
		String str = formDeci(num, 5);

		if (str.contains(".")) {
			int ind = str.indexOf(".");
			return str.length() - ind - 1;
		} else {
			return 0;
		}
	}

	public static double roundToSignificantFigures(double num, int n) {
		// from
		// https://stackoverflow.com/questions/202302/rounding-to-an-arbitrary-number-of-significant-digits
		if (num == 0) {
			return 0;
		}

		final double d = Math.ceil(Math.log10(num < 0 ? -num : num));
		final int power = n - (int) d;

		final double magnitude = Math.pow(10, power);
		final long shifted = Math.round(num * magnitude);
		return shifted / magnitude;
	}

	/**
	 * @param toAddTo counts will be added to this hash
	 * @param toAddFrom hash to add counts from
	 * @return NOTE: couldn't make this generic for {@link Numbers}. Go for it if it's easy
	 */
	public static Hashtable<String, Integer> addHashCounts(Hashtable<String, Integer> toAddTo,
																												 Hashtable<String, Integer> toAddFrom) {
		Set<String> keys = toAddFrom.keySet();
		for (String key : keys) {
			if (!toAddTo.containsKey(key)) {
				toAddTo.put(key, toAddFrom.get(key));
			} else {
				int cur = toAddTo.get(key);
				int curToAdd = toAddFrom.get(key);
				toAddTo.put(key, cur + curToAdd);
			}
		}
		return toAddTo;
	}

	/**
	 * Computes the Levenshtein distance score between two Strings (Copied from
	 * http://rosettacode.org/wiki/Levenshtein_distance#Java)
	 *
	 * @param a First String
	 * @param b Second String
	 * @return Levenshtein difference score
	 */
	public static int distanceLevenshtein(String a, String b) {
		a = a.toLowerCase();
		b = b.toLowerCase();

		// i == 0
		int[] costs = new int[b.length() + 1];
		for (int j = 0; j < costs.length; j++) {
			costs[j] = j;
		}

		for (int i = 1; i <= a.length(); i++) {
			// j == 0; nw = lev(i - 1, j)
			costs[0] = i;
			int nw = i - 1;
			for (int j = 1; j <= b.length(); j++) {
				int cj = Math.min(1 + Math.min(costs[j], costs[j - 1]),
													a.charAt(i - 1) == b.charAt(j - 1) ? nw : nw + 1);
				nw = costs[j];
				costs[j] = cj;
			}
		}
		return costs[b.length()];
	}

	public static void main(String[] args) {
		String[] data = {"kitten", "sitting", "saturday", "sunday", "rosettacode", "raisethysword"};
		for (int i = 0; i < data.length; i += 2) {
			System.out.println("distance(" + data[i] + ", " + data[i + 1] + ") = "
												 + distanceLevenshtein(data[i], data[i + 1]));
		}
		System.exit(1);
	}
}
