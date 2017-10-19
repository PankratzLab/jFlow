package org.genvisis.seq.analysis;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import javax.annotation.Nonnull;

import org.apache.commons.cli.ParseException;
import org.genvisis.CLI;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;

/**
 * Filters markers from a blast.vcf file based on search term occurrence within the value of a given
 * marker token.
 */
public class BlastTable {

	private final String targetToken;
	private final String[] searchTerms;
	private final Logger logger;
	private final String inputFile;
	private final String outputFile;
	private static final String[] MARKER_COLUMNS_TO_INCLUDE = {"ID", "CHROM", "POS"};
	private static final String MARKER_COLUMN_TO_SEARCH = "INFO";

	/**
	 * Constructs a BlastTable instance.
	 * 
	 * @param token String marker token within blast.vcf file. if token exists, it will be searched
	 *        for search terms
	 * @param searchTerms String[] containing terms to search for within the value of given marker
	 *        token
	 * @param inputFile String path to blast.vcf file
	 * @param outputFile String path to text file. this file will be used for marker results after
	 *        filtering
	 * @param logger logger to be used for error and log reporting
	 */
	public BlastTable(@Nonnull String token, @Nonnull String searchTerms[], @Nonnull String inputFile,
										@Nonnull String outputFile, @Nonnull Logger logger) {
		this.targetToken = token;
		this.searchTerms = searchTerms;
		this.logger = logger;
		this.inputFile = inputFile;
		this.outputFile = outputFile;
	}

	/**
	 * Filters markers from a blast.vcf file based on search term occurrence within the value of a
	 * given marker token. Markers whose token contains any of the given search terms will be written
	 * to an output file.
	 * 
	 * @param args
	 * @throws ParseException if a problem is encountered during parsing or a help/usage flag is
	 *         parsed
	 * @throws IOException if a problem is encountered reading from input file or writing output file
	 */
	public static void main(String[] args) throws ParseException, IOException {
		CLI cli = new CLI("BlastTable");
		String inFile = "inputFile";
		cli.addArg(inFile, "Path to input blast.vcf file.", true);
		String outFile = "outputFile";
		cli.addArg(outFile, "Path to output text file.", true);
		String token = "token";
		cli.addArg(token, "Token within blast.vcf info section that will be searched.", true);
		String delim = "::";
		String delimSearchTerms = "searchTerms";
		cli.addArg(delimSearchTerms, "List of intra-token search terms delimited by \"" + delim, true);
		cli.parse(args);
		String[] searchTerms = cli.get(delimSearchTerms).split(delim);
		Logger logger = new Logger();
		BlastTable instance = new BlastTable(cli.get(token), searchTerms, cli.get(inFile),
																				 cli.get(outFile), logger);
		boolean success = instance.createBlastTableFile();
		if (success) {
			logger.report("BlastTable completed successfully.");
		}
	}

	/**
	 * Filters markers from blast.vcf file based on search term occurrence within the value of marker
	 * token. Markers whose token contains any of the given search terms will be written to an output
	 * file.
	 * 
	 * @return a boolean indicating whether the method completed successfully
	 * @throws IOException if a problem is encountered while reading from blast.vcf file or writing
	 *         output file
	 */
	public boolean createBlastTableFile() {
		boolean completedSuccessfully = true;
		try (PrintWriter writer = Files.getAppropriateWriter(outputFile, false);
				 BufferedReader reader = Files.getAppropriateReader(inputFile);) {

			// find line with column labels from input file, and map column headers to column numbers
			Map<String, Integer> mapOfLabelToColumnNumber = new HashMap<>();
			String line;
			boolean foundHeader = false;
			while (!foundHeader && (line = reader.readLine()) != null) {
				if (line.charAt(0) == '#' && line.charAt(1) != '#') { // identify line with column labels
					String[] columns = line.split("\t");
					columns[0] = columns[0].substring(1, columns[0].length()); // remove hash tag from col1
					for (int i = 0; i < columns.length; i++) {
						mapOfLabelToColumnNumber.put(columns[i], i);
					}
					foundHeader = true;
				}
			}

			// check whether all expected column labels are present
			if (hasLabelsMissing((mapOfLabelToColumnNumber.keySet()))) {
				throw new BlastTableException("Canceled blast table creation - expected column labels were not found in blast.vcf file.");
			}

			// write column labels to output file
			String[] outputLabels = new String[MARKER_COLUMNS_TO_INCLUDE.length + searchTerms.length];
			int index = 0;
			for (String label : MARKER_COLUMNS_TO_INCLUDE) {
				outputLabels[index] = label;
				index++;
			}
			for (String label : searchTerms) {
				outputLabels[index] = label;
				index++;
			}
			String outputLine = String.join("\t", outputLabels);
			writer.println(outputLine);

			// read marker data from input file, filter, and write to output file
			while ((line = reader.readLine()) != null) {
				if (line.length() > 1 && line.charAt(0) != '#') {
					String[] columnData = line.split("\t");
					String targetColumn = columnData[mapOfLabelToColumnNumber.get(MARKER_COLUMN_TO_SEARCH)];
					String[] tokens = targetColumn.split(";");
					String token = null;
					boolean foundToken = false;
					for (int i = 0; !foundToken && i < tokens.length; i++) {
						if (tokens[i].length() > targetToken.length()) {
							if ((tokens[i].substring(0, targetToken.length()) + "=").equals(targetToken + "=")) {
								token = tokens[i];
								foundToken = true;
							}
						}
					}
					String[] outputValues = new String[MARKER_COLUMNS_TO_INCLUDE.length + searchTerms.length];
					if (foundToken) {
						String tokenValue = token.substring(targetToken.length() + 1, token.length());
						boolean includeSample = false;
						for (int i = 0; i < searchTerms.length; i++) {
							if (tokenValue.contains(searchTerms[i])) {
								outputValues[MARKER_COLUMNS_TO_INCLUDE.length + i] = token;
								includeSample = true;
							} else {
								outputValues[MARKER_COLUMNS_TO_INCLUDE.length + i] = ".";
							}
						}
						if (includeSample) {
							// transfer desired marker data to output
							for (int i = 0; i < MARKER_COLUMNS_TO_INCLUDE.length; i++) {
								String currentColumnLabel = MARKER_COLUMNS_TO_INCLUDE[i];
								outputValues[i] = columnData[mapOfLabelToColumnNumber.get(currentColumnLabel)];
							}
							writer.println(String.join("\t", outputValues));
						}
					} else {
						throw new BlastTableException("Canceled blast table creation - expected token was not found for every marker in blast.vcf file.");
					}
				}
			}
		} catch (IOException e) {
			logger.reportError(e.getMessage());
			completedSuccessfully = false;
		} catch (BlastTableException e) {
			logger.reportError(e.getMessage());
			completedSuccessfully = false;
		}
		return completedSuccessfully;
	}

	/**
	 * Checks whether expected label values are present in argument set.
	 * 
	 * @param superSet set to check for missing label values
	 * @return true if any expected label values are missing from superset, false if all expected
	 *         values are present
	 */
	private boolean hasLabelsMissing(Set<String> superSet) {
		for (String label : MARKER_COLUMNS_TO_INCLUDE) {
			if (!superSet.contains(label)) {
				return true;
			}
		}
		if (!superSet.contains(MARKER_COLUMN_TO_SEARCH)) {
			return true;
		}
		return false;
	}

	/**
	 * General exception used for internal validation.
	 */
	private class BlastTableException extends Exception {
		private static final long serialVersionUID = 1L; // default UID, no significance

		private BlastTableException(String message) {
			super(message);
		}
	}
}
