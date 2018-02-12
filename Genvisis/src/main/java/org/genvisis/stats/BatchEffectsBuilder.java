package org.genvisis.stats;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import javax.annotation.Nonnull;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class BatchEffectsBuilder {

  private final List<String> batchLabels;
  private final List<String> factorLabels;
  // the batchValuesToInclude list and factorValuesToInclude sublists have parallel indices, each
  // index representing a value (batch value or factor value, respectively) for a particular sample
  private final List<String> batchValuesToInclude;
  // each superlist in factorValuesToInclude represents the sample values for a particular factor.
  private final List<List<Double>> factorValuesToInclude;
  private final Logger logger;

  public BatchEffectsBuilder(Logger logger) {
    batchLabels = new ArrayList<>();
    factorLabels = new ArrayList<>();
    batchValuesToInclude = new ArrayList<>();
    factorValuesToInclude = new ArrayList<>();
    this.logger = logger;
  }

  /**
   * Parses batch and factor data from file. Returns BatchEffects instance containing this data in
   * immutable data structures.
   * 
   * @param batchFilePath String path to file containing labels in first row, sample identifiers in
   *          first column, and batch values in second column.
   * @param factorFilePath String path to file containing labels in first row, sample identifiers in
   *          first column, and factor values in subsequent columns.
   * @return BatchEffects instance containing parsed batch and factor data
   * @throws FileNotFoundException if either String path parameter does not lead to a valid file.
   * @throws ParseException if file with factor data contains fewer than 2 column labels.
   * @throws IOException
   */
  public BatchEffects build(@Nonnull String batchFilePath,
                            @Nonnull String factorFilePath) throws FileNotFoundException,
                                                            ParseException, IOException {
    logger.report("parsing data from files...");

    // parse batch file
    Map<String, String> sampleBatchMap = HashVec.loadFileToHashString(batchFilePath, new int[] {0},
                                                                      new int[] {1}, false, null,
                                                                      true, false);

    // parse factor file
    try (BufferedReader reader = Files.getAppropriateReader(factorFilePath)) {
      String[] rowData;
      int expectedNumColumns = -1;
      // read factor column headers
      String line = reader.readLine();
      if (line != null) {
        rowData = line.split("\t");
        if (rowData.length >= 2) {
          expectedNumColumns = rowData.length;
          for (int i = 1; i < rowData.length; i++) { // this loop starts at 1 to exclude the sample
                                                    // id
                                                    // column
            factorLabels.add(rowData[i]);
            factorValuesToInclude.add(new ArrayList<>()); // each sub-list represents a different
                                                         // factor
          }
        } else {
          throw new ParseException("Error: expected at least 2 column headers in factor data. Found "
                                   + rowData.length + ".", 0);
        }
      }
      // read factor sample data
      long numSampleLinesFromFactorData = 0;
      long numSamplesIncluded = 0;
      while ((line = reader.readLine()) != null) {
        numSampleLinesFromFactorData++;
        rowData = line.split("\t");
        boolean sampleIncluded = addSampleIfValid(rowData, sampleBatchMap, expectedNumColumns);
        if (sampleIncluded) {
          numSamplesIncluded++;
        } else {
          logger.report("line excluded from factor file: " + line);
        }
      }
      logger.report("parsing complete");
      logger.report(numSamplesIncluded + " samples will be included in analysis. "
                    + (numSampleLinesFromFactorData - numSamplesIncluded)
                    + " lines of potential sample data were excluded from the factor file. "
                    + (sampleBatchMap.size() - numSamplesIncluded)
                    + " lines of potential sample data were excluded from the batch file. ");
    }
    return new BatchEffects(this);
  }

  /**
   * If sample passes validation, adds sample values to class-level data structures that will be
   * used for analysis. Sample will fail validation if it is missing any data or if any of its
   * factor values cannot be parsed as a double.
   * 
   * @param factorDataForSample array of factor data. index 0 should be the sample identifier.
   *          subsequent indices should be factor values for that sample.
   * @param sampleBatchMap map of sample identifier to batch value.
   * @param expectedLength expected length of factorDataForSample parameter.
   * @return a boolean value indicating whether the sample was added to class-level data structures
   *         that will be used for analysis.
   */
  private boolean addSampleIfValid(String[] factorDataForSample, Map<String, String> sampleBatchMap,
                                   int expectedLength) {
    boolean includeSample = true;

    if (factorDataForSample.length != expectedLength) {
      includeSample = false;
    }

    // exclude sample if the sample id from the factor file does not match any from the batch file
    // OR if the batch value associated with that sample id is missing
    if (includeSample) {
      String batchValue = sampleBatchMap.get(factorDataForSample[0]);
      if (batchValue == null || ext.isMissingValue(batchValue)) {
        includeSample = false;
      }
    }

    // exclude sample if any of its factor values are missing or do not parse into doubles
    List<Double> parsedValues = new ArrayList<>();
    if (includeSample) {
      for (int i = 0; includeSample && i < factorDataForSample.length; i++) {
        if (ext.isMissingValue(factorDataForSample[i])) {
          includeSample = false;
        }
        if (includeSample && i != 0) {
          try {
            double d = Double.parseDouble(factorDataForSample[i]);
            parsedValues.add(d);
          } catch (NumberFormatException e) {
            includeSample = false;
          }
        }
      }
    }

    // if sample passes validation, add to class-level data structures that will be used for
    // analysis
    if (includeSample) {
      String batchValue = sampleBatchMap.get(factorDataForSample[0]); // index 0 is the sample id
      if (!batchLabels.contains(batchValue)) {
        batchLabels.add(batchValue);
      }
      batchValuesToInclude.add(batchValue);
      for (int i = 0; i < parsedValues.size(); i++) {
        double currentValue = parsedValues.get(i);
        factorValuesToInclude.get(i).add(currentValue);
      }
    }

    // return boolean indication of whether sample was included
    return includeSample;
  }

  /**
   * @return the batchLabels
   */
  List<String> getBatchLabels() {
    return batchLabels;
  }

  /**
   * @return the factorLabels
   */
  List<String> getFactorLabels() {
    return factorLabels;
  }

  /**
   * @return the batchValuesToInclude
   */
  List<String> getBatchValuesToInclude() {
    return batchValuesToInclude;
  }

  /**
   * @return the factorValuesToInclude
   */
  List<List<Double>> getFactorValuesToInclude() {
    return factorValuesToInclude;
  }

  /**
   * @return the logger
   */
  Logger getLogger() {
    return logger;
  }

}
