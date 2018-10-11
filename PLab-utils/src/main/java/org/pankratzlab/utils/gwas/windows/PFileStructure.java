package org.pankratzlab.utils.gwas.windows;

import org.pankratzlab.common.Aliases;

/**
 * Container class representing the minimum structure of a file containing one or more p-values
 */
public class PFileStructure {

  private String pvalFile;
  private String[] markerCols;
  private String[] chrCols;
  private String[] posCols;
  private String[] pvalCols;

  /**
   * Constructs a PFileStructure with default column labels from the {@link Aliases} class.
   *
   * @param pvalFile File to read p-value information from
   */
  public PFileStructure(String pvalFile) {
    this.pvalFile = pvalFile;
    markers(Aliases.MARKER_NAMES);
    chrs(Aliases.CHRS);
    pos(Aliases.POSITIONS);
    pvals(Aliases.PVALUES);
  }

  /**
   * Set the marker name column alias(es) for this file.
   *
   * @return This PFileStructure, for method chaining.
   */
  public PFileStructure markers(String... markerCols) {
    this.markerCols = markerCols;
    return this;
  }

  /**
   * Set the chromosome position column alias(es) for this file.
   *
   * @return This PFileStructure, for method chaining.
   */
  public PFileStructure chrs(String... chrCols) {
    this.chrCols = chrCols;
    return this;
  }

  /**
   * Set the genomic position column alias(es) for this file.
   *
   * @return This PFileStructure, for method chaining.
   */
  public PFileStructure pos(String... posCols) {
    this.posCols = posCols;
    return this;
  }

  /**
   * Set the p-value column name(s) for this file.
   *
   * @return This PFileStructure, for method chaining.
   */
  public PFileStructure pvals(String... pvalCols) {
    this.pvalCols = pvalCols;
    return this;
  }

  /**
   * @return The path to this p-value file
   */
  public String getFile() {
    return pvalFile;
  }

  /**
   * @return The list of marker name column aliases for this p-value file
   */
  public String[] getMarkers() {
    return markerCols;
  }

  /**
   * @return The list of chromosome position column aliases for this p-value file
   */
  public String[] getChrs() {
    return chrCols;
  }

  /**
   * @return The list of genomic position column aliases for this p-value file
   */
  public String[] getPos() {
    return posCols;
  }

  /**
   * @return The list of p-value columns in this p-value file
   */
  public String[] getPvals() {
    return pvalCols;
  }
}
