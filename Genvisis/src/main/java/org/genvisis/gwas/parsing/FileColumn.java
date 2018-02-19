package org.genvisis.gwas.parsing;

public interface FileColumn<T> {

  /**
   * Get the name of the column. Not expected to be unique. Used if writing output and to link
   * columns between files.
   * 
   * @return The column name
   */
  String getName();

  /**
   * Initialize this FileColumn. Repeated calls to this method should be no-ops.
   * 
   * @param FileParser
   */
  void initialize(FileParser parser);

  /**
   * Get this FileColumn's value from a line of data.
   * 
   * @param line {@code String[]} of line data
   * @return {@code <T>}
   * @throws ParseFailureException thrown if there is a
   */
  T getValue(String[] line) throws ParseFailureException;

  public boolean dieOnParseFailure();

}
