package org.genvisis.common.parsing;

public interface FileColumn<T> {

  /**
   * Get the name of the column. This name serves as the unique identifier and must not change after
   * construction.
   * 
   * @return The column name
   */
  String getName();

  /**
   * Get the header label for the column. This will often be the same as {@link #getName()} but will
   * only be called after {@link #initialize(FileParser)} and can therefore use data from the input
   * file where necessary
   * 
   * @return The header for the column when writing to an output file
   */
  String getHeader();

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

  /**
   * Compares the specified object with this {@link FileColumn} for equality. Returns {@code true}
   * if and only if the specified object is also a {@link FileColumn} and both return the same value
   * for {@link #getName()} This definition ensures that the equals method works predictably across
   * different implementations of the {@link FileColumn} interface.
   *
   * @param obj the object to be compared for equality with this {@link FileColumn}
   * @return {@code true} if the specified object is equal to this {@link FileColumn}
   */
  @Override
  boolean equals(Object obj);

  /**
   * Returns the hash code value for this {@link FileColumn}. The hash code of a {@link FileColumn}
   * is defined to be the result of the following calculation: <br>
   * <br>
   * {@code
   *   int hashCode = 31 + ((getName() == null) ? 0 : getName().hashCode());
   * }<br>
   * <br>
   * This ensures that <tt>fileColumn1.equals(fileColumn2)</tt> implies that
   * <tt>fileColumn1.hashCode()==fileColumn1.hashCode()</tt> for any two {@link FileColumn}s,
   * <tt>fileColumn1</tt> and <tt>fileColumn2</tt>, as required by the general contract of
   * {@link Object#hashCode}.
   *
   * @return the hash code value for this {@link FileColumn}
   * @see Object#equals(Object)
   * @see #equals(Object)
   */
  @Override
  int hashCode();

}
