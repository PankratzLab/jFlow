package org.genvisis.gwas.parsing;

/**
 * Factory for creating {@link FileParser} objects.
 */
public class FileParserFactory extends AbstractFileParserFactory {

  protected FileParserFactory(String inputFile, String inputDelim,
                              FileColumn<?>... dataToLoadInOrder) {
    super(inputFile, inputDelim);
    for (FileColumn<?> fc : dataToLoadInOrder) {
      this.dataInOrder.add(fc);
    }
  }

  /**
   * Create a {@link FileParserFactory} to parse the given file and pull out the data requested.
   * <br />
   * <br />
   * {@link FileColumn} objects should <b>not</b> be shared across instances of
   * {@link FileParserFactory}, nor shared with {@link FileLink} objects. Doing so may lead to
   * unexpected behavior (for example, an {@link AliasedFileColumn} shared across instances of
   * {@link FileParserFactory} or {@link FileLink} would redetermine the index of it's aliases for
   * each file). To link {@link FileColumn} instances across objects, give them the same name.
   * <br />
   * <br />
   * The {@code dataToLoadInOrder} argument reflects the desired data, and may not reflect the
   * entirety of what data is loaded in the background. Particularly, {@link ColumnFilter} objects
   * may specify needing other data, which will be loaded and filtered upon, but won't be reflected
   * in the final output data from {@link FileParser#parseToFile(String, String)},
   * {@link FileParser#load()}, or other similar methods.
   * 
   * @param inputFile String full path to file
   * @param dataToLoadInOrder {@link FileColumn}... Desired columns of data
   * @return {@link FileParserFactory}
   */
  public static FileParserFactory setup(String inputFile, FileColumn<?>... dataToLoadInOrder) {
    return setup(inputFile, null, dataToLoadInOrder);
  }

  /**
   * Create a {@link FileParserFactory} to parse the given file and pull out the data requested.
   * <br />
   * <br />
   * {@link FileColumn} objects should <b>not</b> be shared across instances of
   * {@link FileParserFactory}, nor shared with {@link FileLink} objects. Doing so may lead to
   * unexpected behavior (for example, an {@link AliasedFileColumn} shared across instances of
   * {@link FileParserFactory} or {@link FileLink} would redetermine the index of it's aliases for
   * each file). To link {@link FileColumn} instances across objects, give them the same name.
   * <br />
   * <br />
   * The {@code dataToLoadInOrder} argument reflects the desired data, and may not reflect the
   * entirety of what data is loaded in the background. Particularly, {@link ColumnFilter} objects
   * may specify needing other data, which will be loaded and filtered upon, but won't be reflected
   * in the final output data from {@link FileParser#parseToFile(String, String)},
   * {@link FileParser#load()}, or other similar methods.
   * 
   * @param inputFile String full path to file
   * @param inputDelim Delimiter of input file
   * @param dataToLoadInOrder {@link FileColumn}... Desired columns of data
   * @return {@link FileParserFactory}
   */
  public static FileParserFactory setup(String inputFile, String inputDelim,
                                        FileColumn<?>... dataToLoadInOrder) {
    return new FileParserFactory(inputFile, inputDelim, dataToLoadInOrder);
  }

}
