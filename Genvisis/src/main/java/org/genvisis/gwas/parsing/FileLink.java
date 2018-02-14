package org.genvisis.gwas.parsing;

import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;

/**
 * Subclass of {@link AbstractFileParserFactory} to take advantage of {@link FileParser}
 * construction methods.
 */
public class FileLink extends AbstractFileParserFactory {

  private ImmutableSet<FileColumn<?>> keyColumns;
  private List<FileColumn<?>> valueColumns;
  private boolean dieOnMissing = false;
  private boolean dieOnInvalidKey = true;
  Map<ImmutableList<Object>, DataLine> data;

  private FileLink(String inputFile, String inputDelim) {
    super(inputFile, inputDelim);
  }

  /**
   * Get {@link FileColumn} objects this FileLink is keyed on.
   * 
   * @return ImmutableSet&lt;FileColumn&lt;?>>
   */
  public ImmutableSet<FileColumn<?>> getKeys() {
    return ImmutableSet.copyOf(keyColumns);
  }

  /**
   * Set the {@link FileColumn} objects that this FileLink is keyed on.
   * 
   * @param keyColumns FileColumn&lt;?>...
   * @return this
   */
  public FileLink keys(FileColumn<?>... keyColumns) {
    if (this.keyColumns != null) {
      throw new IllegalStateException("Keys have already been defined for linked file "
                                      + this.parser.getInputFile());
    }
    this.keyColumns = ImmutableSet.copyOf(keyColumns);
    for (FileColumn<?> fc : keyColumns) {
      this.parser.dataInOrder.add(fc);
    }
    return this;
  }

  /**
   * @return An unmodifiable view of the value columns.
   */
  public List<FileColumn<?>> getValues() {
    return Collections.unmodifiableList(valueColumns);
  }

  /**
   * Set the {@link FileColumn} objects that this FileLink will load from the linked file.
   * 
   * @param valueColumns FileColumn&lt;?>...
   * @return this
   */
  public FileLink values(FileColumn<?>... valueColumns) {
    if (this.valueColumns != null) {
      throw new IllegalStateException("Values have already been defined for linked file "
                                      + this.parser.getInputFile());
    }
    this.valueColumns = Arrays.asList(valueColumns);
    for (FileColumn<?> fc : valueColumns) {
      this.parser.dataInOrder.add(fc);
    }
    return this;
  }

  /**
   * Tell this FileLink to throw an exception if data for a key is missing.
   * 
   * @return this
   */
  public FileLink dieOnMissing() {
    dieOnMissing = true;
    return this;
  }

  /**
   * Tell this FileLink to NOT throw an exception if it encounters an invalid key when loading data.
   * 
   * @return this
   */
  public FileLink dropInvalidKeys() {
    dieOnInvalidKey = false;
    return this;
  }

  @Override
  /**
   * Builds parser internally with {@code super.build()}, loads the data in the linked file into
   * memory, and then closes the internal parser object. Always returns null;
   */
  public FileParser build() {
    if (data == null) {
      FileParser parser = super.build();
      try {
        data = parser.load(!dieOnInvalidKey,
                           getKeys().toArray(new FileColumn<?>[getKeys().size()]));
        parser.close();
      } catch (IOException e) {
        throw new RuntimeException(e);
      }
    }
    return null;
  }

  /**
   * Create a FileLink object on the given file using the given delimiter to split lines.
   * 
   * @param inputFile String full path to input file
   * @param inputDelim Delimiter to use to split input lines
   * @return new FileLink object
   */
  public static FileLink setup(String inputFile, String inputDelim) {
    return new FileLink(inputFile, inputDelim);
  }

  /**
   * Create a FileLink object on the given file.
   * 
   * @param inputFile String full path to input file
   * @return new FileLink object
   */
  public static FileLink setup(String inputFile) {
    return new FileLink(inputFile, null);
  }

  public DataLine get(List<Object> keyVal) {
    if (data == null) {
      throw new IllegalStateException("FileLink for file " + parser.getInputFile()
                                      + " needs to be built with the build() method.");
    }
    if (!data.containsKey(keyVal) && dieOnMissing) {
      throw new IllegalStateException("Failed to find linked data in file " + parser.getInputFile()
                                      + " for key " + keyVal);
    }
    return data.get(keyVal);
  }

  public String getInputFile() {
    return parser.getInputFile();
  }

}
