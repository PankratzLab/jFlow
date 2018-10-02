package org.pankratzlab.common.parsing;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.pankratzlab.common.Files;

/**
 * Class for creating and configuring {@link FileParser} objects.
 */
public abstract class AbstractFileParserFactory {

  final String inputFile;
  final String inputFileDelim;

  int skippedLines = 0;
  Set<String> skipPrefices;

  boolean failBlank = false;
  boolean failCount = true;

  List<ColumnFilter> filters;
  Map<ColumnFilter, Boolean> filterDeath;

  boolean noHeader = false;
  boolean trimInput = true;

  String parseFailValue = ".";

  final List<FileLink> linkedParsers;
  final Map<FileLink, Map<FileColumn<?>, FileColumn<?>>> linkedFileColumns;
  final List<FileColumn<?>> dataInOrder;
  final List<FileColumn<?>> addlDataToLoad;
  final List<FileColumn<?>> optlDataToLoad;

  protected AbstractFileParserFactory(String inputFile, String inputDelim) {
    if (!Files.exists(inputFile)) {
      throw new IllegalArgumentException(new FileNotFoundException("File missing: " + inputFile));
    }
    this.inputFile = inputFile;
    this.inputFileDelim = inputDelim;

    skipPrefices = new HashSet<>();
    filters = new ArrayList<>();
    filterDeath = new HashMap<>();
    linkedParsers = new ArrayList<>();
    dataInOrder = new ArrayList<>();
    addlDataToLoad = new ArrayList<>();
    optlDataToLoad = new ArrayList<>();
    linkedFileColumns = new HashMap<>();
  }

  /**
   * Add FileColumns to load, iff they exist in the file.
   * 
   * @param columns FileColumn<?> columns to load
   * @return this
   */
  public AbstractFileParserFactory optionalColumns(FileColumn<?>... columns) {
    for (FileColumn<?> fc : columns) {
      this.optlDataToLoad.add(fc);
    }
    return this;
  }

  /**
   * Add a filter to drop lines, specifying whether or not to throw an exception if a value fails.
   * 
   * @param filter {@link ColumnFilter}
   * @param dieIfFail
   * @return this
   */
  public AbstractFileParserFactory filter(ColumnFilter filter, boolean dieIfFail) {
    if (!this.filters.contains(filter)) {
      this.filters.add(filter);
      this.filterDeath.put(filter, dieIfFail);
      for (FileColumn<?> filterCol : filter.getFilterColumns()) {
        if (!dataInOrder.contains(filterCol)) {
          addlDataToLoad.add(filterCol);
        }
      }
    }
    return this;
  }

  /**
   * Add a filter to drop lines, without dying if filtering rejects a line.
   * 
   * @param filter {@link ColumnFilter}
   * @return this
   */
  public AbstractFileParserFactory filter(ColumnFilter filter) {
    return filter(filter, false);
  }

  /**
   * Skip a specified number of lines at the beginning of the file.<br />
   * Multiple calls to this method will overwrite the previous value.
   * 
   * @param num Number of lines to skip
   * @return this
   */
  public AbstractFileParserFactory skipLines(int num) {
    this.skippedLines = num;
    return this;
  }

  /**
   * Skip all lines that start with the specified prefix.<br />
   * Can be called multiple times with different prefices.
   * 
   * @param prefix Line skip prefix
   * @return this
   */
  public AbstractFileParserFactory skipPrefix(String prefix) {
    this.skipPrefices.add(prefix);
    return this;
  }

  /**
   * Ignore lines with an incorrect number of columns. Default behavior is to throw an exception.
   * 
   * @return this
   */
  public AbstractFileParserFactory skipLinesWithIncorrectNumberOfColumns() {
    this.failCount = false;
    return this;
  }

  /**
   * Throw an exception if blank lines are encountered while reading the file. Default behavior is
   * to skip blank lines.
   * 
   * @return this
   */
  public AbstractFileParserFactory failOnBlankLines() {
    this.failBlank = true;
    return this;
  }

  /**
   * Flag the parser to not expect a header. May cause failure if data columns expect headers.
   * 
   * @return this
   */
  public AbstractFileParserFactory noHeader() {
    this.noHeader = true;
    return this;
  }

  /**
   * Flag the parser to not trim input lines.
   * 
   * @return this
   */
  public AbstractFileParserFactory noTrim() {
    this.trimInput = false;
    return this;
  }

  /**
   * Set the value to be used if parsing a value fails.
   * 
   * @param value String
   * @return this
   */
  public AbstractFileParserFactory parseFailValue(String value) {
    this.parseFailValue = value;
    return this;
  }

  /**
   * Check that the key columns in the file link share a name with the data loaded in the primary
   * parser.
   * 
   * @param link Link to check
   */
  private void checkLink(FileLink link) {
    for (FileColumn<?> fc : link.getKeys()) {
      boolean found = false;
      for (FileColumn<?> fc1 : dataInOrder) {
        if (fc1.getName().equals(fc.getName())) {
          found = true;
          break;
        }
      }
      for (FileColumn<?> fc1 : addlDataToLoad) {
        if (fc1.getName().equals(fc.getName())) {
          found = true;
          break;
        }
      }
      if (!found) {
        throw new IllegalArgumentException("Source data must load linked column " + fc.getName());
      }
    }
  }

  /**
   * Link another file to this one with a {@link FileLink}. Key columns, specified by
   * {@link FileLink#keys(FileColumn...)}, must share a name (though not necessarily the same header
   * value) with a column expected to be loaded (either directly or by a filter) in this file.
   * 
   * @param link {@link FileLink}
   * @return this
   */
  public AbstractFileParserFactory link(FileLink link) {
    checkLink(link);

    linkedParsers.add(link);
    linkedFileColumns.put(link, new HashMap<>());
    Set<FileColumn<?>> all = new HashSet<>();
    all.addAll(Collections.unmodifiableList(dataInOrder));
    all.addAll(Collections.unmodifiableList(addlDataToLoad));
    for (FileColumn<?> fc : link.getKeys()) {
      for (FileColumn<?> fc1 : all) {
        if (fc.getName().equals(fc1.getName())) {
          linkedFileColumns.get(link).put(fc, fc1);
          break;
        }
      }
    }
    return this;
  }

  public FileParser build() throws IOException {
    FileParser parser = new FileParser(this);
    parser.loadLinkedData();
    return parser;
  }

}
