package org.genvisis.gwas.parsing;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.ext;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableSet;

public class FileParser implements Iterable<DataLine>, Closeable {

  private final String inputFile;
  private String inputFileDelim;

  private int skippedLines = 0;
  private Set<String> skipPrefices;

  private boolean failBlank = false;
  private boolean failCount = true;

  private List<ColumnFilter> filters;
  private Map<ColumnFilter, Boolean> filterDeath;

  private boolean noHeader = false;
  private boolean trimInput = true;

  private String parseFailValue = ".";

  private List<FileLink> linkedParsers;
  private Map<FileLink, Map<FileColumn<?>, FileColumn<?>>> linkedFileColumns;
  protected List<FileColumn<?>> dataInOrder;
  protected List<FileColumn<?>> addlDataToLoad;
  protected List<FileColumn<?>> optlDataToLoad;
  private List<FileColumn<?>> optlDataFound;

  private BufferedReader reader;
  private boolean opened = false;
  private long lineCount = 0;
  private ImmutableMap<String, Integer> header;

  /**
   * Don't use - create FileParsers with {@link FileParserFactory#setup(String, FileColumn...)}
   * instead.
   * 
   * @param inputFile String full path to file
   * @param inputDelim Input delimiter, or null if the delimiter should be determined from the
   *          header line
   */
  FileParser(AbstractFileParserFactory factory) {
    this.inputFile = factory.inputFile;
    this.inputFileDelim = factory.inputFileDelim;
    skipPrefices = new ImmutableSet.Builder<String>().addAll(factory.skipPrefices).build();
    filters = new ImmutableList.Builder<ColumnFilter>().addAll(factory.filters).build();
    filterDeath = new ImmutableMap.Builder<ColumnFilter, Boolean>().putAll(factory.filterDeath)
                                                                   .build();
    linkedParsers = new ImmutableList.Builder<FileLink>().addAll(factory.linkedParsers).build();
    dataInOrder = new ImmutableList.Builder<FileColumn<?>>().addAll(factory.dataInOrder).build();
    addlDataToLoad = new ImmutableList.Builder<FileColumn<?>>().addAll(factory.addlDataToLoad)
                                                               .build();
    optlDataToLoad = new ImmutableList.Builder<FileColumn<?>>().addAll(factory.optlDataToLoad)
                                                               .build();
    linkedFileColumns = new ImmutableMap.Builder<FileLink, Map<FileColumn<?>, FileColumn<?>>>().putAll(factory.linkedFileColumns)
                                                                                               .build();
  }

  /**
   * Calls {@link FileLink#build()} on each of the linked files.
   */
  void loadLinkedData() {
    for (FileLink fl : linkedParsers) {
      fl.build();
    }
  }

  private List<FileColumn<?>> getOutputColumnsInOrder() {
    List<FileColumn<?>> cols = new ArrayList<>(dataInOrder);
    for (FileLink fl : linkedParsers) {
      for (FileColumn<?> fc : fl.getValues()) {
        if (!cols.contains(fc)) {
          cols.add(fc);
        }
      }
    }
    return cols;
  }

  private void open() throws IOException {
    if (opened) {
      return;
    }
    reader = Files.getAppropriateReader(inputFile);
    String line = null;
    int skip = skippedLines;
    header = null;
    lineCount = 0;
    while ((line = reader.readLine()) != null) {
      lineCount++;
      if (trimInput) {
        line = line.trim();
      }
      boolean skipIt = false;
      if (skipPrefices.size() > 0) {
        for (String skp : skipPrefices) {
          if (line.startsWith(skp)) {
            skipIt = true;
            break;
          }
        }
      }
      if (skipIt) continue;
      if (skip > 0) {
        skip--;
        continue;
      }
      if (inputFileDelim == null) {
        inputFileDelim = ext.determineDelimiter(line);
      }
      if (!noHeader) {
        header = ArrayUtils.immutableIndexMap(line.split(inputFileDelim));
      }
      for (FileColumn<?> fc : dataInOrder) {
        fc.initialize(this);
      }
      for (FileColumn<?> fc : addlDataToLoad) {
        fc.initialize(this);
      }
      optlDataFound = new ArrayList<>();
      for (FileColumn<?> fc : optlDataToLoad) {
        try {
          fc.initialize(this);
          optlDataFound.add(fc);
        } catch (IllegalStateException e) {
          // not in file
        }
      }
      break;
    }

    opened = true;
  }

  /**
   * @return an {@link ImmutableMap} from header Strings to their index in the header with fixed
   *         iteration order by index
   */
  public ImmutableMap<String, Integer> getHeaderMap() {
    return header;
  }

  private DataLine readLine() throws IOException {
    if (reader == null) return null;

    DataLine lineData;
    String[] parts;

    boolean skip = false;
    do {
      skip = false;
      parts = null;
      lineData = null;

      String line = reader.readLine();
      lineCount++;
      if (line == null) break;
      if (trimInput) {
        line = line.trim();
      }
      if (failBlank && line.equals("")) {
        throw new IllegalStateException("Line " + lineCount + " was blank!");
      }
      if (skipPrefices.size() > 0) {
        for (String skp : skipPrefices) {
          if (line.startsWith(skp)) {
            skip = true;
            break;
          }
        }
      }
      if (skip) {
        continue;
      }
      lineData = new DataLine(parseFailValue);
      parts = line.split(inputFileDelim, -1);
      if (failCount && parts.length != header.size()) {
        throw new IllegalStateException("Line " + lineCount + " was the wrong length; expected "
                                        + header.size() + ", found " + parts.length);
      }
      for (FileColumn<?> fc : dataInOrder) {
        lineData.parseOrFail(fc, parts);
      }
      for (FileColumn<?> fc : addlDataToLoad) {
        lineData.parseOrFail(fc, parts);
      }
      for (FileColumn<?> fc : optlDataFound) {
        lineData.parseOrFail(fc, parts);
      }
      for (ColumnFilter fc : filters) {
        if (!fc.filter(lineData)) {
          if (filterDeath.get(fc)) {
            String colNames = "";
            List<FileColumn<?>> cols = fc.getFilterColumns();
            for (int f = 0; f < cols.size(); f++) {
              colNames += cols.get(f).getName();
              if (f < cols.size() - 1) {
                colNames += ", ";
              }
            }
            throw new IllegalStateException("Filter on columns " + colNames + " failed on line "
                                            + lineCount);
          }
          skip = true;
          break;
        }
      }
    } while (skip);

    if (lineData != null) {
      for (FileLink fl : linkedParsers) {
        ImmutableList.Builder<Object> builder = ImmutableList.builder();
        for (FileColumn<?> fc : fl.getKeys()) {
          builder.add(lineData.getString(linkedFileColumns.get(fl).get(fc)));
        }
        ImmutableList<Object> key = builder.build();
        DataLine linkedLine = fl.get(key);
        if (linkedLine == null) {
          // missing data
          for (FileColumn<?> valueCol : fl.getValues()) {
            // check if data already exists
            if (!lineData.has(valueCol)) {
              lineData.fail(valueCol);
            }
          }
        } else {
          for (FileColumn<?> valueCol : fl.getValues()) {
            // check if data doesn't already exist or was a failure
            if (!lineData.hasValid(valueCol) && linkedLine.hasValid(valueCol)) {
              lineData.copyUnsafe(valueCol, linkedLine);
            } else if (lineData.hasValid(valueCol) && linkedLine.hasValid(valueCol)
                       && !linkedLine.getUnsafe(valueCol).equals(lineData.getUnsafe(valueCol))) {
              // if a value already exists and isn't the same nor a parse failure
              throw new IllegalStateException("Different value found for column "
                                              + valueCol.getName() + " with linked key; key="
                                              + key.toString() + ", value1="
                                              + lineData.getString(valueCol) + ", value2="
                                              + linkedLine.getString(valueCol));
            } else {
              lineData.fail(valueCol);
            }
          }
        }
      }
    }

    return lineData;
  }

  /**
   * @return The input file path
   */
  public String getInputFile() {
    return inputFile;
  }

  /**
   * Close the internal {@link BufferedReader}.
   * 
   * @throws IOException
   */
  @Override
  public void close() throws IOException {
    reader.close();
  }

  @Override
  public Iterator<DataLine> iterator() {
    try {
      open();
    } catch (IOException e1) {
      throw new RuntimeException(e1);
    }
    return new Iterator<DataLine>() {

      boolean started = false;
      DataLine currentLine = null;

      @Override
      public DataLine next() {
        DataLine line = null;
        try {
          if (!started) {
            // if we're calling next() without having called hasNext() first
            currentLine = readLine();
            started = true;
          }
          line = currentLine;
          // read the next line, which will return null when we're done
          currentLine = readLine();
        } catch (IOException e) {
          throw new RuntimeException(e);
        }
        return line;
      }

      @Override
      public boolean hasNext() {
        if (!started) {
          /*
           * if we haven't started reading, we don't know if we have or will have any data (due to
           * filters / skipped lines), so load the first line of data.
           */
          try {
            currentLine = readLine();
          } catch (IOException e) {
            throw new RuntimeException(e);
          }
          started = true;
        }
        return currentLine != null;
      }
    };
  }

  /**
   * @return {@link List} of {@link DataLine}s for each line in the input
   * @throws IOException when thrown by {@link #close()}
   */
  public List<DataLine> load() throws IOException {
    List<DataLine> data = new ArrayList<>();
    Iterator<DataLine> iter = iterator();
    while (iter.hasNext()) {
      data.add(iter.next());
    }
    close();
    return data;
  }

  /**
   * @param dropIfKeyFail Set to true to drop any lines where the key value(s) throw a
   *          ParseFailureException; set to false to throw a RuntimeException if invalid keys are
   *          found.
   * @param keyCols Vararg of FileColumns whose String values are to be concatenated together and
   *          used as the key
   * @return {@link Map} from {@link ImmutableList} of keyCols values to {@link DataLine}
   * @throws IOException when thrown by {@link #close()}
   */
  public Map<ImmutableList<Object>, DataLine> load(boolean dropIfKeyFail,
                                                   FileColumn<?>... keyCols) throws IOException {
    Map<ImmutableList<Object>, DataLine> data = new HashMap<>();
    Iterator<DataLine> iter = iterator();
    while (iter.hasNext()) {
      DataLine line = iter.next();
      try {
        data.put(buildParsedKey(line, keyCols), line);
      } catch (ParseFailureException e) {
        if (!dropIfKeyFail) {
          throw new RuntimeException(e);
        }
      }
    }
    close();
    return data;
  }

  /**
   * @param dropIfKeyFail Set to true to drop any lines where the key value throws a
   *          ParseFailureException; set to false to throw a RuntimeException if invalid keys are
   *          found.
   * @param keyCol FileColumn whose String values are to be used as the key
   * @return {@link Map} from value of keyCol to {@link DataLine}
   * @throws IOException when thrown by {@link #close()}
   */
  public <T> Map<T, DataLine> load(boolean dropIfKeyFail, FileColumn<T> keyCol) throws IOException {
    Map<T, DataLine> data = new HashMap<>();
    Iterator<DataLine> iter = iterator();
    while (iter.hasNext()) {
      DataLine line = iter.next();
      try {
        data.put(line.get(keyCol), line);
      } catch (ParseFailureException e) {
        if (!dropIfKeyFail) {
          throw new RuntimeException(e);
        }
      }
    }
    close();
    return data;
  }

  /**
   * Simultaneously parse to file as in {@link #parseToFile(String, String)} and load as in
   * {@link #load()}
   */
  public List<DataLine> parseToFileAndLoad(String outputFile, String outDelim) throws IOException {
    return parseToFileAndLoad(outputFile, outDelim, null);
  }

  public List<DataLine> parseToFileAndLoad(String outputFile, String outDelim,
                                           List<FileColumn<?>> outputOrder) throws IOException {
    try (PrintWriter writer = Files.getAppropriateWriter(outputFile)) {
      List<FileColumn<?>> outputColumns = buildOutputColumns(outputOrder);
      writer.println(buildHeaderLineOut(outputColumns, outDelim));

      List<DataLine> data = new ArrayList<>();
      Iterator<DataLine> iter = iterator();
      while (iter.hasNext()) {
        DataLine line = iter.next();
        data.add(line);
        writer.println(buildLineOut(line, outputColumns, outDelim));
      }
      close();
      return data;
    }
  }

  /**
   * {@link #parseToFileAndLoad(String, String, List, boolean, FileColumn)} with default output
   * order
   */
  public <T> Map<T, DataLine> parseToFileAndLoad(String outputFile, String outDelim,
                                                 boolean dropIfKeyFail,
                                                 FileColumn<T> keyCol) throws IOException {
    return parseToFileAndLoad(outputFile, outDelim, null, dropIfKeyFail, keyCol);
  }

  /**
   * Simultaneously parse to file as in {@link #parseToFile(String, String, List)} and load as in
   * {@link #load(boolean, FileColumn)}
   */
  public <T> Map<T, DataLine> parseToFileAndLoad(String outputFile, String outDelim,
                                                 List<FileColumn<?>> outputOrder,
                                                 boolean dropIfKeyFail,
                                                 FileColumn<T> keyCol) throws IOException {
    try (PrintWriter writer = Files.getAppropriateWriter(outputFile)) {
      List<FileColumn<?>> outputColumns = buildOutputColumns(outputOrder);
      writer.println(buildHeaderLineOut(outputColumns, outDelim));

      Map<T, DataLine> data = new HashMap<>();
      Iterator<DataLine> iter = iterator();
      while (iter.hasNext()) {
        DataLine line = iter.next();
        try {
          data.put(line.get(keyCol), line);
          writer.println(buildLineOut(line, outputColumns, outDelim));
        } catch (ParseFailureException e) {
          if (!dropIfKeyFail) {
            throw new RuntimeException(e);
          }
        }
      }
      close();
      return data;
    }
  }

  /**
   * {@link #parseToFileAndLoad(String, String, List, boolean, FileColumn...)} with default output
   * order
   */
  public Map<ImmutableList<Object>, DataLine> parseToFileAndLoad(String outputFile, String outDelim,
                                                                 boolean dropIfKeyFail,
                                                                 FileColumn<?>... keyCols) throws IOException {
    return parseToFileAndLoad(outputFile, outDelim, null, dropIfKeyFail, keyCols);
  }

  /**
   * Simultaneously parse to file as in {@link #parseToFile(String, String, List)} and load as in
   * {@link #load(boolean, FileColumn...)}
   */
  public Map<ImmutableList<Object>, DataLine> parseToFileAndLoad(String outputFile, String outDelim,
                                                                 List<FileColumn<?>> outputOrder,
                                                                 boolean dropIfKeyFail,
                                                                 FileColumn<?>... keyCols) throws IOException {
    try (PrintWriter writer = Files.getAppropriateWriter(outputFile)) {
      List<FileColumn<?>> outputColumns = buildOutputColumns(outputOrder);
      writer.println(buildHeaderLineOut(outputColumns, outDelim));

      Map<ImmutableList<Object>, DataLine> data = new HashMap<>();
      Iterator<DataLine> iter = iterator();
      while (iter.hasNext()) {
        DataLine line = iter.next();
        try {
          data.put(buildParsedKey(line, keyCols), line);
          writer.println(buildLineOut(line, outputColumns, outDelim));
        } catch (ParseFailureException e) {
          if (!dropIfKeyFail) {
            throw new RuntimeException(e);
          }
        }
      }
      close();
      return data;
    }
  }

  /**
   * {@link #parseToFile(String, String, List)} using default output orser
   */
  public void parseToFile(String outputFile, String outDelim) throws IOException {
    this.parseToFile(outputFile, outDelim, null);
  }

  /**
   * {@link #parseToFile(String, String, boolean, boolean, List)} writing header and not appending
   */
  public void parseToFile(String outputFile, String outDelim,
                          List<FileColumn<?>> outputOrder) throws IOException {
    this.parseToFile(outputFile, outDelim, true, false, outputOrder);
  }

  /**
   * {@link #parseToFile(String, String, boolean, boolean, List)} with default outputOrder
   */
  public void parseToFile(String outputFile, String outDelim, boolean writeHeader,
                          boolean append) throws IOException {
    parseToFile(outputFile, outDelim, writeHeader, append, null);
  }

  /**
   * @param outputFile file to write to
   * @param outDelim delimiter to use in output file
   * @param writeHeader true to include header
   * @param append true to append to outputFile
   * @param outputOrder order to output columns in, any additional output columns will be included
   *          after listed columns
   * @throws IOException when thrown by {@link #close()}
   */
  public void parseToFile(String outputFile, String outDelim, boolean writeHeader, boolean append,
                          List<FileColumn<?>> outputOrder) throws IOException {
    try (PrintWriter writer = Files.getAppropriateWriter(outputFile, append)) {
      List<FileColumn<?>> outputColumns = buildOutputColumns(outputOrder);
      if (writeHeader) {
        writer.println(buildHeaderLineOut(outputColumns, outDelim));
      }

      Iterator<DataLine> iter = iterator();
      while (iter.hasNext()) {
        DataLine line = iter.next();
        writer.println(buildLineOut(line, outputColumns, outDelim));
      }
      close();
    }
  }

  private ImmutableList<Object> buildParsedKey(DataLine line,
                                               FileColumn<?>... keyCols) throws ParseFailureException {
    ImmutableList.Builder<Object> builder = ImmutableList.builder();
    for (FileColumn<?> fc : keyCols) {
      builder.add(line.get(fc));
    }
    return builder.build();
  }

  private List<FileColumn<?>> buildOutputColumns(List<FileColumn<?>> outputOrder) {
    List<FileColumn<?>> outputColumns = outputOrder == null ? getOutputColumnsInOrder()
                                                            : outputOrder;
    if (outputOrder != null) {
      List<FileColumn<?>> allOutput = getOutputColumnsInOrder();
      if (!outputColumns.containsAll(allOutput)) {
        for (FileColumn<?> fc : allOutput) {
          if (!outputColumns.contains(fc)) {
            outputColumns.add(fc);
          }
        }
      }
    }
    return outputColumns;
  }

  private String buildHeaderLineOut(List<FileColumn<?>> outputColumns, String outDelim) {
    StringBuilder lineOut = new StringBuilder();
    for (int i = 0, count = outputColumns.size(); i < count; i++) {
      lineOut.append(outputColumns.get(i).getName());
      if (i < count - 1) {
        lineOut.append(outDelim);
      }
    }
    return lineOut.toString();
  }

  private String buildLineOut(DataLine line, List<FileColumn<?>> outputColumns, String outDelim) {
    StringBuilder lineOut = new StringBuilder();
    for (int i = 0, count = outputColumns.size(); i < count; i++) {
      lineOut.append(line.getString(outputColumns.get(i)));
      if (i < count - 1) {
        lineOut.append(outDelim);
      }
    }
    return lineOut.toString();
  }

}
