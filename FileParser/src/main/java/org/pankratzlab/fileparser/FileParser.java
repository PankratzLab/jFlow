package org.pankratzlab.fileparser;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.Closeable;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.util.zip.ZipInputStream;
import java.util.zip.ZipOutputStream;

import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.xssf.usermodel.XSSFSheet;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.ImmutableSet.Builder;
import com.google.common.collect.Multiset;

public class FileParser implements Iterable<DataLine>, Closeable {

  private final String inputFile;
  private String inputFileDelim;

  private final int skippedLines;
  private final Set<String> skipPrefices;

  private final boolean failBlank;
  private final boolean failCount;

  private final List<ColumnFilter> filters;
  private final Map<ColumnFilter, Boolean> filterDeath;

  private final boolean noHeader;
  private final boolean trimInput;

  private final String parseFailValue;

  private final ImmutableSet<FileLink> linkedParsers;
  private final Map<FileLink, Map<FileColumn<?>, FileColumn<?>>> linkedFileColumns;
  protected final ImmutableSet<FileColumn<?>> dataInOrder;
  protected final ImmutableSet<FileColumn<?>> addlDataToLoad;
  protected final ImmutableSet<FileColumn<?>> optlDataToLoad;
  private ImmutableSet<FileColumn<?>> optlDataFound;

  private BufferedReader reader;
  private long lineCount = 0;
  private ImmutableMap<String, Integer> header;
  private Multiset<ColumnFilter> filterHits;

  /**
   * Don't use - create FileParsers with {@link FileParserFactory#setup(String, FileColumn...)}
   * instead.
   * 
   * @param inputFile String full path to file
   * @param inputDelim Input delimiter, or null if the delimiter should be determined from the
   *          header line
   * @throws IOException
   */
  FileParser(AbstractFileParserFactory factory) throws IOException {
    this.inputFile = factory.inputFile;
    this.inputFileDelim = factory.inputFileDelim;

    this.skippedLines = factory.skippedLines;
    skipPrefices = ImmutableSet.copyOf(factory.skipPrefices);

    this.failBlank = factory.failBlank;
    this.failCount = factory.failCount;

    filters = ImmutableList.copyOf(factory.filters);
    filterDeath = ImmutableMap.copyOf(factory.filterDeath);
    filterHits = HashMultiset.create();

    this.noHeader = factory.noHeader;
    this.trimInput = factory.trimInput;

    this.parseFailValue = factory.parseFailValue;

    linkedParsers = ImmutableSet.copyOf(factory.linkedParsers);
    linkedFileColumns = ImmutableMap.copyOf(factory.linkedFileColumns);
    dataInOrder = ImmutableSet.copyOf(factory.dataInOrder);
    addlDataToLoad = ImmutableSet.copyOf(factory.addlDataToLoad);
    optlDataToLoad = ImmutableSet.copyOf(factory.optlDataToLoad);
    open();
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
    cols.addAll(optlDataFound);
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
    reader = new BufferedReader(getAppropriateInputStream(inputFile));
    for (int i = 0; i < skippedLines; i++) {
      getNextLine();
    }
    if (!noHeader) {
      String headerLine = getNextLine();
      if (headerLine != null) {
        if (inputFileDelim == null) {
          inputFileDelim = determineDelimiter(headerLine);
        }
        String[] headerArray = headerLine.split(inputFileDelim);
        ImmutableMap.Builder<String, Integer> builder = new ImmutableMap.Builder<>();
        for (int i = 0; i < headerArray.length; i++) {
          builder.put(headerArray[i], i);
        }
        header = builder.build();
      }
    } else {
      try (BufferedReader delimReader = new BufferedReader(getAppropriateInputStream(inputFile))) {
        for (int i = 0; i < lineCount; i++)
          delimReader.readLine();
        String firstLine;
        do {
          firstLine = delimReader.readLine();
        } while (firstLine != null && matchesSkipPrefix(firstLine));
        if (firstLine != null) {
          if (trimInput) firstLine = firstLine.trim();
          if (inputFileDelim == null) inputFileDelim = determineDelimiter(firstLine);
          String[] firstLineCols = firstLine.split(inputFileDelim);
          ImmutableMap.Builder<String, Integer> headerBuilder = ImmutableMap.builderWithExpectedSize(firstLineCols.length);
          for (int i = 0; i < firstLineCols.length; i++) {
            headerBuilder.put(String.valueOf(i), i);
          }
          header = headerBuilder.build();
        }
      }
    }

    for (FileColumn<?> fc : dataInOrder) {
      fc.initialize(this);
    }
    for (FileColumn<?> fc : addlDataToLoad) {
      fc.initialize(this);
    }
    Builder<FileColumn<?>> b = ImmutableSet.builder();
    for (FileColumn<?> fc : optlDataToLoad) {
      try {
        fc.initialize(this);
        b.add(fc);
      } catch (IllegalStateException e) {
        // not in file
      }
    }
    optlDataFound = b.build();
  }

  private String getNextLine() throws IOException {
    String line;
    boolean skip;
    do {
      skip = false;
      line = reader.readLine();
      lineCount++;
      if (line == null) break;
      skip = matchesSkipPrefix(line);
      if (skip) {
        continue;
      }
      if (trimInput && line != null) {
        line = line.trim();
      }
    } while (skip);
    return line;
  }

  private boolean matchesSkipPrefix(String line) {
    if (!skipPrefices.isEmpty()) {
      for (String skp : skipPrefices) {
        if (line.startsWith(skp)) {
          return true;
        }
      }
    }
    return false;
  }

  /**
   * @return an {@link ImmutableMap} from header Strings to their index in the header with fixed
   *         iteration order by index
   */
  public ImmutableMap<String, Integer> getHeaderMap() {
    return header;
  }

  public int getFilteredCount(ColumnFilter cf) {
    return filterHits.count(cf);
  }

  private DataLine readLine() throws IOException {
    if (reader == null) return null;

    DataLine lineData;
    String[] parts;

    boolean skip;
    do {
      skip = false;
      parts = null;
      lineData = null;
      String line = getNextLine();
      if (line == null) break;
      if (failBlank && line.equals("")) {
        throw new IllegalStateException("Line " + lineCount + " was blank!");
      }

      lineData = new DataLine(parseFailValue);
      parts = line.split(inputFileDelim, -1);
      if (failCount && parts.length != header.size()) {
        throw new IllegalStateException("Line " + lineCount + " was the wrong length; expected "
                                        + header.size() + ", found " + parts.length);
      }

      for (ColumnFilter cf : filters) {
        for (FileColumn<?> fc : cf.getFilterColumns()) {
          if (!lineData.has(fc)) {
            lineData.parseOrFail(fc, parts);
          }
        }
        if (!cf.filter(lineData)) {
          if (filterDeath.get(cf)) {
            String colNames = "";
            List<FileColumn<?>> cols = cf.getFilterColumns();
            for (int f = 0; f < cols.size(); f++) {
              colNames += cols.get(f).getName();
              if (f < cols.size() - 1) {
                colNames += ", ";
              }
            }
            throw new IllegalStateException("Filter on columns " + colNames + " failed on line "
                                            + lineCount);
          }
          filterHits.add(cf);
          skip = true;
          break;
        }
        if (skip) break;
      }
      if (!skip) {
        for (FileColumn<?> fc : dataInOrder) {
          if (!lineData.has(fc)) {
            lineData.parseOrFail(fc, parts);
          }
        }
        for (FileColumn<?> fc : addlDataToLoad) {
          if (!lineData.has(fc)) {
            lineData.parseOrFail(fc, parts);
          }
        }
        for (FileColumn<?> fc : optlDataFound) {
          if (!lineData.has(fc)) {
            lineData.parseOrFail(fc, parts);
          }
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
    try (PrintWriter writer = openAppropriateWriter(outputFile, false)) {
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
    try (PrintWriter writer = openAppropriateWriter(outputFile, false)) {
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
    try (PrintWriter writer = openAppropriateWriter(outputFile, false)) {
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
  public int parseToFile(String outputFile, String outDelim) throws IOException {
    return this.parseToFile(outputFile, outDelim, null);
  }

  /**
   * {@link #parseToFile(String, String, boolean, boolean, List)} writing header and not appending
   */
  public int parseToFile(String outputFile, String outDelim,
                         List<FileColumn<?>> outputOrder) throws IOException {
    return this.parseToFile(outputFile, outDelim, true, false, outputOrder);
  }

  /**
   * {@link #parseToFile(String, String, boolean, boolean, List)} with default outputOrder
   */
  public int parseToFile(String outputFile, String outDelim, boolean writeHeader,
                         boolean append) throws IOException {
    return parseToFile(outputFile, outDelim, writeHeader, append, null);
  }

  /**
   * @param outputFile file to write to
   * @param outDelim delimiter to use in output file
   * @param writeHeader true to include header
   * @param append true to append to outputFile
   * @param outputOrder order to output columns in, any additional output columns will be included
   *          after listed columns
   * @return number of data lines written (does not include header, if applicable)
   * @throws IOException when thrown by {@link #close()}
   */
  public int parseToFile(String outputFile, String outDelim, boolean writeHeader, boolean append,
                         List<FileColumn<?>> outputOrder) throws IOException {
    int lines = 0;
    try (PrintWriter writer = openAppropriateWriter(outputFile, append)) {
      List<FileColumn<?>> outputColumns = buildOutputColumns(outputOrder);
      if (writeHeader) {
        writer.println(buildHeaderLineOut(outputColumns, outDelim));
      }

      Iterator<DataLine> iter = iterator();
      while (iter.hasNext()) {
        DataLine line = iter.next();
        writer.println(buildLineOut(line, outputColumns, outDelim));
        lines++;
      }
      close();
    }
    return lines;
  }

  public int parseToExcelSheet(XSSFSheet sheet, boolean writeHeader,
                               List<FileColumn<?>> outputOrder) {
    int lines = 0;

    List<FileColumn<?>> outputColumns = buildOutputColumns(outputOrder);

    int rowNum = 0;
    int colNum = 0;
    Row row;
    Cell cell;

    if (writeHeader) {
      row = sheet.createRow(rowNum++);

      for (int i = 0, count = outputColumns.size(); i < count; i++) {
        cell = row.createCell(colNum++);
        cell.setCellValue(outputColumns.get(i).getHeader());
      }
    }

    Iterator<DataLine> iter = iterator();
    while (iter.hasNext()) {
      DataLine line = iter.next();
      colNum = 0;
      row = sheet.createRow(rowNum++);

      for (int i = 0, count = outputColumns.size(); i < count; i++) {
        cell = row.createCell(colNum++);
        String val = line.getString(outputColumns.get(i));
        try {
          if (Double.isNaN(Double.parseDouble(val))) {
            cell.setCellValue(val);
          } else {
            cell.setCellValue(Double.parseDouble(val));
          }
        } catch (NumberFormatException nfe) {
          cell.setCellValue(val);
        }
      }

      lines++;
    }

    return lines;
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
    List<FileColumn<?>> outputColumns;
    if (outputOrder == null) {
      outputColumns = getOutputColumnsInOrder();
    } else {
      outputColumns = new ArrayList<>();
      Set<FileColumn<?>> allOutput = ImmutableSet.copyOf(getOutputColumnsInOrder());
      for (FileColumn<?> fc : outputOrder) {
        if (allOutput.contains(fc)) {
          outputColumns.add(fc);
        }
      }
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
      lineOut.append(outputColumns.get(i).getHeader());
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

  /**
   * @param filename to open for writing
   * @param append true to append to rather than replace the file if it exists
   * @return an appropriate {@link PrintWriter}
   */
  private PrintWriter openAppropriateWriter(String filename,
                                            boolean append) throws FileNotFoundException,
                                                            IOException {
    PrintWriter writer;

    if (filename.endsWith(".gz")) {
      writer = new PrintWriter(new GZIPOutputStream(new FileOutputStream(filename, append)));
    } else if (filename.endsWith(".zip")) {
      writer = new PrintWriter(new ZipOutputStream(new FileOutputStream(filename, append)));
    } else {
      writer = new PrintWriter(new BufferedWriter(new FileWriter(filename, append)));
    }

    return writer;
  }

  private InputStreamReader getAppropriateInputStream(String filename) throws FileNotFoundException,
                                                                       IOException {
    InputStream is;
    if (filename.endsWith(".gz")) {
      is = new GZIPInputStream(new FileInputStream(filename));
    } else if (filename.endsWith(".zip")) {
      is = new ZipInputStream(new FileInputStream(filename));
    } else {
      is = new FileInputStream(filename);
    }
    return new InputStreamReader(is);
  }

  private int countInstancesOf(String str, String pattern) {
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

  private String determineDelimiter(String str) {
    if (str.contains("\t")) {
      return "\t";
    } else if (countInstancesOf(str, ",") > countInstancesOf(str, " ")) {
      return ",";
    } else {
      return "[\\s]+";
    }
  }

}
