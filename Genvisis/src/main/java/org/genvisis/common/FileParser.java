package org.genvisis.common;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;

/**
 * parser to extract a string array based on a separator
 *
 */
public class FileParser extends BufferedReader {

  public static class Builder<T extends Builder<T>> {
    private Logger log = null;
    private String separator = "\t";

    public FileParser build(Reader in) {
      return new FileParser(this, in);
    }

    @SuppressWarnings("unchecked")
    public T log(Logger log) {
      this.log = log;
      return (T) this;
    }

    @SuppressWarnings("unchecked")
    public T separator(String separator) {
      this.separator = separator;
      return (T) this;
    }
  }

  protected String separator;

  protected Logger log;

  protected FileParser(@SuppressWarnings("rawtypes") Builder builder, Reader in) {
    super(in);
    log = builder.log == null ? new Logger() : builder.log;
    separator = builder.separator;
  }

  public String getSeparator() {
    return separator;
  }

  public String[] readLineArray() throws IOException {
    return readLine().trim().split(separator);
  }

}
