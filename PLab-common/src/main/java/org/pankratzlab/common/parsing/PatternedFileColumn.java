package org.pankratzlab.common.parsing;

import java.util.Map;
import java.util.regex.Pattern;

public class PatternedFileColumn extends AbstractFileColumn<String>
                                 implements IndexedFileColumn<String> {

  private final String regex;
  private int matchedIndex;
  private String matchedHeader = null;

  public PatternedFileColumn(String nm, String regex) {
    super(nm);
    this.regex = regex;
  }

  @Override
  public void initialize(FileParser parser) {
    super.initialize(parser);
    if (matchedHeader != null) return;
    Map<String, Integer> headerMap = parser.getHeaderMap();
    if (headerMap == null) {
      throw new IllegalStateException("Column " + getName() + " expected a header for file "
                                      + parser.getInputFile());
    }

    int matched = -1;
    String colname = "";

    for (String key : headerMap.keySet()) {
      if (Pattern.matches(regex, key)) {
        if (matched > -1) {
          // pick the earliest index we encounter
          colname = headerMap.get(key) < matched ? key : colname;
          matched = Math.min(matched, headerMap.get(key));
        } else {
          matched = headerMap.get(key);
          colname = key;
        }
      }
    }

    if (matched == -1) {
      throw new IllegalStateException("Couldn't identify column header for " + getName()
                                      + " in file " + parser.getInputFile());
    }

    this.matchedHeader = colname;
    this.matchedIndex = matched;
  }

  @Override
  public String getValue(String[] line) throws ParseFailureException {
    return line[matchedIndex];
  }

  /**
   * @return Index of column in data based on matched header
   */
  @Override
  public int getIndex() {
    return matchedIndex;
  }

  @Override
  public String getHeader() {
    return matchedHeader;
  }

}
