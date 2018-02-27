package org.genvisis.gwas.parsing;

import java.util.Map;
import java.util.Map.Entry;

/**
 * A {@link FileColumn<String>} subclass that searches the file header for a column name that
 * matches an array of specified aliases. Expects that the file will have a header, and will throw
 * an exception if it doesn't (in which case, {@link ExplicitIndexedFileColumn} should probably be
 * used).
 */
public class AliasedFileColumn extends AbstractFileColumn<String> implements IndexedFileColumn<String> {

  private final Aliases aliases;
  private final boolean aliasHeader;

  private int matchedIndex;
  private String matchedAlias = null;

  /**
   * {@link #AliasedFileColumn(String, Aliases)} with a strategy that fails if it finds multiple
   * matches, and isn't case-sensitive.
   * 
   * @param name String Column name
   * @param aliases String... Aliases to search for.
   */
  public AliasedFileColumn(String name, String... aliases) {
    this(name, new Aliases(aliases));
  }

  /**
   * Create an AliasedFileColumn with the given name that searches for the given aliases, following
   * the {@link Aliases#getStrategy()} value for multiple matches and
   * {@link Aliases#isCaseSensitive()}.
   * 
   * @param name String Column name
   * @param aliases {@link Aliases}
   */
  public AliasedFileColumn(String name, Aliases aliases) {
    super(name, false);
    this.aliases = aliases;
    this.aliasHeader = false;
  }

  /**
   * Create an AliasedFileColumn with the given name that searches for the given aliases, following
   * the {@link Aliases#getStrategy()} value for multiple matches and
   * {@link Aliases#isCaseSensitive()}. The output column header will be the matched input column
   * header and the name will be the first provided alias
   * 
   * @param aliases {@link Aliases}
   */
  public AliasedFileColumn(Aliases aliases) {
    super(aliases.getAliases()[0], false);
    this.aliases = aliases;
    this.aliasHeader = true;
  }

  @Override
  public void initialize(FileParser parser) {
    super.initialize(parser);
    if (matchedAlias != null) return;
    Map<String, Integer> headerMap = parser.getHeaderMap();
    if (headerMap == null) {
      throw new IllegalStateException("Column " + getName() + " expected a header for file "
                                      + parser.getInputFile());
    }

    int matched = aliases.getStrategy().match(headerMap, aliases.getAliases(),
                                              aliases.isCaseSensitive());

    if (matched == -1) {
      throw new IllegalStateException("Couldn't identify column header for " + getName()
                                      + " in file " + parser.getInputFile());
    }
    for (Entry<String, Integer> header : headerMap.entrySet()) {
      if (header.getValue() == matched) {
        this.matchedAlias = header.getKey();
        break;
      }
    }
    this.matchedIndex = matched;
  }

  @Override
  public String getValue(String[] line) throws ParseFailureException {
    return line[matchedIndex];
  }

  @Override
  public String getHeader() {
    if (aliasHeader) return matchedAlias;
    else return super.getHeader();
  }

  /**
   * @return The actual alias found in the file header
   */
  public String getMatchedAlias() {
    return matchedAlias;
  }

  /**
   * @return Index of column in data based on matched header
   */
  @Override
  public int getIndex() {
    return matchedIndex;
  }

}
