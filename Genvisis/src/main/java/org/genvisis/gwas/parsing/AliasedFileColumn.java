package org.genvisis.gwas.parsing;

import java.util.Map;
import java.util.Map.Entry;

/**
 * A {@link FileColumn<String>} subclass that searches the file header for a column name that
 * matches an array of specified aliases. Expects that the file will have a header, and will throw
 * an exception if it doesn't (in which case, {@link IndexedFileColumn} should probably be used).
 */
public class AliasedFileColumn extends AbstractFileColumn<String> {

  private final Aliases aliases;

  private int matchedIndex;
  private String matchedAlias = null;

  /**
   * Create an AliasedFileColumn with the given name that searches for the given aliases, fails if
   * it finds multiple matches, and isn't case-sensitive.
   * 
   * @param name String Column name
   * @param aliases String... Aliases to search for.
   */
  public AliasedFileColumn(String name, String... aliases) {
    super(name, false);
    this.aliases = new Aliases(aliases);
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
  }

  @Override
  public void initialize(Map<String, Integer> headerMap) {
    if (matchedAlias != null) return;
    if (headerMap == null) {
      throw new IllegalStateException("Column " + getName() + " expected a header!");
    }

    int matched = aliases.getStrategy().match(headerMap, aliases.getAliases(),
                                              aliases.isCaseSensitive());

    if (matched == -1) {
      throw new IllegalStateException("Couldn't identify column header for " + getName()
                                      + "; header map: " + headerMap.toString());
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

  /**
   * @return The actual alias found in the file header
   */
  public String getMatchedAlias() {
    return matchedAlias;
  }

  /**
   * @return Index of column in data based on matched header
   */
  public int getMatchedIndex() {
    return matchedIndex;
  }

  @Override
  public int hashCode() {
    final int prime = 31;
    int result = 1;
    result = prime * result + ((aliases == null) ? 0 : aliases.hashCode());
    // can't use mutable properties in hashCode, otherwise the hash could change during use
    return result;
  }

  @Override
  public boolean equals(Object obj) {
    if (this == obj) return true;
    if (obj == null) return false;
    if (getClass() != obj.getClass()) return false;
    AliasedFileColumn other = (AliasedFileColumn) obj;
    if (aliases == null) {
      if (other.aliases != null) return false;
    } else if (!aliases.equals(other.aliases)) return false;
    if (!getName().equals(other.getName())) return false;
    // can't use mutable properties in equals, otherwise equality could change during use
    return true;
  }

}
