package org.pankratzlab.common.parsing;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

public class Aliases {

  public enum MultipleAliasStrategy {
    /**
     * Throws an exception if aliases match against multiple columns.
     */
    FAIL() {

      @Override
      int match(Map<String, Integer> headerMap, String[] aliases, boolean caseSensitive) {
        Set<String> aliasSet = new HashSet<>();
        for (String a : aliases) {
          aliasSet.add(caseSensitive ? a : a.toLowerCase());
        }
        Map<String, String> heads = new HashMap<>();
        for (String head : headerMap.keySet()) {
          heads.put(caseSensitive ? head : head.toLowerCase(), head);
        }
        Set<String> headSet = heads.keySet();
        headSet.retainAll(aliasSet);
        if (heads.size() > 1) {
          // TODO improve error message
          throw new IllegalStateException("Multiple matches found for aliases.");
        }
        for (String h : headSet) {
          return headerMap.get(heads.get(h));
        }
        return -1;
      }
    },
    /**
     * Use the first alias that matches any header value.
     */
    FIRST_ALIAS() {

      @Override
      int match(Map<String, Integer> headerMap, String[] aliases, boolean caseSensitive) {
        Map<String, String> heads = new HashMap<>();
        for (String head : headerMap.keySet()) {
          heads.put(caseSensitive ? head : head.toLowerCase(), head);
        }
        for (int i = 0; i < aliases.length; i++) {
          String key = caseSensitive ? aliases[i] : aliases[i].toLowerCase();
          if (heads.containsKey(key)) {
            return headerMap.get(heads.get(key));
          }
        }
        return -1;
      }
    },
    /**
     * Use the first header value that matches any alias.
     */
    FIRST_MATCH() {

      @Override
      int match(Map<String, Integer> headerMap, String[] aliases, boolean caseSensitive) {
        Set<String> aliasSet = new HashSet<>();
        for (String a : aliases) {
          aliasSet.add(caseSensitive ? a : a.toLowerCase());
        }
        String[] head = new String[headerMap.size()];
        for (Entry<String, Integer> entry : headerMap.entrySet()) {
          head[entry.getValue()] = entry.getKey();
        }
        for (int i = 0; i < head.length; i++) {
          if (aliasSet.contains(caseSensitive ? head[i] : head[i].toLowerCase())) {
            return i;
          }
        }
        return -1;
      }
    };

    abstract int match(Map<String, Integer> headerMap, String[] aliases, boolean caseSensitive);

  }

  private final String[] aliases;
  private final Aliases.MultipleAliasStrategy strategy;
  private final boolean caseSensitive;

  /**
   * Specify the String aliases, default the strategy to {@link MultipleAliasStrategy#FAIL} and
   * non-case-sensitive.
   * 
   * @param aliases String[]
   */
  public Aliases(String[] aliases) {
    this(aliases, MultipleAliasStrategy.FAIL, false);
  }

  /**
   * @param aliases {@code String[]} - List of String aliases for which to look.
   * @param strategy {@link MultipleAliasStrategy} - strategy to use in case of multiple matches.
   * @param caseSensitive {@code boolean} - whether the comparison should be case sensitive or not.
   */
  public Aliases(String[] aliases, Aliases.MultipleAliasStrategy strategy, boolean caseSensitive) {
    this.aliases = aliases;
    this.strategy = strategy;
    this.caseSensitive = caseSensitive;
  }

  public String[] getAliases() {
    return aliases;
  }

  public Aliases.MultipleAliasStrategy getStrategy() {
    return strategy;
  }

  public boolean isCaseSensitive() {
    return caseSensitive;
  }

  @Override
  public int hashCode() {
    final int prime = 31;
    int result = 1;
    result = prime * result + Arrays.hashCode(aliases);
    result = prime * result + (caseSensitive ? 1231 : 1237);
    result = prime * result + ((strategy == null) ? 0 : strategy.hashCode());
    return result;
  }

  @Override
  public boolean equals(Object obj) {
    if (this == obj) return true;
    if (obj == null) return false;
    if (getClass() != obj.getClass()) return false;
    Aliases other = (Aliases) obj;
    if (!Arrays.equals(aliases, other.aliases)) return false;
    if (caseSensitive != other.caseSensitive) return false;
    if (strategy != other.strategy) return false;
    return true;
  }

}
