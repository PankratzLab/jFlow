package org.pankratzlab.common.parsing;

import java.util.Arrays;
import java.util.Collection;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import org.pankratzlab.common.Positions;
import com.google.common.base.Joiner;
import com.google.common.base.Predicates;
import com.google.common.collect.Maps;

/**
 * Commonly used FileColumns. If you end up creating a FileColumn more than once, consider putting
 * it here!
 */
public final class StandardFileColumns {

  private StandardFileColumns() {}

  /**
   * @see #allExcept(String, Collection)
   */
  public static final FileColumn<String> allExcept(String outputDelimeter,
                                                   IndexedFileColumn<?>... excludes) {
    return allExcept(outputDelimeter, Arrays.asList(excludes));
  }

  /**
   * A special {@link FileColumn} that returns all values in a line except the values in the
   * excluded columns.
   * 
   * @param outputDelimeter Delimiter with which to join the values and headers
   * @param excludes {@link IndexedFileColumn}s to exclude
   */
  public static final FileColumn<String> allExcept(String outputDelimeter,
                                                   Collection<IndexedFileColumn<?>> excludes) {
    return new AbstractFileColumn<String>("allExcept_"
                                          + excludes.stream().map(FileColumn::getName)
                                                    .collect(Collectors.joining(outputDelimeter))) {

      private Map<String, Integer> subSetHeaderMap;
      private final String outDelim = outputDelimeter;

      @Override
      public void initialize(FileParser parser) {
        excludes.forEach(e -> e.initialize(parser));
        Map<String, Integer> headerMap = parser.getHeaderMap();
        Set<Integer> excludedIndices = excludes.stream().map(IndexedFileColumn::getIndex)
                                               .collect(Collectors.toSet());
        subSetHeaderMap = Maps.filterValues(headerMap, Predicates.not(excludedIndices::contains));
      }

      @Override
      public String getValue(String[] line) throws ParseFailureException {
        return subSetHeaderMap.values().stream().map(i -> line[i])
                              .collect(Collectors.joining(outDelim));
      }

      @Override
      public String getHeader() {
        return Joiner.on(outDelim).join(subSetHeaderMap.keySet());
      }

    };
  }

  /**
   * Column that computes log10 on a double.
   * 
   * @param name Output name
   * @param val {@link FileColumn&lt;Double&gt;}
   * @return
   */
  public static final FileColumn<Double> log10(final String name, final FileColumn<Double> val) {
    return new CachedFileColumn<Double>(name, val.dieOnParseFailure()) {

      @Override
      public void initialize(FileParser parser) {
        val.initialize(parser);
      }

      @Override
      public Double calculateValue(String[] line) throws ParseFailureException {
        return Math.log10(val.getValue(line));
      }
    };
  }

  /**
   * Creates a non-case-sensitive {@link AliasedFileColumn} around
   * {@link org.pankratzlab.common.Aliases#MARKER_NAMES} that will fail if multiple aliases are found.
   * 
   * @param colName Desired name of column
   * @return AliasedFileColumn
   */
  public static AliasedFileColumn snp(String colName) {
    return new AliasedFileColumn(colName, new Aliases(org.pankratzlab.common.Aliases.MARKER_NAMES));
  }

  /**
   * {@link #chr(String, boolean)} that will not die on parse failure
   */
  public static IndexedFileColumn<Byte> chr(String colName) {
    return chr(colName, false);
  }

  /**
   * Creates a non-case-sensitive {@link NumberWrapperColumn} that uses
   * {@link Positions#chromosomeNumber(String)} to determine the byte value of a chromosome wrapped
   * around an {@link AliasedFileColumn} based on {@link org.pankratzlab.common.Aliases#CHRS} that will
   * fail if multiple aliases are found.
   * 
   * @param colName Desired name of column
   * @param dieOnMissing
   * @return FileColumn<Byte>
   */
  public static IndexedFileColumn<Byte> chr(String colName, boolean dieOnMissing) {
    return new NumberWrapperColumn<>(new AliasedFileColumn(colName,
                                                           new Aliases(org.pankratzlab.common.Aliases.CHRS)),
                                     Positions::chromosomeNumber, dieOnMissing);
  }

  /**
   * Creates a non-case-sensitive {@link IntegerWrapperColumn} wrapped around an
   * {@link AliasedFileColumn} based on {@link org.pankratzlab.common.Aliases#POSITIONS} that will fail
   * if multiple aliases are found.
   * 
   * @param colName Desired name of column
   * @return IntegerWrapperColumn
   */
  public static IndexedFileColumn<Integer> pos(String colName) {
    return pos(colName, false);
  }

  /**
   * Creates a non-case-sensitive {@link IntegerWrapperColumn} wrapped around an
   * {@link AliasedFileColumn} based on {@link org.pankratzlab.common.Aliases#POSITIONS} that will fail
   * if multiple aliases are found.
   * 
   * @param colName Desired name of column
   * @param dieOnMissing
   * @return IntegerWrapperColumn
   */
  public static IndexedFileColumn<Integer> pos(String colName, boolean dieOnMissing) {
    return new IntegerWrapperColumn(new AliasedFileColumn(colName,
                                                          new Aliases(org.pankratzlab.common.Aliases.POSITIONS)),
                                    dieOnMissing);
  }

  /**
   * Creates a non-case-sensitive {@link DoubleWrapperColumn} wrapped around an
   * {@link AliasedFileColumn} based on {@link org.pankratzlab.common.Aliases#PVALUES} that will fail
   * if multiple aliases are found.
   * 
   * @param colName Desired name of column
   * @return DoubleWrapperColumn
   */
  public static IndexedFileColumn<Double> pVal(String colName) {
    return pVal(colName, false);
  }

  /**
   * Creates a non-case-sensitive {@link DoubleWrapperColumn} wrapped around an
   * {@link AliasedFileColumn} based on {@link org.pankratzlab.common.Aliases#PVALUES} that will fail
   * if multiple aliases are found.
   * 
   * @param colName Desired name of column
   * @param dieOnMissing
   * @return DoubleWrapperColumn
   */
  public static IndexedFileColumn<Double> pVal(String colName, boolean dieOnMissing) {
    return new DoubleWrapperColumn(new AliasedFileColumn(colName,
                                                         new Aliases(org.pankratzlab.common.Aliases.PVALUES)),
                                   dieOnMissing);
  }

  /**
   * Creates a non-case-sensitive {@link AliasedFileColumn} around
   * {@link org.pankratzlab.common.Aliases#ALLELES}{@code [0]} that will fail if multiple aliases are
   * found.
   * 
   * @param colName Desired name of column
   * @return AliasedFileColumn
   */
  public static AliasedFileColumn a1(String colName) {
    return new AliasedFileColumn(colName, new Aliases(org.pankratzlab.common.Aliases.ALLELES[0]));
  }

  /**
   * Creates a non-case-sensitive {@link AliasedFileColumn} around
   * {@link org.pankratzlab.common.Aliases#NS} that will fail if multiple aliases are found.
   * 
   * @param colName Desired name of column
   * @return AliasedFileColumn
   */
  public static IndexedFileColumn<Integer> n(String colName) {
    return n(colName, false);
  }

  /**
   * Creates a non-case-sensitive {@link AliasedFileColumn} around
   * {@link org.pankratzlab.common.Aliases#NS} that will fail if multiple aliases are found.
   * 
   * @param colName Desired name of column
   * @param dieOnMissing
   * @return AliasedFileColumn
   */
  public static IndexedFileColumn<Integer> n(String colName, boolean dieOnMissing) {
    return new IntegerWrapperColumn(new AliasedFileColumn(colName,
                                                          new Aliases(org.pankratzlab.common.Aliases.NS)),
                                    dieOnMissing);
  }

  /**
   * Creates a non-case-sensitive {@link AliasedFileColumn} around
   * {@link org.pankratzlab.common.Aliases#ALLELES}{@code [1]} that will fail if multiple aliases are
   * found.
   * 
   * @param colName Desired name of column
   * @return AliasedFileColumn
   */
  public static AliasedFileColumn a2(String colName) {
    return new AliasedFileColumn(colName, new Aliases(org.pankratzlab.common.Aliases.ALLELES[1]));
  }

  /**
   * Creates a non-case-sensitive {@link AliasedFileColumn} around
   * {@link org.pankratzlab.common.Aliases#ALLELE_FREQS} that will fail if multiple aliases are found.
   * 
   * @param colName Desired name of column
   * @return AliasedFileColumn
   */
  public static IndexedFileColumn<Double> alleleFreq(String colName) {
    return alleleFreq(colName, false);
  }

  /**
   * Creates a non-case-sensitive {@link AliasedFileColumn} around
   * {@link org.pankratzlab.common.Aliases#ALLELE_FREQS} that will fail if multiple aliases are found.
   * 
   * @param colName Desired name of column
   * @param dieOnMissing
   * @return AliasedFileColumn
   */
  public static IndexedFileColumn<Double> alleleFreq(String colName, boolean dieOnMissing) {
    return new DoubleWrapperColumn(new AliasedFileColumn(colName,
                                                         new Aliases(org.pankratzlab.common.Aliases.ALLELE_FREQS)),
                                   dieOnMissing);
  }

  /**
   * Creates a non-case-sensitive {@link AliasedFileColumn} around
   * {@link org.pankratzlab.common.Aliases#EFFECTS} that will fail if multiple aliases are found.
   * 
   * @param colName Desired name of column
   * @return FileColumn<Double>
   */
  public static IndexedFileColumn<Double> beta(String colName) {
    return beta(colName, false);
  }

  /**
   * Creates a non-case-sensitive {@link AliasedFileColumn} around
   * {@link org.pankratzlab.common.Aliases#EFFECTS} that will fail if multiple aliases are found.
   * 
   * @param colName Desired name of column
   * @param dieOnMissing
   * @return FileColumn<Double>
   */
  public static IndexedFileColumn<Double> beta(String colName, boolean dieOnMissing) {
    return new DoubleWrapperColumn(new AliasedFileColumn(colName,
                                                         new Aliases(org.pankratzlab.common.Aliases.EFFECTS)),
                                   dieOnMissing);
  }

  /**
   * Creates a non-case-sensitive {@link AliasedFileColumn} around
   * {@link org.pankratzlab.common.Aliases#STD_ERRS} that will fail if multiple aliases are found.
   * 
   * @param colName Desired name of column
   * @return FileColumn<Double>
   */
  public static IndexedFileColumn<Double> stdErr(String colName) {
    return stdErr(colName, false);
  }

  /**
   * Creates a non-case-sensitive {@link AliasedFileColumn} around
   * {@link org.pankratzlab.common.Aliases#STD_ERRS} that will fail if multiple aliases are found.
   * 
   * @param colName Desired name of column
   * @return FileColumn<Double>
   */
  public static IndexedFileColumn<Double> stdErr(String colName, boolean dieOnMissing) {
    return new DoubleWrapperColumn(new AliasedFileColumn(colName,
                                                         new Aliases(org.pankratzlab.common.Aliases.STD_ERRS)),
                                   dieOnMissing);
  }

}
