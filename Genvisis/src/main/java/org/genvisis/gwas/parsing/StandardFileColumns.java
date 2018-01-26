package org.genvisis.gwas.parsing;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.genvisis.common.ArrayUtils;

/**
 * Commonly used FileColumns. If you end up creating a FileColumn more than once, consider putting
 * it here!
 */
public final class StandardFileColumns {
	private StandardFileColumns() {}

	/**
	 * A special {@link FileColumn&lt;String&gt;} that returns all values in a line except the value
	 * in the key column.
	 * 
	 * @param key {@link AliasedFileColumn}
	 * @param outDelim Delimiter with which to join the values
	 * @return {@link FileColumn&lt;String&gt;}
	 */
	public static final FileColumn<String> allButKey(AliasedFileColumn key, String outDelim) {
		return new FileColumn<String>() {

			String outHeader;
			List<Integer> outIndices;

			@Override
			public void initialize(Map<String, Integer> headerMap) {
				key.initialize(headerMap);
				String[] heads = new String[headerMap.size() - 1];
				outIndices = new ArrayList<>();
				int keyInd = key.getMatchedIndex();
				for (Entry<String, Integer> ent : headerMap.entrySet()) {
					if (ent.getValue() == keyInd) {
						continue;
					}
					int v = ent.getValue();
					if (v > keyInd) {
						v--;
					}
					heads[v] = ent.getKey();
					outIndices.add(ent.getValue());
				}
				outHeader = ArrayUtils.toStr(heads, outDelim);
				Collections.sort(outIndices);
			}

			@Override
			public String getValue(String[] line) throws ParseFailureException {
				StringBuilder sb = new StringBuilder();
				for (int i = 0, count = outIndices.size(); i < count; i++) {
					sb.append(line[outIndices.get(i)]);
					if (i < count - 1) {
						sb.append(outDelim);
					}
				}
				return sb.toString();
			}

			@Override
			public String getName() {
				return outHeader;
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
	public static final FileColumn<Double> log10(final String name,
																							 final FileColumn<Double> val) {
		return new CachedFileColumn<Double>(name) {
			@Override
			public void initialize(Map<String, Integer> headerMap) {
				val.initialize(headerMap);
			}

			@Override
			public Double calculateValue(String[] line) throws ParseFailureException {
				return Math.log10(val.getValue(line));
			}

			@Override
			public int hashCode() {
				final int prime = 31;
				int result = 1;
				result = prime * result + ((getName() == null) ? 0 : getName().hashCode());
				result = prime * result + val.hashCode();
				return result;
			}

			@Override
			public boolean equals(Object obj) {
				if (this == obj)
					return true;
				if (obj == null)
					return false;
				if (getClass() != obj.getClass())
					return false;
				FileColumn<?> other;
				if (obj instanceof FileColumn<?>) {
					other = (FileColumn<?>) obj;
				} else {
					return false;
				}
				if (!val.equals(other))
					return false;
				if (getName() == null) {
					if (other.getName() != null)
						return false;
				} else if (!getName().equals(other.getName()))
					return false;
				return true;
			}
		};
	}

	/**
	 * Creates a non-case-sensitive {@link AliasedFileColumn} around
	 * {@link org.genvisis.common.Aliases#MARKER_NAMES} that will fail if multiple aliases are found.
	 * 
	 * @param colName Desired name of column
	 * @return AliasedFileColumn
	 */
	public static AliasedFileColumn snp(String colName) {
		return new AliasedFileColumn(colName, new Aliases(org.genvisis.common.Aliases.MARKER_NAMES));
	}

	/**
	 * Creates a non-case-sensitive {@link IntegerWrapperColumn} wrapped around an
	 * {@link AliasedFileColumn} based on {@link org.genvisis.common.Aliases#CHRS} that will fail if
	 * multiple aliases are found.
	 * 
	 * @param colName Desired name of column
	 * @return IntegerWrapperColumn
	 */
	public static FileColumn<Integer> chr(String colName) {
		return new IntegerWrapperColumn(new AliasedFileColumn(colName,
																													new Aliases(org.genvisis.common.Aliases.CHRS)));
	}

	/**
	 * Creates a non-case-sensitive {@link IntegerWrapperColumn} wrapped around an
	 * {@link AliasedFileColumn} based on {@link org.genvisis.common.Aliases#POSITIONS} that will fail
	 * if multiple aliases are found.
	 * 
	 * @param colName Desired name of column
	 * @return IntegerWrapperColumn
	 */
	public static FileColumn<Integer> pos(String colName) {
		return new IntegerWrapperColumn(new AliasedFileColumn(colName,
																													new Aliases(org.genvisis.common.Aliases.POSITIONS)));
	}

	/**
	 * Creates a non-case-sensitive {@link DoubleWrapperColumn} wrapped around an
	 * {@link AliasedFileColumn} based on {@link org.genvisis.common.Aliases#PVALUES} that will fail
	 * if multiple aliases are found.
	 * 
	 * @param colName Desired name of column
	 * @return DoubleWrapperColumn
	 */
	public static FileColumn<Double> pVal(String colName) {
		return new DoubleWrapperColumn(new AliasedFileColumn(colName,
																												 new Aliases(org.genvisis.common.Aliases.PVALUES)));
	}

	/**
	 * Creates a non-case-sensitive {@link AliasedFileColumn} around
	 * {@link org.genvisis.common.Aliases#ALLELES}{@code [0]} that will fail if multiple aliases are
	 * found.
	 * 
	 * @param colName Desired name of column
	 * @return AliasedFileColumn
	 */
	public static AliasedFileColumn a1(String colName) {
		return new AliasedFileColumn(colName, new Aliases(org.genvisis.common.Aliases.ALLELES[0]));
	}

	/**
	 * Creates a non-case-sensitive {@link AliasedFileColumn} around
	 * {@link org.genvisis.common.Aliases#NS} that will fail if multiple aliases are found.
	 * 
	 * @param colName Desired name of column
	 * @return AliasedFileColumn
	 */
	public static FileColumn<Integer> n(String colName) {
		return new IntegerWrapperColumn(new AliasedFileColumn(colName,
																													new Aliases(org.genvisis.common.Aliases.NS)));
	}

	/**
	 * Creates a non-case-sensitive {@link AliasedFileColumn} around
	 * {@link org.genvisis.common.Aliases#ALLELES}{@code [1]} that will fail if multiple aliases are
	 * found.
	 * 
	 * @param colName Desired name of column
	 * @return AliasedFileColumn
	 */
	public static AliasedFileColumn a2(String colName) {
		return new AliasedFileColumn(colName, new Aliases(org.genvisis.common.Aliases.ALLELES[1]));
	}

	/**
	 * Creates a non-case-sensitive {@link AliasedFileColumn} around
	 * {@link org.genvisis.common.Aliases#ALLELE_FREQS} that will fail if multiple aliases are found.
	 * 
	 * @param colName Desired name of column
	 * @return AliasedFileColumn
	 */
	public static FileColumn<Double> alleleFreq(String colName) {
		return new DoubleWrapperColumn(new AliasedFileColumn(colName,
																												 new Aliases(org.genvisis.common.Aliases.ALLELE_FREQS)));
	}

	/**
	 * Creates a non-case-sensitive {@link AliasedFileColumn} around
	 * {@link org.genvisis.common.Aliases#EFFECTS} that will fail if multiple aliases are found.
	 * 
	 * @param colName Desired name of column
	 * @return FileColumn<Double>
	 */
	public static FileColumn<Double> beta(String colName) {
		return new DoubleWrapperColumn(new AliasedFileColumn(colName,
																												 new Aliases(org.genvisis.common.Aliases.EFFECTS)));
	}

	/**
	 * Creates a non-case-sensitive {@link AliasedFileColumn} around
	 * {@link org.genvisis.common.Aliases#STD_ERRS} that will fail if multiple aliases are found.
	 * 
	 * @param colName Desired name of column
	 * @return FileColumn<Double>
	 */
	public static FileColumn<Double> stdErr(String colName) {
		return new DoubleWrapperColumn(new AliasedFileColumn(colName,
																												 new Aliases(org.genvisis.common.Aliases.STD_ERRS)));
	}

}
