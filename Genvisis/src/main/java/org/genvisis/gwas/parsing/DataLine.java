package org.genvisis.gwas.parsing;

import java.util.HashMap;
import java.util.Map;

/**
 * Line of data as loaded by {@link FileParser}
 */
public class DataLine {

	/**
	 * Default value for a parseFail
	 */
	private Object defaultFailValue;
	private Map<FileColumn<?>, Object> lineValues;
	private Map<FileColumn<?>, Object> failValues;

	/**
	 * Constructor
	 * 
	 * @param defaultFailValue
	 */
	public DataLine(Object defaultFailValue) {
		lineValues = new HashMap<>();
		failValues = new HashMap<>();
		this.defaultFailValue = defaultFailValue;
	}

	/**
	 * Set value
	 * 
	 * @param fc FileColumn
	 * @param value value in line for given column
	 */
	public <T> void put(FileColumn<T> fc, T value) {
		lineValues.put(fc, value);
	}

	/**
	 * Mark the value for this column as a parse failure.
	 * 
	 * @param fc FileColumn
	 */
	public void fail(FileColumn<?> fc) {
		lineValues.put(fc, defaultFailValue);
	}

	/**
	 * Mark the value for this column as a parse failure, setting a specific value as the failure
	 * value.
	 * 
	 * @param fc
	 * @param failValue
	 */
	public void fail(FileColumn<?> fc, Object failValue) {
		lineValues.put(fc, failValue);
		failValues.put(fc, failValue);
	}

	/**
	 * Is there <i>any</i> value for the given FileColumn in this DataLine?<br />
	 * (May be a parse failure value, but at least a value is present.)
	 * 
	 * @param fc
	 * @return
	 */
	public boolean has(FileColumn<?> fc) {
		return lineValues.containsKey(fc);
	}

	/**
	 * Is there a value in this DataLine for the given FileColumn that isn't a parse failure value?
	 * 
	 * @param fc
	 * @return
	 */
	public boolean hasValid(FileColumn<?> fc) {
		boolean hasAtAll = lineValues.containsKey(fc);
		boolean defaultFail = lineValues.get(fc).equals(defaultFailValue);
		boolean nondefaultFail = failValues.containsKey(fc)
														 && failValues.get(fc).equals(lineValues.get(fc));
		return hasAtAll && !defaultFail && !nondefaultFail;
	}

	/**
	 * Get the String representation of the value in this line for the given FileColumn. Calls
	 * {@link Object#toString()}, so if no value is present a {@link NullPointerException} will be
	 * thrown.
	 * 
	 * @param fc
	 * @return
	 */
	public String getString(FileColumn<?> fc) {
		return lineValues.get(fc).toString();
	}

	/**
	 * Get the typed value from this DataLine for a given FileColumn. If the value was a parse failure
	 * value, a {@link ParseFailureException} will be thrown.
	 * 
	 * @param fc
	 * @return
	 * @throws ParseFailureException
	 */
	@SuppressWarnings("unchecked")
	public <T> T get(FileColumn<T> fc) throws ParseFailureException {
		if (defaultFailValue.equals(lineValues.get(fc))
				|| (failValues.containsKey(fc) && failValues.get(fc).equals(lineValues.get(fc)))) {
			throw new ParseFailureException("Failed to successfully parse value for column "
																			+ fc.getName());
		}
		return (T) lineValues.get(fc);
	}

	/**
	 * <b>WARNING: UNSAFE.</b><br />
	 * Should only be used inside an {@code if (hasValid(fc))} block, otherwise may cause
	 * {@link ClassCastException}s. <br />
	 * <br />
	 * Can also be used if the parse failure value for this FileColumn has the same type as the return
	 * value of the {@link FileColumn#getValue(String[])} method.
	 * 
	 * 
	 * @param fc
	 * @return
	 */
	@SuppressWarnings("unchecked")
	public <T> T getUnsafe(FileColumn<T> fc) {
		return (T) lineValues.get(fc);
	}

}
