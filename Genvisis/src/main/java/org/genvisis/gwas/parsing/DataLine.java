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
	 * 
	 * @param fc
	 * @param failValue
	 */
	public void fail(FileColumn<?> fc, Object failValue) {
		lineValues.put(fc, failValue);
		failValues.put(fc, failValue);
	}

	public boolean has(FileColumn<?> fc) {
		return lineValues.containsKey(fc);
	}

	public boolean hasValid(FileColumn<?> fc) {
		boolean hasAtAll = lineValues.containsKey(fc);
		boolean defaultFail = lineValues.get(fc).equals(defaultFailValue);
		boolean nondefaultFail = failValues.containsKey(fc)
														 && failValues.get(fc).equals(lineValues.get(fc));
		return hasAtAll && !defaultFail && !nondefaultFail;
	}

	public String getString(FileColumn<?> fc) {
		return lineValues.get(fc).toString();
	}

	@SuppressWarnings("unchecked")
	public <T> T get(FileColumn<T> fc) throws ParseFailureException {
		if (defaultFailValue.equals(lineValues.get(fc))) {
			throw new ParseFailureException("Failed to successfully parse value for column "
																			+ fc.getName());
		}
		return (T) lineValues.get(fc);
	}

	/**
	 * WARNING: UNSAFE.<br />
	 * Should only be used inside an {@code if (hasValid(fc))} block, otherwise may cause
	 * {@link ClassCastException}s.
	 * 
	 * @param fc
	 * @return
	 */
	@SuppressWarnings("unchecked")
	public <T> T getUnsafe(FileColumn<T> fc) {
		return (T) lineValues.get(fc);
	}

}
