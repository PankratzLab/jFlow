package org.genvisis.common;

import java.util.Vector;

import com.google.common.primitives.Doubles;

public class DoubleVector extends Vector<Double> implements PrimitiveVector {
	
	private static final long serialVersionUID = 7922754862926414680L;

	/**
	 * Create an empty vector
	 */
	public DoubleVector() {
		super();
	}

	/**
	 * Create an empty vector with the specified initial size.
	 */
	public DoubleVector(int initialSize) {
		this(new double[initialSize]);
	}

	/**
	 * Create a vector with the given initial values.
	 */
	public DoubleVector(double[] initialArray) {
		super(Doubles.asList(initialArray));
	}

	/**
	 * Add a given value to this vector if it is not already present.
	 * <p>
	 * FIXME: should use a set instead of creating this method.
	 * </p>
	 * <p>
	 * FIXME: duplicate logic with {@link IntVector} due to primitives being terrible.
	 * </p>
	 *
	 * @param v Value to add if not already present in this vector
	 * @return {@code true} if the value was added.
	 */
	public boolean addIfAbsent(double v) {
		boolean added = false;
		if (!Doubles.contains(Doubles.toArray(this), v)) {
			added = add(v);
		}
		return added;
	}
}
