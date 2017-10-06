package org.genvisis.bgen;


final class BGENBitMath {

	/**
	 * Turn an array of byte values (expected to be bit value flags) into an integer by bitwise-OR-ing
	 * the values and then shifting the result by one position.
	 * 
	 * @param littleEndian Order of importance
	 * @param bitVals boolean vararg of bit flags
	 * @return
	 */
	public static final int bitsToInt(boolean littleEndian, boolean... bitVals) {
		int n = 0;
		for (int i = 0; i < bitVals.length; i++) {
			n |= bitVals[littleEndian ? i : (bitVals.length - 1 - i)] ? 1 : 0;
			if (i < bitVals.length - 1) {
				n <<= 1;
			}
		}
		return n;
	}

	/**
	 * Get the Nth bit from a long value
	 * 
	 * @param v long value
	 * @param p position of desired bit
	 * @return bit value at position <code>p</code> in value <code>v</code>.
	 */
	public static final boolean getBit(long v, int p) {
		return getBit(v, p, 64);
	}

	/**
	 * Get the Nth bit from a byte value
	 * 
	 * @param v byte value
	 * @param p position of desired bit
	 * @return bit value at position <code>p</code> in value <code>v</code>.
	 */
	public static final boolean getBit(byte v, int p) {
		return getBit((long) v, p, 8);
	}

	private static final boolean getBit(long v, int p, int maxBitsInType) {
		if (p < 0 || p > maxBitsInType) {
			throw new IllegalArgumentException("Requested index (" + p
																				 + ") is outside the number of bits in this type [1-"
																				 + maxBitsInType + "]");
		}
		return ((v >> p) & 1) == 1 ? true : false;
	}

	/**
	 * Converts a byte array to an long
	 * 
	 * from http://www.petefreitag.com/item/183.cfm
	 * 
	 * @param b an array of bytes
	 * @return a long representing the value of the byte array
	 */
	public static final long bytesToLong(byte[] b, boolean littleEndian) {
		long l = 0;
		for (int i = 0; i < b.length; i++) {
			l |= b[littleEndian ? (b.length - i - 1) : i] & 0xFF;
			if (i < b.length - 1) {
				l <<= 8;
			}
		}
		return l;
	}

}
