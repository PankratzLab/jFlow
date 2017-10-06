package org.genvisis.bgen;


final class BGENBitMath {

	/**
	 * Turn an array of byte values (expected to be bit value flags) into an integer by bitwise-OR-ing
	 * the values and then shifting the result by one position.
	 * 
	 * @param littleEndian Order of importance
	 * @param bitVals Vararg of bytes, expected to be bits
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
	 * Converts a 4 byte array of unsigned bytes to an long
	 * 
	 * from http://www.petefreitag.com/item/183.cfm
	 * 
	 * @param b an array of 4 unsigned bytes
	 * @return a long representing the unsigned int
	 */
	public static final long unsignedIntToLong(byte[] b, boolean littleEndian) {
		long l = 0;
		l |= b[littleEndian ? 3 : 0] & 0xFF;
		l <<= 8;
		l |= b[littleEndian ? 2 : 1] & 0xFF;
		l <<= 8;
		l |= b[littleEndian ? 1 : 2] & 0xFF;
		l <<= 8;
		l |= b[littleEndian ? 0 : 3] & 0xFF;
		return l;
	}

	/**
	 * Converts a two byte array to an integer
	 * 
	 * @param b a byte array of length 2
	 * @return an int representing the unsigned short
	 */
	public static final int unsignedShortToInt(byte[] b, boolean littleEndian) {
		int i = 0;
		i |= b[littleEndian ? 1 : 0] & 0xFF;
		i <<= 8;
		i |= b[littleEndian ? 0 : 1] & 0xFF;
		return i;
	}

}
