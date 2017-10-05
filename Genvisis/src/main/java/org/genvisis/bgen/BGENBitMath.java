package org.genvisis.bgen;


final class BGENBitMath {

	public static final int bitsToInt(boolean littleEndian, byte... bitVals) {
		int n = 0;
		if (littleEndian) {
			for (int i = 0; i < bitVals.length; i++) {
				n |= bitVals[i];
				if (i < bitVals.length - 1) {
					n <<= 1;
				}
			}
		} else {
			for (int i = bitVals.length - 1; i >= 0; i--) {
				n |= bitVals[i];
				if (i > 0) {
					n <<= 1;
				}
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
	public static final byte getBit(long v, int p) {
		return (byte) ((v >> p) & 1);
	}


	/**
	 * Get the Nth bit from a byte value
	 * 
	 * @param v byte value
	 * @param p position of desired bit
	 * @return bit value at position <code>p</code> in value <code>v</code>.
	 */
	public static final byte getBit(byte v, int p) {
		return (byte) ((v >> p) & 1);
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
