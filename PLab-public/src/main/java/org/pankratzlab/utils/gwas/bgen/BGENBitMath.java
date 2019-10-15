package org.pankratzlab.utils.gwas.bgen;

public final class BGENBitMath {

  /**
   * Turn an array of byte values into a float without using ByteBuffers
   * 
   * @param littleEndian Endian order of the given byte array
   * @param bytes array of four (4) bytes
   * @return float
   * @throws IllegalArgumentException if the given array is not of length 4
   */
  public static float bytesToFloat(boolean littleEndian, byte... bytes) {
    if (bytes.length != 4) {
      throw new IllegalArgumentException("Given array must contain 4 bytes total.");
    }
    return Float.intBitsToFloat(bytes[littleEndian ? 3 : 0] << 24
                                | (bytes[littleEndian ? 2 : 1] & 0xFF) << 16
                                | (bytes[littleEndian ? 1 : 2] & 0xFF) << 8
                                | (bytes[littleEndian ? 0 : 3] & 0xFF));
  }

  public static byte[] floatToBytes(boolean littleEndian, float value) {
    byte[] byts = new byte[4];
    int bits = Float.floatToIntBits(value);
    byts[littleEndian ? 0 : 3] = (byte) (bits & 0xFF);
    byts[littleEndian ? 1 : 2] = (byte) ((bits >> 8) & 0xFF);
    byts[littleEndian ? 2 : 1] = (byte) ((bits >> 16) & 0xFF);
    byts[littleEndian ? 3 : 0] = (byte) ((bits >> 24) & 0xFF);
    return byts;
  }

  /**
   * Turn an array of byte values into a double without using ByteBuffers
   * 
   * @param littleEndian Endian order of the given byte array
   * @param bytes array of eight (8) bytes
   * @return double
   * @throws IllegalArgumentException if the given array is not of length 8
   */
  public static double bytesToDouble(boolean littleEndian, byte... bytes) {
    if (bytes.length != 8) {
      throw new IllegalArgumentException("Given array must contain 8 bytes total.");
    }
    return Double.longBitsToDouble((long) bytes[littleEndian ? 7 : 0] << 56
                                   | (long) (bytes[littleEndian ? 6 : 1] & 0xFF) << 48
                                   | (long) (bytes[littleEndian ? 5 : 2] & 0xFF) << 40
                                   | (long) (bytes[littleEndian ? 4 : 3] & 0xFF) << 32
                                   | (long) (bytes[littleEndian ? 3 : 4] & 0xFF) << 24
                                   | (long) (bytes[littleEndian ? 2 : 5] & 0xFF) << 16
                                   | (long) (bytes[littleEndian ? 1 : 6] & 0xFF) << 8
                                   | (long) (bytes[littleEndian ? 0 : 7] & 0xFF));
  }

  public static byte[] doubleToBytes(boolean littleEndian, double value) {
    byte[] byts = new byte[8];
    long bits = Double.doubleToLongBits(value);
    byts[littleEndian ? 0 : 7] = (byte) (bits & 0xFF);
    byts[littleEndian ? 1 : 6] = (byte) ((bits >> 8) & 0xFF);
    byts[littleEndian ? 2 : 5] = (byte) ((bits >> 16) & 0xFF);
    byts[littleEndian ? 3 : 4] = (byte) ((bits >> 24) & 0xFF);
    byts[littleEndian ? 4 : 3] = (byte) ((bits >> 32) & 0xFF);
    byts[littleEndian ? 5 : 2] = (byte) ((bits >> 40) & 0xFF);
    byts[littleEndian ? 6 : 1] = (byte) ((bits >> 48) & 0xFF);
    byts[littleEndian ? 7 : 0] = (byte) ((bits >> 56) & 0xFF);
    return byts;
  }

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
   * Converts a byte array to an long Adapted from http://www.petefreitag.com/item/183.cfm
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
