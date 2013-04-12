// -Xms1024M -Xmx1024M     or even better: -Xmx15g
package cnv.filesys;

import java.io.*;
import common.Elision;

public class Compression {

	public static final byte[] REDUCED_PRECISION_XY_NAN_BYTES = new byte[] {(byte) 255, (byte) 255};
	public static final byte[] REDUCED_PRECISION_XY_OUT_OF_RANGE_BYTES = new byte[] {(byte) 255, (byte) 254};
	public static final float REDUCED_PRECISION_XY_OUT_OF_RANGE_FLOAT = (float) -1;
	public static final byte[] REDUCED_PRECISION_GCBAF_NAN_BYTES = new byte[] {(byte) 39, (byte) 18};
	public static final byte[] REDUCED_PRECISION_LRR_NAN_BYTES = new byte[] {(byte) 2, (byte) 0, (byte) 0};	//-13.1072
	public static final byte[] REDUCED_PRECISION_LRR_OUT_OF_RANGE_BYTES = new byte[] {(byte) 2, (byte) 0, (byte) 1}; //-13.1071
	public static final float REDUCED_PRECISION_LRR_OUT_OF_RANGE_LRR_FLOAT = (float) -13.1071; //{2, 0, 1}
//	public static final int BYTES_PER_SAMPLE_MARKER = 22;
	public static final int BYTES_PER_SAMPLE_MARKER = 12;
//	public static final int BYTES_PER_SAMPLE_MARKER = 10;


	public static byte[] objToBytes(Object obj) throws IOException{
		ByteArrayOutputStream byteArrayOutputStream = new ByteArrayOutputStream();
		ObjectOutputStream obejctOutputStream = new ObjectOutputStream(byteArrayOutputStream);
		obejctOutputStream.writeObject(obj);
		obejctOutputStream.flush();
		obejctOutputStream.close();
		byteArrayOutputStream.close();
		byte [] data = byteArrayOutputStream.toByteArray();
		return data;
	}

    public static Object bytesToObj(byte[] bytes) throws IOException, ClassNotFoundException {
    	Object obj;
        ByteArrayInputStream byteArrayInStream = new ByteArrayInputStream(bytes);
        ObjectInputStream objInStream = new ObjectInputStream(byteArrayInStream);
        obj = objInStream.readObject();
        objInStream.close();
        byteArrayInStream.close();
        return obj;
    }

    public static Object bytesToObj(byte[] bytes, int startLocation, int length) throws IOException, ClassNotFoundException {
    	Object obj;
        ByteArrayInputStream byteArrayInStream = new ByteArrayInputStream(bytes, startLocation, length);
        ObjectInputStream objInStream = new ObjectInputStream(byteArrayInStream);
        obj = objInStream.readObject();
        objInStream.close();
        byteArrayInStream.close();
        return obj;
    }





	public static byte[] floatToBytes(float data) {
//	    return new byte[] {
//		        (byte)((data >> 24) & 0xff),
//		        (byte)((data >> 16) & 0xff),
//		        (byte)((data >> 8) & 0xff),
//		        (byte)((data >> 0) & 0xff)
//		    };
		
		int intBits;
		intBits = Float.floatToRawIntBits(data);
		return intToBytes(intBits);
	}

	public static void floatToBytes(float value, byte[] array, int startPosition) {
		int data;
		
		data = Float.floatToRawIntBits(value);
		
		array[startPosition+0] = (byte)((data >> 24) & 0xff);
		array[startPosition+1] = (byte)((data >> 16) & 0xff);
		array[startPosition+2] = (byte)((data >> 8) & 0xff);
		array[startPosition+3] = (byte)((data >> 0) & 0xff);
	}
	
	public static float bytesToFloat(byte[] data) {
	    if (data == null || data.length != 4) {
	    	return 0x0;
	    } else {
	    	return Float.intBitsToFloat(bytesToInt(data));
	    }
	}

	public static byte[] intToBytes(int data) {
	    return new byte[] {
		        (byte)((data >> 24) & 0xff),
		        (byte)((data >> 16) & 0xff),
		        (byte)((data >> 8) & 0xff),
		        (byte)(data & 0xff)
		    };
	}

	public static int bytesToInt(byte[] data) {
	    if (data == null || data.length != 4) {
	    	return 0x0;
	    } else {
	    	return bytesToInt(data, 0);
	    }
	}

	public static int bytesToInt(byte[] data, int startPos) {
	    if (data == null || data.length < 4 || startPos > (data.length-4)) {
	    	return 0x0;
	    } else {
//	    	NOTE: type cast not necessary for int
		    return (int)(
		            (0xff & data[startPos]) << 24  |
		            (0xff & data[startPos + 1]) << 16  |
		            (0xff & data[startPos + 2]) << 8   |
		            (0xff & data[startPos + 3])
		            );
	    }
	}

    public static byte[] longToBytes(long v) {
    	return new byte[] {
    			(byte) ((v >>> 56) & 0xFF),
    			(byte) ((v >>> 48) & 0xFF),
    			(byte) ((v >>> 40) & 0xFF),
    			(byte) ((v >>> 32) & 0xFF),
    			(byte) ((v >>> 24) & 0xFF),
    			(byte) ((v >>> 16) & 0xFF),
    			(byte) ((v >>>  8) & 0xFF),
    			(byte) ((v >>>  0) & 0xFF)
    	};
    }


	public static long bytesToLong(byte[] data) {
	    if (data == null || data.length != 8) {
	    	return 0x0;
	    } else {
	    	return bytesToLong(data, 0);
	    }
	}

	public static long bytesToLong(byte[] data, int startPos) {
	    if (data==null || data.length<8 || startPos>(data.length-8)) {
	    	return 0x0;
	    } else {

//		    Note: This block is the wrong code
//		    return (long)(
//		            (0xff & data[startPos]) << 56  |
//		            (0xff & data[startPos+1]) << 48  |
//		            (0xff & data[startPos+2]) << 40	|
//		            (0xff & data[startPos+3]) << 32	|
//		            (0xff & data[startPos+4]) << 24  |
//		            (0xff & data[startPos+5]) << 16  |
//		            (0xff & data[startPos+6]) << 8   |
//		            (0xff & data[startPos+7])
//		            );

//	    	Note: This block is the working code
//		    return (
//		    		(long)(0xff & data[startPos]) << 56 |
//		    		(long)(0xff & data[startPos+1]) << 48 |
//		    		(long)(0xff & data[startPos+2]) << 40 |
//		    		(long)(0xff & data[startPos+3]) << 32 |
//		    		(long)(0xff & data[startPos+4]) << 24 |
//		    		(long)(0xff & data[startPos+5]) << 16 |
//		    		(long)(0xff & data[startPos+6]) << 8 |
//		    		(long)(0xff & data[startPos+7])
//		    		);

	    	return ((long)(bytesToInt(data, startPos)) << 32) + (bytesToInt(data, startPos + 4) & 0xFFFFFFFFL);
	    }
	}





	/**
	 * Converts the float xyValue into byte[]. This is a simple version that does not process the out of range values.
	 * @param xyValue the float to be converted.
	 * @return
	 */
	public static byte[] xyCompress(float xyValue) throws Elision {
		int data;
		if (Float.isNaN(xyValue)) {
			return REDUCED_PRECISION_XY_NAN_BYTES;
		} else if (xyValue>65.533 || xyValue<0) {
			throw new Elision("The value of X or Y (" + xyValue + ") is out of the specified range of [0, 32.767] for compression.");
//			return REDUCED_PRECISION_XY_OUT_OF_RANGE_BYTES;
		} else {
//			Note: Currently, the conversion from float to int is through the operation of *1000. Need to change it to bit conversion to pursue higher accuracy.
			data = (int) Math.round(xyValue * 1000);
			return new byte[] {(byte)((data >> 8) & 0xff), (byte)(data & 0xff)};
		}
	}

	/**
	 * Converts the float xyValue into byte[], which later becomes part of byte[] array.
	 * Out of range values are converted into REDUCED_PRECISION_XY_OUT_OF_RANGE_BYTES, and the program will return false.
	 * Within range values are converted regularly, and the program will return true.
	 * @param xOrY the float to be converted.
	 * @param array the byte[] to hold the output.
	 * @param startPosition the position of byte[] array to hold the output.
	 * @return false means the xOrY is out of the range for compression, while true means the compression is successful. 
	 */
	public static boolean xyCompress(float xyValue, byte[] array, int startPosition) {
		int data;
		if (Float.isNaN(xyValue)) {
			array[startPosition] = REDUCED_PRECISION_XY_NAN_BYTES[0];
			array[startPosition+1] = REDUCED_PRECISION_XY_NAN_BYTES[1];
			return true;
		} else if (xyValue>65.533 || xyValue<0) {
			array[startPosition] = REDUCED_PRECISION_XY_OUT_OF_RANGE_BYTES[0];
			array[startPosition+1] = REDUCED_PRECISION_XY_OUT_OF_RANGE_BYTES[1];
			return false;
		} else {
//			Note: Currently, the conversion from float to int is through the operation of *1000. Need to change it to bit conversion to pursue higher accuracy.
			data = (int) Math.round(xyValue * 1000);
			array[startPosition] = (byte)((data >> 8) & 0xff);
			array[startPosition+1] = (byte)(data & 0xff);
			return true;
		}
	}

	/**
	 * Converts a byte[] into a float X or Y value.
	 * @param data the byte[] to hold the output.
	 * @return the float. 
	 */
	public static float xyDecompress(byte[] data) {
	    if (data == null || data.length != 2) {
	    	return (Float) null;
	    } else if ( data[0]==REDUCED_PRECISION_XY_NAN_BYTES[0] && data[1]==REDUCED_PRECISION_XY_NAN_BYTES[1] ) {
	    	return Float.NaN;
	    } else if ( data[0]==REDUCED_PRECISION_XY_OUT_OF_RANGE_BYTES[0] && data[1]==REDUCED_PRECISION_XY_OUT_OF_RANGE_BYTES[1] ) {
	    	return -1;
	    } else {
//			Note: Currently, the conversion from float to int is through the operation of /1000. Need to change it to bit conversion to pursue higher accuracy.
	    	return (int)(((byte) 0 << 16 & 0xff0000) | (data[0] << 8 & 0xff00) | (0xff & data[1])) / (float)1000;
	    }
	}





	/**
	 * Converts the float GC or BAF value into byte[].
	 * @param gcOrBaf the float to be converted.
	 * @return
	 * @throws Elision 
	 */
	public static byte[] gcBafCompress(float gcOrBaf) throws Elision {
		if (Float.isNaN(gcOrBaf)) {
			return REDUCED_PRECISION_GCBAF_NAN_BYTES;
		} else if (gcOrBaf>1.000 || gcOrBaf<0) {
			throw new Elision("The value of GC or BAF (" + gcOrBaf + ") should be in the range of 0.0000 through 1.0000.");
		} else {
			short data;
//			Note: Currently, the conversion from float to int is through the operation of *10000. Need to change it to bit conversion to pursue higher accuracy.
			data = (short) Math.round(gcOrBaf * 10000);
		    return new byte[] { (byte)( ((byte)0 & 0xc0) | ((data >> 8) & 0x3f) ), (byte)(data & 0xff)};
		}
	}

	/**
	 * Converts the float GC or BAF value into byte[], which later becomes part of byte[] array.
	 * @param gcOrBaf the float to be converted.
	 * @param array the byte[] to hold the output.
	 * @param startPosition the position of byte[] array to hold the output.
	 * @return
	 * @throws Elision 
	 */
	public static void gcBafCompress(float gcOrBaf, byte[] array, int startPosition) throws Elision {
		if (Float.isNaN(gcOrBaf)) {
			array[startPosition] = REDUCED_PRECISION_GCBAF_NAN_BYTES[0];
			array[startPosition+1] = REDUCED_PRECISION_GCBAF_NAN_BYTES[1];
		} else if (gcOrBaf>1.000 || gcOrBaf<0) {
			throw new Elision("The value of GC or BAF (" + gcOrBaf + ") should be in the range of 0.0000 through 1.0000.");
		} else {
			short data;
//			Note: Currently, the conversion from float to int is through the operation of *10000. Need to change it to bit conversion to pursue higher accuracy.
			data = (short) Math.round(gcOrBaf * 10000);
			array[startPosition] = (byte)( ((byte)0 & 0xc0) | ((data >> 8) & 0x3f) );
			array[startPosition+1] = (byte)(data & 0xff);
		}
	}
	
	/**
	 * Converts a byte[] into a float GC or BAF value.
	 * @param data the byte[] to hold the output.
	 * @return the float. 
	 */
	public static float gcBafDecompress(byte[] data) {
	    if (data == null || data.length != 2) {
	    	return (Float) null;
	    } else if (data[0]==REDUCED_PRECISION_GCBAF_NAN_BYTES[0] && data[1]==REDUCED_PRECISION_GCBAF_NAN_BYTES[1]) {
	    	return Float.NaN;
	    } else {
//			Note: Currently, the conversion from float to int is through the operation of /10000. Need to change it to bit conversion to pursue higher accuracy.
	    	return (float)(data[0] << 8 | (0xff & data[1])) / (float)10000;
	    }
	}




	/**
	 * Converts the float LRR value into byte[].
	 * @param lrr the float to be converted.
	 * @return
	 * @throws Elision 
	 */
	public static byte[] lrrCompress(float lrr) throws Elision {
		int data;
		if (Float.isNaN(lrr)) {
			return REDUCED_PRECISION_LRR_NAN_BYTES;
		} else if (lrr<(float)-13.1071 || lrr>(float)13.1071) {
//			return REDUCED_PRECISION_LRR_OUT_OF_RANGE_BYTES;
			throw new Elision("The value of LRR (" + lrr + ") is out of the range of -13.1072 through 13.1069.");
		} else {
//			Note: Currently, the conversion from float to int is through the operation of *10000. Need to change it to bit conversion to pursue higher accuracy.
			data = (int) Math.round(lrr * 10000);
//			return new byte[] {(byte)((data >> 24 & 0x80) | ((data >> 16) & 0x7f)), (byte)((data >> 8) & 0xff), (byte)(data & 0xff)};
			return new byte[] {(byte)((data >> 30 & 0x02) | ((data >> 16) & 0x1)), (byte)((data >> 8) & 0xff), (byte)(data & 0xff)};
		}
	}

	/**
	 * Converts the float LRR value into byte[], which later becomes part of byte[] array.
	 * @param lrr the float to be converted.
	 * @param array the byte[] to hold the output.
	 * @param startPosition the position of byte[] array to hold the output.
	 * @return
	 * @throws Elision 
	 */
	public static byte lrrCompress(float lrr, byte[] array, int startPosition) {
		int data;
		if (Float.isNaN(lrr)) {
			array[startPosition] = REDUCED_PRECISION_LRR_NAN_BYTES[0];
			array[startPosition+1] = REDUCED_PRECISION_LRR_NAN_BYTES[1];
			array[startPosition+2] = REDUCED_PRECISION_LRR_NAN_BYTES[2];
			return 0;
		} else if (lrr<(float)-13.1071 || lrr>(float)13.1071) {
			array[startPosition] = REDUCED_PRECISION_LRR_OUT_OF_RANGE_BYTES[0];
			array[startPosition+1] = REDUCED_PRECISION_LRR_OUT_OF_RANGE_BYTES[1];
			array[startPosition+2] = REDUCED_PRECISION_LRR_OUT_OF_RANGE_BYTES[2];
			return -1;
		} else {
//			Note: Currently, the conversion from float to int is through the operation of *10000. Need to change it to bit conversion to pursue higher accuracy.
			data = (int) Math.round(lrr * 10000);
			array[startPosition] = (byte)((data >> 30 & 0x2) | ((data >> 16) & 0x1));
			array[startPosition+1] = (byte)((data >> 8) & 0xff);
			array[startPosition+2] = (byte)(data & 0xff);
			return 0;
		}
	}
	
	/**
	 * Converts a byte[] into a float LRR value.
	 * If the byte[] represents an out of range value, the program returns REDUCED_PRECISION_LRR_OUT_OF_RANGE_LRR_FLOAT.
	 * @param array the byte[] to hold the output.
	 * @return the float. 
	 */
	public static float lrrDecompress(byte[] data) {
	    if (data == null || data.length != 3) {
	    	return (Float) null;
	    } else if (data[0]==REDUCED_PRECISION_LRR_NAN_BYTES[0] && data[1]==REDUCED_PRECISION_LRR_NAN_BYTES[1] && data[2]==REDUCED_PRECISION_LRR_NAN_BYTES[2]) {
	    	return Float.NaN;
	    } else if (data[0]==REDUCED_PRECISION_LRR_OUT_OF_RANGE_BYTES[0] && data[1]==REDUCED_PRECISION_LRR_OUT_OF_RANGE_BYTES[1] && data[2]==REDUCED_PRECISION_LRR_OUT_OF_RANGE_BYTES[2]) {
	    	return (float) REDUCED_PRECISION_LRR_OUT_OF_RANGE_LRR_FLOAT;
	    } else if (data[0]<2) {
	    	return (int)( 0x0001ffff & (((0x1 & data[0]) << 16) | ((0xff & data[1]) << 8) | (0xff & data[2]))) / (float)10000;
	    } else {
	    	return (int)((0xfffe0000) | ((0x1 & data[0]) << 16) | ((0xff & data[1]) << 8) | (0xff & data[2])) / (float)10000;
	    }
	}








	/**
	 * Converts the AB genotype and forward genotype value into byte.
	 * @param abGenotype
	 * @param forwardGenotype
	 * @return
	 * @throws Elision 
	 */
	public static byte genotypeCompress(byte abGenotype, byte forwardGenotype) {
//	    return (byte) ((byte)((forwardGenotype << 3) & 0xfc) | (abGenotype & 0x03));
		if (abGenotype>=0) {
			return (byte) (((forwardGenotype << 3) & 0xf8) | (abGenotype & 0x03));
		} else {
			return (byte) (((forwardGenotype << 3) & 0xf8) | ((abGenotype & 0x80) >> 5) | ((~abGenotype +1) & 0x03));
		}
	}

	/**
	 * Converts the AB genotype and forward genotype value into byte.
	 * @param abGenotype
	 * @param forwardGenotype
	 * @param array the byte[] to hold the output.
	 * @param startPosition the position of byte[] array to hold the output.
	 * @return
	 * @throws Elision 
	 */
	public static void genotypeCompress(byte abGenotype, byte forwardGenotype, byte[] array, int startPosition) {
//		array[startPosition] = (byte) ((byte)((forwardGenotype << 3) & 0xfc) | (abGenotype & 0x03));
		if (abGenotype>=0) {
			array[startPosition] = (byte) (((forwardGenotype << 3) & 0xf8) | (abGenotype & 0x03));
		} else {
			array[startPosition] = (byte) (((forwardGenotype << 3) & 0xf8) | ((abGenotype & 0x80) >> 5) | ((~abGenotype +1) & 0x03));
		}
	}
	
	/**
	 * Converts a byte into a set of AB genotype and forward genotype values.
	 * @param data the byte for decompression.
	 * @return the float. 
	 */
	public static byte[] genotypeDecompress(byte data) {
		if ((data & 0x04) >= 1) {
			return new byte[] {(byte) (~(data & 0x03) +1), (byte) ((data & 0xf8) >> 3)};
		} else {
			return new byte[] {(byte) (data & 0x03), (byte) ((data & 0xf8) >> 3)};
		}
//    	return new byte[] {(byte) (data & 0x03), (byte) ((data & 0xfc) >> 3)};
//	    if (data == 0) {
//	    	return new byte[] {(byte) -1, (byte) -1};
//	    } else {
//	    	return new byte[] {(byte) (data & 0x03), (byte) ((data & 0xfc) >> 2)};
//	    }
	}






	public static void showBits(byte a) {
		System.out.println((          (a & 0x80)>>7)
							+ "" + ((a & 0x40)>>6)
							+ "" + ((a & 0x20)>>5)
							+ "" + ((a & 0x10)>>4)
							+ "" + ((a & 0x08)>>3)
							+ " " + ((a & 0x04)>>2)
							+ "" + ((a & 0x02)>>1)
							+ "" + ((a & 0x01)));
	}

	public static void showBits(int a) {
//		System.out.println(Integer.toBinaryString(a));
		System.out.println(       ((a & 0x80000000)>>31)
						   + " " + ((a & 0x40000000)>>30)
						   + "" + ((a & 0x20000000)>>29)
						   + "" + ((a & 0x10000000)>>28)
						   + "" + ((a & 0x08000000)>>27)
						   + "" + ((a & 0x04000000)>>26)
						   + "" + ((a & 0x02000000)>>25)
						   + "" + ((a & 0x01000000)>>24)
						   + "" + ((a & 0x00800000)>>23)
						   + " " + ((a & 0x00400000)>>22)
						   + "" + ((a & 0x00200000)>>21)
						   + "" + ((a & 0x00100000)>>20)
						   + "" + ((a & 0x00080000)>>19)
						   + "" + ((a & 0x00040000)>>18)
						   + "" + ((a & 0x00020000)>>17)
						   + "" + ((a & 0x00010000)>>16)
						   + "" + ((a & 0x00008000)>>15)
						   + "" + ((a & 0x00004000)>>14)
						   + "" + ((a & 0x00002000)>>13)
						   + "" + ((a & 0x00001000)>>12)
						   + "" + ((a & 0x00000800)>>11)
						   + "" + ((a & 0x00000400)>>10)
						   + "" + ((a & 0x00000200)>>9)
						   + "" + ((a & 0x00000100)>>8)
						   + "" + ((a & 0x00000080)>>7)
						   + "" + ((a & 0x00000040)>>6)
						   + "" + ((a & 0x00000020)>>5)
						   + "" + ((a & 0x00000010)>>4)
						   + "" + ((a & 0x00000008)>>3)
						   + "" + ((a & 0x00000004)>>2)
						   + "" + ((a & 0x00000002)>>1)
						   + "" + ((a & 0x00000001)));
	}
}
