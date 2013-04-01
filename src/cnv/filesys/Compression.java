// -Xms1024M -Xmx1024M     or even better: -Xmx15g
package cnv.filesys;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.MappedByteBuffer;
import java.nio.channels.FileChannel;
import java.util.*;

import common.Files;


public class Compression {

	public static final byte[] REDUCED_PRECISION_XY_NAN_BYTES = new byte[] {(byte) 255, (byte) 255};
	public static final byte[] REDUCED_PRECISION_XY_OUT_OF_RANGE_BYTES = new byte[] {(byte) 255, (byte) 254};
	public static final float REDUCED_PRECISION_XY_OUT_OF_RANGE_FLOAT = (float) -1;
	public static final byte[] REDUCED_PRECISION_GCBAF_NAN_BYTES = new byte[] {(byte) 39, (byte) 18};
	public static final byte[] REDUCED_PRECISION_LRR_NAN_BYTES = new byte[] {(byte) 2, (byte) 0, (byte) 0};	//-13.1072
	public static final byte[] REDUCED_PRECISION_LRR_OUT_OF_RANGE_BYTES = new byte[] {(byte) 2, (byte) 0, (byte) 1}; //-13.1071
	public static final float REDUCED_PRECISION_LRR_OUT_OF_RANGE_LRR_FLOAT = (float) -13.1071; //{2, 0, 1}
	public static final int BYTES_PER_SAMPLE_MARKER = 22;
	public static final int BYTES_PER_SAMPLE_MARKER_2 = 12;
	public static final int BYTES_PER_SAMPLE_MARKER_10 = 10;


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
		int intBits;
		
		intBits = Float.floatToRawIntBits(data);
				
//	    return new byte[] {
//		        (byte)((data >> 24) & 0xff),
//		        (byte)((data >> 16) & 0xff),
//		        (byte)((data >> 8) & 0xff),
//		        (byte)((data >> 0) & 0xff)
//		    };
		
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
	    if (data == null || data.length != 4) return 0x0;
	    // ----------
	    return bytesToInt(data, 0);
	}

	public static int bytesToInt(byte[] data, int startPos) {
	    if (data == null || data.length < 4 || startPos > (data.length-4)) return 0x0;
	    // ----------
	    return (int)( // NOTE: type cast not necessary for int
	            (0xff & data[startPos]) << 24  |
	            (0xff & data[startPos + 1]) << 16  |
	            (0xff & data[startPos + 2]) << 8   |
	            (0xff & data[startPos + 3])
	            );
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
	    if (data == null || data.length != 8) return 0x0;
	    // ----------
	    return bytesToLong(data, 0);
	}

	public static long bytesToLong(byte[] data, int startPos) {
	    if (data==null || data.length<8 || startPos>(data.length-8)) return 0x0;

	    // This block is the wrong code
//	    return (long)(
//	            (0xff & data[startPos]) << 56  |
//	            (0xff & data[startPos+1]) << 48  |
//	            (0xff & data[startPos+2]) << 40	|
//	            (0xff & data[startPos+3]) << 32	|
//	            (0xff & data[startPos+4]) << 24  |
//	            (0xff & data[startPos+5]) << 16  |
//	            (0xff & data[startPos+6]) << 8   |
//	            (0xff & data[startPos+7])
//	            );

	    // This block is the working code
//	    return (
//	    		(long)(0xff & data[startPos]) << 56 |
//	    		(long)(0xff & data[startPos+1]) << 48 |
//	    		(long)(0xff & data[startPos+2]) << 40 |
//	    		(long)(0xff & data[startPos+3]) << 32 |
//	    		(long)(0xff & data[startPos+4]) << 24 |
//	    		(long)(0xff & data[startPos+5]) << 16 |
//	    		(long)(0xff & data[startPos+6]) << 8 |
//	    		(long)(0xff & data[startPos+7])
//	    		);

	    return ((long)(bytesToInt(data, startPos)) << 32) + (bytesToInt(data, startPos + 4) & 0xFFFFFFFFL);
	}





	public static byte[] reducedPrecisionXYGetBytes(float xy) {
		int data;
		if (xy>32.767 || xy<0) {
			System.err.println("Error - the value of X or Y is over the specified range of [0, 32.767].");
			return null;
		} else {
			data = (int) Math.round(xy * 1000);
			return new byte[] {(byte)((data >> 8) & 0x7f), (byte)(data & 0xff)};
		}

//		short data;
//		data = (short) Math.round(value * 1000);
//	    return new byte[] {
//		        (byte)((data >> 8) & 0xff),
//		        (byte)(data & 0xff)
//		    };

//	    int data;
//		data = Float.floatToRawIntBits(value);// & 0x007fffff;
//		
//		System.out.println(Float.floatToRawIntBits(value)+"\t"+data);
//		return new byte[] {(byte)((data >> 24) & 0xff), (byte)((data >> 16) & 0xff)};
	}

	public static void reducedPrecisionXYGetBytes(float xy, byte[] array, int startPosition) {
		int data;
		if (xy>32.767 || xy<0) {
			System.err.println("Error - the value of X or Y is over the specified range of [0, 32.767].");
		} else {
			data = (int) Math.round(xy * 1000);
			array[startPosition] = (byte)((data >> 8) & 0x7f);
			array[startPosition+1] = (byte)(data & 0xff);
		}
	}
	
	public static byte reducedPrecisionXYGetBytes1(float xy, byte[] array, int startPosition) {
		int data;
		if (xy>32.767 || xy<0) {
			array[startPosition] = (byte)((1 << 7) & 0x80);
			array[startPosition+1] = (byte) 0;
			return (byte)-1;
		} else {
			data = (int) Math.round(xy * 1000);
			array[startPosition] = (byte)((data >> 8) & 0x7f);
			array[startPosition+1] = (byte)(data & 0xff);
			return (byte)0;
		}
	}
	
	public static void reducedPrecisionXYGetBytes(float xy, byte[] array, int startPosition, Vector<Byte> outOfRangeValues) {
		int data;
		byte[] dataArray;
		int location;

		if (xy>=0 && xy<=32.767) {
			data = (int) Math.round(xy * 1000);
			array[startPosition] = (byte)((data >> 8) & 0x7f);
			array[startPosition+1] = (byte)(data & 0xff);
		} else {
			location = outOfRangeValues.size()/4;
			if (location>32767 || xy<0) {
				System.err.println("Error - too many out-of-range X or Y values. Can only accept 32,767 per sample file.");
			} else {
				array[startPosition] = (byte)(((1 << 7) & 0x80) | ((location >> 8) & 0x7f));
				array[startPosition+1] = (byte)(location & 0xff);
				dataArray = floatToBytes(xy);
				outOfRangeValues.add(dataArray[0]);
				outOfRangeValues.add(dataArray[1]);
				outOfRangeValues.add(dataArray[2]);
				outOfRangeValues.add(dataArray[3]);
			}
		}
	}
	
	public static float reducedPrecisionXYGetFloat(byte[] data) {
	    if (data == null || data.length != 2) {
	    	return (Float) null;
	    } else {
	    	return (int)((data[0] << 8 & 0x7f00) | (0xff & data[1])) / (float)1000;
	    }


//	    if (data == null || data.length != 2) {
//	    	return (Float) null;
//	    } else {
////		    return (float)((0xff & data[0]) << 8 | (0xff & data[1])) / (float)1000;
//	    	return (float)(data[0] << 8 | (0xff & data[1])) / (float)1000;
////	    	return (float)(data[0] << 8 | (0xff & data[1]));
//	    }
	}

	public static float reducedPrecisionXYGetFloat1(byte[] data) {
	    if (data == null || data.length != 2) {
	    	return (Float) null;
	    } else if (((data[0] & 0x80) >> 7)!=1) {
	    	return (int)((data[0] << 8 & 0x7f00) | (0xff & data[1])) / (float)1000;
	    } else {
	    	return -1;
	    }
	}

	public static float reducedPrecisionXYGetFloat(byte[] data, float[] outOfRangeValues) {
	    if (data == null || data.length != 2) {
	    	return (Float) null;
	    } else if (((data[0] & 0x80) >> 7)!=1) {
	    	return (int)((data[0] << 8 & 0x7f00) | (0xff & data[1])) / (float)1000;
	    } else {
	    	return outOfRangeValues[(int)((data[0] << 8 & 0x7f00) | (0xff & data[1]))];
	    }
	}
	
	public static byte[] reducedPrecisionXYGetBytes2(float xy) {
		int data;
		if (Float.isNaN(xy)) {
			return REDUCED_PRECISION_XY_NAN_BYTES;
		} else if (xy>65.533 || xy<0) {
			System.err.println("Note: the value of X or Y, " + xy + ", is over the specified range of [0, 65.533]. Thus the byte[] for Out_Of_Range value is returned.");
			return REDUCED_PRECISION_XY_OUT_OF_RANGE_BYTES;
		} else {
			data = (int) Math.round(xy * 1000);
			return new byte[] {(byte)((data >> 8) & 0xff), (byte)(data & 0xff)};
		}
	}

	public static byte reducedPrecisionXYGetBytes2(float xy, byte[] array, int startPosition) {
		int data;
		if (Float.isNaN(xy)) {
			array[startPosition] = REDUCED_PRECISION_XY_NAN_BYTES[0];
			array[startPosition+1] = REDUCED_PRECISION_XY_NAN_BYTES[1];
			return 0;
		} else if (xy>65.533 || xy<0) {
			array[startPosition] = REDUCED_PRECISION_XY_OUT_OF_RANGE_BYTES[0];
			array[startPosition+1] = REDUCED_PRECISION_XY_OUT_OF_RANGE_BYTES[1];
			return (byte)-1;
		} else {
			data = (int) Math.round(xy * 1000);
			array[startPosition] = (byte)((data >> 8) & 0xff);
			array[startPosition+1] = (byte)(data & 0xff);
			return (byte)0;
		}
	}

	public static float reducedPrecisionXYGetFloat2(byte[] data) {
	    if (data == null || data.length != 2) {
	    	return (Float) null;
	    } else if ( data[0]==REDUCED_PRECISION_XY_NAN_BYTES[0] && data[1]==REDUCED_PRECISION_XY_NAN_BYTES[1] ) {
	    	return Float.NaN;
	    } else if ( data[0]==REDUCED_PRECISION_XY_OUT_OF_RANGE_BYTES[0] && data[1]==REDUCED_PRECISION_XY_OUT_OF_RANGE_BYTES[1] ) {
	    	return -1;
	    } else {
	    	return (int)(((byte) 0 << 16 & 0xff0000) | (data[0] << 8 & 0xff00) | (0xff & data[1])) / (float)1000;
	    }
	}





	public static byte[] reducedPrecisionGcBafGetBytes(float gcOrBaf) {
		short data;
		
		data = (short) Math.round(gcOrBaf * 10000);
				
	    return new byte[] {
		        (byte)((data >> 8) & 0xff),
		        (byte)(data & 0xff)
		    };
	}

	public static void reducedPrecisionGcBafGetBytes(float gcOrBaf, byte[] array, int startPosition) {
		short data;
		
		data = (short)  Math.round(gcOrBaf * 10000);
		
		array[startPosition] = (byte)((data >> 8) & 0xff);
		array[startPosition+1] = (byte)(data & 0xff);
	}
	
	public static float reducedPrecisionGcBafGetFloat(byte[] data) {
	    if (data == null || data.length != 2) {
	    	return (Float) null;
	    } else {
//		    return (float)((0xff & data[0]) << 8 | (0xff & data[1])) / (float)1000;
	    	return (float)(data[0] << 8 | (0xff & data[1])) / (float)10000;
	    }
	}

	/**
	 * 
	 * @param gcOrBaf is the float value to be converted into byte[]. It must be in the range of 0.0000 through 1.0000.
	 * @return
	 */
	public static byte[] reducedPrecisionGcBafGetBytes2(float gcOrBaf) {
		if (Float.isNaN(gcOrBaf)) {
			return REDUCED_PRECISION_GCBAF_NAN_BYTES;
		} else if (gcOrBaf>1.000 || gcOrBaf<0) {
			System.out.println("Error: the value of GC or BAF should be in the range of 0.0000 through 1.0000. Please correct your data and re-run the program.");
			return null;
		} else {
			short data;
			data = (short) Math.round(gcOrBaf * 10000);
		    return new byte[] { (byte)( ((byte)0 & 0xc0) | ((data >> 8) & 0x3f) ), (byte)(data & 0xff)};
		}
	}

	public static void reducedPrecisionGcBafGetBytes2(float gcOrBaf, byte[] array, int startPosition) {
		if (Float.isNaN(gcOrBaf)) {
			array[startPosition] = REDUCED_PRECISION_GCBAF_NAN_BYTES[0];
			array[startPosition+1] = REDUCED_PRECISION_GCBAF_NAN_BYTES[1];
		} else if (gcOrBaf>1.000 || gcOrBaf<0) {
			System.out.println("Error: the value of GC or BAF should be in the range of 0.0000 through 1.0000. Please correct your data and re-run the program.");
		} else {
			short data;
			data = (short) Math.round(gcOrBaf * 10000);
			array[startPosition] = (byte)( ((byte)0 & 0xc0) | ((data >> 8) & 0x3f) );
			array[startPosition+1] = (byte)(data & 0xff);
		}
	}
	
	public static float reducedPrecisionGcBafGetFloat2(byte[] data) {
	    if (data == null || data.length != 2) {
	    	return (Float) null;
	    } else if (data[0]==REDUCED_PRECISION_GCBAF_NAN_BYTES[0] && data[1]==REDUCED_PRECISION_GCBAF_NAN_BYTES[1]) {
	    	return Float.NaN;
	    } else {
	    	return (float)(data[0] << 8 | (0xff & data[1])) / (float)10000;
	    }
	}




	public static byte[] reducedPrecisionLrrGetBytes(float lrr) {
//		int data;
//		if (value>=1677.7215 || value<=-2204.8384) {
//			System.out.println("Error - the Log R Ratio value is out of the range of [-204.8384, 204.8383]");
//			return null;
//		} else {
//			data = (int) (value * 10000);
//			return new byte[] {(byte)((data >> 16) & 0xff), (byte)((data >> 8) & 0xff), (byte)(data & 0xff)};
//		}

		int data;
		if (lrr<(float)-838.8608 || lrr>(float)838.8607) {
			System.out.println("Error - the LRR (Log R Ratio value) is out of the range of [-838.8608, 838.8607]");
			return null;
		} else {
			data = (int) Math.round(lrr * 10000);
			return new byte[] {(byte)((data >> 24 & 0x80) | ((data >> 16) & 0x7f)), (byte)((data >> 8) & 0xff), (byte)(data & 0xff)};
		}
	}

	public static void reducedPrecisionLrrGetBytes(float lrr, byte[] array, int startPosition) {

		int data;
		if (lrr<=(float)-838.8608 || lrr>=(float)838.8607) {
			System.out.println("Error - the LRR (Log R Ratio) value is out of the range of [-838.8608, 838.8607]");
		} else {
			data = (int) Math.round(lrr * 10000);
			array[startPosition] = (byte)((data >> 24 & 0x80) | ((data >> 16) & 0x7f));
			array[startPosition+1] = (byte)((data >> 8) & 0xff);
			array[startPosition+2] = (byte)(data & 0xff);
		}
	}
	
	public static float reducedPrecisionLrrGetFloat(byte[] data) {
	    if (data == null || data.length != 3) {
	    	return (Float) null;
	    } else if (data[0]>=0) {
	    	return (int)(((0x7f & data[0]) << 16) | ((0xff & data[1]) << 8) | (0xff & data[2])) / (float)10000;
	    } else {
	    	return (int)((0xff000000) | ((0xff & data[0]) << 16) | ((0xff & data[1]) << 8) | (0xff & data[2])) / (float)10000;
	    }
//	    if (data == null || data.length != 3) {
//	    	return (Float) null;
//	    } else {
//	    	return (int)(((0xff & data[0]) << 16) | ((0xff & data[1]) << 8) | (0xff & data[2])) / (float)10000;
//	    }
	}

	public static byte[] reducedPrecisionLrrGetBytes2(float lrr) {
		int data;
		if (Float.isNaN(lrr)) {
			return REDUCED_PRECISION_LRR_NAN_BYTES;
		} else if (lrr<(float)-13.1071 || lrr>(float)13.1071) {
			System.out.println("Error - the LRR (Log R Ratio value) is out of the range of -13.1072 through 13.1069");
			return REDUCED_PRECISION_LRR_OUT_OF_RANGE_BYTES;
		} else {
			data = (int) Math.round(lrr * 10000);
//			return new byte[] {(byte)((data >> 24 & 0x80) | ((data >> 16) & 0x7f)), (byte)((data >> 8) & 0xff), (byte)(data & 0xff)};
			return new byte[] {(byte)((data >> 30 & 0x02) | ((data >> 16) & 0x1)), (byte)((data >> 8) & 0xff), (byte)(data & 0xff)};
		}
	}

	public static byte reducedPrecisionLrrGetBytes2(float lrr, byte[] array, int startPosition) {
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
			data = (int) Math.round(lrr * 10000);
			array[startPosition] = (byte)((data >> 30 & 0x2) | ((data >> 16) & 0x1));
			array[startPosition+1] = (byte)((data >> 8) & 0xff);
			array[startPosition+2] = (byte)(data & 0xff);
			return 0;
		}
	}
	
	public static float reducedPrecisionLrrGetFloat2(byte[] data) {
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








	public static byte reducedPrecisionGenotypeGetBytes(byte abGenotype, byte forwardGenotype) {
//	    return (byte) ((byte)((forwardGenotype << 3) & 0xfc) | (abGenotype & 0x03));
		if (abGenotype>=0) {
			return (byte) (((forwardGenotype << 3) & 0xf8) | (abGenotype & 0x03));
		} else {
			return (byte) (((forwardGenotype << 3) & 0xf8) | ((abGenotype & 0x80) >> 5) | ((~abGenotype +1) & 0x03));
		}
	}

	public static void reducedPrecisionGenotypeGetBytes(byte abGenotype, byte forwardGenotype, byte[] array, int startPosition) {
//		array[startPosition] = (byte) ((byte)((forwardGenotype << 3) & 0xfc) | (abGenotype & 0x03));
		if (abGenotype>=0) {
			array[startPosition] = (byte) (((forwardGenotype << 3) & 0xf8) | (abGenotype & 0x03));
		} else {
			array[startPosition] = (byte) (((forwardGenotype << 3) & 0xf8) | ((abGenotype & 0x80) >> 5) | ((~abGenotype +1) & 0x03));
		}
	}
	
	public static byte[] reducedPrecisionGenotypeGetTypes(byte data) {
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
