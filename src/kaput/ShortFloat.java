package kaput;

import common.ext;

public class ShortFloat {
	public static short toShort(String str) {
		boolean negative = false;
		int index;
		String num = str;
		short s;

		if (str.equals("NaN")) {
			return Short.MIN_VALUE;
		}
		if (num.startsWith("-")) {
			negative = true;
			num = num.substring(1);
		}
		index = num.indexOf(".");
		if (index==-1) {
			index = num.length();
			if (num.length()>3) {
				System.err.println("Error - cannot parse '"+str+"' to short (too many digits)");
				System.exit(1);
				return -999;
			}
		} else if (index==1&&num.startsWith("0.")) {
			index = 0;
			num = num.substring(2);
		} else {
			num = num.substring(0, index)+num.substring(index+1);
		}
		while (num.length()<4) {
			num += "0";
		}
		if (num.length()>4) {
			String hoo = num.substring(0, 4)+"."+num.charAt(4);
			num = ext.formNum(Math.round(Float.parseFloat(hoo)), 4);
		}

		try {
			s = Short.parseShort((negative?"-":"")+index+num);
		} catch (Exception e) {
			System.err.println("Error - cannot parse '"+str+"' to short:");
			e.printStackTrace();
			System.exit(1);
			return -999;
		}
		return s;
	}

	public static float shortToFloat(short s) {
		boolean negative = false;
		String num = s+"";
		int index;

		if (s==Short.MIN_VALUE) {
			return Float.NaN;
		}
		if (num.startsWith("-")) {
			negative = true;
			num = num.substring(1);
		}
		if (num.length()==5) {
			index = Integer.parseInt(num.substring(0, 1));
			num = num.substring(1);
		} else {
			index = 0;
		}
		return (negative?-1:1)*Float.parseFloat(num)*(float)Math.pow(0.1, 4-index);
	}

	// public static float[] loadBinaryToMatrix(String filename) {
	// BufferedReader reader;
	// PrintWriter writer;
	// String[] line;
	// String temp, trav;
	// Hashtable<String, String> hash = new Hashtable<String, String>();
	// Vector<String> v = new Vector<String>();
	// int count;
	//	
	// float[] myarray = new float[100000];
	// try {
	// File file = new File("data.bin");
	// InputStream is = new FileInputStream(file);
	// DataInputStream dis = new DataInputStream( is );
	// long length = file.length();
	// if (length > Integer.MAX_VALUE) {
	// throw new IOException("File is too large");
	// } else {
	// byte[] bytes = new byte[(int)length];
	// int offset = 0;
	// int numRead = 0;
	// while (offset < bytes.length &&
	// (numRead = is.read(bytes, offset, bytes.length-offset) ) >= 0) {
	// offset += numRead;
	// }
	// if (offset < bytes.length) {
	// throw new IOException("Could not completely read file "+file.getName());
	// }
	// dis.close();
	// is.close();
	// System.out.println("offset="+offset);
	//
	// for (int start = 0; start < offset; start = start + 4) {
	// // myarray[cnt] = arr2float(bytes, start);
	// // cnt++;
	// }
	// }
	// } catch (Exception e) {
	// e.printStackTrace();
	// }
	//
	// return new float[0];
	// }
	//
	// public static float arr2float (byte[] arr, int start) {
	// int i = 0;
	// int len = 4;
	// int cnt = 0;
	// byte[] tmp = new byte[len];
	// for (i = start; i < (start + len); i++) {
	// tmp[cnt] = arr[i];
	// cnt++;
	// }
	// int accum = 0;
	// i = 0;
	// for ( int shiftBy = 0; shiftBy < 32; shiftBy += 8 ) {
	// accum |= ( (long)( tmp[i] & 0xff ) ) << shiftBy;
	// i++;
	// }
	// return Float.intBitsToFloat(accum);
	// }
	//
	// public static void writeFloatMatrix(float[][] matrix, String filename) {
	// DataOutputStream dos;
	//	
	// try {
	// dos = new DataOutputStream(new FileOutputStream(filename));
	// dos.writeInt(matrix.length);
	// dos.writeInt(matrix[0].length);
	// for (int i = 0; i < matrix.length; i++) {
	// if (matrix[i].length != matrix[0].length) {
	// System.err.println("Error - all rows in matrix do not have the same size;
	// failed ot write '"+filename+"'");
	// System.exit(1);
	// }
	// for (int j = 0; j < matrix[i].length; j++) {
	// dos.writeFloat(matrix[i][j]);
	// }
	// }
	// dos.close();
	// } catch (Exception e) {
	// e.printStackTrace();
	// }
	// }
	//
	// public static float[][] readFloatMatrix(String filename) {
	// DataInputStream dis;
	// float[][] matrix = null;
	// int n, m;
	//	
	// try {
	// dis = new DataInputStream(new FileInputStream(filename));
	// n = dis.readInt();
	// m = dis.readInt();
	// matrix = new float[n][m];
	// for (int i = 0; i < n; i++) {
	// for (int j = 0; j < m; j++) {
	// matrix[i][j] = dis.readFloat();
	// }
	// }
	// dis.close();
	// } catch (Exception e) {
	// e.printStackTrace();
	// }
	//	
	// return matrix;
	// }
	//
	// public static void blowme() {
	// BufferedReader reader;
	// PrintWriter writer;
	// String[] line;
	// String temp, trav;
	// Hashtable<String, String> hash = new Hashtable<String, String>();
	// Vector<String> v = new Vector<String>();
	// int count;
	// float[][] matrix;
	//	
	// String filename = "ND13115.ori";
	// try {
	// System.out.println(ext.getTime()+"\tLoading original file");
	// reader = new BufferedReader(new FileReader(filename));
	// matrix = new float[370364][2];
	// count = 0;
	// while (reader.ready()) {
	// line = reader.readLine().trim().split("[\\s]+");
	// matrix[count][0] = Float.parseFloat(line[0]);
	// matrix[count][1] = Float.parseFloat(line[1]);
	// count++;
	// }
	// reader.close();
	//		
	// // System.out.println(ext.getTime()+"\tSaving data to binary stream");
	// // writeFloatMatrix(matrix, filename+".bin");
	// //
	// // matrix = null;
	// //
	// // System.out.println(ext.getTime()+"\tLoading data from binary stream");
	// // matrix = readFloatMatrix(filename+".bin");
	//		
	// System.out.println(ext.getTime()+"\tWriting parsed data to text file");
	// DecimalFormat f = new DecimalFormat("0.####");
	//
	// try {
	// writer = new PrintWriter(new FileWriter(filename+".new"));
	// for (int i = 0; i < matrix.length; i++) {
	// writer.println(f.format(matrix[i][0])+"\t"+f.format(matrix[i][1]));
	// }
	// writer.close();
	// } catch (FileNotFoundException fnfe) {
	// System.err.println("Error: file \""+filename+".new"+"\" not found in
	// current directory");
	// System.exit(1);
	// } catch (IOException ioe) {
	// System.err.println("Error reading file \""+filename+".new"+"\"");
	// System.exit(2);
	// }
	//
	// System.out.println(ext.getTime()+"\tDone!");
	//
	// } catch (FileNotFoundException fnfe) {
	// System.err.println("Error: file \""+filename+"\" not found in current
	// directory");
	// System.exit(1);
	// } catch (IOException ioe) {
	// System.err.println("Error reading file \""+filename+"\"");
	// System.exit(2);
	// }
	//	
	// }

}
