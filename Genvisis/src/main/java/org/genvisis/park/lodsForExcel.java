package org.genvisis.park;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.HashVec;
import org.genvisis.common.IntVector;
import org.genvisis.common.ext;

import com.google.common.primitives.Ints;

public class lodsForExcel {
	public static final String BLANK = "      ";

	public static final String[][] DATA_POINTS = {{"Mapmaker", "Mapmaker	no Dv	Dv"}, {"Dv"},
																								{"Aspex"}, {"Aspex 2pt"},
																								{"Merlin-siblings", "Merlin-sibpairs"},
																								{"Merlin-extended", "Merlin-extended"},
																								{"Genehunter", "Allegro	linear"},
																								{"Dominant", "autosomal dominant"},
																								{"Dom Het", "autosomal dominant het"},
																								{"Recessive", "autosomal recessive"},
																								{"Rec Het", "autosomal recessive het"}};

	public static final String[][] CHR_OFFSETS = {{"D1S468", "4.2"}, {"D2S319", "7.6"},
																								{"D3S1297", "8.3"}, {"D4S412", "4.7"},
																								{"D5S1981", "1.7"}, {"D6S1574", "9.2"},
																								{"D7S531", "5.3"}, {"D8S264", "0.7"},
																								{"D9S288", "9.8"}, {"D10S249", "2.1"},
																								{"D11S4046", "2.8"}, {"D12S352", "0.0"},
																								{"D13S175", "6.0"}, {"D14S261", "6.5"},
																								{"D15S128", "6.1"}, {"D16S423", "10.4"},
																								{"D17S849", "0.6"}, {"D18S59", "0.0"},
																								{"D19S209", "11.0"}, {"D20S117", "2.8"},
																								{"D21S1256", "9.7"}, {"D22S420", "4.1"},
																								{"DXS1060", "15.1"}, {"rs884080", "0.0"},
																								{"rs381726", "1.9"}, {"rs1516337", "1.0"},
																								{"rs963598", "1.4"}, {"rs413666", "0.7"},
																								{"rs719065", "0.0"}, {"rs1881114", "3.5"},
																								{"rs13429", "0.0"}, {"rs1532309", "0.0"},
																								{"rs1476130", "1.2"}, {"rs741737", "0.0"},
																								{"rs476646", "0.0"}, {"rs1838114", "0.0"},
																								{"rs1972373", "0.6"}, {"rs1562203", "0.0"},
																								{"rs8466", "1.1"}, {"rs1609550", "0.1"},
																								{"rs948263", "0.0"}, {"rs1020382", "0.0"},
																								{"rs371791", "0.5"}, {"rs990141", "6.0"},
																								{"rs7288876", "0.0"}, {"rs749706", "11.7"}};

	public class eStruct {
		public int position;

		public String markerName;

		public double[] values;

		public eStruct() {
			markerName = BLANK;
			values = Array.doubleArray(DATA_POINTS.length, -999);
		}
	}

	public lodsForExcel(int[] chrs) throws IOException {
		BufferedReader reader;
		PrintWriter writer;
		PrintWriter all;
		String temp, chrome;
		String[] line;
		eStruct handle;
		String pos;
		String[] poslar;
		Hashtable<String, eStruct> hash = new Hashtable<String, eStruct>();
		int trav, offset;
		int count, num;
		double prev, increment;
		Vector<String> mapFailure = new Vector<String>();
		int numMarkers;
		Vector<String> markers = new Vector<String>();
		double posi, off;

		all = new PrintWriter(new FileWriter("excel-all.xls"));
		all.println("Marker\tPosition\tAdj.Position\t" + Array.toStr(DATA_POINTS, "\t"));

		for (int chr : chrs) {
			chrome = ext.chrome(chr);
			hash.clear();

			reader = new BufferedReader(new FileReader("summary" + chrome + ".out"));
			trav = -1;
			while (reader.ready()) {
				temp = reader.readLine();
				if (temp.equals("")) {
					trav = -1;
				} else if (trav == -1) {
					trav = ext.indexOfAnyStr(temp, DATA_POINTS, false);
					if (trav == -1) {
						System.err.println("Error - failed to index '" + temp + "' as a valid program run");
						System.exit(1);
					}
				} else if (!temp.startsWith("no ")) {
					line = temp.split("[\\s]+");
					offset = 0;
					if (DATA_POINTS[trav][0].startsWith("Aspex")) {
						offset = 1;
					}
					pos = Math.round(Double.valueOf(line[0 + offset]).doubleValue()) + "";
					if (hash.containsKey(pos)) {
						handle = hash.get(pos);
					} else {
						hash.put(pos, handle = new eStruct());
					}
					if (DATA_POINTS[trav][0].startsWith("Aspex")) {
						handle.markerName = line[0];
					}
					handle.values[trav] = Double.parseDouble(line[1 + offset]);
					if (DATA_POINTS[trav][0].startsWith("Mapmaker")) {
						handle.values[trav + 1] = Double.parseDouble(line[1 + 1 + offset]);
					}

				}
			}
			reader.close();

			poslar = HashVec.getNumericKeys(hash);
			for (int j = Integer.parseInt(poslar[0]); j <= Integer.parseInt(poslar[poslar.length
																																							- 1]); j++) {
				if (!hash.containsKey(j + "")) {
					hash.put(j + "", new eStruct());
				}
			}

			try {
				try {
					reader = new BufferedReader(new FileReader("map" + chrome + ".dat"));
				} catch (Exception e) {
					reader = new BufferedReader(new FileReader("../map" + chrome + ".dat"));
				}
				numMarkers = Integer.parseInt(reader.readLine().split("[\\s]+")[0]) - 1;
				markers.removeAllElements();
				for (int j = 0; j < (chr > 22 ? 7 : 6); j++) {
					reader.readLine();
				}
				for (int j = 0; j < numMarkers; j++) {
					markers.add(reader.readLine().split("[\\s]+")[3]);
					reader.readLine();
				}
				reader.readLine();
				line = reader.readLine().split("[\\s]+");
				posi = 0;
				for (int j = 0; j < numMarkers; j++) {
					handle = hash.get(Math.round(posi) + "");
					if (handle.markerName.equals(BLANK)) {
						handle.markerName = markers.elementAt(j);
					} else {
						handle.markerName += "/" + markers.elementAt(j);
					}
					if (j < numMarkers - 1) {
						posi += Double.parseDouble(line[1 + j]);
					}
				}
				reader.close();
			} catch (Exception e) {
				mapFailure.add(chr + "");
			}

			poslar = HashVec.getNumericKeys(hash);
			for (int meth = 0; meth < DATA_POINTS.length; meth++) {
				count = Integer.parseInt(poslar[0]);
				prev = -999;
				while (count <= Integer.parseInt(poslar[poslar.length - 1])) {
					handle = hash.get(count + "");
					if (handle.values[meth] == -999 && prev != -999) {
						num = 0;
						do {
							num++;
							handle = hash.get((count + num) + "");
						} while (count + num <= Integer.parseInt(poslar[poslar.length - 1])
											&& handle.values[meth] == -999);
						if (count + num > Integer.parseInt(poslar[poslar.length - 1])) {
							count = Integer.parseInt(poslar[poslar.length - 1]);
						} else {
							increment = (handle.values[meth] - prev) / (num + 1);
							for (int j = 0; j < num; j++) {
								handle = hash.get(count + "");
								handle.values[meth] = prev + increment;
								if (handle.values[meth] < 0 && handle.values[meth] > -0.000001) {
									handle.values[meth] = 0;
								}
								prev = handle.values[meth];
								count++;
							}
							count--;
						}
					} else if (handle.values[meth] != -999) {
						prev = handle.values[meth];
					}
					count++;
				}
			}

			poslar = HashVec.getNumericKeys(hash);
			try {
				writer = new PrintWriter(new FileWriter("excel" + chrome + ".xls"));
				writer.println("Marker\tPosition\tAdj.Position\t" + Array.toStr(DATA_POINTS, "\t"));

				handle = hash.get(poslar[0]);
				off = -1;
				for (String[] element : CHR_OFFSETS) {
					if (handle.markerName.startsWith(element[0])) {
						off = Double.parseDouble(element[1]);
					}
				}

				for (int j = Integer.parseInt(poslar[0]); j <= Integer.parseInt(poslar[poslar.length
																																								- 1]); j++) {
					handle = hash.get(j + "");
					writer.print(handle.markerName + "\t" + j + "\t" + (off == -1 ? "" : off + j));
					for (int meth = 0; meth < DATA_POINTS.length; meth++) {
						writer.print("\t" + (handle.values[meth] == -999 ? "" : handle.values[meth]));
					}
					writer.println();

					all.print(handle.markerName + "\t" + chr + "::" + j + "\t" + (off == -1 ? "" : off + j));
					for (int meth = 0; meth < DATA_POINTS.length; meth++) {
						all.print("\t" + (handle.values[meth] == -999 ? "" : handle.values[meth]));
					}
					all.println();
				}
				writer.close();
			} catch (Exception e) {
				System.out.println("Error in writing excel" + chrome + ".xls");
				e.printStackTrace();
			}

		}
		all.close();
		if (mapFailure.size() > 0) {
			System.err.println("Map information was not available for the following chromosomes: "
													+ Array.toStr(Array.toStringArray(mapFailure), " "));
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		IntVector iv = new IntVector();

		String usage = "\n"	+ "park.lodsForExcel requires 0+ arguments:\n"
										+ "   numbers of chromosomes to parse (i.e. 1 2 6 13 23 (nothing indicates the autosomes+23 (default)))\n"
										+ "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else {
				try {
					if (Integer.parseInt(arg) < 1 || Integer.parseInt(arg) > 25) {
						throw new NumberFormatException();
					}
					iv.add(Integer.parseInt(arg));
				} catch (NumberFormatException nfe) {
					System.err.println("Error - " + arg + " is not a valid chromosome number");
				}
				numArgs--;
			}
		}

		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			if (iv.size() == 0) {
				new lodsForExcel(Array.toIntArray(Array.stringArraySequence(23, "")));
			} else {
				new lodsForExcel(Ints.toArray(iv));
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
