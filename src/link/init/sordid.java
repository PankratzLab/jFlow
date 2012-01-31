package link.init;

import java.io.*;
import java.util.*;

import link.Markers;

public class sordid {

	public sordid() throws IOException {
		BufferedReader map, readers[];
		PrintWriter writers[];
		StringTokenizer st, markerToke;
		Vector<String> allMarkers = new Vector<String>();
		Vector<String> orderedMarkers = new Vector<String>();
		String[] unorderedMarkers;
		int[][] keys;
		int numMarkers;
		// String[] markers;
		String temp, mark, ID;
		Hashtable<String,String> hash = new Hashtable<String,String>();
		int[] chromosomalKey = new int[24];
		int chrome, allele;
		String[] markerdata;
		Vector<String> chrsPresent = new Vector<String>();
		int[] expectedAlleleCounts;

		int countFiles = 0;
		while ((new File("data"+(countFiles+1)+".txt")).exists()) {
			countFiles++;
		}
		readers = new BufferedReader[countFiles];
		expectedAlleleCounts = new int[countFiles];

		for (int i = 1; i<=23; i++) {
			chromosomalKey[i] = 0; // so we know which chromosomes are in the
			// current file
		}

		for (int count = 0; count<countFiles; count++) {
			try {
				readers[count] = new BufferedReader(new FileReader("data"+(count+1)+".txt"));
				readers[count].readLine();
				st = new StringTokenizer(readers[count].readLine());
				st.nextToken();
				st.nextToken();
				st.nextToken();
				expectedAlleleCounts[count] = st.countTokens()*2;
				while (st.hasMoreTokens()) {
					mark = st.nextToken().toUpperCase(); // each new marker
					allMarkers.add(mark);
					markerToke = new StringTokenizer(mark, "DS");
					if (markerToke.countTokens()!=2) {
						System.err.println("Error parsing marker "+mark);
						System.exit(1);
					}
					temp = markerToke.nextToken();
					if (temp.equals("X")) {
						chrome = 23;
					} else {
						chrome = Integer.valueOf(temp).intValue();
						if (chrome>23) {
							System.err.println("Error parsing marker "+mark);
							System.exit(1);
						}
					}
					if (chromosomalKey[chrome]>0) { // if the chromosome
						// hasn't been seen it's
						// added
						hash.put(chrome+"", hash.get(chrome+"")+" "+mark);
						chromosomalKey[chrome]++; // counts number of markers
						// for each chromosome
					} else {
						chromosomalKey[chrome] = 1;
						hash.put(chrome+"", mark);
						chrsPresent.add(chrome+"");
					}
				}

			} catch (Exception e) {
				System.err.println("Error with marker names in file: data"+(count+1)+".txt");
				e.printStackTrace();
				System.exit(1);
			}
		}

		writers = new PrintWriter[chrsPresent.size()];
		keys = new int[chrsPresent.size()][];

		for (int i = 0; i<chrsPresent.size(); i++) {
			chrome = Integer.valueOf(chrsPresent.elementAt(i)).intValue();
			orderedMarkers.removeAllElements();
			writers[i] = new PrintWriter(new FileWriter("chromosome"+chrome+".dat"));
			keys[i] = new int[chromosomalKey[chrome]];

			writers[i].println("placeholder line");
			writers[i].print("DNA\tFamNo\tIndNo");

			unorderedMarkers = new String[chromosomalKey[chrome]];
			st = new StringTokenizer(hash.get(chrome+""));
			for (int j = 0; j<chromosomalKey[chrome]; j++) {
				unorderedMarkers[j] = st.nextToken();
			}
			Markers.order(unorderedMarkers, true); // the great omnipotent
			// marker database
			map = new BufferedReader(new FileReader("markerMap.dat"));
			map.readLine();
			st = new StringTokenizer(map.readLine());
			orderedMarkers.add(st.nextToken());
			while (map.ready()) {
				map.readLine();
				st = new StringTokenizer(map.readLine());
				orderedMarkers.add(st.nextToken()); // retrieves ordered markers
			}
			map.close();
			numMarkers = orderedMarkers.size();
			for (int j = 0; j<numMarkers; j++) {
				keys[i][j] = allMarkers.indexOf(orderedMarkers.elementAt(j));
				if (keys[i][j]==-1) {
					System.err.println("markers in data file are fubar at marker # "+j+" couldn't find "+orderedMarkers.elementAt(j));
				}
			}
			for (int j = 0; j<keys[i].length; j++) {
				writers[i].print("\t"+allMarkers.elementAt(keys[i][j])+"\t");
			}
			writers[i].println();
		}

		while (readers[0].ready()) {
			markerdata = new String[allMarkers.size()*2];
			allele = 0;
			ID = "null";
			for (int count = 0; count<countFiles; count++) {
				try {
					temp = readers[count].readLine();
					if (temp.startsWith(" ")||temp.startsWith("\t")||temp.equals("")) {
						break;
					}
					st = new StringTokenizer(temp);
					temp = st.nextToken()+"\t"+st.nextToken()+"\t"+st.nextToken();
					if (ID.equals("null")) {
						ID = temp;
					} else if (!ID.equals(temp)) {
						System.err.println("Error - mismatch in data"+chrsPresent.elementAt(count)+".txt");
						System.err.println("      - when expecting "+ID+", saw "+temp);
						System.exit(1);
					}
					if (st.countTokens()!=expectedAlleleCounts[count]) {
						System.err.println("Error - "+temp+" has "+st.countTokens()+" alleles instead of "+expectedAlleleCounts[count]+" in file data"+(count+1)+".txt");
						System.exit(1);
					}
					while (st.hasMoreTokens()) {
						temp = st.nextToken();
						try {
							Integer.valueOf(temp).intValue();
						} catch (NumberFormatException nfe) {
							System.err.println("Replaced a '"+temp+"' with a '0' for data"+(count+1)+".txt");
							temp = "0";
						}
						markerdata[allele++] = temp;
					}
				} catch (Exception e) {
					System.err.println("Error with alleles in file: data"+(count+1)+".txt");
					e.printStackTrace();
					System.exit(1);
				}
			}
			for (int i = 0; i<chrsPresent.size(); i++) {
				chrome = Integer.valueOf(chrsPresent.elementAt(i)).intValue();
				writers[i].print(ID);
				for (int j = 0; j<keys[i].length; j++) {
					writers[i].print("\t"+markerdata[2*keys[i][j]]+"\t"+markerdata[2*keys[i][j]+1]);
				}
				writers[i].println();
			}
		}
		for (int count = 0; count<countFiles; count++) {
			if (readers[count].ready()) {
				System.err.println("Warning - it appears that data"+(count+1)+".txt had lines in addition to those previous.");
			}
			readers[count].close();
		}
		for (int i = 0; i<chrsPresent.size(); i++) {
			writers[i].close();
		}

	}

	public static void main(String[] args) throws IOException {
		if (args.length!=0) {
			System.out.println("I expect no arguments outta you!.");
		}
		try {
			new sordid();
		} catch (Exception e) {
			e.printStackTrace();
			System.err.println("Error in processing.");
		}

	}
}
