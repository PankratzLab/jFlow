package cnv.filesys;

import java.io.*;
import java.util.*;
import common.*;

public class MarkerSet implements Serializable {
	public static final long serialVersionUID = 1L;

	private long fingerprint;
	private String[] markerNames;
	private byte[] chrs;
	private int[] positions;

	public MarkerSet(String[] markerNames, byte[] chrs, int[] positions) {
		if (markerNames.length!=chrs.length||markerNames.length!=positions.length) {
			System.err.println("Error - mismatched number of markers and positions");
			System.exit(1);
		}

		this.markerNames = markerNames;
		this.chrs = chrs;
		this.positions = positions;
		this.fingerprint = fingerprint(markerNames);
	}

	public MarkerSet(String[] rawMarkerNames, byte[] rawChrs, int[] rawPositions, int[] keys) {
		if (rawMarkerNames.length!=rawChrs.length||rawMarkerNames.length!=rawPositions.length||rawMarkerNames.length!=keys.length) {
			System.err.println("Error - mismatched number of markers and positions/keys");
			System.exit(1);
		}

		markerNames = new String[rawMarkerNames.length];
		chrs = new byte[rawChrs.length];
		positions = new int[rawPositions.length];

		for (int i = 0; i<keys.length; i++) {
			markerNames[i] = rawMarkerNames[keys[i]];
			chrs[i] = rawChrs[keys[i]];
			positions[i] = rawPositions[keys[i]];
		}

		fingerprint = fingerprint(markerNames);
	}

	public String[] getMarkerNames() {
		return markerNames;
	}

	public byte[] getChrs() {
		return chrs;
	}

	public int[] getPositions() {
		return positions;
	}

	public int[][] getPositionsByChr() {
		IntVector iv;
		byte chr;
		int[][] positionsByChr;
		boolean done;
		
		positionsByChr = new int[27][0];
		
		chr = 0;
		iv = new IntVector(20000);
		done = false;
		for (int i = 0; !done; i++) {
			if (i==chrs.length || chrs[i] != chr) {
				positionsByChr[chr] = iv.toArray();
				chr = i==chrs.length?0:chrs[i];
				iv = new IntVector(20000);
			}
			if (i==chrs.length) {
				done = true;
			} else {
				iv.add(positions[i]);
			}
        }
		
		return positionsByChr;
	}

	public long getFingerprint() {
		return fingerprint;
	}
	
	public void writeToFile(String filename) {
		PrintWriter writer;

		try {
			writer = new PrintWriter(new FileWriter(filename));
			for (int i = 0; i<markerNames.length; i++) {
				writer.println(markerNames[i]+"\t"+chrs[i]+"\t"+positions[i]);
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing "+filename);
			e.printStackTrace();
		}
	}

	public void serialize(String filename) {
		Files.writeSerial(this, filename);
	}

	public static MarkerSet load(String filename, boolean jar) {
		return (MarkerSet)Files.readSerial(filename, jar, true);
	}

	public static long fingerprint(String[] names) {
		long sum, trav;

		sum = 0;
		for (int i = 0; i<names.length; i++) {
			trav = 0;
			for (int j = 0; j<names[i].length(); j++) {
				trav += names[i].charAt(j);
			}
			sum += trav*(i+1);
		}

		return sum;
	}
	
	public void checkFingerprint(Sample samp) {
		if (samp.getFingerprint() != fingerprint) {
			System.err.println("Error - Sample has a different fingerprint ("+samp.getFingerprint()+") than the MarkerSet ("+fingerprint+")");
		}
	}

	public static void convert(String filename) {
		BufferedReader reader;
		String[] line;
		int count;
		String[] markerNames;
		byte[] chrs;
		int[] positions;

		try {
			reader = new BufferedReader(new FileReader(filename));
			count = 0;
			while (reader.ready()) {
				reader.readLine();
				count++;
			}
			reader.close();

			markerNames = new String[count];
			chrs = new byte[count];
			positions = new int[count];

			reader = new BufferedReader(new FileReader(filename));
			for (int i = 0; i<count; i++) {
				line = reader.readLine().trim().split("[\\s]+");
				markerNames[i] = line[0];
				chrs[i] = (byte)Positions.chromosomeNumber(line[1]);
				positions[i] = Integer.parseInt(line[2]);
			}
			reader.close();

			new MarkerSet(markerNames, chrs, positions).serialize(filename+".ser");
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+filename+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+filename+"\"");
			System.exit(2);
		}
	}

	public static MarkerData[] loadFromList(Project proj, String[] markerNames) {
		Hashtable<String,Vector<String>> hash = new Hashtable<String,Vector<String>>();
		String[] keys, line;
		Vector<String> v;
		long fingerprint, time;
		boolean jar;
		MarkerData[] collection, markerData;
		int index;
		MarkerLookup markerLookup;

		jar = proj.getJarStatus();
		markerLookup = proj.getMarkerLookup();
		fingerprint = proj.getSampleList().getFingerprint();

		time = new Date().getTime();
		for (int i = 0; i<markerNames.length; i++) {
			if (markerLookup.contains(markerNames[i])) {
				line = markerLookup.get(markerNames[i]).split("[\\s]+");
				if (hash.containsKey(line[0])) {
					v = hash.get(line[0]);
				} else {
					hash.put(line[0], v = new Vector<String>());
				}
				v.add(markerNames[i]+"\t"+line[1]);
			} else {
				System.err.println("Error - could not find "+markerNames[0]+" in the lookup table");
			}
		}
		keys = HashVec.getKeys(hash);
		markerData = new MarkerData[markerNames.length];
		for (int i = 0; i<keys.length; i++) {
			v = hash.get(keys[i]);
			collection = MarkerDataCollection.load(proj.getDir(Project.PLOT_DIRECTORY)+keys[i], jar).getCollection();
			for (int j = 0; j<v.size(); j++) {
				line = v.elementAt(j).split("[\\s]+");
				index = ext.indexOfStr(line[0], markerNames);
				if (index==-1) {
					System.err.println("Error - How can this be?");
				} else {
					markerData[index] = collection[Integer.parseInt(line[1])];
					if (markerData[index].getFingerprint()!=fingerprint) {
						System.err.println("Error - ");

					}
				}
			}
		}

		System.out.println("Finished loading MarkerData in "+ext.getTimeElapsed(time));

		return markerData;
	}
}
