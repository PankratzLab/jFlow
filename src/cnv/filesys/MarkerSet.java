package cnv.filesys;

import java.io.*;
import java.util.*;

import javax.swing.JOptionPane;

import common.*;

public class MarkerSet implements Serializable {
	public static final long serialVersionUID = 1L;
	public static final char[] ALLELES = {'A', 'C', 'G', 'T', 'I', 'D'};

	private long fingerprint;
	private String[] markerNames;
	private byte[] chrs;
	private int[] positions;
//	private byte[][] alleles; // two alleles per marker, with the A allele in index 0 and the B allele in index 1, the value is the index in MarkerSet.ALLELES 

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
	
//	public void checkFingerprint(Sample_old samp) {
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
		int[] indices;
		MarkerLookup markerLookup;
		Vector<String> missingMarkers;

		jar = proj.getJarStatus();
		markerLookup = proj.getMarkerLookup();
		fingerprint = proj.getSampleList().getFingerprint();

		missingMarkers = new Vector<String>();
//		time = new Date().getTime();
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
//				System.err.println("Error - could not find "+markerNames[0]+" in the lookup table");
				missingMarkers.add(markerNames[0]);
			}
		}
		if (missingMarkers.size() > 0) {
			JOptionPane.showMessageDialog(null, "Error - the following markers were not found in the MarkerSet: "+Array.toStr(Array.toStringArray(missingMarkers), " "), "Error", JOptionPane.ERROR_MESSAGE);
		}
		
		time = new Date().getTime();
		keys = HashVec.getKeys(hash);
		markerData = new MarkerData[markerNames.length];
		for (int i = 0; i<keys.length; i++) {
			v = hash.get(keys[i]);
			if (!Files.exists(proj.getDir(Project.PLOT_DIRECTORY)+keys[i], jar)) {
				JOptionPane.showMessageDialog(null, "Error - could not load data from '"+proj.getDir(Project.PLOT_DIRECTORY)+keys[i]+"'; because the file could not be found", "Error", JOptionPane.ERROR_MESSAGE);
				return null;
			}
			collection = MarkerDataCollection.load(proj.getDir(Project.PLOT_DIRECTORY)+keys[i], jar).getCollection();
			for (int j = 0; j<v.size(); j++) {
				line = v.elementAt(j).split("[\\s]+");
				indices = ext.indicesOfStr(line[0], markerNames, true, true);
				
				if (indices.length==0) {
					System.err.println("Error - How can this be?");
				}
				for (int k = 0; k<indices.length; k++) {
					try {
						markerData[indices[k]] = collection[Integer.parseInt(line[1])];
					} catch (Exception e) {
						System.err.println("Error - failed to load data for marker '"+line[0]+"' which is in collection "+line[1]+" may need to regenerate the markerLookup file");
					}
					
					if (markerData[indices[k]].getFingerprint()!=fingerprint) {
						System.err.println("Error - mismatched fingerprint after MarkerLookup. Actual in MarkerData: " + markerData[indices[k]].getFingerprint() + ", while expecting: " + fingerprint);
					}					
				}
			}
		}

//		System.out.println("Finished loading MarkerData in "+ext.getTimeElapsed(time));

		return markerData;
	}
	
	public byte[] translateABtoForwardGenotypes(byte[] abGenotypes, char[][] abLookup) {
		byte[] result = new byte[abGenotypes.length];
		String geno;
		
		for (int i=0; i<abGenotypes.length; i++) {
			switch (abGenotypes[i]) {
			case 0:
				geno = abLookup[i][0]+""+abLookup[i][0];
				break;
			case 1:
				geno = abLookup[i][0]+""+abLookup[i][1];
				break;
			case 2:
				geno = abLookup[i][1]+""+abLookup[i][1];
				break;
			case -1:
				geno = Sample.ALLELE_PAIRS[0];
				break;
			default:
				System.err.println("Error - invalid AB genotype: "+abGenotypes[i]);
				geno = null;
			}
//			System.out.println(geno);
			for (byte j=0; j<Sample.ALLELE_PAIRS.length; j++) {
				if (geno.equals(Sample.ALLELE_PAIRS[j])) {
					result[i]=j;
				}
			}
		}
		
		return result;
	}
}
