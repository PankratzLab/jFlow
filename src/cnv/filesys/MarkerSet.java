package cnv.filesys;

import java.io.*;
import java.util.ArrayList;

import common.*;
import filesys.Segment;

public class MarkerSet implements Serializable {
	public static final long serialVersionUID = 1L;
	public static final char[] ALLELES = {'A', 'C', 'G', 'T', 'I', 'D'};
	// TODO these alleles were recently added, should they not be from some place else??

	private long fingerprint;
	private String[] markerNames;
	private byte[] chrs;
	private int[] positions;
//	private char[][] abAlleles;

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

	public int[][] getIndicesByChr() {
		IntVector iv;
		byte chr;
		int[][] indicesByChr;
		boolean done;
		
		indicesByChr = new int[27][0];
		
		chr = 0;
		iv = new IntVector(20000);
		done = false;
		for (int i = 0; !done; i++) {
			if (i==chrs.length || chrs[i] != chr) {
				indicesByChr[chr] = iv.toArray();
				chr = i==chrs.length?0:chrs[i];
				iv = new IntVector(20000);
			}
			if (i==chrs.length) {
				done = true;
			} else {
				iv.add(i);
			}
        }
		
		return indicesByChr;
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

	public int[] getIndicesOfMarkersIn(Segment seg, int[][] indicesByChr, Logger log) {
		return ext.indexLargeFactors(getMarkersIn(seg, indicesByChr), markerNames, true, log, true, false);
	}
	
	
	public String[] getMarkersIn(Segment seg, int[][] indicesByChr) {
		int index = seg.getChr();
		ArrayList<String> markersIn = new ArrayList<String>();
		int[] indices = indicesByChr == null ? getIndicesByChr()[index] : indicesByChr[index];
		for (int i = 0; i < indices.length; i++) {
			int bp = positions[indices[i]];

			if (bp >= seg.getStart() && bp <= seg.getStop()) {
				markersIn.add(markerNames[indices[i]]);
			}
			if (bp > seg.getStop()) {
				break;
			}
		}
		return Array.toStringArray(markersIn);
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
		// TODO remove completely
		return null;
	}

	public static byte[] translateABtoForwardGenotypes(byte[] abGenotypes, char[][] abLookup) {
		byte[] result;
		String geno;
		
		result = Array.byteArray(abGenotypes.length, (byte)-3);
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
