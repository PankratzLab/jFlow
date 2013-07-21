package cnv.var;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import common.ByteVector;
import common.Files;
import common.IntVector;
import common.Positions;
import common.StringVector;
import common.ext;

import filesys.Segment;

public class CNVariant extends Segment {
	public static final long serialVersionUID = 1L;
	public static final String[] PLINK_CNV_HEADER = { "FID", "IID", "CHR", "BP1", "BP2", "TYPE", "SCORE", "SITES" };

	private String familyID;
	private String individualID;
	private int cn;
	private double score;
	private int numMarkers;
	private int source;

	public CNVariant(String[] plinkLine) {
		this(plinkLine, 0);
	}

	public CNVariant(String[] plinkLine, int sourceIndex) {
		familyID = plinkLine[0];
		individualID = plinkLine[1];
		chr = Positions.chromosomeNumber(plinkLine[2]);
		start = Integer.parseInt(plinkLine[3]);
		stop = Integer.parseInt(plinkLine[4]);
		cn = Integer.parseInt(plinkLine[5]);
		score = Double.parseDouble(plinkLine[6]);
		numMarkers = Integer.parseInt(plinkLine[7]);
		source = sourceIndex;
	}

	public CNVariant(String familyID, String individualID, byte chr, int start, int stop, int cn, double score, int numMarkers, int source) {
		this.familyID = familyID;
		this.individualID = individualID;
		this.chr = chr;
		this.start = start;
		this.stop = stop;
		this.cn = cn;
		this.score = score;
		this.numMarkers = numMarkers;
		this.source = source;
	}

	public CNVariant(String ucscLocation) {
		int[] loc = Positions.parseUCSClocation(ucscLocation);
		chr = Positions.chromosomeNumber(loc[0] + "");
		start = loc[1];
		stop = loc[2];
	}

	public void setFamilyID(String famID) {
		familyID = famID;
	}

	public String getFamilyID() {
		return familyID;
	}

	public void setIndividualID(String indID) {
		individualID = indID;
	}

	public String getIndividualID() {
		return individualID;
	}

	public int getSource() { // generally not entered, since loading from a serialized file
		return source;
	}

	public void setSource(int sourceIndex) {
		source = sourceIndex;
	}

	public int getCN() {
		return cn;
	}

	public double getScore() {
		return score;
	}

	public int getNumMarkers() {
		return numMarkers;
	}

	public boolean equalsIncludingIndividual(CNVariant cnv) {
		return cnv.familyID.equals(familyID) && cnv.individualID.equals(individualID) && cnv.chr == chr && cnv.start == start && cnv.stop == stop && cnv.cn == cn;
	}

	public boolean overlapsLocAndIndividual(CNVariant cnv) {
		return familyID.equals(cnv.familyID) && individualID.equals(cnv.individualID) && amountOfOverlapInBasepairs(cnv) > 0;
	}

	public String toPlinkFormat() {
		return familyID + "\t" + individualID + "\t" + chr + "\t" + start + "\t" + stop + "\t" + cn + "\t" + ext.formDeci(score, 2) + "\t" + numMarkers;
	}

	public String getFingerprint() {
		return familyID + "_" + individualID + "_" + chr + "_" + start + "_" + stop + "_" + cn + "_" + numMarkers;
	}

	public static CNVariant[] toArray(Vector<CNVariant> v) {
		CNVariant[] cnvs;

		cnvs = new CNVariant[v == null ? 0 : v.size()];
		for (int i = 0; i < cnvs.length; i++) {
			cnvs[i] = v.elementAt(i);
		}
		return cnvs;
	}

	public static CNVariant[] sortCNVs(CNVariant[] array) {
		CNVariant[] newArray;
		int[] keys;

		keys = quicksort(array);
		newArray = new CNVariant[array.length];
		for (int i = 0; i < keys.length; i++) {
			newArray[i] = array[keys[i]];
		}

		return newArray;
	}

	public static CNVariant[] loadPlinkFile(String filename, boolean jar) {
		return CNVariant.sortCNVs(CNVariant.toArray(loadPlinkFile(filename, null, jar)));
	}

	public static Vector<CNVariant> loadPlinkFile(String filename, Hashtable<String, String> sampleHash, boolean jar) {
		BufferedReader reader;
		Vector<CNVariant> v = null;
		String[] line;

		v = new Vector<CNVariant>();
		try {
			reader = Files.getReader(filename, jar, true, true);

			reader.mark(1000);
			line = reader.readLine().trim().split("[\\s]+");
			if (!line[2].toLowerCase().equals("chr") && Positions.chromosomeNumber(line[2]) != -1) {
				reader.reset();
			}
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				if (sampleHash == null || sampleHash.containsKey(line[0] + "\t" + line[1])) {
					v.add(new CNVariant(line));
				}
			}
			reader.close();

			return v;
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
			fnfe.printStackTrace();
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			ioe.printStackTrace();
		}

		return null;
	}

	// TODO what was this created for? Incomplete
	public static void mergeCNVs(String filename) {
		BufferedReader reader;
		String[] line;
		// PrintWriter writer;
		// Hashtable<String, Vector<CNVariant>> hash;
		// hash = new Hashtable<String, Vector<CNVariant>>();
		StringVector markerNames;
		ByteVector chrs;
		IntVector positions;

		// hash.add(fID+"\t"+iId, Vector<CNVariant>);
		try {
			reader = new BufferedReader(new FileReader(filename));

			markerNames = new StringVector();
			chrs = new ByteVector();
			positions = new IntVector();

			reader.readLine();
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				markerNames.add(line[0]);
				chrs.add((byte) Positions.chromosomeNumber(line[1]));
				positions.add(Integer.parseInt(line[2]));
			}
			reader.close();
			/*
			 * new MarkerSet(markerNames, chrs, positions).serialize(filename+".ser");
			 * 
			 * writer = new PrintWriter(new FileWriter("????")); for (int i = 0; i<temp.length; i++) { writer.println(temp[i].toPlinkFormat()); } writer.close();
			 */

		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			System.exit(2);
		}
	}

	public static void findConsensus(String file1, String file2) {
		PrintWriter writer;
		CNVariant[] list1, list2, consensus;
		Vector<CNVariant> v = new Vector<CNVariant>();

		list1 = loadPlinkFile(file1, false);
		list2 = loadPlinkFile(file2, false);

		for (int i = 0; i < list1.length; i++) {
			for (int j = 0; j < list2.length; j++) {
				if (list1[i].overlapsLocAndIndividual(list2[j])) {
					v.add(list1[i]);
				}
			}
		}
		consensus = sortCNVs(CNVariant.toArray(v));

		try {
			writer = new PrintWriter(new FileWriter(ext.rootOf(file1) + "_" + ext.rootOf(file2) + "_consensus.cnv"));
			for (int i = 0; i < consensus.length; i++) {
				writer.println(consensus[i].toPlinkFormat());
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing consensus");
			e.printStackTrace();
		}
	}

	public static CNVariant[] putInOrder(CNVariant[] array, int[] order) {
		CNVariant[] newArray;

		newArray = new CNVariant[array.length];
		for (int i = 0; i < order.length; i++) {
			newArray[i] = array[order[i]];
		}

		return newArray;
	}

	public static CNVariant[] sort(CNVariant[] array) {
		return putInOrder(array, quicksort(array));
	}

	public CNVariant clone() {
		return new CNVariant(familyID, individualID, chr, start, stop, cn, score, numMarkers, source);
	}

	public static void main(String[] args) {
		String file1 = "C:\\Documents and Settings\\npankrat\\My Documents\\CNV\\penncnv\\again_noGenderProblems\\conf_0kb_5SNP_10.0.cnv";
		String file2 = "C:\\Documents and Settings\\npankrat\\My Documents\\CNV\\quantisnp\\noGenderProblems\\conf_0kb_5SNP_10.0.cnv";
		// String file1 = "Penn_conf_0kb_5SNP_10.0.cnv";
		// String file2 = "Quanti_conf_0kb_5SNP_10.0.cnv";

		try {
			findConsensus(file1, file2);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	// Return the Individual ID when getting as a string
	public String toString() {
		return getIndividualID();
	}
}
