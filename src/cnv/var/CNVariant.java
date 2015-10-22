package cnv.var;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Vector;

import common.ByteVector;
import common.Files;
import common.IntVector;
import common.Logger;
import common.Positions;
import common.Sort;
import common.StringVector;
import common.ext;
import filesys.Segment;

public class CNVariant extends Segment {
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = super.hashCode();
		result = prime * result + cn;
		result = prime * result + ((familyID == null) ? 0 : familyID.hashCode());
		result = prime * result + ((individualID == null) ? 0 : individualID.hashCode());
		result = prime * result + numMarkers;
		long temp;
		temp = Double.doubleToLongBits(score);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		result = prime * result + source;
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (!super.equals(obj))
			return false;
		if (getClass() != obj.getClass())
			return false;
		CNVariant other = (CNVariant) obj;
		if (cn != other.cn)
			return false;
		if (familyID == null) {
			if (other.familyID != null)
				return false;
		} else if (!familyID.equals(other.familyID))
			return false;
		if (individualID == null) {
			if (other.individualID != null)
				return false;
		} else if (!individualID.equals(other.individualID))
			return false;
		if (numMarkers != other.numMarkers)
			return false;
		if (Double.doubleToLongBits(score) != Double.doubleToLongBits(other.score))
			return false;
		if (source != other.source)
			return false;
		return true;
	}

	public static final long serialVersionUID = 1L;
	public static final String[] PLINK_CNV_HEADER = { "FID", "IID", "CHR", "BP1", "BP2", "TYPE", "SCORE", "SITES" };

	private String familyID;
	private String individualID;
	private int cn;
	private double score;
	protected int numMarkers;
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
	
	public boolean overlapsLocAndIndividualSignificantly(CNVariant cnv) {
		boolean overlapsLocAndIndividualSignificantly = familyID.equals(cnv.familyID) && individualID.equals(cnv.individualID) && significantOverlap(cnv);	
		return overlapsLocAndIndividualSignificantly;
	}

	public String toPlinkFormat() {
		return familyID + "\t" + individualID + "\t" + chr + "\t" + start + "\t" + stop + "\t" + cn + "\t" + ext.formDeci(score, 2) + "\t" + numMarkers;
	}

	public String getFingerprint() {
		return familyID + "_" + individualID + "_" + chr + "_" + start + "_" + stop + "_" + cn + "_" + numMarkers;
	}

	public static CNVariant[] toCNVariantArray(Vector<CNVariant> v) {
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
		return CNVariant.sortCNVs(CNVariant.toCNVariantArray(loadPlinkFile(filename, null, jar)));
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
		consensus = sortCNVs(CNVariant.toCNVariantArray(v));

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
	
	
	public enum CONSENSUS_TYPE{
		/**
		 * Check copy number state, 
		 */
		CN_AWARE,/**
		 * Do not check copy number state
		 */
		NOT_CN_AWARE;
	}
	
	public enum OVERLAP_TYPE {
		/**
		 * Must have significant overlap
		 */
		OVERLAP_LOC_AND_INDIVIDUAL_SIGNIFICANTLY, /**
		 * Any overlap
		 */
		OVERLAP_LOC_AND_INDIVIDUAL;
	}

	public static MatchResults findSignificantConsensus(String file1, String file2, CONSENSUS_TYPE cType,OVERLAP_TYPE oType) {
		String output = ext.rootOf(file1) + "_" + ext.rootOf(file2) + "_signif_consensus.cnv";
		return findSignificantConsensus(file1, file2, output, cType, oType);
	}

	public static MatchResults findSignificantConsensus(String file1, String file2,String output, CONSENSUS_TYPE cType,OVERLAP_TYPE oType) {
		PrintWriter writer;
		CNVariant[] list1, list2;
		HashSet<CNVariant> matched1 = new HashSet<CNVariant>();
		HashSet<CNVariant> matched2 = new HashSet<CNVariant>();

		HashSet<CNVariant> unmatched1 = new HashSet<CNVariant>();
		HashSet<CNVariant> unmatched2 = new HashSet<CNVariant>();
		Vector<String> outputLines = new Vector<String>();
		
		list1 = loadPlinkFile(file1, false);
		list2 = loadPlinkFile(file2, false);

		for (int i = 0; i < list1.length; i++) {
			for (int j = 0; j < list2.length; j++) {
				boolean match = false;
				switch(oType){
				case OVERLAP_LOC_AND_INDIVIDUAL:
					match = list1[i].overlapsLocAndIndividual(list2[j]);
					break;
				case OVERLAP_LOC_AND_INDIVIDUAL_SIGNIFICANTLY:
					match = list1[i].overlapsLocAndIndividualSignificantly(list2[j]);
					break;
				default:
					break;
				}
				if (match) {
					switch (cType) {
					case CN_AWARE:
						match = list1[i].getCN() == list2[j].getCN();
						break;
					case NOT_CN_AWARE:
						break;
					default:
						break;

					}
				}
				if (match) {
					matched1.add(list1[i]);
					matched2.add(list2[j]);
				}
			}
		}
		
		for (int i = 0; i < list1.length; i++) {
			if (matched1.contains(list1[i])) {
				outputLines.add(list1[i].toPlinkFormat() + "\t3");
			} else {
				unmatched1.add(list1[i]);
			}
		}

		for (int i = 0; i < list2.length; i++) {
			if (matched2.contains(list2[i])) {
				outputLines.add(list2[i].toPlinkFormat() + "\t3");
			} else {
				unmatched2.add(list2[i]);
			}
		}
		
		Iterator<CNVariant> unmatchedIterator = unmatched1.iterator();
		while(unmatchedIterator.hasNext()) {
			outputLines.add(unmatchedIterator.next().toPlinkFormat() + "\t1");
		}
		unmatchedIterator = unmatched2.iterator();
		while(unmatchedIterator.hasNext()) {
			outputLines.add(unmatchedIterator.next().toPlinkFormat() + "\t2");
		}

		try {
			writer = new PrintWriter(new FileWriter(output));
			for (int i = 0; i < outputLines.size(); i++) {
				writer.println(outputLines.get(i));
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing significant-consensus file");
			e.printStackTrace();
		}
		return new MatchResults(file1, file2, matched1, matched2, unmatched1, unmatched2);
		
	}
	
	public static class MatchResults implements Serializable {
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;
		private String file1;
		private String file2;
		private HashSet<CNVariant> matched1;
		private HashSet<CNVariant> matched2;

		private HashSet<CNVariant> unmatched1;
		private HashSet<CNVariant> unmatched2;

		public MatchResults(String file1, String file2, HashSet<CNVariant> matched1, HashSet<CNVariant> matched2, HashSet<CNVariant> unmatched1, HashSet<CNVariant> unmatched2) {
			super();
			this.file1 = file1;
			this.file2 = file2;
			this.matched1 = matched1;
			this.matched2 = matched2;
			this.unmatched1 = unmatched1;
			this.unmatched2 = unmatched2;
		}

		public String getFile1() {
			return file1;
		}

		public String getFile2() {
			return file2;
		}

		public HashSet<CNVariant> getMatched1() {
			return matched1;
		}

		public HashSet<CNVariant> getMatched2() {
			return matched2;
		}

		public HashSet<CNVariant> getUnmatched1() {
			return unmatched1;
		}

		public HashSet<CNVariant> getUnmatched2() {
			return unmatched2;
		}
		
		public void writeSerial(String fileName){
			Files.writeSerial(this, fileName, true);
		}

		public static MatchResults readSerial(String fileName, Logger log) {
			return (MatchResults) Files.readSerial(fileName, false, log, false, true);
		}
		
		
	}
	
	public static void findSignificantConsensus(String file1, String file2, boolean checkLarger) {
		PrintWriter writer;
		CNVariant[] list1, list2;
		HashSet<CNVariant> matched1 = new HashSet<CNVariant>();
		HashSet<CNVariant> matched2 = new HashSet<CNVariant>();

		HashSet<CNVariant> unmatched1 = new HashSet<CNVariant>();
		HashSet<CNVariant> unmatched2 = new HashSet<CNVariant>();
		Vector<String> outputLines = new Vector<String>();
		
		list1 = loadPlinkFile(file1, false);
		list2 = loadPlinkFile(file2, false);
		
		for (int i = 0; i < list1.length; i++) {
			for (int j = 0; j < list2.length; j++) {

				if (list1[i].familyID.equals(list2[j].familyID) && list1[i].equals(list2[j].individualID) && list1[i].significantOverlap(list2[j], checkLarger)) {
					matched1.add(list1[i]);
					matched2.add(list2[j]);
				}
			}
		}
		
		for (int i = 0; i < list1.length; i++) {
			if (matched1.contains(list1[i])) {
				outputLines.add(list1[i].toPlinkFormat() + "\t3");
			} else {
				unmatched1.add(list1[i]);
			}
		}

		for (int i = 0; i < list2.length; i++) {
			if (matched2.contains(list2[i])) {
				outputLines.add(list2[i].toPlinkFormat() + "\t3");
			} else {
				unmatched2.add(list1[i]);
			}
		}
		
		Iterator<CNVariant> unmatchedIterator = unmatched1.iterator();
		while(unmatchedIterator.hasNext()) {
			outputLines.add(unmatchedIterator.next().toPlinkFormat() + "\t1");
		}
		unmatchedIterator = unmatched2.iterator();
		while(unmatchedIterator.hasNext()) {
			outputLines.add(unmatchedIterator.next().toPlinkFormat() + "\t2");
		}
		

		try {
			writer = new PrintWriter(new FileWriter(ext.rootOf(file1) + "_" + ext.rootOf(file2) + "_signif_consensus.cnv"));
			for (int i = 0; i < outputLines.size(); i++) {
				writer.println(outputLines.get(i));
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing significant-consensus file");
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

	/**
	 * @param array
	 *            sort by the quality scores
	 * @return
	 */
	public static <T extends CNVariant> ArrayList<T> sortByQuality(final T[] array, int direction) {
		double[] scores = new double[array.length];
		for (int i = 0; i < array.length; i++) {
			scores[i] = array[i].getScore();
		}
		int[] order = Sort.quicksort(scores, direction);
		ArrayList<T> sorted = new ArrayList<T>(array.length);
		for (int i = 0; i < order.length; i++) {
			sorted.add(array[order[i]]);
		}
		return sorted;
	}
	
	public CNVariant clone() {
		return new CNVariant(familyID, individualID, chr, start, stop, cn, score, numMarkers, source);
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename1 = "CNVariant1.dat";
		String filename2 = "CNVariant2.dat";
		String logfile = null;
		boolean signif = false;
		
		String usage = 
						"\n" + 
						"cnv.var.CNVariant requires 2+ arguments\n" + 
						"   (1) first CNV filename (i.e. file1=" + filename1 + " (not the default))\n" +
						"   (2) second CNV filename (i.e. file2=" + filename2 + " (not the default))\n" +
						"   (3) (Optional) restrict to significant overlap (i.e. -sig (not the default))\n" + 
						"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file1=")) {
				filename1 = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("file2=")) {
				filename1 = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("-sig")) {
				signif = true;
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			if (signif) {
				findSignificantConsensus(filename1, filename2, CONSENSUS_TYPE.NOT_CN_AWARE,OVERLAP_TYPE.OVERLAP_LOC_AND_INDIVIDUAL_SIGNIFICANTLY);
			} else {
				findConsensus(filename1, filename2);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
//	public static void main(String[] args) {
//		String file1 = "C:\\Documents and Settings\\npankrat\\My Documents\\CNV\\penncnv\\again_noGenderProblems\\conf_0kb_5SNP_10.0.cnv";
//		String file2 = "C:\\Documents and Settings\\npankrat\\My Documents\\CNV\\quantisnp\\noGenderProblems\\conf_0kb_5SNP_10.0.cnv";
//		// String file1 = "Penn_conf_0kb_5SNP_10.0.cnv";
//		// String file2 = "Quanti_conf_0kb_5SNP_10.0.cnv";
//
//		try {
//			findConsensus(file1, file2);
//		} catch (Exception e) {
//			e.printStackTrace();
//		}
//	}

	// Return the Individual ID when getting as a string
	public String toString() {
		return getIndividualID();
	}

	@Override
	public String toAnalysisString() {
		return toPlinkFormat();
	}

	@Override
	public String[] getHeader() {
		return PLINK_CNV_HEADER;
	}

	public static LocusSet<CNVariant> loadLocSet(String cnvFile, Logger log) {
		CNVariant[] cnvs = CNVariant.loadPlinkFile(cnvFile, false);
		LocusSet<CNVariant> cLocusSet = new LocusSet<CNVariant>(cnvs, true, log) {

			/**
		 * 
		 */
			private static final long serialVersionUID = 1L;

		};
		return cLocusSet;

	}

	public static HashSet<String> getUniqueInds(String cnvFile, Logger log) {
		LocusSet<CNVariant> cLocusSet = CNVariant.loadLocSet(cnvFile, log);

		return getUniqueInds(cLocusSet, log);
	}

	public static HashSet<String> getUniqueInds(LocusSet<CNVariant> cLocusSet, Logger log) {
		HashSet<String> uniqueInds = new HashSet<String>();
		for (int i = 0; i < cLocusSet.getLoci().length; i++) {
			CNVariant tmp = cLocusSet.getLoci()[i];
			uniqueInds.add(tmp.getFamilyID() + "\t" + tmp.getIndividualID());
		}
		return uniqueInds;
	}

	/**
	 * Useful when the terms of the cnv are not completely known, such as in CNVCalling
	 *
	 */
	public static class CNVBuilder {
		private byte chr = -1;
		private int start = -1;
		private int stop = -1;
		private String familyID = null;
		private String individualID = null;
		private int cn = -22;
		private double score = -22;
		private int numMarkers = -1;
		private int source = -1;
		private int startIndex = -1;// for tracking external indices, such as positions in project
		private int stopIndex = -1;// for tracking external indices, such as positions in project

		public CNVBuilder() {

		}
		public CNVBuilder(CNVariant cnVariant) {
			this.chr = cnVariant.chr;
			this.start = cnVariant.start;
			this.stop = cnVariant.stop;
			this.familyID = cnVariant.familyID;
			this.individualID = cnVariant.individualID;
			this.cn = cnVariant.cn;
			this.score = cnVariant.score;
			this.numMarkers = cnVariant.numMarkers;
			this.source = cnVariant.source;
		}
		
		public CNVBuilder familyID(String familyID) {
			this.familyID = familyID;
			return this;
		}

		public CNVBuilder individualID(String individualID) {
			this.individualID = individualID;
			return this;
		}

		public CNVBuilder start(int start) {
			this.start = start;
			return this;
		}

		public CNVBuilder stop(int stop) {
			this.stop = stop;
			return this;
		}

		public CNVBuilder startIndex(int startIndex) {
			this.startIndex = startIndex;
			return this;
		}

		public int getStartIndex() {
			return startIndex;
		}
		public CNVBuilder stopIndex(int stopIndex) {
			this.stopIndex = stopIndex;
			return this;
		}

		public CNVBuilder chr(byte chr) {
			this.chr = chr;
			return this;
		}

		public CNVBuilder cn(int cn) {
			this.cn = cn;
			return this;
		}

		public CNVBuilder score(double score) {
			this.score = score;
			return this;
		}

		public CNVBuilder numMarkers(int numMarkers) {
			this.numMarkers = numMarkers;
			return this;
		}

		public CNVBuilder source(int source) {
			this.source = source;
			return this;
		}

		public int getStart() {
			return start;
		}

		public int getStop() {
			return stop;
		}

		public int getCn() {
			return cn;
		}

		public int getNumMarkers() {
			return numMarkers;
		}

		public CNVariant build() {
			return new CNVariant(this);
		}
	}

	private CNVariant(CNVBuilder builder) {
		super(builder.chr, builder.start, builder.stop);
		this.familyID = builder.familyID;
		this.individualID = builder.individualID;
		this.cn = builder.cn;
		this.score = builder.score;
		this.numMarkers = builder.numMarkers;
		this.source = builder.source;
	}
	

}

