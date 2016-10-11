package org.genvisis.filesys;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.Serializable;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.HashVec;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.ext;

public class GeneSet implements Serializable {
	public static final long serialVersionUID = 1L;
	// public static final String REFSEQ_FILE = "refGene.txt";
	// public static final String REFSEQ_FILE = "refGene_hg19.txt";
	public static final String REFSEQ_DB = "RefSeq.genes";
	public static final String REFSEQ_TRACK = "RefSeq.gtrack";
	public static final String REFSEQ_SEGS = "RefSeq_2k.segs";
	public static final String REFSEQ_EXONS = "RefSeq_exons_2k.segs";
	public static final String KNOWN_FILE = "knownGene.txt";
	public static final String KNOWN_LOOKUP_FILE = "kgXref.txt";

	private final GeneData[] set;
	private final long fingerprint;

	public GeneSet(GeneData[] set) {
		this.set = set;
		fingerprint = (long) (Math.random() * Long.MAX_VALUE);
	}

	public GeneSet(Vector<GeneData> setVec) {
		set = new GeneData[setVec.size()];
		for (int i = 0; i < setVec.size(); i++) {
			set[i] = setVec.elementAt(i);
		}

		fingerprint = (long) (Math.random() * Long.MAX_VALUE);
	}

	public long getFingerprint() {
		return fingerprint;
	}

	public GeneData[] getSet() {
		return set;
	}

	public void serialize(String filename) {
		SerializedFiles.writeSerial(this, filename);
	}

	public static GeneSet load(String filename, boolean jar) {
		return (GeneSet) SerializedFiles.readSerial(filename, jar, true);
	}

	public SegmentLists determineGeneSegments(int window) {
		Hashtable<String, Vector<Segment>> hash = new Hashtable<String, Vector<Segment>>();
		Vector<Segment> segs;
		Segment[][] lists;
		int[] chrs;

		for (GeneData element : set) {
			if (hash.containsKey(element.getChr() + "")) {
				segs = hash.get(element.getChr() + "");
			} else {
				hash.put(element.getChr() + "", segs = new Vector<Segment>());
			}
			segs.add(new Segment(	element.getChr(), element.getStart() - window,
														element.getStop() + window));
		}
		chrs = Array.toIntArray(HashVec.getKeys(hash));
		lists = new Segment[Array.max(chrs) + 1][];
		for (int i = 0; i < chrs.length; i++) {
			segs = hash.get(chrs[i] + "");
			Segment.mergeOverlapsAndSort(segs);
			lists[chrs[i]] = Segment.toArray(segs);
		}

		return new SegmentLists(lists);
	}

	public SegmentLists determineExonSegments(int upDownWindow) {
		Hashtable<String, Vector<Segment>> hash = new Hashtable<String, Vector<Segment>>();
		Vector<Segment> segs;
		Segment[][] lists;
		int[] chrs;
		int[][] exonBoundaries;
		GeneData gene;

		for (GeneData element : set) {
			gene = element;
			if (hash.containsKey(gene.getChr() + "")) {
				segs = hash.get(gene.getChr() + "");
			} else {
				hash.put(gene.getChr() + "", segs = new Vector<Segment>());
			}
			segs.add(new Segment(gene.getChr(), gene.getStart() - upDownWindow, gene.getStart()));
			segs.add(new Segment(gene.getChr(), gene.getStop(), gene.getStop() + upDownWindow));
			exonBoundaries = gene.getExonBoundaries();
			for (int[] exonBoundarie : exonBoundaries) {
				segs.add(new Segment(exonBoundarie[0], exonBoundarie[1]));
			}
		}
		chrs = Array.toIntArray(HashVec.getKeys(hash));
		lists = new Segment[Array.max(chrs) + 1][];
		for (int i = 0; i < chrs.length; i++) {
			segs = hash.get(chrs[i] + "");
			Segment.mergeOverlapsAndSort(segs);
			lists[chrs[i]] = Segment.toArray(segs);
		}

		return new SegmentLists(lists);
	}

	public static void parseRefSeqGenes(String sourceFile) {
		BufferedReader reader;
		Hashtable<String, Vector<GeneData>> hash = new Hashtable<String, Vector<GeneData>>();
		Vector<GeneData> v, overlapping, finalList;
		GeneData gene;
		String[] geneNames;
		Vector<String> assessionNumbers;
		byte strand;
		Vector<Segment> exons;
		boolean newlyAdded;
		byte count;
		int[][] exonBoundaries;
		GeneSet geneSet;
		String dir;

		dir = ext.parseDirectoryOfFile(sourceFile);

		System.out.println("Parsing gene info...");
		try {
			reader = new BufferedReader(new FileReader(sourceFile));
			while (reader.ready()) {
				gene = new GeneData(reader.readLine());
				if (gene.isFinalized()) {
					if (hash.containsKey(gene.getGeneName())) {
						v = hash.get(gene.getGeneName());
					} else {
						hash.put(gene.getGeneName(), v = new Vector<GeneData>());
					}
					v.add(gene);
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + sourceFile + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + sourceFile + "\"");
			System.exit(2);
		}

		System.out.println("Collapsing isoforms...");
		geneNames = HashVec.getKeys(hash);
		finalList = new Vector<GeneData>();
		for (int i = 0; i < geneNames.length; i++) {
			v = hash.get(geneNames[i]);
			count = 0;
			while (v.size() > 0) {
				overlapping = new Vector<GeneData>();
				overlapping.add(v.remove(0));
				strand = overlapping.elementAt(0).getStrand();
				newlyAdded = true;
				while (newlyAdded) {
					newlyAdded = false;
					for (int j = 0; j < v.size(); j++) {
						for (int k = 0; k < overlapping.size(); k++) {
							if (v.elementAt(j).overlaps(overlapping.elementAt(k))) {
								overlapping.add(v.remove(j));
								j--;
								k = overlapping.size();
								newlyAdded = true;
							}
						}
					}
				}
				count++; // genes can be in more than one "finished" location, and even be on different
									// chromosomes (though half the time it's an X/Y pairing)

				assessionNumbers = new Vector<String>();
				exons = new Vector<Segment>();
				for (int j = 0; j < overlapping.size(); j++) {
					gene = overlapping.elementAt(j);
					finalList.add(gene); // add all genes/isoforms to list
					assessionNumbers.add(gene.getNcbiAssessionNumbers()[0]);
					if (gene.getStrand() != strand) {
						strand = GeneData.BOTH_STRANDS; // genes can be listed as being on both strands
					}
					exonBoundaries = gene.getExonBoundaries();

					// different isoforms have different combinations of the same superset of exons
					for (int[] exonBoundarie : exonBoundaries) {
						Segment.addIfAbsent(new Segment(exonBoundarie[0], exonBoundarie[1]), exons);
					}
				}

				// different isoforms can splice at different spots in the same exon
				Segment.mergeOverlapsAndSort(exons);

				exonBoundaries = Segment.convertListToSortedBoundaries(exons);
				if (geneNames[i].equals("MIR4444-2")) { // an example of a multi-chr gene
					for (int j2 = 0; j2 < overlapping.size(); j2++) {
						System.out.println(i + "\t" + overlapping.get(j2).getUCSClocation() + "\t" + count);
					}
				}
				finalList.add(new GeneData(	geneNames[i], Array.toStringArray(assessionNumbers),
																		overlapping.elementAt(0).getChr(), true, strand,
																		exonBoundaries[0][0], exonBoundaries[exons.size() - 1][1],
																		exonBoundaries, (count == 1 && v.size() == 0 ? 0 : count),
																		true));
			}
		}

		geneSet = new GeneSet(finalList);
		geneSet.serialize(dir + REFSEQ_DB);
		System.out.println("Parsing gene segments...");
		geneSet.determineGeneSegments(2000).serialize(dir + REFSEQ_SEGS);
		System.out.println("Parsing exon segments...");
		geneSet.determineExonSegments(2000).serialize(dir + REFSEQ_EXONS);
		System.out.println("Parsing gene track...");
		new GeneTrack(dir + REFSEQ_DB).serialize(dir + REFSEQ_TRACK);
	}

	public static void parseUCSCknownGenes() {}

	// public static GeneSet loadRefSeqGenes() {
	// return load(DIRECTORY+REFSEQ_DB, false);
	// }

	public static void main(String[] args) {
		int numArgs = args.length;
		boolean refseq = true;
		boolean known = false;

		String usage = "\n"	+ "filesys.GeneSet requires 0-1 arguments\n"
										+ "   (1) parse RefSeq genes (i.e. -refseq (not the default))\n"
										+ "   (2) parse UCSC known genes (i.e. -known (not the default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("-refseq")) {
				refseq = true;
				numArgs--;
			} else if (arg.startsWith("-known")) {
				known = true;
				numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			if (refseq) {
				// parseRefSeqGenes(Aliases.getPathToFileInReferenceDirectory("gc5Base.txt", true, new
				// Logger()));
				parseRefSeqGenes("N:/statgen/NCBI/testAgain/refGene_hg19.txt");
			}
			if (known) {
				parseUCSCknownGenes();
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
