package seq.qc;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Vector;

import common.Array;
import common.Files;
import common.HashVec;
import common.Logger;
import common.ext;
import filesys.GeneData;
import filesys.GeneTrack;
import filesys.LocusSet;
import filesys.Segment;

/**
 * @author lane0212 Maybe a one-off thing, but this is for summarizing bamQC results across a gene
 */

// coveredMRNA,coveredMRNA-UTRs,numBaits,avgCoverage
public class GeneQC {
	private String bamQCFile;
	private LocusSet<Segment> utrs;
	private LocusSet<GeneData> genes;
	private GeneSummary[] geneSummaries;

	private Logger log;

	public GeneQC(LocusSet<GeneData> genes, String bamQCFile, LocusSet<Segment> utrs, Logger log) {
		super();
		this.genes = genes;
		this.geneSummaries = new GeneSummary[genes.getLoci().length];
		this.bamQCFile = bamQCFile;
		this.utrs = utrs;
		this.log = log;
		initializeSummaries();
	}

	private void initializeSummaries() {
		for (int i = 0; i < genes.getLoci().length; i++) {
			//int totalMrna = 0;
			//int mrnaNoUtrs = 0;
			for (int j = 0; j < genes.getLoci()[i].getExonBoundaries().length; j++) {
				Segment exon = new Segment(genes.getLoci()[i].getChr(), genes.getLoci()[i].getExonBoundaries()[j][0], genes.getLoci()[i].getExonBoundaries()[j][1]);
				//totalMrna += exon.getSize();
				//mrnaNoUtrs += exon.getSize();
				Segment[] utrsOlap = utrs.getOverLappingLoci(exon);
				if (utrsOlap != null) {
					Vector<Segment> mergedUtrs = Segment.toVector(utrsOlap);
					Segment.mergeOverlapsAndSort(mergedUtrs);
					for (int k = 0; k < mergedUtrs.size(); k++) {
					//	mrnaNoUtrs -= exon.getIntersection(mergedUtrs.get(k), log).getSize();
					}
				}
			}
		//	geneSummaries[i] = new GeneSummary(genes.getLoci()[i].getGeneName(), genes.getLoci()[i].getExonBoundaries().length, genes.getLoci()[i].getNcbiAssessionNumbers().length, totalMrna, mrnaNoUtrs);
		}
	}

	public void qcByGene() {
		String output = bamQCFile + ".quickSummary";
		String[] targets = new String[] { "GENE", "NumOtherGenes", "NumExons", "averageCoverage", "averageGC", "numBaitsPerTarget", "MRNA_Overlap", "NonUTRMrnaOverlap" };

		try {
			BufferedReader reader = Files.getAppropriateReader(bamQCFile);
			PrintWriter writer = new PrintWriter(new FileWriter(output));
			writer.println("GENE\tNumExons\tNumOtherGenes\tMRNA_Overlap\tNonUTRMrnaOverlap\t" + Array.toStr(Files.getHeaderOfFile(bamQCFile, log)));

			reader.readLine();
			while (reader.ready()) {
				String[] line = reader.readLine().trim().split("\t");
				Segment seqSeg = new Segment(line[0]);
				int[] gIndices = genes.getOverlappingIndices(seqSeg);
				if (gIndices != null) {
					for (int i = 0; i < gIndices.length; i++) {
						String gene = genes.getLoci()[gIndices[i]].getGeneName();
						int numExons = genes.getLoci()[gIndices[i]].getExonBoundaries().length;
						int numOtherGenes = gIndices.length - 1;
						int mrna = 0;
						int mrnaNonUtr = 0;
						for (int j = 0; j < genes.getLoci()[gIndices[i]].getExonBoundaries().length; j++) {
							Segment curExon = new Segment(genes.getLoci()[gIndices[i]].getChr(), genes.getLoci()[gIndices[i]].getExonBoundaries()[j][0], genes.getLoci()[gIndices[i]].getExonBoundaries()[j][1]);
							if (seqSeg.overlaps(curExon)) {
								mrna += seqSeg.getIntersection(curExon, log).getSize();
								mrnaNonUtr += seqSeg.getIntersection(curExon, log).getSize();
								Segment[] utrsOlap = utrs.getOverLappingLoci(seqSeg.getIntersection(curExon, log));
								if (utrsOlap != null) {
									Vector<Segment> mergedUtrs = Segment.toVector(utrsOlap);
									Segment.mergeOverlapsAndSort(mergedUtrs);
									for (int k = 0; k < mergedUtrs.size(); k++) {
										mrnaNonUtr -= seqSeg.getIntersection(mergedUtrs.get(k), log).getSize();
									}
								}
							}
						}
						if (mrna > 0) {
							writer.println(gene + "\t" + numExons + "\t" + numOtherGenes + "\t" + mrna + "\t" + Math.max(0, mrnaNonUtr) + "\t" + Array.toStr(line));
						}
					}
				}

			}

			reader.close();
			writer.close();

		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + bamQCFile + "\" not found in current directory");
			return;
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + bamQCFile + "\"");
			return;
		}

		String[][] toSumm = HashVec.loadFileToStringMatrix(output, true, ext.indexFactors(targets, Files.getHeaderOfFile(output, log), true, true), false);
		ArrayList<Integer> numMrnaTotal = new ArrayList<Integer>();
		ArrayList<Integer> numMrnaNonUTR = new ArrayList<Integer>();
		Hashtable<String, Integer> index = new Hashtable<String, Integer>();
		// String[] targets = new String[] { "GENE","NumOtherGenes","NumExons", "averageCoverage", "averageGC","numBaitsPerTarget" ,"MRNA_Overlap","NonUTRMrnaOverlap"};

		int curIndex = 0;
		for (int i = 0; i < toSumm.length; i++) {
			if (!index.containsKey(toSumm[i][0])) {
				index.put(toSumm[i][0], curIndex);
				numMrnaTotal.add(Integer.parseInt(toSumm[i][6]));
				numMrnaNonUTR.add(Integer.parseInt(toSumm[i][7]));
				curIndex++;
			} else {
			//	int in = index.get(toSumm[i][0]);
				// numMrnaTotal.get(in)
			}

			System.out.println(Array.toStr(toSumm[i]));
		}

	}

	public GeneSummary[] getGeneSummaries() {
		return geneSummaries;
	}

	private static class GeneSummary {
		//private static final String[] SUMMARY = new String[] { "Gene", "NumExons", "numIsoForms", "TotalMrna", "TotalMrnaMinusUTR", "CoveredMrna", "CoveredMrnaNoUtrs" };
//		private String geneName;
//		private int numExons;
//		private int numIsoForms;
//		private int totalMrna;
//		private int mrnaNoUtrs;
//		private int coveredMrna;
//		private int coveredMrnaNoUTRS;
//		private ArrayList<Segment> utrsAdded;
//
//		public GeneSummary(String geneName, int numExons, int numIsoForms, int totalMrna, int mrnaNoUtrs) {
//			super();
//			this.geneName = geneName;
//			this.numExons = numExons;
//			this.numIsoForms = numIsoForms;
//			this.totalMrna = totalMrna;
//			this.mrnaNoUtrs = mrnaNoUtrs;
//			this.coveredMrna = 0;
//			this.coveredMrnaNoUTRS = 0;
//			this.utrsAdded = new ArrayList<Segment>();
//		}

		

	}

	//
	// public GeneQCResult[] qcByGene() {
	// try {
	// BufferedReader reader = Files.getAppropriateReader(bamQCFile);
	// reader.readLine();
	// String[] targets = new String[] { "UCSC", "averageCoverage", "averageGC" };
	// int[] indices = ext.indexFactors(targets, Files.getHeaderOfFile(bamQCFile, log), true, false);
	// Hashtable<String, GeneQCResult> results = new Hashtable<String, GeneQC.GeneQCResult>();
	// for (int i = 0; i < geneTrack.getGenes().length; i++) {
	// for (int j = 0; j < geneTrack.getGenes()[i].length; j++) {
	// GeneQCResult geneQCResult = new GeneQCResult(geneTrack.getGenes()[i][j].getGeneName(), geneTrack.getGenes()[i][j], geneTrack.getGenes()[i][j].getExonBoundaries().length);
	// int mrnaLen = 0;
	// for (int k = 0; k < geneTrack.getGenes()[i][j].getExonBoundaries().length; k++) {
	//
	// Segment exon = new Segment(geneTrack.getGenes()[i][j].getChr(), geneTrack.getGenes()[i][j].getExonBoundaries()[k][0], geneTrack.getGenes()[i][j].getExonBoundaries()[k][1]);
	// Segment[] utrOlaps = utrs.getOverLappingLoci(exon);
	// ArrayList<Segment> partialExon = new ArrayList<Segment>();
	// if (utrOlaps != null) {
	// for (int l = 0; l < utrOlaps.length; l++) {
	// System.out.println(j + "\t" + exon.getUCSClocation());
	// Segment[] removed = exon.remove(utrOlaps[l], log);
	// if (removed != null) {
	// for (int m = 0; m < removed.length; m++) {
	// partialExon.add(removed[m]);
	// }
	// }
	// }
	// }
	// System.out.println(exon.getUCSClocation());
	// for (int l = 0; l < partialExon.size(); l++) {
	// exon = exon.getUnion(partialExon.get(l), log);
	// System.out.println("Part" + l + "\t" + exon.getUCSClocation());
	// }
	// try {
	// Thread.sleep(1000);
	// } catch (InterruptedException ie) {
	// }
	// geneQCResult.getExonQCResults()[k] = new ExonQCResult(exon, exon.getSize(), k + 1);
	//
	// }
	// }
	// }
	// while (reader.ready()) {
	// String[] line = reader.readLine().trim().split("\t");
	// Segment seqSeg = new Segment(line[0]);
	// GeneData[] geneDatas = geneTrack.getOverlappingGenes(seqSeg);
	// if (geneDatas != null) {
	// for (int i = 0; i < geneDatas.length; i++) {
	//
	// for (int j = 0; j < geneDatas[i].getExonBoundaries().length; j++) {
	// Segment exon =results.get(geneDatas[i].getGeneName()).getExonQCResults()[i].getExonSeg();
	// if(exon.overlaps(seqSeg)){
	//
	// }
	//
	// }
	// }
	// }
	// }
	// reader.close();
	// return null;
	// } catch (FileNotFoundException fnfe) {
	// log.reportError("Error: file \"" + bamQCFile + "\" not found in current directory");
	// return null;
	// } catch (IOException ioe) {
	// log.reportError("Error reading file \"" + bamQCFile + "\"");
	// return null;
	// }
	//
	// }

	public static void test(String bamQCFile, String geneTrackFile, String utr5p, String utr3p, Logger log) {
		LocusSet<Segment> utr5pSegs = LocusSet.loadSegmentSetFromFile(utr5p, 0, 1, 2, 0, true, true, 0, log);
		log.reportTimeInfo("Loaded " + utr5pSegs.getLoci().length + " 5' UTRs");

		LocusSet<Segment> utr3pSegs = LocusSet.loadSegmentSetFromFile(utr3p, 0, 1, 2, 0, true, true, 0, log);
		log.reportTimeInfo("Loaded " + utr3pSegs.getLoci().length + " 3' UTRs");
		LocusSet<Segment> utrs = new LocusSet<Segment>(Array.concatAll(utr5pSegs.getLoci(), utr3pSegs.getLoci()), true, log) {

			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

		};
		ArrayList<GeneData> genes = new ArrayList<GeneData>();

		GeneTrack geneTrack = GeneTrack.load(geneTrackFile, false);

		for (int i = 0; i < geneTrack.getGenes().length; i++) {
			for (int j = 0; j < geneTrack.getGenes()[i].length; j++) {
				genes.add(geneTrack.getGenes()[i][j]);
			}
		}
		LocusSet<GeneData> geneSet = new LocusSet<GeneData>(genes.toArray(new GeneData[genes.size()]), true, log) {

			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

		};

		GeneQC geneQC = new GeneQC(geneSet, bamQCFile, utrs, new Logger(bamQCFile + ".gqc.log"));
		geneQC.qcByGene();
		// String output = bamQCFile + ".gqc.summary";
		// try {
		// PrintWriter writer = new PrintWriter(new FileWriter(output));
		// writer.println(Array.toStr(GeneSummary.SUMMARY));
		// for (int i = 0; i < geneQC.getGeneSummaries().length; i++) {
		// writer.println(geneQC.getGeneSummaries()[i].getSummary());
		// }
		// writer.close();
		// } catch (Exception e) {
		// log.reportError("Error writing to " + output);
		// log.reportException(e);
		// }
	}

	public static void main(String[] args) {
		String bamQCFile = "D:/data/logan/OSv2_seq/geneQC/rrd_bamQC.targets.libraryResults.summary.txt";
		String geneTrackFile = "N:/statgen/NCBI/RefSeq.gtrack";
		String utr5p = "N:/statgen/NCBI/refGene_hg19_5pUTR.txt";
		String utr3p = "N:/statgen/NCBI/refGene_hg19_3pUTR.txt";
		Logger log = new Logger(ext.rootOf(bamQCFile, false) + ".log");
		test(bamQCFile, geneTrackFile, utr5p, utr3p, log);
	}

	// if (gIndices != null) {
	// for (int i = 0; i < gIndices.length; i++) {
	// for (int j = 0; j < genes.getLoci()[gIndices[i]].getExonBoundaries().length; j++) {
	// int curCovTotal = geneSummaries[gIndices[i]].getCoveredMrna();
	// int curCovTotalNoUtr = geneSummaries[gIndices[i]].getCoveredMrnaNoUTRS();
	// Segment curExon = new Segment(genes.getLoci()[gIndices[i]].getChr(), genes.getLoci()[gIndices[i]].getExonBoundaries()[j][0], genes.getLoci()[gIndices[i]].getExonBoundaries()[j][1]);
	// if (seqSeg.overlaps(curExon)) {
	// curCovTotal += seqSeg.getUnion(curExon, log).getSize();
	// curCovTotalNoUtr += seqSeg.getUnion(curExon, log).getSize();
	// Segment[] utrsOlap = utrs.getOverLappingLoci(seqSeg.getUnion(curExon, log));
	// if (utrsOlap != null) {
	// Vector<Segment> mergedUtrs = Segment.toVector(utrsOlap);
	// Segment.mergeOverlapsAndSort(mergedUtrs);
	// for (int k = 0; k < mergedUtrs.size(); k++) {
	// curCovTotalNoUtr -= seqSeg.getUnion(mergedUtrs.get(k), log).getSize();
	// }
	// }
	// geneSummaries[gIndices[i]].setCoveredMrna(curCovTotal);
	// geneSummaries[gIndices[i]].setCoveredMrnaNoUTRS(Math.max(curCovTotalNoUtr, 0));
	// }
	// }
	//
	// }
	// }
}
