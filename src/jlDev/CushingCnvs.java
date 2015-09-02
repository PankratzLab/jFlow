package jlDev;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;

import seq.manage.BedOps;
import seq.qc.Mappability;
import cnv.var.CNVariant;
import cnv.var.LocusSet;
import common.Array;
import common.Logger;
import common.ext;

public class CushingCnvs {
	private static final String[] BASE_HEADER = new String[] { "MAPPABILITY_SCORE", "CALLED_GENE(s)" };

	public static void filter(String mappabilityFile, String cnvFile, String[] cnvRemoveFiles, String geneTrackFile, String callSubsetBed, Logger log) {
		BedOps.verifyBedIndex(mappabilityFile, log);
		// LocusSet<GeneData> gLocusSet = GeneTrack.load(geneTrackFile, false).convertToLocusSet(log);
		LocusSet<CNVariant> cLocusSet = CNVariant.loadLocSet(cnvFile, log);

		SummaryCNV[] summaryCNVs = new SummaryCNV[cLocusSet.getLoci().length];
		for (int i = 0; i < summaryCNVs.length; i++) {
			summaryCNVs[i] = new SummaryCNV(cnvRemoveFiles);
		}
		for (int i = 0; i < cnvRemoveFiles.length; i++) {
			LocusSet<CNVariant> cLocusRemoveSet = CNVariant.loadLocSet(cnvRemoveFiles[i], log);
			int numInds = CNVariant.getUniqueInds(cnvRemoveFiles[i], log).size();
			for (int j = 0; j < cLocusSet.getLoci().length; j++) {
				summaryCNVs[j].getTotalNumInds()[i] = numInds;
				CNVariant currentCNV = cLocusSet.getLoci()[j];
				CNVariant[] overlaps = cLocusRemoveSet.getOverLappingLoci(currentCNV);
				if (overlaps != null && overlaps.length > 0) {
					summaryCNVs[j].getNumOverLaps()[i] += overlaps.length;
					int numSig = 0;
					int numSigCN = 0;
					for (int k = 0; k < overlaps.length; k++) {
						if (currentCNV.significantOverlap(overlaps[k]) && overlaps[k].significantOverlap(currentCNV)) {
							numSig++;
							if (currentCNV.getCN() == overlaps[k].getCN()) {
								numSigCN++;
							}
						}
					}
					summaryCNVs[j].getNumSigOverLaps()[i] += numSig;
					summaryCNVs[j].getNumSigOverLapsSameCopyNumber()[i] += numSigCN;

				}
			}
		}
		Mappability<CNVariant> cnMappability = new Mappability<CNVariant>(cLocusSet, mappabilityFile, callSubsetBed, log);
		cnMappability.computeMappability();

		String output = ext.rootOf(cnvFile, false) + ".qc.summary.txt";
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(output));
			writer.println(Array.toStr(Array.concatAll(CNVariant.PLINK_CNV_HEADER, summaryCNVs[0].getHeader(), BASE_HEADER)));
			for (int i = 0; i < summaryCNVs.length; i++) {
				CNVariant currentCNV = cLocusSet.getLoci()[i];

				writer.print(currentCNV.toPlinkFormat());
				writer.print("\t" + Array.toStr(summaryCNVs[i].getData()));
				writer.print("\t" + cnMappability.getMappabilityResults().get(i).getAverageMapScore() + "\t" + Array.toStr(cnMappability.getMappabilityResults().get(i).getSubsetNames()));
				writer.println();
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + output);
			log.reportException(e);
		}

	}

	private static class SummaryCNV {
		private int[] numOverLaps;
		private int[] numSigOverLaps;
		private int[] numSigOverLapsSameCopyNumber;
		private String[] removeFiles;
		private int[] totalNumInds;

		public SummaryCNV(String[] removeFiles) {
			super();
			this.removeFiles = removeFiles;
			this.numOverLaps = new int[removeFiles.length];
			this.numSigOverLaps = new int[removeFiles.length];
			this.numSigOverLapsSameCopyNumber = new int[removeFiles.length];
			this.totalNumInds = new int[removeFiles.length];
		}

		public String[] getHeader() {
			ArrayList<String> header = new ArrayList<String>();
			for (int i = 0; i < removeFiles.length; i++) {
				String root = ext.removeDirectoryInfo(removeFiles[i]);
				header.add(root + "_" + "NUM_OVERLAP");
				header.add(root + "_" + "NUM_SIG_OVERLAP");
				header.add(root + "_" + "NUM_SIG_OVERLAP_SAME_CN");
				header.add(root + "_" + "NUM_INDS");
			}
			return Array.toStringArray(header);
		}

		public String[] getData() {
			ArrayList<String> data = new ArrayList<String>();
			for (int i = 0; i < removeFiles.length; i++) {
				data.add(numOverLaps[i] + "");
				data.add(numSigOverLaps[i] + "");
				data.add(numSigOverLapsSameCopyNumber[i] + "");
				data.add(totalNumInds[i] + "");

			}
			return Array.toStringArray(data);
		}

		public int[] getNumOverLaps() {
			return numOverLaps;
		}

		public int[] getTotalNumInds() {
			return totalNumInds;
		}

		public int[] getNumSigOverLaps() {
			return numSigOverLaps;
		}

		public int[] getNumSigOverLapsSameCopyNumber() {
			return numSigOverLapsSameCopyNumber;
		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String mappabilityFile = "Mappability.bed";
		String cnvFile = "cnvs.cnv";
		String geneTrackFile = "RefSeq_hg19.gtrack";
		String callSubsetBed = "exons.hg19.bed";
		String[] cnvFreqFiles = new String[] { "file1.cnv,file2.cnv" };
		String usage = "\n" + "one.JL.Mappability requires 0-1 arguments\n";
		usage += "   (1) mappability file (i.e. mapFile=" + mappabilityFile + " (default))\n" + "";
		usage += "   (2) cnv file (i.e. cnvs=" + cnvFile + " (default))\n" + "";
		usage += "   (3) geneTrackFile  (i.e. genes=" + cnvFile + " (default))\n" + "";
		usage += "   (4) call subsetBed  (i.e. callSubset=" + callSubsetBed + " (default))\n" + "";
		usage += "   (5) comma-Delimited list of files to remove  (i.e. cnvFreqFiles=" + Array.toStr(cnvFreqFiles, ",") + " (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("mapFile=")) {
				mappabilityFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("cnvFile=")) {
				cnvFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("genes=")) {
				geneTrackFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("callSubset=")) {
				callSubsetBed = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("cnvFreqFiles=")) {
				cnvFreqFiles = args[i].split("=")[1].split(",");
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
			Logger log = new Logger(ext.rootOf(cnvFile, false) + ".mappability.log");
			filter(mappabilityFile, cnvFile, cnvFreqFiles, geneTrackFile, callSubsetBed, log);
		} catch (Exception e) {

			e.printStackTrace();
		}
	}

}
