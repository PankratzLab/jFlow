package jlDev;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;

import seq.manage.BedOps;
import seq.qc.Mappability;
import seq.qc.Mappability.MappabilityResult;
import stats.Histogram.DynamicAveragingHistogram;
import stats.Rscript.GeomText;
import stats.Rscript.RScatter;
import stats.Rscript.SCATTER_TYPE;
import cnv.filesys.Project;
import cnv.var.CNVariant;
import cnv.var.LocusSet;
import cnv.var.SampleData;
import common.Array;
import common.CmdLine;
import common.Files;
import common.Logger;
import common.ext;
import filesys.Segment.SegmentCompare;

public class CushingCnvs {
	private static final String[] BASE_HEADER = new String[] { "MAPPABILITY_SCORE", "CALLED_GENE(s)" };
	private static final String PLINK = "C:/bin/plink/plink-1.07-dos/plink-1.07-dos/plink.exe";
	private static final String PLINKGENELIST = "C:/bin/plink/glist-hg19.txt";

	public static void filter(Project proj, String mappabilityFile, String[] cnvFiles, String[] cnvRemoveFiles, String geneTrackFile, String callSubsetBed, Logger log) {
		for (int i = 0; i < cnvFiles.length; i++) {

			filter(proj, mappabilityFile, cnvFiles[i], cnvRemoveFiles, geneTrackFile, callSubsetBed, log);
		}

	}

	public static void generatePed(Project proj, String cnvFile1, String cnvFile2) {
		Logger log = new Logger(ext.parseDirectoryOfFile(cnvFile1) + "generatePed.log");
		LocusSet<CNVariant> cnSet1 = CNVariant.loadLocSet(cnvFile1, log);
		LocusSet<CNVariant> cnSet2 = CNVariant.loadLocSet(cnvFile2, log);
		HashSet<String> fidIid1 = CNVariant.getUniqueInds(cnSet1, log);
		HashSet<String> fidIid2 = CNVariant.getUniqueInds(cnSet2, log);
		String dir = ext.parseDirectoryOfFile(cnvFile1) + "plink/";
		new File(dir).mkdirs();
		String outCNVRoot = dir + ext.rootOf(cnvFile1) + "_" + ext.rootOf(cnvFile2);
		String outped = outCNVRoot + ".fam";
		String outCNV = outCNVRoot + ".cnv";

		try {
			PrintWriter writer = new PrintWriter(new FileWriter(outCNV));
			writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER));
			for (int i = 0; i < cnSet1.getLoci().length; i++) {
				writer.println(cnSet1.getLoci()[i].toPlinkFormat());
			}
			for (int i = 0; i < cnSet2.getLoci().length; i++) {
				writer.println(cnSet2.getLoci()[i].toPlinkFormat());
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + outCNV);
			log.reportException(e);
		}

		SampleData sampleData = proj.getSampleData(0, false);

		try {
			PrintWriter writer = new PrintWriter(new FileWriter(outped));

			for (String caseInd : fidIid1) {
				String caseDNA = sampleData.lookup(caseInd)[0];

				if (!sampleData.individualShouldBeExcluded(caseDNA)) {
					if (caseInd.startsWith("CONTROL") || caseInd.startsWith("HapMap")) {
						writer.println(caseInd + "\t0\t0\t" + sampleData.getSexForIndividual(caseDNA) + "\t1");

					} else {
						writer.println(caseInd + "\t0\t0\t" + sampleData.getSexForIndividual(caseDNA) + "\t2");
					}
				}

			}
			for (String controlInd : fidIid2) {
				String controlDNA = sampleData.lookup(controlInd)[0];

				if (!sampleData.individualShouldBeExcluded(controlDNA)) {
					writer.println(controlInd + "\t0\t0\t" + sampleData.getSexForIndividual(controlDNA) + "\t1");
				}
			}

			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + outCNVRoot);
			log.reportException(e);
		}

		ArrayList<String> commandArray = new ArrayList<String>();
		commandArray.add(PLINK);
		commandArray.add("--cnv-list");
		commandArray.add(outCNV);
		commandArray.add("--cnv-make-map");
		commandArray.add("--out");
		commandArray.add(outCNVRoot);
		commandArray.add("--allow-no-sex");

		CmdLine.runCommandWithFileChecks(Array.toStringArray(commandArray), dir, null, null, true, true, false, log);

		commandArray = new ArrayList<String>();
		commandArray.add(PLINK);
		commandArray.add("--cfile");
		commandArray.add(outCNVRoot);
		commandArray.add("--out");
		commandArray.add(outCNVRoot);
		commandArray.add("--allow-no-sex");
		commandArray.add("--mperm");
		commandArray.add("100000");
		
		CmdLine.runCommandWithFileChecks(Array.toStringArray(commandArray), dir, null, null, true, true, true, log);

		
		
		
//		
//		commandArray = new ArrayList<String>();
//		commandArray.add(PLINK);
//		commandArray.add("--cfile");
//		commandArray.add(outCNVRoot);
//		commandArray.add("--out");
//		commandArray.add(outCNVRoot);
//		commandArray.add("--allow-no-sex");
//		commandArray.add("--cnv-indiv-perm");
//		commandArray.add("--mperm");
//		commandArray.add("100");


		//CmdLine.runCommandWithFileChecks(Array.toStringArray(commandArray), dir, null, null, true, true, false, log);
		commandArray = new ArrayList<String>();
		commandArray.add(PLINK);
		commandArray.add("--cfile");
		commandArray.add(outCNVRoot);
		commandArray.add("--out");
		commandArray.add(outCNVRoot);
		commandArray.add("--cnv-count");
		commandArray.add(PLINKGENELIST);
		commandArray.add("--cnv-enrichment-test");

		System.exit(1);
	}

	public static void filter(Project proj, String mappabilityFile, String cnvFile, final String[] cnvRemoveFilesstart, String geneTrackFile, String callSubsetBed, Logger log) {
		BedOps.verifyBedIndex(mappabilityFile, log);
		// LocusSet<GeneData> gLocusSet = GeneTrack.load(geneTrackFile, false).convertToLocusSet(log);
		LocusSet<CNVariant> cLocusSet = CNVariant.loadLocSet(cnvFile, log);
		ArrayList<String> cnvFreqFiles = new ArrayList<String>();
		int numIndss = CNVariant.getUniqueInds(cLocusSet, log).size();
		System.out.println(cnvFile);
		System.out.println(numIndss);

		for (int i = 0; i < cnvRemoveFilesstart.length; i++) {
			if (Files.isDirectory(cnvRemoveFilesstart[i])) {
				String[] tmpsCnvs = Files.list(cnvRemoveFilesstart[i], null, ".cnv", true, false, true);
				for (int j = 0; j < tmpsCnvs.length; j++) {
					cnvFreqFiles.add(tmpsCnvs[j]);
				}
			} else {
				cnvFreqFiles.add(cnvRemoveFilesstart[i]);
			}
		}
		cnvFreqFiles.add(cnvFile);
		String[] cnvRemoveFiles = Array.toStringArray(cnvFreqFiles);
		if (proj != null) {
			for (int i = 0; i < cnvRemoveFiles.length; i++) {
				if (cnvRemoveFiles[i].contains("OSTEO")) {
					generatePed(proj, cnvFile, cnvRemoveFiles[i]);
				}

			}
			return;

		}

		log.reportTimeInfo("found " + cnvRemoveFiles.length + " files to check");
		SummaryCNV[] summaryCNVs = new SummaryCNV[cLocusSet.getLoci().length];
		for (int i = 0; i < summaryCNVs.length; i++) {
			summaryCNVs[i] = new SummaryCNV(cnvRemoveFiles);
		}
		for (int i = 0; i < cnvRemoveFiles.length; i++) {
			log.reportTimeInfo("Checking file " + (i + 1) + " of " + cnvRemoveFiles.length);
			LocusSet<CNVariant> cLocusRemoveSet = null;
			String ser = cnvRemoveFiles[i] + ".ser";
			if (Files.exists(ser)) {
				log.reportTimeInfo("Loading " + ser);
				cLocusRemoveSet = (LocusSet<CNVariant>) LocusSet.readSerialCnvSet(ser, log);
			} else {
				cLocusRemoveSet = CNVariant.loadLocSet(cnvRemoveFiles[i], log);
				cLocusRemoveSet.writeSerial(ser);
			}

			int numInds = CNVariant.getUniqueInds(cLocusRemoveSet, log).size();
			log.reportTimeInfo("Loaded " + cLocusRemoveSet.getLoci().length + "  cnvs from " + cnvRemoveFiles[i] + " with " + numInds + " total samples");

			for (int j = 0; j < cLocusSet.getLoci().length; j++) {
				summaryCNVs[j].getTotalNumInds()[i] = numInds;
				CNVariant currentCNV = cLocusSet.getLoci()[j];
				CNVariant[] overlaps = cLocusRemoveSet.getOverLappingLoci(currentCNV);

				if (overlaps != null && overlaps.length > 0) {
					summaryCNVs[j].getNumOverLaps()[i] += overlaps.length;

					int numSig = 0;
					int numSigCN = 0;
					ArrayList<CNVariant> sigOlaps = new ArrayList<CNVariant>();
					for (int k = 0; k < overlaps.length; k++) {
						if (currentCNV.significantOverlap(overlaps[k], true)) {
							numSig++;
							sigOlaps.add(overlaps[k]);
							if (currentCNV.getCN() == overlaps[k].getCN() || (currentCNV.getCN() < 2 && overlaps[k].getCN() < 2) || (currentCNV.getCN() > 2 && overlaps[k].getCN() > 2)) {
								numSigCN++;
							}
						}
					}
					SegmentCompare segmentCompareAll = currentCNV.new SegmentCompare(overlaps, 0, log);
					segmentCompareAll.compare();
					SegmentCompare segmentCompareSig = currentCNV.new SegmentCompare(sigOlaps.toArray(new CNVariant[sigOlaps.size()]), 0, log);
					segmentCompareSig.compare();

					summaryCNVs[j].getAverageOverlapScore()[i] = segmentCompareAll.getAvgOverlapScore();
					summaryCNVs[j].getNumSigOverLaps()[i] += numSig;
					summaryCNVs[j].getAverageSigOverlapScore()[i] = segmentCompareSig.getAvgOverlapScore();
					summaryCNVs[j].getNumSigOverLapsSameCNDirection()[i] += numSigCN;

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
				writer.print("\t" + cnMappability.getMappabilityResults().get(i).getAverageMapScore() + "\t" + Array.toStr(cnMappability.getMappabilityResults().get(i).getSubsetNames(), "/"));
				writer.println();
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + output);
			log.reportException(e);
		}

		Hashtable<String, Integer> geneCounts = cnMappability.generateInternalGeneCounts();
		int maxNumPer = 0;
		for (String gene : geneCounts.keySet()) {
			if (geneCounts.get(gene) > maxNumPer) {
				maxNumPer = geneCounts.get(gene);
			}
		}
		DynamicAveragingHistogram dynamicAveragingHistogramCNVCentered = new DynamicAveragingHistogram(0, 150, 0);
		String geneCountsFile = ext.rootOf(cnvFile, false) + "geneCounts.txt";
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(geneCountsFile));
			writer.println("GENE\tCNV_COUNT");
			for (String aGene : geneCounts.keySet()) {
				writer.println(aGene + "\t" + geneCounts.get(aGene));
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + geneCountsFile);
			log.reportException(e);
		}

		String outputRawPlot = ext.rootOf(cnvFile, false) + ".qc.summary.plot.txt";

		String[] rawCounts = new String[] { "CNVS_PER_GENE", "MAPPABILITY_SCORE" };
		Hashtable<String, String> problemsAdded = new Hashtable<String, String>();
		ArrayList<GeomText> problemGenes = new ArrayList<GeomText>();
		ArrayList<GeomText> problemGenesSmall = new ArrayList<GeomText>();

		try {
			PrintWriter writer = new PrintWriter(new FileWriter(outputRawPlot));
			writer.println(Array.toStr(rawCounts));

			for (int i = 0; i < cnMappability.getMappabilityResults().size(); i++) {
				MappabilityResult<CNVariant> cnMapp = cnMappability.getMappabilityResults().get(i);
				for (int l = 0; l < cnMapp.getSubsetNames().length; l++) {
					String gene = cnMapp.getSubsetNames()[l];
					int count = geneCounts.get(gene);
					double mapScore = cnMapp.getAverageMapScore();
					dynamicAveragingHistogramCNVCentered.addDataPair(count, mapScore);
					writer.println(count + "\t" + cnMapp.getAverageMapScore());
					if ((gene.startsWith("ZNF") || (gene.startsWith("OR") && !gene.startsWith("ORM") && !gene.startsWith("ORC")) || gene.startsWith("MUC")) && !problemsAdded.containsKey(gene + "_" + count + "_" + mapScore)) {
						GeomText geomText = new GeomText(count, mapScore, 0, "*", 5);
						problemsAdded.put(gene + "_" + count + "_" + mapScore, "DFDSF");
						problemGenes.add(geomText);
						problemGenesSmall.add(new GeomText(count, mapScore, 0, gene, 5));
					}
					if (count > 75 && !problemsAdded.containsKey(gene)) {
						GeomText geomText = new GeomText(count, mapScore, 0, gene, 5);
						problemGenes.add(geomText);
						problemsAdded.put(gene, gene);
						problemGenesSmall.add(geomText);

					}
				}
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + outputRawPlot);
			log.reportException(e);
		}
		if (problemGenes.size() < 50) {
			problemGenes = problemGenesSmall;
		}

		RScatter rscatterNoGeom = new RScatter(outputRawPlot, outputRawPlot + ".no_geom.rscript", ext.removeDirectoryInfo(outputRawPlot) + ".no_geom", outputRawPlot + ".no_geom" + ".jpeg", rawCounts[0], new String[] { rawCounts[1] }, SCATTER_TYPE.POINT, log);
		rscatterNoGeom.setxLabel("CNVS Per Gene");
		rscatterNoGeom.setyRange(new double[] { 0, 1 });
		rscatterNoGeom.setxRange(new double[] { 0, 150 });

		rscatterNoGeom.setyLabel("CNV mappability score");
		rscatterNoGeom.setOverWriteExisting(true);
		rscatterNoGeom.execute();

		RScatter rscatter = new RScatter(outputRawPlot, outputRawPlot + ".rscript", ext.removeDirectoryInfo(outputRawPlot), outputRawPlot + ".jpeg", rawCounts[0], new String[] { rawCounts[1] }, SCATTER_TYPE.POINT, log);
		rscatter.setxLabel("CNVS Per Gene");
		rscatter.setyRange(new double[] { 0, 1 });
		rscatter.setxRange(new double[] { 0, 150 });

		rscatter.setyLabel("CNV mappability score");
		rscatter.setgTexts(problemGenes.toArray(new GeomText[problemGenes.size()]));
		rscatter.setOverWriteExisting(true);
		rscatter.execute();

		String hist = ext.addToRoot(outputRawPlot, ".hist");
		String[] dump = new String[] { rawCounts[0], "NUM_GENES_WITH_COUNT", "AVG_" + rawCounts[1] };
		dynamicAveragingHistogramCNVCentered.average();
		dynamicAveragingHistogramCNVCentered.dump(hist, dump, true);

		RScatter rscatterHist = new RScatter(hist, hist + ".rscript", ext.removeDirectoryInfo(hist), hist + ".jpeg", dump[0], new String[] { dump[2] }, SCATTER_TYPE.POINT, log);
		rscatterHist.setxLabel("CNVS Per Gene");
		rscatterHist.setyRange(new double[] { 0, 1 });
		rscatterHist.setxRange(new double[] { 0, 150 });

		rscatterHist.setyLabel("Average CNV mappability score");
		rscatterHist.setgTexts(problemGenes.toArray(new GeomText[problemGenes.size()]));
		rscatterHist.setOverWriteExisting(true);
		rscatterHist.execute();

		RScatter rscatterHistNoGeom = new RScatter(hist, hist + ".no_geom.rscript", ext.removeDirectoryInfo(hist) + ".no_geom", hist + ".no_geom" + ".jpeg", dump[0], new String[] { dump[2] }, SCATTER_TYPE.POINT, log);
		rscatterHistNoGeom.setxLabel("CNVS Per Gene");
		rscatterHistNoGeom.setyRange(new double[] { 0, 1 });
		rscatterHistNoGeom.setxRange(new double[] { 0, 150 });

		rscatterHistNoGeom.setyLabel("Average CNV mappability score");
		rscatterHistNoGeom.setOverWriteExisting(true);
		rscatterHistNoGeom.execute();

		RScatter rscatterHistNumper = new RScatter(hist, hist + ".numPer.rscript", ext.removeDirectoryInfo(hist) + ".numPer", hist + ".numPer" + ".jpeg", dump[0], new String[] { dump[1] }, SCATTER_TYPE.POINT, log);
		rscatterHistNumper.setxLabel("CNVS Per Gene");
		rscatterHistNumper.setyRange(new double[] { 0, 1 });
		rscatterHistNumper.setxRange(new double[] { 0, 150 });

		rscatterHistNumper.setyLabel("Genes with this number of CNVs");
		rscatterHistNumper.setOverWriteExisting(true);
		rscatterHistNumper.execute();
	}

	private static class SummaryCNV {
		private int[] numOverLaps;
		private int[] numSigOverLaps;
		private int[] numSigOverLapsSameCopyNumber;
		private double[] averageOverlapScore;
		private double[] averageSigOverlapScore;

		private String[] removeFiles;
		private int[] totalNumInds;

		public SummaryCNV(String[] removeFiles) {
			super();
			this.removeFiles = removeFiles;
			this.numOverLaps = new int[removeFiles.length];
			this.averageOverlapScore = new double[removeFiles.length];
			this.numSigOverLaps = new int[removeFiles.length];
			this.averageSigOverlapScore = new double[removeFiles.length];

			this.numSigOverLapsSameCopyNumber = new int[removeFiles.length];
			this.totalNumInds = new int[removeFiles.length];
		}

		public String[] getHeader() {
			ArrayList<String> header = new ArrayList<String>();
			for (int i = 0; i < removeFiles.length; i++) {
				String root = ext.removeDirectoryInfo(removeFiles[i]);
				header.add(root + "_" + "NUM_OVERLAP");
				header.add(root + "_" + "AVG_OVERLAP_SCORE");

				header.add(root + "_" + "NUM_SIG_OVERLAP");
				header.add(root + "_" + "AVG_SIG_OVERLAP_SCORE");

				header.add(root + "_" + "NUM_SIG_OVERLAP_SAME_CN");
				header.add(root + "_" + "NUM_INDS");
			}
			return Array.toStringArray(header);
		}

		public String[] getData() {
			ArrayList<String> data = new ArrayList<String>();
			for (int i = 0; i < removeFiles.length; i++) {
				data.add(numOverLaps[i] + "");
				data.add(averageOverlapScore[i] + "");

				data.add(numSigOverLaps[i] + "");
				data.add(averageSigOverlapScore[i] + "");

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

		public int[] getNumSigOverLapsSameCNDirection() {
			return numSigOverLapsSameCopyNumber;
		}

		public double[] getAverageOverlapScore() {
			return averageOverlapScore;
		}

		public double[] getAverageSigOverlapScore() {
			return averageSigOverlapScore;
		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String mappabilityFile = "Mappability.bed";
		String[] cnvFiles = new String[] { "cnvs1.cnv", "cnvs2.cnv" };
		String geneTrackFile = "RefSeq_hg19.gtrack";
		String callSubsetBed = "exons.hg19.bed";
		String[] cnvFreqFiles = new String[] { "file1.cnv,file2.cnv" };
		String fileName = null;
		String usage = "\n" + "one.JL.Mappability requires 0-1 arguments\n";
		usage += "   (1) mappability file (i.e. mapFile=" + mappabilityFile + " (default))\n" + "";
		usage += "   (2) cnv files (i.e. cnvs=" + Array.toStr(cnvFiles, ",") + " (default))\n" + "";
		usage += "   (3) geneTrackFile  (i.e. genes=" + cnvFiles + " (default))\n" + "";
		usage += "   (4) call subsetBed  (i.e. callSubset=" + callSubsetBed + " (default))\n" + "";
		usage += "   (5) comma-Delimited list of files to remove  (i.e. cnvFreqFiles=" + Array.toStr(cnvFreqFiles, ",") + " (default))\n" + "";
		usage += "   (6) project file  (i.e. proj=" + null + " (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("mapFile=")) {
				mappabilityFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("cnvFile=")) {
				cnvFiles = args[i].split("=")[1].split(",");
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
			} else if (args[i].startsWith("proj=")) {
				fileName = args[i].split("=")[1];
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
			Project proj = null;
			if (fileName != null) {
				proj = new Project(fileName, false);
			}
			Logger log = new Logger(ext.rootOf(cnvFiles[0], false) + ".mappability.log");
			filter(proj, mappabilityFile, cnvFiles, cnvFreqFiles, geneTrackFile, callSubsetBed, log);
		} catch (Exception e) {

			e.printStackTrace();
		}
	}

}
