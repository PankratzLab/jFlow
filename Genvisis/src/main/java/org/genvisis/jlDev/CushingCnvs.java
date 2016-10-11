package org.genvisis.jlDev;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.concurrent.Callable;

import org.genvisis.cnv.analysis.ProjectCNVFiltering;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.Array;
import org.genvisis.common.CNVFilter.FreqFilter;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.WorkerHive;
import org.genvisis.common.ext;
import org.genvisis.filesys.CNVariant;
import org.genvisis.filesys.LocusSet;
import org.genvisis.filesys.LocusSet.TO_STRING_TYPE;
import org.genvisis.filesys.Segment;
import org.genvisis.filesys.Segment.SegmentCompare;
import org.genvisis.seq.manage.BedOps;
import org.genvisis.seq.qc.Mappability;
import org.genvisis.seq.qc.Mappability.MappabilityResult;
import org.genvisis.stats.Histogram.DynamicAveragingHistogram;
import org.genvisis.stats.Rscript.GeomText;
import org.genvisis.stats.Rscript.RScatter;
import org.genvisis.stats.Rscript.SCATTER_TYPE;

public class CushingCnvs {
	private static final String[] BASE_HEADER = new String[] {"MAPPABILITY_SCORE", "CALLED_GENE(s)"};
	// private static final String PLINK = "C:/bin/plink/plink-1.07-dos/plink-1.07-dos/plink.exe";
	private static final String PLINK = "plink";

	// private static final String PLINKGENELIST = "C:/bin/plink/glist-hg19.txt";

	public static void filter(Project proj, String mappabilityFile, String[] cnvFiles,
														String[] cnvRemoveFiles, String geneTrackFile, String callSubsetBed,
														Logger log) {
		for (String cnvFile : cnvFiles) {
			filter(proj, mappabilityFile, cnvFile, cnvRemoveFiles, geneTrackFile, callSubsetBed, log);
		}

	}

	private static class PlinkEmpSeg extends Segment {
		/**
		 *
		 */
		private static final long serialVersionUID = 1L;
		private final double emp1;
		private final double emp2;
		private final String inputFile;
		private static final String[] EMP_HEADER = new String[] {"CHR", "SNP", "EMP1", "EMP2"};
		private static final String[] REPORT_HEADER = new String[] {"CHR", "POS", "EMP1", "EMP2"};

		private PlinkEmpSeg(byte chr, int pos, double emp1, double emp2, String inputFile) {
			super(chr, pos, pos);
			this.emp1 = emp1;
			this.emp2 = emp2;
			this.inputFile = inputFile;
		}

		public double getEmp1() {
			return emp1;
		}

		public double getEmp2() {
			return emp2;
		}

		private String[] getPlinkHeader() {
			String root = ext.rootOf(inputFile);
			String[] header = Array.concatAll(Array.tagOn(REPORT_HEADER, root + "_MIN_EMP", null));
			return header;
		}

		private String[] getData() {
			ArrayList<String> data = new ArrayList<String>();
			data.add(getChr() + "");
			data.add(getStart() + "");
			data.add(emp1 + "");
			data.add(emp2 + "");
			return Array.toStringArray(data);

		}

		private static LocusSet<PlinkEmpSeg> loadPlLocusSet(String filname, Logger log) {
			ArrayList<PlinkEmpSeg> tmp = new ArrayList<PlinkEmpSeg>();
			try {
				BufferedReader reader = Files.getAppropriateReader(filname);
				String[] header = reader.readLine().trim().split("[\\s]+");
				int[] indices = ext.indexFactors(header, EMP_HEADER, true, false);

				while (reader.ready()) {
					String[] line = reader.readLine().trim().split("[\\s]+");
					byte chr = Byte.parseByte(line[indices[0]]);
					String tmpLoc = line[indices[1]];
					tmpLoc = tmpLoc.replaceAll(".*-", "");
					int pos = Integer.parseInt(tmpLoc);
					double emp1 = Double.parseDouble(line[indices[2]]);
					double emp2 = Double.parseDouble(line[indices[3]]);
					tmp.add(new PlinkEmpSeg(chr, pos, emp1, emp2, filname));
				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				log.reportError("Error: file \"" + filname + "\" not found in current directory");
				return null;
			} catch (IOException ioe) {
				log.reportError("Error reading file \"" + filname + "\"");
				return null;
			}
			return new LocusSet<PlinkEmpSeg>(tmp.toArray(new PlinkEmpSeg[tmp.size()]), true, log) {

				/**
				 *
				 */
				private static final long serialVersionUID = 1L;

			};
		}

	}

	private static boolean use(CNVariant cnVariant, LocusSet<CNVariant> setToRemove, Logger log) {
		CNVariant[] overlap = setToRemove.getOverLappingLoci(cnVariant);
		if (overlap == null || overlap.length == 0) {
			return true;
		}
		for (CNVariant element : overlap) {
			double percentOverlap = (double) cnVariant.amountOfOverlapInBasepairs(element)
															/ cnVariant.getSize();
			if (percentOverlap > .5) {// exon is commonly deleted
				return false;
			}
		}
		return true;
	}

	// if (controlFreqFilter < .5 && cnSet1.getLoci()[i].getStop() == 1569278) {
	// CNVariant[] olap = setToRemove.getOverLappingLoci(cnSet1.getLoci()[i]);
	// for (int j = 0; j < setToRemove.getLoci().length; j++) {
	// // System.out.println(setToRemove.getLoci()[j].toPlinkFormat());
	// }
	// for (int j = 0; j < olap.length; j++) {
	// System.out.println(olap[j].toPlinkFormat() + "\t" + cnSet1.getLoci()[i].toPlinkFormat());
	// System.out.println(cnSet1.getLoci()[i].significantOverlap(olap[j], true));
	// }
	// //System.exit(1);
	//
	// }

	public static class PlinkWorker implements Callable<LocusSet<PlinkEmpSeg>> {
		private final Project proj;
		private final String cnvFile1;
		private final String cnvFile2;
		private final double controlFreqFilter;
		private final double conf;

		public PlinkWorker(	Project proj, String cnvFile1, String cnvFile2, double controlFreqFilter,
												double conf) {
			super();
			this.proj = proj;
			this.cnvFile1 = cnvFile1;
			this.cnvFile2 = cnvFile2;
			this.controlFreqFilter = controlFreqFilter;
			this.conf = conf;
		}

		@Override
		public LocusSet<PlinkEmpSeg> call() throws Exception {
			// TODO Auto-generated method stub
			return generatePed(proj, cnvFile1, cnvFile2, controlFreqFilter, conf);
		}
	}

	public static LocusSet<PlinkEmpSeg> generatePed(Project proj, String cnvFile1, String cnvFile2,
																									double controlFreqFilter, double conf) {
		Logger log = new Logger(ext.parseDirectoryOfFile(cnvFile1) + "generatePed.log");
		LocusSet<CNVariant> cnSet1 = CNVariant.loadLocSet(cnvFile1, log);
		LocusSet<CNVariant> cnSet2 = CNVariant.loadLocSet(cnvFile2, log);
		HashSet<String> fidIid1 = CNVariant.getUniqueInds(cnSet1, log);
		HashSet<String> fidIid2 = CNVariant.getUniqueInds(cnSet2, log);
		String dir = ext.parseDirectoryOfFile(cnvFile1) + "plink/";
		new File(dir).mkdirs();
		String outCNVRoot = dir	+ ext.rootOf(cnvFile1) + "_" + ext.rootOf(cnvFile2) + "maxControlFreq"
												+ controlFreqFilter + "_conf" + conf;
		String outped = outCNVRoot + ".fam";
		String outCNV = outCNVRoot + ".cnv";
		String outFilt = outCNVRoot + ".freqFilteredCnvs";
		int numControls = CNVariant.getUniqueInds(cnSet2, proj.getLog()).size();
		int totalLimitedTo = Math.round((float) controlFreqFilter * numControls);
		// public FreqFilter(int totalRequired, int delRequired, int dupRequired, int totalLimitedTo,
		// int delLimitedTo, int dupLimitedTo, double proportionOfProbesThatNeedToPassForFinalInclusion)
		// {

		FreqFilter freqFilter = new FreqFilter(	totalLimitedTo, 0, 0, numControls, numControls,
																						numControls, .9);
		// CNVFilter cnvFilter = new CNVFilter(proj.getLog());
		LocusSet<CNVariant> setToRemove = null;
		setToRemove = ProjectCNVFiltering.filterCNVFile(proj, cnSet2.getLoci(), outFilt, null, false,
																										true, freqFilter, false, true);
		setToRemove.writeRegions(outFilt, TO_STRING_TYPE.REGULAR, true, proj.getLog());

		if (!Files.exists(outFilt)) {
		} else {
		}

		try {
			PrintWriter writer = new PrintWriter(new FileWriter(outCNV));
			writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER));
			int caseRemoved = 0;
			int controlRemoved = 0;
			for (int i = 0; i < cnSet1.getLoci().length; i++) {
				if (use(cnSet1.getLoci()[i], setToRemove, log) && cnSet1.getLoci()[i].getScore() > conf) {

					writer.println(cnSet1.getLoci()[i].toPlinkFormat());
				} else {
					caseRemoved++;
				}
			}
			for (int i = 0; i < cnSet2.getLoci().length; i++) {
				if (use(cnSet2.getLoci()[i], setToRemove, log) && cnSet2.getLoci()[i].getScore() > conf) {
					writer.println(cnSet2.getLoci()[i].toPlinkFormat());
				} else {
					controlRemoved++;
				}
			}
			writer.close();
			log.reportTimeInfo("Removed "	+ caseRemoved + " cnvs of " + cnSet1.getLoci().length
													+ " from cases and " + controlRemoved + " cnvs of "
													+ cnSet2.getLoci().length + " from controls for freq " + controlFreqFilter
													+ " and conf " + conf);
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
					writer.println(controlInd	+ "\t0\t0\t" + sampleData.getSexForIndividual(controlDNA)
													+ "\t1");
				}
			}

			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + outCNVRoot);
			log.reportException(e);
		}

		String mperm = outCNVRoot + ".cnv.summary.mperm";
		String map = outCNVRoot + ".cnv.map";
		new File(map).delete();
		new File(mperm).delete();
		System.out.println(mperm);
		if (!Files.exists(mperm)) {
			ArrayList<String> commandArray = new ArrayList<String>();
			commandArray.add(PLINK);
			commandArray.add("--cnv-list");
			commandArray.add(outCNV);
			commandArray.add("--cnv-make-map");
			commandArray.add("--out");
			commandArray.add(outCNVRoot);
			commandArray.add("--allow-no-sex");
			commandArray.add("--noweb");


			CmdLine.runCommandWithFileChecks(	Array.toStringArray(commandArray), dir, null, null, true,
																				true, false, log);

			commandArray = new ArrayList<String>();
			commandArray.add(PLINK);
			commandArray.add("--cfile");
			commandArray.add(outCNVRoot);
			commandArray.add("--out");
			commandArray.add(outCNVRoot);
			commandArray.add("--allow-no-sex");
			commandArray.add("--mperm");
			commandArray.add("10000");
			commandArray.add("--noweb");


			CmdLine.runCommandWithFileChecks(	Array.toStringArray(commandArray), dir, null, null, true,
																				true, true, log);

			//
			// commandArray = new ArrayList<String>();
			// commandArray.add(PLINK);
			// commandArray.add("--cfile");
			// commandArray.add(outCNVRoot);
			// commandArray.add("--out");
			// commandArray.add(outCNVRoot);
			// commandArray.add("--allow-no-sex");
			// commandArray.add("--cnv-indiv-perm");
			// commandArray.add("--mperm");
			// commandArray.add("100");

			// CmdLine.runCommandWithFileChecks(Array.toStringArray(commandArray), dir, null, null, true,
			// true, false, log);
			// commandArray = new ArrayList<String>();
			// commandArray.add(PLINK);
			// commandArray.add("--cfile");
			// commandArray.add(outCNVRoot);
			// commandArray.add("--out");
			// commandArray.add(outCNVRoot);
			// commandArray.add("--cnv-count");
			// commandArray.add(PLINKGENELIST);
			// commandArray.add("--cnv-enrichment-test");
		}
		return PlinkEmpSeg.loadPlLocusSet(mperm, log);
	}

	public static void filter(Project proj, String mappabilityFile, String cnvFile,
														final String[] cnvRemoveFilesstart, String geneTrackFile,
														String callSubsetBed, Logger log) {
		BedOps.verifyBedIndex(mappabilityFile, log);
		// LocusSet<GeneData> gLocusSet = GeneTrack.load(geneTrackFile, false).convertToLocusSet(log);
		LocusSet<CNVariant> cLocusSet = CNVariant.loadLocSet(cnvFile, log);
		ArrayList<String> cnvFreqFiles = new ArrayList<String>();
		int numIndss = CNVariant.getUniqueInds(cLocusSet, log).size();
		System.out.println(cnvFile);
		System.out.println(numIndss);
		ArrayList<LocusSet<PlinkEmpSeg>> plinkResults = new ArrayList<LocusSet<PlinkEmpSeg>>();
		for (String element : cnvRemoveFilesstart) {
			if (Files.isDirectory(element)) {
				String[] tmpsCnvs = Files.list(element, null, ".cnv", true, false, true);
				for (String tmpsCnv : tmpsCnvs) {
					cnvFreqFiles.add(tmpsCnv);
				}
			} else {
				cnvFreqFiles.add(element);
			}
		}
		cnvFreqFiles.add(cnvFile);
		String[] cnvRemoveFiles = Array.toStringArray(cnvFreqFiles);
		double[] controlFreqFilter = new double[] {0.01, 0.05, 1.1};
		double[] confs = new double[] {10, 20, 0};
		WorkerHive<LocusSet<PlinkEmpSeg>> hive = new WorkerHive<LocusSet<PlinkEmpSeg>>(	6, 10,
																																										proj.getLog());

		if (proj != null) {
			for (String cnvRemoveFile : cnvRemoveFiles) {
				if (cnvRemoveFile.contains("control")) {
					for (double element : controlFreqFilter) {
						for (double conf : confs) {
							PlinkWorker worker = new PlinkWorker(proj, cnvFile, cnvRemoveFile, element, conf);
							hive.addCallable(worker);
						}
					}
				}
			}
			hive.execute(true);
			plinkResults.addAll(hive.getResults());
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
				cLocusRemoveSet = LocusSet.readSerialCnvSet(ser, log);
			} else {
				cLocusRemoveSet = CNVariant.loadLocSet(cnvRemoveFiles[i], log);
				cLocusRemoveSet.writeSerial(ser);
			}

			int numInds = CNVariant.getUniqueInds(cLocusRemoveSet, log).size();
			log.reportTimeInfo("Loaded "	+ cLocusRemoveSet.getLoci().length + "  cnvs from "
													+ cnvRemoveFiles[i] + " with " + numInds + " total samples");

			for (int j = 0; j < cLocusSet.getLoci().length; j++) {
				summaryCNVs[j].getTotalNumInds()[i] = numInds;
				CNVariant currentCNV = cLocusSet.getLoci()[j];
				CNVariant[] overlaps = cLocusRemoveSet.getOverLappingLoci(currentCNV);

				if (overlaps != null && overlaps.length > 0) {
					summaryCNVs[j].getNumOverLaps()[i] += overlaps.length;

					int numSig = 0;
					int numSigCN = 0;
					ArrayList<CNVariant> sigOlaps = new ArrayList<CNVariant>();
					for (CNVariant overlap : overlaps) {
						if (currentCNV.significantOverlap(overlap, true)) {
							numSig++;
							sigOlaps.add(overlap);
							if (currentCNV.getCN() == overlap.getCN()
										|| (currentCNV.getCN() < 2 && overlap.getCN() < 2)
									|| (currentCNV.getCN() > 2 && overlap.getCN() > 2)) {
								numSigCN++;
							}
						}
					}
					SegmentCompare segmentCompareAll = currentCNV.new SegmentCompare(overlaps, 0, log);
					segmentCompareAll.compare();
					SegmentCompare segmentCompareSig = currentCNV.new SegmentCompare(	sigOlaps.toArray(new CNVariant[sigOlaps.size()]),
																																						0, log);
					segmentCompareSig.compare();

					summaryCNVs[j].getAverageOverlapScore()[i] = segmentCompareAll.getAvgOverlapScore();
					summaryCNVs[j].getNumSigOverLaps()[i] += numSig;
					summaryCNVs[j].getAverageSigOverlapScore()[i] = segmentCompareSig.getAvgOverlapScore();
					summaryCNVs[j].getNumSigOverLapsSameCNDirection()[i] += numSigCN;

				}
			}
		}
		Mappability<CNVariant> cnMappability = new Mappability<CNVariant>(cLocusSet, mappabilityFile,
																																			callSubsetBed, log);
		cnMappability.computeMappability();

		String output = ext.rootOf(cnvFile, false) + ".qc.summary.txt";
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(output));
			writer.print("UCSC_LOC\tUCSC_LINK\t"
										+ Array.toStr(Array.concatAll(CNVariant.PLINK_CNV_HEADER,
																									summaryCNVs[0].getHeader(), BASE_HEADER)));
			for (int j = 0; j < plinkResults.size(); j++) {
				writer.print("\t" + Array.toStr(plinkResults.get(j).getLoci()[0].getPlinkHeader()));
			}
			writer.println();
			for (int i = 0; i < summaryCNVs.length; i++) {
				CNVariant currentCNV = cLocusSet.getLoci()[i];
				writer.print(currentCNV.getUCSClocation()	+ "\t" + currentCNV.getUCSCLink("hg19") + "\t"
											+ currentCNV.toPlinkFormat());
				writer.print("\t" + Array.toStr(summaryCNVs[i].getData()));
				writer.print("\t"	+ cnMappability.getMappabilityResults().get(i).getAverageMapScore() + "\t"
											+ Array.toStr(cnMappability.getMappabilityResults().get(i).getSubsetNames(),
																		"/"));
				for (int j = 0; j < plinkResults.size(); j++) {
					PlinkEmpSeg[] overLaps = plinkResults.get(j).getOverLappingLoci(currentCNV);
					if (overLaps == null || overLaps.length == 0 || overLaps.length < 2) {// need to find both
																																								// start and stop
																																								// positions
						writer.print("\t" + Array.toStr(Array.doubleArray(
																															plinkResults.get(j).getLoci()[0]
																																															.getData().length,
																															Double.NaN)));
						log.reportTimeError("Could not find overlapping plink results for "
																+ currentCNV.toPlinkFormat());
					} else {
						// int minEmp1Index = -1;
						int minEmp2Index = -1;
						double minEmp1 = Double.MAX_VALUE;
						double minEmp2 = Double.MAX_VALUE;
						boolean exactMatchStart = false;
						boolean exactMatchStop = false;
						for (int k = 0; k < overLaps.length; k++) {
							if (overLaps[k].getStart() == currentCNV.getStart()
									|| overLaps[k].getStop() == currentCNV.getStop()) {
								if (overLaps[k].getStart() == currentCNV.getStart()) {
									exactMatchStart = true;
								}
								if (overLaps[k].getStart() == currentCNV.getStop()) {
									exactMatchStop = true;
								}
								if (overLaps[k].getEmp1() < minEmp1) {

									minEmp1 = overLaps[k].getEmp1();
									// minEmp1Index = k;
								}
								if (overLaps[k].getEmp2() < minEmp2) {
									minEmp2 = overLaps[k].getEmp1();
									minEmp2Index = k;
								}
							}

						}
						System.out.println(output);
						if (exactMatchStart && exactMatchStop) {
							writer.print("\t" + Array.toStr(overLaps[minEmp2Index].getData()));
						} else {
							writer.print("\t" + Array.toStr(Array.doubleArray(overLaps[0].getData().length,
																																Double.NaN)));

						}
					}
				}

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
		DynamicAveragingHistogram dynamicAveragingHistogramCNVCentered = new DynamicAveragingHistogram(	0,
																																																		150,
																																																		0);
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

		String[] rawCounts = new String[] {"CNVS_PER_GENE", "MAPPABILITY_SCORE"};
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
					if ((gene.startsWith("ZNF")
									|| (gene.startsWith("OR") && !gene.startsWith("ORM") && !gene.startsWith("ORC"))
								|| gene.startsWith("MUC"))
							&& !problemsAdded.containsKey(gene + "_" + count + "_" + mapScore)) {
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

		RScatter rscatterNoGeom = new RScatter(	outputRawPlot, outputRawPlot + ".no_geom.rscript",
																						ext.removeDirectoryInfo(outputRawPlot) + ".no_geom",
																						outputRawPlot + ".no_geom" + ".jpeg", rawCounts[0],
																						new String[] {rawCounts[1]}, SCATTER_TYPE.POINT, log);
		rscatterNoGeom.setxLabel("CNVS Per Gene");
		rscatterNoGeom.setyRange(new double[] {0, 1});
		rscatterNoGeom.setxRange(new double[] {0, 150});

		rscatterNoGeom.setyLabel("CNV mappability score");
		rscatterNoGeom.setOverWriteExisting(true);
		rscatterNoGeom.execute();

		RScatter rscatter = new RScatter(	outputRawPlot, outputRawPlot + ".rscript",
																			ext.removeDirectoryInfo(outputRawPlot),
																			outputRawPlot + ".jpeg", rawCounts[0],
																			new String[] {rawCounts[1]}, SCATTER_TYPE.POINT, log);
		rscatter.setxLabel("CNVS Per Gene");
		rscatter.setyRange(new double[] {0, 1});
		rscatter.setxRange(new double[] {0, 150});

		rscatter.setyLabel("CNV mappability score");
		rscatter.setgTexts(problemGenes.toArray(new GeomText[problemGenes.size()]));
		rscatter.setOverWriteExisting(true);
		rscatter.execute();

		String hist = ext.addToRoot(outputRawPlot, ".hist");
		String[] dump = new String[] {rawCounts[0], "NUM_GENES_WITH_COUNT", "AVG_" + rawCounts[1]};
		dynamicAveragingHistogramCNVCentered.average();
		dynamicAveragingHistogramCNVCentered.dump(hist, dump, true);

		RScatter rscatterHist = new RScatter(	hist, hist + ".rscript", ext.removeDirectoryInfo(hist),
																					hist + ".jpeg", dump[0], new String[] {dump[2]},
																					SCATTER_TYPE.POINT, log);
		rscatterHist.setxLabel("CNVS Per Gene");
		rscatterHist.setyRange(new double[] {0, 1});
		rscatterHist.setxRange(new double[] {0, 150});

		rscatterHist.setyLabel("Average CNV mappability score");
		rscatterHist.setgTexts(problemGenes.toArray(new GeomText[problemGenes.size()]));
		rscatterHist.setOverWriteExisting(true);
		rscatterHist.execute();

		RScatter rscatterHistNoGeom = new RScatter(	hist, hist + ".no_geom.rscript",
																								ext.removeDirectoryInfo(hist) + ".no_geom",
																								hist + ".no_geom" + ".jpeg", dump[0],
																								new String[] {dump[2]}, SCATTER_TYPE.POINT, log);
		rscatterHistNoGeom.setxLabel("CNVS Per Gene");
		rscatterHistNoGeom.setyRange(new double[] {0, 1});
		rscatterHistNoGeom.setxRange(new double[] {0, 150});

		rscatterHistNoGeom.setyLabel("Average CNV mappability score");
		rscatterHistNoGeom.setOverWriteExisting(true);
		rscatterHistNoGeom.execute();

		RScatter rscatterHistNumper = new RScatter(	hist, hist + ".numPer.rscript",
																								ext.removeDirectoryInfo(hist) + ".numPer",
																								hist + ".numPer" + ".jpeg", dump[0],
																								new String[] {dump[1]}, SCATTER_TYPE.POINT, log);
		rscatterHistNumper.setxLabel("CNVS Per Gene");
		rscatterHistNumper.setyRange(new double[] {0, 1});
		rscatterHistNumper.setxRange(new double[] {0, 150});

		rscatterHistNumper.setyLabel("Genes with this number of CNVs");
		rscatterHistNumper.setOverWriteExisting(true);
		rscatterHistNumper.execute();
	}

	private static class SummaryCNV {
		private final int[] numOverLaps;
		private final int[] numSigOverLaps;
		private final int[] numSigOverLapsSameCopyNumber;
		private final double[] averageOverlapScore;
		private final double[] averageSigOverlapScore;

		private final String[] removeFiles;
		private final int[] totalNumInds;

		public SummaryCNV(String[] removeFiles) {
			super();
			this.removeFiles = removeFiles;
			numOverLaps = new int[removeFiles.length];
			averageOverlapScore = new double[removeFiles.length];
			numSigOverLaps = new int[removeFiles.length];
			averageSigOverlapScore = new double[removeFiles.length];

			numSigOverLapsSameCopyNumber = new int[removeFiles.length];
			totalNumInds = new int[removeFiles.length];
		}

		public String[] getHeader() {
			ArrayList<String> header = new ArrayList<String>();
			for (String removeFile : removeFiles) {
				String root = ext.removeDirectoryInfo(removeFile);
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
		String[] cnvFiles = new String[] {"cnvs1.cnv", "cnvs2.cnv"};
		String geneTrackFile = "RefSeq_hg19.gtrack";
		String callSubsetBed = "exons.hg19.bed";
		String[] cnvFreqFiles = new String[] {"file1.cnv,file2.cnv"};
		String fileName = null;
		String usage = "\n" + "one.JL.Mappability requires 0-1 arguments\n";
		usage += "   (1) mappability file (i.e. mapFile=" + mappabilityFile + " (default))\n" + "";
		usage += "   (2) cnv files (i.e. cnvs=" + Array.toStr(cnvFiles, ",") + " (default))\n" + "";
		usage +=
					"   (3) geneTrackFile  (i.e. genes=" + Array.toStr(cnvFiles, ",") + " (default))\n" + "";
		usage += "   (4) call subsetBed  (i.e. callSubset=" + callSubsetBed + " (default))\n" + "";
		usage += "   (5) comma-Delimited list of files to remove  (i.e. cnvFreqFiles="
							+ Array.toStr(cnvFreqFiles, ",") + " (default))\n" + "";
		usage += "   (6) project file  (i.e. proj=" + null + " (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("mapFile=")) {
				mappabilityFile = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("cnvFile=")) {
				cnvFiles = arg.split("=")[1].split(",");
				numArgs--;
			} else if (arg.startsWith("genes=")) {
				geneTrackFile = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("callSubset=")) {
				callSubsetBed = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("cnvFreqFiles=")) {
				cnvFreqFiles = arg.split("=")[1].split(",");
				numArgs--;
			} else if (arg.startsWith("proj=")) {
				fileName = arg.split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + arg);
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
