package org.genvisis.one.ben;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;

import org.genvisis.cnv.analysis.FilterCalls;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.Aliases;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.CNVFilter;
import org.genvisis.common.CNVFilter.CNVFilterPass;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.Positions;
import org.genvisis.common.ext;
import org.genvisis.filesys.CNVariant;
import org.genvisis.filesys.Segment;

public class lab {


	private static void countCNVsForIndividuals(String indivFile, String cnvFile,
																							String outFile) throws IOException {
		Hashtable<String, String> sampleKeyHash = new Hashtable<String, String>();
		BufferedReader reader = new BufferedReader(new FileReader(indivFile));
		String line = null;
		while ((line = reader.readLine()) != null) {
			String[] tmp = line.split("\t");
			sampleKeyHash.put(tmp[0] + "\t" + tmp[1], "");
		}
		reader.close();

		List<CNVariant> cnvs = CNVariant.loadPlinkFile(cnvFile, sampleKeyHash, true, false);

		PrintWriter writer = Files.getAppropriateWriter(outFile);
		for (CNVariant cnv : cnvs) {
			writer.println(cnv.toPlinkFormat());
		}
		writer.flush();
		writer.close();
	}

	private static void idSwap(Project proj, String fileIn) throws IOException {
		BufferedReader reader = new BufferedReader(new FileReader(fileIn));
		String outFile = ext.rootOf(fileIn, false) + ".ids";
		PrintWriter writer = new PrintWriter(new FileWriter(outFile));

		SampleData sampleData = proj.getSampleData(0, false);

		while (reader.ready()) {
			String line = reader.readLine();
			String[] keyVal = sampleData.lookup(line);
			writer.println(ArrayUtils.toStr(keyVal, "\t"));
		}
		writer.flush();
		writer.close();
		reader.close();
	}

	public static void filterCentromeric(String dir, String in, String out,
																			 String markerSetFilenameToBreakUpCentromeres, int build,
																			 Logger log) {
		PrintWriter writer;
		String[] line;
		CNVariant cnv;
		Segment[] centromereMidpoints;
		int[][] centromereBoundaries;
		BufferedReader reader = null;
		FileReader fr = null;

		centromereBoundaries = Positions.determineCentromereBoundariesFromMarkerSet(markerSetFilenameToBreakUpCentromeres,
																																								build, log);
		centromereMidpoints = Positions.computeCentromereMidpoints(centromereBoundaries);

		try {
			fr = new FileReader(dir + in);
			reader = new BufferedReader(fr);
			writer = new PrintWriter(new FileWriter(dir + out));
			writer.println(reader.readLine());
			while (reader.ready()) {
				line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
				cnv = new CNVariant(line);
				if (cnv.overlaps(centromereMidpoints[cnv.getChr()])) {
					writer.println(ArrayUtils.toStr(line));
				}
			}
			fr.close();
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + dir + in + "\" not found in current directory");
			return;
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + dir + in + "\"");
			return;
		} finally {
			if (fr != null) {
				try {
					fr.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}

	public static void breakCentromeric() throws IOException {
		CNVFilter filter = new CNVFilter(null);
		filter.setBreakupCentromeres(true);
		filter.setCentromereBoundariesFromFile("D:/data/gedi_gwas/data/markers.bim");
		filter.computeCentromereMidPoints();

		CNVariant[] centromeric = CNVariant.loadPlinkFile("D:/SIDS and IQ/IQ/merged.cnv", false);

		PrintWriter writer = new PrintWriter(new FileWriter("D:/SIDS and IQ/IQ/merged_split.cnv"));
		writer.println(ArrayUtils.toStr(CNVariant.PLINK_CNV_HEADER, "\t"));

		for (CNVariant cnv : centromeric) {
			CNVFilterPass fp = filter.getCNVFilterPass(cnv);
			if (fp.isCentromeric()) {
				CNVariant[] broken = filter.breakUpCentromere(fp, cnv);
				for (CNVariant newcnv : broken) {
					writer.println(newcnv.toPlinkFormat());
				}
			} else {
				writer.println(cnv.toPlinkFormat());
			}
		}
		writer.flush();
		writer.close();

	}

	private static void breakSeqMetaDataFileIntoChrs() throws IOException {
		String dir = "/panfs/roc/groups/5/pankrat2/shared/skatMeta/snpInfos/exome_chip_v7/";
		String dataFile = "SNPInfo_HumanExome-12v1_rev7b2_slim_wChr.txt";

		String[] outHdr = {"SNP", "CHR", "MapInfo", "sc_exonic", "sc_nonsynSplice", "sc_damaging",
											 "sc_lof", "SKATgene"};

		HashMap<String, PrintWriter> chrWriters = new HashMap<>();
		for (int i = 1; i < 27; i++) {
			PrintWriter writer = Files.getAppropriateWriter(dir + "chr" + i + ".csv");
			writer.println(ArrayUtils.toStr(outHdr, ","));
			chrWriters.put(Positions.CHR_CODES[i], writer);
		}

		BufferedReader reader = Files.getAppropriateReader(dir + dataFile);
		reader.readLine();
		String line;
		while ((line = reader.readLine()) != null) {
			String[] parts = line.split("\t");
			chrWriters.get(parts[1]).println(ArrayUtils.toStr(parts, ","));
		}
		reader.close();

		for (PrintWriter writer : chrWriters.values()) {
			writer.flush();
			writer.close();
		}

	}

	private static void filterNYChoanalCNVs() {
		String filenameOfProblematicRegions = null;
		String individualsToKeepFile = null;
		int commonInOutOrIgnore = FilterCalls.COMMON_IGNORED;
		String markerSetFilenameToBreakUpCentromeres_1 = "/scratch.global/cole0482/ny_choanal/shadow11combo/markerPositions.txt";
		String markerSetFilenameToBreakUpCentromeres_2 = "/scratch.global/cole0482/ny_choanal/shadow12combo/markerPositions.txt";
		int build = 37;
		boolean makeUCSC = false;
		int[] del = new int[] {10, 0};
		int[] dup = new int[] {10, 10};
		int[] number = new int[] {5, 3};
		int score = 10;

		String[][] files = new String[][] { {"/scratch.global/cole0482/ny_choanal/shadow11combo/cnv/",
																				 "23Mgen_merged.cnv", "23_M_filtered.cnv",
																				 markerSetFilenameToBreakUpCentromeres_1},
																			 {"/scratch.global/cole0482/ny_choanal/shadow11combo/cnv/",
																				"23Fgen_merged.cnv", "23_F_filtered.cnv",
																				markerSetFilenameToBreakUpCentromeres_1},
																			 {"/scratch.global/cole0482/ny_choanal/shadow12combo/cnv/",
																				"23Mgen_merged.cnv", "23_M_filtered.cnv",
																				markerSetFilenameToBreakUpCentromeres_2},
																			 {"/scratch.global/cole0482/ny_choanal/shadow12combo/cnv/",
																				"23Fgen_merged.cnv", "23_F_filtered.cnv",
																				markerSetFilenameToBreakUpCentromeres_2},
																			 {"/scratch.global/cole0482/ny_choanal/shadow11combo/cnv/",
																				"24_M_genvisis.cnv", "24_M_filtered.cnv",
																				markerSetFilenameToBreakUpCentromeres_1}};

		for (String[] fileSet : files) {
			FilterCalls.filterCNVs(fileSet[0], fileSet[1], fileSet[2], del, dup, number, score,
														 filenameOfProblematicRegions, commonInOutOrIgnore,
														 individualsToKeepFile, true, fileSet[3], makeUCSC, build,
														 new Logger());
		}
	}


	private static class Alleles {
		String a1;
		String a2;

		public Alleles(String a1, String a2) {
			this.a1 = a1;
			this.a2 = a2;
		}
	}

	public static void famLookup() {
		String famFile = "/home/pankrat2/shared/aric_gw6/ARICGenvisis_CEL_FULL/plinkApril2017/plink.fam";
		String famFileOut = "/home/pankrat2/shared/aric_gw6/ARICGenvisis_CEL_FULL/plinkApril2017/plink_corrected.fam";
		String pedFile = "/home/pankrat2/shared/aric_gw6/ARICGenvisis_CEL_FULL/pedigree_fixed.dat";
		String[][] ped = HashVec.loadFileToStringMatrix(pedFile, false, null, "\t", false, 100000,
																										false);
		String[][] fam = HashVec.loadFileToStringMatrix(famFile, false, null, "[\\s]+", false, 100000,
																										false);
		String dropSampFile = "/home/pankrat2/shared/aric_gw6/ARICGenvisis_CEL_FULL/plinkApril2017/dropSamples.txt";

		HashMap<String, String> dnaToFIDIID = new HashMap<>();
		for (String[] s : ped) {
			dnaToFIDIID.put(s[6], s[0] + "\t" + s[1]);
		}

		String ids;
		StringBuilder sb;
		PrintWriter writer = Files.getAppropriateWriter(famFileOut);
		PrintWriter drops = Files.getAppropriateWriter(dropSampFile);
		for (String[] line : fam) {
			sb = new StringBuilder();
			ids = dnaToFIDIID.get(line[0]);
			if (ids == null) {
				sb.append(line[0]).append("\t").append(line[1]);
				drops.println(line[0]);
			} else {
				sb.append(ids);
			}
			sb.append("\t").append(line[2])
				.append("\t").append(line[3])
				.append("\t").append(line[4])
				.append("\t").append(line[5]);
			writer.println(sb.toString());
		}
		writer.flush();
		writer.close();
		drops.flush();
		drops.close();
	}


	public static void crossRefAndFilter(String f1, String f2, String fOut, boolean not) {
		HashSet<String> f1Data = HashVec.loadFileToHashSet(f1, false);
		HashSet<String> f2Data = HashVec.loadFileToHashSet(f2, false);

		HashSet<String> left = new HashSet<>(f1Data);

		if (not) {
			left.removeAll(f2Data);
		} else {
			left.retainAll(f2Data);
		}

		Files.writeIterable(left, fOut);
	}

	public static void affy6SnpLookup(String file) {
		String affySnpFile = "/home/pankrat2/cole0482/Affy6_SnpList.xln";
		String[][] aff = HashVec.loadFileToStringMatrix(affySnpFile, true, null, "[\\s]+", false,
																										100000, false);
		HashMap<String, String> affRS = new HashMap<>();
		for (String[] line : aff) {
			affRS.put(line[0], line[1]);
		}

		String out = ext.rootOf(file, false) + "_corrected.txt";
		PrintWriter writer = Files.getAppropriateWriter(out);
		String[] mkrs = HashVec.loadFileToStringArray(file, false, new int[] {0}, false);
		int missCnt = 0;
		for (String snp : mkrs) {
			String rs = affRS.get(snp);
			if (rs == null || "---".equals(rs))
				missCnt++;
			writer.println(rs == null || "---".equals(rs) ? snp : rs);
		}
		writer.flush();
		writer.close();

		System.out.println(missCnt + " snps missing an RS number in file: " + file);

	}

	private static HashMap<String, String> loadCallrates(String callrateFile) {
		String[] hdr = Files.getHeaderOfFile(callrateFile, new Logger());
		int[] mkrInds = ext.indexFactors(Aliases.MARKER_NAMES, hdr, false, false);
		int mkrInd = -1;
		for (int m : mkrInds) {
			if (m >= 0) {
				mkrInd = m;
				break;
			}
		}
		int callInd = ext.indexOfStr("CallRate", hdr, false, true);

		if (mkrInd == -1 || callInd == -1) {
			System.err.println("Error - Couldn't find header.");
			return null;
		}

		HashMap<String, String> callrateMap = new HashMap<>();

		String[][] info = HashVec.loadFileToStringMatrix(callrateFile, true,
																										 new int[] {mkrInd, callInd}, false);
		for (String[] line : info) {
			callrateMap.put(line[0], line[1]);
		}

		info = null;
		return callrateMap;
	}

	static class MarkerLookup {
		private final HashMap<String, ArrayList<String>> markerLookup;
		private final HashMap<String, String> reverseLookup;

		public MarkerLookup() {
			markerLookup = new HashMap<>();
			reverseLookup = new HashMap<>();
		}

		public void addMarker(String[] line) {
			ArrayList<String> mkrs = markerLookup.get(line[0]);
			if (mkrs == null) {
				mkrs = new ArrayList<>();
				markerLookup.put(line[0], mkrs);
			}
			for (int i = 1; i < line.length; i++) {
				mkrs.add(line[i]);
				reverseLookup.put(line[i], line[0]);
			}
		}

		public Set<String> getNewMarkers() {
			return markerLookup.keySet();
		}

		public ArrayList<String> getOldMarkers(String s) {
			return markerLookup.get(s);
		}

	}

	public static MarkerLookup loadMarkerInfoFile(String mkrInfoFile) throws IOException {
		BufferedReader reader = Files.getAppropriateReader(mkrInfoFile);
		String line;
		MarkerLookup ml = new MarkerLookup();
		while ((line = reader.readLine()) != null) {
			String[] pts = line.split("[\\s]+");
			ml.addMarker(pts);
		}
		reader.close();
		return ml;
	}

	/**
	 * 
	 * @param mkrInfoFile Duplicates file. Expects the following format:<br/>
	 *        <code>NewMarkerName &lt;...markers to combine into this marker&gt;</code>
	 * @throws IOException
	 */
	public static void markerDuplicateFilter(String mkrInfoFile, String missDropsFile,
																					 String callrateFile, String outFile) throws IOException {

		HashMap<String, String> callrates = loadCallrates(callrateFile);
		if (callrates == null)
			return;
		String[] missDropsArr = HashVec.loadFileToStringArray(missDropsFile, false, null, false);
		HashSet<String> missDrops = new HashSet<>();
		for (String m : missDropsArr) {
			missDrops.add(m);
		}

		MarkerLookup ml = loadMarkerInfoFile(mkrInfoFile);

		HashSet<String> allDrops = new HashSet<>();

		for (String s : ml.getNewMarkers()) {
			ArrayList<String> oldMkrs = ml.getOldMarkers(s);
			if (oldMkrs.size() == 0) {
				continue;
			} else if (oldMkrs.size() == 1) {
				if (missDrops.contains(oldMkrs.get(0))) {
					allDrops.add(oldMkrs.get(0));
				}
			} else {
				boolean droppedAll = true;
				for (String o : oldMkrs) {
					if (!missDrops.contains(o)) {
						droppedAll = false;
						break;
					}
				}
				if (!droppedAll) {
					// if droppedAll == true, all markers will be added through miss_drops, so no worry
					ArrayList<Double> mkrCalls = new ArrayList<>();
					for (String o : oldMkrs) {
						if (missDrops.contains(o)) {
							mkrCalls.add(-1d);
						} else {
							String callStr = callrates.get(o);
							if (callStr == null) {
								System.err.println("Error - no callrate info for marker: " + o);
								mkrCalls.add(-1d);
							} else {
								mkrCalls.add(Double.parseDouble(callStr));
							}
						}
					}
					double maxCall = 0;
					int maxInd = -1;
					for (int i = 0; i < mkrCalls.size(); i++) {
						if (mkrCalls.get(i) > maxCall) {
							maxInd = i;
							maxCall = mkrCalls.get(i);
						}
					}
					for (int i = 0; i < mkrCalls.size(); i++) {
						if (maxInd == i)
							continue;
						allDrops.add(oldMkrs.get(i));
					}
				}
			}
		}

		Files.writeIterable(allDrops, outFile);
	}


	public static void affy6BimLookup() {
		String bimFile = "/home/pankrat2/shared/aric_gw6/ARICGenvisis_CEL_FULL/plinkApril2017/ancestryPipeline/plink.bim";
		String newBimFile = "/home/pankrat2/shared/aric_gw6/ARICGenvisis_CEL_FULL/plinkApril2017/ancestryPipeline/plink_correctedRS.bim";
		String missSnpFile = "/home/pankrat2/shared/aric_gw6/ARICGenvisis_CEL_FULL/plinkApril2017/ancestryPipeline/plink_missingRS.txt";
		String mismatchFile = "/home/pankrat2/shared/aric_gw6/ARICGenvisis_CEL_FULL/plinkApril2017/ancestryPipeline/plink_mismatchAlleles.txt";
		String affySnpFile = "/home/pankrat2/cole0482/Affy6_SnpList.xln";
		String[][] bim = HashVec.loadFileToStringMatrix(bimFile, false, null, "\t", false, 100000,
																										false);
		String[][] aff = HashVec.loadFileToStringMatrix(affySnpFile, true, null, "[\\s]+", false,
																										100000, false);

		System.out.println("Loaded data...");
		System.out.println(bim.length + " lines in .bim file;");
		System.out.println(aff.length + " lines in snp lookup file;");

		HashMap<String, String> affRS = new HashMap<>();
		HashMap<String, Alleles> affMkrs = new HashMap<>();
		for (String[] line : aff) {
			affRS.put(line[0], line[1]);
			affMkrs.put(line[0], new Alleles(line[3], line[4]));
		}

		PrintWriter writer = Files.getAppropriateWriter(newBimFile);
		PrintWriter miss = Files.getAppropriateWriter(missSnpFile);
		PrintWriter mismatch = Files.getAppropriateWriter(mismatchFile);
		mismatch.println("AFFY\tRSID\tA_A1\tA_A2\tRS_A1\tRS_A2");
		for (String[] line : bim) {
			String snp = line[1];
			String rs = affRS.get(snp);
			if (rs == null) {
				miss.println(snp);
				writer.println(ArrayUtils.toStr(line, "\t"));
			} else {
				Alleles all = affMkrs.get(snp);
				if ((all.a1.equals(line[4]) && all.a2.equals(line[5]))
						|| (all.a2.equals(line[4]) && all.a1.equals(line[5]))) {
					if (!rs.equals("---")) {
						line[1] = rs;
					}
					// } else {
					// mismatch.println(snp + "\t" + rs + "\t" + all.a1 + "\t" + all.a2 + "\t" + line[4] +
					// "\t"
					// + line[5]);
				}
				writer.println(ArrayUtils.toStr(line, "\t"));
			}
		}
		writer.flush();
		writer.close();
		miss.flush();
		miss.close();
		mismatch.flush();
		mismatch.close();

		System.out.println("Done!");
	}


	public static void genDupe() {
		String file = "Affy6_SnpList.xln";
		String[][] data = HashVec.loadFileToStringMatrix(file, false, null, "[\\s]+", false, 100000,
																										 false);
		String out = "Affy6_duplicates.txt";
		HashMap<String, ArrayList<String>> map = new HashMap<>();
		for (String[] line : data) {
			if (line[1].equals("---"))
				continue;
			ArrayList<String> list = map.get(line[1]);
			if (list == null) {
				list = new ArrayList<>();
				map.put(line[1], list);
			}
			list.add(line[0]);
		}

		StringBuilder sb;
		PrintWriter writer = Files.getAppropriateWriter(out);
		for (Entry<String, ArrayList<String>> entry : map.entrySet()) {
			if (entry.getValue().size() == 1) {
				continue;
			}
			sb = new StringBuilder();
			sb.append(entry.getKey());
			for (String s : entry.getValue()) {
				sb.append("\t").append(s);
			}
			writer.println(sb.toString());
		}
		writer.flush();
		writer.close();
	}


	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		Project proj;
		String filename = "lab.dat";
		String logfile = null;
		String file = null;

		boolean test = true;
		if (test) {

			// String dir = "/home/pankrat2/shared/aric_gw6/ARICGenvisis_CEL_FULL/plinkApril2017/";
			// String mkrInfoFile = "~/Affy6_duplicates.txt";
			// String missDropsFile = dir + "quality_control/further_analysis_QC/miss_drops.dat";
			// String callrateFile =
			// "/home/pankrat2/shared/aric_gw6/ARICGenvisis_CEL_FULL/ohw_ws_20_ALL1000PCs_gc_corrected_OnTheFly_LRR_035CR_096_SAMP_LRR_TRUE_markerQC.txt";
			// String outFile = dir + "markerDuplicates.out";
			// markerDuplicateFilter(mkrInfoFile, missDropsFile, callrateFile, outFile);


			// genDupe();

			// crossRefAndFilter("/scratch.global/cole0482/testImp/snps/allMarkers.txt",
			// "/scratch.global/cole0482/testImp/snps/miss_drops.dat",
			// "/scratch.global/cole0482/testImp/snps/cleanMarkers.txt", true);
			//
			// affy6SnpLookup("/scratch.global/cole0482/testImp/snps/cleanMarkers.txt");

			// crossRefAndFilter("/scratch.global/cole0482/testImp/snps/allMarkers_corrected.txt",
			// "/scratch.global/cole0482/testImp/snps/miss_drops_corrected.txt",
			// "/scratch.global/cole0482/testImp/snps/cleanMarkers_corrected.txt", true);

			// famLookup();
			affy6BimLookup();

			// String cmd =
			// "java -jar genvisis.jar org.genvisis.imputation.ImputationPipeline"
			// + " proj=projects/poynter.properties"
			// + " ref=/home/pankrat2/shared/bin/ref/1000GP_Phase3_combined.legend.gz"
			// + " chrs=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22"
			// + " plinkDir=plink/"
			// + " outDir=/scratch.global/cole0482/testImp/out/"
			// + " type=PLINK_SHAPE_MINI"
			// + " " ;
			// System.out.println(cmd);


			// proj = new Project("projects/poynter.properties", false);
			// String referenceFile = "/home/pankrat2/shared/bin/ref/1000GP_Phase3_combined.legend.gz";
			// ImputationPipeline ip = new ImputationPipeline(proj, referenceFile);
			// ip.loadDefaultDropFiles(proj.PROJECT_DIRECTORY.getValue() + "plink/");
			// ip.exportToVCF("/scratch.global/cole0482/testImp/output");
			// ip.exportToPlink("/scratch.global/cole0482/testImp/plink");
			// String hapsDir = "/scratch.global/cole0482/testImp/out/";
			// String outDir = "/scratch.global/cole0482/testImp/out/min/";
			//
			// new ImputationImpl.ShapeIt(proj, "/scratch.global/cole0482/testImp/", "plink_chr",
			// hapsDir).createScripts();
			// new ImputationImpl.MiniMac(proj, hapsDir, outDir).createScripts();

			// System.out.println("Username: " + QueueControl.getUserName());
			// System.out.println("Group: " + QueueControl.getCurrentGroup());
			// System.out.println("All Groups: " + QueueControl.getUserGroups().toString());
			// System.out.println();
			// System.out.println("Default Queue: " +
			// QueueControl.findSensibleDefault(QueueControl.parseAllowedQueues(new Logger())).getName());
			//
			// String dir = "F:/temp/filter/";
			// String in = "recodedM.cnv";
			// String out = "recodedM_excl.cnv";
			// String indivFile = "F:/temp/filter/exclude.txt";
			// boolean exclude = true;
			// FilterCalls.filterExclusions(dir, in, out, indivFile, exclude);
			// in = "recodedF.cnv";
			// out = "recodedF_excl.cnv";
			// FilterCalls.filterExclusions(dir, in, out, indivFile, exclude);
			//
			// System.out.println("Removed excluded");
			//
			// CNVFilter cnvF = new CNVFilter(new Logger());
			// cnvF.setProblemRegionsFromFile("F:/temp/filter/problematicRegions_hg19.dat");
			//
			// in = dir + "recodedM_excl.cnv";
			// out = dir + "recodedM_filt.cnv";
			// CNVFilter.filterCNVs(in, out, cnvF, new Logger());
			//
			// in = dir + "recodedF_excl.cnv";
			// out = dir + "recodedF_filt.cnv";
			// CNVFilter.filterCNVs(in, out, cnvF, new Logger());

			// MergeExtractPipeline pipeline = new MergeExtractPipeline();
			// // pipeline.setMarkers(markersFile);
			// pipeline.setRunDirectory("/scratch.global/cole0482/merge/", true);
			// pipeline.setOutputFormat(DosageData.DATABASE_DOSE_FORMAT);
			// pipeline.setOutputFiles(outFile, mapOutFile);
			// pipeline.setRenameMarkers(true);
			// // pipeline.addDataSource("/scratch.global/cole0482/merge/blacks/", "gwas.bed", "gwas.bim",
			// // "gwas.fam");
			// pipeline.addDataSource( "exome", "/scratch.global/cole0482/merge/blacks/", "exome.bed",
			// "exome.bim", "exome.fam");
			// pipeline.addDataSource( "metab", "/scratch.global/cole0482/merge/blacks/", "metab.bed",
			// "metab.bim", "metab.fam");
			// // add more;
			// pipeline.run();


			// String doseFile1 =
			// "/home/pankarne/shared/ARIC_Genomics_Data/GWAS_Chip/1000G/ARIC.whites.impute2/chr3.90069244.95069244.impute2.gz";
			// String mapFile1 =
			// "/home/pankarne/shared/ARIC_Genomics_Data/GWAS_Chip/1000G/ARIC.whites.impute2/chr3.90069244.95069244.impute2_info";
			//
			// String doseFile2 =
			// "/home/pankarne/shared/ARIC_Genomics_Data/GWAS_Chip/1000G/ARIC.whites.impute2/chr3.95069244.100069244.impute2.gz";
			// String mapFile2 =
			// "/home/pankarne/shared/ARIC_Genomics_Data/GWAS_Chip/1000G/ARIC.whites.impute2/chr3.95069244.100069244.impute2_info";
			//
			// String idFile = "/home/pankarne/cole0482/EA.indiv.dup";
			// String outFile = "/scratch.global/cole0482/test.db.xln.gz";
			// String mapOutFile = "/scratch.global/cole0482/mapOut.xln";
			//
			// DosageData dd1 = new DosageData(doseFile1, idFile, mapFile1,
			// DosageData.IMPUTE2_DOSE_FORMAT, null, true, null);
			// DosageData dd2 = new DosageData(doseFile2, idFile, mapFile2,
			// DosageData.IMPUTE2_DOSE_FORMAT, null, true, null);
			// DosageData dd3 = DosageData.combine(dd1, dd2);
			// dd1 = null;
			// dd2 = null;
			// dd1 = DosageData.loadPlinkBinary(dir, plinkRoot);
			// dd2 = DosageData.combine(dd3, dd1);
			// dd1 = null;
			// dd3 = null;
			// dd1 = DosageData.loadPlinkBinary(dir2, plinkRoot2);
			// dd3 = DosageData.combine(dd2, dd1);
			// dd3.writeToFile(outFile, mapOutFile, null, DosageData.DATABASE_DOSE_FORMAT, null);
			// System.out.println("complete!");


			return;
		}


		String usage = "";

		if (numArgs == 0) {
			try {
				countCNVsForIndividuals("D:/data/ny_registry/new_york/stats/puv_ids.txt",
																"D:/data/ny_registry/new_york/stats/recodedM.cnv",
																"D:/data/ny_registry/new_york/stats/puv_cnvs.cnv");
				// testClipboard();
				// BufferedReader reader = new BufferedReader(new
				// FileReader("D:/ForestPlot/Hb_SingleSNP.csv"));
				// String line = reader.readLine();
				// do {
				// System.out.println(line);
				// } while((line = reader.readLine()) != null);

				// filterForMarkers("D:/height/scratch/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq_parsed.xln",
				// "D:/height/scratch/samples_logan.bim");
				// writeRegions();
				// stripMarkerNames();
				// writeExclusions();
				// concatFiles();
				// splitFile();
				// idSwap(new Project("D:/projects/gedi_gwas.properties", false),
				// "D:/data/gedi_gwas/overlap_ids.txt");
				// compareMarkers();
				// filterLDFiles(0.5);
				// formatLDResults();
				// filter();
				// breakCentromeric();
				// filterWrong();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			// mockupGUI();
			return;
		}


		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("proj=")) {
				filename = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("log=")) {
				logfile = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("file=")) {
				file = arg.split("=")[1];
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
			proj = new Project(filename, logfile, false);
			if (file != null) {
				idSwap(proj, file);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
