package org.genvisis.one.ben;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;

import org.genvisis.cnv.analysis.FilterCalls;
import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.manage.MDL;
import org.genvisis.cnv.manage.TransposeData;
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
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.ext;
import org.genvisis.filesys.CNVariant;
import org.genvisis.filesys.DosageData;
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

		List<CNVariant> cnvs = CNVariant.loadPlinkFile(cnvFile, sampleKeyHash, true);

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
		PrintWriter writer = Files.openAppropriateWriter(outFile);

		SampleData sampleData = proj.getSampleData(false);

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
			writer = Files.openAppropriateWriter(dir + out);
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

		CNVariant[] centromeric = CNVariant.loadPlinkFile("D:/SIDS and IQ/IQ/merged.cnv");

		PrintWriter writer = Files.openAppropriateWriter("D:/SIDS and IQ/IQ/merged_split.cnv");
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

		String[][] files = new String[][] {{"/scratch.global/cole0482/ny_choanal/shadow11combo/cnv/",
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

	public static void famRecode() {
		String famFile = "/home/pankrat2/shared/aric_gw6/ARICGenvisis_CEL_FULL/plinkDNA/final/plink.fam";
		String famFileOut = "/home/pankrat2/shared/aric_gw6/ARICGenvisis_CEL_FULL/plinkDNA/final/plink_corrected.fam";
		String lookupFile = "/home/pankrat2/shared/aric_gw6/ARICGenvisis_CEL_FULL/plinkDNA/final/idLookup.txt";

		String[][] ped = HashVec.loadFileToStringMatrix(lookupFile, false, null, "\t", 100000, false);
		String[][] fam = HashVec.loadFileToStringMatrix(famFile, false, null, "[\\s]+", 100000, false);

		HashMap<String, String> lookup = new HashMap<>();
		for (String[] s : ped) {
			lookup.put(s[0], s[1]);
		}

		String id;
		PrintWriter writer = Files.getAppropriateWriter(famFileOut);
		for (String[] line : fam) {
			id = lookup.get(line[0]);
			line[0] = id;
			line[1] = id;
			writer.println(ArrayUtils.toStr(line, "\t"));
		}
		writer.flush();
		writer.close();

	}


	public static void famLookup() {
		String famFile = "/home/pankrat2/shared/aric_gw6/ARICGenvisis_CEL_FULL/plinkApril2017/plink.fam";
		String famFileOut = "/home/pankrat2/shared/aric_gw6/ARICGenvisis_CEL_FULL/plinkApril2017/plink_corrected.fam";
		String pedFile = "/home/pankrat2/shared/aric_gw6/ARICGenvisis_CEL_FULL/pedigree_fixed.dat";
		String[][] ped = HashVec.loadFileToStringMatrix(pedFile, false, null, "\t", 100000, false);
		String[][] fam = HashVec.loadFileToStringMatrix(famFile, false, null, "[\\s]+", 100000, false);
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
		String[][] aff = HashVec.loadFileToStringMatrix(affySnpFile, true, null, "[\\s]+", 100000,
																										false);
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
		int[] mkrInds = ext.indexFactors(Aliases.MARKER_NAMES, hdr, false);
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
																										 new int[] {mkrInd, callInd});
		for (String[] line : info) {
			callrateMap.put(line[0], line[1]);
		}

		info = null;
		return callrateMap;
	}

	public static void affy6BimLookup() {
		// String bimFile =
		// "/home/pankrat2/shared/aric_gw6/ARICGenvisis_CEL_FULL/plinkApril2017/compare/plinkDropFilt.bim";
		// String newBimFile =
		// "/home/pankrat2/shared/aric_gw6/ARICGenvisis_CEL_FULL/plinkApril2017/compare/plinkDropFilt_correctedRS.bim";
		// String missSnpFile =
		// "/home/pankrat2/shared/aric_gw6/ARICGenvisis_CEL_FULL/plinkApril2017/compare/plinkDropFilt_missingRS.txt";
		// String mismatchFile =
		// "/home/pankrat2/shared/aric_gw6/ARICGenvisis_CEL_FULL/plinkApril2017/compare/plinkDropFilt_mismatchAlleles.txt";
		String bimFile = "/home/pankrat2/shared/aric_gw6/ARICGenvisis_CEL_FULL/plinkDNA/final?/plinkNoRSDupe.bim";
		String newBimFile = "/home/pankrat2/shared/aric_gw6/ARICGenvisis_CEL_FULL/plinkDNA/final?/plinkNoRSDupe_correctedRS.bim";
		String missSnpFile = "/home/pankrat2/shared/aric_gw6/ARICGenvisis_CEL_FULL/plinkDNA/final?/plinkNoRSDupe_missingRS.txt";
		String mismatchFile = "/home/pankrat2/shared/aric_gw6/ARICGenvisis_CEL_FULL/plinkDNA/final?/plinkNoRSDupe_mismatchAlleles.txt";
		// String bimFile =
		// "/home/pankrat2/shared/aric_gw6/ARICGenvisis_CEL_FULL/plinkApril2017/ancestryPipeline/plink.bim";
		// String newBimFile =
		// "/home/pankrat2/shared/aric_gw6/ARICGenvisis_CEL_FULL/plinkApril2017/ancestryPipeline/plink_correctedRS.bim";
		// String missSnpFile =
		// "/home/pankrat2/shared/aric_gw6/ARICGenvisis_CEL_FULL/plinkApril2017/ancestryPipeline/plink_missingRS.txt";
		// String mismatchFile =
		// "/home/pankrat2/shared/aric_gw6/ARICGenvisis_CEL_FULL/plinkApril2017/ancestryPipeline/plink_mismatchAlleles.txt";
		// String bimFile = "/scratch.global/cole0482/affy6plink/imputeSrc_EA.bim";
		// String newBimFile = "/scratch.global/cole0482/affy6plink/imputeSrc_EA_corrected.bim";
		// String missSnpFile = "/scratch.global/cole0482/affy6plink/imputeSrc_EA_miss.txt";
		// String mismatchFile = "/scratch.global/cole0482/affy6plink/imputeSrc_EA_mism.txt";
		String affySnpFile = "/home/pankrat2/cole0482/Affy6_SnpList.xln";
		String[][] bim = HashVec.loadFileToStringMatrix(bimFile, false, null, "\t", 100000, false);
		String[][] aff = HashVec.loadFileToStringMatrix(affySnpFile, true, null, "[\\s]+", 100000,
																										false);

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
				// if ((all.a1.equals(line[4]) && all.a2.equals(line[5]))
				// || (all.a2.equals(line[4]) && all.a1.equals(line[5]))) {
				if (!rs.equals("---")) {
					line[1] = rs;
				}
				// } else {
				// mismatch.println(snp + "\t" + rs + "\t" + all.a1 + "\t" + all.a2 + "\t" + line[4] +
				// "\t"
				// + line[5]);
				// }
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
		String[][] data = HashVec.loadFileToStringMatrix(file, false, null, "[\\s]+", 100000, false);
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

	public static void exomeRecode(String bimFile, String newBimFile) {
		// String dir = "/scratch.global/cole0482/affy6plink/back/";
		String dir = "D:/temp/plink/exmToRS/";
		String exomeLookup = dir + "exm_to_rsID_lookup.txt";
		// String bimFile = dir + "exome_EA.bim";
		// String newBimFile = dir + "exome_EA_corrected.bim";
		String[][] exmMkrs = HashVec.loadFileToStringMatrix(exomeLookup, true, null, "\t", 10000,
																												false);
		HashMap<String, String> lookup = new HashMap<>();
		for (String[] mkrs : exmMkrs) {
			if (!".".equals(mkrs[1])) {
				if (lookup.containsKey(mkrs[0]) && !lookup.get(mkrs[0]).equals(mkrs[1])) {
					System.out.println("Duplicate entry: " + mkrs[0] + " -> " + mkrs[1] + " | " + mkrs[0]
														 + " -> " + lookup.get(mkrs[0]));
				}
				lookup.put(mkrs[0], mkrs[1]);
			}
		}

		String[][] bimData = HashVec.loadFileToStringMatrix(bimFile, false, null, "\t", 100000, false);
		PrintWriter writer = Files.getAppropriateWriter(newBimFile);
		for (String[] line : bimData) {
			String mkr = line[1];
			if (lookup.containsKey(mkr)) {
				line[1] = lookup.get(mkr);
			}
			writer.println(ArrayUtils.toStr(line, "\t"));
		}
		writer.flush();
		writer.close();
	}

	public static void filt() {
		String dir = "/home/pankrat2/shared/aric_gw6/ARICGenvisis_CEL_FULL/plinkApril2017/compare/";
		String f2 = dir + "cleanSnps.txt";
		String f1 = dir + "bimSnps.txt";

		String[] f11 = HashVec.loadFileToStringArray(f1, false, null, false);
		String[] f21 = HashVec.loadFileToStringArray(f2, false, null, false);

		HashSet<String> h1 = new HashSet<>();
		for (String s : f11) {
			h1.add(s);
		}

		String o1 = dir + "cleanNotBim.snp";
		PrintWriter writer = Files.getAppropriateWriter(o1);
		for (String s : f21) {
			if (!h1.contains(s)) {
				writer.println(s);
			}
		}
		writer.flush();
		writer.close();
	}

	public static void addPFBToMarkerMetrics() throws IOException {
		String dir = "/scratch.global/cole0482/ny_prev/";
		String fil = "marker_lrr_sd.xln";
		String out = "marker_lrr_sd_pfb.xln";
		String pfb1 = "data/females.pfb";
		String pfb2 = "data/males.pfb";
		String pfb3 = "data/custom.pfb";

		Hashtable<String, String> pfb1Map = HashVec.loadFileToHashString(dir + pfb1, 0, new int[] {3},
																																		 "\t", true);
		Hashtable<String, String> pfb2Map = HashVec.loadFileToHashString(dir + pfb2, 0, new int[] {3},
																																		 "\t", true);
		Hashtable<String, String> pfb3Map = HashVec.loadFileToHashString(dir + pfb3, 0, new int[] {3},
																																		 "\t", true);

		BufferedReader reader = Files.getAppropriateReader(dir + fil);
		PrintWriter writer = Files.getAppropriateWriter(dir + out);
		String line = null;
		String[] pts = null;
		int cnt = 0;
		while ((line = reader.readLine()) != null) {
			pts = line.split("\t", -1);
			pts = ArrayUtils.addStrToArray(cnt == 1 ? pfb1Map.get(pts[0]) : "Female PFB", pts);
			pts = ArrayUtils.addStrToArray(cnt == 1 ? pfb2Map.get(pts[0]) : "Male PFB", pts);
			pts = ArrayUtils.addStrToArray(cnt == 1 ? pfb3Map.get(pts[0]) : "Custom PFB", pts);
			if (cnt == 0) {
				cnt = 1;
			}
			writer.println(ArrayUtils.toStr(pts, "\t"));
		}
		reader.close();
		writer.flush();
		writer.close();
	}

	private static void runMarcotte() {

		String dir = "F:/temp/HB_PLINK/dupeSets/";
		String file = "Marcotte_dupeSet1.ped.in";
		String out = "Marcotte_dupeSet1.ped";

		String[][] data = HashVec.loadFileToStringMatrix(dir + file, false, null, "\t", 3000, false);

		PrintWriter writer = Files.getAppropriateWriter(dir + out);
		for (String[] line : data) {
			StringBuilder sb = new StringBuilder();
			for (int i = 0; i < 6; i++) {
				sb.append(line[i]).append("\t");
			}
			for (int i = 6; i < line.length; i++) {
				sb.append(line[i].charAt(0)).append("\t").append(line[i].charAt(1)).append("\t");
			}
			writer.println(sb.toString());
		}
		writer.flush();
		writer.close();

	}

	private static void run() throws IOException {
		String dir = "F:/testProjectSrc/UKBB_AffyAxiom/00src/";
		String file = "ukb_baf_chr21_v2.txt";

		BufferedReader reader = Files.getAppropriateReader(dir + file);
		String line = null;
		line = reader.readLine();
		int len = line.length();
		reader.close();

		reader = Files.getAppropriateReader(dir + file);
		line = null;
		System.out.println("1st line: " + (len + 1));
		reader.skip(len + 1);
		line = reader.readLine();
		len = line.length();
		System.out.println("2nd line: " + line.length());
		System.out.println("|" + line.substring(0, 10) + "|");
		reader.close();


		InputStreamReader isr = Files.getAppropriateInputStreamReader(dir + file);
		int chr = Integer.MIN_VALUE;
		char[] v;
		long total = 0L;
		int ln = 0;
		while ((chr = isr.read()) != -1) {
			v = Character.toChars(chr);
			for (int i = 0; i < v.length; i++) {
				total++;
				if (v[i] == '\n') {
					ln++;
					System.out.println("line " + ln + " starts @ " + total);
				}
			}
		}
		isr.close();


		reader = Files.getAppropriateReader(dir + file);
		line = null;
		reader.skip(total);
		line = reader.readLine();
		len = line.length();
		System.out.println("2nd line: " + line.length());
		System.out.println("|" + line.substring(0, 10) + "|");
		reader.close();


	}

	private static void testBatching() {
		Set<String> complete;
		String[] listAllSamplesInProj;
		int batchMax;
		Logger log = new Logger();

		complete = new HashSet<>();
		listAllSamplesInProj = new String[] {
																				 "0",
																				 "1",
																				 "2",
																				 "3",
																				 "4",
																				 "5",
																				 "6",
																				 "7",
																				 "8",
																				 "9",
																				 "10",
																				 "11",
																				 "12",
																				 "13",
																				 "14",
																				 "15",
																				 "16",
																				 "17",
																				 "18",
																				 "19",
																				 "20",
																				 "21",
																				 "22",
																				 "23",
		};
		complete.add(listAllSamplesInProj[5]);
		complete.add(listAllSamplesInProj[14]);
		complete.add(listAllSamplesInProj[15]);
		// complete.add(listAllSamplesInProj[16]);
		batchMax = 4;

		int[][] ranges = ArrayUtils.splitUpIntoBinsOfIndices(listAllSamplesInProj, complete,
																												 batchMax, log);
		System.out.println("");
		for (int[] batch : ranges) {
			int batchRange = batch[batch.length - 1] - batch[0] + 1;
			boolean range = batchRange == batchMax || batchRange == batch.length;
			System.out.println("Contigu: " + range + " | " + ArrayUtils.toStr(batch, ", "));
		}
	}

	private static void testRev() {
		Project proj = new Project("D:/projects/FarrarReparse.properties");
		String dir = proj.SAMPLE_DIRECTORY.getValue();
		String[] files = Files.list(dir, Sample.SAMPLE_FILE_EXTENSION);

		for (String f : files) {
			try {
				Sample.loadOutOfRangeValuesFromRandomAccessFile(dir + f);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}

	}

	private static void testRevTran() {
		Project proj = new Project("D:/projects/FarrarReparse.properties");
		TransposeData.reverseTranspose(proj);
		System.out.println("TEST");
		testRev();
	}

	private static void filterFCSList() {
		String filename = "F:/Flow/test3/filter/list.txt";
		String file = "F:/Flow/test3/filter/already.txt";

		String[] samp = HashVec.loadFileToStringArray(filename, false, null, false);
		String[] fnum = HashVec.loadFileToStringArray(file, false, null, false);
		HashSet<String> used = new HashSet<>();
		for (String f : fnum) {
			used.add(f);
		}

		PrintWriter writer = Files.getAppropriateWriter("F:/Flow/test3/filter/listOut.xln");
		for (String s : samp) {
			String temp = s.replace("./", "").replace("_", " ").replace("-", " ").replace("/", " ");
			String[] pts = temp.split(" ");
			String date = null;
			String fn = null;
			for (int i = 0; i < pts.length; i++) {
				if (date == null && pts[i].startsWith("201")) {
					date = pts[i] + "-" + pts[i + 1] + "-" + pts[i + 2];
				}
				if (fn == null && pts[i].startsWith("F")) {
					if (pts[i].length() > 1) {
						if (Character.isDigit(pts[i].charAt(2))) {
							fn = pts[i];
						}
					}
				}
			}
			if (fn != null) {
				writer.println(date + "\t" + fn + "\t"
											 + (used.contains(fn) ? "1" : "0") + "\t" + s);
			}
		}

		writer.flush();
		writer.close();
	}

	private static void vcfTest(String filename, String probAttrName) {
		String dir = ext.parseDirectoryOfFile(filename);
		DosageData dd1 = new DosageData(filename, null, null, false, new Logger());
		DosageData dd2 = DosageData.loadVCF(filename, null, null, null);
		dd1.writeToFile(dir + "pd_var1.db.xln.gz", dir + "pd_var1.xln.gz", null, null, new Logger());
		dd2.writeToFile(dir + "pd_var2.db.xln.gz", dir + "pd_var2.xln.gz", null, null, new Logger());
		System.out.println("Done");
	}

	private static void testUKBBMarkerOutliers() {
		Project proj = new Project("/home/pankrat2/cole0482/projects/UKBioBank.properties");
		String[] files = new File(proj.MARKER_DATA_DIRECTORY.getValue()).list();
		int count = 5;
		Random rand = new Random();
		for (int i = 0; i < files.length && count > 0; i++) {
			String fil = files[rand.nextInt(files.length)];
			if (fil.endsWith(".mdRAF")) {
				count = 0;
				if (!fil.startsWith(proj.MARKER_DATA_DIRECTORY.getValue())) {
					fil = ext.verifyDirFormat(proj.MARKER_DATA_DIRECTORY.getValue()) + fil;
				}
				Hashtable<String, Float> outliers = TransposeData.loadOutliersFromRAF(fil);
				int k = 20;
				for (Entry<String, Float> ent : outliers.entrySet()) {
					System.out.println(ent.getKey() + " | " + ent.getValue());
					if (k-- == 0) {
						return;
					}
				}
			}
		}
	}

	private static void rebuildMarkerOutliers() throws ClassNotFoundException, IOException {
		long t1 = System.nanoTime();
		Project proj = new Project("/home/pankrat2/cole0482/projects/UKBioBank.properties");
		String mkrDir = ext.verifyDirFormat(proj.MARKER_DATA_DIRECTORY.getValue());
		String[] files = new File(mkrDir).list((File f, String e) -> {
			return e.endsWith(MarkerData.MARKER_DATA_FILE_EXTENSION);
		});
		String[] samples = proj.getSamples();
		Map<String, Integer> mkrProjInds = proj.getMarkerIndices();
		Hashtable<String, Float> allOutliers = new Hashtable<>();
		String fil;
		Hashtable<String, Float> mkrOutliers;
		String[] markers;
		String[] pts;
		int ind;
		String samp;
		String type;
		System.out.println("Prepped in " + ext.getTimeElapsedNanos(t1));
		long t2 = System.nanoTime();
		for (String mdraf : files) {
			t1 = System.nanoTime();
			fil = mdraf.startsWith(mkrDir) ? mdraf : mkrDir + mdraf;
			mkrOutliers = TransposeData.loadOutliersFromRAF(fil);
			markers = TransposeData.loadMarkerNamesFromRAF(fil);
			for (Entry<String, Float> outlier : mkrOutliers.entrySet()) {
				pts = outlier.getKey().split("\t");
				ind = mkrProjInds.get(markers[Integer.parseInt(pts[0])]);
				samp = samples[Integer.parseInt(pts[1])];
				type = pts[2];
				allOutliers.put(ind + "\t" + samp + "\t" + type, outlier.getValue());
			}
			System.out.println("Loaded outliers from " + mdraf + " in " + ext.getTimeElapsedNanos(t1));
		}
		SerializedFiles.writeSerial(allOutliers, proj.MARKER_DATA_DIRECTORY.getValue(true, true)
																						 + "outliers.ser");
		System.out.println("Rebuilt outliers in " + ext.getTimeElapsedNanos(t2));
	}

	private static void dumpSingleMDRAFOutliers() {
		String[] files = {
											"/scratch.global/cole0482/UKBB2/project/transposed/markers.12.19398.20890.mdRAF",
											"/scratch.global/cole0482/UKBB2/project/transposed/markers.7.17681.19154.mdRAF",
											"/scratch.global/cole0482/UKBB2/project/transposed/markers.3.5980.7475.mdRAF",
											"/scratch.global/cole0482/UKBB2/project/transposed/markers.4.7415.8898.mdRAF",
		};
		for (String file : files) {
			Hashtable<String, Float> outliers = TransposeData.loadOutliersFromRAF(file);
			PrintWriter writer = Files.getAppropriateWriter(ext.rootOf(file, false) + "_outliers.xln");
			for (Entry<String, Float> out : outliers.entrySet()) {
				writer.println(out.getKey() + "\t" + out.getValue());
			}
			writer.close();
		}
	}

	private static void dumpXValues(Project proj) {
		PrintWriter writer = Files.getAppropriateWriter("/home/pankrat2/cole0482/"
																										+ proj.PROJECT_NAME.getValue()
																										+ "_Xvalues.txt");

		MDL mdl = new MDL(proj, proj.getMarkerSet(), proj.getMarkerNames());
		while (mdl.hasNext()) {
			MarkerData md = mdl.next();
			for (float x : md.getXs()) {
				writer.println(x);
			}
		}
		mdl.shutdown();

		for (String file : new File(proj.MARKER_DATA_DIRECTORY.getValue()).list((e, f) -> {
			return f.endsWith(MarkerData.MARKER_DATA_FILE_EXTENSION);
		})) {
			System.out.println(file);
			for (Entry<String, Float> entry : TransposeData.loadOutliersFromRAF(proj.MARKER_DATA_DIRECTORY.getValue()
																																					+ file)
																										 .entrySet()) {
				if (entry.getKey().endsWith("\tx")) {
					writer.println(entry.getValue());
				}
			}
		}

		writer.close();
	}

	private static void runXYHistogram(Project proj) {
		double scale = proj.XY_SCALE_FACTOR.getValue().doubleValue();
		int binSize = proj.XY_SCALE_FACTOR.getValue().intValue();
		System.out.println("Bin Size: " + binSize);

		HashMap<Integer, AtomicInteger> binXCounts = new HashMap<>();

		boolean[] sampling = new boolean[proj.getMarkerNames().length];
		int every = 5;
		for (int i = 0; i < sampling.length; i++) {
			if (i % every == 0) {
				sampling[i] = true;
			}
		}

		MDL mdl = new MDL(proj, proj.getMarkerSet(), ArrayUtils.subArray(proj.getMarkerNames(),
																																		 sampling));
		while (mdl.hasNext()) {
			MarkerData md = mdl.next();
			for (float x : md.getXs()) {
				if (Float.isNaN(x)) {
					continue;
				}
				int bin = (int) (((double) (x * scale)) / binSize); // correct???
				AtomicInteger cnt = binXCounts.get(bin);
				if (cnt == null) {
					cnt = new AtomicInteger(0);
					binXCounts.put(bin, cnt);
				}
				cnt.incrementAndGet();
			}
		}
		mdl.shutdown();

		for (Entry<Integer, AtomicInteger> entry : binXCounts.entrySet()) {
			System.out.println(entry.getKey() + "\t" + entry.getValue().get());
		}
	}

	static float fromByteArrayBB(byte[] bytes) {
		return ByteBuffer.wrap(bytes).getFloat();
	}

	static void runHRC() {
		String fileDir = "F:/BCX2/sampleFiles/";
		String idFile = "EA.ids_sorted.txt";
		String idLookup = "ids_lookup.txt";
		String[] files = new File(fileDir).list((d, f) -> {
			return f.endsWith(".xln");
		});
		String[] samples = HashVec.loadFileToStringArray(fileDir + idFile, false, null, false);
		Map<String, String> idMap = HashVec.loadFileColumnToMap(fileDir + idLookup, 0, 1, true, null);
		for (String f : files) {
			String out = fileDir + ext.rootOf(f) + ".sample";
			PrintWriter writer = Files.getAppropriateWriter(out);
			Hashtable<String, String> in = HashVec.loadFileToHashString(fileDir + f, 0,
																																	new int[] {2, 3, 1},
																																	" ", true);
			writer.println("ID_1 ID_2 missing PC1 PC2 " + f.split("_")[0]);
			writer.println("0 0 0 C C P");
			for (String s : samples) {
				String v = in.get(idMap.get(s));
				if (v == null) {
					v = "NaN NaN NaN";
				}
				writer.println(s + " " + s + " 0 " + v);
			}
			writer.close();
		}
	}

	static void testWriters() {

		String dir = "C:/mass/";
		int numFiles = 0;
		PrintWriter writer;
		HashMap<String, PrintWriter> writerMap = new HashMap<>();
		boolean temp = true;
		while (temp) {
			String f = dir + ("" + numFiles) + ".out";
			writer = Files.getAppropriateWriter(f, true);
			writerMap.put(f, writer);
			numFiles++;
		}
		for (PrintWriter writer1 : writerMap.values()) {
			writer1.close();
		}
	}

	private static void processAnnotationFiles() throws IOException {
		String dir = "F:/Flow/Annotation/final_1/panel 1 branch/";
		String all = dir + "allAnnots_good.txt";
		HashSet<String> allFiles = new HashSet<>();
		BufferedReader r = Files.getAppropriateReader("F:/Flow/Annotation/final_1/fcsFiles.txt");
		String ln = null;
		while ((ln = r.readLine()) != null) {
			allFiles.add(ln);
		}
		r.close();
		String[] files = new File(dir).list();
		PrintWriter writerAll = Files.getAppropriateWriter(all);
		for (String f : files) {
			PrintWriter writer = Files.getAppropriateWriter(dir + f + ".sh");
			BufferedReader reader = Files.getAppropriateReader(dir + f);
			String line = null;
			HashMap<String, String> datesAndIdentsAndPanels = new HashMap<>();
			while ((line = reader.readLine()) != null) {
				if (line.startsWith("@ANNOT") || "".equals(line.trim())) {
					continue;
				}
				String[] pts = line.split("\\|");

				int p1 = 0;
				if (pts.length == 3 && pts[2].toLowerCase().contains("good")) {
					p1 = pts[1].toLowerCase().contains("p1")
							 || pts[1].toLowerCase().contains("panel_1")
							 || pts[1].toLowerCase().contains("panel 1") ? 1 : 2;

					pts = pts[1].split("/");
					String gate = pts[pts.length - 1].substring(pts[pts.length - 2].length() + 1);
					gate = gate.substring(0, gate.length() - 4);
					String samp = pts[pts.length - 1].substring(0, pts[pts.length - 1].indexOf(".fcs") + 4);
					// writerAll.println(samp);
					pts = samp.split("_");
					String date = pts[0];
					String ident = "";
					for (int i = 0; i < pts.length; i++) {
						if ((pts[i].startsWith("F") && !pts[i].startsWith("FORTESSA"))
								|| pts[i].startsWith("Ctl")
								|| pts[i].startsWith("PBMC")
								|| pts[i].startsWith("RR-")
								|| pts[i].startsWith("ZF-")
								|| pts[i].startsWith("BC-")
								|| pts[i].startsWith("HRS") || pts[i].startsWith("P2-")
								|| pts[i].startsWith("24HR") || pts[i].equals("G") || pts[i].equals("H")) {
							ident = pts[i];
							for (int i1 = i + 1; i1 < pts.length; i1++) {
								ident += "_"
												 + pts[i1].substring(0,
																						 pts[i1].length() - (pts[i1].endsWith(".fcs") ? 4 : 0));
								if (pts[i1].endsWith(".fcs"))
									break;
							}
							break;
						}
					}
					if (ident.equals("")) {
						// System.out.println(samp);
					} else {
						String key = date + "\t" + ident + "\t" + p1 + "\t" + gate;
						datesAndIdentsAndPanels.put(key, samp);
					}
				} else if (pts.length == 3) {
					// System.out.println(pts[2]);
				}
			}
			for (Entry<String, String> d : datesAndIdentsAndPanels.entrySet()) {
				String[] dI = d.getKey().split("\t");
				String fil = "";
				for (String s : allFiles) {
					if (s.contains(dI[0]) && s.contains(dI[1])) {
						if (dI[2].equals("1")) {
							if (s.toLowerCase().contains("p1")
									|| s.toLowerCase().contains("panel_1")
									|| s.toLowerCase().contains("panel 1")) {
								fil = s;
								break;
							}
						} else {
							if (s.toLowerCase().contains("p2")
									|| s.toLowerCase().contains("panel_2")
									|| s.toLowerCase().contains("panel 2")) {
								fil = s;
								break;
							}
						}
					}
				}
				writer.println("cp \"" + fil + "\" ./");
				if (fil.equals("")) {
					System.err.println(d.getKey());
				} else {
					writerAll.println(fil + "\t" + dI[3]);
				}
			}
			writer.close();
			reader.close();
		}
		writerAll.close();
	}

	public static void transposeFile() throws IOException {
		String file = "C:\\Users\\cole0482\\Desktop\\transpose.txt";
		String[][] data = HashVec.loadFileToStringMatrix(file, false, null, "\t", 0, true);
		String[][] data2 = new String[data[0].length][];
		for (int i = 0; i < data[0].length; i++) {
			data2[i] = new String[data.length];
		}

		for (int i = 0; i < data.length; i++) {
			for (int j = 0; j < data[i].length; j++) {
				data2[j][i] = data[i][j];
			}
		}

		for (int i = 0; i < data2.length; i++) {
			System.out.println(ArrayUtils.toStr(data2[i]));
		}
	}

	public static void writeFCSLookup() throws IOException {
		String file1 = "F:/Flow/Annotation/final_1/fcsFiles.txt";
		String file2 = "F:/Flow/Annotation/final_1/fcsLookup.txt";
		BufferedReader reader = Files.getAppropriateReader(file1);
		PrintWriter writer = Files.getAppropriateWriter(file2);
		String line = null;
		while ((line = reader.readLine()) != null) {
			writer.println(line + "\t" + ext.removeDirectoryInfo(line) + "\t"
										 + ext.replaceWithLinuxSafeCharacters(ext.removeDirectoryInfo(line)));
		}
		writer.close();
	}

	public static void main(String[] args) throws IOException, ClassNotFoundException {
		int numArgs = args.length;
		Project proj;
		String filename = "lab.dat";
		String logfile = null;
		String file = null;

		boolean test = true;
		if (test) {

			double[] test2 = (double[]) null;

			// runHRC();
			// QQPlot.main(new String[]
			// {"files=F:/CARDIA 2017/2nd round/results/plots/combined.results"});

			// String[] args1 = {
			// "file=F:/CARDIA 2017/2nd round/results/plots/combined.results"};
			// ManhattanPlot.main(args1);

			// CARDIA2017ResultsProcessor.combineChrXDose("G:/CARDIA_DATA/AA/");
			// CARDIA2017ResultsProcessor.combineChrXInfo("G:/CARDIA_DATA/AA/");

			// byte[] pt5 = ByteBuffer.allocate(4).putFloat(1.5f).array();
			// byte[] pt35 = ByteBuffer.allocate(4).putFloat(65.35f).array();
			// byte[] pt75 = ByteBuffer.allocate(4).putFloat(1578.75f).array();
			// byte[] Opt75 = ByteBuffer.allocate(4).putFloat(-42.75f).array();
			//
			// System.out.println(BGENBitMath.bytesToFloat(true, pt5) + " - " + fromByteArrayBB(pt5));
			// System.out.println(BGENBitMath.bytesToFloat(true, pt35) + " - " + fromByteArrayBB(pt35));
			// System.out.println(BGENBitMath.bytesToFloat(true, pt75) + " - " + fromByteArrayBB(pt75));
			// System.out.println(BGENBitMath.bytesToFloat(true, Opt75) + " - " + fromByteArrayBB(Opt75));

			// String dir = "F:/testProjectSrc/UKBB_AffyAxiom/";
			// UKBBParsingPipeline pipe = new UKBBParsingPipeline();
			// pipe.setSourceDir(dir + "00src/");
			// pipe.setProjectDir(dir + "project/");
			// pipe.setProjectPropertiesDir("D:/projects/");
			// pipe.setFamFile(dir + "ukb1773_l2r_chrY_v2_s488374.fam");
			// pipe.setProjectName("UKBB");
			// pipe.runPipeline();

			// testWriters();
			// processAnnotationFiles();
			// writeFCSLookup();

			// transposeFile();

			// proj = new Project(args[0]);
			// runXYHistogram(proj);

			// dumpSingleMDRAFOutliers();
			// String dir = "/home/pankrat2/shared/aric_gw6/ARICGenvisis_CEL_FULL/plinkApril2017/";
			// String mkrInfoFile = "/home/pankrat2/cole0482/Affy6_duplicates.txt";
			// String missDropsFile = dir + "quality_control/further_analysis_QC/miss_drops.dat";
			// String callrateFile =
			// "/home/pankrat2/shared/aric_gw6/ARICGenvisis_CEL_FULL/ohw_ws_20_ALL1000PCs_gc_corrected_OnTheFly_LRR_035CR_096_SAMP_LRR_TRUE_markerQC.txt";
			// String outFile = dir + "markerDuplicates.out";
			// markerDuplicateFilter(mkrInfoFile, missDropsFile, callrateFile, outFile);

			// addPFBToMarkerMetrics();

			// proj = new Project("/home/pankrat2/cole0482/projects/Poynter1_Shadow.properties", false);
			// System.out.println("Loading PFB for test: " + proj.CUSTOM_PFB_FILENAME.getValue());
			// PFB pfb = PFB.loadPFB(proj, proj.CUSTOM_PFB_FILENAME.getValue());
			// System.out.println(pfb.getPfbs().length + " pfb entries");
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
			// affy6BimLookup();
			// exomeRecode();
			// filt();
			// famRecode();
			// String bimFile;
			// String newBimFile;
			// String dir = "D:/temp/plink/exmToRS/";
			// bimFile = dir + "exome_EA.bim";
			// newBimFile = dir + "exome_EA_corrected.bim";
			// exomeRecode(bimFile, newBimFile);
			// bimFile = dir + "exome_AA.bim";
			// newBimFile = dir + "exome_AA_corrected.bim";
			// exomeRecode(bimFile, newBimFile);

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

			// CNVHelper.generateRegionsFileFromCNVFile("D:/data/ny_registry/new_york/cnvs/prune_belly.cnv");
			// BeastScore.scoreCNVFile(new Project("D:/projects/NY_Registry_Combo_Shadow.properties",
			// false),
			// "D:/data/ny_registry/new_york/cnvs/prune_belly.cnv",
			// true);
			// GeneScorePipeline.preprocessDataFiles(new String[]
			// {"D:/GeneScorePipe/Poynter/SnpInfo_Orig.xln"});

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
			proj = new Project(filename, logfile);
			if (file != null) {
				idSwap(proj, file);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
