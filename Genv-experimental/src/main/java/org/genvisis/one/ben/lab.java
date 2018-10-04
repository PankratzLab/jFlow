package org.genvisis.one.ben;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;
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
import org.genvisis.common.parsing.AliasedFileColumn;
import org.genvisis.common.parsing.FileColumn;
import org.genvisis.common.parsing.FileLink;
import org.genvisis.common.parsing.FileParserFactory;
import org.genvisis.common.parsing.StandardFileColumns;
import org.genvisis.filesys.CNVariant;
import org.genvisis.filesys.DosageData;
import org.genvisis.filesys.Segment;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.HashMultiset;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multiset;
import htsjdk.samtools.seekablestream.SeekableFileStream;
import htsjdk.samtools.util.BlockCompressedFilePointerUtil;
import htsjdk.samtools.util.BlockCompressedInputStream;

public class lab {

  private static void countCNVsForIndividuals(String indivFile, String cnvFile,
                                              String outFile) throws IOException {
    Hashtable<String, String> sampleKeyHash = new Hashtable<>();
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
      sb.append("\t").append(line[2]).append("\t").append(line[3]).append("\t").append(line[4])
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
      if (rs == null || "---".equals(rs)) missCnt++;
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
      if (line[1].equals("---")) continue;
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
    listAllSamplesInProj = new String[] {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                                         "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
                                         "21", "22", "23",};
    complete.add(listAllSamplesInProj[5]);
    complete.add(listAllSamplesInProj[14]);
    complete.add(listAllSamplesInProj[15]);
    // complete.add(listAllSamplesInProj[16]);
    batchMax = 4;

    int[][] ranges = ArrayUtils.splitUpIntoBinsOfIndices(listAllSamplesInProj, complete, batchMax,
                                                         log);
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
        writer.println(date + "\t" + fn + "\t" + (used.contains(fn) ? "1" : "0") + "\t" + s);
      }
    }

    writer.flush();
    writer.close();
  }

  private static void vcfTest(String filename, String probAttrName) {
    String dir = ext.parseDirectoryOfFile(filename);
    DosageData dd1 = new DosageData(filename, null, null, false, new Logger());
    DosageData dd2 = DosageData.loadVCF(filename, null, null, null, null, null, new Logger());
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
    SerializedFiles.writeSerial(allOutliers,
                                proj.MARKER_DATA_DIRECTORY.getValue(true, true) + "outliers.ser");
    System.out.println("Rebuilt outliers in " + ext.getTimeElapsedNanos(t2));
  }

  private static void dumpSingleMDRAFOutliers() {
    String[] files = {"/scratch.global/cole0482/UKBB2/project/transposed/markers.12.19398.20890.mdRAF",
                      "/scratch.global/cole0482/UKBB2/project/transposed/markers.7.17681.19154.mdRAF",
                      "/scratch.global/cole0482/UKBB2/project/transposed/markers.3.5980.7475.mdRAF",
                      "/scratch.global/cole0482/UKBB2/project/transposed/markers.4.7415.8898.mdRAF",};
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

    MDL mdl = new MDL(proj, proj.getMarkerNames());
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

    MDL mdl = new MDL(proj, ArrayUtils.subArray(proj.getMarkerNames(), sampling));
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
                                                                  new int[] {2, 3, 1}, " ", true);
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

  private static void processAnnotationFilesAll() throws IOException {
    String dir = "F:/Flow/Annotation/allMinusManuals/";

    String[] files = new File(dir).list();
    String[] annots = {"abnormal", "bad", "maybe", "not good", "pbs", "unsure"};
    HashSet<String> annotSet = new HashSet<>();
    for (String a : annots) {
      annotSet.add(a);
    }

    HashMap<String, HashMap<String, ArrayList<String>>> fileAnnotMap = new HashMap<>();

    for (String f : files) {
      HashMap<String, ArrayList<String>> annotMap = new HashMap<>();
      fileAnnotMap.put(f, annotMap);
      for (String a : annots) {
        annotMap.put(a, new ArrayList<String>());
      }
      BufferedReader reader = Files.getAppropriateReader(dir + f);
      String line = null;

      while ((line = reader.readLine()) != null) {
        if (line.startsWith("@ANNOT") || "".equals(line.trim())) {
          continue;
        }
        String[] pts = line.split("\\|");
        if (pts.length != 3) {
          continue;
        }
        if (annotSet.contains(pts[2].toLowerCase())) {
          annotMap.get(pts[2].toLowerCase()).add(pts[1]);
        }
      }

      reader.close();
    }

    PrintWriter writer = Files.getAppropriateWriter(dir + "out.xln");
    for (Entry<String, HashMap<String, ArrayList<String>>> fileEntry : fileAnnotMap.entrySet()) {
      String annotFile = fileEntry.getKey();

      for (Entry<String, ArrayList<String>> annotEntry : fileEntry.getValue().entrySet()) {
        String annot = annotEntry.getKey();

        for (String imgFile : annotEntry.getValue()) {
          String[] pathParts = imgFile.split("/");
          String fcs = pathParts[pathParts.length - 2];

          writer.println(annotFile + "\t" + annot + "\t" + fcs + "\t"
                         + pathParts[pathParts.length - 1].substring(fcs.length() + 1));
        }
      }

    }
    writer.close();

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
          p1 = pts[1].toLowerCase().contains("p1") || pts[1].toLowerCase().contains("panel_1")
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
                || pts[i].startsWith("Ctl") || pts[i].startsWith("PBMC") || pts[i].startsWith("RR-")
                || pts[i].startsWith("ZF-") || pts[i].startsWith("BC-") || pts[i].startsWith("HRS")
                || pts[i].startsWith("P2-") || pts[i].startsWith("24HR") || pts[i].equals("G")
                || pts[i].equals("H")) {
              ident = pts[i];
              for (int i1 = i + 1; i1 < pts.length; i1++) {
                ident += "_"
                         + pts[i1].substring(0,
                                             pts[i1].length() - (pts[i1].endsWith(".fcs") ? 4 : 0));
                if (pts[i1].endsWith(".fcs")) break;
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
              if (s.toLowerCase().contains("p1") || s.toLowerCase().contains("panel_1")
                  || s.toLowerCase().contains("panel 1")) {
                fil = s;
                break;
              }
            } else {
              if (s.toLowerCase().contains("p2") || s.toLowerCase().contains("panel_2")
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

  private static void processFreqs() throws IOException {
    String dir = "F:\\Flow\\CBC_processing\\final counts\\freqs\\";

    String p1Lookup = dir + "p1.report.txt";
    String p1Data = dir + "all.p1.cnts_fixed.xln";
    String p1Out = dir + "all.p1.freqs.xln";

    String p2Lookup = dir + "p2.report.txt";
    String p2Data = dir + "all.p2.cnts_fixed.xln";
    String p2Out = dir + "all.p2.freqs.xln";

    String lookup[][];
    BufferedReader reader;
    PrintWriter writer;
    String l;

    reader = Files.getAppropriateReader(p1Data);
    reader.readLine();
    writer = Files.getAppropriateWriter(p1Out);
    lookup = HashVec.loadFileToStringMatrix(p1Lookup, true, new int[] {0, 1, 2});
    l = "Sample";
    for (String[] s : lookup) {
      l += "\t" + s[0];
    }
    writer.println(l);
    l = null;
    while ((l = reader.readLine()) != null) {
      String[] pts = l.split("\t");
      StringBuilder out = new StringBuilder(pts[0]);
      for (String[] ln : lookup) {
        int ind1 = Integer.parseInt(ln[1]);
        int ind2 = Integer.parseInt(ln[2]);
        if (pts[ind1].equals("null")) {
          out.append("\t").append(".");
        } else {
          double cnt1 = Integer.parseInt(pts[ind1]);
          double cnt2 = Integer.parseInt(pts[ind2]);
          if (cnt2 == 0) {
            out.append("\t").append(".");
          } else {
            out.append("\t").append(cnt1 / cnt2);
          }
        }
      }
      writer.println(out.toString());
    }
    writer.close();

    reader = Files.getAppropriateReader(p2Data);
    reader.readLine();
    writer = Files.getAppropriateWriter(p2Out);
    lookup = HashVec.loadFileToStringMatrix(p2Lookup, true, new int[] {0, 1, 2});
    l = "Sample";
    for (String[] s : lookup) {
      l += "\t" + s[0];
    }
    writer.println(l);
    l = null;
    while ((l = reader.readLine()) != null) {
      String[] pts = l.split("\t");
      StringBuilder out = new StringBuilder(pts[0]);
      for (String[] ln : lookup) {
        int ind1 = Integer.parseInt(ln[1]);
        int ind2 = Integer.parseInt(ln[2]);
        if (pts[ind1].equals("null")) {
          out.append("\t").append(".");
        } else {
          double cnt1 = Integer.parseInt(pts[ind1]);
          double cnt2 = Integer.parseInt(pts[ind2]);
          if (cnt2 == 0) {
            out.append("\t").append(".");
          } else {
            out.append("\t").append(cnt1 / cnt2);
          }
        }
      }
      writer.println(out.toString());
    }
    writer.close();

  }

  private static void removeParens() {
    String dir = "F:\\Flow\\CBC_processing\\final counts\\redoEverythingSource\\fixedHeaders\\";
    for (String d : (new File(dir).list())) {
      for (String f : new File(dir + d).list()) {
        String path = dir + d + "\\" + f;
        String[] header = Files.getHeaderOfFile(path, null);
        for (int i = 0; i < header.length; i++) {
          int p = header[i].lastIndexOf(" / ");
          if (p >= 0) {
            header[i] = header[i].substring(p + 3, header[i].length());
          }
          header[i] = header[i].substring(0, header[i].indexOf('(')).trim();
        }
        String finalHdr = "\t" + ArrayUtils.toStr(header, "\t");
        String[] data = HashVec.loadFileToStringArray(path, false, null, false);
        data[0] = finalHdr;
        Files.writeArray(data, path);
      }
    }
  }

  private static void fixBCXResults() {
    String dir = "F:\\BCX2\\results\\";
    String out = "F:\\BCX2\\results\\fixed\\";
    String map = "F:\\BCX2\\results\\rsq\\";

    if (!Files.isWindows()) {
      dir = "/scratch.global/cole0482/CARDIA/HRC/results/";
      out = "/scratch.global/cole0482/CARDIA/HRC/results_fixed/";
      map = "/scratch.global/cardia/HRC/";
    }

    List<FileLink> rsqFiles = new ArrayList<>();
    AliasedFileColumn snpCol = StandardFileColumns.snp("SNP");
    FileColumn<String> rsqCol = StandardFileColumns.allExcept("\t", snpCol);
    for (String f : new File(map).list(new FilenameFilter() {

      @Override
      public boolean accept(File dir, String name) {
        return name.endsWith("rsq.txt");
      }
    })) {
      rsqFiles.add(FileLink.setup(map + f).keys(snpCol).values(rsqCol));
    }

    String[] files = {"CARDIA_BAS_EA_141217_BC.txt.gz", "CARDIA_LYM_EA_141217_BC.txt.gz",
                      "CARDIA_MON_EA_141217_BC.txt.gz", "CARDIA_WBC_EA_141217_BC.txt.gz",
                      "CARDIA_EOS_EA_141217_BC.txt.gz", "CARDIA_MCHC_EA_141217_BC.txt.gz",
                      "CARDIA_NEU_EA_141217_BC.txt.gz", "CARDIA_HCT_EA_141217_BC.txt.gz",
                      "CARDIA_MCH_EA_141217_BC.txt.gz", "CARDIA_PLT_EA_141217_BC.txt.gz",
                      "CARDIA_HGB_EA_141217_BC.txt.gz", "CARDIA_MCV_EA_141217_BC.txt.gz",
                      "CARDIA_RBC_EA_141217_BC.txt.gz"};

    for (String f : files) {
      AliasedFileColumn snpCol1 = StandardFileColumns.snp("SNP");
      FileParserFactory factory = FileParserFactory.setup(dir + f, snpCol1,
                                                          StandardFileColumns.allExcept("\t",
                                                                                        snpCol1));
      for (FileLink fl : rsqFiles) {
        factory.link(fl);
      }
      try {
        factory.build().parseToFile(out + f, "\t");
      } catch (IOException e) {
        e.printStackTrace();
      }
    }

  }

  // private static void fixBCXResultsAlleleFreq() {
  // String dir = "F:\\BCX2\\results\\";
  // String out = "F:\\BCX2\\results\\fixed\\";
  // String map = "F:\\BCX2\\results\\rsq\\";
  //
  // if (!Files.isWindows()) {
  // dir = "/scratch.global/cole0482/CARDIA/HRC/results_fixed/";
  // out = "/scratch.global/cole0482/CARDIA/HRC/results_fixed/freqs/";
  // map = "/scratch.global/cardia/HRC/";
  // }
  //
  // List<FileLink> rsqFiles = new ArrayList<>();
  // AliasedFileColumn snpCol = StandardFileColumns.snp("SNP");
  // FileColumn<String> freqCol = new AliasedFileColumn("EAF", "ALT_frq");
  // for (String f : new File(map).list(new FilenameFilter() {
  // @Override
  // public boolean accept(File dir, String name) {
  // return name.endsWith(".alleleFreqs");
  // }
  // })) {
  // rsqFiles.add(FileLink.setup(map + f).keys(snpCol)
  // .values(freqCol));
  // }
  // for (FileLink fl : rsqFiles) {
  // fl.build();
  // }
  //
  // String[] files = {
  // "CARDIA_BAS_EA_141217_BC.txt.gz",
  // "CARDIA_LYM_EA_141217_BC.txt.gz",
  // "CARDIA_MON_EA_141217_BC.txt.gz",
  // "CARDIA_WBC_EA_141217_BC.txt.gz",
  // "CARDIA_EOS_EA_141217_BC.txt.gz",
  // "CARDIA_MCHC_EA_141217_BC.txt.gz",
  // "CARDIA_NEU_EA_141217_BC.txt.gz",
  // "CARDIA_HCT_EA_141217_BC.txt.gz",
  // "CARDIA_MCH_EA_141217_BC.txt.gz",
  // "CARDIA_PLT_EA_141217_BC.txt.gz",
  // "CARDIA_HGB_EA_141217_BC.txt.gz",
  // "CARDIA_MCV_EA_141217_BC.txt.gz",
  // "CARDIA_RBC_EA_141217_BC.txt.gz"
  // };
  //
  // for (String f : files) {
  // AliasedFileColumn snpCol1 = StandardFileColumns.snp("SNP");
  // FileColumn<String> chr = new AliasedFileColumn("CHR", "CHR");
  // FileColumn<String> POS = new AliasedFileColumn("POS", "POS");
  // FileColumn<String> STRAND = new AliasedFileColumn("STRAND", "STRAND");
  // FileColumn<String> EFFECT_ALLELE = new AliasedFileColumn("EFFECT_ALLELE", "EFFECT_ALLELE");
  // FileColumn<String> OTHER_ALLELE = new AliasedFileColumn("OTHER_ALLELE", "OTHER_ALLELE");
  // FileColumn<String> N = new AliasedFileColumn("N", "N");
  // FileColumn<String> BETA = new AliasedFileColumn("BETA", "BETA");
  // FileColumn<String> SE = new AliasedFileColumn("SE", "SE");
  // FileColumn<String> PVAL = new AliasedFileColumn("PVAL", "PVAL");
  // FileColumn<String> Rsq = new AliasedFileColumn("Rsq", "Rsq");
  // List<FileColumn<?>> order = Lists.newArrayList(snpCol1, chr, POS, STRAND, EFFECT_ALLELE,
  // OTHER_ALLELE, N, freqCol, BETA, SE, PVAL, Rsq);
  // FileParserFactory factory = FileParserFactory.setup(dir + f, snpCol1, chr, POS, STRAND,
  // EFFECT_ALLELE, OTHER_ALLELE, N, BETA, SE,
  // PVAL, Rsq);
  // for (FileLink fl : rsqFiles) {
  // factory.link(fl);
  // }
  // try {
  // factory.build().parseToFile(out + f, "\t", order);
  // } catch (IOException e) {
  // e.printStackTrace();
  // }
  // }
  //
  // }
  //
  // public static void createBCXPlots() {
  // String[] files = {
  // // "CARDIA_BAS_EA_141217_BC.txt.gz",
  // // "CARDIA_LYM_EA_141217_BC.txt.gz",
  // // "CARDIA_MON_EA_141217_BC.txt.gz",
  // "CARDIA_WBC_EA_141217_BC.txt.gz",
  // "CARDIA_EOS_EA_141217_BC.txt.gz",
  // "CARDIA_MCHC_EA_141217_BC.txt.gz",
  // "CARDIA_NEU_EA_141217_BC.txt.gz",
  // "CARDIA_HCT_EA_141217_BC.txt.gz",
  // "CARDIA_MCH_EA_141217_BC.txt.gz",
  // "CARDIA_PLT_EA_141217_BC.txt.gz",
  // "CARDIA_HGB_EA_141217_BC.txt.gz",
  // "CARDIA_MCV_EA_141217_BC.txt.gz",
  // "CARDIA_RBC_EA_141217_BC.txt.gz"
  // };
  // String dir1 = "F:\\BCX2\\results\\";
  // String dir2 = "F:\\BCX2\\results\\fixed\\";
  // String out = "F:\\BCX2\\results\\fixed\\plots\\";
  //
  // for (String f : files) {
  // String f1 = dir1 + f;
  // String o1 = out + f + ".png";
  // String f2 = dir2 + f;
  // String o2 = out + f + ".fixed.png";
  // new Thread(() -> {
  // AFPlot p;
  // p = new AFPlot(null);
  // p.loadFromFile(f1, null);
  // p.waitForData();
  // p.screenshot(o1);
  // p = null;
  // }).start();
  //
  // new Thread(() -> {
  // AFPlot p;
  // p = new AFPlot(null);
  // p.loadFromFile(f2, null);
  // p.waitForData();
  // p.screenshot(o2);
  // p = null;
  // }).start();
  // }
  //
  // }

  public static void checkMissing() {
    Project proj = new Project("/home/pankrat2/cole0482/projects/UKBioBank.properties");
    float[] lrrs;
    int count = 0;
    for (String s : proj.getSamples()) {
      String sF = proj.SAMPLE_DIRECTORY.getValue() + s + ".sampRAF";
      Sample samp = Sample.loadFromRandomAccessFile(sF);
      lrrs = samp.getLRRs();
      if (lrrs == null || lrrs.length == 0 || ArrayUtils.removeNaN(lrrs).length == 0) {
        count++;
      }
      samp = null;
      lrrs = null;
    }
    System.out.println("Found " + count + " samples with missing or NaN LRR data");
  }

  private static final void combineAllWithSelect() throws IOException {
    ImmutableMap.Builder<String, String> b = new ImmutableMap.Builder<>();
    b.put("Comp-BV 605-A \\(CD95\\)", "Comp-BV605-A \\(CD95\\)");
    b.put("Comp-BV 510-A \\(CD28\\)", "Comp-BV510-A \\(CD28\\)");
    b.put("Comp-BB 515-A \\(CD27\\)", "Comp-BB515-A \\(CD27\\)");
    b.put("Comp-BB515-A \\(CD27\\)", "Comp-FITC-A \\(CD27\\)");
    b.put("Comp-BV 421-A \\(CCR7\\)", "Comp-BV421-A \\(CCR7\\)");
    b.put("Comp-BV 711-A \\(CD45RA\\)", "Comp-BV711-A \\(CD45RA\\)");
    b.put("Comp-BUV 395-A \\(CD8\\)", "Comp-BUV396-A \\(CD8\\)");
    b.put("LIVE/DEAD", "L/D");
    b.put("Comp-PE-Cy7 \\(blue\\)-A \\(CD19\\)", "Comp-PE-Cy7-A \\(CD19\\)");
    b.put("Comp-BUV 737-A \\(IgD\\)", "Comp-BUV737-A \\(IgD\\)");

    final ImmutableMap<String, String> dimSwitch = b.build();

    String dir = "/scratch.global/cole0482/fcsVizPipe/r26_TcellSubs_Kmeans_wsp_v8_cleanup/allOuts/";
    //    String dir1 = dir + "outCnts/";
    //    String sel = dir + "nonMan.samp.txt";
    int batch = 6;
    //
    //    HashSet<String> selHash = Sets.newHashSet(HashVec.loadFileToStringArray(sel, false, null, true,
    //                                                                            false, "\t"));
    HashSet<String> fndHash = new HashSet<>();
    //    System.out.println("Loaded " + selHash.size() + " samples to keep");
    //
    //    Map<String, String> xtraFiles = new HashMap<>();
    //    xtraFiles.put("/scratch.global/lanej/flow/manual/panel1_v8_counts/p1.cnts.xln", "V8_manual");
    //    xtraFiles.put("/scratch.global/lanej/flow/manual/panel1_v3_counts/p1.cnts.xln", "V3_manual");
    //    xtraFiles.put("/scratch.global/lanej/flow/manual/kmeans_consolidated_counts/p1.cnts.xln",
    //                  "consol_manual");
    //    xtraFiles.put("/scratch.global/lanej/flow/manual/kmeans_Panel1_bcellsubs_regated_counts/p1.cnts.xln",
    //                  "bcell_manual");

    LinkedHashSet<String> headers = new LinkedHashSet<>();
    for (int i = 0; i < batch + 1; i++) {
      String f = dir + "b" + i + ".cnts.xln";
      if (!Files.exists(f)) continue;
      String[] hdr = Files.getHeaderOfFile(f, null);
      for (String s : hdr) {
        String[] s1 = s.split(" / ");
        String s11 = s1[s1.length - 1];
        for (Entry<String, String> en : dimSwitch.entrySet()) {
          s11 = s11.replaceAll(en.getKey(), en.getValue());
        }
        headers.add(s11);
      }
    }
    //    for (String f : xtraFiles.keySet()) {
    //      String[] hdr = Files.getHeaderOfFile(f, null);
    //      for (String s : hdr) {
    //        String[] s1 = s.split(" / ");
    //        String s11 = s1[s1.length - 1];
    //        for (Entry<String, String> en : dimSwitch.entrySet()) {
    //          s11 = s11.replaceAll(en.getKey(), en.getValue());
    //        }
    //        headers.add(s11);
    //      }
    //    }

    PrintWriter allOut = Files.getAppropriateWriter(dir + "allCounts.xln");
    allOut.print("Sample\tSource");
    for (String h : headers) {
      allOut.print("\t");
      allOut.print(h);
    }
    allOut.println();

    HashSet<String> sampsFound = new HashSet<>();
    //    for (String f : xtraFiles.keySet()) {
    //      if (!Files.exists(f)) continue;
    //      BufferedReader reader = Files.getAppropriateReader(f);
    //      String line = reader.readLine();
    //      String[] hdr = line.split("\t", -1);
    //      String[] fix = new String[hdr.length];
    //      for (int s = 0; s < hdr.length; s++) {
    //        String[] s1 = hdr[s].split(" / ");
    //        String s11 = s1[s1.length - 1];
    //        for (Entry<String, String> en : dimSwitch.entrySet()) {
    //          s11 = s11.replaceAll(en.getKey(), en.getValue());
    //        }
    //        fix[s] = s11;
    //      }
    //
    //      Map<String, List<Integer>> hdrInds = new HashMap<>();
    //      for (int h = 0; h < fix.length; h++) {
    //        if (!hdrInds.containsKey(fix[h])) {
    //          hdrInds.put(fix[h], new ArrayList<>());
    //        }
    //        hdrInds.get(fix[h]).add(h);
    //      }
    //
    //      while ((line = reader.readLine()) != null) {
    //        int ind = line.indexOf("\t");
    //        String samp = line.substring(0, ind);
    //        samp = samp.substring(samp.lastIndexOf('/') + 1);
    //        if (!sampsFound.add(samp)) continue;
    //
    //        String[] sP = line.split("\t");
    //        allOut.print(samp);
    //        allOut.print("\t");
    //        allOut.print(xtraFiles.get(f));
    //        for (String h : headers) {
    //          allOut.print("\t");
    //          if (hdrInds.containsKey(h)) {
    //            boolean p = false;
    //            for (int hI : hdrInds.get(h)) {
    //              if (!sP[hI].equals("null") && !sP[hI].equals("")) {
    //                allOut.print(sP[hI]);
    //                p = true;
    //                break;
    //              }
    //            }
    //            if (!p) {
    //              allOut.print("null");
    //            }
    //          } else {
    //            allOut.print("null");
    //          }
    //        }
    //        allOut.println();
    //      }
    //      reader.close();
    //    }

    for (int i = 0; i < batch + 1; i++) {
      String f = dir + "b" + i + ".cnts.xln";
      //      String f = dir1 + "b" + i + "/p1.cnts.xln";
      if (!Files.exists(f)) continue;
      BufferedReader reader = Files.getAppropriateReader(f);
      String line = reader.readLine();
      String[] hdr = line.split("\t", -1);
      String[] fix = new String[hdr.length];
      for (int s = 0; s < hdr.length; s++) {
        String[] s1 = hdr[s].split(" / ");
        String s11 = s1[s1.length - 1];
        for (Entry<String, String> en : dimSwitch.entrySet()) {
          s11 = s11.replaceAll(en.getKey(), en.getValue());
        }
        fix[s] = s11;
      }

      Map<String, List<Integer>> hdrInds = new HashMap<>();
      for (int h = 0; h < fix.length; h++) {
        if (!hdrInds.containsKey(fix[h])) {
          hdrInds.put(fix[h], new ArrayList<>());
        }
        hdrInds.get(fix[h]).add(h);
      }

      int cnt = 0;
      while ((line = reader.readLine()) != null) {
        int ind = line.indexOf("\t");
        String samp = line.substring(0, ind);
        samp = samp.substring(samp.lastIndexOf('/') + 1);

        //        if (selHash.contains(samp)) {
        if (!fndHash.add(samp)) {
          System.out.println("Duplicate: " + samp);
        }
        if (!sampsFound.add(samp)) continue;
        cnt++;
        String[] sP = line.split("\t");
        allOut.print(samp);
        allOut.print("\tOpenCyto");
        for (String h : headers) {
          allOut.print("\t");
          if (hdrInds.containsKey(h)) {
            boolean p = false;
            for (int hI : hdrInds.get(h)) {
              if (!sP[hI].equals("null") && !sP[hI].equals("")) {
                allOut.print(sP[hI]);
                p = true;
                break;
              }
            }
            if (!p) {
              allOut.print("null");
            }
          } else {
            allOut.print("null");
          }
        }
        allOut.println();
        //        }
      }
      System.out.println("Found " + cnt + " samples in batch " + i);
      reader.close();
    }

    //    System.out.println("Missing " + (selHash.size() - fndHash.size()) + " samples");
    //    selHash.removeAll(fndHash);
    //    Files.writeIterable(selHash, dir + "missing.samp.txt");

    allOut.close();
  }

  static class SNPLine {

    Segment bin;
    String name;
    int dist;
    float mafDist;

  }

  private static void bamSnpComparison(String file, String out, int[] scales) throws IOException {
    BufferedReader reader = Files.getAppropriateReader(file);
    String line = null;
    String[] parts;
    String header = reader.readLine();
    Multiset<Segment> snpCounter = HashMultiset.create();
    Multimap<Segment, SNPLine> binMap = HashMultimap.create();
    Set<Segment> binList = new TreeSet<>();
    while ((line = reader.readLine()) != null) {
      parts = line.split("\t");
      float mafDist = Float.parseFloat(parts[3]);
      Segment seg = new Segment(parts[0]);
      snpCounter.add(seg);
      if (mafDist < .4) {
        SNPLine snp = new SNPLine();
        snp.bin = seg;
        snp.name = parts[1];
        snp.dist = Integer.parseInt(parts[2]);
        snp.mafDist = mafDist;
        binMap.put(snp.bin, snp);
        binList.add(snp.bin);
      }
    }

    System.out.println("Read all variants...");

    PrintWriter writer = Files.getAppropriateWriter(out);

    writer.print(header);
    for (int scale : scales) {
      writer.print("\t" + "MAF_" + scale);
    }
    for (int scale : scales) {
      writer.print("\t" + "SQRT_" + scale);
    }
    writer.println("\tNUM");

    for (Segment bin : binList) {
      Collection<SNPLine> snps = binMap.get(bin);
      for (SNPLine snp : snps) {
        int d = snp.dist;
        writer.print(snp.bin.getUCSClocation() + "\t" + snp.name + "\t" + snp.dist + "\t"
                     + snp.mafDist);
        for (int scale : scales) {
          float m = snp.mafDist * scale;
          writer.print("\t");
          writer.print(m);
        }
        for (int scale : scales) {
          float m = snp.mafDist * scale;
          double dd = Math.sqrt(d * d + m * m);
          writer.print("\t");
          writer.print(dd);
        }
        writer.println("\t" + snpCounter.count(snp.bin));
      }
    }
    writer.close();

    System.out.println("Finished selecting variants!");
  }

  private static void bamSnpSelection(String file, String out, int scale) throws IOException {

    BufferedReader reader = Files.getAppropriateReader(file);
    String line = null;
    String[] parts;
    String header = reader.readLine();
    Multimap<Segment, SNPLine> binMap = HashMultimap.create();
    Set<Segment> binList = new TreeSet<>();
    while ((line = reader.readLine()) != null) {
      parts = line.split("\t");
      SNPLine snp = new SNPLine();
      snp.bin = new Segment(parts[0]);
      snp.name = parts[1];
      snp.dist = Integer.parseInt(parts[2]);
      snp.mafDist = Float.parseFloat(parts[3]);
      if (snp.mafDist < .4) {
        binMap.put(snp.bin, snp);
        binList.add(snp.bin);
      }
    }

    System.out.println("Read all variants...");

    PrintWriter writer = Files.getAppropriateWriter(out);
    writer.println(header + "\t" + "SQRT_DIST");

    for (Segment bin : binList) {
      double closestDist = 10000;
      SNPLine closest = null;
      Collection<SNPLine> snps = binMap.get(bin);
      for (SNPLine snp : snps) {
        int d = snp.dist;
        float m = snp.mafDist * scale;
        double dd = Math.sqrt(d * d + m * m);
        if (dd < closestDist) {
          closestDist = dd;
          closest = snp;
        }
      }
      if (closest != null) {
        writer.println(closest.bin.getUCSClocation() + "\t" + closest.name + "\t" + closest.dist
                       + "\t" + closest.mafDist + "\t" + closestDist);
      } else {
        writer.println(bin.getChromosomeUCSC() + "\t.\t.\t.\t.");
      }
    }
    writer.close();

    System.out.println("Finished selecting variants!");

  }

  private static void buildUKBBLookup(String file) throws IOException {

    SeekableFileStream sfs = new SeekableFileStream(new File(file));
    int line = 0;
    Map<Integer, Long> indexMap = new HashMap<>();

    BlockCompressedInputStream bcis = new BlockCompressedInputStream(sfs);
    while (bcis.readLine() != null) {
      line++;
      if (line % 10000 == 0) {
        System.out.println(line);
      }
      bcis.seek(bcis.getFilePointer());
      indexMap.put(line, sfs.position());
    }
    bcis.close();

    System.out.println("Done with step 1!");

    //    Map<Integer, Long> byteIndexMap = new HashMap<>();
    //    for (int i = 1; i <= indexMap.size(); i++) {
    //      fis = new FileInputStream(file);
    //      bcis = new BlockCompressedInputStream(fis);
    //      bcis.seek(indexMap.get(i));
    //      byteIndexMap.put(i, fis.getChannel().position());
    //      bcis.close();
    //    }

    SerializedFiles.writeSerial(indexMap, ext.rootOf(file, false) + "_lookup.dat");
    System.out.println("Done with step 2!");
  }

  private static void checkUKBBLookup(String file) throws IOException {
    @SuppressWarnings("unchecked")
    Map<Integer, Long> indexMap = (Map<Integer, Long>) SerializedFiles.readSerial(ext.rootOf(file,
                                                                                             false)
                                                                                  + "_lookup.dat");

    SeekableFileStream sfs = new SeekableFileStream(new File(file));
    BlockCompressedInputStream bcis = new BlockCompressedInputStream(sfs);

    System.out.println(indexMap.size() + " lines in index file.");

    for (int i = 0; i < 10; i++) {
      int line = new Random().nextInt(indexMap.size());
      long ind = indexMap.get(line);
      sfs.seek(ind);
      System.out.println("Found " + bcis.readLine().split(" ").length + " in line " + line
                         + ", vfp: " + BlockCompressedFilePointerUtil.asString(ind));
    }
    bcis.close();

  }

  public static void main(String[] args) throws IOException, ClassNotFoundException {
    int numArgs = args.length;
    Project proj;
    String filename = "lab.dat";
    String logfile = null;
    String file = null;

    boolean test = true;
    if (test) {
      // createBCXPlots();
      // processAnnotationFilesAll();

      // processP2Counts_Step1();
      // processP2Counts_Step2();
      // processFreqs();

      // removeParens();
      // processP1Man_1(); // combine data within files
      // processP1Man_2(); // combine files together
      // combine(); // add meta data to file
      // fixBCXResults();
      // fixBCXResultsAlleleFreq();

      //      combineAllWithSelect();
      //      proj = new Project("D:\\projects\\AffyParsingTest.properties");
      //      OutOfRangeValues vals = OutOfRangeValues.construct(proj);
      //      for (String s : proj.getSamples()) {
      //        System.out.println(vals.getSampleOutliersForFile(proj, s).size() + " for " + s);
      //      }

      //      Dump.dumpMdRaf("F:\\testProjectSrc\\Affy1000G_small\\project2\\transposed\\markers.0.mdRAF",
      //                     new int[] {0, 1, 2, 3, 4, 5}, new Logger());
      //      Dump.dumpSampRaf("F:\\testProjectSrc\\Affy1000G_small\\project2\\samples\\"
      //                       + new Project("D:\\projects\\AffyParsingTest.properties").getSamples()[0]
      //                       + ".sampRAF");

      //      bamSnpComparison("G:\\bamTesting\\snpSelection\\All_Variants.xln",
      //                       "G:\\bamTesting\\snpSelection\\selected_scales.xln",
      //                       new int[] {1500, 2000, 5000});
      //      bamSnpSelection("G:\\bamTesting\\snpSelection\\All_Variants.xln",
      //                      "G:\\bamTesting\\snpSelection\\selected_scale2k.xln", 2000);
      //      bamSnpSelection("G:\\bamTesting\\snpSelection\\All_Variants.xln",
      //                      "G:\\bamTesting\\snpSelection\\selected_scale5k.xln", 5000);

      //      proj = new Project("D:/projects/Affy1000G Small3.properties");
      //      SampleList sl = SampleList.load(proj.SAMPLELIST_FILENAME.getValue());
      //      String[] samples = sl.getSamples();
      //      System.out.println();

      //      buildUKBBLookup("F:\\testProjectSrc\\UKBB_AffyAxiom\\00src\\ukb_baf_chr10_v2.txt.gz");
      //      checkUKBBLookup("F:\\testProjectSrc\\UKBB_AffyAxiom\\00src\\ukb_baf_chr10_v2.txt.gz");

      //      BGZipReader reader = new BGZipReader("F:\\testProjectSrc\\UKBB_AffyAxiom\\00src\\ukb_baf_chr10_v2.txt.gz");
      //      BlockCompressedInputStream bgzip = new BlockCompressedInputStream(new File("F:\\testProjectSrc\\UKBB_AffyAxiom\\00src_21gz\\ukb_baf_chr21_v2.txt.gz"));
      //      BufferedReader reader = new BufferedReader(new FileReader("F:\\testProjectSrc\\UKBB_AffyAxiom\\00src_21\\ukb_baf_chr21_v2.txt"));
      //
      //      String txtLn = reader.readLine();
      //      String zipLn = bgzip.readLine();
      //

      BufferedReader reader = Files.getAppropriateReader("/scratch.global/cole0482/UKBB/Axiom_UKB_WCSG.na35.annot.csv.zip");
      int lines = 0;
      while (reader.readLine() != null) {
        lines++;
      }
      System.out.println("Read " + lines
                         + " in /scratch.global/cole0482/UKBB/Axiom_UKB_WCSG.na35.annot.csv.zip");

      // runHRC();
      // QQPlot.main(new String[]
      // {"files=F:/CARDIA 2017/2nd round/results/plots/combined.results"});

      // String[] args1 = {
      // "file=F:/CARDIA 2017/2nd round/results/plots/combined.results"};
      // ManhattanPlot.main(args1);

      // CARDIA2017ResultsProcessor.combineChrXDose("G:/CARDIA_DATA/AA/");
      // CARDIA2017ResultsProcessor.combineChrXInfo("G:/CARDIA_DATA/AA/"); 

      //            String dir = "F:/testProjectSrc/UKBB_AffyAxiom/";
      //      UKBBParsingPipeline pipe = new UKBBParsingPipeline();
      //      pipe.setSourceDir(dir + "00src/");
      //      pipe.setProjectDir(dir + "project_21/");
      //      pipe.setProjectPropertiesDir("D:/projects/");
      //      pipe.setFamFile(dir + "ukb1773_l2r_chrY_v2_s488374.fam");
      //      pipe.setAnnotationCSV(dir + "Axiom_UKB_WCSG.na34.annot.csv");
      //      pipe.setProjectName("UKBB_21");
      //      pipe.runPipeline();

      // testWriters();
      // processAnnotationFiles();1
      // writeFCSLookup();
      //      proj = new Project("D:/projects/gedi_gwas.properties");
      //      MarkerDataLoader.buildOutliersFromMDRAFs(proj);

      //      proj.SAMPLE_DIRECTORY.setValue("G:\\transposeTesting\\sampleFiles\\");

      //      checkMissing();
      //      checkOOR();

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
