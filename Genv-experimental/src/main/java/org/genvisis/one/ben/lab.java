package org.genvisis.one.ben;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.cnv.analysis.FilterCalls;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.manage.PlinkData;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.CNVFilter;
import org.genvisis.common.CNVFilter.CNVFilterPass;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.common.QueueControl;
import org.genvisis.common.QueueControl.JobQueue;
import org.genvisis.common.ext;
import org.genvisis.filesys.CNVariant;
import org.genvisis.filesys.DosageData;
import org.genvisis.filesys.Segment;
import org.genvisis.gwas.MergeExtractPipeline;
import org.genvisis.one.ben.fcs.FCSFileDuplicator;
import org.genvisis.one.ben.fcs.FCSDataLoader;
import org.genvisis.one.ben.fcs.FCSDataLoader.LOAD_STATE;

import com.google.common.io.Closeables;

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

		Vector<CNVariant> cnvs = CNVariant.loadPlinkFile(cnvFile, sampleKeyHash, true, false);

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
				line = reader.readLine().trim().split("[\\s]+");
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



	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		Project proj;
		String filename = "lab.dat";
		String logfile = null;
		String file = null;

		boolean test = true;
		if (test) {


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


