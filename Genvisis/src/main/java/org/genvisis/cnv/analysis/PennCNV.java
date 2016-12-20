package org.genvisis.cnv.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.Vector;

import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.cnv.filesys.Pedigree;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.manage.Resources;
import org.genvisis.cnv.manage.Resources.Resource;
import org.genvisis.cnv.qc.SexChecks;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.PSF;
import org.genvisis.common.Positions;
import org.genvisis.common.SciStringComparator;
import org.genvisis.common.ext;
import org.genvisis.filesys.CNVariant;

public class PennCNV {
	public static final String[] QC_HEADS = {"LRR_mean", "LRR_median", "LRR_SD", "BAF_mean",
	                                         "BAF_median", "BAF_SD", "BAF_DRIFT", "WF", "GCWF"};
	public static final String[] ERRORS = {"large SD for LRR", "drifting BAF values",
	                                       "waviness factor values", "Small-sized CNV calls",
	                                       "NoCall rate"};
	public static final String QC_SUMMARY_FILE = "Sample_QC.xln";
	public static final int MISSING_SCORE = -1;

	public static void batch(Project proj, int numChunks, Vector<String> execList, String pfbFile,
	                         String gcmodelFile, String hmmFile, String scriptSubDir,
	                         String dataSubDir, String resultsSubDir) {
		String commands;
		PrintWriter writer;
		String[] files;
		int step;
		String execDir, dataDir, resultsDir, projDir, pennDir;
		Logger log;

		log = proj.getLog();
		final String runLine = Files.getRunString() + " " + PennCNV.class.getCanonicalName() + " proj="
		                       + new File(proj.getPropertyFilename()).getAbsolutePath();
		projDir = proj.PROJECT_DIRECTORY.getValue();
		execDir = proj.PENNCNV_EXECUTABLE_DIRECTORY.getValue(false, true);
		pennDir = proj.PENNCNV_RESULTS_DIRECTORY.getValue(false, true);
		dataDir = proj.PENNCNV_DATA_DIRECTORY.getValue(false, true) + dataSubDir;
		resultsDir = pennDir + resultsSubDir;

		String penncnvExecutable = execDir + "detect_cnv.pl";
		if (!Files.exists(penncnvExecutable)) {
			log.reportError("WARNING - couldn't find PennCNV executable 'detect_cnv.pl' in given directory: "
			                + execDir);
			if (Files.programExists("detect_cnv.pl")) {
				log.report("PennCNV executable 'detect_cnv.pl' found on the PATH; please check the PENNCNV_EXECUTABLE_DIRECTORY project property.");
			}
		}

		if (pfbFile != null) {
			pfbFile = ext.replaceTilde(pfbFile);
			if (!pfbFile.startsWith("/") && (pfbFile.charAt(1) != ':')) {
				pfbFile = ext.pwd() + pfbFile;
			}
			if (!Files.exists(pfbFile)) {
				log.reportError("Error - pfb file '" + pfbFile + "' does not exist; aborting");
				return;
			}
		}

		if (gcmodelFile != null) {
			gcmodelFile = ext.replaceTilde(gcmodelFile);
			if (!gcmodelFile.startsWith("/") && (gcmodelFile.charAt(1) != ':')) {
				gcmodelFile = ext.pwd() + gcmodelFile;
			}
			if (!Files.exists(gcmodelFile)) {
				log.reportError("Error - gcmodel file '" + gcmodelFile + "' does not exist; aborting");
				return;
			}
		}

		new File(resultsDir).mkdirs();
		new File(dataDir).mkdirs();

		files = new File(dataDir).list(new FilenameFilter() {
			@Override
			public boolean accept(File file, String filename) {
				return file.length() > 1000 && !filename.equals("chrX") && !filename.equals("sexSpecific");
			}
		});

		if (files == null || files.length == 0) {
			log.reportError("Found zero files in " + dataDir);
			log.reportError("Will not proceed");
			return;
		}
		log.report("Found " + files.length + " files in " + dataDir);

		step = (int) Math.ceil((double) files.length / (double) numChunks);
		log.report("Which means the step for " + numChunks + " chunks would be " + step);

		for (int i = 0; i < numChunks; i++) {
			try {
				writer = new PrintWriter(new FileWriter(resultsDir + "list" + (i + 1) + ".txt"));
				for (int j = i * step; j < Math.min(files.length, (i + 1) * step); j++) {
					if (files[j].endsWith(".gz")) {
						writer.println("`gunzip -c " + dataDir + files[j] + "`");
					} else {
						writer.println(dataDir + files[j]);
					}
				}
				writer.close();
			} catch (Exception e) {
				log.reportError("Error writing to " + resultsSubDir + "list" + (i + 1) + ".txt");
				log.reportException(e);
			}
		}

		commands = execDir + "detect_cnv.pl -test -conf -hmm "
		           + (hmmFile == null ? execDir + "lib/hhall.hmm" : hmmFile) + " -pfb "
		           + (pfbFile == null ? execDir + "lib/hhall.hg18.pfb" : pfbFile) + " -gcmodel "
		           + (gcmodelFile == null ? execDir + "lib/hhall.hg18.gcmodel" : gcmodelFile)
		           + " -list " + resultsDir + "list[%0].txt -log " + resultsDir + "[%0].log -out "
		           + resultsDir + "[%0].rawcnv > " + resultsDir + "[%0].out";

		new File(pennDir + scriptSubDir).mkdirs();

		if (execList == null) {
			Files.qsub(pennDir + scriptSubDir + "runPenn", dataDir, numChunks, commands,
			           Matrix.toMatrix(Array.stringArraySequence(numChunks, "")), 2200, 16);
		} else {
			Files.execListAdd(execList, commands, Array.stringArraySequence(numChunks, ""), log);
		}

		Files.writeArray(new String[] {"cd " + projDir,
		                               "cat " + resultsDir + "*.log > " + resultsDir + "penncnv.rawlog",
		                               "cat " + resultsDir + "*.rawcnv > " + resultsDir + "penncnv.rawcnv",
		                               runLine + " rawlog=" + resultsDir + "penncnv.rawlog",
		                               runLine + " rawcnv=" + resultsDir + "penncnv.rawcnv",},
		                 pennDir + scriptSubDir + "assemblePenncnv");
		Files.chmod(pennDir + scriptSubDir + "assemblePenncnv");
	}

	// FIXME need to unify this method with batch
	public static void batchX(Project proj, int numChunks, Vector<String> execList, String pfbFile,
	                          String gcmodelFile, String hmmFile, String scriptSubDir,
	                          String dataSubDir, String resultsSubDir) {
		String commands;
		PrintWriter writer;
		String[] files;
		int step;
		String execDir, dataDir, resultsDir, projDir, pennDir;
		Logger log;

		log = proj.getLog();
		final String runLine = Files.getRunString() + " " + PennCNV.class.getCanonicalName() + " proj="
		                       + new File(proj.getPropertyFilename()).getAbsolutePath();
		projDir = proj.PROJECT_DIRECTORY.getValue();
		execDir = proj.PENNCNV_EXECUTABLE_DIRECTORY.getValue(false, true);
		pennDir = proj.PENNCNV_RESULTS_DIRECTORY.getValue(false, true);
		dataDir = proj.PENNCNV_DATA_DIRECTORY.getValue(false, true) + dataSubDir;
		resultsDir = pennDir + resultsSubDir;

		// if (!Files.exists(proj.getFilename("SAMPLE_DATA_FILENAME", false, false),
		// proj.getJarStatus())) {
		if (!Files.exists(proj.SAMPLE_DATA_FILENAME.getValue(false, false),
		                  proj.JAR_STATUS.getValue())) {
			log.reportError("Error - sample data file " + proj.SAMPLE_DATA_FILENAME.getValue()
			                + " does not exist;");
			return;
		}

		SampleData sampleData = proj.getSampleData(2, false);

		if (sampleData.failedToLoad()) {
			log.reportError("Error - without a sample data file, PennCNV will fail to analyze sex chromosomes");
			return;
		}

		if (pfbFile != null) {
			pfbFile = ext.replaceTilde(pfbFile);
			if (!pfbFile.startsWith("/") && (pfbFile.charAt(1) != ':')) {
				pfbFile = ext.pwd() + pfbFile;
			}
			if (!Files.exists(pfbFile)) {
				log.reportError("Error - pfb file '" + pfbFile + "' does not exist; aborting");
				return;
			}
		}

		if (gcmodelFile != null) {
			gcmodelFile = ext.replaceTilde(gcmodelFile);
			if (!gcmodelFile.startsWith("/") && (gcmodelFile.charAt(1) != ':')) {
				gcmodelFile = ext.pwd() + gcmodelFile;
			}
			if (!Files.exists(gcmodelFile)) {
				log.reportError("Error - gcmodel file '" + gcmodelFile + "' does not exist; aborting");
				return;
			}
		}

		new File(resultsDir).mkdirs();
		new File(dataDir).mkdirs();

		String newGCFile = dataDir + "chrX.gcModel";
		String newPFBFile = dataDir + "chrX.pfb";

		pfbFile = filterFile(proj, pfbFile, newPFBFile, new String[] {"23", "X"});
		gcmodelFile = filterFile(proj, gcmodelFile, newGCFile, new String[] {"23", "X"});

		String sexFileStatus = writeSexFile(proj, sampleData, dataDir, log);
		if (sexFileStatus.length() > 0) {
			log.reportError("Error - " + sexFileStatus);
			return;
		}

		files = new File(dataDir).list(new FilenameFilter() {
			@Override
			public boolean accept(File file, String filename) {
				return file.length() > 1000 && !filename.endsWith(".pfb") && !filename.endsWith(".gcmodel") && !filename.startsWith("sex_file");
			}
		});
		log.report("Found " + files.length + " files");


		if (files == null || files.length == 0) {
			log.reportError("Found zero files in " + dataDir);
			log.reportError("Will not proceed");
			return;
		}
		step = (int) Math.ceil((double) files.length / (double) numChunks);
		log.report("Which means the step for " + numChunks + " chunks would be " + step);

		for (int i = 0; i < numChunks; i++) {
			try {
				writer = new PrintWriter(new FileWriter(resultsDir + "list" + (i + 1) + ".txt"));
				for (int j = i * step; j < Math.min(files.length, (i + 1) * step); j++) {
					if (files[j].endsWith(".gz")) {
						writer.println("`gunzip -c " + dataDir + files[j] + "`");
					} else {
						writer.println(files[j]);
					}
				}
				writer.close();
			} catch (Exception e) {
				log.reportError("Error writing to " + resultsSubDir + "list" + (i + 1) + ".txt");
				log.reportException(e);
			}
		}

		commands = execDir + "detect_cnv.pl -test -conf -hmm "
		           + (hmmFile == null ? execDir + "lib/hhall.hmm" : hmmFile) + " -pfb "
		           + (pfbFile == null ? execDir + "lib/hhall.hg18.pfb" : pfbFile) + " -gcmodel "
		           + (gcmodelFile == null ? execDir + "lib/hhall.hg18.gcmodel" : gcmodelFile)
		           + " -chrx -sexfile " + dataDir + "sex_file.txt -list " + resultsDir
		           + "list[%0].txt -log " + resultsDir + "[%0].log -out " + resultsDir
		           + "[%0].rawcnv > " + resultsDir + "[%0].out";

		new File(pennDir + scriptSubDir).mkdirs();

		if (execList == null) {
			Files.qsub(pennDir + scriptSubDir + "runPennX", dataDir, numChunks, commands,
			           Matrix.toMatrix(Array.stringArraySequence(numChunks, "")), 2200, 16);
		} else {
			Files.execListAdd(execList, commands, Array.stringArraySequence(numChunks, ""), log);
		}
		Files.writeArray(new String[] {"cd " + projDir,
		                               "cat " + resultsDir + "*.log > " + resultsDir
		                                                + "penncnvX.rawlog",
		                               "cat " + resultsDir + "*.rawcnv > " + resultsDir + "penncnvX.rawcnv",
		                               // don't parse warnings; the parseWarnings method isn't written
		                               // to
		                               // parse X-chromosome warnings
		                               runLine + " rawcnv=" + resultsDir + "penncnvX.rawcnv",},
		                 pennDir + scriptSubDir + "assemblePenncnv");
		Files.chmod(pennDir + scriptSubDir + "assemblePenncnv");
	}

	private static String filterFile(Project proj, String fileToFilter, String outputFile,
	                                 String[] chrs) {
		// TODO combine method with filterPFB - literally the same except different names/extensions
		BufferedReader reader = null;
		PrintWriter writer = null;

		try {
			reader = new BufferedReader(new FileReader(fileToFilter));
			writer = new PrintWriter(new FileWriter(outputFile));

			String header;
			String temp;
			String[] line;
			if (reader.ready()) {
				header = reader.readLine();
				writer.println(header);
			}
			while (reader.ready()) {
				temp = reader.readLine();
				line = temp.trim().split("[\\s]+");
				for (String chr : chrs) {
					if (line[1].equals(chr)) {
						writer.println(temp);
					}
				}
			}
		} catch (IOException e) {
			proj.getLog().reportError("Error - filtering failed for file: " + fileToFilter);
			proj.getLog().reportException(e);
			return fileToFilter;
		} finally {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					proj.getLog()
					    .reportError("Error - couldn't properly close file reader for " + fileToFilter);
					proj.getLog().reportException(e);
				}
				reader = null;
			}
			if (writer != null) {
				writer.close();
				writer = null;
			}
		}

		return outputFile;
	}

	private static String writeSexFile(Project proj, SampleData sampleData, String resultsDir,
	                                   Logger log) {
		String sampleDataFile = proj.SAMPLE_DATA_FILENAME.getValue(false, false);
		String[] header = Files.getHeaderOfFile(sampleDataFile, proj.getLog());
		int sexInd = -1;
		for (int i = 0; i < header.length; i++) {
			if (("CLASS=" + SexChecks.EST_SEX_HEADER).equalsIgnoreCase(header[i])) {
				sexInd = i;
				break;
			}
		}
		if (sexInd == -1) {
			return "no estimated sex found in sample data file - please run SexChecks with -check argument to generate the required data";
		}
		Hashtable<String, Vector<String>> sexData = HashVec.loadFileToHashVec(sampleDataFile, 0,
		                                                                      new int[] {sexInd}, "\t",
		                                                                      true, false);
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(resultsDir + "sex_file.txt"));
			for (Map.Entry<String, Vector<String>> lineData : sexData.entrySet()) {
				String estSexStr = lineData.getValue().get(0);
				if (!ext.isMissingValue(estSexStr)) {
					int estSex = Integer.parseInt(estSexStr);
					estSex = SexChecks.KARYOTYPES[estSex].contains("XX") ? 2 : 1;
					writer.println(lineData.getKey() + "\t" + estSex);
				}
			}
			writer.close();
		} catch (IOException e) {
			log.reportException(e);
			return "unable to complete writing of sex_file for PennCNV";
		}
		return "";
	}

	public static void parseWarnings(Project proj, String filename) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, data;
		String temp, sampleID = null;
		Hashtable<String, String[]> hash = new Hashtable<String, String[]>();
		Vector<String> v = new Vector<String>();
		SampleData sampleData;
		int err;
		double lrrSD_cutoff;
		String[] ids;
		long time;
		Logger log;

		time = new Date().getTime();
		log = proj.getLog();
		log.report("Parsing PennCNV warning...");

		sampleData = proj.getSampleData(2, false);
		// lrrSD_cutoff = proj.getDouble(proj.LRRSD_CUTOFF);
		lrrSD_cutoff = proj.LRRSD_CUTOFF.getValue();

		try {
			reader = new BufferedReader(new FileReader(filename));
			while (reader.ready()) {
				temp = reader.readLine();
				temp = translateDerivedSamples(temp, log);
				line = temp.trim().split("[\\s]+");
				try {
					if (temp.contains("quality summary")) {
						sampleID = line[4].substring(line[4].lastIndexOf("/") + 1, line[4].indexOf(":"));
						v.add(sampleID);
						data = new String[QC_HEADS.length + ERRORS.length];
						for (int i = 0; i < ERRORS.length; i++) {
							data[i] = ".";
						}
						if (line.length < QC_HEADS.length + 5) {
							log.reportError("Error - line doesn't have all the expected pieces:");
							log.reportError(temp);
						}
						for (int i = 0; i < QC_HEADS.length; i++) {
							data[ERRORS.length + i] = line[5 + i].split("=")[1];
						}
						hash.put(sampleID, data);
					} else if (temp.startsWith("WARNING")) {
						if (temp.contains("Small-sized CNV calls")) {
							// use old trav
						} else {
							sampleID = line[3].substring(line[3].lastIndexOf("/") + 1);
						}
						data = hash.get(sampleID);
						err = -1;
						for (int i = 0; i < ERRORS.length; i++) {
							if (temp.contains(ERRORS[i])) {
								err = i;
							}
						}
						if (err == -1) {
							log.reportError("Unknown WARNING: " + temp);
						} else {
							data[err] = "1";
						}
					}
				} catch (Exception e) {
					log.reportError("Error with: " + temp);
					log.reportException(e);
				}
			}
			reader.close();

			writer = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue() + QC_SUMMARY_FILE));
			writer.print("Sample\tFID\tIID\tUse_" + ext.formDeci(lrrSD_cutoff, 2));
			for (String element : ERRORS) {
				writer.print("\t" + element);
			}
			for (String element : QC_HEADS) {
				writer.print("\t" + element);
			}
			writer.println();
			Collections.sort(v);
			for (int i = 0; i < v.size(); i++) {
				sampleID = v.get(i);
				data = hash.get(sampleID);
				ids = sampleData.lookup(sampleID);
				writer.print(sampleID + "\t" + (ids == null ? "NotInSampleData\t" + sampleID : ids[1]));
				writer.print("\t" + (data[1].equals("1") || data[2].equals("1")
				                     || Double.parseDouble(data[6]) > lrrSD_cutoff ? "0" : "1"));
				writer.println("\t" + Array.toStr(data));
			}
			writer.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + filename + "\" not found in current directory");
			return;
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + filename + "\"");
			return;
		}
		log.report("Parsed PennCNV warnings in " + ext.getTimeElapsed(time));
	}

	public static String translateDerivedSamples(String str, Logger log) {
		String trav;
		int start, stop;

		start = str.indexOf("`");
		stop = str.lastIndexOf("`");
		if (start == -1) {
			return str;
		}

		trav = str.substring(start + 1, stop);
		if (trav.contains("`")) {
			log.reportError("Error - more than one set of quotes for: " + str);
		}
		if (trav.startsWith("gunzip -c") && trav.endsWith(".gz")) {
			trav = trav.substring(9, trav.length() - 3).trim();
		} else {
			log.reportError("Error - not currently set up to handle the following construction into a sample_ID: "
			                + trav);
		}

		return str.substring(0, start) + trav + str.substring(stop + 1);
	}

	public static void combineResults(Project proj, String[] cnvFiles, String outputFile,
	                                  boolean recode) {
		BufferedReader reader;
		PrintWriter writer;
		Logger log = proj.getLog();

		// TODO check input and output file names for .cnv extension( - error if not? or just
		// warning...?)

		java.util.HashMap<String, java.util.TreeMap<String, java.util.ArrayList<String[]>>> cnvSet =
		                                                                                           new HashMap<String, TreeMap<String, ArrayList<String[]>>>();

		String temp;
		String[] line;
		String key, chr, currFile = null;
		boolean readAll = false;
		try {
			for (String cnvFile : cnvFiles) {
				currFile = cnvFile;
				reader = new BufferedReader(new FileReader(cnvFile));
				if (reader.ready()) {
					// skip header
					reader.readLine();
				}
				while (reader.ready()) {
					temp = reader.readLine();
					line = temp.split("\t");
					key = line[0] + "\t" + line[1];
					// get all CNVs for an individual:
					TreeMap<String, ArrayList<String[]>> chrSets = cnvSet.get(key);
					if (chrSets == null) {
						chrSets = new TreeMap<String, ArrayList<String[]>>();
						cnvSet.put(key, chrSets);
					}
					chr = line[2];
					// get all CNVs for a specific chromosome
					ArrayList<String[]> chrSet = chrSets.get(chr);
					if (chrSet == null) {
						chrSet = new ArrayList<String[]>() {
							private static final long serialVersionUID = 1L;

							@Override
							public boolean add(String[] e) {
								int index = Array.binarySearch(this, e, 0, false);
								super.add(index, e);
								return true;
							}
						};
						chrSets.put(chr, chrSet);
					}
					// add CNV to list
					chrSet.add(Array.subArray(line, 3));
				}
				reader.close();
			}
			readAll = true;
		} catch (FileNotFoundException e) {
			log.reportError("Error: file \"" + currFile + "\" not found in current directory");
			log.reportException(e);
		} catch (IOException e) {
			log.reportException(e);
		}

		if (readAll) {
			try {
				writer = new PrintWriter(new FileWriter(outputFile));
				writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER));
				String FIDIID;
				String cnvChr;
				for (Map.Entry<String, TreeMap<String, ArrayList<String[]>>> sample : cnvSet.entrySet()) {
					FIDIID = sample.getKey();
					for (Map.Entry<String, ArrayList<String[]>> chrSet : sample.getValue().entrySet()) {
						cnvChr = chrSet.getKey();
						if (recode) {
							if ("1".equals(cnvChr)) {
								cnvChr = "23";
							} else if ("2".equals(cnvChr)) {
								cnvChr = "24";
							} else if ("3".equals(cnvChr)) {
								cnvChr = "25";
							} else if ("4".equals(cnvChr)) {
								cnvChr = "26";
							}
						}
						for (String[] cnv : chrSet.getValue()) {
							writer.println(FIDIID + "\t" + cnvChr + "\t" + Array.toStr(cnv, "\t"));
						}
					}
				}
				writer.close();
			} catch (IOException e) {
				log.reportException(e);
			}
		}
	}

	public static void parseResults(Project proj, String filename, boolean denovoOnly) {
		BufferedReader reader;
		PrintWriter writer;
		PrintWriter denovoWriter = null;
		String[] line;
		String temp, trav;
		Vector<String> warnings;
		Hashtable<String, Vector<String>> pedinfo;
		int[] position;
		String score;
		SampleData sampleData;
		String famIndPair;
		Hashtable<String, String> hash;
		String[] ids;
		List<String> inds;
		String[] fams;
		long time;
		int sex;
		Logger log;

		log = proj.getLog();
		log.report("Parsing PennCNV rawcnvs...");
		time = new Date().getTime();

		if (!Files.exists(filename)) {
			log.reportError("Error - could not find file '" + filename + "'");
			return;
		}

		warnings = new Vector<String>();
		sampleData = proj.getSampleData(2, false);
		pedinfo = new Hashtable<String, Vector<String>>();
		Pedigree ped = proj.loadPedigree();
		PrintWriter[] denoValWriter = new PrintWriter[1];
		try {
			reader = new BufferedReader(new FileReader(filename));
			writer = new PrintWriter(new FileWriter(ext.rootOf(filename, false) + ".cnv"));
			writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER));
			hash = new Hashtable<String, String>();
			while (reader.ready()) {
				temp = reader.readLine();
				if (!temp.startsWith("NOTICE:")) {
					temp = translateDerivedSamples(temp, log);
					line = temp.trim().split("[\\s]+");
					position = Positions.parseUCSClocation(line[0]);
					trav = line[4];
					trav = trav.substring(trav.lastIndexOf("/") + 1);
					ids = sampleData.lookup(trav);
					if (ids == null) {
						if (!hash.containsKey(trav)) {
							// log.reportError("Error - '"+trav+"' was not found in
							// "+proj.getFilename(proj.SAMPLE_DATA_FILENAME));
							log.reportError("Error - '" + trav + "' was not found in "
							                + proj.SAMPLE_DATA_FILENAME.getValue());
							hash.put(trav, "");
						}
						famIndPair = trav + "\t" + trav;
					} else {
						famIndPair = ids[1];
					}

					ids = famIndPair.split("\t");
					HashVec.addToHashVec(pedinfo, ids[0], ids[1], true);

					if (line.length < 8 || !line[7].startsWith("conf=")
					    || line[7].toUpperCase().contains("NAN")) {
						score = Integer.toString(MISSING_SCORE);
						if (!warnings.contains(trav) && warnings.size() < 10) {
							log.reportError("Warning - no conf estimates for " + trav);
							warnings.add(trav);
						}
					} else {
						score = ext.formDeci(Double.parseDouble(line[7].substring(5)), 4, true);
					}
					boolean isDenovo = false;
					for (String s : line) {
						if (s.startsWith("statepath=33") || s.startsWith("triostate=33")) {
							isDenovo = true;
						}
					}

					String copynum = line[3].substring(line[3].indexOf("=") + 1);
					String sites = line[1].substring(7);
					StringBuilder lineOut = new StringBuilder(famIndPair).append("\t").append(position[0])
					                                                     .append("\t").append(position[1])
					                                                     .append("\t").append(position[2])
					                                                     .append("\t").append(copynum)
					                                                     .append("\t").append(score)
					                                                     .append("\t").append(sites);
					if (!denovoOnly) {
						writer.println(lineOut.toString());
					}
					if (isDenovo) {
						if (denovoWriter == null) {
							denovoWriter = new PrintWriter(new FileWriter(ext.rootOf(filename, false)
							                                              + "_denovo.cnv"));
							denovoWriter.println(Array.toStr(CNVariant.PLINK_CNV_HEADER));
						}
						denovoWriter.println(lineOut.toString());
						writeValidation(ped, ids, position, copynum, line, filename, denoValWriter, log);
					}
				}
			}
			reader.close();
			writer.close();
			if (denovoWriter != null) {
				denovoWriter.close();
			}
			if (denoValWriter[0] != null) {
				denoValWriter[0].close();
			}

			// FilterCalls.stdFilters(dir, ext.rootOf(filename)+".cnv", MAKE_UCSC_TRACKS);

			writer = new PrintWriter(new FileWriter(ext.rootOf(filename, false) + ".fam"));
			fams = HashVec.getNumericKeys(pedinfo);
			for (String fam : fams) {
				inds = pedinfo.get(fam);
				Collections.sort(inds, new SciStringComparator());
				for (String ind : inds) {
					ids = sampleData.lookup(fam + "\t" + ind);
					if (ids != null) {
						sex = sampleData.getSexForIndividual(ids[0]);
					} else {
						sex = 0;
					}
					int pedIndex = ped == null ? -1 : ped.getIndIndex(fam, ind);
					String fa = pedIndex >= 0 && ped != null ? ped.getFA(pedIndex) : "0";
					String mo = pedIndex >= 0 && ped != null ? ped.getMO(pedIndex) : "0";
					writer.println(fam + "\t" + ind + "\t" + fa + "\t" + mo + "\t" + Math.max(0, sex)
					               + "\t-9");
				}
			}
			writer.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + ext.rootOf(filename, false)
			                + ".cnv\" not found in current directory");
			return;
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + ext.rootOf(filename, false) + "\"");
			return;
		}

		log.report("...finished in " + ext.getTimeElapsed(time));
	}

	private static void writeValidation(Pedigree ped, String[] ids, int[] position, String copynum,
	                                    String[] line, String filename, PrintWriter[] denoValWriter,
	                                    Logger log) {
		int pedIndex = ped.getIndIndex(ids[0], ids[1]);
		if (pedIndex < 0) {
			return;
		}
		int faIndex = ped.getIndexOfFaInIDs(pedIndex);
		int moIndex = ped.getIndexOfMoInIDs(pedIndex);
		if (faIndex < 0 || moIndex < 0) {
			return;
		}

		String cDna = ped.getiDNA(pedIndex);
		String faDna = ped.getiDNA(faIndex);
		String moDna = ped.getiDNA(moIndex);

		String outDir = ext.parseDirectoryOfFile(filename);

		if (denoValWriter[0] == null) {
			try {
				denoValWriter[0] = new PrintWriter(new FileWriter(outDir + "denovoValidation.txt"));
				denoValWriter[0].println("export HMMFILE=");
				denoValWriter[0].println("export PFBFILE=");
				denoValWriter[0].println();
			} catch (IOException e) {
				log.reportException(e);
				return;
			}
		}

		String[] bounds = getBounds(line);
		if (bounds[0] == null || bounds[1] == null) {
			return;
		}

		outDir += "denovo" + File.separator;
		File outFile = new File(outDir);
		if (!outFile.exists()) {
			outFile.mkdirs();
		}

		String childSource = "gunzip -c " + line[4] + ".gz";
		String faSource = childSource.replace(cDna, faDna);
		String moSource = childSource.replace(cDna, moDna);

		if (childSource.contains("sexSpecific")) {
			faSource = faSource.replaceAll("/female/", "/male/");
			moSource = moSource.replaceAll("/male/", "/female/");
		}

		String out = outDir + ids[0] + "_" + ids[1] + "_" + position[0] + "_" + position[1] + "_" + position[2];
		String faFile = out + "_fa.txt";
		String moFile = out + "_mo.txt";
		String childFile = out + "_off.txt";

		StringBuilder extractLine = new StringBuilder(faSource).append(" > ").append(faFile)
		                                                       .append(" && ").append(moSource)
		                                                       .append(" > ").append(moFile)
		                                                       .append(" && ").append(childSource)
		                                                       .append(" > ").append(childFile);

		denoValWriter[0].println(extractLine.toString());

		StringBuilder sb =
		                 new StringBuilder("/home/pankrat2/shared/bin/infer_snp_allele.pl -pfbfile $PFBFILE -hmmfile $HMMFILE").append(" -denovocn ")
		                                                                                                                       .append(copynum)
		                                                                                                                       .append(" -startsnp ")
		                                                                                                                       .append(bounds[0])
		                                                                                                                       .append(" -endsnp ")
		                                                                                                                       .append(bounds[1])
		                                                                                                                       .append(" -outfile ")
		                                                                                                                       .append(out)
		                                                                                                                       .append(".gen  -logfile ")
		                                                                                                                       .append(out)
		                                                                                                                       .append(".log ")
		                                                                                                                       .append(faFile)
		                                                                                                                       .append(" ")
		                                                                                                                       .append(moFile)
		                                                                                                                       .append(" ")
		                                                                                                                       .append(childFile);

		denoValWriter[0].println(sb.toString());
		StringBuilder cleanup = new StringBuilder("rm ").append(faFile).append(" && rm ").append(moFile)
		                                                .append(" && rm ").append(childFile);

		denoValWriter[0].println(cleanup.toString());
	}

	private static String[] getBounds(String[] line) {
		String[] bounds = new String[2];
		for (String s : line) {
			if (s.startsWith("startsnp")) {
				bounds[0] = s.split("=")[1];
			} else if (s.startsWith("endsnp")) {
				bounds[1] = s.split("=")[1];
			}
		}
		return bounds;
	}

	// Available in cnv.Launch
	// public static void fromParameters(String filename, Logger log) {
	// Vector<String> params;
	//
	// params = Files.parseControlFile(filename, "penncnv", new String[]
	// {"proj=/home/npankrat/projects/GEDI.properties", "rawcnv=all.rawcnv", "rawlog=all.log"}, log);
	//
	// if (params != null) {
	// params.add("log=" + log.getFilename());
	// main(Array.toStringArray(params));
	// }
	// }

	/**
	 * Calculate the population BAF (B Allele Frequency based on all the samples available in the)
	 * data. Output is going to be saved on disk. In PennCnv, this file is also called snpFile.
	 *
	 * @param proj The project you are going to run PennCNV on.
	 *
	 *        The output file looks like the the following: Name Chr Position PFB rs1000113 5
	 *        150220269 0.564615751221256 rs1000115 9 112834321 0.565931333264192 rs10001190 4 6335534
	 *        0.5668604380025 rs10002186 4 38517993 0.57141752993563 rs10002743 4 6327482
	 *        0.567557695424774
	 *
	 */
	public static void populationBAF(Project proj) {
		String[] sampleList;
		String output;

		Logger log = proj.getLog();
		String filename = proj.SAMPLE_SUBSET_FILENAME.getValue(true, false);

		if (ext.rootOf(filename) == null || ext.rootOf(filename).equals("")
		    || !Files.exists(filename, proj.JAR_STATUS.getValue())) {
			sampleList = proj.getSampleList().getSamples();
			output = proj.CUSTOM_PFB_FILENAME.getValue(true, false);
		} else if (Files.exists(filename, proj.JAR_STATUS.getValue())) {
			log.report("filename: " + filename);
			sampleList = HashVec.loadFileToStringArray(filename, false, new int[] {0}, false);
			output = proj.PROJECT_DIRECTORY.getValue() + ext.rootOf(filename) + ".pfb";
		} else {
			proj.message("Failed to load \"" + filename + "\"");
			return;
		}

		MarkerSet markerSet = proj.getMarkerSet();
		String[] markerNames = markerSet.getMarkerNames();
		byte[] chrs = markerSet.getChrs();
		int[] positions = markerSet.getPositions();
		double[] bafSum = new double[chrs.length];
		int[] bafCounts = new int[chrs.length];
		int[] genoCounts = new int[chrs.length];
		for (int i = 0; i < sampleList.length; i++) {
			if (i % 100 == 0) {
				log.report("Loading file " + (i + 1) + " of " + sampleList.length);
			}
			Sample samp = proj.getPartialSampleFromRandomAccessFile(sampleList[i], false, false, true,
			                                                        false, true);
			float[] bafs = samp.getBAFs();
			byte[] genotypes = samp.getAB_Genotypes();
			for (int j = 0; j < bafSum.length; j++) {
				if (!Float.isNaN(bafs[j])) {
					bafSum[j] += bafs[j];
					bafCounts[j]++;
					if (genotypes[j] >= 0) {
						genoCounts[j]++;
					}
				}
			}
		}
		double[] bafAverage = new double[chrs.length];
		ArrayList<String> missingGenotypeMarkers = new ArrayList<String>();
		for (int i = 0; i < bafSum.length; i++) {
			boolean cnOnly = proj.getArrayType().isCNOnly(markerNames[i]);
			if (genoCounts[i] != 0 && !cnOnly) {// Since mock genotypes can be present, we demand non-CN
			                                    // only
				bafAverage[i] = bafSum[i] / bafCounts[i];
			} else if (cnOnly) {
				bafAverage[i] = 2;
			} else {
				bafAverage[i] = -1; // This is to more clearly differentiate CN only markers from SNPs
				                    // without callrate

				missingGenotypeMarkers.add(markerNames[i]);
			}
		}


		PSF.checkInterrupted();
		try {

			PrintWriter writer = new PrintWriter(new FileWriter(output));
			writer.println("Name\tChr\tPosition\tPFB");
			for (int i = 0; i < markerNames.length; i++) {
				writer.println(markerNames[i] + "\t"
				               + (chrs[i] < 23 ? chrs[i]
				                               : (chrs[i] == 23 ? "X"
				                                                : (chrs[i] == 24 ? "Y"
				                                                                 : (chrs[i] == 25 ? "XY"
				                                                                                  : (chrs[i] == 26 ? "M"
				                                                                                                   : "Un")))))
				               + "\t" + positions[i] + "\t" + bafAverage[i]);
			}
			writer.close();
			log.report("Population BAF file is now ready at: " + output);
			if (!missingGenotypeMarkers.isEmpty()) {
				String missingGenoFile = ext.addToRoot(output, ".missingGenotypes");
				log.reportTimeInfo(missingGenotypeMarkers.size()
				                   + " markers had missing genotypes and were set to -1 in " + output
				                   + ". These markers can be treated as CN only markers, or removed at your discretion with CNVCaller");
				Files.writeIterable(missingGenotypeMarkers, missingGenoFile);
			}
		} catch (Exception e) {
			log.reportError("Error writing to '" + output + "'");
			log.reportException(e);
		}
	}


	/**
	 * Generate the GCModel file needed by PennCNV software
	 * (http://www.openbioinformatics.org/penncnv/).
	 *
	 * @param proj The project you are going to run PennCNV on.
	 * @param inputGcBaseFullPath The user-supplied genome builds. Positions within each chromosome
	 *        must be sorted by increasing order. For example,
	 *        http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/gc5Base.txt.gz
	 * @param outputGcModelFullPath The name of the GCModel file.
	 * @param numwindow For each SNP, GC Content is calculated for the range of numwindow*5120 before
	 *        its location through numwindow*5120 after. To use the default setting of 100, please
	 *        enter 0.
	 *
	 *        In order to be able to cross-reference with the same feature in PennCNV, this code
	 *        intends to base on cal_gc_snp.pl in PennCNV package. But since the difference between
	 *        Java and Perl, significant structural changes have been made.
	 *
	 *        A sample gcBase file looks like below (There is no header line): 585 chr1 0 5120 chr1.0
	 *        5 1024 0 /gbdb/hg18/wib/gc5Base.wib 0 100 1024 59840 3942400 585 chr1 5120 10240 chr1.1
	 *        5 1024 1024 /gbdb/hg18/wib/gc5Base.wib 0 100 1024 59900 3904400 585 chr1 10240 15360
	 *        chr1.2 5 1024 2048 /gbdb/hg18/wib/gc5Base.wib 0 100 1024 55120 3411200 585 chr1 15360
	 *        20480 chr1.3 5 1024 3072 /gbdb/hg18/wib/gc5Base.wib 0 100 1024 49900 3078800 585 chr1
	 *        20480 25600 chr1.4 5 1024 4096 /gbdb/hg18/wib/gc5Base.wib 0 100 1024 47600 2682400
	 *
	 *        A sample gcModel file (the output) look like this: Name Chr Position GC rs4961 4 2876505
	 *        48.6531211131841 rs3025091 11 102219838 38.1080923507463 rs3092963 3 46371942
	 *        44.4687694340796 rs3825776 15 56534122 40.7894123134328 rs17548390 6 3030512
	 *        45.3604050062189 rs2276302 11 113355350 43.8200598569652
	 *
	 *        snpFile or pfb file (Population B Allele Frequency) is an additional data file that is
	 *        required by the corresponding feature in PennCNV but not here. It looks like the this:
	 *        Name Chr Position PFB rs1000113 5 150220269 0.564615751221256 rs1000115 9 112834321
	 *        0.565931333264192 rs10001190 4 6335534 0.5668604380025 rs10002186 4 38517993
	 *        0.57141752993563 rs10002743 4 6327482 0.567557695424774
	 *
	 */
	public static void gcModel(Project proj, String inputGcBaseFullPath, String outputGcModelFullPath,
	                           int numwindow) {
		MarkerSet markerSet;
		String[] markerNames;
		byte[] chrs;
		int[] positions;

		BufferedReader reader;
		String[] line;
		PrintWriter writer;
		byte curchr;
		int curstart, curend, curcount;
		boolean skip = false;
		double cursum;
		Hashtable<Byte, Byte> seen_chr = new Hashtable<Byte, Byte>();

		int[] snp_count;
		double[] snp_sum;
		byte prechr = -1;
		int prestart = -1;
		int chr_index;
		Logger log;

		log = proj.getLog();

		// generate or load SnpFile (pbf file or Population B Allele Frequency)
		markerSet = proj.getMarkerSet();
		markerNames = markerSet.getMarkerNames();
		chrs = markerSet.getChrs(); // to be used only in the SnpFile (.pdf file) block. But then in the
		                            // outputFile.
		positions = markerSet.getPositions(); // to be used only in the SnpFile (.pdf file) block. But
		                                      // then in the GcFile block.

		// How do we know whether "numwindow==null" ???
		if (numwindow == 0) {
			numwindow = 100;
		}

		snp_count = new int[markerNames.length];
		snp_sum = new double[markerNames.length];
		chr_index = 0;

		// load gcFile
		try {

			PSF.checkInterrupted();

			// If an invalid gc base path was given, try using the resources
			if (!new File(inputGcBaseFullPath).exists()) {
				Resource r = Resources.genome(proj.GENOME_BUILD_VERSION.getValue(), log).getModelBase();
				if (r.isAvailable()) {
					inputGcBaseFullPath = r.getAbsolute();
				}
			}
			reader = Files.getAppropriateReader(inputGcBaseFullPath);
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				curchr = Positions.chromosomeNumber(line[1]);
				curstart = Integer.parseInt(line[2]);
				curend = Integer.parseInt(line[3]);
				curcount = Integer.parseInt(line[11]);
				cursum = Double.parseDouble(line[12]);

				// For each record in gcFile, scan snpFile for the range of +/- numwindow*5120
				if (curchr > 0) {
					if (curchr == prechr) {
						if (curstart < prestart) {
							log.reportError("Error in gcFile: a record in chr" + curchr + " has position "
							                + curstart + ", less then the previous position $prestart");
							reader.close();
							return;
						}
					} else if (seen_chr.containsKey(curchr)) {
						log.reportError("Error in gcFile: rows of the same chromosome must be adjacent. But now chr"
						                + curchr + " occur multiple times in non-continuous segment of the "
						                + inputGcBaseFullPath + ": at " + curchr + ":" + curstart);
						reader.close();
						return;
					} else {
						seen_chr.put(curchr, (byte) 1);

						skip = false;
						if (chrs[chr_index] > curchr) {
							chr_index = 0;
						}
						while (chr_index < (chrs.length) && chrs[chr_index] != curchr) {
							chr_index++;
						}
					}

					if (!(curchr == prechr && skip)) {
						while (chr_index < chrs.length && chrs[chr_index] == curchr
						       && positions[chr_index] < Math.max(curstart - numwindow * 5120, 0)) {
							chr_index++;
						}
						if (chr_index >= chrs.length) {
							chr_index = 0;
							skip = true;
						} else if (chrs[chr_index] != curchr) {
							skip = true;
						} else {
							for (int i = chr_index; i < snp_count.length && chrs[i] == curchr
							                        && positions[i] <= (curend + numwindow * 5120); i++) {
								snp_count[i] += curcount;
								snp_sum[i] += cursum;
							}
						}
					}

					prestart = curstart;
					prechr = curchr;
				}
			}
			reader.close();
		} catch (Exception e) {
			log.reportError("Error reading from '" + inputGcBaseFullPath + "'");
			log.reportException(e);
			return;
		}


		PSF.checkInterrupted();
		// load pfb file or generate it

		// output the result
		try {
			writer = new PrintWriter(new FileWriter(outputGcModelFullPath));
			writer.println("Name\tChr\tPosition\tGC");
			for (int i = 0; i < markerNames.length; i++) {
				// writer.println(markerNames[i]+"\t"+chrs[i]+"\t"+positions[i]+"\t"+(snp_count[i]==0?(snp_sum[i]==0?0:"err"):(snp_sum[i]/snp_count[i])));
				writer.println(markerNames[i] + "\t"
				               + (chrs[i] < 23 ? chrs[i]
				                               : (chrs[i] == 23 ? "X"
				                                                : (chrs[i] == 24 ? "Y"
				                                                                 : (chrs[i] == 25 ? "XY"
				                                                                                  : (chrs[i] == 26 ? "M"
				                                                                                                   : "Un")))))
				               + "\t" + positions[i] + "\t"
				               + (snp_count[i] == 0 ? (snp_sum[i] == 0 ? 0 : "err")
				                                    : (snp_sum[i] / snp_count[i])));
			}
			writer.close();
			log.report("Generated population GC Model " + outputGcModelFullPath);
		} catch (Exception e) {
			log.reportError("Error writing to '" + outputGcModelFullPath + "'");
			log.reportException(e);
		}
	}

	private static String[] getSamplesForTransform(Project proj, boolean excludeExcludeds) {
		if (excludeExcludeds) {
			return Array.subArray(proj.getSamples(), proj.getSamplesToInclude(null));
		} else {
			return proj.getSamples();
		}
	}

	public static void doBatch(Project proj, boolean auto, boolean chrx, boolean sexCent,
	                           boolean transformData, int numChunks, boolean separateQsubFiles,
	                           String pfbFile, String gcmodelFile, String hmmFile,
	                           boolean submitImmed, boolean createCombined, boolean useExcludes,
	                           int threadCount) {
		boolean problem = false;
		Vector<String> execList;
		Logger log = proj.getLog();
		final String runLine = Files.getRunString() + " " + PennCNV.class.getCanonicalName() + " proj="
		                       + new File(proj.getPropertyFilename()).getAbsolutePath();

		String dir = proj.PENNCNV_RESULTS_DIRECTORY.getValue();
		dir += "penn_scripts/";

		if (hmmFile == null || !Files.exists(hmmFile)) {
			hmmFile = Resources.cnv(log).getAllHmm().get();
		}

		if ((pfbFile == null || !Files.exists(pfbFile))
		    && (pfbFile = proj.CUSTOM_PFB_FILENAME.getValue()) == null) {
			System.err.println("Error - could not find " + pfbFile);
			problem = true;
		}
		if ((gcmodelFile == null || !Files.exists(gcmodelFile))
		    && (gcmodelFile = proj.GC_MODEL_FILENAME.getValue()) == null) {
			System.err.println("Error - could not find " + gcmodelFile);
			problem = true;
		}
		if (problem) {
			return;
		}

		if (separateQsubFiles) {
			execList = null;
		} else {
			execList = new Vector<String>();
		}

		String[] samples = getSamplesForTransform(proj, !useExcludes);

		if (auto) {
			if (transformData) {
				log.report("Transforming data for autosomal CNV analysis");
				AnalysisFormats.exportPenncnvSamples(proj, samples, null, null,
				                                     Runtime.getRuntime().availableProcessors());
			}
			log.report("Creating batch scripts for autosomal CNV analysis");
			batch(proj, numChunks, execList, pfbFile, gcmodelFile, hmmFile, "penn_scripts/", "", "");
		}
		if (chrx) {
			MarkerSet ms = proj.getMarkerSet();
			if (ms == null) {
				log.reportError("Error - no marker set available.");
			} else {
				log.report("Transforming data for chromosomal CNV analysis");
				HashSet<String> xMarkers = new HashSet<String>();
				byte[] chrs = ms.getChrs();
				String[] markers = ms.getMarkerNames();
				for (int i = 0; i < chrs.length; i++) {
					if (chrs[i] == 23) {
						xMarkers.add(markers[i]);
					}
				}
				AnalysisFormats.exportPenncnvSamples(proj, samples, xMarkers, "chrX/",
				                                     Runtime.getRuntime().availableProcessors());
			}
			log.report("Creating batch scripts for chromosomal CNV analysis");
			batchX(proj, numChunks, execList, pfbFile, gcmodelFile, hmmFile, "penn_scripts/chrX/",
			       "chrX/", "chrX/");
		}
		if ((auto && chrx) || (auto && createCombined) || (chrx && createCombined)) {
			// write combine script
			String resultsDir = proj.PENNCNV_RESULTS_DIRECTORY.getValue(false, true);
			String outdir = resultsDir + "penn_scripts/";
			new File(outdir).mkdirs();
			String outfile = "combineAutoXCNVs";
			Files.writeArray(new String[] {"cd " + resultsDir,
			                               runLine + " combine=penncnv.cnv,chrX/penncnvX.cnv output=combinedAX.cnv",},
			                 outdir + outfile);
			Files.chmod(outdir + outfile);
		}
		if (sexCent) {
			log.report("Transforming data for 'faked' chromosomal CNV analysis");
			// [males.pfb, females.pfb, sexSpecific.gcModel]

			String[] files = AnalysisFormats.pennCNVSexHackMultiThreaded(proj, gcmodelFile, useExcludes,
			                                                             threadCount);
			// String[] files = AnalysisFormats.pennCNVSexHackSingleThreaded(proj, gcmodelFile);

			log.report("Creating batch scripts for 'faked' chromosomal CNV analysis");
			String scriptDir = "penn_scripts/sexSpecific/";
			batch(proj, numChunks, execList, files[0], files[2], hmmFile, scriptDir + "male/",
			      "sexSpecific/male/", "sexSpecific/male/");
			batch(proj, numChunks, execList, files[1], files[2], hmmFile, scriptDir + "female/",
			      "sexSpecific/female/", "sexSpecific/female/");
			// write combine script
			String resultsDir = proj.PENNCNV_RESULTS_DIRECTORY.getValue(false, true);
			String outdir = resultsDir + "penn_scripts/";
			String outfile = "combineMFCNVs";
			Files.writeArray(new String[] {"cd " + resultsDir,
			                               runLine + " combine=sexSpecific/male/penncnv.cnv output=sexSpecific/male/recodedM.cnv -recode",
			                               runLine + " combine=sexSpecific/female/penncnv.cnv output=sexSpecific/female/recodedF.cnv -recode",
			                               runLine + " combine=sexSpecific/male/recodedM.cnv,sexSpecific/female/recodedF.cnv output=combinedMF.cnv -recode",},
			                 outdir + outfile);
			Files.chmod(outdir + outfile);

			if (auto) {
				outfile = "combineAMFCNVs";
				Files.writeArray(new String[] {"cd " + resultsDir,
				                               runLine + " combine=sexSpecific/male/penncnv.cnv output=sexSpecific/male/recodedM.cnv -recode",
				                               runLine + " combine=sexSpecific/female/penncnv.cnv output=sexSpecific/female/recodedF.cnv -recode",
				                               runLine + " combine=penncnv.cnv,sexSpecific/male/recodedM.cnv,sexSpecific/female/recodedF.cnv output=combinedMF.cnv",},
				                 outdir + outfile);
				Files.chmod(outdir + outfile);
			}

		}

		if (execList != null) {
			Files.qsubExecutor(proj.PROJECT_DIRECTORY.getValue(), execList, null,
			                   proj.PENNCNV_RESULTS_DIRECTORY.getValue() + "runAllPenncnv", 24, 5000, 8);
			log.report("All PennCNV files and scripts have been prepped. The next thing would be to qsub "
			           + proj.PENNCNV_RESULTS_DIRECTORY.getValue() + "runAllPenncnv.pbs");
		}
		List<String> toRun = new ArrayList<String>();
		toRun.add(dir + "assemblePenncnv");
		toRun.add(dir + "chrX/assemblePenncnv");
		if (sexCent) {
			toRun.add(dir + "sexSpecific/female/assemblePenncnv");
			toRun.add(dir + "sexSpecific/male/assemblePenncnv");
			toRun.add(dir + "combineAMFCNVs");
		}
		Files.writeArray(toRun.toArray(new String[toRun.size()]), dir + "parseAllPenncnv");
		Files.chmod(dir + "parseAllPenncnv");

		log.report("Script generation complete. See: " + dir);
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		String rawlog = null;
		String rawcnvs = null;
		int numChunks = 0;
		boolean transformData = true;
		boolean separateQsubs = false;
		boolean auto = true;
		boolean chrx = true;
		boolean sexCent = true;
		Project proj;
		String pfbFile = null;
		String gcmodelFile = null;
		boolean denovoOnly = false;
		boolean parsePFB = false;
		String gc5base = null;
		String logfile = null;
		String[] cnvFiles = null;
		String outputFile = null;
		boolean recode = false;
		boolean submit = false;
		boolean excludes = false;
		String hmmFile = null;
		int numThreads = 1;

		String usage = "\n" + "org.genvisis.cnv.analysis.PennCNV requires 0-1 arguments\n"
		               + "   (0) project properties filename (i.e. proj="
		               + org.genvisis.cnv.Launch.getDefaultDebugProjectFile(false) + " (default))\n"
		               + " AND\n" + "   (1) number of chunks to split everything in to (i.e. chunks="
		               + numChunks + " (default))\n"
		               + "   (2) generate seperate qsub files instead of a single executor chain (i.e. -sepqsub (not the default))\n"
		               + "   (3) generate PennCNV scripts to analyze autosomes (i.e. auto=TRUE (default))\n"
		               + "   (4) generate PennCNV scripts to analyze X Chromosome (i.e. chrx=TRUE (default))\n"
		               + "   (5) recompute centroids of chr23-26 (X, Y, XY, MT) and recode as chr1-4 in subdirectory (i.e. sexSpecificCentroids=TRUE (default))\n"
		               + "   (6) transform sample data into PennCNV data files (i.e. data=TRUE (default))\n"
		               + "   (7) number of threads to use (i.e. threads=" + numThreads + " (default))\n"
		               + "   (8) (optional) use custom pfb file (i.e. pfb=custom.pfb (not the default))\n"
		               + "   (9) (optional) use custom gcmodel file (i.e. gcmodel=custom.gcmodel (not the default))\n"
		               + "   (10) (optional) use an array specific hmm file (i.e. hmm= (no default))\n"
		               +

		               " OR\n"
		               + "   (1) compute file containing project based b allele frequencies for file using parameters in properties file (i.e. -pfb (not the default))\n"
		               + " OR\n"
		               + "   (1) compute a custom gcmodel file for the markers in this project using this file (i.e. gc5base=gc5base.txt (not the default))\n"
		               + " OR\n"
		               + "   (1) parse warnings from log file (i.e. rawlog=final.log (not the default))\n"
		               + " OR\n"
		               + "   (1) raw cnvs to parse (i.e. rawcnv=final.rawcnv (not the default))\n"
		               + "   (2) (optional) parse only de novo variants (i.e. -denovoOnly (not the default))\n"
		               + " OR\n"
		               + "   (1) a comma-separated list of .cnv files to combine together (i.e. combine=/full/path/to/cnv1.cnv,relative/path/to/cnv2.cnv (not the default))\n"
		               + "   (2) full path of the desired output file (i.e. output=/path/to/output/file.cnv (not the default))\n"
		               + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				return;
			} else if (arg.startsWith("proj=")) {
				filename = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("chunks=")) {
				numChunks = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("-sepqsub")) {
				separateQsubs = true;
				numArgs--;
			} else if (arg.startsWith("rawlog=")) {
				rawlog = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("rawcnv=")) {
				rawcnvs = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("-pfb")) {
				parsePFB = true;
				numArgs--;
			} else if (arg.startsWith("gc5base=")) {
				gc5base = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("-denovoOnly")) {
				denovoOnly = true;
				numArgs--;
			} else if (arg.startsWith("pfb=")) {
				pfbFile = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("gcmodel=")) {
				gcmodelFile = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("hmm=")) {
				hmmFile = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("log=")) {
				logfile = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("data=")) {
				transformData = Boolean.parseBoolean(arg.split("=")[1]);
				numArgs--;
			} else if (arg.startsWith("auto=")) {
				auto = Boolean.parseBoolean(arg.split("=")[1]);
				numArgs--;
			} else if (arg.startsWith("chrx=")) {
				chrx = Boolean.parseBoolean(arg.split("=")[1]);
				numArgs--;
			} else if (arg.startsWith("sexSpecificCentroids=")) {
				sexCent = Boolean.parseBoolean(arg.split("=")[1]);
				numArgs--;
			} else if (arg.startsWith("combine=")) {
				cnvFiles = arg.split("=")[1].split(",");
				numArgs--;
			} else if (arg.startsWith("output=")) {
				outputFile = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("-recode")) {
				recode = true;
				numArgs--;
			} else if (arg.startsWith("-submit")) {
				submit = true;
				numArgs--;
			} else if (arg.startsWith("-useExcluded")) {
				excludes = true;
				numArgs--;
			} else if (arg.startsWith("threads=")) {
				numThreads = ext.parseIntArg(arg);
				numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			return;
		}
		try {

			// filename = "C:/workspace/Genvisis/projects/GEDI_exome.properties";
			// parsePFB = true;
			// gc5base = "C:/projects/gcModel/gc5Base.txt";

			// logfile = "penncnv/penncnv.log";
			// rawcnvs = "penncnv/penncnv.rawcnv";

			// filename = "/home/npankrat/projects/GEDI.properties";
			// batch = 60;
			// qsub = true;
			// pfbFile = "gedi.pfb";
			// gcmodelFile = "gedi.gcmodel";
			//
			// batch = 1;
			// filename = "C:/data/FarrarMike/default.properties";
			// qsub = true;
			// pfbFile = "C:/data/FarrarMike/custom.pfb";
			// gcmodelFile = "C:/data/FarrarMike/data/custom.gcmodel";
			// numThreads = 5;

			proj = new Project(filename, logfile, false);
			if (parsePFB) {
				populationBAF(proj);
			}
			if (gc5base != null) {
				gcModel(proj, gc5base, proj.GC_MODEL_FILENAME.getValue(), 100);
			}
			if (numChunks > 0) {
				if (hmmFile == null || !new File(hmmFile).exists()) {
					hmmFile = Resources.cnv(proj.getLog()).getAllHmm().get();
				}
				if (pfbFile == null || !new File(pfbFile).exists()) {
					pfbFile = Resources.cnv(proj.getLog()).genome(proj.GENOME_BUILD_VERSION.getValue())
					                   .getAllPfb().get();
				}
				if (gcmodelFile == null || !new File(pfbFile).exists()) {
					gcmodelFile = Resources.cnv(proj.getLog()).genome(proj.GENOME_BUILD_VERSION.getValue())
					                       .getAllGcmodel().get();
				}
				doBatch(proj, auto, chrx, sexCent, transformData, numChunks, separateQsubs, pfbFile,
				        gcmodelFile, hmmFile, separateQsubs ? submit : false, recode, excludes, numThreads);
			}
			if (rawlog != null) {
				parseWarnings(proj, rawlog);
			}
			if (rawcnvs != null) {
				parseResults(proj, rawcnvs, denovoOnly);
			}

			if (cnvFiles != null && outputFile != null) {
				combineResults(proj, cnvFiles, outputFile, recode);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
