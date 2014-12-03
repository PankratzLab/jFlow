package cnv.analysis;

import java.io.*;
import java.util.*;

import cnv.filesys.MarkerSet;
//import cnv.analysis.FilterCalls;
import cnv.filesys.Project;
import cnv.filesys.Sample;
import cnv.qc.SexChecks;
import cnv.var.CNVariant;
import cnv.var.SampleData;
import common.*;

public class PennCNV {
	public static final String[] QC_HEADS = {"LRR_mean", "LRR_median", "LRR_SD", "BAF_mean", "BAF_median", "BAF_SD", "BAF_DRIFT", "WF", "GCWF"};
	public static final String[] ERRORS = {"large SD for LRR", "drifting BAF values", "waviness factor values", "Small-sized CNV calls", "NoCall rate"};
	public static final String QC_SUMMARY_FILE = "Sample_QC.xln";

	public static void batch(Project proj, int numBatches, boolean qsub, String pfbFile, String gcmodelFile, String scriptSubDir, String dataSubDir, String resultsSubDir) {
		String init, commands;
		PrintWriter writer;
		String[] files;
		int step;
		String execDir, dataDir, resultsDir, projDir, pennDir;
		Logger log;
		
		log = proj.getLog();
		projDir = proj.getProjectDir();
		execDir = proj.getDir(Project.PENNCNV_EXECUTABLE_DIRECTORY);
		pennDir = proj.getDir(Project.PENNCNV_RESULTS_DIRECTORY);
		dataDir = pennDir + proj.getProperty(Project.PENNCNV_DATA_DIRECTORY) + dataSubDir;
		resultsDir = pennDir + resultsSubDir;
		
		if (pfbFile != null) {
			pfbFile = ext.replaceTilde(pfbFile);
			if (!pfbFile.startsWith("/") && (pfbFile.charAt(1) != ':')) {
				pfbFile = ext.pwd() + pfbFile;
			}
			if (!Files.exists(pfbFile)) {
				log.reportError("Error - pfb file '"+pfbFile+"' does not exist; aborting");
				return;
			}
		}

		if (gcmodelFile != null) {
			gcmodelFile = ext.replaceTilde(gcmodelFile);
			if (!gcmodelFile.startsWith("/") && (gcmodelFile.charAt(1) != ':')) {
				gcmodelFile = ext.pwd() + gcmodelFile;
			}
			if (!Files.exists(gcmodelFile)) {
				log.reportError("Error - gcmodel file '"+gcmodelFile+"' does not exist; aborting");
				return;
			}
		}
		
		new File(resultsDir).mkdirs();
		
		files = new File(dataDir).list(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				return file.length()>1000;
			}
		});
		log.report("Found "+files.length+" files");

		step = (int)Math.ceil((double)files.length/(double)numBatches);
		log.report("Which means the step for "+numBatches+" batches would be "+step);
		
		for (int i = 0; i<numBatches; i++) {
			try {
				writer = new PrintWriter(new FileWriter(resultsDir+"list"+(i+1)+".txt"));
				for (int j = i*step; j<Math.min(files.length, (i+1)*step); j++) {
					if (files[j].endsWith(".gz")) {
						writer.println("`gunzip -c "+files[j]+"`");
					} else {
						writer.println(files[j]);
					}
				}
				writer.close();
			} catch (Exception e) {
				log.reportError("Error writing to " + resultsSubDir + "list"+(i+1)+".txt");
				log.reportException(e);
			}
		}

		init = "cd "+dataDir;
		commands = execDir+"detect_cnv.pl -test -conf -hmm "+execDir+"lib/hhall.hmm -pfb "+(pfbFile==null?execDir+"lib/hhall.hg18.pfb":pfbFile)+" -gcmodel "+(gcmodelFile==null?execDir+"lib/hhall.hg18.gcmodel":gcmodelFile)+" -list "+resultsDir+"list[%0].txt -log "+resultsDir+"[%0].log -out "+resultsDir+"[%0].rawcnv > "+resultsDir+"[%0].out";

		if (qsub) {
			Files.qsub(pennDir + scriptSubDir + "runPenn", dataDir, numBatches, commands, Matrix.toMatrix(Array.stringArraySequence(numBatches, "")), 2200, 16);
		} else {
			Files.batchIt(pennDir + scriptSubDir + "penn", init, numBatches, commands, Array.stringArraySequence(numBatches, ""));
		}
		Files.writeList(new String[] {
				"cd " + projDir,
				"cat " + resultsDir + "*.log > " + resultsDir + "penncnv.rawlog",
				"cat " + resultsDir + "*.rawcnv > " + resultsDir + "penncnv.rawcnv",
				"java -cp ~/park.jar cnv.analysis.PennCNV proj="+proj.getPropertyFilename()+" rawlog=" + resultsDir + "penncnv.rawlog",
				"java -cp ~/park.jar cnv.analysis.PennCNV proj="+proj.getPropertyFilename()+" rawcnv=" + resultsDir + "penncnv.rawcnv",
		}, pennDir + scriptSubDir + "assemblePenncnv");
		Files.chmod(pennDir + scriptSubDir + "assemblePenncnv");
	}

	public static void batchX(Project proj, int numBatches, boolean qsub, String pfbFile, String gcmodelFile, String scriptSubDir, String dataSubDir, String resultsSubDir) {
		String init, commands;
		PrintWriter writer;
		String[] files;
		int step;
		String execDir, dataDir, resultsDir, projDir, pennDir;
		Logger log;
		
		log = proj.getLog();
		projDir = proj.getProjectDir();
		execDir = proj.getDir(Project.PENNCNV_EXECUTABLE_DIRECTORY);
		pennDir = proj.getDir(Project.PENNCNV_RESULTS_DIRECTORY);
		dataDir = pennDir + proj.getProperty(Project.PENNCNV_DATA_DIRECTORY) + dataSubDir;
		resultsDir = pennDir + resultsSubDir;
		
		if (!Files.exists(proj.getFilename("SAMPLE_DATA_FILENAME", false, false), proj.getJarStatus())) {
			log.reportError("Error - sample data file " + proj.getProperty("SAMPLE_DATA_FILENAME") + " does not exist;");
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
				log.reportError("Error - pfb file '"+pfbFile+"' does not exist; aborting");
				return;
			}
		}

		if (gcmodelFile != null) {
			gcmodelFile = ext.replaceTilde(gcmodelFile);
			if (!gcmodelFile.startsWith("/") && (gcmodelFile.charAt(1) != ':')) {
				gcmodelFile = ext.pwd() + gcmodelFile;
			}
			if (!Files.exists(gcmodelFile)) {
				log.reportError("Error - gcmodel file '"+gcmodelFile+"' does not exist; aborting");
				return;
			}
		}
		
		new File(resultsDir).mkdirs();
		
		String newGCFile = dataDir + "chrX.gcModel";
		String newPFBFile = dataDir + "chrX.pfb";
		
		pfbFile = filterPFB(proj, pfbFile, newPFBFile, new String[]{"23", "X"});
		gcmodelFile = filterGCModel(proj, gcmodelFile, newGCFile, new String[]{"23", "X"});
		
		String sexFileStatus = writeSexFile(proj, sampleData, dataDir, log);
		if (sexFileStatus.length() > 0) {
			log.reportError("Error - " + sexFileStatus);
			return;
		}
	
		files = new File(dataDir).list(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				return file.length()>1000;
			}
		});
		log.report("Found "+files.length+" files");
	
		step = (int)Math.ceil((double)files.length/(double)numBatches);
		log.report("Which means the step for "+numBatches+" batches would be "+step);
		
		for (int i = 0; i<numBatches; i++) {
			try {
				writer = new PrintWriter(new FileWriter(resultsDir+"list"+(i+1)+".txt"));
				for (int j = i*step; j<Math.min(files.length, (i+1)*step); j++) {
					if (files[j].endsWith(".gz")) {
						writer.println("`gunzip -c "+files[j]+"`");
					} else {
						writer.println(files[j]);
					}
				}
				writer.close();
			} catch (Exception e) {
				log.reportError("Error writing to " + resultsSubDir + "list"+(i+1)+".txt");
				log.reportException(e);
			}
		}
	
		init = "cd " + dataDir;
		commands = execDir + "detect_cnv.pl -test -conf -hmm " + execDir + "lib/hhall.hmm -pfb " + (pfbFile == null ? execDir + "lib/hhall.hg18.pfb" : pfbFile) + " -gcmodel " + (gcmodelFile == null ? execDir + "lib/hhall.hg18.gcmodel" : gcmodelFile) + " -chrx -sexfile " + dataDir + "sex_file.txt -list " + resultsDir + "list[%0].txt -log " + resultsDir + "[%0].log -out " + resultsDir + "[%0].rawcnv > " + resultsDir + "[%0].out";
	
		if (qsub) {
			Files.qsub(pennDir + scriptSubDir + "runPennX", dataDir, numBatches, commands, Matrix.toMatrix(Array.stringArraySequence(numBatches, "")), 2200, 16);
		} else {
			Files.batchIt(pennDir + scriptSubDir + "pennX", init, numBatches, commands, Array.stringArraySequence(numBatches, ""));
		}
		Files.writeList(new String[] {
				"cd " + projDir,
				"cat " + resultsDir + "*.log > " + resultsDir + "penncnvX.rawlog",
				"cat " + resultsDir + "*.rawcnv > " + resultsDir + "penncnvX.rawcnv",
				// don't parse warnings; the parseWarnings method isn't written to parse X-chromosome warnings
				"java -cp ~/park.jar cnv.analysis.PennCNV proj="+proj.getPropertyFilename()+" rawcnv=" + resultsDir + "penncnvX.rawcnv",
		}, pennDir + scriptSubDir + "assemblePenncnvX");
		Files.chmod(pennDir + scriptSubDir + "assemblePenncnvX");
	}

	private static String filterGCModel(Project proj, String gcModelFile, String newGCFile, String[] chrs) {
		// TODO combine method with filterPFB  -  literally the same except different names/extensions
		BufferedReader reader = null;
		PrintWriter writer = null;
		
		try {
			reader = new BufferedReader(new FileReader(gcModelFile));
			writer = new PrintWriter(new FileWriter(newGCFile));

			String header;
			String temp;
			String[] line;
			if (reader.ready()) {
				header = reader.readLine();
				writer.println(header);
			}
			while(reader.ready()) {
				temp = reader.readLine();
				line = temp.trim().split("[\\s]+");
				for (String chr : chrs) {
					if (line[1].equals(chr)) {
						writer.println(temp);
					}
				}
			}
		} catch (IOException e) {
			proj.getLog().reportError("Error - filtering gcModel failed");
			proj.getLog().reportException(e);
			return gcModelFile;
		} finally {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					proj.getLog().reportError("Error - couldn't properly close file reader for " + gcModelFile);
					proj.getLog().reportException(e);
				}
				reader = null;
			}
			if (writer != null) {
				writer.close();
				writer = null;
			}
		}
		
		return newGCFile;
	}

	private static String filterPFB(Project proj, String pfbFile, String newPFBFile, String[] chrs) {
		// TODO combine method with filterGCModel  -  literally the same except different names/extensions
		BufferedReader reader = null;
		PrintWriter writer = null;
		
		try {
			reader = new BufferedReader(new FileReader(pfbFile));
			writer = new PrintWriter(new FileWriter(newPFBFile));

			String header;
			String temp;
			String[] line;
			if (reader.ready()) {
				header = reader.readLine();
				writer.println(header);
			}
			while(reader.ready()) {
				temp = reader.readLine();
				line = temp.trim().split("[\\s]+");
				for (String chr : chrs) {
					if (line[1].equals(chr)) {
						writer.println(temp);
					}
				}
			}
		} catch (IOException e) {
			proj.getLog().reportError("Error - filtering PFB failed");
			proj.getLog().reportException(e);
			return pfbFile;
		} finally {
			if (reader != null) {
				try {
					reader.close();
				} catch (IOException e) {
					proj.getLog().reportError("Error - couldn't properly close file reader for " + pfbFile);
					proj.getLog().reportException(e);
				}
				reader = null;
			}
			if (writer != null) {
				writer.close();
				writer = null;
			}
		}
		
		return newPFBFile;
	}

	private static String writeSexFile(Project proj, SampleData sampleData, String resultsDir, Logger log) {
		String sampleDataFile = proj.getFilename(Project.SAMPLE_DATA_FILENAME, false, false);
		String[] header = Files.getHeaderOfFile(sampleDataFile, proj.getLog());
		int sexInd = -1;
		for (int i = 0; i < header.length; i++) {
			if (("CLASS=" + SexChecks.EST_SEX_HEADER).toUpperCase().equals(header[i].toUpperCase())) {
				sexInd = i;
				break;
			}
		}
		if (sexInd == -1) {
			return "no estimated sex found in sample data file - please run SexChecks with -check argument to generate the required data";
		}
		Hashtable<String, Vector<String>> sexData = HashVec.loadFileToHashVec(sampleDataFile, 0, new int[] { sexInd }, "\t", true, false);
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(resultsDir + "sex_file.txt"));
			for (Map.Entry<String, Vector<String>> lineData : sexData.entrySet()) {
				int estSex = Integer.parseInt(lineData.getValue().get(0));
				estSex = SexChecks.KARYOTYPES[estSex].contains("XX") ? 2 : 1;
				writer.println(lineData.getKey() + "\t" + estSex);
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
		Hashtable<String,String[]> hash = new Hashtable<String,String[]>();
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
		lrrSD_cutoff = proj.getDouble(Project.LRRSD_CUTOFF);
				
		try {
			reader = new BufferedReader(new FileReader(filename));
			while (reader.ready()) {
				temp = reader.readLine();
				temp = translateDerivedSamples(temp, log);
				line = temp.trim().split("[\\s]+");
				try {
					if (temp.contains("quality summary")) {
						sampleID = line[4].substring(line[4].lastIndexOf("/")+1, line[4].indexOf(":"));
						v.add(sampleID);
						data = new String[QC_HEADS.length+ERRORS.length];
						for (int i = 0; i<ERRORS.length; i++) {
							data[i] = ".";
						}
						if (line.length<QC_HEADS.length+5) {
							log.reportError("Error - line doesn't have all the expected pieces:");
							log.reportError(temp);
						}
						for (int i = 0; i<QC_HEADS.length; i++) {
							data[ERRORS.length+i] = line[5+i].split("=")[1];
						}
						hash.put(sampleID, data);
					} else if (temp.startsWith("WARNING")) {
						if (temp.contains("Small-sized CNV calls")) {
							// use old trav
						} else {
							sampleID = line[3].substring(line[3].lastIndexOf("/")+1);
						}
						data = hash.get(sampleID);
						err = -1;
						for (int i = 0; i<ERRORS.length; i++) {
							if (temp.contains(ERRORS[i])) {
								err = i;
							}
						}
						if (err==-1) {
							log.reportError("Unknown WARNING: "+temp);
						} else {
							data[err] = "1";
						}
					}
				} catch (Exception e) {
					log.reportError("Error with: "+temp);
					log.reportException(e);
				}
			}
			reader.close();

			writer = new PrintWriter(new FileWriter(proj.getProjectDir()+QC_SUMMARY_FILE));
			writer.print("Sample\tFID\tIID\tUse_"+ext.formDeci(lrrSD_cutoff, 2));
			for (int i = 0; i<ERRORS.length; i++) {
				writer.print("\t"+ERRORS[i]);
			}
			for (int i = 0; i<QC_HEADS.length; i++) {
				writer.print("\t"+QC_HEADS[i]);
			}
			writer.println();
			int[] keys = Sort.quicksort(Array.toStringArray(v));
			for (int i = 0; i<v.size(); i++) {
				sampleID = v.elementAt(keys[i]);
				data = hash.get(sampleID);
				ids = sampleData.lookup(sampleID);
				writer.print(sampleID+"\t"+(ids==null?"NotInSampleData\t"+sampleID:ids[1]));
				writer.print("\t"+(data[1].equals("1")||data[2].equals("1")||Double.parseDouble(data[6])>lrrSD_cutoff?"0":"1"));
				writer.println("\t"+Array.toStr(data));
			}
			writer.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \""+filename+"\" not found in current directory");
			return;
		} catch (IOException ioe) {
			log.reportError("Error reading file \""+filename+"\"");
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
		
		trav = str.substring(start+1, stop); 
		if (trav.contains("`")) {
			log.reportError("Error - more than one set of quotes for: "+str);
		}
		if (trav.startsWith("gunzip -c") && trav.endsWith(".gz")) {
			trav = trav.substring(9, trav.length()-3).trim();
		} else {
			log.reportError("Error - not currently set up to handle the following construction into a sample_ID: "+trav);
		}
		
		return str.substring(0, start)+trav+str.substring(stop+1);		
	}

	public static void combineResults(Project proj, String[] cnvFiles, String outputFile, boolean recode) {
		BufferedReader reader;
		PrintWriter writer;
		Logger log = proj.getLog();
		
		// TODO check input and output file names for .cnv extension( - error if not? or just warning...?)
		
		java.util.HashMap<String, java.util.TreeMap<String, java.util.ArrayList<String[]>>> cnvSet = new HashMap<String, TreeMap<String, ArrayList<String[]>>>();
		
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
				while(reader.ready()) {
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
			log.reportError("Error: file \""+currFile+"\" not found in current directory");
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
	
	/*
	java -cp ~/park.jar cnv.analysis.PennCNV proj=/home/pankrat2/coleb/projects/NY_Registry_3defects.properties combine=/home/pankrat2/shared/ny_registry/3defects/penncnv/sexSpecific/male/penncnv.cnv,/home/pankrat2/shared/ny_registry/3defects/penncnv/sexSpecific/female/penncnv.cnv output=/home/pankrat2/shared/ny_registry/3defects/penncnv/sexSpecific/combinedMF.cnv -recode
	*/
	
	public static void parseResults(Project proj, String filename, boolean denovoOnly) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String temp, trav;
		Vector<String> warnings;
		Hashtable<String,Vector<String>> pedinfo;
		int[] position;
		String score;
		SampleData sampleData;
		String famIndPair;
		Hashtable<String,String> hash;
		String[] ids, fams, inds;
		long time;
		int sex;
		Logger log;
		
		log = proj.getLog();
		log.report("Parsing PennCNV rawcnvs...");
		time = new Date().getTime();
		
		if (!Files.exists(filename)) {
			log.reportError("Error - could not find file '"+filename+"'");
			return;
		}

		warnings = new Vector<String>();
		sampleData = proj.getSampleData(2, false);
		pedinfo = new Hashtable<String, Vector<String>>();
		try {
			reader = new BufferedReader(new FileReader(filename));
			writer = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+".cnv"));
			writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER));
			hash = new Hashtable<String,String>();
			while (reader.ready()) {
				temp = reader.readLine();
				if (!temp.startsWith("NOTICE:")) {
					temp = translateDerivedSamples(temp, log);
					line = temp.trim().split("[\\s]+");
					position = Positions.parseUCSClocation(line[0]);
					trav = line[4];
					trav = trav.substring(trav.lastIndexOf("/")+1);
					ids = sampleData.lookup(trav);
					if (ids == null) {
						if (!hash.containsKey(trav)) {
							log.reportError("Error - '"+trav+"' was not found in "+proj.getFilename(Project.SAMPLE_DATA_FILENAME));
							hash.put(trav, "");
						}
						famIndPair = trav+"\t"+trav;
					} else {
						famIndPair = ids[1];
					}
					
					ids = famIndPair.split("\t");
					HashVec.addToHashVec(pedinfo, ids[0], ids[1], true);
	
					if (line.length < 8 || !line[7].startsWith("conf=") || line[7].toUpperCase().contains("NAN")) {
						score = "-1";
						if (!warnings.contains(trav) && warnings.size() < 10) {
							log.reportError("Warning - no conf estimates for "+trav);
							warnings.add(trav);
						}
					} else {
						score = ext.formDeci(Double.parseDouble(line[7].substring(5)), 4, true);
					}
					if (!denovoOnly || ext.indexFactors(new String[][] {{"statepath=33"}}, line, false, false, false, false)[0] > 0) {
						writer.println(famIndPair+"\t"+position[0]+"\t"+position[1]+"\t"+position[2]+"\t"+line[3].substring(line[3].indexOf("=")+1)+"\t"+score+"\t"+line[1].substring(7));
					}
				}
			}
			reader.close();
			writer.close();

//			FilterCalls.stdFilters(dir, ext.rootOf(filename)+".cnv", MAKE_UCSC_TRACKS);

			writer = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+".fam"));
			fams = HashVec.getKeys(pedinfo, true, true);
			for (int i = 0; i<fams.length; i++) {
				inds = Sort.putInOrder(Array.toStringArray(pedinfo.get(fams[i])), true);
				for (int j = 0; j < inds.length; j++) {
					ids = sampleData.lookup(fams[i]+"\t"+inds[j]);
					if (ids != null) {
						sex = sampleData.getSexForIndividual(ids[0]);
					} else {
						sex = 0;
					}
					writer.println(fams[i]+"\t"+inds[j]+"\t0\t0\t"+Math.max(0, sex)+"\t-9");
				}
			}
			writer.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \""+ext.rootOf(filename, false)+".cnv\" not found in current directory");
			return;
		} catch (IOException ioe) {
			log.reportError("Error reading file \""+ext.rootOf(filename, false)+"\"");
			return;
		}

		log.report("...finished in " + ext.getTimeElapsed(time));
	}

	// Available in cnv.Launch
//	public static void fromParameters(String filename, Logger log) {
//		Vector<String> params;
//
//		params = Files.parseControlFile(filename, "penncnv", new String[] {"proj=/home/npankrat/projects/GEDI.properties", "rawcnv=all.rawcnv", "rawlog=all.log"}, log);
//
//		if (params != null) {
//			params.add("log=" + log.getFilename());
//			main(Array.toStringArray(params));
//		}
//	}
	
	/**
	 * Calculate the population BAF (B Allele Frequency based on all the samples available in the) data. Output is going to be saved on disk. In PennCnv, this file is also called snpFile.
	 * @param proj The project you are going to run PennCNV on.
	 * 
	 * The output file looks like the the following:
	 * 	Name    	Chr     Position    PFB
	 * 	rs1000113   5       150220269   0.564615751221256
	 * 	rs1000115   9       112834321   0.565931333264192
	 * 	rs10001190  4       6335534		0.5668604380025
	 * 	rs10002186  4       38517993    0.57141752993563
	 * 	rs10002743  4       6327482		0.567557695424774
	 * 
	 */
	public static void populationBAF(Project proj) {
		PrintWriter writer;
		Sample samp;
		String[] sampleList;
		String[] markerNames;
		double[] bafSum;
		int[] bafCounts, genoCounts;
		float[] bafs;
		double[] bafAverage;
//		Hashtable<String, String> samples;
		MarkerSet markerSet;
		byte[] chrs, genotypes;
		int[] positions;
		String filename, output;
		Logger log;

		log = proj.getLog();
		filename = proj.getFilename(Project.SAMPLE_SUBSET_FILENAME, true, false);

		if (ext.rootOf(filename) == null || ext.rootOf(filename).equals("")) {
			sampleList = proj.getSampleList().getSamples();
			output = proj.getProjectDir()+"custom.pfb";
		} else if (Files.exists(filename, proj.getJarStatus())) {
			log.report("filename: "+filename);
			sampleList = HashVec.loadFileToStringArray(filename, false, new int[] {0}, false);
			output = proj.getProjectDir()+ext.rootOf(filename)+".pfb";
		} else {
			proj.message("Failed to load \""+filename+"\"");
			return;
		}

		markerSet = proj.getMarkerSet();
		markerNames = markerSet.getMarkerNames();
		chrs = markerSet.getChrs();
		positions = markerSet.getPositions();
		bafSum = new double[chrs.length];
		bafCounts = new int[chrs.length];
		genoCounts = new int[chrs.length];
		for (int i=0; i<sampleList.length; i++) {
			if (i % 100 == 0) {
				log.report("Loading file "+(i+1)+" of "+sampleList.length);
			}
			samp = proj.getFullSampleFromRandomAccessFile(sampleList[i]);
			bafs = samp.getBAFs();
			genotypes = samp.getAB_Genotypes();
			for (int j=0; j<bafSum.length; j++) {
				if (!Float.isNaN(bafs[j])) {
					bafSum[j] += bafs[j];
					bafCounts[j]++;
					if (genotypes[j] >= 0) {
						genoCounts[j]++;
					}
				}
			}
		}
		bafAverage = new double[chrs.length];
		for (int i=0; i<bafSum.length; i++) {
			if (genoCounts[i]!=0) {
				bafAverage[i] = bafSum[i] / bafCounts[i];
			} else {
				bafAverage[i] = 2;
			}
		}

		try {
			writer = new PrintWriter(new FileWriter(output));
			writer.println("Name\tChr\tPosition\tPFB");
			for (int i = 0; i<markerNames.length; i++) {
//				writer.println(markerNames[i]+"\t"+chrs[i]+"\t"+positions[i]+"\t"+bafAverage[i]);
				writer.println(markerNames[i] + "\t" + (chrs[i]<23?chrs[i]:(chrs[i]==23?"X":(chrs[i]==24?"Y":(chrs[i]==25?"XY":(chrs[i]==26?"M":"Un"))))) + "\t" + positions[i] + "\t" + bafAverage[i]);
			}
			writer.close();
			log.report("Population BAF file is now ready at: " + output);
		} catch (Exception e) {
			log.reportError("Error writing to '" + output + "'");
			log.reportException(e);
		}
	}


	/**
	 * Generate the GCModel file needed by PennCNV software (http://www.openbioinformatics.org/penncnv/).
	 * @param proj The project you are going to run PennCNV on.
	 * @param inputGcBaseFullPath The user-supplied genome builds. Positions within each chromosome must be sorted by increasing order. For example, http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/gc5Base.txt.gz
	 * @param outputGcModelFullPath The name of the GCModel file.
	 * @param numwindow For each SNP, GC Content is calculated for the range of numwindow*5120 before its location through numwindow*5120 after. To use the default setting of 100, please enter 0.
	 * 
	 * In order to be able to cross-reference with the same feature in PennCNV, this code intends to base on cal_gc_snp.pl in PennCNV package. But since the difference between Java and Perl, significant structural changes have been made.
	 * 
	 * A sample gcBase file looks like below (There is no header line):
	 *	585     chr1    0       5120    chr1.0  5       1024    0       /gbdb/hg18/wib/gc5Base.wib      0       100     1024    59840   3942400
	 *	585     chr1    5120    10240   chr1.1  5       1024    1024    /gbdb/hg18/wib/gc5Base.wib      0       100     1024    59900   3904400
	 *	585     chr1    10240   15360   chr1.2  5       1024    2048    /gbdb/hg18/wib/gc5Base.wib      0       100     1024    55120   3411200
	 *	585     chr1    15360   20480   chr1.3  5       1024    3072    /gbdb/hg18/wib/gc5Base.wib      0       100     1024    49900   3078800
	 *	585     chr1    20480   25600   chr1.4  5       1024    4096    /gbdb/hg18/wib/gc5Base.wib      0       100     1024    47600   2682400
	 *
	 * A sample gcModel file (the output) look like this:
	 *	Name    	Chr     Position    GC
	 *	rs4961  	4       2876505 	48.6531211131841
	 *	rs3025091   11      102219838   38.1080923507463
	 *	rs3092963   3       46371942    44.4687694340796
	 *	rs3825776   15      56534122    40.7894123134328
	 *	rs17548390  6       3030512 	45.3604050062189
	 *	rs2276302   11      113355350   43.8200598569652
	 *     
	 * snpFile or pfb file (Population B Allele Frequency) is an additional data file that is required by the corresponding feature in PennCNV but not here. It looks like the this:
	 * 	Name    	Chr     Position    PFB
	 * 	rs1000113   5       150220269   0.564615751221256
	 * 	rs1000115   9       112834321   0.565931333264192
	 * 	rs10001190  4       6335534		0.5668604380025
	 * 	rs10002186  4       38517993    0.57141752993563
	 * 	rs10002743  4       6327482		0.567557695424774
	 * 
	 */
	public static void gcModel(Project proj, String inputGcBaseFullPath, String outputGcModelFullPath, int numwindow) {
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
		chrs = markerSet.getChrs();	//to be used only in the SnpFile (.pdf file) block. But then in the outputFile.
		positions = markerSet.getPositions();	//to be used only in the SnpFile (.pdf file) block. But then in the GcFile block.

		// How do we know whether "numwindow==null" ???
		if (numwindow==0) {
			numwindow = 100;
		}

		snp_count = new int[markerNames.length];
		snp_sum = new double[markerNames.length];
		chr_index=0;

		// load gcFile
		try {
			reader = new BufferedReader(new FileReader(inputGcBaseFullPath));
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				curchr = Positions.chromosomeNumber(line[1]);
				curstart = Integer.parseInt(line[2]);
				curend = Integer.parseInt(line[3]);
				curcount = Integer.parseInt(line[11]);
				cursum = Double.parseDouble(line[12]);
				
				// For each record in gcFile, scan snpFile for the range of +/- numwindow*5120  
				if (curchr>0) {
					if (curchr == prechr) {
						if (curstart < prestart) {
							log.reportError("Error in gcFile: a record in chr"+curchr+" has position "+curstart+", less then the previous position $prestart");
							reader.close();
							return;
						}
					} else if (seen_chr.containsKey(curchr)) {
						log.reportError("Error in gcFile: rows of the same chromosome must be adjacent. But now chr"+curchr+" occur multiple times in non-continuous segment of the "+inputGcBaseFullPath+": at "+curchr+":"+curstart);
						reader.close();
						return;
					} else {
						seen_chr.put(curchr,(byte)1);
	
						skip = false;
						if (chrs[chr_index]>curchr) {
							chr_index=0;
						}
						while (chr_index<(chrs.length) && chrs[chr_index]!=curchr) {
							chr_index ++;
						}
					}
	
					if (!(curchr == prechr && skip)) {
						while (chr_index<chrs.length && chrs[chr_index]==curchr && positions[chr_index]<Math.max(curstart - numwindow*5120, 0)) {
							chr_index ++;
						}
						if (chr_index>=chrs.length) {
							chr_index=0;
							skip=true;
						} else if (chrs[chr_index]!=curchr) {
							skip=true;
						} else {
							for (int i = chr_index; i<snp_count.length && chrs[i]==curchr && positions[i]<=(curend + numwindow*5120); i++) {
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
		
		// load pfb file or generate it
		
		// output the result
		try {
			writer = new PrintWriter(new FileWriter(outputGcModelFullPath));
			writer.println("Name\tChr\tPosition\tGC");
			for (int i = 0; i<markerNames.length; i++) {
//				writer.println(markerNames[i]+"\t"+chrs[i]+"\t"+positions[i]+"\t"+(snp_count[i]==0?(snp_sum[i]==0?0:"err"):(snp_sum[i]/snp_count[i])));
				writer.println(markerNames[i] + "\t" + (chrs[i]<23?chrs[i]:(chrs[i]==23?"X":(chrs[i]==24?"Y":(chrs[i]==25?"XY":(chrs[i]==26?"M":"Un"))))) + "\t"+positions[i] + "\t" + (snp_count[i]==0?(snp_sum[i]==0?0:"err"):(snp_sum[i]/snp_count[i])));
			}
			writer.close();
			log.report("Generated population GC Model "+outputGcModelFullPath);
		} catch (Exception e) {
			log.reportError("Error writing to '" + outputGcModelFullPath + "'");
			log.reportException(e);
		}
	}

	public static void doBatch(Project proj, boolean auto, boolean chrx, boolean sexCent, boolean transformData, int batch, boolean qsub, String pfbFile, String gcmodelFile, boolean submitImmed, boolean createCombined) {
		if (transformData) {
			
			String[] samples = Array.subArray(proj.getSamples(), proj.getSamplesToInclude(null));

			if (auto) {
				proj.getLog().report("Transforming data for autosomal CNV analysis");
				AnalysisFormats.penncnv(proj, samples, null, null);
			}
			if (chrx) {
				MarkerSet ms = proj.getMarkerSet();
				if (ms == null) {
					proj.getLog().reportError("Error - no marker set available.");
				} else {
					proj.getLog().report("Transforming data for chromosomal CNV analysis");
					HashSet<String> xMarkers = new HashSet<String>();
					byte[] chrs = ms.getChrs();
					String[] markers = ms.getMarkerNames();
					for (int i = 0; i < chrs.length; i++) {
						if (chrs[i] == 23) {
							xMarkers.add(markers[i]);
						}
					}
					AnalysisFormats.penncnv(proj, samples, xMarkers, "chrX/");
				}
			}
		}
		if (auto) {
			proj.getLog().report("Creating batch scripts for autosomal CNV analysis");
			batch(proj, batch, qsub, pfbFile, gcmodelFile, "penn_scripts/", "", "");
		}
		if (chrx) {
			proj.getLog().report("Creating batch scripts for chromosomal CNV analysis");
			batchX(proj, batch, qsub, pfbFile, gcmodelFile, "penn_scripts/chrX/", "chrX/", "chrX/");
		}
		if ((auto && chrx) || (auto && createCombined) || (chrx && createCombined) ) {
			// write combine script
			String resultsDir = proj.getDir(Project.PENNCNV_RESULTS_DIRECTORY);
			String outdir = resultsDir + "penn_scripts/";
			String outfile = "combineAutoXCNVs";
			Files.writeList(new String[] {
					"cd " + resultsDir,
					"java -cp ~/park.jar cnv.analysis.PennCNV proj=" + proj.getPropertyFilename() + " combine=penncnv.cnv,chrX/penncnvX.cnv output=combinedAX.cnv",
			}, outdir + outfile);
			Files.chmod(outdir + outfile);
		}
		if (sexCent) {
			proj.getLog().report("Transforming data for 'faked' chromosomal CNV analysis");
			// [males.pfb, females.pfb, sexSpecific.gcModel]
			
			String[] files = AnalysisFormats.pennCNVSexHackMultiThreaded(proj, gcmodelFile);
//			String[] files = AnalysisFormats.pennCNVSexHackSingleThreaded(proj, gcmodelFile);

			proj.getLog().report("Creating batch scripts for 'faked' chromosomal CNV analysis");
			String scriptDir = "penn_scripts/sexSpecific/";
			batch(proj, batch, qsub, files[0], files[2], scriptDir + "male/", "sexSpecific/male/", "sexSpecific/male/");
			batch(proj, batch, qsub, files[1], files[2], scriptDir + "female/", "sexSpecific/female/", "sexSpecific/female/");
			// write combine script
			String resultsDir = proj.getDir(Project.PENNCNV_RESULTS_DIRECTORY);
			String outdir = resultsDir + "penn_scripts/";
			String outfile = "combineMFCNVs";
			Files.writeList(new String[] {
					"cd " + resultsDir,
					"java -cp ~/park.jar cnv.analysis.PennCNV proj=" + proj.getPropertyFilename() + " combine=sexSpecific/male/penncnv.cnv output=sexSpecific/male/recodedM.cnv -recode",
					"java -cp ~/park.jar cnv.analysis.PennCNV proj=" + proj.getPropertyFilename() + " combine=sexSpecific/female/penncnv.cnv output=sexSpecific/female/recodedF.cnv -recode",
					"java -cp ~/park.jar cnv.analysis.PennCNV proj=" + proj.getPropertyFilename() + " combine=sexSpecific/male/recodedM.cnv,sexSpecific/female/recodedF.cnv output=combinedMF.cnv -recode",
			}, outdir + outfile);
			Files.chmod(outdir + outfile);
			
		}
		
		if (qsub && submitImmed) {
			
			// run ./master.runPenn scripts:
			//  - auto
			//  - chrx
			//  - sex/fem
			//  - sex/mal
			
		}
		
		
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		String rawlog = null;
		String rawcnvs = null;
		int batch = 0;
		boolean transformData = true;
		boolean qsub = false;
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
		
		String usage = "\n"+
		"cnv.park.PennCNV requires 0-1 arguments\n"+
		"   (0) project properties filename (i.e. proj="+cnv.Launch.getDefaultDebugProjectFile(false)+" (default))\n"+
		" AND\n" + 
		"   (1) number of batches to do (i.e. batch=12 (not the default))\n"+
		"   (2) generate qsub files instead of batch files (i.e. -qsub (not the default))\n"+
		//"   (3) if -qsub is specified, -submit can be used to automatically start the scripts generated by PennCNV (i.e. -submit (not the default))\n" + 
		"   (3) generate PennCNV scripts to analyze autosomes (i.e. auto=TRUE (default))\n" + 
		"   (4) generate PennCNV scripts to analyze X Chromosome (i.e. chrx=TRUE (default))\n" + 
		"   (5) recompute centroids of chr23-26 (X, Y, XY, MT) and recode as chr1-4 in subdirectory (i.e. sexSpecificCentroids=TRUE (default))\n" + 
		"   (6) transform sample data into PennCNV data files (i.e. data=TRUE (default))\n"+
		"   (7) (optional) use custom pfb file (i.e. pfb=custom.pfb (not the default))\n"+
		"   (8) (optional) use custom gcmodel file (i.e. gcmodel=custom.gcmodel (not the default))\n"+
		" OR\n"+
		"   (1) compute file containing project based b allele frequencies for file using parameters in properties file (i.e. -pfb (not the default))\n"+
		" OR\n"+
		"   (1) compute a custom gcmodel file for the markers in this project using this file (i.e. gc5base=gc5base.txt (not the default))\n"+
		" OR\n"+
		"   (1) parse warnings from log file (i.e. rawlog=final.log (not the default))\n"+
		" OR\n"+
		"   (1) raw cnvs to parse (i.e. rawcnv=final.rawcnv (not the default))\n"+
		"   (2) (optional) parse only de novo variants (i.e. -denovoOnly (not the default))\n"+
		" OR\n" +
		"   (1) a comma-separated list of .cnv files to combine together (i.e. combine=/full/path/to/cnv1.cnv,relative/path/to/cnv2.cnv (not the default))\n" +
		"   (2) full path of the desired output file (i.e. output=/path/to/output/file.cnv (not the default))\n" + 
		"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				return;
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("batch=")) {
				batch = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("-qsub")) {
				qsub = true;
				numArgs--;
			} else if (args[i].startsWith("rawlog=")) {
				rawlog = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("rawcnv=")) {
				rawcnvs = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("-pfb")) {
				parsePFB = true;
				numArgs--;
			} else if (args[i].startsWith("gc5base=")) {
				gc5base = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("-denovoOnly")) {
				denovoOnly = true;
				numArgs--;
			} else if (args[i].startsWith("pfb=")) {
				pfbFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("gcmodel=")) {
				gcmodelFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("data=")) {
				transformData = Boolean.parseBoolean(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("auto=")) {
				auto = Boolean.parseBoolean(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("chrx=")) {
				chrx = Boolean.parseBoolean(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("sexSpecificCentroids=")) {
				sexCent = Boolean.parseBoolean(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("combine=")) {
				cnvFiles = args[i].split("=")[1].split(",");
				numArgs--;
			} else if (args[i].startsWith("output=")) {
				outputFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("-recode")) {
				recode = true;
				numArgs--;
			} else if (args[i].startsWith("-submit")) {
				submit = true;
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			return;
		}
		try {

//			filename = "C:/workspace/Genvisis/projects/GEDI_exome.properties";
//			parsePFB = true;
//			gc5base = "C:/projects/gcModel/gc5Base.txt";

//			logfile = "penncnv/penncnv.log";
//			rawcnvs = "penncnv/penncnv.rawcnv";
			
//			filename = "/home/npankrat/projects/GEDI.properties";
//			batch = 60;
//			qsub = true;
//			pfbFile = "gedi.pfb";
//			gcmodelFile = "gedi.gcmodel";
			
			proj = new Project(filename, logfile, false);
			if (parsePFB) {
				populationBAF(proj);
			}
			if (gc5base != null) {
				gcModel(proj, gc5base, proj.getProjectDir()+"custom.gcmodel", 100);
			}
			if (batch > 0) {
				doBatch(proj, auto, chrx, sexCent, transformData, batch, qsub, pfbFile, gcmodelFile, qsub ? submit : false, recode);
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
