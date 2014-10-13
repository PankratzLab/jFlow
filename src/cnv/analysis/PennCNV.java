package cnv.analysis;

import java.io.*;
import java.util.*;

import cnv.filesys.MarkerSet;
//import cnv.analysis.FilterCalls;
import cnv.filesys.Project;
import cnv.filesys.Sample;
import cnv.var.CNVariant;
import cnv.var.SampleData;
import common.*;

public class PennCNV {
	public static final String[] QC_HEADS = {"LRR_mean", "LRR_median", "LRR_SD", "BAF_mean", "BAF_median", "BAF_SD", "BAF_DRIFT", "WF", "GCWF"};
	public static final String[] ERRORS = {"large SD for LRR", "drifting BAF values", "waviness factor values", "Small-sized CNV calls", "NoCall rate"};
	public static final String QC_SUMMARY_FILE = "Sample_QC.xln";

	public static void batch(Project proj, int numBatches, boolean createLists, boolean qsub, String pfbFile, String gcmodelFile) {
		String init, commands;
		PrintWriter writer;
		String[] files;
		int step;
		String execDir, dataDir, resultsDir;
		Logger log;
		
		log = proj.getLog();
		execDir = proj.getDir(Project.PENNCNV_EXECUTABLE_DIRECTORY);
		dataDir = proj.getDir(Project.PENNCNV_DATA_DIRECTORY);
		resultsDir = proj.getDir(Project.PENNCNV_RESULTS_DIRECTORY);
		
		if (pfbFile != null) {
			pfbFile = ext.replaceTilde(pfbFile);
			if (!pfbFile.startsWith("/")) {
				pfbFile = ext.pwd()+pfbFile;
			}
			if (!Files.exists(pfbFile)) {
				log.reportError("Error - pfb file '"+pfbFile+"' does not exist; aborting");
				return;
			}
		}

		if (gcmodelFile != null) {
			gcmodelFile = ext.replaceTilde(gcmodelFile);
			if (!gcmodelFile.startsWith("/")) {
				gcmodelFile = ext.pwd()+gcmodelFile;
			}
			if (!Files.exists(gcmodelFile)) {
				log.reportError("Error - gcmodel file '"+gcmodelFile+"' does not exist; aborting");
				return;
			}
		}

		files = new File(dataDir).list(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				return file.length()>1000;
			}
		});
		log.report("Found "+files.length+" files");

		step = (int)Math.ceil((double)files.length/(double)numBatches);
		log.report("Which means the step for "+numBatches+" batches would be "+step);
		new File(resultsDir).mkdir();
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
				log.reportError("Error writing to list"+(i+1)+".txt");
				log.reportException(e);
			}
		}

		init = "cd "+dataDir;
		commands = execDir+"detect_cnv.pl -test -conf -hmm "+execDir+"lib/hhall.hmm -pfb "+(pfbFile==null?execDir+"lib/hhall.hg18.pfb":pfbFile)+" -gcmodel "+(gcmodelFile==null?execDir+"lib/hhall.hg18.gcmodel":gcmodelFile)+" -list "+resultsDir+"list[%0].txt -log "+resultsDir+"[%0].log -out "+resultsDir+"[%0].rawcnv > "+resultsDir+"[%0].out";

		if (qsub) {
			Files.qsub("runPenn", dataDir, numBatches, commands, Matrix.toMatrix(Array.stringArraySequence(numBatches, "")), 2200, 8);
		} else {
			Files.batchIt("penn", init, numBatches, commands, Array.stringArraySequence(numBatches, ""));
		}
		Files.writeList(new String[] {
				"cat penncnv/*.log > penncnv.rawlog",
				"cat penncnv/*.rawcnv > penncnv.rawcnv",
				"java -cp ~/park.jar cnv.analysis.PennCNV proj="+proj.getPropertyFilename()+" rawlog=penncnv.rawlog",
				"java -cp ~/park.jar cnv.analysis.PennCNV proj="+proj.getPropertyFilename()+" rawcnv=penncnv.rawcnv",
		}, "assemblePenncnv");
		Files.chmod("assemblePenncnv");
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
			reader = new BufferedReader(new FileReader(proj.getProjectDir()+filename));
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
			log.reportError("Error: file \""+proj.getProjectDir()+filename+"\" not found in current directory");
			return;
		} catch (IOException ioe) {
			log.reportError("Error reading file \""+proj.getProjectDir()+filename+"\"");
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
		
		if (!Files.exists(proj.getProjectDir()+filename)) {
			log.reportError("Error - could not find file '"+proj.getProjectDir()+filename+"'");
			return;
		}

		warnings = new Vector<String>();
		sampleData = proj.getSampleData(2, false);
		pedinfo = new Hashtable<String, Vector<String>>();
		try {
			reader = new BufferedReader(new FileReader(proj.getProjectDir()+filename));
			writer = new PrintWriter(new FileWriter(proj.getProjectDir()+ext.rootOf(filename)+".cnv"));
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
	
					if (line.length<8||!line[7].startsWith("conf=")) {
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

			writer = new PrintWriter(new FileWriter(proj.getProjectDir()+ext.rootOf(filename)+".fam"));
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
			log.reportError("Error: file \""+proj.getProjectDir()+ext.rootOf(filename)+".cnv\" not found in current directory");
			return;
		} catch (IOException ioe) {
			log.reportError("Error reading file \""+proj.getProjectDir()+ext.rootOf(filename)+"\"");
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

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		String rawlog = null;
		String rawcnvs = null;
		int batch = 0;
//		boolean all = false;
		boolean lists = false;
		boolean qsub = false;
		Project proj;
		String pfbFile = null;
		String gcmodelFile = null;
		boolean denovoOnly = false;
		boolean parsePFB = false;
		String gc5base = null;
		String logfile = null;

		String usage = "\n"+
		"cnv.park.PennCNV requires 0-1 arguments\n"+
		"   (0) project properties filename (i.e. proj="+cnv.Launch.getDefaultDebugProjectFile(false)+" (default))\n"+
		"   (1) number of batches to do (i.e. batch=12 (not the default))\n"+
		"   (2) generate qsub files instead of batch files (i.e. -qsub (not the default))\n"+
		"   (3) create lists for batches (i.e. -lists (not the default))\n"+
		"   (4) (optional) use custom pfb file (i.e. pfb=custom.pfb (not the default))\n"+
		"   (5) (optional) use custom gcmodel file (i.e. gcmodel=custom.gcmodel (not the default))\n"+
		" OR\n"+
		"   (1) compute file containing project based b allele frequencies for file using parameters in properties file (i.e. -pfb (not the default))\n"+
		" OR\n"+
		"   (1) compute a custom gcmodel file for the markers in this project using this file (i.e. gc5base=gc5base.txt (not the default))\n"+
		" OR\n"+
		"   (1) parse warnings from log file (i.e. rawlog=final.log (not the default))\n"+
		" OR\n"+
		"   (1) raw cnvs to parse (i.e. rawcnv=final.rawcnv (not the default))\n"+
		"   (2) (optional) parse only de novo variants (i.e. -denovoOnly (not the default))\n"+
//		"   (3) use all individuals (i.e. -all (not the default))\n"+ // ?? does not currently do anything
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
			} else if (args[i].startsWith("-lists")) {
				lists = true;
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
//			} else if (args[i].startsWith("-all")) {
//				all = true;
//				numArgs--;
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
			if (batch>0) {
				batch(proj, batch, lists, qsub, pfbFile, gcmodelFile);
			}
			if (rawlog != null) {
				parseWarnings(proj, rawlog);
			}
			if (rawcnvs != null) {
				parseResults(proj, rawcnvs, denovoOnly);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
