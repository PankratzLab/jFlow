package cnv.analysis;

import java.io.*;
import java.util.*;
//import cnv.analysis.FilterCalls;
import cnv.filesys.Project;
import cnv.var.CNVariant;
import cnv.var.SampleData;
import common.*;

public class PennCNV {
	public static final String[] QC_HEADS = {"LRR_mean", "LRR_median", "LRR_SD", "BAF_mean", "BAF_median", "BAF_SD", "BAF_DRIFT", "WF", "GCWF"};
	public static final String[] ERRORS = {"large SD for LRR", "drifting BAF values", "waviness factor values", "Small-sized CNV calls"};
	public static final String QC_SUMMARY_FILE = "Sample_QC.xln";

	public static void batch(Project proj, int numBatches, boolean createLists, boolean qsub, String pfbFile, String gcmodelFile) {
		String init, commands;
		PrintWriter writer;
		String[] files;
		int step;
		String execDir, dataDir, resultsDir;
		
		
		execDir = proj.getDir(Project.PENNCNV_EXECUTABLE_DIRECTORY);
		dataDir = proj.getDir(Project.PENNCNV_DATA_DIRECTORY);
		resultsDir = proj.getDir(Project.PENNCNV_RESULTS_DIRECTORY);

		files = new File(dataDir).list(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				return file.length()>1000;
			}
		});
		System.out.println("Found "+files.length+" files");

		step = (int)Math.ceil((double)files.length/(double)numBatches);
		System.out.println("Which means the step for "+numBatches+" batches would be "+step);
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
				System.err.println("Error writing to list"+(i+1)+".txt");
				e.printStackTrace();
			}
		}

		init = "cd "+dataDir;
		commands = execDir+"detect_cnv.pl -test -conf -hmm "+execDir+"lib/hhall.hmm -pfb "+(pfbFile==null?execDir+"lib/hhall.hg18.pfb":pfbFile)+" -gcmodel "+(gcmodelFile==null?execDir+"lib/hhall.hg18.gcmodel":gcmodelFile)+" -list ../penncnv/list[%0].txt -log ../penncnv/[%0].log -out ../penncnv/[%0].rawcnv > ../penncnv/[%0].out";

		if (qsub) {
			Files.qsub("runPenn", init+"\n"+commands, Matrix.toMatrix(Array.stringArraySequence(numBatches, "")));
		} else {
			Files.batchIt("penn", init, numBatches, commands, Array.stringArraySequence(numBatches, ""));
		}
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
		
		sampleData = proj.getSampleData(false);
		lrrSD_cutoff = proj.getDouble(Project.LRRSD_CUTOFF);
				
		try {
			reader = new BufferedReader(new FileReader(proj.getProjectDir()+filename));
			while (reader.ready()) {
				temp = reader.readLine();
				temp = translateDerivedSamples(temp);
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
							System.err.println("Error - line doesn't have all the expected pieces:");
							System.err.println(temp);
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
							System.err.println("Unknown WARNING: "+temp);
						} else {
							data[err] = "1";
						}
					}
				} catch (Exception e) {
					System.err.println(temp);
					e.printStackTrace();

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
				writer.print(sampleID+"\t"+sampleData.lookup(sampleID));
				writer.print("\t"+(data[1].equals("1")||data[2].equals("1")||Double.parseDouble(data[6])>lrrSD_cutoff?"0":"1"));
				writer.println("\t"+Array.toStr(data));
			}
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+proj.getProjectDir()+filename+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+proj.getProjectDir()+filename+"\"");
			System.exit(2);
		}
	}
	
	public static String translateDerivedSamples(String str) {
		String trav;
		int start, stop;
		
		start = str.indexOf("`");
		stop = str.lastIndexOf("`");
		if (start == -1) {
			return str;
		}
		
		trav = str.substring(start+1, stop); 
		if (trav.contains("`")) {
			System.err.println("Error - more than one set of quotes for: "+str);
		}
		if (trav.startsWith("gunzip -c") && trav.endsWith(".gz")) {
			trav = trav.substring(9, trav.length()-3).trim();
		} else {
			System.err.println("Error - not currently set up to handle the following construction into a sample_ID: "+trav);
		}
		
		return str.substring(0, start)+trav+str.substring(stop+1);		
	}

	public static void parseResults(Project proj, String filename, boolean all, boolean denovoOnly) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String temp, trav;
		Vector<String> warnings = new Vector<String>();
		Vector<String> pedinfo = new Vector<String>();
		int[] position;
		String score;
		SampleData sampleData;
		String famIndPair;
		Hashtable<String,String> hash;

		sampleData = proj.getSampleData(false);
		try {
			reader = new BufferedReader(new FileReader(proj.getProjectDir()+filename));
			writer = new PrintWriter(new FileWriter(proj.getProjectDir()+ext.rootOf(filename)+".cnv"));
			writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER));
			hash = new Hashtable<String,String>();
			while (reader.ready()) {
				temp = reader.readLine();
				if (!temp.startsWith("NOTICE:")) {
					temp = translateDerivedSamples(temp);
					line = temp.trim().split("[\\s]+");
					position = Positions.parseUCSClocation(line[0]);
					trav = line[4];
					trav = trav.substring(trav.lastIndexOf("/")+1);
					famIndPair = sampleData.lookup(trav);
					if (famIndPair == null) {
						if (!hash.containsKey(trav)) {
							System.err.println("Error - '"+trav+"' was not found in "+proj.getFilename(Project.SAMPLE_DATA_FILENAME));
							hash.put(trav, "");
						}
						famIndPair = trav+"\t"+trav;
					}
	
					if (line.length<8||!line[7].startsWith("conf=")) {
						score = "-1";
						if (!warnings.contains(trav) && warnings.size() < 10) {
							System.err.println("Warning - no conf estimates for "+trav);
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
			for (int i = 0; i<pedinfo.size(); i++) {
				writer.println(pedinfo.elementAt(i));
			}
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+proj.getProjectDir()+ext.rootOf(filename)+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+proj.getProjectDir()+ext.rootOf(filename)+"\"");
			System.exit(2);
		}
	}
	
	public static void createFromParameters(String filename, Logger log) {
		Vector<String> params;

		params = Files.parseControlFile(filename, "penncnv", new String[] { "file=snps.txt", "option1=", "out=finalProduct.out" }, log);

		if (params != null) {
			params.add("log=" + log.getFilename());
			main(Array.toStringArray(params));
		}
	}
	

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = Project.DEFAULT_PROJECT;
//		String logfile = "conf.log";
//		String logfile = "penncnv/penncnv.log";
		String logfile = "";
//		String rawcnvs = "conf.rawcnv";
//		String rawcnvs = "penncnv.rawcnv";
		String rawcnvs = "";
		int batch = 0;
		boolean all = false;
		boolean lists = false;
		boolean qsub = false;
		Project proj;
		String pfbFile = null;
		String gcmodelFile = null;
		boolean denovoOnly = false;

		String usage = "\n"+
		"cnv.park.PennCNV requires 0-1 arguments\n"+
		"   (0) project properties filename (i.e. proj="+filename+" (default))\n"+
		"   (1) number of batches to do (i.e. batch=12 (not the default))\n"+
		"   (2) generate qsub files instead of batch files (i.e. -qsub (not the default))\n"+
		"   (3) create lists for batches (i.e. -lists (not the default))\n"+
		"   (4) log file (i.e. log=final.cnv (not the default))\n"+
		"   (5) raw cnvs to parse (i.e. raw=final.rawcnv (not the default))\n"+
		"   (6) use all indiviudals (i.e. -all (not the default))\n"+
		"   (7) (optional) use custom pfb file (i.e. pfb=custom.pfb (not the default))\n"+
		"   (8) (optional) use custom gcmodel file (i.e. gcmodel=custom.gcmodel (not the default))\n"+
		"   (9) (optional) parse only de novo variants (i.e. -denovoOnly (not the default))\n"+
		"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
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
			} else if (args[i].startsWith("log=")) {
				logfile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("raw=")) {
				rawcnvs = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("-all")) {
				all = true;
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
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			proj = new Project(filename, false);

//			rawcnvs = "penncnv/penncnv.rawcnv";
			
			if (batch>0) {
				batch(proj, batch, lists, qsub, pfbFile, gcmodelFile);
			}
			if (!logfile.equals("")) {
				parseWarnings(proj, logfile);
			}
			if (!rawcnvs.equals("")) {
				parseResults(proj, rawcnvs, all, denovoOnly);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
