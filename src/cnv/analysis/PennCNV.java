package cnv.analysis;

import java.io.*;
import java.util.*;

import javax.swing.JOptionPane;
//import cnv.analysis.FilterCalls;
import cnv.filesys.FullSample;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.var.CNVariant;
import cnv.var.SampleData;
import common.*;

public class PennCNV {
	public static final String[] QC_HEADS = {"LRR_mean", "LRR_median", "LRR_SD", "BAF_mean", "BAF_median", "BAF_SD", "BAF_DRIFT", "WF", "GCWF"};
	public static final String[] ERRORS = {"large SD for LRR", "drifting BAF values", "waviness factor values", "Small-sized CNV calls"};
	public static final String QC_SUMMARY_FILE = "Sample_QC.xln";

	public static void batch(Project proj, int numBatches, boolean createLists) {
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
					writer.println(files[j]);
				}
				writer.close();
			} catch (Exception e) {
				System.err.println("Error writing to list"+(i+1)+".txt");
				e.printStackTrace();
			}
		}

		init = "cd "+dataDir;
		commands = execDir+"detect_cnv.pl -test -conf -hmm "+execDir+"lib/hhall.hmm -pfb "+execDir+"lib/hhall.hg18.pfb -gcmodel "+execDir+"lib/hhall.hg18.gcmodel -list ../penncnv/list[%0].txt -log ../penncnv/[%0].log -out ../penncnv/[%0].rawcnv > ../penncnv/[%0].out";

		Files.batchIt("penn", init, numBatches, commands, Array.stringArraySequence(numBatches, ""));
	}

	public static void parseWarnings(Project proj, String filename) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, data;
		String temp, trav = null;
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
				line = temp.trim().split("[\\s]+");
				try {
					if (temp.contains("quality summary")) {
						trav = line[4].substring(line[4].lastIndexOf("/")+1, line[4].indexOf(":"));
						v.add(trav);
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
						hash.put(trav, data);
					} else if (temp.startsWith("WARNING")) {
						if (temp.contains("Small-sized CNV calls")) {
							// use old trav
						} else {
							trav = line[3].substring(line[3].lastIndexOf("/")+1);
						}
						data = hash.get(trav);
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
				trav = v.elementAt(keys[i]);
				data = hash.get(trav);
				writer.print(trav+"\t"+sampleData.lookup(trav));
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

	public static void parseResults(Project proj, String filename) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String temp, trav;
		Vector<String> warnings = new Vector<String>();
		Vector<String> pedinfo = new Vector<String>();
		int[] position;
		String score;
		SampleData sampleData;
		String ind;
		Hashtable<String,String> hash;

		sampleData = proj.getSampleData(false);
		try {
			reader = new BufferedReader(new FileReader(proj.getProjectDir()+filename));
			writer = new PrintWriter(new FileWriter(proj.getProjectDir()+ext.rootOf(filename)+".cnv"));
			writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER));
			hash = new Hashtable<String,String>();
			while (reader.ready()) {
				temp = reader.readLine();
				line = temp.trim().split("[\\s]+");
				position = Positions.parseUCSClocation(line[0]);
				trav = line[4];
				trav = trav.substring(trav.lastIndexOf("/")+1);
				ind = sampleData.lookup(trav);
				if (ind == null) {
					if (!hash.containsKey(trav)) {
						System.err.println("Error - '"+trav+"' was not found in "+proj.getFilename(Project.SAMPLE_DATA_FILENAME));
						hash.put(trav, "");
					}
					ind = trav+"\t"+trav;
				}

				if (line.length<8||!line[7].startsWith("conf=")) {
					score = "-1";
					if (!warnings.contains(trav)) {
						System.err.println("Warning - no conf estimates for "+trav);
						warnings.add(trav);
					}
				} else {
					score = ext.formDeci(Double.parseDouble(line[7].substring(5)), 4, true);
				}
				writer.println(ind+"\t"+position[0]+"\t"+position[1]+"\t"+position[2]+"\t"+line[3].substring(line[3].indexOf("=")+1)+"\t"+score+"\t"+line[1].substring(7));
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

	public static void populationBAF(Project proj) {
		PrintWriter writer;
		FullSample samp;
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
		
		filename = proj.getFilename(Project.SAMPLE_SUBSET_FILENAME, true, false);
		
		if (ext.rootOf(filename) == null || ext.rootOf(filename).equals("")) {
			sampleList = proj.getSampleList().getSamples();
			output = proj.getProjectDir()+"custom.pfb";
		} else if (Files.exists(filename, proj.getJarStatus())) {
			System.out.print("filename: "+filename);
			sampleList = HashVec.loadFileToStringArray(filename, false, new int[] {0}, false);
			output = proj.getProjectDir()+ext.rootOf(filename)+".pfb";
		} else {
			JOptionPane.showMessageDialog(null, "Failed to load \""+filename+"\"", "Error", JOptionPane.ERROR_MESSAGE);
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
			samp = proj.getFullSample(sampleList[i]);
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
				writer.println(markerNames[i]+"\t"+chrs[i]+"\t"+positions[i]+"\t"+bafAverage[i]);
			}
			writer.close();
			System.out.println("Population BAF file is now ready at: "+output);
		} catch (Exception e) {
			System.err.println("Error writing to '" + output + "'");
			e.printStackTrace();
		}
	}

	public static void gcModel(Project proj, String gcFile, String outputFile, int numwindow) {
		MarkerSet markerSet;
		String[] markerNames;
		byte[] chrs;
		int[] positions;

		BufferedReader reader;
		String[] line;
		double[] gc;
		PrintWriter writer;
		byte curchr;
		int curstart, curend, curcount;
		double cursum;
		Hashtable<Byte, Byte> seen_chr = new Hashtable<Byte, Byte>();

		int[] snp_count;
		double[] snp_sum;
		byte prechr     = -1;
		int prestart    = -1;
		int preend      = 0;
		int windowcount = 0;
		int windowsum   = 0;

		int chr_index;
//		byte[] seen_chr;
//		int count= 0;
		
//		alldata[];						//to be used in both the SnpFile (.pdf file) block and the GcFile block.
		
		// generate or load SnpFile (pbf file or Population B Allele Frequency)
		markerSet = proj.getMarkerSet();
		markerNames = markerSet.getMarkerNames();
		chrs = markerSet.getChrs();	//to be used only in the SnpFile (.pdf file) block. But then in the outputFile.
		positions = markerSet.getPositions();	//to be used only in the SnpFile (.pdf file) block. But then in the GcFile block.
//		while (SnpFile.length) {
//			push @{$alldata{$curchr}}, [$curpos, $curname];
//			chrs[i] = curchr;
//			positions[i] = curpos;
//		}

		// How do we know whether "numwindow==null" ???
		if (numwindow==0) {
			numwindow = 100;
		}

		snp_count = new int[markerNames.length];
		snp_sum = new double[markerNames.length];
		chr_index=0;
		// load gcFile
		try {
			reader = new BufferedReader(new FileReader(gcFile));
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				curchr = Positions.chromosomeNumber(line[1]);
				curstart = Integer.parseInt(line[2]);
				curend = Integer.parseInt(line[3]);
				curcount = Integer.parseInt(line[11]);
				cursum = Double.parseDouble(line[12]);
				
				if (curchr > 0) { //skip the random chromosomes!
					//skip these annotations for 1Mb sliding windows (since we only care about 5kb window)
//					m/1K.wib/ and next;
//					curdata = alldata[curchr];
//					if (curdata) {
					
					if (curchr == prechr) {
						if (curstart < prestart) {
							System.err.println("Error in gcFile: a record in chr"+curchr+" has position "+curstart+", less then the previous position $prestart");
							System.exit(1);
						}
					} else if (seen_chr.containsKey(curchr)) {
							System.err.println("Error in gcFile: rows of the same chromosome must be adjacent. But now chr"+curchr+" occur multiple times in non-continuous segment of the "+gcFile+": at "+curchr+":"+curstart);
							System.exit(1);
					} else {
						seen_chr.put(curchr,(byte)1);
					}
					
//					chr_index[curchr]= 0;
//					(curdata-1);

					if (chrs[chr_index]==curchr) {
						while (positions[chr_index]<Math.max(curstart - numwindow*5120, 0)) {
							chr_index ++;						
						}
					} else if (chrs[chr_index]<curchr) {
						while (chr_index<chrs.length && chrs[chr_index]!=curchr) {
							chr_index ++;
						}
						if (chr_index==chrs.length) {
							chr_index --;
						}
					} else {
						chr_index = 0;
						while (chr_index<chrs.length && chrs[chr_index]!=curchr) {
							chr_index ++;
						}
						if (chr_index==chrs.length) {
							chr_index --;
						}
					}
					for (int i = chr_index; i<snp_count.length && chrs[i]==curchr && positions[i]<=(curend + numwindow*5120); i++) {
						snp_count[i] += curcount;
						snp_sum[i] += cursum;
					}
//					}
//					count++;
				}
				prestart = curstart;
				prechr = curchr;
				}
			reader.close();
		} catch (Exception e) {
			System.err.println("Error reading from '" + gcFile + "'");
			e.printStackTrace();
			System.exit(0);
		}
		
		// load pfb file or generate it
		
		// output the result
		try {
			writer = new PrintWriter(new FileWriter(outputFile));
			writer.println("Name\tChr\tPosition\tGC");
			for (int i = 0; i<markerNames.length; i++) {
				writer.println(markerNames[i]+"\t"+chrs[i]+"\t"+positions[i]+"\t"+(snp_count[i]==0?(snp_sum[i]==0?0:"err"):(snp_sum[i]/snp_count[i])));
			}
			writer.close();
			System.out.println("Population BAF file is now ready at: "+outputFile);
		} catch (Exception e) {
			System.err.println("Error writing to '" + outputFile + "'");
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = Project.DEFAULT_PROJECT;
//		String logfile = "conf.log";
		String logfile = "penncnv/penncnv.log";
//		String logfile = "";
//		String rawcnvs = "conf.rawcnv";
//		String rawcnvs = "penncnv.rawcnv";
		String rawcnvs = null;
		int batch = 0;
		boolean lists = false;
		Project proj;
		boolean parsePFB = false;

		String usage = "\n"+
		"cnv.analysis.PennCNV requires 0-1 arguments\n"+
		"   (1) project properties filename (i.e. proj="+filename+" (default))\n"+
		"   (2) number of batches to do (i.e. batch=12 (not the default))\n"+
		"   (3) create lists for batches (i.e. -lists (not the default))\n"+
		" OR\n"+
		"   (2) compute populationo frequence of b allele file (using paramters in properties file) (i.e. -pfb (not the default))\n"+
		" OR\n"+
		"   (1) parse warnings from log file (i.e. log=final.log (not the default))\n"+
		" OR\n"+
		"   (1) raw cnvs to parse (i.e. raw=final.rawcnv (not the default))\n"+
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
			} else if (args[i].startsWith("log=")) {
				logfile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("raw=")) {
				rawcnvs = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("-pfb")) {
				parsePFB = true;
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			proj = new Project(filename, false);
			
			if (batch>0) {
				batch(proj, batch, lists);
			} else if (parsePFB) {
				populationBAF(proj);
			} else if (logfile != null) {
				parseWarnings(proj, logfile);
			} else if (rawcnvs != null) {
				parseResults(proj, rawcnvs);
			}

		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
