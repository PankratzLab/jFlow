// quantile normalize

package cnv.analysis;

import java.io.*;
import java.util.Hashtable;
import java.util.Vector;

import cnv.filesys.*;
import common.Array;
import common.Files;
import common.HashVec;
import common.Matrix;
import common.ext;
import filesys.Segment;

public class AnalysisFormats implements Runnable {
	public static final String[] PROGRAM_OPTIONS = {"PennCNV", "QuantiSNP"};
	public static final String[] OUTPUT_DIRECTORIES = {"results/PennCNV/", "results/QuantiSNP/"};
	public static final int PENN_CNV = 1;
	public static final int QUANTISNP = 2;
	// public static final int EM_ITERATIONS = 10;
	public static final int EM_ITERATIONS = 25;
	private Project proj;
	private String[] samples;
	private int program;
	private Hashtable<String,String> hash;

	public AnalysisFormats(Project proj, String[] samples, int program, Hashtable<String,String> hash) {
		this.proj = proj;
		this.samples = samples;
		this.program = program;
		this.hash = hash;
	}

	public void run() {
		switch (program) {
		case PENN_CNV:
			penncnv(proj, samples, hash);
			break;
		case QUANTISNP:
			quantisnp(proj, samples, hash);
			break;
		default:
			System.err.println("Error - invalid program option: "+program);
			break;
		}

	}

	public static void penncnv(Project proj, String[] samples, Hashtable<String,String> hash) {
		PrintWriter writer;
		String[] markerNames = proj.getMarkerNames();
		Sample samp;
		float[] lrrs, bafs;
		byte[] genotypes;
		boolean jar, gzip;
		String dir;

		dir = proj.getDir(Project.PENNCNV_DATA_DIRECTORY);
		new File(dir).mkdirs();
		jar = proj.getJarStatus();
		gzip = proj.getBoolean(Project.PENNCNV_GZIP_YESNO);
		for (int i = 0; i<samples.length; i++) {
			System.out.println(ext.getTime()+"\tTransforming "+(i+1)+" of "+samples.length);
			if (Files.exists(proj.getDir(Project.SAMPLE_DIRECTORY) + samples[i] + Sample.SAMPLE_DATA_FILE_EXTENSION, jar)) {
//				samp = Sample.loadFromRandomAccessFile(proj.getDir(Project.SAMPLE_DIRECTORY) + samples[i] + Sample.SAMPLE_DATA_FILE_EXTENSION, false, false, true, true, false, jar);
				samp = Sample.loadFromRandomAccessFile(proj.getDir(Project.SAMPLE_DIRECTORY) + samples[i] + Sample.SAMPLE_DATA_FILE_EXTENSION, false, false, true, true, true, jar);
			} else {
//				System.err.println("Error - the " + Sample.SAMPLE_DATA_FILE_EXTENSION + " file is not found for the following item in sampleList " + samples[i]);
				System.err.println("Error - the " + samples[i] + Sample.SAMPLE_DATA_FILE_EXTENSION + " is not found.");
				return;
			}
//			samp = proj.getPartialSampleFromRandomAccessFile(samples[i]);
			lrrs = samp.getLRRs();
			bafs = samp.getBAFs();
			genotypes = samp.getAB_Genotypes();

			try {
				writer = Files.getAppropriateWriter(dir+samples[i]+(gzip?".gz":""));
				writer.println("Name\t"+samples[i]+".GType\t"+samples[i]+".Log R Ratio\t"+samples[i]+".B Allele Freq");
				for (int j = 0; j<markerNames.length; j++) {
					if (hash == null || hash.containsKey(markerNames[j])) {
						writer.println(markerNames[j]+"\t"+(genotypes[j]==-1?"NC":Sample.AB_PAIRS[genotypes[j]])+"\t"+lrrs[j]+"\t"+bafs[j]);
					}
				}
				writer.close();
			} catch (Exception e) {
				System.err.println("Error writing PennCNV data for "+samples[i]);
				e.printStackTrace();
			}
		}
	}

	public static void quantisnp(Project proj, String[] samples, Hashtable<String,String> hash) {
		PrintWriter writer;
		MarkerSet set;
		String[] markerNames = proj.getMarkerNames();
		Sample samp;
		float[] lrrs, bafs;
		byte[] chrs;
		int[] positions;
		
		set = proj.getMarkerSet();
		chrs = set.getChrs();
		positions = set.getPositions();

		new File(proj.getProjectDir()+"quanti_data/").mkdirs();
		for (int i = 0; i<samples.length; i++) {
			System.out.println(ext.getTime()+"\tTransforming "+(i+1)+" of "+samples.length);
			samp = proj.getPartialSampleFromRandomAccessFile(samples[i]);
			set.checkFingerprint(samp);
			lrrs = samp.getLRRs();
			bafs = samp.getBAFs();
			
			try {
				writer = new PrintWriter(new FileWriter(proj.getProjectDir()+"quanti_data/"+samples[i]));
				writer.println("Name\tChr\tPosition\t"+samples[i]+".Log R Ratio\t"+samples[i]+".B Allele Freq");
				for (int j = 0; j<markerNames.length; j++) {
					if (hash == null || hash.containsKey(markerNames[j])) {
						writer.println(markerNames[j]+"\t"+chrs[j]+"\t"+positions[j]+"\t"+lrrs[j]+"\t"+bafs[j]);
					}
				}
				writer.close();
			} catch (Exception e) {
				System.err.println("Error writing QuantiSNP data for "+samples[i]);
				e.printStackTrace();
			}
		}
	}
	
	public static void batchQuantiSNP(Project proj, int numBatches) {
		Vector<String[]> v = new Vector<String[]>();
		Hashtable<String,String> genders;
		String[] inputs, outputs;
		String commands, gender;

		genders = HashVec.loadFileToHashString(proj.getFilename(Project.SAMPLE_DATA_FILENAME), "DNA", new String[] {"CLASS=Gender"}, "");
		
		inputs = new File(proj.getProjectDir()+"quanti_data/").list(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				return filename.endsWith(".qs");
			}
		});
		
		if (inputs == null) {
			System.err.println("Error - QuantiSNP inputs files have not yet been created");
		}

		outputs = new File(proj.getProjectDir()+"results/QuantiSNP/").list(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				return filename.endsWith("_output.out");
			}
		});
		
		if (outputs == null) {
			System.out.println("Found "+inputs.length+" samples; creating output directory for QuantiSNP in results/QuantiSNP/");
			outputs = new String[0];
		} else {
			System.out.println("Found "+inputs.length+" samples, as well as results for "+outputs.length+" that have been done (not necessarily the same ones)");
		}
		
		for (int i = 0; i<inputs.length; i++) {
			if (ext.indexOfStr(ext.rootOf(inputs[i])+"_output.out", outputs)==-1) {
				if (genders.containsKey(ext.rootOf(inputs[i]))) {
					gender = genders.get(ext.rootOf(inputs[i]));
					if (gender.equals("M") || gender.equals("1")) {
						gender = "male";
					} else if (gender.equals("F") || gender.equals("2")) {
						gender = "female";
					} else {
						System.err.println("Error - '"+gender+"' is not a valid gender (expecting M/F or 1/2)");
					}
				} else {
					System.err.println("Error - no gender found for subject '"+ext.rootOf(inputs[i])+"'");
					gender = null;
				}					

				v.add(new String[] {ext.rootOf(inputs[i]), gender});
			}
		}

		System.out.println("Made "+numBatches+" batch files that will take care of the "+v.size()+" files yet to parse");

//		commands = "quantisnp.exe --config ../windows/config.dat  --emiters "+EM_ITERATIONS+" --Lsetting 2000000 --maxcopy 3 --printRS --doGCcorrect --gcdir ../gc/b36/ --output "+OUTPUT_DIRECTORIES[1]+"[%0].out --gender [%1]--input-files ../source/[%0].qs 300\n\n";
		commands = "quantisnp --output "+OUTPUT_DIRECTORIES[1]+"[%0].out --gender [%1] --input-files ../source/[%0].qs 300\n\n";
		Files.batchIt("batch", null, numBatches, commands, Matrix.toStringArrays(v));
	}	

	public static void launch(Project proj, int program, String markers, int numThreads) {
		Vector<Vector<String>> sampleLists = new Vector<Vector<String>>();
		String[] samples = proj.getSamples();
		Thread[] threads;
		Hashtable<String,String> hash;

		for (int i = 0; i<numThreads; i++) {
			sampleLists.add(new Vector<String>());
		}
		for (int i = 0; i<samples.length; i++) {
			sampleLists.elementAt(i%numThreads).add(samples[i]);
		}
		if (markers == null) {
			hash = null;
		} else {
			hash = HashVec.loadFileToHashNull(proj.getProjectDir()+markers, false);
		}
		threads = new Thread[numThreads];
		for (int i = 0; i<numThreads; i++) {
			threads[i] = new Thread(new AnalysisFormats(proj, Array.toStringArray(sampleLists.elementAt(i)), program, hash));
			threads[i].start();
			try {
				Thread.sleep(100L);
			} catch (InterruptedException ex) {}
		}
	}
	
	public static void filter(Project proj, String regions, String list, String outfile) {
		PrintWriter writer;
		MarkerSet markers;
		Segment[] segs;
		String[] markerNames;
		byte[] chrs;
		int[] positions;
		int countFromList, countInRegions, countOverlap;
		Hashtable<String,String> hash;
		
		if (outfile == null) {
			System.err.println("Error - outfile is defined as null; need to provide a filename before results can be filtered");
			return;
		}
		
		if (regions.equals("")) {
			segs = new Segment[0];
		} else {
			segs = Segment.loadUCSCregions(proj.getProjectDir()+regions, false);
		}

		if (list.equals("")) {
			hash = new Hashtable<String,String>();
		} else {
			hash = HashVec.loadFileToHashNull(proj.getProjectDir()+list, false);
		}

		markers = proj.getMarkerSet();
		markerNames = markers.getMarkerNames();
		chrs = markers.getChrs();
		positions = markers.getPositions();
		
		try {
	        writer = new PrintWriter(new FileWriter(proj.getProjectDir()+outfile));
	        countFromList = countInRegions = countOverlap = 0;
			for (int i = 0; i<markerNames.length; i++) {
				if (segs.length > 0 && Segment.overlapsAny(new Segment(chrs[i], positions[i], positions[i]), segs)) {
					countInRegions++;
					if (hash.containsKey(markerNames[i])) {
						countFromList++;
						countOverlap++;
					}
				} else if (hash.containsKey(markerNames[i])) {
					countFromList++;
				} else {
					writer.println(markerNames[i]);
				}
	        }
			System.out.println("Started off with "+chrs.length+" markers in the dataset");
			System.out.println("   "+countFromList+" of "+hash.size()+" markers on the list were removed ("+ext.formDeci((countFromList-countOverlap)/(double)chrs.length*100, 2, true)+"% of total)");
			System.out.println("   "+countInRegions+" were found within the list of regions and removed ("+ext.formDeci((countInRegions-countOverlap)/(double)chrs.length*100, 2, true)+"% of total)");
			System.out.println("   "+countOverlap+" overlap in filtering criteria ("+ext.formDeci(countOverlap/(double)chrs.length*100, 2, true)+"% of total)");
			System.out.println("Leaving behind "+(chrs.length-countFromList-countInRegions+countOverlap)+" in final marker list ("+ext.formDeci((chrs.length-countFromList-countInRegions+countOverlap)/(double)chrs.length*100, 2, true)+"% of total)");
	        writer.close();
        } catch (Exception e) {
	        System.err.println("Error writing to "+outfile);
	        e.printStackTrace();
        }		
	}

	public static void main(String[] args) {
		int numArgs = args.length;
//		String filename = Project.DEFAULT_PROJECT;
		String filename = "C:/workspace/Genvisis/projects/OSv2.properties";
		int numThreads = 6;
		int program = PENN_CNV;
		String filterRegions = "";
		String filterList = "";
		String markers = null;
		Project proj;

		String usage = "\n"+
		"filesys.AnalysisFormats requires 0-1 arguments\n"+
		"   (1) project file (i.e. proj="+filename+" (default))\n"+
		"   (2) number of threads to use (i.e. threads="+numThreads+" (default))\n"+
		"   (3) filter markers out within specified regions (i.e. filterRegions=problematicRegions.dat (not the default))\n"+
		"   (4) filter markers out from list (i.e. filterList=drops.dat (not the default))\n"+
		"   (5) input/output file of final list of markers to use (all markers if null) (i.e. markers="+markers+" (default))\n"+
		"   (6) program option (i.e. program="+program+" (default))\n";
		for (int i = 0; i<PROGRAM_OPTIONS.length; i++) {
			usage += "           "+(i+1)+" = "+PROGRAM_OPTIONS[i]+"\n";
		}

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("threads=")) {
				numThreads = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("program=")) {
				program = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("filterRegions=")) {
				filterRegions = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("filterList=")) {
				filterList = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("markers=")) {
				markers = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}

//		filterRegions = "data/problematicRegions.dat";
//		filterList = "data/drops.dat";
//		markers = "finalMarkerList.dat";
		try {
			proj = new Project(filename, false);
			if (!filterRegions.equals("") || !filterList.equals("")) {
				filter(proj, filterRegions, filterList, markers);
			} else {
				launch(proj, program, markers, numThreads);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
