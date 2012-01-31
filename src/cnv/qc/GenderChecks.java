// -Xms1024M -Xmx1024M
package cnv.qc;

import java.io.*;
import java.util.*;

import cnv.filesys.*;
import cnv.var.SampleData;
import common.*;
import stats.*;

public class GenderChecks {
	public static final String[] SAMPLE_FIELDS = {"DNA", "IID", "CLASS=Gender"};
	public static final String[] SNP_FIELDS = {"Sample", "X", "Y", "X Raw", "Y Raw", "Theta", "R", "B Allele Freq", "Log R Ratio", "AlleleCount"};
	public static final String WIN_PLOTS_DIR = "C:\\Documents and Settings\\npankrat\\My Documents\\jProjects\\CNVisualize\\plots\\scratch\\";
	public static final String WIN_LOOKUP_DIR = "C:\\Documents and Settings\\npankrat\\My Documents\\jProjects\\CNVisualize\\lookup\\";
	// public static final String LINUX_PLOTS_DIR = "/home/npankrat/gwas/dists/";
	// public static final String LINUX_LOOKUP_DIR = "/home/npankrat/gwas/dists/lookup/";
	public static final String LINUX_PLOTS_DIR = "/home/genanal/genetics_anals/analysis/plots/";
	public static final String LINUX_LOOKUP_DIR = "/home/genanal/genetics_anals/analysis/plots/lookup/";
	public static final String RESULTS_DIR = "results/genderChecks/";

	public static void sexCheck(Project proj) {
		PrintWriter writer;
		boolean[][] sexlinked;
		int numXs, numYs, numX_15_85;
		double lrrsX, lrrsY;
		long time;
		MarkerSet markerSet;
		byte[] chrs;
		String[] samples;
		Sample samp;
		float[] lrrs, bafs;
		SampleData sampleData;
		String trav;

		time = new Date().getTime();
		markerSet = proj.getMarkerSet();
		chrs = markerSet.getChrs();

		sexlinked = new boolean[chrs.length][2];
		for (int i = 0; i<chrs.length; i++) {
			if (chrs[i]==23) {
				sexlinked[i][0] = true;
			}
			if (chrs[i]==24) {
				sexlinked[i][1] = true;
			}
		}
		
		sampleData = proj.getSampleData(false);

		System.out.println("Took "+ext.getTimeElapsed(time)+" to hash samples");
		time = new Date().getTime();

		samples = proj.getSamples();

		try {
			writer = new PrintWriter(new FileWriter(proj.getProjectDir()+"sexCheck.xln"));
			writer.println("Sample\tID\tGender\tMean X LRR\tNum X markers\tNum X 15-85\t% 15-85\tMean Y LRR\tNum Y markers");
			for (int i = 0; i<samples.length; i++) {
				numXs = numYs = numX_15_85 = 0;
				lrrsX = lrrsY = 0;

				samp = proj.getSample(samples[i]);
				if (samp==null) {
					System.err.println("Error - could not load sample: "+samples[i]);
					System.exit(1);
				}
				if (markerSet.getFingerprint()!=samp.getFingerprint()) {
					System.err.println("Error - mismatched MarkerSet fingerprints for sample "+samples[i]);
					System.exit(1);
				}
				lrrs = samp.getLRRs();
				bafs = samp.getBAFs();

				for (int j = 0; j<lrrs.length; j++) {
					if (sexlinked[j][0]&&!Double.isNaN(lrrs[j])) {
						lrrsX += lrrs[j];
						numXs++;
						if (bafs[j] > 0.15 && bafs[j] < 0.85) {
							numX_15_85++;
						}
					}
					if (sexlinked[j][1]&&!Double.isNaN(lrrs[j])) {
						lrrsY += lrrs[j];
						numYs++;
					}
				}

				if (i%100 == 0) {
					System.out.println("parsed "+samples[i]+" ("+(i+1)+" of "+samples.length+")");
				}
				trav = sampleData.lookup(samples[i]);
				if (trav == null) {
					System.err.println("Error - no data for sample '"+samples[i]+"'");
					writer.print(samples[i]+"\t"+".\t.");
				} else {
					writer.print(samples[i]+"\t"+trav);
				}
				writer.println("\t"+ext.formDeci(lrrsX/numXs, 7)+"\t"+numXs+"\t"+numX_15_85+"\t"+ext.formDeci((double)numX_15_85/(double)numXs, 4)+"\t"+ext.formDeci(lrrsY/numYs, 7)+"\t"+numYs);
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to sexCheck.xln");
			e.printStackTrace();
		}
		System.out.println("Took "+ext.getTimeElapsed(time)+" to parse "+samples.length+" samples");
	}

	public static void markerByMarker(Project proj, String filename, String plots_dir) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		int count;
		String[] snpNames, samples;
		int[] fieldIndices;
		int[] genders;
		Vector<double[]> xys, raw_xys, thetars, lrrbafs;
		Vector<String> intensityDeps;
		LogisticRegression lr;
		String trav, output;
		boolean satisfied = false;
		int sampIndex;

		trav = ext.rootOf(filename);
		System.err.println("Lookup: "+trav);

		while (!satisfied&&trav.length()>0) {
			try {
				Integer.parseInt(trav.charAt(trav.length()-1)+"");
				satisfied = true;
			} catch (NumberFormatException nfe) {
				trav = trav.substring(0, trav.length()-1);
			}
		}
		if (trav.length()==0) {
			System.err.println("Error - expecting a filename that starts with a number (not '"+ext.rootOf(filename)+"')");
			System.exit(1);
		}

		snpNames = Array.toStringArray(HashVec.loadFileToVec(filename, false, false, true, true));
		// does this do what it's supposed to? The proj.getFilename used to be a ManageFiles.getFilename. Is this the same? 
		samples =  Array.toStringArray(HashVec.loadFileToVec(plots_dir+trav+"/"+proj.getFilename(snpNames[0]), true, true, false));
//		samples = null;

		genders = Array.intArray(samples.length, -1);
		try {
			reader = new BufferedReader(new FileReader(proj.getFilename(Project.SAMPLE_DATA_FILENAME)));
			fieldIndices = ext.indexFactors(SAMPLE_FIELDS, reader.readLine().trim().split("\t", -1), false, true);
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				sampIndex = ext.indexOfStr(line[fieldIndices[0]], samples);
				if (sampIndex!=-1) {
					genders[sampIndex] = Integer.parseInt(line[fieldIndices[2]]);
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+proj.getFilename(Project.SAMPLE_DATA_FILENAME)+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+proj.getFilename(Project.SAMPLE_DATA_FILENAME)+"\"");
			System.exit(2);
		}

		if (!new File(RESULTS_DIR).exists()) {
			new File(RESULTS_DIR).mkdirs();
		}
		try {
			writer = new PrintWriter(new FileWriter(RESULTS_DIR+ext.rootOf(filename)+"_genderChecks.xln"));
			writer.println("SNP\tX abs(T)\tY abs(T)\tRaw X abs(T)\tRaw Y abs(T)\tTheta abs(T)\tR abs(T)\tBAF abs(T)\tLRR abs(T)\tX p\tY p\tXY r2\tRaw X p\tRaw Y p\tRaw XY r2\tTheta p\tR p\tTheta/R r2\tBAF p\tLRR p\tBAF/LRR r2");
			for (int i = 0; i<snpNames.length; i++) {
				System.err.println(snpNames[i]);
				intensityDeps = new Vector<String>();
				xys = new Vector<double[]>();
				raw_xys = new Vector<double[]>();
				thetars = new Vector<double[]>();
				lrrbafs = new Vector<double[]>();

				try {
					reader = new BufferedReader(new FileReader(plots_dir+trav+"/"+snpNames[i]));
					fieldIndices = ext.indexFactors(SNP_FIELDS, reader.readLine().trim().split("\t", -1), false, true);
					count = 0;
					while (reader.ready()) {
						line = reader.readLine().trim().split("[\\s]+");
						if (!line[fieldIndices[0]].equals(samples[count])) {
							System.err.println("Error - out of sync on marker "+snpNames[i]);
							System.exit(1);
						}
						if (genders[count]!=-1) {
							if (!line[fieldIndices[8]].equals("NaN")) {
								intensityDeps.add(genders[count]+"");
								xys.add(new double[] {Double.parseDouble(line[fieldIndices[1]]), Double.parseDouble(line[fieldIndices[2]])});
								raw_xys.add(new double[] {Double.parseDouble(line[fieldIndices[3]]), Double.parseDouble(line[fieldIndices[4]])});
								thetars.add(new double[] {Double.parseDouble(line[fieldIndices[5]]), Double.parseDouble(line[fieldIndices[6]])});
								lrrbafs.add(new double[] {Double.parseDouble(line[fieldIndices[7]]), Double.parseDouble(line[fieldIndices[8]])});
							}
						}
						count++;
					}
					reader.close();
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \""+plots_dir+trav+"/"+snpNames[i]+"\" not found in current directory");
					System.exit(1);
				} catch (IOException ioe) {
					System.err.println("Error reading file \""+plots_dir+trav+"/"+snpNames[i]+"\"");
					System.exit(2);
				}

				line = new String[30];
				output = snpNames[i];

				if (intensityDeps.size()==0) {
					System.err.println("Warning - no data for marker "+snpNames[i]);
					output += "\t.\t.\t.\t.\t.\t.\t.\t.";
					line[12] = "NoData";
					line[13] = "NoData";
					line[14] = "NoData";
				} else {
					output += "\t"+Math.abs(new Ttest(Array.toIntArray(Array.toStringArray(intensityDeps)), Matrix.extractColumn(Matrix.toDoubleArrays(xys), 0)).getT());
					output += "\t"+Math.abs(new Ttest(Array.toIntArray(Array.toStringArray(intensityDeps)), Matrix.extractColumn(Matrix.toDoubleArrays(xys), 1)).getT());
					output += "\t"+Math.abs(new Ttest(Array.toIntArray(Array.toStringArray(intensityDeps)), Matrix.extractColumn(Matrix.toDoubleArrays(raw_xys), 0)).getT());
					output += "\t"+Math.abs(new Ttest(Array.toIntArray(Array.toStringArray(intensityDeps)), Matrix.extractColumn(Matrix.toDoubleArrays(raw_xys), 1)).getT());
					output += "\t"+Math.abs(new Ttest(Array.toIntArray(Array.toStringArray(intensityDeps)), Matrix.extractColumn(Matrix.toDoubleArrays(thetars), 0)).getT());
					output += "\t"+Math.abs(new Ttest(Array.toIntArray(Array.toStringArray(intensityDeps)), Matrix.extractColumn(Matrix.toDoubleArrays(thetars), 1)).getT());
					output += "\t"+Math.abs(new Ttest(Array.toIntArray(Array.toStringArray(intensityDeps)), Matrix.extractColumn(Matrix.toDoubleArrays(lrrbafs), 0)).getT());
					output += "\t"+Math.abs(new Ttest(Array.toIntArray(Array.toStringArray(intensityDeps)), Matrix.extractColumn(Matrix.toDoubleArrays(lrrbafs), 1)).getT());
				}

				lr = null;
				try {
					lr = new LogisticRegression(intensityDeps, xys);
					output += "\t"+lr.getWalds()[1]+"\t"+lr.getWalds()[2]+"\t"+(lr.getRsquare()<0?".":lr.getRsquare());
					lr.dumpData(snpNames[i]+"_xys.xln");
				} catch (Exception e) {
					output += "\t.\t.\t.";
				}

				try {
					lr = new LogisticRegression(intensityDeps, raw_xys);
					output += "\t"+lr.getWalds()[1]+"\t"+lr.getWalds()[2]+"\t"+(lr.getRsquare()<0?".":lr.getRsquare());
					lr.dumpData(snpNames[i]+"_raw_xys.xln");
				} catch (Exception e) {
					output += "\t.\t.\t.";
				}

				try {
					lr = new LogisticRegression(intensityDeps, thetars);
					output += "\t"+lr.getWalds()[1]+"\t"+lr.getWalds()[2]+"\t"+(lr.getRsquare()<0?".":lr.getRsquare());
					lr.dumpData(snpNames[i]+"_thetar.xln");
				} catch (Exception e) {
					output += "\t.\t.\t.";
				}

				try {
					lr = new LogisticRegression(intensityDeps, lrrbafs);
					output += "\t"+lr.getWalds()[1]+"\t"+lr.getWalds()[2]+"\t"+(lr.getRsquare()<0?".":lr.getRsquare());
					// lr.dumpData(snpNames[i]+"_lrr_baf.xln");
				} catch (Exception e) {
					output += "\t.\t.\t.";
				}

				writer.println(output);
				writer.flush();
			}

			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing results");
			e.printStackTrace();
		}

	}

	public static void batch(String lookup_dir, int numBatches) {
		PrintWriter writer;

		String[] filenames = new File(lookup_dir).list(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				return filename.endsWith(".txt");
			}
		});

		if (filenames==null) {
			System.err.println("Error - directory not found: "+lookup_dir);

		}
		System.out.println("Found "+filenames.length+" files to run");

		// String commands = "vis scatter.GenderChecks file="+lookup_dir+"[%0]";
		String commands = "java -cp /home/genanal/vis.jar scatter.GenderChecks file="+lookup_dir+"[%0]";

		Files.batchIt("gender", null, numBatches, commands, filenames);
		try {
			writer = new PrintWriter(new FileWriter("masterBatch"));
			for (int i = 0; i<numBatches; i++) {
				writer.println("nohup ./gender."+(i+1)+" > "+(i+1)+".out &");
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error witing masterBatch");
			e.printStackTrace();
		}

	}

	public static void parse(Project proj) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String trav;
		Hashtable<String,String> hash = new Hashtable<String,String>();

		File[] filenames = new File(RESULTS_DIR).listFiles(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				return filename.endsWith("_genderChecks.xln");
			}
		});

		if (filenames==null) {
			System.err.println("Error - directory not found: "+RESULTS_DIR);

		} else {
			System.out.println("Found results for "+filenames.length+" lookup files");
		}

		hash = HashVec.loadFileToHashString(proj.getFilename(Project.MARKERSET_FILENAME), 0, new int[] {1, 2}, "\t", true);

		try {
			writer = new PrintWriter(new FileWriter("GenderChecks.xln"));
			for (int i = 0; i<filenames.length; i++) {
				System.out.println((i+1)+" of "+filenames.length);
				try {
					reader = new BufferedReader(new FileReader(filenames[i]));
					if (i==0) {
						line = reader.readLine().trim().split("\t");
						writer.println(line[0]+"\tChr\tPosition\t"+Array.toStr(Array.subArray(line, 1, line.length)));
					} else {
						reader.readLine();
					}
					while (reader.ready()) {
						line = reader.readLine().trim().split("[\\s]+");
						trav = hash.containsKey(line[0])?hash.get(line[0]):".\t.";
						writer.println(line[0]+"\t"+trav+"\t"+Array.toStr(Array.subArray(line, 1, line.length)));
						hash.remove(line[0]);
					}
					reader.close();
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \""+filenames[i].getName()+"\" not found in current directory");
					System.exit(1);
				} catch (IOException ioe) {
					System.err.println("Error reading file \""+filenames[i].getName()+"\"");
					System.exit(2);
				}
			}
			line = HashVec.getKeys(hash);
			for (int i = 0; i<line.length; i++) {
				writer.println(line[i]+"\t"+hash.get(line[i]));
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing results");
			e.printStackTrace();
		}
	}

	public static void dropMarkers(String allMarkers, String markersToDrop) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		Vector<String> v = new Vector<String>();

		v = HashVec.loadFileToVec(markersToDrop, false, true, true, true);

		try {
			reader = new BufferedReader(new FileReader(allMarkers));
			writer = new PrintWriter(new FileWriter(ext.rootOf(allMarkers)+"_dropped.out"));
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				if (!v.contains(line[0])) {
					writer.println(Array.toStr(line));
				}
			}
			writer.close();
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+allMarkers+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+allMarkers+"\"");
			System.exit(2);
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		// String filename = "C:\\Documents and Settings\\npankrat\\My Documents\\jProjects\\CNVisualize\\1.txt";
		// String filename = "";
		boolean check = false;
		int batch = 0;
		boolean parse = false;
		String plotsDir = System.getProperty("os.name").startsWith("Windows")?WIN_PLOTS_DIR:LINUX_PLOTS_DIR;
		String lookupDir = System.getProperty("os.name").startsWith("Windows")?WIN_LOOKUP_DIR:LINUX_LOOKUP_DIR;
		String markersToDrop = "data/drops.dat";
		String allMarkers = "data/markerListWithIndices.dat";
		boolean drop = false;
		Project proj;
		String filename = Project.DEFAULT_PROJECT;

		String usage = "\\n"+
		"qc.GenderChecks requires 0-1 arguments\n"+
		"   (1) project file (i.e. proj="+filename+" (default))\n"+
		"   (2) check sex of indiviudals (i.e. -check (not the default))\n"+
		"   (3) marker lookup file (i.e. file=1a.txt (not the default))\n"+
		"   (4) batch all (i.e. batch=4 (not the default))\n"+
		"   (5) parse all results (i.e. -parse (not the default))\n"+
		"   (6) drop markers (i.e. -drop (not the default))\n"+
		"   (7) file with all markers (i.e. all="+allMarkers+" (default file))\n"+
		"   (8) list of bad markers (i.e. drop="+markersToDrop+" (default file))\n"+
		"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("-check")) {
				check = true;
				numArgs--;
			} else if (args[i].startsWith("batch=")) {
				batch = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("-parse")) {
				parse = true ;
				numArgs--;
			}
		}

		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}

		check = true;
		try {
			proj = new Project(filename, false);

			if (check) {
				sexCheck(proj);
			} else if (batch>0) {
				batch(lookupDir, batch);
			} else if (parse) {
				parse(proj);
			} else if (drop) {
				dropMarkers(allMarkers, markersToDrop);
			} else {
				markerByMarker(proj, null, plotsDir);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}