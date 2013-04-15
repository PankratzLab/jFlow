// -Xms1024M -Xmx1024M
package cnv.qc;

import java.io.*;
import java.util.*;

import cnv.filesys.*;
import cnv.var.SampleData;
import common.*;
import stats.*;

public class SexChecks {
	public static final String[] SEX_HEADER = {"Sample", "FID", "IID", "Sex", "Estimated Sex", "Mean X LRR", "Num X markers", "Num X 10-90", "% 10-90", "Mean Y LRR", "Num Y markers"};
	public static final String[] SAMPLE_FIELDS = {"DNA", "IID", "CLASS=Gender"};
	public static final String[] SNP_FIELDS = {"Sample", "X", "Y", "X Raw", "Y Raw", "Theta", "R", "B Allele Freq", "Log R Ratio", "AlleleCount"};
	public static final String RESULTS_DIR = "results/genderChecks/";
	public static final float NUM_SD_FOR_MALE_OUTLIERS = 3.0f;
	public static final float NUM_SD_FOR_FEMALE_OUTLIERS = 3.0f;

	static MarkerSet markerSet;
	static SampleData sampleData;
	static Sample samp;
	static String[] samples;
	static int[] numXs;
	static int[] numYs;
	static int[] numX_10_90;
	static double[] lrrsX;
	static double[] lrrsY;

	public static void sexCheck(Project proj) {
		long time;
		time = new Date().getTime();

		markerSet = proj.getMarkerSet();
		sampleData = proj.getSampleData(2, false);

		System.out.println("Took "+ext.getTimeElapsed(time)+" to hash samples");
		time = new Date().getTime();

		samples = proj.getSamples();
		lrrCounts(samples, proj, sexlinked(markerSet.getChrs()));
		writeToFile(proj, estimateSex(samples));

		System.out.println("Took "+ext.getTimeElapsed(time)+" to parse "+samples.length+" samples");
	}

	public static void markerByMarker(Project proj) {
		PrintWriter writer;
		String[] samples;
		Vector<double[]> xys, baflrrs; 
		float[] xs, ys, lrrs, bafs;
		Vector<String> intensityDeps;
		LogisticRegression lr;
		String output;
		MarkerData[] markerData;
        String[] files;
        SampleData sampleData;
        int[] sexes;
        
        sampleData = proj.getSampleData(2, false);
        samples = proj.getSamples();
        sexes = new int[samples.length];
        for (int i = 0; i < samples.length; i++) {
        	sexes[i] = sampleData.getSexForIndividual(samples[i]);
		}
 		
		files = Files.list(proj.getDir(Project.PLOT_DIRECTORY), ".scat", false);

		try {
			writer = new PrintWriter(new FileWriter(proj.getDir(Project.RESULTS_DIRECTORY)+"markerGenderChecks.xln"));
			writer.println("SNP\tX abs(T)\tY abs(T)\tBAF abs(T)\tLRR abs(T)\tX p\tY p\tXY r2\tBAF p\tLRR p\tBAF/LRR r2");
			for (int i=0; i<files.length; i++) {
				markerData = MarkerDataCollection.load(proj.getDir(Project.PLOT_DIRECTORY)+files[i], proj.getJarStatus()).getCollection();
				for (int j = 0; j < markerData.length; j++) {
					output = markerData[j].getMarkerName();
					
					xs = markerData[j].getXs();
					ys = markerData[j].getYs();
					bafs = markerData[j].getBAFs();
					lrrs = markerData[j].getLRRs();
					
					intensityDeps = new Vector<String>();
					xys = new Vector<double[]>();
					baflrrs = new Vector<double[]>();
					for (int s = 0; s < samples.length; s++) {
						if (ext.isValidDouble(lrrs[s]+"")) {
							intensityDeps.add(sexes[s]+"");
							xys.add(new double[] {xs[s], ys[s]});
							baflrrs.add(new double[] {bafs[s], lrrs[s]});
						}
					}
					
					
					if (intensityDeps.size()==0) {
						System.err.println("Warning - no data for marker "+markerData[j].getMarkerName());
						output += "\t.\t.\t.\t.\t.\t.\t.\t.";
					} else {
						output += "\t"+Math.abs(new Ttest(Array.toIntArray(Array.toStringArray(intensityDeps)), Matrix.extractColumn(Matrix.toDoubleArrays(xys), 0)).getPvalue());
						output += "\t"+Math.abs(new Ttest(Array.toIntArray(Array.toStringArray(intensityDeps)), Matrix.extractColumn(Matrix.toDoubleArrays(xys), 1)).getPvalue());
						output += "\t"+Math.abs(new Ttest(Array.toIntArray(Array.toStringArray(intensityDeps)), Matrix.extractColumn(Matrix.toDoubleArrays(baflrrs), 0)).getPvalue());
						output += "\t"+Math.abs(new Ttest(Array.toIntArray(Array.toStringArray(intensityDeps)), Matrix.extractColumn(Matrix.toDoubleArrays(baflrrs), 1)).getPvalue());
					}
	
					lr = null;
					try {
						lr = new LogisticRegression(intensityDeps, xys);
						output += "\t"+lr.getSigs()[1]+"\t"+lr.getSigs()[2]+"\t"+(lr.getRsquare()<0?".":lr.getRsquare());
					} catch (Exception e) {
						output += "\t.\t.\t.";
					}
					try {
						lr = new LogisticRegression(intensityDeps, baflrrs);
						output += "\t"+lr.getSigs()[1]+"\t"+lr.getSigs()[2]+"\t"+(lr.getRsquare()<0?".":lr.getRsquare());
					} catch (Exception e) {
						output += "\t.\t.\t.";
					}
	
					writer.println(output);
					writer.flush();
				}
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing results");
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

	private static boolean[][] sexlinked(byte[] chrs) {
		boolean[][] sexlinked;
		sexlinked = new boolean[chrs.length][2];
		for (int i = 0; i<chrs.length; i++) {
			if (chrs[i]==23) {
				sexlinked[i][0] = true;
			}
			if (chrs[i]==24) {
				sexlinked[i][1] = true;
			}
		}
		return sexlinked;
	}

	private static void lrrCounts(String[] samples, Project proj, boolean[][] sexlinked) {
		float[] lrrs, bafs;
		numXs = new int[samples.length];
		numYs = new int[samples.length];
		numX_10_90 = new int[samples.length];
		lrrsX = new double[samples.length];
		lrrsY = new double[samples.length];

		for (int i = 0; i<samples.length; i++) {
			numXs[i] = numYs[i] = numX_10_90[i] = 0;
			lrrsX[i] = lrrsY[i] = 0;

//			samp = proj.getSample(samples[i]);
			samp = proj.getPartialSampleFromRandomAccessFile(samples[i]);
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
					lrrsX[i] += lrrs[j];
					numXs[i]++;
//					if (bafs[j] > 0.15 && bafs[j] < 0.85) {
					if (bafs[j] > 0.10 && bafs[j] < 0.9) {
						numX_10_90[i]++;
					}
				}
				if (sexlinked[j][1]&&!Double.isNaN(lrrs[j])) {
					lrrsY[i] += lrrs[j];
					numYs[i]++;
				}
			}

			if (i%100 == 0) {
				System.out.println("parsed "+samples[i]+" ("+(i+1)+" of "+samples.length+")");
			}
		}
	}


	public static byte[] estimateSex(String[] samples) {
//		the result will be used for color code for the points;
//		mean for M and F;
//		standard deviation for M and F;

		byte[] result = new byte[samples.length];
		double putativeMaleMeanY, putativeFemaleMeanY;
		double[] values;
		double maleMeanX, maleStdDevX, femaleMeanX, femaleStdDevX;
		int numMales, numFemales, putativeSex;
		
		
		values = new double[samples.length];
		for (int i = 0; i < values.length; i++) {
			values[i] = lrrsY[i]/numYs[i];
		}
		if (Array.isBimodal(values, 0.01, 100)) {
			
		}

		putativeMaleMeanY=0;
		putativeFemaleMeanY=0;
		numMales = numFemales = 0;
		for (int i=0; i<samples.length; i++) {
			putativeSex = sampleData.getSexForIndividual(samples[i]);
			switch (putativeSex) {
			case 1:
				putativeMaleMeanY += (lrrsY[i]/numYs[i]);
				numMales++;
				break;
			case 2:
				putativeFemaleMeanY += (lrrsY[i]/numYs[i]);
				numFemales++;
				break;
			default:
					
			}
		}
		putativeMaleMeanY=putativeMaleMeanY/(double)numMales;
		putativeFemaleMeanY=putativeFemaleMeanY/(double)numFemales;

		DoubleVector males, females;
		males = new DoubleVector();
		females = new DoubleVector();
		for (int i=0; i<samples.length; i++) {
			if (Math.abs(lrrsY[i]/numYs[i] - putativeMaleMeanY) < Math.abs(lrrsY[i]/numYs[i] - putativeFemaleMeanY)) {
				males.add((lrrsX[i]/numXs[i]));
			} else {
				females.add((lrrsX[i]/numXs[i]));
			}
		}
		
		values = males.toArray();
		maleMeanX = Array.mean(values);
		maleStdDevX = Array.stdev(values);
		values = females.toArray();
		femaleMeanX = Array.mean(values);
		femaleStdDevX = Array.stdev(values);

		for (int i=0; i<samples.length; i++) {
			if (Math.abs(lrrsY[i]/numYs[i] - putativeMaleMeanY) < Math.abs(lrrsY[i]/numYs[i] - putativeFemaleMeanY)) {
				if (lrrsX[i]/numXs[i] > (maleMeanX + NUM_SD_FOR_MALE_OUTLIERS*maleStdDevX)) {
					if (samp.hasBimodalBAF((byte)23, 0, Integer.MAX_VALUE)) {
						result[i] = 4; // mosaic Klinefelter
					} else {
						result[i] = 3; // full Klinefelter
					}
				} else {
					result[i] = 1; // normal male
				}
			} else {
				if (lrrsX[i]/numXs[i] > (femaleMeanX + NUM_SD_FOR_FEMALE_OUTLIERS*femaleStdDevX) && samp.hasBimodalBAF((byte)23, 0, Integer.MAX_VALUE)) {
					result[i] = 5; // Triple X syndrome
				} else if (lrrsX[i]/numXs[i] < (femaleMeanX - NUM_SD_FOR_FEMALE_OUTLIERS*femaleStdDevX) && (double)numX_10_90[i]/(double)numXs[i] < 0.15) {
					result[i] = 6; // Turner
				} else if (lrrsX[i]/numXs[i] < (femaleMeanX - NUM_SD_FOR_FEMALE_OUTLIERS*femaleStdDevX) && samp.hasBimodalBAF((byte)23, 0, Integer.MAX_VALUE)) {
					result[i] = 7; // mosaic Turner
				} else {
					result[i] = 2; // normal female
				}
			}
		}

		return result;
	}

	public double[] medianXYBySex(double[][] data, byte[] sexes) {
//		return new double[] {Arrays.sort(data[][0]).[data.length/2], Arrays.sort(data[][1]).[data.length/2]}
		ArrayList<Double> maleX, maleY, femaleX, femaleY;
		maleX	= new ArrayList<Double>();
		maleY	= new ArrayList<Double>();
		femaleX	= new ArrayList<Double>(); 
		femaleY	= new ArrayList<Double>();
		
		for (int i=0; i<data.length; i++) {
			if (sexes[i]==1) {
				maleX.add(data[i][0]);
				maleY.add(data[i][1]);
			} else {
				femaleX.add(data[i][0]);
				femaleY.add(data[i][1]);
			}
		}
		Collections.sort(maleX);
		Collections.sort(maleY);
		Collections.sort(femaleX);
		Collections.sort(femaleY);
		return new double[] {maleX.get(maleX.size()/2), maleY.get(maleY.size()/2), femaleX.get(femaleX.size()/2), femaleY.get(femaleY.size()/2)};
	}

	public static void writeToFile (Project proj, byte[] estimatedSex) {
		PrintWriter writer;
		String famIndPair;

		try {
			writer = new PrintWriter(new FileWriter(proj.getProjectDir()+"sexCheck.xln"));
			writer.println(Array.toStr(SEX_HEADER));
			for (int i = 0; i<samples.length; i++) {
				famIndPair = sampleData.lookup(samples[i])[1];
				if (famIndPair == null) {
					System.err.println("Error - no data for sample '"+samples[i]+"'");
					writer.print(samples[i]+"\t"+".\t.\t-9\t-9");
				} else {
					writer.print(samples[i]+"\t"+famIndPair+"\t"+sampleData.getSexForIndividual(samples[i])+"\t"+estimatedSex[i]);
				}
				writer.println("\t"+ext.formDeci(lrrsX[i]/numXs[i], 7)+"\t"+numXs[i]+"\t"+numX_10_90[i]+"\t"+ext.formDeci((double)numX_10_90[i]/(double)numXs[i], 4)+"\t"+ext.formDeci(lrrsY[i]/numYs[i], 7)+"\t"+numYs[i]);
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to sexCheck.xln");
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		boolean check = false;
		boolean parse = false;
		String markersToDrop = "data/drops.dat";
		String allMarkers = "data/markerListWithIndices.dat";
		boolean drop = false;
		Project proj;
		String filename = Project.DEFAULT_PROJECT;

		String usage = "\\n"+
		"qc.GenderChecks requires 0-1 arguments\n"+
		"   (1) project file (i.e. proj="+filename+" (default))\n"+
		"   (2) check sex of indiviudals (i.e. -check (not the default))\n"+
		" OR\n"+
		"   (2) parse all results (i.e. -parse (not the default))\n"+
		"   (3) drop markers (i.e. -drop (not the default))\n"+
		"   (4) file with all markers (i.e. all="+allMarkers+" (default file))\n"+
		"   (5) list of bad markers (i.e. drop="+markersToDrop+" (default file))\n"+
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
			} else if (args[i].startsWith("-parse")) {
				parse = true ;
				numArgs--;
			}
		}

		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}

//		check = true;
		try {
			proj = new Project(filename, false);

			if (check) {
				sexCheck(proj);
			} else if (parse) {
				parse(proj);
			} else if (drop) {
				dropMarkers(allMarkers, markersToDrop);
			} else {
				markerByMarker(proj);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}