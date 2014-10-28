// -Xms1024M -Xmx1024M
package cnv.qc;

import java.io.*;
import java.util.*;

import cnv.filesys.*;
import cnv.manage.MarkerDataLoader;
import cnv.var.SampleData;
import common.*;
import stats.*;

public class SexChecks {
	public static final String EST_SEX_HEADER = "Estimated Sex;1=Male;2=Female;3=Klinefelter;4=Mosaic Klinefelter;5=Triple X;6=Turner;7=Mosaic Turner";
	public static final String[] SEX_HEADER = {"Sample", "FID", "IID", "Sex", EST_SEX_HEADER, "Mean X LRR", "Num X markers", "Num X 10-90", "% 10-90", "Mean Y LRR", "Num Y markers"};
	public static final String[] SAMPLE_FIELDS = {"DNA", "IID", "CLASS=Gender"};
	public static final String[] SNP_FIELDS = {"Sample", "X", "Y", "X Raw", "Y Raw", "Theta", "R", "B Allele Freq", "Log R Ratio", "AlleleCount"};
	public static final String RESULTS_DIR = "genderChecks/";
	public static final float NUM_SD_FOR_MALE_OUTLIERS = 5.0f;
	public static final float NUM_SD_FOR_FEMALE_OUTLIERS = 5.0f;

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
		Logger log;
		
		log = proj.getLog();
		time = new Date().getTime();

		markerSet = proj.getMarkerSet();
		sampleData = proj.getSampleData(2, false);
		if (sampleData.failedToLoad()) {
			log.reportError("Error - without a sample data file, sexChecks will fail");
			return;
		}

		log.report("Took "+ext.getTimeElapsed(time)+" to hash samples");
		time = new Date().getTime();

		samples = proj.getSamples();
		lrrCounts(samples, proj, sexlinked(markerSet.getChrs()));
		byte[] estSex = estimateSex(samples);
		writeToFile(proj, estSex);
		
		String[] classes = sampleData.getClasses();
		int sexInd = -1;
		for (int i = 0; i < classes.length; i++) {
			if (SexChecks.EST_SEX_HEADER.equals(classes[i])) {
				sexInd = i;
				break;
			}
		}
		if (sexInd == -1) {
			Hashtable<String, String> linkData = new Hashtable<String, String>();
			for (int i = 0; i < samples.length; i++) {
				linkData.put(samples[i], "" + estSex[i]);
			}
			if (!sampleData.addData(linkData, "DNA", new String[] {"CLASS=" + EST_SEX_HEADER}, ".", "", log)) {
				log.reportError("Error - failed to write Estimated Sex to sample data file");
			}
		} else {
			log.report("Warning - sample data already contains estimated sex; will not process data into sample data file."); 
		}
		
		
		log.report("Took "+ext.getTimeElapsed(time)+" to parse "+samples.length+" samples");
	}

	public static void markerByMarker(Project proj) {
		PrintWriter writer;
		String[] samples;
		Vector<double[]> xys, baflrrs; 
		float[] xs, ys, lrrs, bafs;
		Vector<String> intensityDeps;
		LogisticRegression lr;
		String output;
        SampleData sampleData;
        int[] sexes;
    	MarkerDataLoader markerDataLoader;
		MarkerData markerData;
		String[] markerNames;
		long time;
		Logger log;
        
		log = proj.getLog();
        sampleData = proj.getSampleData(2, false);
        samples = proj.getSamples();
        sexes = new int[samples.length];
        for (int i = 0; i < samples.length; i++) {
        	sexes[i] = sampleData.getSexForIndividual(samples[i]);
		}
 		

		try {
			writer = new PrintWriter(new FileWriter(proj.getDir(Project.RESULTS_DIRECTORY)+"markerGenderChecks.xln"));
			writer.println("SNP\tX abs(T)\tY abs(T)\tBAF abs(T)\tLRR abs(T)\tX p\tY p\tXY r2\tBAF p\tLRR p\tBAF/LRR r2");
			
	        time = new Date().getTime();
	        markerNames = proj.getMarkerNames();
			markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markerNames);
	        for (int i = 0; i < markerNames.length; i++) {
	        	markerData = markerDataLoader.requestMarkerData(i);
	        	if (i % 100 == 0) {
	        		log.report(ext.getTime()+"\tMarker "+i+" of "+markerNames.length);
	        	}

				output = markerData.getMarkerName();
				
				xs = markerData.getXs();
				ys = markerData.getYs();
				bafs = markerData.getBAFs();
				lrrs = markerData.getLRRs();
				
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
					log.reportError("Warning - no data for marker "+markerData.getMarkerName());
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
				markerDataLoader.releaseIndex(i);
			}
	        log.reportError("Finished in " + ext.getTimeElapsed(time));
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing results");
			log.reportException(e);
		}

	}

	public static void parse(Project proj) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String trav;
		Hashtable<String,String> hash = new Hashtable<String,String>();
		Logger log;

		log = proj.getLog();
		
		File[] filenames = new File(proj.getDir(Project.RESULTS_DIRECTORY)+RESULTS_DIR).listFiles(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				return filename.endsWith("_genderChecks.xln");
			}
		});

		if (filenames==null) {
			log.reportError("Error - directory not found: "+proj.getDir(Project.RESULTS_DIRECTORY)+RESULTS_DIR);

		} else {
			log.report("Found results for "+filenames.length+" lookup files");
		}

		hash = HashVec.loadFileToHashString(proj.getFilename(Project.MARKERSET_FILENAME), 0, new int[] {1, 2}, "\t", true);

		try {
			writer = new PrintWriter(new FileWriter("GenderChecks.xln"));
			for (int i = 0; i<filenames.length; i++) {
				log.report((i+1)+" of "+filenames.length);
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
					log.reportError("Error: file \""+filenames[i].getName()+"\" not found in current directory");
					writer.close();
					return;
				} catch (IOException ioe) {
					log.reportError("Error reading file \""+filenames[i].getName()+"\"");
					writer.close();
					return;
				}
			}
			line = HashVec.getKeys(hash);
			for (int i = 0; i<line.length; i++) {
				writer.println(line[i]+"\t"+hash.get(line[i]));
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing results");
			log.reportException(e);
		}
	}

	public static void dropMarkers(String allMarkers, String markersToDrop) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		HashSet<String> hashSet;

		hashSet = HashVec.loadFileToHashSet(markersToDrop, false);

		try {
			reader = new BufferedReader(new FileReader(allMarkers));
			writer = new PrintWriter(new FileWriter(ext.rootOf(allMarkers)+"_dropped.out"));
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				if (!hashSet.contains(line[0])) {
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
		Logger log;
		
		log = proj.getLog();
		
		numXs = new int[samples.length];
		numYs = new int[samples.length];
		numX_10_90 = new int[samples.length];
		lrrsX = new double[samples.length];
		lrrsY = new double[samples.length];

		for (int i = 0; i<samples.length; i++) {
			numXs[i] = numYs[i] = numX_10_90[i] = 0;
			lrrsX[i] = lrrsY[i] = 0;

			samp = proj.getPartialSampleFromRandomAccessFile(samples[i]);
			if (samp == null) {
				log.reportError("Error - could not load sample: " + samples[i]);
				return;
			}
			if (markerSet.getFingerprint() != samp.getFingerprint()) {
				log.reportError("Error - mismatched MarkerSet fingerprints for sample " + samples[i]);
				return;
			}
			lrrs = samp.getLRRs();
			bafs = samp.getBAFs();

			for (int j = 0; j < lrrs.length; j++) {
				if (sexlinked[j][0] && !Double.isNaN(lrrs[j])) {
					lrrsX[i] += lrrs[j];
					numXs[i]++;
					if (bafs[j] > 0.10 && bafs[j] < 0.9) {
						numX_10_90[i]++;
					}
				}
				if (sexlinked[j][1] && !Double.isNaN(lrrs[j])) {
					lrrsY[i] += lrrs[j];
					numYs[i]++;
				}
			}

			if (i % 100 == 0) {
				log.report("parsed "+samples[i]+" ("+(i+1)+" of "+samples.length+")");
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
		// TODO not sure where we left off with this...
		if (Array.isBimodal(values, 0.01, 100)) {
			
		}

		putativeMaleMeanY=0;
		putativeFemaleMeanY=0;
		numMales = numFemales = 0;
		for (int i=0; i<samples.length; i++) {
			putativeSex = sampleData.getSexForIndividual(samples[i]);
			if (putativeSex == -1) {
				if (lrrsY[i] / numYs[i] < -0.75) {
					putativeSex = 2;
				} else {
					putativeSex = 1;
				}
			}
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
		putativeMaleMeanY = putativeMaleMeanY / (double) numMales;
		putativeFemaleMeanY = putativeFemaleMeanY / (double) numFemales;

		DoubleVector males, females;
		males = new DoubleVector();
		females = new DoubleVector();
		for (int i = 0; i < samples.length; i++) {
			if (Math.abs(lrrsY[i] / numYs[i] - putativeMaleMeanY) < Math.abs(lrrsY[i] / numYs[i] - putativeFemaleMeanY)) {
				males.add((lrrsX[i] / numXs[i]));
			} else {
				females.add((lrrsX[i] / numXs[i]));
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

	private static void writeToFile (Project proj, byte[] estimatedSex) {
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
	
	public static void identifyPseudoautosomalBreakpoints(Project proj) {
		PrintWriter writer;
		String[] samples;
		float[] lrrs;
		MarkerData markerData;
        SampleData sampleData;
        int[] sexes;
        byte[] abGenotypes;
        String markerName;
        ClusterFilterCollection clusterFilterCollection;
        float gcThreshold;
        long time;
        DoubleVector[] values; // sex
        MarkerDataLoader markerDataLoader;
        String[] markerList;
        String line, eol;
		MarkerSet markerSet;
        String[] markerNames;
        boolean[] sexChrs;
        byte[] chrs;
        int[][] genotypeCounts;
        boolean[] samplesToExclude;
        
        if (System.getProperty("os.name").startsWith("Windows")) {
        	eol = "\r\n";
		} else {
			eol = "\n";
		}
        
        sampleData = proj.getSampleData(2, false);
        samplesToExclude = proj.getSamplesToExclude();
        samples = proj.getSamples();
        sexes = new int[samples.length];
        for (int i = 0; i < samples.length; i++) {
        	sexes[i] = Math.max(0, sampleData.getSexForIndividual(samples[i]));
		}
        
        markerSet = proj.getMarkerSet();
        markerNames = markerSet.getMarkerNames();
        chrs = markerSet.getChrs();
        sexChrs = new boolean[chrs.length];
        for (int i = 0; i < chrs.length; i++) {
        	sexChrs[i] = chrs[i]>=23;
		}
        markerList = Array.subArray(markerNames, sexChrs);
 		
        clusterFilterCollection = proj.getClusterFilterCollection();
        gcThreshold = Float.parseFloat(proj.getProperty(Project.GC_THRESHOLD));

        try {
			writer = new PrintWriter(new FileWriter(proj.getDir(Project.RESULTS_DIRECTORY, true, false)+"pseudoautosomalSearch.xln"));
			writer.println("SNP\tChr\tPosition\tmLRR_M\tmLRR_F\thet_M\thet_F\tmiss_M\tmiss_F");
			
			markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markerList);
			time = new Date().getTime();
			line = "";
			for (int i = 0; i < markerList.length; i++) {
				markerData = markerDataLoader.requestMarkerData(i);

				markerName = markerData.getMarkerName();
				lrrs = markerData.getLRRs();
				abGenotypes = markerData.getAbGenotypesAfterFilters(clusterFilterCollection, markerName, gcThreshold);
				
				genotypeCounts = new int[2][4]; // sex, genotype
				values = new DoubleVector[2]; // sex
				values[0] = new DoubleVector();
				values[1] = new DoubleVector();
				for (int s = 0; s < samples.length; s++) {
					if (ext.isValidDouble(lrrs[s]+"") && !samplesToExclude[s]) {
						if (sexes[s] == 1 || sexes[s] == 2) {
							values[sexes[s]-1].add(lrrs[s]);
							genotypeCounts[sexes[s]-1][abGenotypes[s]+1]++;
						}
					}
				}

				line += markerName +"\t"+ markerData.getChr() +"\t"+ markerData.getPosition();
				if (values[0].size() > 0) {
					line += "\t"+Array.mean(values[0].toArray());
				} else {
					line += "\t.";
				}
				if (values[1].size() > 0) {
					line += "\t"+Array.mean(values[1].toArray());
				} else {
					line += "\t.";
				}
				if (genotypeCounts[0][1]+genotypeCounts[0][2]+genotypeCounts[0][3] > 0) {
					line += "\t"+(double)genotypeCounts[0][2]/(double)(genotypeCounts[0][1]+genotypeCounts[0][2]+genotypeCounts[0][3]);
				} else {
					line += "\t.";
				}
				if (genotypeCounts[1][1]+genotypeCounts[1][2]+genotypeCounts[1][3] > 0) {
					line += "\t"+(double)genotypeCounts[1][2]/(double)(genotypeCounts[1][1]+genotypeCounts[1][2]+genotypeCounts[1][3]);
				} else {
					line += "\t.";
				}
				line += "\t"+genotypeCounts[0][0];
				line += "\t"+genotypeCounts[1][0];
				line += eol;
				
				if (line.length() > 25000) {
					writer.print(line);
					writer.flush();
					line = "";
				}
				markerDataLoader.releaseIndex(i);
			}
			writer.print(line);
			System.out.println("Finished analyzing "+markerList.length+" in "+ext.getTimeElapsed(time));

			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing results");
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
		String filename = null;
		boolean par = false;

		String usage = "\\n"+
		"qc.GenderChecks requires 0-1 arguments\n"+
		"   (1) project properties filename (i.e. proj="+cnv.Launch.getDefaultDebugProjectFile(false)+" (default))\n"+
		"   (2) check sex of indiviudals (i.e. -check (not the default))\n"+
		" OR\n"+
		"   (2) parse all results (i.e. -parse (not the default))\n"+
		"   (3) drop markers (i.e. -drop (not the default))\n"+
		"   (4) file with all markers (i.e. all="+allMarkers+" (default file))\n"+
		"   (5) list of bad markers (i.e. drop="+markersToDrop+" (default file))\n"+
		" OR\n"+
		"   (2) check sex chromosomes for pseudoautosomal regions (i.e. -PARcheck (not the default))\n"+
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
			} else if (args[i].startsWith("-drop")) {
				drop = true ;
				numArgs--;
			} else if (args[i].startsWith("all=")) {
				allMarkers = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("drop=")) {
				markersToDrop = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("-PARcheck")) {
				par = true ;
				numArgs--;
			}
		}

		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}

//		check = true;
//		par = true;
//		filename = "D:/home/npankrat/projects/GEDI_exomeRAF.properties";
		try {
			proj = new Project(filename, false);
			
			if (check) {
				sexCheck(proj);
			} else if (parse) {
				parse(proj);
			} else if (par) {
				identifyPseudoautosomalBreakpoints(proj);
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