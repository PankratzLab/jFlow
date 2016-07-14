package cnv.analysis.pca;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Hashtable;

import cnv.filesys.Project;
import cnv.var.SampleData;
import common.Array;
import common.Files;
import common.HashVec;
import common.Logger;
import common.ext;

public class PCImputeRace {
	public static final String[] RACES = {"White", "African American", "Hispanic", "Asian"};
	public static final String[] STEP_PCS_HEADER = {"FID", "IID", "PC1", "PC2"};
	public static final String[] CORRECTED_PCS_HEADER = {"FID", "IID", "PC1", "PC2", "%African", "%Asian", "%White"};
	public static final String[] IMPUTED_RACE_SAMPLE_DATA_HEADERS = new String[] {"Class=ImputedRace;1=White;2=African American;3=Hispanic;4=Asian", "% African", "% Asian", "% European"};
	
	private Project proj;
	
	private String[] fidiids;
	
	private double[] pc1;
	private double[] pc2;
	
	private int[] europeanSeeds;
	private int[] africanSeeds;
	private int[] asianSeeds;
	
	private Logger log;
	
	
	
	/**
	 * @param fidiids
	 * @param pc1
	 * @param pc2
	 * @param europeanSeeds
	 * @param africanSeeds
	 * @param asianSeeds
	 * @param log
	 */
	public PCImputeRace(Project proj, String[] fidiids, double[] pc1, double[] pc2, int[] europeanSeeds, int[] africanSeeds,
			int[] asianSeeds, Logger log) {
		super();
		this.proj = proj;
		this.fidiids = fidiids;
		this.pc1 = pc1;
		this.pc2 = pc2;
		this.europeanSeeds = europeanSeeds;
		this.africanSeeds = africanSeeds;
		this.asianSeeds = asianSeeds;
		this.log = log;
	}

	public void correctPCsToRace(String outFile) {
		PrintWriter writer;
		if (pc1.length != pc2.length) {
			log.reportTimeError("PC1 and PC2 arrays must be the same length, 1 entry for each sample");
			return;
		}
		
		log.report("Checking for PC1 and PC2 predictions of African and Asian");
		
		if (!maybeSwapPCs()) return;
		
		writer = Files.getAppropriateWriter(ext.rootOf(outFile, false) + "_step0.mds");
		writer.println(Array.toStr(STEP_PCS_HEADER));
		
		for (int i = 0; i < pc1.length; i++) {
			writer.println(fidiids[i] + "\t" + pc1[i] + "\t" + pc2[i] + "\t" + (1.0 - (pc1[i] + pc2[i])));
		}
		writer.close();
		
		log.report("Setting European Seeds to below 0,0");
		
		
		double maxEuropeanPC1 = Double.NEGATIVE_INFINITY;
		double maxEuropeanPC2 = Double.NEGATIVE_INFINITY;
		
		for (int i : europeanSeeds) {
			if (pc1[i] > maxEuropeanPC1) {
				maxEuropeanPC1 = pc1[i];
			}
			if (pc2[i] > maxEuropeanPC2) {
				maxEuropeanPC2 = pc2[i];
			}
		}
		
		pc1 = Array.minus(pc1, maxEuropeanPC1);
		pc2 = Array.minus(pc2, maxEuropeanPC2);
		
		writer = Files.getAppropriateWriter(ext.rootOf(outFile, false) + "_step1.mds");
		writer.println(Array.toStr(STEP_PCS_HEADER));
		
		for (int i = 0; i < pc1.length; i++) {
			writer.println(fidiids[i] + "\t" + pc1[i] + "\t" + pc2[i] + "\t" + (1.0 - (pc1[i] + pc2[i])));
		}
		writer.close();
		
		log.report("Converting to Polar Coordinates");
		
		double[][] polar = cartesianToPolar(pc1, pc2);
		
		log.report("Forcing 90 degree difference between Asians and Africans");
		
		double maxAfricanTheta = Double.NEGATIVE_INFINITY;
		for (int i : africanSeeds) {
			if (polar[i][0] > maxAfricanTheta) {
				maxAfricanTheta = polar[i][0];
			}
		}
		
		double minAsianTheta = Double.MAX_VALUE;
		for (int i : asianSeeds) {
			if (polar[i][0] < minAsianTheta) {
				minAsianTheta = polar[i][0];
			}
		}
		
		for (int i = 0; i < polar.length; i++) {
			double theta = polar[i][0];
			if (theta <= minAsianTheta && theta >= maxAfricanTheta) {
				// If in target range, scale such that minAsianTheta becomes pi/2 (90) and maxAfricanTheta becomes 0 (0)
				polar[i][0] = ((theta - maxAfricanTheta) / (minAsianTheta - maxAfricanTheta)) * (Math.PI / 2.0);
			} else {
				// If out of target range, scale opposite such that minAsianTheta becomes pi/2 (90) and maxAfricanTheta becomes 2*pi (360)
				if (theta < maxAfricanTheta) {
					// Values from pi (180) to 2*pi (360) are expressed from -pi (-180) to 0 (0), express all values under maxAfricanTheta as >pi (180) for scaling
					theta = 2.0 * Math.PI + theta;
				}
				polar[i][0] = ((theta - minAsianTheta) / (2.0 * Math.PI + maxAfricanTheta - minAsianTheta )) * (3.0 * Math.PI / 2.0) + (Math.PI / 2.0) ;
			}
		}
		
		log.report("Converting back to Cartesian");
		
		double[][] rect = polarToCartesian(polar);
		
		pc1 = rect[0];
		pc2 = rect[1];
		
		writer = Files.getAppropriateWriter(ext.rootOf(outFile, false) + "_step2.mds");
		writer.println(Array.toStr(STEP_PCS_HEADER));
		
		for (int i = 0; i < pc1.length; i++) {
			writer.println(fidiids[i] + "\t" + pc1[i] + "\t" + pc2[i] + "\t" + (1.0 - (pc1[i] + pc2[i])));
		}
		writer.close();
		
		log.report("Forcing scale where African and Asian are at 1.0");
		
		double minAfricanPC1 = Double.MAX_VALUE;
		for (int i : africanSeeds) {
			if (pc1[i] < minAfricanPC1) {
				minAfricanPC1 = pc1[i];
			}
		}
		
		double minAsianPC2 = Double.MAX_VALUE;
		for (int i : asianSeeds) {
			if (pc2[i] < minAsianPC2) {
				minAsianPC2 = pc2[i];
			}
		}
		
		for (int i = 0; i < pc1.length; i++) {
			pc1[i] = pc1[i] / minAfricanPC1;
			pc2[i] = pc2[i] / minAsianPC2;
		}
		double[] pctAfrican = new double[pc1.length];
		double[] pctAsian = new double[pc1.length];
		double[] pctEuropean = new double[pc1.length];
		int[] imputedRace = new int[pc1.length];
		
		for (int i = 0; i < pc1.length; i++) {
			pctAfrican[i] = Math.min(Math.max(pc1[i], 0.0),1.0);
			pctAsian[i] = Math.min(Math.max(pc2[i], 0.0),1.0);
			pctEuropean[i] = 1.0 - (pctAfrican[i] + pctAsian[i]);
			if (pctEuropean[i] > 0.95) {
				imputedRace[i] = 1;
			} else if (pctAsian[i] > 0.999) {
				imputedRace[i] = 4;
			} else if (pctAfrican[i] / Math.max(pctAsian[i], 0.015) > 10) {
				imputedRace[i] = 2;
			} else {
				imputedRace[i] = 3;
			}
		}
		
		log.report("Writing Results");
		
		SampleData sampleData = proj.getSampleData(0, false);
		
		Hashtable<String, String> dataToAdd = new Hashtable<String, String>();
		for (int i = 0; i < pc1.length; i++) {
			dataToAdd.put(sampleData.lookup(fidiids[i])[0], imputedRace[i] + "\t" +
															pctAfrican[i] + "\t" +
															pctAsian[i] + "\t" +
															pctEuropean[i]);
		}
		
		
		sampleData.addData(dataToAdd, "DNA", IMPUTED_RACE_SAMPLE_DATA_HEADERS, ".", "\t", log);
		
		writer = Files.getAppropriateWriter(outFile);
		writer.println(Array.toStr(CORRECTED_PCS_HEADER));
		
		for (int i = 0; i < pc1.length; i++) {
			writer.println(fidiids[i] + "\t" + pc1[i] + "\t" + pc2[i] + "\t" + pctAfrican[i] + "\t" + pctAsian[i] + "\t" + pctEuropean[i]);
		}
		writer.close();
		
	}
	
	private boolean maybeSwapPCs() {
		double[] europeanPC1 = new double[europeanSeeds.length], europeanPC2 = new double[europeanSeeds.length];
		double[] africanPC1 = new double[africanSeeds.length], africanPC2 = new double[africanSeeds.length];
		double[] asianPC1 = new double[asianSeeds.length], asianPC2 = new double[asianSeeds.length];
		
		for (int i = 0; i < europeanSeeds.length; i++) {
			europeanPC1[i] = pc1[europeanSeeds[i]];
			europeanPC2[i] = pc2[europeanSeeds[i]];
		}
		
		for (int i = 0; i < africanSeeds.length; i++) {
			africanPC1[i] = pc1[africanSeeds[i]];
			africanPC2[i] = pc2[africanSeeds[i]];
		}
		
		for (int i = 0; i < asianSeeds.length; i++) {
			asianPC1[i] = pc1[asianSeeds[i]];
			asianPC2[i] = pc2[asianSeeds[i]];
		}
		
		double europeanMeanPC1 = Array.mean(europeanPC1, true);
		double europeanMeanPC2 = Array.mean(europeanPC2, true);
		double africanMeanPC1 = Array.mean(africanPC1, true);
		double africanMeanPC2 = Array.mean(africanPC2, true);
		double asianMeanPC1 = Array.mean(asianPC1, true);
		double asianMeanPC2 = Array.mean(asianPC2, true);
		
		if (Math.abs(africanMeanPC1) > Math.abs(asianMeanPC1) 
				&& Math.abs(asianMeanPC2) > Math.abs(africanMeanPC2)
				&& Math.abs(africanMeanPC1) > Math.abs(europeanMeanPC1) 
				&& Math.abs(asianMeanPC2) > Math.abs(europeanMeanPC2)) {
			// PC1 = African, PC2 = Asian
		} else if (Math.abs(asianMeanPC1) > Math.abs(africanMeanPC1) 
				&& Math.abs(africanMeanPC2) > Math.abs(asianMeanPC2) 
				&& Math.abs(asianMeanPC1) > Math.abs(europeanMeanPC1) 
				&& Math.abs(africanMeanPC2) > Math.abs(europeanMeanPC2)) {
			// PC1 = Asian, PC2 = African
			double[] temp;
			
			temp = pc1;
			pc1 = pc2;
			pc2 = temp;
			
			double tempMean;
			
			tempMean = europeanMeanPC1;
			europeanMeanPC1 = europeanMeanPC2;
			europeanMeanPC2 = tempMean;
			
			tempMean = africanMeanPC1;
			africanMeanPC1 = africanMeanPC2;
			africanMeanPC2 = tempMean;
			
			tempMean = asianMeanPC1;
			asianMeanPC1 = asianMeanPC2;
			asianMeanPC2 = tempMean;
		} else {
			log.reportError("PC1 and PC2 do not appear to predict African and Asian, race cannot be imputed");
			return false;
		}
		
		if (africanMeanPC1 < 0) {
			pc1 = Array.multiply(pc1, -1.0);

		}
		
		if (asianMeanPC2 < 0) {
			pc2 = Array.multiply(pc2, -1.0);
		}
		return true;
	}
	
	public static double[][] cartesianToPolar(double[] x, double[] y) {
		double[][] polar = new double[x.length][2];
		for (int i = 0 ; i < x.length; i++) {
			polar[i][0] = Math.atan2(y[i], x[i]);
			polar[i][1] = Math.hypot(x[i], y[i]);
		}
		return polar;
	}
	
	public static double[][] polarToCartesian(double[][] polar) {
		double[] x = new double[polar.length];
		double[] y = new double[polar.length];
		
		for (int i = 0; i < polar.length; i++) {
			double theta = polar[i][0];
			double r = polar[i][1];
			x[i] = r * Math.cos(theta);
			y[i] = r * Math.sin(theta);
		}
		
		return new double[][] {x, y};
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		Project proj = null;
		String inFile = null;

		String usage = "\n" + "cnv.analysis.pca.PCImputeRace requires 2 arguments\n";
		usage += "   (1) Project Filename (i.e. proj=" + inFile + " (default))\n" + "";
		usage += "   (2) Input Filename (i.e. inFile=" + inFile + " (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				proj = new Project(args[i].split("=")[1], false);
				numArgs--;
			} else if (args[i].startsWith("inFile=")) {
				inFile = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			String[] input = HashVec.loadFileToStringArray(inFile, true, null, false);
			String[] fidiids = new String[input.length];
			double[] pc1 = new double[input.length];
			double[] pc2 = new double[input.length];
			ArrayList<Integer> europeans = new ArrayList<Integer>();
			ArrayList<Integer> africans = new ArrayList<Integer>();
			ArrayList<Integer> asians = new ArrayList<Integer>();
			
			for (int i = 0;  i < input.length; i++) {
				String[] line = input[i].split("\t");
				fidiids[i] = line[0] + "\t" + line[1];
				pc1[i] = Double.parseDouble(line[2]);
				pc2[i] = Double.parseDouble(line[3]);
				try {
					int race = Integer.parseInt(line[4]);
					switch (race) {
					case 1:
						europeans.add(i);
						break;
					case 2:
						africans.add(i);
						break;
					case 3:
						asians.add(i);
						break;
					default:
						break;
					}
				} catch (NumberFormatException nfe) { }
			}

			PCImputeRace raceChecker = new PCImputeRace(proj,
					fidiids,
					pc1,
					pc2, 
					Array.toIntArray(europeans), 
					Array.toIntArray(africans),
					Array.toIntArray(asians), 
					new Logger());
			raceChecker.correctPCsToRace(ext.rootOf(inFile, false) + "_Corrected_PCS.mds");
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
}
