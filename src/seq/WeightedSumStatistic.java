package seq;

import java.io.*;

import common.*;
import filesys.*;
import stats.*;

public class WeightedSumStatistic {
	public static final int NUM_REPS = 10000;

	private Pedfile pedfile;
	private SnpMarkerSet map;
	private String[][] ids;
	private String[] markerNames;
	private double[] freqs;
	private double[] controlFreqs;
	private double[] weights;
	private double[] scores;
	private char[][] alleles;
	private int[] genotypeCounts;
	private int[][] markerCounts;
	private byte[] affectionStatus;
	private double[] trait;
	private boolean binary;
	private double stat;
	private double[] perms;
	private double sig;
	private int[] subset;
	private int[] valid_subset;
	private double mafThreshold;

	public WeightedSumStatistic(String plinkRoot) {
		pedfile = new Pedfile(plinkRoot+".ped");
		map = new SnpMarkerSet(plinkRoot+".map");
		ids = pedfile.getFamilyStructure().getIDs();
		markerNames = map.getMarkerNames();
		freqs = pedfile.getFrequencies();
		alleles = pedfile.getAlleles();
		genotypeCounts = pedfile.getGenotypeCounts();
		markerCounts = pedfile.getMarkerCounts();
		affectionStatus = pedfile.getFamilyStructure().getAffections();
		controlFreqs = getControlFreqs();
		subset = Array.intArray(markerNames.length);
		valid_subset = subset;
		trait = Array.toDoubleArray(affectionStatus);
		binary = true;
		mafThreshold = 1;
		reset();
	}
	
	public String[][] getIDs() {
		return ids;
	}

	public String[] getMarkerNames() {
		return Array.subArray(markerNames, valid_subset);
	}

	public double[] getFreqs() {
		return freqs;
	}

	public double[] getControlFreqs() {
		int count;
		
		if (controlFreqs == null) {
			controlFreqs = new double[freqs.length];
			for (int i = 0; i < controlFreqs.length; i++) {
				count = 0;
				for (int j = 0; j < markerCounts.length; j++) {
					if (markerCounts[j][i] != -1) {
						controlFreqs[i] += markerCounts[j][i];
						count++;
					}
				}
				controlFreqs[i] /= count;
			}
		}
		
		return controlFreqs;
	}

	public double[] getWeights() {
		if (weights == null) {
			weights = computeWeights(affectionStatus, markerCounts, valid_subset);
		}
		return weights;
	}

	public double[] getScores() {
		weights = getWeights();
		if (scores == null) {
	        scores = computeScores(markerCounts, weights, valid_subset);
		}
		return scores;
	}

	public double getStat() {
		if (stat == Double.MIN_VALUE) {
			scores = getScores();
			if (binary) {
				stat = sumOfRanks(trait, scores);
			} else {
				stat = Correlation.Spearman(new double[][] {trait, scores})[0];
//				stat = Correlation.Pearson(new double[][] {trait, scores})[0];
			}
		}
		
		return stat;
	}

	public double getSig() {
		double mean, stdev, z;
		
		perms = getPermutations();
		if (sig == Double.MAX_VALUE) {
	        mean = Array.mean(perms);
	        stdev = Array.stdev(perms);
	        z = (stat - mean)/stdev;
	        
	        sig = ProbDist.NormDist(Math.abs(z)*-1);
		}
		
		return sig;
	}

	public double[] getPermutations() {
		int[] order;
		double[] permWeights, permScores;
		
		if (perms == null) {
			perms = new double[NUM_REPS];
	        for (int i = 0; i<NUM_REPS; i++) {
	        	order = Array.random(affectionStatus.length);
	        	permWeights = computeWeights(Sort.putInOrder(affectionStatus, order), markerCounts, valid_subset);
	        	permScores = computeScores(markerCounts, permWeights, valid_subset);
				if (binary) {
					perms[i] = sumOfRanks(Sort.putInOrder(trait, order), permScores);
				} else {
					perms[i] = Correlation.Spearman(new double[][] {Sort.putInOrder(trait, order), permScores})[0];
//					perms[i] = Correlation.Pearson(new double[][] {Sort.putInOrder(trait, order), permScores})[0];
				}	        	
	        }
		}
		
		return perms;
	}
	
	public double getEmpiricalSig() {
        int count;
        
        stat = getStat();
        if ((stat+"").equalsIgnoreCase("NaN")) {
        	return 1.1;
        }
        
        count = 0;
        perms = getPermutations();
        for (int i = 0; i<perms.length; i++) {
        	if (binary) {
            	if (perms[i] <= stat) {
            		count++;
            	}
        	} else {
            	if (perms[i] >= stat) {
            		count++;
            	}
        	}
        }

        return (double)(count+1)/(double)(NUM_REPS+1);
	}

	
	public char[][] getAlleles() {
		return alleles;
	}

	public int[] getGenotypeCounts() {
		return genotypeCounts;
	}

	public int[][] getMarkerCounts() {
		return markerCounts;
	}

	public void setThresholdMAF(double mafThreshold) {
		this.mafThreshold = mafThreshold;
		validateSubset();
		reset();
	}
	
	public void setAffectionStatus(byte[] affectionStatus) {
		this.affectionStatus = affectionStatus;
		reset();
	}
	
	public void setSubset(String[] markers) {
		this.subset = ext.indexFactors(markers, markerNames, false, true);
		validateSubset();
		reset();
	}
	
	public void setTrait(double[] trait) {
		int type;
		
		this.trait = trait;
		type = Array.determineType(trait, true);
		if (type == -1) {
			System.exit(1);
		}
		binary = type==0;
		stat = Double.MIN_VALUE;
		sig = Double.MAX_VALUE;
	}
	
	public void validateSubset() {
		boolean[] passes;
		int count;
		
		controlFreqs = getControlFreqs();
		if (mafThreshold < 1) {
			passes = new boolean[subset.length];
			for (int i = 0; i < subset.length; i++) {
				passes[i] = controlFreqs[subset[i]] < mafThreshold;
//				passes[i] = freqs[subset[i]] < mafThreshold;
			}
			valid_subset = new int[Array.booleanArraySum(passes)];
			if (valid_subset.length < subset.length) {
				count = 0;
				for (int i = 0; i < subset.length; i++) {
					if (passes[i]) {
						valid_subset[count] = subset[i];
						count++;
					}
				}
				subset = valid_subset;
			}			
		} else {
			valid_subset = subset;
		}
	}

	public void reset() {
		weights = null;
		scores = null;
		perms = null;
		stat = Double.MIN_VALUE;
		sig = Double.MAX_VALUE;
	}
	
	public static void demo(String dir, String root) {
        PrintWriter writer;
        int count;
		double actualX;
		double[] repXs;
		
        double mean, stdev, z;

        WeightedSumStatistic wss;
    	String[] markerNames;
    	double[] freqs;
    	double[] weights;
    	double[] scores;
    	char[][] alleles;
    	int[] genotypeCounts;
    	int[][] markerCounts;

    	wss = new WeightedSumStatistic(dir+root);
    	
    	markerNames = wss.getMarkerNames();
    	freqs = wss.getFreqs();
    	weights = wss.getWeights();
    	scores = wss.getScores();
    	alleles = wss.getAlleles();
    	markerCounts = wss.getMarkerCounts();
    	genotypeCounts = wss.getGenotypeCounts();
				
		try {
	        writer = new PrintWriter(new FileWriter(dir+root+"_weights.xln"));
	        writer.println("SNP\tA1\tA2\tMAF\tNCHROBS\tweight");
	        for (int i = 0; i<markerNames.length; i++) {
	        	writer.println(markerNames[i]+"\t"+alleles[i][0]+"\t"+alleles[i][1]+"\t"+freqs[i]+"\t"+genotypeCounts[i]+"\t"+weights[i]);
            }
	        writer.close();
        } catch (Exception e) {
	        System.err.println("Error writing to "+dir+root+"_weights.xln");
	        e.printStackTrace();
        }
        scores = wss.getScores();

        actualX = wss.getStat();
		try {
	        writer = new PrintWriter(new FileWriter(dir+root+"_scores.xln"));
	        writer.println("Score\t"+Array.toStr(markerNames));
	        writer.println("Weights:\t"+Array.toStr(weights));
	        for (int i = 0; i<scores.length; i++) {
	        	writer.println(scores[i]+"\t"+Array.toStr(markerCounts[i]));
            }
	        writer.close();
        } catch (Exception e) {
	        System.err.println("Error writing to "+dir+root+"_scores.xln");
	        e.printStackTrace();
        }
        
        System.out.println("N SNPs: "+markerNames.length);
        System.out.println("Actual X: "+actualX);

        repXs = wss.getPermutations();
        
        mean = Array.mean(repXs);
        stdev = Array.stdev(repXs);
        z = (actualX - mean)/stdev;
        
        System.out.println("Mean null X: "+mean);
        System.out.println("Stdev null X: "+stdev);
        System.out.println("Z: "+ext.formDeci(z, 3));
        System.out.println("p-value for z-score: "+ext.prettyP(ProbDist.NormDist(Math.abs(z)*-1)));
        count = 0;
        for (int i = 0; i<repXs.length; i++) {
        	if (repXs[i] < actualX) {
        		count++;
        	}
        }
        System.out.println("empirical p-value: "+ext.prettyP((double)count/(double)NUM_REPS));

        System.out.println(markerNames.length+"\t"+ext.formDeci(z, 3)+"\t"+ext.prettyP(ProbDist.NormDist(Math.abs(z)*-1))+"\t"+ext.prettyP((double)count/(double)NUM_REPS));

        try {
	        writer = new PrintWriter(new FileWriter("logfile2.out", true));
	        writer.println(markerNames.length+"\t"+ext.formDeci(z, 3)+"\t"+ext.prettyP(ProbDist.NormDist(Math.abs(z)*-1))+"\t"+ext.prettyP((double)count/(double)NUM_REPS));
	        writer.println();
	        writer.close();
        } catch (Exception e) {
	        System.err.println("Error writing to "+"logfile.out");
	        e.printStackTrace();
        }
	}
	
//	public static double computeX(byte[] affectionStatus, int[][] markerCounts) {
//        double[] weights, scores;
//        
//		weights = computeWSS_weights(affectionStatus, markerCounts);
//        scores = computeWSS_scores(markerCounts, weights);
//
//        return computeWSS_sumOfRanks(affectionStatus, scores);
//	}
	
	public static double sumOfRanks(double[] affectionStatus, double[] scores) {
        int[] ranks;
        double sum;
        
        ranks = Sort.ranks(scores, Sort.DESCENDING);
        sum = 0;
        for (int i = 0; i<ranks.length; i++) {
        	if (affectionStatus[i] == 2) {
        		sum += ranks[i];
        	}
        }
        
		return sum;
	}
	
	public static double[] computeScores(int[][] markerCounts, double[] weights, int[] subset) {
		double[] scores;
		
		scores = new double[markerCounts.length];
		for (int i = 0; i<markerCounts.length; i++) {
			for (int j = 0; j<subset.length; j++) {
				if (markerCounts[i][subset[j]] != -1) {
					scores[i] += (double)markerCounts[i][subset[j]]/weights[j];
				}
            }
        }
		
		return scores;
	}

	public static double[] computeWeights(byte[] affectionStatus, int[][] markerCounts, int[] subset) {
		double[] weights;
		int numGenotyped, numVariantsInUnaffecteds, numUnaffecteds;
		
		weights = new double[subset.length];
		for (int i = 0; i<subset.length; i++) {
			numGenotyped = numVariantsInUnaffecteds = numUnaffecteds = 0;
			for (int j = 0; j<affectionStatus.length; j++) {
				if (markerCounts[j][subset[i]] >= 0) {
					numGenotyped++;
					if (affectionStatus[j] == 1) {
						numUnaffecteds++;
						numVariantsInUnaffecteds += markerCounts[j][subset[i]];
					}
				}
            }
			weights[i] = computeWSS_weight(numGenotyped, numVariantsInUnaffecteds, numUnaffecteds);
        }
		
		return weights;
	}
	
	public static double computeWSS_weight(int numGenotyped, int numVariantsInUnaffecteds, int numUnaffecteds) {
		double q;
		
		q = ((double)numVariantsInUnaffecteds+1)/((double)numUnaffecteds + 2);
		if (q > 1) {
//			System.err.println("Error - weight is not accurate and can be undefined if there are more variants than controls (i.e. control MAF > 0.50)");
		}

		return Math.sqrt(numGenotyped*q*(1-q));
	}
	
	public static void batch(String modelsFile, String filename) {
        PrintWriter writer;
        String[] line;
        String[] models;
        
        models = HashVec.loadFileToStringArray(modelsFile, false, new int[] {0,1}, false);        
        try {
	        writer = new PrintWriter(new FileWriter(filename));
	        for (int i = 0; i<models.length; i++) {
	        	line = models[i].split("[\\s]+");
		        writer.println("plink --bfile allThree --recode --keep "+line[0]+" --extract "+line[1]);
		        writer.println("java -cp /home/npankrat/" + common.PSF.Java.GENVISIS + " seq.WeightedSumStatistic set="+line[1]);
		        writer.println();
            }
	        
	        writer.close();
        } catch (Exception e) {
	        System.err.println("Error writing to "+filename);
	        e.printStackTrace();
        }

	}
	
	public static void main(String[] args) {
	    int numArgs = args.length;
//	    String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\UMN\\Myron\\ExcisionPathway\\RareVariants\\";
//	    String dir = "D:\\RareVariants\\";
	    String dir = "";
	    String root = "plink";
//	    String root = "func_plink";
	    String models = "models.dat";

	    String usage = "\n"+
	    "seq.WeightedSumStatistic requires 0-1 arguments\n"+
	    "   (1) filename (i.e. file="+root+" (default))\n"+
	    "";

	    for (int i = 0; i<args.length; i++) {
		    if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
			    System.err.println(usage);
			    System.exit(1);
		    } else if (args[i].startsWith("root=")) {
			    root = args[i].split("=")[1];
			    numArgs--;
		    } else if (args[i].startsWith("set=")) {
		    	try {
			    	PrintWriter writer = new PrintWriter(new FileWriter("logfile2.out", true));
			        writer.println(args[i].split("=")[1]);
			        writer.close();
				    numArgs--;
				} catch (Exception e) {
					System.err.println("Error writing");
					e.printStackTrace();
				}
		    }
	    }
	    if (numArgs!=0) {
		    System.err.println(usage);
		    System.exit(1);
	    }
	    try {
	    	if (new File("models.dat").exists()) {
	    		batch(models, "batchModels.bat");
	    	} else {
	    		demo(dir, root);
	    	}
	    } catch (Exception e) {
		    e.printStackTrace();
	    }
    }
}
