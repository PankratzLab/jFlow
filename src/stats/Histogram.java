package stats;

import java.io.*;
import common.*;

public class Histogram {
	public final static int DEFAULT_WINDOW = 7;
	public final static float DEFAULT_PEAK_THRESHOLD = 0.02f;

	private int[] counts;
	private int sumTotal;
	private float[] smoothed;
	private double min;
	private double max;
	private int sigfigs;
	private int window;
	private float peakThreshold;
	
	public Histogram(double[] array, double min, double max, int sigfigs) {
		this(array, array.length, min, max, sigfigs);
	}
	
	public Histogram(double[] array, int arraySize, double min, double max, int sigfigs) {
		this.min = min;
		this.max = max;
		this.sigfigs = sigfigs;
		this.window = DEFAULT_WINDOW;
		this.peakThreshold = DEFAULT_PEAK_THRESHOLD;
		this.sumTotal = array.length;

		counts = new int[(int)((max-min)*Math.pow(10, sigfigs))+1];

		for (int i = 0; i<arraySize; i++) {
			if (Double.isNaN(array[i])) {
			
			} else if (array[i] < min) {
				counts[0]++;
			} else if (array[i] > max) {
				counts[counts.length-1]++;
			} else {
				counts[(int)(array[i]*Math.pow(10, sigfigs)-min*Math.pow(10, sigfigs))]++;
			}
        }
	}

	/**
	 * When the data to be used is not known before hand (many many data points from a .bam file, etc...) This is generally the constructor used for {@link DynamicHistogram}
	 */
	private Histogram(double min, double max, int sigfigs) {
		this.min = min;
		this.max = max;
		this.sigfigs = sigfigs;
		this.window = DEFAULT_WINDOW;
		this.peakThreshold = DEFAULT_PEAK_THRESHOLD;
		this.counts = new int[(int) ((max - min) * Math.pow(10, sigfigs)) + 1];
		this.sumTotal = 0;
	}
	
	public Histogram(float[] array, float min, float max, int sigfigs) {
		this(array, array.length, min, max, sigfigs);
	}

	public Histogram(float[] array, int arraySize, float min, float max, int sigfigs) {
		this.min = min;
		this.max = max;
		this.sigfigs = sigfigs;
		this.window = DEFAULT_WINDOW;
		this.peakThreshold = DEFAULT_PEAK_THRESHOLD;
		this.sumTotal = array.length;

		counts = new int[(int)((max-min)*Math.pow(10, sigfigs))+1];

		for (int i = 0; i<arraySize; i++) {
			if (Float.isNaN(array[i])) {

			} else if (array[i] < min) {
				counts[0]++;
			} else if (array[i] > max) {
				counts[counts.length-1]++;
			} else {
				counts[(int)(array[i]*Math.pow(10, sigfigs)-min*Math.pow(10, sigfigs))]++;
			}
        }
	}

	public double[] getBins() {
		double[] bins = new double[counts.length];
		int count = 0;
		double d = min, step = Math.pow(10, -1*sigfigs);
				
		while (count < counts.length && d<=max) {
			bins[count++] = d;
			d += step;
		}
		
		return bins;
	}

	public int[] getCounts() {
		return counts;
	}
	
	public float[] getSmoothed() {
		smooth();
		return smoothed;
	}
	
	public void setWindow(int window) {
		this.window = window;
	}
	
	public void setPeakThreshold(float peakThreshold) {
		this.peakThreshold = peakThreshold;
	}
	
	public void smooth() {
		float localSum, win;
		int ahead = (window-1)/2;
		int behind = (window-1) - ahead;
		
		if (smoothed != null) {
			return;
		}
		
		smoothed = new float[counts.length];
			
		localSum = 0;
		for (int i = 0; i<ahead; i++) {
			localSum += counts[i];
		}
		for (int i = 0; i<behind+1; i++) {
			localSum += counts[i+ahead];
			smoothed[i] = localSum/(float)(i+1+ahead);
		}
		win = (float)window;
		for (int i = behind+1; i<counts.length-ahead; i++) {
			localSum += counts[i+ahead] - counts[i-behind-1];
			smoothed[i] = localSum/win;
        }

		for (int i = counts.length-ahead; i<counts.length; i++) {
			localSum -= counts[i-behind-1];
			smoothed[i] = localSum/(float)(counts.length-(i-behind));
		}
		
	}
	
	public int getMaxIndex() {
		int max = -1, maxIndex = -1;
		
		for (int i = 0; i<counts.length; i++) {
			if (counts[i] > max) {
				max = counts[i];
				maxIndex = i;
			}
        }
		
		return maxIndex;
	}
	
	public double getMaxBin() {
		return min+getMaxIndex()*Math.pow(10, -1*sigfigs);
	}
	
	public int[] getLocalMaxima() {
		IntVector localMaxima = new IntVector();
		float sum;
		int lastMaxima;
		
		smooth();
		sum = (float)sumTotal;
		lastMaxima = -999;
		
		if (smoothed[0]/sum > peakThreshold && smoothed[0] > smoothed[1]) {
			localMaxima.add(0);
		} 
		for (int i = 1; i<smoothed.length-1; i++) {
			if (smoothed[i]/sum > peakThreshold && smoothed[i] >= smoothed[i-1] && smoothed[i] >= smoothed[i+1] && i - lastMaxima > 5) {
				localMaxima.add(i);
				lastMaxima = i;
			} 
        }
		if (smoothed[smoothed.length-1]/sum > peakThreshold && smoothed[smoothed.length-1] > smoothed[smoothed.length-2]) {
			localMaxima.add(smoothed.length-1);
		} 
		
		return localMaxima.toArray();
	}
	
	public int getNumPeaks() {
		return getLocalMaxima().length;
	}
	
	public void dump(String filename) {
        PrintWriter writer;
        double[] bins;
        float[] smoothd = getSmoothed();
        
		try {
			bins = getBins();
			writer = new PrintWriter(new FileWriter(filename));
			writer.println("Bin\tCounts");
			for (int i = 0; i<bins.length; i++) {
				writer.println(ext.formDeci(bins[i], sigfigs)+"\t"+counts[i]+"\t"+ext.formDeci(smoothd[i], 5));
            }
	        writer.close();
        } catch (Exception e) {
	        System.err.println("Error writing to "+filename);
	        e.printStackTrace();
        }
    }
	
	public double getMin() {
		return min;
	}

	public double getMax() {
		return max;
	}

	public int getSigfigs() {
		return sigfigs;
	}
	
	public void addToSumTotal(){
		sumTotal++;
	}

	/**
	 * Extension of {@link Histogram} that allows data to be added on the fly, useful when there are many millions of data points that do not need to be kept in memory
	 *
	 */
	public static class DynamicHistogram extends Histogram {

		public DynamicHistogram(double min, double max, int sigfigs) {
			super(min, max, sigfigs);
		}

		public void addDataPointToHistogram(double data) {
			addToSumTotal();
			if (Double.isNaN(data)) {

			} else if (data < getMin()) {
				getCounts()[0]++;
			} else if (data > getMax()) {
				getCounts()[getCounts().length - 1]++;
			} else {
				getCounts()[(int) (data * Math.pow(10, getSigfigs()) - getMin() * Math.pow(10, getSigfigs()))]++;
			}
		}

		/**
		 * If you would like to initialize several histograms to the same value
		 */
		public static DynamicHistogram[] initHistograms(int numHistograms, double min, double max, int sigFigs) {
			DynamicHistogram[] dynamicHistograms = new DynamicHistogram[numHistograms];
			for (int i = 0; i < dynamicHistograms.length; i++) {
				dynamicHistograms[i] = new DynamicHistogram(min, max, sigFigs);
			}
			return dynamicHistograms;
		}
	}

	public static void main(String[] args) {
	    int numArgs = args.length;
	    String filename = "Histogram.dat";

	    String usage = "\n"+"stats.Histogram requires 0-1 arguments\n"+"   (1) filename (i.e. file="+filename+" (default))\n"+"";

	    for (int i = 0; i<args.length; i++) {
		    if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
			    System.err.println(usage);
			    System.exit(1);
		    } else if (args[i].startsWith("file=")) {
			    filename = args[i].split("=")[1];
			    numArgs--;
		    }
	    }
	    if (numArgs!=0) {
		    System.err.println(usage);
		    System.exit(1);
	    }
	    try {
	    	float[] test = new float[] {0, 0, 0, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 8};
		    new Histogram(test, 0, 8, 0).dump("test.xln");
	    } catch (Exception e) {
		    e.printStackTrace();
	    }
    }
}
