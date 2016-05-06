package stats;

import java.io.*;

import common.*;

public class Histogram implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	public final static int DEFAULT_WINDOW = 7;
	public final static float DEFAULT_PEAK_THRESHOLD = 0.02f;
	
			// These correspond to steps of:    0.10     0.20                 0.25              0.50
			// To add more options, use POWER(10, -1*desired_value) and add it here
	public final static double[] EXTRA_STEPS = {0.0, 0.698970004336019, 0.602059991327962, 0.301029995663981};

	protected int[] counts;
	private int sumTotal;
	private float[] smoothed;
	protected double min;
	protected double max;
	protected int sigfigs;
	protected int extrastep;
	private int window;
	private float peakThreshold;
	private String title;
	
	public Histogram(double[] array) {
		this.window = DEFAULT_WINDOW;
		this.peakThreshold = DEFAULT_PEAK_THRESHOLD;
		this.sumTotal = array.length;
		
		min = Array.min(array);
		max = Array.max(array);

		sigfigs = 0;
		extrastep = 0;
		
		while ((max - min) / determineStep() < (double)array.length/50.0) {
			extrastep++;
			if (extrastep == EXTRA_STEPS.length) {
				sigfigs++;
				extrastep = 0;
			}
		}
		System.out.println("optimal step is "+sigfigs+" + "+EXTRA_STEPS[extrastep]+" = "+(sigfigs+EXTRA_STEPS[extrastep])+ " log10 transformed = "+determineStep());
		
		computeCounts(array);
	}

	public String getTitle() {
		return title;
	}

	public void setTitle(String title) {
		this.title = title;
	}

	public int getExtrastep() {
		return extrastep;
	}
	
	public Histogram(double[] array, int sigfigs, int extrastep) {
	    this(array, Array.min(array), Array.max(array), sigfigs, extrastep);
	}

	public Histogram(double[] array, double min, double max, int sigfigs, int extrastep) {
		this.min = min;
		this.max = max;
		this.sigfigs = sigfigs;
		this.extrastep = extrastep;
		this.window = DEFAULT_WINDOW;
		this.peakThreshold = DEFAULT_PEAK_THRESHOLD;
		this.sumTotal = array.length;

		computeCounts(array);
	}
	
	private void computeCounts(double[] array) {
		double start;
		double halfStep;
		
		start = determineStart();
		halfStep = determineStep()/2;
		int step = (int)((max-start)*Math.pow(10, sigfigs+EXTRA_STEPS[extrastep]))+2;
		counts = new int[step];
		
		System.out.println("min: "+min+"; max: "+max+"; start: "+start+"; step: "+determineStep(sigfigs, extrastep)+"; # bins: "+counts.length);


		for (int i = 0; i<array.length; i++) {
			if (Double.isNaN(array[i])) {
			
			} else if (array[i] < start) {
				counts[0]++;
			} else if (array[i] > max) {
				counts[counts.length-1]++;
			} else {
				int num = (int)((array[i]+halfStep)*Math.pow(10, sigfigs+EXTRA_STEPS[extrastep])-start*Math.pow(10, sigfigs+EXTRA_STEPS[extrastep]));
				if (num >= counts.length) {
					System.err.println("Error - trying to place value "+array[i]+" in position "+num+", but we only went up to "+counts.length);
				}
				counts[num]++;
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
		this.extrastep = 0;
		this.counts = new int[(int) ((max - min) * Math.pow(10, sigfigs)) + 1];
		this.sumTotal = 0;

		this.counts = new int[(int)((max-determineStart())*Math.pow(10, sigfigs))+1]; // does not currently use EXTRA_STEPS[extrastep]
	}
	
	public Histogram(float[] array, float min, float max, int sigfigs) {
		this(Array.toDoubleArray(array), min, max, sigfigs, 0);
	}
	
	protected Histogram() {
        // TODO Auto-generated constructor stub
    }

    public double determineStep() {
		return determineStep(sigfigs, extrastep);
	}
	
	public static double determineStep(int sigfigs, int extrastep) {
		return Double.parseDouble(ext.formDeci(Math.pow(10, -1*(sigfigs+EXTRA_STEPS[extrastep])), sigfigs+2));
	}
	
	public static int[] reverseStep(double step) {
		int[] trav, best;
		double diff, closestDiff;
		boolean done;
		
		best = null;
		done = false;
		trav = new int[] {0, EXTRA_STEPS.length-1};
		closestDiff = Double.POSITIVE_INFINITY;
		while (!done) {
			double travStep = determineStep(trav[0], trav[1]);
			System.out.println(trav[0]+","+trav[1]+"\t"+travStep);
			diff = Math.abs(determineStep(trav[0], trav[1])-step);
			if (diff < closestDiff) {
				best = Array.clone(trav);
				closestDiff = diff;
			} else if (diff > closestDiff) {
				done = true;
			}
			trav[1]--;
			if (trav[1] == 0) {
				trav[0]++;
			}
			if (trav[1] < 0) {
				trav[1] = EXTRA_STEPS.length-1;
			}
		}

		return best;
	}

	public double determineStart() {
		double d, step;

		d = 0;
		step = determineStep();
		if (d > min) {
			while (d > min) {
				d -= step;
			}
		} else {
			while (d + step < min) {
				d += step;
			}
		}

		return d;
	}

	public double[] getBins() {
		return getBins(false);
	}
	
	public double[] getBins(boolean verbose) {
		double[] bins = new double[counts.length];
		int count = 0;
		double d = min, step = Double.parseDouble(ext.formDeci(Math.pow(10, -1*(sigfigs+EXTRA_STEPS[extrastep])), sigfigs+1));
		
		d = determineStart();
		
		if (verbose) {
			System.out.println("min: " + min + "; max: " + max + "; start: " + d + "; step: " + step);
		}
		while (count < counts.length && d<=max) {
			bins[count++] = d;
			d += step;
		}
		
		return bins;
	}
	
	public String[] getBinsInString() {
		double[] bins;
		String[] result;
		
		bins = getBins();
		result = new String[bins.length];
		for (int i = 0; i < result.length; i++) {
			result[i] = ext.formDeci(bins[i], sigfigs+1+(extrastep==2?1:0));	
		}
		
		return result;
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
		return min+getMaxIndex()*Math.pow(10, -1*sigfigs+EXTRA_STEPS[extrastep]);
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
				writer.println(ext.formDeci(bins[i], sigfigs+(extrastep==2?1:0))+"\t"+counts[i]+"\t"+ext.formDeci(smoothd[i], 5));
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
	
	public String getSummary() {
        String[] bins;
        float[] smoothd;
		String result;
        
		bins = getBinsInString();
		smoothd = getSmoothed();
		
		result = "Bin\tCounts\tAvgOver"+window+"bins\r\n";
		for (int i = 0; i<bins.length; i++) {
			result += bins[i]+"\t"+counts[i]+"\t"+ext.formDeci(smoothd[i], 5)+"\r\n";
        }

		return result;
	}
	
	public void addToSumTotal(){
		sumTotal++;
	}

	/**
	 * Extension of {@link Histogram} that allows data to be added on the fly, useful when there are many millions of data points that do not need to be kept in memory
	 *
	 */
	public static class DynamicHistogram extends Histogram implements Serializable {

		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

		public DynamicHistogram(double min, double max, int sigfigs) {
			super(min, max, sigfigs);
		}

		public int addDataPointToHistogram(double data) {
			addToSumTotal();
			int indexAdded = -1;
			if (Double.isNaN(data)) {

			} else if (data < getMin()) {
				indexAdded = 0;
			} else if (data > getMax()) {
				indexAdded = getCounts().length - 1;
			} else {
				indexAdded = (int) ((data + determineStep() / 2) * Math.pow(10, getSigfigs()) - determineStart() * Math.pow(10, getSigfigs()));

			}
			if (indexAdded >= 0) {
				getCounts()[indexAdded]++;
			}

			return indexAdded;
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

		/**
		 * @param histograms
		 *            must be the same length
		 */
		public static DumpResult dumpToSameFile(DynamicHistogram[] histograms, String[] titles, String output,boolean proportion, Logger log) {
			double max = -1 * Double.MAX_VALUE;
			double min = Double.MAX_VALUE;
			if (histograms.length != titles.length) {
				log.reportTimeError("Titles and histogram array must be the same length");
				return null;
			}
			try {
		
				double[] totals = new double[histograms.length];
				for (int i = 0; i < totals.length; i++) {
					totals[i] = proportion ? Array.sum(histograms[i].getCounts()) : 1;
				}
				PrintWriter writer = new PrintWriter(new FileWriter(output));
				String[] cumulative = Array.tagOn(titles, "Cumulative_", null);
				String[] allTitles = new String[titles.length * 2];
				int index = 0;
				for (int i = 0; i < titles.length; i++) {
					allTitles[index] = titles[i];
					index++;
					allTitles[index] = cumulative[i];
					index++;
				}
				int[] cuCounts = new int[histograms.length];
				writer.println("Bin\t" + Array.toStr(allTitles));
				for (int i = 0; i < histograms[0].getBins().length; i++) {
					writer.print(histograms[0].getBins()[i]);
					for (int j = 0; j < histograms.length; j++) {
						double val =  (proportion ? histograms[j].getCounts()[i] / totals[j] : histograms[j].getCounts()[i]);
						cuCounts[j] += val;
						writer.print("\t" + val + "\t" + cuCounts[j]);
						if(val>max){
							max=val;
						}else if(val<min){
							min =val;
						}
					}
					writer.println();
				}
				writer.close();
			} catch (Exception e) {
				log.reportError("Error writing to " + output);
				log.reportException(e);
			}
			return new DumpResult(max, min);
		}
	}


	
	/**
	 * @author lane0212<br>
	 *         A histogram that also tracks the average of a separate paired metric...sort of like an "averageIf" in excel
	 */
	public static class DynamicAveragingHistogram extends DynamicHistogram {
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;
		public static final String DUMP_COUNT = "_Count";
		public static final String DUMP_PROP = "_Proportion";
		public static final String DUMP_AVG = "_Average";
		private double[] averages;

		public DynamicAveragingHistogram(double min, double max, int sigfigs) {
			super(min, max, sigfigs);
			this.averages = new double[getCounts().length];
		}

		public double[] getAverages() {

			return averages;
		}

		public void setAverages(double[] averages) {
			this.averages = averages;
		}

		/**
		 * @param dataBin
		 *            this determines which bin of the histogram will be tallied
		 * @param dataToAverage
		 *            this data point will be added to the the bins average
		 */
		public void addDataPair(double dataBin, double dataToAverage) {
			if (!Double.isNaN(dataToAverage)) {
				int binIndex = addDataPointToHistogram(dataBin);
				if (binIndex >= 0) {
					averages[binIndex] += dataToAverage;
				}
			}
		}

		public void average() {
			for (int i = 0; i < getCounts().length; i++) {
				if (getCounts()[i] > 0) {
					averages[i] = averages[i] / getCounts()[i];
				}
			}
		}

		/**
		 * @param filename to dump to 
		 * @param titles
		 *            alternate labels for Bin, Counts, Average in that order
		 * @param countsGtZero
		 *            only dump counts greater than 0
		 */
		public void dump(String filename, String[] titles, boolean countsGtZero) {

			try {
				double[] bins = getBins();

				PrintWriter writer = new PrintWriter(new FileWriter(filename));
				if (titles != null) {
					writer.println(Array.toStr(titles));
				} else {
					writer.println("Bin\tCounts\tAverage");
				}
				for (int i = 0; i < bins.length; i++) {
					if (!countsGtZero || getCounts()[i] > 0)
						writer.println(getBins()[i] + "\t" + getCounts()[i] + "\t" + averages[i]);
				}
				writer.close();
			} catch (Exception e) {
				System.err.println("Error writing to " + filename);
				e.printStackTrace();
			}
		}

		/**
		 * @param histograms
		 *            must be the same length
		 */
		public static void dumpToSameFile(DynamicAveragingHistogram[] histograms, String[] titles, String output, Logger log) {

			if (histograms.length != titles.length) {
				log.reportTimeError("Titles and histogram array must be the same length");
				return;
			}
			try {
				String[] finalTitles = Array.concatAll(Array.tagOn(titles, null, DUMP_COUNT), Array.tagOn(titles, null, DUMP_PROP), Array.tagOn(titles, null, DUMP_AVG));
				PrintWriter writer = new PrintWriter(new FileWriter(output));
				writer.println("Bin\t" + Array.toStr(finalTitles));
				double[] totals = new double[histograms.length];
				for (int i = 0; i < totals.length; i++) {
					totals[i] = Array.sum(histograms[i].getCounts());
				}
				for (int i = 0; i < histograms[0].getBins().length; i++) {
					double currentBin = histograms[0].getBins()[i];
					writer.print(currentBin);
					for (int j = 0; j < histograms.length; j++) {
						if (histograms[j].getBins()[i] != currentBin) {
							String error = "mismatched bins between histograms";
							log.reportTimeError(error);
							writer.close();
							throw new IllegalStateException(error);
						}
						double val = histograms[j].getCounts()[i];
						writer.print("\t" + val);
					}
					for (int j = 0; j < histograms.length; j++) {
						if (histograms[j].getBins()[i] != currentBin) {
							String error = "mismatched bins between histograms";
							log.reportTimeError(error);
							writer.close();
							throw new IllegalStateException(error);
						}
						double val = (double) histograms[j].getCounts()[i] / totals[j];
						writer.print("\t" + val);
					}
					for (int j = 0; j < histograms.length; j++) {
						if (histograms[j].getBins()[i] != currentBin) {
							String error = "mismatched bins between histograms";
							log.reportTimeError(error);
							writer.close();
							throw new IllegalStateException(error);
						}
						double val = histograms[j].getAverages()[i];
						writer.print("\t" + val);
					}
					writer.println();
				}
				writer.close();
			} catch (Exception e) {
				log.reportError("Error writing to " + output);
				log.reportException(e);
			}
		}
	}

	public static class DumpResult {
		private double maxCount;
		private double minCount;
		public DumpResult(double maxCount, double minCount) {
			super();
			this.maxCount = maxCount;
			this.minCount = minCount;
		}
		public double getMaxCount() {
			return maxCount;
		}
		public double getMinCount() {
			return minCount;
		}

	}

	public static void main(String[] args) {
	    int numArgs = args.length;
	    String filename = "Histogram.dat";

	    String usage = "\n"+
	    "stats.Histogram requires 0-1 arguments\n"+
	    "   (1) filename (i.e. file="+filename+" (default))\n"+
	    "";

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
//	    	float[] test = new float[] {0, 0, 0, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 8};
//		    new Histogram(test, 0, 8, 0).dump("test.xln");
	    	double[] test = new double[] {1.490768942, -1.094579018, 1.645354986, 0.013646307, -0.695236359, 1.262092982, -0.825076523, 1.634672757, 1.647866972, -0.271143259, -1.0760067, 1.23669905, -0.531958778, -0.049129586, 0.392974736, -0.272928712, 0.710031172, 1.602224773, -1.11794446, -2.128879598, -1.835109853, -0.377314082, -0.160942974, 0.055953153, 0.965978934, 1.204199968, -0.804212092, 0.290432333, 0.09123259, 1.796191396, 0.297993442, 0.651277312, 0.934207425, 1.211702172, -0.636430058, -1.414079671, -0.148228087, 0.351882719, 1.864783488, 0.074212634, -0.211102404, 0.271863083, -1.665709565, 0.08254936, -0.841179783, 0.126467344, -0.384886687, -0.77316704, 0.332620843, -1.647107979, 0.34094433, 1.556280176, -0.083923434, 0.432382122, 0.433788501, -0.316962568, -2.406639498, -0.085522764, -0.186979043, 0.278039403, 0.126011001, 1.468510171, 0.571926828, -0.058868504, -1.315359829, 0.035885671, 0.32233104, -0.459557597, 0.469290603, 1.137721779, 2.179248884, 1.498077933, 1.393277066, 0.273138359, 0.410165892, -0.365423104, 1.020732945, -0.686632488, 0.701253837, -1.19700993, -0.250149509, 0.416293419, 0.061132087, 0.831890244, -0.839656801, -0.475777277, -0.682231088, 0.105960557, 0.312311005, -0.512195559, -0.280520721, -0.746578294, 1.678569634, -0.505597422, 0.680584223, 0.14884956, 0.038804831, 0.457528379, -2.254511154, -0.098324035, -0.087044587, 0.857642068, -0.074809563, -1.606020686, 0.916085602, 0.721052954, -0.299804932, 0.969867318, 1.357378844, 0.219138391, 0.459802779, -0.382949053, 0.156849559, -2.080530596, -0.904422482, -0.505037665, -0.663657581, -0.119361591, 0.3754383, 0.95891091, 1.209612682, 0.590890203, 1.097682398, -0.823872618, 0.010057052, -1.225818035, -0.817688854, -1.378813822, 0.199843251, 0.673475559, 0.042591797, -0.446316491, 0.544995381, -2.337249932, -0.217460085, 0.547698366, -1.027430897, 1.937712788, 1.147399294, 1.262548634, -0.466981718, 1.921706631, 0.893803596, 0.636545065, -1.160241934, 1.453718881, -0.787537426, -0.660558523, 0.566338958, -1.412135759, -0.910148719, -0.734923863, -0.594608853, 2.267016554, -0.320478167, -2.143256633, -0.382620474, -1.907681047, -0.097689199, 0.623604101, -0.824959314, -0.824257624, 0.430352141, -0.420152784, 0.698001374, 0.35434313, -0.486679341, -0.278826248, -1.604419186, 1.382443656, 0.798609757, -0.164586565, 2.193192195, 0.704808497, 0.737897087, -0.622126937, -0.256323853, 0.30546939, 0.750274258, 0.815253016, -1.428402249, -0.47757267, 0.659788015, -1.11378197, -0.902233795, 0.255005482, -1.874420792, -0.964488721, 0.317077431, -0.043456686, -1.246411895, -0.784492981, -0.157103248, -0.518998751, -0.898316996, 0.160897403, -0.909944206, -0.463023208, 0.40005037, 0.273804892};
		    System.out.println(new Histogram(test).getSummary());
	    } catch (Exception e) {
		    e.printStackTrace();
	    }
    }
}
