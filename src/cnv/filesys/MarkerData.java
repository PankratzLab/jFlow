package cnv.filesys;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;

import cnv.analysis.CentroidCompute;
import cnv.analysis.pca.PrincipalComponentsIntensity;
import cnv.analysis.pca.PrincipalComponentsResiduals;
import cnv.manage.MarkerDataLoader;
import cnv.var.SampleData;
import stats.Correlation;
import common.AlleleFreq;
import common.Array;
import common.DoubleVector;
import common.Logger;

public class MarkerData implements Serializable {
	public static final long serialVersionUID = 1L;
	public static final String MARKER_DATA_FILE_EXTENSION = ".mdRAF";
	//Use for more alternate correction types
	//public static final String[][] TYPES = { { "X Raw", "Y Raw" }, { "X", "Y" }, { "Theta", "R" }, { "B Allele Freq", "Log R Ratio" }, { "BAF2", "LRR2" }, { "missThetaX", "missThetaY" },{ "twostageX", "twostageY" } ,{ "NstageX", "NstageY" },{ "N+1stageXResidual", "N+1stageYResidual" }};

	// TODO remove X Raw / Y Raw from the entire project
	
	private String markerName;
	private byte chr;
	private int position; 
	private long fingerprint;
	private float[] gcs;
	private float[] xRaws;
	private float[] yRaws;
	private float[] xs;
	private float[] ys;
	private float[] thetas;
	private float[] rs;
	private float[] lrrs;
	private float[] bafs;
	private byte[] abGenotypes;
	private byte[] forwardGenotypes;

	public MarkerData(String markerName, byte chr, int position, long fingerprint, float[] gcs, float[] xRaws, float[] yRaws, float[] xs, float[] ys, float[] thetas, float[] rs, float[] bafs, float[] lrrs, byte[] abGenotypes, byte[] forwardGenotypes) {
		this.markerName = markerName;
		this.chr = chr;
		this.position = position;
		this.fingerprint = fingerprint;
		this.gcs = gcs;
		this.xRaws = xRaws;
		this.yRaws = yRaws;
		this.xs = xs;
		this.ys = ys;
		this.thetas = thetas;
		this.rs = rs;
		this.bafs = bafs;
		this.lrrs = lrrs;
		this.abGenotypes = abGenotypes;
		this.forwardGenotypes = forwardGenotypes;
	}

	public float[][] getDatapoints(int type) {
		return getDatapoints(type, null, null, false, 1, 0, null, true, null, 0, 0, 0, 0, 1, false, new Logger());
	}

	public float[][] getDatapoints(int type, int[] sampleSex, boolean[] samplesToUse, boolean intensityOnly, double missingnessThreshold, double confThreshold, ClusterFilterCollection clusterFilterCollection, boolean medianCenter, PrincipalComponentsResiduals pcResids, int numComponents, int nstage, double residStandardDeviationFilter, double correctionRatio, int numThreads, boolean correctedData, Logger log) {
		switch (type) {
//		case 0:
//			return new float[][] { xRaws, yRaws };
		case 0:
			if (correctedData) {
				return getCorrectedIntesity(sampleSex, samplesToUse, missingnessThreshold, confThreshold, clusterFilterCollection, medianCenter, pcResids, numComponents, 2, nstage, residStandardDeviationFilter, correctionRatio, PrincipalComponentsIntensity.XY_RETURN, numThreads, log);
			} else {
				return new float[][] { xs, ys };
			}
		case 1:
			// return new float[][] {thetas, rs};
			if (correctedData) {
				return getCorrectedIntesity(sampleSex, samplesToUse, missingnessThreshold, confThreshold, clusterFilterCollection, medianCenter, pcResids, numComponents, 2, nstage, residStandardDeviationFilter, correctionRatio, PrincipalComponentsIntensity.THETA_R_RETURN, numThreads, log);
			} else {
				return new float[][] { getThetas(), getRs() };
			}
		case 2:
			if (correctedData) {
				return getCorrectedIntesity(sampleSex, samplesToUse, missingnessThreshold, confThreshold, clusterFilterCollection, medianCenter, pcResids, numComponents, 2, nstage, residStandardDeviationFilter, correctionRatio, PrincipalComponentsIntensity.BAF_LRR_RETURN, numThreads, log);
			} else {
				return new float[][] { bafs, lrrs };
			}
//		case 3:
//			return getRecomputedLRR_BAF(sampleSex, samplesToUse, intensityOnly, missingnessThreshold, confThreshold, clusterFilterCollection, medianCenter, false, log);
//			// case 5:
//			// return getCorrectedIntesity(sampleSex, missingnessThreshold, confThreshold, clusterFilterCollection, medianCenter, pcResids, numComponents, 0, 0, 0, false, log);
//			// case 6:
//			// return getCorrectedIntesity(sampleSex, missingnessThreshold, confThreshold, clusterFilterCollection, medianCenter, pcResids, numComponents, 1, 0, 0, false, log);
//			// case 7:
//			// return getCorrectedIntesity(sampleSex, missingnessThreshold, confThreshold, clusterFilterCollection, medianCenter, pcResids, numComponents, 2, nstage, 0, false, log);
//		//new case 5
//		case 4:
//			return getCorrectedIntesity(sampleSex, missingnessThreshold, confThreshold, clusterFilterCollection, medianCenter, pcResids, numComponents, 2, nstage, residStandardDeviationFilter, false,numThreads, log);
//		case 2:
//			return getCorrectedIntesity(sampleSex, missingnessThreshold, confThreshold, clusterFilterCollection, medianCenter, pcResids, numComponents, 2, nstage, residStandardDeviationFilter, true, numThreads,log);

		default:
			log.reportError("Error - invalid plot type");
			return null;
		}
	}

	/**
	 * Warning - the behavior of this method will likely be volatile for some time as the correction methods are improved
	 */
	public float[][] getCorrectedIntesity(int[] sampleSex, boolean[] samplesToUse, double missingnessThreshold, double confThreshold, ClusterFilterCollection clusterFilterCollection, boolean medianCenter, PrincipalComponentsResiduals pcResids, int numComponents, int correctionType, int nStage, double residStandardDeviationFilter,double correctionRatio, String typeToReturn, int numThreads, Logger log) {
		if (pcResids == null || numComponents == 0) {
			if (pcResids == null && numComponents > 0) {
				log.report("Info - an intensity PCA file was not found as specified by INTENSITY_PC_FILENAME property");
			}
			if (typeToReturn.equals(PrincipalComponentsIntensity.BAF_LRR_RETURN)) {
				// TODO, currently getRecomputedLRR_BAF may differ from PrincipalComponentsIntensity due to the intensity only flag...at false
				return getRecomputedLRR_BAF(sampleSex, samplesToUse, false, missingnessThreshold, confThreshold, clusterFilterCollection, medianCenter, false, log);
			}
			if (typeToReturn.equals(PrincipalComponentsIntensity.XY_RETURN)) {
				return new float[][] { xs, ys };
			}
			if (typeToReturn.equals(PrincipalComponentsIntensity.THETA_R_RETURN)) {
				return new float[][] { getThetas(), getRs() };
			} else {
				log.reportError("Error - invalid correction type");
				return null;
			}
		} else {
			numThreads = Math.min(numThreads, 6);// currently can only utilize 6
			PrincipalComponentsIntensity pcIntensity = new PrincipalComponentsIntensity(pcResids, this, true, sampleSex, pcResids.getProj().getSamplesToInclude(null, false), missingnessThreshold, confThreshold, clusterFilterCollection, medianCenter, (numComponents > 130 ? true : false), correctionType, nStage, residStandardDeviationFilter, correctionRatio, numThreads, false, null);
			pcIntensity.correctXYAt(numComponents);
			// This will display the genotypes after correction in scatter plot for testing, note that you will lose the original
			// setAbGenotypes(pcIntensity.getCentroidCompute().getClustGenotypes());
			return pcIntensity.getCorrectedIntensity(typeToReturn, true);
		}
	}

	/**
	 * Warning - the behavior of this method will likely be volatile for some time as the correction methods are improved
	 */
	public byte[] getAlternateGenotypes(int[] sampleSex, boolean[] samplesToUse, double missingnessThreshold, double confThreshold, ClusterFilterCollection clusterFilterCollection, boolean medianCenter, PrincipalComponentsResiduals pcResids, int numComponents, int correctionType, int nStage, double residStandardDeviationFilter, double correctionRatio, boolean correctLRR, int numThreads, Logger log) {
		if (pcResids == null || numComponents == 0) {
			if (pcResids == null && numComponents > 0) {
				log.report("Info - an intensity PCA file was not found as specified by INTENSITY_PC_FILENAME property");
			}
			return abGenotypes;
		} else {
			PrincipalComponentsIntensity pcIntensity = new PrincipalComponentsIntensity(pcResids, this, true, sampleSex, samplesToUse, missingnessThreshold, confThreshold, clusterFilterCollection, medianCenter, (numComponents > 130 ? true : false), correctionType, nStage, residStandardDeviationFilter, correctionRatio, numThreads, false, null);
			pcIntensity.correctXYAt(numComponents);
			return pcIntensity.getCentroidCompute().getClustGenotypes();
		}
	}

	public void setAbGenotypes(byte[] abGenotypes) {
		this.abGenotypes = abGenotypes;
	}

	public String getMarkerName() {
		return markerName;
	}

	public byte getChr() {
		return chr;
	}

	public int getPosition() {
		return position;
	}
	
	public long getFingerprint() {
		return fingerprint;
	}

	public float[] getGCs() {
		return gcs;
	}

	public float[] getX_Raws() {
		return xRaws;
	}

	public float[] getY_Raws() {
		return yRaws;
	}

	public float[] getXs() {
		return xs;
	}

	public float[] getYs() {
		return ys;
	}

	public float[] getThetas() {
		if (thetas == null) {
			thetas = new float[xs.length];
			for (int i = 0; i<xs.length; i++) {
				thetas[i] = Centroids.calcTheta(xs[i], ys[i]);
            }
		}
		return thetas;
	}

	public float[] getRs() {
		if (rs == null) {
			rs = new float[xs.length];
			for (int i = 0; i<xs.length; i++) {
				rs[i] = Centroids.calcR(xs[i], ys[i]);
            }
		}
		return rs;
	}

	public void recompute(float[][] centroids) {
		float[] thetas, rs;

		thetas = getThetas();
		rs = getRs();
		
		bafs = new float[xs.length];
		lrrs = new float[xs.length];
		for (int i = 0; i<xs.length; i++) {
			bafs[i] = Centroids.calcBAF(thetas[i], centroids);
			lrrs[i] = Centroids.calcLRR(thetas[i], rs[i], centroids);
        }
	}

	public float[] getBAFs() {
		return bafs;
	}

	public float[] getLRRs() {
		return lrrs;
	}

	public byte[] getAbGenotypes() {
		// TODO if null then compute from cluster centers
		return abGenotypes;
	}

	public int[] getAB_GenotypeCount() {
		int[] result;
		result = new int[3];
		for (int i=0; i<abGenotypes.length; i++) {
			if(abGenotypes[i]!=-1) {
				result[abGenotypes[i]]++;
			}
		}
		return result;
	}
	
	/**
	 * See {@link CentroidCompute#Centroid} for further usage
	 * 
	 * @param LRRonly
	 *            only return the lrr values. Warning the float[][] returned will have bafs=new float[0]. To get the LRRs, use getRecomputedLRR_BAF(...arguments...)[1];
	 * @return float[][] organized as float[0] = new BAFs, float[1]= new LRRs
	 */
	public float[][] getRecomputedLRR_BAF(int[] sampleSex, boolean[] samplesToUse, boolean intensityOnly, double missingnessThreshold, double confThreshold, ClusterFilterCollection clusterFilterCollection, boolean medianCenter, boolean LRRonly, Logger log) {
		CentroidCompute cent = getCentroid(sampleSex, samplesToUse, intensityOnly, missingnessThreshold, confThreshold, clusterFilterCollection, medianCenter, log);
		cent.computeCentroid();
		if (cent.failed()) {
			return new float[][] { bafs, lrrs };
		} else {
			return new float[][] { (LRRonly ? new float[0] : cent.getRecomputedBAF()), cent.getRecomputedLRR() };
		}
	}

	/**
	 * See {@link CentroidCompute#Centroid} for further usage
	 */
	public CentroidCompute getCentroid(int[] sampleSex, boolean[] samplesToUse, boolean intensityOnly, double missingnessThreshold, double confThreshold, ClusterFilterCollection clusterFilterCollection, boolean medianCenter, Logger log) {
		return new CentroidCompute(this, sampleSex, samplesToUse, intensityOnly, missingnessThreshold, confThreshold, clusterFilterCollection, medianCenter, log);
	}

	public boolean[] getHighlightStatus(ClusterFilter clusterFilter) {
		byte[] original;
		boolean[] result;
		float[] realX;
		float[] realY;
		
		original = getAbGenotypes(); 
		result = new boolean[original.length];
		for (int i=0; i<original.length; i++) {
			result[i] = false;
		}
		result=new boolean[original.length];
			// get appropriate data (X/Y Theta/R LRR/BAF)
			switch(clusterFilter.getPlotType()) {
			case 0:
				realX = getXs();
				realY = getYs();
				break;
			case 1:
				realX = getThetas();
				realY = getRs();
				break;
			case 2:
				realX = getBAFs();
				realY = getLRRs();
				break;
			default:
				realX = getXs();
				realY = getYs();
//			case 0:
//				realX = getX_Raws();
//				realY = getY_Raws();
//				break;
//			case 1:
//				realX = getXs();
//				realY = getYs();
//				break;
//			case 2:
//				realX = getThetas();
//				realY = getRs();
//				break;
//			case 3:
//				realX = getBAFs();
//				realY = getLRRs();
//				break;
//			default:
//				realX = getXs();
//				realY = getYs();
			}
			// iterate through all samples
			for (int j=0; j<result.length; j++) {
				if (realX[j]>=clusterFilter.getXMin()
						&& realY[j]>=clusterFilter.getYMin()
						&& realX[j]<=clusterFilter.getXMax()
						&& realY[j]<=clusterFilter.getYMax()) {
					result[j]=true;
				}
			}
		return result;		
		
	}
	
	public byte[] getAbGenotypesAfterFilters(ClusterFilterCollection clusterFilterCollection, String markerName, float gcThreshold, Logger log) {
		byte[] result, original;
		float[] realX;
		float[] realY;
		float clusterFilterXMin;
		float clusterFilterYMin;
		float clusterFilterXMax;
		float clusterFilterYMax;
		ArrayList<ClusterFilter> clusterFilters;
		ClusterFilter clusterFilter;
		int counter;
		
		if (clusterFilterCollection == null) {
			return getAbGenotypes();
		}
		
		
		original = getAbGenotypes();
		result = new byte[original.length];
		if (gcThreshold > 1 || gcThreshold < 0) {
			log.reportError("Error - Invalid GC threshold: " + gcThreshold + ", expecting a decimal number between 0 and 1. Use 0 to include everything.");
			return null;
		} else if (gcThreshold == 0 || gcs == null) {
			for (int i=0; i<original.length; i++) {
				result[i] = original[i];
			}
		} else {
			for (int i=0; i<original.length; i++) {
				if (gcs[i] < gcThreshold) {
					result[i]=(byte)-1;
				} else {
					result[i] = original[i];
				}
			}
		}
		
		clusterFilters = clusterFilterCollection.getClusterFilters(markerName);
		for (int i=0; clusterFilters != null && i < clusterFilters.size(); i++) {
			// get appropriate data (X/Y Theta/R LRR/BAF)
			clusterFilter = clusterFilters.get(i);
			switch(clusterFilter.getPlotType()) {
			case 0:
//				realX = getX_Raws();
//				realY = getY_Raws();
//				break;
//			case 1:
				realX = getXs();
				realY = getYs();
				break;
			case 1:
				realX = getThetas();
				realY = getRs();
				break;
			case 2:
				realX = getBAFs();
				realY = getLRRs();
				break;
			default:
				realX = getXs();
				realY = getYs();
			}

			clusterFilterXMin = clusterFilter.getXMin();
			clusterFilterYMin = clusterFilter.getYMin();
			clusterFilterXMax = clusterFilter.getXMax();
			clusterFilterYMax = clusterFilter.getYMax();

			try {
				counter = 0;
				for (int j=0; j<result.length; j++) {
					if (realX[j] >= clusterFilterXMin && realY[j] >= clusterFilterYMin && realX[j] <= clusterFilterXMax && realY[j] <= clusterFilterYMax) {
						result[j] = clusterFilter.getCluterGenotype();
						counter ++;
					}
				}
				if (counter == 0) {
					clusterFilterCollection.deleteClusterFilter(markerName, (byte) i);
					i--;
				}
			} catch (Exception e) {
				if (realX == null || realY == null) {
					log.reportError("Error - values for marker '"+markerName+"' were null; now returning null, which will probably create a new null pointer exception ");
				}
				log.reportException(e);
				return null;
			}
		}
		return result;		
	}

	public byte[] getForwardGenotypes() {
		return forwardGenotypes;
	}
	
	public double[] compareLRRs(float[][] centroids, Logger log) {
		double[] originalLRRs, compLRRs;
		double error;
		int count;
		
		error = 0;
		count = 0;
		originalLRRs = new double[lrrs.length];
		compLRRs = new double[lrrs.length];
		for (int i = 0; i<lrrs.length; i++) {
			if (!Float.isNaN(lrrs[i])) {
				originalLRRs[count] = lrrs[i];
				compLRRs[count] = Centroids.calcLRR(Centroids.calcTheta(xs[i], ys[i]), Centroids.calcR(xs[i], ys[i]), centroids);
				if (Double.isNaN(compLRRs[count])) {
					log.reportError("Error - compLRR is invalid ("+compLRRs[count]+") where oriLRR is not ("+lrrs[count]+")");
				} else {
					error += Math.abs(compLRRs[count]-originalLRRs[count]);
					count++;
					if (Double.isNaN(error) || Double.isInfinite(error)) {
						log.reportError("Started with index "+i+", compLRR of '"+compLRRs[count]+"', and oriLRR of '"+originalLRRs[count]+"'");
						return new double[] {-999,-999};
					}
				}
			}
        }
		
//		System.out.println("error="+error+" count="+count);
		return new double[] {Correlation.Pearson(Array.subArray(originalLRRs, 0, count), Array.subArray(compLRRs, 0, count))[0], error/count};
	}

	public void dump(SampleData sampleData, String filename, String[] samples, boolean includeMarkerName, Logger log) {
		PrintWriter writer;
		boolean hasExcludedIndividuals;
		
		if ((xs != null && samples != null && samples.length!=xs.length) || (lrrs != null && samples != null && samples.length!=lrrs.length)) {
			log.reportError("Error - Number of samples (n="+samples.length+") does not match up with the number of LRRs/BAFs/etc (n="+lrrs.length+")");
			return;
        }
		
		thetas = getThetas();
		rs = getRs();
		
		hasExcludedIndividuals = sampleData != null && sampleData.hasExcludedIndividuals();

		try {
        	writer = new PrintWriter(new FileWriter(filename));
        	writer.println( (includeMarkerName? "Marker\t" : "")
        					+ (samples==null? "SampleIndex" : "SampleId")
							+ (gcs==null? "" : "\tGC")
							+ (xRaws==null? "" : "\tRaw X")
							+ (yRaws==null? "" : "\tRaw Y")
							+ (xs==null? "" : "\tX")
							+ (ys==null? "" : "\tY")
							+ (thetas==null? "" : "\tTheta")
							+ (rs==null? "" : "\tR")
							+ (lrrs==null? "" : "\tLRR")
							+ (bafs==null? "" : "\tBAF")
							+ (abGenotypes==null? "" : "\tAB_Genotypes")
							+ (forwardGenotypes==null? "" : "\tForward_Genotypes")
							+ (hasExcludedIndividuals? "\tExclude_Sample" : "")
							);
        	for (int i = 0; i<xs.length; i++) {
        		writer.println(   (includeMarkerName? markerName + "\t" : "") 
        						+ (samples !=null? samples[i] : i) 
    	        				+ (gcs != null? "\t" + gcs[i] : "")
    	        				+ (xRaws != null? "\t" + xRaws[i] : "")
    	        				+ (yRaws != null? "\t" + yRaws[i] : "")
    	        				+ (xs != null? "\t" + xs[i] : "")
    	        				+ (ys != null? "\t" + ys[i]: "")
    	        				+ (thetas != null? "\t" + thetas[i] : "")
    	        				+ (rs != null? "\t" + rs[i] : "")
    	        				+ (lrrs != null? "\t" + lrrs[i] : "")
    	        				+ (bafs != null? "\t" + bafs[i] : "")
    	        				+ (abGenotypes != null? "\t" + abGenotypes[i] : "")
    	        				+ (forwardGenotypes != null? "\t" + Sample.ALLELE_PAIRS[forwardGenotypes[i]] : "")
    	        				+ (hasExcludedIndividuals? "\t" + (sampleData.individualShouldBeExcluded(samples[i])?1:0) : "")
    	        				);
        	}
        	writer.close();
        } catch (Exception e) {
        	log.reportError("Error writing "+filename);
        	log.reportException(e);
        }
	}

	public int detectCNP() {
		return detectCNP(0.05, 0.15);
	}

	public int detectCNP(double proportionOfLastPeakRequiredForNewLocalMinima, double proportionOfGlobalMaxRequiredForLocalMaxima) {
		float[][] clusterCenters;
		DoubleVector x = new DoubleVector();
		DoubleVector y = new DoubleVector();
		float ab_bb_xMidpoint, aa_ab_yMidpoint; 
		
		clusterCenters = Centroids.computeClusterCenters(this, null, 0.50);
////		if (clusterCenters[0].length == 0 || clusterCenters[1].length == 0 || clusterCenters[2].length == 0) {
////			return -1;
////		}
//		if (clusterCenters[0]!=null && clusterCenters[0].length!=0 && clusterCenters[1]!=null && clusterCenters[1].length!=0 && (clusterCenters[2]==null || clusterCenters[2].length==0)) {
////			ab_bb_xMidpoint = 0;	//should this be maximum?
//			ab_bb_xMidpoint = 1;
//			aa_ab_yMidpoint = Array.mean(new float[] {clusterCenters[0][1], clusterCenters[1][1]});
//		} else if ((clusterCenters[0]==null || clusterCenters[0].length==0) && clusterCenters[1]!=null && clusterCenters[1].length!= 0 && clusterCenters[2]!=null && clusterCenters[2].length!=0) {
//			ab_bb_xMidpoint = Array.mean(new float[] {clusterCenters[1][0], clusterCenters[2][0]});
////			aa_ab_yMidpoint = 0;	//should this be maximum?
//			aa_ab_yMidpoint = 1;
//		} else {
//			ab_bb_xMidpoint = Array.mean(new float[] {clusterCenters[1][0], clusterCenters[2][0]});
//			aa_ab_yMidpoint = Array.mean(new float[] {clusterCenters[0][1], clusterCenters[1][1]});
//		}
		
		if (clusterCenters[0] == null) {
			clusterCenters[0] = new float[] {1, 0};
		}

		if (clusterCenters[1] == null) {
			clusterCenters[1] = new float[] {1, 1};
		}

		if (clusterCenters[2] == null) {
			clusterCenters[2] = new float[] {0, 1};
		}		
		
		ab_bb_xMidpoint = Array.mean(new float[] {clusterCenters[1][0], clusterCenters[2][0]});
		aa_ab_yMidpoint = Array.mean(new float[] {clusterCenters[0][1], clusterCenters[1][1]});
		

		x = new DoubleVector();
		y = new DoubleVector();

		for (int i=0; i<xs.length; i++) {
			if (xs[i] < ab_bb_xMidpoint) {
				x.add(xs[i]);
			}
			if (ys[i] < aa_ab_yMidpoint) {
				y.add(ys[i]);
			}
		}

//		if (Array.isBimodal(x.toArray()) || Array.isBimodal(y.toArray())) {
//			return 1;
//		} else {
//			return 0;
//		}
//		System.out.print(x.size()+"\t"+y.size()+"\t");

		return (x.size()>0?Array.getLocalModes(x.toArray(), proportionOfLastPeakRequiredForNewLocalMinima, proportionOfGlobalMaxRequiredForLocalMaxima).length:0) + (y.size()>0?Array.getLocalModes(y.toArray(), proportionOfLastPeakRequiredForNewLocalMinima, proportionOfGlobalMaxRequiredForLocalMaxima).length:0);
	}

	public void setGC(float gc, int i) {
		gcs[i] = gc;
	}

	public void setX(float x, int i) {
		xs[i] = x;
	}

	public void setY(float y, int i) {
		ys[i] = y;
	}

	public void setBaf(float baf, int i) {
		bafs[i] = baf;
	}

	public void setLrr(float lrr, int i) {
		lrrs[i] = lrr;
	}

	public char[] getAB_AlleleMappings() {
		byte[] abGenotypes, forwardGenotypes;
		char[] mapping;
		
		mapping = new char[] {'0','0'};
		abGenotypes = getAbGenotypes();
		forwardGenotypes = getForwardGenotypes();

		for (int i = 0; i < abGenotypes.length; i++) {
			if (mapping[0] == '0' && (abGenotypes[i] == 0 || abGenotypes[i] == 1)) {
				mapping[0] = Sample.ALLELE_PAIRS[forwardGenotypes[i]].charAt(0);
			} else if (mapping[1] == '0' && (abGenotypes[i] == 1 || abGenotypes[i] == 2)) {
				mapping[1] = Sample.ALLELE_PAIRS[forwardGenotypes[i]].charAt(1);
			} else if (mapping[0] != '0' && mapping[1] != '0'){
				return mapping;
			}
		}
		
		return mapping;
	}

	// samplesToBeUsed, sex, and clusterFilterCollection can be null
	// however if sex is null then chrX genotype counts will be inaccurate and show deviation from Hardy-Weinberg equilibrium
	public int[] getGenotypeCounts(boolean[] samplesToBeUsed, String[] sex, ClusterFilterCollection clusterFilterCollection, float gcThreshold, Logger log) {
		int[] genotypeCounts = new int[3];
		String sexSpecific;
		byte[] genoytpes;
		
		if (clusterFilterCollection == null) {
			genoytpes = getAbGenotypes();
		} else {
			genoytpes = getAbGenotypesAfterFilters(clusterFilterCollection, markerName, gcThreshold, log);
		}
		// 0=non-specific, 1=male, 2=female
		if (chr == 23) {
			sexSpecific = "2"; // this is not technically correct, but there is no place to add the A/null and B/null male genotypes in this class; will add a new getMAF for chrX
		} else if (chr == 24) {
			sexSpecific = "1";
		} else {
			sexSpecific = "0";
		}
		for (int i = 0; i < genoytpes.length; i++) {
			if (genoytpes[i] >= 0 && (samplesToBeUsed == null || samplesToBeUsed[i]) && (sex == null || sexSpecific.equals("0") || sex[i].equals(sexSpecific))) {
				genotypeCounts[genoytpes[i]]++;
			}
		}
		return genotypeCounts;
	}

	// samplesToBeUsed, sex, and clusterFilters can be null
	// however if sex is null then chrX frequencies will be slightly inaccurate
	public double getMAF(boolean[] samplesToBeUsed, String[] sex, ClusterFilterCollection clusterFilters, float gcThreshold, Logger log) {
		double freqB;
		
		freqB = getFrequencyOfB(samplesToBeUsed, sex, clusterFilters, gcThreshold, log);
		if (freqB > 0.5) {
			return 1-freqB;
		} else {
			return freqB;
		}
	}

	// samplesToBeUsed, sex, and clusterFilters can be null
	// however if sex is null then chrX frequencies will be slightly inaccurate
	public double getFrequencyOfB(boolean[] samplesToBeUsed, String[] sex, ClusterFilterCollection clusterFilterCollection, float gcThreshold, Logger log) {
		int[] alleleCounts = new int[2];
		byte[] genoytpes;

		if (chr == 23 && sex != null) {
			if (clusterFilterCollection == null) {
				genoytpes = getAbGenotypes();
			} else {
				genoytpes = getAbGenotypesAfterFilters(clusterFilterCollection, markerName, gcThreshold, log);
			}

			alleleCounts = new int[2];
			for (int i = 0; i < genoytpes.length; i++) {
				if (genoytpes[i] >= 0 && (samplesToBeUsed == null || samplesToBeUsed[i])) {
					if (sex[i].equals("2")) {
						switch (genoytpes[i]) {
						case 0:
							alleleCounts[0] += 2;
							break;
						case 1:
							alleleCounts[0]++;
							alleleCounts[1]++;
							break;
						case 2:
							alleleCounts[1] += 2;
							break;
						default:
						}
					} else if (sex[i].equals("1")) {
						switch (genoytpes[i]) {
						case 0:
							alleleCounts[0]++;
							break;
						case 2:
							alleleCounts[1]++;
							break;
						default:
						}
					}
				}
			}
			return (double)alleleCounts[1] / (double)(alleleCounts[0]+alleleCounts[1]);
			
		} else {
			return AlleleFreq.calcFrequency(getGenotypeCounts(samplesToBeUsed, sex, clusterFilterCollection, gcThreshold, log));
		}
	}
	
	public byte[] compressMarker() {
		byte[] result, temp;

		result = new byte[0];
		if (gcs!=null) {
			temp = float2Byte(gcs);
			System.arraycopy(temp, 0, result, result.length, temp.length);
		}
		if (xRaws!=null) {
			temp = float2Byte(xRaws);
			System.arraycopy(temp, 0, result, result.length, temp.length);
		}
		if (yRaws!=null) {
			temp = float2Byte(yRaws);
			System.arraycopy(temp, 0, result, result.length, temp.length);
		}
		if (xs!=null) {
			temp = float2Byte(xs);
			System.arraycopy(temp, 0, result, result.length, temp.length);
		}
		if (ys!=null) {
			temp = float2Byte(ys);
			System.arraycopy(temp, 0, result, result.length, temp.length);
		}
		if (thetas!=null) {
			temp = float2Byte(thetas);
			System.arraycopy(temp, 0, result, result.length, temp.length);
		}
		if (rs!=null) {
			temp = float2Byte(rs);
			System.arraycopy(temp, 0, result, result.length, temp.length);
		}
		if (lrrs!=null) {
			temp = float2Byte(lrrs);
			System.arraycopy(temp, 0, result, result.length, temp.length);
		}
		if (bafs!=null) {
			temp = float2Byte(bafs);
			System.arraycopy(temp, 0, result, result.length, temp.length);
		}
		if (abGenotypes!=null) {
			System.arraycopy(abGenotypes, 0, result, result.length, abGenotypes.length);
		}
		if (forwardGenotypes!=null) {
			System.arraycopy(forwardGenotypes, 0, result, result.length, forwardGenotypes.length);
		}

		return result;
	}

	public static final byte[] float2Byte(float[] inData) {
		int j=0;
		int length=inData.length;
		byte[] outData=new byte[length*4];
		for (int i=0;i<length;i++) {
	    	int data=Float.floatToIntBits(inData[i]);
	    	outData[j++]=(byte)(data>>>24);
	    	outData[j++]=(byte)(data>>>16);
	    	outData[j++]=(byte)(data>>>8);
	    	outData[j++]=(byte)(data>>>0);
	    }
		return outData;
	}
	
	public static MarkerData[] loadMarkerData(Project proj, String[] markers) {
		if (markers.length > 500) {
			proj.getLog().reportTimeWarning("This method is generally used for loading a small subset of markers in the same thread, currently loading " + markers.length + " markers");
		}
		MarkerData[] markerDatas = new MarkerData[markers.length];
		MarkerDataLoader markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markers);
		for (int i = 0; i < markerDatas.length; i++) {
			markerDatas[i] = markerDataLoader.requestMarkerData(i);
			markerDataLoader.releaseIndex(i);
		}
		return markerDatas;
	}
}
