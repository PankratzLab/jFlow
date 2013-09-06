package cnv.filesys;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;

import stats.Correlation;
import common.Array;
import common.DoubleVector;

public class MarkerData implements Serializable {
	public static final long serialVersionUID = 1L;
	public static final String MARKER_DATA_FILE_EXTENSION = ".mdRAF";
	public static final String[][] TYPES = { {"X Raw", "Y Raw"}, {"X", "Y"}, {"Theta", "R"}, {"B Allele Freq", "Log R Ratio"}};
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
	private String[] alleleMappings;

	public MarkerData(String markerName, byte chr, int position, long fingerprint, float[] gcs, float[] xRaws, float[] yRaws, float[] xs, float[] ys, float[] thetas, float[] rs, float[] bafs, float[] lrrs, byte[] abGenotypes, String[] alleleMappings) {
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
		this.alleleMappings = alleleMappings;
	}

	public float[][] getDatapoints(int type) {
		switch (type) {
		case 0:
			return new float[][] {xRaws, yRaws};
		case 1:
			return new float[][] {xs, ys};
		case 2:
//			return new float[][] {thetas, rs};
			return new float[][] {getThetas(), getRs()};
		case 3:
			return new float[][] {bafs, lrrs};
		default:
			System.err.println("Error - invalid plot type");
			return null;
		}
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

	public byte[] getAB_Genotypes() {
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

	public boolean[] getHighlightStatus(ClusterFilter clusterFilter) {
		byte[] original;
		boolean[] result;
		float[] realX;
		float[] realY;
		
		original = getAB_Genotypes(); 
		result = new boolean[original.length];
		for (int i=0; i<original.length; i++) {
			result[i] = false;
		}
		result=new boolean[original.length];
			// get appropriate data (X/Y Theta/R LRR/BAF)
			switch(clusterFilter.getPlotType()) {
			case 0:
				realX = getX_Raws();
				realY = getY_Raws();
				break;
			case 1:
				realX = getXs();
				realY = getYs();
				break;
			case 2:
				realX = getThetas();
				realY = getRs();
				break;
			case 3:
				realX = getBAFs();
				realY = getLRRs();
				break;
			default:
				realX = getXs();
				realY = getYs();
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
	
	public byte[] getAbGenotypesAfterFilters(ClusterFilterCollection clusterFilterCollection, String markerName, float gcThreshold) {
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
			return getAB_Genotypes();
		}
		
		
		original = getAB_Genotypes();
		result = new byte[original.length];
		if (gcThreshold > 1 || gcThreshold < 0) {
			System.err.println("Error - Invalid GC threshold: " + gcThreshold + ", expecting a decimal number between 0 and 1. Use 0 to include everything.");
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
				realX = getX_Raws();
				realY = getY_Raws();
				break;
			case 1:
				realX = getXs();
				realY = getYs();
				break;
			case 2:
				realX = getThetas();
				realY = getRs();
				break;
			case 3:
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
		}
		return result;		
	}

	public String[] getAlleleMappings() {
		return alleleMappings;
	}
	
	public double[] compareLRRs(float[][] centroids) {
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
					System.err.println("Error - compLRR is invalid ("+compLRRs[count]+") where oriLRR is not ("+lrrs[count]+")");
				} else {
					error += Math.abs(compLRRs[count]-originalLRRs[count]);
					count++;
					if (Double.isNaN(error) || Double.isInfinite(error)) {
						System.err.println("Started with index "+i+", compLRR of '"+compLRRs[count]+"', and oriLRR of '"+originalLRRs[count]+"'");
						return new double[] {-999,-999};
					}
				}
			}
        }
		
//		System.out.println("error="+error+" count="+count);
		return new double[] {Correlation.Pearson(Array.subArray(originalLRRs, 0, count), Array.subArray(compLRRs, 0, count))[0], error/count};
	}
	

	public void dump(String filename, String[] samples) {
		dump(filename, samples, false);
	}

	public void dump(String filename, String[] samples, boolean includeMarkerName) {
		PrintWriter writer;

		if ((xs != null && samples != null && samples.length!=xs.length) || (lrrs != null && samples != null && samples.length!=lrrs.length)) {
			System.err.println("Error - Number of samples (n="+samples.length+") does not match up with the number of LRRs/BAFs/etc (n="+lrrs.length+")");
			return;
        }

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
							+ (alleleMappings==null? "" : "\tForward_Genotypes"));
        	dump(writer, samples, includeMarkerName);
        	writer.close();
        } catch (Exception e) {
        	System.err.println("Error writing "+filename);
        	e.printStackTrace();
        }
	}

	public void dump(PrintWriter writer, String[] samples, boolean includeMarkerName) {
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
	        				+ (alleleMappings != null? "\t" + alleleMappings[i] : ""));
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

	public void setAbGenotype(byte abGenotype, int i) {
		abGenotypes[i] = abGenotype;
	}

	public void setAlleleMapping(String alleleMapping, int i) {
		alleleMappings[i] = alleleMapping;
	}

//	public byte[] detectClusters() {
//		return detectClusters(0.01, 3);
//	}

//	public byte[] detectClusters(double epsilon, int minPoints) {
//		byte currentClusterLabel;
//		byte[] clusterLabels;
//		IntVector neighborPoints;
//		IntVector neighborPointsExt;
//		boolean exit;
//		ArrayList<ArrayList<Integer>> clusters;
//		
//		currentClusterLabel = 1;
//		clusterLabels = new byte[xs.length];
//		for (int i=0; i<xs.length; i++) {
//			if (clusterLabels[i]!=0) {
//				neighborPoints = new IntVector();
//				for (int j=0; j<xs.length; j++) {
//					if (Distance.euclidean(new double[] {xs[i], ys[i]}, new double[] {xs[j], ys[j]}) <= epsilon) {
//						neighborPoints.add(j);
//					}
//				}
//	
//				if (neighborPoints.size() < minPoints) {
//					//TODO border point vs. noise point
//				} else {
//					clusterLabels[i]=currentClusterLabel;
//					exit = false;
//					while (!exit) {
//						neighborPointsExt = new IntVector();
//						for (int j=i; j<neighborPoints.size(); j++) {
//							if (clusterLabels[j]!=0) {
//								clusterLabels[j]=currentClusterLabel;
//								for (int k=0; k<xs.length; k++) {
//									if (!neighborPoints.contains(k) && !neighborPointsExt.contains(k) && Distance.euclidean(new double[] {xs[j], ys[j]}, new double[] {xs[k], ys[k]}) <= epsilon) {
//										neighborPointsExt.add(k);
//										clusterLabels[k]=currentClusterLabel;
//									}
//								}
//							}
//							if (neighborPointsExt.size()>0) {
//								neighborPoints=neighborPointsExt;
//								exit = false;
//							} else {
//								exit = true;
//							}
//						}
//					}
//				}
//				currentClusterLabel ++;
//			}
//		}
//
//		return clusterLabels;
//	}

//	public byte[] detectClusters(double epsilon, int minPoints) {
//		byte currentClusterLabel;
//		byte[] clusterLabels;
//		IntVector neighborPoints;
//		IntVector neighborPointsExt;
//		boolean exit;
//		
//		currentClusterLabel = 1;
//		clusterLabels = new byte[xs.length];
//		for (int i=0; i<xs.length; i++) {
//			if (clusterLabels[i]!=0) {
//				neighborPoints = new IntVector();
//				for (int j=0; j<xs.length; j++) {
//					if (Distance.euclidean(new double[] {xs[i], ys[i]}, new double[] {xs[j], ys[j]}) <= epsilon) {
//						neighborPoints.add(j);
//					}
//				}
//	
//				if (neighborPoints.size() < minPoints) {
//					//TODO border point vs. noise point
//				} else {
//					clusterLabels[i]=currentClusterLabel;
//					exit = false;
//					while (!exit) {
//						neighborPointsExt = new IntVector();
//						for (int j=i; j<neighborPoints.size(); j++) {
//							if (clusterLabels[j]!=0) {
//								clusterLabels[j]=currentClusterLabel;
//								for (int k=0; k<xs.length; k++) {
//									if (!neighborPoints.contains(k) && !neighborPointsExt.contains(k) && Distance.euclidean(new double[] {xs[j], ys[j]}, new double[] {xs[k], ys[k]}) <= epsilon) {
//										neighborPointsExt.add(k);
//										clusterLabels[k]=currentClusterLabel;
//									}
//								}
//							}
//							if (neighborPointsExt.size()>0) {
//								neighborPoints=neighborPointsExt;
//								exit = false;
//							} else {
//								exit = true;
//							}
//						}
//					}
//				}
//				currentClusterLabel ++;
//			}
//		}
//
//		return clusterLabels;
//	}

}
