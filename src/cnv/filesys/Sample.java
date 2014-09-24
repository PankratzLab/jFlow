package cnv.filesys;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Date;
import java.util.Enumeration;
import java.util.Hashtable;

import common.Array;
import common.DoubleVector;
import common.Elision;
import common.Files;
import common.HashVec;
import common.Logger;
import common.ext;

public class Sample implements Serializable {
	public static final long serialVersionUID = 1L;
	public static final String[][] DATA_FIELDS = {{"GC Score", "GCscore", "confidence"}, {"X Raw"}, {"Y Raw"}, {"X", "Xvalue", "Log Ratio", "intensity_1"}, {"Y", "Yvalue", "Strength", "intensity_2"}, {"Theta"}, {"R"}, {"B Allele Freq"}, {"Log R Ratio"}};
	public static final String[][] GENOTYPE_FIELDS = {{"Allele1 - Forward", "Allele1", "genotype1", "Allele1 - Top"}, {"Allele2 - Forward", "Allele2", "genotype2", "Allele2 - Top"}, {"Allele1 - AB"}, {"Allele2 - AB"}}; // , {"Forward Strand Base Calls"}, {"Call Codes"}
	public static final String[] ALL_STANDARD_GENOTYPE_FIELDS = {"Allele1 - AB", "Allele2 - AB", "Allele1 - Forward", "Allele2 - Forward", "Allele1 - Top", "Allele2 - Top", "Allele1 - Design", "Allele2 - Design"};
	public static final String[] ALLELE_PAIRS = {"--", "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT", "DD", "DI", "II", "ID"};
	public static final String[] ALT_NULL = {"-", "0"};
	public static final String[] ALT_NULLS = {"--", "00", "---", "NoCall", "NC"};
	public static final String[] AB_PAIRS = {"AA", "AB", "BB"};
	public static final String SAMPLE_DATA_FILE_EXTENSION = ".sampRAF";
//	public static final byte PARAMETER_SECTION_BYTES = 13;
	public static final byte PARAMETER_SECTION_BYTES = 17;
	public static final byte PARAMETER_SECTION_NUMMARKERS_LOCATION = 0;
	public static final byte PARAMETER_SECTION_NUMMARKERS_LENGTH = 4;
	public static final byte PARAMETER_SECTION_NULLSTAT_LOCATION = 4;
	public static final byte PARAMETER_SECTION_NULLSTAT_LENGTH = 1;
	public static final byte PARAMETER_SECTION_OUTLIERSECTIONLENGTH_LOCATION = 5;
	public static final byte PARAMETER_SECTION_OUTLIERSECTIONLENGTH_LENGTH = 4;
//	public static final byte PARAMETER_SECTION_FINGPRNT_LOC = 5;
	public static final byte PARAMETER_SECTION_FINGPRNT_LOCATION = 9;
	public static final byte PARAMETER_SECTION_FINGPRNT_LENGTH = 8;
	public static final int MAX_ROWS_PER_WRITE_OPERATION = 500;
	public static final byte NULLSTATUS_GC_LOCATION = 0;
	public static final byte NULLSTATUS_X_LOCATION = 1;
	public static final byte NULLSTATUS_Y_LOCATION = 2;
	public static final byte NULLSTATUS_BAF_LOCATION = 3;
	public static final byte NULLSTATUS_LRR_LOCATION = 4;
	public static final byte NULLSTATUS_ABGENOTYPE_LOCATION = 5;
	public static final byte NULLSTATUS_FOWARDGENOTYPE_LOCATION = 6;
	public static final byte NULLSTATUS_CAN_XY_BE_NEGATIVE_LOCATION = 7;

	private byte[] abGenotypes;
	private byte[] forwardGenotypes;
	private float[] xs;
	private float[] ys;
	private float[] thetas;
	private float[] rs;
	private float[] lrrs;
	private float[] bafs;
	private float[] gcs;
	private long fingerprint;
	private String sampleName;
	private byte nullStatus;
	private boolean canXYBeNegative;

	public Sample(String sampleName, long fingerprint, float[] gcs, float[] xs, float[] ys, float[] bafs, float[] lrrs, byte[] forwardGenotypes, byte[] abGenotypes, boolean canXYBeNegative) {
		this.sampleName = sampleName;
		this.fingerprint = fingerprint;
		this.gcs = gcs;
		this.xs = xs;
		this.ys = ys;
		this.bafs = bafs;
		this.lrrs = lrrs;
		this.forwardGenotypes = forwardGenotypes;
		this.abGenotypes = abGenotypes;
		this.canXYBeNegative = canXYBeNegative;
		updateNullStatus();
		// note for future expansion: 8th bit should indicate a second nullStatus byte will be used, to ensure backward compatibility
	}

	public Sample(String sampleName, long fingerprint, float[][] data, byte[][] genotypes, boolean canXYBeNegative) {
		this.sampleName = sampleName;
		this.fingerprint = fingerprint;
		this.gcs = data[0];
		this.xs = data[3];
		this.ys = data[4];
		this.thetas = null;
		this.rs = null;
		this.bafs = data[7];
		this.lrrs = data[8];
		this.forwardGenotypes = genotypes[0];
		this.abGenotypes = genotypes[1];
		this.canXYBeNegative = canXYBeNegative;
		updateNullStatus();
	}

	public void updateNullStatus() {
		nullStatus = 0;
		if (gcs==null) {
			nullStatus = (byte) (nullStatus | (1 << NULLSTATUS_GC_LOCATION));
		}
		if (xs==null) {
			nullStatus = (byte) (nullStatus | (1 << NULLSTATUS_X_LOCATION));
		}
		if (ys==null) {
			nullStatus = (byte) (nullStatus | (1 << NULLSTATUS_Y_LOCATION));
		}
		if (bafs==null) {
			nullStatus = (byte) (nullStatus | (1 << NULLSTATUS_BAF_LOCATION));
		}
		if (lrrs==null) {
			nullStatus = (byte) (nullStatus | (1 << NULLSTATUS_LRR_LOCATION));
		}
		if (abGenotypes==null) {
			nullStatus = (byte) (nullStatus | (1 << NULLSTATUS_ABGENOTYPE_LOCATION));
		}
		if (forwardGenotypes==null) {
			nullStatus = (byte) (nullStatus | (1 << NULLSTATUS_FOWARDGENOTYPE_LOCATION));
		}
		if (canXYBeNegative) {
			nullStatus = (byte) (nullStatus | (1 << NULLSTATUS_CAN_XY_BE_NEGATIVE_LOCATION));
		}
	}

	public static byte updateNullStatus(float[] gcs, float[] xs, float[] ys, float[] bafs, float[] lrrs, byte[] abGenotypes, byte[] forwardGenotypes, boolean canXYBeNegative) {
		byte nullStatus;

		nullStatus = 0;
		if (gcs==null) {
			nullStatus = (byte) (nullStatus | (1 << Sample.NULLSTATUS_GC_LOCATION));
		}
		if (xs==null) {
			nullStatus = (byte) (nullStatus | (1 << Sample.NULLSTATUS_X_LOCATION));
		}
		if (ys==null) {
			nullStatus = (byte) (nullStatus | (1 << Sample.NULLSTATUS_Y_LOCATION));
		}
		if (bafs==null) {
			nullStatus = (byte) (nullStatus | (1 << Sample.NULLSTATUS_BAF_LOCATION));
		}
		if (lrrs==null) {
			nullStatus = (byte) (nullStatus | (1 << Sample.NULLSTATUS_LRR_LOCATION));
		}
		if (abGenotypes==null) {
			nullStatus = (byte) (nullStatus | (1 << Sample.NULLSTATUS_ABGENOTYPE_LOCATION));
		}
		if (forwardGenotypes==null) {
			nullStatus = (byte) (nullStatus | (1 << Sample.NULLSTATUS_FOWARDGENOTYPE_LOCATION));
		}
		if (canXYBeNegative) {
			nullStatus = (byte) (nullStatus | (1 << Sample.NULLSTATUS_CAN_XY_BE_NEGATIVE_LOCATION));
		}

		return nullStatus;
	}

	public byte getNBytesPerSampleMarker() {
		byte nBytesPerSampleMarker;
		
		nBytesPerSampleMarker = (byte) Compression.BYTES_PER_SAMPLE_MARKER;
		if (gcs==null) {
			nBytesPerSampleMarker -= Compression.REDUCED_PRECISION_GCBAF_NUM_BYTES;
		}
		if (xs==null) {
			nBytesPerSampleMarker -= Compression.REDUCED_PRECISION_XY_NUM_BYTES;
		}
		if (ys==null) {
			nBytesPerSampleMarker -= Compression.REDUCED_PRECISION_XY_NUM_BYTES;
		}
		if (bafs==null) {
			nBytesPerSampleMarker -= Compression.REDUCED_PRECISION_GCBAF_NUM_BYTES;
		}
		if (lrrs==null) {
			nBytesPerSampleMarker -= Compression.REDUCED_PRECISION_LRR_NUM_BYTES;
		}
		if (abGenotypes==null && forwardGenotypes==null) {
			nBytesPerSampleMarker -= Compression.REDUCED_PRECISION_ABFORWARD_GENOTYPE_NUM_BYTES;
		}

		return nBytesPerSampleMarker;
	}

	public static byte getNBytesPerSampleMarker(byte nullStatus) {
		byte nBytesPerSampleMarker;
		
		nBytesPerSampleMarker = (byte) Compression.BYTES_PER_SAMPLE_MARKER;
		if (isGcNull(nullStatus)) {
			nBytesPerSampleMarker -= Compression.REDUCED_PRECISION_GCBAF_NUM_BYTES;
		}
		if (isXOrYNull(nullStatus)) {
			nBytesPerSampleMarker -= Compression.REDUCED_PRECISION_XY_NUM_BYTES;
		}
		if (((nullStatus >> NULLSTATUS_Y_LOCATION) & 0x01) == 1) {
			nBytesPerSampleMarker -= Compression.REDUCED_PRECISION_XY_NUM_BYTES;
		}
		if (isBafOrLrrNull(nullStatus)) {
			nBytesPerSampleMarker -= Compression.REDUCED_PRECISION_GCBAF_NUM_BYTES;
		}
		if (((nullStatus >> NULLSTATUS_LRR_LOCATION) & 0x01) == 1) {
			nBytesPerSampleMarker -= Compression.REDUCED_PRECISION_LRR_NUM_BYTES;
		}
		if (isAbOrForwardGenotypeNull(nullStatus)) {
			nBytesPerSampleMarker -= Compression.REDUCED_PRECISION_ABFORWARD_GENOTYPE_NUM_BYTES;
		}

		return nBytesPerSampleMarker;
	}

	public static boolean isNegativeXOrYAllowed (byte nullStatus) {
		return ((nullStatus >> NULLSTATUS_CAN_XY_BE_NEGATIVE_LOCATION) & 0x01) == 1;
	}

	public static boolean isGcNull(byte nullStatus) {
		return ((nullStatus >> NULLSTATUS_GC_LOCATION) & 0x01) == 1;
	}

	public static boolean isXOrYNull(byte nullStatus) {
		return (((nullStatus >> NULLSTATUS_X_LOCATION) & 0x01) == 1 || ((nullStatus >> NULLSTATUS_Y_LOCATION) & 0x01) == 1);
	}

	public static boolean isXNull(byte nullStatus) {
		return ((nullStatus >> NULLSTATUS_X_LOCATION) & 0x01) == 1;
	}

	public static boolean isYNull(byte nullStatus) {
		return ((nullStatus >> NULLSTATUS_Y_LOCATION) & 0x01) == 1;
	}

	public static boolean isBafOrLrrNull(byte nullStatus) {
		return (((nullStatus >> NULLSTATUS_BAF_LOCATION) & 0x01) == 1 || ((nullStatus >> NULLSTATUS_LRR_LOCATION) & 0x01) == 1);
	}

	public static boolean isBafNull(byte nullStatus) {
		return ((nullStatus >> NULLSTATUS_BAF_LOCATION) & 0x01) == 1;
	}

	public static boolean isLrrNull(byte nullStatus) {
		return ((nullStatus >> NULLSTATUS_LRR_LOCATION) & 0x01) == 1;
	}

	public static boolean isAbOrForwardGenotypeNull(byte nullStatus) {
		return (((nullStatus >> NULLSTATUS_ABGENOTYPE_LOCATION) & 0x01) == 1 || ((nullStatus >> NULLSTATUS_FOWARDGENOTYPE_LOCATION) & 0x01) == 1);
	}

	public float[][] getAllData() {
		return new float[][] {gcs, xs, ys, thetas, rs, bafs, lrrs};
	}

	public byte[][] getAllGenotypes() {
		return new byte[][] {forwardGenotypes, abGenotypes};
	}

	public String getSampleName() {
		return sampleName;
	}

	public long getFingerprint() {
		return fingerprint;
	}

	public float[] getGCs() {
		return gcs;
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

	public float[] getBAFs() {
		return bafs;
	}

	public float[] getBAFs(float[][][] centroids) {
		float[] thetas, bafs;

		thetas = getThetas();
		bafs = new float[xs.length];
		for (int i = 0; i<xs.length; i++) {
			bafs[i] = Centroids.calcBAF(thetas[i], centroids[i]);
        }

		return bafs;
	}

	public float[] getLRRs() {
		return lrrs;
	}

	public float[] getLRRs(float[][][] centroids) {
		float[] thetas, rs, lrrs;

		thetas = getThetas();
		rs = getRs();
		lrrs = new float[xs.length];
		for (int i = 0; i<xs.length; i++) {
			lrrs[i] = Centroids.calcLRR(thetas[i], rs[i], centroids[i]);
        }

		return lrrs;
	}

	public byte[] getForwardGenotypes() {
		return forwardGenotypes;
	}

	public byte[] getForwardGenotypes(float gcThreshold) {
		byte[] result = new byte[forwardGenotypes.length];
		
		for (int i = 0; i < result.length; i++) {
			if (gcs[i] <= gcThreshold) {
				result[i] = (byte)0;
			} else {
				result[i] = forwardGenotypes[i];
			}
		}

		return result;
	}
	
	public boolean getCanXYBeNegative() {
		return canXYBeNegative;
	}

	public byte[] getAB_GenotypesAfterFilters(String[] markerNames, ClusterFilterCollection clusterFilterCollection, float gcThreshold) {
		byte[] result;
		float realX;
		float realY;
		ClusterFilter clusterFilter;
		ArrayList<ClusterFilter> clusterFilterArray;

		getThetas();
		getRs();

		if (markerNames.length != abGenotypes.length) {
			System.err.println("Error - need to pass the full set of markerNames to getAB_GenotypesAfterFilters");
			return null;
		}
	
		result = new byte[abGenotypes.length];
		for (int i = 0; i < abGenotypes.length; i ++) {
			if (gcs[i]<gcThreshold) {
				result[i]=(byte)-1;
			} else {
				result[i]=abGenotypes[i];
			}

			clusterFilterArray = clusterFilterCollection.getClusterFilters(markerNames[i]);
			if (clusterFilterArray != null) {
				for (int j=0; j<clusterFilterArray.size(); j++) {
					clusterFilter = clusterFilterArray.get(j);
					switch(clusterFilter.getPlotType()) {
					case 0:
						realX = -9;
						realY = -9;
						System.err.println("Error - PlotType cannot be 0 for ClusterFilter #"+(j+1)+" for marker '"+markerNames[i]+"' as we've done away with raw Xs and Ys");
						break;
					case 1:
						realX = xs[i];
						realY = ys[i];
						break;
					case 2:
						realX = thetas[i];
						realY = rs[i];
						break;
					case 3:
						realX = bafs[i];
						realY = lrrs[i];
						break;
					default:
						System.err.println("Error - invalid PlotType in ClusterFilter #"+(j+1)+" for marker '" + i + "'");
						realX = -9;
						realY = -9;
					}
					if (realX>=clusterFilter.getXMin() && realY>=clusterFilter.getYMin() && realX<=clusterFilter.getXMax() && realY<=clusterFilter.getYMax()) {
						result[i]=clusterFilter.getCluterGenotype();
					}
				}
			}
		}
		return result;
	}
	
	public byte[] getAB_GenotypesAfterFilters(String[] targetMarkers, int[] indicesOfTargetMarkersInProj, ClusterFilterCollection clusterFilterCollection, float gcThreshold) {
		byte[] result;
		float realX;
		float realY;
		ClusterFilter clusterFilter;
		ArrayList<ClusterFilter> clusterFilterArray;

		if (indicesOfTargetMarkersInProj == null) {
			return getAB_GenotypesAfterFilters(targetMarkers, clusterFilterCollection, gcThreshold);
		} else {
			result = new byte[targetMarkers.length];
			for (int i=0; i < targetMarkers.length; i++) {
				if (gcs[indicesOfTargetMarkersInProj[i]] < gcThreshold) {
					result[i] = (byte) -1;
				} else {
					result[i] = abGenotypes[indicesOfTargetMarkersInProj[i]];
				}
	
				clusterFilterArray = clusterFilterCollection.getClusterFilters(targetMarkers[i]);
				if (clusterFilterArray != null) {
					for (int j=0; j < clusterFilterArray.size(); j++) {
						clusterFilter = clusterFilterArray.get(j);
						switch(clusterFilter.getPlotType()) {
						case 0:
							realX = -9;
							realY = -9;
							System.err.println("Error - PlotType cannot be 0 for ClusterFilter #" + (j + 1) + " for marker '" + targetMarkers[i] + "' as we've done away with raw Xs and Ys");
							break;
						case 1:
							realX = xs[indicesOfTargetMarkersInProj[i]];
							realY = ys[indicesOfTargetMarkersInProj[i]];
							break;
						case 2:
							realX = thetas[indicesOfTargetMarkersInProj[i]];
							realY = rs[indicesOfTargetMarkersInProj[i]];
							break;
						case 3:
							realX = bafs[indicesOfTargetMarkersInProj[i]];
							realY = lrrs[indicesOfTargetMarkersInProj[i]];
							break;
						default:
							System.err.println("Error - invalid PlotType in ClusterFilter #" + (j + 1) + " for marker '" + targetMarkers[i] + "'");
							realX = -9;
							realY = -9;
						}
						if (realX>=clusterFilter.getXMin() && realY>=clusterFilter.getYMin() && realX<=clusterFilter.getXMax() && realY<=clusterFilter.getYMax()) {
							result[i] = clusterFilter.getCluterGenotype();
						}
					}
				}
			}
			return result;
		}
	}
	
	public byte[] getAB_Genotypes() {
		return abGenotypes;
	}

	public byte[] getAB_Genotypes(int[] indicesOfTargetMarkersInProj) {
		byte[] result;
		
		if (indicesOfTargetMarkersInProj == null) {
			return abGenotypes;

		} else {
			result = new byte[indicesOfTargetMarkersInProj.length];
			for (int i = 0; i < result.length; i++) {
				result[i] = abGenotypes[indicesOfTargetMarkersInProj[i]];
			}
			return result;
		}
	}

	public void setAB_Genotypes(byte[] abGenotypes) {
		this.abGenotypes = abGenotypes;
		updateNullStatus();
	}

	public void writeToFile(String[] markerNames, String filename) {
		PrintWriter writer;

		if (markerNames.length!=lrrs.length) {
			System.err.println("Error - MarkerNames (n="+markerNames.length+") do not match up with the number of LRRs/BAFs/etc (n="+lrrs.length+")");
			System.exit(1);
		}
		try {
			writer = new PrintWriter(new FileWriter(filename));
			writer.println("SNP\tGC Score\tX\tY\tTheta\tR\tLRR\tBAF\tGenotypes");
			for (int i = 0; i<markerNames.length; i++) {
				writer.println(markerNames[i]+"\t"+gcs[i]+"\t"+xs[i]+"\t"+ys[i]+"\t"+thetas[i]+"\t"+rs[i]+"\t"+lrrs[i]+"\t"+bafs[i]+"\t"+ALLELE_PAIRS[forwardGenotypes[i]]);
				
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing "+filename);
			e.printStackTrace();
		}
	}

	public float[][][] loadCentroids(Project proj) {
		Centroids centroids;
		
		centroids = Centroids.load(proj.getFilename(Project.CUSTOM_CENTROIDS_FILENAME), proj.getJarStatus());
		if (centroids.getFingerprint() != getFingerprint()) {
			System.err.println("Error - mismatched fingerprint for "+sampleName);
		}
		return centroids.getCentroids();
	}
	
	public void compareCalculationsFile(Project proj, String[] markerNames, String filename) {
		PrintWriter writer;
		float[][][] centroids;
		float[] compBAFs, compLRRs;

		centroids = loadCentroids(proj);
		compBAFs = getBAFs(centroids);
		compLRRs = getLRRs(centroids);
		
		if (markerNames.length!=lrrs.length) {
			System.err.println("Error - MarkerNames (n="+markerNames.length+") do not match up with the number of LRRs/BAFs/etc (n="+lrrs.length+")");
			System.exit(1);
		}
		try {
			writer = new PrintWriter(new FileWriter(filename));
			writer.println("SNP\tX\tY\tTheta\tR\tcompTheta\tcompR\tLRR\tBAF\tcompLRR\tcompBAF");
			for (int i = 0; i<markerNames.length; i++) {
				writer.println(markerNames[i]+"\t"+xs[i]+"\t"+ys[i]+"\t"+thetas[i]+"\t"+rs[i]+"\t"+bafs[i]+"\t"+lrrs[i]+"\t"+compBAFs[i]+"\t"+compLRRs[i]);
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing "+filename);
			e.printStackTrace();
		}
	}

	public void compLRRs(Project proj) {
		float[] compLRRs, diffs;

		compLRRs = getLRRs(loadCentroids(proj));
		diffs = new float[lrrs.length];
		for (int i = 0; i<lrrs.length; i++) {
			diffs[i] = lrrs[i]-compLRRs[i];
		}
		Files.writeSerial(diffs, proj.getProjectDir()+"comps/"+sampleName+".comp");
	}

	public void serialize(String filename) {
		Files.writeSerial(this, filename);
	}

	/**
	 * Save the instance of Sample to hard drive in Random Access File format.
	 * @param filename
	 */
	public void saveToRandomAccessFile(String filename) {
		saveToRandomAccessFile(filename, null, null);
	}
		
	public void saveToRandomAccessFile(String filename, Hashtable<String, Float> allOutliers, String sampleName) {
		File fileTmp;
		RandomAccessFile rafFile;
		Hashtable<String, Float> outOfRangeValuesEachSample;
		byte[] outOfRangeValuesWriteBuffer;
		int bytesRemained;
		byte[] writeBuffer = null;
		byte[] parameters;
		byte[] temp;
		int writeBufferIndex;
		byte bytesPerSampleMarker;
		boolean isOutlier;

		fileTmp = new File(filename);
		if (new File(ext.parseDirectoryOfFile(filename)).getFreeSpace() <= ((long)xs.length * Compression.BYTES_PER_SAMPLE_MARKER)) {
			System.err.println("Not enough space (available: "+ext.prettyUpSize(new File(ext.parseDirectoryOfFile(filename)).getFreeSpace(), 1)+") for all the new data to be created (required: "+ext.prettyUpSize(((long)xs.length * Compression.BYTES_PER_SAMPLE_MARKER), 1)+").");
			return;
		}
		if (fileTmp.exists()) {
			fileTmp.delete();
		}

		bytesPerSampleMarker = getNBytesPerSampleMarker();
		bytesRemained = xs.length *  bytesPerSampleMarker;
		long time = new Date().getTime();
		outOfRangeValuesEachSample = new Hashtable<String, Float>();
		try {
			rafFile = new RandomAccessFile(filename, "rw");

//			rafFile.writeInt(xs.length);
//			rafFile.writeByte(nullStatus);
//			rafFile.writeLong(fingerprint);
			parameters = new byte[PARAMETER_SECTION_BYTES];
			temp = Compression.intToBytes(xs.length);
			for (int i=0; i<temp.length; i++) {
				parameters[PARAMETER_SECTION_NUMMARKERS_LOCATION + i] = temp[i];
			}

			parameters[PARAMETER_SECTION_NULLSTAT_LOCATION] = nullStatus;

			temp = Compression.longToBytes(fingerprint);
			for (int i=0; i<temp.length; i++) {
				parameters[PARAMETER_SECTION_FINGPRNT_LOCATION + i] = temp[i];
			}
			rafFile.write(parameters);

			writeBufferIndex = 0;
			for (int j = 0; j<xs.length; j++) {
				if (writeBufferIndex == 0) {
					writeBuffer = new byte[ Math.min(Integer.MAX_VALUE, bytesRemained) ];
				}
				if (gcs != null) {
					try {
						Compression.gcBafCompress(gcs[j], writeBuffer, writeBufferIndex);
					} catch (Elision e) {
						System.err.println("Error - problem saving sampleRAF file '"+filename+"' for sample "+sampleName+" since the marker in index "+j+" has a GC value of "+gcs[j]+"; in the past this has happened if the final report file was truncated and did not contain all of the markers");
					}
					writeBufferIndex += Compression.REDUCED_PRECISION_GCBAF_NUM_BYTES;
				}
				if (xs != null) {
					if (canXYBeNegative) {
						isOutlier = ! Compression.xyCompressAllowNegative(xs[j], writeBuffer, writeBufferIndex);
					} else {
						if (writeBufferIndex == 2726570) {
							System.out.println();
						}
						isOutlier = ! Compression.xyCompressPositiveOnly(xs[j], writeBuffer, writeBufferIndex);
					}
					if (isOutlier) {
						outOfRangeValuesEachSample.put(j + "\tx", xs[j]);
						if (allOutliers != null) {
							allOutliers.put(j + "\t" + sampleName + "\tx", xs[j]);
						}
					}
					writeBufferIndex += Compression.REDUCED_PRECISION_XY_NUM_BYTES;
				}
				if (ys != null) {
					if (canXYBeNegative) {
						isOutlier = ! Compression.xyCompressAllowNegative(ys[j], writeBuffer, writeBufferIndex);
					} else {
						isOutlier = ! Compression.xyCompressPositiveOnly(ys[j], writeBuffer, writeBufferIndex);
					}
					if (isOutlier) {
						outOfRangeValuesEachSample.put(j + "\ty", ys[j]);
						if (allOutliers != null) {
							allOutliers.put(j + "\t" + sampleName + "\ty", ys[j]);
						}
					}
					writeBufferIndex += Compression.REDUCED_PRECISION_XY_NUM_BYTES;
				}
				if (bafs != null) {
					Compression.gcBafCompress(bafs[j], writeBuffer, writeBufferIndex);
					writeBufferIndex += Compression.REDUCED_PRECISION_GCBAF_NUM_BYTES;
				}
				if (lrrs != null) {
					if (Compression.lrrCompress(lrrs[j], writeBuffer, writeBufferIndex) == -1) {
						outOfRangeValuesEachSample.put(j + "\tlrr", lrrs[j]);
						if (allOutliers != null) {
							allOutliers.put(j + "\t" + sampleName + "\tlrr", lrrs[j]);
						}
					}
					writeBufferIndex += Compression.REDUCED_PRECISION_LRR_NUM_BYTES;
				}
				if (abGenotypes != null || forwardGenotypes != null) {
					writeBuffer[writeBufferIndex] = Compression.genotypeCompress(abGenotypes == null?-1:abGenotypes[j], forwardGenotypes == null?0:forwardGenotypes[j]);
					writeBufferIndex += Compression.REDUCED_PRECISION_ABFORWARD_GENOTYPE_NUM_BYTES;
				}
				bytesRemained -= bytesPerSampleMarker;
				if ( writeBufferIndex>=writeBuffer.length ) {
					rafFile.write(writeBuffer);
					writeBufferIndex = 0;
				}
			}

			if (outOfRangeValuesEachSample!=null && outOfRangeValuesEachSample.size()>0) {
				outOfRangeValuesWriteBuffer = Compression.objToBytes(outOfRangeValuesEachSample);
				rafFile.write(outOfRangeValuesWriteBuffer);
				rafFile.seek(PARAMETER_SECTION_OUTLIERSECTIONLENGTH_LOCATION);
				rafFile.writeInt(outOfRangeValuesWriteBuffer.length);
			}

			rafFile.close();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (Elision e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		time = new Date().getTime() - time;
//		System.out.println("Random Access File. ---- Finished writing to all the samples in " + (time/60000) + " min " + ((time%60000)/1000) + " sec");
	}

	public static Sample loadFromSerialized(String filename, boolean jar) {
		return (Sample)Files.readSerial(filename, jar, true);
	}

	public static Sample loadFromRandomAccessFile(String filename, boolean jar) {
		return loadFromRandomAccessFile(filename, true, true, true, true, true, jar);
	}

	/**
	 * Load data from Random Access Files organized by samples.
	 */
	@SuppressWarnings("unchecked")
	public static Sample loadFromRandomAccessFile(String filename, boolean loadGC, boolean loadXY, boolean loadBAF, boolean loadLRR, boolean loadAbOrForwardGenotypes, boolean jar) {
		Sample result = null;
		int numMarkers;
		RandomAccessFile file;
		String sampleName;
		long fingerPrint;
        float[] gcs = null, xs = null, ys = null, lrrs = null, bafs = null;
        byte[] abGenotypes = null, fwdGenotypes = null, readBuffer;
        byte[] genoTypeTmp;
        byte[] temp;
        int index, indexStart;
        int numBytesOfOutOfRangeValues;
        Hashtable<String, Float> outOfRangeValues = null;
        int outlierSectionLocation;
        byte nullStatus;
        byte bytesPerSampleMarker;

        sampleName = ext.rootOf(filename);
		try {
			file = new RandomAccessFile(filename, "r");
			readBuffer = new byte[(int) file.length()];	//numMarkers * BYTES_PER_SAMPLE_MARKER
			file.read(readBuffer);
			file.close();

//			numMarkers = Compression.bytesToInt(new byte[]{readBuffer[0], readBuffer[1], readBuffer[2], readBuffer[3]});
			temp = new byte[PARAMETER_SECTION_NUMMARKERS_LENGTH];
			for (int i=0; i<temp.length; i++) {
				temp[i] = readBuffer[PARAMETER_SECTION_NUMMARKERS_LOCATION + i];
			}
			numMarkers = Compression.bytesToInt(temp);

//			nullStatus = readBuffer[4];
//			temp = new byte[PARAMETER_SECTION_NULLSTAT_LEN];
//			for (int i=0; i<temp.length; i++) {
//				temp[i] = readBuffer[PARAMETER_SECTION_NULLSTAT_LOC + i];
//			}
			nullStatus = readBuffer[PARAMETER_SECTION_NULLSTAT_LOCATION];

			temp = new byte[PARAMETER_SECTION_OUTLIERSECTIONLENGTH_LENGTH];
			for (int i=0; i<temp.length; i++) {
				temp[i] = readBuffer[PARAMETER_SECTION_OUTLIERSECTIONLENGTH_LOCATION + i];
			}
			numBytesOfOutOfRangeValues = Compression.bytesToInt(temp);

//			fingerPrint = Compression.bytesToLong(new byte[]{readBuffer[5], readBuffer[6], readBuffer[7], readBuffer[8], readBuffer[9], readBuffer[10], readBuffer[11], readBuffer[12]});
			temp = new byte[PARAMETER_SECTION_FINGPRNT_LENGTH];
			for (int i=0; i<temp.length; i++) {
				temp[i] = readBuffer[PARAMETER_SECTION_FINGPRNT_LOCATION + i];
			}
			fingerPrint = Compression.bytesToLong(temp);

//			bytesPerSampleMarker = (byte) (Compression.BYTES_PER_SAMPLE_MARKER - (nullStatus & 0x01) - (nullStatus >>1 & 0x01) - (nullStatus >>2 & 0x01) - (nullStatus >>3 & 0x01) - (nullStatus >>4 & 0x01) - (nullStatus >>5 & 0x01) - (nullStatus >>6 & 0x01));
			bytesPerSampleMarker = getNBytesPerSampleMarker(nullStatus);

//			numBytesOfOutOfRangeValues = Compression.bytesToInt(new byte[]{readBuffer[outlierSectionLocation], readBuffer[outlierSectionLocation+1], readBuffer[outlierSectionLocation+2], readBuffer[outlierSectionLocation+3]});
			if (numBytesOfOutOfRangeValues>0) {
				outlierSectionLocation = PARAMETER_SECTION_BYTES + numMarkers * bytesPerSampleMarker;
				outOfRangeValues = (Hashtable<String, Float>) Compression.bytesToObj(readBuffer, outlierSectionLocation, numBytesOfOutOfRangeValues);
			}

			indexStart = PARAMETER_SECTION_BYTES;
			index = indexStart;
			if (! isGcNull(nullStatus)) {
				if (loadGC) {
					gcs = new float[numMarkers];
					for (int j=0; j<numMarkers; j++) {
						gcs[j] = Compression.gcBafDecompress(new byte[] {readBuffer[index], readBuffer[index + 1]});
						index += bytesPerSampleMarker;
					}
				}
				indexStart += 2;
			}
			index = indexStart;
			if (! isXNull(nullStatus)) {
				if (loadXY) {
					xs = new float[numMarkers];
					for (int j=0; j<numMarkers; j++) {
						xs[j] = Compression.xyDecompressPositiveOnly(new byte[] {readBuffer[index], readBuffer[index + 1]});
						if (xs[j]==Compression.REDUCED_PRECISION_XY_OUT_OF_RANGE_FLAG_FLOAT) {
							xs[j] = outOfRangeValues.get(j + "\tx");
						}
						index += bytesPerSampleMarker;
					}
				}
				indexStart += 2;
			}
			index = indexStart;
			if (! isYNull(nullStatus)) {
				if (loadXY) {
					ys = new float[numMarkers];
					for (int j=0; j<numMarkers; j++) {
						ys[j] = Compression.xyDecompressPositiveOnly(new byte[] {readBuffer[index], readBuffer[index + 1]});
						if (ys[j]==Compression.REDUCED_PRECISION_XY_OUT_OF_RANGE_FLAG_FLOAT) {
							ys[j] = outOfRangeValues.get(j + "\ty");
						}
						index += bytesPerSampleMarker;
					}
					
				}
				indexStart += 2;
			}
			index = indexStart;
			if (! isBafNull(nullStatus)) {
				if (loadBAF) {
					bafs = new float[numMarkers];
					for (int j=0; j<numMarkers; j++) {
						bafs[j] = Compression.gcBafDecompress(new byte[] {readBuffer[index], readBuffer[index + 1]});
						index += bytesPerSampleMarker;
					}
				}
				indexStart += 2;
			}
			index = indexStart;
			if (! isLrrNull(nullStatus)) {
				if (loadLRR) {
					lrrs = new float[numMarkers];
					for (int j=0; j<numMarkers; j++) {
						lrrs[j] = Compression.lrrDecompress(new byte[] {readBuffer[index], readBuffer[index + 1], readBuffer[index + 2]});
						if (lrrs[j] == Compression.REDUCED_PRECISION_LRR_OUT_OF_RANGE_LRR_FLAG_FLOAT) {
							lrrs[j] = outOfRangeValues.get(j + "\tlrr");
						}
						index += bytesPerSampleMarker;
					}
				}
				indexStart += 3;
			}
			index = indexStart;
			if ((!isAbOrForwardGenotypeNull(nullStatus)) && loadAbOrForwardGenotypes) {
				abGenotypes = new byte[numMarkers];
				fwdGenotypes = new byte[numMarkers];
				for (int j=0; j<numMarkers; j++) {
					genoTypeTmp = Compression.genotypeDecompress(readBuffer[index]);
					abGenotypes[j] = genoTypeTmp[0];
					fwdGenotypes[j] = genoTypeTmp[1];
					index += bytesPerSampleMarker;
				}
			}
			result = new Sample(sampleName, fingerPrint, gcs, xs, ys, bafs, lrrs, fwdGenotypes, abGenotypes, isNegativeXOrYAllowed(nullStatus));
		} catch (IOException e) {
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
		}
		return result;
	}

//	@SuppressWarnings("unchecked")
//	public static void loadFromRandomAccessFileWithoutDecompress(RandomAccessFile sampleFile, byte[] readBuffer, int indexOfCurrentSample, int indexOfFirstMarkerToLoad, byte bytesPerSampleMarker, int numMarkersInProj, Hashtable<String, Float> allOutliers) {
//		byte[] outliersBuffer;
//		Hashtable<String, Float> sampleOutlierHash;
//		Enumeration<String> keys;
//		String currentKey;
//		int outlierSectionSize = 0;
//
//		try {
//	    	if (allOutliers != null) {
//	    			sampleFile.seek(Sample.PARAMETER_SECTION_OUTLIERSECTIONLENGTH_LOC);
//	    			outlierSectionSize = sampleFile.readInt();
//	    	}
//
//	    	sampleFile.seek(Sample.PARAMETER_SECTION_BYTES + indexOfFirstMarkerToLoad * bytesPerSampleMarker);
//	    	sampleFile.read(readBuffer);
//
//    		if (outlierSectionSize > 0) {
//		    	sampleFile.seek(Sample.PARAMETER_SECTION_BYTES + numMarkersInProj * bytesPerSampleMarker);
//	    		outliersBuffer = new byte[outlierSectionSize];
//	    		sampleFile.read(outliersBuffer);
//	    		sampleOutlierHash = (Hashtable<String, Float>) Compression.bytesToObj(outliersBuffer);
//	    		keys = sampleOutlierHash.keys();
//	    		while (keys.hasMoreElements()) {
//	    			currentKey = keys.nextElement();
//	    			allOutliers.put(indexOfCurrentSample + "\t" + currentKey, sampleOutlierHash.get(currentKey));
//	    		}
//    		}
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		} catch (ClassNotFoundException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
//	}

	@SuppressWarnings("unchecked")
	public static void loadFromRandomAccessFileWithoutDecompress(RandomAccessFile sampleFile, byte[] readBuffer, boolean seekOrLoadWholeFile, int indexOfCurrentSample, int indexOfFirstMarkerToLoad, byte bytesPerSampleMarker, int numMarkersInProj, Hashtable<String, Float> allOutliers, Logger log) {
		byte[] outliersBuffer;
		Hashtable<String, Float> sampleOutlierHash;
		Enumeration<String> keys;
		String currentKey;
		int outlierSectionSize = 0;
		byte[] readBufferLocal;
		int pointer;
		long seekPointer;

		try {
			if (seekOrLoadWholeFile) {
		    	if (allOutliers != null) {
		    		sampleFile.seek(Sample.PARAMETER_SECTION_OUTLIERSECTIONLENGTH_LOCATION);
		    		outlierSectionSize = sampleFile.readInt();
		    	}

		    	seekPointer = Sample.PARAMETER_SECTION_BYTES + indexOfFirstMarkerToLoad * bytesPerSampleMarker;
		    	if (seekPointer != sampleFile.getFilePointer()) {
		    		sampleFile.seek(seekPointer);
		    	}
				sampleFile.read(readBuffer);

				if (outlierSectionSize > 0) {
			    	sampleFile.seek(Sample.PARAMETER_SECTION_BYTES + numMarkersInProj * bytesPerSampleMarker);
		    		outliersBuffer = new byte[outlierSectionSize];
		    		sampleFile.read(outliersBuffer);
		    		sampleOutlierHash = (Hashtable<String, Float>) Compression.bytesToObj(outliersBuffer);
		    		keys = sampleOutlierHash.keys();
		    		while (keys.hasMoreElements()) {
		    			currentKey = keys.nextElement();
		    			allOutliers.put(indexOfCurrentSample + "\t" + currentKey, sampleOutlierHash.get(currentKey));
		    		}
	    		}

			} else {
				readBufferLocal = new byte[(int) sampleFile.length()];
				sampleFile.seek(0);
				sampleFile.read(readBufferLocal);
		    	pointer = Sample.PARAMETER_SECTION_BYTES + indexOfFirstMarkerToLoad * bytesPerSampleMarker;
		    	for (int i=0; i<readBuffer.length; i++) {
		    		readBuffer[i] = readBufferLocal[pointer];
		    		pointer ++;
		    	}

				if (allOutliers != null) {
		    		pointer = Sample.PARAMETER_SECTION_OUTLIERSECTIONLENGTH_LOCATION;
	    			outlierSectionSize = Compression.bytesToInt(readBufferLocal, pointer);
		    		if (outlierSectionSize > 0) {
		    			pointer = Sample.PARAMETER_SECTION_BYTES + numMarkersInProj * bytesPerSampleMarker;
			    		sampleOutlierHash = (Hashtable<String, Float>) Compression.bytesToObj(readBufferLocal, pointer, outlierSectionSize);
			    		keys = sampleOutlierHash.keys();
			    		while (keys.hasMoreElements()) {
			    			currentKey = keys.nextElement();
			    			allOutliers.put(indexOfCurrentSample + "\t" + currentKey, sampleOutlierHash.get(currentKey));
			    		}
		    		}
		    	}
			}
		} catch (IOException ioe) {
			log.reportError("Error reading from a "+Sample.SAMPLE_DATA_FILE_EXTENSION+" file");
			log.reportException(ioe);
		} catch (ClassNotFoundException cnfe) {
			log.reportError("Error reading from a "+Sample.SAMPLE_DATA_FILE_EXTENSION+" file");
			log.reportException(cnfe);
		}
	}


	/*
	 * 
	 */
	@SuppressWarnings("unchecked")
	public static void loadSampleFileWithoutDecompress(String sampleFileName, byte[] outputBuffer, Hashtable<String, Float> allOutliers, int indexOfCurrentSample, int indexOfFirstMarkerToLoad, byte bytesPerSampleMarker, int numMarkersInProj) {
		RandomAccessFile sampleFile;
		byte[] outliersBuffer;
		Hashtable<String, Float> sampleOutlierHash;
		Enumeration<String> keys;
		String currentKey;
		int outlierSectionSize = 0;
		long seekPointer;

		try {
			sampleFile = new RandomAccessFile(sampleFileName, "r");
	    	if (allOutliers != null) {
	    		sampleFile.seek(Sample.PARAMETER_SECTION_OUTLIERSECTIONLENGTH_LOCATION);
	    		outlierSectionSize = sampleFile.readInt();
	    	}

	    	seekPointer = Sample.PARAMETER_SECTION_BYTES + indexOfFirstMarkerToLoad * bytesPerSampleMarker;
	    	if (seekPointer != sampleFile.getFilePointer()) {
	    		sampleFile.seek(seekPointer);
	    	}
			sampleFile.read(outputBuffer);

			if (outlierSectionSize > 0) {
		    	sampleFile.seek(Sample.PARAMETER_SECTION_BYTES + numMarkersInProj * bytesPerSampleMarker);
	    		outliersBuffer = new byte[outlierSectionSize];
	    		sampleFile.read(outliersBuffer);
	    		sampleOutlierHash = (Hashtable<String, Float>) Compression.bytesToObj(outliersBuffer);
	    		keys = sampleOutlierHash.keys();
	    		while (keys.hasMoreElements()) {
	    			currentKey = keys.nextElement();
	    			allOutliers.put(indexOfCurrentSample + "\t" + currentKey, sampleOutlierHash.get(currentKey));
	    		}
    		}
		} catch (IOException e) {
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
		} catch (ArrayIndexOutOfBoundsException e) {
			e.printStackTrace();
		}
	}


	public static byte getNullstatusFromRandomAccessFile(String filename, boolean jar) {
		byte nullStatusOfTheFile = Byte.MIN_VALUE;
		RandomAccessFile sampleFile;

		try {
			sampleFile = new RandomAccessFile(filename, "r");
			sampleFile.readInt();
			nullStatusOfTheFile = sampleFile.readByte();
			sampleFile.close();
		} catch (FileNotFoundException e1) {
			e1.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return nullStatusOfTheFile;
	}

	@SuppressWarnings("unchecked")
	public static Hashtable<String, Float> loadOutOfRangeValues(RandomAccessFile file, long outlierSectionLocation) throws Exception {
        Hashtable<String, Float> outOfRangeValues = null;
        int numBytesOfOutOfRangeValues;
        byte[] readBuffer;
		
		file.seek(outlierSectionLocation);
		numBytesOfOutOfRangeValues = file.readInt();
		if (numBytesOfOutOfRangeValues>0) {
			readBuffer = new byte[numBytesOfOutOfRangeValues];
			file.read(readBuffer);
			outOfRangeValues = (Hashtable<String, Float>) Compression.bytesToObj(readBuffer);
		}
		
		return outOfRangeValues;
	}

	// TODO if Outlier class is created, then move these two methods there
	public static Hashtable<String, Float> loadOutOfRangeValuesFromSerializable(Project proj) throws Exception {
		return loadOutOfRangeValuesFromSerializable(proj.getDir(Project.SAMPLE_DIRECTORY, true) + "outliers.ser");
	}

	@SuppressWarnings("unchecked")
	public static Hashtable<String, Float> loadOutOfRangeValuesFromSerializable(String filenameOfSerializableOutOfRangeFile) throws Exception {
		return (Hashtable<String, Float>) Files.readSerial(filenameOfSerializableOutOfRangeFile);
	}

	public static void dumpOutOfRangeValues(String outOfRangeFileName, String outFileName, boolean includeKeysOnly) {
		Hashtable<String, Float> outOfRangeValues;

		if (outOfRangeFileName==null || outFileName==null) {
			System.err.println("Error - outOfRangeFileName and outFileName cannot be null.");
		}
		try {
			outOfRangeValues = loadOutOfRangeValuesFromSerializable(outOfRangeFileName);
			dumpOutOfRangeValues(outOfRangeValues, outFileName, includeKeysOnly);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void dumpOutOfRangeValues(Hashtable<String, Float> outOfRangeValues, String outFileName, boolean includeKeysOnly) {
		String[] keys;
		PrintWriter outfile;

		try {
			keys = HashVec.getKeys(outOfRangeValues);
			outfile = new PrintWriter(new FileOutputStream(outFileName));
			if (! includeKeysOnly) {
				for (int i = 0; i < keys.length; i ++) {
					outfile.println(keys[i] + "\t" + outOfRangeValues.get(keys[i]));
				}
			} else {
				for (int i = 0; i < keys.length; i ++) {
					outfile.println(keys[i]);
				}
			}
			outfile.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}

	/**
	 * Load data from Random Access Files organized by samples.
	 */
	public static void testLoadTime(Project proj) {
		String[] samplesProj;
		RandomAccessFile file;
		byte[] readBuffer;
		long time1;

		samplesProj = proj.getSamples();
		time1 = new Date().getTime();
		try {
//			for (int i=0; i<samplesProj.length; i++) {
			for (int i=0; i<100; i++) {
				file = new RandomAccessFile(proj.getDir(Project.SAMPLE_DIRECTORY, true) + samplesProj[i] + SAMPLE_DATA_FILE_EXTENSION, "r");
				readBuffer = new byte[(int) file.length()];	//numMarkers * BYTES_PER_SAMPLE_MARKER
				file.read(readBuffer);
				file.close();
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		System.out.println("Finished loadToMarkerData from RAF in: "+(new Date().getTime() - time1)+" ms");
	}

	public boolean hasBimodalBAF(byte chr, int startPosition, int endPosition) {
		DoubleVector dv = new DoubleVector();
		for (int i=0; i<bafs.length; i++) {
			if (bafs[i]>0.15 && bafs[i]<0.85) {
				dv.add(bafs[i]);
			}
		}
		return Array.isBimodal(dv.toArray());
	}

	public static void main(String[] args) {
		Project proj = new Project(cnv.Launch.getDefaultDebugProjectFile(true), false);
		String[] samples = proj.getSamples();
		Sample samp;
		
		for (int i = 0; i<samples.length; i++) {
			samp = proj.getFullSampleFromRandomAccessFile(samples[i]);
			samp.compareCalculationsFile(proj, proj.getMarkerNames(), proj.getProjectDir()+samples[i]+"_comp.xln");
        }

		//tests
//		testLoadTime(new Project("C:/workspace/Genvisis/projects/GEDI_exome.properties", false));
//		System.out.println(Compression.bytesToLong(Compression.longToBytes(15351532491l), 0));
    }

}
