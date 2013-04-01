package cnv.filesys;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Date;
import java.util.Hashtable;
import java.util.TreeMap;


import common.Array;
import common.DoubleVector;
import common.Files;
import common.ext;

public class Sample implements Serializable {
	public static final long serialVersionUID = 1L;
	public static final String[][] DATA_FIELDS = {{"GC Score", "GCscore"}, {"X Raw"}, {"Y Raw"}, {"X", "Xvalue"}, {"Y", "Yvalue"}, {"Theta"}, {"R"}, {"B Allele Freq"}, {"Log R Ratio"}};
	public static final String[][] GENOTYPE_FIELDS = {{"Allele1 - Forward", "Allele1"}, {"Allele2 - Forward", "Allele2"}, {"Allele1 - AB"}, {"Allele2 - AB"}};
	public static final String[] ALL_STANDARD_GENOTYPE_FIELDS = {"Allele1 - AB", "Allele2 - AB", "Allele1 - Forward", "Allele2 - Forward", "Allele1 - Top", "Allele2 - Top", "Allele1 - Design", "Allele2 - Design"};
	public static final String[] ALLELE_PAIRS = {"--", "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT", "DD", "DI", "II", "ID"};
	public static final String[] ALT_NULL = {"-", "0"};
	public static final String[] ALT_NULLS = {"--", "00"};
	public static final String[] AB_PAIRS = {"AA", "AB", "BB"};
	public static final String SAMPLE_DATA_FILE_EXTENSION = ".sampRAF";
	public static final byte PARAMETER_SECTION_BYTES = 13;
	public static final byte PARAMETER_SECTION_NUMMARK_LOC = 0;
	public static final byte PARAMETER_SECTION_NUMMARK_LEN = 4;
	public static final byte PARAMETER_SECTION_NULLSTAT_LOC = 4;
	public static final byte PARAMETER_SECTION_NULLSTAT_LEN = 1;
	public static final byte PARAMETER_SECTION_FINGPRNT_LOC = 5;
	public static final byte PARAMETER_SECTION_FINGPRNT_LEN = 8;
	public static final int MAX_ROWS_PER_WRITE_OPERATION = 500;
	public static final byte NULLSTATUS_GC_LOCATION = 0;
	public static final byte NULLSTATUS_X_LOCATION = 1;
	public static final byte NULLSTATUS_Y_LOCATION = 2;
	public static final byte NULLSTATUS_BAF_LOCATION = 3;
	public static final byte NULLSTATUS_LRR_LOCATION = 4;
	public static final byte NULLSTATUS_ABGENOTYPE_LOCATION = 5;
	public static final byte NULLSTATUS_FOWARDGENOTYPE_LOCATION = 6;

	private byte[] abGenotypes;
	private byte[] forwardGenotypes;
	private float[] xs;
	private float[] ys;
	private float[] thetas;
	private float[] rs;
	private float[] lrrs;
	private float[] bafs;
	private float[] gcs;
//	private float[] xRaws;
//	private float[] yRaws;
	private long fingerprint;
	private String sampleName;
	private byte nullStatus;

//	public Sample(String sampleName, long fingerprint, float[] gcs, float[] xs, float[] ys, float[] bafs, float[] lrrs, byte[] forwardGenotypes, byte[] abGenotypes) {
//		this.sampleName = sampleName;
//		this.fingerprint = fingerprint;
//		this.gcs = gcs;
////		this.xRaws = xRaws;
////		this.yRaws = yRaws;
//		this.xs = xs;
//		this.ys = ys;
////		this.rs = rs;
////		this.thetas = thetas;
//		this.bafs = bafs;
//		this.lrrs = lrrs;
//		this.forwardGenotypes = forwardGenotypes;
//		this.abGenotypes = abGenotypes;
//	}

	public Sample(String sampleName, long fingerprint, float[] gcs, float[] xs, float[] ys, float[] bafs, float[] lrrs, byte[] forwardGenotypes, byte[] abGenotypes) {
		this.sampleName = sampleName;
		this.fingerprint = fingerprint;
		this.gcs = gcs;
		this.xs = xs;
		this.ys = ys;
		this.bafs = bafs;
		this.lrrs = lrrs;
		this.forwardGenotypes = forwardGenotypes;
		this.abGenotypes = abGenotypes;
		updateNullStatus();
		// note for future expansion: 8th bit should indicate a second nullStatus byte will be used, to ensure backward compatibility
	}

//	public Sample(String sampleName, long fingerprint, float[] gcs, float[] xRaws, float[] yRaws, float[] xs, float[] ys, float[] rs, float[] thetas, float[] bafs, float[] lrrs, byte[] forwardGenotypes, byte[] abGenotypes) {
//		this.sampleName = sampleName;
//		this.fingerprint = fingerprint;
//		this.gcs = gcs;
////		this.xRaws = xRaws;
////		this.yRaws = yRaws;
//		this.xs = xs;
//		this.ys = ys;
////		this.rs = rs;
////		this.thetas = thetas;
////		this.rs = rs;
//		this.bafs = bafs;
//		this.lrrs = lrrs;
//		this.forwardGenotypes = forwardGenotypes;
//		this.abGenotypes = abGenotypes;
//	}

//	public FullSample(String sampleName, long fingerprint, float[][] data, byte[][] genotypes) {
//		this.sampleName = sampleName;
//		this.fingerprint = fingerprint;
//		this.gcs = data[0];
//		this.xRaws = data[1];
//		this.yRaws = data[2];
//		this.xs = data[3];
//		this.ys = data[4];
//		this.thetas = data[5];
//		this.rs = data[6];
//		this.bafs = data[7];
//		this.lrrs = data[8];
//		this.forwardGenotypes = genotypes[0];
//		this.abGenotypes = genotypes[1];
//	}
//

	public Sample(String sampleName, long fingerprint, float[][] data, byte[][] genotypes) {
		this.sampleName = sampleName;
		this.fingerprint = fingerprint;
		this.gcs = data[0];
//		this.xRaws = null;
//		this.yRaws = null;
		this.xs = data[3];
		this.ys = data[4];
//		this.thetas = null;
//		this.rs = null;
		this.bafs = data[7];
		this.lrrs = data[8];
		this.forwardGenotypes = genotypes[0];
		this.abGenotypes = genotypes[1];
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
		if (lrrs==null) {
			nullStatus = (byte) (nullStatus | (1 << NULLSTATUS_LRR_LOCATION));
		}
		if (bafs==null) {
			nullStatus = (byte) (nullStatus | (1 << NULLSTATUS_BAF_LOCATION));
		}
		if (abGenotypes==null) {
			nullStatus = (byte) (nullStatus | (1 << NULLSTATUS_ABGENOTYPE_LOCATION));
		}
		if (forwardGenotypes==null) {
			nullStatus = (byte) (nullStatus | (1 << NULLSTATUS_FOWARDGENOTYPE_LOCATION));
		}
	}

	public float[][] getAllData() {
//		return new float[][] {gcs, xRaws, yRaws, xs, ys, thetas, rs, bafs, lrrs};
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

	/*
	public float[] getX_Raws() {
		return xRaws;
	}

	public float[] getY_Raws() {
		return yRaws;
	}
	*/

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
	

	public byte[] getAB_GenotypesAfterFilters(String[] markerNames, ClusterFilterCollection clusterFilterCollection, float gcThreshold) {
		byte[] result;
		float realX;
		float realY;
		ClusterFilter clusterFilter;
		ArrayList<ClusterFilter> clusterFilterArray;

		result = new byte[abGenotypes.length];
		for (int i=0; i<markerNames.length; i++) {
			if (getGCs()[i]<gcThreshold) {
				result[i]=(byte)-1;
			} else {
				result[i]=abGenotypes[i];
			}
			
			clusterFilterArray = clusterFilterCollection.getClusterFilters(markerNames[i]);
			if (clusterFilterArray != null) {
				for (int j=0; j<clusterFilterArray.size(); j++) {
					clusterFilter = clusterFilterArray.get(j);
					switch(clusterFilter.getPlotType()) {
//					case 0:
//						realX = xRaws[i];
//						realY = yRaws[i];
//						break;
//					case 1:
//						realX = xs[i];
//						realY = ys[i];
//						break;
//					case 2:
//						realX = thetas[i];
//						realY = rs[i];
//						break;
//					case 3:
//						realX = bafs[i];
//						realY = lrrs[i];
//						break;
					case 0:
						realX = xs[i];
						realY = ys[i];
						break;
					case 1:
						realX = thetas[i];
						realY = rs[i];
						break;
					case 2:
						realX = bafs[i];
						realY = lrrs[i];
						break;
					default:
						System.err.println("Error - invalid PlotType in ClusterFilter #"+(j+1)+" for marker '"+markerNames[i]+"'");
						realX = -9;
						realY = -9;
					}
					if (realX>=clusterFilter.getXMin() && realY>=clusterFilter.getYMin() && realX<=clusterFilter.getXMax() && realY<=clusterFilter.getYMax()) {
						result[i]=clusterFilter.getNewGenotype();
					}
				}
			}
		}
		return result;
	}
	
	public byte[] getAB_Genotypes() {
		return abGenotypes;
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
//			writer.println("SNP\tGC Score\tRaw X\tRaw Y\tX\tY\tTheta\tR\tLRR\tBAF\tGenotypes");
			writer.println("SNP\tGC Score\tX\tY\tLRR\tBAF\tGenotypes");
			for (int i = 0; i<markerNames.length; i++) {
//				writer.println(markerNames[i]+"\t"+gcs[i]+"\t"+xRaws[i]+"\t"+yRaws[i]+"\t"+xs[i]+"\t"+ys[i]+"\t"+thetas[i]+"\t"+rs[i]+"\t"+lrrs[i]+"\t"+bafs[i]+"\t"+ALLELE_PAIRS[forwardGenotypes[i]]);
				writer.println(markerNames[i]+"\t"+gcs[i]+"\t"+"\t"+xs[i]+"\t"+ys[i]+"\t"+lrrs[i]+"\t"+bafs[i]+"\t"+ALLELE_PAIRS[forwardGenotypes[i]]);
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
//			for (int i = 0; i<1000; i++) {
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

//	/**
//	 * Save the data to Random Access File.
//	 * @param filename
//	 * @param allMarkersInProj
//	 */
//	public void saveToRandomAccessFile3(String filename) {
//		File fileTmp;
//		RandomAccessFile rafFile;
//		Hashtable<String, Float> outOfRangeValuesEachSample;
//		byte[] outOfRangeValuesWriteBuffer;
//		int bytesRemained;
//		byte[] writeBuffer = null;
//		int writeBufferIndex;
//		byte bytesPerSampMark;
//		
//		fileTmp = new File(filename);
////		if (fileTmp.getFreeSpace() <= ((long)allMarkersInProj.length * Compression.BYTES_PER_SAMPLE_MARKER_2)) {
//		if (new File(ext.parseDirectoryOfFile(filename)).getFreeSpace() <= ((long)xs.length * Compression.BYTES_PER_SAMPLE_MARKER_2)) {
//			System.err.println("Not enough space (available: "+ext.prettyUpSize(new File(ext.parseDirectoryOfFile(filename)).getFreeSpace(), 1)+") for all the new data to be created (required: "+ext.prettyUpSize(((long)xs.length * Compression.BYTES_PER_SAMPLE_MARKER_2), 1)+").");
//			return;
//		}
//		if (fileTmp.exists()) {
//			fileTmp.delete();
////			DataInRafByMarker.backupOlderRafs(filename, ".fsampRaf");
//		}
//
//		bytesPerSampMark = (byte) (Compression.BYTES_PER_SAMPLE_MARKER_2 - (nullStatus & 0x01) - (nullStatus >>1 & 0x01) - (nullStatus >>2 & 0x01) - (nullStatus >>3 & 0x01) - (nullStatus >>4 & 0x01) - (nullStatus >>5 & 0x01) - (nullStatus >>6 & 0x01));
//		bytesRemained = xs.length *  bytesPerSampMark;
//		long time = new Date().getTime();
//		outOfRangeValuesEachSample = new Hashtable<String, Float>();
//		try {
//			rafFile = new RandomAccessFile(filename, "rw");
////			rafFile.writeInt(allMarkersInProj.length);
//			rafFile.writeInt(xs.length);
//			rafFile.writeByte(nullStatus);
//			rafFile.writeLong(fingerprint);
//			writeBufferIndex = 0;
//			writeBuffer = new byte[ Math.min(Integer.MAX_VALUE, bytesRemained) ];
//			if (writeBufferIndex == 0) {
//			}
//			if (gcs != null) {
//				for (int j = 0; j<xs.length; j++) {
//					Compression.reducedPrecisionGcBafGetBytes2(gcs[j], writeBuffer, writeBufferIndex);
//					writeBufferIndex += 2;
//				}
//			}
//			if (xs != null) {
//				if (Compression.reducedPrecisionXYGetBytes2(xs[j], writeBuffer, writeBufferIndex) == -1) {
////						outOfRangeValuesEachSample.put(sampleName+"\t"+allMarkersInProj[j]+"\tx", xs[j]);
//					outOfRangeValuesEachSample.put(j + "\tx", xs[j]);
//				}
//				writeBufferIndex += 2;
//			}
//			if (ys != null) {
//				if (Compression.reducedPrecisionXYGetBytes2(ys[j], writeBuffer, writeBufferIndex) == -1) {
////						outOfRangeValuesEachSample.put(sampleName+"\t"+allMarkersInProj[j]+"\ty", ys[j]);
//					outOfRangeValuesEachSample.put(j + "\ty", ys[j]);
//				}
//				writeBufferIndex += 2;
//			}
//			if (bafs != null) {
//				Compression.reducedPrecisionGcBafGetBytes2(bafs[j], writeBuffer, writeBufferIndex);
//				writeBufferIndex += 2;
//			}
//			if (lrrs != null) {
//				if (j==3942) {
//					System.out.println();
//				}
//				if (Compression.reducedPrecisionLrrGetBytes2(lrrs[j], writeBuffer, writeBufferIndex) == -1) {
////						outOfRangeValuesEachSample.put(sampleName+"\t"+allMarkersInProj[j]+"\tlrr", lrrs[j]);
//					outOfRangeValuesEachSample.put(j + "\tlrr", lrrs[j]);
//				}
//				writeBufferIndex += 3;
//			}
//			if (abGenotypes != null || forwardGenotypes != null) {
//				writeBuffer[writeBufferIndex] = Compression.reducedPrecisionGenotypeGetBytes(abGenotypes == null?-1:abGenotypes[j], forwardGenotypes == null?0:forwardGenotypes[j]);
//				writeBufferIndex += 1;
//			}
////				bytesRemained -= Compression.BYTES_PER_SAMPLE_MARKER_2;
//			bytesRemained -= bytesPerSampMark;
//			if ( writeBufferIndex>=writeBuffer.length ) {
//				rafFile.write(writeBuffer);
//				writeBufferIndex = 0;
//			}
//			if (outOfRangeValuesEachSample==null || outOfRangeValuesEachSample.size()==0) {
//				rafFile.writeInt(0);
//			} else {
//				outOfRangeValuesWriteBuffer = Compression.objToBytes(outOfRangeValuesEachSample);
//				rafFile.writeInt(outOfRangeValuesWriteBuffer.length);
//				rafFile.write(outOfRangeValuesWriteBuffer);
//			}
//			rafFile.close();
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//		time = new Date().getTime() - time;
////		System.out.println("Random Access File. ---- Finished writing to all the samples in " + (time/60000) + " min " + ((time%60000)/1000) + " sec");
//	}

	/**
	 * Save the data to Random Access File.
	 * @param filename
	 * @param allMarkersInProj
	 */
	public void saveToRandomAccessFile2(String filename) {
		File fileTmp;
		RandomAccessFile rafFile;
		Hashtable<String, Float> outOfRangeValuesEachSample;
		byte[] outOfRangeValuesWriteBuffer;
		int bytesRemained;
		byte[] writeBuffer = null;
		int writeBufferIndex;
		byte bytesPerSampMark;
		
		fileTmp = new File(filename);
//		if (fileTmp.getFreeSpace() <= ((long)allMarkersInProj.length * Compression.BYTES_PER_SAMPLE_MARKER_2)) {
		if (new File(ext.parseDirectoryOfFile(filename)).getFreeSpace() <= ((long)xs.length * Compression.BYTES_PER_SAMPLE_MARKER_2)) {
			System.err.println("Not enough space (available: "+ext.prettyUpSize(new File(ext.parseDirectoryOfFile(filename)).getFreeSpace(), 1)+") for all the new data to be created (required: "+ext.prettyUpSize(((long)xs.length * Compression.BYTES_PER_SAMPLE_MARKER_2), 1)+").");
			return;
		}
		if (fileTmp.exists()) {
			fileTmp.delete();
//			DataInRafByMarker.backupOlderRafs(filename, ".fsampRaf");
		}

		bytesPerSampMark = (byte) (Compression.BYTES_PER_SAMPLE_MARKER_2 - (nullStatus & 0x01) - (nullStatus >>1 & 0x01) - (nullStatus >>2 & 0x01) - (nullStatus >>3 & 0x01) - (nullStatus >>4 & 0x01) - (nullStatus >>5 & 0x01) - (nullStatus >>6 & 0x01));
		bytesRemained = xs.length *  bytesPerSampMark;
		long time = new Date().getTime();
		outOfRangeValuesEachSample = new Hashtable<String, Float>();
		try {
			rafFile = new RandomAccessFile(filename, "rw");
//			rafFile.writeInt(allMarkersInProj.length);
			rafFile.writeInt(xs.length);
			rafFile.writeByte(nullStatus);
			rafFile.writeLong(fingerprint);
			writeBufferIndex = 0;
			for (int j = 0; j<xs.length; j++) {
				if (writeBufferIndex == 0) {
					writeBuffer = new byte[ Math.min(Integer.MAX_VALUE, bytesRemained) ];
				}
				if (gcs != null) {
					Compression.reducedPrecisionGcBafGetBytes2(gcs[j], writeBuffer, writeBufferIndex);
					writeBufferIndex += 2;
				}
				if (xs != null) {
					if (Compression.reducedPrecisionXYGetBytes2(xs[j], writeBuffer, writeBufferIndex) == -1) {
//						outOfRangeValuesEachSample.put(sampleName+"\t"+allMarkersInProj[j]+"\tx", xs[j]);
						outOfRangeValuesEachSample.put(j + "\tx", xs[j]);
					}
					writeBufferIndex += 2;
				}
				if (ys != null) {
					if (Compression.reducedPrecisionXYGetBytes2(ys[j], writeBuffer, writeBufferIndex) == -1) {
//						outOfRangeValuesEachSample.put(sampleName+"\t"+allMarkersInProj[j]+"\ty", ys[j]);
						outOfRangeValuesEachSample.put(j + "\ty", ys[j]);
					}
					writeBufferIndex += 2;
				}
				if (bafs != null) {
					Compression.reducedPrecisionGcBafGetBytes2(bafs[j], writeBuffer, writeBufferIndex);
					writeBufferIndex += 2;
				}
				if (lrrs != null) {
					if (Compression.reducedPrecisionLrrGetBytes2(lrrs[j], writeBuffer, writeBufferIndex) == -1) {
//						outOfRangeValuesEachSample.put(sampleName+"\t"+allMarkersInProj[j]+"\tlrr", lrrs[j]);
						outOfRangeValuesEachSample.put(j + "\tlrr", lrrs[j]);
					}
					writeBufferIndex += 3;
				}
				if (abGenotypes != null || forwardGenotypes != null) {
					writeBuffer[writeBufferIndex] = Compression.reducedPrecisionGenotypeGetBytes(abGenotypes == null?-1:abGenotypes[j], forwardGenotypes == null?0:forwardGenotypes[j]);
					writeBufferIndex += 1;
				}
//				bytesRemained -= Compression.BYTES_PER_SAMPLE_MARKER_2;
				bytesRemained -= bytesPerSampMark;
				if ( writeBufferIndex>=writeBuffer.length ) {
					rafFile.write(writeBuffer);
					writeBufferIndex = 0;
				}
			}
			if (outOfRangeValuesEachSample==null || outOfRangeValuesEachSample.size()==0) {
				rafFile.writeInt(0);
			} else {
				outOfRangeValuesWriteBuffer = Compression.objToBytes(outOfRangeValuesEachSample);
				rafFile.writeInt(outOfRangeValuesWriteBuffer.length);
				rafFile.write(outOfRangeValuesWriteBuffer);
			}
			rafFile.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		time = new Date().getTime() - time;
//		System.out.println("Random Access File. ---- Finished writing to all the samples in " + (time/60000) + " min " + ((time%60000)/1000) + " sec");
	}

	/**
	 * Save the instance of Sample to hard drive in Random Access File format.
	 * @param filename
	 * @param allMarkersInProj
	 */
//	public void saveToRandomAccessFile(String filename, String[] allMarkersInProj) {
	public void saveToRandomAccessFile(String filename) {
		File fileTmp;
		RandomAccessFile rafFile;
		Hashtable<String, Float> outOfRangeValuesEachSample;
		byte[] outOfRangeValuesWriteBuffer;
		int bytesRemained;
		byte[] writeBuffer = null;
		int writeBufferIndex;
		
		fileTmp = new File(filename);
//		if (fileTmp.getFreeSpace() <= ((long)allMarkersInProj.length * Compression.BYTES_PER_SAMPLE_MARKER_2)) {
		if (new File(ext.parseDirectoryOfFile(filename)).getFreeSpace() <= ((long)xs.length * Compression.BYTES_PER_SAMPLE_MARKER_2)) {
			System.err.println("Not enough space (available: "+ext.prettyUpSize(new File(ext.parseDirectoryOfFile(filename)).getFreeSpace(), 1)+") for all the new data to be created (required: "+ext.prettyUpSize(((long)xs.length * Compression.BYTES_PER_SAMPLE_MARKER_2), 1)+").");
			return;
		}
		if (fileTmp.exists()) {
			fileTmp.delete();
//			DataInRafByMarker.backupOlderRafs(filename, ".fsampRaf");
		}

//		bytesRemained = allMarkersInProj.length * Compression.BYTES_PER_SAMPLE_MARKER_2;
		bytesRemained = xs.length * Compression.BYTES_PER_SAMPLE_MARKER_2;
		long time = new Date().getTime();
		outOfRangeValuesEachSample = new Hashtable<String, Float>();
		try {
			rafFile = new RandomAccessFile(filename, "rw");
//			rafFile.writeInt(allMarkersInProj.length);
			rafFile.writeInt(xs.length);
			rafFile.writeLong(fingerprint);
			writeBufferIndex = 0;
			for (int j = 0; j<xs.length; j++) {
				if (writeBufferIndex == 0) {
					writeBuffer = new byte[ Math.min(Integer.MAX_VALUE, bytesRemained) ];
				}
				Compression.reducedPrecisionGcBafGetBytes2(gcs[j], writeBuffer, writeBufferIndex);
				writeBufferIndex += 2;
				if (Compression.reducedPrecisionXYGetBytes2(xs[j], writeBuffer, writeBufferIndex) == -1) {
//					outOfRangeValuesEachSample.put(sampleName+"\t"+allMarkersInProj[j]+"\tx", xs[j]);
					outOfRangeValuesEachSample.put(sampleName+"\t" + j + "\tx", xs[j]);
				}
				writeBufferIndex += 2;
				if (Compression.reducedPrecisionXYGetBytes2(ys[j], writeBuffer, writeBufferIndex) == -1) {
//					outOfRangeValuesEachSample.put(sampleName+"\t"+allMarkersInProj[j]+"\ty", ys[j]);
					outOfRangeValuesEachSample.put(sampleName + "\t" + j + "\ty", ys[j]);
				}
				writeBufferIndex += 2;
				Compression.reducedPrecisionGcBafGetBytes2(bafs[j], writeBuffer, writeBufferIndex);
				writeBufferIndex += 2;
				if (Compression.reducedPrecisionLrrGetBytes2(lrrs[j], writeBuffer, writeBufferIndex) == -1) {
//					outOfRangeValuesEachSample.put(sampleName+"\t"+allMarkersInProj[j]+"\tlrr", lrrs[j]);
					outOfRangeValuesEachSample.put(sampleName + "\t" + j + "\tlrr", lrrs[j]);
				}
				writeBufferIndex += 3;
				writeBuffer[writeBufferIndex] = Compression.reducedPrecisionGenotypeGetBytes(abGenotypes == null?-1:abGenotypes[j], forwardGenotypes == null?0:forwardGenotypes[j]);
				writeBufferIndex += 1;
				bytesRemained -= Compression.BYTES_PER_SAMPLE_MARKER_2;
				if ( writeBufferIndex>=writeBuffer.length ) {
					rafFile.write(writeBuffer);
					writeBufferIndex = 0;
				}
			}
			outOfRangeValuesWriteBuffer = Compression.objToBytes(outOfRangeValuesEachSample);
			if (outOfRangeValuesWriteBuffer==null || outOfRangeValuesWriteBuffer.length==0) {
				rafFile.writeInt(0);
			} else {
				rafFile.writeInt(outOfRangeValuesWriteBuffer.length);
				rafFile.write(outOfRangeValuesWriteBuffer);
			}
			rafFile.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		time = new Date().getTime() - time;
		System.out.println("Random Access File. ---- Finished writing to all the samples in " + (time/60000) + " min " + ((time%60000)/1000) + " sec");
	}

//	public Sample_old convertToSample() {
//		return new Sample_old(fingerprint, lrrs, bafs, abGenotypes);
//	}

	public static Sample loadFromSerialized(String filename, boolean jar) {
		return (Sample)Files.readSerial(filename, jar, true);
	}
	
//	public static FullSample loadFromSerialized(String filename, boolean jar) {
//		return (FullSample)Files.readSerial(filename, jar, true);
//	}
	
//	public static Sample loadFromRandomAccessFile3(String filename, boolean jar) {
//		return loadFromRandomAccessFile2(filename, true, true, true, true, true, jar);
//	}
	
//	/**
//	 * Load data from Random Access Files organized by samples.
//	 */
//	public static Sample loadFromRandomAccessFile3(String filename, boolean loadGC, boolean loadXY, boolean loadBAF, boolean loadLRR, boolean loadAbOrForwardGenotypes, boolean jar) {
//		Sample result = null;
//		int numMarkers;
//		RandomAccessFile file;
//		String sampleName;
//		long fingerPrint;
//        float[] gcs = null, xs = null, ys = null, lrrs = null, bafs = null;
//        byte[] abGenotypes = null, fwdGenotypes = null, readBuffer;
//        byte[] genoTypeTmp;
//        int index, indexStart;
//        int numBytesOfOutOfRangeValues;
//        Hashtable<String, Float> outOfRangeValues = null;
//        int outlierSectionLocation;
//        byte nullStatus;
//        byte bytesPerSampMark;
////      long time;
//
////		What would happen if not all the markers have exactly the same samples? In SQL this will result in shift the rows and will affect the matching of the data.
//        sampleName = ext.rootOf(filename);
////		long time1 = new Date().getTime();
//		try {
////			time = new Date().getTime();
//			file = new RandomAccessFile(filename, "r");
//			readBuffer = new byte[(int) file.length()];	//numMarkers * BYTES_PER_SAMPLE_MARKER_2
//			file.read(readBuffer);
//			file.close();
//			numMarkers = Compression.bytesToInt(new byte[]{readBuffer[0], readBuffer[1], readBuffer[2], readBuffer[3]});
//			nullStatus = readBuffer[4];
//			fingerPrint = Compression.bytesToLong(new byte[]{readBuffer[5], readBuffer[6], readBuffer[7], readBuffer[8], readBuffer[9], readBuffer[10], readBuffer[11], readBuffer[12]});
//			bytesPerSampMark = (byte) (Compression.BYTES_PER_SAMPLE_MARKER_2 - (nullStatus & 0x01) - (nullStatus >>1 & 0x01) - (nullStatus >>2 & 0x01) - (nullStatus >>3 & 0x01) - (nullStatus >>4 & 0x01) - (nullStatus >>5 & 0x01) - (nullStatus >>6 & 0x01));
//			outlierSectionLocation = PARAMETER_SECTION_BYTES + numMarkers * bytesPerSampMark;
//			numBytesOfOutOfRangeValues = Compression.bytesToInt(new byte[]{readBuffer[outlierSectionLocation], readBuffer[outlierSectionLocation+1], readBuffer[outlierSectionLocation+2], readBuffer[outlierSectionLocation+3]});
//			if (numBytesOfOutOfRangeValues>0) {
//				outOfRangeValues = (Hashtable<String, Float>) Compression.bytesToObj(readBuffer, outlierSectionLocation+4, outlierSectionLocation + numBytesOfOutOfRangeValues);
//			}
////			System.out.print("\t"+(new Date().getTime() - time));
//
//			indexStart = PARAMETER_SECTION_BYTES;
//			index = indexStart;
////			time = new Date().getTime();
//			if (((nullStatus>>NULLSTATUS_GC_LOCATION) & 0x01) != 1) {
//				if (loadGC) {
//					gcs = new float[numMarkers];
//					for (int j=0; j<numMarkers; j++) {
//						gcs[j] = Compression.reducedPrecisionGcBafGetFloat2(new byte[] {readBuffer[index], readBuffer[index + 1]});
//						index += bytesPerSampMark;
//					}
//				}
//				indexStart += 2;
//			}
//			index = indexStart;
//			if (((nullStatus>>NULLSTATUS_X_LOCATION) & 0x01) != 1) {
//				if (loadXY) {
//					xs = new float[numMarkers];
//					for (int j=0; j<numMarkers; j++) {
//						xs[j] = Compression.reducedPrecisionXYGetFloat2(new byte[] {readBuffer[index], readBuffer[index + 1]});
//						if (xs[j]==-1) {
////							xs[j] = outOfRangeValues.get(sampleName+"\t"+allMarkersProj[j]+"\tx");
//							xs[j] = outOfRangeValues.get(j + "\tx");
//						}
//						index += bytesPerSampMark;
//					}
//				}
//				indexStart += 2;
//			}
//			index = indexStart;
//			if (((nullStatus>>NULLSTATUS_Y_LOCATION) & 0x01) != 1) {
//				if (loadXY) {
//					ys = new float[numMarkers];
//					for (int j=0; j<numMarkers; j++) {
//						ys[j] = Compression.reducedPrecisionXYGetFloat2(new byte[] {readBuffer[index], readBuffer[index + 1]});
//						if (ys[j]==-1) {
////							ys[j] = outOfRangeValues.get(sampleName+"\t"+allMarkersProj[j]+"\ty");
//							ys[j] = outOfRangeValues.get(j + "\ty");
//						}
//						index += bytesPerSampMark;
//					}
//					
//				}
//				indexStart += 2;
//			}
//			index = indexStart;
//			if (((nullStatus>>NULLSTATUS_BAF_LOCATION) & 0x01) != 1) {
//				if (loadBAF) {
//					bafs = new float[numMarkers];
//					for (int j=0; j<numMarkers; j++) {
//						bafs[j] = Compression.reducedPrecisionGcBafGetFloat2(new byte[] {readBuffer[index], readBuffer[index + 1]});
//						index += bytesPerSampMark;
//					}
//				}
//				indexStart += 2;
//			}
//			index = indexStart;
//			if (((nullStatus>>NULLSTATUS_LRR_LOCATION) & 0x01) != 1) {
//				if (loadLRR) {
//					lrrs = new float[numMarkers];
//					for (int j=0; j<numMarkers; j++) {
//						lrrs[j] = Compression.reducedPrecisionLrrGetFloat2(new byte[] {readBuffer[index], readBuffer[index + 1], readBuffer[index + 2]});
//						if (lrrs[j] == Compression.REDUCED_PRECISION_LRR_OUT_OF_RANGE_LRR_FLOAT) {
////							lrrs[j] = outOfRangeValues.get(sampleName+"\t"+allMarkersProj[j]+"\tlrr");
//							lrrs[j] = outOfRangeValues.get(j + "\tlrr");
//						}
//						index += bytesPerSampMark;
//					}
//				}
//				indexStart += 3;
//			}
//			index = indexStart;
//			if ((((nullStatus>>NULLSTATUS_ABGENOTYPE_LOCATION) & 0x01) != 1 || ((nullStatus>>NULLSTATUS_FOWARDGENOTYPE_LOCATION) & 0x01) != 1) && loadAbOrForwardGenotypes) {
//				abGenotypes = new byte[numMarkers];
//				fwdGenotypes = new byte[numMarkers];
//				for (int j=0; j<numMarkers; j++) {
//					genoTypeTmp = Compression.reducedPrecisionGenotypeGetTypes(readBuffer[index]);
//					abGenotypes[j] = genoTypeTmp[0];
//					fwdGenotypes[j] = genoTypeTmp[1];
//					index += bytesPerSampMark;
//				}
//			}
////			System.out.print("\t"+(new Date().getTime() - time));
////			time = new Date().getTime();
//			result = new Sample(sampleName, fingerPrint, gcs, xs, ys, bafs, lrrs, fwdGenotypes, abGenotypes);
////			System.out.println("\t"+(new Date().getTime() - time));
//		} catch (IOException e) {
//			e.printStackTrace();
//		} catch (ClassNotFoundException e) {
//			e.printStackTrace();
//		}
////		System.out.println("Finished loadToMarkerData from RAF in: "+(new Date().getTime() - time1)+" ms");
//		return result;
//	}

//	public static byte[] loadByteStreamFromRandomAccessFile3(String filename, boolean jar, int beginMarkerIndex, int endMarkerIndex) {
//		int numMarkers;
//		RandomAccessFile file;
//		byte[] parameterSection = null, readBuffer = null;
//        int numBytesOfOutOfRangeValues;
//        int outlierSectionLocation;
//        byte nullStatus;
//        byte bytesPerSampMark;
//
//        ext.rootOf(filename);
//		try {
//			parameterSection = new byte[PARAMETER_SECTION_BYTES];	//numMarkers * BYTES_PER_SAMPLE_MARKER_2
//			file = new RandomAccessFile(filename, "r");
//			file.read(parameterSection);
//			numMarkers = Compression.bytesToInt(new byte[]{parameterSection[0], parameterSection[1], parameterSection[2], parameterSection[3]});
////			numMarkers = endMarkerIndex - beginMarkerIndex + 1;
//			nullStatus = parameterSection[4];
//			Compression.bytesToLong(new byte[]{parameterSection[5], parameterSection[6], parameterSection[7], parameterSection[8], parameterSection[9], parameterSection[10], parameterSection[11], parameterSection[12]});
//			bytesPerSampMark = (byte) (Compression.BYTES_PER_SAMPLE_MARKER_2 - (nullStatus & 0x01) - (nullStatus >>1 & 0x01) - (nullStatus >>2 & 0x01) - (nullStatus >>3 & 0x01) - (nullStatus >>4 & 0x01) - (nullStatus >>5 & 0x01) - (nullStatus >>6 & 0x01));
//			outlierSectionLocation = PARAMETER_SECTION_BYTES + numMarkers * bytesPerSampMark;
//			numBytesOfOutOfRangeValues = Compression.bytesToInt(new byte[]{parameterSection[outlierSectionLocation], parameterSection[outlierSectionLocation+1], parameterSection[outlierSectionLocation+2], parameterSection[outlierSectionLocation+3]});
//			readBuffer = new byte[PARAMETER_SECTION_BYTES + (endMarkerIndex - beginMarkerIndex + 1) * bytesPerSampMark + numBytesOfOutOfRangeValues];
//			file.seek(pos);
//			file.read(readBuffer, PARAMETER_SECTION_BYTES+beginMarkerIndex*bytesPerSampMark, PARAMETER_SECTION_BYTES+endMarkerIndex*bytesPerSampMark);
//			file.close();
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
////		System.out.println("Finished loadToMarkerData from RAF in: "+(new Date().getTime() - time1)+" ms");
//		return parameterSection;
//	}

	public static Sample loadFromRandomAccessFile2(String filename, boolean jar) {
		return loadFromRandomAccessFile2(filename, true, true, true, true, true, jar);
	}

	/**
	 * Load data from Random Access Files organized by samples.
	 */
	public static Sample loadFromRandomAccessFile2(String filename, boolean loadGC, boolean loadXY, boolean loadBAF, boolean loadLRR, boolean loadAbOrForwardGenotypes, boolean jar) {
		Sample result = null;
		int numMarkers;
		RandomAccessFile file;
		String sampleName;
		long fingerPrint;
        float[] gcs = null, xs = null, ys = null, lrrs = null, bafs = null;
        byte[] abGenotypes = null, fwdGenotypes = null, readBuffer;
        byte[] genoTypeTmp;
        int index, indexStart;
        int numBytesOfOutOfRangeValues;
        Hashtable<String, Float> outOfRangeValues = null;
        int outlierSectionLocation;
        byte nullStatus;
        byte bytesPerSampMark;
//      long time;

//		What would happen if not all the markers have exactly the same samples? In SQL this will result in shift the rows and will affect the matching of the data.
        sampleName = ext.rootOf(filename);
//		long time1 = new Date().getTime();
		try {
//			time = new Date().getTime();
			file = new RandomAccessFile(filename, "r");
			readBuffer = new byte[(int) file.length()];	//numMarkers * BYTES_PER_SAMPLE_MARKER_2
			file.read(readBuffer);
			file.close();
			numMarkers = Compression.bytesToInt(new byte[]{readBuffer[0], readBuffer[1], readBuffer[2], readBuffer[3]});
			nullStatus = readBuffer[4];
			fingerPrint = Compression.bytesToLong(new byte[]{readBuffer[5], readBuffer[6], readBuffer[7], readBuffer[8], readBuffer[9], readBuffer[10], readBuffer[11], readBuffer[12]});
			bytesPerSampMark = (byte) (Compression.BYTES_PER_SAMPLE_MARKER_2 - (nullStatus & 0x01) - (nullStatus >>1 & 0x01) - (nullStatus >>2 & 0x01) - (nullStatus >>3 & 0x01) - (nullStatus >>4 & 0x01) - (nullStatus >>5 & 0x01) - (nullStatus >>6 & 0x01));
			outlierSectionLocation = PARAMETER_SECTION_BYTES + numMarkers * bytesPerSampMark;
			numBytesOfOutOfRangeValues = Compression.bytesToInt(new byte[]{readBuffer[outlierSectionLocation], readBuffer[outlierSectionLocation+1], readBuffer[outlierSectionLocation+2], readBuffer[outlierSectionLocation+3]});
			if (numBytesOfOutOfRangeValues>0) {
				outOfRangeValues = (Hashtable<String, Float>) Compression.bytesToObj(readBuffer, outlierSectionLocation+4, numBytesOfOutOfRangeValues);
			}
//			System.out.print("\t"+(new Date().getTime() - time));

			indexStart = PARAMETER_SECTION_BYTES;
			index = indexStart;
//			time = new Date().getTime();
			if (((nullStatus>>NULLSTATUS_GC_LOCATION) & 0x01) != 1) {
				if (loadGC) {
					gcs = new float[numMarkers];
					for (int j=0; j<numMarkers; j++) {
						gcs[j] = Compression.reducedPrecisionGcBafGetFloat2(new byte[] {readBuffer[index], readBuffer[index + 1]});
						index += bytesPerSampMark;
					}
				}
				indexStart += 2;
			}
			index = indexStart;
			if (((nullStatus>>NULLSTATUS_X_LOCATION) & 0x01) != 1) {
				if (loadXY) {
					xs = new float[numMarkers];
					for (int j=0; j<numMarkers; j++) {
						xs[j] = Compression.reducedPrecisionXYGetFloat2(new byte[] {readBuffer[index], readBuffer[index + 1]});
						if (xs[j]==Compression.REDUCED_PRECISION_XY_OUT_OF_RANGE_FLOAT) {
//							xs[j] = outOfRangeValues.get(sampleName+"\t"+allMarkersProj[j]+"\tx");
							xs[j] = outOfRangeValues.get(j + "\tx");
						}
						index += bytesPerSampMark;
					}
				}
				indexStart += 2;
			}
			index = indexStart;
			if (((nullStatus>>NULLSTATUS_Y_LOCATION) & 0x01) != 1) {
				if (loadXY) {
					ys = new float[numMarkers];
					for (int j=0; j<numMarkers; j++) {
						ys[j] = Compression.reducedPrecisionXYGetFloat2(new byte[] {readBuffer[index], readBuffer[index + 1]});
						if (ys[j]==Compression.REDUCED_PRECISION_XY_OUT_OF_RANGE_FLOAT) {
//							ys[j] = outOfRangeValues.get(sampleName+"\t"+allMarkersProj[j]+"\ty");
							ys[j] = outOfRangeValues.get(j + "\ty");
						}
						index += bytesPerSampMark;
					}
					
				}
				indexStart += 2;
			}
			index = indexStart;
			if (((nullStatus>>NULLSTATUS_BAF_LOCATION) & 0x01) != 1) {
				if (loadBAF) {
					bafs = new float[numMarkers];
					for (int j=0; j<numMarkers; j++) {
						bafs[j] = Compression.reducedPrecisionGcBafGetFloat2(new byte[] {readBuffer[index], readBuffer[index + 1]});
						index += bytesPerSampMark;
					}
				}
				indexStart += 2;
			}
			index = indexStart;
			if (((nullStatus>>NULLSTATUS_LRR_LOCATION) & 0x01) != 1) {
				if (loadLRR) {
					lrrs = new float[numMarkers];
					for (int j=0; j<numMarkers; j++) {
						lrrs[j] = Compression.reducedPrecisionLrrGetFloat2(new byte[] {readBuffer[index], readBuffer[index + 1], readBuffer[index + 2]});
						if (lrrs[j] == Compression.REDUCED_PRECISION_LRR_OUT_OF_RANGE_LRR_FLOAT) {
//							lrrs[j] = outOfRangeValues.get(sampleName+"\t"+allMarkersProj[j]+"\tlrr");
							lrrs[j] = outOfRangeValues.get(j + "\tlrr");
						}
						index += bytesPerSampMark;
					}
				}
				indexStart += 3;
			}
			index = indexStart;
			if ((((nullStatus>>NULLSTATUS_ABGENOTYPE_LOCATION) & 0x01) != 1 || ((nullStatus>>NULLSTATUS_FOWARDGENOTYPE_LOCATION) & 0x01) != 1) && loadAbOrForwardGenotypes) {
				abGenotypes = new byte[numMarkers];
				fwdGenotypes = new byte[numMarkers];
				for (int j=0; j<numMarkers; j++) {
					genoTypeTmp = Compression.reducedPrecisionGenotypeGetTypes(readBuffer[index]);
					abGenotypes[j] = genoTypeTmp[0];
					fwdGenotypes[j] = genoTypeTmp[1];
					index += bytesPerSampMark;
				}
			}
//			System.out.print("\t"+(new Date().getTime() - time));
//			time = new Date().getTime();
			result = new Sample(sampleName, fingerPrint, gcs, xs, ys, bafs, lrrs, fwdGenotypes, abGenotypes);
//			System.out.println("\t"+(new Date().getTime() - time));
		} catch (IOException e) {
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
		}
//		System.out.println("Finished loadToMarkerData from RAF in: "+(new Date().getTime() - time1)+" ms");
		return result;
	}

	/**
	 * Load data from Random Access Files organized by samples.
	 */
//	public static Sample loadFromRandomAccessFile(String filename, boolean jar, String sampleName, String[] allMarkersProj, long fingerPrint) {
	public static Sample loadFromRandomAccessFile(String filename, boolean jar) {
		Sample result = null;
		int numMarkers;
		RandomAccessFile file;
		String sampleName;
		long fingerPrint;
        float[] gcs, xs, ys, lrrs, bafs;
        byte[] abGenotypes, fwdGenotypes, readBuffer;
        byte[] genoTypeTmp;
        int index;
        int numBytesOfOutOfRangeValues;
        Hashtable<String, Float> outOfRangeValues = null;
        int outlierSectionLocation;
//      long time;

//		What would happen if not all the markers have exactly the same samples? In SQL this will result in shift the rows and will affect the matching of the data.
        sampleName = ext.rootOf(filename);
//		long time1 = new Date().getTime();
		try {
//			time = new Date().getTime();
			file = new RandomAccessFile(filename, "r");
			readBuffer = new byte[(int) file.length()];	//numMarkers * BYTES_PER_SAMPLE_MARKER_2
			file.read(readBuffer);
			file.close();
			numMarkers = Compression.bytesToInt(new byte[]{readBuffer[0], readBuffer[1], readBuffer[2], readBuffer[3]});
			fingerPrint = Compression.bytesToLong(new byte[]{readBuffer[4], readBuffer[5], readBuffer[6], readBuffer[7], readBuffer[8], readBuffer[9], readBuffer[10], readBuffer[11]});
			outlierSectionLocation = 12 + numMarkers * Compression.BYTES_PER_SAMPLE_MARKER_2;
			numBytesOfOutOfRangeValues = Compression.bytesToInt(new byte[]{readBuffer[outlierSectionLocation], readBuffer[outlierSectionLocation+1], readBuffer[outlierSectionLocation+2], readBuffer[outlierSectionLocation+3]});
			if (numBytesOfOutOfRangeValues>0) {
				outOfRangeValues = (Hashtable<String, Float>) Compression.bytesToObj(readBuffer, outlierSectionLocation+4, numBytesOfOutOfRangeValues);
			}
//			System.out.print("\t"+(new Date().getTime() - time));

			index = 12;
			gcs = new float[numMarkers];
			xs = new float[numMarkers];
			ys = new float[numMarkers];
			lrrs = new float[numMarkers];
			bafs = new float[numMarkers];
			abGenotypes = new byte[numMarkers];
			fwdGenotypes = new byte[numMarkers];
//			time = new Date().getTime();
			for (int j=0; j<numMarkers; j++) {
				gcs[j] = Compression.reducedPrecisionGcBafGetFloat2(new byte[] {readBuffer[index], readBuffer[index + 1]});
				xs[j] = Compression.reducedPrecisionXYGetFloat2(new byte[] {readBuffer[index + 2], readBuffer[index + 3]});
				if (xs[j]==-1) {
//					xs[j] = outOfRangeValues.get(sampleName+"\t"+allMarkersProj[j]+"\tx");
					xs[j] = outOfRangeValues.get(sampleName + "\t" + j + "\tx");
				}
				ys[j] = Compression.reducedPrecisionXYGetFloat2(new byte[] {readBuffer[index + 4], readBuffer[index + 5]});
				if (ys[j]==-1) {
//					ys[j] = outOfRangeValues.get(sampleName+"\t"+allMarkersProj[j]+"\ty");
					ys[j] = outOfRangeValues.get(sampleName + "\t" + j + "\ty");
				}
				bafs[j] = Compression.reducedPrecisionGcBafGetFloat2(new byte[] {readBuffer[index + 6], readBuffer[index + 7]});
				lrrs[j] = Compression.reducedPrecisionLrrGetFloat2(new byte[] {readBuffer[index + 8], readBuffer[index + 9], readBuffer[index + 10]});
				if (lrrs[j] == Compression.REDUCED_PRECISION_LRR_OUT_OF_RANGE_LRR_FLOAT) {
//					lrrs[j] = outOfRangeValues.get(sampleName+"\t"+allMarkersProj[j]+"\tlrr");
					lrrs[j] = outOfRangeValues.get(sampleName + "\t" + j + "\tlrr");
				}
				genoTypeTmp = Compression.reducedPrecisionGenotypeGetTypes(readBuffer[index + 11]);
				abGenotypes[j] = genoTypeTmp[0];
				fwdGenotypes[j] = genoTypeTmp[1];
				index += Compression.BYTES_PER_SAMPLE_MARKER_2;
			}
//			System.out.print("\t"+(new Date().getTime() - time));
//			time = new Date().getTime();
			result = new Sample(sampleName, fingerPrint, gcs, xs, ys, bafs, lrrs, fwdGenotypes, abGenotypes);
//			System.out.println("\t"+(new Date().getTime() - time));
		} catch (IOException e) {
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
		}
//		System.out.println("Finished loadToMarkerData from RAF in: "+(new Date().getTime() - time1)+" ms");
		return result;
	}

	/**
	 * Load data from Random Access Files organized by samples.
	 */
	public static Sample loadFromRandomAccessFileOnlyLrrBafGenotype2(String filename, boolean jar) {
		return loadFromRandomAccessFile2(filename, false, false, true, true, false, jar);
	}

//	public static Sample loadFromRandomAccessFileOnlyLrrBafGenotype(String filename, boolean jar, String sampleName, String[] allMarkersProj, long fingerPrint) {
	public static Sample loadFromRandomAccessFileOnlyLrrBafGenotype(String filename, boolean jar) {
		Sample result = null;
		int numMarkers;
		RandomAccessFile file;
		String sampleName;
		long fingerPrint;
        float[] lrrs, bafs;
        byte[] abGenotypes, readBuffer;
        int index;
        int numBytesOfOutOfRangeValues;
        Hashtable<String, Float> outOfRangeValues = null;
        int outlierSectionLocation;
//      long time;

//		What would happen if not all the markers have exactly the same samples? In SQL this will result in shift the rows and will affect the matching of the data.
        sampleName = ext.rootOf(filename);
//		long time1 = new Date().getTime();
		try {
//			time = new Date().getTime();
			file = new RandomAccessFile(filename, "r");
			readBuffer = new byte[(int) file.length()];	//numMarkers * BYTES_PER_SAMPLE_MARKER_2
			file.read(readBuffer);
			file.close();
			numMarkers = Compression.bytesToInt(new byte[]{readBuffer[0], readBuffer[1], readBuffer[2], readBuffer[3]});
			fingerPrint = Compression.bytesToLong(new byte[]{readBuffer[4], readBuffer[5], readBuffer[6], readBuffer[7], readBuffer[8], readBuffer[9], readBuffer[10], readBuffer[11]});
			outlierSectionLocation = 12 + numMarkers * Compression.BYTES_PER_SAMPLE_MARKER_2;
			numBytesOfOutOfRangeValues = Compression.bytesToInt(new byte[]{readBuffer[outlierSectionLocation], readBuffer[outlierSectionLocation+1], readBuffer[outlierSectionLocation+2], readBuffer[outlierSectionLocation+3]});
			if (numBytesOfOutOfRangeValues>0) {
				outOfRangeValues = (Hashtable<String, Float>) Compression.bytesToObj(readBuffer, outlierSectionLocation+4, numBytesOfOutOfRangeValues);
			}
//			System.out.print("\t"+(new Date().getTime() - time));

			index = 12;
			lrrs = new float[numMarkers];
			bafs = new float[numMarkers];
			abGenotypes = new byte[numMarkers];
//			time = new Date().getTime();
			for (int j=0; j<numMarkers; j++) {
				bafs[j] = Compression.reducedPrecisionGcBafGetFloat2(new byte[] {readBuffer[index + 6], readBuffer[index + 7]});
				lrrs[j] = Compression.reducedPrecisionLrrGetFloat2(new byte[] {readBuffer[index + 8], readBuffer[index + 9], readBuffer[index + 10]});
				if (lrrs[j] == Compression.REDUCED_PRECISION_LRR_OUT_OF_RANGE_LRR_FLOAT) {
//					lrrs[j] = outOfRangeValues.get(sampleName+"\t"+allMarkersProj[j]+"\tlrr");
					lrrs[j] = outOfRangeValues.get(sampleName + "\t" + j + "\tlrr");
				}
				abGenotypes[j] = Compression.reducedPrecisionGenotypeGetTypes(readBuffer[index + 11])[0];
				index += Compression.BYTES_PER_SAMPLE_MARKER_2;
			}
//			System.out.print("\t"+(new Date().getTime() - time));
//			time = new Date().getTime();
			result = new Sample(sampleName, fingerPrint, null, null, null, bafs, lrrs, null, abGenotypes);
//			System.out.println("\t"+(new Date().getTime() - time));
		} catch (IOException e) {
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
		}
//		System.out.println("Finished loadToMarkerData from RAF in: "+(new Date().getTime() - time1)+" ms");
		return result;
	}

	
	public static MarkerData[] loadMarkerDataFromRandomAccessFile2(Project proj, String[] markerNames, boolean jar) {
		MarkerSet markerSet;
		String[] markerNamesProj, samplesProj;
		byte[] chrsProj;
		int[] positionsProj;
		int[] indeciesOfmarkersToLoad;

		samplesProj = proj.getSamples();
		markerNamesProj = proj.getMarkerNames();
		markerSet = proj.getMarkerSet();
		chrsProj = markerSet.getChrs();
		positionsProj = markerSet.getPositions();

		indeciesOfmarkersToLoad = new int[markerNames.length];
		for (int i=0; i<markerNames.length; i++) {
			for (int j=0; j<markerNamesProj.length; j++) {
				if (markerNames[i].equals(markerNamesProj[j])) {
					indeciesOfmarkersToLoad[i] = j;
				}
			}
		}
		//We don't need to sort the markerNames by the order they appear in the project.
		
		return loadMarkerDataFromRandomAccessFile2(markerNamesProj, chrsProj, positionsProj, samplesProj, indeciesOfmarkersToLoad, proj.getDir(Project.SAMPLE_DIRECTORY, true), true, true, true, true, true, true);
//		return loadMarkerDataFromRandomAccessFile3(markerNamesProj, chrsProj, positionsProj, samplesProj, indeciesOfmarkersToLoad, proj.getDir(Project.SAMPLE_DIRECTORY, true), true, true, true, true, true, true);
	}
	
	public static MarkerData[] loadMarkerDataFromRandomAccessFile2(String[] markerNamesProj, byte[] chrsProj, int[] positionsProj, String[] samplesProj, int[] markerIndeciesInProj, String sampleDataFileDir) {
		return loadMarkerDataFromRandomAccessFile2(markerNamesProj, chrsProj, positionsProj, samplesProj, markerIndeciesInProj, sampleDataFileDir, true, true, true, true, true, true);
	}
	
//	/**
//	 * Load MarkerData from Random Access Files organized by samples.
//	 */
//	public static MarkerData[] loadMarkerDataFromRandomAccessFile2(String[] markerNamesProj, byte[] chrsProj, int[] positionsProj, String[] samplesProj, int[] indeciesInProjOfMarkersToLoad, String sampleDataFileDir, boolean loadGC, boolean loadXY, boolean loadBAF, boolean loadLRR, boolean loadAbGenotype) {
//		MarkerData[] result = null;
//		int numMarkers, numSamplesProj;
//		RandomAccessFile file;
//		long fingerPrint;
//        float[][] gcs = null, xs = null, ys = null, lrrs = null, bafs = null;
//        byte[][] abGenotypes = null;
//        String[][] alleleMappings = null;
//        byte[] genoTypeTmp, readBuffer;
//        int index, indexStart;
//        int numBytesOfOutOfRangeValues;
//        Hashtable<String, Float> outOfRangeValues = null;
//        int outlierSectionLocation;
//        byte nullStatus;
//        byte bytesPerSampMark;
//        long fingerprintForMarkerData;
////      long time;
//
////		What would happen if not all the markers have exactly the same samples? In SQL this will result in shift the rows and will affect the matching of the data.
//        numSamplesProj = samplesProj.length;
//        numMarkers = indeciesInProjOfMarkersToLoad.length;
//		fingerprintForMarkerData = MarkerSet.fingerprint(samplesProj);
//
//		long time1 = new Date().getTime();
//        for (int i=0; i<numSamplesProj; i++) {
//    		try {
//    			file = new RandomAccessFile(sampleDataFileDir + samplesProj[i] + Sample.SAMPLE_DATA_FILE_EXTENSION, "r");
//    			readBuffer = new byte[(int) file.length()];	//numMarkers * BYTES_PER_SAMPLE_MARKER_2
//    			file.read(readBuffer);
//    			file.close();
////    			numMarkersProj = Compression.bytesToInt(new byte[]{readBuffer[0], readBuffer[1], readBuffer[2], readBuffer[3]});
//    			nullStatus = readBuffer[4];
//    			fingerPrint = Compression.bytesToLong(new byte[]{readBuffer[5], readBuffer[6], readBuffer[7], readBuffer[8], readBuffer[9], readBuffer[10], readBuffer[11], readBuffer[12]});
//    			bytesPerSampMark = (byte) (Compression.BYTES_PER_SAMPLE_MARKER_2 - (nullStatus & 0x01) - (nullStatus >>1 & 0x01) - (nullStatus >>2 & 0x01) - (nullStatus >>3 & 0x01) - (nullStatus >>4 & 0x01) - (nullStatus >>5 & 0x01) - (nullStatus >>6 & 0x01));
//    			outlierSectionLocation = PARAMETER_SECTION_BYTES + bytesPerSampMark * Compression.bytesToInt(new byte[]{readBuffer[0], readBuffer[1], readBuffer[2], readBuffer[3]});
//    			numBytesOfOutOfRangeValues = Compression.bytesToInt(new byte[]{readBuffer[outlierSectionLocation], readBuffer[outlierSectionLocation+1], readBuffer[outlierSectionLocation+2], readBuffer[outlierSectionLocation+3]});
//    			if (numBytesOfOutOfRangeValues>0) {
//    				outOfRangeValues = (Hashtable<String, Float>) Compression.bytesToObj(readBuffer, outlierSectionLocation+4, outlierSectionLocation + numBytesOfOutOfRangeValues);
//    			}
//
//    			indexStart = PARAMETER_SECTION_BYTES;
//    			index = indexStart;
////    			time = new Date().getTime();
//    			if (((nullStatus>>NULLSTATUS_GC_LOCATION) & 0x01) != 1) {
//    				if (loadGC) {
//						gcs = new float[numMarkers][numSamplesProj];
//						for (int j=0; j<numMarkers; j++) {
//							index = indexStart + indeciesInProjOfMarkersToLoad[j] * bytesPerSampMark;
//							gcs[j][i] = Compression.reducedPrecisionGcBafGetFloat2(new byte[] {readBuffer[index], readBuffer[index + 1]});
//						}
//    				}
//    				indexStart += 2;
//    			}
//    			index = indexStart;
//    			if (((nullStatus>>NULLSTATUS_X_LOCATION) & 0x01) != 1) {
//    				if (loadXY) {
//						xs = new float[numMarkers][numSamplesProj];
//						for (int j=0; j<numMarkers; j++) {
//							index = indexStart + indeciesInProjOfMarkersToLoad[j] * bytesPerSampMark;
//							xs[j][i] = Compression.reducedPrecisionXYGetFloat2(new byte[] {readBuffer[index], readBuffer[index + 1]});
//							if (xs[j][i]==Compression.REDUCED_PRECISION_XY_OUT_OF_RANGE_FLOAT) {
//								xs[j][i] = outOfRangeValues.get(indeciesInProjOfMarkersToLoad[j] + "\tx");
//							}
//						}
//    				}
//    				indexStart += 2;
//    			}
//    			index = indexStart;
//    			if (((nullStatus>>NULLSTATUS_Y_LOCATION) & 0x01) != 1) {
//    				if (loadXY) {
//						ys = new float[numMarkers][numSamplesProj];
//						for (int j=0; j<numMarkers; j++) {
//							index = indexStart + indeciesInProjOfMarkersToLoad[j] * bytesPerSampMark;
//							ys[j][i] = Compression.reducedPrecisionXYGetFloat2(new byte[] {readBuffer[index], readBuffer[index + 1]});
//							if (ys[j][i]==Compression.REDUCED_PRECISION_XY_OUT_OF_RANGE_FLOAT) {
//								ys[j][i] = outOfRangeValues.get(indeciesInProjOfMarkersToLoad[j] + "\ty");
//							}
//						}
//    				}
//    				indexStart += 2;
//    			}
//    			index = indexStart;
//    			if (((nullStatus>>NULLSTATUS_BAF_LOCATION) & 0x01) != 1) {
//    				if (loadBAF) {
//						bafs = new float[numMarkers][numSamplesProj];
//						for (int j=0; j<numMarkers; j++) {
//							index = indexStart + indeciesInProjOfMarkersToLoad[j] * bytesPerSampMark;
//							bafs[j][i] = Compression.reducedPrecisionGcBafGetFloat2(new byte[] {readBuffer[index], readBuffer[index + 1]});
//						}
//    				}
//    				indexStart += 2;
//    			}
//    			index = indexStart;
//    			if (((nullStatus>>NULLSTATUS_LRR_LOCATION) & 0x01) != 1) {
//    				if (loadLRR) {
//						lrrs = new float[numMarkers][numSamplesProj];
//						for (int j=0; j<numMarkers; j++) {
//							index = indexStart + indeciesInProjOfMarkersToLoad[j] * bytesPerSampMark;
//							lrrs[j][i] = Compression.reducedPrecisionLrrGetFloat2(new byte[] {readBuffer[index], readBuffer[index + 1], readBuffer[index + 2]});
//							if (lrrs[j][i] == Compression.REDUCED_PRECISION_LRR_OUT_OF_RANGE_LRR_FLOAT) {
//								lrrs[j][i] = outOfRangeValues.get(indeciesInProjOfMarkersToLoad[j] + "\tlrr");
//							}
//						}
//    				}
//    				indexStart += 3;
//    			}
//    			index = indexStart;
//    			if ((((nullStatus>>NULLSTATUS_ABGENOTYPE_LOCATION) & 0x01) != 1 || ((nullStatus>>NULLSTATUS_FOWARDGENOTYPE_LOCATION) & 0x01) != 1) & loadAbGenotype) {
//    				abGenotypes = new byte[numMarkers][numSamplesProj];
//    				alleleMappings = new String[numMarkers][numSamplesProj];
//    				for (int j=0; j<numMarkers; j++) {
//						index = indexStart + indeciesInProjOfMarkersToLoad[j] * bytesPerSampMark;
//    					genoTypeTmp = Compression.reducedPrecisionGenotypeGetTypes(readBuffer[index]);
//    					abGenotypes[j][i] = genoTypeTmp[0];
//    					alleleMappings[j][i] = ALLELE_PAIRS[genoTypeTmp[1]];
//    				}
//    			}
//
//    		} catch (IOException e) {
//    			e.printStackTrace();
//    		} catch (ClassNotFoundException e) {
//    			e.printStackTrace();
//    		}
//        }
//        for (int i=0; i<numMarkers; i++) {
//    	        result[i] = new MarkerData(markerNamesProj[indeciesInProjOfMarkersToLoad[i]], chrsProj[indeciesInProjOfMarkersToLoad[i]], positionsProj[indeciesInProjOfMarkersToLoad[i]], fingerprintForMarkerData, gcs[i], new float[] {0}, new float[] {0}, xs[i], ys[i], new float[] {0}, new float[] {0}, bafs[i], lrrs[i], abGenotypes[i], alleleMappings[i]);
//        }
//        System.out.println("Total time loading MarkerData array from Sample Data Files: " + ((new Date().getTime() - time1)/1000) + "sec.");
//		return result;
//	}

	/**
	 * Load MarkerData from Random Access Files organized by samples.
	 */
	public static MarkerData[] loadMarkerDataFromRandomAccessFile2(String[] markerNamesProj, byte[] chrsProj, int[] positionsProj, String[] samplesProj, int[] indeciesInProjOfMarkersToLoad, String sampleDataFileDir, boolean loadGC, boolean loadXY, boolean loadBAF, boolean loadLRR, boolean loadAbGenotype, boolean loadAlleleMapping) {
		MarkerData[] result = null;
		int numMarkers, numSamplesProj;
		RandomAccessFile file;
        float[][] gcs = null, xs = null, ys = null, lrrs = null, bafs = null;
        byte[][] abGenotypes = null;
        String[][] alleleMappings = null;
        byte[] readBuffer;
        byte nullStatus;
        long fingerprintForMarkerData;
        long time1, time2;

//		What would happen if not all the markers have exactly the same samples? In SQL this will result in shift the rows and will affect the matching of the data.
        numSamplesProj = samplesProj.length;
        numMarkers = indeciesInProjOfMarkersToLoad.length;

		time1 = new Date().getTime();
		try {
			file = new RandomAccessFile(sampleDataFileDir + samplesProj[0] + Sample.SAMPLE_DATA_FILE_EXTENSION, "r");
			readBuffer = new byte[PARAMETER_SECTION_BYTES];	//numMarkers * BYTES_PER_SAMPLE_MARKER_2
			file.read(readBuffer);
			file.close();
//			fingerPrint = Compression.bytesToLong(new byte[]{readBuffer[5], readBuffer[6], readBuffer[7], readBuffer[8], readBuffer[9], readBuffer[10], readBuffer[11], readBuffer[12]});
			nullStatus = readBuffer[4];
			if (((nullStatus>>NULLSTATUS_GC_LOCATION) & 0x01) != 1 && loadGC) {
				gcs = new float[numMarkers][numSamplesProj];
			}
			if (((nullStatus>>NULLSTATUS_X_LOCATION) & 0x01) != 1 && loadXY) {
				xs = new float[numMarkers][numSamplesProj];
			}
			if (((nullStatus>>NULLSTATUS_Y_LOCATION) & 0x01) != 1 && loadXY) {
				ys = new float[numMarkers][numSamplesProj];
			}
			if (((nullStatus>>NULLSTATUS_BAF_LOCATION) & 0x01) != 1 && loadBAF) {
				bafs = new float[numMarkers][numSamplesProj];
			}
			if (((nullStatus>>NULLSTATUS_LRR_LOCATION) & 0x01) != 1 && loadLRR) {
				lrrs = new float[numMarkers][numSamplesProj];
			}
			if (((nullStatus>>NULLSTATUS_ABGENOTYPE_LOCATION) & 0x01) != 1 && loadAbGenotype) {
				abGenotypes = new byte[numMarkers][numSamplesProj];
			}
			if (((nullStatus>>NULLSTATUS_FOWARDGENOTYPE_LOCATION) & 0x01) != 1 & loadAlleleMapping) {
				alleleMappings = new String[numMarkers][numSamplesProj];
			}
    		time2 = new Date().getTime();
	        for (int i=0; i<numSamplesProj; i++) {
	        	loadMarkersByReadWholeRaf(samplesProj[i], i, nullStatus, indeciesInProjOfMarkersToLoad, sampleDataFileDir, gcs, xs, ys, lrrs, bafs, abGenotypes, alleleMappings);
//	        	loadMarkersByReadPartRaf(samplesProj[i], i, nullStatus, indeciesInProjOfMarkersToLoad, sampleDataFileDir, gcs, xs, ys, lrrs, bafs, abGenotypes, alleleMappings);
	        }
	        System.out.println("Loaded " + numMarkers + " markers in " + (new Date().getTime() - time2) + "ms.");
	        result = new MarkerData[numMarkers];
			fingerprintForMarkerData = MarkerSet.fingerprint(samplesProj);
	        for (int i=0; i<numMarkers; i++) {
	    	        result[i] = new MarkerData(markerNamesProj[indeciesInProjOfMarkersToLoad[i]], chrsProj[indeciesInProjOfMarkersToLoad[i]], positionsProj[indeciesInProjOfMarkersToLoad[i]], fingerprintForMarkerData, gcs[i], new float[] {0}, new float[] {0}, xs[i], ys[i], new float[] {0}, new float[] {0}, bafs[i], lrrs[i], abGenotypes[i], alleleMappings[i]);
	        }
		} catch (IOException e) {
			e.printStackTrace();
		}
        System.out.println("Total time loading MarkerData array from Sample Data Files: " + ((new Date().getTime() - time1)/1000) + "sec.");
		return result;
	}

	public static void loadMarkersByReadWholeRaf(String sampleToLoad, int i, byte nullStatus, int[] markersToLoadIndeciesInProj, String sampleDataFileDir, float[][] gcs, float[][] xs, float[][] ys, float[][] lrrs, float[][] bafs, byte[][] abGenotypes, String[][] alleleMappings) {
		int numMarkers;
		RandomAccessFile file;
		byte[] genoTypeTmp, readBuffer;
        int index, indexStart;
        int numBytesOfOutOfRangeValues;
        Hashtable<String, Float> outOfRangeValues = null;
        int outlierSectionLocation;
        byte bytesPerSampMark;
        long fingerprintInSampleDataFile;

        numMarkers = markersToLoadIndeciesInProj.length;

		try {
			file = new RandomAccessFile(sampleDataFileDir + sampleToLoad + Sample.SAMPLE_DATA_FILE_EXTENSION, "r");
			readBuffer = new byte[(int) file.length()];	//numMarkers * BYTES_PER_SAMPLE_MARKER_2
			file.read(readBuffer);
			file.close();

			fingerprintInSampleDataFile = Compression.bytesToLong(new byte[]{readBuffer[5], readBuffer[6], readBuffer[7], readBuffer[8], readBuffer[9], readBuffer[10], readBuffer[11], readBuffer[12]});
			bytesPerSampMark = (byte) (Compression.BYTES_PER_SAMPLE_MARKER_2 - (nullStatus & 0x01) - (nullStatus >>1 & 0x01) - (nullStatus >>2 & 0x01) - (nullStatus >>3 & 0x01) - (nullStatus >>4 & 0x01) - (nullStatus >>5 & 0x01) - (nullStatus >>6 & 0x01));
			outlierSectionLocation = PARAMETER_SECTION_BYTES + bytesPerSampMark * Compression.bytesToInt(new byte[]{readBuffer[0], readBuffer[1], readBuffer[2], readBuffer[3]});
			numBytesOfOutOfRangeValues = Compression.bytesToInt(new byte[]{readBuffer[outlierSectionLocation], readBuffer[outlierSectionLocation+1], readBuffer[outlierSectionLocation+2], readBuffer[outlierSectionLocation+3]});
			if (numBytesOfOutOfRangeValues>0) {
				outOfRangeValues = (Hashtable<String, Float>) Compression.bytesToObj(readBuffer, outlierSectionLocation+4, numBytesOfOutOfRangeValues);
			}

			indexStart = PARAMETER_SECTION_BYTES;
			index = indexStart;
			if (((nullStatus>>NULLSTATUS_GC_LOCATION) & 0x01) != 1) {
				if (gcs != null) {
					for (int j=0; j<numMarkers; j++) {
						index = indexStart + markersToLoadIndeciesInProj[j] * bytesPerSampMark;
						gcs[j][i] = Compression.reducedPrecisionGcBafGetFloat2(new byte[] {readBuffer[index], readBuffer[index + 1]});
					}
				}
				indexStart += 2;
			}
			index = indexStart;
			if (((nullStatus>>NULLSTATUS_X_LOCATION) & 0x01) != 1) {
				if (xs != null) {
					for (int j=0; j<numMarkers; j++) {
						index = indexStart + markersToLoadIndeciesInProj[j] * bytesPerSampMark;
						xs[j][i] = Compression.reducedPrecisionXYGetFloat2(new byte[] {readBuffer[index], readBuffer[index + 1]});
						if (xs[j][i]==Compression.REDUCED_PRECISION_XY_OUT_OF_RANGE_FLOAT) {
							xs[j][i] = outOfRangeValues.get(markersToLoadIndeciesInProj[j] + "\tx");
						}
					}
				}
				indexStart += 2;
			}
			index = indexStart;
			if (((nullStatus>>NULLSTATUS_Y_LOCATION) & 0x01) != 1) {
				if (ys != null) {
					for (int j=0; j<numMarkers; j++) {
						index = indexStart + markersToLoadIndeciesInProj[j] * bytesPerSampMark;
						ys[j][i] = Compression.reducedPrecisionXYGetFloat2(new byte[] {readBuffer[index], readBuffer[index + 1]});
						if (ys[j][i]==Compression.REDUCED_PRECISION_XY_OUT_OF_RANGE_FLOAT) {
							ys[j][i] = outOfRangeValues.get(markersToLoadIndeciesInProj[j] + "\ty");
						}
					}
				}
				indexStart += 2;
			}
			index = indexStart;
			if (((nullStatus>>NULLSTATUS_BAF_LOCATION) & 0x01) != 1) {
				if (bafs != null) {
					for (int j=0; j<numMarkers; j++) {
						index = indexStart + markersToLoadIndeciesInProj[j] * bytesPerSampMark;
						bafs[j][i] = Compression.reducedPrecisionGcBafGetFloat2(new byte[] {readBuffer[index], readBuffer[index + 1]});
					}
				}
				indexStart += 2;
			}
			index = indexStart;
			if (((nullStatus>>NULLSTATUS_LRR_LOCATION) & 0x01) != 1) {
				if (lrrs != null) {
					for (int j=0; j<numMarkers; j++) {
						try {
							index = indexStart + markersToLoadIndeciesInProj[j] * bytesPerSampMark;
							lrrs[j][i] = Compression.reducedPrecisionLrrGetFloat2(new byte[] {readBuffer[index], readBuffer[index + 1], readBuffer[index + 2]});
							if (lrrs[j][i] == Compression.REDUCED_PRECISION_LRR_OUT_OF_RANGE_LRR_FLOAT) {
								lrrs[j][i] = outOfRangeValues.get(markersToLoadIndeciesInProj[j] + "\tlrr");
							}
						} catch (Exception e) {
							System.out.println("hi: ");
						}
					}
				}
				indexStart += 3;
			}
			index = indexStart;
			if ((((nullStatus>>NULLSTATUS_ABGENOTYPE_LOCATION) & 0x01) != 1 || ((nullStatus>>NULLSTATUS_FOWARDGENOTYPE_LOCATION) & 0x01) != 1)) {
				if ((abGenotypes != null) || (alleleMappings != null)) {
					for (int j=0; j<numMarkers; j++) {
						index = indexStart + markersToLoadIndeciesInProj[j] * bytesPerSampMark;
						genoTypeTmp = Compression.reducedPrecisionGenotypeGetTypes(readBuffer[index]);
						if (abGenotypes != null) {
							abGenotypes[j][i] = genoTypeTmp[0];
						}
						if (alleleMappings != null) {
							alleleMappings[j][i] = ALLELE_PAIRS[genoTypeTmp[1]];
						}
					}
				}
			}

		} catch (IOException e) {
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Load MarkerData from Random Access Files organized by samples.
	 */
	public static MarkerData[] loadMarkerDataFromRandomAccessFile3(String[] markerNamesProj, byte[] chrsProj, int[] positionsProj, String[] samplesProj, int[] markersToLoadProjInd, String sampleDataFileDir, boolean loadGC, boolean loadXY, boolean loadBAF, boolean loadLRR, boolean loadAbGenotype, boolean loadAlleleMapping) {
		MarkerData[] result = null;
		int numMarkersToLoad, numSamplesProj;
		RandomAccessFile[] files;
		RandomAccessFile currentFile;
        float[][] gcs = null, xs = null, ys = null, lrrs = null, bafs = null;
        byte[][] abGenotypes = null;
        String[][] alleleMappings = null;
        byte[] readBuffer;
        byte nullStatus;
        long fingerprintForMarkerData;
        int bytesPerSampMark;
        long seekLocation;
        long outlierSectionLocations;
        Hashtable<String, Float>[] outOfRangeValues;
        TreeMap<Integer, Integer> sorting;
        int[] markersToLoadIndProjSorted;
        long time1, time2;

//		What would happen if not all the markers have exactly the same samples? In SQL this will result in shift the rows and will affect the matching of the data.
        numMarkersToLoad = markersToLoadProjInd.length;
        if (markersToLoadProjInd==null || numMarkersToLoad==0) {
        	return null;
        }

        sorting = new TreeMap<Integer, Integer>();
        for (int i=0; i<numMarkersToLoad; i++) {
        	sorting.put(markersToLoadProjInd[i], i);
        }
        markersToLoadIndProjSorted = new int[numMarkersToLoad];
        for (int i=0; i<markersToLoadProjInd.length; i++) {
        	markersToLoadIndProjSorted[i] = sorting.get(markersToLoadProjInd[i]);
        }
        
        numSamplesProj = samplesProj.length;
        files = new RandomAccessFile[numSamplesProj];
        outOfRangeValues = new Hashtable[numSamplesProj];

		time1 = new Date().getTime();
		try {
    		time2 = new Date().getTime();
			files[0] = new RandomAccessFile(sampleDataFileDir + samplesProj[0] + Sample.SAMPLE_DATA_FILE_EXTENSION, "r");
			readBuffer = new byte[PARAMETER_SECTION_BYTES];	//numMarkers * BYTES_PER_SAMPLE_MARKER_2
			files[0].read(readBuffer);
//			fingerPrint = Compression.bytesToLong(new byte[]{readBuffer[5], readBuffer[6], readBuffer[7], readBuffer[8], readBuffer[9], readBuffer[10], readBuffer[11], readBuffer[12]});
			nullStatus = readBuffer[4];
			bytesPerSampMark = (byte) (Compression.BYTES_PER_SAMPLE_MARKER_2 - (nullStatus & 0x01) - (nullStatus >>1 & 0x01) - (nullStatus >>2 & 0x01) - (nullStatus >>3 & 0x01) - (nullStatus >>4 & 0x01) - (nullStatus >>5 & 0x01) - (nullStatus >>6 & 0x01));
			outlierSectionLocations = PARAMETER_SECTION_BYTES + bytesPerSampMark * Compression.bytesToInt(readBuffer, PARAMETER_SECTION_NUMMARK_LOC);
	        System.out.println("Read in parameters in " + (new Date().getTime() - time2) + "ms.");

    		time2 = new Date().getTime();
			if (((nullStatus>>NULLSTATUS_GC_LOCATION) & 0x01) != 1 && loadGC) {
				gcs = new float[numMarkersToLoad][numSamplesProj];
			}
			if (((nullStatus>>NULLSTATUS_X_LOCATION) & 0x01) != 1 && loadXY) {
				xs = new float[numMarkersToLoad][numSamplesProj];
			}
			if (((nullStatus>>NULLSTATUS_Y_LOCATION) & 0x01) != 1 && loadXY) {
				ys = new float[numMarkersToLoad][numSamplesProj];
			}
			if (((nullStatus>>NULLSTATUS_BAF_LOCATION) & 0x01) != 1 && loadBAF) {
				bafs = new float[numMarkersToLoad][numSamplesProj];
			}
			if (((nullStatus>>NULLSTATUS_LRR_LOCATION) & 0x01) != 1 && loadLRR) {
				lrrs = new float[numMarkersToLoad][numSamplesProj];
			}
			if (((nullStatus>>NULLSTATUS_ABGENOTYPE_LOCATION) & 0x01) != 1 && loadAbGenotype) {
				abGenotypes = new byte[numMarkersToLoad][numSamplesProj];
			}
			if (((nullStatus>>NULLSTATUS_FOWARDGENOTYPE_LOCATION) & 0x01) != 1 & loadAlleleMapping) {
				alleleMappings = new String[numMarkersToLoad][numSamplesProj];
			}
	        System.out.println("Initialized variables in " + (new Date().getTime() - time2) + "ms.");





//    		time2 = new Date().getTime();
//			for (int i=1; i<numSamplesProj; i++) {
//				files[i] = new RandomAccessFile(sampleDataFileDir + samplesProj[i] + Sample.SAMPLE_DATA_FILE_EXTENSION, "r");
//    		}
//	        System.out.println("Opened all sample Files in " + (new Date().getTime() - time2) + "ms.");
//
//    		time2 = new Date().getTime();
//	        for (int i=0; i<numMarkersToLoad; i++) {
//	        	seekLocation = (long)PARAMETER_SECTION_BYTES + (long)markersToLoadProjInd[i] * (long)bytesPerSampMark;
//	        	loadMarkersByReadPartRaf3(files, nullStatus, seekLocation, outlierSectionLocations, outOfRangeValues, markersToLoadProjInd[markersToLoadIndProjSorted[i]], markersToLoadIndProjSorted[i], bytesPerSampMark, gcs, xs, ys, lrrs, bafs, abGenotypes, alleleMappings);
//	        }
//	        System.out.println("Loaded " + numMarkersToLoad + " markers in " + (new Date().getTime() - time2) + "ms.");
//
//	        for (int i=0; i<numSamplesProj; i++) {
//				files[i].close();
//    		}

	        for (int i=0; i<numMarkersToLoad; i++) {
	        	seekLocation = (long)PARAMETER_SECTION_BYTES + (long)markersToLoadProjInd[i] * (long)bytesPerSampMark;
	        	time2 = new Date().getTime();
				for (int j=1; j<numSamplesProj; j++) {
					currentFile = new RandomAccessFile(sampleDataFileDir + samplesProj[j] + Sample.SAMPLE_DATA_FILE_EXTENSION, "r");
		        	loadMarkersByReadPartRaf3a(currentFile, j, nullStatus, seekLocation, outlierSectionLocations, outOfRangeValues, markersToLoadProjInd[markersToLoadIndProjSorted[i]], markersToLoadIndProjSorted[i], bytesPerSampMark, gcs, xs, ys, lrrs, bafs, abGenotypes, alleleMappings);
		        	currentFile.close();
				}
				System.out.println("Loaded marker " + i + "in " + (new Date().getTime() - time2) + "ms.");
	        }





	        result = new MarkerData[numMarkersToLoad];
			fingerprintForMarkerData = MarkerSet.fingerprint(samplesProj);
	        for (int i=0; i<numMarkersToLoad; i++) {
	    	        result[i] = new MarkerData(markerNamesProj[markersToLoadProjInd[i]], chrsProj[markersToLoadProjInd[i]], positionsProj[markersToLoadProjInd[i]], fingerprintForMarkerData, gcs==null?null:gcs[i], new float[] {0}, new float[] {0}, xs==null?null:xs[i], ys==null?null:ys[i], new float[] {0}, new float[] {0}, bafs==null?null:bafs[i], lrrs==null?null:lrrs[i], abGenotypes==null?null:abGenotypes[i], alleleMappings==null?null:alleleMappings[i]);
	        }
		} catch (IOException e) {
			e.printStackTrace();
		}
        System.out.println("Total time loading MarkerData array from Sample Data Files: " + ((new Date().getTime() - time1)/1000) + "sec.");
		return result;
	}

	public static void loadMarkersByReadPartRaf3(RandomAccessFile[] files, byte nullStatus, long seekLocation, long outlierSectionLocations, Hashtable<String, Float>[] outOfRangeValues, int markersIndexProj, int markerIndexToLoad, int bytesPerSampMark, float[][] gcs, float[][] xs, float[][] ys, float[][] lrrs, float[][] bafs, byte[][] abGenotypes, String[][] alleleMappings) {
		byte[] genoTypeTmp, readBuffer;
        int index, indexStart;

		try {
			for (int i=0; i<files.length; i++) {
				readBuffer = new byte[bytesPerSampMark];
				files[i].seek(seekLocation);
				files[i].read(readBuffer);

				indexStart = 0;
				index = indexStart;
				if (((nullStatus>>NULLSTATUS_GC_LOCATION) & 0x01) != 1) {
					if (gcs != null) {
						gcs[markerIndexToLoad][i] = Compression.reducedPrecisionGcBafGetFloat2(new byte[] {readBuffer[index], readBuffer[index + 1]});
					}
					indexStart += 2;
				}
				index = indexStart;
				if (((nullStatus>>NULLSTATUS_X_LOCATION) & 0x01) != 1) {
					if (xs != null) {
						xs[markerIndexToLoad][i] = Compression.reducedPrecisionXYGetFloat2(new byte[] {readBuffer[index], readBuffer[index + 1]});
						if (xs[markerIndexToLoad][i]==Compression.REDUCED_PRECISION_XY_OUT_OF_RANGE_FLOAT) {
							if (outOfRangeValues[i] == null) {
								outOfRangeValues[i] = loadOutOfRangeValues(files[i], outlierSectionLocations);
							}
							xs[markerIndexToLoad][i] = outOfRangeValues[i].get(markersIndexProj + "\tx");
						}
					}
					indexStart += 2;
				}
				index = indexStart;
				if (((nullStatus>>NULLSTATUS_Y_LOCATION) & 0x01) != 1) {
					if (ys != null) {
						ys[markerIndexToLoad][i] = Compression.reducedPrecisionXYGetFloat2(new byte[] {readBuffer[index], readBuffer[index + 1]});
						if (ys[markerIndexToLoad][i]==Compression.REDUCED_PRECISION_XY_OUT_OF_RANGE_FLOAT) {
							if (outOfRangeValues[i] == null) {
								outOfRangeValues[i] = loadOutOfRangeValues(files[i], outlierSectionLocations);
							}
							ys[markerIndexToLoad][i] = outOfRangeValues[i].get(markersIndexProj + "\ty");
						}
					}
					indexStart += 2;
				}
				index = indexStart;
				if (((nullStatus>>NULLSTATUS_BAF_LOCATION) & 0x01) != 1) {
					if (bafs != null) {
							bafs[markerIndexToLoad][i] = Compression.reducedPrecisionGcBafGetFloat2(new byte[] {readBuffer[index], readBuffer[index + 1]});
						}
					indexStart += 2;
				}
				index = indexStart;
				if (((nullStatus>>NULLSTATUS_LRR_LOCATION) & 0x01) != 1) {
					if (lrrs != null) {
						lrrs[markerIndexToLoad][i] = Compression.reducedPrecisionLrrGetFloat2(new byte[] {readBuffer[index], readBuffer[index + 1], readBuffer[index + 2]});
						if (lrrs[markerIndexToLoad][i] == Compression.REDUCED_PRECISION_LRR_OUT_OF_RANGE_LRR_FLOAT) {
							if (outOfRangeValues[i] == null) {
								outOfRangeValues[i] = loadOutOfRangeValues(files[i], outlierSectionLocations);
							}
							lrrs[markerIndexToLoad][i] = outOfRangeValues[i].get(markersIndexProj + "\tlrr");
						}
					}
					indexStart += 3;
				}
				index = indexStart;
				if ((((nullStatus>>NULLSTATUS_ABGENOTYPE_LOCATION) & 0x01) != 1 || ((nullStatus>>NULLSTATUS_FOWARDGENOTYPE_LOCATION) & 0x01) != 1)) {
					if ((abGenotypes != null) || (alleleMappings != null)) {
						genoTypeTmp = Compression.reducedPrecisionGenotypeGetTypes(readBuffer[index]);
						if (abGenotypes != null) {
							abGenotypes[markerIndexToLoad][i] = genoTypeTmp[0];
						}
						if (alleleMappings != null) {
							alleleMappings[markerIndexToLoad][i] = ALLELE_PAIRS[genoTypeTmp[1]];
						}
					}
				}
	
			}

		} catch (IOException e) {
			e.printStackTrace();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static void loadMarkersByReadPartRaf3a(RandomAccessFile file, int sampleIndex, byte nullStatus, long seekLocation, long outlierSectionLocations, Hashtable<String, Float>[] outOfRangeValues, int markersIndexProj, int markerIndexToLoad, int bytesPerSampMark, float[][] gcs, float[][] xs, float[][] ys, float[][] lrrs, float[][] bafs, byte[][] abGenotypes, String[][] alleleMappings) {
		byte[] genoTypeTmp, readBuffer;
        int index, indexStart;

		try {
			readBuffer = new byte[bytesPerSampMark];
			file.seek(seekLocation);
			file.read(readBuffer);

			indexStart = 0;
			index = indexStart;
			if (((nullStatus>>NULLSTATUS_GC_LOCATION) & 0x01) != 1) {
				if (gcs != null) {
					gcs[markerIndexToLoad][sampleIndex] = Compression.reducedPrecisionGcBafGetFloat2(new byte[] {readBuffer[index], readBuffer[index + 1]});
				}
				indexStart += 2;
			}
			index = indexStart;
			if (((nullStatus>>NULLSTATUS_X_LOCATION) & 0x01) != 1) {
				if (xs != null) {
					xs[markerIndexToLoad][sampleIndex] = Compression.reducedPrecisionXYGetFloat2(new byte[] {readBuffer[index], readBuffer[index + 1]});
					if (xs[markerIndexToLoad][sampleIndex]==Compression.REDUCED_PRECISION_XY_OUT_OF_RANGE_FLOAT) {
						if (outOfRangeValues[sampleIndex] == null) {
							outOfRangeValues[sampleIndex] = loadOutOfRangeValues(file, outlierSectionLocations);
						}
						xs[markerIndexToLoad][sampleIndex] = outOfRangeValues[sampleIndex].get(markersIndexProj + "\tx");
					}
				}
				indexStart += 2;
			}
			index = indexStart;
			if (((nullStatus>>NULLSTATUS_Y_LOCATION) & 0x01) != 1) {
				if (ys != null) {
					ys[markerIndexToLoad][sampleIndex] = Compression.reducedPrecisionXYGetFloat2(new byte[] {readBuffer[index], readBuffer[index + 1]});
					if (ys[markerIndexToLoad][sampleIndex]==Compression.REDUCED_PRECISION_XY_OUT_OF_RANGE_FLOAT) {
						if (outOfRangeValues[sampleIndex] == null) {
							outOfRangeValues[sampleIndex] = loadOutOfRangeValues(file, outlierSectionLocations);
						}
						ys[markerIndexToLoad][sampleIndex] = outOfRangeValues[sampleIndex].get(markersIndexProj + "\ty");
					}
				}
				indexStart += 2;
			}
			index = indexStart;
			if (((nullStatus>>NULLSTATUS_BAF_LOCATION) & 0x01) != 1) {
				if (bafs != null) {
						bafs[markerIndexToLoad][sampleIndex] = Compression.reducedPrecisionGcBafGetFloat2(new byte[] {readBuffer[index], readBuffer[index + 1]});
					}
				indexStart += 2;
			}
			index = indexStart;
			if (((nullStatus>>NULLSTATUS_LRR_LOCATION) & 0x01) != 1) {
				if (lrrs != null) {
					lrrs[markerIndexToLoad][sampleIndex] = Compression.reducedPrecisionLrrGetFloat2(new byte[] {readBuffer[index], readBuffer[index + 1], readBuffer[index + 2]});
					if (lrrs[markerIndexToLoad][sampleIndex] == Compression.REDUCED_PRECISION_LRR_OUT_OF_RANGE_LRR_FLOAT) {
						if (outOfRangeValues[sampleIndex] == null) {
							outOfRangeValues[sampleIndex] = loadOutOfRangeValues(file, outlierSectionLocations);
						}
						lrrs[markerIndexToLoad][sampleIndex] = outOfRangeValues[sampleIndex].get(markersIndexProj + "\tlrr");
					}
				}
				indexStart += 3;
			}
			index = indexStart;
			if ((((nullStatus>>NULLSTATUS_ABGENOTYPE_LOCATION) & 0x01) != 1 || ((nullStatus>>NULLSTATUS_FOWARDGENOTYPE_LOCATION) & 0x01) != 1)) {
				if ((abGenotypes != null) || (alleleMappings != null)) {
					genoTypeTmp = Compression.reducedPrecisionGenotypeGetTypes(readBuffer[index]);
					if (abGenotypes != null) {
						abGenotypes[markerIndexToLoad][sampleIndex] = genoTypeTmp[0];
					}
					if (alleleMappings != null) {
						alleleMappings[markerIndexToLoad][sampleIndex] = ALLELE_PAIRS[genoTypeTmp[1]];
					}
				}
			}

		} catch (IOException e) {
			e.printStackTrace();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
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

	/**
	 * Load data from Random Access Files organized by samples.
	 */
	public static void testLoadTime(Project proj) {
		String[] samplesProj;
		RandomAccessFile file;
		byte[] readBuffer;
		long time1, time2;

		samplesProj = proj.getSamples();
		time1 = new Date().getTime();
		try {
//			for (int i=0; i<samplesProj.length; i++) {
			for (int i=0; i<100; i++) {
				file = new RandomAccessFile(proj.getDir(Project.SAMPLE_DIRECTORY, true) + samplesProj[i] + SAMPLE_DATA_FILE_EXTENSION, "r");
				readBuffer = new byte[(int) file.length()];	//numMarkers * BYTES_PER_SAMPLE_MARKER_2
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

	//	public static Sample loadFromRandomAccessFile3(String filename, boolean jar) {
	//		return loadFromRandomAccessFile2(filename, true, true, true, true, true, jar);
	//	}
		
	//	/**
	//	 * Load data from Random Access Files organized by samples.
	//	 */
	//	public static Sample loadFromRandomAccessFile3(String filename, boolean loadGC, boolean loadXY, boolean loadBAF, boolean loadLRR, boolean loadAbOrForwardGenotypes, boolean jar) {
	//		Sample result = null;
	//		int numMarkers;
	//		RandomAccessFile file;
	//		String sampleName;
	//		long fingerPrint;
	//        float[] gcs = null, xs = null, ys = null, lrrs = null, bafs = null;
	//        byte[] abGenotypes = null, fwdGenotypes = null, readBuffer;
	//        byte[] genoTypeTmp;
	//        int index, indexStart;
	//        int numBytesOfOutOfRangeValues;
	//        Hashtable<String, Float> outOfRangeValues = null;
	//        int outlierSectionLocation;
	//        byte nullStatus;
	//        byte bytesPerSampMark;
	////      long time;
	//
	////		What would happen if not all the markers have exactly the same samples? In SQL this will result in shift the rows and will affect the matching of the data.
	//        sampleName = ext.rootOf(filename);
	////		long time1 = new Date().getTime();
	//		try {
	////			time = new Date().getTime();
	//			file = new RandomAccessFile(filename, "r");
	//			readBuffer = new byte[(int) file.length()];	//numMarkers * BYTES_PER_SAMPLE_MARKER_2
	//			file.read(readBuffer);
	//			file.close();
	//			numMarkers = Compression.bytesToInt(new byte[]{readBuffer[0], readBuffer[1], readBuffer[2], readBuffer[3]});
	//			nullStatus = readBuffer[4];
	//			fingerPrint = Compression.bytesToLong(new byte[]{readBuffer[5], readBuffer[6], readBuffer[7], readBuffer[8], readBuffer[9], readBuffer[10], readBuffer[11], readBuffer[12]});
	//			bytesPerSampMark = (byte) (Compression.BYTES_PER_SAMPLE_MARKER_2 - (nullStatus & 0x01) - (nullStatus >>1 & 0x01) - (nullStatus >>2 & 0x01) - (nullStatus >>3 & 0x01) - (nullStatus >>4 & 0x01) - (nullStatus >>5 & 0x01) - (nullStatus >>6 & 0x01));
	//			outlierSectionLocation = PARAMETER_SECTION_BYTES + numMarkers * bytesPerSampMark;
	//			numBytesOfOutOfRangeValues = Compression.bytesToInt(new byte[]{readBuffer[outlierSectionLocation], readBuffer[outlierSectionLocation+1], readBuffer[outlierSectionLocation+2], readBuffer[outlierSectionLocation+3]});
	//			if (numBytesOfOutOfRangeValues>0) {
	//				outOfRangeValues = (Hashtable<String, Float>) Compression.bytesToObj(readBuffer, outlierSectionLocation+4, outlierSectionLocation + numBytesOfOutOfRangeValues);
	//			}
	////			System.out.print("\t"+(new Date().getTime() - time));
	//
	//			indexStart = PARAMETER_SECTION_BYTES;
	//			index = indexStart;
	////			time = new Date().getTime();
	//			if (((nullStatus>>NULLSTATUS_GC_LOCATION) & 0x01) != 1) {
	//				if (loadGC) {
	//					gcs = new float[numMarkers];
	//					for (int j=0; j<numMarkers; j++) {
	//						gcs[j] = Compression.reducedPrecisionGcBafGetFloat2(new byte[] {readBuffer[index], readBuffer[index + 1]});
	//						index += bytesPerSampMark;
	//					}
	//				}
	//				indexStart += 2;
	//			}
	//			index = indexStart;
	//			if (((nullStatus>>NULLSTATUS_X_LOCATION) & 0x01) != 1) {
	//				if (loadXY) {
	//					xs = new float[numMarkers];
	//					for (int j=0; j<numMarkers; j++) {
	//						xs[j] = Compression.reducedPrecisionXYGetFloat2(new byte[] {readBuffer[index], readBuffer[index + 1]});
	//						if (xs[j]==-1) {
	////							xs[j] = outOfRangeValues.get(sampleName+"\t"+allMarkersProj[j]+"\tx");
	//							xs[j] = outOfRangeValues.get(j + "\tx");
	//						}
	//						index += bytesPerSampMark;
	//					}
	//				}
	//				indexStart += 2;
	//			}
	//			index = indexStart;
	//			if (((nullStatus>>NULLSTATUS_Y_LOCATION) & 0x01) != 1) {
	//				if (loadXY) {
	//					ys = new float[numMarkers];
	//					for (int j=0; j<numMarkers; j++) {
	//						ys[j] = Compression.reducedPrecisionXYGetFloat2(new byte[] {readBuffer[index], readBuffer[index + 1]});
	//						if (ys[j]==-1) {
	////							ys[j] = outOfRangeValues.get(sampleName+"\t"+allMarkersProj[j]+"\ty");
	//							ys[j] = outOfRangeValues.get(j + "\ty");
	//						}
	//						index += bytesPerSampMark;
	//					}
	//					
	//				}
	//				indexStart += 2;
	//			}
	//			index = indexStart;
	//			if (((nullStatus>>NULLSTATUS_BAF_LOCATION) & 0x01) != 1) {
	//				if (loadBAF) {
	//					bafs = new float[numMarkers];
	//					for (int j=0; j<numMarkers; j++) {
	//						bafs[j] = Compression.reducedPrecisionGcBafGetFloat2(new byte[] {readBuffer[index], readBuffer[index + 1]});
	//						index += bytesPerSampMark;
	//					}
	//				}
	//				indexStart += 2;
	//			}
	//			index = indexStart;
	//			if (((nullStatus>>NULLSTATUS_LRR_LOCATION) & 0x01) != 1) {
	//				if (loadLRR) {
	//					lrrs = new float[numMarkers];
	//					for (int j=0; j<numMarkers; j++) {
	//						lrrs[j] = Compression.reducedPrecisionLrrGetFloat2(new byte[] {readBuffer[index], readBuffer[index + 1], readBuffer[index + 2]});
	//						if (lrrs[j] == Compression.REDUCED_PRECISION_LRR_OUT_OF_RANGE_LRR_FLOAT) {
	////							lrrs[j] = outOfRangeValues.get(sampleName+"\t"+allMarkersProj[j]+"\tlrr");
	//							lrrs[j] = outOfRangeValues.get(j + "\tlrr");
	//						}
	//						index += bytesPerSampMark;
	//					}
	//				}
	//				indexStart += 3;
	//			}
	//			index = indexStart;
	//			if ((((nullStatus>>NULLSTATUS_ABGENOTYPE_LOCATION) & 0x01) != 1 || ((nullStatus>>NULLSTATUS_FOWARDGENOTYPE_LOCATION) & 0x01) != 1) && loadAbOrForwardGenotypes) {
	//				abGenotypes = new byte[numMarkers];
	//				fwdGenotypes = new byte[numMarkers];
	//				for (int j=0; j<numMarkers; j++) {
	//					genoTypeTmp = Compression.reducedPrecisionGenotypeGetTypes(readBuffer[index]);
	//					abGenotypes[j] = genoTypeTmp[0];
	//					fwdGenotypes[j] = genoTypeTmp[1];
	//					index += bytesPerSampMark;
	//				}
	//			}
	////			System.out.print("\t"+(new Date().getTime() - time));
	////			time = new Date().getTime();
	//			result = new Sample(sampleName, fingerPrint, gcs, xs, ys, bafs, lrrs, fwdGenotypes, abGenotypes);
	////			System.out.println("\t"+(new Date().getTime() - time));
	//		} catch (IOException e) {
	//			e.printStackTrace();
	//		} catch (ClassNotFoundException e) {
	//			e.printStackTrace();
	//		}
	////		System.out.println("Finished loadToMarkerData from RAF in: "+(new Date().getTime() - time1)+" ms");
	//		return result;
	//	}
	
	//	public static byte[] loadByteStreamFromRandomAccessFile3(String filename, boolean jar, int beginMarkerIndex, int endMarkerIndex) {
	//		int numMarkers;
	//		RandomAccessFile file;
	//		byte[] parameterSection = null, readBuffer = null;
	//        int numBytesOfOutOfRangeValues;
	//        int outlierSectionLocation;
	//        byte nullStatus;
	//        byte bytesPerSampMark;
	//
	//        ext.rootOf(filename);
	//		try {
	//			parameterSection = new byte[PARAMETER_SECTION_BYTES];	//numMarkers * BYTES_PER_SAMPLE_MARKER_2
	//			file = new RandomAccessFile(filename, "r");
	//			file.read(parameterSection);
	//			numMarkers = Compression.bytesToInt(new byte[]{parameterSection[0], parameterSection[1], parameterSection[2], parameterSection[3]});
	////			numMarkers = endMarkerIndex - beginMarkerIndex + 1;
	//			nullStatus = parameterSection[4];
	//			Compression.bytesToLong(new byte[]{parameterSection[5], parameterSection[6], parameterSection[7], parameterSection[8], parameterSection[9], parameterSection[10], parameterSection[11], parameterSection[12]});
	//			bytesPerSampMark = (byte) (Compression.BYTES_PER_SAMPLE_MARKER_2 - (nullStatus & 0x01) - (nullStatus >>1 & 0x01) - (nullStatus >>2 & 0x01) - (nullStatus >>3 & 0x01) - (nullStatus >>4 & 0x01) - (nullStatus >>5 & 0x01) - (nullStatus >>6 & 0x01));
	//			outlierSectionLocation = PARAMETER_SECTION_BYTES + numMarkers * bytesPerSampMark;
	//			numBytesOfOutOfRangeValues = Compression.bytesToInt(new byte[]{parameterSection[outlierSectionLocation], parameterSection[outlierSectionLocation+1], parameterSection[outlierSectionLocation+2], parameterSection[outlierSectionLocation+3]});
	//			readBuffer = new byte[PARAMETER_SECTION_BYTES + (endMarkerIndex - beginMarkerIndex + 1) * bytesPerSampMark + numBytesOfOutOfRangeValues];
	//			file.seek(pos);
	//			file.read(readBuffer, PARAMETER_SECTION_BYTES+beginMarkerIndex*bytesPerSampMark, PARAMETER_SECTION_BYTES+endMarkerIndex*bytesPerSampMark);
	//			file.close();
	//		} catch (IOException e) {
	//			e.printStackTrace();
	//		}
	////		System.out.println("Finished loadToMarkerData from RAF in: "+(new Date().getTime() - time1)+" ms");
	//		return parameterSection;
	//	}
	
	public static void main(String[] args) {
		Project proj = new Project(Project.DEFAULT_CURRENT, false);
		String[] samples = proj.getSamples();
		Sample fsamp;
		
		for (int i = 0; i<samples.length; i++) {
			fsamp = proj.getFullSampleFromRandomAccessFile(samples[i]);
			fsamp.compareCalculationsFile(proj, proj.getMarkerNames(), proj.getProjectDir()+samples[i]+"_comp.xln");
        }

		//tests
//		testLoadTime(new Project("C:/workspace/Genvisis/projects/GEDI_exome.properties", false));
//		System.out.println(Compression.bytesToLong(Compression.longToBytes(15351532491l), 0));
    }

}
