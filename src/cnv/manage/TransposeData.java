// -Xms1024M -Xmx1024M     or even better: -Xmx15g
package cnv.manage;

import java.io.*;
import java.util.*;

import cnv.filesys.Compression;
import cnv.filesys.MarkerData;
import cnv.filesys.Sample;
import cnv.filesys.MarkerLookup;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import common.*;

// TODO whole file needs clean up
public class TransposeData {
	public static final int DEFAULT_MARKERS_PER_FILE = 210000;
	public static final byte MARKDATA_PARAMETER_TOTAL_LEN = 21;
	public static final byte MARKDATA_NUMSAMPS_START = 0;
	public static final byte MARKDATA_NUMSAMPS_LEN = 4;
	public static final byte MARKDATA_NUMMARKS_START = 4;
	public static final byte MARKDATA_NUMMARKS_LEN = 4;
	public static final byte MARKDATA_NULLSTATUS_START = 8;
	public static final byte MARKDATA_NULLSTATUS_LEN = 1;
	public static final byte MARKDATA_FINGERPRINT_START = 9;
	public static final byte MARKDATA_FINGERPRINT_LEN = 8;
	public static final byte MARKDATA_MARKNAMELEN_START = 17;
	public static final byte MARKDATA_MARKNAMELEN_LEN = 4;
	public static final byte MARKDATA_MARKNAME_START = 21;

	@SuppressWarnings("unchecked")
	public static void transposeData(Project proj, long maxMarkerFileSize, boolean keepAllSampleFilesOpen) {
		String[] allSamplesInProj;
		String[] allMarkersInProj;
		int numBytesPerMarker;
		int numMarkersInWriterBuffer;
		int numMarkersLastReadRoundLastBuffer = 0;
        int numMarkersPerBufferChunk;
		int maxNumMarkersPerFile;
//		int[] numMarkersInEachFile;
		int numChunksInWriteBuffer;
		int numBufferChunksPerFile;
		int numMarkerFiles;
		int numMarkersLastReadRound;
		int timesEachSampleFileBeRead;
		String[] markerFilenames;
		int markerFileIndex;
		int markerFileBufferChunkIndex;
		int firstMarkerInCurrentFile;
		int markerIndex;
//		int writeBufferSeekStep;
		RandomAccessFile[] sampleFiles = null;
		RandomAccessFile sampleFile;
        RandomAccessFile markerFile;
        byte[][] markerDataWriteBuffer = null;
        byte[] markerDataWriteBufferParameter;
        byte[] sampleDataReadBuffer;
        int indexInSampleDataReadBuffer;
        int locationInWriteBufferChunk;
		Hashtable<String,String> markerLookupHash = new Hashtable<String,String>();
		byte[][] markersInEachFile;
		String[] markersInEachFile1;
		byte[] markerNamesByteStream;
//		byte[] sampFileMarkerNameReadBuffer;
		int numSampFileOutliers;
		byte[] sampFileOutliersReadBuffer = null;
//		float[] outOfRangeValuesInSampleFile;
		Hashtable<String, Float> outOfRangeValuesInSampleFile = null;
//		int indexOutOfRangeValuesInSampleFile;
		Hashtable<String, Float>[] markFileWriteBufferOutliers;
		Hashtable<String, Float> markFileWriteBufferOutliersAdj;
		byte[] markFileWriteBufferOutliersBytes;
//		int indexMarkFileOutliers;
		Hashtable<String, Float> allOutliers;
		byte nullStatus = 0;
        byte bytesPerSampMark;
        int numMarkersCurrentLoop;
		long timer1, timer2, timer3, timer4, timer5, timer6;

		allSamplesInProj = proj.getSamples();
//		String[] samplesTemp = proj.getSamples();
//		allSamplesInProj = new String[10];
//		for (int i=0; i<10; i++) {
//			allSamplesInProj[i] = samplesTemp[i];
//		}
		allMarkersInProj = proj.getMarkerNames();
		try {
			sampleFile = new RandomAccessFile(proj.getDir(Project.SAMPLE_DIRECTORY, true) + allSamplesInProj[0] + Sample.SAMPLE_DATA_FILE_EXTENSION, "r");
			sampleFile.readInt();
			nullStatus = sampleFile.readByte();
		} catch (FileNotFoundException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
			return;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		bytesPerSampMark = (byte) (Compression.BYTES_PER_SAMPLE_MARKER - (nullStatus & 0x01) - (nullStatus >>1 & 0x01) - (nullStatus >>2 & 0x01) - (nullStatus >>3 & 0x01) - (nullStatus >>4 & 0x01) - (nullStatus >>5 & 0x01) - (nullStatus >>6 & 0x01));
		numBytesPerMarker = allSamplesInProj.length * bytesPerSampMark;
//		writeBufferSeekStep = numBytesPerMarker - BYTES_PER_SAMPLE_MARKER;
		if (new File(proj.getProjectDir()).getFreeSpace() <= (allSamplesInProj.length * (long)allMarkersInProj.length * bytesPerSampMark)) {
			System.err.println("Not enough space (available: "+ext.prettyUpSize(new File(proj.getProjectDir()).getFreeSpace(), 1)+") for all the new data to be created (required: "+ext.prettyUpSize(new File(proj.getProjectDir()).getFreeSpace(), 1)+").");
			return;
		}


		if (maxMarkerFileSize == 0) {
			maxMarkerFileSize = Integer.MAX_VALUE;
		}
		numMarkersInWriterBuffer = Math.min((int)((0.65*(double)Runtime.getRuntime().maxMemory())/numBytesPerMarker), allMarkersInProj.length);
//		numMarkersPerBufferChunk = Math.min((int) Integer.MAX_VALUE/numBytesPerMarker, numMarkersInWriterBuffer);
		if (maxMarkerFileSize > Integer.MAX_VALUE) {
			numMarkersPerBufferChunk = (int) Math.min((double) Integer.MAX_VALUE/numBytesPerMarker, (double) numMarkersInWriterBuffer);
		} else {
			numMarkersPerBufferChunk = (int) Math.min((double) maxMarkerFileSize/numBytesPerMarker, (double) numMarkersInWriterBuffer);
		}
		numChunksInWriteBuffer = numMarkersInWriterBuffer/numMarkersPerBufferChunk;
		numMarkersInWriterBuffer = numChunksInWriteBuffer * numMarkersPerBufferChunk;
		timesEachSampleFileBeRead = (int) Math.ceil((double)allMarkersInProj.length / (double)numMarkersInWriterBuffer);
		numMarkersLastReadRound = allMarkersInProj.length % numMarkersInWriterBuffer;

		maxNumMarkersPerFile = (int) Math.min((double)maxMarkerFileSize/(double)numBytesPerMarker, allMarkersInProj.length);
		//TODO remove here
//		numBufferChunksPerFile = maxNumMarkersPerFile / numMarkersPerBufferChunk;
		numBufferChunksPerFile = (int) Math.ceil((double)maxNumMarkersPerFile / (double)numMarkersPerBufferChunk);
		maxNumMarkersPerFile = numMarkersPerBufferChunk * numBufferChunksPerFile;
		numMarkerFiles = (int) Math.ceil((double)allMarkersInProj.length / (double)maxNumMarkersPerFile);
		markerFilenames = new String[numMarkerFiles];
		markersInEachFile = new byte[numMarkerFiles][];
//		numMarkersInEachFile = new int[numMarkerFiles];
		backupOlderRafs(proj.getDir(Project.MARKER_DATA_DIRECTORY, true, new Logger(), false), MarkerData.MARKER_DATA_FILE_EXTENSION);
		markerIndex = 0;
		for (int i=0; i<numMarkerFiles; i++) {
			markerFilenames[i]=proj.getDir(Project.MARKER_DATA_DIRECTORY)+"markers." + i + MarkerData.MARKER_DATA_FILE_EXTENSION;
			numMarkersCurrentLoop = Math.min(maxNumMarkersPerFile, allMarkersInProj.length - markerIndex);
			markersInEachFile1 = new String[numMarkersCurrentLoop];
//			numMarkersInEachFile[i] =; 
			for (int j=0; j<numMarkersCurrentLoop; j++) {
				markerLookupHash.put(allMarkersInProj[markerIndex], "markers." + i + MarkerData.MARKER_DATA_FILE_EXTENSION + "\t" + j);
				markersInEachFile1[j] = allMarkersInProj[markerIndex];
				markerIndex ++;
			}
			try {
				markersInEachFile[i] = Compression.objToBytes(markersInEachFile1);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		new MarkerLookup(markerLookupHash).serialize(proj.getFilename(Project.MARKERLOOKUP_FILENAME, false, false));
		System.out.println("Total markers in the project:\t"+allMarkersInProj.length
						  +"\nTotal samples in the project:\t"+allSamplesInProj.length
						  +"\nMaximum memory size available:\t"+Runtime.getRuntime().maxMemory()/1024/1024/1024+"."+(Runtime.getRuntime().maxMemory()/1024/1024/10 - Runtime.getRuntime().maxMemory()/1024/1024/1024*102)+" gb"
						  +"\nNumber of markers in write buffer at a time:\t"+numMarkersInWriterBuffer
						  +"\nNumber of markers per write buffer Chunk:\t"+numMarkersPerBufferChunk
						  +"\nMaximum file size of marker files:\t"+maxMarkerFileSize/1024/1024/1024+"."+((int)(maxMarkerFileSize/1024/1024/10.24) - (int)(maxMarkerFileSize/1024/1024/1024*102.4))+" gb"
						  +"\nNumber of markers per marker file:\t"+maxNumMarkersPerFile
						  +"\nNumber of write buffer chunks per marker file:\t"+numBufferChunksPerFile
						  +"\nNumber of marker files:\t"+numMarkerFiles);

		timer1 = new Date().getTime();
		timer3 = 0;
		firstMarkerInCurrentFile = 0;
		sampleDataReadBuffer = new byte[numMarkersPerBufferChunk * bytesPerSampMark];
		markerDataWriteBuffer = new byte[numChunksInWriteBuffer][numMarkersPerBufferChunk * numBytesPerMarker];
//		markFileOutliersWriteBufferVec = new Vector[numChunksInWriteBuffer];
		markFileWriteBufferOutliers = new Hashtable[numChunksInWriteBuffer];
		allOutliers = new Hashtable<String, Float>();
		try {
			if (keepAllSampleFilesOpen) {
				sampleFiles = new RandomAccessFile[allSamplesInProj.length];
				for (int i=0; i<allSamplesInProj.length; i++) {
					sampleFiles[i] = new RandomAccessFile(proj.getDir(Project.SAMPLE_DIRECTORY, true) + allSamplesInProj[i] + Sample.SAMPLE_DATA_FILE_EXTENSION, "r");
				}
			}

			markerFileIndex = 0;

			markerNamesByteStream = Compression.objToBytes(markersInEachFile[markerFileIndex]);
			markerDataWriteBufferParameter = new byte[MARKDATA_PARAMETER_TOTAL_LEN + markerNamesByteStream.length];
			System.arraycopy(Compression.intToBytes(allSamplesInProj.length), 0, markerDataWriteBufferParameter, MARKDATA_NUMSAMPS_START, MARKDATA_NUMSAMPS_LEN);
			System.arraycopy(Compression.intToBytes(maxNumMarkersPerFile), 0, markerDataWriteBufferParameter, MARKDATA_NUMMARKS_START, MARKDATA_NUMMARKS_LEN);
			markerDataWriteBufferParameter[MARKDATA_NULLSTATUS_START] = nullStatus;
			System.arraycopy(Compression.longToBytes(MarkerSet.fingerprint(allSamplesInProj)), 0, markerDataWriteBufferParameter, MARKDATA_FINGERPRINT_START, MARKDATA_FINGERPRINT_LEN);
			System.arraycopy(Compression.intToBytes(markerNamesByteStream.length), 0, markerDataWriteBufferParameter, MARKDATA_MARKNAMELEN_START, MARKDATA_MARKNAMELEN_LEN);
			System.arraycopy(markerNamesByteStream, 0, markerDataWriteBufferParameter, MARKDATA_MARKNAME_START, markerNamesByteStream.length);
			markerFile = new RandomAccessFile(markerFilenames[markerFileIndex], "rw");
			markerFile.write(markerDataWriteBufferParameter);

//			markerFile = new RandomAccessFile(markerFilenames[markerFileIndex], "rw");
////			sampFileMarkerNameReadBuffer = markersInEachFile[markerFileIndex].getBytes();
//			markerFile.writeInt(allSamplesInProj.length);
//			markerFile.writeInt(maxNumMarkersPerFile);
//			markerFile.writeByte(nullStatus);
//			markerFile.writeLong(MarkerSet.fingerprint(allSamplesInProj));
////			markerFile.writeInt(sampFileMarkerNameReadBuffer.length);
////			markerFile.write(sampFileMarkerNameReadBuffer);
////			markerFile.writeInt(markersInEachFile[markerFileIndex].length);
////			markerFile.write(markersInEachFile[markerFileIndex]);
//			markerFile.writeInt(markerNamesByteStream.length);
//			markerFile.write(markerNamesByteStream);

//			markFileOutliersWriteBufferVecAdj = new Vector<Byte>();
			markFileWriteBufferOutliersAdj = new Hashtable<String, Float>();
			markerFileBufferChunkIndex = 0;
			for(int i=0; i<timesEachSampleFileBeRead; i++) {
				timer2 = new Date().getTime();
				timer5 = 0;
				timer6 = 0;
				if ((i+1)==timesEachSampleFileBeRead && numMarkersLastReadRound!=0) {
					numChunksInWriteBuffer = (int) Math.ceil((double)numMarkersLastReadRound / (double)numMarkersPerBufferChunk);
//					markFileOutliersWriteBufferVec = new Vector[numChunksInWriteBuffer];
					markFileWriteBufferOutliers = new Hashtable[numChunksInWriteBuffer];
					if ( numChunksInWriteBuffer>1 ) {
						markerDataWriteBuffer = null;
						markerDataWriteBuffer = new byte[numChunksInWriteBuffer][numMarkersPerBufferChunk * numBytesPerMarker];
						numMarkersLastReadRoundLastBuffer = numMarkersLastReadRound % numMarkersPerBufferChunk;
						if ( numMarkersLastReadRoundLastBuffer!=0 ) {
							markerDataWriteBuffer[numChunksInWriteBuffer-1] = new byte[numMarkersLastReadRoundLastBuffer * numBytesPerMarker];
						}
					} else {
						markerDataWriteBuffer = null;
						markerDataWriteBuffer = new byte[numChunksInWriteBuffer][numMarkersLastReadRound * numBytesPerMarker];
					}
				}
//				for (int j=0; j<markFileOutliersWriteBufferVec.length; j++) {
//					markFileOutliersWriteBufferVec[j] = new Vector<Byte>();
//				}
				for (int j=0; j<markFileWriteBufferOutliers.length; j++) {
					markFileWriteBufferOutliers[j] = new Hashtable<String, Float>();
				}
				for (int j=0; j<allSamplesInProj.length; j++) {
					if (!keepAllSampleFilesOpen) {
						sampleFile = new RandomAccessFile(proj.getDir(Project.SAMPLE_DIRECTORY, true) + allSamplesInProj[j] + Sample.SAMPLE_DATA_FILE_EXTENSION, "r"); //TODO
					} else {
						sampleFile = sampleFiles[j];
					}
					markerIndex = firstMarkerInCurrentFile;
				    for (int k=0; k < numChunksInWriteBuffer; k++) {
				    	timer4 = new Date().getTime();
//				    	sampleFiles[j].read(sampleDataReadBuffer);
				    	sampleFile.seek(Sample.PARAMETER_SECTION_BYTES + markerIndex * bytesPerSampMark);
				    	sampleFile.read(sampleDataReadBuffer);
				    	// Read in the over range value array
				    	sampleFile.seek(Sample.PARAMETER_SECTION_BYTES + allMarkersInProj.length * bytesPerSampMark);
				    	numSampFileOutliers = sampleFile.readInt();
				    	if (numSampFileOutliers > 0) {
				    		sampFileOutliersReadBuffer = new byte[numSampFileOutliers];
//				    		outOfRangeValuesInSampleFile = new float[numOutOfRangeValuesInSampleFile];
				    		sampleFile.read(sampFileOutliersReadBuffer);
//				    		outOfRangeValuesInSampleFile = new Hashtable<String, Float>();
				    		outOfRangeValuesInSampleFile = (Hashtable<String, Float>) Compression.bytesToObj(sampFileOutliersReadBuffer);
				    	}
				    	indexInSampleDataReadBuffer = 0;
				    	timer5 += (new Date().getTime() - timer4);
				    	timer4 = new Date().getTime();
						locationInWriteBufferChunk = j * bytesPerSampMark;
				    	for (int l=0; l<numMarkersPerBufferChunk && markerIndex<allMarkersInProj.length; l++) {
				    		for (int m=0; m < bytesPerSampMark; m++) {
			    				markerDataWriteBuffer[k][locationInWriteBufferChunk + m] = sampleDataReadBuffer[indexInSampleDataReadBuffer];
				    			if (m==3 && (sampleDataReadBuffer[indexInSampleDataReadBuffer-1] == Compression.REDUCED_PRECISION_XY_OUT_OF_RANGE_BYTES[0]
				    					     && sampleDataReadBuffer[indexInSampleDataReadBuffer] == Compression.REDUCED_PRECISION_XY_OUT_OF_RANGE_BYTES[1])) {
				    				allOutliers.put(markerIndex + "\t" + j +"\tx", outOfRangeValuesInSampleFile.get(markerIndex + "\tx"));
				    				markFileWriteBufferOutliers[k].put(markerIndex + "\t" + j + "\tx", outOfRangeValuesInSampleFile.get(markerIndex + "\tx"));
				    			} else if (m==5 && (sampleDataReadBuffer[indexInSampleDataReadBuffer-1] == Compression.REDUCED_PRECISION_XY_OUT_OF_RANGE_BYTES[0]
			    					                && sampleDataReadBuffer[indexInSampleDataReadBuffer] == Compression.REDUCED_PRECISION_XY_OUT_OF_RANGE_BYTES[1])) {
				    				allOutliers.put(markerIndex + "\t" + j + "\ty", outOfRangeValuesInSampleFile.get(markerIndex + "\ty"));
				    				markFileWriteBufferOutliers[k].put(markerIndex + "\t" + j + "\ty", outOfRangeValuesInSampleFile.get(markerIndex + "\ty"));
				    			} else if (m==10 && (sampleDataReadBuffer[indexInSampleDataReadBuffer-2] == Compression.REDUCED_PRECISION_LRR_OUT_OF_RANGE_BYTES[0]
				    								 && sampleDataReadBuffer[indexInSampleDataReadBuffer-1] == Compression.REDUCED_PRECISION_LRR_OUT_OF_RANGE_BYTES[1]
				    								 && sampleDataReadBuffer[indexInSampleDataReadBuffer] == Compression.REDUCED_PRECISION_LRR_OUT_OF_RANGE_BYTES[2])) {
				    				allOutliers.put(markerIndex + "\t" + j + "\tlrr", outOfRangeValuesInSampleFile.get(markerIndex + "\tlrr"));
				    				markFileWriteBufferOutliers[k].put(markerIndex + "\t" + j + "\tlrr", outOfRangeValuesInSampleFile.get(markerIndex + "\tlrr"));
				    			}
//				    			if ((m==2 && (0x80 & sampleDataReadBuffer[indexInSampleDataReadBuffer])>0) || (m==4 && (0x80 & sampleDataReadBuffer[indexInSampleDataReadBuffer])>0)) {
//				    				allOutliers.put(allSamplesInProj[j]+"\t"+allMarkersInProj[markerIndex]+"\t"+(m==2?"x":"y"), outOfRangeValuesInSampleFile.get(allSamplesInProj[j]+"\t"+allMarkersInProj[markerIndex]+"\t"+(m==2?"x":"y")));
//				    				markFileWriteBufferOutliers[k].put(allSamplesInProj[j]+"\t"+allMarkersInProj[markerIndex]+"\t"+(m==2?"x":"y"), outOfRangeValuesInSampleFile.get(allSamplesInProj[j]+"\t"+allMarkersInProj[markerIndex]+"\t"+(m==2?"x":"y")));
//				    			}
			    				indexInSampleDataReadBuffer ++;
					    	}
				    		locationInWriteBufferChunk += numBytesPerMarker;
							markerIndex ++;
				    	}
				    	timer6 += (new Date().getTime() - timer4);
					}
				    if (!keepAllSampleFilesOpen) {
				    	sampleFile.close();
				    }
				}
				timer2 = new Date().getTime() - timer2;
				System.out.println("Load in sample files - round " + i + ": " + timer2/1000 + " secs = "+timer5/1000+" (file.read) + "+timer6/1000+" (array shuffle)");
				for (int j=0; j<markerDataWriteBuffer.length; j++) {
			    	timer4 = new Date().getTime();
					markerFile.write(markerDataWriteBuffer[j]);
					markFileWriteBufferOutliersAdj.putAll(markFileWriteBufferOutliers[j]);
			    	timer3 += (new Date().getTime() - timer4);
					markerFileBufferChunkIndex ++;
					if (markerFileBufferChunkIndex >= numBufferChunksPerFile) {
						if (markFileWriteBufferOutliersAdj == null || markFileWriteBufferOutliersAdj.size() == 0) {
							markerFile.writeInt(0);
						} else {
							markFileWriteBufferOutliersBytes = Compression.objToBytes(markFileWriteBufferOutliersAdj);
							markerFile.writeInt(markFileWriteBufferOutliersBytes.length);
							markerFile.write(markFileWriteBufferOutliersBytes);
						}
						markerFile.close();
						System.out.println("Write marker file " + markerFileIndex + ":\t"+ timer3/1000 + " sec."); // \toutlier #: " + (markFileOUtliersWriteBufferArray.length/4));
						markerFileBufferChunkIndex = 0;
						markerFileIndex ++;
						markFileWriteBufferOutliersAdj = new Hashtable<String, Float>();
						if ((i+1) != timesEachSampleFileBeRead || (j+1) != markerDataWriteBuffer.length) {
//							markerFile = new RandomAccessFile(markerFilenames[markerFileIndex], "rw");
////							sampFileMarkerNameReadBuffer = markersInEachFile[markerFileIndex].getBytes();
//							markerFile.writeInt(allSamplesInProj.length);
//							if ((markerFileIndex+1) == numMarkerFiles) {
//								markerFile.writeInt(allMarkersInProj.length % maxNumMarkersPerFile);
//							} else {
//								markerFile.writeInt(maxNumMarkersPerFile);
//							}
//							markerFile.writeByte(nullStatus);
//							markerFile.writeLong(MarkerSet.fingerprint(allSamplesInProj));
////							markerFile.writeInt(sampFileMarkerNameReadBuffer.length);
////							markerFile.write(sampFileMarkerNameReadBuffer);
//							markerFile.writeInt(markersInEachFile[markerFileIndex].length);
//							markerFile.write(markersInEachFile[markerFileIndex]);

							markerNamesByteStream = Compression.objToBytes(markersInEachFile[markerFileIndex]);
							markerDataWriteBufferParameter = new byte[MARKDATA_PARAMETER_TOTAL_LEN + markerNamesByteStream.length];
							System.arraycopy(Compression.intToBytes(allSamplesInProj.length), 0, markerDataWriteBufferParameter, MARKDATA_NUMSAMPS_START, MARKDATA_NUMSAMPS_LEN);
							System.arraycopy(Compression.intToBytes((markerFileIndex+1)==numMarkerFiles? (allMarkersInProj.length % maxNumMarkersPerFile) : maxNumMarkersPerFile), 0, markerDataWriteBufferParameter, MARKDATA_NUMMARKS_START, MARKDATA_NUMMARKS_LEN);
							markerDataWriteBufferParameter[MARKDATA_NULLSTATUS_START] = nullStatus;
							System.arraycopy(Compression.longToBytes(MarkerSet.fingerprint(allSamplesInProj)), 0, markerDataWriteBufferParameter, MARKDATA_FINGERPRINT_START, MARKDATA_FINGERPRINT_LEN);
							System.arraycopy(Compression.intToBytes(markerNamesByteStream.length), 0, markerDataWriteBufferParameter, MARKDATA_MARKNAMELEN_START, MARKDATA_MARKNAMELEN_LEN);
							System.arraycopy(markerNamesByteStream, 0, markerDataWriteBufferParameter, MARKDATA_MARKNAME_START, markerNamesByteStream.length);
							markerFile = new RandomAccessFile(markerFilenames[markerFileIndex], "rw");
							markerFile.write(markerDataWriteBufferParameter);

//							markFileOutliersWriteBufferVecAdj = new Vector<Byte>();
						}
						timer3 = 0;
					}
				}
				firstMarkerInCurrentFile += numMarkersInWriterBuffer;
			}
			if (markerFileBufferChunkIndex != 0) {
				if (markFileWriteBufferOutliersAdj == null || markFileWriteBufferOutliersAdj.size() == 0) {
					markerFile.writeInt(0);
				} else {
					markFileWriteBufferOutliersBytes = Compression.objToBytes(markFileWriteBufferOutliersAdj);
					markerFile.writeInt(markFileWriteBufferOutliersBytes.length);
					markerFile.write(markFileWriteBufferOutliersBytes);
				}
				markerFile.close();
			}
			if (allOutliers != null && allOutliers.size() != 0) {
				Files.writeSerial(allOutliers, proj.getDir(Project.MARKER_DATA_DIRECTORY) + "outliers" + MarkerData.MARKER_DATA_FILE_EXTENSION);
			}
			System.out.println("Write marker file " + markerFileIndex + ":\t" + timer3/1000 + " sec."); // \toutlier #: " + (markFileOUtliersWriteBufferArray.length/4));
			
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		timer1 = (new Date().getTime() - timer1);
		System.out.println("******\nTotal Time: "+ (timer1/3600000) +" hrs "+ (timer1/60000 - 60*(timer1/3600000)) +" mins "+ (timer1/1000 - 60*(timer1/60000)) +" secs.");
	}

	/*
	 * Delete existing files of a specific file name extension.
	 * @param directory	the folder where the existing files are
	 * @param fileNameExtension	the file name extension, including the ".". For example, ".scat", ".xls", ".csv", and etc. 
	 */
	public static void deleteOlderRafs(String directory, final String fileNameExtesion) {
		File[] files;

		// List of files to be backed up.
		files = new File(directory).listFiles(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				return filename.endsWith(fileNameExtesion);
			}
		});

		// Move files to backup folder
		for (int i = 0; i<files.length; i++) {
			files[i].delete();
		}
		
		if (files.length >= 0) {
			System.out.println("Older version of the data in " + directory + " has been deleted from the hard drive.");
		}
	}

	/*
	 * Move existing files of a specific file name extension to a sub-folder named "Backup.?"
	 * @param directory	the folder where the existing files are
	 * @param fileNameExtension	the file name extension, including the ".". For example, ".scat", ".xls", ".csv", and etc. 
	 */
	public static void backupOlderRafs(String directory, final String fileNameExtesion) {
		File[] files;
		int count;

		// Check existing backup folders.
		count = 0;
		do {
			count++;
		} while (new File(directory+"Backup."+count).exists());

		// List of files to be backed up.
		files = new File(directory).listFiles(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				return filename.endsWith(fileNameExtesion);
			}
		});

		if (files.length>0) {
			// Create a new backup folder.
			new File(directory+"Backup."+count+"/").mkdirs();
			
			// Move files to backup folder
			for (int i = 0; i<files.length; i++) {
				files[i].renameTo(new File(directory + "Backup." + count + "/" + files[i].getName()));
			}

			System.out.println("Older version of the data in " + directory + " has been moved to " + directory + "Backup." + count + "/ To save disk space, please manually delete these files.");
		}
	}

	/**
	 * Create MarkerLookup from MarkerData (RAF and half-precision format)
	 * @param proj
	 */
	// TODO refactor without the number 
	public static void createMarkerLookup2(Project proj) {
		Hashtable<String,String> hash = new Hashtable<String,String>();
		String[] files, markerNames;
		byte[] readBuffer;
		RandomAccessFile currentFile;
		long time;

		time = new Date().getTime();
		System.out.println("Creating MarkerLookup file");
		files = new File(proj.getDir(Project.MARKER_DATA_DIRECTORY)).list(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				return filename.endsWith(MarkerData.MARKER_DATA_FILE_EXTENSION);
			}
		});
		if (files==null) {
			System.err.println("Error - failed to create MarkerLookup -- directory does not exist: "+proj.getDir(Project.MARKER_DATA_DIRECTORY));
		} else if (files.length==0) {
			System.err.println("Error - failed to create MarkerLookup -- no "+MarkerData.MARKER_DATA_FILE_EXTENSION+" files available");
		} else {
			for (int i = 0; i<files.length; i++) {
				try {
					currentFile = new RandomAccessFile(proj.getDir(Project.MARKER_DATA_DIRECTORY)+files[i], "r");
					readBuffer = new byte[currentFile.readInt()];
					currentFile.read(readBuffer);
					markerNames = new String(readBuffer).split("\t");
					for (int j = 0; j<markerNames.length; j++) {
						hash.put(markerNames[j], files[i]+"\t"+j);
					}
				} catch (FileNotFoundException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}

			new MarkerLookup(hash).serialize(proj.getFilename(Project.MARKERLOOKUP_FILENAME));

			System.out.println("Created MarkerLookup in "+ext.getTimeElapsed(time));
		}
	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		Project proj;
		String filename = Project.DEFAULT_PROJECT;
		// String demo = "flagged_results.txt";
		boolean transpose = false;
		int maxFileSize = 0;
		boolean lookup = false;

		String usage = "\n"+
		"filesys.ExtractPlots requires 0-1 arguments\n"+
		"   (1) project file (i.e. proj="+filename+" (default))\n"+
		"   (2) transpose data (i.e. -transpose ("+(transpose?"":"not the ")+"default))\n"+
		"   (3) maximum size of each file in bytes (i.e. max="+maxFileSize+" (default))\n"+
		"  OR:\n"+
		"   (7) create marker lookup table (i.e. -lookup ("+(lookup?"":"not the ")+"default))\n"+
		"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("-transpose")) {
				transpose = true;
				numArgs--;
			} else if (args[i].startsWith("max=")) {
				maxFileSize = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("-lookup")) {
				lookup = true;
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}

		proj = new Project(filename, false);

		if (!proj.getDir(Project.SOURCE_DIRECTORY).equals("")&&!new File(proj.getDir(Project.SOURCE_DIRECTORY)).exists()) {
			System.err.println("Error - the project source location is invalid: "+proj.getDir(Project.SOURCE_DIRECTORY));
			return;
		}

//		transpose = true;
//		lookup = true;
		try {
			if (transpose) {
				transposeData(proj, maxFileSize, true);
			} else if (lookup) {
				createMarkerLookup2(proj);
			} else {
				System.err.println("Error - no selection was made");
				System.out.println(usage);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
