// -Xms1024M -Xmx1024M     or even better: -Xmx15g
package cnv.manage;

import java.io.*;
import java.text.SimpleDateFormat;
import java.util.*;

import cnv.LaunchProperties;
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

	/**
	 * Transpose data from sample based structure into Marker based structure.
	 * The input is a set of files based on the sample oriented structure, and
	 * the output is a set of data files based on the marker oriented structure.
	 * @param proj The Genivis project containing the data to be transposed
	 * @param maxMarkerFileSize The maximum size in mb of a single output
	 * file. If set to 0, the default value 2000 is taken on.
	 * @param keepAllSampleFilesOpen To keep all the files or just a single one
	 * open at a time.
	 * @param log The log file. If set as null, the default file
	 * Genvisis_TransposeData_yyyyMMdd_HHmmss.log will be used.
	 */
	@SuppressWarnings("unchecked")
	public static void transposeData(Project proj, long maxMarkerFileSize, boolean keepAllSampleFilesOpen, Logger log) {
		String[] allSamplesInProj;
		String[] allMarkersInProj;
		int numBytesPerMarker;
		int numMarkersWrtBuffer;
		int numMarkersLastReadRoundLastBuffer = 0;
        int numMarksPerBufferChunk;
		int maxNumMarkersPerFile;
//		int[] numMarkersInEachFile;
		int numChunksWriteBuffer;
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
//		byte[] markerNamesByteStream;
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
		byte nullStatus = Byte.MIN_VALUE;
        byte bytesPerSampMark;
        int numMarkersCurrentLoop;
        byte backupCount;
		long timer1, timer2, timer3, timer4, timer5, timer6;
		boolean done, safe;
		long memoryReserve;
		int oomeLoops;

		if (log == null) {
			log = new Logger(proj.getDir(Project.MARKER_DATA_DIRECTORY) + "Genvisis_" + (new SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date())) + ".log");
		}
		log.report("Transposing data of the project " + proj.getProjectDir());

		allSamplesInProj = proj.getSamples();
		allMarkersInProj = proj.getMarkerNames();
		nullStatus = Sample.getNullstatusFromRandomAccessFile(proj.getDir(Project.SAMPLE_DIRECTORY, true) + allSamplesInProj[0] + Sample.SAMPLE_DATA_FILE_EXTENSION, false);
		bytesPerSampMark = (byte) (Compression.BYTES_PER_SAMPLE_MARKER - (nullStatus & 0x01) - (nullStatus >>1 & 0x01) - (nullStatus >>2 & 0x01) - (nullStatus >>3 & 0x01) - (nullStatus >>4 & 0x01) - (nullStatus >>5 & 0x01) - (nullStatus >>6 & 0x01));
		numBytesPerMarker = allSamplesInProj.length * bytesPerSampMark;
		if (new File(proj.getProjectDir()).getFreeSpace() <= (allSamplesInProj.length * (long)allMarkersInProj.length * bytesPerSampMark)) {
//			System.err.println("Not enough space (available: " + ext.prettyUpSize(new File(proj.getProjectDir()).getFreeSpace(), 1)+") for all the new data to be created (required: "+ext.prettyUpSize(new File(proj.getProjectDir()).getFreeSpace(), 1)+").");
			log.reportError("Not enough disk space for all the new data to be created. Available: " + ext.prettyUpSize(new File(proj.getProjectDir()).getFreeSpace(), 1) + "; Required: " + ext.prettyUpSize(new File(proj.getProjectDir()).getFreeSpace(), 1) + ").");
			return;
		}
		if (maxMarkerFileSize == 0) {
			maxMarkerFileSize = Integer.MAX_VALUE;
		}

		done = false;
		safe = false;
		oomeLoops = 0;
		timer1 = new Date().getTime();
		memoryReserve = (long)((double)Runtime.getRuntime().maxMemory() * 0.25);
		while (!done) {
			try {
		//		numMarkersInWriterBuffer = Math.min((int)((0.65*(double)Runtime.getRuntime().maxMemory())/numBytesPerMarker), allMarkersInProj.length);
				numMarkersWrtBuffer = Math.min((int)(((double)Runtime.getRuntime().maxMemory() - memoryReserve) / numBytesPerMarker), allMarkersInProj.length);
		//		numMarkersPerBufferChunk = Math.min((int) Integer.MAX_VALUE/numBytesPerMarker, numMarkersInWriterBuffer);
				if (maxMarkerFileSize > Integer.MAX_VALUE) {
					numMarksPerBufferChunk = (int) Math.min((double) Integer.MAX_VALUE/numBytesPerMarker, (double) numMarkersWrtBuffer);
				} else {
					numMarksPerBufferChunk = (int) Math.min((double) maxMarkerFileSize/numBytesPerMarker, (double) numMarkersWrtBuffer);
				}
				numChunksWriteBuffer = numMarkersWrtBuffer/numMarksPerBufferChunk;
				numMarkersWrtBuffer = numChunksWriteBuffer * numMarksPerBufferChunk;
				timesEachSampleFileBeRead = (int) Math.ceil((double)allMarkersInProj.length / (double)numMarkersWrtBuffer);
				numMarkersLastReadRound = allMarkersInProj.length % numMarkersWrtBuffer;
		
				maxNumMarkersPerFile = (int) Math.min((double)maxMarkerFileSize/(double)numBytesPerMarker, allMarkersInProj.length);
				//TODO remove here
		//		numBufferChunksPerFile = maxNumMarkersPerFile / numMarkersPerBufferChunk;
				numBufferChunksPerFile = (int) Math.ceil((double)maxNumMarkersPerFile / (double)numMarksPerBufferChunk);
				maxNumMarkersPerFile = numMarksPerBufferChunk * numBufferChunksPerFile;
				numMarkerFiles = (int) Math.ceil((double)allMarkersInProj.length / (double)maxNumMarkersPerFile);
				markerFilenames = new String[numMarkerFiles];
				markersInEachFile = new byte[numMarkerFiles][];
		//		numMarkersInEachFile = new int[numMarkerFiles];
				backupCount = backupOlderRafs(proj.getDir(Project.MARKER_DATA_DIRECTORY, true, new Logger(), false), new String[] {MarkerData.MARKER_DATA_FILE_EXTENSION, "outliers.ser"}, true);
				backupOlderRafs(proj.getFilename(Project.MARKERLOOKUP_FILENAME, false, false), backupCount);
				markerIndex = 0;
				for (int i=0; i<numMarkerFiles; i++) {
					markerFilenames[i] = proj.getDir(Project.MARKER_DATA_DIRECTORY) + "markers." + i + MarkerData.MARKER_DATA_FILE_EXTENSION;
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
		
				timer3 = 0;
				firstMarkerInCurrentFile = 0;
		//		sampleDataReadBuffer = new byte[numMarkersPerBufferChunk * bytesPerSampMark];
				sampleDataReadBuffer = new byte[numMarkersWrtBuffer * bytesPerSampMark];
				
//				log.report("Before WriteBuffer\t"+ext.reportMemoryUsage());
//				log.report("numChunksInWriteBuffer="+numChunksWriteBuffer+"\tnumMarkersPerBufferChunk="+numMarksPerBufferChunk+"\tnumBytesPerMarker="+numBytesPerMarker);
//				log.report(numChunksWriteBuffer+"*"+numMarksPerBufferChunk+"*"+numBytesPerMarker+"="+(numChunksWriteBuffer*numMarksPerBufferChunk*numBytesPerMarker)+" ="+ext.prettyUpSize((long)numChunksWriteBuffer*(long)numMarksPerBufferChunk*(long)numBytesPerMarker, 2));
				markerDataWriteBuffer = new byte[numChunksWriteBuffer][numMarksPerBufferChunk * numBytesPerMarker];
		//		markFileOutliersWriteBufferVec = new Vector[numChunksInWriteBuffer];
				markFileWriteBufferOutliers = new Hashtable[numChunksWriteBuffer];
				allOutliers = new Hashtable<String, Float>();
//				log.report("After WriteBuffer\t"+ext.reportMemoryUsage());
				
				if (!safe) {
					safe = true;
					throw new OutOfMemoryError(); // one more round to be on the safe side
				}
				log.report("Memory optimization took "+oomeLoops+" round(s) ("+ext.prettyUpSize(memoryReserve, 1)+" was kept in reserve)");
				
				new MarkerLookup(markerLookupHash).serialize(proj.getFilename(Project.MARKERLOOKUP_FILENAME, false, false));
				log.report("Markers in the project =\t" + allMarkersInProj.length
								  + "\nSamples in the project =\t" + allSamplesInProj.length
								  + "\nMemory available =\t" + Runtime.getRuntime().maxMemory()/1024/1024/1024 + "." + (Runtime.getRuntime().maxMemory()/1024/1024/10 - Runtime.getRuntime().maxMemory()/1024/1024/1024*102) + " gb"
								  + "\nMarkers / write buffer =\t" + numMarkersWrtBuffer
								  + "\nMarkers / write buffer chunk =\t" + numMarksPerBufferChunk
								  + "\nWrite buffer chunks =\t" + numChunksWriteBuffer
								  + "\nMarker file size <=\t" + maxMarkerFileSize/1024/1024/1024 + "." + ((int)(maxMarkerFileSize/1024/1024/10.24) - (int)(maxMarkerFileSize/1024/1024/1024*102.4)) + " gb"
								  + "\nMarkers / marker file =\t" + maxNumMarkersPerFile
								  + "\nMarker files =\t" + numMarkerFiles
								  + "\nWrite buffer chunks / marker file =\t" + numBufferChunksPerFile);
				
				try {
					if (keepAllSampleFilesOpen) {
						sampleFiles = new RandomAccessFile[allSamplesInProj.length];
						for (int i=0; i<allSamplesInProj.length; i++) {
							sampleFiles[i] = new RandomAccessFile(proj.getDir(Project.SAMPLE_DIRECTORY, true) + allSamplesInProj[i] + Sample.SAMPLE_DATA_FILE_EXTENSION, "r");
						}
					}
		
					markerFileIndex = 0;
		
		//			markerNamesByteStream = Compression.objToBytes(markersInEachFile[markerFileIndex]);
		//			markerDataWriteBufferParameter = new byte[MARKDATA_PARAMETER_TOTAL_LEN + markerNamesByteStream.length];
					markerDataWriteBufferParameter = new byte[MARKDATA_PARAMETER_TOTAL_LEN + markersInEachFile[markerFileIndex].length];
					System.arraycopy(Compression.intToBytes(allSamplesInProj.length), 0, markerDataWriteBufferParameter, MARKDATA_NUMSAMPS_START, MARKDATA_NUMSAMPS_LEN);
					System.arraycopy(Compression.intToBytes(maxNumMarkersPerFile), 0, markerDataWriteBufferParameter, MARKDATA_NUMMARKS_START, MARKDATA_NUMMARKS_LEN);
					markerDataWriteBufferParameter[MARKDATA_NULLSTATUS_START] = nullStatus;
					System.arraycopy(Compression.longToBytes(MarkerSet.fingerprint(allSamplesInProj)), 0, markerDataWriteBufferParameter, MARKDATA_FINGERPRINT_START, MARKDATA_FINGERPRINT_LEN);
					System.arraycopy(Compression.intToBytes(markersInEachFile[markerFileIndex].length), 0, markerDataWriteBufferParameter, MARKDATA_MARKNAMELEN_START, MARKDATA_MARKNAMELEN_LEN);
					System.arraycopy(markersInEachFile[markerFileIndex], 0, markerDataWriteBufferParameter, MARKDATA_MARKNAME_START, markersInEachFile[markerFileIndex].length);
					markerFile = new RandomAccessFile(markerFilenames[markerFileIndex], "rw");
					markerFile.write(markerDataWriteBufferParameter);
		
		//			markFileOutliersWriteBufferVecAdj = new Vector<Byte>();
					markFileWriteBufferOutliersAdj = new Hashtable<String, Float>();
					markerFileBufferChunkIndex = 0;
					for(int i=0; i<timesEachSampleFileBeRead; i++) {
						log.report("Iteration "+i+"\t"+ext.reportMemoryUsage());
						timer2 = new Date().getTime();
						timer5 = 0;
						timer6 = 0;
						if ((i+1)==timesEachSampleFileBeRead && numMarkersLastReadRound!=0) {
							numChunksWriteBuffer = (int) Math.ceil((double)numMarkersLastReadRound / (double)numMarksPerBufferChunk);
		//					markFileOutliersWriteBufferVec = new Vector[numChunksInWriteBuffer];
							markFileWriteBufferOutliers = new Hashtable[numChunksWriteBuffer];
							if ( numChunksWriteBuffer>1 ) {
								markerDataWriteBuffer = null;
								markerDataWriteBuffer = new byte[numChunksWriteBuffer][numMarksPerBufferChunk * numBytesPerMarker];
								numMarkersLastReadRoundLastBuffer = numMarkersLastReadRound % numMarksPerBufferChunk;
								if ( numMarkersLastReadRoundLastBuffer!=0 ) {
									markerDataWriteBuffer[numChunksWriteBuffer-1] = new byte[numMarkersLastReadRoundLastBuffer * numBytesPerMarker];
								}
							} else {
								markerDataWriteBuffer = null;
								markerDataWriteBuffer = new byte[numChunksWriteBuffer][numMarkersLastReadRound * numBytesPerMarker];
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
					    	sampleFile.seek(Sample.PARAMETER_SECTION_BYTES + markerIndex * bytesPerSampMark);
					    	sampleFile.read(sampleDataReadBuffer);
					    	sampleFile.seek(Sample.PARAMETER_SECTION_BYTES + allMarkersInProj.length * bytesPerSampMark);
					    	numSampFileOutliers = sampleFile.readInt();
					    	if (numSampFileOutliers > 0) {
					    		sampFileOutliersReadBuffer = new byte[numSampFileOutliers];
		//			    		outOfRangeValuesInSampleFile = new float[numOutOfRangeValuesInSampleFile];
					    		sampleFile.read(sampFileOutliersReadBuffer);
		//			    		outOfRangeValuesInSampleFile = new Hashtable<String, Float>();
					    		outOfRangeValuesInSampleFile = (Hashtable<String, Float>) Compression.bytesToObj(sampFileOutliersReadBuffer);
					    	}
					    	indexInSampleDataReadBuffer = 0;
						    for (int k=0; k < numChunksWriteBuffer; k++) {
						    	timer4 = new Date().getTime();
						    	timer5 += (new Date().getTime() - timer4);
						    	timer4 = new Date().getTime();
								locationInWriteBufferChunk = j * bytesPerSampMark;
						    	for (int l=0; l<numMarksPerBufferChunk && markerIndex<allMarkersInProj.length; l++) {
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
						log.report("Load in sample files - round " + i + ": " + timer2/1000 + " secs = "+timer5/1000+" (file.read) + "+timer6/1000+" (array shuffle)");
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
								log.report("Write marker file " + markerFileIndex + ":\t"+ timer3/1000 + " sec."); // \toutlier #: " + (markFileOUtliersWriteBufferArray.length/4));
								markerFileBufferChunkIndex = 0;
								markerFileIndex ++;
								markFileWriteBufferOutliersAdj = new Hashtable<String, Float>();
								if ((i+1) != timesEachSampleFileBeRead || (j+1) != markerDataWriteBuffer.length) {
									markerDataWriteBufferParameter = new byte[MARKDATA_PARAMETER_TOTAL_LEN + markersInEachFile[markerFileIndex].length];
									System.arraycopy(Compression.intToBytes(allSamplesInProj.length), 0, markerDataWriteBufferParameter, MARKDATA_NUMSAMPS_START, MARKDATA_NUMSAMPS_LEN);
									System.arraycopy(Compression.intToBytes((markerFileIndex+1)==numMarkerFiles? (allMarkersInProj.length % maxNumMarkersPerFile) : maxNumMarkersPerFile), 0, markerDataWriteBufferParameter, MARKDATA_NUMMARKS_START, MARKDATA_NUMMARKS_LEN);
									markerDataWriteBufferParameter[MARKDATA_NULLSTATUS_START] = nullStatus;
									System.arraycopy(Compression.longToBytes(MarkerSet.fingerprint(allSamplesInProj)), 0, markerDataWriteBufferParameter, MARKDATA_FINGERPRINT_START, MARKDATA_FINGERPRINT_LEN);
									System.arraycopy(Compression.intToBytes(markersInEachFile[markerFileIndex].length), 0, markerDataWriteBufferParameter, MARKDATA_MARKNAMELEN_START, MARKDATA_MARKNAMELEN_LEN);
									System.arraycopy(markersInEachFile[markerFileIndex], 0, markerDataWriteBufferParameter, MARKDATA_MARKNAME_START, markersInEachFile[markerFileIndex].length);
									markerFile = new RandomAccessFile(markerFilenames[markerFileIndex], "rw");
									markerFile.write(markerDataWriteBufferParameter);
								}
								timer3 = 0;
							}
						}
						firstMarkerInCurrentFile += numMarkersWrtBuffer;
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
		//				Files.writeSerial(allOutliers, proj.getDir(Project.MARKER_DATA_DIRECTORY) + "outliers" + MarkerData.MARKER_DATA_FILE_EXTENSION);
						Files.writeSerial(allOutliers, proj.getDir(Project.MARKER_DATA_DIRECTORY) + "outliers.ser");
					}
					log.report("Write marker file " + markerFileIndex + ":\t" + timer3/1000 + " sec."); // \toutlier #: " + (markFileOUtliersWriteBufferArray.length/4));
					
				} catch (FileNotFoundException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (IOException e) {
					e.printStackTrace();
				} catch (ClassNotFoundException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			} catch (OutOfMemoryError oome) {
				memoryReserve *= 1.1;
				oomeLoops++;
			}
			
		} 
		timer1 = (new Date().getTime() - timer1);
		log.report("Finished transposing data. Total Time used: "+ (timer1/3600000) +" hrs "+ (timer1/60000 - 60*(timer1/3600000)) +" mins "+ (timer1/1000 - 60*(timer1/60000)) +" secs.\n\n");
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

	public static void backupOlderRafs(String directory, final String nameExtOrSuffixOfFilesToBackup) {
		backupOlderRafs(directory, new String[] {nameExtOrSuffixOfFilesToBackup}, true);
	}

	public static void backupOlderRafs(String directory, final String nameExtOrSuffixOfFilesToBackup, boolean backupToFolderOrRename) {
		backupOlderRafs(directory, new String[] {nameExtOrSuffixOfFilesToBackup}, backupToFolderOrRename);
	}

	public static void backupOlderRafs(String fileFullPath, byte i) {
//		new File(fileFullPath).renameTo(new File(ext.parseDirectoryOfFile(fileFullPath) + ext.rootOf(fileFullPath) + "_" + (new SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date())) + "." + fileFullPath.substring(fileFullPath.lastIndexOf(".") + 1, fileFullPath.length())));
		new File(fileFullPath).renameTo(new File(ext.parseDirectoryOfFile(fileFullPath) + ext.rootOf(fileFullPath) + "_Backup" + i + "." + fileFullPath.substring(fileFullPath.lastIndexOf(".") + 1, fileFullPath.length())));
}

	/*
	 * Move existing files of a specific file name extension to a sub-folder named "Backup.?"
	 * @param directory	the folder where the existing files are
	 * @param fileNameExtension	the file name extension, including the ".". For example, ".scat", ".xls", ".csv", and etc. 
	 */
	public static byte backupOlderRafs(String directory, final String[] nameExtOrSuffixOfFilesToBackup, boolean backupToFolderOrRename) {
		File[] files;
		byte count;
		String filename;

		// List of files to be backed up.
		files = new File(directory).listFiles(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				boolean result = false;
				for (int i=0; i<nameExtOrSuffixOfFilesToBackup.length; i++) {
					if (filename.endsWith(nameExtOrSuffixOfFilesToBackup[i])) {
						result = true;
						break;
					}
				}
				return result;
			}
		});

		count = 0;
		if (files.length>0) {
			// Create a new backup folder.
			if (backupToFolderOrRename) {
				do {
					count++;
				} while (new File(directory+"Backup."+count).exists());
				new File(directory+"Backup."+count+"/").mkdirs();
				
				// Move files to backup folder
				for (int i = 0; i<files.length; i++) {
					files[i].renameTo(new File(directory + "Backup." + count + "/" + files[i].getName()));
				}
	
				System.out.println("Older version of the data in " + directory + " has been moved to " + directory + "Backup." + count + "/ To save disk space, please manually delete these files.");
			} else {
				for (int i=0; i<files.length; i++) {
					filename = files[i].getName();
					files[i].renameTo(new File(files[i].getParent() + ext.rootOf(filename) + "_" + (new SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date())) + "." + filename.substring(filename.lastIndexOf(".") + 1, filename.length())));
				}
			}
		}
		return count;
	}

//	/**
//	 * Create MarkerLookup when creating new .
//	 * @param proj
//	 */
//	public static void createMarkerLookup(Project proj) {
//		for (int i=0; i<numMarkerFiles; i++) {
//			markerFilenames[i] = proj.getDir(Project.MARKER_DATA_DIRECTORY) + "markers." + i + MarkerData.MARKER_DATA_FILE_EXTENSION;
//			numMarkersCurrentLoop = Math.min(maxNumMarkersPerFile, allMarkersInProj.length - markerIndex);
//			markersInEachFile1 = new String[numMarkersCurrentLoop];
////			numMarkersInEachFile[i] =; 
//			for (int j=0; j<numMarkersCurrentLoop; j++) {
//				markerLookupHash.put(allMarkersInProj[markerIndex], "markers." + i + MarkerData.MARKER_DATA_FILE_EXTENSION + "\t" + j);
//				markersInEachFile1[j] = allMarkersInProj[markerIndex];
//				markerIndex ++;
//			}
//			try {
//				markersInEachFile[i] = Compression.objToBytes(markersInEachFile1);
//			} catch (IOException e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			}
//		}
//	}

	
	/**
	 * Create MarkerLookup from transposed data files in RAF and half-precision format
	 * @param proj
	 */
	public static void recreateMarkerLookup(Project proj) {
		Hashtable<String,String> hash = new Hashtable<String,String>();
		String[] files, markerNames;
		byte[] readBuffer;
		RandomAccessFile currentFile;
		long time;

		time = new Date().getTime();
		System.out.println("Creating MarkerLookup file");
		files = new File(proj.getDir(Project.MARKER_DATA_DIRECTORY)).list(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				return filename.endsWith(MarkerData.MARKER_DATA_FILE_EXTENSION) && filename.startsWith("marker");
			}
		});
		if (files==null) {
			System.err.println("Error - failed to create MarkerLookup -- directory does not exist: "+proj.getDir(Project.MARKER_DATA_DIRECTORY));
		} else if (files.length==0) {
			System.err.println("Error - failed to create MarkerLookup -- no "+MarkerData.MARKER_DATA_FILE_EXTENSION+" files available");
		} else {
			for (int i = 0; i<files.length; i++) {
				try {
					currentFile = new RandomAccessFile(proj.getDir(Project.MARKER_DATA_DIRECTORY) + files[i], "r");
					currentFile.seek(MARKDATA_MARKNAMELEN_START);
					readBuffer = new byte[currentFile.readInt()];
					currentFile.read(readBuffer);
					markerNames = (String[]) Compression.bytesToObj(readBuffer);
					for (int j = 0; j<markerNames.length; j++) {
						hash.put(markerNames[j], files[i]+"\t"+j);
//						markerLookupHash.put(allMarkersInProj[markerIndex], "markers." + i + MarkerData.MARKER_DATA_FILE_EXTENSION + "\t" + j);
					}
				} catch (FileNotFoundException e) {
					e.printStackTrace();
				} catch (Exception e) {
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
//		String filename = Project.DEFAULT_PROJECT;
		String filename = "C:/workspace/Genvisis/projects/gedi_exome.properties";
		// String demo = "flagged_results.txt";
//		boolean transpose = false;
		boolean transpose = true;
		int maxFileSize = 0;
		boolean lookup = false;
//		boolean lookup = true;

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
				transposeData(proj, maxFileSize, false, null);
			} else if (lookup) {
				recreateMarkerLookup(proj);
			} else {
				System.err.println("Error - no selection was made");
				System.out.println(usage);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
