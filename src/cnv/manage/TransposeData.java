// -Xms1024M -Xmx1024M     or even better: -Xmx15g
package cnv.manage;

import java.io.*;
import java.text.SimpleDateFormat;
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
		String[] allSampleNamesInProj;
		String[] allMarkerNamesInProj;
		int numBytesPerMarker;
		int numMarkersWrtBuffer;
		int numMarkersLastReadRoundLastBuffer = 0;
        int numMarksPerBufferChunk;
		int maxNumMarkersPerFile;
		int numWrBufferChunks;
		int numBufferChunksEachMarkerFile;
		int numMarkerFiles;
		int numMarkersLastReadRound;
		int timesEachSampleFileBeRead;
		String[] markerFilenames;
		int markerFileIndex;
		int indexWrBufferChunk;
		int numBufferChunksNeededCurrentMarkFile;
		int counterMarkerFileBufferChunks;
		int indexFirstMarkerCurrentIteration;
		int markerIndex;
		boolean isFileClosed;
		RandomAccessFile[] sampleFiles = null;
		RandomAccessFile sampleFile;
        RandomAccessFile markerFile = null;
        byte[][] writeBuffer = null;
        byte[] markerDataWriteBufferParameter;
        byte[] readBuffer;
		Hashtable<String,String> markerLookupHash = new Hashtable<String,String>();
		byte[][] markersInEachFile;
		String[] markersInEachFile1;
		Vector<Hashtable<String, Float>> markFileWriteBufferOutliers;
		byte[] markFileWriteBufferOutliersBytes;
		Hashtable<String, Float> allOutliers;
		byte nullStatus = Byte.MIN_VALUE;
        byte bytesPerSampMark;
        int numMarkersCurrentLoop;
        byte backupCount;
		boolean done, safe;
		long memoryReserve;
		int oomeLoops;
		long fingerPrint;
		long timerOverAll, timerLoadFiles, timerTransposeMemory, timerWriteFiles, timerTmp;
		SimpleDateFormat timeFormat;

		if (log == null) {
			log = new Logger(proj.getDir(Project.MARKER_DATA_DIRECTORY) + "Genvisis_" + (new SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date())) + ".log");
		}
		log.report("Transposing data of the project " + proj.getProjectDir());

		timeFormat = new SimpleDateFormat("HH:mm:ss.SSS");
		timeFormat.setTimeZone(TimeZone.getTimeZone("GMT"));
//		timeFormat = new SimpleDateFormat("'00':mm:ss.SSS");

		allSampleNamesInProj = proj.getSamples();
		allMarkerNamesInProj = proj.getMarkerNames();
		fingerPrint = MarkerSet.fingerprint(allSampleNamesInProj);
		nullStatus = Sample.getNullstatusFromRandomAccessFile(proj.getDir(Project.SAMPLE_DIRECTORY, true) + allSampleNamesInProj[0] + Sample.SAMPLE_DATA_FILE_EXTENSION, false);

		bytesPerSampMark = (byte) (Compression.BYTES_PER_SAMPLE_MARKER - (nullStatus & 0x01) - (nullStatus >>1 & 0x01) - (nullStatus >>2 & 0x01) - (nullStatus >>3 & 0x01) - (nullStatus >>4 & 0x01) - (nullStatus >>5 & 0x01) - (nullStatus >>6 & 0x01));
		numBytesPerMarker = allSampleNamesInProj.length * bytesPerSampMark;
		if (new File(proj.getProjectDir()).getFreeSpace() <= (allSampleNamesInProj.length * (long)allMarkerNamesInProj.length * bytesPerSampMark)) {
			log.reportError("Not enough disk space for all the new data to be created. Available: " + ext.prettyUpSize(new File(proj.getProjectDir()).getFreeSpace(), 1) + "; Required: " + ext.prettyUpSize(new File(proj.getProjectDir()).getFreeSpace(), 1) + ").");
			return;
		}
		if (maxMarkerFileSize == 0) {
			maxMarkerFileSize = Integer.MAX_VALUE;
		}

		done = false;
		safe = false;
		oomeLoops = 0;
		timerOverAll = new Date().getTime();
		memoryReserve = (long)((double)Runtime.getRuntime().maxMemory() * 0.25);
		while (!done) {
			try {
				numMarkersWrtBuffer = Math.min((int)(((double)Runtime.getRuntime().maxMemory() - memoryReserve) / numBytesPerMarker), allMarkerNamesInProj.length);
				if (maxMarkerFileSize > Integer.MAX_VALUE) {
					numMarksPerBufferChunk = (int) Math.min((double) Integer.MAX_VALUE / numBytesPerMarker, (double) numMarkersWrtBuffer);
				} else {
					numMarksPerBufferChunk = (int) Math.min((double) maxMarkerFileSize / numBytesPerMarker, (double) numMarkersWrtBuffer);
				}
				numWrBufferChunks = numMarkersWrtBuffer/numMarksPerBufferChunk;
				numMarkersWrtBuffer = numWrBufferChunks * numMarksPerBufferChunk;
				timesEachSampleFileBeRead = (int) Math.ceil((double)allMarkerNamesInProj.length / (double)numMarkersWrtBuffer);
				numMarkersLastReadRound = allMarkerNamesInProj.length % numMarkersWrtBuffer;
		
				maxNumMarkersPerFile = (int) Math.min((double)maxMarkerFileSize / (double)numBytesPerMarker, allMarkerNamesInProj.length);
		//		numBufferChunksPerFile = maxNumMarkersPerFile / numMarkersPerBufferChunk;
				numBufferChunksEachMarkerFile = (int) Math.ceil((double)maxNumMarkersPerFile / (double)numMarksPerBufferChunk);
				maxNumMarkersPerFile = numMarksPerBufferChunk * numBufferChunksEachMarkerFile;
				numMarkerFiles = (int) Math.ceil((double)allMarkerNamesInProj.length / (double)maxNumMarkersPerFile);
				counterMarkerFileBufferChunks = (int) Math.ceil((double)allMarkerNamesInProj.length / (double)numMarksPerBufferChunk);
				markerFilenames = new String[numMarkerFiles];
				markersInEachFile = new byte[numMarkerFiles][];
//				numMarkersInEachFile = new int[numMarkerFiles];
				backupCount = backupOlderFiles(proj.getDir(Project.MARKER_DATA_DIRECTORY, true, new Logger(), false), new String[] {MarkerData.MARKER_DATA_FILE_EXTENSION, "outliers.ser"}, true);
				backupOlderFile(proj.getFilename(Project.MARKERLOOKUP_FILENAME, false, false), backupCount);
				log.report( "--\nTotal Markers =\t" + allMarkerNamesInProj.length
						  + "\nTotal Samples =\t" + allSampleNamesInProj.length
						  + "\nMemory available =\t" + Runtime.getRuntime().maxMemory()/1024/1024/1024 + "." + (Runtime.getRuntime().maxMemory()/1024/1024/10 - Runtime.getRuntime().maxMemory()/1024/1024/1024*102) + " gb"
						  + "\nMarkers / writeBuffer =\t" + numMarkersWrtBuffer
						  + "\nMarkers / writeBufferChunk =\t" + numMarksPerBufferChunk
						  + "\nWriteBuffer chunks =\t" + numWrBufferChunks
						  + "\nMarkerFile size <=\t" + maxMarkerFileSize/1024/1024/1024 + "." + ((int)(maxMarkerFileSize/1024/1024/10.24) - (int)(maxMarkerFileSize/1024/1024/1024*102.4)) + " gb"
						  + "\nMarkers / markerFile =\t" + maxNumMarkersPerFile
						  + "\nMarkerFiles =\t" + numMarkerFiles
						  + "\nWriteBufferChunks / markerFile =\t" + numBufferChunksEachMarkerFile);


				markerIndex = 0;
				for (int i=0; i<numMarkerFiles; i++) {
					markerFilenames[i] = proj.getDir(Project.MARKER_DATA_DIRECTORY) + "markers." + i + MarkerData.MARKER_DATA_FILE_EXTENSION;
					numMarkersCurrentLoop = Math.min(maxNumMarkersPerFile, allMarkerNamesInProj.length - markerIndex);
					markersInEachFile1 = new String[numMarkersCurrentLoop];
					for (int j=0; j<numMarkersCurrentLoop; j++) {
						markerLookupHash.put(allMarkerNamesInProj[markerIndex], "markers." + i + MarkerData.MARKER_DATA_FILE_EXTENSION + "\t" + j);
						markersInEachFile1[j] = allMarkerNamesInProj[markerIndex];
						markerIndex ++;
					}
					markersInEachFile[i] = Compression.objToBytes(markersInEachFile1);
				}
		
				indexFirstMarkerCurrentIteration = 0;
				readBuffer = new byte[numMarkersWrtBuffer * bytesPerSampMark];
				
				writeBuffer = new byte[numWrBufferChunks][numMarksPerBufferChunk * numBytesPerMarker];
//				markFileOutliersWriteBufferVec = new Vector[numChunksInWriteBuffer];
//				log.report("After WriteBuffer\t"+ext.reportMemoryUsage());
				
				if (!safe) {
					safe = true;
					throw new OutOfMemoryError(); // one more round to be on the safe side
				}
				log.report("Memory optimization took "+oomeLoops+" round(s) ("+ext.prettyUpSize(memoryReserve, 1)+" was kept in reserve)");
				
				new MarkerLookup(markerLookupHash).serialize(proj.getFilename(Project.MARKERLOOKUP_FILENAME, false, false));

//				timerLoadFiles = 0;
				timerTmp = new Date().getTime();
				if (keepAllSampleFilesOpen) {
					sampleFiles = new RandomAccessFile[allSampleNamesInProj.length];
					for (int i=0; i<allSampleNamesInProj.length; i++) {
						sampleFiles[i] = new RandomAccessFile(proj.getDir(Project.SAMPLE_DIRECTORY, true) + allSampleNamesInProj[i] + Sample.SAMPLE_DATA_FILE_EXTENSION, "r");
					}
				}
//				timerLoadFiles += (new Date().getTime() - timerTmp);
	
				markerFileIndex = 0;
				indexWrBufferChunk = 0;
				numBufferChunksNeededCurrentMarkFile = numBufferChunksEachMarkerFile;
				isFileClosed = true;
	
//				allOutliers = new Hashtable<String, Float>();
//				markFileWriteBufferOutliersAdj = new Hashtable<String, Float>();
				if (new File(proj.getDir(Project.SAMPLE_DIRECTORY, true) + "outliers.ser").exists()) {
					allOutliers = (Hashtable<String, Float>) Files.readSerial(proj.getDir(Project.SAMPLE_DIRECTORY, true) + "outliers.ser");
				} else {
					allOutliers = new Hashtable<String, Float>();
				}
				markFileWriteBufferOutliers = getOutlierHashForEachFile(allOutliers, numMarkerFiles, maxNumMarkersPerFile, allSampleNamesInProj);

				timerWriteFiles = 0;
				log.report("--\ni\tLoad\tTranpose\tWrite");
				for(int i=0; i<timesEachSampleFileBeRead; i++) {
//					log.report(i + " - " + ext.reportMemoryUsage());
					if ((i+1)==timesEachSampleFileBeRead && numMarkersLastReadRound != 0) {
						readBuffer = new byte[numMarkersLastReadRound * bytesPerSampMark];
						numWrBufferChunks = (int) Math.ceil((double)numMarkersLastReadRound / (double)numMarksPerBufferChunk);
//						markFileWriteBufferOutliers = new Hashtable[numWrBufferChunks];
						if ( numWrBufferChunks>1 ) {
							writeBuffer = null;
							writeBuffer = new byte[numWrBufferChunks][numMarksPerBufferChunk * numBytesPerMarker];
							numMarkersLastReadRoundLastBuffer = numMarkersLastReadRound % numMarksPerBufferChunk;
							if ( numMarkersLastReadRoundLastBuffer!=0 ) {
								writeBuffer[numWrBufferChunks-1] = new byte[numMarkersLastReadRoundLastBuffer * numBytesPerMarker];
							}
						} else {
							writeBuffer = null;
							writeBuffer = new byte[numWrBufferChunks][numMarkersLastReadRound * numBytesPerMarker];
						}
					}
//					for (int j=0; j<markFileWriteBufferOutliers.length; j++) {
//						markFileWriteBufferOutliers[j] = new Hashtable<String, Float>();
//					}

					timerLoadFiles = 0;
					timerTransposeMemory = 0;
					for (int j=0; j<allSampleNamesInProj.length; j++) {
						timerTmp = new Date().getTime();
						if (! keepAllSampleFilesOpen) {
							sampleFile = new RandomAccessFile(proj.getDir(Project.SAMPLE_DIRECTORY, true) + allSampleNamesInProj[j] + Sample.SAMPLE_DATA_FILE_EXTENSION, "r");
						} else {
							sampleFile = sampleFiles[j];
						}

//						if (i == 0) {
//							loadSampleFileToWriteBuffer1(sampleFile, readBuffer, j, indexFirstMarkerCurrentIteration, bytesPerSampMark, allMarkersInProj.length, allOutliers);
//						} else {
//							loadSampleFileToWriteBuffer1(sampleFile, readBuffer, j, indexFirstMarkerCurrentIteration, bytesPerSampMark, allMarkersInProj.length, null);
//						}
						Sample.loadFromRandomAccessFileWithoutDecompress(sampleFile, readBuffer, true, j, indexFirstMarkerCurrentIteration, bytesPerSampMark, allMarkerNamesInProj.length, null);
//						Sample.loadFromRandomAccessFileWithoutDecompress(sampleFile, readBuffer, ! keepAllSampleFilesOpen, j, indexFirstMarkerCurrentIteration, bytesPerSampMark, allMarkerNamesInProj.length, null);

						timerLoadFiles += (new Date().getTime() - timerTmp);
						timerTmp = new Date().getTime();

						transposeBuffer(writeBuffer, readBuffer, bytesPerSampMark, indexFirstMarkerCurrentIteration, j, allSampleNamesInProj.length);

						timerTransposeMemory += (new Date().getTime() - timerTmp);

					    if (!keepAllSampleFilesOpen) {
					    	sampleFile.close();
					    }
					}
					log.report(i + "\t" + timeFormat.format(timerLoadFiles) + "\t" + timeFormat.format(timerTransposeMemory), false, true);

					while (numBufferChunksNeededCurrentMarkFile > 0 && (indexWrBufferChunk + numBufferChunksNeededCurrentMarkFile) <= writeBuffer.length) {
//						if (numBufferChunksNeededCurrentMarkFile == numBufferChunksEachMarkerFile) {
						if (markFileWriteBufferOutliers == null || markFileWriteBufferOutliers.elementAt(markerFileIndex).size() == 0) {
							markFileWriteBufferOutliersBytes = new byte[0];
						} else {
							markFileWriteBufferOutliersBytes = Compression.objToBytes(markFileWriteBufferOutliers.elementAt(markerFileIndex));
						}

						timerTmp = new Date().getTime();

						if (isFileClosed) {
							markerDataWriteBufferParameter = getWriteBufferParameterSection(allSampleNamesInProj.length, markerFileIndex == (numMarkerFiles - 1)? allMarkerNamesInProj.length % maxNumMarkersPerFile : maxNumMarkersPerFile, nullStatus, fingerPrint, markersInEachFile[markerFileIndex]);
							markerFile = new RandomAccessFile(markerFilenames[markerFileIndex], "rw");
//							writeBufferToRAF(markersInEachFile, indexWrBufferChunk, indexWrBufferChunk + numBufferChunksNeeded - 1, markerFile, markerDataWriteBufferParameter, markFileWriteBufferOutliersBytes);
							writeBufferToRAF(writeBuffer, indexWrBufferChunk, indexWrBufferChunk + numBufferChunksNeededCurrentMarkFile - 1, markerFile, markerDataWriteBufferParameter, markFileWriteBufferOutliersBytes);
							markerFile.close();
							counterMarkerFileBufferChunks -= numBufferChunksNeededCurrentMarkFile;
							indexWrBufferChunk += numBufferChunksNeededCurrentMarkFile;
							numBufferChunksNeededCurrentMarkFile = Math.min(numBufferChunksEachMarkerFile, counterMarkerFileBufferChunks);

							timerWriteFiles += (new Date().getTime() - timerTmp);
							log.report("\t" + timeFormat.format(timerWriteFiles), false, true);
							timerWriteFiles = 0;
						} else {
							writeBufferToRAF(writeBuffer, indexWrBufferChunk, indexWrBufferChunk + numBufferChunksNeededCurrentMarkFile - 1, markerFile, null, markFileWriteBufferOutliersBytes);
							markerFile.close();
							counterMarkerFileBufferChunks -= numBufferChunksNeededCurrentMarkFile;
							indexWrBufferChunk += numBufferChunksNeededCurrentMarkFile;
							numBufferChunksNeededCurrentMarkFile = Math.min(numBufferChunksEachMarkerFile, counterMarkerFileBufferChunks);
							isFileClosed = true;

						}
						markerFileIndex ++;
						timerWriteFiles += (new Date().getTime() - timerTmp);
						log.report("\t" + timeFormat.format(timerWriteFiles), false, true);
						timerWriteFiles = 0;
					}
					if (numBufferChunksNeededCurrentMarkFile > 0 && indexWrBufferChunk < writeBuffer.length) {
						timerTmp = new Date().getTime();

						counterMarkerFileBufferChunks -= (writeBuffer.length - indexWrBufferChunk);
						if (isFileClosed) {
							markerDataWriteBufferParameter = getWriteBufferParameterSection(allSampleNamesInProj.length, markerFileIndex == (numMarkerFiles - 1)? allMarkerNamesInProj.length % maxNumMarkersPerFile : maxNumMarkersPerFile, nullStatus, fingerPrint, markersInEachFile[markerFileIndex]);
							markerFile = new RandomAccessFile(markerFilenames[markerFileIndex], "rw");
							writeBufferToRAF(writeBuffer, indexWrBufferChunk, writeBuffer.length - 1, markerFile, markerDataWriteBufferParameter, null);
							isFileClosed = false;
							numBufferChunksNeededCurrentMarkFile -= (writeBuffer.length - indexWrBufferChunk);
						} else if ((indexWrBufferChunk + numBufferChunksNeededCurrentMarkFile) > writeBuffer.length) {
							writeBufferToRAF(writeBuffer, indexWrBufferChunk, writeBuffer.length - 1, markerFile, null, null);
							numBufferChunksNeededCurrentMarkFile -= (writeBuffer.length - indexWrBufferChunk);
						} else {
							if (markFileWriteBufferOutliers == null || markFileWriteBufferOutliers.elementAt(markerFileIndex).size() == 0) {
								markFileWriteBufferOutliersBytes = new byte[0];
							} else {
								markFileWriteBufferOutliersBytes = Compression.objToBytes(markFileWriteBufferOutliers.elementAt(markerFileIndex));
							}
							writeBufferToRAF(writeBuffer, indexWrBufferChunk, indexWrBufferChunk + numBufferChunksNeededCurrentMarkFile - 1, markerFile, null, markFileWriteBufferOutliersBytes);
							markerFile.close();
							isFileClosed = true;
							markerFileIndex ++;
							numBufferChunksNeededCurrentMarkFile = Math.min(numBufferChunksEachMarkerFile, counterMarkerFileBufferChunks);
						}

//						numBufferChunksNeededCurrentMarkFile = numBufferChunksEachMarkerFile + indexWrBufferChunk - writeBuffer.length;

						timerWriteFiles += (new Date().getTime() - timerTmp);
					}
					indexWrBufferChunk = 0;

					indexFirstMarkerCurrentIteration += numMarkersWrtBuffer;
					log.report("");
				}

				if (allOutliers != null && allOutliers.size() != 0) {
					Files.writeSerial(allOutliers, proj.getDir(Project.MARKER_DATA_DIRECTORY) + "outliers.ser");
				}

				done = true;

			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			} catch (OutOfMemoryError oome) {
				memoryReserve *= 1.1;
				oomeLoops++;
			}

		} 
		timerOverAll = (new Date().getTime() - timerOverAll);
		log.report("--\nFinished transposing data. Total Time used: "+ timeFormat.format(timerOverAll));
	}


	/**
	 * Rebuild MarkerLookup from transposed data files (.mkRAF).
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


	private static void transposeBuffer(byte[][] writeBuffer, byte[] readBuffer, byte bytesPerSampMark, int indexOfFirstMarkerInBuffer, int indexOfCurrentSample, int numSamplesInProj) {
		int numMarkersInChunk;
		int indexInReadBuffer;
		int indexInChunk;
		int step;

		indexInReadBuffer = 0;
		step = numSamplesInProj * bytesPerSampMark;
//		indexInChunk = indexOfCurrentSample * bytesPerSampMark;
		for (int i=0; i < writeBuffer.length; i++) {
			numMarkersInChunk = writeBuffer[i].length / bytesPerSampMark / numSamplesInProj;
			indexInChunk = indexOfCurrentSample * bytesPerSampMark;
//			if (i>0) {
//				indexInChunk -= writeBuffer[i-1].length;
//			}

			for (int j=0; j < numMarkersInChunk; j++) {
				for (int k=0; k < bytesPerSampMark; k++) {
					try {
						writeBuffer[i][indexInChunk + k] = readBuffer[indexInReadBuffer];
					} catch (ArrayIndexOutOfBoundsException e) {
						System.out.println("j:" + j + "\tnumMarkersInChunk:" + numMarkersInChunk);
						System.out.println("i:" + i + "\tj:"+j+"\tk:"+k+"\tindexInChunk: "+indexInChunk);
						System.out.println("writebuffer size: " + writeBuffer.length + " , "+writeBuffer[i].length+"\treadBuffer size: "+readBuffer.length);
						e.printStackTrace();
						System.exit(1);
					}
	    			indexInReadBuffer ++;
	    		}
				indexInChunk += step;
			}
		}

//    	indexInSampleDataReadBuffer = 0;
//	    for (int k=0; k < numChunksWriteBuffer; k++) {
//	    	timer4 = new Date().getTime();
//	    	timer5 += (new Date().getTime() - timer4);
//	    	timer4 = new Date().getTime();
//			locationInWriteBufferChunk = j * bytesPerSampMark;
//	    	for (int l=0; l<numMarksPerBufferChunk && markerIndex<allMarkersInProj.length; l++) {
//	    		for (int m=0; m < bytesPerSampMark; m++) {
//    				markerDataWriteBuffer[k][locationInWriteBufferChunk + m] = sampleDataReadBuffer[indexInSampleDataReadBuffer];
//	    			if (m==3 && (sampleDataReadBuffer[indexInSampleDataReadBuffer-1] == Compression.REDUCED_PRECISION_XY_OUT_OF_RANGE_BYTES[0]
//	    					     && sampleDataReadBuffer[indexInSampleDataReadBuffer] == Compression.REDUCED_PRECISION_XY_OUT_OF_RANGE_BYTES[1])) {
//	    				allOutliers.put(markerIndex + "\t" + j +"\tx", outOfRangeValuesInSampleFile.get(markerIndex + "\tx"));
//	    				markFileWriteBufferOutliers[k].put(markerIndex + "\t" + j + "\tx", outOfRangeValuesInSampleFile.get(markerIndex + "\tx"));
//	    			} else if (m==5 && (sampleDataReadBuffer[indexInSampleDataReadBuffer-1] == Compression.REDUCED_PRECISION_XY_OUT_OF_RANGE_BYTES[0]
//    					                && sampleDataReadBuffer[indexInSampleDataReadBuffer] == Compression.REDUCED_PRECISION_XY_OUT_OF_RANGE_BYTES[1])) {
//	    				allOutliers.put(markerIndex + "\t" + j + "\ty", outOfRangeValuesInSampleFile.get(markerIndex + "\ty"));
//	    				markFileWriteBufferOutliers[k].put(markerIndex + "\t" + j + "\ty", outOfRangeValuesInSampleFile.get(markerIndex + "\ty"));
//	    			} else if (m==10 && (sampleDataReadBuffer[indexInSampleDataReadBuffer-2] == Compression.REDUCED_PRECISION_LRR_OUT_OF_RANGE_BYTES[0]
//	    								 && sampleDataReadBuffer[indexInSampleDataReadBuffer-1] == Compression.REDUCED_PRECISION_LRR_OUT_OF_RANGE_BYTES[1]
//	    								 && sampleDataReadBuffer[indexInSampleDataReadBuffer] == Compression.REDUCED_PRECISION_LRR_OUT_OF_RANGE_BYTES[2])) {
//	    				allOutliers.put(markerIndex + "\t" + j + "\tlrr", outOfRangeValuesInSampleFile.get(markerIndex + "\tlrr"));
//	    				markFileWriteBufferOutliers[k].put(markerIndex + "\t" + j + "\tlrr", outOfRangeValuesInSampleFile.get(markerIndex + "\tlrr"));
//	    			}
////				    			if ((m==2 && (0x80 & sampleDataReadBuffer[indexInSampleDataReadBuffer])>0) || (m==4 && (0x80 & sampleDataReadBuffer[indexInSampleDataReadBuffer])>0)) {
////				    				allOutliers.put(allSamplesInProj[j]+"\t"+allMarkersInProj[markerIndex]+"\t"+(m==2?"x":"y"), outOfRangeValuesInSampleFile.get(allSamplesInProj[j]+"\t"+allMarkersInProj[markerIndex]+"\t"+(m==2?"x":"y")));
////				    				markFileWriteBufferOutliers[k].put(allSamplesInProj[j]+"\t"+allMarkersInProj[markerIndex]+"\t"+(m==2?"x":"y"), outOfRangeValuesInSampleFile.get(allSamplesInProj[j]+"\t"+allMarkersInProj[markerIndex]+"\t"+(m==2?"x":"y")));
////				    			}
//    				indexInSampleDataReadBuffer ++;
//		    	}
//	    		locationInWriteBufferChunk += numBytesPerMarker;
//				markerIndex ++;
//	    	}
//	    	timer6 += (new Date().getTime() - timer4);
//		}

	}

	public static void writeBufferToRAF(int numBufferChunksNeededCurrentMarkFile, int indexWrBufferChunk, byte[][] writeBuffer, Vector<Hashtable<String, Float>> markFileWriteBufferOutliers, boolean isFileClosed, int markerFileIndex, int numMarkerFiles, String[] allSampleNamesInProj, int counterMarkerFileBufferChunks, String[] allMarkerNamesInProj, String[] markerFilenames, int numBufferChunksEachMarkerFile) {
//		Logger log;
//		long timerTmp, timerWriteFiles;
//		byte[] markFileWriteBufferOutliersBytes, markerDataWriteBufferParameter;
//		RandomAccessFile markerFile;
//		SimpleDateFormat timeFormat;
//		
//		while (numBufferChunksNeededCurrentMarkFile > 0 && (indexWrBufferChunk + numBufferChunksNeededCurrentMarkFile) <= writeBuffer.length) {
////			if (numBufferChunksNeededCurrentMarkFile == numBufferChunksEachMarkerFile) {
//			if (markFileWriteBufferOutliers == null || markFileWriteBufferOutliers.elementAt(markerFileIndex).size() == 0) {
//				markFileWriteBufferOutliersBytes = new byte[0];
//			} else {
//				markFileWriteBufferOutliersBytes = Compression.objToBytes(markFileWriteBufferOutliers.elementAt(markerFileIndex));
//			}
//
//			timerTmp = new Date().getTime();
//
//			if (isFileClosed) {
//				markerDataWriteBufferParameter = getWriteBufferParameterSection(allSampleNamesInProj.length, markerFileIndex == (numMarkerFiles - 1)? allMarkerNamesInProj.length % maxNumMarkersPerFile : maxNumMarkersPerFile, nullStatus, fingerPrint, markersInEachFile[markerFileIndex]);
//				markerFile = new RandomAccessFile(markerFilenames[markerFileIndex], "rw");
////				writeBufferToRAF(markersInEachFile, indexWrBufferChunk, indexWrBufferChunk + numBufferChunksNeeded - 1, markerFile, markerDataWriteBufferParameter, markFileWriteBufferOutliersBytes);
//				writeBufferToRAF(writeBuffer, indexWrBufferChunk, indexWrBufferChunk + numBufferChunksNeededCurrentMarkFile - 1, markerFile, markerDataWriteBufferParameter, markFileWriteBufferOutliersBytes);
//				markerFile.close();
//				counterMarkerFileBufferChunks -= numBufferChunksNeededCurrentMarkFile;
//				indexWrBufferChunk += numBufferChunksNeededCurrentMarkFile;
//				numBufferChunksNeededCurrentMarkFile = Math.min(numBufferChunksEachMarkerFile, counterMarkerFileBufferChunks);
//
//				timerWriteFiles += (new Date().getTime() - timerTmp);
//				log.report("\t" + timeFormat.format(timerWriteFiles), false, true);
//				timerWriteFiles = 0;
//			} else {
//				writeBufferToRAF(writeBuffer, indexWrBufferChunk, indexWrBufferChunk + numBufferChunksNeededCurrentMarkFile - 1, markerFile, null, markFileWriteBufferOutliersBytes);
//				markerFile.close();
//				counterMarkerFileBufferChunks -= numBufferChunksNeededCurrentMarkFile;
//				indexWrBufferChunk += numBufferChunksNeededCurrentMarkFile;
//				numBufferChunksNeededCurrentMarkFile = Math.min(numBufferChunksEachMarkerFile, counterMarkerFileBufferChunks);
//				isFileClosed = true;
//
//			}
//			markerFileIndex ++;
//			timerWriteFiles += (new Date().getTime() - timerTmp);
//			log.report("\t" + timeFormat.format(timerWriteFiles), false, true);
//			timerWriteFiles = 0;
//		}
//		if (numBufferChunksNeededCurrentMarkFile > 0 && indexWrBufferChunk < writeBuffer.length) {
//			timerTmp = new Date().getTime();
//
//			counterMarkerFileBufferChunks -= (writeBuffer.length - indexWrBufferChunk);
//			if (isFileClosed) {
//				markerDataWriteBufferParameter = getWriteBufferParameterSection(allSampleNamesInProj.length, markerFileIndex == (numMarkerFiles - 1)? allMarkerNamesInProj.length % maxNumMarkersPerFile : maxNumMarkersPerFile, nullStatus, fingerPrint, markersInEachFile[markerFileIndex]);
//				markerFile = new RandomAccessFile(markerFilenames[markerFileIndex], "rw");
//				writeBufferToRAF(writeBuffer, indexWrBufferChunk, writeBuffer.length - 1, markerFile, markerDataWriteBufferParameter, null);
//				isFileClosed = false;
//				numBufferChunksNeededCurrentMarkFile -= (writeBuffer.length - indexWrBufferChunk);
//			} else if ((indexWrBufferChunk + numBufferChunksNeededCurrentMarkFile) > writeBuffer.length) {
//				writeBufferToRAF(writeBuffer, indexWrBufferChunk, writeBuffer.length - 1, markerFile, null, null);
//				numBufferChunksNeededCurrentMarkFile -= (writeBuffer.length - indexWrBufferChunk);
//			} else {
//				if (markFileWriteBufferOutliers == null || markFileWriteBufferOutliers.elementAt(markerFileIndex).size() == 0) {
//					markFileWriteBufferOutliersBytes = new byte[0];
//				} else {
//					markFileWriteBufferOutliersBytes = Compression.objToBytes(markFileWriteBufferOutliers.elementAt(markerFileIndex));
//				}
//				writeBufferToRAF(writeBuffer, indexWrBufferChunk, indexWrBufferChunk + numBufferChunksNeededCurrentMarkFile - 1, markerFile, null, markFileWriteBufferOutliersBytes);
//				markerFile.close();
//				isFileClosed = true;
//				markerFileIndex ++;
//				numBufferChunksNeededCurrentMarkFile = Math.min(numBufferChunksEachMarkerFile, counterMarkerFileBufferChunks);
//			}
//
////			numBufferChunksNeededCurrentMarkFile = numBufferChunksEachMarkerFile + indexWrBufferChunk - writeBuffer.length;
//
//			timerWriteFiles += (new Date().getTime() - timerTmp);
//		}
//		indexWrBufferChunk = 0;
//
//		indexFirstMarkerCurrentIteration += numMarkersWrtBuffer;
//		log.report("");
	}

	public static void writeBufferToRAF(byte[][] buffer, String fileName, byte[] head, byte[] tail) {
		writeBufferToRAF(buffer, 0, buffer.length, fileName, head, tail);
	}

	public static void writeBufferToRAF(byte[][] buffer, int start, int end, String fileName, byte[] head, byte[] tail) {
		RandomAccessFile markerFile;
		try {
			markerFile = new RandomAccessFile(fileName, "rw");
			writeBufferToRAF(buffer, start, end, markerFile, head, tail);
			markerFile.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
		
	/*
	 * Write buffer to Random Access File
	 * tail == null means this file is not yet ended;
	 * teil.length = 0 means this file is ended, but there is no content in the tail.
	 */
	public static void writeBufferToRAF(byte[][] buffer, int indexOfStart, int indexOfEnd, RandomAccessFile markerFile, byte[] head, byte[] tail) {
		
		if (buffer==null || indexOfStart<0 || indexOfEnd>=buffer.length || indexOfEnd<indexOfStart) {
			System.err.println("Transpose Data encoutered the following error: buffer be null, or start index of buffer is negative, or end index is less than the start index, or end index is over the buffer size.");
			System.exit(1);
		}
		
		try {
			if (head != null && head.length != 0) {
					markerFile.write(head);
			}
			for (int i=indexOfStart; i<=indexOfEnd; i++) {
	//	    	timer4 = new Date().getTime();
				markerFile.write(buffer[i]);
			}
			if (tail!=null) {
				markerFile.writeInt(tail.length);
				if (tail.length!=0) {
					markerFile.write(tail);
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private static Vector<Hashtable<String, Float>> getOutlierHashForEachFile(Hashtable<String, Float> allOutliers, int numMarkerFiles, int numMarkersInEachFile, String[] sampleNames) {
		Vector<Hashtable<String, Float>> result;
		Enumeration<String> keys;
		String key;
		String[] line;
		int sampleIndex = -1;
		
		result = new Vector<Hashtable<String,Float>>(numMarkerFiles);
		for (int j=0; j<numMarkerFiles; j++) {
			result.add(new Hashtable<String, Float>(allOutliers.size()));
		}

		keys = allOutliers.keys();
		while (keys.hasMoreElements()) {
			key = keys.nextElement();
			line = key.split("\t");
			for (int i=0; i<sampleNames.length; i++) {
				if (sampleNames[i].equals(line[1])) {
					sampleIndex = i;
					break;
				}
			}
			if (sampleIndex == -1) {
				System.err.println("Error - Cannot find the sample " + line[1] + "in the project's sample list.");
			}
			result.elementAt(Integer.parseInt(line[0]) / numMarkersInEachFile).put((Integer.parseInt(line[0]) % numMarkersInEachFile) + "\t" + sampleIndex, allOutliers.get(key));
		}

		return result;
	}

//	private static byte[] getWriteBufferParameterSection(int numSampsInProj, int numMarkersInCurrentFile, byte nullStatus, long fingerPrint, String[] currentFileMarkerNames) {
//		byte[] markerNamesBytes;
//
//		try {
//			markerNamesBytes = Compression.objToBytes(currentFileMarkerNames);
//		} catch (IOException e) {
//			e.printStackTrace();
//			return null;
//		}
//
//		return getWriteBufferParameterSection(numSampsInProj, numMarkersInCurrentFile, nullStatus, fingerPrint, markerNamesBytes);
//	}
//
	private static byte[] getWriteBufferParameterSection(int numSampsInProj, int numMarkersInCurrentFile, byte nullStatus, long fingerPrint, byte[] currentFileMarkerNamesInBytes) {
		byte[] markerFileHead;

		markerFileHead = new byte[MARKDATA_PARAMETER_TOTAL_LEN + currentFileMarkerNamesInBytes.length];
		System.arraycopy(Compression.intToBytes(numSampsInProj), 0, markerFileHead, MARKDATA_NUMSAMPS_START, MARKDATA_NUMSAMPS_LEN);
		System.arraycopy(Compression.intToBytes(numMarkersInCurrentFile), 0, markerFileHead, MARKDATA_NUMMARKS_START, MARKDATA_NUMMARKS_LEN);
		markerFileHead[MARKDATA_NULLSTATUS_START] = nullStatus;
		System.arraycopy(Compression.longToBytes(fingerPrint), 0, markerFileHead, MARKDATA_FINGERPRINT_START, MARKDATA_FINGERPRINT_LEN);
		System.arraycopy(Compression.intToBytes(currentFileMarkerNamesInBytes.length), 0, markerFileHead, MARKDATA_MARKNAMELEN_START, MARKDATA_MARKNAMELEN_LEN);
		System.arraycopy(currentFileMarkerNamesInBytes, 0, markerFileHead, MARKDATA_MARKNAME_START, currentFileMarkerNamesInBytes.length);
	
		return markerFileHead;
	}

	/*
	 * Delete existing files of a specific file name extension.
	 * @param directory	the folder where the existing files are
	 * @param fileNameExtension	the file name extension, including the ".". For example, ".scat", ".xls", ".csv", and etc. 
	 */
//	public static void deleteOlderRafs(String directory, final String fileNameExtesion) {
//		deleteOlderRafs(directory, new String[] {fileNameExtesion});
//	}

//	public static void deleteOlderRafs(String directory, final String[] fileNameExtesions) {
//		deleteOlderRafs(directory, null, fileNameExtesions);
//	}

//	public static void deleteOlderRafs(String directory, final String[] fileNamePrefixes, final String[] fileNameExtesions) {
//		deleteOlderRafs(directory, fileNamePrefixes, fileNameExtesions, false);
//	}

//	public static void deleteOlderRafs(String directory, final String[] fileNamePrefixes, final String[] fileNameExtesions, final boolean prefixAndSuffix) {
//		deleteOlderRafs(directory, fileNamePrefixes, fileNameExtesions, prefixAndSuffix, null);
//	}

	public static void deleteOlderRafs(String directory, final String[] fileNamePrefixes, final String[] fileNameExtesions, final boolean prefixAndSuffix, final String[] excludes) {
		File[] files;

		// List of files to be backed up.
		files = new File(directory).listFiles(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				boolean result = false;
				filename = ext.removeDirectoryInfo(filename);
				if (prefixAndSuffix) {
					if (fileNamePrefixes.length != fileNameExtesions.length) {
						System.err.println("Array fileNamePrefixes and Array fileNameExtesions must have the same number of elements.");
						System.exit(1);
					}
					for (int i=0; !result && fileNameExtesions!= null && i<fileNameExtesions.length; i++) {
						if (filename.startsWith(fileNamePrefixes[i]) && filename.endsWith(fileNameExtesions[i])) {
							result = true;
						}
					}
				} else {
					for (int i=0; !result && fileNameExtesions!= null && i<fileNameExtesions.length; i++) {
						if (filename.endsWith(fileNameExtesions[i])) {
							result = true;
						}
					}
					for (int i=0; !result && fileNamePrefixes!= null && i<fileNamePrefixes.length; i++) {
						if (filename.startsWith(fileNamePrefixes[i])) {
							result = true;
						}
					}
				}
				if (result) {
					for (int i=0; i<excludes.length; i++) {
						if (ext.removeDirectoryInfo(filename).equals(excludes[i])) {
							result = false;
							break;
						}
					}
				}
				return result;
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
//	public static void backupOlderRafs(String directory, final String nameExtOrSuffixOfFilesToBackup) {
//		backupOlderFiles(directory, new String[] {nameExtOrSuffixOfFilesToBackup}, true);
//	}

//	public static void backupOlderRafs(String directory, final String nameExtOrSuffixOfFilesToBackup, boolean backupToFolderOrRename) {
//		backupOlderFiles(directory, new String[] {nameExtOrSuffixOfFilesToBackup}, backupToFolderOrRename);
//	}

	public static byte backupOlderFiles(String directory, final String[] nameExtOrSuffixOfFilesToBackup, boolean backupToFolderOrRename) {
		return backupOlderFiles(directory, null, nameExtOrSuffixOfFilesToBackup, backupToFolderOrRename);
	}

	public static byte backupOlderFiles(String directory, final String[] namePrefixOfFilesToBackup, final String[] nameExtOrSuffixOfFilesToBackup, boolean backupToFolderOrRename) {
		File[] files;
		byte count;
		String filename;

		// List of files to be backed up.
		files = new File(directory).listFiles(new FilenameFilter() {
			public boolean accept(File file, String filename) {
				boolean result = false;
				filename = ext.removeDirectoryInfo(filename);
				for (int i=0; !result && nameExtOrSuffixOfFilesToBackup!= null && i<nameExtOrSuffixOfFilesToBackup.length; i++) {
					if (filename.endsWith(nameExtOrSuffixOfFilesToBackup[i])) {
						result = true;
					}
				}
				for (int i=0; !result && namePrefixOfFilesToBackup!= null && i<namePrefixOfFilesToBackup.length; i++) {
					if (filename.startsWith(namePrefixOfFilesToBackup[i])) {
						result = true;
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

	public static void backupOlderFile(String fileFullPath, byte i) {
//		new File(fileFullPath).renameTo(new File(ext.parseDirectoryOfFile(fileFullPath) + ext.rootOf(fileFullPath) + "_" + (new SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date())) + "." + fileFullPath.substring(fileFullPath.lastIndexOf(".") + 1, fileFullPath.length())));
		new File(fileFullPath).renameTo(new File(ext.parseDirectoryOfFile(fileFullPath) + ext.rootOf(fileFullPath) + "_Backup" + i + "." + fileFullPath.substring(fileFullPath.lastIndexOf(".") + 1, fileFullPath.length())));
	}

	@SuppressWarnings("unchecked")
	public static MarkerData[] loadFromRAF(String markerFilename, int[] indeciesOfMarkersToLoad) {
		MarkerData[] result = null;
		RandomAccessFile file;
        int numBytesPerMarker;
        float[] gcs = null;
        float[] xs = null;
        float[] ys = null;
        float[] bafs = null;
        float[] lrrs = null;
        byte[] abGenotypes = null;
        String[] alleleMappings = null;
        byte[] genotypeTmp;
        long seekLocation;
        byte[][] readBuffer;
        int indexReadBuffer;
        byte[] parameters;
        int markernamesSectionLength;
        Hashtable<String, Float> outOfRangeValues = null;
        byte nullStatus = 0;
        byte bytesPerSampMark = 0;
        int indexStart;
        int numSamples;
        int numMarkersInThisFile;
        long fingerPrint;
        String[] markerNames;
        int lengthOfOutOfRangeHashtable;
//        boolean loadGC, loadXY, loadBAF, loadLRR;

        try {
			file = new RandomAccessFile(markerFilename, "r");
        	parameters = new byte[TransposeData.MARKDATA_PARAMETER_TOTAL_LEN];
			file.read(parameters);
			numSamples = Compression.bytesToInt(parameters, TransposeData.MARKDATA_NUMSAMPS_START);
			numMarkersInThisFile = Compression.bytesToInt(parameters, TransposeData.MARKDATA_NUMMARKS_START);
			nullStatus = parameters[TransposeData.MARKDATA_NULLSTATUS_START];
//			loadGC = ((nullStatus >>6 & 0x01) == 1);
//			loadXY = ((nullStatus >>5 & 0x01) == 1 || (nullStatus >>4 & 0x01) == 1);
//			loadBAF = ((nullStatus >>3 & 0x01) == 1);
//			loadLRR = ((nullStatus >>2 & 0x01) == 1);
			bytesPerSampMark = (byte) (Compression.BYTES_PER_SAMPLE_MARKER - (nullStatus & 0x01) - (nullStatus >>1 & 0x01) - (nullStatus >>2 & 0x01) - (nullStatus >>3 & 0x01) - (nullStatus >>4 & 0x01) - (nullStatus >>5 & 0x01) - (nullStatus >>6 & 0x01));
			numBytesPerMarker = bytesPerSampMark * numSamples;
			fingerPrint = Compression.bytesToLong(parameters, MARKDATA_FINGERPRINT_START);
			markernamesSectionLength = Compression.bytesToInt(parameters, MARKDATA_MARKNAMELEN_START);

			parameters = new byte[markernamesSectionLength];
			file.read(parameters);
			markerNames = (String[]) Compression.bytesToObj(parameters);

			result = new MarkerData[indeciesOfMarkersToLoad.length];
			readBuffer = new byte[indeciesOfMarkersToLoad.length][numBytesPerMarker];

//	        for (int i=indexStartMarker; i<indexEndMarker; i++) {
	        for (int i=0; i<indeciesOfMarkersToLoad.length; i++) {
		        //if(indeciesInFile[i-1]+1==indeciesInFile[i]) {No need to seek}
	        	if (indeciesOfMarkersToLoad[i] < 0 || indeciesOfMarkersToLoad[i] >= numMarkersInThisFile) {
					System.err.println("Skipped the marker index " + indeciesOfMarkersToLoad[i] + ", because it is out of range.");
	        	} else {
	        		seekLocation = (long)TransposeData.MARKDATA_PARAMETER_TOTAL_LEN + (long)markernamesSectionLength + indeciesOfMarkersToLoad[i] * (long)numBytesPerMarker;
					file.seek(seekLocation);
					file.read(readBuffer[i]);
	        	}
			}

	        // TODO this is read every time, wouldn't it be faster to use the serialized version?
	        file.seek((long)TransposeData.MARKDATA_PARAMETER_TOTAL_LEN + (long)markernamesSectionLength + (long) numMarkersInThisFile * (long)numBytesPerMarker);
			lengthOfOutOfRangeHashtable = file.readInt();
			if (lengthOfOutOfRangeHashtable>0) {
				parameters = new byte[lengthOfOutOfRangeHashtable];
				file.read(parameters);
				outOfRangeValues = (Hashtable<String, Float>) Compression.bytesToObj(parameters);
			}

			file.close();

	        for (int i=0; i<indeciesOfMarkersToLoad.length; i++) {
				indexReadBuffer = 0;
	
				indexStart = 0;
				indexReadBuffer = indexStart;
				if (((nullStatus>>Sample.NULLSTATUS_GC_LOCATION) & 0x01) != 1) {
						gcs = new float[numSamples];
						for (int j=0; j<numSamples; j++) {
							gcs[j] = Compression.gcBafDecompress(new byte[] {readBuffer[i][indexReadBuffer], readBuffer[i][indexReadBuffer + 1]});
							indexReadBuffer += bytesPerSampMark;
						}
					indexStart += 2;
				}
				indexReadBuffer = indexStart;
				if (((nullStatus>>Sample.NULLSTATUS_X_LOCATION) & 0x01) != 1) {
						xs = new float[numSamples];
						for (int j=0; j<numSamples; j++) {
							xs[j] = Compression.xyDecompress(new byte[] {readBuffer[i][indexReadBuffer], readBuffer[i][indexReadBuffer + 1]});
							if (xs[j]==Compression.REDUCED_PRECISION_XY_OUT_OF_RANGE_FLOAT) {
	//							xs[j] = outOfRangeValues.get(sampleName+"\t"+allMarkersProj[j]+"\tx");
								xs[j] = outOfRangeValues.get(i + "\t" + j + "\tx");
							}
							indexReadBuffer += bytesPerSampMark;
						}
					indexStart += 2;
				}
				indexReadBuffer = indexStart;
				if (((nullStatus>>Sample.NULLSTATUS_Y_LOCATION) & 0x01) != 1) {
						ys = new float[numSamples];
						for (int j=0; j<numSamples; j++) {
							ys[j] = Compression.xyDecompress(new byte[] {readBuffer[i][indexReadBuffer], readBuffer[i][indexReadBuffer + 1]});
							if (ys[j]==Compression.REDUCED_PRECISION_XY_OUT_OF_RANGE_FLOAT) {
	//							ys[j] = outOfRangeValues.get(sampleName+"\t"+allMarkersProj[j]+"\ty");
								ys[j] = outOfRangeValues.get(i + "\t" + j + "\ty");
							}
							indexReadBuffer += bytesPerSampMark;
						}
	
					indexStart += 2;
				}
				indexReadBuffer = indexStart;
				if (((nullStatus>>Sample.NULLSTATUS_BAF_LOCATION) & 0x01) != 1) {
						bafs = new float[numSamples];
						for (int j=0; j<numSamples; j++) {
							bafs[j] = Compression.gcBafDecompress(new byte[] {readBuffer[i][indexReadBuffer], readBuffer[i][indexReadBuffer + 1]});
							indexReadBuffer += bytesPerSampMark;
						}
					indexStart += 2;
				}
				indexReadBuffer = indexStart;
				if (((nullStatus>>Sample.NULLSTATUS_LRR_LOCATION) & 0x01) != 1) {
						lrrs = new float[numSamples];
						for (int j=0; j<numSamples; j++) {
							lrrs[j] = Compression.lrrDecompress(new byte[] {readBuffer[i][indexReadBuffer], readBuffer[i][indexReadBuffer + 1], readBuffer[i][indexReadBuffer + 2]});
							if (lrrs[j] == Compression.REDUCED_PRECISION_LRR_OUT_OF_RANGE_LRR_FLOAT) {
	//							lrrs[j] = outOfRangeValues.get(sampleName+"\t"+allMarkersProj[j]+"\tlrr");
								lrrs[j] = outOfRangeValues.get(i + "\t" + j + "\tlrr");
							}
							indexReadBuffer += bytesPerSampMark;
						}
					indexStart += 3;
				}
				indexReadBuffer = indexStart;
				if (((nullStatus>>Sample.NULLSTATUS_ABGENOTYPE_LOCATION) & 0x01) != 1 || ((nullStatus>>Sample.NULLSTATUS_FOWARDGENOTYPE_LOCATION) & 0x01) != 1) {
					abGenotypes = new byte[numSamples];
					alleleMappings = new String[numSamples];
					for (int j=0; j<numSamples; j++) {
						genotypeTmp = Compression.genotypeDecompress(readBuffer[i][indexReadBuffer]);
						abGenotypes[j] = genotypeTmp[0];
						alleleMappings[j] = Sample.ALLELE_PAIRS[genotypeTmp[1]];
						indexReadBuffer += bytesPerSampMark;
					}
				}
		        result[i] = new MarkerData(markerNames[i], (byte)0, 0, fingerPrint, gcs, null, null, xs, ys, null, null, bafs, lrrs, abGenotypes, alleleMappings);
	        }
		} catch (FileNotFoundException e) {
			System.err.println("Error - could not find RAF marker file '"+markerFilename+"'");
			e.printStackTrace();
		} catch (IOException e) {
			System.err.println("Error reading RAF marker file '"+markerFilename+"'");
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
		}
		return result;
	}




	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		Project proj;
//		String filename = Project.DEFAULT_PROJECT;
		String filename = "C:/workspace/Genvisis/projects/gedi_exome.properties";
//		String filename = "C:/workspace/Genvisis/projects/practice.properties";
//		String demo = "flagged_results.txt";
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
				transposeData(proj, maxFileSize, true, null);
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
