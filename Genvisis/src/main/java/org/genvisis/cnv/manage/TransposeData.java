// -Xms1024M -Xmx1024M or even better: -Xmx15g
package org.genvisis.cnv.manage;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TimeZone;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import javax.swing.JOptionPane;

import org.genvisis.cnv.filesys.Compression;
import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.MarkerLookup;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.ext;

// TODO whole file needs clean up
public class TransposeData {
	public static final int DEFAULT_MARKERS_PER_FILE = 210000;
	public static final byte MARKERDATA_PARAMETER_TOTAL_LEN = 21;
	public static final byte MARKERDATA_NUMSAMPLES_START = 0;
	public static final byte MARKERDATA_NUMSAMPLES_LEN = 4;
	public static final byte MARKERDATA_NUMMARKERS_START = 4;
	public static final byte MARKERDATA_NUMMARKERS_LEN = 4;
	public static final byte MARKERDATA_NULLSTATUS_START = 8;
	public static final byte MARKERDATA_NULLSTATUS_LEN = 1;
	public static final byte MARKERDATA_FINGERPRINT_START = 9;
	public static final byte MARKERDATA_FINGERPRINT_LEN = 8;
	public static final byte MARKERDATA_MARKERNAMELEN_START = 17;
	public static final byte MARKERDATA_MARKERNAMELEN_LEN = 4;
	public static final byte MARKERDATA_MARKERNAME_START = 21;

	private static final String TEMP_SUCCESS_FILE = "success.temp";

	/**
	 * Transpose data from sample based structure into Marker based structure. The input is a set of
	 * files based on the sample oriented structure, and the output is a set of data files based on
	 * the marker oriented structure.
	 *
	 * @param proj The Genvisis project containing the data to be transposed
	 * @param markerFileSizeSuggested The maximum size in mb of a single output file. If set to 0, the
	 *        default value 2000 is taken on.
	 * @param keepAllSampleFilesOpen To keep all the files or just a single one open at a time.
	 * @param log The log file. If set as null, the default file
	 *        Genvisis_TransposeData_yyyyMMdd_HHmmss.log will be used.
	 */
	@SuppressWarnings("unchecked")
	public static void transposeData(Project proj, long markerFileSizeSuggested,
																	 boolean keepAllSampleFilesOpen) {
		boolean isFileClosed;
		boolean done; // safe;
		byte nullStatus = 0;
		byte numBytesPerSampleMarker;
		byte backupCount;
		int numBytes_Mark;
		int numMarkers_WriteBuffer;
		int numMarkers_Chunk;
		int numMarkers_File;
		int numChunks_WriteBuffer;
		int numChunks_File;
		int numFiles;
		int numMarkers_LastRound;
		int numRounds_LoadSampFile;
		int markerFileIndex;
		int numBufferChunksNeededCurrentMarkFile;
		int countTotalChunks_MarkerFile;
		int countTotalChunks_writeBuffer;
		int indexFirstMarkerCurrentIteration;
		int markerIndex;
		int numMarkersCurrentLoop;
		int index_WriteBufferEnd;
		int numChunks_RemainedInWriteBuffer;
		long fingerPrint;
		long timerOverAll, timerLoadFiles, timerTransposeMemory, timerWriteFiles, timerTmp;
		byte[] markFileParameterSection;
		byte[] readBuffer;
		byte[] markFileOutliersBytes = null;
		byte[][] writeBuffer = null;
		byte[][] markersInEachFile;
		int[] writeBufferSizes;
		int[] parameters;
		String[] allSampleNamesInProj;
		String[] allMarkerNamesInProj;
		String[] markersInEachFile1;
		String[] markerFilenames;
		Hashtable<String, String> markerLookup = new Hashtable<String, String>();
		Hashtable<String, Float>[] markFileOutliers;
		Hashtable<String, Float> allOutliers;
		SimpleDateFormat timeFormat;
		RandomAccessFile sampleFile;
		RandomAccessFile[] sampleFiles = null;
		Logger log;

		log = proj.getLog();
		log.report("Transposing data for the project in " + proj.PROJECT_DIRECTORY.getValue());

		timeFormat = new SimpleDateFormat("HH:mm:ss.SSS");
		timeFormat.setTimeZone(TimeZone.getTimeZone("GMT"));
		allSampleNamesInProj = proj.getSamples();
		allMarkerNamesInProj = proj.getMarkerNames();
		fingerPrint = org.genvisis.cnv.filesys.MarkerSet.fingerprint(allSampleNamesInProj);
		if (!Files.exists(proj.SAMPLE_DIRECTORY.getValue(true, true) + allSampleNamesInProj[0]
											+ Sample.SAMPLE_FILE_EXTENSION)) {
			log.reportError("Could not locate file: " + proj.SAMPLE_DIRECTORY.getValue(true, true)
											+ allSampleNamesInProj[0] + Sample.SAMPLE_FILE_EXTENSION
											+ "; aborting transposeData");
			return;
		}
		nullStatus = Sample.getNullstatusFromRandomAccessFile(proj.SAMPLE_DIRECTORY.getValue(true, true)
																													+ allSampleNamesInProj[0]
																													+ Sample.SAMPLE_FILE_EXTENSION);
		numBytesPerSampleMarker = Sample.getNBytesPerSampleMarker(nullStatus);
		numBytes_Mark = allSampleNamesInProj.length * numBytesPerSampleMarker;
		if (new File(proj.PROJECT_DIRECTORY.getValue()).getFreeSpace() <= (allSampleNamesInProj.length
																																			 * (long) allMarkerNamesInProj.length
																																			 * numBytesPerSampleMarker)) {
			log.reportError("Not enough disk space for all the new data to be created. Available: "
											+ ext.prettyUpSize(new File(proj.PROJECT_DIRECTORY.getValue()).getFreeSpace(),
																				 1)
											+ "; Required: "
											+ ext.prettyUpSize(new File(proj.PROJECT_DIRECTORY.getValue()).getFreeSpace(),
																				 1)
											+ ").");
			return;
		}
		if (markerFileSizeSuggested <= 0) {
			markerFileSizeSuggested = Integer.MAX_VALUE;
		}

		done = false;
		timerOverAll = new Date().getTime();
		numMarkers_WriteBuffer = Math.min(getOptimalNumSamplesBasedOnHeapSpace(-1, numBytes_Mark),
																			allMarkerNamesInProj.length);
		while (!done) {
			PSF.checkInterrupted();
			try {
				numMarkers_File = (int) Math.min((double) markerFileSizeSuggested / (double) numBytes_Mark,
																				 allMarkerNamesInProj.length);
				parameters = getOptimizedFileAndBufferSize(numMarkers_WriteBuffer,
																									 allMarkerNamesInProj.length, numBytes_Mark,
																									 numMarkers_File);
				numMarkers_File = parameters[0];
				numMarkers_WriteBuffer = parameters[1];
				numChunks_WriteBuffer = parameters[2];
				numMarkers_Chunk = parameters[3];
				numChunks_File = parameters[4];
				markerFileSizeSuggested = (long) numMarkers_File * (long) numBytes_Mark;
				numRounds_LoadSampFile = (int) Math.ceil((double) allMarkerNamesInProj.length
																								 / (double) numMarkers_WriteBuffer);
				numMarkers_LastRound = allMarkerNamesInProj.length % numMarkers_WriteBuffer;
				numFiles = (int) Math.ceil((double) allMarkerNamesInProj.length / (double) numMarkers_File);
				countTotalChunks_MarkerFile = (int) Math.ceil((double) allMarkerNamesInProj.length
																											/ (double) numMarkers_Chunk);
				countTotalChunks_writeBuffer = countTotalChunks_MarkerFile;
				markerFilenames = new String[numFiles];
				markersInEachFile = new byte[numFiles][];
				backupCount = backupOlderFiles(proj.MARKER_DATA_DIRECTORY.getValue(true, false),
																			 new String[] {MarkerData.MARKER_DATA_FILE_EXTENSION,
																										 "outliers.ser"},
																			 true);
				backupOlderFile(proj.MARKERLOOKUP_FILENAME.getValue(false, false), backupCount);
				log.report("--\nProject:\t" + ext.addCommas(allMarkerNamesInProj.length) + " markers\t"
									 + ext.addCommas(allSampleNamesInProj.length) + " samples" + "\nHeapSpace:\t"
									 + ext.prettyUpSize(Runtime.getRuntime().maxMemory(), 1) + " max"
									 + "\nwriteBuffer:\t" + ext.addCommas(numMarkers_WriteBuffer) + " markers\t"
									 + numChunks_WriteBuffer + " chunks\t" + ext.addCommas(numMarkers_Chunk)
									 + " Markers/chunk\t"
									 + ext.formDeci((double) numMarkers_WriteBuffer * numBytes_Mark
																	/ Runtime.getRuntime().maxMemory() * 100, 1)
									 + "% heap efficiency" + "\nMarkerFile:\t" + ext.addCommas(numMarkers_File)
									 + " markers\t" + numChunks_File + " chunks\t"
									 + markerFileSizeSuggested / 1024 / 1024 / 1024 + "."
									 + ((int) (markerFileSizeSuggested / 1024 / 1024 / 10.24)
											- (int) (markerFileSizeSuggested / 1024 / 1024 / 1024 * 102.4))
									 + " gb\t" + numFiles + " files");


				markerIndex = 0;
				for (int i = 0; i < numFiles; i++) {
					markerFilenames[i] = proj.MARKER_DATA_DIRECTORY.getValue(false, true) + "markers." + i
															 + MarkerData.MARKER_DATA_FILE_EXTENSION;
					numMarkersCurrentLoop = Math.min(numMarkers_File,
																					 allMarkerNamesInProj.length - markerIndex);
					markersInEachFile1 = new String[numMarkersCurrentLoop];
					for (int j = 0; j < numMarkersCurrentLoop; j++) {
						markerLookup.put(allMarkerNamesInProj[markerIndex],
														 "markers." + i + MarkerData.MARKER_DATA_FILE_EXTENSION + "\t" + j);
						markersInEachFile1[j] = allMarkerNamesInProj[markerIndex];
						markerIndex++;
					}
					markersInEachFile[i] = Compression.objToBytes(markersInEachFile1);
				}
				PSF.checkInterrupted();
				new MarkerLookup(markerLookup).serialize(proj.MARKERLOOKUP_FILENAME.getValue(false, false));


				PSF.checkInterrupted();
				indexFirstMarkerCurrentIteration = 0;
				readBuffer = new byte[numMarkers_WriteBuffer * numBytesPerSampleMarker];
				writeBuffer = new byte[numChunks_WriteBuffer][numMarkers_Chunk * numBytes_Mark];
				writeBufferSizes = new int[numChunks_WriteBuffer];
				for (int i = 0; i < writeBufferSizes.length; i++) {
					writeBufferSizes[i] = writeBuffer[i].length;
				}

				timerTmp = new Date().getTime();
				if (keepAllSampleFilesOpen) {
					sampleFiles = new RandomAccessFile[allSampleNamesInProj.length];
					for (int i = 0; i < allSampleNamesInProj.length; i++) {
						PSF.checkInterrupted();
						try {
							sampleFiles[i] = new RandomAccessFile(proj.SAMPLE_DIRECTORY.getValue(true, true)
																										+ allSampleNamesInProj[i]
																										+ Sample.SAMPLE_FILE_EXTENSION, "r");
						} catch (FileNotFoundException fnfe) {
							log.reportError("Error - file not found: "
															+ proj.SAMPLE_DIRECTORY.getValue(true, true)
															+ allSampleNamesInProj[i]
															+ Sample.SAMPLE_FILE_EXTENSION);
							log.reportError("        if you get this error and the file does exist, then likely your operating system (especially common on a linux platform) is not allowing this many files to be open at once; rerun without using that option at the command line");
							log.reportError("        Transpose aborted");
							return;
						}
					}
				}

				if (new File(proj.SAMPLE_DIRECTORY.getValue(true, true) + "outliers.ser").exists()) {
					allOutliers = (Hashtable<String, Float>) SerializedFiles.readSerial(proj.SAMPLE_DIRECTORY.getValue(true,
																																																						 true)
																																							+ "outliers.ser");
				} else {
					allOutliers = new Hashtable<String, Float>();
				}
				markFileOutliers = getOutlierHashForEachMdRafFile(allOutliers, numFiles, numMarkers_File,
																													allSampleNamesInProj);

				PSF.checkInterrupted();

				markerFileIndex = 0;
				numBufferChunksNeededCurrentMarkFile = 0;
				isFileClosed = true;

				timerWriteFiles = 0;
				log.report("--\ni (<" + numRounds_LoadSampFile + ")\tLoad\tTranspose\tWrite");
				for (int i = 0; i < numRounds_LoadSampFile; i++) {

					PSF.checkInterrupted();
					if ((i + 1) == numRounds_LoadSampFile && numMarkers_LastRound != 0) {
						numChunks_WriteBuffer = (int) Math.ceil((double) numMarkers_LastRound
																										/ (double) numMarkers_Chunk);
						for (int j = numChunks_WriteBuffer; j < writeBuffer.length; j++) {
							writeBufferSizes[j] = 0;
						}
						writeBufferSizes[numChunks_WriteBuffer - 1] = (numMarkers_LastRound % numMarkers_Chunk)
																													* numBytes_Mark;
					}

					// --- Step 1 --- Load sample files into buffer
					timerLoadFiles = 0;
					timerTransposeMemory = 0;
					for (int j = 0; j < allSampleNamesInProj.length; j++) {

						PSF.checkInterrupted();
						timerTmp = new Date().getTime();
						if (!keepAllSampleFilesOpen) {
							sampleFile = new RandomAccessFile(proj.SAMPLE_DIRECTORY.getValue(true, true)
																								+ allSampleNamesInProj[j]
																								+ Sample.SAMPLE_FILE_EXTENSION, "r");
						} else {
							sampleFile = sampleFiles[j];
						}

						Sample.loadFromRandomAccessFileWithoutDecompress(sampleFile, readBuffer, true, j,
																														 indexFirstMarkerCurrentIteration,
																														 numBytesPerSampleMarker,
																														 allMarkerNamesInProj.length, null,
																														 log);

						PSF.checkInterrupted();

						timerLoadFiles += (new Date().getTime() - timerTmp);
						timerTmp = new Date().getTime();

						transposeBuffer(writeBuffer, writeBufferSizes, readBuffer, numBytesPerSampleMarker,
														indexFirstMarkerCurrentIteration, j, allSampleNamesInProj.length);

						PSF.checkInterrupted();

						timerTransposeMemory += (new Date().getTime() - timerTmp);

						if (!keepAllSampleFilesOpen) {
							sampleFile.close();
						}
					}
					log.report(i + "\t" + timeFormat.format(timerLoadFiles) + "\t"
										 + timeFormat.format(timerTransposeMemory), false, true);


					// --- Step 2 --- Dump write buffer to marker files
					numChunks_RemainedInWriteBuffer = Math.min(countTotalChunks_writeBuffer,
																										 numChunks_WriteBuffer);
					countTotalChunks_writeBuffer -= numChunks_RemainedInWriteBuffer;
					for (int j = 0; j < numChunks_RemainedInWriteBuffer; j++) {

						PSF.checkInterrupted();
						if (isFileClosed) {
							numBufferChunksNeededCurrentMarkFile = Math.min(numChunks_File,
																															countTotalChunks_MarkerFile);
							countTotalChunks_MarkerFile -= numBufferChunksNeededCurrentMarkFile;
							if (markerFileIndex == (numFiles - 1)
									&& allMarkerNamesInProj.length % numMarkers_File > 0) {
								numMarkers_File = allMarkerNamesInProj.length % numMarkers_File;
							}
							if (i == 1) {
								System.out.println("");
							}
							markFileParameterSection = getParameterSectionForMdRaf(allSampleNamesInProj.length,
																																		 numMarkers_File,
																																		 nullStatus,
																																		 fingerPrint,
																																		 markersInEachFile[markerFileIndex]);
							timerWriteFiles = 0;
						} else {
							markFileParameterSection = null;
						}

						if ((j + numBufferChunksNeededCurrentMarkFile) > numChunks_RemainedInWriteBuffer) {
							index_WriteBufferEnd = numChunks_RemainedInWriteBuffer - 1;
							markFileOutliersBytes = null;
							isFileClosed = false;
							numBufferChunksNeededCurrentMarkFile -= (numChunks_RemainedInWriteBuffer - j);
						} else {
							index_WriteBufferEnd = j + numBufferChunksNeededCurrentMarkFile - 1;
							if (markFileOutliers == null || markFileOutliers[markerFileIndex].size() == 0) {
								markFileOutliersBytes = new byte[0];
							} else {
								markFileOutliersBytes = Compression.objToBytes(markFileOutliers[markerFileIndex]);
							}
							isFileClosed = true;
						}

						timerTmp = new Date().getTime();
						writeBufferToRAF(writeBuffer, writeBufferSizes, j, index_WriteBufferEnd,
														 markerFilenames[markerFileIndex], markFileParameterSection,
														 markFileOutliersBytes, true);

						PSF.checkInterrupted();

						j = index_WriteBufferEnd;
						if (isFileClosed) {
							markerFileIndex++;
						}
						timerWriteFiles += (new Date().getTime() - timerTmp);
						log.report("\t" + timeFormat.format(timerWriteFiles), false, true);

					}

					indexFirstMarkerCurrentIteration += numMarkers_WriteBuffer;
					log.report("");
				}

				if (allOutliers != null && allOutliers.size() != 0) {
					SerializedFiles.writeSerial(allOutliers, proj.MARKER_DATA_DIRECTORY.getValue(false, true)
																									 + "outliers.ser");
				}

				done = true;

			} catch (FileNotFoundException e) {
				System.err.println(ext.getTime() + "\tFileNotFoundException");
				e.printStackTrace();
			} catch (IOException e) {
				System.err.println(ext.getTime() + "\tIOException");
				e.printStackTrace();
			} catch (OutOfMemoryError oome) {
				numMarkers_WriteBuffer = getOptimalNumSamplesBasedOnHeapSpace(numMarkers_WriteBuffer, -1);
				deleteOlderRafs(proj.MARKER_DATA_DIRECTORY.getValue(true, false), null,
												new String[] {MarkerData.MARKER_DATA_FILE_EXTENSION, "outliers.ser"},
												false,
												null);
			} catch (Exception e) {
				e.printStackTrace();
			}

		}
		log.report("--\nFinished transposing data. Total Time used: "
							 + ext.getTimeElapsed(timerOverAll)
							 + "\n");
	}


	public static int[] getOptimizedFileAndBufferSize(int numMarkers_WriteBuffer,
																										int numMarkers_Proj,
																										int numBytes_Mark, int numMarkers_File) {
		int numMarkers1_File;
		int numMarkers1_WriteBuffer;
		int numMarkers1_Chunk;
		int numChunks1_WriteBuffer;
		int numChunks1_File;

		int numRounds_SampFile;
		int numRounds_MarkFile;
		int numMarkers_WriteBuffer_Min;
		int numChunks_WriteBuffer_Max;
		int numChunks_WriteBuffer_Min;
		int numMarkers_Chunk_Max;
		int numMarkers_Chunk_Min;

		if ((Integer.MAX_VALUE / numBytes_Mark) < numMarkers_WriteBuffer) {
			numRounds_SampFile = (int) Math.ceil((double) numMarkers_Proj / numMarkers_WriteBuffer);
			numMarkers_WriteBuffer_Min = (int) Math.ceil((double) numMarkers_Proj / numRounds_SampFile);

			numMarkers_Chunk_Max = Integer.MAX_VALUE / numBytes_Mark;
			numRounds_MarkFile = (int) Math.ceil((double) numMarkers_Proj / numMarkers_Chunk_Max);
			numMarkers_Chunk_Min = numMarkers_Proj / numRounds_MarkFile;

			numChunks_WriteBuffer_Max = numMarkers_WriteBuffer / numMarkers_Chunk_Min;
			numChunks_WriteBuffer_Min = numMarkers_WriteBuffer_Min / numMarkers_Chunk_Max;

			while (numChunks_WriteBuffer_Max == numChunks_WriteBuffer_Min
						 && numMarkers_WriteBuffer % numMarkers_Chunk_Min != numMarkers_WriteBuffer_Min
																																 % numMarkers_Chunk_Max) {
				numRounds_MarkFile++;
				numMarkers_Chunk_Min = numMarkers_Proj / numRounds_MarkFile;
				numChunks_WriteBuffer_Max = numMarkers_WriteBuffer / numMarkers_Chunk_Min;
			}

			numChunks1_WriteBuffer = numChunks_WriteBuffer_Max;
			if ((numChunks1_WriteBuffer * numMarkers_Chunk_Min) >= numMarkers_WriteBuffer_Min) {
				numMarkers1_Chunk = numMarkers_Chunk_Min;
			} else {
				numMarkers1_Chunk = numMarkers_WriteBuffer_Min / numChunks1_WriteBuffer;
			}

		} else {
			numChunks1_WriteBuffer = 1;
			numMarkers1_Chunk = numMarkers_WriteBuffer;
		}

		numMarkers1_WriteBuffer = numChunks1_WriteBuffer * numMarkers1_Chunk;
		numChunks1_File = (int) Math.round((double) numMarkers_File / numMarkers1_Chunk);
		numMarkers1_File = numChunks1_File * numMarkers1_Chunk;

		return new int[] {numMarkers1_File, numMarkers1_WriteBuffer, numChunks1_WriteBuffer,
											numMarkers1_Chunk, numChunks1_File};
	}


	/**
	 * Rebuild MarkerLookup from transposed data files (.mkRAF).
	 *
	 * @param proj
	 */
	public static void recreateMarkerLookup(Project proj) {
		Hashtable<String, String> hash = new Hashtable<String, String>();
		String[] files, markerNames;
		byte[] readBuffer;
		RandomAccessFile currentFile;
		long time;

		time = new Date().getTime();
		System.out.println("Creating MarkerLookup file");
		files = new File(proj.MARKER_DATA_DIRECTORY.getValue(false, true)).list(new FilenameFilter() {
			@Override
			public boolean accept(File file, String filename) {
				return filename.endsWith(MarkerData.MARKER_DATA_FILE_EXTENSION)
							 && filename.startsWith("marker");
			}
		});
		if (files == null) {
			System.err.println("Error - failed to create MarkerLookup -- marker data directory does not exist: "
												 + proj.MARKER_DATA_DIRECTORY.getValue(false, true));
			System.err.println("      - Did you forget to transpose the data?");
		} else if (files.length == 0) {
			System.err.println("Error - failed to create MarkerLookup -- no "
												 + MarkerData.MARKER_DATA_FILE_EXTENSION + " files available");
		} else {
			for (int i = 0; i < files.length; i++) {
				try {
					proj.getLog().report((i + 1) + " of " + files.length);
					currentFile = new RandomAccessFile(proj.MARKER_DATA_DIRECTORY.getValue(false, true)
																						 + files[i], "r");
					currentFile.seek(MARKERDATA_MARKERNAMELEN_START);
					readBuffer = new byte[currentFile.readInt()];
					currentFile.read(readBuffer);
					markerNames = (String[]) Compression.bytesToObj(readBuffer);
					for (int j = 0; j < markerNames.length; j++) {
						hash.put(markerNames[j], files[i] + "\t" + j);
					}
				} catch (FileNotFoundException e) {
					e.printStackTrace();
				} catch (Exception e) {
					e.printStackTrace();
				}
			}

			new MarkerLookup(hash).serialize(proj.MARKERLOOKUP_FILENAME.getValue(true, false));

			System.out.println("Created MarkerLookup in " + ext.getTimeElapsed(time));
		}
	}


	private static void transposeBuffer(byte[][] writeBuffer, int[] writeBufferSizes,
																			byte[] readBuffer, byte bytesPerSampleMarker,
																			int indexOfFirstMarkerInBuffer, int indexOfCurrentSample,
																			int numSamplesInProj) {
		int numMarkersInChunk;
		int indexInReadBuffer;
		int indexInChunk;
		int step;

		indexInReadBuffer = 0;
		step = numSamplesInProj * bytesPerSampleMarker;
		for (int i = 0; i < writeBuffer.length; i++) {
			numMarkersInChunk = writeBufferSizes[i] / bytesPerSampleMarker / numSamplesInProj;
			indexInChunk = indexOfCurrentSample * bytesPerSampleMarker;

			for (int j = 0; j < numMarkersInChunk; j++) {
				for (int k = 0; k < bytesPerSampleMarker; k++) {
					try {
						writeBuffer[i][indexInChunk + k] = readBuffer[indexInReadBuffer];
					} catch (ArrayIndexOutOfBoundsException e) {
						System.out.println("j:" + j + "\tnumMarkersInChunk:" + numMarkersInChunk);
						System.out.println("i:" + i + "\tj:" + j + "\tk:" + k + "\tindexInChunk: "
															 + indexInChunk);
						System.out.println("writebuffer size: " + writeBuffer.length + " , "
															 + writeBuffer[i].length + "\treadBuffer size: " + readBuffer.length);
						e.printStackTrace();
						System.exit(1);
					}
					indexInReadBuffer++;
				}
				indexInChunk += step;
			}
		}
	}

	public static void reverseTransposeStreaming(Project proj) {
		long fingerprintForMarkers, fingerprintForSamples;

		Hashtable<String, Float>[] sampRafFileOutliers;
		Hashtable<String, Float> allOutliers;
		SimpleDateFormat timeFormat;
		String filename;
		Logger log;
		String successFile;
		Set<String> successes;
		// ---
		String smp;

		log = proj.getLog();
		log.report("Reverse transposing data for the project in " + proj.PROJECT_DIRECTORY.getValue());

		timeFormat = new SimpleDateFormat("HH:mm:ss.SSS");
		timeFormat.setTimeZone(TimeZone.getTimeZone("GMT"));

		successFile = proj.MARKER_DATA_DIRECTORY.getValue() + TEMP_SUCCESS_FILE;
		successes = Files.exists(successFile) ? HashVec.loadFileToHashSet(successFile, false)
																					: new HashSet<>();
		if (successes.size() > 0) {
			log.reportTime("Skipping " + successes.size() + " already existing samples.");
		}

		final String[] listOfAllSamplesInProj = proj.getSamples();
		final String[] listOfAllMarkersInProj = proj.getMarkerNames();

		@SuppressWarnings("deprecation")
		long f = org.genvisis.cnv.filesys.MarkerSet.fingerprint(listOfAllSamplesInProj);
		fingerprintForSamples = f;
		fingerprintForMarkers = proj.getMarkerSet().getFingerprint();

		// compute null status (see GenvisisTechDoc.docx for explanation)
		final byte nullStatus = getNullstatusFromRandomAccessFile(proj.MARKER_DATA_DIRECTORY.getValue(false,
																																																	true)
																															+ proj.getMarkerLookup()
																																		.getFirstMarkerDataRafFilename(),
																															false);
		int numBytesPerSampleMarker = Sample.getNBytesPerSampleMarker(nullStatus);
		// numBytesPerSamp = listOfAllMarkersInProj.length * numBytesPerSampleMarker;

		allOutliers = loadMarkerDataOutliers(proj);
		sampRafFileOutliers = getOutlierHashForEachSampleRafFile(allOutliers,
																														 listOfAllSamplesInProj);

		final MarkerLookup mkrLookup = proj.getMarkerLookup();

		// initialize all sample files:
		boolean[] load = new boolean[listOfAllSamplesInProj.length];
		for (int s = 0; s < listOfAllSamplesInProj.length; s++) {
			smp = listOfAllSamplesInProj[s];
			filename = proj.SAMPLE_DIRECTORY.getValue(true, true) + smp + Sample.SAMPLE_FILE_EXTENSION;
			load[s] = !successes.contains(filename);
		}

		// TODO tweak these parameters
		int maxThreads = 20;
		maxThreads = Math.min(Runtime.getRuntime().availableProcessors(), maxThreads);
		int threadBuffer = (int) ((Runtime.getRuntime().maxMemory() * 0.75) / maxThreads);

		int numSamplesPerThread = (listOfAllSamplesInProj.length / (2 * maxThreads));

		log.reportTime("Reverse-transposing " + (listOfAllSamplesInProj.length - successes.size())
									 + " samples (" + successes.size() + " already parsed) @ " + numSamplesPerThread
									 + " samples per thread with "
									 + maxThreads + " threads.");

		ExecutorService transposeService = Executors.newFixedThreadPool(maxThreads);

		PrintWriter successWriter = Files.getAppropriateWriter(successFile, Files.exists(successFile));
		for (int s = 0; s < listOfAllSamplesInProj.length; s += numSamplesPerThread) {
			final int s1 = s;
			final int s2 = Math.min(listOfAllSamplesInProj.length, s1 + numSamplesPerThread);
			Runnable runner = new Runnable() {
				@Override
				public void run() {
					int numMarkersCurrentFile;
					String smp;
					String mkr;
					String markerFile;
					String filename;
					String markerLookupEntry;
					byte[] sampOutlrs;
					RandomAccessFile markerRAF;
					byte[] parameterReadBuffer;
					long time;
					long time2;
					String openMarkerFile;
					try {
						Map<Integer, String> sampFiles = new HashMap<>();
						time = System.nanoTime();
						for (int i = s1; i < s2; i++) {
							smp = listOfAllSamplesInProj[i];
							filename = proj.SAMPLE_DIRECTORY.getValue(true, true) + smp
												 + Sample.SAMPLE_FILE_EXTENSION;
							if (!load[i]) {
								continue;
							}
							sampFiles.put(i, filename);

							FileOutputStream fos = new FileOutputStream(filename, false);
							// write header
							if (sampRafFileOutliers == null || sampRafFileOutliers[i].size() == 0) {
								sampOutlrs = new byte[0];
							} else {
								sampOutlrs = Compression.objToBytes(sampRafFileOutliers[i]);
							}
							fos.write(getParameterSectionForSampRaf(listOfAllMarkersInProj.length,
																											nullStatus,
																											sampOutlrs.length,
																											fingerprintForMarkers));
							fos.flush();
							fos.close();
							sampOutlrs = null;
						}
						log.reportTime("Wrote headers for samples " + s1 + " to " + s2 + " in "
													 + ext.getTimeElapsedNanos(time));

						// read markers
						parameterReadBuffer = new byte[TransposeData.MARKERDATA_PARAMETER_TOTAL_LEN];
						openMarkerFile = mkrLookup.get(listOfAllMarkersInProj[0])
																			.split(PSF.Regex.GREEDY_WHITESPACE)[0];
						// open first marker file
						markerRAF = new RandomAccessFile(
																						 proj.MARKER_DATA_DIRECTORY.getValue() + openMarkerFile,
																						 "r");
						markerRAF.read(parameterReadBuffer);

						time2 = System.nanoTime();
						// iterate through all markers
						for (int m = 0; m < listOfAllMarkersInProj.length; m += 0) {
							mkr = listOfAllMarkersInProj[m];
							markerLookupEntry = mkrLookup.get(mkr);
							markerFile = markerLookupEntry.split(PSF.Regex.GREEDY_WHITESPACE)[0];
							if (!openMarkerFile.equals(markerFile)) {
								// close current file and open next
								markerRAF.close();
								log.reportTime("Parsed " + openMarkerFile + " in "
															 + ext.getTimeElapsedNanos(time2));
								openMarkerFile = markerFile;
								time2 = System.nanoTime();
								markerRAF = new RandomAccessFile(
																								 proj.MARKER_DATA_DIRECTORY.getValue()
																								 + openMarkerFile,
																								 "r");
								markerRAF.read(parameterReadBuffer);
							}
							numMarkersCurrentFile = Compression.bytesToInt(parameterReadBuffer,
																														 TransposeData.MARKERDATA_NUMMARKERS_START);
							int numMkrsBuffer = Math.min(numMarkersCurrentFile,
																					 threadBuffer / (sampFiles.size()
																													 * numBytesPerSampleMarker));

							MarkerDataLoader.loadAndWriteFromMarkerDataRafWithoutDecompressRange(markerRAF,
																																										sampFiles,
																																										parameterReadBuffer,
																																										ArrayUtils.arrayOfIndices(numMarkersCurrentFile),
																																										numMkrsBuffer,
																																										fingerprintForSamples,
																																										log);
							m += numMarkersCurrentFile;
						}
						markerRAF.close();
						for (Entry<Integer, String> sampFile : sampFiles.entrySet()) {
							synchronized (successWriter) {
								successWriter.println(sampFile.getValue());
								successWriter.flush();
							}

							log.reportTime("Parsed markerData for sample " + sampFile.getKey().intValue()
														 + " in "
														 + ext.getTimeElapsedNanos(time));
						}
					} catch (FileNotFoundException e) {
						log.reportException(e);
					} catch (IOException e) {
						throw new RuntimeException(e);
					}
				}
			};
			transposeService.submit(runner);
		}

		transposeService.shutdown();
		try {
			transposeService.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
		} catch (InterruptedException e) {
			/**/
		}

		successWriter.close();
	}

	public static void reverseTranspose(Project proj) {
		boolean done;
		byte nullStatus = 0;
		byte numBytesPerSampleMarker;
		int numBytes_PerSamp;
		int numSamples_WriteBuffer;
		int indexFirstMarkerCurrentIteration;
		long fingerprintForMarkers, fingerprintForSamples;
		long timerOverAll, timerLoadFiles, timerTransposeMemory, timerWrite, timerTmp;
		byte[] markFileParameterSection;
		byte[][] readBuffer;
		byte[] markFileOutliersBytes = null;
		byte[][] writeBuffer = null;
		String logTemp;
		String[] listOfAllSamplesInProj;
		String[] listOfAllMarkersInProj;
		Hashtable<String, Float>[] sampRafFileOutliers;
		Hashtable<String, Float> allOutliers;
		String[] markerDataRafFilenames;
		SimpleDateFormat timeFormat;
		String markerFile;
		String filename;
		Logger log;
		Set<String> successes;
		int[][] batchIndices;
		int[] batch;
		boolean range;

		log = proj.getLog();
		log.report("Reverse transposing data for the project in " + proj.PROJECT_DIRECTORY.getValue());

		timeFormat = new SimpleDateFormat("HH:mm:ss.SSS");
		timeFormat.setTimeZone(TimeZone.getTimeZone("GMT"));
		listOfAllSamplesInProj = proj.getSamples();
		listOfAllMarkersInProj = proj.getMarkerNames();
		@SuppressWarnings("deprecation")
		long f = org.genvisis.cnv.filesys.MarkerSet.fingerprint(listOfAllSamplesInProj);
		fingerprintForSamples = f;
		fingerprintForMarkers = proj.getMarkerSet().getFingerprint();
		nullStatus = getNullstatusFromRandomAccessFile(proj.MARKER_DATA_DIRECTORY.getValue(false, true)
																									 + proj.getMarkerLookup()
																												 .getFirstMarkerDataRafFilename(),
																									 false);
		numBytesPerSampleMarker = Sample.getNBytesPerSampleMarker(nullStatus);
		if (!checkFreeSpaceOrLog(proj, listOfAllSamplesInProj.length, listOfAllMarkersInProj.length,
														 numBytesPerSampleMarker)) {
			return;
		}
		numBytes_PerSamp = listOfAllMarkersInProj.length * numBytesPerSampleMarker;

		done = false;
		timerOverAll = new Date().getTime();
		numSamples_WriteBuffer = Math.min(getOptimalNumSamplesBasedOnHeapSpace(-1, numBytes_PerSamp),
																			listOfAllSamplesInProj.length);

		String successFile = proj.MARKER_DATA_DIRECTORY.getValue() + TEMP_SUCCESS_FILE;
		successes = Files.exists(successFile) ? HashVec.loadFileToHashSet(successFile, false)
																					: new HashSet<>();
		PrintWriter writer = Files.getAppropriateWriter(successFile, true);

		while (!done) { // used to recover from OutOfMemory ONLY
			try {
				batchIndices = ArrayUtils.splitUpIntoBinsOfIndices(listOfAllSamplesInProj, successes,
																													 numSamples_WriteBuffer, log);

				writeStatMsg(log, listOfAllSamplesInProj.length, listOfAllMarkersInProj.length,
										 numSamples_WriteBuffer, numBytes_PerSamp);

				writeBuffer = new byte[numSamples_WriteBuffer][numBytes_PerSamp];

				timerTmp = new Date().getTime();
				allOutliers = loadMarkerDataOutliers(proj);
				sampRafFileOutliers = getOutlierHashForEachSampleRafFile(allOutliers,
																																 listOfAllSamplesInProj);

				markerDataRafFilenames = proj.getMarkerLookup().getMarkerDataRafFilenames();
				log.reportTime("--");
				log.reportTime("i (<" + batchIndices.length + ")\tLoad\tTranspose\tWrite");
				for (int i = 0; i < batchIndices.length; i++) {
					logTemp = "";
					batch = batchIndices[i];
					int batchRange = batch[batch.length - 1] - batch[0] + 1;
					range = batchRange == numSamples_WriteBuffer || batchRange == batch.length;

					// --- Step 1 --- Load mdRaf files into buffer, and transpose the buffer
					timerLoadFiles = 0;
					timerTransposeMemory = 0;
					timerWrite = 0;
					indexFirstMarkerCurrentIteration = 0;
					for (String markerDataRafFilename : markerDataRafFilenames) {
						timerTmp = new Date().getTime();
						markerFile = proj.MARKER_DATA_DIRECTORY.getValue(true, true) + markerDataRafFilename;
						if (range) {
							readBuffer = MarkerDataLoader.loadFromMarkerDataRafWithoutDecompressRange(markerFile,
																																												null,
																																												batch[0],
																																												batch.length,
																																												fingerprintForSamples,
																																												log);
						} else {
							readBuffer = MarkerDataLoader.loadFromMarkerDataRafWithoutDecompressIndices(markerFile,
																																													null,
																																													batch,
																																													fingerprintForSamples,
																																													log);
						}
						timerLoadFiles += (new Date().getTime() - timerTmp);

						timerTmp = new Date().getTime();
						reverseTransposeBuffer(readBuffer, writeBuffer, numBytesPerSampleMarker,
																	 indexFirstMarkerCurrentIteration);
						timerTransposeMemory += (new Date().getTime() - timerTmp);
						indexFirstMarkerCurrentIteration += readBuffer.length;
					}
					logTemp += (i + "\t" + timeFormat.format(timerLoadFiles) + "\t"
											+ timeFormat.format(timerTransposeMemory));

					// --- Step 2 --- Dump write buffer to marker files
					for (int j = 0; j < batch.length; j++) {
						timerTmp = new Date().getTime();
						if (sampRafFileOutliers == null || sampRafFileOutliers[batch[j]].size() == 0) {
							markFileOutliersBytes = new byte[0];
						} else {
							markFileOutliersBytes = Compression.objToBytes(sampRafFileOutliers[batch[j]]);
						}
						markFileParameterSection = getParameterSectionForSampRaf(listOfAllMarkersInProj.length,
																																		 nullStatus,
																																		 markFileOutliersBytes.length,
																																		 fingerprintForMarkers);

						filename = proj.SAMPLE_DIRECTORY.getValue(true, true)
											 + listOfAllSamplesInProj[batch[j]]
											 + Sample.SAMPLE_FILE_EXTENSION;
						boolean write = true;
						if (Files.exists(filename)) {
							int lenByt = Files.getSize(filename);
							if (lenByt == numBytes_PerSamp + Sample.PARAMETER_SECTION_BYTES) {
								write = false;
							}
						}
						if (write) {
							writeBufferToRAF(writeBuffer, null, j, j, filename,
															 markFileParameterSection, markFileOutliersBytes, false);
							writer.println(filename);
							timerWrite += (new Date().getTime() - timerTmp);
						}
					}
					logTemp += "\t" + timeFormat.format(timerWrite);

					log.reportTime(logTemp);
				}

				if (allOutliers != null && allOutliers.size() != 0) {
					SerializedFiles.writeSerial(allOutliers,
																			proj.SAMPLE_DIRECTORY.getValue(false, true) + "outliers.ser");
				}

				done = true;
				writer.close();

			} catch (FileNotFoundException e) {
				System.err.println(ext.getTime() + "\tFileNotFoundException");
				e.printStackTrace();
				proj.message("Major error - file not found:\n" + e.getMessage()
										 + "\nPlease check log for more info.", "Error - File Missing!",
										 JOptionPane.ERROR_MESSAGE);
				System.exit(1);
			} catch (IOException e) {
				System.err.println(ext.getTime() + "\tIOException");
				e.printStackTrace();
				proj.message("Major error:\n" + e.getMessage() + "\nPlease check log for more info.",
										 "Error!", JOptionPane.ERROR_MESSAGE);
				System.exit(1);
			} catch (OutOfMemoryError oome) {
				log.reportTimeWarning("Ran out of memory!  Trying again...");
			} catch (Exception e) {
				e.printStackTrace();
				proj.message("Major error:\n" + e.getMessage() + "\nPlease check log for more info.",
										 "Error!", JOptionPane.ERROR_MESSAGE);
				System.exit(1);
			}

		}
		timerOverAll = (new Date().getTime() - timerOverAll);
		log.report("--\nFinished reversely transposing data. Total Time used: "
							 + timeFormat.format(timerOverAll));
		new File(successFile).delete();
	}

	@SuppressWarnings("unchecked")
	private static Hashtable<String, Float> loadMarkerDataOutliers(Project proj) {
		if (new File(proj.MARKER_DATA_DIRECTORY.getValue(true, true) + "outliers.ser").exists()) {
			return (Hashtable<String, Float>) SerializedFiles.readSerial(proj.MARKER_DATA_DIRECTORY.getValue(true,
																																																			 true)
																																	 + "outliers.ser");
		} else {
			return new Hashtable<String, Float>();
		}
	}


	private static void writeStatMsg(Logger log, int numSamp, int numMark, long sampBuffer,
																	 long numBytes_PerSamp) {
		log.reportTime("--");
		log.reportTime("Project:\t" + numMark + " markers\t" + numSamp + " samples");
		log.reportTime("HeapSpace:\t" + ext.prettyUpSize(Runtime.getRuntime().maxMemory(), 1) + " max");
		log.reportTime("writeBuffer:\t" + sampBuffer + " samples\t"
									 + ext.formDeci((double) sampBuffer * numBytes_PerSamp
																	/ Runtime.getRuntime().maxMemory() * 100, 1)
									 + "% heap efficiency");
		log.reportTime("Bytes per Sample:\t" + numBytes_PerSamp);
	}

	private static boolean checkFreeSpaceOrLog(Project proj, int numSamp, int numMark,
																						 int numBytesPerSampleMarker) {
		Logger log = proj.getLog();
		long freeSpace = new File(proj.PROJECT_DIRECTORY.getValue()).getFreeSpace();
		long reqrSpace = (numSamp * numMark * numBytesPerSampleMarker);
		if (freeSpace <= reqrSpace) {
			log.reportError("Not enough disk space for all the new data to be created. Available: "
											+ ext.prettyUpSize(freeSpace, 1)
											+ "; Required: "
											+ ext.prettyUpSize(reqrSpace, 1)
											+ ").");
			return false;
		}
		return true;
	}


	private static void reverseTransposeBuffer(byte[][] input_markerLeadBuffer,
																						 byte[][] output_sampleLeadBuffer,
																						 byte bytesPerSampleMarker,
																						 int indexOfFirstMarkerInBuffer) {
		int numSamples;
		// int numMarkers;
		int indexSample;
		int indexMarker;
		int indexByte;

		numSamples = input_markerLeadBuffer[0].length / bytesPerSampleMarker;
		indexMarker = indexOfFirstMarkerInBuffer * bytesPerSampleMarker;
		for (int i = 0; i < input_markerLeadBuffer.length; i++) {
			indexSample = 0;
			indexByte = 0;
			for (int j = 0; j < input_markerLeadBuffer[i].length; j++) {
				if (indexByte == bytesPerSampleMarker) {
					indexByte = 0;
					indexSample++;
				}
				try {
					output_sampleLeadBuffer[indexSample][indexMarker
																							 + indexByte] = input_markerLeadBuffer[i][j];
				} catch (ArrayIndexOutOfBoundsException e) {
					System.out.println("j:" + indexByte + "\tnumMarkersInChunk:" + numSamples);
					System.out.println("i:" + i + "\tj:" + indexByte + "\tk:" + indexByte
														 + "\tindexInChunk: "
														 + indexMarker);
					System.out.println("writebuffer size: " + output_sampleLeadBuffer.length + " , "
														 + output_sampleLeadBuffer[i].length + "\treadBuffer size: "
														 + input_markerLeadBuffer.length);
					e.printStackTrace();
					System.exit(1);
				}
				indexByte++;
			}
			indexMarker += bytesPerSampleMarker;
		}
	}

	public static void writeBufferToRAF(byte[][] buffer, int[] bufferLength, int start, int end,
																			String fileName, byte[] head, byte[] tail,
																			boolean writeTailLength) {
		BufferedOutputStream markerFile;
		try {
			markerFile = new BufferedOutputStream(new FileOutputStream(fileName,
																																 head == null ? true : false));
			writeBufferToRAF(buffer, bufferLength, start, end, markerFile, head, tail, writeTailLength);
			markerFile.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

		// ObjectOutputStream markerFile;
		// try {
		// markerFile = new ObjectOutputStream(new FileOutputStream(fileName, head==null? true :
		// false));
		// writeBufferToRAF(buffer, bufferLength, start, end, markerFile, head, tail);
		// markerFile.close();
		// } catch (IOException e) {
		// e.printStackTrace();
		// }


		// RandomAccessFile markerFile;
		// try {
		// markerFile = new RandomAccessFile(fileName, "rw");
		// if (head != null && head.length != 0) {
		// markerFile.seek(markerFile.length());
		// }
		// writeBufferToRAF(buffer, bufferLength, start, end, markerFile, head, tail);
		// markerFile.close();
		// } catch (IOException e) {
		// e.printStackTrace();
		// }
	}


	public static void writeBufferToRAF(byte[][] buffer, int[] bufferLength, int indexOfStart,
																			int indexOfEnd, BufferedOutputStream markerFile, byte[] head,
																			byte[] tail, boolean writeTailLength) {
		if (buffer == null || indexOfStart < 0 || indexOfEnd >= buffer.length
				|| indexOfEnd < indexOfStart) {
			System.err.println("\nTranspose Data encoutered the following error: buffer is null, or start index of buffer is negative, or end index is less than the start index, or end index is over the buffer size.");
			System.exit(1);
		}

		try {
			if (head != null && head.length != 0) {
				markerFile.write(head);
			}
			if (bufferLength != null) {
				for (int i = 0; i < bufferLength.length; i++) {
					if (bufferLength[i] == 0) {
						if (indexOfEnd >= i) {
							indexOfEnd = i;
						}
						break;
					}
				}
			}
			for (int i = indexOfStart; i <= indexOfEnd; i++) {
				if (bufferLength != null) {
					markerFile.write(buffer[i], 0, bufferLength[i]);
				} else {
					markerFile.write(buffer[i]);
				}
			}
			if (tail != null && tail.length != 0) {
				if (writeTailLength) {
					markerFile.write(Compression.intToBytes(tail.length));
				}
				if (tail.length != 0) {
					markerFile.write(tail);
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}


	/*
	 * Write buffer to Random Access File tail == null means this file is not yet ended; teil.length =
	 * 0 means this file is ended, but there is no content in the tail.
	 */
	public static void writeBufferToRAF(byte[][] buffer, int[] bufferLength, int indexOfStart,
																			int indexOfEnd, RandomAccessFile markerFile, byte[] head,
																			byte[] tail) {
		if (buffer == null || indexOfStart < 0 || indexOfEnd >= buffer.length
				|| indexOfEnd < indexOfStart) {
			System.err.println("\nTranspose Data encoutered the following error: buffer be null, or start index of buffer is negative, or end index is less than the start index, or end index is over the buffer size.");
			System.exit(1);
		}

		try {
			if (head != null && head.length != 0) {
				markerFile.write(head);
			}
			for (int i = indexOfStart; i <= indexOfEnd && bufferLength[i] != 0; i++) {
				markerFile.write(buffer[i], 0, bufferLength[i]);
			}
			if (tail != null) {
				markerFile.writeInt(tail.length);
				if (tail.length != 0) {
					markerFile.write(tail);
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	@SuppressWarnings("unchecked")
	private static Hashtable<String, Float>[] getOutlierHashForEachMdRafFile(Hashtable<String, Float> allOutliers,
																																					 int numMarkerFiles,
																																					 int numMarkersInEachFile,
																																					 String[] sampleNames) {
		Hashtable<String, Float>[] result;
		Enumeration<String> keys;
		String key;
		String[] line;
		int sampleIndex = -1;

		result = new Hashtable[numMarkerFiles];
		for (int i = 0; i < numMarkerFiles; i++) {
			result[i] = new Hashtable<String, Float>();
		}

		keys = allOutliers.keys();

		while (keys.hasMoreElements()) {
			key = keys.nextElement();
			line = key.split("\t");
			for (int i = 0; i < sampleNames.length; i++) {
				if (sampleNames[i].equals(line[1])) {
					sampleIndex = i;
					break;
				}
			}
			if (sampleIndex == -1) {
				System.err.println("Error - Cannot find the sample " + line[1]
													 + " in the project's sample list.");
			}
			result[Integer.parseInt(line[0]) / numMarkersInEachFile].put(
																																	 (Integer.parseInt(line[0])
																																		% numMarkersInEachFile)
																																	 + "\t" + sampleIndex + "\t"
																																	 + line[2],
																																	 allOutliers.get(key));
		}

		return result;
	}

	@SuppressWarnings("unchecked")
	private static Hashtable<String, Float>[] getOutlierHashForEachSampleRafFile(Hashtable<String, Float> allOutliers,
																																							 String[] sampleNamesWholeProj) {
		Hashtable<String, Float>[] result;
		Enumeration<String> keys;
		String key;
		String[] line;

		HashMap<String, Integer> sampleIndexMap = new HashMap<>();
		result = new Hashtable[sampleNamesWholeProj.length];
		for (int i = 0; i < result.length; i++) {
			result[i] = new Hashtable<String, Float>();
			sampleIndexMap.put(sampleNamesWholeProj[i], i);
		}

		keys = allOutliers.keys();
		while (keys.hasMoreElements()) {
			key = keys.nextElement();
			line = key.split("\t");
			if (sampleIndexMap.containsKey(line[1])) {
				result[sampleIndexMap.get(line[1])].put(Integer.parseInt(line[0]) + "\t" + line[2],
																								allOutliers.get(key));
			}
		}

		return result;
	}

	public static byte[] getParameterSectionForMdRaf(int numSamplesInProj,
																									 int numMarkersInCurrentFile, byte nullStatus,
																									 long fingerPrint,
																									 byte[] currentFileMarkerNamesInumBytes) {
		byte[] markerFileHead;

		markerFileHead = new byte[MARKERDATA_PARAMETER_TOTAL_LEN
															+ currentFileMarkerNamesInumBytes.length];
		System.arraycopy(Compression.intToBytes(numSamplesInProj), 0, markerFileHead,
										 MARKERDATA_NUMSAMPLES_START, MARKERDATA_NUMSAMPLES_LEN);
		System.arraycopy(Compression.intToBytes(numMarkersInCurrentFile), 0, markerFileHead,
										 MARKERDATA_NUMMARKERS_START, MARKERDATA_NUMMARKERS_LEN);
		markerFileHead[MARKERDATA_NULLSTATUS_START] = nullStatus;
		System.arraycopy(Compression.longToBytes(fingerPrint), 0, markerFileHead,
										 MARKERDATA_FINGERPRINT_START, MARKERDATA_FINGERPRINT_LEN);
		System.arraycopy(Compression.intToBytes(currentFileMarkerNamesInumBytes.length), 0,
										 markerFileHead, MARKERDATA_MARKERNAMELEN_START, MARKERDATA_MARKERNAMELEN_LEN);
		System.arraycopy(currentFileMarkerNamesInumBytes, 0, markerFileHead,
										 MARKERDATA_MARKERNAME_START, currentFileMarkerNamesInumBytes.length);

		return markerFileHead;
	}

	public static byte[] getParameterSectionForSampRaf(int numMarkersInProj, byte nullStatus,
																										 int numBytesOfOutlierHashtableInCurrentFile,
																										 long fingerprintForMarkers) {
		byte[] markerFileHead;

		markerFileHead = new byte[Sample.PARAMETER_SECTION_BYTES];
		System.arraycopy(Compression.intToBytes(numMarkersInProj), 0, markerFileHead,
										 Sample.PARAMETER_SECTION_NUMMARKERS_LOCATION,
										 Sample.PARAMETER_SECTION_NUMMARKERS_LENGTH);
		markerFileHead[Sample.PARAMETER_SECTION_NULLSTAT_LOCATION] = nullStatus;
		System.arraycopy(Compression.intToBytes(numBytesOfOutlierHashtableInCurrentFile), 0,
										 markerFileHead, Sample.PARAMETER_SECTION_OUTLIERSECTIONLENGTH_LOCATION,
										 Sample.PARAMETER_SECTION_OUTLIERSECTIONLENGTH_LENGTH);
		System.arraycopy(Compression.longToBytes(fingerprintForMarkers), 0, markerFileHead,
										 Sample.PARAMETER_SECTION_FINGPRNT_LOCATION,
										 Sample.PARAMETER_SECTION_FINGPRNT_LENGTH);

		return markerFileHead;
	}

	public static int getOptimalNumSamplesBasedOnHeapSpace(int currentNumSamplesInHeapSpace,
																												 int numBytesPerSample) {
		if (currentNumSamplesInHeapSpace <= 0) {
			return (int) ((Runtime.getRuntime().maxMemory() * 0.75) / numBytesPerSample);
		} else {
			return currentNumSamplesInHeapSpace *= 0.9;
		}
	}


	public static byte getNullstatusFromRandomAccessFile(String filename, boolean isSampRafOrMdRaf) {
		byte nullStatusOfTheFile = Byte.MIN_VALUE;
		RandomAccessFile rafFile;

		try {
			rafFile = new RandomAccessFile(filename, "r");
			if (isSampRafOrMdRaf) {
				rafFile.readInt();
				nullStatusOfTheFile = rafFile.readByte();
			} else {
				rafFile.readInt();
				rafFile.readInt();
				nullStatusOfTheFile = rafFile.readByte();
			}
			rafFile.close();
		} catch (FileNotFoundException e1) {
			e1.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return nullStatusOfTheFile;
	}


	public static void deleteOlderRafs(String directory, final String[] fileNamePrefixes,
																		 final String[] fileNameExtesions,
																		 final boolean prefixAndSuffix, final String[] excludes) {
		File[] files;

		// List of files to be backed up.
		files = new File(directory).listFiles(new FilenameFilter() {
			@Override
			public boolean accept(File file, String filename) {
				boolean result = false;
				filename = ext.removeDirectoryInfo(filename);
				if (prefixAndSuffix) {
					if (fileNamePrefixes.length != fileNameExtesions.length) {
						System.err.println("Array fileNamePrefixes and Array fileNameExtesions must have the same number of elements.");
						System.exit(1);
					}
					for (int i = 0; !result && fileNameExtesions != null
													&& i < fileNameExtesions.length; i++) {
						if (filename.startsWith(fileNamePrefixes[i])
								&& filename.endsWith(fileNameExtesions[i])) {
							result = true;
						}
					}
				} else {
					for (int i = 0; !result && fileNameExtesions != null
													&& i < fileNameExtesions.length; i++) {
						if (filename.endsWith(fileNameExtesions[i])) {
							result = true;
						}
					}
					for (int i = 0; !result && fileNamePrefixes != null && i < fileNamePrefixes.length; i++) {
						if (filename.startsWith(fileNamePrefixes[i])) {
							result = true;
						}
					}
				}
				if (result && excludes != null) {
					for (String exclude : excludes) {
						if (ext.removeDirectoryInfo(filename).equals(exclude)) {
							result = false;
							break;
						}
					}
				}
				return result;
			}
		});

		// Move files to backup folder
		for (File file : files) {
			file.delete();
		}

		if (files.length > 0) {
			System.out.println("Older version of the data in " + directory
												 + " has been deleted from the hard drive.");
		}
	}

	public static byte backupOlderFiles(String directory,
																			final String[] nameExtOrSuffixOfFilesToBackup,
																			boolean backupToFolderOrRename) {
		return backupOlderFiles(directory, null, nameExtOrSuffixOfFilesToBackup,
														backupToFolderOrRename);
	}

	public static byte backupOlderFiles(String directory, final String[] namePrefixOfFilesToBackup,
																			final String[] nameExtOrSuffixOfFilesToBackup,
																			boolean backupToFolderOrRename) {
		File[] files;
		byte count;
		String filename;

		// List of files to be backed up.
		files = new File(directory).listFiles(new FilenameFilter() {
			@Override
			public boolean accept(File file, String filename) {
				boolean result = false;
				filename = ext.removeDirectoryInfo(filename);
				for (int i = 0; !result && nameExtOrSuffixOfFilesToBackup != null
												&& i < nameExtOrSuffixOfFilesToBackup.length; i++) {
					if (filename.endsWith(nameExtOrSuffixOfFilesToBackup[i])) {
						result = true;
					}
				}
				for (int i = 0; !result && namePrefixOfFilesToBackup != null
												&& i < namePrefixOfFilesToBackup.length; i++) {
					if (filename.startsWith(namePrefixOfFilesToBackup[i])) {
						result = true;
					}
				}
				return result;
			}
		});

		count = 0;
		if (files.length > 0) {
			// Create a new backup folder.
			if (backupToFolderOrRename) {
				do {
					count++;
				} while (new File(directory + "Backup." + count).exists());
				new File(directory + "Backup." + count + "/").mkdirs();

				// Move files to backup folder
				for (File file : files) {
					file.renameTo(new File(directory + "Backup." + count + "/" + file.getName()));
				}

				System.out.println("Older version of the data in " + directory + " has been moved to "
													 + directory + "Backup." + count
													 + "/ To save disk space, please manually delete these files.");
			} else {
				for (File file : files) {
					filename = file.getName();
					file.renameTo(new File(file.getParent() + ext.rootOf(filename) + "_"
																 + (new SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date()))
																 + "." + filename.substring(filename.lastIndexOf(".") + 1,
																														filename.length())));
				}
			}
		}
		return count;
	}

	public static void backupOlderFile(String fileFullPath, byte i) {
		new File(fileFullPath).renameTo(new File(ext.parseDirectoryOfFile(fileFullPath)
																						 + ext.rootOf(fileFullPath) + "_Backup" + i + "."
																						 + fileFullPath.substring(fileFullPath.lastIndexOf(".")
																																			+ 1, fileFullPath.length())));
	}

	@SuppressWarnings("unchecked")
	public static Hashtable<String, Float> loadOutliersFromRAF(String mdraf) {
		RandomAccessFile file;
		int numBytesPerMarker;
		byte[] parameters;
		int markernamesSectionLength;
		Hashtable<String, Float> outOfRangeValues = null;
		byte nullStatus = 0;
		byte numBytesPerSampleMarker = 0;
		int nSamples;
		int numMarkersInThisFile;
		int lengthOfOutOfRangeHashtable;

		try {
			file = new RandomAccessFile(mdraf, "r");
			parameters = new byte[TransposeData.MARKERDATA_PARAMETER_TOTAL_LEN];
			file.read(parameters);
			nSamples = Compression.bytesToInt(parameters, TransposeData.MARKERDATA_NUMSAMPLES_START);
			numMarkersInThisFile = Compression.bytesToInt(parameters,
																										TransposeData.MARKERDATA_NUMMARKERS_START);
			nullStatus = parameters[TransposeData.MARKERDATA_NULLSTATUS_START];
			numBytesPerSampleMarker = Sample.getNBytesPerSampleMarker(nullStatus);
			numBytesPerMarker = numBytesPerSampleMarker * nSamples;
			markernamesSectionLength = Compression.bytesToInt(parameters, MARKERDATA_MARKERNAMELEN_START);

			parameters = new byte[markernamesSectionLength];
			file.read(parameters);

			file.seek((long) TransposeData.MARKERDATA_PARAMETER_TOTAL_LEN
								+ (long) markernamesSectionLength
								+ (long) numMarkersInThisFile * (long) numBytesPerMarker);
			lengthOfOutOfRangeHashtable = file.readInt();
			if (lengthOfOutOfRangeHashtable > 0) {
				parameters = new byte[lengthOfOutOfRangeHashtable];
				file.read(parameters);
				outOfRangeValues = (Hashtable<String, Float>) Compression.bytesToObj(parameters);
			}
			file.close();

		} catch (FileNotFoundException e) {
			System.err.println("Error - could not find RAF marker file '" + mdraf + "'");
			e.printStackTrace();
		} catch (IOException e) {
			System.err.println("Error reading RAF marker file '" + mdraf + "'");
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
		}

		return outOfRangeValues;
	}

	public static String[] loadMarkerNamesFromRAF(String markerFilename) throws IOException,
																																			 ClassNotFoundException {
		RandomAccessFile file;
		byte[] parameters;
		file = new RandomAccessFile(markerFilename, "r");
		parameters = new byte[TransposeData.MARKERDATA_PARAMETER_TOTAL_LEN];
		file.read(parameters);
		int markernamesSectionLength = Compression.bytesToInt(parameters,
																													MARKERDATA_MARKERNAMELEN_START);

		parameters = new byte[markernamesSectionLength];
		file.read(parameters);
		String[] markerNames = (String[]) Compression.bytesToObj(parameters);
		file.close();
		return markerNames;
	}

	@SuppressWarnings("unchecked")
	public static MarkerData[] loadFromRAF(String markerFilename, int[] targertMarkIndicesInFile) {
		MarkerData[] result = null;
		RandomAccessFile file;
		int numBytesPerMarker;
		float[] gcs = null;
		float[] xs = null;
		float[] ys = null;
		float[] bafs = null;
		float[] lrrs = null;
		byte[] abGenotypes = null;
		byte[] forwardGenotypes = null;
		byte[] genotypeTmp;
		long seekLocation;
		byte[][] readBuffer;
		int indexReadBuffer;
		byte[] parameters;
		int markernamesSectionLength;
		Hashtable<String, Float> outOfRangeValues = null;
		byte nullStatus = 0;
		byte numBytesPerSampleMarker = 0;
		int indexStart;
		int nSamples;
		int numMarkersInThisFile;
		long fingerPrint;
		String[] markerNames;
		int lengthOfOutOfRangeHashtable;
		boolean isGcNull, isXNull, isYNull, isBafNull, isLrrNull, isGenotypeNull, isNegativeXOrYAllowed;

		try {
			file = new RandomAccessFile(markerFilename, "r");
			parameters = new byte[TransposeData.MARKERDATA_PARAMETER_TOTAL_LEN];
			file.read(parameters);
			nSamples = Compression.bytesToInt(parameters, TransposeData.MARKERDATA_NUMSAMPLES_START);
			numMarkersInThisFile = Compression.bytesToInt(parameters,
																										TransposeData.MARKERDATA_NUMMARKERS_START);
			nullStatus = parameters[TransposeData.MARKERDATA_NULLSTATUS_START];
			isGcNull = Sample.isGcNull(nullStatus);
			isXNull = Sample.isXNull(nullStatus);
			isYNull = Sample.isYNull(nullStatus);
			isBafNull = Sample.isBafNull(nullStatus);
			isLrrNull = Sample.isLrrNull(nullStatus);
			isGenotypeNull = Sample.isAbAndForwardGenotypeNull(nullStatus);
			isNegativeXOrYAllowed = Sample.isNegativeXOrYAllowed(nullStatus);
			numBytesPerSampleMarker = Sample.getNBytesPerSampleMarker(nullStatus);
			numBytesPerMarker = numBytesPerSampleMarker * nSamples;
			fingerPrint = Compression.bytesToLong(parameters, MARKERDATA_FINGERPRINT_START);
			markernamesSectionLength = Compression.bytesToInt(parameters, MARKERDATA_MARKERNAMELEN_START);

			parameters = new byte[markernamesSectionLength];
			file.read(parameters);
			markerNames = (String[]) Compression.bytesToObj(parameters);

			result = new MarkerData[targertMarkIndicesInFile.length];
			readBuffer = new byte[targertMarkIndicesInFile.length][numBytesPerMarker];

			for (int i = 0; i < targertMarkIndicesInFile.length; i++) {
				if (targertMarkIndicesInFile[i] < 0
						|| targertMarkIndicesInFile[i] >= numMarkersInThisFile) {
					System.err.println("Skipped the marker index " + targertMarkIndicesInFile[i]
														 + ", because it is out of range.");
				} else {
					seekLocation = (long) TransposeData.MARKERDATA_PARAMETER_TOTAL_LEN
												 + (long) markernamesSectionLength
												 + targertMarkIndicesInFile[i] * (long) numBytesPerMarker;
					file.seek(seekLocation);
					file.read(readBuffer[i]);
				}
			}

			// TODO this is read every time, wouldn't it be faster to use the serialized version?
			file.seek((long) TransposeData.MARKERDATA_PARAMETER_TOTAL_LEN
								+ (long) markernamesSectionLength
								+ (long) numMarkersInThisFile * (long) numBytesPerMarker);
			lengthOfOutOfRangeHashtable = file.readInt();
			if (lengthOfOutOfRangeHashtable > 0) {
				parameters = new byte[lengthOfOutOfRangeHashtable];
				file.read(parameters);
				outOfRangeValues = (Hashtable<String, Float>) Compression.bytesToObj(parameters);
			}
			file.close();

			for (int i = 0; i < targertMarkIndicesInFile.length; i++) {
				indexReadBuffer = 0;

				indexStart = 0;
				indexReadBuffer = indexStart;
				if (!isGcNull) {
					gcs = new float[nSamples];
					for (int j = 0; j < nSamples; j++) {
						gcs[j] = Compression.gcBafDecompress(new byte[] {readBuffer[i][indexReadBuffer],
																														 readBuffer[i][indexReadBuffer + 1]});
						indexReadBuffer += numBytesPerSampleMarker;
					}
					indexStart += 2;
				}
				indexReadBuffer = indexStart;
				if (!isXNull) {
					xs = new float[nSamples];
					for (int j = 0; j < nSamples; j++) {
						if (isNegativeXOrYAllowed) {
							xs[j] = Compression.xyDecompressAllowNegative(new byte[] {
																																				readBuffer[i][indexReadBuffer],
																																				readBuffer[i][indexReadBuffer
																																											+ 1]});
						} else {
							xs[j] = Compression.xyDecompressPositiveOnly(new byte[] {
																																			 readBuffer[i][indexReadBuffer],
																																			 readBuffer[i][indexReadBuffer
																																										 + 1]});
						}
						if (xs[j] == Compression.REDUCED_PRECISION_XY_OUT_OF_RANGE_FLAG_FLOAT) {
							xs[j] = outOfRangeValues.get(targertMarkIndicesInFile[i] + "\t" + j + "\tx");
						}
						indexReadBuffer += numBytesPerSampleMarker;
					}
					indexStart += 2;
				}
				indexReadBuffer = indexStart;
				if (!isYNull) {
					ys = new float[nSamples];
					for (int j = 0; j < nSamples; j++) {
						if (isNegativeXOrYAllowed) {
							ys[j] = Compression.xyDecompressAllowNegative(new byte[] {
																																				readBuffer[i][indexReadBuffer],
																																				readBuffer[i][indexReadBuffer
																																											+ 1]});
						} else {
							ys[j] = Compression.xyDecompressPositiveOnly(new byte[] {
																																			 readBuffer[i][indexReadBuffer],
																																			 readBuffer[i][indexReadBuffer
																																										 + 1]});
						}
						if (ys[j] == Compression.REDUCED_PRECISION_XY_OUT_OF_RANGE_FLAG_FLOAT) {
							ys[j] = outOfRangeValues.get(targertMarkIndicesInFile[i] + "\t" + j + "\ty");
						}
						indexReadBuffer += numBytesPerSampleMarker;
					}

					indexStart += 2;
				}
				indexReadBuffer = indexStart;
				if (!isBafNull) {
					bafs = new float[nSamples];
					for (int j = 0; j < nSamples; j++) {
						bafs[j] = Compression.gcBafDecompress(new byte[] {readBuffer[i][indexReadBuffer],
																															readBuffer[i][indexReadBuffer + 1]});
						indexReadBuffer += numBytesPerSampleMarker;
					}
					indexStart += 2;
				}
				indexReadBuffer = indexStart;
				if (!isLrrNull) {
					lrrs = new float[nSamples];
					for (int j = 0; j < nSamples; j++) {
						lrrs[j] = Compression.lrrDecompress(new byte[] {readBuffer[i][indexReadBuffer],
																														readBuffer[i][indexReadBuffer + 1],
																														readBuffer[i][indexReadBuffer + 2]});
						if (lrrs[j] == Compression.REDUCED_PRECISION_LRR_OUT_OF_RANGE_LRR_FLAG_FLOAT) {
							lrrs[j] = outOfRangeValues.get(targertMarkIndicesInFile[i] + "\t" + j + "\tlrr");
						}
						indexReadBuffer += numBytesPerSampleMarker;
					}
					indexStart += 3;
				}
				indexReadBuffer = indexStart;
				if (!isGenotypeNull) {
					abGenotypes = new byte[nSamples];
					forwardGenotypes = new byte[nSamples];
					for (int j = 0; j < nSamples; j++) {
						genotypeTmp = Compression.genotypeDecompress(readBuffer[i][indexReadBuffer]);
						abGenotypes[j] = genotypeTmp[0];
						forwardGenotypes[j] = genotypeTmp[1];
						indexReadBuffer += numBytesPerSampleMarker;
					}
				}
				result[i] = new MarkerData(markerNames[i], (byte) 0, 0, fingerPrint, gcs, null, null, xs,
																	 ys, null, null, bafs, lrrs, abGenotypes, forwardGenotypes);
			}
		} catch (FileNotFoundException e) {
			System.err.println("Error - could not find RAF marker file '" + markerFilename + "'");
			e.printStackTrace();
		} catch (IOException e) {
			System.err.println("Error reading RAF marker file '" + markerFilename + "'");
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
		}
		return result;
	}



	public static void showHeap(String caller) {
		System.out.println(caller + "\tmaxMem: "
											 + ext.prettyUpSize(Runtime.getRuntime().maxMemory(), 1)
											 + "\tfreeMem: " + ext.prettyUpSize(Runtime.getRuntime().freeMemory(), 1)
											 + "\ttotalMem: " + ext.prettyUpSize(Runtime.getRuntime().totalMemory(), 1));
	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		Project proj;
		String filename = null;
		boolean transpose = false;
		boolean reversetranspose = false;
		boolean keepFilesOpen = false;
		int maxFileSize = 0;
		boolean lookup = false;

		String usage = "\n" + "TransposeData requires 0-1 arguments\n"
									 + "   (1) project properties filename (i.e. proj="
									 + org.genvisis.cnv.Launch.getDefaultDebugProjectFile(false) + " (default))\n"
									 + "   (2) transpose data (i.e. -transpose (" + (transpose ? "" : "not the ")
									 + "default))\n" + "   (3) keep all files open at once (i.e. -keepFilesOpen ("
									 + (keepFilesOpen ? "" : "not the ")
									 + "default; not recommended usually allowed on linux servers))\n"
									 + "   (4) maximum size of each file in bytes (i.e. max=" + maxFileSize
									 + " (default))\n" + "  OR:\n" + "   (1) project file (i.e. proj=" + filename
									 + " (default))\n" + "   (5) reversetranspose data (i.e. -reversetranspose ("
									 + (reversetranspose ? "" : "not the ") + "default))\n" + "  OR:\n"
									 + "   (1) project file (i.e. proj=" + filename + " (default))\n"
									 + "   (6) create marker lookup table (i.e. -lookup ("
									 + (lookup ? "" : "not the ") + "default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("proj=")) {
				filename = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("-transpose")) {
				transpose = true;
				numArgs--;
			} else if (arg.startsWith("-keepFilesOpen")) {
				keepFilesOpen = true;
				numArgs--;
			} else if (arg.startsWith("max=")) {
				maxFileSize = Integer.parseInt(arg.split("=")[1]);
				numArgs--;
			} else if (arg.startsWith("-lookup")) {
				lookup = true;
				numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}

		proj = new Project(filename);

		// transpose = true;
		// lookup = true;
		try {
			if (transpose) {
				transposeData(proj, maxFileSize, keepFilesOpen);
			} else if (reversetranspose) {
				reverseTranspose(proj);
			} else if (lookup) {
				recreateMarkerLookup(proj);
			} else {
				System.err.println("Error - must spefify one of the following: -transpose, or -lookup. Refer to the following usage for details.");
				System.out.println(usage);
			}
		} catch (Exception e) {
			System.err.println(ext.getTime() + "\tExcepted");
			e.printStackTrace();
		}
	}
}
