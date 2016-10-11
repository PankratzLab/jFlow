package org.genvisis.cnv.manage;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.Arrays;
import java.util.Date;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.cnv.filesys.Compression;
import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.MarkerLookup;
import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.common.Array;
import org.genvisis.common.CountHash;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.ext;

public class MarkerDataLoader implements Runnable {
	private Project proj;
	private String[] markerNames;
	private byte[] chrs;
	private int[] positions;
	private MarkerData[] markerData;
	private boolean[] loaded;

	private Hashtable<String, String> filenames;
	private Hashtable<String, Vector<String>> hash;
	private MarkerLookup markerLookup;
	private int currentIndexBeingLoaded;
	private int currentDirection;
	private boolean initiated;
	private boolean killed;
	private boolean killComplete;
	private long sampleFingerprint;
	private int readAheadLimit;
	private int numberCurrentlyLoaded;
	// private boolean plinkFormat;
	// private String plinkFileRoot;
	// private int[] plinkSampleIndices;
	private Thread thread;
	private Logger log;
	private CountHash finalWaitTimeCounts;
	private HashSet<Integer> waitTimesSeen;

	public MarkerDataLoader(Project proj, String[] markerNames, int amountToLoadAtOnceInMB) {
		Hashtable<String, Integer> markerHash;
		Vector<String> missingMarkers, v;
		String[] markerNamesProj;
		int[] positionsProj;
		byte[] chrsProj;
		String[] line;
		int index;

		this.proj = proj;
		log = proj.getLog();
		this.markerNames = markerNames;

		if (markerNames == null) {
			log.reportError("The list of markers for MarkerDataLoader to load was null");
			return;
		}

		if (markerNames.length == 0) {
			log.reportError("The list of markers for MarkerDataLoader to load was empty (n=0)");
			killed = true;
			return;
		}

		chrs = new byte[markerNames.length];
		positions = new int[markerNames.length];
		markerData = new MarkerData[markerNames.length];
		loaded = new boolean[markerNames.length];
		sampleFingerprint = proj.getSampleList().getFingerprint();
		numberCurrentlyLoaded = 0;
		initiated = false;
		waitTimesSeen = new HashSet<Integer>();
		finalWaitTimeCounts = new CountHash();

		markerNamesProj = proj.getMarkerNames();
		chrsProj = proj.getMarkerSet().getChrs();
		positionsProj = proj.getMarkerSet().getPositions();
		markerHash = HashVec.loadToHashIndices(markerNames, log);
		for (int i = 0; i < markerNamesProj.length; i++) {
			if (markerHash.containsKey(markerNamesProj[i])) {
				index = markerHash.get(markerNamesProj[i]);
				chrs[index] = chrsProj[i];
				positions[index] = positionsProj[i];
				break;
			}
		}

		// if (plinkFileRoot != null) {
		// this.plinkFormat = true;
		// this.plinkFileRoot = plinkFileRoot;
		// markerLookup = PlinkData.parseMarkerLookup(plinkFileRoot); // Note: markerLookup from
		// PlinkData contains alleles, unlike its standard counterpart
		// plinkSampleIndices = PlinkData.parseSampleIndicesForProject(proj, plinkFileRoot);
		// } else {
		// this.plinkFormat = false;
		// this.plinkFileRoot = null;
		markerLookup = proj.getMarkerLookup();
		// }

		if (amountToLoadAtOnceInMB <= 0) {
			amountToLoadAtOnceInMB = (int) ((double) Runtime.getRuntime().maxMemory()	/ 1024 / 1024
																			* 0.80);
			log.report("80% of max memory available ("
										+ ext.prettyUpSize(Runtime.getRuntime().maxMemory(), 1)
									+ ") will be used by MarkerDataLoader: ");
		}
		readAheadLimit = -1;

		missingMarkers = new Vector<String>();
		filenames = new Hashtable<String, String>();
		hash = new Hashtable<String, Vector<String>>();
		for (String markerName : markerNames) {
			if (markerLookup.contains(markerName)) {
				line = markerLookup.get(markerName).split("[\\s]+");
				if (hash.containsKey(line[0])) {
					v = hash.get(line[0]);
				} else {
					hash.put(line[0], v = new Vector<String>(100000));
					filenames.put(line[0], "");
				}
				// if (plinkFormat) {
				// v.add(markerNames[i]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]);
				// } else {
				v.add(markerName + "\t" + line[1]);
				// }

				if (readAheadLimit == -1) {
					readAheadLimit = (int) Math.floor((double) amountToLoadAtOnceInMB	* 1024 * 1024
																						/ determineNBytesPerMarker(proj.MARKER_DATA_DIRECTORY.getValue(	false,
																																																						true)
																																				+ line[0], log));
					if (readAheadLimit > markerNames.length) {
						log.report("Read ahead limit was computed to be greater than the number of markers; all markers will be loaded into memory at once");
					} else {
						log.report("Read ahead limit was computed to be "	+ ext.addCommas(readAheadLimit)
												+ " markers at a time");
					}
				}
			} else {
				missingMarkers.add(markerName);
			}
		}
		if (missingMarkers.size() > 0) {
			proj.message("Error - the following markers were not found in the MarkerSet: "
										+ Array.toStr(Array.toStringArray(missingMarkers), " "));
		}
		currentIndexBeingLoaded = 0;
		currentDirection = +1;
		killed = false;
	}

	// public void setChrsPositions(byte[] chrs, int[] positions) {
	// this.chrs = chrs;
	// this.positions = positions;
	// }

	public static int determineNBytesPerMarker(String filename, Logger log) {
		RandomAccessFile file;
		byte[] parameterReadBuffer;
		byte nullStatus;
		int nSampObserved, nBytesPerSampleMarker;

		parameterReadBuffer = new byte[TransposeData.MARKERDATA_PARAMETER_TOTAL_LEN];
		try {
			file = new RandomAccessFile(filename, "r");
			file.read(parameterReadBuffer);
			nullStatus = parameterReadBuffer[TransposeData.MARKERDATA_NULLSTATUS_START];
			// bytesPerSampleMarker = Compression.BYTES_PER_SAMPLE_MARKER - (nullStatus & 0x01) -
			// (nullStatus >>1 & 0x01) - (nullStatus >>2 & 0x01) - (nullStatus >>3 & 0x01) - (nullStatus
			// >>4 & 0x01) - (nullStatus >>5 & 0x01) - (nullStatus >>6 & 0x01);
			nBytesPerSampleMarker = Sample.getNBytesPerSampleMarker(nullStatus);
			nSampObserved = Compression.bytesToInt(	parameterReadBuffer,
																							TransposeData.MARKERDATA_NUMSAMPLES_START);
			file.close();
			return nSampObserved * nBytesPerSampleMarker;

		} catch (FileNotFoundException e) {
			log.reportError("Error - could not find RAF marker file '" + filename + "'");
			e.printStackTrace();
		} catch (IOException e) {
			log.reportError("Error reading RAF marker file '" + filename + "'");
			e.printStackTrace();
		}

		return -1;
	}

	// public boolean hasIndexBeenLoaded(int index) {
	// return markerData[index] != null;
	// }

	public void requestIndexBeTheNextFileToLoad(int indexRequested) {
		if (indexRequested < currentIndexBeingLoaded) {
			currentDirection = -1;
		} else {
			currentDirection = +1;
		}
		currentIndexBeingLoaded = indexRequested;
	}

	public void releaseIndex(int index) {
		if (markerData[index] == null) {
			log.report("Index " + index + " is not currently active; already null");
		} else {
			markerData[index] = null;
			numberCurrentlyLoaded--;
		}
	}

	public void kill() {
		killed = true;
		if (!thread.isAlive()) {
			killComplete = true;
		}
	}

	public boolean killComplete() {
		return killComplete;
	}

	public Thread getThread() {
		return thread;
	}

	public MarkerData getMarkerData(int markerIndex) {
		if (!initiated) {
			log.reportError("Error - cannot getMarkerData before the loader has been launched. Check to see if the MarkerDataLoader constructor was used instead of one of the loadMarkerDataFromListIn* methods.");
			return null;
		}
		if (loaded[markerIndex]) {
			return markerData[markerIndex];
		} else {
			return null;
		}
	}

	public MarkerData requestMarkerData(int markerIndex) {
		MarkerData markerData;
		int count;

		if (!initiated) {
			log.reportError("Error - cannot requestMarkerData before the loader has been launched. Check to see if the MarkerDataLoader constructor was used instead of one of the loadMarkerDataFromListIn* methods.");
			return null;
		}

		count = 0;
		while ((markerData = getMarkerData(markerIndex)) == null) {
			requestIndexBeTheNextFileToLoad(markerIndex);
			try {
				Thread.sleep(250);
				count++;
			} catch (InterruptedException ie) {
			}
			if (count > 8 && count % 8 == 0) {
				if (!waitTimesSeen.contains(count)) {
					log.reportError("Have been waiting on markerDataLoader to load "
													+ markerNames[markerIndex] + " for " + (count / 4) + " seconds");
					waitTimesSeen.add(count);
				}
			}
		}

		finalWaitTimeCounts.add(count + "");

		return markerData;
	}

	public void reportWaitTimes() {
		String[] values;
		int[] counts;

		finalWaitTimeCounts.sortValuesAlphanumerically();
		values = finalWaitTimeCounts.getValues();
		counts = finalWaitTimeCounts.getCounts();

		log.report("Final wait times for markerDataLoader:");
		log.report("All markers:\t" + markerNames.length);
		for (int i = 0; i < values.length; i++) {
			log.report(ext.formDeci(Double.parseDouble(values[i]) / 4.0, 2, true)	+ " seconds:\t"
									+ counts[i]);
		}
	}

	@Override
	@SuppressWarnings("unchecked")
	public void run() {
		MarkerData[] collection;
		String[] line;
		Vector<String> v, allRemaining;
		String filename;
		int[] markerIndicesInProj = null;
		int[] markerIndicesInFile;
		int[][] markerIndicesInSelection;
		long fingerprint;
		int count;
		MarkerSet markerSet;
		String[] allMarkersInProj;
		byte[] allChrsInProj;
		int[] allPosInProj;
		String[] allSampsInProj;
		long time, subtime;
		int maxPerCycle;
		Hashtable<String, Float> outlierHash;
		// String[][] plinkMarkerAlleles;

		if (killed) {
			return;
		}

		initiate();
		// maxPerCycle = proj.getInt(proj.MAX_MARKERS_LOADED_PER_CYCLE);
		maxPerCycle = proj.MAX_MARKERS_LOADED_PER_CYCLE.getValue();

		fingerprint = proj.getSampleList().getFingerprint();

		markerSet = proj.getMarkerSet();
		allMarkersInProj = markerSet.getMarkerNames();
		allChrsInProj = markerSet.getChrs();
		allPosInProj = markerSet.getPositions();
		allSampsInProj = proj.getSamples();
		if (new File(proj.MARKER_DATA_DIRECTORY.getValue(false, true) + "outliers.ser").exists()) {
			outlierHash = (Hashtable<String, Float>) SerializedFiles.readSerial(
																																					proj.MARKER_DATA_DIRECTORY.getValue(false,
																																																							true)
																																					+ "outliers.ser");
		} else {
			outlierHash = new Hashtable<String, Float>();
		}
		time = subtime = new Date().getTime();
		count = 0;
		// log.report("filenames.size(): " + filenames.size());
		while (filenames.size() > 0) {
			while (loaded[currentIndexBeingLoaded]) {
				if (currentIndexBeingLoaded + 1 == markerNames.length) {
					currentDirection = -1;
				}
				if (currentIndexBeingLoaded - 1 == 0) {
					currentDirection = +1;
				}
				currentIndexBeingLoaded += currentDirection;
			}
			filename = markerLookup.get(markerNames[currentIndexBeingLoaded]).split("\t")[0];
			if (!Files.exists(proj.MARKER_DATA_DIRECTORY.getValue(false, true)	+ filename,
												proj.JAR_STATUS.getValue())) {
				proj.message("Error - could not load data from '"
												+ proj.MARKER_DATA_DIRECTORY.getValue(false, true) + filename
											+ "'; because the file could not be found");
				return;
			}

			allRemaining = hash.get(filename);

			v = new Vector<String>();
			for (int i = 0; allRemaining.size() > 0 && i < maxPerCycle; i++) {
				v.addElement(allRemaining.remove(0));
			}
			log.report("Loaded up "	+ v.size() + " from MarkerData file:" + filename + "; "
									+ allRemaining.size() + " remaining for that file; "
									+ ext.getTimeElapsed(subtime), true, true, 8);
			subtime = new Date().getTime();
			if (allRemaining.size() == 0) {
				filenames.remove(filename);
			}
			markerIndicesInFile = new int[v.size()];
			markerIndicesInProj = new int[v.size()];
			markerIndicesInSelection = new int[v.size()][];
			// plinkMarkerAlleles = new String[v.size()][];
			for (int j = 0; j < v.size() && !killed; j++) {
				line = v.elementAt(j).split("[\\s]+");
				markerIndicesInProj[j] = ext.indexOfStr(line[0], allMarkersInProj, true, true); // modified
																																												// here to
																																												// fix the
																																												// bug
				markerIndicesInFile[j] = Integer.parseInt(line[1]);
				markerIndicesInSelection[j] = ext.indicesOfStr(line[0], markerNames, true, true);
				if (markerIndicesInSelection[j].length > 1) {
					log.report("FYI, marker "	+ line[0] + " was requested "
											+ markerIndicesInSelection[j].length + " times");
				}
				// if (plinkFormat) {
				// plinkMarkerAlleles[j] = new String[] {line[2], line[3]};
				// }
			}

			// if (plinkFormat) {
			// collection = PlinkData.loadBedUsingRAF(allMarkersInProj, allChrsInProj, allPosInProj,
			// allSampsInProj, proj.MARKER_DATA_DIRECTORY.getValue(false, true) + plinkFileRoot + ".bed",
			// markerIndicesInProj, markerIndicesInFile, sampleFingerprint, plinkSampleIndices);
			// collection = PlinkData.loadPedUsingRAF(allMarkersInProj, allChrsInProj, allPosInProj,
			// allSampsInProj, proj.getDir(proj.MARKER_DATA_DIRECTORY) + plinkFileRoot + ".bed",
			// markerIndicesInProj, markerIndicesInFile, sampleFingerprint, plinkSampleIndices);
			// } else {
			collection = loadFromRAF(	allMarkersInProj, allChrsInProj, allPosInProj, allSampsInProj,
																proj.MARKER_DATA_DIRECTORY.getValue(false, true) + filename,
																markerIndicesInProj, markerIndicesInFile, sampleFingerprint,
																outlierHash, log);
			// }

			for (int k = 0; k < markerIndicesInProj.length && !killed; k++) {
				for (int i = 0; i < markerIndicesInSelection[k].length; i++) {
					markerData[markerIndicesInSelection[k][i]] = collection[k];
					loaded[markerIndicesInSelection[k][i]] = true;
					count++;
					numberCurrentlyLoaded++;
					if (markerData[markerIndicesInSelection[k][i]].getFingerprint() != fingerprint) {
						log.reportError("Error - mismatched fingerprint after MarkerLookup. Actual in MarkerData: "
															+ markerData[markerIndicesInSelection[k][i]].getFingerprint()
														+ ", while expecting: " + fingerprint);
					}
				}
			}

			while (!killed && numberCurrentlyLoaded > readAheadLimit) {
				try {
					Thread.sleep(500);
				} catch (InterruptedException ie) {
				}
				log.report("Currently loading index "	+ currentIndexBeingLoaded
										+ " but readAheadLimit has been reached (" + numberCurrentlyLoaded + " over "
										+ readAheadLimit + "); so waiting");
			}
			if (killed) {
				filenames = new Hashtable<String, String>();
			}

		}
		if (killed) {
			log.report("MarkerDataLoader killed");
			markerNames = null;
			loaded = null;
			markerData = null;
			System.gc();
			killComplete = true;
		} else {
			log.report("Independent thread has finished loading "	+ count + " markers to MarkerData[] in "
									+ ext.getTimeElapsed(time));
		}

	}

	public Hashtable<String, Vector<String>> getBatches() {
		return hash;
	}


	@SuppressWarnings("unchecked")
	public static Hashtable<String, Float> loadOutliers(Project proj) {
		if (new File(proj.MARKER_DATA_DIRECTORY.getValue(false, true) + "outliers.ser").exists()) {
			return (Hashtable<String, Float>) SerializedFiles.readSerial(proj.MARKER_DATA_DIRECTORY.getValue(	false,
																																																				true)
																																		+ "outliers.ser");
		} else {
			return new Hashtable<String, Float>();
		}
	}

	public static MarkerData[] loadFromRAF(	String[] allMarkersInProj, byte[] allChrsProj,
																					int[] allPositionsInProj, String[] allSampsInProj,
																					String currentMarkFilename, int[] markerIndicesInProj,
																					int[] markerIndcesInFile, long sampleFingerprint,
																					Hashtable<String, Float> outOfRangeValues, Logger log) {
		return loadFromRAF(	allMarkersInProj, allChrsProj, allPositionsInProj, allSampsInProj,
												currentMarkFilename, markerIndicesInProj, markerIndcesInFile, true, true,
												true, true, true, sampleFingerprint, outOfRangeValues, log);
	}


	/**
	 * Load MarkerData of selected markers from a single Random Access File of 12-byte
	 * (half-precision) format that has incorporated the marker names. Load marker data (RAF) approach
	 * 7b.
	 *
	 * @param proj
	 * @param markerNamesOfInterest
	 * @return
	 */
	public static MarkerData[] loadFromRAF(	String[] allMarkersInProj, byte[] allChrsInProj,
																					int[] allPosInProj, String[] allSampsInProj,
																					String currentMarkFilename, int[] markersIndicesInProj,
																					int[] markersIndicesInFile, boolean loadGC,
																					boolean loadXY, boolean loadBAF, boolean loadLRR,
																					boolean loadAbGenotype, long sampleFingerprint,
																					Hashtable<String, Float> outOfRangeValues, Logger log) {
		MarkerData[] result;
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
		byte[][] readBuffer = null;
		int indexReadBuffer;
		byte[] parameterReadBuffer;
		long fingerprint;
		int numBytesMarkernamesSection;
		// Hashtable<String, Float> outOfRangeValues = null;
		int numSamplesProj;
		byte nullStatus = 0;
		byte bytesPerSampleMarker = 0;
		int indexStart;
		int numSamplesObserved;
		boolean isGcNull, isXNull, isYNull, isBafNull, isLrrNull, isGenotypeNull, isNegativeXYAllowed;

		fingerprint = -1;
		numSamplesProj = allSampsInProj.length;
		result = new MarkerData[markersIndicesInFile.length];
		parameterReadBuffer = new byte[TransposeData.MARKERDATA_PARAMETER_TOTAL_LEN];
		try {
			file = new RandomAccessFile(currentMarkFilename, "r");
			file.read(parameterReadBuffer);
			// numMarkers = Compression.bytesToInt(parameterReadBuffer,
			// TransposeData.MARKERDATA_NUMMARKS_START);
			nullStatus = parameterReadBuffer[TransposeData.MARKERDATA_NULLSTATUS_START];
			// bytesPerSampleMarker = (byte) (Compression.BYTES_PER_SAMPLE_MARKER - (nullStatus & 0x01) -
			// (nullStatus >>1 & 0x01) - (nullStatus >>2 & 0x01) - (nullStatus >>3 & 0x01) - (nullStatus
			// >>4 & 0x01) - (nullStatus >>5 & 0x01) - (nullStatus >>6 & 0x01));
			bytesPerSampleMarker = Sample.getNBytesPerSampleMarker(nullStatus);
			numBytesPerMarker = bytesPerSampleMarker * numSamplesProj;
			readBuffer = new byte[markersIndicesInFile.length][numBytesPerMarker];
			numSamplesObserved = Compression.bytesToInt(parameterReadBuffer,
																									TransposeData.MARKERDATA_NUMSAMPLES_START);
			if (numSamplesObserved != numSamplesProj) {
				log.reportError("Error - mismatched number of samples between sample list (n="
													+ numSamplesProj + ") and file '" + currentMarkFilename + "' (n="
												+ numSamplesObserved + ")");
				System.exit(1);
			}
			fingerprint = Compression.bytesToLong(parameterReadBuffer,
																						TransposeData.MARKERDATA_FINGERPRINT_START);
			if (fingerprint != sampleFingerprint) {
				log.reportError("Error - mismatched sample fingerprints between sample list and file '"
												+ currentMarkFilename + "'");
				System.exit(1);
			}

			numBytesMarkernamesSection = Compression.bytesToInt(parameterReadBuffer,
																													TransposeData.MARKERDATA_MARKERNAMELEN_START);
			// TODO to optimize here. Adjacent markers can be read in at once.
			for (int i = 0; i < markersIndicesInFile.length; i++) {
				seekLocation = (long) TransposeData.MARKERDATA_PARAMETER_TOTAL_LEN
													+ (long) numBytesMarkernamesSection
												+ markersIndicesInFile[i] * (long) numBytesPerMarker;
				file.seek(seekLocation);
				file.read(readBuffer[i]);
			}

			file.close();
		} catch (FileNotFoundException e) {
			log.reportError("Error - could not find RAF marker file '" + currentMarkFilename + "'");
			log.reportException(e);
		} catch (IOException ioe) {
			log.reportError("Error reading RAF marker file '" + currentMarkFilename + "'");
			log.reportException(ioe);
		}

		isGcNull = Sample.isGcNull(nullStatus);
		isXNull = Sample.isXNull(nullStatus);
		isYNull = Sample.isYNull(nullStatus);
		isBafNull = Sample.isBafNull(nullStatus);
		isLrrNull = Sample.isLrrNull(nullStatus);
		isGenotypeNull = Sample.isAbAndForwardGenotypeNull(nullStatus);
		isNegativeXYAllowed = Sample.isNegativeXOrYAllowed(nullStatus);

		for (int i = 0; i < markersIndicesInFile.length; i++) {
			indexReadBuffer = 0;

			indexStart = 0;
			indexReadBuffer = indexStart;
			// time = new Date().getTime();
			if (!isGcNull) {
				if (loadGC) {
					gcs = new float[numSamplesProj];
					for (int j = 0; j < numSamplesProj; j++) {
						gcs[j] = Compression.gcBafDecompress(new byte[] {	readBuffer[i][indexReadBuffer],
																															readBuffer[i][indexReadBuffer + 1]});
						indexReadBuffer += bytesPerSampleMarker;
					}
				}
				indexStart += 2;
			}
			indexReadBuffer = indexStart;
			if (!isXNull) {
				if (loadXY) {
					xs = new float[numSamplesProj];
					for (int j = 0; j < numSamplesProj; j++) {
						if (isNegativeXYAllowed) {
							xs[j] =
										Compression.xyDecompressAllowNegative(new byte[] {readBuffer[i][indexReadBuffer],
																																			readBuffer[i][indexReadBuffer
																																										+ 1]});
						} else {
							xs[j] =
										Compression.xyDecompressPositiveOnly(new byte[] {	readBuffer[i][indexReadBuffer],
																																			readBuffer[i][indexReadBuffer
																																										+ 1]});
						}
						if (xs[j] == Compression.REDUCED_PRECISION_XY_OUT_OF_RANGE_FLAG_FLOAT) {
							// xs[j] = outOfRangeValues.get(sampleName+"\t"+allMarkersProj[j]+"\tx");
							Float f = outOfRangeValues.get(markersIndicesInProj[i]	+ "\t" + allSampsInProj[j]
																							+ "\tx");
							if (f == null) {
								// System.out.println(outOfRangeValues.keySet());
								if (allSampsInProj == null	|| allMarkersInProj == null
										|| markersIndicesInProj == null) {
									log.reportError("Error - value for X is an outlier for the marker that is at index "
																		+ i
																	+ " in the targetMarkers (if subset) or in all markers (if all), for the sample with the index of "
																	+ j
																	+ "; however, the original value is not present in the outlier hash.  Check if outliers.ser exists and re-run.");
								} else {
									log.reportError("Error - value for X is an outlier [Marker: "
																		+ allMarkersInProj[markersIndicesInProj[i]] + ", Sample: "
																	+ allSampsInProj[j]
																	+ "], original value not present in outlier hash.  Check if outliers.ser exists and re-run.");
								}
								return null; // TODO better return value than null?
							}
							xs[j] = f.floatValue();
						}
						indexReadBuffer += bytesPerSampleMarker;
					}
				}
				indexStart += 2;
			}
			indexReadBuffer = indexStart;
			if (!isYNull) {
				if (loadXY) {
					ys = new float[numSamplesProj];
					for (int j = 0; j < numSamplesProj; j++) {
						if (isNegativeXYAllowed) {
							ys[j] =
										Compression.xyDecompressAllowNegative(new byte[] {readBuffer[i][indexReadBuffer],
																																			readBuffer[i][indexReadBuffer
																																										+ 1]});
						} else {
							ys[j] =
										Compression.xyDecompressPositiveOnly(new byte[] {	readBuffer[i][indexReadBuffer],
																																			readBuffer[i][indexReadBuffer
																																										+ 1]});
						}
						if (ys[j] == Compression.REDUCED_PRECISION_XY_OUT_OF_RANGE_FLAG_FLOAT) {
							// ys[j] = outOfRangeValues.get(sampleName+"\t"+allMarkersProj[j]+"\ty");
							Float f = outOfRangeValues.get(markersIndicesInProj[i]	+ "\t" + allSampsInProj[j]
																							+ "\ty");
							if (f == null) {
								if (allSampsInProj == null	|| allMarkersInProj == null
										|| markersIndicesInProj == null) {
									log.reportError("Error - value for Y is an outlier for the marker that is at index "
																		+ i
																	+ " in the targetMarkers (if subset) or in all markers (if all), for the sample with the index of "
																	+ j
																	+ "; however, the original value is not present in the outlier hash.  Check if outliers.ser exists and re-run.");
								} else {
									log.reportError("Error - value for Y is an outlier [Marker: "
																		+ allMarkersInProj[markersIndicesInProj[i]] + ", Sample: "
																	+ allSampsInProj[j]
																	+ "], original value not present in outlier hash.  Check if outliers.ser exists and re-run.");
								}
								return null; // TODO better return value than null?
							}
							ys[j] = f.floatValue();
						}
						indexReadBuffer += bytesPerSampleMarker;
					}

				}
				indexStart += 2;
			}
			indexReadBuffer = indexStart;
			if (!isBafNull) {
				if (loadBAF) {
					bafs = new float[numSamplesProj];
					for (int j = 0; j < numSamplesProj; j++) {
						bafs[j] = Compression.gcBafDecompress(new byte[] {readBuffer[i][indexReadBuffer],
																															readBuffer[i][indexReadBuffer + 1]});
						indexReadBuffer += bytesPerSampleMarker;
					}
				}
				indexStart += 2;
			}
			indexReadBuffer = indexStart;
			if (!isLrrNull) {
				if (loadLRR) {
					lrrs = new float[numSamplesProj];
					for (int j = 0; j < numSamplesProj; j++) {
						lrrs[j] = Compression.lrrDecompress(new byte[] {readBuffer[i][indexReadBuffer],
																														readBuffer[i][indexReadBuffer + 1],
																														readBuffer[i][indexReadBuffer + 2]});
						if (lrrs[j] == Compression.REDUCED_PRECISION_LRR_OUT_OF_RANGE_LRR_FLAG_FLOAT) {
							// lrrs[j] = outOfRangeValues.get(sampleName+"\t"+allMarkersProj[j]+"\tlrr");
							Float f = outOfRangeValues.get(markersIndicesInProj[i]	+ "\t" + allSampsInProj[j]
																							+ "\tlrr");
							if (f == null) {
								if (allSampsInProj == null	|| allMarkersInProj == null
										|| markersIndicesInProj == null) {
									log.reportError("Error - value for LRR is an outlier for the marker that is at index "
																		+ i
																	+ " in the targetMarkers (if subset) or in all markers (if all), for the sample with the index of "
																	+ j
																	+ "; however, the original value is not present in the outlier hash.  Check if outliers.ser exists and re-run.");
								} else {
									log.reportError("Error - value for LRR is an outlier [Marker: "
																		+ allMarkersInProj[markersIndicesInProj[i]] + ", Sample: "
																	+ allSampsInProj[j]
																	+ "], original value not present in outlier hash.  Check if outliers.ser exists and re-run.");
								}
								return null; // TODO better return value than null?
							}
							lrrs[j] = f.floatValue();
						}
						indexReadBuffer += bytesPerSampleMarker;
					}
				}
				indexStart += 3;
			}
			indexReadBuffer = indexStart;
			if (!isGenotypeNull && loadAbGenotype) {
				abGenotypes = new byte[numSamplesProj];
				forwardGenotypes = new byte[numSamplesProj];
				for (int j = 0; j < numSamplesProj; j++) {
					genotypeTmp = Compression.genotypeDecompress(readBuffer[i][indexReadBuffer]);
					abGenotypes[j] = genotypeTmp[0];
					forwardGenotypes[j] = genotypeTmp[1];
					indexReadBuffer += bytesPerSampleMarker;
				}
			}
			result[i] = new MarkerData(
																	allMarkersInProj == null	? null
																														: allMarkersInProj[markersIndicesInProj[i]],
																	allChrsInProj == null	? (byte) -1
																												: allChrsInProj[markersIndicesInProj[i]],
																	allPosInProj == null ? -1 : allPosInProj[markersIndicesInProj[i]],
																	fingerprint, gcs, null, null, xs, ys, null, null, bafs, lrrs,
																	abGenotypes, forwardGenotypes);
		}
		return result;
	}

	public static byte[][] loadFromMarkerDataRafWithoutDecompress(String markerDataRafFileName,
																																int[] indicesOfMarkersInTheFileToLoad,
																																int indexOfFirstSampleToLoad,
																																int numOfSamplesToLoad,
																																long fingerprintShouldBe,
																																Logger log) {
		byte[][] result = null;
		RandomAccessFile file;
		int numBytesPerMarker;
		long seekLocation;
		byte[] parameterReadBuffer;
		long fingerprintActual;
		int numBytesMarkernamesSection;
		int numSamplesProj;
		byte numBytesPerSampleMarker;
		int numMarkersCurrentFile;

		parameterReadBuffer = new byte[TransposeData.MARKERDATA_PARAMETER_TOTAL_LEN];

		try {
			file = new RandomAccessFile(markerDataRafFileName, "r");
			file.read(parameterReadBuffer);
			numSamplesProj = Compression.bytesToInt(parameterReadBuffer,
																							TransposeData.MARKERDATA_NUMSAMPLES_START);
			numMarkersCurrentFile = Compression.bytesToInt(	parameterReadBuffer,
																											TransposeData.MARKERDATA_NUMMARKERS_START);
			fingerprintActual = Compression.bytesToLong(parameterReadBuffer,
																									TransposeData.MARKERDATA_FINGERPRINT_START);
			if (fingerprintActual != fingerprintShouldBe) {
				log.reportError("Error - mismatched sample fingerprints between sample list and file '"
												+ markerDataRafFileName + "'");
				System.exit(1);
			}
			numBytesPerSampleMarker =
															Sample.getNBytesPerSampleMarker(parameterReadBuffer[TransposeData.MARKERDATA_NULLSTATUS_START]);
			numBytesPerMarker = numBytesPerSampleMarker * numSamplesProj;
			if (numOfSamplesToLoad <= 0 || numOfSamplesToLoad > numSamplesProj) {
				numOfSamplesToLoad = numSamplesProj;
			}
			numBytesMarkernamesSection = Compression.bytesToInt(parameterReadBuffer,
																													TransposeData.MARKERDATA_MARKERNAMELEN_START);
			if (indicesOfMarkersInTheFileToLoad != null) {
				Arrays.sort(indicesOfMarkersInTheFileToLoad);
				result = new byte[indicesOfMarkersInTheFileToLoad.length][numBytesPerSampleMarker
																																	* numOfSamplesToLoad];
				for (int i = 0; i < indicesOfMarkersInTheFileToLoad.length; i++) {
					seekLocation = (long) TransposeData.MARKERDATA_PARAMETER_TOTAL_LEN
														+ (long) numBytesMarkernamesSection
													+ indicesOfMarkersInTheFileToLoad[i] * (long) numBytesPerMarker
													+ (long) indexOfFirstSampleToLoad * numBytesPerSampleMarker;
					file.seek(seekLocation);
					file.read(result[i]);
				}
			} else {
				result = new byte[numMarkersCurrentFile][numBytesPerSampleMarker * numOfSamplesToLoad];
				for (int i = 0; i < numMarkersCurrentFile; i++) {
					seekLocation = (long) TransposeData.MARKERDATA_PARAMETER_TOTAL_LEN
														+ (long) numBytesMarkernamesSection + i * (long) numBytesPerMarker
													+ (long) indexOfFirstSampleToLoad * numBytesPerSampleMarker;
					file.seek(seekLocation);
					file.read(result[i]);
				}
			}
			file.close();
		} catch (FileNotFoundException e) {
			log.reportError("Error - could not find RAF marker file '" + markerDataRafFileName + "'");
			log.reportException(e);
		} catch (IOException e) {
			log.reportError("Error reading RAF marker file '" + markerDataRafFileName + "'");
			log.reportException(e);
		}

		return result;
	}

	public static MarkerDataLoader loadMarkerDataFromListInSeparateThread(Project proj,
																																				String[] markerList) {
		MarkerDataLoader markerDataLoader;
		Thread thread2;
		int amountToLoadAtOnceInMB;

		proj.getLog().report("Marker data is loading in an independent thread.");
		// amountToLoadAtOnceInMB = proj.getInt(proj.MAX_MEMORY_USED_TO_LOAD_MARKER_DATA);
		amountToLoadAtOnceInMB = proj.MAX_MEMORY_USED_TO_LOAD_MARKER_DATA.getValue();
		markerDataLoader = new MarkerDataLoader(proj, markerList, amountToLoadAtOnceInMB);
		if (markerDataLoader.isKilled()) {
			return null;
		}
		markerDataLoader.initiate();
		thread2 = new Thread(markerDataLoader);
		thread2.start();
		markerDataLoader.registerThread(thread2);

		return markerDataLoader;
	}

	private void initiate() {
		initiated = true;
	}

	public boolean isKilled() {
		return killed;
	}

	private void registerThread(Thread thread) {
		this.thread = thread;
	}

	public static MarkerDataLoader loadMarkerDataFromListInSameThread(Project proj,
																																		String[] markerList) {
		MarkerDataLoader markerDataLoader;
		int amountToLoadAtOnceInMB;

		proj.getLog().report("Marker data is loading in the same thread.");
		// amountToLoadAtOnceInMB = proj.getInt(proj.MAX_MEMORY_USED_TO_LOAD_MARKER_DATA);
		amountToLoadAtOnceInMB = proj.MAX_MEMORY_USED_TO_LOAD_MARKER_DATA.getValue();
		markerDataLoader = new MarkerDataLoader(proj, markerList, amountToLoadAtOnceInMB);
		markerDataLoader.run();

		return markerDataLoader;
	}

	public static int getnumBytes_Mark(byte nullStatus, int numberOfSamples) {
		return numberOfSamples * Sample.getNBytesPerSampleMarker(nullStatus);
	}
}
