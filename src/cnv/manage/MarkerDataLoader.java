package cnv.manage;

import javax.swing.JOptionPane;

import common.Array;
import common.Files;
import common.HashVec;
import common.ext;

import cnv.filesys.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.*;

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
	private boolean killed;
	private long sampleFingerprint;
	private int readAheadLimit;
	private int numberCurrentlyLoaded;
	private boolean plinkFormat;
	private String plinkFileRoot;
	private int[] plinkSampleIndices;

	public MarkerDataLoader(Project proj, String[] markerNames, int amountToLoadAtOnceInMB) {
		this(proj, markerNames, amountToLoadAtOnceInMB, null);	
	}
	
	public MarkerDataLoader(Project proj, String[] markerNames, int amountToLoadAtOnceInMB, String plinkFileRoot) {
		this.proj = proj;
		this.markerNames = markerNames;
		this.chrs = new byte[markerNames.length];
		this.positions = new int[markerNames.length];
		markerData = new MarkerData[markerNames.length];
		loaded = new boolean[markerNames.length];
		sampleFingerprint = proj.getSampleList().getFingerprint();
		numberCurrentlyLoaded = 0;

		Hashtable<String, Integer> markerHash;
		Vector<String> missingMarkers, v;
		String[] markerNamesProj;
		int[] positionsProj;
		byte[] chrsProj;
		String[] line;
		int index;

		markerNamesProj = proj.getMarkerNames();
		chrsProj = proj.getMarkerSet().getChrs();
		positionsProj = proj.getMarkerSet().getPositions();
		markerHash = HashVec.loadFileToHashIndex(markerNames);
		for (int i = 0; i<markerNamesProj.length; i++) {
			if (markerHash.containsKey(markerNamesProj[i])) {
				index = markerHash.get(markerNamesProj[i]);
				chrs[index] = chrsProj[i];
				positions[index] = positionsProj[i];
				break;
			}
		}

		if (plinkFileRoot != null) {
			this.plinkFormat = true;
			this.plinkFileRoot = plinkFileRoot;
			markerLookup = PlinkData.parseMarkerLookup(plinkFileRoot); // Note: markerLookup from PlinkData contains alleles, unlike its standard counterpart 
			plinkSampleIndices = PlinkData.parseSampleIndicesForProject(proj, plinkFileRoot, proj.getLog());
		} else {
			this.plinkFormat = false;
			this.plinkFileRoot = null;
			markerLookup = proj.getMarkerLookup();
		}

		if (amountToLoadAtOnceInMB <= 0) {
			amountToLoadAtOnceInMB = (int)((double)Runtime.getRuntime().maxMemory()/1024/1024*0.80);
			System.out.println("80% of max memory available ("+ext.prettyUpSize(Runtime.getRuntime().maxMemory(), 1)+") will be used by MarkerDataLoader: ");
		}
		readAheadLimit = -1;

		missingMarkers = new Vector<String>();
		filenames = new Hashtable<String, String>();
		hash = new Hashtable<String, Vector<String>>();
		for (int i = 0; i < markerNames.length; i++) {
			if (markerLookup.contains(markerNames[i])) {
				line = markerLookup.get(markerNames[i]).split("[\\s]+");
				if (hash.containsKey(line[0])) {
					v = hash.get(line[0]);
				} else {
					hash.put(line[0], v = new Vector<String>());
					filenames.put(line[0], "");
				}
				if (plinkFormat) {
					v.add(markerNames[i]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]);
				} else {
					v.add(markerNames[i]+"\t"+line[1]);
				}

				if (readAheadLimit == -1) {
					readAheadLimit = (int)Math.floor((double)amountToLoadAtOnceInMB *1024*1024 / (double)determineBytesPerMarkerData(proj.getDir(Project.MARKER_DATA_DIRECTORY)+line[0]));
					System.out.println("Read ahead limit was computed to be "+readAheadLimit+" markers at a time");
				}
			} else {
				missingMarkers.add(markerNames[0]);
			}
		}
		if (missingMarkers.size() > 0) {
			JOptionPane.showMessageDialog(null, "Error - the following markers were not found in the MarkerSet: "+Array.toStr(Array.toStringArray(missingMarkers), " "), "Error", JOptionPane.ERROR_MESSAGE);
		}
		currentIndexBeingLoaded = 0;
		currentDirection = +1;
		killed = false;
	}

//	public void setChrsPositions(byte[] chrs, int[] positions) {
//		this.chrs = chrs;
//		this.positions = positions;
//	}

	public static int determineBytesPerMarkerData(String filename) {
		RandomAccessFile file;
		byte[] parameterReadBuffer;
		byte nullStatus;
		int numSamplesObserved, bytesPerSampMark;

		parameterReadBuffer = new byte[TransposeData.MARKDATA_PARAMETER_TOTAL_LEN];
        try {
			file = new RandomAccessFile(filename, "r");
			file.read(parameterReadBuffer);
			nullStatus = parameterReadBuffer[TransposeData.MARKDATA_NULLSTATUS_START];
			bytesPerSampMark = Compression.BYTES_PER_SAMPLE_MARKER - (nullStatus & 0x01) - (nullStatus >>1 & 0x01) - (nullStatus >>2 & 0x01) - (nullStatus >>3 & 0x01) - (nullStatus >>4 & 0x01) - (nullStatus >>5 & 0x01) - (nullStatus >>6 & 0x01);
			numSamplesObserved = Compression.bytesToInt(parameterReadBuffer, TransposeData.MARKDATA_NUMSAMPS_START);
			file.close();
			return numSamplesObserved*bytesPerSampMark;
		} catch (FileNotFoundException e) {
			System.err.println("Error - could not find RAF marker file '"+filename+"'");
			e.printStackTrace();
		} catch (IOException e) {
			System.err.println("Error reading RAF marker file '"+filename+"'");
			e.printStackTrace();
		}

        return -1;
	}

//	public boolean hasIndexBeenLoaded(int index) {
//		return markerData[index] != null;
//	}

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
			System.out.println("Index "+index+" is not currently active; already null");
		} else {
			markerData[index] = null;
			numberCurrentlyLoaded--;
		}
	}

	public void kill() {
		killed = true;
	}

	public MarkerData getMarkerData(int markerIndex) {
		if (loaded[markerIndex]) {
			return markerData[markerIndex];
		} else {
			return null;
		}
	}

	public MarkerData requestMarkerData(int markerIndex) {
		MarkerData markerData;
		int count;

		count = 0;
		while((markerData = getMarkerData(markerIndex)) == null) {
			requestIndexBeTheNextFileToLoad(markerIndex);
			try {
				Thread.sleep(250);
				count++;
			} catch (InterruptedException ie) {
			}
			if (count > 8 && count % 8 == 0) {
				System.err.println("Error - have been waiting on markerDataLoader to load "+markerNames[markerIndex]+" for "+(count/4)+" seconds");
			}
		}

		return markerData;
	}

	@SuppressWarnings("unchecked")
	public void run() {
		MarkerData[] collection;
		String[] line;
		Vector<String> v, allRemaining;
		String filename;
		int[] marksIndicesInProj = null;
		int[] marksIndicesInFile;
		int[][] markerIndicesInSelection;
		long fingerprint;
		int count;
		MarkerSet markerSet;
		String[] allMarkersInProj;
		byte[] allChrsInProj;
		int[] allPosInProj;
		String[] allSampsInProj;
		long time;
		int maxPerCycle;
		Hashtable<String, Float> outlierHash;
		String[][] plinkMarkerAlleles;
		
		maxPerCycle = proj.getInt(Project.MAX_MARKERS_LOADED_PER_CYCLE);

		fingerprint = proj.getSampleList().getFingerprint();

		System.out.println("Marker data is loading in an independent thread.");
		markerSet = proj.getMarkerSet();
		allMarkersInProj = markerSet.getMarkerNames();
		allChrsInProj = markerSet.getChrs();
		allPosInProj = markerSet.getPositions();
		allSampsInProj = proj.getSamples();
		if (new File(proj.getDir(Project.MARKER_DATA_DIRECTORY) + "outliers.ser").exists()) {
			outlierHash = (Hashtable<String, Float>) Files.readSerial(proj.getDir(Project.MARKER_DATA_DIRECTORY) + "outliers.ser");
		} else {
			outlierHash = new Hashtable<String, Float>();
		}
		time = new Date().getTime();
		count = 0;
		while (filenames.size() > 0) {
			while (loaded[currentIndexBeingLoaded]) {
				if (currentIndexBeingLoaded+1 == markerNames.length) {
					currentDirection = -1;
				}
				if (currentIndexBeingLoaded-1 == 0) {
					currentDirection = +1;
				}
				currentIndexBeingLoaded += currentDirection;
			}
			filename = markerLookup.get(markerNames[currentIndexBeingLoaded]).split("\t")[0];
//			System.err.print(".");
			if (!Files.exists(proj.getDir(Project.MARKER_DATA_DIRECTORY)+filename, proj.getJarStatus())) {
				JOptionPane.showMessageDialog(null, "Error - could not load data from '"+proj.getDir(Project.MARKER_DATA_DIRECTORY)+filename+"'; because the file could not be found", "Error", JOptionPane.ERROR_MESSAGE);
				return;
			}

			allRemaining = hash.get(filename);

			v = new Vector<String>();
			for (int i = 0; allRemaining.size() >0 && i < maxPerCycle; i++) {
				v.addElement(allRemaining.remove(0));
			}
			System.out.println("Loaded up "+v.size()+" from MarkerData file:"+filename+"; "+allRemaining.size()+" remaining for that file; "+ext.getTimeElapsed(time));
			time = new Date().getTime();
			if (allRemaining.size() == 0) {
				filenames.remove(filename);
			}
			marksIndicesInFile = new int[v.size()];
			marksIndicesInProj = new int[v.size()];
			markerIndicesInSelection = new int[v.size()][];
			plinkMarkerAlleles = new String[v.size()][];
			for (int j = 0; j<v.size() && !killed; j++) {
				line = v.elementAt(j).split("[\\s]+");
				marksIndicesInProj[j] = ext.indexOfStr(line[0], allMarkersInProj, true, true);	//modified here to fix the bug
				marksIndicesInFile[j] = Integer.parseInt(line[1]);
				markerIndicesInSelection[j] = ext.indicesOfStr(line[0], markerNames, true, true);
				if (markerIndicesInSelection[j].length > 1) {
					System.out.println("FYI, marker "+line[0]+" was requested "+markerIndicesInSelection[j].length+" times");
				}
				if (plinkFormat) {
					plinkMarkerAlleles[j] = new String[] {line[2], line[3]};
				}
			}

			if (plinkFormat) {
				collection = PlinkData.loadBedUsingRAF(allMarkersInProj, allChrsInProj, allPosInProj, allSampsInProj, proj.getDir(Project.MARKER_DATA_DIRECTORY) + plinkFileRoot + ".bed", marksIndicesInProj, marksIndicesInFile, sampleFingerprint, plinkSampleIndices);
//				collection = PlinkData.loadPedUsingRAF(allMarkersInProj, allChrsInProj, allPosInProj, allSampsInProj, proj.getDir(Project.MARKER_DATA_DIRECTORY) + plinkFileRoot + ".bed", marksIndicesInProj, marksIndicesInFile, sampleFingerprint, plinkSampleIndices);
			} else {
				collection = loadFromRAF(allMarkersInProj, allChrsInProj, allPosInProj, allSampsInProj, proj.getDir(Project.MARKER_DATA_DIRECTORY)+filename, marksIndicesInProj, marksIndicesInFile, sampleFingerprint, outlierHash);
			}

			for (int k = 0; k < marksIndicesInProj.length && !killed; k++) {
				for (int i = 0; i < markerIndicesInSelection[k].length; i++) {
					markerData[markerIndicesInSelection[k][i]] = collection[k];
					loaded[markerIndicesInSelection[k][i]] = true;
					count++;
					numberCurrentlyLoaded++;
					if (markerData[markerIndicesInSelection[k][i]].getFingerprint()!=fingerprint) {
						System.err.println("Error - mismatched fingerprint after MarkerLookup. Actual in MarkerData: " + markerData[markerIndicesInSelection[k][i]].getFingerprint() + ", while expecting: " + fingerprint);
					}
				}
			}

			while (!killed && numberCurrentlyLoaded > readAheadLimit) {
				try {
					Thread.sleep(500);
				} catch (InterruptedException ie) {}
				System.out.println("Currently loading index "+currentIndexBeingLoaded+" but readAheadLimit has been reached ("+numberCurrentlyLoaded+" over "+readAheadLimit+"); so waiting");
			}
			if (killed) {
				filenames = new Hashtable<String, String>();
			}

		}
		if (killed) {
			System.out.println("Killed");
			markerNames = null;
			loaded = null;
			markerData = null;
			System.gc();
		} else {
			System.out.println("Independent thread has finished loading "+count+" markers to MarkerData[] in "+ ext.getTimeElapsed(time));
		}

	}

	public static MarkerData[] loadFromRAF(String[] markerNamesProj, byte[] chrsProj, int[] positionsProj, String[] samplesProj, String markerFilename, int[] markerIndeciesInProj, int[] markerIndeciesInFile, long sampleFingerprint, Hashtable<String, Float> outOfRangeValues) {
		return loadFromRAF(markerNamesProj, chrsProj, positionsProj, samplesProj, markerFilename, markerIndeciesInProj, markerIndeciesInFile, true, true, true, true, true, sampleFingerprint, outOfRangeValues);
	}


	/**
	 * Load MarkerData of selected markers from a single Random Access File of 12-byte (half-precision) format that has incorporated the marker names. Load marker data (RAF) approach 7b.
	 * @param proj
	 * @param markerNamesOfInterest
	 * @return
	 */
	public static MarkerData[] loadFromRAF(String[] allMarksInProj, byte[] allChrsInProj, int[] allPosInProj, String[] allSampsInProj, String markerFilename, int[] markersIndeciesInProj, int[] markersIndicesInFile, boolean loadGC, boolean loadXY, boolean loadBAF, boolean loadLRR, boolean loadAbGenotype, long sampleFingerprint, Hashtable<String, Float> outOfRangeValues) {
		MarkerData[] result;
		RandomAccessFile file;
        int numBytesPerMarker;
        float[] gcs = null;
        float[] xs = null;
        float[] ys = null;
        float[] bafs = null;
        float[] lrrs = null;
        byte[] abGenotypes;
        String[] alleleMappings;
        byte[] genotypeTmp;
        long seekLocation;
        byte[][] readBuffer = null;
        int indexReadBuffer;
        byte[] parameterReadBuffer;
        long fingerprint;
        int numBytesMarkernamesSection;
//        Hashtable<String, Float> outOfRangeValues = null;
        int numSamplesProj;
        byte nullStatus = 0;
        byte bytesPerSampMark = 0;
        int indexStart;
        int numSamplesObserved;

        fingerprint = -1;
        numSamplesProj = allSampsInProj.length;
		result = new MarkerData[markersIndicesInFile.length];
		parameterReadBuffer = new byte[TransposeData.MARKDATA_PARAMETER_TOTAL_LEN];
        try {
			file = new RandomAccessFile(markerFilename, "r");
			file.read(parameterReadBuffer);
//			numMarkers = Compression.bytesToInt(parameterReadBuffer, TransposeData.MARKDATA_NUMMARKS_START);
			nullStatus = parameterReadBuffer[TransposeData.MARKDATA_NULLSTATUS_START];
			bytesPerSampMark = (byte) (Compression.BYTES_PER_SAMPLE_MARKER - (nullStatus & 0x01) - (nullStatus >>1 & 0x01) - (nullStatus >>2 & 0x01) - (nullStatus >>3 & 0x01) - (nullStatus >>4 & 0x01) - (nullStatus >>5 & 0x01) - (nullStatus >>6 & 0x01));
			numBytesPerMarker = bytesPerSampMark * numSamplesProj;
			readBuffer = new byte[markersIndeciesInProj.length][numBytesPerMarker];
			numSamplesObserved = Compression.bytesToInt(parameterReadBuffer, TransposeData.MARKDATA_NUMSAMPS_START);
			if (numSamplesObserved != numSamplesProj) {
				System.err.println("Error - mismatched number of samples between sample list (n="+numSamplesProj+") and file '"+markerFilename+"' (n="+numSamplesObserved+")");
			}
			fingerprint = Compression.bytesToLong(parameterReadBuffer, TransposeData.MARKDATA_FINGERPRINT_START);
			if (fingerprint != sampleFingerprint) {
				System.err.println("Error - mismatched sample fingerprints between sample list and file '"+markerFilename+"'");
			}

			numBytesMarkernamesSection = Compression.bytesToInt(parameterReadBuffer, TransposeData.MARKDATA_MARKNAMELEN_START);
			//TODO to optimize here. Adjacent markers can be read in at once.
	        for (int i=0; i<markersIndicesInFile.length; i++) {
		        //if(indeciesInFile[i-1]+1==indeciesInFile[i]) {No need to seek}
		        seekLocation = (long)TransposeData.MARKDATA_PARAMETER_TOTAL_LEN + (long)numBytesMarkernamesSection + markersIndicesInFile[i] * (long)numBytesPerMarker;
				file.seek(seekLocation);
				file.read(readBuffer[i]);
			}

	        // TODO this is read every time, wouldn't it be faster to use the serialized version?
	        // Read in the Out of Range Value array
//	        file.seek((long)TransposeData.MARKDATA_PARAMETER_TOTAL_LEN + (long)numBytesMarkernamesSection + (long) numMarkersInThisFile * (long)numBytesPerMarker);
////	        System.out.println("number of markers in this file: "+numMarkersInThisFile);
//			lengthOfOutOfRangeHashtable = file.readInt();
//			if (lengthOfOutOfRangeHashtable>0) {
//				outOfRangeValuesReadBuffer = new byte[lengthOfOutOfRangeHashtable];
//				file.read(outOfRangeValuesReadBuffer);
//				outOfRangeValues = (Hashtable<String, Float>)Compression.bytesToObj(outOfRangeValuesReadBuffer);
//			}

			file.close();
		} catch (FileNotFoundException e) {
			System.err.println("Error - could not find RAF marker file '"+markerFilename+"'");
			e.printStackTrace();
		} catch (IOException e) {
			System.err.println("Error reading RAF marker file '"+markerFilename+"'");
			e.printStackTrace();
		}

        for (int i=0; i<markersIndicesInFile.length; i++) {
	        gcs = new float[numSamplesProj];
	        xs = new float[numSamplesProj];
	        ys = new float[numSamplesProj];
	        bafs = new float[numSamplesProj];
	        lrrs = new float[numSamplesProj];
	        abGenotypes = new byte[numSamplesProj] ;
	        alleleMappings = new String[numSamplesProj];
			indexReadBuffer = 0;

			indexStart = 0;
			indexReadBuffer = indexStart;
//			time = new Date().getTime();
			if (((nullStatus>>Sample.NULLSTATUS_GC_LOCATION) & 0x01) != 1) {
				if (loadGC) {
					gcs = new float[numSamplesProj];
					for (int j=0; j<numSamplesProj; j++) {
						gcs[j] = Compression.gcBafDecompress(new byte[] {readBuffer[i][indexReadBuffer], readBuffer[i][indexReadBuffer + 1]});
						indexReadBuffer += bytesPerSampMark;
					}
				}
				indexStart += 2;
			}
			indexReadBuffer = indexStart;
			if (((nullStatus>>Sample.NULLSTATUS_X_LOCATION) & 0x01) != 1) {
				if (loadXY) {
					xs = new float[numSamplesProj];
					for (int j=0; j<numSamplesProj; j++) {
						xs[j] = Compression.xyDecompress(new byte[] {readBuffer[i][indexReadBuffer], readBuffer[i][indexReadBuffer + 1]});
						if (xs[j]==Compression.REDUCED_PRECISION_XY_OUT_OF_RANGE_FLOAT) {
//							xs[j] = outOfRangeValues.get(sampleName+"\t"+allMarkersProj[j]+"\tx");
							xs[j] = outOfRangeValues.get(markersIndeciesInProj[i] + "\t" + allSampsInProj[j] + "\tx");
						}
						indexReadBuffer += bytesPerSampMark;
					}
				}
				indexStart += 2;
			}
			indexReadBuffer = indexStart;
			if (((nullStatus>>Sample.NULLSTATUS_Y_LOCATION) & 0x01) != 1) {
				if (loadXY) {
					ys = new float[numSamplesProj];
					for (int j=0; j<numSamplesProj; j++) {
						ys[j] = Compression.xyDecompress(new byte[] {readBuffer[i][indexReadBuffer], readBuffer[i][indexReadBuffer + 1]});
						if (ys[j]==Compression.REDUCED_PRECISION_XY_OUT_OF_RANGE_FLOAT) {
//							ys[j] = outOfRangeValues.get(sampleName+"\t"+allMarkersProj[j]+"\ty");
							ys[j] = outOfRangeValues.get(markersIndeciesInProj[i] + "\t" + allSampsInProj[j] + "\ty");
						}
						indexReadBuffer += bytesPerSampMark;
					}

				}
				indexStart += 2;
			}
			indexReadBuffer = indexStart;
			if (((nullStatus>>Sample.NULLSTATUS_BAF_LOCATION) & 0x01) != 1) {
				if (loadBAF) {
					bafs = new float[numSamplesProj];
					for (int j=0; j<numSamplesProj; j++) {
						bafs[j] = Compression.gcBafDecompress(new byte[] {readBuffer[i][indexReadBuffer], readBuffer[i][indexReadBuffer + 1]});
						indexReadBuffer += bytesPerSampMark;
					}
				}
				indexStart += 2;
			}
			indexReadBuffer = indexStart;
			if (((nullStatus>>Sample.NULLSTATUS_LRR_LOCATION) & 0x01) != 1) {
				if (loadLRR) {
					lrrs = new float[numSamplesProj];
					for (int j=0; j<numSamplesProj; j++) {
						lrrs[j] = Compression.lrrDecompress(new byte[] {readBuffer[i][indexReadBuffer], readBuffer[i][indexReadBuffer + 1], readBuffer[i][indexReadBuffer + 2]});
						if (lrrs[j] == Compression.REDUCED_PRECISION_LRR_OUT_OF_RANGE_LRR_FLOAT) {
//							lrrs[j] = outOfRangeValues.get(sampleName+"\t"+allMarkersProj[j]+"\tlrr");
							lrrs[j] = outOfRangeValues.get(markersIndeciesInProj[i] + "\t" + allSampsInProj[j] + "\tlrr");
						}
						indexReadBuffer += bytesPerSampMark;
					}
				}
				indexStart += 3;
			}
			indexReadBuffer = indexStart;
			if ((((nullStatus>>Sample.NULLSTATUS_ABGENOTYPE_LOCATION) & 0x01) != 1 || ((nullStatus>>Sample.NULLSTATUS_FOWARDGENOTYPE_LOCATION) & 0x01) != 1) && loadAbGenotype) {
				abGenotypes = new byte[numSamplesProj];
				for (int j=0; j<numSamplesProj; j++) {
					genotypeTmp = Compression.genotypeDecompress(readBuffer[i][indexReadBuffer]);
					abGenotypes[j] = genotypeTmp[0];
					if (genotypeTmp[1] >= Sample.ALLELE_PAIRS.length) {
						System.err.println("Error - invalid allelePair designation ("+genotypeTmp[1]+") as there are only "+Sample.ALLELE_PAIRS.length+" that are defined");
						alleleMappings[j] = "UU";
					} else {
						alleleMappings[j] = Sample.ALLELE_PAIRS[genotypeTmp[1]];
					}
					indexReadBuffer += bytesPerSampMark;
				}
			}
	        result[i] = new MarkerData(allMarksInProj[markersIndeciesInProj[i]], allChrsInProj[markersIndeciesInProj[i]], allPosInProj[markersIndeciesInProj[i]], fingerprint, gcs, null, null, xs, ys, null, null, bafs, lrrs, abGenotypes, alleleMappings);
        }
		return result;
	}

	public static MarkerDataLoader loadMarkerDataFromList(Project proj, String[] markerList) {
		MarkerDataLoader markerDataLoader;
		Thread thread2;
		int amountToLoadAtOnceInMB;
		
		amountToLoadAtOnceInMB = proj.getInt(Project.MAX_MEMORY_USED_TO_LOAD_MARKER_DATA);
		markerDataLoader= new MarkerDataLoader(proj, markerList, amountToLoadAtOnceInMB);
		thread2 = new Thread(markerDataLoader);
		thread2.start();

		return markerDataLoader;
	}
}
