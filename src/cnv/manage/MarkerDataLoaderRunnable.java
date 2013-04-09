package cnv.manage;

import javax.swing.JOptionPane;

import common.Array;
import common.Files;
import common.HashVec;
import common.ext;

import cnv.filesys.*;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.*;

public class MarkerDataLoaderRunnable implements Runnable {
	public static final int MAX_PER_CYCLE = 100;
	
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
	
	public MarkerDataLoaderRunnable(Project proj, String[] markerNames) {
		this.proj = proj;
		this.markerNames = markerNames;
		this.chrs = new byte[markerNames.length];
		this.positions = new int[markerNames.length];
		markerData = new MarkerData[markerNames.length];
		loaded = new boolean[markerNames.length];

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
		
		markerLookup = proj.getMarkerLookup();
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
				v.add(markerNames[i]+"\t"+line[1]);
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
	
	public MarkerDataLoaderRunnable(Project proj, String[] markerNames, byte[] chrs, int[] positions) {
		this(proj, markerNames);
		this.chrs = chrs;
		this.positions = positions;
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

	public void loadData() {
		MarkerData[] collection;
		String[] line;
		Vector<String> v, allRemaining;
		String filename;
		int[] markerIndicesInProj = null;
		int[] markerIndicesInFile;
		int[] markerIndicesInSelection;
		long fingerprint;
		int count;
		MarkerSet markerSet;
		String[] markerNamesInProj;
		byte[] chrsInProj;
		int[] positionsInProj;
		String[] samplesNamesProj;
		long time, start;
		
		fingerprint = proj.getSampleList().getFingerprint();

		System.out.println("Marker data is loading in an independent thread.");
		markerSet = proj.getMarkerSet();
		markerNamesInProj = markerSet.getMarkerNames();
		chrsInProj = markerSet.getChrs();
		positionsInProj = markerSet.getPositions();
		samplesNamesProj = proj.getSamples();
		start = time = new Date().getTime();
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
			if (!Files.exists(proj.getDir(Project.PLOT_DIRECTORY)+filename, proj.getJarStatus())) {
				JOptionPane.showMessageDialog(null, "Error - could not load data from '"+proj.getDir(Project.PLOT_DIRECTORY)+filename+"'; because the file could not be found", "Error", JOptionPane.ERROR_MESSAGE);
				return;
			}
			
			allRemaining = hash.get(filename);
			
			v = new Vector<String>();
			for (int i = 0; allRemaining.size() >0 && i < MAX_PER_CYCLE; i++) {
				v.addElement(allRemaining.remove(0));
			}
			System.out.println("Loaded up "+v.size()+" from marker Data file:"+filename+"; "+allRemaining.size()+" remaining for that file; "+ext.getTimeElapsed(time));
			time = new Date().getTime();
			if (allRemaining.size() == 0) {
				filenames.remove(filename);
			}
			markerIndicesInFile = new int[v.size()];
			markerIndicesInProj = new int[v.size()];	//TODO modified here to fix the bug
			markerIndicesInSelection = new int[v.size()];
			for (int j = 0; j<v.size() && !killed; j++) {
				line = v.elementAt(j).split("[\\s]+");
				markerIndicesInProj[j] = ext.indexOfStr(line[0], markerNamesInProj, true, true);	//modified here to fix the bug
				markerIndicesInFile[j] = Integer.parseInt(line[1]);
				markerIndicesInSelection[j] = ext.indexOfStr(line[0], markerNames, true, true);
			}

			collection = loadFromRAF9(markerNamesInProj, chrsInProj, positionsInProj, samplesNamesProj, proj.getDir(Project.PLOT_DIRECTORY)+filename, markerIndicesInProj, markerIndicesInFile);

//			for (int j = 0; j<v.size() && !killed; j++) {
				for (int k = 0; k<markerIndicesInFile.length && !killed; k++) {
					markerData[markerIndicesInSelection[k]] = collection[k];
					loaded[markerIndicesInSelection[k]] = true;
					count++;
					if (markerData[markerIndicesInSelection[k]].getFingerprint()!=fingerprint) {
						System.err.println("Error - mismatched fingerprint after MarkerLookup. Actual in MarkerData: " + markerData[markerIndicesInSelection[k]].getFingerprint() + ", while expecting: " + fingerprint);
					}					
				}
//			}

			if (killed) {
				filenames = new Hashtable<String, String>();
			}

		}
		if (killed) {
			System.out.println("Killed");
			// need to clean up
			markerNames = null;
			loaded = null;
			markerData = null;
			System.gc();
		} else {
			System.out.println("Independent thread has finished loading "+count+" markers to MarkerData[] in "+ ext.getTimeElapsed(time));
		}
		
	}

	public static MarkerData[] loadFromRAF9(String[] markerNamesProj, byte[] chrsProj, int[] positionsProj, String[] samplesProj, String markerFilename, int[] markerIndeciesInProj, int[] markerIndeciesInFile) {
		return loadFromRAF9(markerNamesProj, chrsProj, positionsProj, samplesProj, markerFilename, markerIndeciesInProj, markerIndeciesInFile, true, true, true, true, true);
	}


	/**
	 * Load MarkerData of selected markers from a single Random Access File of 12-byte (half-precision) format that has incorporated the marker names. Load marker data (RAF) approach 7b.
	 * @param proj
	 * @param markerNamesOfInterest
	 * @return
	 */
	@SuppressWarnings("unchecked")
	public static MarkerData[] loadFromRAF9(String[] markerNamesProj, byte[] chrsProj, int[] positionsProj, String[] samplesProj, String markerFilename, int[] markerIndeciesInProj, int[] markerIndicesInFile, boolean loadGC, boolean loadXY, boolean loadBAF, boolean loadLRR, boolean loadAbGenotype) {
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
//	    ByteBuffer readBuffer2;
        int indexReadBuffer;
        byte[] parameterReadBuffer;
        int numMarkersInThisFile;
        long fingerprint = 0;
        int numBytesMarkernamesSection;
        int lengthOfOutOfRangeHashtable;
        byte[] outOfRangeValuesReadBuffer;
        Hashtable<String, Float> outOfRangeValues = null;
        int numSamplesProj;
        byte nullStatus = 0;
        byte bytesPerSampMark = 0;
        int indexStart;

        numSamplesProj = samplesProj.length;
		result = new MarkerData[markerIndicesInFile.length];
		parameterReadBuffer = new byte[TransposeData.MARKDATA_PARAMETER_TOTAL_LEN];
        try {
			file = new RandomAccessFile(markerFilename, "r");
			file.read(parameterReadBuffer);
//			numMarkersInThisFile = Compression.bytesToInt(new byte[] {parameterReadBuffer[4], parameterReadBuffer[5], parameterReadBuffer[6], parameterReadBuffer[7]});
			numMarkersInThisFile = Compression.bytesToInt(parameterReadBuffer, TransposeData.MARKDATA_NUMMARKS_START);
			nullStatus = parameterReadBuffer[TransposeData.MARKDATA_NULLSTATUS_START];
			bytesPerSampMark = (byte) (Compression.BYTES_PER_SAMPLE_MARKER - (nullStatus & 0x01) - (nullStatus >>1 & 0x01) - (nullStatus >>2 & 0x01) - (nullStatus >>3 & 0x01) - (nullStatus >>4 & 0x01) - (nullStatus >>5 & 0x01) - (nullStatus >>6 & 0x01));
			numBytesPerMarker = bytesPerSampMark * numSamplesProj;
			readBuffer = new byte[markerIndeciesInProj.length][numBytesPerMarker];
//			fingerprint = Compression.bytesToLong(new byte[] {parameterReadBuffer[9], parameterReadBuffer[10], parameterReadBuffer[11], parameterReadBuffer[12], parameterReadBuffer[13], parameterReadBuffer[14], parameterReadBuffer[15], parameterReadBuffer[16]});
//			numBytesMarkernamesSection = Compression.bytesToInt(new byte[] {parameterReadBuffer[17], parameterReadBuffer[18], parameterReadBuffer[19], parameterReadBuffer[20]});
			fingerprint = Compression.bytesToLong(parameterReadBuffer, TransposeData.MARKDATA_FINGERPRINT_START);
//			fingerprint = Compression.bytesToLong(Compression.longToBytes(15351532491l), 0);
			numBytesMarkernamesSection = Compression.bytesToInt(parameterReadBuffer, TransposeData.MARKDATA_MARKNAMELEN_START);
			//TODO to optimize here. Adjacent markers can be read in at once.
	        for (int i=0; i<markerIndicesInFile.length; i++) {
		        //if(indeciesInFile[i-1]+1==indeciesInFile[i]) {No need to seek}
		        seekLocation = (long)TransposeData.MARKDATA_PARAMETER_TOTAL_LEN + (long)numBytesMarkernamesSection + markerIndicesInFile[i] * (long)numBytesPerMarker;
				file.seek(seekLocation);
				file.read(readBuffer[i]);
			}
	        
	        // Read in the Out of Range Value array
	        file.seek((long)TransposeData.MARKDATA_PARAMETER_TOTAL_LEN + (long)numBytesMarkernamesSection + (long) numMarkersInThisFile * (long)numBytesPerMarker);
//	        System.out.println("number of markers in this file: "+numMarkersInThisFile);
			lengthOfOutOfRangeHashtable = file.readInt();
			if (lengthOfOutOfRangeHashtable>0) {
				outOfRangeValuesReadBuffer = new byte[lengthOfOutOfRangeHashtable];
				file.read(outOfRangeValuesReadBuffer);
				outOfRangeValues = (Hashtable<String, Float>)Compression.bytesToObj(outOfRangeValuesReadBuffer);
			}
			file.close();
		} catch (FileNotFoundException e) {
			System.err.println("Error - could not find RAF marker file '"+markerFilename+"'");
			e.printStackTrace();
		} catch (IOException e) {
			System.err.println("Error reading RAF marker file '"+markerFilename+"'");
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			System.err.println("Error - could not convert the last section of the RAF marker file '"+markerFilename+"' into an hashtable of outliers");
			e.printStackTrace();
		}

        for (int i=0; i<markerIndicesInFile.length; i++) {
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
						gcs[j] = Compression.reducedPrecisionGcBafGetFloat2(new byte[] {readBuffer[i][indexReadBuffer], readBuffer[i][indexReadBuffer + 1]});
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
						xs[j] = Compression.reducedPrecisionXYGetFloat2(new byte[] {readBuffer[i][indexReadBuffer], readBuffer[i][indexReadBuffer + 1]});
						if (xs[j]==Compression.REDUCED_PRECISION_XY_OUT_OF_RANGE_FLOAT) {
//							xs[j] = outOfRangeValues.get(sampleName+"\t"+allMarkersProj[j]+"\tx");
							xs[j] = outOfRangeValues.get(markerIndeciesInProj[i] + "\t" + j + "\tx");
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
						ys[j] = Compression.reducedPrecisionXYGetFloat2(new byte[] {readBuffer[i][indexReadBuffer], readBuffer[i][indexReadBuffer + 1]});
						if (ys[j]==Compression.REDUCED_PRECISION_XY_OUT_OF_RANGE_FLOAT) {
//							ys[j] = outOfRangeValues.get(sampleName+"\t"+allMarkersProj[j]+"\ty");
							ys[j] = outOfRangeValues.get(markerIndeciesInProj[i] + "\t" + j + "\ty");
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
						bafs[j] = Compression.reducedPrecisionGcBafGetFloat2(new byte[] {readBuffer[i][indexReadBuffer], readBuffer[i][indexReadBuffer + 1]});
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
						lrrs[j] = Compression.reducedPrecisionLrrGetFloat2(new byte[] {readBuffer[i][indexReadBuffer], readBuffer[i][indexReadBuffer + 1], readBuffer[i][indexReadBuffer + 2]});
						if (lrrs[j] == Compression.REDUCED_PRECISION_LRR_OUT_OF_RANGE_LRR_FLOAT) {
//							lrrs[j] = outOfRangeValues.get(sampleName+"\t"+allMarkersProj[j]+"\tlrr");
							lrrs[j] = outOfRangeValues.get(markerIndeciesInProj[i] + "\t" + j + "\tlrr");
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
					genotypeTmp = Compression.reducedPrecisionGenotypeGetTypes(readBuffer[i][indexReadBuffer]);
					abGenotypes[j] = genotypeTmp[0];
					alleleMappings[j] = Sample.ALLELE_PAIRS[genotypeTmp[1]];
					indexReadBuffer += bytesPerSampMark;
				}
			}
	        result[i] = new MarkerData(markerNamesProj[markerIndeciesInProj[i]], chrsProj[markerIndeciesInProj[i]], positionsProj[markerIndeciesInProj[i]], fingerprint, gcs, new float[] {0}, new float[] {0}, xs, ys, new float[] {0}, new float[] {0}, bafs, lrrs, abGenotypes, alleleMappings);
        }
		return result;
	}

	public static MarkerDataLoaderRunnable loadMarkerDataFromList(Project proj, String[] markerList) {
		MarkerDataLoaderRunnable markerDataLoader;
		Thread thread2;
		
		markerDataLoader= new MarkerDataLoaderRunnable(proj, markerList);
		thread2 = new Thread(markerDataLoader);
		thread2.start();
		
		return markerDataLoader;
	}
	
	
	@Override
	public void run() {
//		launch third thread
		loadData();
//		while (true) {
//			loadData();
//			synchronized(this) {
//				while (true) wait();
//			}
//		}
	}
	
}
