package cnv.manage;

import javax.swing.JOptionPane;
import javax.swing.SwingWorker;

import common.Array;
import common.Files;
import common.ext;

import cnv.filesys.*;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.util.*;

public class MarkerDataLoaderRunnable implements Runnable {
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
//	private int indexRequested;
//	private boolean loadFirstFileOnly;
	private boolean killed;
	
	public MarkerDataLoaderRunnable(Project proj, String[] markerNames) {
		this.proj = proj;
		this.markerNames = markerNames;
		this.chrs = new byte[markerNames.length];
		this.positions = new int[markerNames.length];
		markerData = new MarkerData[markerNames.length];
		loaded = new boolean[markerNames.length];

		String[] markerNamesProj = proj.getMarkerNames();
		byte[] chrsProj = proj.getMarkerSet().getChrs();
		int[] positionsProj = proj.getMarkerSet().getPositions();
		for (int i = 0; i<markerNames.length; i++) {
			for (int j = 0; j<markerNamesProj.length; j++) {
				if (markerNames[i].equals(markerNamesProj[j])) {
					chrs[i] = chrsProj[j];
					positions[i] = positionsProj[j];
					break;
				}
			}
		}
		Vector<String> missingMarkers, v;
		String[] line;
		
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
//				System.err.println("Error - could not find "+markerNames[0]+" in the lookup table");
				missingMarkers.add(markerNames[0]);
			}
		}
		if (missingMarkers.size() > 0) {
			JOptionPane.showMessageDialog(null, "Error - the following markers were not found in the MarkerSet: "+Array.toStr(Array.toStringArray(missingMarkers), " "), "Error", JOptionPane.ERROR_MESSAGE);
		}
		currentIndexBeingLoaded = 0;
		currentDirection = +1;
//		loadFirstFileOnly = true;
		killed = false;
	}
	
	public MarkerDataLoaderRunnable(Project proj, String[] markerNames, byte[] chrs, int[] positions) {
		this(proj, markerNames);
		this.chrs = chrs;
		this.positions = positions;
	}
	
//	public void setLoadOnlyFirstFile() {
//		loadFirstFileOnly = true;
//	}
	
//	public boolean hasIndexBeenLoaded(int index) {
//		return markerData[index] != null;
//	}
//	
//	public MarkerData[] getMarkerData() {
//		return markerData;
//	}
	
	public void requestIndexBeTheNextFileToLoad(int indexRequested) {
		if (indexRequested < currentIndexBeingLoaded) {
			currentDirection = -1;
		} else {
			currentDirection = +1;
		}
		currentIndexBeingLoaded = indexRequested;

//		fileRequested = markerLookup.get(markerNames[index]).split("\t")[0];
//		System.err.println("Next file will be "+fileRequested);
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
	
//	public void loadMarkerData(int markerIndex) {
////		try {
////			wait(5000);
////		} catch (InterruptedException e) {
////			// TODO Auto-generated catch block
////			e.printStackTrace();
////		}
//		markerData[markerIndex] = MarkerSet.loadFromList(proj, new String[] {markerNames[markerIndex]})[0];
//		loaded[markerIndex] = true;
//	}

////	public MarkerData[] loadData() {
//	public void loadData() {
//		MarkerData[] collection;
//		String[] line;
//		Vector<String> v;
//		String filename;
//		int[] indices;
//		long fingerprint;
//		int count;
//		
//		System.out.println("doingInBackground");
//		
//		fingerprint = proj.getSampleList().getFingerprint();
//
//		System.out.println("SwingWorker is running.");
//		long time = new Date().getTime();
//		count = 0;
//		while (filenames.size() > 0) {
////			if (currentIndexBeingLoaded < 0 || currentIndexBeingLoaded == markerNames.length) {
////				System.err.println("Error - before");
////			}
//			while (loaded[currentIndexBeingLoaded]) {
//				if (currentIndexBeingLoaded+1 == markerNames.length) {
//					currentDirection = -1;
//				}
//				if (currentIndexBeingLoaded-1 == 0) {
//					currentDirection = +1;
//				}
//				currentIndexBeingLoaded += currentDirection;
////				if (currentIndexBeingLoaded < 0 || currentIndexBeingLoaded == markerNames.length) {
////					System.err.println("Error - within");
////				}
//			}
//			filename = markerLookup.get(markerNames[currentIndexBeingLoaded]).split("\t")[0];
//			filenames.remove(filename);
////			System.err.println("Currently loading "+indexOfCurrentFile);
//			System.err.print(".");
//			if (!Files.exists(proj.getDir(Project.PLOT_DIRECTORY)+filename, proj.getJarStatus())) {
//				JOptionPane.showMessageDialog(null, "Error - could not load data from '"+proj.getDir(Project.PLOT_DIRECTORY)+filename+"'; because the file could not be found", "Error", JOptionPane.ERROR_MESSAGE);
////				return null;
//				return;
//			}
//			
//			//TODO
//			collection = MarkerDataCollection.load(proj.getDir(Project.PLOT_DIRECTORY)+filename, proj.getJarStatus()).getCollection();
//			v = hash.get(filename);
//			for (int j = 0; j<v.size() && !killed; j++) {
//				line = v.elementAt(j).split("[\\s]+");
//				indices = ext.indicesOfStr(line[0], markerNames, true, true);
//				
//				for (int k = 0; k<indices.length; k++) {
//					try {
//						markerData[indices[k]] = collection[Integer.parseInt(line[1])];
//						loaded[indices[k]] = true;
//						count++;
//					} catch (Exception e) {
//						System.err.println("Error - failed to load data for marker '"+line[0]+"' which is in collection "+line[1]+" may need to regenerate the markerLookup file");
//					}
//					
//					if (markerData[indices[k]].getFingerprint()!=fingerprint) {
//						System.err.println("Error - mismatched fingerprint after MarkerLookup");
//					}					
//				}
//			}
//
//			if (killed) {
//				filenames = new Hashtable<String, String>();
//			}
////			indexOfCurrentFile++;
//			
////			if (loadFirstFileOnly) {
////				cancel(true);
////				loadFirstFileOnly = false;
////			}
//
//		}
//		if (killed) {
//			System.out.println("Killed");
//			// need to clean up
//			markerNames = null;
//			loaded = null;
//			markerData = null;
//			System.gc();
//		} else {
//			System.out.println("SwingWorker has finished loading "+count+" additional markers to MarkerData[] in "+ ext.getTimeElapsed(time));
//		}
//		
////		return markerData;
//	}

	public void loadData() {
		MarkerData[] collection;
		String[] line;
		Vector<String> v;
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
		
		System.out.println("doingInBackground");
		
		fingerprint = proj.getSampleList().getFingerprint();

		System.out.println("SwingWorker is running.");
		markerSet = proj.getMarkerSet();
		markerNamesInProj = markerSet.getMarkerNames();
		chrsInProj = markerSet.getChrs();
		positionsInProj = markerSet.getPositions();
		samplesNamesProj = proj.getSamples();
		long time = new Date().getTime();
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
			System.out.println("loading marker Data file:"+filename);
			filenames.remove(filename);
			System.err.print(".");
			if (!Files.exists(proj.getDir(Project.PLOT_DIRECTORY)+filename, proj.getJarStatus())) {
				JOptionPane.showMessageDialog(null, "Error - could not load data from '"+proj.getDir(Project.PLOT_DIRECTORY)+filename+"'; because the file could not be found", "Error", JOptionPane.ERROR_MESSAGE);
				return;
			}
			
			v = hash.get(filename);
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

			for (int j = 0; j<v.size() && !killed; j++) {
				for (int k = 0; k<markerIndicesInFile.length; k++) {
					markerData[markerIndicesInSelection[k]] = collection[k];
					loaded[markerIndicesInSelection[k]] = true;
					count++;
					if (markerData[markerIndicesInSelection[k]].getFingerprint()!=fingerprint) {
						System.err.println("Error - mismatched fingerprint after MarkerLookup. Actual in MarkerData: " + markerData[markerIndicesInSelection[k]].getFingerprint() + ", while expecting: " + fingerprint);
					}					
				}
			}

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
			System.out.println("SwingWorker has finished loading "+count+" additional markers to MarkerData[] in "+ ext.getTimeElapsed(time));
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
	public static MarkerData[] loadFromRAF9(String[] markerNamesProj, byte[] chrsProj, int[] positionsProj, String[] samplesProj, String markerFilename, int[] markerIndeciesInProj, int[] markerIndeciesInFile, boolean loadGC, boolean loadXY, boolean loadBAF, boolean loadLRR, boolean loadAbGenotype) {
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
		result = new MarkerData[markerIndeciesInFile.length];
		parameterReadBuffer = new byte[TransposeData.MARKDATA_PARAMETER_TOTAL_LEN];
        try {
			file = new RandomAccessFile(markerFilename, "r");
			file.read(parameterReadBuffer);
//			numMarkersInThisFile = Compression.bytesToInt(new byte[] {parameterReadBuffer[4], parameterReadBuffer[5], parameterReadBuffer[6], parameterReadBuffer[7]});
			numMarkersInThisFile = Compression.bytesToInt(parameterReadBuffer, TransposeData.MARKDATA_NUMMARKS_START);
			nullStatus = parameterReadBuffer[TransposeData.MARKDATA_NULLSTATUS_START];
			bytesPerSampMark = (byte) (Compression.BYTES_PER_SAMPLE_MARKER_2 - (nullStatus & 0x01) - (nullStatus >>1 & 0x01) - (nullStatus >>2 & 0x01) - (nullStatus >>3 & 0x01) - (nullStatus >>4 & 0x01) - (nullStatus >>5 & 0x01) - (nullStatus >>6 & 0x01));
			numBytesPerMarker = bytesPerSampMark * numSamplesProj;
			readBuffer = new byte[markerIndeciesInProj.length][numBytesPerMarker];
//			fingerprint = Compression.bytesToLong(new byte[] {parameterReadBuffer[9], parameterReadBuffer[10], parameterReadBuffer[11], parameterReadBuffer[12], parameterReadBuffer[13], parameterReadBuffer[14], parameterReadBuffer[15], parameterReadBuffer[16]});
//			numBytesMarkernamesSection = Compression.bytesToInt(new byte[] {parameterReadBuffer[17], parameterReadBuffer[18], parameterReadBuffer[19], parameterReadBuffer[20]});
			fingerprint = Compression.bytesToLong(parameterReadBuffer, TransposeData.MARKDATA_FINGERPRINT_START);
//			fingerprint = Compression.bytesToLong(Compression.longToBytes(15351532491l), 0);
			numBytesMarkernamesSection = Compression.bytesToInt(parameterReadBuffer, TransposeData.MARKDATA_MARKNAMELEN_START);
			//TODO to optimize here. Adjacent markers can be read in at once.
	        for (int i=0; i<markerIndeciesInFile.length; i++) {
		        //if(indeciesInFile[i-1]+1==indeciesInFile[i]) {No need to seek}
		        seekLocation = (long)TransposeData.MARKDATA_PARAMETER_TOTAL_LEN + (long)numBytesMarkernamesSection + markerIndeciesInFile[i] * (long)numBytesPerMarker;
				file.seek(seekLocation);
				file.read(readBuffer[i]);
			}
	        
	        // Read in the Out of Range Value array
	        file.seek((long)TransposeData.MARKDATA_PARAMETER_TOTAL_LEN + (long)numBytesMarkernamesSection + (long) numMarkersInThisFile * (long)numBytesPerMarker);
	        System.out.println("number of markers in this file: "+numMarkersInThisFile);
			lengthOfOutOfRangeHashtable = file.readInt();
			if (lengthOfOutOfRangeHashtable>0) {
				outOfRangeValuesReadBuffer = new byte[lengthOfOutOfRangeHashtable];
				file.read(outOfRangeValuesReadBuffer);
				outOfRangeValues = (Hashtable<String, Float>) Compression.bytesToObj(outOfRangeValuesReadBuffer);
			}
			file.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

        for (int i=0; i<markerIndeciesInFile.length; i++) {
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
