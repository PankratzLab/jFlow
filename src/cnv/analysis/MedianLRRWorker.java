package cnv.analysis;

import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.Hashtable;
import java.util.Vector;
import java.util.concurrent.CancellationException;
import java.util.concurrent.ExecutionException;

import javax.swing.JOptionPane;
import javax.swing.JProgressBar;
import javax.swing.SwingWorker;

import cnv.filesys.MarkerData;
import cnv.filesys.MarkerSet;
import cnv.filesys.MeanLRRset;
import cnv.filesys.Project;
import cnv.filesys.Sample;
import cnv.filesys.SampleList;
import cnv.manage.MarkerDataLoader;
import cnv.manage.Transforms;
import cnv.var.CNVariant;
import common.Array;
import common.Files;
import common.HashVec;
import common.IntVector;
import common.Logger;
import common.Positions;
import common.ext;
import filesys.Segment;

public class MedianLRRWorker extends SwingWorker<Integer, Integer> {
	private static final String[] MEDIAN_EXTENSIONS = { ".mlrr", ".log" };
	private static final String MARKER_REGION_PREFIX = "probeset_id";
	private static final String MARKER_REGION_DELIMITER = ";";
	private static final String MARKER_REGION_REGEX = "[%0]";
	private static final int MARKER_REGION_START_OF_MARKERS = 2;
	private static final double[] QUANTILES = { 0.5 };
	private static final String[] MEDIAN_WORKER_JOBS = { "Parsing and intitializing Regions", "Computing Median Log R Ratios for " + MARKER_REGION_REGEX, "Waiting for data to Load for Region " + MARKER_REGION_REGEX, "Creating Output Files" };
	public static final String[] CLASSES_TO_DUMP = { "IID" };
	private Project proj;
	private String[] input;
	private String outputBase;
	private int transformationType;
	private int scope;
	private Logger computelog;
	private String[] samples, markerNames;
	private byte[] chrs;
	private int[] positions;
	private Hashtable<String, String> hash;
	private MarkerSet markerSet;
	private SampleList sampleList;
	private int[] processTracker;
	private boolean[] transChrs;

	private JProgressBar progressBar;

	public MedianLRRWorker(Project proj, String[] input, int transformationType, int scope, String outputBase, JProgressBar jProgressBar) {
		this.proj = proj;
		this.input = input;
		this.outputBase = outputBase;
		this.transformationType = transformationType;
		this.scope = scope;
		this.progressBar = jProgressBar;
		addPropertyChangeListener(new PropertyChangeListener() {
			public void propertyChange(PropertyChangeEvent evt) {
				if ("progress".equals(evt.getPropertyName())) {
					progressBar.setValue((Integer) evt.getNewValue());
				}
			}
		});

		this.computelog = new Logger(proj.getProjectDir() + outputBase + ".log");
		this.markerSet = proj.getMarkerSet();
		this.markerNames = markerSet.getMarkerNames();
		this.chrs = markerSet.getChrs();
		this.positions = markerSet.getPositions();
		this.hash = proj.getFilteredHash();
		this.sampleList = proj.getSampleList();
		this.samples = sampleList.getSamples();
		// total ,processed;
		this.processTracker = new int[2];
		this.transChrs = Array.booleanArray(27, false);
		progressBar.setValue(0);
		setProgress(0);
	}

	protected Integer doInBackground() throws Exception {
		System.out.println(proj.getProjectDir() + outputBase);
		newJob(MEDIAN_WORKER_JOBS[0]);
		MarkerRegion[] markerRegions = parseRegions();
		process(processTracker[0]);
		if (processTracker[0] < 1) {
			String Error = "Error - dataset did not contain markers in the Region(s) of interest";
			computelog.reportError(Error);
			warnAndCancel(Error);
		} else {
			process(0);
			float[][] regionMedianValues;
			if (transformationType == 0) {
				newJob(ext.replaceAllWith(MEDIAN_WORKER_JOBS[1], "[%" + 0 + "]", Transforms.TRANFORMATIONS[transformationType]));
				assignMarkerProgress(markerRegions.length);
				regionMedianValues = getRawValueMedianForRegions(markerRegions);
				newJob(MEDIAN_WORKER_JOBS[3]);
				process(processTracker[0] - 1);
				printResults(regionMedianValues, markerRegions);
			} else {
				JOptionPane.showMessageDialog(null, "Sorry, this method is not yet supported...");
			}
			Thread.sleep(1000);
		}
		return 42;
	}

	// total,processed,updateat,tracker,scale;

	private void assignMarkerProgress(int extraFactor) {
		progressBar.setMaximum(processTracker[0]);

	}

	private void assignSampleProgress() {
		progressBar.setMaximum(samples.length);
	}

	// update progressBar
	protected void process(Integer chunk) {
		progressBar.setValue(chunk);
	}

	// set ProgressBar Name
	protected void newJob(String job) {
		progressBar.setString(job);
	}

	protected void warnAndCancel(String message) {
		JOptionPane.showMessageDialog(null, message);
		this.cancel(true);
	}

	protected void done() {
		try {
			get();
			JOptionPane.showMessageDialog(null, "Log R Ratio Summarization Complete");
			setProgress(100);
			progressBar.setValue(100);
			progressBar.setStringPainted(false);
			progressBar.setVisible(false);
		} catch (ExecutionException e) {
			computelog.reportError("Error - Could not Compute Median Log R Ratio Values");
			e.printStackTrace();
		} catch (InterruptedException e) {
			computelog.reportError("Error - Median Log R Ratio Computation was Interupted ");
		} catch (CancellationException cce) {
			computelog.reportError("Error - Cancelling Median Log R Ratio Computation ");

		}
	}

	// First Get Markers For Region if Region, else region =markers input;

	private MarkerRegion[] parseRegions() {
		MarkerRegion[] regions = new MarkerRegion[input.length];
		if (input != null) {
			for (int i = 0; i < input.length; i++) {
				process((int) Math.round((float) (i + 1 / input.length)));
				if (input[i].split(MARKER_REGION_DELIMITER)[0].equals(MARKER_REGION_PREFIX)) {
					String[] regionMarkers = input[i].split(MARKER_REGION_DELIMITER);
					if (regionMarkers.length <= 1) {
						String Error = "Error - markers must be " + MARKER_REGION_DELIMITER + " delimited for Marker Region " + input[i];
						computelog.reportError(Error);
						warnAndCancel(Error);
					}
					String regionName = regionMarkers[1];
					regions[i] = new MarkerRegion(regionName, i, 1);
					regions[i].addMarkers(regionMarkers);

				} else if (input[i].startsWith("chr") && checkUCSCRegions(input[i])) {
					regions[i] = new MarkerRegion(input[i], i, 0);
					regions[i].addUCSCMarkers(input[i]);
				} else {
					String Error = "Error - Improper Formatting for Input " + input[i] + "\n Input must one per line of \n(1) UCSC Formatted\n or\n(2) " + MARKER_REGION_PREFIX + MARKER_REGION_DELIMITER + "(Your Region Name)" + MARKER_REGION_DELIMITER + "marker names(" + MARKER_REGION_DELIMITER + ") delimited";
					computelog.reportError(Error);
					warnAndCancel(Error);
				}
			}
		} else {
			String Error = "Error - Zero Regions Were Entered";
			computelog.reportError(Error);
			warnAndCancel(Error);
		}
		return regions;
	}

	private float[][] getRawValueMedianForRegions(MarkerRegion[] markerRegions) {
		float[][] rawMediansRegions = new float[markerRegions.length][samples.length];
		for (int i = 0; i < markerRegions.length; i++) {
			rawMediansRegions[i] = getRawValueMedianForRegion(markerRegions[i]);
		}
		return rawMediansRegions;
	}

	// For each region load lrrs from markerData if and assign medians depending on anlayis

	private void printResults(float[][] regionMedianValues, MarkerRegion[] markerRegions) {
		printRegionMarkers(markerRegions);
		printMedianLRRs(regionMedianValues, markerRegions);
	}

	private void printMedianLRRs(float[][] regionMedianValues, MarkerRegion[] markerRegions) {
		String output = proj.getProjectDir() + "LRR_MEDIAN_" + ext.rootOf(outputBase) + ".xln";
		Hashtable<String, String> hashSamps = HashVec.loadFileToHashString(proj.getFilename(Project.SAMPLE_DATA_FILENAME), "DNA", CLASSES_TO_DUMP, "\t");
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(output));
			writer.print("Sample");
			for (int i = 0; i < CLASSES_TO_DUMP.length; i++) {
				writer.print("\t" + CLASSES_TO_DUMP[i].substring(CLASSES_TO_DUMP[i].lastIndexOf("=") + 1));
			}
			for (int i = 0; i < markerRegions.length; i++) {
				writer.print("\t" + markerRegions[i].getRegionName() + "_" + ext.replaceWithLinuxSafeCharacters(outputBase, false));
			}
			writer.println();
			for (int i = 0; i < samples.length; i++) {
				writer.print(samples[i] + "\t" + (hashSamps.containsKey(samples[i]) ? hashSamps.get(samples[i]) : Array.stringArray(CLASSES_TO_DUMP.length, ".")));
				for (int j = 0; j < markerRegions.length; j++) {
					writer.print("\t" + regionMedianValues[j][i]);
				}
				writer.println();
			}
			writer.close();
		} catch (Exception e) {
			computelog.reportError("Error writing to Medain Log R Ratios to " + output);
			e.printStackTrace();
		}
	}

	private void printRegionMarkers(MarkerRegion[] markerRegions) {
		String output = proj.getProjectDir() + "MarkersIn_" + ext.rootOf(outputBase) + ".xln";
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(output));
			for (int i = 0; i < markerRegions.length; i++) {
				String[] regionMarkers = markerRegions[i].returnMarkers();
				if (markerRegions[i].getRegionType() == 0) {
					Segment seg = new Segment(markerRegions[i].getRegionName());
					writer.println(seg.getUCSClocation() + "\t" + regionMarkers.length + "\t" + seg.getChr() + "\t" + seg.getStart() + "\t" + seg.getStop() + "\t" + Array.toStr(regionMarkers));
				} else if (markerRegions[i].getRegionType() == 1) {
					writer.println(markerRegions[i].getRegionName() + "\t" + regionMarkers.length + "\tNA\tNA\tNA\t" + Array.toStr(regionMarkers));
				} else {
					computelog.reportError("Error - Unknonwn Region Type " + markerRegions[i].getRegionType());
				}
			}
			writer.close();
		} catch (Exception e) {
			computelog.reportError("Error writing the list of marker names within the regions to " + output);
			e.printStackTrace();
		}
	}

	// use this to eliminate NANs
	private static float[] getMedianLrrs(float[][] lrrs) {
		float[] medians = new float[lrrs[0].length];
		ArrayList<ArrayList<Float>> sampleLRRS = new ArrayList<ArrayList<Float>>();
		for (int i = 0; i < lrrs[0].length; i++) {
			sampleLRRS.add(new ArrayList<Float>());
		}
		for (int i = 0; i < lrrs.length; i++) {
			for (int k = 0; k < lrrs[i].length; k++) {
				// Not Including NaNs
				// if marker i, sample k
				if (Float.isNaN(lrrs[i][k])) {
					continue;
				} else {
					// add to sample k array list, marker i , sample k value
					sampleLRRS.get(k).add((lrrs[i][k]));
				}
			}
		}
		for (int i = 0; i < sampleLRRS.size(); i++) {
			medians[i] = Array.quants(toFloatArray(sampleLRRS.get(i)), QUANTILES)[0];
		}
		return medians;
	}

	private static float[] toFloatArray(ArrayList<Float> al) {
		float[] d = new float[al.size()];
		for (int i = 0; i < al.size(); i++) {
			d[i] = al.get(i);
		}
		return d;
	}

	private boolean checkUCSCRegions(String region) {
		boolean valid = false;
		int[] newLocation = Positions.parseUCSClocation(region);
		if ((newLocation == null || newLocation.length != 3 || newLocation[0] < 0) || (newLocation[1] < 0) || (newLocation[2] < 0)) {
			valid = false;
		} else {
			valid = true;
		}
		return valid;
	}

	private void assignMarkerIndices(MarkerRegion[] markerRegions) {
		for (int i = 0; i < markerNames.length; i++) {
			for (int j = 0; j < markerRegions.length; j++) {
				// UCSC indices already Assigned;
				if (markerRegions[i].getRegionType() == 1) {
					markerRegions[i].assignIndex(i, markerNames[i]);
				}
			}
		}
	}

	private void assignChr(MarkerRegion[] markerRegions) {
		for (int i = 0; i < markerRegions.length; i++) {
			String[] regionMarkers = markerRegions[i].returnMarkers();
			for (int k = 0; k < regionMarkers.length; k++) {
				if (!markerRegions[i].getIndex().containsKey(regionMarkers[k])) {
					computelog.reportError("Error - could not find " + regionMarkers[k] + " in markerLookup file , was this probeset analyzed?");
				} else {
					int index = markerRegions[i].getIndex().get(regionMarkers[k]);
					transChrs[(int) chrs[index]] = true;
				}
			}
		}

	}

	private class MarkerRegion {
		private String regionName;
		private int regionID;
		private int regionType; // 0 = UCSC ,1 =markerOnly
		private ArrayList<Integer> markerIndex;
		private ArrayList<String> markersInRegion;
		private Hashtable<String, Boolean> inRegion;
		private Hashtable<String, Integer> index;

		public Hashtable<String, Integer> getIndex() {
			return index;
		}

		public MarkerRegion(String regionName, int regionID, int regionType) {
			this.regionName = regionName;
			this.regionID = regionID;
			this.regionType = regionType;
			this.markerIndex = new ArrayList<Integer>();
			this.markersInRegion = new ArrayList<String>();
			this.inRegion = new Hashtable<String, Boolean>();
			this.index = new Hashtable<String, Integer>();
		}

		public void assignPositions() {
			for (int i = 0; i < markersInRegion.size(); i++) {
				if (index.containsKey(markersInRegion.get(i))) {
					markerIndex.add(index.get(markersInRegion.get(i)));
				} else {
					computelog.reportError("Error -could not find marker " + markersInRegion.get(i) + " in Marker Lookup");
				}
			}
		}

		public void assignIndex(int anindex, String markerName) {
			if (inRegion.containsKey(markerName)) {
				index.put(markerName, anindex);
			}
		}

		public String getRegionName() {
			return regionName;
		}

		public int getRegionType() {
			return regionType;
		}

		public void addMarker(String marker) {
			markersInRegion.add(marker);
		}

		public void addMarkerIndex(int position) {
			markerIndex.add(position);
		}

		public void markerInRegion(String marker) {
			inRegion.put(marker, true);
		}

		public void checkNumMarkers() {
			if (markersInRegion.size() < 1) {
				warnAndCancel("Error - All markers were filtered out of region " + regionName);
			}
		}

		public void addMarkers(String[] input) {
			for (int i = MARKER_REGION_START_OF_MARKERS; i < input.length; i++) {
				if (hash.containsKey(input[i])) {
					computelog.report(input[i] + " was filtered out");
				} else {
					processTracker[0]++;
					addMarker(input[i]);
					markerInRegion(input[i]);
				}
			}
			checkNumMarkers();
		}

		public void addUCSCMarkers(String UCSCLine) {
			Segment seg = new Segment(UCSCLine);
			for (int i = 0; i < positions.length; i++) {
				if (chrs[i] == seg.getChr() && positions[i] >= seg.getStart() && positions[i] <= seg.getStop()) {
					if (hash.containsKey(markerNames[i])) {
						computelog.report(markerNames[i] + " was filtered out");
					} else {
						processTracker[0]++;
						addMarker(markerNames[i]);
						markerInRegion(markerNames[i]);
						addMarkerIndex(i);
					}
				}
			}
			checkNumMarkers();
		}

		public String[] returnMarkers() {
			return markersInRegion.toArray(new String[markersInRegion.size()]);
		}
	}

	private float[] getRawValueMedianForRegion(MarkerRegion markerRegion) {
		String[] regionMarkers = markerRegion.returnMarkers();
		float[][] sampleLrrs = new float[regionMarkers.length][samples.length];
		float[] lrrs;
		newJob(ext.replaceAllWith(MEDIAN_WORKER_JOBS[2], "[%" + 0 + "]", markerRegion.getRegionName()));
		MarkerDataLoader markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSameThread(proj, regionMarkers, computelog);
		newJob(ext.replaceAllWith(MEDIAN_WORKER_JOBS[1], "[%" + 0 + "]", markerRegion.getRegionName()));

		for (int i = 0; i < regionMarkers.length; i++) {
			MarkerData markerData = markerDataLoader.requestMarkerData(i);
			processTracker[1]++;
			if (markerData.getFingerprint() != sampleList.getFingerprint()) {
				String Error = "Error - mismatched fingerprint for " + markerData.getMarkerName();
				computelog.reportError(Error);
				warnAndCancel(Error);
			}
			lrrs = markerData.getLRRs();
			markerDataLoader.releaseIndex(i);
			for (int j = 0; j < samples.length; j++) {

				try {
					// marker;samples...
					sampleLrrs[i][j] = lrrs[j];
				} catch (ArrayIndexOutOfBoundsException aioobe) {
					computelog.report("" + i + "\t" + j);
					System.exit(1);
				}
			}

		}
		process(processTracker[1]);
		return getMedianLrrs(sampleLrrs);
	}

	// assign

	private void getNormalizedMedianForRegions(MarkerRegion[] markerRegions) {
		assignMarkerIndices(markerRegions);
		// chromosome
		if (scope == 0) {
			assignChr(markerRegions);
		}

	}

	// private float[] getNormalizedMedianForRegion(MarkerRegion markerRegion) {
	// String[] regionMarkers = markerRegion.returnMarkers();
	// float[][] sampleLrrs = new float[regionMarkers.length][samples.length];
	// for (int i = 0; i < samples.length; i++) {
	// if (i % 100 == 0) {
	// computelog.report((i + 1) + " of " + samples.length + " (" + ext.getTimeElapsed(time) + ")");
	// time = new Date().getTime();
	// }
	// Sample samp = proj.getPartialSampleFromRandomAccessFile(samples[i]);
	// for (int trans = 0; trans < Transforms.TRANSFORMATION_TYPES.length; trans++) {
	// lrrs = samp.getLRRs();
	// if (trans > 0) {
	// lrrs = Transforms.transform(lrrs, trans, markerChrIndices, transChrs);
	// }
	// for (int j = 0; j < regions.length; j++) {
	// sum = 0;
	// count = 0;
	// for (int k = 0; k < indices[j].length; k++) {
	// if (!Double.isNaN(lrrs[indices[j][k]])) {
	// sum += lrrs[indices[j][k]];
	// count++;
	// }
	// }
	// data[j][i][trans] = (float) (sum / (double) count);
	// }
	// }
	// }
	// return getMedianLrrs(sampleLrrs);
	// }

	public static void createFilesFromFullSample(Project proj, String regionsFile, Logger log) {
		PrintWriter writer;
		MarkerSet markerSet;
		SampleList sampleList;
		Sample samp;
		String[] samples, markerNames;
		byte[] chrs;
		int[] positions;
		Segment[] regions;
		IntVector[] components;
		int[][] indices;
		float[] lrrs;
		float[][][] data;
		double sum;
		int count;
		long time;
		int[] numberOfMarkers;
		Hashtable<String, String> hash;
		int[][] markerChrIndices;
		boolean[] transChrs;
		transChrs = Array.booleanArray(27, false);
		log.report("Computing LRR values from FullSample. While this is slower, it allows for additional columns containing normalized values.");

		markerSet = proj.getMarkerSet();
		markerNames = markerSet.getMarkerNames();
		chrs = markerSet.getChrs();
		positions = markerSet.getPositions();

		hash = proj.getFilteredHash();

		regions = CNVariant.loadUCSCregions(proj.getProjectDir() + regionsFile, false);
		components = IntVector.newIntVectors(regions.length);
		for (int i = 0; i < positions.length; i++) {
			for (int j = 0; j < regions.length; j++) {
				if (chrs[i] == regions[j].getChr() && positions[i] >= regions[j].getStart() && positions[i] <= regions[j].getStop()) {
					if (hash.containsKey(markerNames[i])) {
						log.report(markerNames[i] + " was filtered out");
					} else {
						components[j].add(i);
					}
				}
			}
		}
		indices = IntVector.toIntMatrix(components);
		numberOfMarkers = new int[regions.length];
		for (int i = 0; i < regions.length; i++) {
			numberOfMarkers[i] = indices[i].length;
		}

		transChrs = Array.booleanArray(27, false);
		try {
			writer = new PrintWriter(new FileWriter(proj.getProjectDir() + "MarkersIn_" + ext.rootOf(regionsFile) + ".xln"));
			for (int i = 0; i < regions.length; i++) {
				writer.print(regions[i].getUCSClocation() + "\t" + numberOfMarkers[i] + "\t" + regions[i].getChr() + "\t" + regions[i].getStart() + "\t" + regions[i].getStop());
				transChrs[regions[i].getChr()] = true;
				for (int j = 0; j < indices[i].length; j++) {
					writer.print("\t" + markerNames[indices[i][j]]);
				}
				writer.println();
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing the list of marker names within the regions");
			e.printStackTrace();
		}

		markerChrIndices = markerSet.getIndicesByChr();
		sampleList = proj.getSampleList();
		samples = sampleList.getSamples();
		data = new float[regions.length][samples.length][Transforms.TRANSFORMATION_TYPES.length];
		log.report("Computing mean Log R ratios for:");
		time = new Date().getTime();
		for (int i = 0; i < samples.length; i++) {
			if (i % 100 == 0) {
				log.report((i + 1) + " of " + samples.length + " (" + ext.getTimeElapsed(time) + ")");
				time = new Date().getTime();
			}
			samp = proj.getPartialSampleFromRandomAccessFile(samples[i]);
			for (int trans = 0; trans < Transforms.TRANSFORMATION_TYPES.length; trans++) {
				lrrs = samp.getLRRs();
				if (trans > 0) {
					lrrs = Transforms.transform(lrrs, trans, markerChrIndices, transChrs);
				}
				for (int j = 0; j < regions.length; j++) {
					sum = 0;
					count = 0;
					for (int k = 0; k < indices[j].length; k++) {
						if (!Double.isNaN(lrrs[indices[j][k]])) {
							sum += lrrs[indices[j][k]];
							count++;
						}
					}
					data[j][i][trans] = (float) (sum / (double) count);
				}
			}
		}

		new MeanLRRset(sampleList.getFingerprint(), regions, numberOfMarkers, data, Transforms.TRANFORMATIONS).serialize(proj.getProjectDir() + ext.rootOf(regionsFile) + ".mlrr");
	}

	public static void dump(Project proj, String[] phenotypes, String mlrrSetFile, String regionToDumpOrNullForAll, int transformationToUse, Logger log) {
		PrintWriter writer;
		MeanLRRset mlrrSet;
		float[][][] data;
		SampleList sampleList;
		String[] samples;
		Segment[] regions;
		Hashtable<String, String> hash;
		int index;
		String[] transformations;

		if (!Files.exists(mlrrSetFile)) {
			log.report("Error - mlrr dataset file '" + mlrrSetFile + "' was never created");
			return;
		}

		hash = HashVec.loadFileToHashString(proj.getFilename(Project.SAMPLE_DATA_FILENAME), "DNA", phenotypes, "\t");
		sampleList = proj.getSampleList();
		samples = sampleList.getSamples();

		mlrrSet = MeanLRRset.load(mlrrSetFile, false);
		data = mlrrSet.getData();
		regions = mlrrSet.getRegions();
		if (mlrrSet.getSampleFingerprint() != sampleList.getFingerprint()) {
			log.reportError("Error - the SampleList fingerprint for the MeanLRRset (" + mlrrSet.getSampleFingerprint() + ") does not match the Project's SampleList fingerprint (" + sampleList.getFingerprint() + ")");
			return;
		}

		if (regionToDumpOrNullForAll != null) {
			regionToDumpOrNullForAll = ext.replaceAllWith(regionToDumpOrNullForAll, new String[][] { { ",", "" } });
		}

		try {
			index = -1;
			if (regionToDumpOrNullForAll == null) {
				writer = new PrintWriter(new FileWriter(proj.getProjectDir() + ext.rootOf(mlrrSetFile) + "_dump.xln"));
			} else {
				writer = new PrintWriter(new FileWriter(proj.getProjectDir() + ext.replaceAllWith(regionToDumpOrNullForAll, ":", "_") + ".xln"));

				for (int i = 0; i < regions.length; i++) {
					if (regionToDumpOrNullForAll.equals(regions[i].getUCSClocation())) {
						index = i;
					}
				}
				if (index == -1) {
					log.reportError("Error - Region flagged to be dumped (" + regionToDumpOrNullForAll + ") was not found in the MeanLRRset's regions list");
					return;
				}

			}
			writer.print("Sample");
			for (int i = 0; i < phenotypes.length; i++) {
				writer.print("\t" + phenotypes[i].substring(phenotypes[i].lastIndexOf("=") + 1));
			}
			transformations = mlrrSet.getTransformations();
			for (int i = 0; i < regions.length; i++) {
				if (transformationToUse == -1) {
					for (int j = 0; j < transformations.length; j++) {
						writer.print("\t" + regions[i].getUCSClocation() + "_" + ext.replaceWithLinuxSafeCharacters(transformations[j], false));
					}
				} else {
					writer.print("\t" + regions[i].getUCSClocation() + "_" + ext.replaceWithLinuxSafeCharacters(transformations[transformationToUse], false));
				}
			}
			writer.println();
			for (int i = 0; i < samples.length; i++) {
				writer.print(samples[i] + "\t" + (hash.containsKey(samples[i]) ? hash.get(samples[i]) : Array.stringArray(phenotypes.length, ".")));
				if (regionToDumpOrNullForAll == null) {
					for (int j = 0; j < regions.length; j++) {
						if (transformationToUse == -1) {
							writer.print("\t" + Array.toStr(data[j][i]));
						} else {
							writer.print("\t" + data[j][i][transformationToUse]);
						}
					}
				} else {
					writer.print("\t" + Array.toStr(data[index][i]));
				}
				writer.println();
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + ext.rootOf(mlrrSetFile) + ".xln");
			e.printStackTrace();
		}
	}

	public static void createFilesFromMarkerData(Project proj, String regionsFile, Logger log) {
		PrintWriter writer;
		MarkerSet markerSet;
		SampleList sampleList;
		String[] samples, markerNames;
		byte[] chrs;
		int[] positions;
		Segment[] regions;
		Hashtable<String, Vector<String>> components;
		float[] lrrs;
		float[][][] data;
		int[] counts;
		int[] numberOfMarkers;
		Hashtable<String, String> hash;
		MarkerDataLoader markerDataLoader;
		MarkerData markerData;

		log.report("Computing LRR values from MarkerData. This is faster, but does not allow for additional columns containing normalized values.");

		markerSet = proj.getMarkerSet();
		markerNames = markerSet.getMarkerNames();
		chrs = markerSet.getChrs();
		positions = markerSet.getPositions();

		hash = proj.getFilteredHash();

		regions = CNVariant.loadUCSCregions(proj.getProjectDir() + regionsFile, false);
		if (regions == null) {
			return;
		}
		components = new Hashtable<String, Vector<String>>();
		for (int i = 0; i < positions.length; i++) {
			for (int j = 0; j < regions.length; j++) {
				if (chrs[i] == regions[j].getChr() && positions[i] >= regions[j].getStart() && positions[i] <= regions[j].getStop()) {
					if (hash.containsKey(markerNames[i])) {
						log.report(markerNames[i] + " was filtered out");
					} else {
						HashVec.addToHashVec(components, j + "", markerNames[i], false);
					}
				}
			}
		}

		sampleList = proj.getSampleList();
		samples = sampleList.getSamples();
		numberOfMarkers = new int[regions.length];
		data = new float[regions.length][samples.length][1]; // only mean will be computed; no normalization
		try {
			writer = new PrintWriter(new FileWriter(proj.getProjectDir() + "MarkersIn_" + ext.rootOf(regionsFile) + ".xln"));
			for (int i = 0; i < regions.length; i++) {
				markerNames = Array.toStringArray(components.get(i + ""));
				numberOfMarkers[i] = markerNames.length;
				writer.println(regions[i].getUCSClocation() + "\t" + markerNames.length + "\t" + regions[i].getChr() + "\t" + regions[i].getStart() + "\t" + regions[i].getStop() + "\t" + Array.toStr(markerNames));
				log.report("Computing mean Log R ratios for: " + regions[i].getUCSClocation());
				markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markerNames, log);
				counts = new int[samples.length];
				for (int j = 0; j < markerNames.length; j++) {
					markerData = markerDataLoader.requestMarkerData(j);
					if (markerData.getFingerprint() != sampleList.getFingerprint()) {
						log.reportError("Error - mismatched fingerprint for " + markerData.getMarkerName());
					}

					lrrs = markerData.getLRRs();
					for (int k = 0; k < samples.length; k++) {
						if (!ext.isMissingValue(lrrs[k] + "")) {
							data[i][k][0] += lrrs[k];
							counts[k]++;
						}
					}
				}
				for (int j = 0; j < samples.length; j++) {
					data[i][j][0] /= (float) counts[j];
				}
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing the list of marker names within the regions");
			e.printStackTrace();
		}

		new MeanLRRset(sampleList.getFingerprint(), regions, numberOfMarkers, data, new String[] { Transforms.TRANFORMATIONS[0] }).serialize(proj.getProjectDir() + ext.rootOf(regionsFile) + ".mlrr");
	}
}
