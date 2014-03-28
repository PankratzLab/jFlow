package cnv.analysis;

import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
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
	private static final String[] FILE_PREFIXES = { "LRR_MEDIAN_", "MarkersIn_" };
	private static final String[] FILE_EXT = { ".xln" };
	private static final String MARKER_REGION_PREFIX = "probeset_id";
	private static final String MARKER_REGION_DELIMITER = ";";
	private static final String MARKER_REGION_REGEX = "[%0]";
	private static final int MARKER_REGION_START_OF_MARKERS = 2;
	private static final double[] QUANTILES = { 0.5 };
	private static final String[] MEDIAN_WORKER_JOBS = { "Parsing and intitializing Regions", "Computing Median Log R Ratios for " + MARKER_REGION_REGEX, "Waiting for data to Load for Region " + MARKER_REGION_REGEX, "Creating Output Files" };
	private static final String[] CLASSES_TO_DUMP = { "IID" };
	private static final String[] MARKER_REGION_RESULTS_SUFFIX = { "MEDIAN", "MAD" };
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

	public static boolean checkExists(Project proj, String outputBase) {
		boolean exists = false;
		for (int i = 0; i < FILE_PREFIXES.length; i++) {
			for (int j = 0; j < FILE_EXT.length; j++) {
				if (Files.exists(proj.getProjectDir() + FILE_PREFIXES[i] + outputBase + FILE_EXT[j])) {
					exists = true;
				}
			}
		}
		return exists;
	}

	// for access from command line
	public static void computeMedianLrrs(Project proj, String regionFileName, int transfromationType, int scope, String outputBase, Logger log) {
		String[] input = readToArray(proj.getProjectDir() + regionFileName, log);
		MedianLRRWorker medianLRRWorker = new MedianLRRWorker(proj, input, transfromationType, scope, outputBase, null, log);
		log.report("Starting job for " + input.length + " regions");
		medianLRRWorker.execute();

	}

	public MedianLRRWorker(Project proj, String[] input, int transformationType, int scope, String outputBase, JProgressBar jProgressBar, Logger log) {
		this.proj = proj;
		this.input = input;
		this.outputBase = outputBase;
		this.transformationType = transformationType;
		this.scope = scope;
		this.progressBar = jProgressBar == null ? new JProgressBar() : jProgressBar;
		addPropertyChangeListener(new PropertyChangeListener() {
			public void propertyChange(PropertyChangeEvent evt) {
				if ("progress".equals(evt.getPropertyName())) {
					progressBar.setValue((Integer) evt.getNewValue());
				}
			}
		});

		this.computelog = log == null ? new Logger(proj.getProjectDir() + outputBase + ".log") : log;
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

	}

	protected Integer doInBackground() throws Exception {
		System.out.println(proj.getProjectDir() + outputBase);
		progressBar.setValue(0);
		setProgress(0);
		newJob(MEDIAN_WORKER_JOBS[0]);
		MarkerRegion[] markerRegions = parseRegions();
		process(processTracker[0]);
		if (processTracker[0] < 1) {
			String Error = "Error - did not find any markers in the Region(s) of interest";
			computelog.reportError(Error);
			warnAndCancel(Error);
		} else {
			String job;
			process(0);
			RegionResults regionResults;
			if (transformationType == 0) {
				job = ext.replaceAllWith(MEDIAN_WORKER_JOBS[1], "[%" + 0 + "]", Transforms.TRANFORMATIONS[transformationType]);
				newJob(job);
				computelog.report(job);
				assignMarkerProgress();
				regionResults = getRawValueMedianForRegions(markerRegions);
				newJob(MEDIAN_WORKER_JOBS[3]);
				process(processTracker[0] - 1);
				printResults(regionResults, markerRegions);
			} else {
				assignSampleProgress();
				job = ext.replaceAllWith(MEDIAN_WORKER_JOBS[1], "[%" + 0 + "]", Transforms.TRANFORMATIONS[transformationType] + " " + Transforms.SCOPES[scope]);
				computelog.report(job);
				// String[] smallSamples = { samples[0] };
				// samples = smallSamples;
				newJob(job);
				regionResults = getNormalizedMedianForRegions(markerRegions);
				printResults(regionResults, markerRegions);
			}
			Thread.sleep(1000);
		}
		return 42;
	}

	// total,processed,updateat,tracker,scale;

	private void assignMarkerProgress() {
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
			e.printStackTrace();
		} catch (CancellationException cce) {
			computelog.reportError("Error - Cancelling Median Log R Ratio Computation ");
			cce.printStackTrace();

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
			String Error = "Error - Zero Regions Were Entered!";
			computelog.reportError(Error);
			warnAndCancel(Error);
		}
		return regions;
	}

	// For each region load lrrs from markerData and compute

	private RegionResults getRawValueMedianForRegions(MarkerRegion[] markerRegions) {
		float[][] rawMediansRegions = new float[markerRegions.length][samples.length];
		float[][] rawMADRegions = new float[markerRegions.length][samples.length];

		for (int i = 0; i < markerRegions.length; i++) {
			// stores median ,MAD for each sample
			float[][] results = getRawValueResultsForRegion(markerRegions[i]);
			rawMediansRegions[i] = results[0];
			rawMADRegions[i] = results[1];
		}
		return new RegionResults(rawMediansRegions, rawMADRegions);
	}

	private float[][] getRawValueResultsForRegion(MarkerRegion markerRegion) {
		String[] regionMarkers = markerRegion.returnMarkers();
		float[][] sampleLrrs = new float[regionMarkers.length][samples.length];
		newJob(ext.replaceAllWith(MEDIAN_WORKER_JOBS[2], "[%" + 0 + "]", markerRegion.getRegionName()));
		MarkerDataLoader markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, regionMarkers, computelog);
		newJob(ext.replaceAllWith(MEDIAN_WORKER_JOBS[1], "[%" + 0 + "]", markerRegion.getRegionName()));
		for (int i = 0; i < regionMarkers.length; i++) {
			MarkerData markerData = markerDataLoader.requestMarkerData(i);
			processTracker[1]++;
			if (markerData.getFingerprint() != sampleList.getFingerprint()) {
				String Error = "Error - mismatched fingerprint for " + markerData.getMarkerName();
				computelog.reportError(Error);
				warnAndCancel(Error);
			}
			float[] lrrs = markerData.getLRRs();
			markerDataLoader.releaseIndex(i);
			for (int j = 0; j < samples.length; j++) {
				try {
					// marker;samples...
					sampleLrrs[i][j] = lrrs[j];
				} catch (ArrayIndexOutOfBoundsException aioobe) {
					computelog.report("" + i + "\t" + j);
					aioobe.printStackTrace();
					System.exit(1);
				}
			}
		}
		process(processTracker[1]);
		return getSampleResults(sampleLrrs);
	}

	private float[][] getSampleResults(float[][] lrrs) {
		// store median, MAD
		float[][] results = new float[2][lrrs[0].length];
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
			results[0][i] = Array.quants(toFloatArray(sampleLRRS.get(i)), QUANTILES)[0];
			results[1][i] = getMAD(toFloatArray(sampleLRRS.get(i)), results[0][i]);
		}
		return results;
	}

	private void printResults(RegionResults regionResults, MarkerRegion[] markerRegions) {
		printRegionMarkers(markerRegions);
		printMedianLRRs(regionResults, markerRegions);
	}

	// print median values
	private void printMedianLRRs(RegionResults regionResults, MarkerRegion[] markerRegions) {
		String output = proj.getProjectDir() + FILE_PREFIXES[0] + outputBase + FILE_EXT[0];
		Hashtable<String, String> hashSamps = HashVec.loadFileToHashString(proj.getFilename(Project.SAMPLE_DATA_FILENAME), "DNA", CLASSES_TO_DUMP, "\t");
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(output));
			writer.print("Sample");
			for (int i = 0; i < CLASSES_TO_DUMP.length; i++) {
				writer.print("\t" + CLASSES_TO_DUMP[i].substring(CLASSES_TO_DUMP[i].lastIndexOf("=") + 1));
			}
			for (int i = 0; i < markerRegions.length; i++) {
				for (int j = 0; j < MARKER_REGION_RESULTS_SUFFIX.length; j++) {
					writer.print("\t" + MARKER_REGION_RESULTS_SUFFIX[j] + "_" + markerRegions[i].getRegionName());
				}
			}
			writer.println();
			for (int i = 0; i < samples.length; i++) {
				writer.print(samples[i] + "\t" + (hashSamps.containsKey(samples[i]) ? hashSamps.get(samples[i]) : Array.stringArray(CLASSES_TO_DUMP.length, ".")));
				for (int j = 0; j < markerRegions.length; j++) {
					writer.print("\t" + regionResults.getMedianAt(j, i) + "\t" + regionResults.getMADAt(j, i));
				}
				writer.println();
			}
			writer.close();
		} catch (Exception e) {
			computelog.reportError("Error writing to Medain Log R Ratios to " + output);
			e.printStackTrace();
		}
	}

	// print markers in region
	private void printRegionMarkers(MarkerRegion[] markerRegions) {
		String output = proj.getProjectDir() + FILE_PREFIXES[1] + outputBase + FILE_EXT[0];
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
				if (markerRegions[j].getRegionType() == 1) {
					markerRegions[j].assignIndex(i, markerNames[i]);
				}
			}
		}
		for (int i = 0; i < markerRegions.length; i++) {
			markerRegions[i].assignPositions();
		}
	}

	// only if transforming by chromosome
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

	// for storing region information
	private class MarkerRegion {
		private String regionName;
		// private int regionID;
		private int regionType; // 0 = UCSC ,1 =markerOnly
		private ArrayList<Integer> markerIndex;
		private ArrayList<String> markersInRegion;
		private Hashtable<String, Boolean> inRegion;
		private Hashtable<String, Integer> index;

		private Hashtable<String, Integer> getIndex() {
			return index;
		}

		private MarkerRegion(String regionName, int regionID, int regionType) {
			this.regionName = regionName;
			// this.regionID = regionID;
			this.regionType = regionType;
			// stores marker indices
			this.markerIndex = new ArrayList<Integer>();
			// stores markers in the region
			this.markersInRegion = new ArrayList<String>();
			// for assigning in region
			this.inRegion = new Hashtable<String, Boolean>();
			// for assigning index
			this.index = new Hashtable<String, Integer>();
		}

		public ArrayList<Integer> getMarkerIndex() {
			return markerIndex;
		}

		private void assignPositions() {
			for (int i = 0; i < markersInRegion.size(); i++) {
				if (index.containsKey(markersInRegion.get(i))) {
					markerIndex.add(index.get(markersInRegion.get(i)));
				} else {
					computelog.reportError("Error -could not find marker " + markersInRegion.get(i) + " in Marker Lookup");
				}
			}
		}

		private void assignIndex(int anindex, String markerName) {
			if (inRegion.containsKey(markerName)) {
				index.put(markerName, anindex);
			}
		}

		private String getRegionName() {
			return regionName;
		}

		private int getRegionType() {
			return regionType;
		}

		private void addMarker(String marker) {
			markersInRegion.add(marker);
		}

		private void addMarkerIndex(int position) {
			markerIndex.add(position);
		}

		private void markerInRegion(String marker) {
			inRegion.put(marker, true);
		}

		private void checkNumMarkers() {
			if (markersInRegion.size() < 1) {
				String Error = "Error - All markers were filtered out of region " + regionName;
				computelog.reportError(Error);
				warnAndCancel(Error);
			}
		}

		private void addMarkers(String[] input) {
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

		private void addUCSCMarkers(String UCSCLine) {
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
						assignIndex(i, markerNames[i]);
					}
				}
			}
			checkNumMarkers();
		}

		private String[] returnMarkers() {
			return markersInRegion.toArray(new String[markersInRegion.size()]);
		}
	}

	private class RegionResults {
		// stored as region, sample
		private float[][] regionMedianValues;
		private float[][] regionMADValues;

		public RegionResults(float[][] regionMedianValues, float[][] regionMADValues) {
			super();
			this.regionMedianValues = regionMedianValues;
			this.regionMADValues = regionMADValues;
		}

		public float getMedianAt(int region, int sample) {
			return regionMedianValues[region][sample];
		}

		public float getMADAt(int region, int sample) {
			return regionMADValues[region][sample];
		}

	}

	private RegionResults getNormalizedMedianForRegions(MarkerRegion[] markerRegions) {
		int[][] indices;
		float[][] regionMedianValues = new float[markerRegions.length][samples.length];
		float[][] regionMADValues = new float[markerRegions.length][samples.length];
		assignMarkerIndices(markerRegions);
		computelog.report("Computing mean Log R ratios for:");
		// norm by genome
		if (scope == 1) {
			indices = new int[][] { Array.intArray(proj.getPartialSampleFromRandomAccessFile(samples[0]).getLRRs().length) };
		} else {
			// only normalize chrs with markers in regions
			indices = markerSet.getIndicesByChr();
			assignChr(markerRegions);
		}
		long time = new Date().getTime();
		for (int i = 0; i < samples.length; i++) {
			newJob(ext.replaceAllWith(MEDIAN_WORKER_JOBS[1], "[%" + 0 + "]", samples[i]));
			process(i + 1);
			if (i % 100 == 0) {
				time = new Date().getTime();
				computelog.report((i + 1) + " of " + samples.length + " (" + ext.getTimeElapsed(time) + ")");
			}
			float[] lrrs = proj.getPartialSampleFromRandomAccessFile(samples[i]).getLRRs();
			// genome
			if (scope == 1) {
				lrrs = Transforms.transform(lrrs, transformationType, false, markerSet);
			}
			// default to norm by chromosome
			else {
				//markerSet.getIndicesByChr()
				lrrs = Transforms.transform(lrrs, transformationType, indices, transChrs);
			}
			for (int j = 0; j < markerRegions.length; j++) {
				ArrayList<Integer> regionPositions = markerRegions[j].getMarkerIndex();
				ArrayList<Float> regionLrrs = new ArrayList<Float>();
				for (int k = 0; k < regionPositions.size(); k++) {
					if (Float.isNaN(lrrs[regionPositions.get(k)])) {
						continue;
					} else {
						regionLrrs.add(lrrs[regionPositions.get(k)]);
					}
				}

				regionMedianValues[j][i] = Array.quants(toFloatArray(regionLrrs), QUANTILES)[0];
				regionMADValues[j][i] = getMAD(toFloatArray(regionLrrs), regionMedianValues[j][i]);
			}
		}
		return new RegionResults(regionMedianValues, regionMADValues);
	}

	private float getMAD(float[] regionLrrs, float median) {
		float[] diffs = new float[regionLrrs.length];
		for (int i = 0; i < regionLrrs.length; i++) {
			diffs[i] = Math.abs((regionLrrs[i] - median));
		}
		return Array.quants(diffs, QUANTILES)[0];
	}

	private static String[] readToArray(String filename, Logger log) {

		ArrayList<String> lines = new ArrayList<String>();
		try {
			FileReader fileReader = new FileReader(filename);
			BufferedReader reader = new BufferedReader(fileReader);
			while (reader.ready()) {
				lines.add(reader.readLine().trim());
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error unable to find" + filename);

		} catch (IOException e) {
			log.reportError("Error reading file " + filename);
			e.printStackTrace();
		}
		return lines.toArray(new String[lines.size()]);
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = Project.DEFAULT_PROJECT;
		String regionFileName = "cnps.txt";
		String headless = "true";
		int transformationType = 0;
		int scope = 0;
		long time;
		String outputBase = Transforms.TRANFORMATIONS[transformationType];
		String logfile = outputBase + ".log";
		Logger log;
		Project proj;

		String usage = "cnv.analysis.MedianLRRWorker requires 2 arguments\n" + "" + "   (1) project file (i.e. proj=" + filename + " (default))\n" + "   (2) filename of the regions (one per line) in UCSC format (chr8:25129632-25130278) \n" + "       OR:\n" + "       formatted as \"" + MARKER_REGION_PREFIX + MARKER_REGION_DELIMITER + "(Your Region Name)" + MARKER_REGION_DELIMITER + "marker name 1" + MARKER_REGION_DELIMITER + "marker name 2...\"" + MARKER_REGION_DELIMITER + "\n" + "       (i.e. regions=" + regionFileName + "(default))\n" + "       OPTIONAL:\n" + "   (3) transformation type (i.e. transform=0 (default, " + Transforms.TRANFORMATIONS[transformationType] + ")) \n" + "       transformations are: " + Array.toStr(Transforms.TRANFORMATIONS) + "\n" + "   (4) scope of transformation (i.e. scope=0 (default))\n" + "       scopes are: " + Array.toStr(Transforms.SCOPES) + "\n" + "   (5) base name of the output files (i.e out=" + outputBase + " (default))\n" + "   (6) name of the log file (i.e. log=" + logfile + "\n" + "   (7) run program in headless mode to quiet gui errors when X11 forwarding\n is un-available (i.e. headless=true (default));" + "";
		// String usage= "cnv.analysis.MedianLRRWorker requires 2 arguments\n"+"" +
		// "   (1) project file (i.e. proj="+filename+" (default))\n"+
		// "   (2) filename of the regions (one per line) in UCSC format (chr8:25129632-25130278) \n"+
		// "       OR:\n"+
		// "       formatted as \"" + MARKER_REGION_PREFIX + MARKER_REGION_DELIMITER + "(Your Region Name)" + MARKER_REGION_DELIMITER + "marker name 1" + MARKER_REGION_DELIMITER + "marker name 2...\""+ MARKER_REGION_DELIMITER +"\n"+
		// "       (i.e. regions="+regionFileName+"(default))\n"+
		// "       OPTIONAL:\n"+
		// "   (3) transformation type (i.e. transform=0 (default, "+Transforms.TRANFORMATIONS[transformationType]+")) \n"+
		// "       transformations are: "+ Array.toStr(Transforms.TRANFORMATIONS)+"\n"+
		// "   (4) scope of transformation (i.e. scope=0 (default))\n"+
		// "       scopes are: "+Array.toStr(Transforms.SCOPES)+"\n"+
		// "   (5) base name of the output files (i.e out="+outputBase+" (default))\n"+
		// "   (6) name of the log file (i.e. log="+logfile+"\n"+
		// "   (7) run program in headless mode to quiet gui errors when X11 forwarding\n is un-available (i.e. headless=true (default));"+
		// "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				return;
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("regions=")) {
				regionFileName = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("transform=")) {
				transformationType = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("scope=")) {
				scope = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("out=")) {
				outputBase = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("headless=")) {
				headless = ext.parseStringArg(args[i], null);
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			return;
		}
		try {
			log = new Logger(logfile);
			time = new Date().getTime();
			proj = new Project(filename, false);
			System.setProperty("java.awt.headless", headless);
			MedianLRRWorker medianLRRWorker = new MedianLRRWorker(proj, readToArray(proj.getProjectDir() + regionFileName, log), transformationType, scope, outputBase, null, log);
			medianLRRWorker.execute();
			while (!medianLRRWorker.isDone()) {
				Thread.sleep(100);
			}
			log.report("Finished in " + ext.getTimeElapsed(time));
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

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
