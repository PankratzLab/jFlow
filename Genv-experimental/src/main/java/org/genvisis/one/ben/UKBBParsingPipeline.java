package org.genvisis.one.ben;

import htsjdk.variant.variantcontext.Allele;

import java.io.BufferedReader;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Timer;
import java.util.TimerTask;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;

import org.apache.commons.collections4.iterators.ArrayIterator;
import org.genvisis.CLI;
import org.genvisis.cnv.filesys.Compression;
import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.MarkerDetailSet;
import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.genvisis.cnv.filesys.MarkerLookup;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.filesys.SampleList;
import org.genvisis.cnv.manage.Markers;
import org.genvisis.cnv.manage.PlinkData;
import org.genvisis.cnv.manage.Resources.GENOME_BUILD;
import org.genvisis.cnv.manage.TransposeData;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Elision;
import org.genvisis.common.Files;
import org.genvisis.common.GenomicPosition;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.Positions;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.ext;

public class UKBBParsingPipeline {

	int maxMDRAFThreadsPerChr = 8;
	int maxChrThreads = 4;
	int maxMarkersPerMDRAF = 1500;
	int logEveryNMins = 5;
	double scaleFactor = 100d;

	Logger log = new Logger();

	String sourceDir;
	String projDir;
	String propFileDir;
	String projName;

	// AB fields in RAF files are only 2 bits and can't hold full genotype values
	final boolean importGenoAsAB = false;

	Project proj;
	HashMap<Integer, FileSet> fileSets;
	long fingerprint = 0L;
	String famFile;
	String annotFile;
	String[][] famData;

	public void setProjectName(String projName2) {
		projName = projName2;
	}

	public void setProjectPropertiesDir(String propFileDir2) {
		propFileDir = propFileDir2;
	}

	public void setProjectDir(String projDir2) {
		projDir = projDir2;
	}

	public void setSourceDir(String sourceDir2) {
		sourceDir = sourceDir2;
	}

	public void setFamFile(String famFile2) {
		famFile = famFile2;
	}

	public void setAnnotationCSV(String csv) {
		annotFile = csv;
	}

	private void setScaleFactor(double scale) {
		scaleFactor = scale;
	}

	private void setMarkersPerMDRAF(int mkrsPerFile) {
		maxMarkersPerMDRAF = mkrsPerFile;
	}

	private void setThreadsPerChromosomeCount(int threadsPerChr) {
		maxMDRAFThreadsPerChr = threadsPerChr;
	}

	private void setChromosomeThreadCount(int chrThreads) {
		maxChrThreads = chrThreads;
	}

	public void runPipeline() {
		createProject();
		createSampleList();
		discoverFilesets();
		parseMarkerData();
		writeLookup();
		createMarkerPositions();
		writeMarkerSet();
		writeMarkerDetailSet();
		try {
			writeOutliers();
		} catch (ClassNotFoundException | IOException e) {
			log.reportError("Problem occurred while loading outliers from marker files. Attempting to continue...");
			log.reportException(e);
		}
		createSampRAFsFromMDRAFs();
		createPED();
		createSampleData();
	}

	protected void createProject() {
		String propFile = ext.verifyDirFormat(propFileDir)
											+ ext.replaceWithLinuxSafeCharacters(projName) + ".properties";
		if (!Files.exists(propFile)) {
			Files.write((new Project()).PROJECT_NAME.getName() + "=" + projName, propFile);
			proj = new Project(propFile);
			proj.PROJECT_NAME.setValue(projName);
			proj.PROJECT_DIRECTORY.setValue(projDir);
			proj.SOURCE_DIRECTORY.setValue(sourceDir);
			proj.XY_SCALE_FACTOR.setValue(scaleFactor);
			proj.TARGET_MARKERS_FILENAMES.setValue(new String[] {});
			proj.SOURCE_FILENAME_EXTENSION.setValue("NULL");
			proj.ID_HEADER.setValue("NULL");
			proj.SOURCE_FILENAME_EXTENSION.setValue("NULL");
			proj.GENOME_BUILD_VERSION.setValue(GENOME_BUILD.HG19);

			proj.ARRAY_TYPE.setValue(ARRAY.AFFY_AXIOM);

			proj.saveProperties();
			log.reportTime("Created project properties file: " + propFile);
		} else {
			log.reportTime("Project properties file already exists; skipping creation.");
			proj = new Project(propFile);
		}
	}

	protected void createSampleList() {
		famData = HashVec.loadFileToStringMatrix(famFile, false, new int[] {0, 1, 2, 3,
																																				4, 5});
		String[] allSamples = Matrix.extractColumn(famData, 1);
		@SuppressWarnings("deprecation")
		long f = org.genvisis.cnv.filesys.MarkerSet.fingerprint(allSamples);
		fingerprint = f;
		if (!Files.exists(proj.SAMPLELIST_FILENAME.getValue())) {
			SampleList sl = new SampleList(allSamples);
			sl.serialize(proj.SAMPLELIST_FILENAME.getValue());
			sl = null;
			allSamples = null;
			log.reportTime("Created sample list file: " + proj.SAMPLELIST_FILENAME.getValue());
		} else {
			log.reportTime("Project sample list file already exists; skipping creation.");
		}
	}

	protected void discoverFilesets() {
		List<File> allSrcFiles = discoverAllSourceFiles(sourceDir, log);

		fileSets = new HashMap<>();

		for (File srcF : allSrcFiles) {
			String name = srcF.getName();
			String[] parts = name.split("_");
			int chr = Positions.chromosomeNumber(parts[2]);
			FileSet fs = fileSets.get(chr);
			if (fs == null) {
				fs = new FileSet();
				fs.log = this.log;
				fs.chr = chr;
				fs.famFile = famFile;
				fileSets.put(chr, fs);
			}

			switch (parts[1]) {
				case "l2r":
					if (fs.lrrFile != null) {
						log.reportError("Identified assumed duplicate " + "LRR" + " files: {" + fs.lrrFile
														+ " | " + srcF.getAbsolutePath() + "}.  Parsing has failed.");
						System.exit(1);
					} else {
						fs.lrrFile = srcF.getAbsolutePath();
					}
					break;
				case "baf":
					if (fs.bafFile != null) {
						log.reportError("Identified assumed duplicate " + "BAF" + " files: {" + fs.bafFile
														+ " | " + srcF.getAbsolutePath() + "}.  Parsing has failed.");
						System.exit(1);
					} else {
						fs.bafFile = srcF.getAbsolutePath();
					}
					break;
				case "con":
					if (fs.conFile != null) {
						log.reportError("Identified assumed duplicate " + "CON" + " files: {" + fs.conFile
														+ " | " + srcF.getAbsolutePath() + "}.  Parsing has failed.");
						System.exit(1);
					} else {
						fs.conFile = srcF.getAbsolutePath();
					}
					break;
				case "int":
					if (fs.intFile != null) {
						log.reportError("Identified assumed duplicate " + "INT" + " files: {" + fs.intFile
														+ " | " + srcF.getAbsolutePath() + "}.  Parsing has failed.");
						System.exit(1);
					} else {
						fs.intFile = srcF.getAbsolutePath();
					}
					break;
				case "snp":
					if (fs.bimFile != null) {
						log.reportError("Identified assumed duplicate " + "SNP" + " files: {" + fs.bimFile
														+ " | " + srcF.getAbsolutePath() + "}.  Parsing has failed.");
						System.exit(1);
					} else {
						fs.bimFile = srcF.getAbsolutePath();
					}
					break;
				case "cal":
					if (fs.bedFile != null) {
						log.reportError("Identified assumed duplicate " + "CAL" + " files: {" + fs.bedFile
														+ " | " + srcF.getAbsolutePath() + "}.  Parsing has failed.");
						System.exit(1);
					} else {
						fs.bedFile = srcF.getAbsolutePath();
					}
					break;
				default:
					log.reportTimeWarning("Parsing NOT YET IMPLEMENTED for the following data file: "
																+ srcF.getAbsolutePath());
					break;
			}
		}

		ArrayList<Integer> missingChrs = new ArrayList<>();
		ArrayList<Integer> incompleteChrs = new ArrayList<>();
		for (int i = 0; i < 27; i++) {
			if (!fileSets.containsKey(i)) {
				missingChrs.add(i);
			} else if (!fileSets.get(i).isComplete()) {
				incompleteChrs.add(i);
				fileSets.remove(i);
			}
		}

		if (!missingChrs.isEmpty()) {
			log.reportTimeWarning("Missing " + missingChrs.size() + " chrs: "
														+ ArrayUtils.toStr(missingChrs, ", "));
		}
		if (!incompleteChrs.isEmpty()) {
			log.reportTimeWarning("Removed chrs with incomplete files: "
														+ ArrayUtils.toStr(incompleteChrs, ", "));
		}
	}

	protected void parseMarkerData() {
		long time = System.nanoTime();

		ExecutorService executor = Executors.newFixedThreadPool(Math.min(fileSets.size(),
																																		 maxChrThreads));

		final MarkerLogger mkrLog = new MarkerLogger(log, logEveryNMins);
		CountDownLatch cdl = new CountDownLatch(fileSets.size());

		for (Integer v : fileSets.keySet()) {
			final FileSet fs = fileSets.get(v);
			Runnable runn = new Runnable() {
				@Override
				public void run() {
					try {
						fs.buildLookups(log);
						createMDRAF(fs, mkrLog);
						cdl.countDown();
					} catch (IOException | Elision e) {
						log.reportError("FAILED TO PARSE DATA FILE FOR CHROMOSOME: " + fs.chr);
						e.printStackTrace();
						System.exit(1);
					} catch (Exception e) {
						log.reportError("FAILED TO PARSE DATA FILE FOR CHROMOSOME: " + fs.chr);
						e.printStackTrace();
						System.exit(1);
					}
				}
			};
			executor.submit(runn);
		}
		log.reportTime("Queued MarkerData parsing for " + fileSets.size() + " chrs.");

		executor.shutdown();
		try {
			cdl.await(Long.MAX_VALUE, TimeUnit.DAYS);
		} catch (InterruptedException e) {
			/**/
		}
		mkrLog.stopLogging();

		log.reportTime("Completed parsing MarkerData in " + ext.getTimeElapsedNanos(time));
	}

	protected void writeLookup() {
		if (!Files.exists(proj.MARKERLOOKUP_FILENAME.getValue())) {
			Hashtable<String, String> lookup = new Hashtable<>();
			for (FileSet fs : fileSets.values()) {
				String mdRAF;
				List<String[]> mkrBins = fs.mkrBins;
				int start = 0;
				int end = 0;
				for (String[] bin : mkrBins) {
					end = start + bin.length;
					mdRAF = getMDRAFName(fs.chr, start, end);
					for (int m = 0; m < bin.length; m++) {
						lookup.put(bin[m], mdRAF + "\t" + m);
					}
					start = end;
				}
			}
			new MarkerLookup(lookup).serialize(proj.MARKERLOOKUP_FILENAME.getValue());
			log.reportTime("Created marker lookup file: " + proj.MARKERLOOKUP_FILENAME.getValue());
		} else {
			log.reportTime("Project marker lookup file already exists; skipping creation.");
		}
	}

	protected void createMarkerPositions() {
		String file = proj.MARKER_POSITION_FILENAME.getValue();
		if (!Files.exists(file)) {
			long time = System.nanoTime();
			PrintWriter writer = Files.getAppropriateWriter(file);
			writer.println("Marker\tChr\tPosition");

			for (int i = 0; i < 27; i++) {
				if (!fileSets.containsKey(i)) {
					continue;
				}
				FileSet fs = fileSets.get(i);

				try {
					BufferedReader reader = Files.getAppropriateReader(fs.bimFile);

					String line = null;
					String[] parts;

					while ((line = reader.readLine()) != null) {
						parts = line.split("\t", -1);
						writer.println(parts[1] + "\t" + parts[0] + "\t" + parts[3]);
					}
					reader.close();
				} catch (IOException e) {
					log.reportError("ERROR - Major failure when create markerPositions file.");
					log.reportException(e);
					System.exit(1);
				}
			}
			writer.flush();
			writer.close();

			log.reportTime("Completed markerPositions: " + proj.MARKER_POSITION_FILENAME.getValue()
										 + " in "
										 + ext.getTimeElapsedNanos(time));
		} else {
			log.reportTime("Project markerPositions file already exists; skipping creation.");
		}
	}

	protected void writeMarkerSet() {
		if (!Files.exists(proj.MARKERSET_FILENAME.getValue())) {
			long t1 = System.nanoTime();
			String[] mkrs = HashVec.loadFileToStringArray(proj.MARKER_POSITION_FILENAME.getValue(), true,
																										new int[] {0}, false);
			Markers.orderMarkers(mkrs, proj.MARKER_POSITION_FILENAME.getValue(),
													 proj.MARKERSET_FILENAME.getValue(), log);
			mkrs = null;
			log.reportTime("Completed markerSet file in " + ext.getTimeElapsedNanos(t1));
		} else {
			log.reportTime("Project marker set file already exists; skipping creation.");
		}
	}

	protected void writeMarkerDetailSet() {
		if (!Files.exists(proj.MARKER_DETAILS_FILENAME.getValue())) {
			if (null == annotFile || "".equals(annotFile)) {
				log.reportError("Affymetrix annotation file not found.  To call CNVs, A/B alleles will be needed and determining these requires an Affymetrix annotation file.");
				return;
			}
			Iterable<Marker> markerIterable;
			try {
				markerIterable = loadMarkers();
				MarkerDetailSet mds = new MarkerDetailSet(markerIterable);
				mds.serialize(proj.MARKER_DETAILS_FILENAME.getValue());
			} catch (IOException e) {
				e.printStackTrace();
			}
		} else {
			log.report("Project marker detail set file already exists; skipping creation.");
		}
	}

	private Iterable<Marker> loadMarkers() throws IOException {
		String[] mkrs = HashVec.loadFileToStringArray(proj.MARKER_POSITION_FILENAME.getValue(), true,
																									new int[] {0}, false);
		String[] headerFactors = {"Affy SNP ID", "dbSNP RS ID", "Chromosome",
															"Physical Position", "Allele A", "Allele B", "Ref Allele"};
		int[] hdrInds = null;
		BufferedReader reader = Files.getAppropriateReader(annotFile);
		HashMap<String, Integer> markerIndices = new HashMap<>();
		ArrayList<Integer> allInds = new ArrayList<>(mkrs.length);
		for (int i = 0; i < mkrs.length; i++) {
			markerIndices.put(mkrs[i], i);
			allInds.add(i);
		}
		final Marker[] markers = new Marker[mkrs.length];
		String line = null;
		boolean foundHeader = false;
		String[] parts;
		String mkrName;
		int inAnnotNotProj = 0;
		int ind;
		while ((line = reader.readLine()) != null) {
			if (line.charAt(0) == '#')
				continue;
			if (!foundHeader) {
				foundHeader = true; // first nonComment line
				hdrInds = ext.indexFactors(headerFactors, ext.splitCommasIntelligently(line, true, log),
																	 false);
				continue;
			}
			parts = ext.splitCommasIntelligently(line, true, log);
			mkrName = isMissing(parts[hdrInds[1]]) ? parts[hdrInds[0]] : parts[hdrInds[1]];

			if (!markerIndices.containsKey(mkrName)) {
				if (markerIndices.containsKey(parts[hdrInds[0]])) {
					mkrName = parts[hdrInds[0]];
					ind = markerIndices.get(mkrName);
				} else {
					inAnnotNotProj++;
					ind = -1;
				}
			} else {
				ind = markerIndices.get(mkrName);
			}
			if (ind != -1) {
				allInds.remove(Integer.valueOf(ind));
				markers[markerIndices.get(mkrName)] = new Marker(
																												 mkrName,
																												 new GenomicPosition(
																																						 !isMissing(parts[hdrInds[2]])
																																																					? Positions.chromosomeNumber(parts[hdrInds[2]],
																																																																			 log)
																																																					: (byte) 0,
																																						 !isMissing(parts[hdrInds[3]])
																																																					? Integer.parseInt(parts[hdrInds[3]])
																																																					: 0),
																												 !isMissing(parts[hdrInds[4]])
																																											? Allele.create(parts[hdrInds[4]],
																																																			parts[hdrInds[4]] == parts[hdrInds[6]])
																																											: Allele.NO_CALL,
																												 !isMissing(parts[hdrInds[5]])
																																											? Allele.create(parts[hdrInds[5]],
																																																			parts[hdrInds[5]] == parts[hdrInds[6]])
																																											: Allele.NO_CALL);
			}
		}
		if (inAnnotNotProj > 0) {
			log.reportTimeWarning("Found "
														+ inAnnotNotProj
														+ " markers in the annotation file that were not present in the project data.");
		}
		if (allInds.size() > 0) {
			log.reportError(allInds.size()
											+ " markers were not found in the annotation file! Missing marker indices: "
											+ allInds.toString());
		}
		reader.close();
		markerIndices = null;
		mkrs = null;
		return new Iterable<MarkerDetailSet.Marker>() {
			@Override
			public Iterator<Marker> iterator() {
				return new ArrayIterator<Marker>(markers);
			}
		};
	}

	private boolean isMissing(String value) {
		return value.equals("---") || value.equals("\"---\"") || value.equals("-")
					 || value.equals("\"-\"");
	}

	protected void writeOutliers() throws ClassNotFoundException, IOException {
		Hashtable<String, Float> allOutliers = new Hashtable<>();

		String[] mdRAFs = new File(proj.MARKER_DATA_DIRECTORY.getValue()).list((File f, String n) -> {
			return n.endsWith(MarkerData.MARKER_DATA_FILE_EXTENSION);
		});

		Map<String, Integer> mkrInds = proj.getMarkerIndices();
		String[] samples = proj.getSamples();
		for (String mdRAF : mdRAFs) {
			String[] fileMkrs = TransposeData.loadMarkerNamesFromRAF(proj.MARKER_DATA_DIRECTORY.getValue()
																															 + mdRAF);
			Hashtable<String, Float> outliers = TransposeData.loadOutliersFromRAF(proj.MARKER_DATA_DIRECTORY.getValue()
																																						+ mdRAF);
			for (Entry<String, Float> outlier : outliers.entrySet()) {
				String[] pts = outlier.getKey().split("\t");
				int mkrInd = mkrInds.get(fileMkrs[Integer.parseInt(pts[0])]);
				int sampInd = Integer.parseInt(pts[1]);
				allOutliers.put(mkrInd + "\t" + samples[sampInd] + "\t" + pts[2], outlier.getValue());
			}
		}

		SerializedFiles.writeSerial(allOutliers, proj.MARKER_DATA_DIRECTORY.getValue(true, true)
																						 + "outliers.ser");
	}

	protected void createSampRAFsFromMDRAFs() {
		if (!checkIfSampRAFsExist()) {
			TransposeData.reverseTranspose(proj);
		}
	}

	private boolean checkIfSampRAFsExist() {
		String[] listOfAllSamplesInProj = proj.getSamples();
		boolean allExist = true;
		for (String sample : listOfAllSamplesInProj) {
			String file = proj.SAMPLE_DIRECTORY.getValue(false, true) + sample
										+ Sample.SAMPLE_FILE_EXTENSION;
			if (!Files.exists(file)) {
				allExist = false;
				break;
			}
		}
		return allExist;
	}

	protected void createPED() {
		if (!Files.exists(proj.PEDIGREE_FILENAME.getValue())) {
			long t1 = System.nanoTime();
			PrintWriter writer = Files.getAppropriateWriter(proj.PEDIGREE_FILENAME.getValue());
			for (String[] ind : famData) {
				StringBuilder sb = new StringBuilder(ind[0]).append("\t");
				sb.append(ind[1]).append("\t");
				sb.append(ind[2]).append("\t");
				sb.append(ind[3]).append("\t");
				sb.append(ind[4]).append("\t");
				sb.append("0").append("\t");
				sb.append(ind[1]);
				writer.println(sb.toString());
			}
			writer.flush();
			writer.close();
			log.reportTime("Created ped file in " + ext.getTimeElapsedNanos(t1));
		} else {
			log.reportTime("Project sample list file already exists; skipping creation.");
		}
	}

	protected void createSampleData() {
		if (!Files.exists(proj.SAMPLE_DATA_FILENAME.getValue())) {
			try {
				SampleData.createSampleData(proj.PEDIGREE_FILENAME.getValue(), null, proj);
			} catch (Elision e) {
				e.printStackTrace();
				System.exit(1);
			}
		} else {
			log.reportTime("Project SampleData file already exists; skipping creation.");
		}
	}

	private byte getNullStatus() {
		return Sample.updateNullStatus(new float[0], new float[0],
																	 new float[0], new float[0],
																	 new float[0], importGenoAsAB ? new byte[0] : null,
																	 importGenoAsAB ? null : new byte[0], true);
	}


	private RandomAccessFile openMDRAF(String filename, int nInd, byte nullStatus,
																		 long fingerprint,
																		 String[] mkrNames) throws IOException {
		byte[] mkrBytes = Compression.objToBytes(mkrNames);
		byte[] mdRAFHeader = TransposeData.getParameterSectionForMdRaf(nInd,
																																	 mkrNames.length,
																																	 nullStatus,
																																	 fingerprint,
																																	 mkrBytes);
		mkrBytes = null;

		String file = proj.MARKER_DATA_DIRECTORY.getValue(true, true) + filename;
		if (Files.exists(file)) {
			new File(file).delete();
		}

		RandomAccessFile mdRAF = new RandomAccessFile(file, "rw");
		mdRAF.write(mdRAFHeader);
		mdRAFHeader = null;

		return mdRAF;
	}

	private void createMDRAF(FileSet fs, MarkerLogger mkrLog) throws IOException, Elision {
		long t1 = System.nanoTime();
		final String[][] bimData = HashVec.loadFileToStringMatrix(fs.bimFile, false, new int[] {0, 1,
																																														2, 3,
																																														4, 5});
		log.reportTime("Loaded chr" + fs.chr + " .bim data in " + ext.getTimeElapsedNanos(t1));
		final int nInd = famData.length;

		final byte nullStatus = getNullStatus();

		String[] mkrNames = Matrix.extractColumn(bimData, 1);
		fs.mkrBins = ArrayUtils.splitUpArray(mkrNames, (mkrNames.length / maxMarkersPerMDRAF) + 1, log);
		String[] existing = readWritten();
		ArrayList<Runnable> runners = new ArrayList<Runnable>();

		final CountDownLatch cdl = new CountDownLatch(fs.mkrBins.size());
		int start = 0;
		int end = 0;
		for (int i = 0; i < fs.mkrBins.size(); i++) {
			end = start + fs.mkrBins.get(i).length;
			final int s = start;
			final int ii = i;
			String mdRAFName = getMDRAFName(fs.chr, start, end);
			if (ext.indexOfStr(mdRAFName, existing) == -1) {
				runners.add(new Runnable() {
					@Override
					public void run() {
						String mdRAF;
						try {
							long t1 = System.nanoTime();
							mdRAF = write(fs, s, fs.mkrBins.get(ii), nullStatus, nInd, bimData, mkrLog);
							log.reportTime("Wrote " + mdRAF + " in " + ext.getTimeElapsedNanos(t1));
							success(mdRAF);
							cdl.countDown();
						} catch (Elision e) {
							e.printStackTrace();
						} catch (IOException e) {
							e.printStackTrace();
						} catch (Exception e) {
							e.printStackTrace();
						}
					}
				});
			} else {
				cdl.countDown();
			}
			start = end;
		}

		if (runners.size() > 0) {
			ExecutorService executor = Executors.newFixedThreadPool(Math.min(runners.size(),
																																			 maxMDRAFThreadsPerChr));
			for (Runnable r : runners) {
				executor.submit(r);
			}
			executor.shutdown();
			try {
				cdl.await(Long.MAX_VALUE, TimeUnit.DAYS);
			} catch (InterruptedException e) {
				/**/
			}

		}
	}

	private synchronized String[] readWritten() {
		String temp = proj.PROJECT_DIRECTORY.getValue() + "markersWritten.txt";
		return Files.exists(temp) ? HashVec.loadFileToStringArray(temp, false, null, false)
														 : new String[0];
	}

	private synchronized void success(String mkrFile) {
		String temp = proj.PROJECT_DIRECTORY.getValue() + "markersWritten.txt";
		Files.appendStringToFile(temp, mkrFile);
	}

	private String getMDRAFName(int chr, int start, int end) {
		return "markers." + chr + "." + start + "."
					 + end + MarkerData.MARKER_DATA_FILE_EXTENSION;
	}

	private String write(FileSet fs, int startBatchInd, String[] mkrNames, byte nullStatus, int nInd,
											 String[][] bimData, MarkerLogger markerLogger) throws Elision, IOException {

		FileReaderSet frs = new FileReaderSet();
		frs.open(fs);

		byte[] magicBytes = new byte[3];
		frs.bedIn.read(magicBytes);
		if (magicBytes[2] == 0) {
			log.reportError("Error - .bed file is sample-dominant: " + fs.bedFile);
			System.exit(1);
		}

		if (startBatchInd > 0) {
			frs.readAhead(fs, 0, startBatchInd);
		}

		RandomAccessFile mdRAF;

		Hashtable<String, Float> outOfRangeTable = new Hashtable<>();
		String mdRAFName = getMDRAFName(fs.chr, startBatchInd, (startBatchInd + mkrNames.length));

		mdRAF = openMDRAF(mdRAFName, famData.length, nullStatus, fingerprint, mkrNames);

		int bedBlockSize = (int) Math.ceil(nInd / 4.0);
		int binBlockSize = nInd * 8;

		int bytesPerSamp = Sample.getNBytesPerSampleMarker(getNullStatus());
		int markerBlockSize = nInd * bytesPerSamp;
		byte[] mkrBuff;
		String[] parts;
		double log2 = Math.log(2);
		for (int i = startBatchInd; i < startBatchInd + mkrNames.length; i++) {
			markerLogger.logMarker(fs.chr, startBatchInd, mkrNames.length, i);

			mkrBuff = new byte[markerBlockSize];
			int buffInd = 0;

			parts = frs.conIn.readLine().split(" ");
			if (parts.length != nInd) {
				log.reportError("Mismatched # of confidence scores {fnd: " + parts.length + ", exp: "
												+ nInd + "}.  File: " + fs.conFile);
				System.exit(1);
			}
			float gcV;
			for (String con : parts) {
				gcV = 1 - Float.parseFloat(con);
				if (gcV == 2) { // i.e. 1 - (-1)
					gcV = Float.NaN;
				}
				Compression.gcBafCompress(gcV, mkrBuff, buffInd);
				buffInd += Compression.REDUCED_PRECISION_GCBAF_NUM_BYTES;
			}
			parts = null;

			byte[] intensBytes = new byte[binBlockSize];
			frs.binIn.seek(((long) i) * ((long) binBlockSize));
			frs.binIn.read(intensBytes);

			float a, b, x, y;
			boolean oor;
			for (int bitInd = 0, binInd = 0; bitInd < nInd; bitInd++) {
				byte[] intA = {intensBytes[binInd++], intensBytes[binInd++], intensBytes[binInd++],
											 intensBytes[binInd++]};
				byte[] intB = {intensBytes[binInd++], intensBytes[binInd++], intensBytes[binInd++],
											 intensBytes[binInd++]};

				a = ByteBuffer.wrap(intA).order(ByteOrder.LITTLE_ENDIAN).getFloat();
				b = ByteBuffer.wrap(intB).order(ByteOrder.LITTLE_ENDIAN).getFloat();

				x = (float) (Math.log(a / b) / log2);
				y = (float) ((Math.log(a * b) / log2) / 2);
				// x = DoubleMath.log2((a / b), RoundingMode.HALF_UP);
				// y = DoubleMath.log2((a * b), RoundingMode.HALF_UP) / 2;
				// x = (float) (Math.log(x) / Math.log(2));
				// y = (float) (Math.log(y) / Math.log(2));

				oor = !Compression.xyCompressAllowNegative(x == -1
																													? Float.NaN
																													: (float) (x),
																									 mkrBuff, buffInd);
				if (oor) {
					outOfRangeTable.put((i - startBatchInd) + "\t" + bitInd + "\tx", x);
				}

				oor = !Compression.xyCompressAllowNegative(y == -1
																													? Float.NaN
																													: (float) (y),
																									 mkrBuff,
																									 buffInd
																											 + (nInd * Compression.REDUCED_PRECISION_XY_NUM_BYTES));
				if (oor) {
					outOfRangeTable.put((i - startBatchInd) + "\t" + bitInd + "\ty", y);
				}

				buffInd += Compression.REDUCED_PRECISION_XY_NUM_BYTES;
			}
			intensBytes = null;

			parts = frs.bafIn.readLine().split(" ");
			if (parts.length != nInd) {
				log.reportError("Mismatched # of BAF values {fnd: " + parts.length + ", exp: "
												+ nInd + "}.  File: " + fs.bafFile);
				System.exit(1);
			}
			float baf = Float.NaN;
			for (String bafStr : parts) {
				try {
					baf = Float.parseFloat(bafStr);
				} catch (NumberFormatException e) {
					if (ext.isMissingValue(bafStr)) {
						baf = Float.NaN;
					} else {
						log.reportError("Malformed BAF value: " + bafStr + " for chr " + fs.chr);
						System.exit(1);
					}
				}
				Compression.gcBafCompress(baf, mkrBuff, buffInd);
				buffInd += Compression.REDUCED_PRECISION_GCBAF_NUM_BYTES;
			}
			parts = null;

			parts = frs.lrrIn.readLine().split(" ");
			if (parts.length != nInd) {
				log.reportError("Mismatched # of L2R values {fnd: " + parts.length + ", exp: "
												+ nInd + "}.  File: " + fs.lrrFile);
				System.exit(1);
			}
			float lrr = Float.NaN;
			for (int samp = 0; samp < parts.length; samp++) {
				try {
					lrr = Float.parseFloat(parts[samp]);
				} catch (NumberFormatException e) {
					if (ext.isMissingValue(parts[samp])) {
						lrr = Float.NaN;
					} else {
						log.reportError("Malformed LRR value: " + parts[samp] + " for chr " + fs.chr);
						System.exit(1);
					}
				}
				oor = -1 == Compression.lrrCompress(lrr, mkrBuff, buffInd);
				buffInd += Compression.REDUCED_PRECISION_LRR_NUM_BYTES;
				if (oor) {
					outOfRangeTable.put((i - startBatchInd) + "\t" + samp + "\tlrr", lrr);
				}
			}
			parts = null;

			byte[] markerBytes = new byte[bedBlockSize];
			frs.bedIn.seek(3L + ((long) i) * ((long) bedBlockSize)); // ALWAYS seek forward
			frs.bedIn.read(markerBytes);
			for (int bitInd = 0; bitInd < markerBytes.length; bitInd++) {
				byte bedByte = markerBytes[bitInd];
				byte[] genotypes = PlinkData.decodeBedByte(bedByte);

				for (int g = 0; g < genotypes.length; g++) {
					int idInd = bitInd * 4 + g;
					if (idInd == -1 || idInd >= nInd) {
						continue;
					}
					mkrBuff[buffInd] = Compression.genotypeCompress(importGenoAsAB ? genotypes[g] : (byte) -1,
																													importGenoAsAB ? (byte) 0 : genotypes[g]);
					buffInd += Compression.REDUCED_PRECISION_ABFORWARD_GENOTYPE_NUM_BYTES;
				}
				genotypes = null;
			}
			markerBytes = null;

			mdRAF.write(mkrBuff);
			mkrBuff = null;
		}

		log.reportTime("Writing " + outOfRangeTable.size() + " outliers to " + mdRAFName);
		byte[] oorBytes = Compression.objToBytes(outOfRangeTable);
		mdRAF.write(Compression.intToBytes(oorBytes.length));
		mdRAF.write(oorBytes);
		log.reportTime("Finished writing outliers to " + mdRAFName);

		mdRAF.close();

		frs.close();
		frs = null;

		oorBytes = null;
		outOfRangeTable = null;
		return mdRAFName;
	}

	private static void addDirs(File src, List<File> subDirs) {
		for (File f : src.listFiles(File::isDirectory)) {
			if (!subDirs.contains(f)) {
				subDirs.add(f);
				addDirs(f, subDirs);
			}
		}
	}

	private static List<File> discoverAllSourceFiles(String dir, Logger log) {
		File srcDir = new File(dir);
		ArrayList<File> subDirs = new ArrayList<>();
		subDirs.add(srcDir);
		addDirs(srcDir, subDirs);

		log.reportTime("Found " + subDirs.size() + " source directories to scan.");

		/**
		 * Delim = "_"; <br />
		 * At least 3 tokens; <br />
		 * First token must be "ukb"; <br />
		 * Second token must be from FileSet.EXP; <br />
		 * Third token must start with "chr"; <br />
		 * Will not accept any .fam file; <br />
		 */
		FilenameFilter ff = new FilenameFilter() {
			@Override
			public boolean accept(File dir, String name) {
				String[] parts = name.split("_");
				if (parts.length <= 2 || parts.length > 4) {
					return false;
				}
				if (name.endsWith(".fam")) {
					return false;
				}
				if (!parts[0].equals("ukb")) {
					return false;
				}
				if (ext.indexOfStr(parts[1], FileSet.EXP) == -1) {
					return false;
				}
				if (!parts[2].startsWith("chr")) {
					return false;
				}
				return true;
			}
		};

		ArrayList<File> allSrcFiles = new ArrayList<>();
		for (File f : subDirs) {
			for (File s : f.listFiles(ff)) {
				allSrcFiles.add(s);
			}
		}

		log.reportTime("Found " + allSrcFiles.size() + " potential source data files.");
		return allSrcFiles;
	}

	private class MarkerLogger {

		public MarkerLogger(Logger log, int minutes) {
			mkrMap = new ConcurrentHashMap<>();
			keyMap = new ConcurrentHashMap<>();
			this.log = log;
			timer = new Timer();
			startLogging(minutes);
		}

		Logger log;
		ConcurrentMap<String, Integer> mkrMap;
		ConcurrentMap<String, int[]> keyMap;

		public void logMarker(int chr, int bStrt, int bLen, int currMkr) {
			String keyStr = chr + "\t" + bStrt + "\t" + bLen;
			synchronized (keyStr) {
				Integer exists = mkrMap.put(keyStr, currMkr);
				if (exists == null) {
					keyMap.put(keyStr, new int[] {chr, bStrt, bLen});
				} else if (exists > currMkr) {
					mkrMap.put(keyStr, exists);
				}
			}
		}

		private void startLogging(int minutes) {
			// wait two minutes, then every N minutes
			timer.schedule(logging, 1000 * 60 * 2, 1000 * 60 * minutes);
		}

		public void stopLogging() {
			timer.cancel();
		}

		Timer timer;
		TimerTask logging = new TimerTask() {
			@Override
			public void run() {
				if (keyMap.size() == 0)
					return;
				long t1 = System.nanoTime();
				List<String> sortedKeys = keyMap.keySet().stream().sorted((o1, o2) -> {
					String[] pts1 = o1.split("\t");
					String[] pts2 = o2.split("\t");
					for (int i = 0; i < pts1.length; i++) {
						int comp1 = Integer.parseInt(pts1[i]);
						int comp2 = Integer.parseInt(pts2[i]);
						int comp = Integer.compare(comp1, comp2);
						if (comp != 0) {
							return comp;
						}
					}
					return 0;
				}).collect(Collectors.toList());

				int max = 0;
				HashMap<Integer, ArrayList<String>> chrKeys = new HashMap<>();
				ArrayList<Integer> keys = new ArrayList<>();
				for (String k : sortedKeys) {
					int[] bin = keyMap.get(k);
					int chr = bin[0];
					double pct = (mkrMap.get(k) - bin[1]) / (double) bin[2];
					// String pctStr = (mkrMap.get(k) - bin[1]) + "/" + bin[2];
					int v = (int) (100 * pct);
					if (!chrKeys.containsKey(chr)) {
						chrKeys.put(chr, new ArrayList<>());
						keys.add(chr);
					}
					// chrKeys.get(chr).add(pct < 0.993 ? v + "% (" + pctStr + ")" : "C");
					chrKeys.get(chr).add(pct < 0.993 ? v + "%" : "C");
					// chrKeys.get(chr).add((mkrMap.get(k) - bin[1]) + "/" + bin[2]);
				}
				for (ArrayList<String> cK : chrKeys.values()) {
					if (cK.size() > max) {
						max = cK.size();
					}
				}
				StringBuilder header = new StringBuilder("CHR");
				for (Integer c : keys) {
					header.append("\t").append(c);
				}

				log.reportTime(header.toString());
				for (int i = 0; i < max; i++) {
					StringBuilder line = new StringBuilder("BIN_").append(i);
					for (Integer c : keys) {
						line.append("\t");
						if (chrKeys.get(c).size() <= i) {
							line.append("-");
						} else {
							line.append(chrKeys.get(c).get(i));
						}
					}
					log.reportTime(line.toString());
				}
				log.reportTime("(wrote to log in " + ext.getTimeElapsedNanos(t1));
				log.reportTime("");

			}
		};

	}

	private class FileReaderSet {
		RandomAccessFile bedIn;
		RandomAccessFile binIn;
		BufferedReader conIn;
		BufferedReader lrrIn;
		BufferedReader bafIn;

		public void open(FileSet fs) throws IOException {
			bedIn = new RandomAccessFile(fs.bedFile, "r");
			binIn = new RandomAccessFile(fs.intFile, "r");
			conIn = new BufferedReader(Files.getAppropriateInputStreamReader(fs.conFile), 4000000);
			lrrIn = new BufferedReader(Files.getAppropriateInputStreamReader(fs.lrrFile), 4000000);
			bafIn = new BufferedReader(Files.getAppropriateInputStreamReader(fs.bafFile), 2000000);
		}

		public void readAhead(FileSet fs, int start, int binSize) throws IOException {
			int line = start + binSize;
			conIn.skip(fs.conInds.get(line) - (start == 0 ? 0 : fs.conInds.get(start)));
			lrrIn.skip(fs.lrrInds.get(line) - (start == 0 ? 0 : fs.lrrInds.get(start)));
			bafIn.skip(fs.bafInds.get(line) - (start == 0 ? 0 : fs.bafInds.get(start)));
		}

		public void close() throws IOException {
			bedIn.close();
			binIn.close();
			conIn.close();
			lrrIn.close();
			bafIn.close();

			bedIn = null;
			binIn = null;
			conIn = null;
			lrrIn = null;
			bafIn = null;
		}
	}

	static class FileSet {
		int chr;

		String bimFile;
		String bedFile;
		String famFile;

		String bafFile;
		String lrrFile;
		String intFile;
		String conFile;

		/** genotype info: */
		String bgiFile;
		String bgenFile;

		HashMap<Integer, Long> conInds;
		HashMap<Integer, Long> lrrInds;
		HashMap<Integer, Long> bafInds;

		List<String[]> mkrBins;

		Logger log;

		/**
		 * Expected files
		 */
		final static String[] EXP = {
																 "l2r",
																 "baf",
																 "con",
																 "int",
																 "snp",
																 "cal",
		};

		public void buildLookups(Logger log) throws IOException {
			long t1 = System.nanoTime();
			ArrayList<Runnable> runners = new ArrayList<>();
			CountDownLatch cdl = new CountDownLatch(3);
			runners.add(new Runnable() {
				@Override
				public void run() {
					try {
						conInds = buildLookup(conFile, log);
					} catch (IOException e) {
						e.printStackTrace();
						System.exit(1);
					}
					cdl.countDown();
				}
			});
			runners.add(new Runnable() {
				@Override
				public void run() {
					try {
						lrrInds = buildLookup(lrrFile, log);
					} catch (IOException e) {
						e.printStackTrace();
						System.exit(1);
					}
					cdl.countDown();
				}
			});
			runners.add(new Runnable() {
				@Override
				public void run() {
					try {
						bafInds = buildLookup(bafFile, log);
					} catch (IOException e) {
						e.printStackTrace();
						System.exit(1);
					}
					cdl.countDown();
				}
			});


			ExecutorService executor = Executors.newFixedThreadPool(Math.min(runners.size(), 3));
			for (Runnable r : runners) {
				executor.submit(r);
			}
			executor.shutdown();
			try {
				cdl.await(Long.MAX_VALUE, TimeUnit.DAYS);
			} catch (InterruptedException e) {
				/**/
			}
			log.reportTime("Built chr" + chr + " lookups in " + ext.getTimeElapsedNanos(t1));
		}

		@SuppressWarnings("unchecked")
		private HashMap<Integer, Long> buildLookup(String file, Logger log) throws IOException {
			long t1 = System.nanoTime();
			String lookupFile = ext.rootOf(file, false) + "_lookup.dat";
			if (Files.exists(lookupFile)) {
				HashMap<Integer, Long> indMap = (HashMap<Integer, Long>) SerializedFiles.readSerial(lookupFile);
				return indMap;
			}
			HashMap<Integer, Long> lineIndices = new HashMap<>();

			InputStreamReader isr = Files.getAppropriateInputStreamReader(file);
			int chr = Integer.MIN_VALUE;
			char[] v;
			long total = 0L;
			int line = 0;
			while ((chr = isr.read()) != -1) {
				v = Character.toChars(chr);
				if (v.length == 1) {
					total++;
					if (v[0] == '\n') {
						line++;
						lineIndices.put(line, total);
					}
				} else {
					for (char element : v) {
						total++;
						if (element == '\n') {
							line++;
							lineIndices.put(line, total);
						}
					}
				}
			}
			isr.close();

			SerializedFiles.writeSerial(lineIndices, lookupFile);
			log.reportTime("Built {" + ext.rootOf(lookupFile) + "} in " + ext.getTimeElapsedNanos(t1));
			return lineIndices;
		}

		boolean isComplete() {
			boolean a = Files.exists(bimFile);
			boolean b = Files.exists(bedFile);
			boolean c = Files.exists(famFile);
			boolean d = Files.exists(bafFile);
			boolean e = Files.exists(lrrFile);
			boolean f = Files.exists(intFile);
			boolean g = Files.exists(conFile);
			if (!a) {
				log.reportTimeWarning("Missing BIM file for chr " + chr);
			}
			if (!b) {
				log.reportTimeWarning("Missing BED file for chr " + chr);
			}
			if (!c) {
				log.reportTimeWarning("Missing FAM file for chr " + chr);
			}
			if (!d) {
				log.reportTimeWarning("Missing BAF file for chr " + chr);
			}
			if (!e) {
				log.reportTimeWarning("Missing LRR file for chr " + chr);
			}
			if (!f) {
				log.reportTimeWarning("Missing INT file for chr " + chr);
			}
			if (!g) {
				log.reportTimeWarning("Missing CON file for chr " + chr);
			}
			return a && b && c && d && e && f && g;
		}
	}

	private static final String ARG_SRC_DIR = "source=";
	private static final String ARG_PROJ_DIR = "projDir=";
	private static final String ARG_PROP_DIR = "propDir=";
	private static final String ARG_PROJ_NAME = "projName=";
	private static final String ARG_FAM_FILE = "fam=";
	private static final String ARG_ANNOT_CSV = "annot=";

	private static final String DESC_SRC_DIR = "Directory of UK BioBank source files";
	private static final String DESC_PROJ_DIR = "Directory in which to create parsed files";
	private static final String DESC_PROP_DIR = "Directory in which to create project properties file";
	private static final String DESC_PROJ_NAME = "Project name";
	private static final String DESC_FAM_FILE = "UKBioBank-provided, PLINK-formatted .fam file";
	private static final String DESC_ANNOT_CSV = "Affymetrix annotation .csv file";


	public static void main(String[] args) {
		String sourceDir = "";
		String projDir = "";
		String propFileDir = "";
		String projName = "";
		String famFile = "";
		String annotFile = "";
		int chrThreads;
		int threadsPerChr;
		int mkrsPerFile;
		double scale;

		CLI cli = new CLI(UKBBParsingPipeline.class);

		cli.addArg(ARG_SRC_DIR, DESC_SRC_DIR, true);
		cli.addArg(ARG_PROJ_DIR, DESC_PROJ_DIR, true);
		cli.addArg(ARG_PROP_DIR, DESC_PROP_DIR, true);
		cli.addArgWithDefault(ARG_PROJ_NAME, DESC_PROJ_NAME, "UKBioBank");
		cli.addArg(ARG_FAM_FILE, DESC_FAM_FILE, true);
		cli.addArg(ARG_ANNOT_CSV, DESC_ANNOT_CSV, false);
		cli.addArg("chrThreads", "Number of chromosome threads", false);
		cli.addArg("threadsPerChr", "Number of threads PER chromosome", false);
		cli.addArg("markersPerFile", "Number of markers per mdRAF file", false);
		cli.addArg("scaleFactor", "Scaling factor for X/Y", false);

		cli.parseWithExit(args);

		sourceDir = cli.get(ARG_SRC_DIR);
		projDir = cli.get(ARG_PROJ_DIR);
		propFileDir = cli.get(ARG_PROP_DIR);
		projName = cli.get(ARG_PROJ_NAME);
		famFile = cli.get(ARG_FAM_FILE);
		if (cli.has(ARG_ANNOT_CSV)) {
			annotFile = cli.get(ARG_ANNOT_CSV);
		}
		chrThreads = cli.has("chrThreads") ? cli.getI("chrThreads") : -1;
		threadsPerChr = cli.has("threadsPerChr") ? cli.getI("threadsPerChr") : -1;
		mkrsPerFile = cli.has("mkrsPerFile") ? cli.getI("mkrsPerFile") : -1;
		scale = cli.has("scaleFactor") ? cli.getD("scaleFactor") : -1;

		UKBBParsingPipeline parser = new UKBBParsingPipeline();
		parser.setSourceDir(sourceDir);
		parser.setProjectDir(projDir);
		parser.setProjectPropertiesDir(propFileDir);
		parser.setProjectName(projName);
		parser.setFamFile(famFile);
		if (!"".equals(annotFile)) {
			parser.setAnnotationCSV(annotFile);
		}
		if (chrThreads != -1) {
			parser.setChromosomeThreadCount(chrThreads);
		}
		if (threadsPerChr != -1) {
			parser.setThreadsPerChromosomeCount(threadsPerChr);
		}
		if (mkrsPerFile != -1) {
			parser.setMarkersPerMDRAF(mkrsPerFile);
		}
		if (scale > 0) {
			parser.setScaleFactor(scale);
		}
		parser.runPipeline();
	}

}
