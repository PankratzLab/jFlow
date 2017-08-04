package org.genvisis.one.ben;

import java.io.BufferedReader;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.genvisis.CLI;
import org.genvisis.cnv.filesys.Compression;
import org.genvisis.cnv.filesys.MarkerData;
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
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.Positions;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.ext;

public class UKBBParsingPipeline {

	Logger log = new Logger();

	String sourceDir;
	String projDir;
	String propFileDir;
	String projName;

	Project proj;
	HashMap<Integer, FileSet> fileSets;
	long fingerprint = 0L;
	String famFile;
	String[][] famData;
	// MkrName + SampName + type -> Value || Write as MkrProjInd + SampName + type
	ConcurrentHashMap<String, Float> oorValues;

	protected void createSampleData() {
		try {
			SampleData.createSampleData(proj.PEDIGREE_FILENAME.getValue(), null, proj);
		} catch (Elision e) {
			e.printStackTrace();
			System.exit(1);
		}
	}

	protected void createPED() {
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
		log.report("Created ped file in " + TimeUnit.NANOSECONDS.toMillis(System.nanoTime() - t1));
	}

	protected void createSampRAFsFromMDRAFs() {
		TransposeData.reverseTranspose(proj);
	}

	protected void writeLookup() {
		Hashtable<String, String> lookup = new Hashtable<>();
		for (FileSet fs : fileSets.values()) {

			int binSize = 2500;

			String mdRAF;
			String[] mkrs = HashVec.loadFileToStringArray(fs.bimFile, false, new int[] {1}, false);
			for (int i = 0; i < mkrs.length; i++) {
				int bin = i / binSize;
				mdRAF = "markers." + fs.chr + "." + bin + "." + (bin + binSize)
								+ MarkerData.MARKER_DATA_FILE_EXTENSION;
				lookup.put(mkrs[i], mdRAF);
			}
			mkrs = null;
		}
		new MarkerLookup(lookup).serialize(proj.MARKERLOOKUP_FILENAME.getValue());
	}

	protected void writeMarkerSet() {
		String[] mkrs = HashVec.loadFileToStringArray(proj.MARKER_POSITION_FILENAME.getValue(), true,
																									new int[] {0}, false);
		Markers.orderMarkers(mkrs, proj.MARKER_POSITION_FILENAME.getValue(),
												 proj.MARKERSET_FILENAME.getValue(), log);
	}

	protected void writeOutliers() {
		Hashtable<String, Float> allOutliers = new Hashtable<>();
		Map<String, Integer> mkrInds = proj.getMarkerIndices();
		for (Entry<String, Float> ent : oorValues.entrySet()) {
			String[] pts = ent.getKey().split("\t");
			allOutliers.put(mkrInds.get(pts[0]) + "\t" + pts[1] + "\t" + pts[0], ent.getValue());
		}
		oorValues = null;

		SerializedFiles.writeSerial(allOutliers, proj.MARKER_DATA_DIRECTORY.getValue(true, true)
																						 + "outliers.ser");
	}

	private byte[] compressMarkerData(MarkerData md, int mkrIndLocal,
																		Hashtable<String, Float> oorTable) throws Elision {
		int nInd = md.getGCs().length;
		int bytesPerSamp = Sample.getNBytesPerSampleMarker(getNullStatus());
		int markerBlockSize = nInd * bytesPerSamp;

		byte[] mkrBuff = new byte[markerBlockSize];
		int buffInd = 0;

		// GCs
		for (int i = 0; i < nInd; i++) {
			Compression.gcBafCompress(md.getGCs()[i], mkrBuff, buffInd);
			buffInd += Compression.REDUCED_PRECISION_GCBAF_NUM_BYTES;
		}

		// Xs
		for (int i = 0; i < nInd; i++) {
			boolean oor = !Compression.xyCompressAllowNegative(md.getXs()[i], mkrBuff, buffInd);
			buffInd += Compression.REDUCED_PRECISION_XY_NUM_BYTES;
			if (oor) {
				oorTable.put(mkrIndLocal + "\t" + i + "\tx", md.getXs()[i]);
			}
		}

		// Ys
		for (int i = 0; i < nInd; i++) {
			boolean oor = !Compression.xyCompressAllowNegative(md.getYs()[i], mkrBuff, buffInd);
			buffInd += Compression.REDUCED_PRECISION_XY_NUM_BYTES;
			if (oor) {
				oorTable.put(mkrIndLocal + "\t" + i + "\ty", md.getYs()[i]);
			}
		}

		// BAFs
		for (int i = 0; i < nInd; i++) {
			Compression.gcBafCompress(md.getBAFs()[i], mkrBuff, buffInd);
			buffInd += Compression.REDUCED_PRECISION_GCBAF_NUM_BYTES;
		}

		// LRRs
		for (int i = 0; i < nInd; i++) {
			boolean oor = -1 == Compression.lrrCompress(md.getLRRs()[i], mkrBuff, buffInd);
			buffInd += Compression.REDUCED_PRECISION_LRR_NUM_BYTES;
			if (oor) {
				oorTable.put(mkrIndLocal + "\t" + i + "\tlrr", md.getLRRs()[i]);
			}
		}

		// Genotypes
		for (int i = 0; i < nInd; i++) {
			mkrBuff[buffInd] = Compression.genotypeCompress(md.getAbGenotypes() == null
																																								 ? -1
																																								 : md.getAbGenotypes()[i],
																											md.getForwardGenotypes() == null
																																											? 0
																																											: md.getForwardGenotypes()[i]);
			buffInd += Compression.REDUCED_PRECISION_ABFORWARD_GENOTYPE_NUM_BYTES;
		}

		return mkrBuff;
	}

	private byte getNullStatus() {
		return Sample.updateNullStatus(new float[0], new float[0],
																	 new float[0], new float[0],
																	 new float[0], new byte[0],
																	 null, true);
	}


	private RandomAccessFile openMDRAF(String filename, int nInd, int nMkr, byte nullStatus,
																		 long fingerprint,
																		 String[] mkrNames) throws IOException {
		byte[] mkrBytes = Compression.objToBytes(mkrNames);
		byte[] mdRAFHeader = TransposeData.getParameterSectionForMdRaf(nInd,
																																	 nMkr,
																																	 nullStatus,
																																	 fingerprint,
																																	 mkrBytes);
		mkrBytes = null;

		RandomAccessFile mdRAF = new RandomAccessFile(proj.MARKER_DATA_DIRECTORY.getValue(true, true)
																									+ filename, "rw");
		mdRAF.write(mdRAFHeader);
		mdRAFHeader = null;

		return mdRAF;
	}

	private void createMDRAF(FileSet fs) throws IOException, Elision {
		String[][] bimData = HashVec.loadFileToStringMatrix(fs.bimFile, false, new int[] {0, 1, 2, 3,
																																											4, 5}, false);

		int nInd = famData.length;

		byte nullStatus = getNullStatus();

		FileReaderSet frs = new FileReaderSet();
		frs.bedIn = new RandomAccessFile(fs.bedFile, "r");
		frs.binIn = new RandomAccessFile(fs.intFile, "r");
		frs.conIn = new BufferedReader(Files.getAppropriateInputStreamReader(fs.conFile), 4000000);
		frs.lrrIn = new BufferedReader(Files.getAppropriateInputStreamReader(fs.lrrFile), 4000000);
		frs.bafIn = new BufferedReader(Files.getAppropriateInputStreamReader(fs.bafFile), 2000000);

		byte[] magicBytes = new byte[3];
		frs.bedIn.read(magicBytes);
		if (magicBytes[2] == 0) {
			log.reportError("Error - .bed file is sample-dominant: " + fs.bedFile);
			System.exit(1);
		}

		int binSize = 2500;
		String[] mkrNames = Matrix.extractColumn(bimData, 1);
		List<String[]> mkrBins = ArrayUtils.splitUpArray(mkrNames, (mkrNames.length / binSize) + 1, log);
		String[] existing = readWritten();

		for (int i = 0; i < mkrBins.size(); i++) {
			String mdRAFName = "markers." + fs.chr + "." + (i * binSize) + "."
												 + ((i * binSize) + mkrBins.get(i).length)
												 + MarkerData.MARKER_DATA_FILE_EXTENSION;
			if (ext.indexOfStr(mdRAFName, existing) == -1) {
				String mdRAF = write(fs, i * binSize, mkrBins.get(i), nullStatus, nInd, bimData,
														 frs);
				success(mdRAF);
			} else {
				frs.readAhead(fs, i, binSize);
			}
		}

		frs.bedIn.close();
		frs.binIn.close();
		frs.conIn.close();
		frs.lrrIn.close();
		frs.bafIn.close();

		frs.bedIn = null;
		frs.binIn = null;
		frs.conIn = null;
		frs.lrrIn = null;
		frs.bafIn = null;
		bimData = null;
		frs = null;

		log.reportTime("Wrote MDRAF for Chr " + fs.chr);

		System.gc();
	}

	private String[] readWritten() {
		String temp = proj.PROJECT_DIRECTORY.getValue() + "markersWritten.txt";
		return Files.exists(temp) ? HashVec.loadFileToStringArray(temp, false, null, false)
														 : new String[0];
	}

	private void success(String mkrFile) {
		String temp = proj.PROJECT_DIRECTORY.getValue() + "markersWritten.txt";
		Files.appendStringToFile(temp, mkrFile);
	}

	private static class FileReaderSet {
		RandomAccessFile bedIn;
		RandomAccessFile binIn;
		BufferedReader conIn;
		BufferedReader lrrIn;
		BufferedReader bafIn;

		public void readAhead(FileSet fs, int start, int binSize) throws IOException {
			int line = start + binSize;
			conIn.skip(fs.conInds.get(line) - fs.conInds.get(start));
			lrrIn.skip(fs.lrrInds.get(line) - fs.lrrInds.get(start));
			bafIn.skip(fs.bafInds.get(line) - fs.bafInds.get(start));
		}
	}

	private String write(FileSet fs, int startBatchInd, String[] mkrNames, byte nullStatus, int nInd,
											 String[][] bimData, FileReaderSet frs) throws Elision, IOException {
		RandomAccessFile mdRAF;

		Hashtable<String, Float> outOfRangeTable = new Hashtable<>();
		String mdRAFName = "markers." + fs.chr + "." + startBatchInd + "."
											 + (startBatchInd + mkrNames.length) + MarkerData.MARKER_DATA_FILE_EXTENSION;

		mdRAF = openMDRAF(mdRAFName, famData.length, mkrNames.length, nullStatus, fingerprint, mkrNames);

		int bedBlockSize = (int) Math.ceil(nInd / 4.0);
		int binBlockSize = nInd * 8;

		int logPer = mkrNames.length / 10;
		String conLine;
		String[] cons;
		float[] gcs;
		String lrrLine;
		String[] lrrStrs;
		float[] lrrs;
		String bafLine;
		String[] bafStrs;
		float[] bafs;
		for (int i = startBatchInd; i < startBatchInd + mkrNames.length; i++) {
			if (i > 0 && i % logPer == 0) {
				log.reportTime("Processed " + i + " markers for chr " + fs.chr);
			}

			String[] bimLine = bimData[i];
			frs.bedIn.seek(3L + ((long) i) * ((long) bedBlockSize)); // ALWAYS seek forward
			frs.binIn.seek(((long) i) * ((long) binBlockSize));

			byte[] markerBytes = new byte[bedBlockSize];
			byte[] intensBytes = new byte[binBlockSize];
			byte[] sampGeno = new byte[nInd];
			float[] sampAs = new float[nInd];
			float[] sampBs = new float[nInd];
			frs.bedIn.read(markerBytes);
			frs.binIn.read(intensBytes);

			for (int bitInd = 0, binInd = 0; bitInd < markerBytes.length; bitInd++) {
				byte bedByte = markerBytes[bitInd];
				byte[] genotypes = PlinkData.decodeBedByte(bedByte);

				for (int g = 0; g < genotypes.length; g++) {
					int idInd = bitInd * 4 + g;
					if (idInd == -1 || idInd >= sampGeno.length) {
						continue;
					}
					sampGeno[idInd] = genotypes[g];
				}
				genotypes = null;

				byte[] intA = {intensBytes[binInd++], intensBytes[binInd++], intensBytes[binInd++],
											 intensBytes[binInd++]};
				byte[] intB = {intensBytes[binInd++], intensBytes[binInd++], intensBytes[binInd++],
											 intensBytes[binInd++]};
				sampAs[binInd / 8] = ByteBuffer.wrap(intA).order(ByteOrder.LITTLE_ENDIAN).getFloat();
				sampBs[binInd / 8] = ByteBuffer.wrap(intB).order(ByteOrder.LITTLE_ENDIAN).getFloat();
			}
			markerBytes = null;
			intensBytes = null;

			conLine = frs.conIn.readLine();
			cons = conLine.split(" ");
			gcs = new float[nInd];
			if (cons.length != nInd) {
				log.reportError("Mismatched # of confidence scores {fnd: " + cons.length + ", exp: "
												+ nInd + "}.  File: " + fs.conFile);
				System.exit(1);
			}
			for (int k = 0; k < cons.length; k++) {
				gcs[k] = 1 - Float.parseFloat(cons[k]);
				if (gcs[k] == 2) { // i.e. 1 - (-1)
					gcs[k] = Float.NaN;
				}
			}
			conLine = null;
			cons = null;

			lrrLine = frs.lrrIn.readLine();
			lrrStrs = lrrLine.split(" ");
			lrrs = new float[nInd];

			if (lrrStrs.length != nInd) {
				log.reportError("Mismatched # of L2R values {fnd: " + lrrStrs.length + ", exp: "
												+ nInd + "}.  File: " + fs.lrrFile);
				System.exit(1);
			}
			for (int k = 0; k < lrrStrs.length; k++) {
				try {
					lrrs[k] = Float.parseFloat(lrrStrs[k]);
				} catch (NumberFormatException e) {
					if (ext.isMissingValue(lrrStrs[k])) {
						lrrs[k] = Float.NaN;
					} else {
						log.reportError("Malformed LRR value: " + lrrStrs[k] + " for chr " + fs.chr);
						System.exit(1);
					}
				}
			}
			lrrLine = null;
			lrrStrs = null;

			bafLine = frs.bafIn.readLine();
			bafStrs = bafLine.split(" ");
			bafs = new float[nInd];
			if (bafStrs.length != nInd) {
				log.reportError("Mismatched # of BAF values {fnd: " + bafStrs.length + ", exp: "
												+ nInd + "}.  File: " + fs.bafFile);
				System.exit(1);
			}
			for (int k = 0; k < bafStrs.length; k++) {
				try {
					bafs[k] = Float.parseFloat(bafStrs[k]);
				} catch (NumberFormatException e) {
					if (ext.isMissingValue(bafStrs[k])) {
						bafs[k] = Float.NaN;
					} else {
						log.reportError("Malformed BAF value: " + bafStrs[k] + " for chr " + fs.chr);
						System.exit(1);
					}
				}
			}
			bafLine = null;
			bafStrs = null;

			MarkerData md = new MarkerData(bimLine[1], (byte) Integer.parseInt(bimLine[0]),
																		 Integer.parseInt(bimLine[3]), fingerprint,
																		 gcs, null, null, sampAs, sampBs, null, null, bafs, lrrs,
																		 sampGeno, null);
			byte[] writeData = compressMarkerData(md, i - startBatchInd, outOfRangeTable);

			mdRAF.write(writeData);

			writeData = null;
			md = null;
			gcs = null;
			sampAs = null;
			sampBs = null;
			bafs = null;
			lrrs = null;
			sampGeno = null;

			if (i % 500 == 0) {
				System.out.print("Cleaning up...");
				long t1 = System.nanoTime();
				System.gc();
				System.out.println("done in " + TimeUnit.NANOSECONDS.toMillis(System.nanoTime() - t1)
													 + " millis");
			}
		}

		byte[] oorBytes = Compression.objToBytes(outOfRangeTable);
		mdRAF.write(Compression.intToBytes(oorBytes.length));
		mdRAF.write(oorBytes);

		mdRAF.close();

		for (Entry<String, Float> entry : outOfRangeTable.entrySet()) {
			String[] pts = entry.getKey().split("\t");
			int mkrInd = Integer.parseInt(pts[0]);
			int sampInd = Integer.parseInt(pts[1]);
			oorValues.put(bimData[startBatchInd + mkrInd][1] + "\t" + famData[sampInd][1] + "\t" + pts[2],
										entry.getValue());
		}

		oorBytes = null;
		outOfRangeTable = null;
		return mdRAFName;
	}

	protected void parseMarkerData() {
		long time = System.nanoTime();
		oorValues = new ConcurrentHashMap<String, Float>();

		int avail = Runtime.getRuntime().availableProcessors();
		ExecutorService executor = Executors.newFixedThreadPool(Math.min(fileSets.size(),
																																		 Math.max(1, avail - 1)));

		CountDownLatch cdl = new CountDownLatch(fileSets.size());

		for (Integer v : fileSets.keySet()) {
			final FileSet fs = fileSets.get(v);
			final Thread parentThread = Thread.currentThread();
			Runnable runn = new Runnable() {
				@Override
				public void run() {
					try {
						fs.readLookups();
						createMDRAF(fs);
						cdl.countDown();
						parentThread.interrupt();
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
			log.reportTime("Queued MarkerData parsing for chr " + fs.chr);
		}

		executor.shutdown();
		try {
			cdl.await(Long.MAX_VALUE, TimeUnit.DAYS);
		} catch (InterruptedException e) {
			/**/
		}

		long elapsed = TimeUnit.NANOSECONDS.toSeconds(System.nanoTime() - time);
		log.report("Completed parsing MarkerData in " + elapsed + " seconds.");
	}

	private void completedMarkerPositions(long time) {
		long elapsed = TimeUnit.NANOSECONDS.toSeconds(System.nanoTime() - time);
		log.report("Completed markerPositions: " + proj.MARKER_POSITION_FILENAME.getValue() + " in "
							 + elapsed + " seconds.");
	}

	protected void createMarkerPositions() {
		String file = proj.MARKER_POSITION_FILENAME.getValue();
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

		completedMarkerPositions(time);
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

		final static String[] EXP = {
																 "l2r",
																 "baf",
																 "con",
																 "int",
																 "snp",
																 "cal",
		};

		public void buildLookups() throws IOException {
			conInds = buildLookup(conFile, 4000000);
			lrrInds = buildLookup(lrrFile, 4000000);
			bafInds = buildLookup(bafFile, 2000000);
		}

		@SuppressWarnings("unchecked")
		public void readLookups() {
			if (conInds == null) {
				conInds = (HashMap<Integer, Long>) SerializedFiles.readSerial(ext.rootOf(conFile, false)
																																			+ "_lookup.dat");
			}
			if (lrrInds == null) {
				lrrInds = (HashMap<Integer, Long>) SerializedFiles.readSerial(ext.rootOf(lrrFile, false)
																																			+ "_lookup.dat");
			}
			if (bafInds == null) {
				bafInds = (HashMap<Integer, Long>) SerializedFiles.readSerial(ext.rootOf(bafFile, false)
																																			+ "_lookup.dat");
			}
		}

		@SuppressWarnings("unchecked")
		private HashMap<Integer, Long> buildLookup(String file, int buffer) throws IOException {
			if (Files.exists(ext.rootOf(file, false) + "_lookup.dat"))
				return (HashMap<Integer, Long>) SerializedFiles.readSerial(ext.rootOf(file, false)
																																	 + "_lookup.dat");
			long t1 = System.nanoTime();
			BufferedReader reader = new BufferedReader(Files.getAppropriateInputStreamReader(file),
																								 buffer);
			HashMap<Integer, Long> lineIndices = new HashMap<>();
			int ind = 1;
			String line;
			long total = 0L;
			while ((line = reader.readLine()) != null) {
				total += line.length();
				lineIndices.put(ind, total);
				ind++;
			}
			reader.close();
			SerializedFiles.writeSerial(lineIndices, ext.rootOf(file, false) + "_lookup.dat");
			System.out.println("Built lookup for file {" + file + "} in "
												 + TimeUnit.NANOSECONDS.toSeconds(System.nanoTime() - t1) + " seconds");
			return lineIndices;
		}

		boolean isComplete() {
			return Files.exists(bimFile) && Files.exists(bedFile) && Files.exists(famFile)
						 && Files.exists(bafFile)
						 && Files.exists(lrrFile) && Files.exists(intFile) && Files.exists(conFile);
		}

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
				if (parts.length <= 2)
					return false;
				if (name.endsWith(".fam"))
					return false;
				if (!parts[0].equals("ukb"))
					return false;
				if (ext.indexOfStr(parts[1], FileSet.EXP) == -1)
					return false;
				if (!parts[2].startsWith("chr"))
					return false;
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
		for (int i = 0; i < 27; i++) {
			if (!fileSets.containsKey(i)) {
				missingChrs.add(i);
			} else if (!fileSets.get(i).isComplete()) {
				log.reportTimeWarning("Removing chr " + fileSets.get(i).chr + " - incomplete file set.");
				fileSets.remove(i);
			} else {
				try {
					fileSets.get(i).buildLookups();
				} catch (IOException e) {
					e.printStackTrace();
					System.exit(1);
				}
			}
		}
		if (missingChrs.size() > 0) {
			log.reportTimeWarning("Missing " + missingChrs.size() + " chrs: "
														+ ArrayUtils.toStr(missingChrs, ", "));
		}
	}

	protected void createSampleList() {
		famData = HashVec.loadFileToStringMatrix(famFile, false, new int[] {0, 1, 2, 3,
																																				4, 5}, false);

		String[] allSamples = Matrix.extractColumn(famData, 1);
		@SuppressWarnings("deprecation")
		long f = org.genvisis.cnv.filesys.MarkerSet.fingerprint(allSamples);
		fingerprint = f;
		SampleList sl = new SampleList(allSamples);
		sl.serialize(proj.SAMPLELIST_FILENAME.getValue());
		sl = null;
		allSamples = null;
	}

	protected void createProject() {
		String propFile = ext.verifyDirFormat(propFileDir)
											+ ext.replaceWithLinuxSafeCharacters(projName) + ".properties";
		Files.write((new Project()).PROJECT_NAME.getName() + "=" + projName, propFile);
		proj = new Project(propFile, false);
		proj.PROJECT_NAME.setValue(projName);
		proj.PROJECT_DIRECTORY.setValue(projDir);
		proj.SOURCE_DIRECTORY.setValue(sourceDir);
		proj.XY_SCALE_FACTOR.setValue(1d);
		proj.TARGET_MARKERS_FILENAMES.setValue(new String[] {});
		proj.SOURCE_FILENAME_EXTENSION.setValue("NULL");
		proj.ID_HEADER.setValue("NULL");
		proj.SOURCE_FILENAME_EXTENSION.setValue("NULL");
		proj.GENOME_BUILD_VERSION.setValue(GENOME_BUILD.HG19);

		proj.ARRAY_TYPE.setValue(ARRAY.AFFY_AXIOM);

		proj.saveProperties();
	}

	public void runPipeline() {
		createProject();
		createSampleList();
		discoverFilesets();
		createMarkerPositions();
		parseMarkerData();
		writeLookup();
		writeMarkerSet();
		writeOutliers();
		createSampRAFsFromMDRAFs();
		createPED();
		createSampleData();
	}

	private void setProjectName(String projName2) {
		this.projName = projName2;
	}

	private void setProjectPropertiesDir(String propFileDir2) {
		this.propFileDir = propFileDir2;
	}

	private void setProjectDir(String projDir2) {
		this.projDir = projDir2;
	}

	private void setSourceDir(String sourceDir2) {
		this.sourceDir = sourceDir2;
	}

	private void setFamFile(String famFile2) {
		this.famFile = famFile2;
	}

	public static void run() {
		String dir = "F:/testProjectSrc/UKBB_AffyAxiom/";
		String sourceDir = dir + "00src/";
		String projDir = dir + "project/";
		String propFileDir = "D:/projects/";
		String projName = "UKBioBank";
		String famFile = dir + "ukb1773_l2r_chrY_v2_s488374.fam";

		if (!Files.isWindows()) {
			dir = "/scratch.global/cole0482/UKBB/";
			sourceDir = dir + "00src/";
			projDir = dir + "project/";
			propFileDir = "/home/pankrat2/cole0482/projects/";
			famFile = dir + "ukb1773_l2r_chrY_v2_s488374.fam";
		}

		UKBBParsingPipeline parser = new UKBBParsingPipeline();
		parser.setSourceDir(sourceDir);
		parser.setProjectDir(projDir);
		parser.setProjectPropertiesDir(propFileDir);
		parser.setProjectName(projName);
		parser.setFamFile(famFile);
		parser.runPipeline();
	}

	private static final String ARG_SRC_DIR = "source=";
	private static final String ARG_PROJ_DIR = "projDir=";
	private static final String ARG_PROP_DIR = "propFile=";
	private static final String ARG_PROJ_NAME = "projName=";
	private static final String ARG_FAM_FILE = "fam=";

	private static final String DESC_SRC_DIR = "Directory of UK BioBank source files";
	private static final String DESC_PROJ_DIR = "Directory in which to create parsed files";
	private static final String DESC_PROP_DIR = "Directory in which to create project properties file";
	private static final String DESC_PROJ_NAME = "Project name";
	private static final String DESC_FAM_FILE = "UKBioBank-provided, PLINK-formatted .fam file";

	public static void main(String[] args) {
		String sourceDir = "";
		String projDir = "";
		String propFileDir = "";
		String projName = "";
		String famFile = "";

		boolean testing = true; // TODO REMOVE
		if (testing) {
			run();
			return;
		}

		CLI cli = new CLI(UKBBParsingPipeline.class);

		cli.addArg(ARG_SRC_DIR, DESC_SRC_DIR, true);
		cli.addArg(ARG_PROJ_DIR, DESC_PROJ_DIR, true);
		cli.addArg(ARG_PROP_DIR, DESC_PROP_DIR, true);
		cli.addArgWithDefault(ARG_PROJ_NAME, DESC_PROJ_NAME, "UKBioBank");
		cli.addArg(ARG_FAM_FILE, DESC_FAM_FILE, true);

		cli.parseWithExit(args);

		sourceDir = cli.get(ARG_SRC_DIR);
		projDir = cli.get(ARG_PROJ_DIR);
		propFileDir = cli.get(ARG_PROP_DIR);
		projName = cli.get(ARG_PROJ_NAME);
		famFile = cli.get(ARG_FAM_FILE);

		UKBBParsingPipeline parser = new UKBBParsingPipeline();
		parser.setSourceDir(sourceDir);
		parser.setProjectDir(projDir);
		parser.setProjectPropertiesDir(propFileDir);
		parser.setProjectName(projName);
		parser.setFamFile(famFile);
		parser.runPipeline();
	}

}
