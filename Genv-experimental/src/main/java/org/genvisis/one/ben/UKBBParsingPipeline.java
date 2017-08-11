package org.genvisis.one.ben;

import java.io.BufferedReader;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.math.RoundingMode;
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

	private static final int MAX_MDRAF_THREADS = 3;
	private static final int MAX_CHR_THREADS = 8;

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

	public void runPipeline() {
		createProject();
		createSampleList();
		discoverFilesets();
		parseMarkerData();
		writeLookup();
		createMarkerPositions();
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
		log.reportTime("Created ped file in " + TimeUnit.NANOSECONDS.toMillis(System.nanoTime() - t1));
	}

	protected void createSampRAFsFromMDRAFs() {
		TransposeData.reverseTranspose(proj);
	}

	protected void writeLookup() {
		Hashtable<String, String> lookup = new Hashtable<>();
		for (FileSet fs : fileSets.values()) {
			String mdRAF;
			List<String[]> mkrBins = fs.mkrBins;
			int start = 0;
			int end = 0;
			for (String[] bin : mkrBins) {
				end = start + bin.length;
				mdRAF = getMDRAFName(fs.chr, start, end);
				for (String s : bin) {
					lookup.put(s, mdRAF);
				}
				start = end;
			}
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

	private byte getNullStatus() {
		return Sample.updateNullStatus(new float[0], new float[0],
																	 new float[0], new float[0],
																	 new float[0], null,
																	 new byte[0], true);
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

		RandomAccessFile mdRAF = new RandomAccessFile(proj.MARKER_DATA_DIRECTORY.getValue(true, true)
																									+ filename, "rw");
		mdRAF.write(mdRAFHeader);
		mdRAFHeader = null;

		return mdRAF;
	}

	private void createMDRAF(FileSet fs) throws IOException, Elision {
		final String[][] bimData = HashVec.loadFileToStringMatrix(fs.bimFile, false, new int[] {0, 1, 2, 3,
																																											4, 5}, false);

		final int nInd = famData.length;

		final byte nullStatus = getNullStatus();

		String[] mkrNames = Matrix.extractColumn(bimData, 1);
		fs.mkrBins = ArrayUtils.splitUpArray(mkrNames, (mkrNames.length / 2500) + 1, log);
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
							mdRAF = write(fs, s, fs.mkrBins.get(ii), nullStatus, nInd, bimData);
							log.reportTime("Wrote " + mdRAF + " in " + TimeUnit.NANOSECONDS.toHours(System.nanoTime() - t1) + " hours.");
							success(mdRAF);
							cdl.countDown();
						} catch (Elision e) {
							e.printStackTrace();
						} catch (IOException e) {
							e.printStackTrace();
						}
					}
				});
			} else {
				runners.add(new Runnable() {
					@Override
					public void run() {
						Hashtable<String, Float> outliers = TransposeData.loadOutliersFromRAF(proj.MARKER_DATA_DIRECTORY.getValue() + mdRAFName);
						for (Entry<String, Float> outlier : outliers.entrySet()) {
							String[] pts = outlier.getKey().split("\t");
							int mkrInd = Integer.parseInt(pts[0]);
							int sampInd = Integer.parseInt(pts[1]);
							oorValues.put(bimData[s + mkrInd][1] + "\t" + famData[sampInd][1] + "\t" + pts[2],
														outlier.getValue());
						}
						cdl.countDown();
					}
				});
			}
			start = end;
		}

		if (runners.size() > 0) {
  		ExecutorService executor = Executors.newFixedThreadPool(Math.min(runners.size(), MAX_MDRAF_THREADS));
  		for (Runnable r : runners) {
  			executor.submit(r);
  		}
  		executor.shutdown();
  		try {
  			cdl.await(Long.MAX_VALUE, TimeUnit.DAYS);
  		} catch (InterruptedException e) {
  			/**/
  		}
  
  		System.gc();
		}
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

		public void open(FileSet fs) throws IOException {
			this.bedIn = new RandomAccessFile(fs.bedFile, "r");
			this.binIn = new RandomAccessFile(fs.intFile, "r");
			this.conIn = new BufferedReader(Files.getAppropriateInputStreamReader(fs.conFile), 4000000);
			this.lrrIn = new BufferedReader(Files.getAppropriateInputStreamReader(fs.lrrFile), 4000000);
			this.bafIn = new BufferedReader(Files.getAppropriateInputStreamReader(fs.bafFile), 2000000);
		}
		
		public void readAhead(FileSet fs, int start, int binSize) throws IOException {
			int line = start + binSize;
			conIn.skip(fs.conInds.get(line) - fs.conInds.get(start));
			lrrIn.skip(fs.lrrInds.get(line) - fs.lrrInds.get(start));
			bafIn.skip(fs.bafInds.get(line) - fs.bafInds.get(start));
		}

		public void close() throws IOException {
			this.bedIn.close();
			this.binIn.close();
			this.conIn.close();
			this.lrrIn.close();
			this.bafIn.close();

			this.bedIn = null;
			this.binIn = null;
			this.conIn = null;
			this.lrrIn = null;
			this.bafIn = null;
		}
	}

	private String getMDRAFName(int chr, int start, int end) {
		return "markers." + chr + "." + start + "."
					 + end + MarkerData.MARKER_DATA_FILE_EXTENSION;
	}

	private String write(FileSet fs, int startBatchInd, String[] mkrNames, byte nullStatus, int nInd,
											 String[][] bimData) throws Elision, IOException {

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

		int logPer = mkrNames.length / 10;
		String conLine;
		String[] cons;
		String lrrLine;
		String[] lrrStrs;
		String bafLine;
		String[] bafStrs;
		for (int i = startBatchInd; i < startBatchInd + mkrNames.length; i++) {
			if (i > 0 && i % logPer == 0) {
				log.reportTime("Processed " + (10 * (i / logPer)) + "% of markers for chr " + fs.chr);
			}

			int bytesPerSamp = Sample.getNBytesPerSampleMarker(getNullStatus());
			int markerBlockSize = nInd * bytesPerSamp;

			byte[] mkrBuff = new byte[markerBlockSize];
			int buffInd = 0;

			conLine = frs.conIn.readLine();
			cons = conLine.split(" ");
			conLine = null;
			if (cons.length != nInd) {
				log.reportError("Mismatched # of confidence scores {fnd: " + cons.length + ", exp: "
												+ nInd + "}.  File: " + fs.conFile);
				System.exit(1);
			}
			float gcV;
			for (int k = 0; k < cons.length; k++) {
				gcV = 1 - Float.parseFloat(cons[k]);
				if (gcV == 2) { // i.e. 1 - (-1)
					gcV = Float.NaN;
				}
				Compression.gcBafCompress(gcV, mkrBuff, buffInd);
				buffInd += Compression.REDUCED_PRECISION_GCBAF_NUM_BYTES;
			}
			cons = null;

			byte[] intensBytes = new byte[binBlockSize];
			frs.binIn.seek(((long) i) * ((long) binBlockSize));
			frs.binIn.read(intensBytes);

			float x, y;
			boolean oor;
			for (int bitInd = 0, binInd = 0; bitInd < nInd; bitInd++) {
				byte[] intA = {intensBytes[binInd++], intensBytes[binInd++], intensBytes[binInd++],
											 intensBytes[binInd++]};
				byte[] intB = {intensBytes[binInd++], intensBytes[binInd++], intensBytes[binInd++],
											 intensBytes[binInd++]};

				x = ByteBuffer.wrap(intA).order(ByteOrder.LITTLE_ENDIAN).getFloat();
				y = ByteBuffer.wrap(intB).order(ByteOrder.LITTLE_ENDIAN).getFloat();
//				x = LongMath.log2((long) (x / y), RoundingMode.HALF_UP);
//				y = LongMath.log2((long) (x * y), RoundingMode.HALF_UP) / 2;
//				x = (float) (Math.log(x) / Math.log(2));
//				y = (float) (Math.log(y) / Math.log(2));
				
				oor = !Compression.xyCompressAllowNegative(x, mkrBuff, buffInd);
				if (oor) {
					outOfRangeTable.put((i - startBatchInd) + "\t" + i + "\tx", x);
				}

				oor = !Compression.xyCompressAllowNegative(y,
																									 mkrBuff,
																									 buffInd
																											 + (nInd * Compression.REDUCED_PRECISION_XY_NUM_BYTES));
				if (oor) {
					outOfRangeTable.put((i - startBatchInd) + "\t" + i + "\ty", y);
				}

				buffInd += Compression.REDUCED_PRECISION_XY_NUM_BYTES;
			}
			intensBytes = null;

			bafLine = frs.bafIn.readLine();
			bafStrs = bafLine.split(" ");
			bafLine = null;
			if (bafStrs.length != nInd) {
				log.reportError("Mismatched # of BAF values {fnd: " + bafStrs.length + ", exp: "
												+ nInd + "}.  File: " + fs.bafFile);
				System.exit(1);
			}
			float baf = Float.NaN;
			for (int k = 0; k < bafStrs.length; k++) {
				try {
					baf = Float.parseFloat(bafStrs[k]);
				} catch (NumberFormatException e) {
					if (ext.isMissingValue(bafStrs[k])) {
						baf = Float.NaN;
					} else {
						log.reportError("Malformed BAF value: " + bafStrs[k] + " for chr " + fs.chr);
						System.exit(1);
					}
				}
				Compression.gcBafCompress(baf, mkrBuff, buffInd);
				buffInd += Compression.REDUCED_PRECISION_GCBAF_NUM_BYTES;
			}
			bafStrs = null;

			lrrLine = frs.lrrIn.readLine();
			lrrStrs = lrrLine.split(" ");
			lrrLine = null;
			if (lrrStrs.length != nInd) {
				log.reportError("Mismatched # of L2R values {fnd: " + lrrStrs.length + ", exp: "
												+ nInd + "}.  File: " + fs.lrrFile);
				System.exit(1);
			}
			float lrr = Float.NaN;
			for (int k = 0; k < lrrStrs.length; k++) {
				try {
					lrr = Float.parseFloat(lrrStrs[k]);
				} catch (NumberFormatException e) {
					if (ext.isMissingValue(lrrStrs[k])) {
						lrr = Float.NaN;
					} else {
						log.reportError("Malformed LRR value: " + lrrStrs[k] + " for chr " + fs.chr);
						System.exit(1);
					}
				}
				oor = -1 == Compression.lrrCompress(lrr, mkrBuff, buffInd);
				buffInd += Compression.REDUCED_PRECISION_LRR_NUM_BYTES;
				if (oor) {
					outOfRangeTable.put((i - startBatchInd) + "\t" + i + "\tlrr", lrr);
				}
			}
			lrrStrs = null;

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
//					mkrBuff[buffInd] = Compression.genotypeCompress((byte) -1, genotypes[g]);
					mkrBuff[buffInd] = Compression.genotypeCompress(genotypes[g], (byte) 0);
					buffInd += Compression.REDUCED_PRECISION_ABFORWARD_GENOTYPE_NUM_BYTES;
				}
				genotypes = null;
			}
			markerBytes = null;

			mdRAF.write(mkrBuff);
			mkrBuff = null;

		}

		byte[] oorBytes = Compression.objToBytes(outOfRangeTable);
		mdRAF.write(Compression.intToBytes(oorBytes.length));
		mdRAF.write(oorBytes);

		mdRAF.close();
		
		frs.close();
		frs = null;
		
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

		ExecutorService executor = Executors.newFixedThreadPool(Math.min(fileSets.size(), MAX_CHR_THREADS));

		CountDownLatch cdl = new CountDownLatch(fileSets.size());

		for (Integer v : fileSets.keySet()) {
			final FileSet fs = fileSets.get(v);
			Runnable runn = new Runnable() {
				@Override
				public void run() {
					try {
						fs.buildLookups(log);
						createMDRAF(fs);
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

		long elapsed = TimeUnit.NANOSECONDS.toDays(System.nanoTime() - time);
		log.reportTime("Completed parsing MarkerData in " + elapsed + " days.");
	}

	private void completedMarkerPositions(long time) {
		long elapsed = TimeUnit.NANOSECONDS.toSeconds(System.nanoTime() - time);
		log.reportTime("Completed markerPositions: " + proj.MARKER_POSITION_FILENAME.getValue() + " in "
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

		List<String[]> mkrBins;

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
			log.reportTime("Built chr" + chr + " lookups in " + TimeUnit.NANOSECONDS.toMinutes(System.nanoTime() - t1) + " minutes.");
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
					for (int i = 0; i < v.length; i++) {
						total++;
						if (v[i] == '\n') {
							line++;
							lineIndices.put(line, total);
						}
					}
				}
			}
			isr.close();

			SerializedFiles.writeSerial(lineIndices, lookupFile);
			log.reportTime("Built {" + ext.rootOf(lookupFile) + "} in "
												 + TimeUnit.NANOSECONDS.toMinutes(System.nanoTime() - t1) + " minutes");
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
				if (parts.length <= 2 || parts.length > 4)
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
			log.reportTimeWarning("Removed chrs with incomplete files: " + ArrayUtils.toStr(incompleteChrs, ", "));
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

	public static void run() {
		String dir = "F:/testProjectSrc/UKBB_AffyAxiom/";
		String sourceDir = dir + "00src/";
		String projDir = dir + "project/";
		String propFileDir = "D:/projects/";
		String projName = "UKBioBank";
		String famFile = dir + "ukb1773_l2r_chrY_v2_s488374.fam";

		if (!Files.isWindows()) {
			dir = "/scratch.global/cole0482/UKBB/";
//			sourceDir = dir + "00src/";
			sourceDir = "/scratch.global/bb/all/";
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
