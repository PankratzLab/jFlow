package org.genvisis.one.ben.fcs.auto;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Queue;
import java.util.Random;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.one.ben.fcs.auto.proc.ProcessorFactory;
import org.genvisis.one.ben.fcs.auto.proc.SampleProcessor;
import org.genvisis.one.ben.fcs.gating.Workbench.SampleNode;

public class SamplingPipeline {

	// private static final String CSV_DIR = "/panfs/roc/groups/15/thyagara/shared/HRS/UPLOAD ANALYZED
	// DATA/";
	// private static final String WSP_DIR = "/panfs/roc/groups/15/thyagara/shared/HRS/UPLOAD WSP/";
	// private static final String FCS_DIR = "/panfs/roc/groups/15/thyagara/shared/HRS/UPLOAD FCS HRS
	// samples/";
	// private static final String OUTLIERS_FILE = "/panfs/roc/groups/15/thyagara/shared/HRS/UPLOAD
	// FCS HRS samples/outliers.txt";
	// private static final String OUT_DIR = "/scratch.global/cole0482/FCS/";

	private static final String CSV_EXT = ".csv";
	private static final FilenameFilter CSV_FILTER = new FilenameFilter() {
		@Override
		public boolean accept(File dir, String name) {
			return name.endsWith(CSV_EXT);
		}
	};


	private static final String FCS_EXT = ".fcs";
	private static final FilenameFilter FCS_FILTER = new FilenameFilter() {
		@Override
		public boolean accept(File dir, String name) {
			return name.endsWith(FCS_EXT);
		}
	};

	private static final int OUTLIER_SAMPLE_EVERY = 0; // no outlier sampling // every fifth sampling
																										 // will be from an outlier, if possible

	WSPLoader wspLoader;

	Logger log = new Logger();

	HashMap<String, AnalysisData> p1d;
	HashMap<String, AnalysisData> p2d;
	HashMap<String, String> fileToPathMap1;
	HashMap<String, String> fileToPathMap2;
	HashSet<String> outliersP1;
	HashSet<String> outliersP2;
	ArrayList<String> p1Sampling;
	ArrayList<String> p2Sampling;

	boolean analysisLoaded = false;
	boolean wspLoaded = false;
	boolean fcsPathsLoaded = false;
	boolean outliersLoaded = false;
	boolean filesSampled = false;

	boolean isFinished = false;

	private double samplingPct = 1;

	final String csvDir;
	final String wspDir;
	final String fcsDir;
	final String outliersFile;
	final String outDir;
	final HashSet<String> highPriority;
	final HashSet<String> lowPriority;
	final int panelToRun;
	final ProcessorFactory<? extends SampleProcessor> processorFactory;

	public SamplingPipeline(double sampPct, String csvDir, String wspDir, String fcsDir,
													String outliersFile, String outDir, int panel, String[] priorityFiles,
													ProcessorFactory<? extends SampleProcessor> processorFactory) {
		this.samplingPct = sampPct;
		this.csvDir = csvDir;
		this.wspDir = wspDir;
		this.fcsDir = fcsDir;
		this.outliersFile = outliersFile;
		this.outDir = outDir;
		this.processorFactory = processorFactory;
		this.panelToRun = panel;
		this.highPriority = priorityFiles != null && priorityFiles[0] != null
																																					? HashVec.loadFileToHashSet(priorityFiles[0],
																																																			false)
																																					: null;
		this.lowPriority = priorityFiles != null && priorityFiles.length > 1
											 && priorityFiles[1] != null
																									 ? HashVec.loadFileToHashSet(priorityFiles[1],
																																							 false)
																									 : null;
		p1d = new HashMap<>();
		p2d = new HashMap<>();
		fileToPathMap1 = new HashMap<>();
		fileToPathMap2 = new HashMap<>();
		outliersP1 = new HashSet<>();
		outliersP2 = new HashSet<>();
		p1Sampling = new ArrayList<>();
		p2Sampling = new ArrayList<>();
	}

	class AnalysisData {
		String fcsFilename;
		HashMap<String, Double> data = new HashMap<>();
	}

	public void run() {
		loadFCSFilePaths();
		loadAnalysisData();
		loadOutliers();
		loadWSPInfo();
		doFileSampling();
		doDataSampling();
		setFinished();
	}

	int checkPanel(String lwr) {
		int panel = 0;
		if (lwr.contains("panel 1") || lwr.contains("p1") || lwr.contains("panel one")) {
			panel = 1;
		} else if (lwr.contains("panel 2") || lwr.contains("p2") || lwr.contains("panel two")) {
			panel = 2;
		}
		if (panel > 0 && (panelToRun == -1 || panelToRun == panel)) {
			return panel;
		}
		return 0;
	}

	private void loadFCS(File dir) {
		if (!dir.canRead()) {
			log.reportError("Can't read directory: " + dir.getAbsolutePath());
			return;
		}
		String[] fcsFiles = dir.list(FCS_FILTER);
		if (fcsFiles == null) {
			log.reportError("No FCS files found in directory: " + dir.getAbsolutePath());
		} else {
			for (String s : fcsFiles) {
				int pnl = checkPanel(s.toLowerCase());
				if (pnl == 1) {
					fileToPathMap1.put(s, new File(ext.verifyDirFormat(dir.getAbsolutePath())
																				 + s).getAbsolutePath());
				} else if (pnl == 2) {
					fileToPathMap2.put(s, new File(ext.verifyDirFormat(dir.getAbsolutePath())
																				 + s).getAbsolutePath());
				}
			}
		}
		File[] sub = dir.listFiles();
		for (File f : sub) {
			if (f.isDirectory()) {
				loadFCS(f);
			}
		}
	}


	public void loadFCSFilePaths() {
		File dir = new File(fcsDir);
		loadFCS(dir);
		fcsPathsLoaded = true;
	}

	public void loadAnalysisData() {
		String[] analysisFiles = csvDir == null ? new String[0] : (new File(csvDir)).list(CSV_FILTER);

		ArrayList<String> p1 = new ArrayList<>();
		ArrayList<String> p2 = new ArrayList<>();

		for (String s : analysisFiles) {
			String lwr = s.toLowerCase();
			int pnl = checkPanel(lwr);
			if (pnl == 1) {
				p1.add(s);
			} else if (pnl == 2) {
				p2.add(s);
			}
		}

		for (String f : p1) {
			if (!f.startsWith(csvDir)) {
				f = csvDir + f;
			}
			String[][] strData = HashVec.loadFileToStringMatrix(f, false, null);
			for (int i = 1; i < strData.length; i++) {
				String[] line = strData[i];
				if (!"sd".equalsIgnoreCase(line[0]) && !"mean".equalsIgnoreCase(line[0])) {
					AnalysisData ad = new AnalysisData();
					ad.fcsFilename = line[0];
					for (int j = 1; j < line.length; j++) {
						if (!"".equals(strData[0][j])) {
							ad.data.put(strData[0][j], "".equals(line[j]) ? Double.NaN : Double.valueOf(line[j]));
						}
					}
					p1d.put(ad.fcsFilename, ad);
				}
			}
		}
		for (String f : p2) {
			if (!f.startsWith(csvDir)) {
				f = csvDir + f;
			}
			String[][] strData = HashVec.loadFileToStringMatrix(f, false, null);
			for (int i = 1; i < strData.length; i++) {
				String[] line = strData[i];
				if (!"sd".equalsIgnoreCase(line[0]) && !"mean".equalsIgnoreCase(line[0])) {
					AnalysisData ad = new AnalysisData();
					ad.fcsFilename = line[0];
					for (int j = 1; j < line.length; j++) {
						ad.data.put(strData[0][j], Double.valueOf(line[j]));
					}
					p2d.put(ad.fcsFilename, ad);
				}
			}
		}

		log.reportTime("Loaded analysis data for " + p1d.size() + " panel 1 samples and " + p2d.size()
									 + " panel 2 samples");
		analysisLoaded = true;
	}

	public void loadWSPInfo() {
		wspLoader = new WSPLoader();
		boolean success = wspLoader.loadWorkspaces(wspDir);
		if (success) {
			wspLoaded = true;
		} else {
			log.reportError("Loading WSP data failed - please check log and try again.");
		}
	}

	public void loadOutliers() {
		if (outliersFile != null) {
			String[] out = HashVec.loadFileToStringArray(outliersFile, false, null, false);
			for (String s : out) {
				int pnl = checkPanel(s.toLowerCase());
				if (pnl == 1) {
					outliersP1.add(s);
				} else if (pnl == 2) {
					outliersP2.add(s);
				}
			}
		}
		outliersLoaded = true;
	}

	protected boolean checkReady() {
		return analysisLoaded && wspLoaded && fcsPathsLoaded && outliersLoaded;
	}

	protected boolean checkFinished() {
		return isFinished;
	}

	private void setFinished() {
		isFinished = true;
	}

	public void doFileSampling() {
		if (!checkReady()) {
			return;
		}

		ArrayList<String> p1 = new ArrayList<>(fileToPathMap1.keySet());
		ArrayList<String> p2 = new ArrayList<>(fileToPathMap2.keySet());
		if (!p1d.isEmpty()) {
			p1.retainAll(p1d.keySet());
		}
		if (!p2d.isEmpty()) {
			p2.retainAll(p2d.keySet());
		}

		ArrayList<String> p1O = new ArrayList<>(outliersP1);
		ArrayList<String> p2O = new ArrayList<>(outliersP2);
		p1O.retainAll(p1);
		p2O.retainAll(p2);

		p1.removeAll(p1O);
		p2.removeAll(p2O);

		int p1Cnt = p1.isEmpty() ? 0 : Math.max(1, (int) (p1.size() * samplingPct));
		int p2Cnt = p2.isEmpty() ? 0 : Math.max(1, (int) (p2.size() * samplingPct));

		log.reportTime("Sampling " + p1Cnt + " panel 1 files and " + p2Cnt + " panel 2 files.");
		if (samplingPct > .95) {
			p1Sampling.addAll(p1);
			p2Sampling.addAll(p2);
		} else {
			Random rand = new Random();
			for (int i = 0; i < p1Cnt; i++) {
				if (!p1O.isEmpty() && OUTLIER_SAMPLE_EVERY > 0
						&& (i > 0 && i % OUTLIER_SAMPLE_EVERY == 0)) {
					p1Sampling.add(p1O.get(rand.nextInt(p1O.size())));
				} else {
					p1Sampling.add(p1.get(rand.nextInt(p1.size())));
				}
			}
			for (int i = 0; i < p2Cnt; i++) {
				if (!p2O.isEmpty() && OUTLIER_SAMPLE_EVERY > 0
						&& (i > 0 && i % OUTLIER_SAMPLE_EVERY == 0)) {
					p2Sampling.add(p2O.get(rand.nextInt(p2O.size())));
				} else {
					p2Sampling.add(p2.get(rand.nextInt(p2.size())));
				}
			}
		}

		new File(outDir).mkdirs();
		PrintWriter writer;
		if (p1Sampling.size() > 0) {
			writer = Files.getAppropriateWriter(outDir + "p1.files.txt");
			for (String s : p1Sampling) {
				writer.println(s);
			}
			writer.flush();
			writer.close();
		}
		if (p2Sampling.size() > 0) {
			writer = Files.getAppropriateWriter(outDir + "p2.files.txt");
			for (String s : p2Sampling) {
				writer.println(s);
			}
			writer.flush();
			writer.close();
		}

		filesSampled = true;
	}

	public void doDataSampling() {
		if (!checkReady() && !filesSampled) {
			return;
		}
		int proc = Math.max(1, Runtime.getRuntime().availableProcessors() / 2);
		final ThreadPoolExecutor threadPool1 = new ThreadPoolExecutor(
																																	proc,
																																	proc,
																																	0L,
																																	TimeUnit.MILLISECONDS,
																																	new LinkedBlockingQueue<Runnable>());
		final ThreadPoolExecutor threadPool2 = new ThreadPoolExecutor(
																																	proc,
																																	proc,
																																	0L,
																																	TimeUnit.MILLISECONDS,
																																	new LinkedBlockingQueue<Runnable>());
		final ConcurrentLinkedQueue<SampleNode> p1Queue = new ConcurrentLinkedQueue<>();
		final ConcurrentLinkedQueue<SampleNode> p2Queue = new ConcurrentLinkedQueue<>();

		PrintWriter sampWspMatch1 = Files.getAppropriateWriter(outDir + "/matchFile_p1.xln");
		PrintWriter sampWspMatch2 = Files.getAppropriateWriter(outDir + "/matchFile_p2.xln");

		if (highPriority != null) {
			for (String s1 : p1Sampling) {
				if (highPriority.contains(ext.replaceWithLinuxSafeCharacters(ext.removeDirectoryInfo(s1)))) {
					SampleNode sn = wspLoader.getPanel1Nodes().get(s1);
					if (sn != null) {
						sn.fcsFile = fileToPathMap1.get(s1);
						sampWspMatch1.println(s1 + "\t" + sn.wspFile + "\t" + sn.fcsFile);
						p1Queue.add(sn);
					} else {
						log.reportError("Couldn't find WSP node for panel 1 fcs file: " + s1);
					}
					break;
				}
			}
			for (String s1 : p2Sampling) {
				if (highPriority.contains(ext.replaceWithLinuxSafeCharacters(ext.removeDirectoryInfo(s1)))) {
					SampleNode sn = wspLoader.panel2Nodes.get(s1);
					if (sn != null) {
						sn.fcsFile = fileToPathMap2.get(s1);
						sampWspMatch2.println(s1 + "\t" + sn.wspFile + "\t" + sn.fcsFile);
						p2Queue.add(sn);
					} else {
						log.reportError("Couldn't find WSP node for panel 2 fcs file: " + s1);
					}
					break;
				}
			}
		}

		for (String s : p1Sampling) {
			if (highPriority != null
					&& highPriority.contains(ext.replaceWithLinuxSafeCharacters(ext.removeDirectoryInfo(s)))) {
				continue;
			}
			if (lowPriority != null
					&& lowPriority.contains(ext.replaceWithLinuxSafeCharacters(ext.removeDirectoryInfo(s)))) {
				continue;
			}
			SampleNode sn = wspLoader.getPanel1Nodes().get(s);
			if (sn != null) {
				sn.fcsFile = fileToPathMap1.get(s);
				sampWspMatch1.println(s + "\t" + sn.wspFile + "\t" + sn.fcsFile);
				p1Queue.add(sn);
			} else {
				log.reportError("Couldn't find WSP node for panel 1 fcs file: " + s);
			}
		}
		for (String s : p2Sampling) {
			if (highPriority != null
					&& highPriority.contains(ext.replaceWithLinuxSafeCharacters(ext.removeDirectoryInfo(s)))) {
				continue;
			}
			if (lowPriority != null
					&& lowPriority.contains(ext.replaceWithLinuxSafeCharacters(ext.removeDirectoryInfo(s)))) {
				continue;
			}
			SampleNode sn = wspLoader.panel2Nodes.get(s);
			if (sn != null) {
				sn.fcsFile = fileToPathMap2.get(s);
				sampWspMatch2.println(s + "\t" + sn.wspFile + "\t" + sn.fcsFile);
				p2Queue.add(sn);
			} else {
				log.reportError("Couldn't find WSP node for panel 2 fcs file: " + s);
			}
		}

		if (lowPriority != null) {
			for (String s : p1Sampling) {
				if (lowPriority.contains(ext.replaceWithLinuxSafeCharacters(ext.removeDirectoryInfo(s)))) {
					SampleNode sn = wspLoader.getPanel1Nodes().get(s);
					if (sn != null) {
						sn.fcsFile = fileToPathMap1.get(s);
						sampWspMatch1.println(s + "\t" + sn.wspFile + "\t" + sn.fcsFile);
						p1Queue.add(sn);
					} else {
						log.reportError("Couldn't find WSP node for panel 1 fcs file: " + s);
					}
				}
			}
			for (String s : p2Sampling) {
				if (lowPriority.contains(ext.replaceWithLinuxSafeCharacters(ext.removeDirectoryInfo(s)))) {
					SampleNode sn = wspLoader.panel2Nodes.get(s);
					if (sn != null) {
						sn.fcsFile = fileToPathMap2.get(s);
						sampWspMatch2.println(s + "\t" + sn.wspFile + "\t" + sn.fcsFile);
						p2Queue.add(sn);
					} else {
						log.reportError("Couldn't find WSP node for panel 2 fcs file: " + s);
					}
				}
			}
		}

		sampWspMatch1.close();
		sampWspMatch2.close();

		log.reportTime("Finished matching sample files to workspace files.");

		AbstractPipelineRunnable p1Run = new AbstractPipelineRunnable(threadPool1, processorFactory,
																																	p1Queue, log);
		AbstractPipelineRunnable p2Run = new AbstractPipelineRunnable(threadPool2, processorFactory,
																																	p2Queue, log);

		ExecutorService serve = Executors.newFixedThreadPool(2);
		serve.submit(p1Run);
		serve.submit(p2Run);

		serve.shutdown();
		try {
			serve.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
		} catch (InterruptedException e) {
		}

		processorFactory.cleanup(p1Run);
		processorFactory.cleanup(p2Run);
	}
}


class AbstractPipelineRunnable implements Runnable {
	final ThreadPoolExecutor threadPool;
	final ProcessorFactory<? extends SampleProcessor> processorFactory;
	final Queue<SampleNode> queue;
	final Logger log;
	final int myIndex;
	static int sharedIndex = 1;

	public <T extends SampleProcessor> AbstractPipelineRunnable(ThreadPoolExecutor threadPool,
																															final ProcessorFactory<? extends SampleProcessor> factory,
																															Queue<SampleNode> queue, Logger log) {
		this.threadPool = threadPool;
		this.processorFactory = factory;
		this.queue = queue;
		this.log = log;
		myIndex = sharedIndex;
		sharedIndex++;
	}

	@Override
	public void run() {
		while (!threadPool.isShutdown()) {
			final SampleNode sn = queue.poll();
			if (sn == null) {
				if (threadPool.getQueue().isEmpty() && threadPool.getActiveCount() == 0) {
					threadPool.shutdown();
				}
				try {
					Thread.sleep(100);
				} catch (InterruptedException e) {
				}
				continue;
			}
			log.reportTime("Submitting processing job for " + sn.fcsFile);
			threadPool.execute(() -> {
				try {
					processorFactory.createProcessor(this, myIndex).processSample(sn, log);
				} catch (IOException e) {
					log.reportException(e);
				} catch (OutOfMemoryError e) {
					// cleanup and re-add to queue
					System.gc();
					queue.add(sn);
				}
			});
			try {
				Thread.sleep(500);
			} catch (InterruptedException e) {
			}
		}
		try {
			threadPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
		} catch (InterruptedException e1) {
		}
	}

}
