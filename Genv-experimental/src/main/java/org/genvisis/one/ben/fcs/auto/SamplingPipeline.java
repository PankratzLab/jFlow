package org.genvisis.one.ben.fcs.auto;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Queue;
import java.util.Random;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.one.ben.fcs.gating.Workbench.SampleNode;
import org.genvisis.one.ben.fcs.sub.EMInitializer;

public class SamplingPipeline {
	
//	private static final String CSV_DIR = "/panfs/roc/groups/15/thyagara/shared/HRS/UPLOAD ANALYZED DATA/";
//	private static final String WSP_DIR = "/panfs/roc/groups/15/thyagara/shared/HRS/UPLOAD WSP/";
//	private static final String FCS_DIR = "/panfs/roc/groups/15/thyagara/shared/HRS/UPLOAD FCS HRS samples/";
//	private static final String OUTLIERS_FILE = "/panfs/roc/groups/15/thyagara/shared/HRS/UPLOAD FCS HRS samples/outliers.txt"; 
//	private static final String OUT_DIR = "/scratch.global/cole0482/FCS/";
	
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
	
	private static final double SAMPLING_PCT = .25;
	private static final int OUTLIER_SAMPLE_EVERY = 5; // every fifth sampling will be from an outlier, if possible 
	
	private static final String P1_DATA_ROOT = "p1.";
	private static final String P2_DATA_ROOT = "p2.";
	
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
	
	final String csvDir;
	final String wspDir;
	final String fcsDir;
	final String outliersFile;
	final String outDir;
	
	public SamplingPipeline(String csvDir, String wspDir, String fcsDir, String outliersFile, String outDir) {
		this.csvDir = csvDir;
		this.wspDir = wspDir;
		this.fcsDir = fcsDir;
		this.outliersFile = outliersFile;
		this.outDir = outDir;
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
	}
	
	int checkPanel(String lwr) {
		if (lwr.contains("panel 1") || lwr.contains("p1")) {
			return 1;
		} else if (lwr.contains("panel 2") || lwr.contains("p2")) {
			return 2;
		}
		return 0;
	}
	
	private void loadFCS(File dir) {
		String[] fcsFiles = dir.list(FCS_FILTER);
		for (String s : fcsFiles) {
			int pnl = checkPanel(s.toLowerCase());
			if (pnl == 1) {
				fileToPathMap1.put(s, new File(ext.verifyDirFormat(dir.getAbsolutePath()) + s).getAbsolutePath());
			} else if (pnl == 2) {
				fileToPathMap2.put(s, new File(ext.verifyDirFormat(dir.getAbsolutePath()) + s).getAbsolutePath());
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
		String[] analysisFiles = (new File(csvDir)).list(CSV_FILTER);
		
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
			String[][] strData = HashVec.loadFileToStringMatrix(f, false, null, false);
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
			String[][] strData = HashVec.loadFileToStringMatrix(f, false, null, false);
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
		
		log.reportTime("Loaded analysis data for " + p1d.size() + " panel 1 samples and " + p2d.size() + " panel 2 samples");
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
	
	private boolean checkReady() {
		return analysisLoaded && wspLoaded && fcsPathsLoaded && outliersLoaded;
	}
	
	public void doFileSampling() {
		if (!checkReady()) {
			return;
		}

		ArrayList<String> p1 = new ArrayList<>(fileToPathMap1.keySet());
		ArrayList<String> p2 = new ArrayList<>(fileToPathMap2.keySet());
		p1.retainAll(p1d.keySet());
		p2.retainAll(p2d.keySet());
		
		ArrayList<String> p1O = new ArrayList<>(outliersP1);
		ArrayList<String> p2O = new ArrayList<>(outliersP2);
		p1O.retainAll(p1);
		p2O.retainAll(p2);
		
		int p1Cnt = p1.isEmpty() ? 0 : Math.max(1, (int) (p1.size() * SAMPLING_PCT));
		int p2Cnt = p2.isEmpty() ? 0 : Math.max(1, (int) (p2.size() * SAMPLING_PCT));
		
		log.reportTime("Sampling " + p1Cnt + " panel 1 files and " + p2Cnt + " panel 2 files.");
		
		Random rand = new Random();
		for (int i = 0; i < p1Cnt; i++) {
			if (i > 0 && i % OUTLIER_SAMPLE_EVERY == 0) {
				p1Sampling.add(p1O.get(rand.nextInt(p1O.size())));
			} else {
				p1Sampling.add(p1.get(rand.nextInt(p1.size())));
			}
		}
		for (int i = 0; i < p2Cnt; i++) {
			if (i > 0 && i % OUTLIER_SAMPLE_EVERY == 0) {
				p2Sampling.add(p2O.get(rand.nextInt(p2O.size())));
			} else {
				p2Sampling.add(p2.get(rand.nextInt(p2.size())));
			}
		}
		
		PrintWriter writer = Files.getAppropriateWriter(outDir + "p1.files.txt");
		for (String s : p1Sampling) {
			writer.println(s);
		}
		writer.flush();
		writer.close();
		writer = Files.getAppropriateWriter(outDir + "p2.files.txt");
		for (String s : p2Sampling) {
			writer.println(s);
		}
		writer.flush();
		writer.close();
		
		filesSampled = true;
	}
	
	public void doDataSampling() {
		if (!checkReady() && !filesSampled) {
			return;
		}
	  int proc = Runtime.getRuntime().availableProcessors() / 2;
	  final ThreadPoolExecutor threadPool1 = new ThreadPoolExecutor(proc, proc, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<Runnable>());
	  final ThreadPoolExecutor threadPool2 = new ThreadPoolExecutor(proc, proc, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<Runnable>());
	  final ConcurrentLinkedQueue<SampleNode> p1Queue = new ConcurrentLinkedQueue<>();
	  final ConcurrentLinkedQueue<SampleNode> p2Queue = new ConcurrentLinkedQueue<>();
//	  Class<? extends SampleProcessor> processorClass = LeafDataSampler.class;
	  Class<? extends SampleProcessor> processorClass = PercentageWriter.class;
	  
	  for (String s : p1Sampling) {
	  	SampleNode sn = wspLoader.panel1Nodes.get(s);
	  	if (sn != null) {
	  		sn.fcsFile = fileToPathMap1.get(s);
	  		p1Queue.add(sn);
	  	} else {
	  		log.reportError("Couldn't find WSP node for panel 1 fcs file: " + s);
	  	}
	  }
	  for (String s : p2Sampling) {
	  	SampleNode sn = wspLoader.panel2Nodes.get(s);
	  	if (sn != null) {
	  		sn.fcsFile = fileToPathMap2.get(s);
	  		p2Queue.add(sn);
	  	} else {
	  		log.reportError("Couldn't find WSP node for panel 2 fcs file: " + s);
	  	}
	  }
  	
	  PipelineRunnable p1Run = new PipelineRunnable(outDir + P1_DATA_ROOT, threadPool1, processorClass, p1Queue, log);
	  PipelineRunnable p2Run = new PipelineRunnable(outDir + P2_DATA_ROOT, threadPool2, processorClass, p2Queue, log);
	  Thread p1T = new Thread(p1Run);
	  Thread p2T = new Thread(p2Run);
	  p1T.start();
	  p2T.start();
	  
	  while(p1T.isAlive() || p2T.isAlive()) {
	    try {
	      Thread.sleep(1000);
	    } catch (InterruptedException e) { /**/ }
	  }
	  
	  p1Run.cleanup();
//	  p2Run.cleanup();
	}	
}

class PipelineRunnable implements Runnable {
	final ThreadPoolExecutor threadPool;
	final Class<? extends SampleProcessor> processorClass;
	final Queue<SampleNode> queue;
	final Logger log;
	final List<String> params;
	final int SAMPLING = 500;
	final Map<String, PrintWriter> gateWriters;
	final Map<String, AtomicInteger> writeCounts;
	final String outputFileRoot;
	final Map<String, Map<String, Double>> resultMap;
	
	public <T extends SampleProcessor> PipelineRunnable(String outputFileRoot, ThreadPoolExecutor threadPool, Class<T> processorClass, Queue<SampleNode> queue, Logger log) {
		this.threadPool = threadPool;
		this.processorClass = processorClass;
		this.queue = queue;
		this.log = log;
		this.params = new ArrayList<>();
    for (String s : EMInitializer.DATA_COLUMNS) {
      params.add(s);
    }
    this.outputFileRoot = outputFileRoot;
    gateWriters = new ConcurrentHashMap<>();
    writeCounts = new ConcurrentHashMap<>();
    resultMap = new ConcurrentHashMap<>();
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
        } catch (InterruptedException e) {}
        continue;
      }
      Runnable run = new Runnable() {
        @Override
        public void run() {
          try {
//          	processorClass.getConstructor(int.class, String.class, Map.class, Map.class, List.class).newInstance(SAMPLING, outputFileRoot, gateWriters, writeCounts, params).processSample(sn, log);
            processorClass.getConstructor(Map.class).newInstance(resultMap).processSample(sn, log);
          } catch (IOException e) {
            log.reportError(e.getMessage());
            // do not re-add to queue
          } catch (OutOfMemoryError e) {
            System.gc();
            queue.add(sn);
            // cleanup and re-add to queue
          } catch (InstantiationException | IllegalAccessException | IllegalArgumentException | InvocationTargetException | NoSuchMethodException | SecurityException e) {
          	log.reportException(e);
					} 
        }
      };
      log.reportTime("Submitting processing job for " + sn.fcsFile);
      threadPool.execute(run);
	    try {
	      Thread.sleep(500);
	    } catch (InterruptedException e) {}
    }
	  try {
	    threadPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
	  } catch (InterruptedException e1) {}
  }
  
  private static final String[] order = {
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A)",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H)",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A)",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / B cells (CD3- CD19+) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / B cells (CD3- CD19+) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / IgD+ memory Bcells (CD27+) (Comp-BUV 737-A (IgD) v Comp-BB515-A (CD27))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / B cells (CD3- CD19+) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / IgD- memory Bcells (CD27+) (Comp-BUV 737-A (IgD) v Comp-BB515-A (CD27))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / B cells (CD3- CD19+) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / naive Bcells (CD27- IgD+) (Comp-BUV 737-A (IgD) v Comp-BB515-A (CD27))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / activated cytotoxic Tcells (CD8+ HLA-DR+) (Comp-PE-CF594-A (HLA-DR) v Comp-BUV 395-A (CD8))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / central memory cytotoxic Tcells (CCR7+ , CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / central memory cytotoxic Tcells (CCR7+ , CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / central memory cytotoxic Tcells (CD95+ CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / central memory cytotoxic Tcells (CCR7+ , CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / effector memory cytotoxic Tcells (CD95+ CD28-) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / central memory cytotoxic Tcells (CCR7+ , CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / naive cytotoxic Tcells (CD95- CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / central memory cytotoxic Tcells (CD95+ CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector cytotoxic Tcells  (CCR7-  CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector cytotoxic Tcells  (CCR7-  CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / central memory cytotoxic Tcells (CD95+ CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector cytotoxic Tcells  (CCR7-  CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / cytotoxic Tcells CD27- , CD28+ (Comp-BB515-A (CD27) v Comp-BV 510-A (CD28))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector cytotoxic Tcells  (CCR7-  CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / effector memory cytotoxic Tcells (CD95+ CD28-) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector cytotoxic Tcells  (CCR7-  CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / naive cytotoxic Tcells (CD95- CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector cytotoxic Tcells  (CCR7-  CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / pE cytotoxic Tcells (CD27-  CD28-) (Comp-BB515-A (CD27) v Comp-BV 510-A (CD28))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector cytotoxic Tcells  (CCR7-  CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / pE1 cytotoxic Tcells (CD27+  CD28+) (Comp-BB515-A (CD27) v Comp-BV 510-A (CD28))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector cytotoxic Tcells  (CCR7-  CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / pE2 cytotoxic Tcells (CD27+ , CD28-) (Comp-BB515-A (CD27) v Comp-BV 510-A (CD28))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector memory cytotoxic Tcells (CCR7- , CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector memory cytotoxic Tcells (CCR7- , CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / central memory cytotoxic Tcells (CD95+ CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector memory cytotoxic Tcells (CCR7- , CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / effector memory cytotoxic Tcells (CD95+ CD28-) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector memory cytotoxic Tcells (CCR7- , CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / EM1 cytotoxic Tcells (CD27+  CD28+) (Comp-BB515-A (CD27) v Comp-BV 510-A (CD28))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector memory cytotoxic Tcells (CCR7- , CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / EM2 cytotoxic Tcells (CD27+  CD28-) (Comp-BB515-A (CD27) v Comp-BV 510-A (CD28))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector memory cytotoxic Tcells (CCR7- , CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / EM3 cytotoxic Tcells (CD27-  CD28-) (Comp-BB515-A (CD27) v Comp-BV 510-A (CD28))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector memory cytotoxic Tcells (CCR7- , CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / EM4 cytotoxic Tcells (CD27-  CD28+) (Comp-BB515-A (CD27) v Comp-BV 510-A (CD28))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector memory cytotoxic Tcells (CCR7- , CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / naive cytotoxic Tcells (CD95- CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector memory cytotoxic Tcells (CD95+ CD28-) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / naive cytotoxic Tcells (CCR7+ , CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / naive cytotoxic Tcells (CCR7+ , CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / central memory cytotoxic Tcells (CD95+ CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / naive cytotoxic Tcells (CCR7+ , CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / effector memory cytotoxic Tcells (CD95+ CD28-) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / naive cytotoxic Tcells (CCR7+ , CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / naive cytotoxic Tcells (CD95- CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / cytotoxic Tcells-CD8+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / naive cytotoxic Tcells (CD95- CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / activated helper Tcells (CD4+ HLA-DR+) (Comp-PE-CF594-A (HLA-DR) v Comp-APC-Cy7-A (CD4))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / central memory helper Tcells (CCR7+ CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / central memory helper Tcells (CCR7+ CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / central memory helper Tcells (CD95+, CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / central memory helper Tcells (CCR7+ CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / effector memory helper Tcells (CD95+, CD28-) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / central memory helper Tcells (CCR7+ CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / naive helper Tcells (CD95-, CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / central memory helper Tcells (CD95+, CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector helper Tcells (CCR7- CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector helper Tcells (CCR7- CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / central memory helper Tcells (CD95+, CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector helper Tcells (CCR7- CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / effector memory helper Tcells (CD95+, CD28-) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector helper Tcells (CCR7- CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / naive helper Tcells (CD95-, CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector memory helper Tcells (CCR7- CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector memory helper Tcells (CCR7- CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / central memory helper Tcells (CD95+, CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector memory helper Tcells (CCR7- CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / effector memory helper Tcells (CD95+, CD28-) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector memory helper Tcells (CCR7- CD45RA-) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / naive helper Tcells (CD95-, CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / effector memory helper Tcells (CD95+, CD28-) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / naive helper Tcells (CCR7+ CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / naive helper Tcells (CCR7+ CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / central memory helper Tcells (CD95+, CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / naive helper Tcells (CCR7+ CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / effector memory helper Tcells (CD95+, CD28-) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / naive helper Tcells (CCR7+ CD45RA+) (Comp-BV 421-A (CCR7) v Comp-BV 711-A (CD45RA)) / naive helper Tcells (CD95-, CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
"Lymphocytes (SSC-A v FSC-A) (FSC-A v SSC-A) / Single Cells (FSC-H v FSC-W) (FSC-W v FSC-H) / Live cells (PE-) (Comp-PE-A (L/D) v SSC-A) / Tcells (CD3+ CD19-) (Comp-APC-A (CD3) v Comp-PE-Cy7-A (CD19)) / Helper Tcells-CD4+ (Comp-APC-Cy7-A (CD4) v Comp-BUV 395-A (CD8)) / naive helper Tcells (CD95-, CD28+) (Comp-BV 605-A (CD95) v Comp-BV 510-A (CD28))",
  };
  
  public void cleanup() {
//  	for (Entry<String, PrintWriter> entry : gateWriters.entrySet()) {
//  		entry.getValue().flush();
//  		entry.getValue().close();
//  		log.reportTime("Wrote " + writeCounts.get(entry.getKey()).get() + " for gate " + entry.getKey());
//  	}
  	PrintWriter writer = Files.getAppropriateWriter("F:/Flow/gatingTest/p1.pcts.xln");
  	boolean header = false;
  	for (Entry<String, Map<String, Double>> entry : resultMap.entrySet()) {
  		if (!header) {
	  		StringBuilder sb = new StringBuilder();
	  		for (String s : order) {
	  			sb.append("\t").append(s);
	  		}
	  		writer.println(sb.toString());
	  		header = true;
  		}
  		StringBuilder sb = new StringBuilder(entry.getKey());
  		for (String s : order) {
  			sb.append("\t").append(entry.getValue().get(s));
  		}
  		writer.println(sb.toString());
  	}
  	writer.flush();
  	writer.close();
  }
  
}
