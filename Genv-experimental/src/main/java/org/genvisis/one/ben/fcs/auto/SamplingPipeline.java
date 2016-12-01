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
import java.util.Queue;
import java.util.Random;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.one.ben.fcs.gating.Workbench.SampleNode;
import org.genvisis.one.ben.fcs.sub.EMInitializer;

public class SamplingPipeline {
	
	private static final String CSV_DIR = "/panfs/roc/groups/15/thyagara/shared/HRS/UPLOAD ANALYZED DATA/";
	private static final String CSV_EXT = ".csv";
	private static final FilenameFilter CSV_FILTER = new FilenameFilter() {
		@Override
		public boolean accept(File dir, String name) {
			return name.endsWith(CSV_EXT);
		}
	};
	
	private static final String WSP_DIR = "/panfs/roc/groups/15/thyagara/shared/HRS/UPLOAD WSP/";
	
	private static final String FCS_DIR = "/panfs/roc/groups/15/thyagara/shared/HRS/UPLOAD FCS HRS samples/";
	private static final String FCS_EXT = ".fcs";
	private static final FilenameFilter FCS_FILTER = new FilenameFilter() {
		@Override
		public boolean accept(File dir, String name) {
			return name.endsWith(FCS_EXT);
		}
	};
	
	private static final String OUTLIERS_FILE = "/panfs/roc/groups/15/thyagara/shared/HRS/UPLOAD FCS HRS samples/outliers.txt"; 
	private static final double SAMPLING_PCT = .5;
	private static final int OUTLIER_SAMPLE_EVERY = 5; // every fifth sampling will be from an outlier, if possible 
	
	private static final String OUT_DIR = "/scratch.global/cole0482/FCS/";
	private static final String P1_DATA_ROOT = OUT_DIR + "p1.";
	private static final String P2_DATA_ROOT = OUT_DIR + "p2.";
	
	WSPLoader wspLoader;

	Logger log = new Logger();
	
	HashMap<String, AnalysisData> p1d;
	HashMap<String, AnalysisData> p2d;
	ArrayList<String> p1fcs;
	ArrayList<String> p2fcs;
	HashMap<String, String> fileToPathMap;
	HashSet<String> outliersP1;
	HashSet<String> outliersP2;
	ArrayList<String> p1Sampling;
	ArrayList<String> p2Sampling;

	boolean analysisLoaded = false;
	boolean wspLoaded = false;
	boolean fcsPathsLoaded = false;
	boolean outliersLoaded = false;
	boolean filesSampled = false;
	
	private SamplingPipeline() {
		p1d = new HashMap<>();
		p2d = new HashMap<>();
		p1fcs = new ArrayList<>();
		p2fcs = new ArrayList<>();
		fileToPathMap = new HashMap<>();
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
	
	public void loadFCSFilePaths() {
		File[] sub = new File(FCS_DIR).listFiles();
		for (File f : sub) {
			if (f.isDirectory()) {
				String[] fcsFiles = f.list(FCS_FILTER);
				for (String s : fcsFiles) {
					fileToPathMap.put(s, new File(s).getAbsolutePath());
				}
			}
		}
		fcsPathsLoaded = true;
	}
	
	public void loadAnalysisData() {
		String[] analysisFiles = (new File(CSV_DIR)).list(CSV_FILTER);
		
		ArrayList<String> p1 = new ArrayList<>();
		ArrayList<String> p2 = new ArrayList<>();
		
		for (String s : analysisFiles) {
			String lwr = s.toLowerCase();
			if (lwr.contains("panel 1") || lwr.contains("p1")) {
				p1.add(s);
			} else if (lwr.contains("panel 2") || lwr.contains("p2")) {
				p2.add(s);
			}
		}
		
		for (String f : p1) {
			if (!f.startsWith(CSV_DIR)) {
				f = CSV_DIR + f;
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
					p1d.put(ad.fcsFilename, ad);
				}
			}
		}
		for (String f : p2) {
			if (!f.startsWith(CSV_DIR)) {
				f = CSV_DIR + f;
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
		
		p1fcs.addAll(p1d.keySet());
		p2fcs.addAll(p2d.keySet());
		
		log.reportTime("Loaded analysis data for " + p1fcs.size() + " panel 1 samples and " + p2fcs.size() + " panel 2 samples");
		analysisLoaded = true;
	}
	
	public void loadWSPInfo() {
		wspLoader = new WSPLoader();
		boolean success = wspLoader.loadWorkspaces(WSP_DIR);
		if (success) {
			wspLoaded = true;
		} else {
			log.reportError("Loading WSP data failed - please check log and try again.");
		}
	}
	
	public void loadOutliers() {
		String[] out = HashVec.loadFileToStringArray(OUTLIERS_FILE, false, null, false);
		for (String s : out) {
			if (s.toLowerCase().contains("panel 1") || s.toLowerCase().contains("p1")) {
				outliersP1.add(s);
			} else if (s.toLowerCase().contains("panel 2") || s.toLowerCase().contains("p2")) {
				outliersP2.add(s);
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

		ArrayList<String> p1 = new ArrayList<>(p1fcs);
		ArrayList<String> p2 = new ArrayList<>(p2fcs);
		p1.retainAll(p1d.keySet());
		p2.retainAll(p2d.keySet());
		
		ArrayList<String> p1O = new ArrayList<>(outliersP1);
		ArrayList<String> p2O = new ArrayList<>(outliersP2);
		p1O.retainAll(p1);
		p2O.retainAll(p2);
		
		int p1Cnt = (int) (p1.size() * SAMPLING_PCT);
		int p2Cnt = (int) (p1.size() * SAMPLING_PCT);
		
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
		
		PrintWriter writer = Files.getAppropriateWriter(OUT_DIR + "p1.files.txt");
		for (String s : p1Sampling) {
			writer.println(s);
		}
		writer.flush();
		writer.close();
		writer = Files.getAppropriateWriter(OUT_DIR + "p2.files.txt");
		for (String s : p1Sampling) {
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
	  Class<? extends SampleProcessor> processorClass = LeafDataSampler.class;
	  
	  for (String s : p1Sampling) {
	  	SampleNode sn = wspLoader.panel1Nodes.get(s);
	  	if (sn != null) {
	  		p1Queue.add(sn);
	  	} else {
	  		log.reportError("Couldn't find WSP node for panel 1 fcs file: " + s);
	  	}
	  }
	  for (String s : p2Sampling) {
	  	SampleNode sn = wspLoader.panel2Nodes.get(s);
	  	if (sn != null) {
	  		p2Queue.add(sn);
	  	} else {
	  		log.reportError("Couldn't find WSP node for panel 2 fcs file: " + s);
	  	}
	  }
	  
	  PipelineRunnable p1Run = new PipelineRunnable(P1_DATA_ROOT, threadPool1, processorClass, p1Queue, log);
	  PipelineRunnable p2Run = new PipelineRunnable(P2_DATA_ROOT, threadPool2, processorClass, p2Queue, log);
	  Thread p1T = new Thread(p1Run);
	  Thread p2T = new Thread(p2Run);
	  p1T.start();
	  p2T.start();
	  
	  while(p1T.isAlive() || p2T.isAlive()) {
	    try {
	      Thread.sleep(1000);
	    } catch (InterruptedException e) {}
	  }
	  
	  p1Run.cleanup();
	  p2Run.cleanup();
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
	final String outputFileRoot;
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
              processorClass.getConstructor(WSPLoader.class).newInstance(SAMPLING, outputFileRoot, gateWriters, params).processSample(sn, log);
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
  
  public void cleanup() {
  	for (PrintWriter writer : gateWriters.values()) {
  		writer.flush();
  		writer.close();
  	}
  }
  
}
