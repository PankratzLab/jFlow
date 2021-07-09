package org.genvisis.fcs.auto;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

import org.genvisis.fcs.auto.proc.ProcessorFactory;
import org.genvisis.fcs.auto.proc.SampleProcessor;
import org.genvisis.fcs.gating.Workbench.SampleNode;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;

import com.google.common.collect.HashMultiset;

public class SamplingPipeline {

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

  class PanelData {
    final Panel panel;
    final Map<String, AnalysisData> data;
    final Map<String, String> filePathMap;
    final Set<String> outliers;
    final List<String> sampling;

    public PanelData(Panel p) {
      this.panel = p;
      data = new HashMap<>();
      filePathMap = new HashMap<>();
      outliers = new HashSet<>();
      sampling = new ArrayList<>();
    }
  }

  Map<Panel, PanelData> panelData;

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
  final List<Panel> panelsToRun;
  final ProcessorFactory<? extends SampleProcessor> processorFactory;

  public SamplingPipeline(double sampPct, String csvDir, String wspDir, String fcsDir,
                          String outliersFile, String outDir, List<Panel> panels,
                          String[] priorityFiles,
                          ProcessorFactory<? extends SampleProcessor> processorFactory) {
    this.samplingPct = sampPct;
    this.csvDir = csvDir;
    this.wspDir = wspDir;
    this.fcsDir = fcsDir;
    this.outliersFile = outliersFile;
    this.outDir = outDir;
    this.processorFactory = processorFactory;
    this.panelsToRun = panels;
    this.highPriority = priorityFiles != null
                        && priorityFiles[0] != null ? HashVec.loadFileToHashSet(priorityFiles[0], false) : null;
    this.lowPriority = priorityFiles != null && priorityFiles.length > 1
                       && priorityFiles[1] != null ? HashVec.loadFileToHashSet(priorityFiles[1], false) : null;
    this.panelData = new HashMap<>();
    for (Panel p : panels) {
      this.panelData.put(p, new PanelData(p));
    }
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

  Panel checkPanel(String lwr) {
    for (Panel p : panelsToRun) {
      for (String a : p.aliases) {
        if (lwr.contains(a)) {
          return p;
        }
      }
    }
    return null;
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
        if (!new File(ext.verifyDirFormat(dir.getAbsolutePath()) + s).canRead()) {
          log.reportError("Can't read file: " + ext.verifyDirFormat(dir.getAbsolutePath()) + s);
          continue;
        }
        Panel pnl = checkPanel(s.toLowerCase());
        if (panelData.containsKey(pnl)) {
          panelData.get(pnl).filePathMap.put(s, new File(ext.verifyDirFormat(dir.getAbsolutePath())
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

    Map<String, Panel> panelMap = new HashMap<>();

    for (String s : analysisFiles) {
      String lwr = s.toLowerCase();
      Panel pnl = checkPanel(lwr);
      if (pnl != null) {
        panelMap.put(s, pnl);
      } else {
        log.reportWarning("Couldn't identify panel for analysis file " + s);
      }
    }

    HashMultiset<Panel> panelCounts = HashMultiset.create();
    for (String f : panelMap.keySet()) {
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
          panelData.get(panelMap.get(f)).data.put(ad.fcsFilename, ad);
          panelCounts.add(panelMap.get(f));
        }
      }
    }

    if (panelCounts.size() > 0) {
      StringBuilder logMsg = new StringBuilder("Loaded analysis data for ");
      int count = 1;
      for (Panel p : panelCounts.elementSet()) {
        logMsg.append(panelCounts.count(p)).append(" ").append(p.getName()).append(" samples");
        if (count < panelCounts.elementSet().size() - 1) {
          logMsg.append(" and ");
        }
        count++;
      }
      log.reportTime(logMsg.toString());
    }
    analysisLoaded = true;
  }

  public void loadWSPInfo() {
    wspLoader = new WSPLoader(panelsToRun.toArray(new Panel[panelsToRun.size()]));
    int success = wspLoader.loadWorkspaces(wspDir);
    if (success > 0) {
      log.reportTime("Loaded WSP data from " + success + " files.");
      wspLoaded = true;
    } else {
      log.reportError("Loading WSP data failed - please check log and try again.");
    }
  }

  public void loadOutliers() {
    if (outliersFile != null) {
      String[] out = HashVec.loadFileToStringArray(outliersFile, false, null, false);
      for (String s : out) {
        Panel pnl = checkPanel(s.toLowerCase());
        if (pnl != null) {
          panelData.get(pnl).outliers.add(s);
        } else {
          log.reportWarning("Couldn't identify panel for outlier " + s);
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

    for (Panel p : panelsToRun) {
      PanelData pd = panelData.get(p);
      List<String> files = new ArrayList<>(pd.filePathMap.keySet());
      if (!pd.data.isEmpty()) {
        files.retainAll(pd.data.keySet());
      }
      List<String> outliers = new ArrayList<>(pd.outliers);
      outliers.retainAll(files);
      files.removeAll(outliers);
      int cnt = files.isEmpty() ? 0 : (int) Math.max(1, (int) files.size() * samplingPct);
      log.reportTime("Sampling " + cnt + " " + p.getName() + " files");
      if (samplingPct > .95) {
        pd.sampling.addAll(files);
      } else {
        Random rand = new Random();
        for (int i = 0; i < cnt; i++) {
          if (!outliers.isEmpty() && OUTLIER_SAMPLE_EVERY > 0
              && (i > 0 && i % OUTLIER_SAMPLE_EVERY == 0)) {
            pd.sampling.add(outliers.get(rand.nextInt(outliers.size())));
          } else {
            pd.sampling.add(files.get(rand.nextInt(files.size())));
          }
        }
      }
    }

    new File(outDir).mkdirs();
    PrintWriter writer;
    for (Panel p : panelsToRun) {
      if (panelData.get(p).sampling.size() > 0) {
        writer = Files.getAppropriateWriter(outDir + ext.replaceWithLinuxSafeCharacters(p.getName())
                                            + ".files.txt");
        for (String s : panelData.get(p).sampling) {
          writer.println(s);
        }
        writer.flush();
        writer.close();
      }
    }

    filesSampled = true;
  }

  public void doDataSampling() {
    if (!checkReady() && !filesSampled) {
      return;
    }

    final Map<Panel, ConcurrentLinkedQueue<SampleNode>> panelQueues = new HashMap<>();
    for (Panel p : panelsToRun) {
      panelQueues.put(p, new ConcurrentLinkedQueue<>());
    }

    Map<Panel, PrintWriter> sampWspMatch = new HashMap<>();
    for (Panel p : panelsToRun) {
      sampWspMatch.put(p,
                       Files.getAppropriateWriter(outDir + "/matchFile_"
                                                  + ext.replaceWithLinuxSafeCharacters(p.getName())
                                                  + ".xln"));
    }

    if (highPriority != null) {
      for (Panel p : panelsToRun) {
        for (String s : panelData.get(p).sampling) {
          if (highPriority.contains(ext.replaceWithLinuxSafeCharacters(ext.removeDirectoryInfo(s)))) {
            SampleNode sn = wspLoader.getPanelNodes(p).get(s);
            if (sn != null) {
              sn.fcsFile = panelData.get(p).filePathMap.get(s);
              if (panelQueues.get(p).contains(sn)) {
                log.report("Duplicate " + p.getName() + " sampleNode: " + sn.fcsFile + " | "
                           + sn.wspFile + " | " + sn.id);
                continue;
              }
              sampWspMatch.get(p).println(s + "\t" + sn.wspFile + "\t" + sn.fcsFile);
              panelQueues.get(p).add(sn);
            }
            break;
          }
        }
      }
    }

    for (Panel p : panelsToRun) {
      for (String s : panelData.get(p).sampling) {
        if (highPriority != null
            && highPriority.contains(ext.replaceWithLinuxSafeCharacters(ext.removeDirectoryInfo(s)))) {
          continue;
        }
        if (lowPriority != null
            && lowPriority.contains(ext.replaceWithLinuxSafeCharacters(ext.removeDirectoryInfo(s)))) {
          continue;
        }
        SampleNode sn = wspLoader.getPanelNodes(p).get(s);
        if (sn != null) {
          sn.fcsFile = panelData.get(p).filePathMap.get(s);
          if (panelQueues.get(p).contains(sn)) {
            log.report("Duplicate " + p.getName() + " sampleNode: " + sn.fcsFile + " | "
                       + sn.wspFile + " | " + sn.id);
            continue;
          }
          sampWspMatch.get(p).println(s + "\t" + sn.wspFile + "\t" + sn.fcsFile);
          panelQueues.get(p).add(sn);
        } else {
          log.reportWarning("Couldn't find WSP match for sample: " + s);
        }
      }
    }

    if (lowPriority != null) {
      for (Panel p : panelsToRun) {
        for (String s : panelData.get(p).sampling) {
          if (lowPriority.contains(ext.replaceWithLinuxSafeCharacters(ext.removeDirectoryInfo(s)))) {
            SampleNode sn = wspLoader.getPanelNodes(p).get(s);
            if (sn != null) {
              sn.fcsFile = panelData.get(p).filePathMap.get(s);
              if (panelQueues.get(p).contains(sn)) {
                log.report("Duplicate " + p.getName() + " sampleNode: " + sn.fcsFile + " | "
                           + sn.wspFile + " | " + sn.id);
                continue;
              }
              sampWspMatch.get(p).println(s + "\t" + sn.wspFile + "\t" + sn.fcsFile);
              panelQueues.get(p).add(sn);
            }
          }
        }
      }
    }

    for (PrintWriter pw : sampWspMatch.values()) {
      pw.close();
    }

    log.reportTime("Finished matching sample files to workspace files.");

    ExecutorService serve = Executors.newFixedThreadPool(panelsToRun.size());
    int proc = Math.max(1, Runtime.getRuntime().availableProcessors() / panelsToRun.size());
    List<Runnable> runn = new ArrayList<>();
    for (Panel p : panelsToRun) {
      final ThreadPoolExecutor threadPool = new ThreadPoolExecutor(proc, proc, 0L,
                                                                   TimeUnit.MILLISECONDS,
                                                                   new LinkedBlockingQueue<Runnable>());
      AbstractPipelineRunnable run = new AbstractPipelineRunnable(threadPool, processorFactory,
                                                                  panelQueues.get(p), log);
      runn.add(run);
    }
    for (Runnable run : runn) {
      serve.submit(run);
    }
    serve.shutdown();
    try {
      serve.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
    } catch (InterruptedException e) {}

    for (Runnable run : runn) {
      processorFactory.cleanup(run);
    }
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
        } catch (InterruptedException e) {}
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
      } catch (InterruptedException e) {}
    }
    try {
      threadPool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
    } catch (InterruptedException e1) {}
  }

}
