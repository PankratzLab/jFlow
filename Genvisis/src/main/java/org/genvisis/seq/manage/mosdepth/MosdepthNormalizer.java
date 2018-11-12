package org.genvisis.seq.manage.mosdepth;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import org.genvisis.CLI;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.common.SerializedFiles;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.BEDFileReader;
import htsjdk.tribble.bed.BEDFeature;

public class MosdepthNormalizer {

  public static String FILE_EXT = ".norm.gz";
  public static String INDEX_EXT = ".idx";

  Logger log = new Logger();
  int numThreads;

  String mosSrcDir;
  String mosSrcExt;
  String[] mosdepthFiles;

  String outputDir;

  public MosdepthNormalizer(String dir, String ext, String outputDir, int threads, Logger log) {
    this.mosSrcDir = org.genvisis.common.ext.verifyDirFormat(dir);
    this.mosSrcExt = ext;
    this.mosdepthFiles = Files.list(dir, ext);
    this.outputDir = org.genvisis.common.ext.verifyDirFormat(outputDir);
    this.numThreads = threads;
    this.log = log;
  }

  public void run() {
    long t1 = System.nanoTime();
    ExecutorService exec = Executors.newFixedThreadPool(numThreads);
    for (String s : mosdepthFiles) {
      exec.submit(new MosdepthNormWorker(mosSrcDir + s, outputDir + getOutputFileName(s)));
    }
    exec.shutdown();
    try {
      exec.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
    } catch (InterruptedException e) {
      log.reportException(e);
    }
    log.reportTime("Normalized " + mosdepthFiles.length + " mosdepth files in "
                   + ext.getTimeElapsedNanos(t1));
  }

  private static String getOutputFileName(String input) {
    return ext.rootOf(input, true) + FILE_EXT;
  }

  private static String getIndexFileName(String output) {
    return ext.rootOf(output, false) + INDEX_EXT;
  }

  class NormalizedMosdepthReader implements AutoCloseable {

    private static final int BYTES_PER = 25;

    String file;
    RandomAccessFile raf;
    Map<Segment, Long> index;

    public NormalizedMosdepthReader(String filename) throws IOException {
      this.file = filename;
      this.raf = new RandomAccessFile(this.file, "r");
      buildIndex();
    }

    @SuppressWarnings("unchecked")
    private void buildIndex() throws IOException {
      String idxFile = getIndexFileName(file);
      if (Files.exists(idxFile)) {
        this.index = (Map<Segment, Long>) SerializedFiles.readSerial(idxFile);
      } else {
        rebuildIndex();
      }
    }

    private void rebuildIndex() throws IOException {
      this.raf.seek(0);
      this.index = new HashMap<Segment, Long>();
      boolean con = true;
      long ind = 0;
      while (con) {
        this.index.put(new Segment(this.raf.readByte(), (int) this.raf.readLong(),
                                   (int) this.raf.readLong()),
                       ind++);
        this.raf.skipBytes(8); // skip the double value 
        if (this.raf.getFilePointer() == this.raf.length()) {
          con = false;
        }
      }
      SerializedFiles.writeSerial(index, getIndexFileName(file));
      this.raf.seek(0);
    }

    public double get(Segment seg) throws IOException {
      if (index.containsKey(seg)) {
        long ind = index.get(seg);
        long pos = BYTES_PER * ind;
        if (pos != raf.getFilePointer()) {
          raf.seek(pos);
        }
        byte chr = raf.readByte();
        long stt = raf.readLong();
        long stp = raf.readLong();
        double d = raf.readDouble();
        if (chr != seg.getChr()) {
          throw new RuntimeException(); // TODO error msgl;
        }
        if (((int) stt) != seg.getStart()) {
          throw new RuntimeException(); // TODO error msg
        }
        if (((int) stp) != seg.getStop()) {
          throw new RuntimeException(); // TODO error msg
        }
        return d;
      } else {
        throw new RuntimeException(); // TODO error msg
      }
    }

    @Override
    public void close() throws Exception {
      this.raf.close();
      this.index = null;
    }

  }

  class MosdepthNormWorker implements Runnable {

    private String in;
    private String out;

    public MosdepthNormWorker(String inputFile, String outputFile) {
      this.in = inputFile;
      this.out = outputFile;
      System.out.println("Parsing " + in + " to " + out);
    }

    @Override
    public void run() {
      try {
        RandomAccessFile raf = null;
        new File(ext.parseDirectoryOfFile(out)).mkdirs();
        try {
          raf = new RandomAccessFile(out, "rw");
        } catch (FileNotFoundException e) {
          // TODO Auto-generated catch block
          e.printStackTrace();
        }
        BEDFileReader reader = new BEDFileReader(in, false);
        Map<Segment, Long> segIndexMap = new HashMap<>();
        long ind = 0;

        int sz = Files.countLines(in, 0);
        double median = ArrayUtils.median(reader.iterator().stream().map(bf -> bf.getScore()), sz);

        for (int i = 0; i <= 27; i++) {
          List<BEDFeature> iter = reader.query(Positions.chromosomeNumberInverse(i), 0,
                                               Integer.MAX_VALUE)
                                        .toList();
          for (BEDFeature bf : iter) {
            float scor = bf.getScore();
            double newscor = scor / median;
            try {
              raf.write((byte) i);
              raf.writeLong((long) bf.getStart());
              raf.writeLong((long) bf.getEnd());
              raf.writeDouble(newscor);
              segIndexMap.put(new Segment(bf.getContig(), bf.getStart(), bf.getEnd()), ind);
            } catch (IOException e) {
              // TODO Auto-generated catch block
              e.printStackTrace();
            }
            ind++;
          }
        }

        SerializedFiles.writeSerial(segIndexMap, getIndexFileName(out));
        reader.close();
        try {
          raf.close();
        } catch (IOException e) {
          // TODO Auto-generated catch block
          e.printStackTrace();
        }
      } catch (Exception e) {
        e.printStackTrace();
      }
    }

  }

  public static void main(String[] args) {
    CLI cli = new CLI(MosdepthNormalizer.class);
    cli.addArg(CLI.ARG_INDIR, CLI.DESC_INDIR);
    cli.addArg("ext", "Mosdepth output file extension");
    cli.addArg(CLI.ARG_OUTDIR, CLI.DESC_OUTDIR);
    cli.addArg(CLI.ARG_THREADS, CLI.DESC_THREADS);

    cli.parseWithExit(args);

    new MosdepthNormalizer(cli.get(CLI.ARG_INDIR), cli.get("ext"), cli.get(CLI.ARG_OUTDIR),
                           cli.getI(CLI.ARG_THREADS), new Logger()).run();
  }

}
