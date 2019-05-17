package org.genvisis.one.ben;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.ext;
import org.pankratzlab.utils.gwas.bgen.BGENReader;
import org.pankratzlab.utils.gwas.bgen.BGENTools;

public class BGENIndexer {

  // /**
  // * Read and dump, without saving, record metadata
  // *
  // * @param mapFileOut File to which to dump metadata values
  // *
  // * @throws IOException
  // */
  // public static void dumpMap(BGENReader reader, String mapFileOut) throws IOException {
  // long skip = 0;
  // try {
  // skip = Files.exists(ext.rootOf(mapFileOut, false) + ".done")
  // ? Long.parseLong(HashVec.loadFileToStringArray(ext.rootOf(mapFileOut,
  // false)
  // + ".done",
  // false,
  // null,
  // false)[0])
  // : 0;
  // } catch (NumberFormatException e) {
  // skip = 0;
  // }
  // PrintWriter writer = new PrintWriter(new GZIPOutputStream(new FileOutputStream(mapFileOut,
  // false)));
  // Iterator<BGENRecord> mapIter = new BGENIterators.BGENIterator(reader, reader.raf, false, false,
  // false, skip);
  // long written = 0;
  // StringBuilder sb;
  // while (mapIter.hasNext()) {
  // BGENRecordMetaData md = mapIter.next().metaData;
  // sb = new StringBuilder();
  // sb.append(md.id).append("\t");
  // sb.append(md.rsId).append("\t");
  // sb.append(md.chr).append("\t");
  // sb.append(md.pos).append("\t");
  //
  // for (int i = 0; i < md.alleles.length; i++) {
  // sb.append(md.alleles[i]);
  // if (i < md.alleles.length - 1) {
  // sb.append(",");
  // }
  // }
  // sb.append("\t");
  //
  // sb.append(md.N).append("\t");
  // sb.append(md.ptrByt).append("\t");
  // sb.append(md.lenByt).append("\t");
  // sb.append(md.blockLength);
  // writer.println(sb.toString());
  // written++;
  // if (written > 0 && written % 2500 == 0) {
  // writer.flush();
  // Files.write(Long.toString(written), ext.rootOf(mapFileOut, false) + ".done");
  // }
  // }
  // writer.close();
  // }

  private static void indexOld() throws IOException {
    String dir = "/scratch.global/bb/all/EGAD00010001225/001/";
    String[] files = new File(dir).list((File f, String s) -> {
      return s.endsWith(".bgen");
    });
    for (String file : files) {
      final String outFile = dir + ext.rootOf(file, true) + ".bgm.gz";
      if (!Files.exists(outFile)) {
        System.out.println("Submitting mapping job for " + file);
        BGENReader reader = BGENReader.open(dir + file, true);
        BGENTools.serializeMapInfo(reader, outFile);
        reader.close();
        System.out.println("Finished with " + file);
      }
    }
  }

  private static void indexNew2() throws InterruptedException {
    String dir = "/scratch.global/bb/all/EGAD00010001225/001/";
    String[] files = new File(dir).list((File f, String s) -> {
      return s.endsWith(".bgen");
    });
    ArrayList<Runnable> runners = new ArrayList<>();
    for (String file : files) {
      final String outFile = dir + ext.rootOf(file, true) + ".bgm.gz";
      boolean skip = Files.exists(outFile);
      if (!skip) {
        runners.add(new Runnable() {

          @Override
          public void run() {
            try {
              long t1 = System.nanoTime();
              System.out.println("Submitting mapping job for " + file);
              BGENReader reader = BGENReader.open(dir + file, true);
              System.out.println("Read map info for " + file + " in " + ext.getTimeElapsedNanos(t1)
                                 + "... serializing now ...");
              t1 = System.nanoTime();
              BGENTools.serializeMapInfo(reader, outFile);
              reader.close();
              System.out.println("Finished serializing " + file + " in "
                                 + ext.getTimeElapsedNanos(t1));
            } catch (IOException e) {
              e.printStackTrace();
            }
          }
        });
      }
    }
    int proc = Math.min(runners.size(), Runtime.getRuntime().availableProcessors());
    ExecutorService exec = Executors.newFixedThreadPool(proc);
    for (Runnable r : runners) {
      exec.submit(r);
    }
    exec.shutdown();
    exec.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
  }

  // private static void indexNew1() throws InterruptedException {
  // String dir = "/scratch.global/bb/all/EGAD00010001225/001/";
  // String[] files = new File(dir).list((File f, String s) -> {
  // return s.endsWith(".bgen");
  // });
  // ArrayList<Runnable> runners = new ArrayList<>();
  // for (String file : files) {
  // final String outFile = dir + ext.rootOf(file, true) + ".bgm";
  // if (!(Files.exists(outFile) || Files.exists(outFile + ".temp.gz"))) {
  // runners.add(new Runnable() {
  // @Override
  // public void run() {
  // System.out.println("Submitting mapping job for " + file);
  // BGENReader reader;
  // try {
  // reader = BGENReader.open(dir + file, false);
  // reader.dumpMap(outFile + ".temp.gz");
  // reader.close();
  // System.out.println("Finished with " + file);
  // } catch (IOException e) {
  // e.printStackTrace();
  // }
  // }
  // });
  // }
  // }
  // int proc = Math.min(runners.size(), Runtime.getRuntime().availableProcessors());
  // ExecutorService exec = Executors.newFixedThreadPool(proc);
  // for (Runnable r : runners) {
  // exec.submit(r);
  // }
  // exec.shutdown();
  // exec.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
  // }

  public static void main(String[] args) throws InterruptedException, IOException {
    // indexOld();
    // indexNew1();
    indexNew2();
  }

}
