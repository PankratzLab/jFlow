package org.genvisis.cnv.seq.manage;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.atomic.AtomicInteger;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.seq.manage.BamSample.NORMALIZATON_METHOD;
import org.genvisis.seq.ReferenceGenome;
import org.genvisis.seq.SeqVariables.ASSEMBLY_NAME;
import org.genvisis.seq.manage.BamOps;
import org.genvisis.seq.manage.BamPile;
import org.genvisis.seq.manage.BamSegPileUp;
import org.genvisis.seq.manage.BamSegPileUp.BamPileResult;
import org.genvisis.seq.qc.FilterNGS;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.WorkerTrain.AbstractProducer;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.filesys.Segment;
import com.google.common.collect.Sets;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class PileupProducer extends AbstractProducer<BamPileResult> {

  private final AtomicInteger index;
  private final String[] bamFiles;
  private final ASSEMBLY_NAME aName;
  private final String serDir;
  private final Logger log;
  private final Segment[] pileSegs;
  private final FilterNGS filterNGS;
  private final ReferenceGenome referenceGenome;
  private final Map<String, QueryInterval[]> qiMap;
  private final Project proj;
  private final boolean writeSer;
  private final boolean compileSamp;
  private final NORMALIZATON_METHOD normMethod;

  public PileupProducer(Project proj, String[] bamFiles, String serDir,
                        ReferenceGenome referenceGenome, FilterNGS filterNGS, Segment[] pileSegs,
                        ASSEMBLY_NAME aName, boolean writeSerial, boolean compileSample,
                        NORMALIZATON_METHOD normMethod, Logger log) {
    super();
    this.proj = proj;
    this.writeSer = writeSerial;
    this.compileSamp = compileSample;
    this.normMethod = normMethod;
    this.index = new AtomicInteger(0);
    this.bamFiles = bamFiles;
    this.aName = aName;
    this.serDir = serDir;
    this.referenceGenome = referenceGenome;
    this.log = log;
    this.pileSegs = pileSegs;
    this.filterNGS = filterNGS;
    new File(serDir).mkdirs();
    qiMap = new HashMap<String, QueryInterval[]>();
    setupQIs();
  }

  private void setupQIs() {
    long t1 = System.nanoTime();
    Map<SAMSequenceDictionary, QueryInterval[]> dicts = new HashMap<>();
    for (String bam : bamFiles) {
      try (SamReader reader = BamOps.getDefaultReader(bam, ValidationStringency.STRICT,
                                                      Sets.immutableEnumSet(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES))) {
        SAMSequenceDictionary sd = reader.getFileHeader().getSequenceDictionary();
        QueryInterval[] qi = null;
        for (SAMSequenceDictionary d : dicts.keySet()) {
          boolean same = false;
          try {
            sd.assertSameDictionary(d);
            same = true;
          } catch (AssertionError e) {
            // not the same
          }
          if (same) {
            qi = dicts.get(d);
            break;
          }
        }
        if (qi == null) {
          qi = BamOps.convertSegsToQI(pileSegs, reader.getFileHeader(), 0, true, aName.addChr(),
                                      log);
          dicts.put(sd, qi);
        }
        qiMap.put(bam, qi);
      } catch (IOException e) {
        throw new IllegalStateException("Couldn't read header of BAM file: " + bam, e);
      }
    }
    dicts = null;
    log.reportTime("Constructed query intervals in " + ext.getTimeElapsedNanos(t1));
  }

  @Override
  public boolean hasNext() {
    return index.get() < bamFiles.length;
  }

  @Override
  public Callable<BamPileResult> next() {
    final String bamFile = bamFiles[index.getAndIncrement()];
    return new Callable<BamSegPileUp.BamPileResult>() {

      @Override
      public BamPileResult call() throws Exception {
        return processBamFile(proj, bamFile, serDir, referenceGenome, pileSegs, qiMap.get(bamFile),
                              filterNGS, aName, log, writeSer, compileSamp, normMethod);
      }
    };
  }

  @Override
  public void shutdown() {

  }

  public static BamPileResult processBamFile(Project proj, String bamFile, String serDir,
                                             ReferenceGenome referenceGenome, Segment[] pileSegs,
                                             QueryInterval[] qi, FilterNGS filterNGS,
                                             ASSEMBLY_NAME aName, Logger log, boolean writeSerial,
                                             boolean compileSample,
                                             NORMALIZATON_METHOD normMethod) throws Exception {
    String ser = serDir + ext.rootOf(bamFile) + ".ser";
    if (!Files.exists(ser)) {
      BamOps.verifyIndex(bamFile, log);
      BamSegPileUp bamSegPileUp = new BamSegPileUp(bamFile, referenceGenome, pileSegs, qi,
                                                   filterNGS, aName, log);
      BamPile[] bamPilesFinal = bamSegPileUp.pileup();
      int numMiss = 0;
      for (BamPile element : bamPilesFinal) {
        numMiss += element.getNumBasesWithMismatch();
      }
      log.reportTimeInfo(bamFile + " had " + numMiss + " mismatched bases");
      if (writeSerial) {
        BamPile.writeSerial(bamPilesFinal, ser);
      }
      if (compileSample) {
        BamSample bamSample = new BamSample(proj, bamFile, bamPilesFinal, normMethod);
        bamSample.writeSample(proj.getMarkerSet().getFingerprint());
      }
    }
    return new BamPileResult(bamFile, ser);
  }
}
