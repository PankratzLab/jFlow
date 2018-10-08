package org.genvisis.one.JL.topMed;

import java.io.File;
import java.util.HashSet;
import org.genvisis.seq.manage.VCFOps;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.WorkerHive;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.CLI;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

/**
 * Geared toward selecting matching rsIDs from multiple *chr*.vcf.gz files, and creating single
 * merged and sorted .vcf
 */
public class CollapseToRsIds {

  private static class ExtractResult {

    private final String outVcf;
    private int foundTotal;
    private int scannedTotal;

    /**
     * @param outVcf
     * @param foundTotal
     * @param scannedTotal
     */
    private ExtractResult(String outVcf, int foundTotal, int scannedTotal) {
      super();
      this.outVcf = outVcf;
      this.foundTotal = foundTotal;
      this.scannedTotal = scannedTotal;
    }
  }

  private static void run(String inputDir, String rsIDFile, String outDir, String preChr,
                          String postChr, int threads) {
    new File(outDir).mkdirs();

    Logger log = new Logger(outDir + "merge.log");

    VCFFileReader readerTmp = new VCFFileReader(new File(Files.listFullPaths(inputDir,
                                                                             VCFOps.VCF_EXTENSIONS.GZIP_VCF.getLiteral())[0]),
                                                false);
    HashSet<String> rsIds = HashVec.loadFileToHashSet(rsIDFile, false);
    log.reportTimeInfo("Loaded " + rsIds.size() + " ids from " + rsIDFile);

    VCFHeader header = readerTmp.getFileHeader();
    SAMSequenceDictionary newDictionary = readerTmp.getFileHeader().getSequenceDictionary();
    readerTmp.close();
    WorkerHive<ExtractResult> hive = new WorkerHive<>(threads, 10, log);
    int numProcessing = 0;
    for (SAMSequenceRecord record : newDictionary.getSequences()) {
      String vcf = inputDir + preChr + record.getSequenceName() + postChr;
      if (!Files.exists(vcf)) {
        log.reportTimeWarning("skipping " + vcf);
      } else if (!Files.exists(vcf + ".tbi")) {
        log.reportTimeWarning("skipping " + vcf + " due to missing index file");
      } else {
        log.reportTimeInfo("Submitting vcf " + vcf);
        hive.addCallable(() -> extract(outDir, preChr, postChr, log, rsIds, header, newDictionary,
                                       record, vcf));
        numProcessing++;

      }

    }
    log.reportTimeInfo("Processing total of " + numProcessing + " vcf files");
    hive.execute(true);
    hive.getResults();

    log.reportTimeInfo("Finished scanning files");
    String outVcf = outDir + preChr + ext.rootOf(rsIDFile) + postChr;
    VariantContextWriter writer = VCFOps.initWriter(outVcf, VCFOps.DEFUALT_WRITER_OPTIONS,
                                                    newDictionary);
    writer.writeHeader(header);
    for (SAMSequenceRecord record : newDictionary.getSequences()) {
      String vcf = getOut(outDir, preChr, postChr, record);
      log.reportTimeInfo("Processing vcf " + vcf);
      if (!Files.exists(vcf)) {
        log.reportTimeWarning("skipping " + vcf);
      } else {
        VCFFileReader reader = new VCFFileReader(new File(vcf));
        for (VariantContext vc : reader) {
          writer.add(vc);
        }
        reader.close();
      }
    }
    writer.close();
  }

  private static ExtractResult extract(String outDir, String preChr, String postChr, Logger log,
                                       HashSet<String> rsIds, VCFHeader header,
                                       SAMSequenceDictionary newDictionary,
                                       SAMSequenceRecord record, String vcf) {
    int foundTotal = 0;
    int scannedTotal = 0;
    String outVcf = getOut(outDir, preChr, postChr, record);
    if (!Files.exists(outVcf)) {
      VariantContextWriter writer = VCFOps.initWriter(outVcf, VCFOps.DEFUALT_WRITER_OPTIONS,
                                                      newDictionary);
      VCFFileReader reader = new VCFFileReader(new File(vcf), false);
      writer.writeHeader(header);
      for (VariantContext vc : reader) {
        scannedTotal++;
        if (rsIds.contains(vc.getID())) {
          writer.add(vc);
          foundTotal++;
        }
        if (scannedTotal % 10000 == 0) {
          log.reportTimeInfo("Currently scanning file " + vcf);
          log.reportTimeInfo("Scanned " + scannedTotal + " variants, found " + foundTotal
                             + " variants");

        }
      }
      reader.close();
      Files.write("scanned", vcf + ".scanned");
      log.reportTimeInfo("Scanned " + scannedTotal + " variants, found " + foundTotal
                         + " variants");
      writer.close();
    }
    return new ExtractResult(outVcf, foundTotal, scannedTotal);
  }

  private static String getOut(String outDir, String preChr, String postChr,
                               SAMSequenceRecord record) {
    return outDir + preChr + record.getSequenceName() + ".aims" + postChr;
  }

  public static void main(String[] args) {
    CLI c = new CLI(CollapseToRsIds.class);
    c.addArg(CLI.ARG_INDIR, "directory of .bcf files", true);
    c.addArg("rsIDFile", "file containing rsIDs", true);
    c.addArg(CLI.ARG_OUTDIR, CLI.DESC_OUTDIR, true);

    c.parseWithExit(args);

    run(c.get(CLI.ARG_INDIR), c.get("rsIDFile"), c.get(CLI.ARG_OUTDIR), "freeze.6a.",
        ".pass_and_fail.gtonly.minDP0.vcf.gz", 24);
  }
}
