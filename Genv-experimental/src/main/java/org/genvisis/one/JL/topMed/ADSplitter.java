
package org.genvisis.one.JL.topMed;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.manage.dictionary.SAMSequenceDictionarySplitter;
import org.genvisis.seq.manage.dictionary.SAMSequenceDictionarySplitter.Query;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.CLI;
import org.pankratzlab.shared.filesys.Positions;
import org.pankratzlab.shared.qsub.Qsub;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.vcf.VCFFileReader;

public class ADSplitter {

  public static void main(String[] args) {
    CLI c = new CLI(PopAfSplitter.class);
    c.addArg(CLI.ARG_VCF, "vcf file");
    c.addArg("crams", "file of .cram file", true);
    c.addArg(CLI.ARG_OUTDIR, CLI.DESC_OUTDIR, false);
    c.addArg("batches", "number of batches", "150");
    c.addArg(CLI.ARG_REFERENCE_GENOME, "reference genome", "ref.fa");
    c.addArg(CLI.ARG_THREADS, "threads for GATK", "6");

    c.addArg("gatk", "gatk");

    c.parseWithExit(args);

    String outDir = c.get(CLI.ARG_OUTDIR);
    new File(outDir).mkdirs();
    int memoryInMb = 32000;
    int wallTime = 64;
    int numBatches = c.getI("batches");
    String vcf = c.get(CLI.ARG_VCF);

    String[] allCrams = HashVec.loadFileToStringArray(c.get("crams"), false, new int[] {0}, false);

    VCFFileReader reader = new VCFFileReader(new File(vcf));

    reader.close();
    List<SAMSequenceRecord> recordsToUse = new ArrayList<>();

    for (SAMSequenceRecord record : reader.getFileHeader().getSequenceDictionary().getSequences()) {
      if (Positions.chromosomeNumber(record.getSequenceName()) > 0) {
        recordsToUse.add(record);
      }
    }

    SAMSequenceDictionary samSequenceDictionary = new SAMSequenceDictionary(recordsToUse);

    List<Query> queries = SAMSequenceDictionarySplitter.parseToQueries(samSequenceDictionary,
                                                                       numBatches, new Logger());

    List<String> baseCommands = new ArrayList<>();
    baseCommands.add("java");
    baseCommands.add("-Xmx" + memoryInMb + "m");
    baseCommands.add("-jar ");
    baseCommands.add(c.get("gatk"));
    baseCommands.add("-T");
    baseCommands.add("VariantAnnotator");
    baseCommands.add("-R");
    baseCommands.add(c.get(CLI.ARG_REFERENCE_GENOME));
    baseCommands.add("-V");
    baseCommands.add(vcf);
    baseCommands.add("-A");
    baseCommands.add("DepthPerAlleleBySample");
    baseCommands.add("-nt");
    baseCommands.add(Integer.toString(c.getI(CLI.ARG_THREADS)));
    for (String cram : allCrams) {
      baseCommands.add("-I");
      baseCommands.add(cram);
    }

    for (Query q : queries) {
      List<String> currentCmd = new ArrayList<>(baseCommands);

      String outPut = outDir + VCFOps.getAppropriateRoot(vcf, true) + "_" + q.getChrom() + "_"
                      + q.getStart() + "-" + q.getEnd() + ".pbs";
      String outvcf = outDir + VCFOps.getAppropriateRoot(vcf, true) + "_" + q.getChrom() + "_"
                      + q.getStart() + "-" + q.getEnd() + ".vcf.gz";

      currentCmd.add("-L");
      currentCmd.add(q.getChrom() + ":" + q.getStart() + "-" + q.getEnd());
      currentCmd.add("-o");
      currentCmd.add(outvcf);
      Qsub.qsub(outPut, ArrayUtils.toStr(currentCmd, " "), memoryInMb, wallTime,
                c.getI(CLI.ARG_THREADS));
    }
  }

}
