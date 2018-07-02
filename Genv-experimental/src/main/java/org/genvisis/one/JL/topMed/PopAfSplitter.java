package org.genvisis.one.JL.topMed;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.genvisis.CLI;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.qsub.Qsub;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.manage.dictionary.SAMSequenceDictionarySplitter;
import org.genvisis.seq.manage.dictionary.SAMSequenceDictionarySplitter.Query;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.vcf.VCFFileReader;

public class PopAfSplitter {

  public static void main(String[] args) {
    CLI c = new CLI(PopAfSplitter.class);
    c.addArg("vcfDir", "directory containing .vcf.gz files");
    c.addArg("popFile",
             "File containing a header with \"IID\" and \"Population\", every sample in this file must be"
                        + " in the vcf. Allele frequency will be computed for each population",
             true);
    c.addArg(CLI.ARG_OUTDIR, CLI.DESC_OUTDIR, false);
    c.addArg("batches", "number of batches", "150");
    c.addArg("genvisis", "genvisis");

    c.parseWithExit(args);

    String outDir = c.get(CLI.ARG_OUTDIR);
    new File(outDir).mkdirs();
    int memoryInMb = 12000;
    int wallTime = 64;
    int numBatches = c.getI("batches");

    String[] allFiles = Files.listFullPaths(c.get("vcfDir"), ".vcf.gz");

    String vcf = allFiles[0];

    //    String java = "";

    //    freeze.5b.chr21.pass

    Map<String, String> filechrmap = new HashMap<>();
    for (int i = 0; i < 24; i++) {

      String chr = Positions.getChromosomeUCSC(i, true);
      for (String file : allFiles) {
        if (file.contains("freeze.5b." + chr + ".pass")
            || file.contains("freeze.5b." + chr + ".aims..pass")) {
          if (filechrmap.containsKey(chr)) {
            throw new IllegalArgumentException("multiple");
          } else {
            filechrmap.put(chr, file);
          }
        }
      }
    }
    VCFFileReader reader = new VCFFileReader(new File(vcf));

    reader.close();
    List<SAMSequenceRecord> recordsToUse = new ArrayList<>();

    for (SAMSequenceRecord record : reader.getFileHeader().getSequenceDictionary().getSequences()) {
      if (filechrmap.containsKey(record.getSequenceName())) {
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
    baseCommands.add(c.get("genvisis"));

    baseCommands.add("seq.analysis.freq.VCFPopulationFrequency");
    baseCommands.add(CLI.ARG_OUTDIR + "=" + outDir);
    baseCommands.add("popFile=" + c.get("popFile"));
    baseCommands.add(CLI.ARG_THREADS + "=1");

    for (Query q : queries) {
      List<String> currentCmd = new ArrayList<>(baseCommands);
      if (!filechrmap.containsKey(q.getChrom())) {
        throw new IllegalArgumentException("invalid chrom " + q.getChrom());
      }
      String vcfCurrent = filechrmap.get(q.getChrom());
      String outPut = outDir + VCFOps.getAppropriateRoot(vcfCurrent, true) + "_" + q.getChrom()
                      + "_" + q.getStart() + "-" + q.getEnd() + ".pbs";
      currentCmd.add(CLI.ARG_VCF + "=" + vcfCurrent);
      currentCmd.add("region=" + q.getChrom() + ":" + q.getStart() + "-" + q.getEnd());

      Qsub.qsub(outPut, ArrayUtils.toStr(currentCmd, " "), memoryInMb, wallTime, 1);
    }
  }

}
