package org.genvisis.seq.manage;

import java.io.File;
import java.util.ArrayList;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.WorkerHive;
import org.genvisis.common.ext;
import org.genvisis.seq.analysis.Blast;
import org.genvisis.seq.analysis.Blast.BlastWorker;
import org.genvisis.seq.analysis.Blast.FastaEntry;

import htsjdk.tribble.annotation.Strand;

public class Adapter {
  public static final String BARCODE = "barcode";

  /**
   * blast a sequence (typically soft clipped) against a blast database of adapter content
   */
  public static String[] blast(int blastWordSize, ArrayList<Adapter> adapters, String[] sequences,
                               String outputDir, String root, int numThreads, Logger log) {
    new File(outputDir).mkdirs();
    String db = outputDir + "adapter.db.fa";
    String containsFile = db + ".adapterContent";
    String[] names = getAllNames(adapters);

    if (!Files.exists(containsFile) || !Files.exists(db)
        || !ext.containsAll(names, HashVec.loadFileToStringArray(containsFile, false, new int[] {0},
                                                                 true))) {
      String contains = "";
      String adapterFasta = "";
      for (Adapter adapter : adapters) {
        contains += adapter.getName() + "\n";
        adapterFasta += ">" + adapter.getAdaptorSequence() + "_" + adapter.getName() + "\n";
        adapterFasta += adapter.getAdaptorSequence() + "\n";
        // adapterFasta += ">" + adapter.getName() + "_Reverse_Strand\n";
        // adapterFasta += adapter.getAdaptorRSSequence() + "\n";
      }
      if (Files.exists(db)) {
        for (String db_EXT : Blast.DB_EXTs) {
          if (Files.exists(db + db_EXT)) {
            new File(db + db_EXT).delete();
          }
        }
      }
      Files.write(adapterFasta, db);
      Files.write(contains, containsFile);
    }
    Blast blast = new Blast(db, blastWordSize, 0, log, true, true);
    blast.setEvalue(10000);
    FastaEntry[] fastaEntries = new FastaEntry[sequences.length];
    for (int i = 0; i < sequences.length; i++) {
      fastaEntries[i] = new FastaEntry(sequences[i] + "_softClippedSequence", sequences[i]);
    }
    ArrayList<FastaEntry[]> splits = Array.splitUpArray(fastaEntries, numThreads, log);
    ArrayList<BlastWorker> workers = new ArrayList<Blast.BlastWorker>();
    String[] tmps = new String[splits.size()];
    if (fastaEntries != null && fastaEntries.length > 0) {
      for (int i = 0; i < splits.size(); i++) {

        String tmp = outputDir + root + "adapterBlast_" + i + ".gz";
        workers.add(new BlastWorker(blast, splits.get(i), tmp));
        tmps[i] = tmp;
      }
    }
    WorkerHive<Blast.BlastResultsSummary[]> hive =
        new WorkerHive<Blast.BlastResultsSummary[]>(numThreads, 10, log);
    hive.addCallables(workers.toArray(new BlastWorker[workers.size()]));
    hive.setReportEvery(1);
    hive.execute(true);
    return tmps;
  }

  public static Adapter getAgilentQXTAdapter() {
    String seq = "CTGTCTCTTGATCACA";
    String name = "AgilentQXT_Adapter";
    String description = "Adapter sequence for Agilent QXT protocol";
    return new Adapter(seq, name, description);
  }

  private static String[] getAllNames(ArrayList<Adapter> adapters) {
    String[] names = new String[adapters.size()];
    for (int i = 0; i < names.length; i++) {
      names[i] = adapters.get(i).getName();
    }
    return names;

  }

  public static Adapter[] getBarcodes(String[] barCodes) {
    Adapter[] adapters = new Adapter[barCodes.length];
    for (int i = 0; i < barCodes.length; i++) {
      adapters[i] = new Adapter(barCodes[i], BARCODE, BARCODE);
    }
    return adapters;
  }

  public static ArrayList<Adapter> getCurrentAdapters(Adapter[] extras) {
    ArrayList<Adapter> adapters = new ArrayList<Adapter>();
    adapters.add(getAgilentQXTAdapter());
    adapters.add(getIlluminaTopAdapter());
    adapters.add(getIlluminaBottomAdapter());
    for (Adapter extra : extras) {
      adapters.add(extra);
    }
    return adapters;
  }

  public static Adapter getIlluminaBottomAdapter() {
    String seq = "GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG";
    String name = "IlluminaBottom_Adapter";
    String description = "Adapter sequence for Illumina Bottom adapter";
    return new Adapter(seq, name, description);
  }

  public static Adapter getIlluminaTopAdapter() {
    String seq = "ACACTCTTTCCCTACACGACGCTCTTCCGATC";
    String name = "IlluminaTop_Adapter";
    String description = "Adapter sequence for Illumina Top adapter";
    return new Adapter(seq, name, description);
  }

  private final String adaptorSequence;

  private final String adaptorRSSequence;// reverse strand

  private final String name;

  private final String description;

  public Adapter(String adaptorSequence, String name, String description) {
    super();
    this.adaptorSequence = adaptorSequence;
    adaptorRSSequence = StrandOps.flipsIfNeeded(adaptorSequence, Strand.NEGATIVE, false);
    this.name = name;
    this.description = description;
  }

  public String getAdaptorRSSequence() {
    return adaptorRSSequence;
  }

  public String getAdaptorSequence() {
    return adaptorSequence;
  }

  public String getDescription() {
    return description;
  }

  public String getName() {
    return name;
  }

}
