package org.genvisis.one.JL.topMed;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Date;
import java.util.List;
import java.util.StringJoiner;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.seq.manage.BamOps;
import com.google.common.math.Stats;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class CramStats {

  public static void main(String[] args) {
    String dir = args[1];
    String ref = args[2];
    String out = args[3];

    System.setProperty("samjdk.reference_fasta", ref);

    new File(out).mkdirs();
    String[] crams = Files.listFullPaths(dir, ".cram");
    System.out.println("found " + crams.length + " crams");

    String outFile = out + "summary.txt";
    PrintWriter writer = Files.getAppropriateWriter(outFile);

    writer.println("Sample\tSEQUENCING_CENTER\tDATE_RUN_PRODUCED\tPLATFORM\tLIBRARY\tDESCRIPTION\tAVERAGE_READ_LENGTH\tMAX_READ_LENGTH");

    for (String cram : crams) {

      SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault();
      samReaderFactory.validationStringency(ValidationStringency.STRICT);

      SamReader reader = samReaderFactory.open(new File(cram));

      SAMFileHeader header = reader.getFileHeader();
      Stats rl = BamOps.estimateReadSize(cram, new Logger());
      List<SAMReadGroupRecord> rgs = header.getReadGroups();
      for (SAMReadGroupRecord rg : rgs) {
        StringJoiner joiner = new StringJoiner("\t");
        joiner.add(rg.getSample());
        joiner.add(rg.getSequencingCenter());
        joiner.add(rg.getRunDate() == null ? new Date().toString() : rg.getRunDate().toString());
        joiner.add(rg.getPlatform());
        joiner.add(rg.getLibrary());
        joiner.add(rg.getDescription());
        joiner.add(Double.toString(rl.mean()));
        joiner.add(Double.toString(rl.max()));

        writer.println(joiner.toString());

      }

      try {
        reader.close();
      } catch (IOException e) {
        e.printStackTrace();
      }
    }
    writer.close();

  }
}
