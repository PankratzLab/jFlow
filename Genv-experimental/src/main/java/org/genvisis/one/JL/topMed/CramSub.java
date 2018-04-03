package org.genvisis.one.JL.topMed;

import java.io.File;
import java.io.IOException;
import org.genvisis.common.Files;
import org.genvisis.common.ext;
import org.genvisis.seq.manage.SamRecordOps;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

/**
 * make smaller .crams for testing
 */
public class CramSub {

  public static void main(String[] args) {
    String dir = args[1];
    String ref = args[2];
    String out = args[3];

    System.setProperty("samjdk.reference_fasta", ref);

    new File(out).mkdirs();
    String[] crams = Files.listFullPaths(dir, ".cram");
    System.out.println("found " + crams.length + " crams");
    for (String cram : crams) {
      System.out.println(cram);
      SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault();
      samReaderFactory.validationStringency(ValidationStringency.STRICT);

      SamReader reader = samReaderFactory.open(new File(cram));

      String outFile = out + ext.addToRoot(ext.removeDirectoryInfo(cram), ".sub");
      SAMFileWriter sAMFileWriter = new SAMFileWriterFactory().setCreateIndex(true)
                                                              .makeSAMOrBAMWriter(reader.getFileHeader(),
                                                                                  true,
                                                                                  new File(outFile));
      SAMRecordIterator sIterator = reader.query("chr1", 10000000, 20000000, true);

      int num = 0;
      while (sIterator.hasNext()) {
        num++;

        SAMRecord record = sIterator.next();
        if (num % 10000 == 0) {
          System.out.println(num + "\t" + SamRecordOps.getDisplayLoc(record));
        }
        sAMFileWriter.addAlignment(record);
      }
      System.out.println(num);
      try {
        reader.close();
      } catch (IOException e) {
        e.printStackTrace();
      }
      sAMFileWriter.close();
    }

  }

}
