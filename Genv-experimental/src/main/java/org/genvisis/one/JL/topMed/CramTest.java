package org.genvisis.one.JL.topMed;

import java.io.File;
import java.io.IOException;

import org.genvisis.seq.manage.SamRecordOps;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

/**
 * Test of .cram reading/writing etc
 */
public class CramTest {

  public static void main(String[] args) {
    String cram = "my.cram";
    String ref = "hs38DH.fa";
    String out = "testOut.cram";
    // Setting system property is nice, allows anywhere to have access to the reference fasta
    System.setProperty("samjdk.reference_fasta", ref);

    // SamReader reader = BamOps.getDefaultReader(cram, ValidationStringency.LENIENT);
    SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault();
    samReaderFactory.validationStringency(ValidationStringency.STRICT);
    // samReaderFactory.referenceSequence(new File(ref));
    // for (Option option : options) {
    // samReaderFactory.enable(option);
    // }
    SamReader reader = samReaderFactory.open(new File(cram));

    SAMFileWriter sAMFileWriter = new SAMFileWriterFactory().setCreateIndex(true)
                                                            .makeSAMOrBAMWriter(reader.getFileHeader(),
                                                                                true,
                                                                                new File(out));
    // SAMRecordIterator sIterator = reader.iterator();
    SAMRecordIterator sIterator = reader.query("chr1", 1, 10000000, true);

    int num = 0;
    while (sIterator.hasNext()) {
      num++;

      SAMRecord record = sIterator.next();
      if (num % 10000 == 0) {
        System.out.println(num + "\t" + SamRecordOps.getDisplayLoc(record));
      }
      // System.out.println(record.toString());
      sAMFileWriter.addAlignment(record);
      // System.out.println(SamRecordOps.getDisplayLoc(record));
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
