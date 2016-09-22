/**
 * 
 */
package org.genvisis.seq.manage.mtdna;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.seq.manage.VCFOps;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;

/**
 * Class for managing the GenBank vcfs
 *
 */
public class GenBankMtDNA {


  /**
   * Trim whitespace from genbank vcf
   * 
   * @param input
   * @param output
   * @param log
   */
  public static void fixGenBankVCF(String input, String output, Logger log) {
    try {
      BufferedReader reader = Files.getAppropriateReader(input);
      PrintWriter writer = Files.getAppropriateWriter(output + ".tmp");
      while (reader.ready()) {
        String line = reader.readLine();
        if (line.startsWith("#")) {
          writer.println(line);
        } else {
          String[] data = line.trim().split("\t");
          // Info is index 7
          data[7] = data[7].replaceAll(" ", "");
          writer.println(Array.toStr(data));
        }
      }
      reader.close();
      writer.close();

    } catch (FileNotFoundException e) {
      log.reportException(e);

    } catch (IOException e) {
      log.reportException(e);

    }
    VCFFileReader reader = new VCFFileReader(new File(output + ".tmp"), false);
    new File(ext.parseDirectoryOfFile(output)).mkdirs();


    VariantContextWriter writer = VCFOps.initWriter(output, VCFOps.DEFUALT_WRITER_OPTIONS,
                                                    reader.getFileHeader().getSequenceDictionary());

    writer.writeHeader(reader.getFileHeader());
    CloseableIterator<VariantContext> vcIter = reader.iterator();
    while (vcIter.hasNext()) {
      try {
        VariantContextBuilder builder = new VariantContextBuilder(vcIter.next());
        writer.add(builder.make());
      } catch (Exception e) {
        log.reportTimeError("could not parse vc ,skipping");
        log.reportException(e);

      }
    }
    reader.close();
    writer.close();

  }

}
