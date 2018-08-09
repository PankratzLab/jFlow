package org.genvisis.seq.manage;

import java.io.File;
import java.io.PrintWriter;
import org.apache.commons.collections4.BidiMap;
import org.apache.commons.collections4.bidimap.DualHashBidiMap;
import org.genvisis.cnv.manage.Resources.GENOME_BUILD;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.common.ext;
import org.genvisis.filesys.LocusSet;
import org.genvisis.filesys.Segment;
import com.google.common.collect.Sets;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

public class NGSBinSNPSelector {

  private static final boolean DEBUG = true;

  LocusSet<Segment> bins;
  String snpVCF;

  public void run() {
    PrintWriter writerDEBUG = DEBUG ? Files.getAppropriateWriter("./All_Variants.xln") : null;
    PrintWriter writer = Files.getAppropriateWriter("./selectedVariants.xln");

    VCFHeader header;
    try (VCFFileReader reader = new VCFFileReader(new File(snpVCF),
                                                  Files.exists(snpVCF + ".tbi"))) {
      header = reader.getFileHeader();

      // discover which chrs are in this vcf file
      BidiMap<String, Integer> contigMap = new DualHashBidiMap<>();
      header.getContigLines().forEach(vch -> {
        // ensure parsability
        contigMap.put(vch.getID(), (int) Positions.chromosomeNumber(vch.getID()));
      });

      for (Segment bin : bins.getLoci()) {
        // no variants in file for this bin's chromosome  
        if (!contigMap.containsValue((int) bin.getChr())) {
          System.err.println("Error - bin contig " + bin.getChr() + " was not present in the VCF!");
          continue;
        }

        int mid = bin.getStart() + ((bin.getStop() - (bin.getStart() + 1)) / 2);
        //  example:      mid = 1 + ((1000 - (1 + 1)) / 2)
        //                    = 1 + 998 / 2
        //                    = 500
        //------------------
        // Dropping the '+1' results in 500.5, which relies on the properties of 
        // integer casting to drop to the true midpoint of 500. 
        CloseableIterator<VariantContext> iter = reader.query("" + bin.getChr(), bin.getStart(),
                                                              bin.getStop());
        VariantContext selected = null;
        double lowestScale = Double.MAX_VALUE;
        int count = 0;
        while (iter.hasNext()) {
          count++;
          VariantContext vc = iter.next();

          int dist = Math.abs(vc.getStart() - mid);
          String af = vc.getAttributes().get("AF").toString();
          if (af.contains(",")) {
            // multi-allelic
            continue;
          }
          double afD = Double.parseDouble(af);
          double afDist = Math.abs(afD - .5);
          double scale = Math.sqrt(dist * dist + afDist * afDist);

          String id = vc.getID();
          if (ext.isMissingValue(id)) {
            id = vc.getContig() + ":" + vc.getStart();
          }
          if (DEBUG) {
            writerDEBUG.println(bin.getUCSClocation() + "\t" + id + "\t" + dist + "\t" + afDist
                                + "\t" + scale);
          }

          if (scale < lowestScale) {
            selected = vc;
            lowestScale = scale;
          }
        }

        if (selected != null) {
          String id = selected.getID();
          if (ext.isMissingValue(id)) {
            id = selected.getContig() + ":" + selected.getStart();
          }
          int dist = Math.abs(selected.getStart() - mid);
          double afDist = Math.abs(selected.getAttributeAsDouble("AF", Double.NaN) - .5);
          double scale = Math.sqrt(dist * dist + afDist * afDist);

          writer.println(bin.getUCSClocation() + "\t" + id + "\t" + dist + "\t" + afDist + "\t"
                         + scale + "\t" + count);
        } else {
          writer.println(bin.getUCSClocation() + "\t.\t.\t.\t." + "\t" + count);
        }
      }

    }

    writer.close();
    if (DEBUG) {
      writerDEBUG.close();
    }

  }

  public static void main(String[] args) {
    NGSBinSNPSelector selector = new NGSBinSNPSelector();
    selector.snpVCF = "/home/spectorl/shared/Ewing_WGS/VQSR/ES_recalibrated_snps_indels.vcf.gz";
    selector.bins = new ReferenceGenome(GENOME_BUILD.HG19,
                                        new Logger()).getBins(1000,
                                                              Sets.newHashSet(1, 2, 3, 4, 5, 6, 7,
                                                                              8, 9, 10, 11, 12, 13,
                                                                              14, 15, 16, 17, 18,
                                                                              19, 20, 21, 22, 23,
                                                                              24, 25));
    selector.run();
  }

}
