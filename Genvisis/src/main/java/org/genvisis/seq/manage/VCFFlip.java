package org.genvisis.seq.manage;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;

/**
 * @author lane0212 Flip a vcf's calls from alt to ref
 */
public class VCFFlip {
  // VCFOps.copyHeader(reader, writer, null, HEADER_COPY_TYPE.FULL_COPY, log);

  private static VCFFileReader reader;

  public static void flip(String inputVCF, boolean justGenoFlip, Logger log) {
    reader = new VCFFileReader(new File(inputVCF), true);
    String output = inputVCF;
    if (inputVCF.endsWith(".vcf.gz")) {
      output = output.replaceAll(".vcf.gz", ".flipped.vcf.gz");
    } else {
      output = ext.addToRoot(inputVCF, ".flipped");
    }
    VariantContextWriter writer = VCFOps.initWriter(output, null,
                                                    reader.getFileHeader().getSequenceDictionary());

    VariantContextWriterBuilder builderError = new VariantContextWriterBuilder().setOutputFile(inputVCF
                                                                                               + ".errors.vcf.gz");
    builderError.setBuffer(2);
    builderError.setReferenceDictionary(reader.getFileHeader().getSequenceDictionary());
    VariantContextWriter writerError = builderError.build();

    VCFHeader vcfHeader = new VCFHeader(reader.getFileHeader());
    VCFFormatHeaderLine vcfFormatHeaderLine = new VCFFormatHeaderLine("GT", 1,
                                                                      VCFHeaderLineType.String,
                                                                      "Ref to alt flipped genotypes from "
                                                                                                + inputVCF);
    vcfHeader.addMetaDataLine(vcfFormatHeaderLine);
    writer.writeHeader(vcfHeader);
    writerError.writeHeader(vcfHeader);

    int skipped = 0;
    int count = 0;
    for (VariantContext vc : reader) {

      if (vc.isBiallelic()) {
        count++;
        VariantContextBuilder builder = new VariantContextBuilder(vc);

        GenotypesContext genotypesContext = vc.getGenotypes();
        GenotypesContext newGenotypes = GenotypesContext.create();

        if (!justGenoFlip) {
          List<Allele> alleles = vc.getAlleles();
          List<Allele> newAlleles = new ArrayList<Allele>();
          Allele newRef = null;
          Allele newAlt = null;

          for (Allele a : alleles) {
            if (a.isReference()) {
              if (newAlt != null) {
                log.reportError("Multiple new alts....");
                return;
              }
              newAlt = Allele.create(a.getBases(), false);
              newAlleles.add(Allele.create(a.getBases(), false));
            } else {
              if (newRef != null) {
                log.reportError("Multiple new refs...." + newRef.getBaseString());
              }
              newRef = Allele.create(a.getBases(), true);

              newAlleles.add(Allele.create(a.getBases(), true));
            }
          }
          if (newRef.getBases().length != newAlt.getBases().length) {
            builder.start(vc.getStart());
            builder.stop(vc.getStart() + newRef.getBases().length - 1);

          }
          builder.alleles(newAlleles);

          ArrayList<Allele> newHomRef = new ArrayList<Allele>();
          newHomRef.add(newRef);
          newHomRef.add(newRef);
          ArrayList<Allele> newHomAlt = new ArrayList<Allele>();
          newHomAlt.add(newAlt);
          newHomAlt.add(newAlt);

          ArrayList<Allele> newHet = new ArrayList<Allele>();
          newHet.add(newRef);
          newHet.add(newAlt);

          for (Genotype g : genotypesContext) {
            GenotypeBuilder gbuilder = new GenotypeBuilder(g);
            if (g.isHomRef()) {
              gbuilder.alleles(newHomAlt);
            } else if (g.isHomVar()) {
              gbuilder.alleles(newHomRef);
            } else if (g.isHet()) {
              gbuilder.alleles(newHet);
            }
            newGenotypes.add(gbuilder.make());
          }
        } else {
          builder.alleles(vc.getAlleles());
          ArrayList<Allele> newHomRef = new ArrayList<Allele>();
          newHomRef.add(vc.getReference());
          newHomRef.add(vc.getReference());
          ArrayList<Allele> newHomAlt = new ArrayList<Allele>();
          newHomAlt.add(vc.getAlternateAlleles().get(0));
          newHomAlt.add(vc.getAlternateAlleles().get(0));

          ArrayList<Allele> newHet = new ArrayList<Allele>();
          newHet.add(vc.getReference());
          newHet.add(vc.getAlternateAlleles().get(0));

          for (Genotype g : genotypesContext) {
            GenotypeBuilder gbuilder = new GenotypeBuilder(g);
            if (g.isHomRef()) {
              gbuilder.alleles(newHomAlt);
            } else if (g.isHomVar()) {
              gbuilder.alleles(newHomRef);
            } else if (g.isHet()) {
              gbuilder.alleles(newHet);
            }
            newGenotypes.add(gbuilder.make());
          }
        }

        builder.genotypes(newGenotypes);
        try {
          if (vc.isBiallelic()) {
            VariantContext flippedVC = builder.make();
            writer.add(flippedVC);
          } else {
            writerError.add(vc);
          }

        } catch (IllegalArgumentException ile) {
          log.reportException(ile);
          skipped++;
          log.reportError("Could not flip variant context " + vc.toStringWithoutGenotypes()
                          + ", reverting to original context...");
          writerError.add(vc);
          // writer.add(vc);
        }
        if (count % 10000 == 0) {
          log.reportTimeInfo("Flipped " + count + " variants");
        }

      } else {
        writerError.add(vc);
      }
    }
    writer.close();
    writerError.close();
    if (skipped > 0) {
      log.reportTimeWarning("Some variants were invalid, we skipped " + skipped);
    }
    reader.close();
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String vcf = "VCFFlip.dat";
    String logfile = null;
    Logger log;

    String usage = "\n" + "seq.manage.VCFFlip requires 0-1 arguments\n"
                   + "   (1) vcf file (i.e. vcf=" + vcf + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("vcf=")) {
        vcf = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("log=")) {
        logfile = arg.split("=")[1];
        numArgs--;
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      log = new Logger(logfile);
      flip(vcf, true, log);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

}
