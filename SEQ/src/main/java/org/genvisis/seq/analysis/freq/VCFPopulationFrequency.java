/**
 * 
 */
package org.genvisis.seq.analysis.freq;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.StringJoiner;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.manage.VCFOps.VcfPopulation;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.WorkerHive;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.CLI;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;

/**
 * Class to compute population specific allele frequencies from a vcf
 */
public class VCFPopulationFrequency {

  /**
  * 
  */
  private static final String REGION = "region";
  /**
   * 
   */
  private static final String POP_FILE = "popFile";
  /**
   * 
   */
  private static final String VCF_DIR = "vcfDir";

  private static String getPopFreqKey(String population) {
    return "AF_" + population;
  }

  private static String getPopNKey(String population) {
    return "N_Alleles_" + population;
  }

  /**
   * Generate {@link VCFInfoHeaderLine}s of "AF_" and "N_Alleles_" for each population in the
   * {@link VcfPopulation}
   */
  private static List<VCFInfoHeaderLine> getNewInfoLines(VcfPopulation vpop) {
    List<VCFInfoHeaderLine> vcfInfos = new ArrayList<>();
    for (String population : vpop.getSubPop().keySet()) {
      VCFInfoHeaderLine infoF = new VCFInfoHeaderLine(getPopFreqKey(population), 1,
                                                      VCFHeaderLineType.String,
                                                      "Alternate allele frequency of population "
                                                                                + population
                                                                                + ". If multiple alternate alleles exist,"
                                                                                + " this will be a comma delimited list in the same order of the alternate alleles in the .vcf");
      vcfInfos.add(infoF);
      VCFInfoHeaderLine infoN = new VCFInfoHeaderLine(getPopNKey(population), 1,
                                                      VCFHeaderLineType.String,
                                                      "Number of total alleles for population "
                                                                                + population
                                                                                + " (total number of population samples is  "
                                                                                + vpop.getSubPop()
                                                                                      .get(population)
                                                                                      .size()
                                                                                + ")");
      vcfInfos.add(infoN);

    }

    return vcfInfos;
  }

  /**
   * Add the {@link VCFHeader} to the {@link VariantContextWriter}
   * 
   * @param vpop two new {@link VCFInfoHeaderLine}s will be added for each population (from
   *          {@link VCFPopulationFrequency#getNewInfoLines(VcfPopulation)}
   * @param header the base header to use
   * @param writer {@link VariantContextWriter} to write header to
   */
  private static void addHeader(VcfPopulation vpop, VCFHeader header, VariantContextWriter writer) {
    List<VCFInfoHeaderLine> vcfInfos = getNewInfoLines(vpop);
    for (VCFInfoHeaderLine vcfInfo : vcfInfos) {
      header.addMetaDataLine(vcfInfo);
    }
    writer.writeHeader(new VCFHeader(header.getMetaDataInInputOrder(), new HashSet<String>()));
  }

  /**
   * Initialize {@link VariantContextWriter} to site only mode
   * 
   * @param reader will use this reader's {@link SAMSequenceDictionary}
   * @param outputVcf desired output file
   * @return {@link VariantContextWriter}
   */
  private static VariantContextWriter buildVcfWriter(VCFFileReader reader, String outputVcf) {
    VariantContextWriterBuilder builder = new VariantContextWriterBuilder().setOutputFile(outputVcf);
    if (reader.getFileHeader().getSequenceDictionary() != null) {
      builder.setReferenceDictionary(reader.getFileHeader().getSequenceDictionary());
    }
    builder.setOption(Options.DO_NOT_WRITE_GENOTYPES);

    return builder.build();
  }

  /**
   * ensure all samples in the vpop are present in the reader's header
   * 
   * @param vpop {@link VcfPopulation}
   * @param reader {@link VCFFileReader}
   * @param log
   * @return true if error
   */
  private static boolean validateSamples(VcfPopulation vpop, VCFFileReader reader, Logger log) {
    ArrayList<String> samples = reader.getFileHeader().getSampleNamesInOrder();
    HashSet<String> popInds = vpop.getAllPopInds();
    boolean error = false;
    for (String sample : popInds) {
      if (!samples.contains(sample)) {
        error = true;
        log.reportError("Sample " + sample + "was not found in the .vcf");

      }
    }
    return error;
  }

  /**
   * compute population specific allele frequencies for a vcf
   * 
   * @param vcf input vcf file to compute population specific frequencies on
   * @param popFile population definition file
   * @param outDir output directory
   * @param region if not null, only this region will be scanned
   */
  private static boolean run(String vcf, String popFile, String outDir, String region) {

    //    validate input files exist
    new File(outDir).mkdirs();
    Logger log = new Logger(outDir + "popFreq.log");
    if (!Files.exists(popFile)) {
      log.reportError(popFile + " does not exist");
      return false;
    }
    if (!Files.exists(vcf)) {
      log.reportError(vcf + " does not exist");
      return false;
    }

    //    Load populations of interest
    VcfPopulation vpop = VcfPopulation.load(popFile, POPULATION_TYPE.AF, log);
    vpop.report();

    //    initialize reader and validate samples
    VCFFileReader reader = new VCFFileReader(new File(vcf), false);

    boolean error = validateSamples(vpop, reader, log);
    if (error) {
      reader.close();
      return false;
    }

    //    Initialize writer and add updated header 
    String outputVcf = outDir + VCFOps.getAppropriateRoot(vcf, true) + "_"
                       + ext.rootOf(popFile, true)
                       + (region == null ? ""
                                         : "_region_"
                                           + ext.replaceWithLinuxSafeCharacters(region, true))
                       + "_AF.vcf.gz";

    if (!(Files.exists(outputVcf) && Files.exists(outputVcf + ".tbi"))) {
      log.reportTimeInfo("Creating site only .vcf at " + outputVcf);
      VariantContextWriter writer = buildVcfWriter(reader, outputVcf);
      addHeader(vpop, reader.getFileHeader(), writer);

      //    Compute population specific AFs for all variants
      int numVariantsScanned = 0;

      CloseableIterator<VariantContext> iter;
      if (region == null) {
        iter = reader.iterator();
      } else {
        String[] split = region.split(":");
        String chr = split[0];
        int start = Integer.parseInt(split[1].split("-")[0]);
        int stop = Integer.parseInt(split[1].split("-")[1]);
        iter = reader.query(chr, start, stop);
        log.reportTime("restricting pop AF to " + region);
      }

      while (iter.hasNext()) {
        VariantContext vc = iter.next();
        VariantContextBuilder vcBuilder = new VariantContextBuilder(vc);
        List<Allele> alleles = vc.getAlternateAlleles();

        for (String population : vpop.getSubPop().keySet()) {
          //        number of total alleles
          int nTotal = vc.getCalledChrCount(vpop.getSubPop().get(population));
          StringJoiner afJoiner = new StringJoiner(",");

          for (Allele a : alleles) {
            //  number of this allele
            int nAllele = vc.getCalledChrCount(a, vpop.getSubPop().get(population));

            double afA = 0;
            if (nTotal > 0) {
              afA = (double) nAllele / nTotal;
            }
            afJoiner.add(Double.toString(afA));

          }
          vcBuilder.attribute(getPopFreqKey(population), afJoiner.toString());
          vcBuilder.attribute(getPopNKey(population), Integer.toString(nTotal));
        }
        writer.add(vcBuilder.make());
        numVariantsScanned++;
        if (numVariantsScanned % 10000 == 0) {
          log.reportTimeInfo("Computed population specific AF for " + numVariantsScanned
                             + " variants");
        }
      }
      reader.close();
      writer.close();
    } else {
      log.reportTimeWarning("Output " + outputVcf + " already exists, skipping");
    }
    return true;
  }

  public static void runDir(String vcfDir, String popFile, String outDir, String region,
                            int threads) {
    new File(outDir).mkdirs();
    Logger log = new Logger(outDir + "popFreq.log");

    String[] vcfs = Files.listFullPaths(vcfDir, VCFOps.VCF_EXTENSIONS.GZIP_VCF.getLiteral());
    log.reportTimeInfo("found " + vcfs.length + " " + VCFOps.VCF_EXTENSIONS.GZIP_VCF.getLiteral()
                       + "files in " + vcfDir);

    WorkerHive<Boolean> hive = new WorkerHive<>(threads, 10, log);
    for (final String vcf : vcfs) {
      hive.addCallable(() -> run(vcf, popFile, outDir, region));
    }
    hive.execute(true);

  }

  public static void main(String[] args) {
    CLI c = new CLI(VCFPopulationFrequency.class);
    c.addArg(CLI.ARG_VCF, CLI.DESC_VCF, false);
    c.addArg(POP_FILE,
             "File containing a header with \"IID\" and \"Population\", every sample in this file must be"
                       + " in the vcf. Allele frequency will be computed for each population",
             true);
    c.addArg(CLI.ARG_OUTDIR, CLI.DESC_OUTDIR, false);
    c.addArg(VCF_DIR,
             "a directory containing .vcf.gz files, if flagged, the threads option will be valid and "
                      + CLI.ARG_VCF + " will not be used",
             false);
    c.addArg(CLI.ARG_THREADS, "number of threads, only if " + VCF_DIR + " option is flagged",
             false);

    c.addArg(REGION, "only compute population specific allele frequencies in this UCSC region ",
             false);

    c.parseWithExit(args);

    if (c.has(VCF_DIR)) {
      runDir(c.get(VCF_DIR), c.get(POP_FILE), c.get(CLI.ARG_OUTDIR), c.get(REGION),
             c.getI(CLI.ARG_THREADS));
    } else {
      run(c.get(CLI.ARG_VCF), c.get(POP_FILE), c.get(CLI.ARG_OUTDIR), c.get(REGION));
    }

  }

}
