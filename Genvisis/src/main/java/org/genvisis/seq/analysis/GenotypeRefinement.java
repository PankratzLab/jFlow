package org.genvisis.seq.analysis;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Set;
import java.util.concurrent.Callable;
import org.genvisis.CLI;
import org.genvisis.cnv.manage.Resources.GENOME_BUILD;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.Positions;
import org.genvisis.common.WorkerTrain;
import org.genvisis.common.WorkerTrain.Producer;
import org.genvisis.common.ext;
import org.genvisis.filesys.Pedfile;
import org.genvisis.seq.analysis.GATK.GenotypeRefiner;
import org.genvisis.seq.analysis.GATK.Mutect;
import org.genvisis.seq.manage.ReferenceGenome;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.manage.VCFOps.ChrSplitResults;
import org.genvisis.seq.manage.VCFOps.VcfPopulation;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import org.genvisis.seq.manage.VCOps;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

/**
 * Methods to wrap the GATK genotype refinement pipeline
 * https://gatkforums.broadinstitute.org/gatk/discussion/4723/genotype-refinement-workflow
 */
public class GenotypeRefinement {

  /**
   * 
   */
  private static final String LO_CONF_DE_NOVO = "loConfDeNovo";
  /**
   * 
   */
  private static final String HI_CONF_DE_NOVO = "hiConfDeNovo";

  private GenotypeRefinement() {

  }

  // https://software.broadinstitute.org/gatk/documentation/article.php?id=4726
  // https://software.broadinstitute.org/gatk/documentation/article.php?id=4727

  // Use 1000G_phase3_v4_20130502.sites.vcf as supporting VCF for genotype refinement
  // http://gatkforums.broadinstitute.org/gatk/discussion/5071/calculategenotypeposteriors-error

  private static void refineGenotypes(String vcf, String ped, String outputDir, GATK gatk,
                                      int threads) {
    new File(outputDir).mkdirs();
    Logger log = new Logger(outputDir + "gtRefine.log");
    String[] sampsinPed = HashVec.loadFileToStringArray(ped, false, Pedfile.COMMENT_INDICATORS,
                                                        new int[] {1}, true);
    String[] vcfSamps = VCFOps.getSamplesInFile(vcf);
    if (!ext.containsAll(sampsinPed, vcfSamps)) {
      throw new IllegalArgumentException("All pedigree samples must be in vcf file");
    }
    log.reportTimeInfo(vcfSamps.length + " samples in VCF");
    log.reportTimeInfo(sampsinPed.length + " samples in ped");

    String analysisVCF = vcf;
    if (vcfSamps.length > sampsinPed.length) {
      analysisVCF = subToPedSamples(vcf, ped, outputDir, log, sampsinPed);
    }

    ChrSplitResults[] splits = VCFOps.splitByChrs(analysisVCF, outputDir, 24, true, log);

    RefinementProducer producer = new RefinementProducer(gatk, splits, outputDir, ped, log);

    List<GenotypeRefiner> refinedResults = new ArrayList<>();
    try (WorkerTrain<GenotypeRefiner> train = new WorkerTrain<>(producer, threads, 10, log)) {
      while (train.hasNext()) {
        GenotypeRefiner tmp = train.next();
        if (!tmp.isFail()) {
          refinedResults.add(tmp);
        } else {
          log.reportError("Could not refine genotypes for " + tmp.toString());
        }
      }
    }
    log.reportTimeWarning("Starting hacky portion to get quick look at denovo variants");

    // combineResults(refinedResults, analysisVCF, outputDir, log);
    // Later, when the original vcf is already annotated, we can just simply combine results from
    // contig denovo vcfs using combineResults
    // Just don't want to annotate every variant right now
    List<String> filtDenovos = filterAndCombineDenovos(refinedResults, analysisVCF, outputDir, log);

    for (String filtDenovoVCF : filtDenovos) {
      GATK_Genotyper.annotateOnlyWithDefualtLocations(filtDenovoVCF, null,
                                                      PSF.Ext.DEFAULT_MEMORY_MB, true, false, log);
    }
  }

  /**
   * This method is untested, but is the basic replacement for
   * {@link GenotypeRefinement#filterAndCombineDenovos}
   * 
   * @param refinedResults this assumes that each {@link GenotypeRefiner} has been successfully run
   *          (has result vcf)
   * @param originalVCF for output naming
   * @param outputDir
   * @param log
   */
  private static void combineResults(List<GenotypeRefiner> refinedResults, String originalVCF,
                                     String outputDir, Logger log) {
    String refinedVCF = outputDir + VCFOps.getAppropriateRoot(originalVCF, true)
                        + ".refined.vcf.gz";
    if (!Files.exists(refinedVCF) || !Files.exists(VCFOps.getIndex(refinedVCF))) {

      VCFFileReader tmp = new VCFFileReader(new File(refinedResults.get(0).getDenovoVCF()));

      VCFHeader header = tmp.getFileHeader(); // header from results should contain new refinement
                                             // annotations

      SAMSequenceDictionary dict = header.getSequenceDictionary();
      VariantContextWriter outputWriter = VCFOps.initWriter(refinedVCF,
                                                            VCFOps.DEFUALT_WRITER_OPTIONS, dict);
      outputWriter.writeHeader(header);

      tmp.close();
      int total = 0;
      for (SAMSequenceRecord record : dict.getSequences()) { // should iterate in sorted
        // order
        for (GenotypeRefiner refiner : refinedResults) {
          String denovoVCF = refiner.getDenovoVCF();

          String contig = VCFOps.getFirstContig(denovoVCF, log);
          if (Positions.chromosomeNumber(contig) == Positions.chromosomeNumber(record.getSequenceName())) { // so
            // efficient
            log.reportTimeInfo("Matching " + ext.removeDirectoryInfo(denovoVCF) + " to contig "
                               + record.getSequenceName());
            VCFFileReader reader = new VCFFileReader(new File(denovoVCF), true);

            for (VariantContext vc : reader) {
              total++;
              if (!vc.getContig().equals(contig)) {
                reader.close();
                outputWriter.close();
                throw new IllegalArgumentException("Method designed to combine vcfs split by contigs");
              }
              outputWriter.add(vc);

              if (total % 1000000 == 0) {
                log.reportTimeInfo("Reading file " + ext.removeDirectoryInfo(denovoVCF)
                                   + ", written  " + total + " variants ");
              }
            }
            reader.close();
          }
        }
      }
      outputWriter.close();
    } else {
      log.reportFileExists(refinedVCF);
    }
  }

  /**
   * @param refinedResults
   * @param outputVCF
   * @param outputDir Extracts potential denovos from refined results, and converts contigs to hg19
   *          ...for now
   */
  private static List<String> filterAndCombineDenovos(List<GenotypeRefiner> refinedResults,
                                                      String originalVCF, String outputDir,
                                                      Logger log) {
    String hiConfVCF = outputDir + VCFOps.getAppropriateRoot(originalVCF, true)
                       + ".hiConfDenovo.vcf.gz";
    String hiLowConfVCF = outputDir + VCFOps.getAppropriateRoot(originalVCF, true)
                          + ".hiLowConfDenovo.vcf.gz";
    List<String> filtDenovos = new ArrayList<>();
    filtDenovos.add(hiConfVCF);
    filtDenovos.add(hiLowConfVCF);
    if (!Files.exists(hiConfVCF) && !Files.exists(hiLowConfVCF)) {
      SAMSequenceDictionary newDictionary = new ReferenceGenome(GENOME_BUILD.HG19,
                                                                log).getDictionary();
      VariantContextWriter hiLowConfwriter = VCFOps.initWriter(hiLowConfVCF,
                                                               VCFOps.DEFUALT_WRITER_OPTIONS,
                                                               newDictionary);
      VariantContextWriter hiConfwriter = VCFOps.initWriter(hiConfVCF,
                                                            VCFOps.DEFUALT_WRITER_OPTIONS,
                                                            newDictionary);
      int potentialDenovos = 0;
      int total = 0;
      boolean writeHeader = true;
      for (SAMSequenceRecord record : newDictionary.getSequences()) { // should iterate in sorted
                                                                     // order
        for (GenotypeRefiner refiner : refinedResults) {
          String denovoVCF = refiner.getDenovoVCF();

          String contig = VCFOps.getFirstContig(denovoVCF, log);
          if (Positions.chromosomeNumber(contig) == Positions.chromosomeNumber(record.getSequenceName())) { // so
                                                                                                           // efficient
            log.reportTimeInfo("Matching " + ext.removeDirectoryInfo(denovoVCF) + " to contig "
                               + record.getSequenceName());
            VCFFileReader reader = new VCFFileReader(new File(denovoVCF), true);
            if (writeHeader) {
              VCFHeader header = reader.getFileHeader();
              header.setSequenceDictionary(newDictionary);
              hiConfwriter.writeHeader(header);
              hiLowConfwriter.writeHeader(header);
            }
            writeHeader = false;
            for (VariantContext vc : reader) {
              total++;
              if (!vc.getContig().equals(contig)) {
                reader.close();
                hiConfwriter.close();
                hiLowConfwriter.close();
                throw new IllegalArgumentException("Method designed to combine contig split vcfs");
              } else if (!vc.isFiltered() && isPotentialDenovo(vc)) { // Lots of tranche garbage
                VariantContextBuilder builder = new VariantContextBuilder(vc);
                builder.chr(record.getSequenceName());
                VariantContext vcWrite = builder.make();
                hiLowConfwriter.add(vcWrite);
                if (isHighConfPotentialDenovo(vcWrite)) {
                  hiConfwriter.add(vcWrite);
                }
                potentialDenovos++;
              }
              if (total % 100000 == 0) {
                log.reportTimeInfo("Reading file " + ext.removeDirectoryInfo(denovoVCF) + ", found "
                                   + potentialDenovos + " potential denovos from " + total
                                   + " total variants scanned");
              }
            }
            reader.close();
          }
        }
      }
      hiConfwriter.close();
      hiLowConfwriter.close();
      log.reportTimeInfo("found " + potentialDenovos + " potential denovos in total from " + total
                         + " variants scanned");
    }
    return filtDenovos;
  }

  private static boolean isHighConfPotentialDenovo(VariantContext vc) {
    String[] confDenovos = VCOps.getAnnotationsFor(new String[] {HI_CONF_DE_NOVO}, vc, ".");
    return !".".equals(confDenovos[0]);
  }

  private static boolean isPotentialDenovo(VariantContext vc) {
    String[] confDenovos = VCOps.getAnnotationsFor(new String[] {HI_CONF_DE_NOVO, LO_CONF_DE_NOVO},
                                                   vc, ".");
    return !".".equals(confDenovos[0]) || !".".equals(confDenovos[1]);
  }

  private static class RefinementProducer implements Producer<GenotypeRefiner> {

    private ChrSplitResults[] splits;
    private GATK gatk;
    private String outputDir;
    private String ped;
    private Logger log;
    private int index;

    public RefinementProducer(GATK gatk, ChrSplitResults[] splits, String outputDir, String ped,
                              Logger log) {
      super();
      this.gatk = gatk;
      this.splits = splits;
      this.outputDir = outputDir;
      this.ped = ped;
      this.log = log;
      this.index = 0;

    }

    /*
     * (non-Javadoc)
     * @see java.util.Iterator#hasNext()
     */
    @Override
    public boolean hasNext() {
      return index < splits.length;
    }

    /*
     * (non-Javadoc)
     * @see java.util.Iterator#next()
     */
    @Override
    public Callable<GenotypeRefiner> next() {
      if (index >= splits.length) {
        throw new NoSuchElementException();
      }
      final ChrSplitResults current = splits[index];

      index++;
      return () -> gatk.refineGenotypes(current.getOutputVCF(), ped, outputDir, log);
    }

    /*
     * (non-Javadoc)
     * @see org.genvisis.common.WorkerTrain.Producer#shutdown()
     */
    @Override
    public void shutdown() {
      //
    }

  }

  /**
   * Subset the vcf to only samples in the ped file. I _think_ I was hitting GATK errors if a sample
   * was not represented
   * 
   * @param vcf
   * @param ped
   * @param outputDir
   * @param log
   * @param sampsinPed
   * @return
   */
  private static String subToPedSamples(String vcf, String ped, String outputDir, Logger log,
                                        String[] sampsinPed) {
    String analysisVCF;
    Set<String> sub = new HashSet<>();
    Map<String, Set<String>> pop = new HashMap<>();
    pop.put(ext.rootOf(ped), sub);
    for (String element : sampsinPed) {
      pop.get(ext.rootOf(ped)).add(element);
    }
    VcfPopulation vpop = new VcfPopulation(pop, pop, POPULATION_TYPE.ANY, log);
    String vpopTmp = outputDir + ext.rootOf(ped) + ".vpop";
    vpop.dump(vpopTmp);
    log.reportTimeWarning("Subsetting vcf to samples found in " + ped);
    analysisVCF = VcfPopulation.splitVcfByPopulation(vcf, vpopTmp, true, false, false, false,
                                                     log)[0];
    return analysisVCF;
  }

  /**
   * @param args
   */
  public static void main(String[] args) {
    CLI c = new CLI(GenotypeRefinement.class);
    c.addArg(CLI.ARG_OUTDIR, CLI.DESC_OUTDIR);
    c.addArg(CLI.ARG_REFERENCE_GENOME, CLI.DESC_REFERENCE_GENOME);
    c.addArg(CLI.ARG_VCF, CLI.DESC_VCF);
    c.addArg("gatk", "gatk location");
    c.addArg("supportSnps", "supporting snps for priors");
    c.addArg("ped", "pedigree file");
    c.addArg(CLI.ARG_THREADS, CLI.DESC_THREADS);

    c.parseWithExit(args);

    new File(c.get(CLI.ARG_OUTDIR)).mkdirs();
    Logger log = new Logger(c.get(CLI.ARG_OUTDIR) + "TN.log");
    GATK gatk = new Mutect(c.get("gatk"), c.get(CLI.ARG_REFERENCE_GENOME),
                           PSF.Ext.DEFAULT_MEMORY_MB, null, null, null, null, true, false, log);
    gatk.setSupportingSnps(c.get("supportSnps"));
    refineGenotypes(c.get(CLI.ARG_VCF), c.get("ped"), c.get(CLI.ARG_OUTDIR), gatk,
                    c.getI(CLI.ARG_THREADS));

  }
}
