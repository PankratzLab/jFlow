package org.genvisis.seq.analysis;

import org.genvisis.seq.analysis.GATK_Genotyper.ANNOTATION_BUILD;
import org.genvisis.seq.manage.VCFOps;
import org.pankratzlab.common.CmdLine;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.PSF;

/**
 * Class for automating the annotation of a vcf using ANNOVAR
 */
public class ANNOVAR {

  public static final String ANNOVAR_COMMAND = "annovar=";
  public static final String TABLE_ANNOVAR = "table_annovar.pl";
  private static final String PROTOCOL = "-protocol";
  private static final String DEFAULT_HG19_PROTOCOLS = "refGene,cytoBand,genomicSuperDups,esp6500si_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,snp138,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,popfreq_max_20150413,popfreq_all_20150413,esp6500siv2_all,esp6500siv2_aa,esp6500siv2_ea,cosmic70,dbnsfp30a,mitimpact24";
  private static final String DEFAULT_HG19_OPERATIONS = "g,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f";
  private static final String DEFAULT_HG38_PROTOCOLS = "refGene,knownGene,ensGene,cytoBand,genomicSuperDups,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,avsnp150,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_eas,1000g2015aug_eur,esp6500siv2_all,esp6500siv2_aa,esp6500siv2_ea,cosmic70,dbnsfp33a,clinvar_20180603,exac03nontcga,intervar_20180118";
  private static final String DEFAULT_HG38_OPERATIONS = "g,g,g,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f";

  private static final String OPERATION = "-operation";
  private static final String REMOVE = "-remove";
  private static final String DEFUALT_ANNOVAR_DB = "humandb/";
  private static final String BUILD_VERSION = "-buildver";
  private static final String OUT = "-out";
  private static final String NA_STRING = "-nastring";
  private static final String DEFAULT_NA_STRING = ".";
  private static final String VCF_INPUT = "-vcfinput";
  private static final String MULTI_ANNO = "multianno.vcf";
  private static final String THREAD = "--thread";

  private final String annovarLocation;
  private String annodvarDB;
  private final boolean fail, verbose, overWriteExistingOutput;
  private final Logger log;

  public ANNOVAR(String annovarLocation, boolean verbose, boolean overWriteExistingOutput,
                 Logger log) {
    super();
    this.annovarLocation = annovarLocation;
    this.verbose = verbose;
    this.overWriteExistingOutput = overWriteExistingOutput;
    this.log = log;
    fail = !verify();
  }

  private boolean verify() {
    boolean verify = true;
    if (annodvarDB == null) {
      // we don't care for now
    }
    if (!Files.exists(annovarLocation)) {
      verify = false;
      log.reportError("Warning - could not find Annovar directory  " + annovarLocation);
    }
    if (!Files.exists(annovarLocation + TABLE_ANNOVAR)) {
      verify = false;
      log.reportError("Warning - could not find " + annovarLocation + TABLE_ANNOVAR);
    }
    if (!Files.exists(annovarLocation + DEFUALT_ANNOVAR_DB)) {
      verify = false;
      log.reportError("Warning - could not find " + annovarLocation + DEFUALT_ANNOVAR_DB);
    }
    return verify;
  }

  public String getAnnovarLocation() {
    return annovarLocation;
  }

  public boolean isFail() {
    return fail;
  }

  public AnnovarResults AnnovarAVCF(String inputVCF, ANNOTATION_BUILD build, int numthreads,
                                    Logger log) {
    AnnovarResults annovarResults = new AnnovarResults(inputVCF, build);
    annovarResults.parse();
    boolean progress = AnnovarAVCF(annovarResults.getInputVCF(), annovarResults.getOutput(),
                                   annovarResults.getOutputVCF(), annovarResults.getBuild(),
                                   numthreads, log);
    annovarResults.setFail(!progress);
    if (!progress) {
      log.reportError("The annovar has failed :( \n This is most likely caused by not having the required database files. Please see \"http://www.openbioinformatics.org/annovar/annovar_startup.html\" for commands to download the necessary files");
    }
    return annovarResults;
  }

  private boolean AnnovarAVCF(String inputVCF, String outputBase, String outputVCF,
                              ANNOTATION_BUILD build, int numthreads, Logger log) {
    boolean progress = !fail;
    if (progress) {
      // THREAD, numthreads + ""

      String protocols;
      String operations;
      switch (build) {
        case HG19:
          protocols = DEFAULT_HG19_PROTOCOLS;
          operations = DEFAULT_HG19_OPERATIONS;
          break;
        case HG38:
          protocols = DEFAULT_HG38_PROTOCOLS;
          operations = DEFAULT_HG38_OPERATIONS;
          break;
        default:
          throw new IllegalArgumentException("Invalid build for ANNOVAR, " + build);

      }

      String[] command = new String[] {PSF.Cmd.PERL, annovarLocation + TABLE_ANNOVAR, inputVCF,
                                       annovarLocation + DEFUALT_ANNOVAR_DB, BUILD_VERSION,
                                       build.getAnnovarBuild(), OUT, outputBase, REMOVE, PROTOCOL,
                                       protocols, OPERATION, operations, NA_STRING,
                                       DEFAULT_NA_STRING, VCF_INPUT, THREAD,
                                       Integer.toString(numthreads)};
      progress = CmdLine.runCommandWithFileChecks(command, "", new String[] {inputVCF},
                                                  new String[] {outputVCF}, verbose,
                                                  overWriteExistingOutput, false, log);
    }
    return progress;
  }

  public static class AnnovarResults {

    private final String inputVCF;
    private final ANNOTATION_BUILD build;
    private String output;
    private String outputVCF;
    private boolean fail;
    private Logger log;

    public AnnovarResults(String inputVCF, ANNOTATION_BUILD build) {
      super();
      this.inputVCF = inputVCF;
      this.build = build;
      fail = false;
    }

    public void parse() {
      output = VCFOps.getAppropriateRoot(inputVCF, false);
      outputVCF = output + "." + build.getAnnovarBuild() + "_" + MULTI_ANNO;
    }

    public String getInputVCF() {
      return inputVCF;
    }

    public ANNOTATION_BUILD getBuild() {
      return build;
    }

    public boolean isFail() {
      return fail;
    }

    public void setFail(boolean fail) {
      this.fail = fail;
    }

    public String getOutput() {
      return output;
    }

    public String getOutputVCF() {
      return outputVCF;
    }

    public Logger getLog() {
      return log;
    }

  }

}
