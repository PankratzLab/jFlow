package org.genvisis.seq.analysis;

import java.util.ArrayList;

import org.genvisis.common.Array;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.seq.manage.VCFOps;

/**
 * @author lane0212 Class for annotating a vcf using a local version of oncotator
 *         -http://www.broadinstitute.org/oncotator/
 */
public class Oncotator {
  public static class Oncotator_Analysis {
    private static final String ONCO_TAG = ".onco";
    private final String inputVcf;
    private final String outputVcf;
    private boolean fail;

    public Oncotator_Analysis(String inputVcf) {
      super();
      this.inputVcf = inputVcf;
      outputVcf = VCFOps.getAppropriateRoot(inputVcf, false) + ONCO_TAG + ".vcf";
      fail = false;
    }

    public String getInputVcf() {
      return inputVcf;
    }

    public String getOutputVcf() {
      return outputVcf;
    }

    public boolean isFail() {
      return fail;
    }

    public void setFail(boolean fail) {
      this.fail = fail;
    }

  }

  private static final String ONCOTATOR = "oncotator";
  private static final String DB_DIR = "--db-dir";
  private static final String I = "-i";

  private static final String O = "-o";

  private static final String VCF = "VCF";

  private static boolean verify(String oncoDBLoc, Logger log) {
    boolean ver = CmdLine.run(ONCOTATOR, "");
    if (!ver) {
      log.reportTimeError(ONCOTATOR + " was not found on the system path, cannot run ");
    } else {
      ver = Files.exists(oncoDBLoc);
      if (!ver) {
        log.reportTimeError(
            ONCOTATOR + " db directory " + oncoDBLoc + " did not exist, cannot run ");
      }
    }
    return ver;
  }

  private final String oncoDBLoc;
  private final String genomeBuild;
  private final Logger log;

  private final boolean fail;

  private final boolean overWriteExistingOutput;

  /**
   * @param oncoDBLoc this directory is required to exist
   * @param genomeBuild hg19, hg18, etc
   * @param overWriteExistingOutput if true, existing output will be over written
   * @param log
   */
  public Oncotator(String oncoDBLoc, String genomeBuild, boolean overWriteExistingOutput,
      Logger log) {
    super();
    this.oncoDBLoc = oncoDBLoc;
    this.genomeBuild = genomeBuild;
    this.overWriteExistingOutput = overWriteExistingOutput;
    this.log = log;
    fail = !verify(oncoDBLoc, log);
  }

  private boolean oncoAnnoate(Oncotator_Analysis oncotator_Analysis) {
    boolean anno = true;
    String[] inputFiles = new String[] {oncotator_Analysis.getInputVcf()};
    String[] outputFiles = new String[] {oncotator_Analysis.getOutputVcf()};

    ArrayList<String> command = new ArrayList<String>();
    command.add(ONCOTATOR);
    command.add(I);
    command.add(VCF);
    command.add(DB_DIR);
    command.add(oncoDBLoc);
    command.add(O);
    command.add(VCF);
    command.add(oncotator_Analysis.getInputVcf());
    command.add(oncotator_Analysis.getOutputVcf());
    command.add(genomeBuild);
    anno = CmdLine.runCommandWithFileChecks(Array.toStringArray(command), "", inputFiles,
        outputFiles, true, overWriteExistingOutput, false, log);

    return anno;
  }

  public Oncotator_Analysis oncoAnnotate(String vcf) {
    Oncotator_Analysis oncotator_Analysis = new Oncotator_Analysis(vcf);
    if (!fail) {
      boolean annot = oncoAnnoate(oncotator_Analysis);
      oncotator_Analysis.setFail(!annot);
    } else {
      oncotator_Analysis.setFail(true);
    }
    return oncotator_Analysis;
  }

}
