package org.genvisis.cnv.workflow.steps;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.EnumSet;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.cnv.workflow.StepBuilder;
import org.genvisis.cnv.workflow.Variables;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.parsing.AliasedFileColumn;
import org.pankratzlab.common.parsing.DataLine;
import org.pankratzlab.common.parsing.FileColumn;
import org.pankratzlab.common.parsing.FileParser;
import org.pankratzlab.common.parsing.FileParserFactory;
import org.pankratzlab.common.parsing.ParseFailureException;
import org.pankratzlab.common.parsing.StandardFileColumns;

public class AxiomManifestParsingStep extends Step {

  public static final String AXIOM_MANIFEST_DESC = "An Affymetrix Axiom Manifest file";
  public static final String AXIOM_EXAMPLE_MANIFEST = "Axiom_tx_v1.na35.annot.csv";
  public static final String NAME = "Parse Axiom Manifest";
  public static final String DESC = "";
  public static final String RSLOOKUP_SUFFIX = ".rsLookup.txt";

  public static AxiomManifestParsingStep create(Project proj) {
    String defaultFile = Files.exists(proj.SNP_DATA_FILE.getValue()) ? proj.SNP_DATA_FILE.getValue()
                                                                     : AXIOM_EXAMPLE_MANIFEST;
    final Requirement<File> manifestReq = new Requirement.FileRequirement("manifestFile",
                                                                          AXIOM_MANIFEST_DESC,
                                                                          new File(defaultFile));
    return new AxiomManifestParsingStep(proj, manifestReq);
  }

  final Project proj;
  final Requirement<File> manifestReq;

  private AxiomManifestParsingStep(Project proj, Requirement<File> manifestReq) {
    super(NAME, DESC, RequirementSetBuilder.and().add(manifestReq),
          EnumSet.noneOf(Requirement.Flag.class));
    this.proj = proj;
    this.manifestReq = manifestReq;
  }

  @Override
  public void setNecessaryPreRunProperties(Variables variables) {
    proj.SNP_DATA_FILE.setValue(variables.getString(manifestReq));
  }

  @Override
  public void run(Variables variables) {
    FileColumn<String> colPsID = new AliasedFileColumn("Probe Set ID", "Probe Set ID");
    FileColumn<String> colAffyID = new AliasedFileColumn("Affy SNP ID", "Affy SNP ID");
    FileColumn<String> colRSID = new AliasedFileColumn("RS ID", "dbSNP RS ID");
    FileColumn<Byte> colChr = StandardFileColumns.chr("Chromosome");
    FileColumn<Integer> colPos = StandardFileColumns.pos("Position");

    try {
      FileParser parser = FileParserFactory.setup(variables.get(manifestReq).getPath(), colPsID,
                                                  colAffyID, colRSID, colChr, colPos)
                                           .skipPrefix("#").build();
      PrintWriter mkrPosWriter = Files.getAppropriateWriter(proj.MARKER_POSITION_FILENAME.getValue());
      PrintWriter lookupWriter = Files.getAppropriateWriter(proj.PROJECT_DIRECTORY.getValue()
                                                            + ext.rootOf(variables.get(manifestReq)
                                                                                  .getPath(),
                                                                         true)
                                                            + AxiomManifestParsingStep.RSLOOKUP_SUFFIX);
      for (DataLine dl : parser) {
        try {
          String probeset = dl.get(colPsID);
          String affyId = dl.get(colAffyID);
          String rsId = dl.get(colRSID);
          Byte chr = dl.get(colChr);
          Integer pos = dl.get(colPos);

          mkrPosWriter.println(probeset + "\t" + chr + "\t" + pos);
          lookupWriter.println(probeset + "\t" + (rsId.equals("---") ? affyId : rsId));
        } catch (ParseFailureException e) {
          e.printStackTrace();
        }
      }
      mkrPosWriter.close();
      lookupWriter.close();
    } catch (IOException e) {
      e.printStackTrace();
    }

  }

  @Override
  public boolean checkIfOutputExists(Variables variables) {
    return Files.exists(proj.MARKER_POSITION_FILENAME.getValue(false, false));
  }

  @Override
  public String getCommandLine(Variables variables) {
    return getStepCommandLine(proj, variables);
  }

  public static void main(String[] args) {
    Project proj = Step.parseProject(args);
    StepBuilder sb = new StepBuilder(proj);
    Step step = sb.generateAxiomManifestParsingStep();
    Variables variables = step.parseArguments(args);
    Step.run(proj, step, variables);
  }

}
