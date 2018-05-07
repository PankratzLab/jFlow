package org.genvisis.cnv.workflow.steps;

import java.util.EnumSet;
import java.util.List;
import java.util.Map;
import org.genvisis.CLI;
import org.genvisis.cnv.filesys.ABLookup;
import org.genvisis.cnv.filesys.ABLookup.ABSource;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.workflow.Requirement;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.genvisis.cnv.workflow.Step;
import org.genvisis.common.Files;
import org.genvisis.common.ext;
import com.google.common.base.Joiner;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;

public class ABLookupStep extends Step {

  public static ABLookupStep create(Step parseSamplesStep) {
    Requirement parseSamplesStepReq = new Requirement.StepRequirement(parseSamplesStep);
    return new ABLookupStep(parseSamplesStepReq);
  }

  private ABLookupStep(Requirement parseSamplesReq) {
    super("Generate AB Lookup File", "", RequirementSetBuilder.and().add(parseSamplesReq),
          EnumSet.of(Requirement.Flag.RUNTIME));
  }

  @Override
  public void setNecessaryPreRunProperties(Project proj,
                                           Map<Step, Map<Requirement, String>> variables) {
    // Nothing to do here
  }

  @Override
  public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
    String filename = proj.PROJECT_DIRECTORY.getValue()
                      + ext.addToRoot(ABLookup.DEFAULT_AB_FILE, "_parsed");
    ABLookup.parseABLookup(proj, ABSource.VCF, filename);

    if (ABLookup.fillInMissingAlleles(proj, filename, proj.getLocationOfSNP_Map(true), false)) {
      ABLookup.applyABLookupToFullSampleFiles(proj, filename);
    } else {
      throw new RuntimeException("Failed to fill in missing alleles - please check log for more info.");
    }
  }

  @Override
  public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
    String filename = proj.PROJECT_DIRECTORY.getValue()
                      + ext.addToRoot(ABLookup.DEFAULT_AB_FILE, "_parsed");
    String projFile = proj.getPropertyFilename();
    String mapFile = proj.getLocationOfSNP_Map(true);

    List<String> baseCommand = ImmutableList.of(Files.getRunString(), ABLookup.class.getName(),
                                                CLI.formCmdLineArg(CLI.ARG_PROJ, projFile));
    List<String> commandVcf = Lists.newArrayList(baseCommand);
    commandVcf.add(CLI.formCmdLineArg(CLI.ARG_OUTFILE, filename));
    commandVcf.add(CLI.formCmdLineFlag(ABLookup.FLAGS_VCF));

    List<String> commandPartial = Lists.newArrayList(baseCommand);
    commandPartial.add(CLI.formCmdLineArg(ABLookup.ARGS_PARTAB, filename));
    commandPartial.add(CLI.formCmdLineArg(ABLookup.ARGS_MAP, mapFile));

    List<String> commandProp = Lists.newArrayList(ImmutableList.of(Files.getRunString(),
                                                                   Project.class.getName(),
                                                                   CLI.formCmdLineArg(CLI.ARG_PROJ,
                                                                                      projFile)));
    commandProp.add(CLI.formCmdLineArg(proj.AB_LOOKUP_FILENAME.getName(), filename));

    List<String> commandApply = Lists.newArrayList(baseCommand);
    commandApply.add(CLI.formCmdLineFlag(ABLookup.FLAGS_APPLYAB));

    StringBuilder cmd = new StringBuilder();

    cmd.append(Joiner.on(" ").join(commandVcf)).append("\n");
    cmd.append(Joiner.on(" ").join(commandPartial)).append("\n");
    cmd.append(Joiner.on(" ").join(commandProp)).append("\n");
    cmd.append(Joiner.on(" ").join(commandApply));

    return cmd.toString();
  }

  @Override
  public boolean checkIfOutputExists(Project proj, Map<Step, Map<Requirement, String>> variables) {
    return Files.exists(proj.AB_LOOKUP_FILENAME.getValue(false, false));
  }
}
