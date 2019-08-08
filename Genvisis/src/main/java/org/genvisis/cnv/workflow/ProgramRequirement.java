package org.genvisis.cnv.workflow;

import java.util.Map;
import java.util.Set;

import org.genvisis.cnv.workflow.Requirement.ParsedRequirement;
import org.pankratzlab.common.Files;

public class ProgramRequirement extends ParsedRequirement<String> {

  public ProgramRequirement(String key, String description, String defaultValue) {
    super(key, description, Requirement.RequirementInputType.FILE, defaultValue);
  }

  @Override
  public boolean checkRequirement(String arg, Set<Step> stepSelections,
                                  Map<Step, Variables> variables) {
    return Files.programExists(arg);
  }

  @Override
  public String parseValue(String raw) {
    return raw;
  }

}
