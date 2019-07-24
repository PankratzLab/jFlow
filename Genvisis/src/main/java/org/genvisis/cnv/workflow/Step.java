package org.genvisis.cnv.workflow;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.cli.ParseException;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.gui.GenvisisWorkflowGUI;
import org.genvisis.cnv.workflow.Requirement.ParsedRequirement;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
import org.pankratzlab.common.CLI;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.gui.Task;

import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

public abstract class Step {

  public static final Step EMPTY_STEP = new Step("", "", RequirementSetBuilder.and(),
                                                 new ArrayList<>()) {

    @Override
    public void setNecessaryPreRunProperties(Variables variables) {
      // noop
    }

    @Override
    public void run(Variables variables) {
      // noop
    }

    @Override
    public String getCommandLine(Variables variables) {
      return "";
    }

    @Override
    public boolean checkIfOutputExists(Variables variables) {
      return true;
    }
  };

  public static enum FINAL_CODE {
    COMPLETE("Complete"), FAILED("Failed"), CANCELLED("Cancelled");

    private String message;

    FINAL_CODE(String msg) {
      message = msg;
    }

    public String getMessage() {
      return message;
    }
  }

  private final String name;
  private final String desc;
  private final RequirementSet requirements;
  private final Set<Step> dependentSteps; // Not included in equality to prevent infinite recursion
  private final Set<Requirement.Flag> stepFlags;
  private final CLI cli;

  /**
   * @param name displayed in the workflow
   * @param desc description to display as alt-text
   * @param requirements 2-d array of {@link Requirement}s, where one {@code Requirement} in each
   *          2nd dimension array must be met for the {@code Step} to run
   * @param flags {@link Requirement.Flag}s for the {@code Step}, {@code Flag.MULTITHREADED} is
   *          included by default if {@link GenvisisWorkflow.#getNumThreadsReq()} is a requirement
   * @param priority determines order in the workflow
   */
  public Step(String name, String desc, RequirementSet requirements,
              Collection<Requirement.Flag> flags) {
    this.name = name;
    this.desc = desc;
    this.requirements = requirements;
    this.stepFlags = Sets.immutableEnumSet(flags);
    this.dependentSteps = new HashSet<>();
    for (Requirement<?> req : requirements.getFlatRequirementsList()) {
      if (req instanceof Requirement.StepRequirement && req != null) {
        Step requiredStep = ((Requirement.StepRequirement) req).getRequiredStep();
        if (requiredStep != null) {
          requiredStep.addDependent(this);
        }
      }
    }
    cli = new CLI(this.getClass());
    cli.addArg(CLI.ARG_PROJ, CLI.DESC_PROJ, CLI.EXAMPLE_PROJ);
    for (Requirement<?> r : this.getRequirements().getFlatRequirementsList()) {
      if (r instanceof ParsedRequirement) {
        cli.addArg(((ParsedRequirement<?>) r).getKey(), r.getDescription());
      }
    }

  }

  private void addDependent(Step step) {
    this.dependentSteps.add(step);
  }

  public Set<Step> getSelfAndDependents() {
    ImmutableSet.Builder<Step> dependentStepsBuilder = ImmutableSet.builder();
    dependentStepsBuilder.add(this);
    for (Step step : dependentSteps) {
      dependentStepsBuilder.addAll(step.getSelfAndDependents());
    }
    return dependentStepsBuilder.build();
  }

  public String getName() {
    return this.name;
  }

  public String getDescription() {
    return this.desc;
  }

  /**
   * Set any Project Property values, if necessary.
   * 
   * @param proj
   * @param variables Map of Requirement to String value for this Step only
   */
  public abstract void setNecessaryPreRunProperties(Variables variables);

  /**
   * Run this Step
   * 
   * @param proj
   * @param variables Map of Requirement to String value for this Step only
   */
  public abstract void run(Variables variables);

  /**
   * Used to cancel a step
   * 
   * @param proj {@link Project}
   */
  public void gracefulDeath(Project proj) {
    cleanupAfterFailure(proj);
    return;
  }

  /**
   * Removes incomplete files from a failed or cancelled execution
   * 
   * @param proj {@link Project}
   */
  public void cleanupAfterFailure(Project proj) {
    return;
  }

  /**
   * Check if the requirements for this step are satisfied. These requirements possibly include
   * other steps, which requires the full Step-to-variables map.
   * 
   * @param proj Project
   * @param stepSelections Set of selected steps
   * @param variables Full map of each selected step to Requirement values
   * @return
   */
  public boolean hasRequirements(Set<Step> stepSelections, Map<Step, Variables> variables) {
    if (variables.get(this) == null) {
      return false;
    }
    return this.requirements.satisfiesRequirements(this, stepSelections, variables);
  }

  /**
   * @return An array of {@link Requirement}s. At least one element of each subarray must be met to
   *         satisfy the step pre-requisites - effectively this means elements of the first array
   *         are AND'd together, while elements of the second array are OR'd.
   */
  public RequirementSet getRequirements() {
    return requirements;
  }

  /**
   * Check if the output from this step already exists
   * 
   * @param proj
   * @param variables Map of Requirement to String value for this Step only
   * @return
   */
  public abstract boolean checkIfOutputExists(Variables variables);

  /**
   * Get the command line invocation for this Step, based on the applied variables.
   * 
   * @param proj
   * @param variables Map of Requirement to String value for this Step only
   * @return
   */
  public abstract String getCommandLine(Variables variables);

  public String getStepCommandLine(Project proj, Variables variables) {
    Map<String, String> args = Maps.newHashMap();
    args.put(CLI.ARG_PROJ, proj.getPropertyFilename());
    for (Requirement<?> r : this.getRequirements().getFlatRequirementsList()) {
      if (r instanceof ParsedRequirement && variables.hasValid(r)) {
        args.put(((ParsedRequirement<?>) r).getKey(), variables.getString(r));
      }
    }
    return Files.getRunString() + " " + CLI.formCmdLine(this.getClass(), args);
  }

  public Variables parseArguments(String[] args) {
    try {
      cli.parse(args);
    } catch (ParseException e) {
      throw new IllegalArgumentException(e);
    }
    Variables vars = new Variables();
    for (Requirement<?> r : this.getRequirements().getFlatRequirementsList()) {
      if (r instanceof ParsedRequirement) {
        vars.parseOrFail(r, cli.get(((ParsedRequirement<?>) r).getKey()));
      }
    }
    return vars;
  }

  @SuppressWarnings("unchecked")
  public Variables getDefaultRequirementValues() {
    Variables varMap = new Variables();
    for (@SuppressWarnings("rawtypes")
    Requirement r : this.requirements.getFlatRequirementsList()) {
      varMap.put(r, r.getDefaultValue());
    }
    return varMap;
  }

  public Collection<Requirement.Flag> getFlags() {
    return stepFlags;
  }

  @Override
  public int hashCode() {
    final int prime = 31;
    int result = 1;
    result = prime * result + ((desc == null) ? 0 : desc.hashCode());
    result = prime * result + ((name == null) ? 0 : name.hashCode());
    result = prime * result + requirements.hashCode();
    result = prime * result + ((stepFlags == null) ? 0 : stepFlags.hashCode());
    return result;
  }

  @Override
  public boolean equals(Object obj) {
    if (this == obj) return true;
    if (obj == null) return false;
    if (getClass() != obj.getClass()) return false;
    Step other = (Step) obj;
    if (desc == null) {
      if (other.desc != null) return false;
    } else if (!desc.equals(other.desc)) return false;
    if (name == null) {
      if (other.name != null) return false;
    } else if (!name.equals(other.name)) return false;
    if (!requirements.equals(other.requirements)) return false;
    if (stepFlags == null) {
      if (other.stepFlags != null) return false;
    } else if (!stepFlags.equals(other.stepFlags)) return false;
    return true;
  }

  public Task<Void, Void> createTask(GenvisisWorkflowGUI gui, Variables variables,
                                     List<Step> selectedSteps) {
    StepTask st = new StepTask(gui, this, selectedSteps, variables);
    return st;
  }

  public static Project parseProject(String[] args) {
    Project proj = null;
    for (String a : args) {
      if (a.startsWith(CLI.ARG_PROJ)) {
        proj = new Project(a.split("=")[1]);
      }
    }
    if (proj == null) {
      throw new IllegalArgumentException("Error - no project properties argument found.");
    }
    return proj;
  }

}
