package org.genvisis.cnv.workflow;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.gui.GenvisisWorkflowGUI;
import org.genvisis.common.gui.Task;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Sets;

public abstract class Step {

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
  private final Set<Step> relatedSteps; // Not included in equality to prevent infinite recursion
  private final Set<Requirement.Flag> stepFlags;

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
    ImmutableSet.Builder<Step> relatedStepsBuilder = ImmutableSet.builder();
    relatedStepsBuilder.add(this);
    for (Requirement<?> req : requirements.getFlatRequirementsList()) {
      if (req instanceof Requirement.StepRequirement && req != null) {
        Step requiredStep = ((Requirement.StepRequirement) req).getRequiredStep();
        if (requiredStep != null) {
          relatedStepsBuilder.add(requiredStep);
          relatedStepsBuilder.addAll(requiredStep.getRelatedSteps());
        }
      }
    }
    this.relatedSteps = relatedStepsBuilder.build();
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

  @SuppressWarnings("unchecked")
  public Variables getDefaultRequirementValues() {
    Variables varMap = new Variables();
    for (@SuppressWarnings("rawtypes")
    Requirement r : this.requirements.getFlatRequirementsList()) {
      varMap.put(r, r.getDefaultValue());
    }
    return varMap;
  }

  /**
   * @return A {@link Collection} of the complete network of {@link Step)s related to this {@code
   *         Step) - including this {@code Step), direct and transitive dependencies.
   */
  public Collection<Step> getRelatedSteps() {
    return relatedSteps;
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

}
