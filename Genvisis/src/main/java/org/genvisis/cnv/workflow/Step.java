package org.genvisis.cnv.workflow;

import java.util.Collection;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.gui.GenvisisWorkflowGUI;
import org.genvisis.common.gui.Task;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Sets;

public abstract class Step implements Comparable<Step> {

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

  private String name;
  private String desc;
  private RequirementSet requirements;
  private final Set<Step> relatedSteps; // Not included in equality to prevent infinite recursion
  private Set<Requirement.Flag> stepFlags;
  private final double priority;

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
              Collection<Requirement.Flag> flags, double priority) {
    this.name = name;
    this.desc = desc;
    this.requirements = requirements;
    this.stepFlags = EnumSet.copyOf(flags);
    ImmutableSet.Builder<Step> relatedStepsBuilder = ImmutableSet.builder();
    relatedStepsBuilder.add(this);
    for (Requirement req : requirements.getFlatRequirementsList()) {
      if (req instanceof Requirement.StepRequirement && req != null) {
        Step requiredStep = ((Requirement.StepRequirement) req).getRequiredStep();
        if (requiredStep != null) {
          relatedStepsBuilder.add(requiredStep);
          relatedStepsBuilder.addAll(requiredStep.getRelatedSteps());
        }
      }
    }
    this.relatedSteps = relatedStepsBuilder.build();
    this.stepFlags = Sets.immutableEnumSet(flags);
    this.priority = priority;
  }

  public String getName() {
    return this.name;
  }

  public String getDescription() {
    return this.desc;
  }

  public abstract void setNecessaryPreRunProperties(Project proj,
                                                    Map<Step, Map<Requirement, String>> variables);

  public abstract void run(Project proj, Map<Step, Map<Requirement, String>> variables);

  /**
   * Used to cancel a step
   * 
   * @param proj
   */
  public void gracefulDeath(Project proj) {
    return;
  }

  public boolean hasRequirements(Project proj, Set<Step> stepSelections,
                                 Map<Step, Map<Requirement, String>> variables) {
    if (variables.get(this) == null) {
      return false;
    }
    return this.requirements.satisfiesRequirements(proj, this, stepSelections, variables);
  }

  /**
   * @return An array of {@link Requirement}s. At least one element of each subarray must be met to
   *         satisfy the step pre-requisites - effectively this means elements of the first array
   *         are AND'd together, while elements of the second array are OR'd.
   */
  public RequirementSet getRequirements() {
    return requirements;
  }

  public abstract boolean checkIfOutputExists(Project proj,
                                              Map<Step, Map<Requirement, String>> variables);

  public abstract String getCommandLine(Project proj,
                                        Map<Step, Map<Requirement, String>> variables);

  public Map<Requirement, String> getDefaultRequirementValues() {
    Map<Requirement, String> varMap = new HashMap<>();
    for (Requirement r : this.requirements.getFlatRequirementsList()) {
      varMap.put(r, r.getDefaultValue().toString());
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

  public double getPriority() {
    return priority;
  }

  @Override
  public int compareTo(Step o) {
    // Preferably, just compare on priority. Otherwise, compare hash to prevent collisions
    int priorityCmp = Double.compare(getPriority(), o.getPriority());
    if (priorityCmp != 0) return priorityCmp;
    return Integer.compare(hashCode(), o.hashCode());
  }

  @Override
  public int hashCode() {
    final int prime = 31;
    int result = 1;
    result = prime * result + ((desc == null) ? 0 : desc.hashCode());
    result = prime * result + ((name == null) ? 0 : name.hashCode());
    long temp;
    temp = Double.doubleToLongBits(priority);
    result = prime * result + (int) (temp ^ (temp >>> 32));
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
    if (Double.doubleToLongBits(priority) != Double.doubleToLongBits(other.priority)) return false;
    if (!requirements.equals(other.requirements)) return false;
    if (stepFlags == null) {
      if (other.stepFlags != null) return false;
    } else if (!stepFlags.equals(other.stepFlags)) return false;
    return true;
  }

  public Task<Void, Void> createTask(GenvisisWorkflowGUI gui, Project proj,
                                     Map<Step, Map<Requirement, String>> variables,
                                     Set<Step> selectedSteps) {
    StepTask st = new StepTask(gui, this, proj, selectedSteps, variables);
    return st;
  }

}
