package org.genvisis.cnv.workflow;

import java.util.List;
import java.util.Map;
import java.util.Set;
import org.genvisis.cnv.filesys.Project;

public abstract class RequirementSet {

  public static class RequirementSetBuilder {

    private RequirementSetBuilder() {}

    public static RequirementSet or() {
      return new OrRequirementSet();
    }

    public static RequirementSet and() {
      return new AndRequirementSet();
    }

  }

  protected List<Requirement<?>> reqs = new java.util.ArrayList<>();
  protected List<RequirementSet> reqSets = new java.util.ArrayList<>();

  public RequirementSet add(Requirement<?> r) {
    // UnitaryRequirementSet is used to preserve the order of added requirements
    reqSets.add(new UnitaryRequirementSet().add(r));
    return this;
  }

  public RequirementSet add(RequirementSet rs) {
    reqSets.add(rs);
    return this;
  }

  public List<Requirement<?>> getFlatRequirementsList() {
    List<Requirement<?>> reqList = new java.util.ArrayList<>();
    for (Requirement<?> r : reqs) {
      reqList.add(r);
    }
    for (RequirementSet rs : reqSets) {
      reqList.addAll(rs.getFlatRequirementsList());
    }
    return reqList;
  }

  public int size() {
    return getFlatRequirementsList().size();
  }

  abstract boolean satisfiesRequirements(Project proj, Step step, Set<Step> stepSelections,
                                         Map<Step, Variables> variables);

  @Override
  public int hashCode() {
    final int prime = 31;
    int result = 1;
    result = prime * result + ((reqSets == null) ? 0 : reqSets.hashCode());
    result = prime * result + ((reqs == null) ? 0 : reqs.hashCode());
    return result;
  }

  @Override
  public boolean equals(Object obj) {
    if (this == obj) return true;
    if (obj == null) return false;
    if (getClass() != obj.getClass()) return false;
    RequirementSet other = (RequirementSet) obj;
    if (reqSets == null) {
      if (other.reqSets != null) return false;
    } else if (!reqSets.equals(other.reqSets)) return false;
    if (reqs == null) {
      if (other.reqs != null) return false;
    } else if (!reqs.equals(other.reqs)) return false;
    return true;
  }

  private static final class UnitaryRequirementSet extends RequirementSet {

    private UnitaryRequirementSet() {}

    @Override
    public RequirementSet add(Requirement<?> r) {
      if (this.reqs.size() > 0) {
        throw new RuntimeException("UnitaryRequirementSet can only take one Requirement object.");
      }
      this.reqs.add(r);
      return this;
    }

    @Override
    public String getJoinString() {
      return "";
    }

    @Override
    boolean satisfiesRequirements(Project proj, Step step, Set<Step> stepSelections,
                                  Map<Step, Variables> variables) {
      return this.reqs.get(0).checkRequirement(proj,
                                               variables.get(step).get(this.reqs.get(0)).toString(),
                                               stepSelections, variables);
    }

    @Override
    public List<RequirementSet> getRequirementSets() {
      return new java.util.ArrayList<>();
    }

  }

  public static final class OrRequirementSet extends RequirementSet {

    private OrRequirementSet() {}

    @Override
    boolean satisfiesRequirements(Project proj, Step step, Set<Step> stepSelections,
                                  Map<Step, Variables> variables) {
      for (Requirement<?> r : reqs) {
        if (r.checkRequirement(proj, variables.get(step).get(r).toString(), stepSelections,
                               variables)) {
          return true;
        }
      }
      for (RequirementSet rs : reqSets) {
        if (rs.satisfiesRequirements(proj, step, stepSelections, variables)) {
          return true;
        }
      }
      return false;
    }

    @Override
    public String getJoinString() {
      return "OR";
    }

    @Override
    public List<Requirement<?>> getRequirements() {
      return new java.util.ArrayList<>();
    }

  }

  public static final class AndRequirementSet extends RequirementSet {

    private AndRequirementSet() {}

    @Override
    boolean satisfiesRequirements(Project proj, Step step, Set<Step> stepSelections,
                                  Map<Step, Variables> variables) {
      for (Requirement<?> r : reqs) {
        if (!r.checkRequirement(proj, variables.get(step).get(r).toString(), stepSelections,
                                variables)) {
          return false;
        }
      }
      for (RequirementSet rs : reqSets) {
        if (!rs.satisfiesRequirements(proj, step, stepSelections, variables)) {
          return false;
        }
      }
      return true;
    }

    @Override
    public String getJoinString() {
      return "AND";
    }

    @Override
    public List<Requirement<?>> getRequirements() {
      return new java.util.ArrayList<>();
    }
  }

  public List<Requirement<?>> getRequirements() {
    return reqs;
  }

  public List<RequirementSet> getRequirementSets() {
    return reqSets;
  }

  public abstract String getJoinString();

}
