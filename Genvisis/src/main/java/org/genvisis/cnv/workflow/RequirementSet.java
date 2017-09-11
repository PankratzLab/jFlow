package org.genvisis.cnv.workflow;

import java.util.List;
import java.util.Map;
import java.util.Set;

public abstract class RequirementSet {

	public static class OrRequirementSet extends RequirementSet {
		@Override
		boolean satisfiesRequirements(String arg, Set<Step> stepSelections,
																	Map<Step, Map<Requirement, String>> variables) {
			for (Requirement r : reqs) {
				if (r.checkRequirement(arg, stepSelections, variables)) {
					return true;
				}
			}
			for (RequirementSet rs : reqSets) {
				if (rs.satisfiesRequirements(arg, stepSelections, variables)) {
					return true;
				}
			}
			return false;
		}
	}

	public static class AndRequirementSet extends RequirementSet {
		@Override
		boolean satisfiesRequirements(String arg, Set<Step> stepSelections,
																	Map<Step, Map<Requirement, String>> variables) {
			for (Requirement r : reqs) {
				if (!r.checkRequirement(arg, stepSelections, variables)) {
					return false;
				}
			}
			for (RequirementSet rs : reqSets) {
				if (!rs.satisfiesRequirements(arg, stepSelections, variables)) {
					return false;
				}
			}
			return true;
		}
	}

	List<Requirement> reqs;
	List<RequirementSet> reqSets;

	void stubUsageForPedigree() {
		Requirement r1 = null, r211 = null, r212 = null, r22 = null, r31 = null, r32 = null;
		RequirementSet overall = new RequirementSet.OrRequirementSet();
		RequirementSet rs2 = new RequirementSet.AndRequirementSet(), rs21 = new RequirementSet.OrRequirementSet(), rs3 = new RequirementSet.AndRequirementSet();

		overall.addRequirement(r1);

		rs21.addRequirement(r211);
		rs21.addRequirement(r212);
		rs2.addRequirement(rs21);
		rs2.addRequirement(r22);
		overall.addRequirement(rs2);

		rs3.addRequirement(r31);
		rs3.addRequirement(r32);
		rs3.addRequirement(rs21);
		overall.addRequirement(rs3);
	}

	void addRequirement(Requirement r) {
		reqs.add(r);
	}

	void addRequirement(RequirementSet rs) {
		reqSets.add(rs);
	}

	abstract boolean satisfiesRequirements(String arg, Set<Step> stepSelections,
																				 Map<Step, Map<Requirement, String>> variables);
}
