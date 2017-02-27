package org.genvisis.cnv.manage;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.Files;

public abstract class GenvisisWorkflowStep implements Comparable<GenvisisWorkflowStep> {
	public abstract static class Requirement {
		private final String description;
		private final RequirementInputType type;
		private final Object defaultValue;

		/**
		 * @param description
		 * @param type
		 */
		protected Requirement(String description, RequirementInputType type) {
			this(description, type, null);
		}

		/**
		 * @param description
		 * @param type
		 * @param defaultValue
		 */
		protected Requirement(String description, RequirementInputType type, Object defaultValue) {
			super();
			this.description = description;
			this.type = type;
			this.defaultValue = defaultValue != null ? defaultValue : "";
		}

		protected static int checkIntArgOrNeg1(String val) {
			int valInt = -1;
			try {
				valInt = Integer.parseInt(val);
			} catch (NumberFormatException e) {
				// leave as -1
			}
			return valInt;
		}

		protected double checkDoubleArgOrNeg1(String val) {
			double valDou = -1.0;
			try {
				valDou = Double.parseDouble(val);
			} catch (NumberFormatException e) {
				// leave as -1.0
			}
			return valDou;
		}

		public String getDescription() {
			return description;
		}

		public RequirementInputType getType() {
			return type;
		}

		public Object getDefaultValue() {
			return defaultValue;
		}


		public abstract boolean checkRequirement(String arg, Set<GenvisisWorkflowStep> stepSelections,
																						 Map<GenvisisWorkflowStep, Map<Requirement, String>> variables);

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + ((defaultValue == null) ? 0 : defaultValue.hashCode());
			result = prime * result + ((description == null) ? 0 : description.hashCode());
			result = prime * result + ((type == null) ? 0 : type.hashCode());
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			Requirement other = (Requirement) obj;
			if (defaultValue == null) {
				if (other.defaultValue != null)
					return false;
			} else if (!defaultValue.equals(other.defaultValue))
				return false;
			if (description == null) {
				if (other.description != null)
					return false;
			} else if (!description.equals(other.description))
				return false;
			if (type != other.type)
				return false;
			return true;
		}

	}

	public static class StepRequirement extends Requirement {
		private final GenvisisWorkflowStep requiredStep;

		public StepRequirement(GenvisisWorkflowStep requiredStep) {
			super(stepReqMessage(requiredStep), RequirementInputType.NONE);
			this.requiredStep = requiredStep;
		}

		@Override
		public boolean checkRequirement(String arg, Set<GenvisisWorkflowStep> stepSelections,
																		Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
			return stepSelections.contains(requiredStep) || requiredStep.checkIfOutputExists(variables);
		}

		private static String stepReqMessage(GenvisisWorkflowStep requiredStep) {
			return "[" + requiredStep.getName()
						 + "] step must have been run already or must be selected and valid";
		}


	}

	public static class FileRequirement extends Requirement {

		public FileRequirement(String description, String defaultValue) {
			super(description, RequirementInputType.FILE, defaultValue);
		}

		@Override
		public boolean checkRequirement(String arg, Set<GenvisisWorkflowStep> stepSelections,
																		Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
			return Files.exists(arg);
		}
	}

	public static class OutputFileRequirement extends FileRequirement {

		public OutputFileRequirement(String description, String defaultValue) {
			super(description, defaultValue);
		}

		@Override
		public boolean checkRequirement(String arg, Set<GenvisisWorkflowStep> stepSelections,
																		Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
			return !Files.exists(arg);
		}
	}

	public static class OptionalFileRequirement extends FileRequirement {
		public OptionalFileRequirement(String description, String defaultValue) {
			super(description, defaultValue);
		}

		@Override
		public boolean checkRequirement(String arg, Set<GenvisisWorkflowStep> stepSelections,
																		Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
			return "".equals(arg) || Files.exists(arg);
		}
	}

	public static class BoolRequirement extends Requirement {

		protected BoolRequirement(String description, boolean defaultValue) {
			super(description, RequirementInputType.BOOL, defaultValue);
		}

		@Override
		public boolean checkRequirement(String arg, Set<GenvisisWorkflowStep> stepSelections,
																		Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
			return Boolean.parseBoolean(arg);
		}

	}

	public static class OptionalBoolRequirement extends BoolRequirement {
		protected OptionalBoolRequirement(String description, boolean defaultValue) {
			super(description, defaultValue);
		}

		@Override
		public boolean checkRequirement(String arg, Set<GenvisisWorkflowStep> stepSelections,
																		Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
			return true;
		}
	}

	public static class DoubleRequirement extends Requirement {

		private final double min;
		private final double max;

		protected DoubleRequirement(String description, double defaultValue, double min, double max) {
			super(description, RequirementInputType.NUMBER, defaultValue);
			this.min = min;
			this.max = max;
		}

		@Override
		public boolean checkRequirement(String arg, Set<GenvisisWorkflowStep> stepSelections,
																		Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
			double value;
			try {
				value = Double.parseDouble(arg);
			} catch (NumberFormatException e) {
				return false;
			}
			return value >= min && value <= max;
		}

	}

	public static class IntRequirement extends Requirement {

		private final int min;
		private final int max;

		protected IntRequirement(String description, int defaultValue, int min, int max) {
			super(description, RequirementInputType.NUMBER, defaultValue);
			this.min = min;
			this.max = max;
		}

		@Override
		public boolean checkRequirement(String arg, Set<GenvisisWorkflowStep> stepSelections,
																		Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
			int value;
			try {
				value = Integer.parseInt(arg);
			} catch (NumberFormatException e) {
				return false;
			}
			return value >= min && value <= max;
		}

	}

	public static class EnumRequirement extends Requirement {

		protected EnumRequirement(String description, Enum<?> defaultValue) {
			super(description, RequirementInputType.ENUM, defaultValue);
		}

		@Override
		public boolean checkRequirement(String arg, Set<GenvisisWorkflowStep> stepSelections,
																		Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
			return true;
		}

	}

	public enum RequirementInputType {
		NONE, FILE, DIR, STRING, NUMBER, BOOL, ENUM
	}

	public enum FLAG {
		MEMORY, RUNTIME, MULTITHREADED
	}

	private String name;
	private String desc;
	private GenvisisWorkflowStep.Requirement[][] requirements;
	private boolean failed = false;
	private ArrayList<String> failReasons = new ArrayList<String>();
	private final Set<GenvisisWorkflowStep> relatedSteps;
	private final Set<GenvisisWorkflowStep.FLAG> stepFlags;
	private final double priority;

	GenvisisWorkflowStep(String name, String desc, GenvisisWorkflowStep.Requirement[][] requirements,
											 GenvisisWorkflowStep[] relatedSteps,
											 GenvisisWorkflowStep.FLAG[] flags, double priority) {
		this.name = name;
		this.desc = desc;
		this.requirements = requirements;
		final Set<GenvisisWorkflowStep> steps = new HashSet<GenvisisWorkflowStep>();
		steps.add(this);
		if (relatedSteps != null) {
			for (final GenvisisWorkflowStep s : relatedSteps) {
				steps.add(s);
				steps.addAll(s.getRelatedSteps());
			}
		}
		this.relatedSteps = Collections.unmodifiableSet(steps);
		final Set<GenvisisWorkflowStep.FLAG> flgs = new HashSet<GenvisisWorkflowStep.FLAG>();
		if (flags != null) {
			for (final GenvisisWorkflowStep.FLAG f : flags) {
				flgs.add(f);
			}
		}
		this.stepFlags = Collections.unmodifiableSet(flgs);
		this.priority = priority;
	}

	public String getName() {
		return this.name;
	}

	public String getDescription() {
		return this.desc;
	}

	public boolean getFailed() {
		return failed;
	}

	protected void setFailed(String reason) {
		failed = true;
		failReasons.add(reason);
	}

	public List<String> getFailureMessages() {
		return failReasons;
	}

	public abstract void setNecessaryPreRunProperties(Project proj,
																										Map<GenvisisWorkflowStep, Map<GenvisisWorkflowStep.Requirement, String>> variables);

	public abstract void run(Project proj,
													 Map<GenvisisWorkflowStep, Map<GenvisisWorkflowStep.Requirement, String>> variables);

	public void gracefulDeath(Project proj) {
		return;
	}

	public boolean hasRequirements(Project proj, Set<GenvisisWorkflowStep> stepSelections,
																 Map<GenvisisWorkflowStep, Map<GenvisisWorkflowStep.Requirement, String>> variables) {
		if (variables.get(this) == null) {
			return false;
		}
		for (GenvisisWorkflowStep.Requirement[] group : getRequirements()) {
			boolean groupMet = false;
			for (GenvisisWorkflowStep.Requirement req : group) {
				if (req.checkRequirement(variables.get(this).get(req), stepSelections, variables)) {
					groupMet = true;
					break;
				}
			}
			if (!groupMet)
				return false;
		}
		return true;
	}

	/**
	 * @return An array of {@link GenvisisWorkflowStep.Requirement}s. At least one element of each
	 *         subarray must be met to satisfy the step pre-requisites - effectively this means
	 *         elements of the first array are AND'd together, while elements of the second array are
	 *         OR'd.
	 */
	public GenvisisWorkflowStep.Requirement[][] getRequirements() {
		// TODO unify requirement names, AND/OR structure, input types and default values to avoid
		// maintaining these parallel arrays
		return requirements;
	}

	public abstract boolean checkIfOutputExists(Map<GenvisisWorkflowStep, Map<GenvisisWorkflowStep.Requirement, String>> variables);

	public void resetRun() {
		failed = false;
		failReasons.clear();
	}

	public abstract String getCommandLine(Project proj,
																				Map<GenvisisWorkflowStep, Map<GenvisisWorkflowStep.Requirement, String>> variables);

	/**
	 * @return A {@link Collection} of the complete network of {@link GenvisisWorkflowStep}s related
	 *         to this {@code GenvisisWorkflowStep} - including this {@code GenvisisWorkflowStep},
	 *         direct and transitive dependencies.
	 */
	public Collection<GenvisisWorkflowStep> getRelatedSteps() {
		return relatedSteps;
	}

	public Collection<GenvisisWorkflowStep.FLAG> getFlags() {
		return stepFlags;
	}



	public double getPriority() {
		return priority;
	}

	@Override
	public int compareTo(GenvisisWorkflowStep o) {
		return Double.compare(priority, o.priority);
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((desc == null) ? 0 : desc.hashCode());
		result = prime * result + ((failReasons == null) ? 0 : failReasons.hashCode());
		result = prime * result + (failed ? 1231 : 1237);
		result = prime * result + ((name == null) ? 0 : name.hashCode());
		long temp;
		temp = Double.doubleToLongBits(priority);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		result = prime * result + Arrays.deepHashCode(requirements);
		result = prime * result + ((stepFlags == null) ? 0 : stepFlags.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		GenvisisWorkflowStep other = (GenvisisWorkflowStep) obj;
		if (desc == null) {
			if (other.desc != null)
				return false;
		} else if (!desc.equals(other.desc))
			return false;
		if (failReasons == null) {
			if (other.failReasons != null)
				return false;
		} else if (!failReasons.equals(other.failReasons))
			return false;
		if (failed != other.failed)
			return false;
		if (name == null) {
			if (other.name != null)
				return false;
		} else if (!name.equals(other.name))
			return false;
		if (Double.doubleToLongBits(priority) != Double.doubleToLongBits(other.priority))
			return false;
		if (!Arrays.deepEquals(requirements, other.requirements))
			return false;
		if (stepFlags == null) {
			if (other.stepFlags != null)
				return false;
		} else if (!stepFlags.equals(other.stepFlags))
			return false;
		return true;
	}
}
