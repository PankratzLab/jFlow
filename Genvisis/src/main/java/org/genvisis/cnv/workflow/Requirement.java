package org.genvisis.cnv.workflow;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.genvisis.cnv.manage.Resources.Resource;
import org.genvisis.common.Files;
import org.genvisis.gwas.MarkerQC;
import org.genvisis.stats.Maths;

import com.google.common.base.Joiner;
import com.google.common.base.Splitter;

public abstract class Requirement {

	public static class StepRequirement extends Requirement {
		private final Step requiredStep;

		public StepRequirement(Step requiredStep) {
			super(stepReqMessage(requiredStep), Requirement.RequirementInputType.NONE);
			this.requiredStep = requiredStep;
		}

		@Override
		public boolean checkRequirement(String arg, Set<Step> stepSelections,
																		Map<Step, Map<Requirement, String>> variables) {
			return stepSelections.contains(requiredStep) || requiredStep.checkIfOutputExists(variables);
		}

		public Step getRequiredStep() {
			return requiredStep;
		}

		private static String stepReqMessage(Step requiredStep) {
			String msg = "[" + (requiredStep == null ? "" : requiredStep.getName())
									 + "] step must have been run already or must be selected";
			return msg;
		}


	}

	public static class FileRequirement extends Requirement {

		public FileRequirement(String description, String defaultValue) {
			super(description, Requirement.RequirementInputType.FILE, defaultValue);
		}

		@Override
		public boolean checkRequirement(String arg, Set<Step> stepSelections,
																		Map<Step, Map<Requirement, String>> variables) {
			return exists(arg);
		}

		protected boolean exists(String arg) {
			return Files.exists(arg);
		}
	}

	public static class OutputFileRequirement extends FileRequirement {

		public OutputFileRequirement(String description, String defaultValue) {
			super(description, defaultValue);
		}

		@Override
		public boolean checkRequirement(String arg, Set<Step> stepSelections,
																		Map<Step, Map<Requirement, String>> variables) {
			return !exists(arg);
		}
	}

	public static class OptionalFileRequirement extends FileRequirement {
		public OptionalFileRequirement(String description, String defaultValue) {
			super(description, defaultValue);
		}

		@Override
		public boolean checkRequirement(String arg, Set<Step> stepSelections,
																		Map<Step, Map<Requirement, String>> variables) {
			return "".equals(arg) || super.checkRequirement(arg, stepSelections, variables);
		}
	}

	public static class ResourceRequirement extends Requirement {

		private final Resource resource;

		/**
		 * 
		 * @param resourceDescription a description of the {@link Resource} required. This will be used
		 *        as part of the full description describing the download/use of the resource
		 * @param resource the {@link Resource} required
		 */
		public ResourceRequirement(String resourceDescription, Resource resource) {
			super(generateRequirementDescription(resourceDescription, resource),
						Requirement.RequirementInputType.NONE);
			this.resource = resource;
		}

		private static String generateRequirementDescription(String resourceDescription,
																												 Resource resource) {
			String requirementDescription;
			if (resource.isLocallyAvailable()) {
				requirementDescription = "Use locally available " + resourceDescription;
			} else {
				requirementDescription = "Download remotely available " + resourceDescription;
			}
			return requirementDescription;
		}

		@Override
		public boolean checkRequirement(String arg, Set<Step> stepSelections,
																		Map<Step, Map<Requirement, String>> variables) {
			return true;
		}

		public Resource getResource() {
			return resource;
		}

	}

	public static class BoolRequirement extends Requirement {

		public BoolRequirement(String description, boolean defaultValue) {
			super(description, Requirement.RequirementInputType.BOOL, defaultValue);
		}

		@Override
		public boolean checkRequirement(String arg, Set<Step> stepSelections,
																		Map<Step, Map<Requirement, String>> variables) {
			return Boolean.parseBoolean(arg);
		}

	}

	public static class OptionalBoolRequirement extends BoolRequirement {
		public OptionalBoolRequirement(String description, boolean defaultValue) {
			super(description, defaultValue);
		}

		@Override
		public boolean checkRequirement(String arg, Set<Step> stepSelections,
																		Map<Step, Map<Requirement, String>> variables) {
			return true;
		}
	}

	public static class DoubleRequirement extends Requirement {

		private final double min;
		private final double max;

		public DoubleRequirement(String description, double defaultValue, double min, double max) {
			super(description, Requirement.RequirementInputType.NUMBER, defaultValue);
			this.min = min;
			this.max = max;
		}

		@Override
		public boolean checkRequirement(String arg, Set<Step> stepSelections,
																		Map<Step, Map<Requirement, String>> variables) {
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

		public IntRequirement(String description, int defaultValue, int min, int max) {
			super(description, Requirement.RequirementInputType.NUMBER, defaultValue);
			this.min = min;
			this.max = max;
		}

		@Override
		public boolean checkRequirement(String arg, Set<Step> stepSelections,
																		Map<Step, Map<Requirement, String>> variables) {
			int value;
			try {
				value = Integer.parseInt(arg);
			} catch (NumberFormatException e) {
				return false;
			}
			return value >= min && value <= max;
		}

	}

	public static class PosIntRequirement extends IntRequirement {

		public PosIntRequirement(String description, int defaultValue) {
			super(description, defaultValue, 1, Integer.MAX_VALUE);
		}

	}

	public static class ListSelectionRequirement extends Requirement {

		private static final char SELECTION_LIST_DELIM = ',';
		private static final Joiner SELECTION_LIST_JOINER = Joiner.on(SELECTION_LIST_DELIM);
		private static final Splitter SELECTION_LIST_SPLITTER = Splitter.on(SELECTION_LIST_DELIM);

		private final Collection<String> options;
		private final boolean allowNone;

		public ListSelectionRequirement(String description, Collection<String> options,
																		Collection<String> defaultOptions, boolean allowNone) {
			super(description, Requirement.RequirementInputType.LISTSELECTION, defaultOptions);
			if (!options.containsAll(defaultOptions))
				throw new IllegalArgumentException("All defaultOptions are not in options");
			this.options = options;
			this.allowNone = allowNone;
		}

		@Override
		public boolean checkRequirement(String arg, Set<Step> stepSelections,
																		Map<Step, Map<Requirement, String>> variables) {
			return allowNone || !arg.isEmpty();
		}

		public Collection<String> getOptions() {
			return options;
		}

		@SuppressWarnings("unchecked")
		public Collection<String> getDefaultOptions() {
			return (Collection<String>) getDefaultValue();
		}

		public static String createArgValString(Iterable<?> selections) {
			return SELECTION_LIST_JOINER.join(selections);
		}

		public static List<String> parseArgValString(String arg) {
			return SELECTION_LIST_SPLITTER.splitToList(arg);
		}

	}

	public static class EnumRequirement extends Requirement {

		public EnumRequirement(String description, Enum<?> defaultValue) {
			super(description, RequirementInputType.ENUM, defaultValue);
		}

		@Override
		public boolean checkRequirement(String arg, Set<Step> stepSelections,
																		Map<Step, Map<Requirement, String>> variables) {
			return true;
		}

	}

	public static class ThresholdRequirement extends Requirement {

		public ThresholdRequirement(String description, String defaultValue) {
			super(description, RequirementInputType.STRING, defaultValue);
		}

		@Override
		public boolean checkRequirement(String arg, Set<Step> stepSelections,
																		Map<Step, Map<Requirement, String>> variables) {
			Maths.COMPARISON op = MarkerQC.findOperator(arg);
			if (op == null)
				return false;
			try {
				Double.parseDouble(arg.substring(op.getSymbol().length()));
			} catch (NumberFormatException nfe) {
				return false;
			}
			return true;
		}

	}

	public enum RequirementInputType {
		NONE, FILE, DIR, STRING, NUMBER, BOOL, ENUM, LISTSELECTION
	}

	public enum Flag {
		MEMORY, RUNTIME, MULTITHREADED
	}

	private final String description;
	private final Requirement.RequirementInputType type;
	private final Object defaultValue;

	/**
	 * @param description
	 * @param type
	 */
	public Requirement(String description, Requirement.RequirementInputType type) {
		this(description, type, null);
	}

	/**
	 * @param description
	 * @param type
	 * @param defaultValue
	 */
	public Requirement(String description, Requirement.RequirementInputType type,
										 Object defaultValue) {
		super();
		this.description = description;
		this.type = type;
		this.defaultValue = defaultValue != null ? defaultValue : "";
	}

	public static int checkIntArgOrNeg1(String val) {
		int valInt = -1;
		try {
			valInt = Integer.parseInt(val);
		} catch (NumberFormatException e) {
			// leave as -1
		}
		return valInt;
	}

	public static double checkDoubleArgOrNeg1(String val) {
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

	public Requirement.RequirementInputType getType() {
		return type;
	}

	public Object getDefaultValue() {
		return defaultValue;
	}


	public abstract boolean checkRequirement(String arg, Set<Step> stepSelections,
																					 Map<Step, Map<Requirement, String>> variables);

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
