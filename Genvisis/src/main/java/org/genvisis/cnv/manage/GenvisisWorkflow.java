package org.genvisis.cnv.manage;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.EnumSet;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;

import org.genvisis.CLI;
import org.genvisis.cnv.Launch;
import org.genvisis.cnv.analysis.AnalysisFormats;
import org.genvisis.cnv.analysis.Mosaicism;
import org.genvisis.cnv.analysis.pca.PCA;
import org.genvisis.cnv.analysis.pca.PCAPrep;
import org.genvisis.cnv.analysis.pca.PCImputeRace;
import org.genvisis.cnv.analysis.pca.PCImputeRace.RACE;
import org.genvisis.cnv.analysis.pca.PrincipalComponentsIntensity.CHROMOSOME_X_STRATEGY;
import org.genvisis.cnv.analysis.pca.PrincipalComponentsIntensity.CORRECTION_TYPE;
import org.genvisis.cnv.filesys.ABLookup;
import org.genvisis.cnv.filesys.ABLookup.ABSource;
import org.genvisis.cnv.filesys.Centroids;
import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.Pedigree;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.gui.GenvisisWorkflowGUI;
import org.genvisis.cnv.hmm.CNVCaller;
import org.genvisis.cnv.hmm.CNVCaller.PFB_MANAGEMENT_TYPE;
import org.genvisis.cnv.prop.Property;
import org.genvisis.cnv.qc.GcAdjustor;
import org.genvisis.cnv.qc.LrrSd;
import org.genvisis.cnv.qc.MarkerBlast;
import org.genvisis.cnv.qc.MarkerBlastQC;
import org.genvisis.cnv.qc.MarkerMetrics;
import org.genvisis.cnv.qc.SampleQC;
import org.genvisis.cnv.util.Java6Helper;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Elision;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.ext;
import org.genvisis.gwas.Ancestry;
import org.genvisis.gwas.FurtherAnalysisQc;
import org.genvisis.gwas.MarkerQC;
import org.genvisis.gwas.MarkerQC.QC_METRIC;
import org.genvisis.gwas.PlinkMendelianChecker;
import org.genvisis.gwas.Qc;
import org.genvisis.gwas.RelationAncestryQc;
import org.genvisis.stats.Maths;

import com.google.common.base.Joiner;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

public class GenvisisWorkflow {

	public abstract class Step implements Comparable<Step> {
		private String name;
		private String desc;
		private Requirement[][] requirements;
		private boolean failed = false;
		private ArrayList<String> failReasons = new ArrayList<String>();
		private final Set<Step> relatedSteps; // Not included in equality to prevent infinite recursion
		private Set<Flag> stepFlags;
		private final double priority;

		/**
		 * 
		 * @param name displayed in the workflow
		 * @param desc description to display as alt-text
		 * @param requirements 2-d array of {@link Requirement}s, where one {@code Requirement} in each
		 *        2nd dimension array must be met for the {@code Step} to run
		 * @param flags {@link Flag}s for the {@code Step}, {@code Flag.MULTITHREADED} is included by
		 *        default if {@link GenvisisWorkflow.#getNumThreadsReq()} is a requirement
		 * @param priority determines order in the workflow
		 */
		public Step(String name, String desc, Requirement[][] requirements, Collection<Flag> flags,
								double priority) {
			this.name = name;
			this.desc = desc;
			this.requirements = requirements;
			final Set<Step> buildRelatedSteps = new HashSet<Step>();
			this.stepFlags = EnumSet.copyOf(flags);
			for (Requirement[] group : requirements) {
				for (Requirement req : group) {
					if (req instanceof StepRequirement) {
						Step requiredStep = ((StepRequirement) req).getRequiredStep();
						buildRelatedSteps.add(requiredStep);
						buildRelatedSteps.addAll(requiredStep.getRelatedSteps());
					} else if (req == getNumThreadsReq()) {
						stepFlags.add(Flag.MULTITHREADED);
					}
				}
			}
			this.relatedSteps = Collections.unmodifiableSet(buildRelatedSteps);
			this.stepFlags = Sets.immutableEnumSet(flags);
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

		public void setFailed(String reason) {
			failed = true;
			failReasons.add(reason);
		}

		public List<String> getFailureMessages() {
			return failReasons;
		}

		public abstract void setNecessaryPreRunProperties(Project proj,
																											Map<Step, Map<Requirement, String>> variables);

		public abstract void run(Project proj, Map<Step, Map<Requirement, String>> variables);

		public void gracefulDeath(Project proj) {
			return;
		}

		public boolean hasRequirements(Project proj, Set<Step> stepSelections,
																	 Map<Step, Map<Requirement, String>> variables) {
			if (variables.get(this) == null) {
				return false;
			}
			for (Requirement[] group : getRequirements()) {
				boolean groupMet = false;
				for (Requirement req : group) {
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
		 * @return An array of {@link Requirement}s. At least one element of each subarray must be met
		 *         to satisfy the step pre-requisites - effectively this means elements of the first
		 *         array are AND'd together, while elements of the second array are OR'd.
		 */
		public Requirement[][] getRequirements() {
			return requirements;
		}

		public abstract boolean checkIfOutputExists(Map<Step, Map<Requirement, String>> variables);

		public void resetRun() {
			failed = false;
			failReasons.clear();
		}

		public abstract String getCommandLine(Project proj,
																					Map<Step, Map<Requirement, String>> variables);

		/**
		 * @return A {@link Collection} of the complete network of {@link Step)s related to this
		 *         {@code Step) - including this {@code Step), direct and transitive dependencies.
		 */
		public Collection<Step> getRelatedSteps() {
			return relatedSteps;
		}

		public Collection<Flag> getFlags() {
			return stepFlags;
		}

		public double getPriority() {
			return priority;
		}

		@Override
		public int compareTo(Step o) {
			// Preferably, just compare on priority. Otherwise, compare hash to prevent collisions
			int priorityCmp = Double.compare(getPriority(), o.getPriority());
			if (priorityCmp != 0)
				return priorityCmp;
			return Java6Helper.compare(hashCode(), o.hashCode());
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
			Step other = (Step) obj;
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


	public abstract static class Requirement {
		private final String description;
		private final RequirementInputType type;
		private final Object defaultValue;

		/**
		 * @param description
		 * @param type
		 */
		public Requirement(String description, RequirementInputType type) {
			this(description, type, null);
		}

		/**
		 * @param description
		 * @param type
		 * @param defaultValue
		 */
		public Requirement(String description, RequirementInputType type, Object defaultValue) {
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

		public double checkDoubleArgOrNeg1(String val) {
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

	public static class StepRequirement extends Requirement {
		private final Step requiredStep;

		public StepRequirement(Step requiredStep) {
			super(stepReqMessage(requiredStep), RequirementInputType.NONE);
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
			return "[" + requiredStep.getName()
						 + "] step must have been run already or must be selected";
		}


	}

	public static class FileRequirement extends Requirement {

		public FileRequirement(String description, String defaultValue) {
			super(description, RequirementInputType.FILE, defaultValue);
		}

		@Override
		public boolean checkRequirement(String arg, Set<Step> stepSelections,
																		Map<Step, Map<Requirement, String>> variables) {
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
			return !Files.exists(arg);
		}
	}

	public static class OptionalFileRequirement extends FileRequirement {
		public OptionalFileRequirement(String description, String defaultValue) {
			super(description, defaultValue);
		}

		@Override
		public boolean checkRequirement(String arg, Set<Step> stepSelections,
																		Map<Step, Map<Requirement, String>> variables) {
			return "".equals(arg) || Files.exists(arg);
		}
	}

	public static class BoolRequirement extends Requirement {

		public BoolRequirement(String description, boolean defaultValue) {
			super(description, RequirementInputType.BOOL, defaultValue);
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
			super(description, RequirementInputType.NUMBER, defaultValue);
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
			super(description, RequirementInputType.NUMBER, defaultValue);
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
			Maths.OPERATOR op = MarkerQC.findOperator(arg);
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
		NONE, FILE, DIR, STRING, NUMBER, BOOL, ENUM
	}

	public enum Flag {
		MEMORY, RUNTIME, MULTITHREADED
	}

	private static final String PROJ_PROP_UPDATE_STR = " org.genvisis.cnv.filesys.Project proj=";
	private static final String PLINK_SUBDIR = "plink/";
	private static final String PLINKROOT = "plink";
	final Project proj;
	private final SortedSet<Step> steps;
	private final Requirement numThreadsReq;
	Logger log;
	private final Launch launch;

	public GenvisisWorkflow(Project project, Launch launch) {
		if (project == null) {
			throw new IllegalArgumentException(this.getClass().getName()
																				 + " cannot be constructed with a null "
																				 + Project.class.getName());
		}
		proj = project;
		log = project.getLog();
		this.launch = launch;
		numThreadsReq = new PosIntRequirement("Number of Threads to Use", proj.NUM_THREADS.getValue());

		steps = Collections.unmodifiableSortedSet(generateSteps());
	}

	public void showDialogAndRun() {
		GenvisisWorkflowGUI gui;
		gui = new GenvisisWorkflowGUI(proj, launch, steps);
		if (!gui.getCancelled()) {
			gui.setModal(true);
			gui.setVisible(true);

			if (gui.getCancelled()) {
				return;
			}
		}
	}

	public Requirement getNumThreadsReq() {
		return numThreadsReq;
	}

	/**
	 * Helper class to minimize manual bookkeeping when instantiating steps. Each
	 * {@code generateXXXXStep} method should use the {@link #priority()} method to get its priorty,
	 * and call {@link #register(Step)} on the constructed step.
	 * <p>
	 * TODO: to reduce the risk of coding mistakes, convert the priority and register methods to
	 * intrinsic functions of the steps themselves
	 * </p>
	 */
	private class StepBuilder {
		private SortedSet<Step> buildSteps;
		private double p;

		public StepBuilder() {
			buildSteps = Sets.newTreeSet();
			p = 0.0;
		}

		/**
		 * @return All steps {@link #register(Step)}ed by this step builder thus far
		 */
		public SortedSet<Step> getSortedSteps() {
			return buildSteps;
		}

		/**
		 * @return The next step priority
		 */
		private double priority() {
			return ++p;
		}

		/**
		 * Register the given step in the list returned by {@link #getSortedSteps()}
		 */
		private Step register(Step s) {
			buildSteps.add(s);
			return s;
		}

		private Step generateMarkerPositionsStep() {
			final Requirement manifestReq = new FileRequirement("(Preferred) An Illumina Manifest file.",
																													"");
			final Requirement snpMapReq = new FileRequirement("An Illumina SNP_map file.",
																												proj.getLocationOfSNP_Map(false));
			return register(new Step("Create Marker Positions (if not already exists)", "",
															 new Requirement[][] {{manifestReq, snpMapReq}, {getNumThreadsReq()}},
															 EnumSet.noneOf(Flag.class), priority()) {

				@Override
				public void setNecessaryPreRunProperties(Project proj,
																								 Map<Step, Map<Requirement, String>> variables) {
					int numThreads = resolveThreads(variables.get(this).get(getNumThreadsReq()));
					maybeSetProjNumThreads(numThreads);
				}

				@Override
				public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
					proj.getLog().report("Generating marker positions file");
					String manifestFile = variables.get(this).get(manifestReq);
					String snpMapFile = variables.get(this).get(snpMapReq);
					if (Files.exists(manifestFile)) {
						proj.getLog().report("BLASTing probes from " + manifestFile);
						MarkerBlast.blastEm(proj, manifestFile, MarkerBlast.FILE_SEQUENCE_TYPE.MANIFEST_FILE,
																proj.NUM_THREADS.getValue());
					} else {
						proj.getLog().report("Generating marker positions file from " + snpMapFile);
						Markers.generateMarkerPositions(proj, snpMapFile);
					}
				}

				@Override
				public boolean checkIfOutputExists(Map<Step, Map<Requirement, String>> variables) {
					return Files.exists(proj.MARKER_POSITION_FILENAME.getValue(false, false));
				}

				@Override
				public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
					String manifestFile = variables.get(this).get(manifestReq);
					String snpMapFile = variables.get(this).get(snpMapReq);
					int numThreads = resolveThreads(variables.get(this).get(numThreadsReq));
					String projFile = proj.getPropertyFilename();
					if (Files.exists(manifestFile)) {
						List<String> command = ImmutableList.of(Files.getRunString(),
																										MarkerBlast.class.getName(),
																										"fileSeq=" + manifestFile,
																										"proj=" + proj.getPropertyFilename(),
																										PSF.Ext.NUM_THREADS_COMMAND + numThreads);
						return Joiner.on(" ").join(command);
					} else {
						return Files.getRunString() + " cnv.manage.Markers proj=" + projFile + " snps="
									 + snpMapFile;
					}
				}
			});
		}

		private Step generateIlluminaParseSamplesStep(final Step markerPositionsStep) {
			final Requirement markerPositionsStepReq = new StepRequirement(markerPositionsStep);

			final Requirement markerPositionsReq = new FileRequirement("Parsed markerPositions file must already exist.",
																																 proj.MARKER_POSITION_FILENAME.getValue(false,
																																																				false));

			return register(new Step("Parse Illumina Sample Files", "",
															 new Requirement[][] {{markerPositionsStepReq, markerPositionsReq},
																										{getNumThreadsReq()}},
															 EnumSet.of(Flag.MEMORY, Flag.RUNTIME, Flag.MEMORY),
															 priority()) {

				@Override
				public void setNecessaryPreRunProperties(Project proj,
																								 Map<Step, Map<Requirement, String>> variables) {
					String projFile = proj.MARKER_POSITION_FILENAME.getValue(false, false);
					String mkrFile = variables.get(this).get(markerPositionsReq);
					mkrFile = ext.verifyDirFormat(mkrFile);
					mkrFile = mkrFile.substring(0, mkrFile.length() - 1);
					if (!mkrFile.equals(projFile)) {
						proj.MARKER_POSITION_FILENAME.setValue(mkrFile);
					}
					int numThreads = resolveThreads(variables.get(this).get(getNumThreadsReq()));
					maybeSetProjNumThreads(numThreads);
				}

				@Override
				public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
					int numThreads = proj.NUM_THREADS.getValue();
					proj.getLog().report("Parsing sample files");
					int retCode = org.genvisis.cnv.manage.SourceFileParser.createFiles(proj, numThreads);
					switch (retCode) {
						case 0:
							setFailed("Operation failure, please check log for more information.");
							break;
						case 1:
						case 6:
						default:
							break;
					}
				}

				@Override
				public boolean checkIfOutputExists(Map<Step, Map<Requirement, String>> variables) {
					String sampleDirectory = proj.SAMPLE_DIRECTORY.getValue(false, false);
					boolean mkrSetFile = Files.exists(proj.MARKERSET_FILENAME.getValue(false, false));
					boolean returnValue = mkrSetFile;
					returnValue = returnValue && Files.exists(sampleDirectory);
					returnValue = returnValue && Files.list(sampleDirectory, Sample.SAMPLE_FILE_EXTENSION,
																									false).length > 0;
					returnValue = returnValue && proj.getSampleList() != null;
					returnValue = returnValue && proj.getSampleList().getSamples().length > 0;
					return returnValue;
				}

				@Override
				public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
					String projPropFile = proj.getPropertyFilename();
					StringBuilder kvCmd = new StringBuilder(Files.getRunString()).append(PROJ_PROP_UPDATE_STR)
																																			 .append(projPropFile);
					StringBuilder kvPairs = new StringBuilder();
					String projFile = proj.MARKER_POSITION_FILENAME.getValue(false, false);
					String mkrFile = variables.get(this).get(markerPositionsReq);
					mkrFile = ext.verifyDirFormat(mkrFile);
					mkrFile = mkrFile.substring(0, mkrFile.length() - 1);
					if (!mkrFile.equals(projFile)) {
						kvPairs.append(" MARKER_POSITION_FILENAME=").append(mkrFile);
					}
					int numThreads = resolveThreads(variables.get(this).get(getNumThreadsReq()));
					if (numThreads != proj.NUM_THREADS.getValue()) {
						kvPairs.append(" ").append(proj.NUM_THREADS.getName()).append("=").append(numThreads);
					}
					StringBuilder command = new StringBuilder();
					if (kvPairs.length() != 0) {
						command.append(kvCmd).append(kvPairs).append("\n");
					}
					command.append(Files.getRunString()).append(" cnv.manage.SourceFileParser proj=")
								 .append(projPropFile).append(" ").append(PSF.Ext.NUM_THREADS_COMMAND)
								 .append(numThreads);
					return command.toString();
				}

			});
		}

		private Step generateAffyParseSamplesStep() {

			final Requirement markerPositionsReq = new Requirement("markerPositions file must already exist.",
																														 RequirementInputType.FILE,
																														 proj.MARKER_POSITION_FILENAME.getValue(false,
																																																		false)) {

				@Override
				public boolean checkRequirement(String arg, Set<Step> stepSelections,
																				Map<Step, Map<Requirement, String>> variables) {
					return Files.exists(arg);
				}
			};

			return register(new Step("Parse Sample Files", "",
															 new Requirement[][] {{markerPositionsReq}, {getNumThreadsReq()}},
															 EnumSet.of(Flag.MEMORY, Flag.RUNTIME, Flag.MEMORY),
															 priority()) {

				@Override
				public void setNecessaryPreRunProperties(Project proj,
																								 Map<Step, Map<Requirement, String>> variables) {
					String projFile = proj.MARKER_POSITION_FILENAME.getValue(false, false);
					String mkrFile = variables.get(this).get(markerPositionsReq);
					mkrFile = ext.verifyDirFormat(mkrFile);
					mkrFile = mkrFile.substring(0, mkrFile.length() - 1);
					if (!mkrFile.equals(projFile)) {
						proj.MARKER_POSITION_FILENAME.setValue(mkrFile);
					}
					int numThreads = resolveThreads(variables.get(this).get(getNumThreadsReq()));
					maybeSetProjNumThreads(numThreads);
				}

				@Override
				public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
					int numThreads = resolveThreads(variables.get(this).get(getNumThreadsReq()));
					proj.getLog().report("Parsing sample files");
					int retCode = org.genvisis.cnv.manage.SourceFileParser.createFiles(proj, numThreads);
					switch (retCode) {
						case 0:
							setFailed("Operation failure, please check log for more information.");
							break;
						case 6:
						case 1:
						default:
							break;
					}
				}

				@Override
				public boolean checkIfOutputExists(Map<Step, Map<Requirement, String>> variables) {
					String sampleDirectory = proj.SAMPLE_DIRECTORY.getValue(false, false);
					boolean mkrSetFile = Files.exists(proj.MARKERSET_FILENAME.getValue(false, false));
					mkrSetFile = mkrSetFile && Files.exists(sampleDirectory);
					mkrSetFile = mkrSetFile
											 && Files.list(sampleDirectory, Sample.SAMPLE_FILE_EXTENSION,
																		 false).length > 0;
					mkrSetFile = mkrSetFile && proj.getSampleList() != null;
					mkrSetFile = mkrSetFile && proj.getSampleList().getSamples().length > 0;
					return mkrSetFile;
				}

				@Override
				public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
					String projPropFile = proj.getPropertyFilename();
					StringBuilder kvCmd = new StringBuilder(Files.getRunString()).append(PROJ_PROP_UPDATE_STR)
																																			 .append(projPropFile);
					StringBuilder kvPairs = new StringBuilder();
					String projFile = proj.MARKER_POSITION_FILENAME.getValue(false, false);
					String mkrFile = variables.get(this).get(markerPositionsReq);
					mkrFile = ext.verifyDirFormat(mkrFile);
					mkrFile = mkrFile.substring(0, mkrFile.length() - 1);
					if (!mkrFile.equals(projFile)) {
						kvPairs.append(" MARKER_POSITION_FILENAME=").append(mkrFile);
					}

					int numThreads = resolveThreads(variables.get(this).get(getNumThreadsReq()));
					if (numThreads != proj.NUM_THREADS.getValue()) {
						kvPairs.append(" ").append(proj.NUM_THREADS.getName()).append("=").append(numThreads);
					}
					StringBuilder command = new StringBuilder();
					if (kvPairs.length() != 0) {
						command.append(kvCmd).append(kvPairs).append("\n");
					}
					command.append(Files.getRunString()).append(" cnv.manage.SourceFileParser proj=")
								 .append(projPropFile).append(" ").append(PSF.Ext.NUM_THREADS_COMMAND)
								 .append(numThreads);
					return command.toString();
				}

			});
		}

		private Step generateCreateSampleDataStep(final Step parseSamplesStep) {
			final Requirement parseSamplesStepReq = new StepRequirement(parseSamplesStep);

			final Requirement createMinimalSampleDataReq = new BoolRequirement("Create a minimal SampleData.txt file from sample files",
																																				 true);

			final String pedPreset = proj.PEDIGREE_FILENAME.getValue();

			final Requirement pedigreeReq = new Requirement("Either a Pedigree.dat file, or any file with a header containing all of the following elements (in any order):  \""
																											+ ArrayUtils.toStr(MitoPipeline.PED_INPUT,
																																				 ", ")
																											+ "\"", RequirementInputType.FILE,
																											pedPreset) {

				@Override
				public boolean checkRequirement(String arg, Set<Step> stepSelections,
																				Map<Step, Map<Requirement, String>> variables) {
					return Files.exists(arg);
				}
			};

			// check for SampleMap only if we haven't found a pedigree
			final String sampMapPreset = Files.exists(pedPreset) ? null : getLocationOfSampleMap();

			final Requirement sampMapReq = new Requirement("A Sample_Map.csv file, with at least two columns having headers \""
																										 + MitoPipeline.SAMPLEMAP_INPUT[1] + "\" and \""
																										 + MitoPipeline.SAMPLEMAP_INPUT[2] + "\"",
																										 RequirementInputType.FILE, sampMapPreset) {

				@Override
				public boolean checkRequirement(String arg, Set<Step> stepSelections,
																				Map<Step, Map<Requirement, String>> variables) {
					return Files.exists(arg);
				}
			};
			return register(new Step("Create SampleData.txt File", "",
															 new Requirement[][] {{parseSamplesStepReq},
																										{createMinimalSampleDataReq, pedigreeReq,
																										 sampMapReq}},
															 EnumSet.noneOf(Flag.class), priority()) {

				@Override
				public void setNecessaryPreRunProperties(Project proj,
																								 Map<Step, Map<Requirement, String>> variables) {
					// Nothing to do
				}

				@Override
				public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
					Boolean minimal = Boolean.parseBoolean(variables.get(this)
																													.get(createMinimalSampleDataReq));
					String pedFile = minimal ? null : variables.get(this).get(pedigreeReq);
					String sampleMapCsv = minimal ? null : variables.get(this).get(sampMapReq);

					proj.getLog().report("Creating SampleData.txt");
					try {
						int retStat = SampleData.createSampleData(pedFile, sampleMapCsv, proj);
						if (retStat == -1) {
							setFailed("SampleData already exists - please delete and try again.");
							return;
						}
					} catch (Elision e) {
						setFailed(e.getMessage());
						return;
					}
				}

				@Override
				public boolean checkIfOutputExists(Map<Step, Map<Requirement, String>> variables) {
					return Files.exists(proj.SAMPLE_DATA_FILENAME.getValue(false, false));
				}

				@Override
				public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
					String projPropFile = proj.getPropertyFilename();
					Boolean minimal = Boolean.parseBoolean(variables.get(this)
																													.get(createMinimalSampleDataReq));
					String pedFile = minimal ? "" : variables.get(this).get(pedigreeReq);
					String sampleMapCsv = minimal ? "" : variables.get(this).get(sampMapReq);
					StringBuilder cmd = new StringBuilder();
					cmd.append(Files.getRunString()).append(" cnv.var.SampleData proj=").append(projPropFile);
					if (!"".equals(pedFile)) {
						cmd.append(" ped=").append(pedFile);
					}
					if (!"".equals(sampleMapCsv)) {
						cmd.append(" sampleMap=").append(sampleMapCsv);
					}
					return cmd.toString();
				}

			});
		}

		private Step generateTransposeStep(final Step parseSamplesStep) {
			return register(new Step("Transpose Data into Marker-Dominant Files", "",
															 new Requirement[][] {{new StepRequirement(parseSamplesStep)}},
															 EnumSet.of(Flag.MEMORY), priority()) {

				@Override
				public void setNecessaryPreRunProperties(Project proj,
																								 Map<Step, Map<Requirement, String>> variables) {
					// Nothing to do here
				}

				@Override
				public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
					proj.getLog().report("Transposing data");
					TransposeData.transposeData(proj, 2000000000, false); // compact if no LRR was provided
				}

				@Override
				public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
					String projPropFile = proj.getPropertyFilename();
					StringBuilder cmd = new StringBuilder();
					return cmd.append(Files.getRunString())
										.append(" cnv.manage.TransposeData -transpose proj="
														+ projPropFile + " max=" + 2000000000)
										.toString();
				}

				@Override
				public boolean checkIfOutputExists(Map<Step, Map<Requirement, String>> variables) {
					return Files.list(proj.MARKER_DATA_DIRECTORY.getValue(false, false),
														MarkerData.MARKER_DATA_FILE_EXTENSION, false).length > 0;
				}
			});
		}

		private Step generateGCModelStep() {
			final Requirement gcBaseReq = new FileRequirement("A GC Base file must exist.",
																												Resources.genome(proj.GENOME_BUILD_VERSION.getValue(),
																																				 proj.getLog())
																																 .getModelBase().get());
			final Requirement gcModelOutputReq = new OutputFileRequirement("GCModel output file must be specified.",
																																		 proj.GC_MODEL_FILENAME.getValue());

			return register(new Step("Compute GCMODEL File", "",
															 new Requirement[][] {{gcBaseReq}, {gcModelOutputReq}},
															 EnumSet.noneOf(Flag.class), priority()) {

				@Override
				public void setNecessaryPreRunProperties(Project proj,
																								 Map<Step, Map<Requirement, String>> variables) {
					String setGCOutputFile = proj.GC_MODEL_FILENAME.getValue();
					String gcOutputFile = variables.get(this).get(gcModelOutputReq);
					if (!ext.verifyDirFormat(setGCOutputFile).equals(gcOutputFile)) {
						proj.GC_MODEL_FILENAME.setValue(gcOutputFile);
					}
				}

				@Override
				public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
					String gcBaseFile = variables.get(this).get(gcBaseReq);
					String gcOutputFile = variables.get(this).get(gcModelOutputReq);
					org.genvisis.cnv.analysis.PennCNV.gcModel(proj, gcBaseFile, gcOutputFile, 100);
				}

				@Override
				public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
					String kvCmd = "";

					String setGCOutputFile = proj.GC_MODEL_FILENAME.getValue();
					String gcOutputFile = variables.get(this).get(gcModelOutputReq);
					if (!ext.verifyDirFormat(setGCOutputFile).equals(gcOutputFile)) {
						kvCmd += " GC_MODEL_FILENAME=" + gcOutputFile;
					}

					String projPropFile = proj.getPropertyFilename();
					StringBuilder cmd = new StringBuilder();
					if (kvCmd.length() > 0) {
						cmd.append(Files.getRunString()).append(PROJ_PROP_UPDATE_STR + projPropFile)
							 .append(kvCmd)
							 .append("\n");
					}
					String gcBaseFile = variables.get(this).get(gcBaseReq);
					return cmd.append(Files.getRunString())
										.append(" cnv.analysis.PennCNV proj=" + proj.getPropertyFilename() + " log="
														+ proj.getLog().getFilename() + " gc5base=" + gcBaseFile)
										.toString();
				}

				@Override
				public boolean checkIfOutputExists(Map<Step, Map<Requirement, String>> variables) {
					String gcOutputFile = variables.get(this).get(gcModelOutputReq);
					return Files.exists(gcOutputFile);
				}

			});
		}

		private Step generateSampleQCStep(final Step parseSamplesStep) {
			final Requirement parseSamplesStepReq = new StepRequirement(parseSamplesStep);

			return register(new Step("Run Sample QC Metrics", "",
															 new Requirement[][] {{parseSamplesStepReq}, {getNumThreadsReq()},},
															 EnumSet.noneOf(Flag.class), priority()) {

				@Override
				public void setNecessaryPreRunProperties(Project proj,
																								 Map<Step, Map<Requirement, String>> variables) {
					int numThreads = resolveThreads(variables.get(this).get(getNumThreadsReq()));
					maybeSetProjNumThreads(numThreads);
				}

				@Override
				public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
					proj.getLog().report("Running LrrSd");
					int numThreads = proj.NUM_THREADS.getValue();
					LrrSd.init(proj, null, null, numThreads, false);
				}

				@Override
				public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
					int numThreads = resolveThreads(variables.get(this).get(getNumThreadsReq()));
					String projPropFile = proj.getPropertyFilename();
					StringBuilder cmd = new StringBuilder();
					cmd.append(Files.getRunString()).append(" cnv.qc.LrrSd").append(" proj=")
						 .append(projPropFile).append(" ").append(PSF.Ext.NUM_THREADS_COMMAND)
						 .append(numThreads)
						 .append(" projectMarkers=TRUE");
					return cmd.toString();
				}

				@Override
				public boolean checkIfOutputExists(Map<Step, Map<Requirement, String>> variables) {
					return Files.exists(proj.SAMPLE_QC_FILENAME.getValue(false, false));
				}
			});
		}

		private Step generateMarkerQCStep(final Step parseSamplesStep) {
			final Requirement parseSamplesStepReq = new StepRequirement(parseSamplesStep);

			final Requirement exportAllReq = new OptionalBoolRequirement("Export all markers in project.",
																																	 true);

			final Requirement targetMarkersReq = new FileRequirement("A targetMarkers files listing the markers to QC.",
																															 proj.TARGET_MARKERS_FILENAMES.getValue()[0]);

			return register(new Step("Run Marker QC Metrics", "",
															 new Requirement[][] {{parseSamplesStepReq},
																										{exportAllReq, targetMarkersReq},
																										{getNumThreadsReq()}},
															 EnumSet.noneOf(Flag.class), priority()) {

				@Override
				public void setNecessaryPreRunProperties(Project proj,
																								 Map<Step, Map<Requirement, String>> variables) {
					// Nothing to do here
				}

				@Override
				public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
					boolean allMarkers = Boolean.parseBoolean(variables.get(this).get(exportAllReq));
					String tgtFile = allMarkers ? null : variables.get(this).get(targetMarkersReq);
					int numThreads = resolveThreads(variables.get(this).get(getNumThreadsReq()));
					MarkerMetrics.fullQC(proj, null, tgtFile, true, numThreads);
				}

				@Override
				public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
					boolean allMarkers = Boolean.parseBoolean(variables.get(this).get(exportAllReq));
					String tgtFile = variables.get(this).get(targetMarkersReq);
					int numThreads = resolveThreads(variables.get(this).get(getNumThreadsReq()));
					return Files.getRunString() + " cnv.qc.MarkerMetrics -fullQC" + " proj="
								 + proj.getPropertyFilename() + (allMarkers ? "" : " markers=" + tgtFile) + " "
								 + PSF.Ext.NUM_THREADS_COMMAND + numThreads;
				}

				@Override
				public boolean checkIfOutputExists(Map<Step, Map<Requirement, String>> variables) {
					String markerMetricsFile = proj.MARKER_METRICS_FILENAME.getValue();
					return Files.exists(markerMetricsFile);
				}

			});
		}

		private Step generateSexChecksStep(final Step parseSamplesStep,
																			 final Step sampleDataStep, final Step transposeStep,
																			 final Step sampleQCStep) {
			final Requirement parseSamplesStepReq = new StepRequirement(parseSamplesStep);
			final Requirement sampleDataStepReq = new StepRequirement(sampleDataStep);
			final Requirement transposeStepReq = new StepRequirement(transposeStep);
			final Requirement sampleQCStepReq = new StepRequirement(sampleQCStep);
			final Requirement addToSampleDataReq = new OptionalBoolRequirement("Add Estimated Sex to Sample Data",
																																				 true);

			final Requirement noCrossHybeReq = new BoolRequirement("Use only X and Y chromosome R values to identify sex discriminating markers",
																														 false);

			final Requirement oneHittersReq = new FileRequirement("List of markers that do not cross hybridize",
																														MarkerBlastQC.defaultOneHitWondersFilename(proj.BLAST_ANNOTATION_FILENAME.getValue()));

			final Requirement blastVCFReq = new FileRequirement("BLAST annotation VCF to generate list of markers that do not cross hybridize from",
																													proj.BLAST_ANNOTATION_FILENAME.getValue());

			return register(new Step("Run Sex Checks", "",
															 new Requirement[][] {{parseSamplesStepReq}, {sampleDataStepReq},
																										{transposeStepReq}, {sampleQCStepReq},
																										{addToSampleDataReq},
																										{noCrossHybeReq, oneHittersReq, blastVCFReq}},
															 EnumSet.noneOf(Flag.class), priority()) {

				@Override
				public void setNecessaryPreRunProperties(Project proj,
																								 Map<Step, Map<Requirement, String>> variables) {
					// Nothing to do here
				}

				@Override
				public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
					proj.getLog().report("Running SexCheck");
					boolean addToSampleData = Boolean.parseBoolean(variables.get(this)
																																	.get(addToSampleDataReq));
					String discriminatingMarkersFile;
					if (Boolean.parseBoolean(variables.get(this).get(noCrossHybeReq))) {
						discriminatingMarkersFile = null;
					} else {
						discriminatingMarkersFile = variables.get(this).get(oneHittersReq);
						if (!Files.exists(discriminatingMarkersFile)) {
							MarkerBlastQC.getOneHitWonders(proj, variables.get(this).get(blastVCFReq),
																						 discriminatingMarkersFile, 0.8, proj.getLog());
						}
					}
					org.genvisis.cnv.qc.SexChecks.sexCheck(proj, addToSampleData, discriminatingMarkersFile);
				}

				@Override
				public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
					String projPropFile = proj.getPropertyFilename();
					StringBuilder cmd = new StringBuilder();
					boolean addToSampleData = Boolean.parseBoolean(variables.get(this)
																																	.get(addToSampleDataReq));
					String discriminatingMarkersFile;
					if (Boolean.parseBoolean(variables.get(this).get(noCrossHybeReq))) {
						discriminatingMarkersFile = null;
					} else {
						discriminatingMarkersFile = variables.get(this).get(oneHittersReq);
						if (!Files.exists(discriminatingMarkersFile)) {
							cmd.append(Files.getRunString())
								 .append(" cnv.qc.MarkerBlastQC proj=" + projPropFile + " blastVCF="
												 + variables.get(this).get(blastVCFReq))
								 .append("\n");
						}
					}
					return cmd.append(Files.getRunString())
										.append(" cnv.qc.SexChecks -check proj=" + projPropFile).toString()
								 + (discriminatingMarkersFile == null ? ""
																											: " useMarkers=" + discriminatingMarkersFile)
								 + (addToSampleData ? "" : " -skipSampleData");
				}

				@Override
				public boolean checkIfOutputExists(Map<Step, Map<Requirement, String>> variables) {
					return Files.exists(proj.SEXCHECK_RESULTS_FILENAME.getValue());
				}

			});
		}

		private Step generatePlinkExportStep(final Step parseSamplesStep) {
			final Requirement parseSamplesStepReq = new StepRequirement(parseSamplesStep);
			final Requirement pedigreeRequirement = new FileRequirement("A pedigree.dat file must exist.",
																																	proj.PEDIGREE_FILENAME.getValue(false,
																																																	false));
			final Requirement createPedigreeRequirement = new BoolRequirement("Create a minimal pedigree.dat file [will pull information from SexChecks step results].",
																																				false);

			return register(new Step("Create PLINK Files", "",
															 new Requirement[][] {{parseSamplesStepReq},
																										{pedigreeRequirement,
																										 createPedigreeRequirement}},
															 EnumSet.of(Flag.MEMORY), priority()) {

				@Override
				public void setNecessaryPreRunProperties(Project proj,
																								 Map<Step, Map<Requirement, String>> variables) {
					if (!Boolean.parseBoolean(variables.get(this).get(createPedigreeRequirement))) {
						String projPedFile = proj.PEDIGREE_FILENAME.getValue(false, false);
						String pedFile = variables.get(this).get(pedigreeRequirement);
						if (!pedFile.equals(projPedFile)) {
							proj.PEDIGREE_FILENAME.setValue(pedFile);
						}
					}
				}

				@Override
				public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
					if (Boolean.parseBoolean(variables.get(this).get(createPedigreeRequirement))) {
						proj.getLog().report("Creating Pedigree File");
						Pedigree.build(proj, null, null, false);
					}
					if (!Files.exists(proj.PEDIGREE_FILENAME.getValue())) {
						setFailed("Creation of Pedigree file in [Create/Run PLINK Files] step failed.");
						return;
					}

					proj.getLog().report("Running PLINK");

					boolean create = PlinkData.saveGenvisisToPlinkBedSet(proj, PLINK_SUBDIR + PLINKROOT,
																															 null, null, -1);
					if (!create) {
						setFailed("Creation of initial PLINK files failed.");
					}
					proj.PLINK_DIR_FILEROOTS.addValue(proj.PROJECT_DIRECTORY.getValue() + PLINK_SUBDIR
																						+ PLINKROOT);
				}

				@Override
				public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
					String kvCmd = "";

					if (!Boolean.parseBoolean(variables.get(this).get(createPedigreeRequirement))) {
						String projPedFile = proj.PEDIGREE_FILENAME.getValue(false, false);
						String pedFile = variables.get(this).get(pedigreeRequirement);
						if (!pedFile.equals(projPedFile)) {
							kvCmd += " PEDIGREE_FILENAME=" + pedFile;
						}
					}

					String projPropFile = proj.getPropertyFilename();
					StringBuilder cmd = new StringBuilder();
					if (kvCmd.length() > 0) {
						cmd.append(Files.getRunString()).append(PROJ_PROP_UPDATE_STR).append(projPropFile)
							 .append(kvCmd).append("\n");
					}
					if (Boolean.parseBoolean(variables.get(this).get(createPedigreeRequirement))) {
						cmd.append(Files.getRunString()).append(" cnv.filesys.Pedigree proj=")
							 .append(projPropFile).append("\n");
					}
					cmd.append(Files.getRunString())
						 .append(" cnv.manage.PlinkData -genvisisToBed plinkdata=" + PLINK_SUBDIR + PLINKROOT
										 + " gcthreshold=-1 proj=")
						 .append(proj.getPropertyFilename());
					return cmd.toString();
				}

				@Override
				public boolean checkIfOutputExists(Map<Step, Map<Requirement, String>> variables) {
					boolean plinkFilesExist = Files.checkAllFiles(getPlinkDir(),
																												PSF.Plink.getPlinkBedBimFamSet(PLINKROOT),
																												false, log);
					boolean pedGenerated = Boolean.parseBoolean(variables.get(this)
																															 .get(createPedigreeRequirement));
					boolean pedCheck = pedGenerated ? Files.exists(proj.PEDIGREE_FILENAME.getValue()) : true;
					return plinkFilesExist && pedCheck;
				}

			});
		}

		private Step generateGwasQCStep(Step plinkExportStep) {
			final Requirement plinkExportStepReq = new StepRequirement(plinkExportStep);
			String defaultCallrate;
			switch (proj.getArrayType()) {
				case AFFY_GW6:
				case AFFY_GW6_CN:
					defaultCallrate = MarkerQC.DEFAULT_AFFY_CALLRATE_THRESHOLD;
					break;
				case ILLUMINA:
					defaultCallrate = MarkerQC.DEFAULT_ILLUMINA_CALLRATE_THRESHOLD;
					break;
				default:
					throw new IllegalArgumentException("Invalid " + proj.getArrayType().getClass().getName()
																						 + ": " + proj.getArrayType().toString());
			}
			final Requirement callrateReq = new ThresholdRequirement(QC_METRIC.CALLRATE.getUserDescription(),
																															 defaultCallrate);

			return register(new Step("Run GWAS QC", "",
															 new Requirement[][] {{plinkExportStepReq}, {callrateReq}},
															 EnumSet.noneOf(Flag.class), priority()) {

				@Override
				public void setNecessaryPreRunProperties(Project proj,
																								 Map<Step, Map<Requirement, String>> variables) {
					// not needed for step
				}

				@Override
				public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
					String dir = getPlinkDir();
					Map<QC_METRIC, String> markerQCThresholds = Maps.newEnumMap(RelationAncestryQc.DEFAULT_QC_METRIC_THRESHOLDS);
					markerQCThresholds.put(QC_METRIC.CALLRATE, variables.get(this).get(callrateReq));
					new RelationAncestryQc(dir, PLINKROOT, markerQCThresholds, log).run(false);
					RelationAncestryQc.fullGamut(dir, PLINKROOT, false, proj.getLog());
					if (new File(dir + Qc.QC_SUBDIR + RelationAncestryQc.GENOME_DIR + PLINKROOT
											 + ".genome").exists()) {
						proj.GENOME_CLUSTER_FILENAME.setValue(dir + Qc.QC_SUBDIR + RelationAncestryQc.GENOME_DIR
																									+ PLINKROOT + ".genome");
						proj.saveProperties();
					}
					new PlinkMendelianChecker(proj).run();
				}

				@Override
				public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
					Map<Requirement, String> stepVars = variables.get(this);

					String dir = getPlinkDir();
					Map<QC_METRIC, String> markerQCThresholds = Maps.newEnumMap(RelationAncestryQc.DEFAULT_QC_METRIC_THRESHOLDS);
					markerQCThresholds.put(QC_METRIC.CALLRATE, variables.get(this).get(callrateReq));


					List<String> commandChunks = Lists.newArrayList();
					commandChunks.add(Files.getRunString());
					commandChunks.add(RelationAncestryQc.class.getName());
					commandChunks.add(CLI.formCmdLineArg(CLI.ARG_INDIR, getPlinkDir()));
					commandChunks.add(CLI.formCmdLineArg(CLI.ARG_PLINKROOT, PLINKROOT));
					commandChunks.add(CLI.formCmdLineArg(RelationAncestryQc.ARGS_KEEPGENOME, "false"));
					commandChunks.add(CLI.formCmdLineArg(QC_METRIC.CALLRATE.getKey(),
																							 stepVars.get(callrateReq)));
					commandChunks.add("\n" + Files.getRunString());
					commandChunks.add(PROJ_PROP_UPDATE_STR + proj.getPropertyFilename());
					commandChunks.add(proj.GENOME_CLUSTER_FILENAME.getName() + "=" + dir + Qc.QC_SUBDIR
														+ RelationAncestryQc.GENOME_DIR + PLINKROOT + ".genome");
					commandChunks.add("\n" + Files.getRunString());
					commandChunks.add(PlinkMendelianChecker.class.getName());
					commandChunks.add("proj=" + proj.getPropertyFilename());
					return Joiner.on(' ').join(commandChunks);
				}

				@Override
				public boolean checkIfOutputExists(Map<Step, Map<Requirement, String>> variables) {
					String dir = getPlinkDir();
					for (int i = 0; i < org.genvisis.gwas.RelationAncestryQc.FOLDERS_CREATED.length; i++) {
						for (int j = 0; j < org.genvisis.gwas.RelationAncestryQc.FILES_CREATED[i].length; j++) {
							if (!Files.exists(dir + org.genvisis.gwas.RelationAncestryQc.FOLDERS_CREATED[i]
																+ org.genvisis.gwas.RelationAncestryQc.FILES_CREATED[i][j])) {
								return false;
							}
						}
					}
					return Files.checkAllFiles(PlinkMendelianChecker.parseOutputDirectory(proj),
																		 PlinkMendelianChecker.OUTPUTS, false, log);
				}
			});
		}

		private Step generateAncestryStep(final Step gwasQCStep) {
			final Requirement gwasQCStepReq = new StepRequirement(gwasQCStep);
			final Requirement putativeWhitesReq = new FileRequirement("File with FID/IID pairs of putative white samples",
																																"");

			return register(new Step("Run Ancestry Checks", "",
															 new Requirement[][] {{gwasQCStepReq}, {putativeWhitesReq}},
															 EnumSet.noneOf(Flag.class), priority()) {

				@Override
				public void setNecessaryPreRunProperties(Project proj,
																								 Map<Step, Map<Requirement, String>> variables) {
					// not needed for step
				}

				@Override
				public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
					String putativeWhites = variables.get(this).get(putativeWhitesReq);
					String hapMapPlinkRoot = Resources.hapMap(log).getUnambiguousHapMapFounders().get();
					int hapMapDotIndex = hapMapPlinkRoot.lastIndexOf('.');
					if (hapMapDotIndex > 0 && PSF.Plink.getPlinkBedBimFamSet("")
																						 .contains(hapMapPlinkRoot.substring(hapMapDotIndex))) {
						hapMapPlinkRoot = hapMapPlinkRoot.substring(0, hapMapDotIndex);
					}
					String ancestryDir = getAncestryDir();
					Ancestry.runPipeline(ancestryDir, putativeWhites, hapMapPlinkRoot, proj,
															 new Logger(ancestryDir + "ancestry.log"));
				}

				@Override
				public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
					String putativeWhites = variables.get(this).get(putativeWhitesReq);
					String hapMapPlinkRoot = Resources.hapMap(log).getUnambiguousHapMapFounders().get();
					int hapMapDotIndex = hapMapPlinkRoot.lastIndexOf('.');
					if (hapMapDotIndex > 0 && PSF.Plink.getPlinkBedBimFamSet("")
																						 .contains(hapMapPlinkRoot.substring(hapMapDotIndex))) {
						hapMapPlinkRoot = hapMapPlinkRoot.substring(0, hapMapDotIndex);
					}

					String ancestryDir = getAncestryDir();
					String command = Files.getRunString() + " gwas.Ancestry -runPipeline dir=" + ancestryDir;
					command += " putativeWhites=" + putativeWhites;
					command += " proj=" + proj.getPropertyFilename();
					command += " hapMapPlinkRoot=" + hapMapPlinkRoot;
					command += " log=" + ancestryDir + "ancestry.log";
					return command;
				}

				@Override
				public boolean checkIfOutputExists(Map<Step, Map<Requirement, String>> variables) {
					String ancestryDir = getAncestryDir();
					return Files.exists(ancestryDir + Ancestry.RACE_FREQS_FILENAME)
								 && Files.exists(ancestryDir + Ancestry.RACE_IMPUTATIONAS_FILENAME);
				}
			});
		}

		private Step generateFurtherAnalysisQCStep(Step plinkExportStep,
																							 Step gwasQCStep,
																							 Step ancestryStep) {
			final Requirement plinkExportStepReq = new StepRequirement(plinkExportStep);
			final Requirement gwasQCStepReq = new StepRequirement(gwasQCStep);
			final Requirement ancestryStepReq = new StepRequirement(ancestryStep);
			final Requirement unrelatedsFileReq = new FileRequirement("File with list of unrelated FID/IID pairs to use for marker QC",
																																"");
			final Requirement europeansFilesReq = new FileRequirement("File with list of European samples to use for Hardy-Weinberg equilibrium tests",
																																"");

			List<Requirement[]> requirementsList = Lists.newArrayList();
			requirementsList.add(new Requirement[] {plinkExportStepReq});
			requirementsList.add(new Requirement[] {gwasQCStepReq, unrelatedsFileReq});
			requirementsList.add(new Requirement[] {ancestryStepReq, europeansFilesReq});
			final Map<QC_METRIC, Requirement> metricRequirements = Maps.newEnumMap(QC_METRIC.class);
			for (QC_METRIC metric : QC_METRIC.values()) {
				Map<QC_METRIC, String> defaultThresholds = FurtherAnalysisQc.getDefaultMarkerQCThresholds(proj.getArrayType());
				String defaultVal = defaultThresholds.get(metric);
				final Requirement metricReq = new ThresholdRequirement(metric.getUserDescription(),
																															 defaultVal);
				requirementsList.add(new Requirement[] {metricReq});
				metricRequirements.put(metric, metricReq);
			}

			return register(new Step("Run Further Analysis QC", "",
															 requirementsList.toArray(new Requirement[requirementsList.size()][]),
															 EnumSet.noneOf(Flag.class), priority()) {

				@Override
				public void setNecessaryPreRunProperties(Project proj,
																								 Map<Step, Map<Requirement, String>> variables) {
					// not needed for step
				}

				@Override
				public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
					Map<Requirement, String> stepVars = variables.get(this);

					String unrelatedsFile = resolveUnrelatedsFile(stepVars);

					String europeansFile = resolveEuropeansFile(stepVars);

					Map<QC_METRIC, String> markerQCThresholds = Maps.newEnumMap(QC_METRIC.class);
					for (QC_METRIC metric : QC_METRIC.values()) {
						Requirement req = metricRequirements.get(metric);
						markerQCThresholds.put(metric, stepVars.get(req));
					}
					new FurtherAnalysisQc(getPlinkDir(), PLINKROOT, markerQCThresholds, unrelatedsFile,
																europeansFile, log).runFurtherAnalysisQC();
				}

				@Override
				public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
					Map<Requirement, String> stepVars = variables.get(this);

					String unrelatedsFile = resolveUnrelatedsFile(stepVars);

					String europeansFile = resolveEuropeansFile(stepVars);

					List<String> commandChunks = Lists.newArrayList();
					commandChunks.add(Files.getRunString());
					commandChunks.add(FurtherAnalysisQc.class.getName());
					commandChunks.add(CLI.formCmdLineArg(FurtherAnalysisQc.ARG_UNRELATEDS, unrelatedsFile));
					commandChunks.add(CLI.formCmdLineArg(FurtherAnalysisQc.ARG_EUROPEANS, europeansFile));
					commandChunks.add(CLI.formCmdLineArg(CLI.ARG_INDIR, getPlinkDir()));
					commandChunks.add(CLI.formCmdLineArg(CLI.ARG_PLINKROOT, PLINKROOT));
					for (QC_METRIC metric : QC_METRIC.values()) {
						Requirement req = metricRequirements.get(metric);
						commandChunks.add(CLI.formCmdLineArg(metric.getKey(), stepVars.get(req)));
					}
					return Joiner.on(' ').join(commandChunks);
				}

				@Override
				public boolean checkIfOutputExists(Map<Step, Map<Requirement, String>> variables) {
					String dir = getPlinkDir() + Qc.QC_SUBDIR + FurtherAnalysisQc.FURTHER_ANALYSIS_DIR;
					String qcdPlinkroot = PLINKROOT + FurtherAnalysisQc.FURTHER_ANALYSIS_QC_PLINK_SUFFIX;
					return PSF.Plink.bedBimFamExist(dir + qcdPlinkroot);
				}

				private String resolveUnrelatedsFile(Map<Requirement, String> stepVars) {
					String unrelatedsFile = stepVars.get(unrelatedsFileReq);
					if (!Files.exists(unrelatedsFile)) {
						unrelatedsFile = getAncestryDir() + RelationAncestryQc.UNRELATEDS_FILENAME;
					}
					return unrelatedsFile;
				}

				private String resolveEuropeansFile(Map<Requirement, String> stepVars) {
					String europeansFile = stepVars.get(europeansFilesReq);
					if (europeansFile == null || "".equals(europeansFile)) {
						String raceImputationFilename = getAncestryDir() + Ancestry.RACE_IMPUTATIONAS_FILENAME;
						europeansFile = PCImputeRace.formRaceListFilename(RACE.WHITE, raceImputationFilename);
					}
					return europeansFile;
				}

			});
		}

		private Step generateMosaicArmsStep(final Step parseSamplesStep) {
			final Requirement parseSamplesStepReq = new StepRequirement(parseSamplesStep);
			return register(new Step("Create Mosaic Arms File", "",
															 new Requirement[][] {{parseSamplesStepReq}, {getNumThreadsReq()}},
															 EnumSet.noneOf(Flag.class), priority()) {

				@Override
				public void setNecessaryPreRunProperties(Project proj,
																								 Map<Step, Map<Requirement, String>> variables) {

					int numThreads = resolveThreads(variables.get(this).get(getNumThreadsReq()));
					maybeSetProjNumThreads(numThreads);
				}

				@Override
				public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
					Mosaicism.findOutliers(proj);
				}

				@Override
				public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
					String kvCmd = "";


					int numThreads = resolveThreads(variables.get(this).get(getNumThreadsReq()));
					if (numThreads != proj.NUM_THREADS.getValue()) {
						kvCmd += " " + proj.NUM_THREADS.getName() + "=" + numThreads;
					}

					String projPropFile = proj.getPropertyFilename();
					StringBuilder cmd = new StringBuilder();
					if (kvCmd.length() > 0) {
						cmd.append(Files.getRunString()).append(PROJ_PROP_UPDATE_STR + projPropFile)
							 .append(kvCmd)
							 .append("\n");
					}
					return cmd.append(Files.getRunString())
										.append(" cnv.analysis.Mosaicism proj=" + proj.getPropertyFilename())
										.toString();
				}

				@Override
				public boolean checkIfOutputExists(Map<Step, Map<Requirement, String>> variables) {
					return Files.exists(proj.RESULTS_DIRECTORY.getValue(false, false) + "Mosaicism.xln");
				}
			});
		}

		private Step generateAnnotateSampleDataStep(final Step sampleQCStep,
																								final Step createSampleDataStep,
																								final Step gwasQCStep) {
			final Requirement sampleQCStepReq = new StepRequirement(sampleQCStep);
			final Requirement createSampleDataStepReq = new StepRequirement(createSampleDataStep);
			final Requirement skipIDingDuplicatesReq = new BoolRequirement("Skip identifying duplicates",
																																		 false);
			final Requirement gwasQCStepReq = new StepRequirement(gwasQCStep);
			final Requirement notGcCorrectedLrrSdReq = new BoolRequirement("Do not use GC corrected LRR SD?",
																																		 false);
			final Requirement gcCorrectedLrrSdReq = new Requirement("GC Corrected LRR SD must exist in Sample QC File",
																															RequirementInputType.NONE) {

				@Override
				public boolean checkRequirement(String arg, Set<Step> stepSelections,
																				Map<Step, Map<Requirement, String>> variables) {
					String sampleQCFile = proj.SAMPLE_QC_FILENAME.getValue();
					return Files.exists(sampleQCFile)
								 && ext.indexOfStr("LRR_SD_Post_Correction",
																	 Files.getHeaderOfFile(sampleQCFile, proj.getLog())) != -1;
				}

			};
			final Requirement lrrSdThresholdReq = new DoubleRequirement("LRR SD Threshold",
																																	proj.LRRSD_CUTOFF.getValue(),
																																	proj.LRRSD_CUTOFF.getMinValue(),
																																	proj.LRRSD_CUTOFF.getMaxValue());

			final Requirement callrateThresholdReq = new DoubleRequirement("Callrate Threshold",
																																		 proj.SAMPLE_CALLRATE_THRESHOLD.getValue(),
																																		 proj.SAMPLE_CALLRATE_THRESHOLD.getMinValue(),
																																		 proj.SAMPLE_CALLRATE_THRESHOLD.getMaxValue());
			final Requirement numQReq = new PosIntRequirement("Number of Quantiles to Generate", 10);
			final Requirement replaceFIDIIDReq = new OptionalBoolRequirement("Replace FID and IID with data from Pedigree",
																																			 false);
			return register(new Step("Annotate Sample Data File", "",
															 new Requirement[][] {{sampleQCStepReq}, {createSampleDataStepReq},
																										{skipIDingDuplicatesReq, gwasQCStepReq},
																										{notGcCorrectedLrrSdReq, gcCorrectedLrrSdReq},
																										{lrrSdThresholdReq}, {callrateThresholdReq},
																										{numQReq},
																										{replaceFIDIIDReq}},
															 EnumSet.noneOf(Flag.class), priority()) {

				@Override
				public void setNecessaryPreRunProperties(Project proj,
																								 Map<Step, Map<Requirement, String>> variables) {
					double projLrrSdThreshold = proj.LRRSD_CUTOFF.getValue();
					double lrrSdThreshold = Double.parseDouble(variables.get(this).get(lrrSdThresholdReq));
					double projCallrateThreshold = proj.SAMPLE_CALLRATE_THRESHOLD.getValue();
					double callrateThreshold = Double.parseDouble(variables.get(this)
																																 .get(callrateThresholdReq));

					if (projLrrSdThreshold != lrrSdThreshold) {
						proj.LRRSD_CUTOFF.setValue(lrrSdThreshold);
					}
					if (projCallrateThreshold != callrateThreshold) {
						proj.SAMPLE_CALLRATE_THRESHOLD.setValue(callrateThreshold);
					}
				}

				@Override
				public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
					boolean checkDuplicates = !Boolean.parseBoolean(variables.get(this)
																																	 .get(skipIDingDuplicatesReq));
					String duplicatesSetFile = null;
					if (checkDuplicates) {
						duplicatesSetFile = getPlinkDir() + Qc.QC_SUBDIR + RelationAncestryQc.GENOME_DIR
																+ PLINKROOT
																+ ".genome_duplicatesSet.dat";
					}
					boolean gcCorrectedLrrSd = !Boolean.parseBoolean(variables.get(this)
																																		.get(notGcCorrectedLrrSdReq));
					int numQ = Integer.parseInt(variables.get(this).get(numQReq));
					boolean correctFidIids = Boolean.parseBoolean(variables.get(this).get(replaceFIDIIDReq));
					SampleQC.parseAndAddToSampleData(proj, numQ, 0, false, gcCorrectedLrrSd,
																					 duplicatesSetFile,
																					 correctFidIids);
				}

				@Override
				public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {

					double projLrrSdThreshold = proj.LRRSD_CUTOFF.getValue();
					double lrrSdThreshold = Double.parseDouble(variables.get(this).get(lrrSdThresholdReq));
					double projCallrateThreshold = proj.SAMPLE_CALLRATE_THRESHOLD.getValue();
					double callrateThreshold = Double.parseDouble(variables.get(this)
																																 .get(callrateThresholdReq));

					String projPropFile = proj.getPropertyFilename();

					boolean checkDuplicates = !Boolean.parseBoolean(variables.get(this)
																																	 .get(skipIDingDuplicatesReq));
					String duplicatesSetFile = null;
					if (checkDuplicates) {
						duplicatesSetFile = getPlinkDir() + Qc.QC_SUBDIR + RelationAncestryQc.GENOME_DIR
																+ PLINKROOT
																+ ".genome_duplicatesSet.dat";
					}
					boolean gcCorrectedLrrSd = !Boolean.parseBoolean(variables.get(this)
																																		.get(notGcCorrectedLrrSdReq));
					int numQ = Integer.parseInt(variables.get(this).get(numQReq));
					boolean correctFidIids = Boolean.parseBoolean(variables.get(this).get(replaceFIDIIDReq));

					String kvCmd = "";

					if (projLrrSdThreshold != lrrSdThreshold) {
						kvCmd += " LRRSD_CUTOFF=" + lrrSdThreshold;
					}
					if (projCallrateThreshold != callrateThreshold) {
						kvCmd += " SAMPLE_CALLRATE_THRESHOLD=" + callrateThreshold;
					}

					StringBuilder cmd = new StringBuilder();
					if (kvCmd.length() > 0) {
						cmd.append(Files.getRunString()).append(PROJ_PROP_UPDATE_STR + projPropFile)
							 .append(kvCmd)
							 .append("\n");
					}
					cmd.append(Files.getRunString())
						 .append(" cnv.qc.SampleQC proj=" + projPropFile + " numQ=" + numQ
										 + " justQuantiles=false" + " gcCorrectedLrrSd=" + gcCorrectedLrrSd
										 + " duplicatesSetFile=" + duplicatesSetFile + " correctFidIids="
										 + correctFidIids);
					return cmd.toString();
				}

				@Override
				public boolean checkIfOutputExists(Map<Step, Map<Requirement, String>> variables) {
					String sampleDataFile = proj.SAMPLE_DATA_FILENAME.getValue();
					if (!Files.exists(sampleDataFile)) {
						return false;
					}
					boolean checkDuplicates = !Boolean.parseBoolean(variables.get(this)
																																	 .get(skipIDingDuplicatesReq));
					String[] header = Files.getHeaderOfFile(sampleDataFile, proj.getLog());
					if (checkDuplicates && ext.indexOfStr("DuplicateId", header, false, true) == -1) {
						return false;
					}
					String[] reqHdr = {"Class=Exclude", "ExcludeNote", "Use", "UseNote", "Use_cnv",
														 "Use_cnvNote"};
					int[] facts = ext.indexFactors(reqHdr, header, false, false);
					for (int i : facts) {
						if (i == -1) {
							return false;
						}
					}
					return true;
				}

			});
		}

		private Step generateCreatePCsStep(Step transposeStep) {
			// FIXME http://genvisis.org/MitoPipeline/#illumina_marker_lists has illumina markers.. this
			// should be linked to, or
			// these steps split or something...
			final Requirement transposeStepReq = new StepRequirement(transposeStep);
			final Requirement medianMarkersReq = new FileRequirement("MedianMarkers file must exist.",
																															 "");
			final Requirement lrrSdThresholdReq = new DoubleRequirement("LRR SD threshold to filter samples.",
																																	proj.LRRSD_CUTOFF.getValue(),
																																	proj.LRRSD_CUTOFF.getMinValue(),
																																	proj.LRRSD_CUTOFF.getMaxValue());
			final Requirement callrateThresholdReq = new DoubleRequirement("Call rate threshold to filter markers.",
																																		 MitoPipeline.DEFAULT_MKR_CALLRATE_FILTER,
																																		 0.0, 1.0);
			final Requirement qcPassingOnlyReq = new OptionalBoolRequirement("Compute PCs with samples passing QC only",
																																			 true);
			final Requirement imputeNaNs = new OptionalBoolRequirement("Impute mean value for NaN", true);
			final Requirement recomputeLrrPCMarkersReq = new OptionalBoolRequirement("Should recompute Log-R ratio for PC markers?",
																																							 true);
			final Requirement recomputeLrrMedianMarkersReq = new OptionalBoolRequirement("Should recompute Log-R ratio for median markers?",
																																									 true);
			final Requirement homozygousOnlyReq = new OptionalBoolRequirement("Homozygous only?", true);
			final Requirement gcRegressionDistanceReq = new PosIntRequirement("Regression distance for the GC adjustment",
																																				GcAdjustor.DEFAULT_REGRESSION_DISTANCE[0]);
			final Requirement pcSelectionSamplesReq = new OptionalFileRequirement("A file listing a subset of samples (DNA ID) to use for determining optimal PC selection, typically a list of unrelated and single race samples. If a list is not provided, only samples passing sample qc thresholds will be used.",
																																						"");
			final Requirement externalBetaFileReq = new OptionalFileRequirement("An external beta file to optimize PC selection.",
																																					"");

			return register(new Step("Create Principal Components File and Mitochondrial Copy-Number Estimates File",
															 "",
															 new Requirement[][] {{transposeStepReq}, {medianMarkersReq},
																										{lrrSdThresholdReq}, {callrateThresholdReq},
																										{qcPassingOnlyReq}, {imputeNaNs},
																										{recomputeLrrPCMarkersReq},
																										{recomputeLrrMedianMarkersReq},
																										{homozygousOnlyReq}, {gcRegressionDistanceReq},
																										{getNumThreadsReq()}, {pcSelectionSamplesReq},
																										{externalBetaFileReq}},
															 EnumSet.noneOf(Flag.class), priority()) {

				@Override
				public void setNecessaryPreRunProperties(Project proj,
																								 Map<Step, Map<Requirement, String>> variables) {
					double sampleLRRSdFilter = Double.parseDouble(variables.get(this).get(lrrSdThresholdReq));
					if (sampleLRRSdFilter < 0) {
						switch (proj.ARRAY_TYPE.getValue()) {
							case AFFY_GW6:
							case AFFY_GW6_CN:
								proj.LRRSD_CUTOFF.setValue(0.35);
								proj.getLog()
										.reportTimeInfo("Setting " + proj.LRRSD_CUTOFF.getName()
																		+ " to default 0.35 for array " + proj.ARRAY_TYPE.getValue());
								break;
							case ILLUMINA:
								proj.LRRSD_CUTOFF.setValue(0.30);
								proj.getLog()
										.reportTimeInfo("Setting " + proj.LRRSD_CUTOFF.getName()
																		+ " to default 0.30 for array " + proj.ARRAY_TYPE.getValue());
								break;
							default:
								throw new IllegalArgumentException("Invalid Array type");
						}
					} else {
						proj.LRRSD_CUTOFF.setValue(sampleLRRSdFilter);
					}
				}

				@Override
				public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
					String medianMarkers = variables.get(this).get(medianMarkersReq);
					double markerCallRateFilter = Double.parseDouble(variables.get(this)
																																		.get(callrateThresholdReq));
					// FIXME: This gcCorrect assignment was carried over from the old indexed version but
					// appears incorrect
					boolean gcCorrect = Boolean.parseBoolean(variables.get(this).get(qcPassingOnlyReq));
					boolean imputeMeanForNaN = Boolean.parseBoolean(variables.get(this).get(imputeNaNs));
					boolean recomputeLRRPCs = Boolean.parseBoolean(variables.get(this)
																																	.get(recomputeLrrPCMarkersReq));
					boolean recomputeLRRMedian = Boolean.parseBoolean(variables.get(this)
																																		 .get(recomputeLrrMedianMarkersReq));
					boolean homozygousOnly = Boolean.parseBoolean(variables.get(this).get(homozygousOnlyReq));
					int bpGcModel = GcAdjustor.GcModel.DEFAULT_GC_MODEL_BIN_FASTA;
					int regressionDistance = Integer.parseInt(variables.get(this)
																														 .get(gcRegressionDistanceReq));
					int numComponents = MitoPipeline.DEFAULT_NUM_COMPONENTS;
					int numThreads = resolveThreads(variables.get(this).get(getNumThreadsReq()));
					String outputBase = MitoPipeline.FILE_BASE;

					String betaOptFile = variables.get(this).get(pcSelectionSamplesReq);
					String betaFile = variables.get(this).get(externalBetaFileReq);

					boolean markerQC = true;
					double[] pvalOpt = MitoPipeline.DEFAULT_PVAL_OPTS;
					String pedFile = null;
					String useFile = null;
					boolean sampLrr = true;
					boolean plot = false;
					int retCode = PCAPrep.prepPCA(proj, numThreads, outputBase, markerQC,
																				markerCallRateFilter,
																				useFile, proj.getSampleList(), proj.getLog());
					if (retCode == PCAPrep.SUCCESS_CODE) {
						MitoPipeline.estimateMtDNACN(proj, numThreads, medianMarkers, numComponents, outputBase,
																				 homozygousOnly, markerCallRateFilter, betaOptFile, pedFile,
																				 recomputeLRRPCs, recomputeLRRMedian, sampLrr,
																				 imputeMeanForNaN, gcCorrect, bpGcModel, regressionDistance,
																				 proj.GENOME_BUILD_VERSION.getValue(), pvalOpt, betaFile,
																				 plot, false, proj.getLog());
					} else {
						setFailed(PCAPrep.errorMessage(retCode));
					}
				}

				@Override
				public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
					String medianMarkers = variables.get(this).get(medianMarkersReq);
					double lrrSD = Double.parseDouble(variables.get(this).get(lrrSdThresholdReq));
					double markerCallRateFilter = Double.parseDouble(variables.get(this)
																																		.get(callrateThresholdReq));
					// FIXME: This gcCorrect assignment was carried over from the old indexed version but
					// appears incorrect
					boolean gcCorrect = Boolean.parseBoolean(variables.get(this).get(qcPassingOnlyReq));
					boolean imputeMeanForNaN = Boolean.parseBoolean(variables.get(this).get(imputeNaNs));
					boolean recomputeLRRPCs = Boolean.parseBoolean(variables.get(this)
																																	.get(recomputeLrrPCMarkersReq));
					boolean recomputeLRRMedian = Boolean.parseBoolean(variables.get(this)
																																		 .get(recomputeLrrMedianMarkersReq));
					boolean homozygousOnly = Boolean.parseBoolean(variables.get(this).get(homozygousOnlyReq));
					int bpGcModel = GcAdjustor.GcModel.DEFAULT_GC_MODEL_BIN_FASTA;
					int regressionDistance = Integer.parseInt(variables.get(this)
																														 .get(gcRegressionDistanceReq));
					int numComponents = MitoPipeline.DEFAULT_NUM_COMPONENTS;
					int numThreads = resolveThreads(variables.get(this).get(getNumThreadsReq()));
					String outputBase = MitoPipeline.FILE_BASE;

					String betaOptFile = variables.get(this).get(pcSelectionSamplesReq);
					String betaFile = variables.get(this).get(externalBetaFileReq);
					boolean sampLrr = true;


					String projPropFile = proj.getPropertyFilename();
					StringBuilder cmd = new StringBuilder();
					cmd.append(Files.getRunString()).append(" org.genvisis.cnv.manage.MitoPipeline")
						 .append(" proj=").append(projPropFile).append(" mitochondrialMarkers=")
						 .append(medianMarkers).append(" numComponents=").append(numComponents)
						 .append(" imputeMeanForNaN=").append(imputeMeanForNaN).append(" recomputeLRR_PCs=")
						 .append(recomputeLRRPCs).append(" recomputeLRR_Median=").append(recomputeLRRMedian)
						 .append(" gcCorrect=").append(gcCorrect).append(" bpGcModel=").append(bpGcModel)
						 .append(" LRRSD=").append(lrrSD).append(" markerCallRate=")
						 .append(markerCallRateFilter)
						 .append(" regressionDistance=").append(regressionDistance).append(" sampLRR=")
						 .append(sampLrr).append(" ").append(PSF.Ext.NUM_THREADS_COMMAND).append(numThreads)
						 .append(" log=").append(proj.getLog().getFilename()).append(" output=")
						 .append(outputBase);
					if (!"".equals(betaOptFile)) {
						cmd.append(" ").append(MitoPipeline.PC_OPT_FILE).append("=").append(betaOptFile);
					}
					if (!"".equals(betaFile)) {
						cmd.append(" betas=").append(betaFile);
					}
					if (!homozygousOnly) {
						cmd.append(" -allCalls ");
					}

					cmd.append(" -SkipProjectCreationWithLongUndocumentedFlag ");

					return cmd.toString();
				}

				@Override
				public boolean checkIfOutputExists(Map<Step, Map<Requirement, String>> variables) {
					String outputBase = proj.PROJECT_DIRECTORY.getValue() + MitoPipeline.FILE_BASE;
					String finalReport = outputBase + PCA.FILE_EXTs[0];// PrincipalComponentsResiduals.MT_REPORT_EXT[0];
					// boolean mkrFiles = true;
					// for (String file :
					// PrincipalComponentsResiduals.MT_REPORT_MARKERS_USED) {
					// if (!Files.exists(outputBase + file)) {
					// mkrFiles = false;
					// break;
					// }
					// }
					return Files.exists(finalReport) /* && mkrFiles */;
				}
			});
		}

		private Step generatePFBStep(final Step parseSamplesStep) {
			final Requirement parseSamplesStepReq = new StepRequirement(parseSamplesStep);
			final Requirement sampleSubsetReq = new FileRequirement("A Sample subset file must exist.",
																															proj.SAMPLE_SUBSET_FILENAME.getValue());
			String defaultOutputFile;
			if (Files.exists(proj.SAMPLE_SUBSET_FILENAME.getValue())) {
				defaultOutputFile = ext.rootOf(proj.SAMPLE_SUBSET_FILENAME.getValue()) + ".pfb";
			} else {
				defaultOutputFile = proj.CUSTOM_PFB_FILENAME.getValue();
			}
			final Requirement outputFileReq = new OutputFileRequirement("PFB (population BAF) output file must be specified.",
																																	defaultOutputFile);
			return register(new Step("Compute Population BAF files", "",
															 new Requirement[][] {{parseSamplesStepReq, sampleSubsetReq},
																										{outputFileReq}},
															 EnumSet.noneOf(Flag.class), priority()) {

				@Override
				public void setNecessaryPreRunProperties(Project proj,
																								 Map<Step, Map<Requirement, String>> variables) {
					String setSubSampFile = proj.SAMPLE_SUBSET_FILENAME.getValue();
					String subSampFile = variables.get(this).get(sampleSubsetReq);
					String setPFBFile = proj.CUSTOM_PFB_FILENAME.getValue();
					String pfbOutputFile = variables.get(this).get(outputFileReq);

					if (!ext.verifyDirFormat(setSubSampFile).equals(subSampFile)) {
						proj.SAMPLE_SUBSET_FILENAME.setValue(subSampFile);
					}
					if (!ext.verifyDirFormat(setPFBFile).equals(pfbOutputFile)) {
						proj.CUSTOM_PFB_FILENAME.setValue(pfbOutputFile);
					}
				}

				@Override
				public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
					org.genvisis.cnv.analysis.PennCNV.populationBAF(proj);
				}

				@Override
				public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
					String kvCmd = "";

					String setSubSampFile = proj.SAMPLE_SUBSET_FILENAME.getValue();
					String subSampFile = variables.get(this).get(sampleSubsetReq);
					String setPFBFile = proj.CUSTOM_PFB_FILENAME.getValue();
					String pfbOutputFile = variables.get(this).get(outputFileReq);

					if (!ext.verifyDirFormat(setSubSampFile).equals(subSampFile)) {
						kvCmd += " SAMPLE_SUBSET_FILENAME=" + subSampFile;
					}
					if (!ext.verifyDirFormat(setPFBFile).equals(pfbOutputFile)) {
						kvCmd += " CUSTOM_PFB_FILENAME=" + pfbOutputFile;
					}

					String projPropFile = proj.getPropertyFilename();
					StringBuilder cmd = new StringBuilder();
					if (kvCmd.length() > 0) {
						cmd.append(Files.getRunString()).append(PROJ_PROP_UPDATE_STR + projPropFile)
							 .append(kvCmd)
							 .append("\n");
					}
					return cmd.append(Files.getRunString())
										.append(" cnv.analysis.PennCNV -pfb proj=" + proj.getPropertyFilename()
														+ " log="
														+ proj.getLog().getFilename())
										.toString();
				}

				@Override
				public boolean checkIfOutputExists(Map<Step, Map<Requirement, String>> variables) {
					String subSampFile = variables.get(this).get(sampleSubsetReq);
					String pfbOutputFile = variables.get(this).get(outputFileReq);
					return Files.exists(pfbOutputFile) || Files.exists(ext.rootOf(subSampFile) + ".pfb");
				}
			});
		}

		private Step generateSexCentroidsStep(Step gcModelStep) {
			final Requirement gcModelStepReq = new StepRequirement(gcModelStep);
			final Requirement gcModelFileReq = new FileRequirement("Full GC Model File.",
																														 proj.GC_MODEL_FILENAME.getValue());
			return register(new Step("Create Sex-Specific Centroids; Filter PFB and GCMODEL Files", "",
															 new Requirement[][] {{gcModelStepReq, gcModelFileReq},
																										{getNumThreadsReq()},},
															 EnumSet.of(Flag.RUNTIME),
															 priority()) {

				@Override
				public void setNecessaryPreRunProperties(Project proj,
																								 Map<Step, Map<Requirement, String>> variables) {

					int numThreads = resolveThreads(variables.get(this).get(getNumThreadsReq()));
					maybeSetProjNumThreads(numThreads);
				}

				@Override
				public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
					String malePFB;
					String femalePFB;
					String centFilePathM;
					String centFilePathF;
					String newGCFile;
					String outputDir = proj.DATA_DIRECTORY.getValue();
					newGCFile = outputDir + "sexSpecific.gcModel";
					malePFB = outputDir + "males.pfb";
					femalePFB = outputDir + "females.pfb";
					centFilePathM = outputDir + "sexSpecific_Male.cent";
					centFilePathF = outputDir + "sexSpecific_Female.cent";

					int numThreads = resolveThreads(variables.get(this).get(getNumThreadsReq()));
					String gcModelFile = variables.get(this).get(gcModelFileReq);
					Centroids.computeSexSpecificCentroids(proj, new String[] {malePFB, femalePFB},
																								new String[] {centFilePathM, centFilePathF},
																								numThreads);

					AnalysisFormats.filterSexSpecificGCModel(proj, gcModelFile, newGCFile);
				}

				@Override
				public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {

					int numThreads = resolveThreads(variables.get(this).get(getNumThreadsReq()));
					String mainCmd = Files.getRunString() + " cnv.filesys.Centroids proj="
													 + proj.getPropertyFilename() + " -sexSpecific "
													 + PSF.Ext.NUM_THREADS_COMMAND + numThreads;
					String gcModelFile = variables.get(this).get(gcModelFileReq);
					String gcCmd = Files.getRunString() + " cnv.analysis.AnalysisFormats proj="
												 + proj.getPropertyFilename() + " gcmodel=" + gcModelFile;
					return mainCmd + "\n" + gcCmd;
				}

				@Override
				public boolean checkIfOutputExists(Map<Step, Map<Requirement, String>> variables) {
					String malePFB;
					String femalePFB;
					String centFilePathM;
					String centFilePathF;
					String newGCFile;
					String outputDir = proj.DATA_DIRECTORY.getValue();
					malePFB = outputDir + "males.pfb";
					femalePFB = outputDir + "females.pfb";
					centFilePathM = outputDir + "sexSpecific_Male.cent";
					centFilePathF = outputDir + "sexSpecific_Female.cent";
					newGCFile = outputDir + "sexSpecific.gcModel";
					boolean exists = Files.exists(malePFB);
					exists = exists && Files.exists(femalePFB);
					exists = exists && Files.exists(centFilePathM);
					exists = exists && Files.exists(centFilePathF);
					exists = exists && Files.exists(newGCFile);
					return exists;
				}

			});
		}

		private Step generateCNVStep(Step pfbStep, Step gcModelStep) {
			final Requirement hmmFile = new FileRequirement("Hidden Markov Model File Must Exist",
																											proj.HMM_FILENAME.getValue());
			final Requirement pfbStepReq = new StepRequirement(pfbStep);
			final Requirement pfbFileReq = new FileRequirement("PFB File Must Exist",
																												 proj.CUSTOM_PFB_FILENAME.getValue());
			final Requirement gcModelStepReq = new StepRequirement(gcModelStep);
			final Requirement gcModelFileReq = new FileRequirement("GCMODEL File Must Exist",
																														 proj.GC_MODEL_FILENAME.getValue());
			final Requirement outputFileReq = new OutputFileRequirement("Output filename.",
																																	"cnvs/genvisis.cnv") {
				@Override
				public boolean checkRequirement(String arg, Set<Step> stepSelections,
																				Map<Step, Map<Requirement, String>> variables) {
					return super.checkRequirement(proj.PROJECT_DIRECTORY.getValue() + arg, stepSelections,
																				variables);
				}
			};

			return register(new Step("Call CNVs", "",
															 new Requirement[][] {{hmmFile}, {pfbStepReq, pfbFileReq},
																										{gcModelStepReq, gcModelFileReq},
																										{getNumThreadsReq()},
																										{outputFileReq}},
															 EnumSet.of(Flag.MEMORY),
															 priority()) {

				@Override
				public void setNecessaryPreRunProperties(Project proj,
																								 Map<Step, Map<Requirement, String>> variables) {
					String hmmP = proj.HMM_FILENAME.getValue();
					String hmmG = variables.get(this).get(hmmFile);
					if (!hmmP.equals(hmmG)) {
						proj.HMM_FILENAME.setValue(hmmG);
					}
					String pfbP = proj.CUSTOM_PFB_FILENAME.getValue();
					String pfbG = variables.get(this).get(pfbFileReq);
					if (!pfbP.equals(pfbG)) {
						proj.CUSTOM_PFB_FILENAME.setValue(pfbG);
					}
					String gcmP = proj.GC_MODEL_FILENAME.getValue();
					String gcmG = variables.get(this).get(gcModelFileReq);
					if (!gcmP.equals(gcmG)) {
						proj.GC_MODEL_FILENAME.setValue(gcmG);
					}
					int numThreads = resolveThreads(variables.get(this).get(getNumThreadsReq()));
					maybeSetProjNumThreads(numThreads);
				}

				@Override
				public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
					int numThreads = resolveThreads(variables.get(this).get(getNumThreadsReq()));
					maybeSetProjNumThreads(numThreads);
					String output = variables.get(this).get(outputFileReq); // gets PROJ_DIR prepended, so NOT
																																	// ABSOLUTE
					(new File(ext.parseDirectoryOfFile(proj.PROJECT_DIRECTORY.getValue() + output))).mkdirs();
					CNVCaller.callAutosomalCNVs(proj, output, proj.getSamples(), null, null,
																			CNVCaller.DEFAULT_MIN_SITES, CNVCaller.DEFAULT_MIN_CONF,
																			PFB_MANAGEMENT_TYPE.PENNCNV_DEFAULT, numThreads, 1);
					proj.CNV_FILENAMES.addValue(proj.PROJECT_DIRECTORY.getValue() + output);
					proj.saveProperties(new Property[] {proj.CNV_FILENAMES});
				}

				@Override
				public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
					String kvCmd = "";

					String hmmP = proj.HMM_FILENAME.getValue();
					String hmmG = variables.get(this).get(hmmFile);
					if (!hmmP.equals(hmmG)) {
						kvCmd += " HMM_FILENAME=" + hmmG;
					}
					String pfbP = proj.CUSTOM_PFB_FILENAME.getValue();
					String pfbG = variables.get(this).get(pfbFileReq);
					if (!pfbP.equals(pfbG)) {
						kvCmd += " CUSTOM_PFB_FILENAME=" + pfbG;
					}
					String gcmP = proj.GC_MODEL_FILENAME.getValue();
					String gcmG = variables.get(this).get(gcModelFileReq);
					if (!gcmP.equals(gcmG)) {
						kvCmd += " GC_MODEL_FILENAME=" + gcmG;
					}

					int numThreads = resolveThreads(variables.get(this).get(getNumThreadsReq()));
					if (numThreads != proj.NUM_THREADS.getValue()) {
						kvCmd += " " + proj.NUM_THREADS.getName() + "=" + numThreads;
					}
					String projPropFile = proj.getPropertyFilename();
					StringBuilder cmd = new StringBuilder();
					if (kvCmd.length() > 0) {
						cmd.append(Files.getRunString()).append(PROJ_PROP_UPDATE_STR + projPropFile)
							 .append(kvCmd)
							 .append("\n");
					}
					return cmd.append(Files.getRunString()).append(" cnv.hmm.CNVCaller proj=" + projPropFile)
										.append(" out=" + variables.get(this).get(outputFileReq)).append(" ")
										.append(PSF.Ext.NUM_THREADS_COMMAND).append(numThreads).toString();
				}

				@Override
				public boolean checkIfOutputExists(Map<Step, Map<Requirement, String>> variables) {
					String output = variables.get(this).get(outputFileReq);
					return Files.exists(proj.PROJECT_DIRECTORY.getValue() + output);
				}

			});
		}

		private Step generateShadowStep(final Step parseSamplesStep) {
			final Requirement parseSamplesStepReq = new StepRequirement(parseSamplesStep);
			final Requirement numPCsReq = new PosIntRequirement("Number of principal components for correction.",
																													MitoPipeline.DEFAULT_NUM_COMPONENTS);
			final Requirement outputBaseReq = new OutputFileRequirement("Output file path (relative to project directory) and baseName for principal components correction files",
																																	MitoPipeline.FILE_BASE) {
				@Override
				public boolean checkRequirement(String arg, Set<Step> stepSelections,
																				Map<Step, Map<Requirement, String>> variables) {
					String outputBase = proj.PROJECT_DIRECTORY.getValue() + arg;
					String finalReport = outputBase + PCA.FILE_EXTs[0];
					return super.checkRequirement(finalReport, stepSelections, variables);
				}
			};
			final Requirement callrateReq = new DoubleRequirement("Call-rate filter for determining high-quality markers",
																														MitoPipeline.DEFAULT_MKR_CALLRATE_FILTER,
																														0.0, 1.0);
			final Requirement recomputeLrrReq = new OptionalBoolRequirement("Re-compute Log-R Ratio values? (usually false if LRRs already exist)",
																																			false);
			final Requirement tempDirReq = new OptionalFileRequirement("Temporary directory for intermediate files (which tend to be very large)",
																																 "");
			final Requirement correctionStrategyReq = new EnumRequirement("Correction Type",
																																		CORRECTION_TYPE.XY);
			final Requirement sexChromosomeStrategyReq = new EnumRequirement("Sex Chromosome Strategy",
																																			 CHROMOSOME_X_STRATEGY.BIOLOGICAL);

			return register(new Step("Create 'Shadow' Project", "",
															 new Requirement[][] {{parseSamplesStepReq}, {numPCsReq},
																										{outputBaseReq},
																										{callrateReq}, {recomputeLrrReq}, {tempDirReq},
																										{correctionStrategyReq},
																										{sexChromosomeStrategyReq},
																										{getNumThreadsReq()},},
															 EnumSet.of(Flag.MEMORY, Flag.RUNTIME), priority()) {

				@Override
				public void setNecessaryPreRunProperties(Project proj,
																								 Map<Step, Map<Requirement, String>> variables) {
					// not needed for step
				}

				@Override
				public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
					int numComponents = Integer.parseInt(variables.get(this).get(numPCsReq));
					String outputBase = variables.get(this).get(outputBaseReq);
					double markerCallRateFilter = Double.parseDouble(variables.get(this).get(callrateReq));
					boolean recomputeLRRPCs = Boolean.parseBoolean(variables.get(this).get(recomputeLrrReq));
					String tmpDir = variables.get(this).get(tempDirReq);
					if ("".equals(tmpDir.trim())) {
						tmpDir = null;
					}
					CORRECTION_TYPE type = CORRECTION_TYPE.valueOf(variables.get(this)
																																	.get(correctionStrategyReq));
					CHROMOSOME_X_STRATEGY strategy = CHROMOSOME_X_STRATEGY.valueOf(variables.get(this)
																																									.get(sexChromosomeStrategyReq));

					int totalThreads = resolveThreads(variables.get(this).get(getNumThreadsReq()));
					String retMsg = PRoCtOR.shadow(proj, tmpDir, outputBase, markerCallRateFilter,
																				 recomputeLRRPCs, type, strategy, numComponents,
																				 totalThreads);
					if (!"".equals(retMsg)) {
						setFailed(retMsg);
					}
				}

				@Override
				public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
					int numComponents = Integer.parseInt(variables.get(this).get(numPCsReq));
					String outputBase = variables.get(this).get(outputBaseReq);
					double markerCallRateFilter = Double.parseDouble(variables.get(this).get(callrateReq));
					boolean recomputeLRRPCs = Boolean.parseBoolean(variables.get(this).get(recomputeLrrReq));
					String tmpDir = variables.get(this).get(tempDirReq);
					if ("".equals(tmpDir.trim())) {
						tmpDir = null;
					}
					String correctionType = variables.get(this).get(correctionStrategyReq);
					String strategy = variables.get(this).get(sexChromosomeStrategyReq);

					int totalThreads = resolveThreads(variables.get(this).get(getNumThreadsReq()));

					String projPropFile = proj.getPropertyFilename();
					StringBuilder cmd = new StringBuilder();
					cmd.append(Files.getRunString()).append(" org.genvisis.cnv.manage.PRoCtOR")
						 .append(" proj=")
						 .append(projPropFile).append(" numComponents=").append(numComponents)
						 .append(" outputBase=").append(outputBase).append(" callrate=")
						 .append(markerCallRateFilter).append(" recomputeLRR=").append(recomputeLRRPCs)
						 .append(" type=").append(correctionType).append(" sexStrategy=").append(strategy)
						 .append(" numThreads=").append(totalThreads);
					if (tmpDir != null) {
						cmd.append(" tmp=").append(tmpDir);
					}

					return cmd.toString();
				}

				@Override
				public boolean checkIfOutputExists(Map<Step, Map<Requirement, String>> variables) {
					String outputBase = proj.PROJECT_DIRECTORY.getValue()
															+ variables.get(this).get(outputBaseReq);
					String finalReport = outputBase + PCA.FILE_EXTs[0];
					return Files.exists(finalReport);
				}

			});

		}

		private Step generateABLookupStep(final Step parseSamplesStep) {
			final Requirement parseSamplesStepReq = new StepRequirement(parseSamplesStep);
			return register(new Step("Generate AB Lookup File", "",
															 new Requirement[][] {{parseSamplesStepReq}},
															 EnumSet.of(Flag.RUNTIME), priority()) {

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

					if (ABLookup.fillInMissingAlleles(proj, filename, proj.getLocationOfSNP_Map(true),
																						false)) {
						ABLookup.applyABLookupToFullSampleFiles(proj, filename);
					} else {
						setFailed("Failed to fill in missing alleles - please check log for more info.");
					}
				}

				@Override
				public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
					String filename = proj.PROJECT_DIRECTORY.getValue()
														+ ext.addToRoot(ABLookup.DEFAULT_AB_FILE, "_parsed");
					String projFile = proj.getPropertyFilename();
					String mapFile = proj.getLocationOfSNP_Map(true);

					List<String> baseCommand = ImmutableList.of(Files.getRunString(),
																											ABLookup.class.getName(),
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
				public boolean checkIfOutputExists(Map<Step, Map<Requirement, String>> variables) {
					return Files.exists(proj.AB_LOOKUP_FILENAME.getValue(false, false));
				}
			});
		}

	}

	private SortedSet<Step> generateSteps() {
		StepBuilder sb = new StepBuilder();

		Step parseSamplesStep;
		if (proj.getArrayType() == ARRAY.AFFY_GW6 || proj.getArrayType() == ARRAY.AFFY_GW6_CN) {
			parseSamplesStep = sb.generateAffyParseSamplesStep();
		} else {
			Step markerPositions = sb.generateMarkerPositionsStep();
			parseSamplesStep = sb.generateIlluminaParseSamplesStep(markerPositions);
		}
		Step createSampleDataStep = sb.generateCreateSampleDataStep(parseSamplesStep);
		Step transposeStep = sb.generateTransposeStep(parseSamplesStep);
		Step gcModelStep = sb.generateGCModelStep();
		Step sampleQCStep = sb.generateSampleQCStep(parseSamplesStep);
		sb.generateMarkerQCStep(parseSamplesStep);
		sb.generateSexChecksStep(parseSamplesStep, createSampleDataStep,
														 transposeStep, sampleQCStep);
		sb.generateABLookupStep(parseSamplesStep);
		Step plinkExportStep = sb.generatePlinkExportStep(parseSamplesStep);
		Step gwasQCStep = sb.generateGwasQCStep(plinkExportStep);
		Step ancestryStep = sb.generateAncestryStep(gwasQCStep);
		sb.generateFurtherAnalysisQCStep(plinkExportStep,
																		 gwasQCStep, ancestryStep);
		sb.generateMosaicArmsStep(parseSamplesStep);
		sb.generateAnnotateSampleDataStep(sampleQCStep,
																			createSampleDataStep,
																			gwasQCStep);
		sb.generateCreatePCsStep(transposeStep);
		Step pfbStep = sb.generatePFBStep(parseSamplesStep);
		sb.generateSexCentroidsStep(gcModelStep);
		sb.generateCNVStep(pfbStep, gcModelStep);
		sb.generateShadowStep(parseSamplesStep);

		return sb.getSortedSteps();
	}

	private String getLocationOfSampleMap() {
		String filename;

		String projDir = proj.PROJECT_DIRECTORY.getValue();
		String snpMap = "Sample_Map.csv";
		String snpMapGz = "Sample_Map.csv.gz";
		if (Files.exists(projDir + snpMap)) {
			filename = projDir + snpMap;
		} else if (Files.exists(projDir + snpMapGz)) {
			filename = projDir + snpMapGz;
		} else {
			String srcDir = proj.SOURCE_DIRECTORY.getValue();
			if (Files.exists(srcDir + snpMap)) {
				filename = srcDir + snpMap;
			} else if (Files.exists(srcDir + snpMapGz)) {
				filename = srcDir + snpMapGz;
			} else {
				return null;
			}
		}
		return filename;
	}

	private String getPlinkDir() {
		return ext.verifyDirFormat(proj.PROJECT_DIRECTORY.getValue() + PLINK_SUBDIR);
	}

	private String getAncestryDir() {
		return getPlinkDir() + Qc.QC_SUBDIR + RelationAncestryQc.ANCESTRY_DIR;
	}

	private int resolveThreads(String arg) {
		int numThreads = Requirement.checkIntArgOrNeg1(arg);
		if (numThreads <= 0) {
			numThreads = proj.NUM_THREADS.getValue();
		}
		return numThreads;
	}

	private void maybeSetProjNumThreads(int numThreads) {
		if (numThreads != proj.NUM_THREADS.getValue()) {
			proj.NUM_THREADS.setValue(numThreads);
		}
	}

}
