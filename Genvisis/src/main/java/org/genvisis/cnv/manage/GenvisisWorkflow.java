package org.genvisis.cnv.manage;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
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
import org.genvisis.cnv.qc.MarkerBlast.FILE_SEQUENCE_TYPE;
import org.genvisis.cnv.qc.MarkerBlastQC;
import org.genvisis.cnv.qc.MarkerMetrics;
import org.genvisis.cnv.qc.SampleQC;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Elision;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.ext;
import org.genvisis.gwas.Ancestry;
import org.genvisis.gwas.PlinkMendelianChecker;
import org.genvisis.gwas.Qc;
import org.genvisis.gwas.RelationAncestryQc;

import com.google.common.base.Joiner;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

public class GenvisisWorkflow {

	public abstract static class Step implements Comparable<Step> {
		private String name;
		private String desc;
		private GenvisisWorkflow.Requirement[][] requirements;
		private boolean failed = false;
		private ArrayList<String> failReasons = new ArrayList<String>();
		private final Set<Step> relatedSteps;
		private final Set<GenvisisWorkflow.FLAG> stepFlags;
		private final double priority;

		Step(String name, String desc,
				 GenvisisWorkflow.Requirement[][] requirements,
				 Step[] relatedSteps,
				 GenvisisWorkflow.FLAG[] flags, double priority) {
			this.name = name;
			this.desc = desc;
			this.requirements = requirements;
			final Set<Step> steps = new HashSet<Step>();
			steps.add(this);
			if (relatedSteps != null) {
				for (final Step s : relatedSteps) {
					steps.add(s);
					steps.addAll(s.getRelatedSteps());
				}
			}
			this.relatedSteps = Collections.unmodifiableSet(steps);
			final Set<GenvisisWorkflow.FLAG> flgs = new HashSet<GenvisisWorkflow.FLAG>();
			if (flags != null) {
				for (final GenvisisWorkflow.FLAG f : flags) {
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
																											Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables);

		public abstract void run(Project proj,
														 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables);

		public void gracefulDeath(Project proj) {
			return;
		}

		public boolean hasRequirements(Project proj, Set<Step> stepSelections,
																	 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
			if (variables.get(this) == null) {
				return false;
			}
			for (GenvisisWorkflow.Requirement[] group : getRequirements()) {
				boolean groupMet = false;
				for (GenvisisWorkflow.Requirement req : group) {
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
		 * @return An array of {@link GenvisisWorkflow.Requirement}s. At least one element of each
		 *         subarray must be met to satisfy the step pre-requisites - effectively this means
		 *         elements of the first array are AND'd together, while elements of the second array
		 *         are OR'd.
		 */
		public GenvisisWorkflow.Requirement[][] getRequirements() {
			// TODO unify requirement names, AND/OR structure, input types and default values to avoid
			// maintaining these parallel arrays
			return requirements;
		}

		public abstract boolean checkIfOutputExists(Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables);

		public void resetRun() {
			failed = false;
			failReasons.clear();
		}

		public abstract String getCommandLine(Project proj,
																					Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables);

		/**
		 * @return A {@link Collection} of the complete network of {@link Step}s related to this
		 *         {@code Step} - including this {@code Step}, direct and transitive dependencies.
		 */
		public Collection<Step> getRelatedSteps() {
			return relatedSteps;
		}

		public Collection<GenvisisWorkflow.FLAG> getFlags() {
			return stepFlags;
		}



		public double getPriority() {
			return priority;
		}

		@Override
		public int compareTo(Step o) {
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

		private static String stepReqMessage(Step requiredStep) {
			return "[" + requiredStep.getName()
						 + "] step must have been run already or must be selected and valid";
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

		protected BoolRequirement(String description, boolean defaultValue) {
			super(description, RequirementInputType.BOOL, defaultValue);
		}

		@Override
		public boolean checkRequirement(String arg, Set<Step> stepSelections,
																		Map<Step, Map<Requirement, String>> variables) {
			return Boolean.parseBoolean(arg);
		}

	}

	public static class OptionalBoolRequirement extends BoolRequirement {
		protected OptionalBoolRequirement(String description, boolean defaultValue) {
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

		protected DoubleRequirement(String description, double defaultValue, double min, double max) {
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

		protected IntRequirement(String description, int defaultValue, int min, int max) {
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

	public static class EnumRequirement extends Requirement {

		protected EnumRequirement(String description, Enum<?> defaultValue) {
			super(description, RequirementInputType.ENUM, defaultValue);
		}

		@Override
		public boolean checkRequirement(String arg, Set<Step> stepSelections,
																		Map<Step, Map<Requirement, String>> variables) {
			return true;
		}

	}

	public enum RequirementInputType {
		NONE, FILE, DIR, STRING, NUMBER, BOOL, ENUM
	}

	public enum FLAG {
		MEMORY, RUNTIME, MULTITHREADED
	}

	private static final String PROJ_PROP_UPDATE_STR = " org.genvisis.cnv.filesys.Project proj=";
	final Project proj;
	private final SortedSet<Step> steps;
	private final GenvisisWorkflow.Requirement numThreadsReq;
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
		numThreadsReq = new GenvisisWorkflow.Requirement("Number of Threads to Use",
																										 GenvisisWorkflow.RequirementInputType.NUMBER,
																										 proj.NUM_THREADS.getValue()) {

			@Override
			public boolean checkRequirement(String arg, Set<Step> stepSelections,
																			Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				int numThreads = checkIntArgOrNeg1(arg);
				return numThreads > 0;
			}
		};

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

	private SortedSet<Step> generateSteps() {
		SortedSet<Step> buildSteps = Sets.newTreeSet();
		double priority = 1.0;
		Step parseSamplesStep;
		if (proj.getArrayType() == ARRAY.AFFY_GW6 || proj.getArrayType() == ARRAY.AFFY_GW6_CN) {
			parseSamplesStep = generateAffyParseSamplesStep(priority++);
		} else {
			Step markerPositions = generateMarkerPositionsStep(priority++);
			buildSteps.add(markerPositions);
			parseSamplesStep = generateIlluminaParseSamplesStep(priority++, markerPositions);
		}
		Step createSampleDataStep = generateCreateSampleDataStep(priority++,
																														 parseSamplesStep);
		Step transposeStep = generateTransposeStep(priority++, parseSamplesStep);
		Step gcModelStep = generateGCModelStep(priority++);
		Step sampleQCStep = generateSampleQCStep(priority++, parseSamplesStep);
		Step markerQCStep = generateMarkerQCStep(priority++, parseSamplesStep);
		Step sexChecksStep = generateSexChecksStep(priority++, parseSamplesStep,
																							 createSampleDataStep, transposeStep,
																							 sampleQCStep);
		Step abLookupStep = generateABLookupStep(priority++, parseSamplesStep);
		Step plinkExportStep = generatePlinkExportStep(priority++, parseSamplesStep);
		Step gwasQCStep = generateGwasQCStep(priority++, plinkExportStep);
		Step mosaicArmsStep = generateMosaicArmsStep(priority++, parseSamplesStep);
		Step annotateSampleDataStep = generateAnnotateSampleDataStep(priority++,
																																 sampleQCStep,
																																 createSampleDataStep,
																																 gwasQCStep);
		Step createPCsStep = generateCreatePCsStep(priority++, transposeStep);
		Step pfbStep = generatePFBStep(priority++, parseSamplesStep);
		Step sexCentroidsStep = generateSexCentroidsStep(priority++, gcModelStep);
		Step cnvStep = generateCNVStep(priority++, pfbStep, gcModelStep);
		Step shadowStep = generateShadowStep(priority++, parseSamplesStep);

		buildSteps.add(parseSamplesStep);
		buildSteps.add(createSampleDataStep);
		buildSteps.add(gcModelStep);
		buildSteps.add(sampleQCStep);
		buildSteps.add(markerQCStep);
		buildSteps.add(sexChecksStep);
		buildSteps.add(abLookupStep);
		buildSteps.add(plinkExportStep);
		buildSteps.add(gwasQCStep);
		buildSteps.add(mosaicArmsStep);
		buildSteps.add(annotateSampleDataStep);
		buildSteps.add(createPCsStep);
		buildSteps.add(pfbStep);
		buildSteps.add(sexCentroidsStep);
		buildSteps.add(cnvStep);
		buildSteps.add(shadowStep);

		return buildSteps;
	}

	private Step generateMarkerPositionsStep(final double priority) {
		final GenvisisWorkflow.Requirement snpMapReq = new GenvisisWorkflow.FileRequirement("An Illumina SNP_map file.",
																																												proj.getLocationOfSNP_Map(false));
		final GenvisisWorkflow.Requirement manifestReq = new GenvisisWorkflow.FileRequirement("An Illumina Manifest file.",
																																													proj.getLocationOfSNP_Map(false));
		return new Step("Create Marker Positions (if not already exists)", "",
										new GenvisisWorkflow.Requirement[][] {{snpMapReq}, {manifestReq}}, null, null,
										priority) {

			@Override
			public void setNecessaryPreRunProperties(Project proj,
																							 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				// not needed for step
			}

			@Override
			public void run(Project proj,
											Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				proj.getLog().report("Generating marker positions file");
				String snpMap = variables.get(this).get(snpMapReq);
				String manifest = variables.get(this).get(manifestReq);
				if (Files.exists(snpMap)) {
					org.genvisis.cnv.manage.Markers.generateMarkerPositions(proj, snpMap);
				} else if (Files.exists(manifest)) {
					org.genvisis.cnv.manage.Markers.extractMarkerPositionsFromManifest(manifest,
																																						 ARRAY.ILLUMINA,
																																						 proj.GENOME_BUILD_VERSION.getValue(),
																																						 FILE_SEQUENCE_TYPE.MANIFEST_FILE,
																																						 proj.MARKER_POSITION_FILENAME.getValue(false,
																																																										false),
																																						 proj.getLog());
				}
			}

			@Override
			public boolean checkIfOutputExists(Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				return Files.exists(proj.MARKER_POSITION_FILENAME.getValue(false, false));
			}

			@Override
			public String getCommandLine(Project proj,
																	 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				String projFile = proj.getPropertyFilename();
				String snpMap = variables.get(this).get(snpMapReq);
				String manifest = variables.get(this).get(manifestReq);
				String baseCommand = Files.getRunString() + " cnv.manage.Markers proj=" + projFile;
				if (Files.exists(snpMap)) {
					return baseCommand + " snps=" + snpMap;
				} else {
					return baseCommand + " snps=" + manifest + " -manifest";
				}
			}

		};
	}

	private Step generateIlluminaParseSamplesStep(double priority,
																								final Step markerPositionsStep) {
		final GenvisisWorkflow.Requirement markerPositionsStepReq = new GenvisisWorkflow.StepRequirement(markerPositionsStep);

		final GenvisisWorkflow.Requirement markerPositionsReq = new GenvisisWorkflow.FileRequirement("Parsed markerPositions file must already exist.",
																																																 proj.MARKER_POSITION_FILENAME.getValue(false,
																																																																				false));

		return new Step("Parse Illumina Sample Files", "",
										new GenvisisWorkflow.Requirement[][] {{markerPositionsStepReq,
																													 markerPositionsReq},
																													{numThreadsReq}},
										new Step[] {markerPositionsStep},
										new GenvisisWorkflow.FLAG[] {GenvisisWorkflow.FLAG.MEMORY,
																								 GenvisisWorkflow.FLAG.RUNTIME,
																								 GenvisisWorkflow.FLAG.MEMORY},
										priority) {

			@Override
			public void setNecessaryPreRunProperties(Project proj,
																							 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				String projFile = proj.MARKER_POSITION_FILENAME.getValue(false, false);
				String mkrFile = variables.get(this).get(markerPositionsReq);
				mkrFile = ext.verifyDirFormat(mkrFile);
				mkrFile = mkrFile.substring(0, mkrFile.length() - 1);
				if (!mkrFile.equals(projFile)) {
					proj.MARKER_POSITION_FILENAME.setValue(mkrFile);
				}
				int numThreads = resolveThreads(variables.get(this).get(numThreadsReq));
				maybeSetProjNumThreads(numThreads);
			}

			@Override
			public void run(Project proj,
											Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
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
			public boolean checkIfOutputExists(Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				String sampleDirectory = proj.SAMPLE_DIRECTORY.getValue(false, false);
				boolean mkrSetFile = Files.exists(proj.MARKERSET_FILENAME.getValue(false, false));
				boolean returnValue = mkrSetFile;
				returnValue = returnValue && Files.exists(sampleDirectory);
				returnValue = returnValue
											&& Files.list(sampleDirectory, Sample.SAMPLE_FILE_EXTENSION,
																		false).length > 0;
				returnValue = returnValue && proj.getSampleList() != null;
				returnValue = returnValue && proj.getSampleList().getSamples().length > 0;
				return returnValue;
			}

			@Override
			public String getCommandLine(Project proj,
																	 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
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
				int numThreads = resolveThreads(variables.get(this).get(numThreadsReq));
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

		};
	}

	private Step generateAffyParseSamplesStep(double priority) {

		final GenvisisWorkflow.Requirement markerPositionsReq = new GenvisisWorkflow.Requirement("markerPositions file must already exist.",
																																														 GenvisisWorkflow.RequirementInputType.FILE,
																																														 proj.MARKER_POSITION_FILENAME.getValue(false,
																																																																		false)) {

			@Override
			public boolean checkRequirement(String arg, Set<Step> stepSelections,
																			Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				return Files.exists(arg);
			}
		};

		return new Step("Parse Sample Files", "",
										new GenvisisWorkflow.Requirement[][] {{markerPositionsReq},
																													{numThreadsReq}},
										null,
										new GenvisisWorkflow.FLAG[] {GenvisisWorkflow.FLAG.MEMORY,
																								 GenvisisWorkflow.FLAG.RUNTIME,
																								 GenvisisWorkflow.FLAG.MEMORY},
										priority) {

			@Override
			public void setNecessaryPreRunProperties(Project proj,
																							 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				String projFile = proj.MARKER_POSITION_FILENAME.getValue(false, false);
				String mkrFile = variables.get(this).get(markerPositionsReq);
				mkrFile = ext.verifyDirFormat(mkrFile);
				mkrFile = mkrFile.substring(0, mkrFile.length() - 1);
				if (!mkrFile.equals(projFile)) {
					proj.MARKER_POSITION_FILENAME.setValue(mkrFile);
				}
				int numThreads = resolveThreads(variables.get(this).get(numThreadsReq));
				maybeSetProjNumThreads(numThreads);
			}

			@Override
			public void run(Project proj,
											Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				int numThreads = resolveThreads(variables.get(this).get(numThreadsReq));
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
			public boolean checkIfOutputExists(Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				String sampleDirectory = proj.SAMPLE_DIRECTORY.getValue(false, false);
				boolean mkrSetFile = Files.exists(proj.MARKERSET_FILENAME.getValue(false, false));
				mkrSetFile = mkrSetFile && Files.exists(sampleDirectory);
				mkrSetFile = mkrSetFile
										 && Files.list(sampleDirectory, Sample.SAMPLE_FILE_EXTENSION, false).length > 0;
				mkrSetFile = mkrSetFile && proj.getSampleList() != null;
				mkrSetFile = mkrSetFile && proj.getSampleList().getSamples().length > 0;
				return mkrSetFile;
			}

			@Override
			public String getCommandLine(Project proj,
																	 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
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

				int numThreads = resolveThreads(variables.get(this).get(numThreadsReq));
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

		};
	}

	private Step generateCreateSampleDataStep(double priority,
																						final Step parseSamplesStep) {
		final GenvisisWorkflow.Requirement parseSamplesStepReq = new GenvisisWorkflow.StepRequirement(parseSamplesStep);

		final GenvisisWorkflow.Requirement createMinimalSampleDataReq = new GenvisisWorkflow.BoolRequirement("Create a minimal SampleData.txt file from sample files",
																																																				 true);

		final String pedPreset = proj.PEDIGREE_FILENAME.getValue();

		final GenvisisWorkflow.Requirement pedigreeReq = new GenvisisWorkflow.Requirement("Either a Pedigree.dat file, or any file with a header containing all of the following elements (in any order):  \""
																																											+ ArrayUtils.toStr(MitoPipeline.PED_INPUT,
																																																				 ", ")
																																											+ "\"",
																																											GenvisisWorkflow.RequirementInputType.FILE,
																																											pedPreset) {

			@Override
			public boolean checkRequirement(String arg, Set<Step> stepSelections,
																			Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				return Files.exists(arg);
			}
		};

		// check for SampleMap only if we haven't found a pedigree
		final String sampMapPreset = Files.exists(pedPreset) ? null : getLocationOfSampleMap();

		final GenvisisWorkflow.Requirement sampMapReq = new GenvisisWorkflow.Requirement("A Sample_Map.csv file, with at least two columns having headers \""
																																										 + MitoPipeline.SAMPLEMAP_INPUT[1]
																																										 + "\" and \""
																																										 + MitoPipeline.SAMPLEMAP_INPUT[2]
																																										 + "\"",
																																										 GenvisisWorkflow.RequirementInputType.FILE,
																																										 sampMapPreset) {

			@Override
			public boolean checkRequirement(String arg, Set<Step> stepSelections,
																			Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				return Files.exists(arg);
			}
		};
		return new Step("Create SampleData.txt File", "",
										new GenvisisWorkflow.Requirement[][] {{parseSamplesStepReq},
																													{createMinimalSampleDataReq, pedigreeReq,
																													 sampMapReq}},
										new Step[] {parseSamplesStep},
										new GenvisisWorkflow.FLAG[] {},
										priority) {

			@Override
			public void setNecessaryPreRunProperties(Project proj,
																							 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				// Nothing to do
			}

			@Override
			public void run(Project proj,
											Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				Boolean minimal = Boolean.parseBoolean(variables.get(this).get(createMinimalSampleDataReq));
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
			public boolean checkIfOutputExists(Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				return Files.exists(proj.SAMPLE_DATA_FILENAME.getValue(false, false));
			}

			@Override
			public String getCommandLine(Project proj,
																	 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				String projPropFile = proj.getPropertyFilename();
				Boolean minimal = Boolean.parseBoolean(variables.get(this).get(createMinimalSampleDataReq));
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

		};
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

	private Step generateTransposeStep(double priority,
																		 final Step parseSamplesStep) {
		return new Step("Transpose Data into Marker-Dominant Files", "",
										new GenvisisWorkflow.Requirement[][] {{new GenvisisWorkflow.StepRequirement(parseSamplesStep)}},
										new Step[] {parseSamplesStep},
										new GenvisisWorkflow.FLAG[] {GenvisisWorkflow.FLAG.MEMORY}, priority) {

			@Override
			public void setNecessaryPreRunProperties(Project proj,
																							 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				// Nothing to do here
			}

			@Override
			public void run(Project proj,
											Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				proj.getLog().report("Transposing data");
				TransposeData.transposeData(proj, 2000000000, false); // compact if no LRR was provided
			}

			@Override
			public String getCommandLine(Project proj,
																	 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				String projPropFile = proj.getPropertyFilename();
				StringBuilder cmd = new StringBuilder();
				return cmd.append(Files.getRunString()).append(" cnv.manage.TransposeData -transpose proj="
																											 + projPropFile + " max=" + 2000000000)
									.toString();
			}

			@Override
			public boolean checkIfOutputExists(Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				return Files.list(proj.MARKER_DATA_DIRECTORY.getValue(false, false),
													MarkerData.MARKER_DATA_FILE_EXTENSION, false).length > 0;
			}
		};
	}

	private Step generateGCModelStep(double priority) {
		final GenvisisWorkflow.Requirement gcBaseReq = new GenvisisWorkflow.FileRequirement("A GC Base file must exist.",
																																												Resources.genome(proj.GENOME_BUILD_VERSION.getValue(),
																																																				 proj.getLog())
																																																 .getModelBase()
																																																 .get());
		final GenvisisWorkflow.Requirement gcModelOutputReq = new GenvisisWorkflow.OutputFileRequirement("GCModel output file must be specified.",
																																																		 proj.GC_MODEL_FILENAME.getValue());

		return new Step("Compute GCMODEL File", "",
										new GenvisisWorkflow.Requirement[][] {{gcBaseReq},
																													{gcModelOutputReq}},
										null, new GenvisisWorkflow.FLAG[] {}, priority) {

			@Override
			public void setNecessaryPreRunProperties(Project proj,
																							 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				String setGCOutputFile = proj.GC_MODEL_FILENAME.getValue();
				String gcOutputFile = variables.get(this).get(gcModelOutputReq);
				if (!ext.verifyDirFormat(setGCOutputFile).equals(gcOutputFile)) {
					proj.GC_MODEL_FILENAME.setValue(gcOutputFile);
				}
			}

			@Override
			public void run(Project proj,
											Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				String gcBaseFile = variables.get(this).get(gcBaseReq);
				String gcOutputFile = variables.get(this).get(gcModelOutputReq);
				org.genvisis.cnv.analysis.PennCNV.gcModel(proj, gcBaseFile, gcOutputFile, 100);
			}

			@Override
			public String getCommandLine(Project proj,
																	 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				String kvCmd = "";

				String setGCOutputFile = proj.GC_MODEL_FILENAME.getValue();
				String gcOutputFile = variables.get(this).get(gcModelOutputReq);
				if (!ext.verifyDirFormat(setGCOutputFile).equals(gcOutputFile)) {
					kvCmd += " GC_MODEL_FILENAME=" + gcOutputFile;
				}

				String projPropFile = proj.getPropertyFilename();
				StringBuilder cmd = new StringBuilder();
				if (kvCmd.length() > 0) {
					cmd.append(Files.getRunString()).append(PROJ_PROP_UPDATE_STR + projPropFile).append(kvCmd)
						 .append("\n");
				}
				String gcBaseFile = variables.get(this).get(gcBaseReq);
				return cmd.append(Files.getRunString())
									.append(" cnv.analysis.PennCNV proj=" + proj.getPropertyFilename() + " log="
													+ proj.getLog().getFilename() + " gc5base=" + gcBaseFile)
									.toString();
			}

			@Override
			public boolean checkIfOutputExists(Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				String gcOutputFile = variables.get(this).get(gcModelOutputReq);
				return Files.exists(gcOutputFile);
			}

		};
	}

	private Step generateSampleQCStep(double priority,
																		final Step parseSamplesStep) {
		final GenvisisWorkflow.Requirement parseSamplesStepReq = new GenvisisWorkflow.StepRequirement(parseSamplesStep);

		return new Step("Run Sample QC Metrics", "",
										new GenvisisWorkflow.Requirement[][] {{parseSamplesStepReq},
																													{numThreadsReq},},
										new Step[] {parseSamplesStep},
										new GenvisisWorkflow.FLAG[] {GenvisisWorkflow.FLAG.MULTITHREADED}, priority) {

			@Override
			public void setNecessaryPreRunProperties(Project proj,
																							 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				int numThreads = resolveThreads(variables.get(this).get(numThreadsReq));
				maybeSetProjNumThreads(numThreads);
			}

			@Override
			public void run(Project proj,
											Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				proj.getLog().report("Running LrrSd");
				int numThreads = proj.NUM_THREADS.getValue();
				LrrSd.init(proj, null, null, numThreads, false);
			}

			@Override
			public String getCommandLine(Project proj,
																	 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				int numThreads = resolveThreads(variables.get(this).get(numThreadsReq));
				String projPropFile = proj.getPropertyFilename();
				StringBuilder cmd = new StringBuilder();
				cmd.append(Files.getRunString()).append(" cnv.qc.LrrSd").append(" proj=")
					 .append(projPropFile)
					 .append(" ").append(PSF.Ext.NUM_THREADS_COMMAND).append(numThreads)
					 .append(" projectMarkers=TRUE");
				return cmd.toString();
			}

			@Override
			public boolean checkIfOutputExists(Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				return Files.exists(proj.SAMPLE_QC_FILENAME.getValue(false, false));
			}
		};
	}

	private Step generateMarkerQCStep(double priority,
																		final Step parseSamplesStep) {
		final GenvisisWorkflow.Requirement parseSamplesStepReq = new GenvisisWorkflow.StepRequirement(parseSamplesStep);

		final GenvisisWorkflow.Requirement exportAllReq = new GenvisisWorkflow.OptionalBoolRequirement("Export all markers in project.",
																																																	 true);

		final GenvisisWorkflow.Requirement targetMarkersReq = new GenvisisWorkflow.FileRequirement("A targetMarkers files listing the markers to QC.",
																																															 proj.TARGET_MARKERS_FILENAMES.getValue()[0]);

		return new Step("Run Marker QC Metrics", "",
										new GenvisisWorkflow.Requirement[][] {{parseSamplesStepReq},
																													{exportAllReq, targetMarkersReq},
																													{numThreadsReq}},
										new Step[] {parseSamplesStep},
										new GenvisisWorkflow.FLAG[] {GenvisisWorkflow.FLAG.MULTITHREADED},
										priority) {

			@Override
			public void setNecessaryPreRunProperties(Project proj,
																							 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				// Nothing to do here
			}

			@Override
			public void run(Project proj,
											Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				boolean allMarkers = Boolean.parseBoolean(variables.get(this).get(exportAllReq));
				String tgtFile = allMarkers ? null : variables.get(this).get(targetMarkersReq);
				int numThreads = resolveThreads(variables.get(this).get(numThreadsReq));
				MarkerMetrics.fullQC(proj, null, tgtFile, true, numThreads);
			}

			@Override
			public String getCommandLine(Project proj,
																	 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				boolean allMarkers = Boolean.parseBoolean(variables.get(this).get(exportAllReq));
				String tgtFile = variables.get(this).get(targetMarkersReq);
				int numThreads = resolveThreads(variables.get(this).get(numThreadsReq));
				return Files.getRunString() + " cnv.qc.MarkerMetrics -fullQC" + " proj="
							 + proj.getPropertyFilename() + (allMarkers ? "" : " markers=" + tgtFile) + " "
							 + PSF.Ext.NUM_THREADS_COMMAND + numThreads;
			}

			@Override
			public boolean checkIfOutputExists(Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				String markerMetricsFile = proj.MARKER_METRICS_FILENAME.getValue();
				return Files.exists(markerMetricsFile);
			}

		};
	}

	private Step generateSexChecksStep(double priority,
																		 final Step parseSamplesStep,
																		 final Step sampleDataStep,
																		 final Step transposeStep,
																		 final Step sampleQCStep) {
		final GenvisisWorkflow.Requirement parseSamplesStepReq = new GenvisisWorkflow.StepRequirement(parseSamplesStep);
		final GenvisisWorkflow.Requirement sampleDataStepReq = new GenvisisWorkflow.StepRequirement(sampleDataStep);
		final GenvisisWorkflow.Requirement transposeStepReq = new GenvisisWorkflow.StepRequirement(transposeStep);
		final GenvisisWorkflow.Requirement sampleQCStepReq = new GenvisisWorkflow.StepRequirement(sampleQCStep);
		final GenvisisWorkflow.Requirement addToSampleDataReq = new GenvisisWorkflow.OptionalBoolRequirement("Add Estimated Sex to Sample Data",
																																																				 true);

		final GenvisisWorkflow.Requirement noCrossHybeReq = new GenvisisWorkflow.BoolRequirement("Use only X and Y chromosome R values to identify sex discriminating markers",
																																														 false);

		final GenvisisWorkflow.Requirement oneHittersReq = new GenvisisWorkflow.FileRequirement("List of markers that do not cross hybridize",
																																														MarkerBlastQC.defaultOneHitWondersFilename(proj.BLAST_ANNOTATION_FILENAME.getValue()));

		final GenvisisWorkflow.Requirement blastVCFReq = new GenvisisWorkflow.FileRequirement("BLAST annotation VCF to generate list of markers that do not cross hybridize from",
																																													proj.BLAST_ANNOTATION_FILENAME.getValue());

		return new Step("Run Sex Checks", "",
										new GenvisisWorkflow.Requirement[][] {{parseSamplesStepReq},
																													{sampleDataStepReq},
																													{transposeStepReq},
																													{sampleQCStepReq},
																													{addToSampleDataReq},
																													{noCrossHybeReq,
																													 oneHittersReq,
																													 blastVCFReq}},
										new Step[] {parseSamplesStep, sampleDataStep,
																transposeStep, sampleQCStep},
										new GenvisisWorkflow.FLAG[] {},
										priority) {

			@Override
			public void setNecessaryPreRunProperties(Project proj,
																							 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				// Nothing to do here
			}

			@Override
			public void run(Project proj,
											Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				proj.getLog().report("Running SexCheck");
				boolean addToSampleData = Boolean.parseBoolean(variables.get(this).get(addToSampleDataReq));
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
			public String getCommandLine(Project proj,
																	 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				String projPropFile = proj.getPropertyFilename();
				StringBuilder cmd = new StringBuilder();
				boolean addToSampleData = Boolean.parseBoolean(variables.get(this).get(addToSampleDataReq));
				String discriminatingMarkersFile;
				if (Boolean.parseBoolean(variables.get(this).get(noCrossHybeReq))) {
					discriminatingMarkersFile = null;
				} else {
					discriminatingMarkersFile = variables.get(this).get(oneHittersReq);
					if (!Files.exists(discriminatingMarkersFile)) {
						cmd.append(Files.getRunString()).append(" cnv.qc.MarkerBlastQC proj=" + projPropFile
																										+ " blastVCF="
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
			public boolean checkIfOutputExists(Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				return Files.exists(proj.SEXCHECK_RESULTS_FILENAME.getValue());
			}

		};
	}

	private Step generateABLookupStep(double priority,
																		final Step parseSamplesStep) {
		final GenvisisWorkflow.Requirement parseSamplesStepReq = new GenvisisWorkflow.StepRequirement(parseSamplesStep);
		return new Step("Generate AB Lookup File", "",
										new GenvisisWorkflow.Requirement[][] {{parseSamplesStepReq}},
										new Step[] {parseSamplesStep},
										new GenvisisWorkflow.FLAG[] {GenvisisWorkflow.FLAG.RUNTIME}, priority) {

			@Override
			public void setNecessaryPreRunProperties(Project proj,
																							 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				// Nothing to do here
			}

			@Override
			public void run(Project proj,
											Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				String filename = proj.PROJECT_DIRECTORY.getValue()
													+ ext.addToRoot(ABLookup.DEFAULT_AB_FILE, "_parsed");
				ABLookup.parseABLookup(proj, ABSource.VCF, filename);

				if (ABLookup.fillInMissingAlleles(proj, filename, proj.getLocationOfSNP_Map(true), false)) {
					ABLookup.applyABLookupToFullSampleFiles(proj, filename);
				} else {
					setFailed("Failed to fill in missing alleles - please check log for more info.");
				}
			}

			@Override
			public String getCommandLine(Project proj,
																	 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
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
			public boolean checkIfOutputExists(Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				return Files.exists(proj.AB_LOOKUP_FILENAME.getValue(false, false));
			}
		};

	}

	private Step generatePlinkExportStep(double priority,
																			 final Step parseSamplesStep) {
		final GenvisisWorkflow.Requirement parseSamplesStepReq = new GenvisisWorkflow.StepRequirement(parseSamplesStep);
		final GenvisisWorkflow.Requirement pedigreeRequirement = new GenvisisWorkflow.FileRequirement("A pedigree.dat file must exist.",
																																																	proj.PEDIGREE_FILENAME.getValue(false,
																																																																	false));
		final GenvisisWorkflow.Requirement createPedigreeRequirement = new GenvisisWorkflow.BoolRequirement("Create a minimal pedigree.dat file [will pull information from SexChecks step results].",
																																																				false);

		return new Step("Create PLINK Files", "",
										new GenvisisWorkflow.Requirement[][] {{parseSamplesStepReq},
																													{pedigreeRequirement,
																													 createPedigreeRequirement}},
										new Step[] {parseSamplesStep},
										new GenvisisWorkflow.FLAG[] {GenvisisWorkflow.FLAG.MEMORY}, priority) {

			@Override
			public void setNecessaryPreRunProperties(Project proj,
																							 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				if (!Boolean.parseBoolean(variables.get(this).get(createPedigreeRequirement))) {
					String projPedFile = proj.PEDIGREE_FILENAME.getValue(false, false);
					String pedFile = variables.get(this).get(pedigreeRequirement);
					if (!pedFile.equals(projPedFile)) {
						proj.PEDIGREE_FILENAME.setValue(pedFile);
					}
				}
			}

			@Override
			public void run(Project proj,
											Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				if (Boolean.parseBoolean(variables.get(this).get(createPedigreeRequirement))) {
					proj.getLog().report("Creating Pedigree File");
					Pedigree.build(proj, null, null, false);
				}
				if (!Files.exists(proj.PEDIGREE_FILENAME.getValue())) {
					setFailed("Creation of Pedigree file in [Create/Run PLINK Files] step failed.");
					return;
				}

				proj.getLog().report("Running PLINK");

				boolean create = PlinkData.saveGenvisisToPlinkBedSet(proj, "plink/plink", null, null, -1,
																														 true);
				if (!create) {
					setFailed("Creation of initial PLINK files failed.");
				}
				proj.PLINK_DIR_FILEROOTS.addValue(proj.PROJECT_DIRECTORY.getValue() + "plink/plink");
			}

			@Override
			public String getCommandLine(Project proj,
																	 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
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
						 .append(projPropFile)
						 .append("\n");
				}
				cmd.append(Files.getRunString())
					 .append(" cnv.manage.PlinkData -genvisisToBed plinkdata=plink/plink gcthreshold=-1 proj=")
					 .append(proj.getPropertyFilename());
				return cmd.toString();
			}

			@Override
			public boolean checkIfOutputExists(Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				boolean plinkFilesExist = Files.checkAllFiles(proj.PROJECT_DIRECTORY.getValue() + "plink/",
																											PSF.Plink.getPlinkBedBimFamSet("plink"),
																											false, log);
				boolean pedGenerated = Boolean.parseBoolean(variables.get(this)
																														 .get(createPedigreeRequirement));
				boolean pedCheck = pedGenerated ? Files.exists(proj.PEDIGREE_FILENAME.getValue()) : true;
				return plinkFilesExist && pedCheck;
			}

		};
	}

	private Step generateGwasQCStep(double priority,
																	Step plinkExportStep) {
		// TODO: Move Ancestry to its own step
		final GenvisisWorkflow.Requirement plinkExportStepReq = new GenvisisWorkflow.StepRequirement(plinkExportStep);
		final GenvisisWorkflow.Requirement genomeForRelatedsReq = new GenvisisWorkflow.OptionalBoolRequirement("Keep genome info for unrelateds only",
																																																					 false);
		final GenvisisWorkflow.Requirement skipAncestryReq = new GenvisisWorkflow.BoolRequirement("Skip ancestry checks",
																																															false);
		final GenvisisWorkflow.Requirement putativeWhitesReq = new GenvisisWorkflow.FileRequirement("File with FID/IID pairs of putative white samples",
																																																"");
		final GenvisisWorkflow.Requirement hapMapFoundersReq = new GenvisisWorkflow.FileRequirement("PLINK root of HapMap founders",
																																																Ancestry.DEFAULT_HAPMAP_PLINKROOT) {
			@Override
			public boolean checkRequirement(String arg, Set<Step> stepSelections,
																			Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				String hapMapPlinkRoot = arg;
				int dotIndex = hapMapPlinkRoot.lastIndexOf('.');
				if (dotIndex > 0
						&& PSF.Plink.getPlinkBedBimFamSet("").contains(hapMapPlinkRoot.substring(dotIndex))) {
					hapMapPlinkRoot = hapMapPlinkRoot.substring(0, dotIndex);
				}
				return Files.checkAllFiles("", PSF.Plink.getPlinkBedBimFamSet(hapMapPlinkRoot), false, log);
			}
		};

		return new Step("Run GWAS QC", "",
										new GenvisisWorkflow.Requirement[][] {{plinkExportStepReq},
																													{genomeForRelatedsReq},
																													{skipAncestryReq,
																													 putativeWhitesReq,
																													 hapMapFoundersReq}},
										new Step[] {plinkExportStep}, new GenvisisWorkflow.FLAG[] {},
										priority) {

			@Override
			public void setNecessaryPreRunProperties(Project proj,
																							 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				// not needed for step
			}

			@Override
			public void run(Project proj,
											Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				String dir = getPlinkDir();
				boolean keepUnrelatedsOnly = Boolean.parseBoolean(variables.get(this)
																																	 .get(genomeForRelatedsReq));
				boolean skipAncestry = Boolean.parseBoolean(variables.get(this).get(skipAncestryReq));
				String putativeWhites = variables.get(this).get(putativeWhitesReq);
				String hapMapPlinkRoot = variables.get(this).get(hapMapFoundersReq);
				int hapMapDotIndex = hapMapPlinkRoot.lastIndexOf('.');
				if (hapMapDotIndex > 0 && PSF.Plink.getPlinkBedBimFamSet("")
																					 .contains(hapMapPlinkRoot.substring(hapMapDotIndex))) {
					hapMapPlinkRoot = hapMapPlinkRoot.substring(0, hapMapDotIndex);
				}
				RelationAncestryQc.fullGamut(dir, null, keepUnrelatedsOnly, proj.getLog());
				if (new File(dir + Qc.QC_DIR + RelationAncestryQc.GENOME_DIR + "plink.genome").exists()) {
					proj.GENOME_CLUSTER_FILENAME.setValue(dir + Qc.QC_DIR + RelationAncestryQc.GENOME_DIR
																								+ "plink.genome");
					proj.saveProperties();
				}
				if (!keepUnrelatedsOnly) {
					new PlinkMendelianChecker(proj).run();
				}
				if (!skipAncestry) {
					String ancestryDir = dir + Qc.QC_DIR + RelationAncestryQc.ANCESTRY_DIR;
					Ancestry.runPipeline(ancestryDir, putativeWhites, hapMapPlinkRoot, proj,
															 new Logger(ancestryDir + "ancestry.log"));
				}
			}

			@Override
			public String getCommandLine(Project proj,
																	 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				String dir = getPlinkDir();
				boolean keepUnrelatedsOnly = Boolean.parseBoolean(variables.get(this)
																																	 .get(genomeForRelatedsReq));
				boolean skipAncestry = Boolean.parseBoolean(variables.get(this).get(skipAncestryReq));
				String putativeWhites = variables.get(this).get(putativeWhitesReq);
				String hapMapPlinkRoot = variables.get(this).get(hapMapFoundersReq);
				int hapMapDotIndex = hapMapPlinkRoot.lastIndexOf('.');
				if (hapMapDotIndex > 0 && PSF.Plink.getPlinkBedBimFamSet("")
																					 .contains(hapMapPlinkRoot.substring(hapMapDotIndex))) {
					hapMapPlinkRoot = hapMapPlinkRoot.substring(0, hapMapDotIndex);
				}

				String command = Files.getRunString() + " gwas.Qc dir=" + dir
												 + " keepGenomeInfoForRelatedsOnly=" + keepUnrelatedsOnly;
				command += "\n";
				command += Files.getRunString() + " " + PROJ_PROP_UPDATE_STR + proj.getPropertyFilename();
				command += " " + proj.GENOME_CLUSTER_FILENAME.getName() + "=" + dir
									 + Qc.QC_DIR + RelationAncestryQc.GENOME_DIR + "plink.genome";
				if (!keepUnrelatedsOnly) {
					command += "\n";
					command += Files.getRunString() + " org.genvisis.gwas.PlinkMendelianChecker proj="
										 + proj.getPropertyFilename();
				}
				if (!skipAncestry) {
					String ancestryDir = dir + Qc.QC_DIR + RelationAncestryQc.ANCESTRY_DIR;
					command += "\n";
					command += Files.getRunString() + " gwas.Ancestry -runPipeline dir=" + ancestryDir;
					command += " putativeWhites=" + putativeWhites;
					command += " proj=" + proj.getPropertyFilename();
					command += " hapMapPlinkRoot=" + hapMapPlinkRoot;
					command += " log=" + ancestryDir + "ancestry.log";
				}
				return command;
			}

			@Override
			public boolean checkIfOutputExists(Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				String dir = getPlinkDir();
				boolean allExist = true;
				boolean skipAncestry = Boolean.parseBoolean(variables.get(this).get(skipAncestryReq));
				folders: for (int i = 0; i < org.genvisis.gwas.RelationAncestryQc.FOLDERS_CREATED.length; i++) {
					for (int j = 0; j < org.genvisis.gwas.RelationAncestryQc.FILES_CREATED[i].length; j++) {
						if (!Files.exists(dir + org.genvisis.gwas.RelationAncestryQc.FOLDERS_CREATED[i]
															+ org.genvisis.gwas.RelationAncestryQc.FILES_CREATED[i][j])) {
							allExist = false;
							break folders;
						}
					}
				}
				if (!skipAncestry
						&& (!Files.exists(dir + Qc.QC_DIR + RelationAncestryQc.ANCESTRY_DIR + "freqsByRace.xln")
								|| !Files.exists(dir + Qc.QC_DIR + RelationAncestryQc.ANCESTRY_DIR
																 + "raceImputations.mds"))) {
					allExist = false;
				}
				return allExist;
			}

		};
	}

	private Step generateMosaicArmsStep(double priority,
																			final Step parseSamplesStep) {
		final GenvisisWorkflow.Requirement parseSamplesStepReq = new GenvisisWorkflow.StepRequirement(parseSamplesStep);
		return new Step("Create Mosaic Arms File", "",
										new GenvisisWorkflow.Requirement[][] {{parseSamplesStepReq},
																													{numThreadsReq}},
										new Step[] {parseSamplesStep},
										new GenvisisWorkflow.FLAG[] {GenvisisWorkflow.FLAG.MULTITHREADED}, priority) {

			@Override
			public void setNecessaryPreRunProperties(Project proj,
																							 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {

				int numThreads = resolveThreads(variables.get(this).get(numThreadsReq));
				maybeSetProjNumThreads(numThreads);
			}

			@Override
			public void run(Project proj,
											Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				Mosaicism.findOutliers(proj);
			}

			@Override
			public String getCommandLine(Project proj,
																	 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				String kvCmd = "";


				int numThreads = resolveThreads(variables.get(this).get(numThreadsReq));
				if (numThreads != proj.NUM_THREADS.getValue()) {
					kvCmd += " " + proj.NUM_THREADS.getName() + "=" + numThreads;
				}

				String projPropFile = proj.getPropertyFilename();
				StringBuilder cmd = new StringBuilder();
				if (kvCmd.length() > 0) {
					cmd.append(Files.getRunString()).append(PROJ_PROP_UPDATE_STR + projPropFile).append(kvCmd)
						 .append("\n");
				}
				return cmd.append(Files.getRunString())
									.append(" cnv.analysis.Mosaicism proj=" + proj.getPropertyFilename()).toString();
			}

			@Override
			public boolean checkIfOutputExists(Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				return Files.exists(proj.RESULTS_DIRECTORY.getValue(false, false) + "Mosaicism.xln");
			}
		};
	}

	private Step generateAnnotateSampleDataStep(double priority,
																							final Step sampleQCStep,
																							final Step createSampleDataStep,
																							final Step gwasQCStep) {
		final GenvisisWorkflow.Requirement sampleQCStepReq = new GenvisisWorkflow.StepRequirement(sampleQCStep);
		final GenvisisWorkflow.Requirement createSampleDataStepReq = new GenvisisWorkflow.StepRequirement(createSampleDataStep);
		final GenvisisWorkflow.Requirement skipIDingDuplicatesReq = new GenvisisWorkflow.BoolRequirement("Skip identifying duplicates",
																																																		 false);
		final GenvisisWorkflow.Requirement gwasQCStepReq = new GenvisisWorkflow.StepRequirement(gwasQCStep);
		final GenvisisWorkflow.Requirement notGcCorrectedLrrSdReq = new GenvisisWorkflow.BoolRequirement("Do not use GC corrected LRR SD?",
																																																		 false);
		final GenvisisWorkflow.Requirement gcCorrectedLrrSdReq = new GenvisisWorkflow.Requirement("GC Corrected LRR SD must exist in Sample QC File",
																																															GenvisisWorkflow.RequirementInputType.NONE) {

			@Override
			public boolean checkRequirement(String arg, Set<Step> stepSelections,
																			Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				String sampleQCFile = proj.SAMPLE_QC_FILENAME.getValue();
				return Files.exists(sampleQCFile)
							 && ext.indexOfStr("LRR_SD_Post_Correction",
																 Files.getHeaderOfFile(sampleQCFile, proj.getLog())) != -1;
			}

		};
		final GenvisisWorkflow.Requirement lrrSdThresholdReq = new GenvisisWorkflow.DoubleRequirement("LRR SD Threshold",
																																																	proj.LRRSD_CUTOFF.getValue(),
																																																	proj.LRRSD_CUTOFF.getMinValue(),
																																																	proj.LRRSD_CUTOFF.getMaxValue());

		final GenvisisWorkflow.Requirement callrateThresholdReq = new GenvisisWorkflow.DoubleRequirement("Callrate Threshold",
																																																		 proj.SAMPLE_CALLRATE_THRESHOLD.getValue(),
																																																		 proj.SAMPLE_CALLRATE_THRESHOLD.getMinValue(),
																																																		 proj.SAMPLE_CALLRATE_THRESHOLD.getMaxValue());
		final GenvisisWorkflow.Requirement numQReq = new GenvisisWorkflow.IntRequirement("Number of Quantiles to Generate",
																																										 10, 1,
																																										 Integer.MAX_VALUE);
		final GenvisisWorkflow.Requirement replaceFIDIIDReq = new GenvisisWorkflow.OptionalBoolRequirement("Replace FID and IID with data from Pedigree",
																																																			 false);
		return new Step("Annotate Sample Data File", "",
										new GenvisisWorkflow.Requirement[][] {{sampleQCStepReq},
																													{createSampleDataStepReq},
																													{skipIDingDuplicatesReq, gwasQCStepReq},
																													{notGcCorrectedLrrSdReq,
																													 gcCorrectedLrrSdReq},
																													{lrrSdThresholdReq},
																													{callrateThresholdReq},
																													{numQReq},
																													{replaceFIDIIDReq}},
										new Step[] {sampleQCStep,
																createSampleDataStep,
																gwasQCStep},
										new GenvisisWorkflow.FLAG[] {}, priority) {

			@Override
			public void setNecessaryPreRunProperties(Project proj,
																							 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
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
			public void run(Project proj,
											Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				boolean checkDuplicates = !Boolean.parseBoolean(variables.get(this)
																																 .get(skipIDingDuplicatesReq));
				String duplicatesSetFile = null;
				if (checkDuplicates) {
					String dir = "plink/";
					duplicatesSetFile = proj.PROJECT_DIRECTORY.getValue() + dir
															+ Qc.QC_DIR + RelationAncestryQc.GENOME_DIR
															+ "plink.genome_duplicatesSet.dat";
				}
				boolean gcCorrectedLrrSd = !Boolean.parseBoolean(variables.get(this)
																																	.get(notGcCorrectedLrrSdReq));
				int numQ = Integer.parseInt(variables.get(this).get(numQReq));
				boolean correctFidIids = Boolean.parseBoolean(variables.get(this).get(replaceFIDIIDReq));
				SampleQC.parseAndAddToSampleData(proj, numQ, 0, false, gcCorrectedLrrSd, duplicatesSetFile,
																				 correctFidIids);
			}

			@Override
			public String getCommandLine(Project proj,
																	 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {

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
					String dir = "plink/";
					duplicatesSetFile = proj.PROJECT_DIRECTORY.getValue() + dir
															+ Qc.QC_DIR + RelationAncestryQc.GENOME_DIR
															+ "plink.genome_duplicatesSet.dat";
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
					cmd.append(Files.getRunString()).append(PROJ_PROP_UPDATE_STR + projPropFile).append(kvCmd)
						 .append("\n");
				}
				cmd.append(Files.getRunString())
					 .append(" cnv.qc.SampleQC proj=" + projPropFile + " numQ=" + numQ
									 + " justQuantiles=false"
									 + " gcCorrectedLrrSd=" + gcCorrectedLrrSd + " duplicatesSetFile="
									 + duplicatesSetFile + " correctFidIids=" + correctFidIids);
				return cmd.toString();
			}

			@Override
			public boolean checkIfOutputExists(Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
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

		};
	}

	private Step generateCreatePCsStep(double priority,
																		 Step transposeStep) {
		// FIXME http://genvisis.org/MitoPipeline/#illumina_marker_lists has illumina markers.. this
		// should be linked to, or
		// these steps split or something...
		final GenvisisWorkflow.Requirement transposeStepReq = new GenvisisWorkflow.StepRequirement(transposeStep);
		final GenvisisWorkflow.Requirement medianMarkersReq = new GenvisisWorkflow.FileRequirement("MedianMarkers file must exist.",
																																															 "");
		final GenvisisWorkflow.Requirement lrrSdThresholdReq = new GenvisisWorkflow.DoubleRequirement("LRR SD threshold to filter samples.",
																																																	proj.LRRSD_CUTOFF.getValue(),
																																																	proj.LRRSD_CUTOFF.getMinValue(),
																																																	proj.LRRSD_CUTOFF.getMaxValue());
		final GenvisisWorkflow.Requirement callrateThresholdReq = new GenvisisWorkflow.DoubleRequirement("Call rate threshold to filter markers.",
																																																		 MitoPipeline.DEFAULT_MKR_CALLRATE_FILTER,
																																																		 0.0,
																																																		 1.0);
		final GenvisisWorkflow.Requirement qcPassingOnlyReq = new GenvisisWorkflow.OptionalBoolRequirement("Compute PCs with samples passing QC only",
																																																			 true);
		final GenvisisWorkflow.Requirement imputeNaNs = new GenvisisWorkflow.OptionalBoolRequirement("Impute mean value for NaN",
																																																 true);
		final GenvisisWorkflow.Requirement recomputeLrrPCMarkersReq = new GenvisisWorkflow.OptionalBoolRequirement("Should recompute Log-R ratio for PC markers?",
																																																							 true);
		final GenvisisWorkflow.Requirement recomputeLrrMedianMarkersReq = new GenvisisWorkflow.OptionalBoolRequirement("Should recompute Log-R ratio for median markers?",
																																																									 true);
		final GenvisisWorkflow.Requirement homozygousOnlyReq = new GenvisisWorkflow.OptionalBoolRequirement("Homozygous only?",
																																																				true);
		final GenvisisWorkflow.Requirement gcRegressionDistanceReq = new GenvisisWorkflow.IntRequirement("Regression distance for the GC adjustment",
																																																		 GcAdjustor.DEFAULT_REGRESSION_DISTANCE[0],
																																																		 1,
																																																		 Integer.MAX_VALUE);
		final GenvisisWorkflow.Requirement pcSelectionSamplesReq = new GenvisisWorkflow.OptionalFileRequirement("A file listing a subset of samples (DNA ID) to use for determining optimal PC selection, typically a list of unrelated and single race samples. If a list is not provided, only samples passing sample qc thresholds will be used.",
																																																						"");
		final GenvisisWorkflow.Requirement externalBetaFileReq = new GenvisisWorkflow.OptionalFileRequirement("An external beta file to optimize PC selection.",
																																																					"");

		return new Step("Create Principal Components File and Mitochondrial Copy-Number Estimates File",
										"",
										new GenvisisWorkflow.Requirement[][] {{transposeStepReq},
																													{medianMarkersReq},
																													{lrrSdThresholdReq},
																													{callrateThresholdReq},
																													{qcPassingOnlyReq},
																													{imputeNaNs},
																													{recomputeLrrPCMarkersReq},
																													{recomputeLrrMedianMarkersReq},
																													{homozygousOnlyReq},
																													{gcRegressionDistanceReq},
																													{numThreadsReq},
																													{pcSelectionSamplesReq},
																													{externalBetaFileReq}},
										new Step[] {transposeStep},
										new GenvisisWorkflow.FLAG[] {GenvisisWorkflow.FLAG.MULTITHREADED}, priority) {

			@Override
			public void setNecessaryPreRunProperties(Project proj,
																							 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
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
			public void run(Project proj,
											Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
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
				int regressionDistance = Integer.parseInt(variables.get(this).get(gcRegressionDistanceReq));
				int numComponents = MitoPipeline.DEFAULT_NUM_COMPONENTS;
				int numThreads = resolveThreads(variables.get(this).get(numThreadsReq));
				String outputBase = MitoPipeline.FILE_BASE;

				String betaOptFile = variables.get(this).get(pcSelectionSamplesReq);
				String betaFile = variables.get(this).get(externalBetaFileReq);

				boolean markerQC = true;
				double[] pvalOpt = MitoPipeline.DEFAULT_PVAL_OPTS;
				String pedFile = null;
				String useFile = null;
				boolean sampLrr = true;
				boolean plot = false;
				int retCode = PCAPrep.prepPCA(proj, numThreads, outputBase, markerQC, markerCallRateFilter,
																			useFile, proj.getSampleList(), proj.getLog());
				if (retCode == PCAPrep.SUCCESS_CODE) {
					MitoPipeline.estimateMtDNACN(proj, numThreads, medianMarkers, numComponents, outputBase,
																			 homozygousOnly, markerCallRateFilter, betaOptFile, pedFile,
																			 recomputeLRRPCs, recomputeLRRMedian, sampLrr,
																			 imputeMeanForNaN,
																			 gcCorrect, bpGcModel, regressionDistance,
																			 proj.GENOME_BUILD_VERSION.getValue(), pvalOpt, betaFile,
																			 plot,
																			 false, proj.getLog());
				} else {
					setFailed(PCAPrep.errorMessage(retCode));
				}
			}

			@Override
			public String getCommandLine(Project proj,
																	 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
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
				int regressionDistance = Integer.parseInt(variables.get(this).get(gcRegressionDistanceReq));
				int numComponents = MitoPipeline.DEFAULT_NUM_COMPONENTS;
				int numThreads = resolveThreads(variables.get(this).get(numThreadsReq));
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
					 .append(" LRRSD=").append(lrrSD).append(" markerCallRate=").append(markerCallRateFilter)
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
			public boolean checkIfOutputExists(Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
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
		};
	}

	private Step generatePFBStep(double priority,
															 final Step parseSamplesStep) {
		final GenvisisWorkflow.Requirement parseSamplesStepReq = new GenvisisWorkflow.StepRequirement(parseSamplesStep);
		final GenvisisWorkflow.Requirement sampleSubsetReq = new GenvisisWorkflow.FileRequirement("A Sample subset file must exist.",
																																															proj.SAMPLE_SUBSET_FILENAME.getValue());
		String defaultOutputFile;
		if (Files.exists(proj.SAMPLE_SUBSET_FILENAME.getValue())) {
			defaultOutputFile = ext.rootOf(proj.SAMPLE_SUBSET_FILENAME.getValue()) + ".pfb";
		} else {
			defaultOutputFile = proj.CUSTOM_PFB_FILENAME.getValue();
		}
		final GenvisisWorkflow.Requirement outputFileReq = new GenvisisWorkflow.OutputFileRequirement("PFB (population BAF) output file must be specified.",
																																																	defaultOutputFile);
		return new Step("Compute Population BAF files", "",
										new GenvisisWorkflow.Requirement[][] {{parseSamplesStepReq,
																													 sampleSubsetReq},
																													{outputFileReq}},
										new Step[] {parseSamplesStep},
										new GenvisisWorkflow.FLAG[] {}, priority) {

			@Override
			public void setNecessaryPreRunProperties(Project proj,
																							 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
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
			public void run(Project proj,
											Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				org.genvisis.cnv.analysis.PennCNV.populationBAF(proj);
			}

			@Override
			public String getCommandLine(Project proj,
																	 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
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
					cmd.append(Files.getRunString()).append(PROJ_PROP_UPDATE_STR + projPropFile).append(kvCmd)
						 .append("\n");
				}
				return cmd.append(Files.getRunString())
									.append(" cnv.analysis.PennCNV -pfb proj=" + proj.getPropertyFilename() + " log="
													+ proj.getLog().getFilename())
									.toString();
			}

			@Override
			public boolean checkIfOutputExists(Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				String subSampFile = variables.get(this).get(sampleSubsetReq);
				String pfbOutputFile = variables.get(this).get(outputFileReq);
				return Files.exists(pfbOutputFile) || Files.exists(ext.rootOf(subSampFile) + ".pfb");
			}
		};
	}

	private Step generateSexCentroidsStep(double priority,
																				Step gcModelStep) {
		final GenvisisWorkflow.Requirement gcModelStepReq = new GenvisisWorkflow.StepRequirement(gcModelStep);
		final GenvisisWorkflow.Requirement gcModelFileReq = new GenvisisWorkflow.FileRequirement("Full GC Model File.",
																																														 proj.GC_MODEL_FILENAME.getValue());
		return new Step("Create Sex-Specific Centroids; Filter PFB and GCMODEL Files",
										"",
										new GenvisisWorkflow.Requirement[][] {{gcModelStepReq,
																													 gcModelFileReq},
																													{numThreadsReq},},
										new Step[] {gcModelStep},
										new GenvisisWorkflow.FLAG[] {GenvisisWorkflow.FLAG.RUNTIME,
																								 GenvisisWorkflow.FLAG.MULTITHREADED,},
										priority) {

			@Override
			public void setNecessaryPreRunProperties(Project proj,
																							 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {

				int numThreads = resolveThreads(variables.get(this).get(numThreadsReq));
				maybeSetProjNumThreads(numThreads);
			}

			@Override
			public void run(Project proj,
											Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
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

				int numThreads = resolveThreads(variables.get(this).get(numThreadsReq));
				String gcModelFile = variables.get(this).get(gcModelFileReq);
				Centroids.computeSexSpecificCentroids(proj, new String[] {malePFB, femalePFB},
																							new String[] {centFilePathM, centFilePathF},
																							numThreads);

				AnalysisFormats.filterSexSpecificGCModel(proj, gcModelFile, newGCFile);
			}

			@Override
			public String getCommandLine(Project proj,
																	 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {

				int numThreads = resolveThreads(variables.get(this).get(numThreadsReq));
				String mainCmd = Files.getRunString() + " cnv.filesys.Centroids proj="
												 + proj.getPropertyFilename() + " -sexSpecific "
												 + PSF.Ext.NUM_THREADS_COMMAND
												 + numThreads;
				String gcModelFile = variables.get(this).get(gcModelFileReq);
				String gcCmd = Files.getRunString() + " cnv.analysis.AnalysisFormats proj="
											 + proj.getPropertyFilename() + " gcmodel=" + gcModelFile;
				return mainCmd + "\n" + gcCmd;
			}

			@Override
			public boolean checkIfOutputExists(Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
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
		};
	}

	private Step generateCNVStep(double priority, Step pfbStep,
															 Step gcModelStep) {
		final GenvisisWorkflow.Requirement hmmFile = new GenvisisWorkflow.FileRequirement("Hidden Markov Model File Must Exist",
																																											proj.HMM_FILENAME.getValue());
		final GenvisisWorkflow.Requirement pfbStepReq = new GenvisisWorkflow.StepRequirement(pfbStep);
		final GenvisisWorkflow.Requirement pfbFileReq = new GenvisisWorkflow.FileRequirement("PFB File Must Exist",
																																												 proj.CUSTOM_PFB_FILENAME.getValue());
		final GenvisisWorkflow.Requirement gcModelStepReq = new GenvisisWorkflow.StepRequirement(gcModelStep);
		final GenvisisWorkflow.Requirement gcModelFileReq = new GenvisisWorkflow.FileRequirement("GCMODEL File Must Exist",
																																														 proj.GC_MODEL_FILENAME.getValue());
		final GenvisisWorkflow.Requirement outputFileReq = new GenvisisWorkflow.OutputFileRequirement("Output filename.",
																																																	"cnvs/genvisis.cnv") {
			@Override
			public boolean checkRequirement(String arg, Set<Step> stepSelections,
																			Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				return super.checkRequirement(proj.PROJECT_DIRECTORY.getValue() + arg, stepSelections,
																			variables);
			}
		};

		return new Step("Call CNVs", "",
										new GenvisisWorkflow.Requirement[][] {{hmmFile},
																													{pfbStepReq,
																													 pfbFileReq},
																													{gcModelStepReq,
																													 gcModelFileReq},
																													{numThreadsReq},
																													{outputFileReq}},
										new Step[] {pfbStep, gcModelStep},
										new GenvisisWorkflow.FLAG[] {GenvisisWorkflow.FLAG.MEMORY,
																								 GenvisisWorkflow.FLAG.MULTITHREADED,},
										priority) {

			@Override
			public void setNecessaryPreRunProperties(Project proj,
																							 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
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
				int numThreads = resolveThreads(variables.get(this).get(numThreadsReq));
				maybeSetProjNumThreads(numThreads);
			}

			@Override
			public void run(Project proj,
											Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				int numThreads = resolveThreads(variables.get(this).get(numThreadsReq));
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
			public String getCommandLine(Project proj,
																	 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
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

				int numThreads = resolveThreads(variables.get(this).get(numThreadsReq));
				if (numThreads != proj.NUM_THREADS.getValue()) {
					kvCmd += " " + proj.NUM_THREADS.getName() + "=" + numThreads;
				}
				String projPropFile = proj.getPropertyFilename();
				StringBuilder cmd = new StringBuilder();
				if (kvCmd.length() > 0) {
					cmd.append(Files.getRunString()).append(PROJ_PROP_UPDATE_STR + projPropFile).append(kvCmd)
						 .append("\n");
				}
				return cmd.append(Files.getRunString()).append(" cnv.hmm.CNVCaller proj=" + projPropFile)
									.append(" out=" + variables.get(this).get(outputFileReq)).append(" ")
									.append(PSF.Ext.NUM_THREADS_COMMAND).append(numThreads).toString();
			}

			@Override
			public boolean checkIfOutputExists(Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				String output = variables.get(this).get(outputFileReq);
				return Files.exists(proj.PROJECT_DIRECTORY.getValue() + output);
			}
		};
	}

	private Step generateShadowStep(double priority,
																	final Step parseSamplesStep) {
		final GenvisisWorkflow.Requirement parseSamplesStepReq = new GenvisisWorkflow.StepRequirement(parseSamplesStep);
		final GenvisisWorkflow.Requirement numPCsReq = new GenvisisWorkflow.IntRequirement("Number of principal components for correction.",
																																											 MitoPipeline.DEFAULT_NUM_COMPONENTS,
																																											 1,
																																											 Integer.MAX_VALUE);
		final GenvisisWorkflow.Requirement outputBaseReq = new GenvisisWorkflow.OutputFileRequirement("Output file path (relative to project directory) and baseName for principal components correction files",
																																																	MitoPipeline.FILE_BASE) {
			@Override
			public boolean checkRequirement(String arg, Set<Step> stepSelections,
																			Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				String outputBase = proj.PROJECT_DIRECTORY.getValue() + arg;
				String finalReport = outputBase + PCA.FILE_EXTs[0];
				return super.checkRequirement(finalReport, stepSelections, variables);
			}
		};
		final GenvisisWorkflow.Requirement callrateReq = new GenvisisWorkflow.DoubleRequirement("Call-rate filter for determining high-quality markers",
																																														MitoPipeline.DEFAULT_MKR_CALLRATE_FILTER,
																																														0.0,
																																														1.0);
		final GenvisisWorkflow.Requirement recomputeLrrReq = new GenvisisWorkflow.OptionalBoolRequirement("Re-compute Log-R Ratio values? (usually false if LRRs already exist)",
																																																			false);
		final GenvisisWorkflow.Requirement tempDirReq = new GenvisisWorkflow.OptionalFileRequirement("Temporary directory for intermediate files (which tend to be very large)",
																																																 "");
		final GenvisisWorkflow.Requirement correctionStrategyReq = new GenvisisWorkflow.EnumRequirement("Correction Type",
																																																		CORRECTION_TYPE.XY);
		final GenvisisWorkflow.Requirement sexChromosomeStrategyReq = new GenvisisWorkflow.EnumRequirement("Sex Chromosome Strategy",
																																																			 CHROMOSOME_X_STRATEGY.BIOLOGICAL);

		return new Step("Create 'Shadow' Project", "",
										new GenvisisWorkflow.Requirement[][] {{parseSamplesStepReq},
																													{numPCsReq},
																													{outputBaseReq},
																													{callrateReq},
																													{recomputeLrrReq},
																													{tempDirReq},
																													{correctionStrategyReq},
																													{sexChromosomeStrategyReq},
																													{numThreadsReq},},
										new Step[] {parseSamplesStep},
										new GenvisisWorkflow.FLAG[] {GenvisisWorkflow.FLAG.MEMORY,
																								 GenvisisWorkflow.FLAG.RUNTIME},
										priority) {

			@Override
			public void setNecessaryPreRunProperties(Project proj,
																							 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				// not needed for step
			}

			@Override
			public void run(Project proj,
											Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
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

				int totalThreads = resolveThreads(variables.get(this).get(numThreadsReq));
				String retMsg = PRoCtOR.shadow(proj, tmpDir, outputBase, markerCallRateFilter,
																			 recomputeLRRPCs, type, strategy, numComponents,
																			 totalThreads);
				if (!"".equals(retMsg)) {
					setFailed(retMsg);
				}
			}

			@Override
			public String getCommandLine(Project proj,
																	 Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
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

				int totalThreads = resolveThreads(variables.get(this).get(numThreadsReq));

				String projPropFile = proj.getPropertyFilename();
				StringBuilder cmd = new StringBuilder();
				cmd.append(Files.getRunString()).append(" org.genvisis.cnv.manage.PRoCtOR").append(" proj=")
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
			public boolean checkIfOutputExists(Map<Step, Map<GenvisisWorkflow.Requirement, String>> variables) {
				String outputBase = proj.PROJECT_DIRECTORY.getValue()
														+ variables.get(this).get(outputBaseReq);
				String finalReport = outputBase + PCA.FILE_EXTs[0];
				return Files.exists(finalReport);
			}
		};

	}

	private String getPlinkDir() {
		String dir = "plink/";
		dir = ext.verifyDirFormat(proj.PROJECT_DIRECTORY.getValue() + dir);
		return dir;

	}

	private int resolveThreads(String arg) {
		int numThreads = GenvisisWorkflow.Requirement.checkIntArgOrNeg1(arg);
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
