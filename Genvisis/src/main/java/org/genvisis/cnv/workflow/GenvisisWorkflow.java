package org.genvisis.cnv.workflow;

import java.io.File;
import java.util.Collections;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.StringJoiner;

import org.genvisis.CLI;
import org.genvisis.cnv.Launch;
import org.genvisis.cnv.analysis.Mosaicism;
import org.genvisis.cnv.analysis.pca.PCA;
import org.genvisis.cnv.analysis.pca.PCAPrep;
import org.genvisis.cnv.analysis.pca.PCImputeRace;
import org.genvisis.cnv.analysis.pca.PCImputeRace.RACE;
import org.genvisis.cnv.analysis.pca.PrincipalComponentsCompute.PRE_PROCESSING_METHOD;
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
import org.genvisis.cnv.hmm.CNVCaller.CALLING_SCOPE;
import org.genvisis.cnv.hmm.CNVCaller.PFB_MANAGEMENT_TYPE;
import org.genvisis.cnv.manage.MitoPipeline;
import org.genvisis.cnv.manage.PRoCtOR;
import org.genvisis.cnv.manage.PlinkData;
import org.genvisis.cnv.manage.Resources;
import org.genvisis.cnv.manage.TransposeData;
import org.genvisis.cnv.prop.Property;
import org.genvisis.cnv.qc.AffyMarkerBlast;
import org.genvisis.cnv.qc.GcAdjustor;
import org.genvisis.cnv.qc.IlluminaMarkerBlast;
import org.genvisis.cnv.qc.LrrSd;
import org.genvisis.cnv.qc.MarkerBlast;
import org.genvisis.cnv.qc.MarkerBlast.FILE_SEQUENCE_TYPE;
import org.genvisis.cnv.qc.MarkerBlastQC;
import org.genvisis.cnv.qc.MarkerMetrics;
import org.genvisis.cnv.qc.SampleQC;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.cnv.workflow.RequirementSet.RequirementSetBuilder;
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
import org.genvisis.qsub.Qsub;

import com.google.common.base.Joiner;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

public class GenvisisWorkflow {

	private static final String NUM_THREADS_DESC = "Number of Threads to Use";
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
		numThreadsReq = new Requirement.PosIntRequirement(NUM_THREADS_DESC, proj.NUM_THREADS.getValue());

		steps = Collections.unmodifiableSortedSet(generateSteps(!project.IS_PC_CORRECTED_PROJECT.getValue()));
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

		private Step generateIlluminaMarkerPositionsStep() {
			final Requirement snpMapReq = new Requirement.FileRequirement(
																																		"An Illumina SNP_map file.",
																																		proj.getLocationOfSNP_Map(false));
			final Requirement manifestReq = new Requirement.FileRequirement(
																																			"An Illumina Manifest file.",
																																			proj.getLocationOfSNP_Map(false));

			final RequirementSet reqSet = RequirementSetBuilder.and()
																												 .add(RequirementSetBuilder.or()
																																									 .add(snpMapReq)
																																									 .add(manifestReq));
			return register(new Step("Create Marker Positions (if not already exists)", "", reqSet,
															 EnumSet.noneOf(Requirement.Flag.class),
															 priority()) {

				@Override
				public void setNecessaryPreRunProperties(Project proj,
																								 Map<Step, Map<Requirement, String>> variables) {
					// not needed for step
				}

				@Override
				public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
					proj.getLog().report("Generating marker positions file");
					String snpMap = variables.get(this).get(snpMapReq);
					String manifest = variables.get(this).get(manifestReq);
					if (Files.exists(snpMap)) {
						org.genvisis.cnv.manage.Markers.generateMarkerPositions(proj, snpMap);
					} else if (Files.exists(manifest)) {
						MarkerBlast.extractMarkerPositionsFromManifest(manifest,
																													 ARRAY.ILLUMINA,
																													 FILE_SEQUENCE_TYPE.MANIFEST_FILE,
																													 proj.MARKER_POSITION_FILENAME.getValue(false,
																																																	false),
																													 Files.determineDelimiter(manifest, log),
																													 log);
					}
				}

				@Override
				public boolean checkIfOutputExists(Map<Step, Map<Requirement, String>> variables) {
					return Files.exists(proj.MARKER_POSITION_FILENAME.getValue(false, false));
				}

				@Override
				public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
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

			});
		}

		private Step generateIlluminaMarkerBlastAnnotationStep(final Step parseSamplesStep) {
			final Requirement parseSamplesStepReq = new Requirement.StepRequirement(parseSamplesStep);
			final Requirement manifestFileReq = new Requirement.FileRequirement(
																																					ext.capitalizeFirst(IlluminaMarkerBlast.DESC_MANIFEST),
																																					IlluminaMarkerBlast.EXAMPLE_MANIFEST);

			final RequirementSet reqSet = RequirementSetBuilder.and()
																												 .add(parseSamplesStepReq)
																												 .add(manifestFileReq)
																												 .add(getNumThreadsReq());

			return register(new Step("Run Marker BLAST Annotation", "", reqSet,
															 EnumSet.of(Requirement.Flag.MEMORY, Requirement.Flag.RUNTIME,
																					Requirement.Flag.MULTITHREADED),
															 priority()) {

				@Override
				public void setNecessaryPreRunProperties(Project proj,
																								 Map<Step, Map<Requirement, String>> variables) {
					// Not necessary for this step

				}

				@Override
				public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
					String manifestFile = variables.get(this).get(manifestFileReq);
					int numThreads = resolveThreads(variables.get(this).get(getNumThreadsReq()));
					new IlluminaMarkerBlast(proj, numThreads, manifestFile).blastEm();
				}

				@Override
				public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
					String manifestFile = variables.get(this).get(manifestFileReq);
					int numThreads = resolveThreads(variables.get(this).get(getNumThreadsReq()));
					ImmutableMap.Builder<String, String> argsBuilder = ImmutableMap.builder();
					argsBuilder.put(CLI.ARG_PROJ, proj.getPropertyFilename());
					argsBuilder.put(IlluminaMarkerBlast.ARG_MANIFEST, manifestFile);
					argsBuilder.put(CLI.ARG_THREADS, String.valueOf(numThreads));
					return Files.getRunString() + " "
								 + CLI.formCmdLine(IlluminaMarkerBlast.class, argsBuilder.build());
				}

				@Override
				public boolean checkIfOutputExists(Map<Step, Map<Requirement, String>> variables) {
					return Files.exists(proj.BLAST_ANNOTATION_FILENAME.getValue());
				}
			});
		}

		private Step generateAffyMarkerBlastAnnotationStep(final Step parseSamplesStep) {
			final Requirement parseSamplesStepReq = new Requirement.StepRequirement(parseSamplesStep);

			final Requirement probeFileReq = new Requirement.FileRequirement(
																																			 ext.capitalizeFirst(AffyMarkerBlast.DESC_PROBE_FILE),
																																			 AffyMarkerBlast.EXAMPLE_PROBE_FILE);
			final Requirement annotFileReq = new Requirement.FileRequirement(
																																			 ext.capitalizeFirst(AffyMarkerBlast.DESC_ANNOT_FILE),
																																			 AffyMarkerBlast.EXAMPLE_ANNOT_FILE);

			final RequirementSet reqSet = RequirementSetBuilder.and().add(parseSamplesStepReq)
																												 .add(probeFileReq).add(annotFileReq)
																												 .add(getNumThreadsReq());

			return register(new Step("Run Marker BLAST Annotation", "", reqSet,
															 EnumSet.of(Requirement.Flag.MEMORY, Requirement.Flag.RUNTIME,
																					Requirement.Flag.MULTITHREADED),
															 priority()) {

				@Override
				public void setNecessaryPreRunProperties(Project proj,
																								 Map<Step, Map<Requirement, String>> variables) {
					// Not necessary for this step

				}

				@Override
				public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
					String annotFile = variables.get(this).get(annotFileReq);
					String probeFile = variables.get(this).get(probeFileReq);
					int numThreads = resolveThreads(variables.get(this).get(getNumThreadsReq()));
					new AffyMarkerBlast(proj, numThreads, probeFile, annotFile).blastEm();
				}

				@Override
				public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
					String annotFile = variables.get(this).get(annotFileReq);
					String probeFile = variables.get(this).get(probeFileReq);
					int numThreads = resolveThreads(variables.get(this).get(getNumThreadsReq()));
					ImmutableMap.Builder<String, String> argsBuilder = ImmutableMap.builder();
					argsBuilder.put(CLI.ARG_PROJ, proj.getPropertyFilename());
					argsBuilder.put(AffyMarkerBlast.ARG_PROBE_FILE, probeFile);
					argsBuilder.put(AffyMarkerBlast.ARG_ANNOT_FILE, annotFile);
					argsBuilder.put(CLI.ARG_THREADS, String.valueOf(numThreads));
					return Files.getRunString() + " "
								 + CLI.formCmdLine(AffyMarkerBlast.class, argsBuilder.build());
				}

				@Override
				public boolean checkIfOutputExists(Map<Step, Map<Requirement, String>> variables) {
					return Files.exists(proj.BLAST_ANNOTATION_FILENAME.getValue());
				}
			});
		}

		private Step generateParseSamplesStep() {
			return generateParseSamplesStep(null);
		}

		private Step generateParseSamplesStep(final Step markerPositionsStep) {

			final Requirement markerPositionsReq = new Requirement.FileRequirement(
																																						 "Marker Positions file must already exist.",
																																						 proj.MARKER_POSITION_FILENAME.getValue(false,
																																																										false));

			final RequirementSet reqSet = RequirementSetBuilder.and();
			if (markerPositionsStep == null) {
				reqSet.add(markerPositionsReq).add(getNumThreadsReq());
			} else {
				final Requirement markerPositionsStepReq = new Requirement.StepRequirement(
																																									 markerPositionsStep);
				reqSet.add(RequirementSetBuilder.or().add(markerPositionsReq).add(markerPositionsStepReq))
							.add(getNumThreadsReq());
			}

			return register(new Step("Parse Sample Files", "", reqSet,
															 EnumSet.of(Requirement.Flag.MEMORY, Requirement.Flag.RUNTIME,
																					Requirement.Flag.MULTITHREADED), priority()) {

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
							throw new RuntimeException(
																				 "Operation failure, please check log for more information.");
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
					returnValue = returnValue && proj.getSampleList() != null;
					returnValue = returnValue && Files.exists(sampleDirectory);

					int numSamples = proj.getSampleList().getSamples().length;
					returnValue = returnValue && numSamples > 0;
					returnValue = returnValue
												&& Files.countFiles(sampleDirectory, Sample.SAMPLE_FILE_EXTENSION) == numSamples;
					// checking the validity / completeness of each sample would be a Good Thing, but too
					// costly time-wise for larger projects
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

		private Step generateCreateSampleDataStep(final Step parseSamplesStep) {
			final Requirement parseSamplesStepReq = new Requirement.StepRequirement(parseSamplesStep);

			final Requirement createMinimalSampleDataReq = new Requirement.BoolRequirement(
																																										 "Create a minimal SampleData.txt file from sample files",
																																										 true);

			final String pedPreset = proj.PEDIGREE_FILENAME.getValue();

			final Requirement pedigreeReq = new Requirement(
																											"Either a Pedigree.dat file, or any file with a header containing all of the following elements (in any order):  \""
																													+ ArrayUtils.toStr(MitoPipeline.PED_INPUT,
																																						 ", ")
																													+ "\"",
																											Requirement.RequirementInputType.FILE,
																											pedPreset) {

				@Override
				public boolean checkRequirement(String arg, Set<Step> stepSelections,
																				Map<Step, Map<Requirement, String>> variables) {
					return Files.exists(arg);
				}
			};

			// check for SampleMap only if we haven't found a pedigree
			final String sampMapPreset = Files.exists(pedPreset) ? null : getLocationOfSampleMap();

			final Requirement sampMapReq = new Requirement(
																										 "A Sample_Map.csv file, with at least two columns having headers \""
																												 + MitoPipeline.SAMPLEMAP_INPUT[1]
																												 + "\" and \""
																												 + MitoPipeline.SAMPLEMAP_INPUT[2] + "\"",
																										 Requirement.RequirementInputType.FILE,
																										 sampMapPreset) {

				@Override
				public boolean checkRequirement(String arg, Set<Step> stepSelections,
																				Map<Step, Map<Requirement, String>> variables) {
					return Files.exists(arg);
				}
			};
			final RequirementSet reqSet = RequirementSetBuilder.and()
																												 .add(parseSamplesStepReq)
																												 .add(RequirementSetBuilder.or()
																																									 .add(createMinimalSampleDataReq)
																																									 .add(pedigreeReq)
																																									 .add(sampMapReq));
			return register(new Step("Create SampleData.txt File", "", reqSet,
															 EnumSet.noneOf(Requirement.Flag.class), priority()) {

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
							throw new RuntimeException("SampleData already exists - please delete and try again.");
						}
					} catch (Elision e) {
						throw new RuntimeException(e.getMessage());
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
			final RequirementSet reqSet = RequirementSetBuilder.and()
																												 .add(new Requirement.StepRequirement(
																																															parseSamplesStep));
			return register(new Step("Transpose Data into Marker-Dominant Files", "", reqSet,
															 EnumSet.of(Requirement.Flag.MEMORY), priority()) {

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
					return Files.countFiles(proj.MARKER_DATA_DIRECTORY.getValue(false, false),
																	MarkerData.MARKER_DATA_FILE_EXTENSION) > 0;
				}
			});
		}

		private Step generateGCModelStep() {
			final Requirement.ResourceRequirement gcBaseResourceReq = new Requirement.ResourceRequirement(
																																																		"GC Base file",
																																																		Resources.genome(proj.GENOME_BUILD_VERSION.getValue(),
																																																										 proj.getLog())
																																																						 .getModelBase());
			final Requirement gcModelOutputReq = new Requirement.OutputFileRequirement(
																																								 "GCModel output file must be specified.",
																																								 proj.GC_MODEL_FILENAME.getValue());

			final RequirementSet reqSet = RequirementSetBuilder.and().add(gcBaseResourceReq)
																												 .add(gcModelOutputReq);
			return register(new Step("Compute GCMODEL File", "", reqSet,
															 EnumSet.noneOf(Requirement.Flag.class), priority()) {

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
					String gcBaseFile = gcBaseResourceReq.getResource().getAbsolute();
					String gcOutputFile = variables.get(this).get(gcModelOutputReq);
					org.genvisis.cnv.analysis.PennCNV.gcModel(proj, gcBaseFile, gcOutputFile, 100);
				}

				@Override
				public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
					String kvCmd = "";

					String setGCOutputFile = proj.GC_MODEL_FILENAME.getValue();
					String gcOutputFile = variables == null ? null
																								 : variables.get(this).get(gcModelOutputReq);
					if (gcOutputFile != null && !ext.verifyDirFormat(setGCOutputFile).equals(gcOutputFile)) {
						kvCmd += " GC_MODEL_FILENAME=" + gcOutputFile;
					}

					String projPropFile = proj.getPropertyFilename();
					StringBuilder cmd = new StringBuilder();
					if (kvCmd.length() > 0) {
						cmd.append(Files.getRunString()).append(PROJ_PROP_UPDATE_STR + projPropFile)
							 .append(kvCmd)
							 .append("\n");
					}
					String gcBaseFile = gcBaseResourceReq.getResource().getAbsolute();
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
			final Requirement parseSamplesStepReq = new Requirement.StepRequirement(parseSamplesStep);
			final RequirementSet reqSet = RequirementSetBuilder.and().add(parseSamplesStepReq)
																												 .add(getNumThreadsReq());
			return register(new Step("Run Sample QC Metrics", "", reqSet,
															 EnumSet.of(Requirement.Flag.MULTITHREADED), priority()) {

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
			final Requirement parseSamplesStepReq = new Requirement.StepRequirement(parseSamplesStep);

			final Requirement exportAllReq = new Requirement.OptionalBoolRequirement(
																																							 "Export all markers in project.",
																																							 true);

			String[] tgtMkrFiles = proj.TARGET_MARKERS_FILENAMES.getValue();
			final Requirement targetMarkersReq = new Requirement.FileRequirement(
																																					 "A targetMarkers files listing the markers to QC.",
																																					 tgtMkrFiles != null
																																							 && tgtMkrFiles.length >= 1
																																																				 ? tgtMkrFiles[0]
																																																				 : "");
			final Set<String> sampleDataHeaders;
			if (Files.exists(proj.SAMPLE_DATA_FILENAME.getValue()) && proj.getSampleData(false) != null) {
				sampleDataHeaders = proj.getSampleData(false).getMetaHeaders();
			} else {
				sampleDataHeaders = Sets.newHashSet();
			}
			final Set<String> defaultBatchHeaders = Sets.intersection(sampleDataHeaders,
																																MarkerMetrics.DEFAULT_SAMPLE_DATA_BATCH_HEADERS);
			final Requirement.ListSelectionRequirement batchHeadersReq = new Requirement.ListSelectionRequirement(
																																																						"SampleData column headers to use as batches for batch effects calculations",
																																																						sampleDataHeaders,
																																																						defaultBatchHeaders,
																																																						true);

			final RequirementSet reqSet = RequirementSetBuilder.and()
																												 .add(parseSamplesStepReq)
																												 .add(RequirementSetBuilder.or()
																																									 .add(exportAllReq)
																																									 .add(targetMarkersReq))
																												 .add(batchHeadersReq)
																												 .add(getNumThreadsReq());

			return register(new Step("Run Marker QC Metrics", "", reqSet,
															 EnumSet.of(Requirement.Flag.MULTITHREADED), priority()) {

				@Override
				public void setNecessaryPreRunProperties(Project proj,
																								 Map<Step, Map<Requirement, String>> variables) {
					// Nothing to do here
				}

				@Override
				public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
					boolean allMarkers = Boolean.parseBoolean(variables.get(this).get(exportAllReq));
					String tgtFile = allMarkers ? null : variables.get(this).get(targetMarkersReq);
					boolean[] samplesToExclude = proj.getSamplesToExclude();
					int numThreads = resolveThreads(variables.get(this).get(getNumThreadsReq()));
					Set<String> batchHeaders = ImmutableSet.copyOf(Requirement.ListSelectionRequirement.parseArgValString(variables.get(this)
																																																												 .get(batchHeadersReq)));
					MarkerMetrics.fullQC(proj, samplesToExclude, tgtFile, true, batchHeaders, numThreads);
				}

				@Override
				public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
					boolean allMarkers = Boolean.parseBoolean(variables.get(this).get(exportAllReq));
					String tgtFile = variables.get(this).get(targetMarkersReq);
					int numThreads = resolveThreads(variables.get(this).get(getNumThreadsReq()));
					List<String> batchHeaders = Requirement.ListSelectionRequirement.parseArgValString(variables.get(this)
																																																			.get(batchHeadersReq));
					StringJoiner args = new StringJoiner(" ");
					args.add(Files.getRunString());
					args.add(MarkerMetrics.class.getCanonicalName());
					args.add("-fullQC");
					args.add("proj=" + proj.getPropertyFilename());
					if (!allMarkers)
						args.add("markers=" + tgtFile);
					String batchHeadersArg = String.join(",", batchHeaders);
					args.add(MarkerMetrics.BATCH_HEADERS_ARG + "=" + batchHeadersArg);
					args.add(PSF.Ext.NUM_THREADS_COMMAND + numThreads);
					return args.toString();
				}

				@Override
				public boolean checkIfOutputExists(Map<Step, Map<Requirement, String>> variables) {
					String markerMetricsFile = proj.MARKER_METRICS_FILENAME.getValue();
					return Files.exists(markerMetricsFile);
				}

			});
		}

		private Step generateSexChecksStep(final Step parseSamplesStep, final Step markerBlastStep,
																			 final Step sampleDataStep, final Step transposeStep,
																			 final Step sampleQCStep) {
			final Requirement parseSamplesStepReq = new Requirement.StepRequirement(parseSamplesStep);
			final Requirement sampleDataStepReq = new Requirement.StepRequirement(sampleDataStep);
			final Requirement transposeStepReq = new Requirement.StepRequirement(transposeStep);
			final Requirement sampleQCStepReq = new Requirement.StepRequirement(sampleQCStep);
			final Requirement addToSampleDataReq = new Requirement.OptionalBoolRequirement(
																																										 "Add Estimated Sex to Sample Data",
																																										 true);

			final Requirement oneHittersReq = new Requirement.FileRequirement(
																																				"List of markers that do not cross hybridize",
																																				MarkerBlastQC.defaultOneHitWondersFilename(proj.BLAST_ANNOTATION_FILENAME.getValue()));
			final Requirement markerBlastStepReq = new Requirement.StepRequirement(markerBlastStep);
			final Requirement noCrossHybeReq = new Requirement.BoolRequirement(
																																				 "Use only X and Y chromosome R values to identify sex discriminating markers",
																																				 false);

			final RequirementSet reqSet = RequirementSetBuilder.and()
																												 .add(parseSamplesStepReq)
																												 .add(sampleDataStepReq)
																												 .add(transposeStepReq)
																												 .add(sampleQCStepReq)
																												 .add(RequirementSetBuilder.or()
																																									 .add(oneHittersReq)
																																									 .add(markerBlastStepReq)
																																									 .add(noCrossHybeReq));

			return register(new Step("Run Sex Checks", "", reqSet,
															 EnumSet.noneOf(Requirement.Flag.class), priority()) {

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
							MarkerBlastQC.getOneHitWonders(proj, proj.BLAST_ANNOTATION_FILENAME.getValue(),
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
												 + proj.BLAST_ANNOTATION_FILENAME.getValue())
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
			final Requirement parseSamplesStepReq = new Requirement.StepRequirement(parseSamplesStep);
			final Requirement pedigreeRequirement = new Requirement.FileRequirement(
																																							"A pedigree.dat file must exist.",
																																							proj.PEDIGREE_FILENAME.getValue(false,
																																																							false));
			final Requirement createPedigreeRequirement = new Requirement.BoolRequirement(
																																										"Create a minimal pedigree.dat file [will pull information from SexChecks step results].",
																																										false);

			final RequirementSet reqSet = RequirementSetBuilder.and()
																												 .add(parseSamplesStepReq)
																												 .add(RequirementSetBuilder.or()
																																									 .add(pedigreeRequirement)
																																									 .add(createPedigreeRequirement));


			return register(new Step("Create PLINK Files", "", reqSet,
															 EnumSet.of(Requirement.Flag.MEMORY), priority()) {

				@Override
				public Map<Requirement, String> getDefaultRequirementValues() {
					Map<Requirement, String> varMap = super.getDefaultRequirementValues();
					if (!Files.exists(proj.PEDIGREE_FILENAME.getValue())) {
						// if no pedigree, default to creating a minimal one
						varMap.put(createPedigreeRequirement, Boolean.TRUE.toString());
					}
					return varMap;
				}

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
						throw new RuntimeException(
																			 "Creation of Pedigree file in [Create/Run PLINK Files] step failed.");
					}

					proj.getLog().report("Running PLINK");

					boolean create = PlinkData.saveGenvisisToPlinkBedSet(proj, PLINK_SUBDIR + PLINKROOT,
																															 null, null);
					if (!create) {
						throw new RuntimeException("Creation of initial PLINK files failed.");
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
										 + " proj=")
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
			final Requirement plinkExportStepReq = new Requirement.StepRequirement(plinkExportStep);
			String defaultCallrate;
			switch (proj.getArrayType()) {
				case AFFY_GW6:
				case AFFY_GW6_CN:
					defaultCallrate = MarkerQC.DEFAULT_AFFY_CALLRATE_THRESHOLD;
					break;
				case AFFY_AXIOM:
				case ILLUMINA:
					defaultCallrate = MarkerQC.DEFAULT_ILLUMINA_CALLRATE_THRESHOLD;
					break;
				default:
					throw new IllegalArgumentException("Invalid " + proj.getArrayType().getClass().getName()
																						 + ": " + proj.getArrayType().toString());
			}
			final Requirement callrateReq = new Requirement.ThresholdRequirement(
																																					 QC_METRIC.CALLRATE.getUserDescription(),
																																					 defaultCallrate);
			final RequirementSet reqSet = RequirementSetBuilder.and().add(plinkExportStepReq)
																												 .add(callrateReq);

			return register(new Step("Run GWAS QC", "", reqSet,
															 EnumSet.noneOf(Requirement.Flag.class), priority()) {

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
					if (new File(dir + Qc.QC_SUBDIR + RelationAncestryQc.GENOME_DIR + PLINKROOT
											 + ".genome").exists()) {
						proj.GENOME_CLUSTER_FILENAME.setValue(dir + Qc.QC_SUBDIR
																									+ RelationAncestryQc.GENOME_DIR
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
			final Requirement gwasQCStepReq = new Requirement.StepRequirement(gwasQCStep);
			final Requirement putativeWhitesReq = new Requirement.FileRequirement(
																																						"File with FID/IID pairs of putative white samples",
																																						"");
			final Requirement.ResourceRequirement hapMapFoundersReq = new Requirement.ResourceRequirement(
																																																		"PLINK root of HapMap founders",
																																																		Resources.hapMap(log)
																																																						 .getUnambiguousHapMapFounders());

			final RequirementSet reqSet = RequirementSetBuilder.and()
																												 .add(gwasQCStepReq).add(putativeWhitesReq)
																												 .add(hapMapFoundersReq);

			return register(new Step("Run Ancestry Checks", "", reqSet,
															 EnumSet.noneOf(Requirement.Flag.class), priority()) {

				@Override
				public void setNecessaryPreRunProperties(Project proj,
																								 Map<Step, Map<Requirement, String>> variables) {
					// not needed for step
				}

				@Override
				public void run(Project proj, Map<Step, Map<Requirement, String>> variables) {
					String putativeWhites = variables.get(this).get(putativeWhitesReq);
					String hapMapPlinkRoot = hapMapFoundersReq.getResource().getAbsolute();
					String ancestryDir = getAncestryDir();
					Ancestry.runPipeline(ancestryDir, putativeWhites, hapMapPlinkRoot, proj,
															 new Logger(ancestryDir + "ancestry.log"));
				}

				@Override
				public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
					String putativeWhites = variables.get(this).get(putativeWhitesReq);
					String hapMapPlinkRoot = hapMapFoundersReq.getResource().getAbsolute();

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
			final Requirement plinkExportStepReq = new Requirement.StepRequirement(plinkExportStep);
			final Requirement gwasQCStepReq = new Requirement.StepRequirement(gwasQCStep);
			final Requirement ancestryStepReq = new Requirement.StepRequirement(ancestryStep);
			final Requirement unrelatedsFileReq = new Requirement.FileRequirement(
																																						"File with list of unrelated FID/IID pairs to use for marker QC",
																																						"");
			final Requirement europeansFilesReq = new Requirement.FileRequirement(
																																						"File with list of European samples to use for Hardy-Weinberg equilibrium tests",
																																						"");

			final RequirementSet reqSet = RequirementSetBuilder.and()
																												 .add(plinkExportStepReq)
																												 .add(RequirementSetBuilder.or()
																																									 .add(gwasQCStepReq)
																																									 .add(unrelatedsFileReq))
																												 .add(RequirementSetBuilder.or()
																																									 .add(ancestryStepReq)
																																									 .add(europeansFilesReq));
			final Map<QC_METRIC, Requirement> metricRequirements = Maps.newEnumMap(QC_METRIC.class);
			for (QC_METRIC metric : QC_METRIC.values()) {
				Map<QC_METRIC, String> defaultThresholds = FurtherAnalysisQc.getDefaultMarkerQCThresholds(proj.getArrayType());
				String defaultVal = defaultThresholds.get(metric);
				final Requirement metricReq = new Requirement.ThresholdRequirement(
																																					 metric.getUserDescription(),
																																					 defaultVal);
				reqSet.add(metricReq);
				metricRequirements.put(metric, metricReq);
			}

			return register(new Step("Run Further Analysis QC", "", reqSet,
															 EnumSet.noneOf(Requirement.Flag.class), priority()) {

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
					return PSF.Plink.bedBimFamExist(dir + qcdPlinkroot)
								 && Files.exists(dir + FurtherAnalysisQc.SAMPLE_QC_DROPS, false)
								 && Files.exists(dir + FurtherAnalysisQc.MARKER_QC_DROPS, false);
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
			final Requirement parseSamplesStepReq = new Requirement.StepRequirement(parseSamplesStep);
			final RequirementSet reqSet = RequirementSetBuilder.and().add(parseSamplesStepReq)
																												 .add(getNumThreadsReq());
			return register(new Step("Create Mosaic Arms File", "", reqSet,
															 EnumSet.of(Requirement.Flag.MULTITHREADED), priority()) {

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
			final Requirement sampleQCStepReq = new Requirement.StepRequirement(sampleQCStep);
			final Requirement createSampleDataStepReq = new Requirement.StepRequirement(
																																									createSampleDataStep);
			final Requirement skipIDingDuplicatesReq = new Requirement.BoolRequirement(
																																								 "Skip identifying duplicates",
																																								 false);
			final Requirement gwasQCStepReq = new Requirement.StepRequirement(gwasQCStep);
			final Requirement notGcCorrectedLrrSdReq = new Requirement.BoolRequirement(
																																								 "Do not use GC corrected LRR SD?",
																																								 false);
			final Requirement gcCorrectedLrrSdReq = new Requirement(
																															"GC Corrected LRR SD must exist in Sample QC File",
																															Requirement.RequirementInputType.NONE) {

				@Override
				public boolean checkRequirement(String arg, Set<Step> stepSelections,
																				Map<Step, Map<Requirement, String>> variables) {
					String sampleQCFile = proj.SAMPLE_QC_FILENAME.getValue();
					return Files.exists(sampleQCFile)
								 && ext.indexOfStr("LRR_SD_Post_Correction",
																	 Files.getHeaderOfFile(sampleQCFile, proj.getLog())) != -1;
				}

			};
			final Requirement lrrSdThresholdReq = new Requirement.DoubleRequirement(
																																							"LRR SD Threshold",
																																							proj.LRRSD_CUTOFF.getValue(),
																																							proj.LRRSD_CUTOFF.getMinValue(),
																																							proj.LRRSD_CUTOFF.getMaxValue());

			final Requirement callrateThresholdReq = new Requirement.DoubleRequirement(
																																								 "Callrate Threshold",
																																								 proj.SAMPLE_CALLRATE_THRESHOLD.getValue(),
																																								 proj.SAMPLE_CALLRATE_THRESHOLD.getMinValue(),
																																								 proj.SAMPLE_CALLRATE_THRESHOLD.getMaxValue());
			final Requirement numQReq = new Requirement.PosIntRequirement(
																																		"Number of Quantiles to Generate",
																																		10);
			final Requirement replaceFIDIIDReq = new Requirement.OptionalBoolRequirement(
																																									 "Replace FID and IID with data from Pedigree",
																																									 false);

			final RequirementSet reqSet = RequirementSetBuilder.and()
																												 .add(sampleQCStepReq)
																												 .add(createSampleDataStepReq)
																												 .add(RequirementSetBuilder.or()
																																									 .add(skipIDingDuplicatesReq)
																																									 .add(gwasQCStepReq))
																												 .add(RequirementSetBuilder.or()
																																									 .add(notGcCorrectedLrrSdReq)
																																									 .add(gcCorrectedLrrSdReq))
																												 .add(lrrSdThresholdReq)
																												 .add(callrateThresholdReq)
																												 .add(numQReq)
																												 .add(replaceFIDIIDReq);

			return register(new Step("Annotate Sample Data File", "", reqSet,
															 EnumSet.noneOf(Requirement.Flag.class), priority()) {

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
					if (checkDuplicates
							&& ext.indexOfStr(SampleQC.DUPLICATE_ID_HEADER, header, false, true) == -1) {
						return false;
					}
					String[] reqHdr = {SampleQC.EXCLUDE_HEADER, "ExcludeNote", "Use", "UseNote", "Use_cnv",
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

		private Step generateMitoCNEstimateStep(Step transposeStep) {
			// FIXME http://genvisis.org/MitoPipeline/#illumina_marker_lists has illumina markers.. this
			// should be linked to, or
			// these steps split or something...
			final Requirement transposeStepReq = new Requirement.StepRequirement(transposeStep);
			final Requirement medianMarkersReq = new Requirement.FileRequirement(
																																					 "MedianMarkers file must exist.",
																																					 "");
			final Requirement lrrSdThresholdReq = new Requirement.DoubleRequirement(
																																							"LRR SD threshold to filter samples.",
																																							proj.LRRSD_CUTOFF.getValue(),
																																							proj.LRRSD_CUTOFF.getMinValue(),
																																							proj.LRRSD_CUTOFF.getMaxValue());
			final Requirement callrateThresholdReq = new Requirement.DoubleRequirement(
																																								 "Call rate threshold to filter markers.",
																																								 MitoPipeline.DEFAULT_MKR_CALLRATE_FILTER,
																																								 0.0, 1.0);
			final Requirement qcPassingOnlyReq = new Requirement.OptionalBoolRequirement(
																																									 "Compute PCs with samples passing QC only",
																																									 true);
			final Requirement imputeNaNs = new Requirement.OptionalBoolRequirement(
																																						 "Impute mean value for NaN",
																																						 true);
			final Requirement recomputeLrrPCMarkersReq = new Requirement.OptionalBoolRequirement(
																																													 "Should recompute Log-R ratio for PC markers?",
																																													 true);
			final Requirement recomputeLrrMedianMarkersReq = new Requirement.OptionalBoolRequirement(
																																															 "Should recompute Log-R ratio for median markers?",
																																															 true);
			final Requirement homozygousOnlyReq = new Requirement.OptionalBoolRequirement(
																																										"Homozygous only?",
																																										true);
			final Requirement gcRegressionDistanceReq = new Requirement.PosIntRequirement(
																																										"Regression distance for the GC adjustment",
																																										GcAdjustor.DEFAULT_REGRESSION_DISTANCE[0]);
			final Requirement pcSelectionSamplesReq = new Requirement.OptionalFileRequirement(
																																												"A file listing a subset of samples (DNA ID) to use for determining optimal PC selection, typically a list of unrelated and single race samples. If a list is not provided, only samples passing sample qc thresholds will be used.",
																																												"");
			final Requirement externalBetaFileReq = new Requirement.OptionalFileRequirement(
																																											"An external beta file to optimize PC selection.",
																																											"");

			final RequirementSet reqSet = RequirementSetBuilder.and()
																												 .add(transposeStepReq)
																												 .add(medianMarkersReq)
																												 .add(lrrSdThresholdReq)
																												 .add(callrateThresholdReq)
																												 .add(qcPassingOnlyReq)
																												 .add(imputeNaNs)
																												 .add(recomputeLrrPCMarkersReq)
																												 .add(recomputeLrrMedianMarkersReq)
																												 .add(homozygousOnlyReq)
																												 .add(gcRegressionDistanceReq)
																												 .add(getNumThreadsReq())
																												 .add(pcSelectionSamplesReq)
																												 .add(externalBetaFileReq);

			return register(new Step("Create Mitochondrial Copy-Number Estimates File",
															 "", reqSet, EnumSet.of(Requirement.Flag.MULTITHREADED), priority()) {

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
						MitoPipeline.estimateMtDNACN(proj, numThreads, medianMarkers, numComponents,
																				 outputBase,
																				 homozygousOnly, markerCallRateFilter, betaOptFile,
																				 pedFile,
																				 recomputeLRRPCs, recomputeLRRMedian, sampLrr,
																				 imputeMeanForNaN, gcCorrect, bpGcModel,
																				 regressionDistance,
																				 proj.GENOME_BUILD_VERSION.getValue(), pvalOpt, betaFile,
																				 plot, false, PRE_PROCESSING_METHOD.NONE, proj.getLog());
					} else {
						throw new RuntimeException(PCAPrep.errorMessage(retCode));
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
			final Requirement parseSamplesStepReq = new Requirement.StepRequirement(parseSamplesStep);
			final Requirement sampleSubsetReq = new Requirement.FileRequirement(
																																					"A Sample subset file must exist.",
																																					proj.SAMPLE_SUBSET_FILENAME.getValue());
			String defaultOutputFile;
			if (Files.exists(proj.SAMPLE_SUBSET_FILENAME.getValue())) {
				defaultOutputFile = ext.rootOf(proj.SAMPLE_SUBSET_FILENAME.getValue()) + ".pfb";
			} else {
				defaultOutputFile = proj.CUSTOM_PFB_FILENAME.getValue();
			}
			final Requirement outputFileReq = new Requirement.OutputFileRequirement(
																																							"PFB (population BAF) output file must be specified.",
																																							defaultOutputFile);

			final RequirementSet reqSet = RequirementSetBuilder.and()
																												 .add(RequirementSetBuilder.or()
																																									 .add(parseSamplesStepReq)
																																									 .add(sampleSubsetReq))
																												 .add(outputFileReq);

			return register(new Step("Compute Population BAF files", "", reqSet,
															 EnumSet.noneOf(Requirement.Flag.class), priority()) {

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
					String subSampFile = variables == null ? null : variables.get(this).get(sampleSubsetReq);
					String setPFBFile = proj.CUSTOM_PFB_FILENAME.getValue();
					String pfbOutputFile = variables == null ? null : variables.get(this).get(outputFileReq);

					if (subSampFile != null && !ext.verifyDirFormat(setSubSampFile).equals(subSampFile)) {
						kvCmd += " SAMPLE_SUBSET_FILENAME=" + subSampFile;
					}
					if (pfbOutputFile != null && !ext.verifyDirFormat(setPFBFile).equals(pfbOutputFile)) {
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

		private Step generateSexCentroidsStep() {
			final RequirementSet reqSet = RequirementSetBuilder.and().add(getNumThreadsReq());
			return register(new Step(
															 "Create Sex-Specific Centroids; Filter PFB file",
															 "",
															 reqSet,
															 EnumSet.of(Requirement.Flag.RUNTIME, Requirement.Flag.MULTITHREADED),
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
					String outputDir = proj.DATA_DIRECTORY.getValue();
					malePFB = outputDir + "males.pfb";
					femalePFB = outputDir + "females.pfb";
					centFilePathM = outputDir + "sexSpecific_Male.cent";
					centFilePathF = outputDir + "sexSpecific_Female.cent";

					int numThreads = resolveThreads(variables.get(this).get(getNumThreadsReq()));
					Centroids.computeSexSpecificCentroids(proj, new String[] {malePFB, femalePFB},
																								new String[] {centFilePathM, centFilePathF},
																								numThreads);

				}

				@Override
				public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {

					int numThreads = resolveThreads(variables == null ? "-1"
																													 : variables.get(this)
																																			.get(getNumThreadsReq()));
					String mainCmd = Files.getRunString() + " cnv.filesys.Centroids proj="
													 + proj.getPropertyFilename() + " -sexSpecific "
													 + PSF.Ext.NUM_THREADS_COMMAND + numThreads;
					return mainCmd;
				}

				@Override
				public boolean checkIfOutputExists(Map<Step, Map<Requirement, String>> variables) {
					String malePFB;
					String femalePFB;
					String centFilePathM;
					String centFilePathF;
					String outputDir = proj.DATA_DIRECTORY.getValue();
					malePFB = outputDir + "males.pfb";
					femalePFB = outputDir + "females.pfb";
					centFilePathM = outputDir + "sexSpecific_Male.cent";
					centFilePathF = outputDir + "sexSpecific_Female.cent";
					boolean exists = Files.exists(malePFB);
					exists = exists && Files.exists(femalePFB);
					exists = exists && Files.exists(centFilePathM);
					exists = exists && Files.exists(centFilePathF);
					return exists;
				}

			});
		}

		private Step generateCNVStep(Step pfbStep, Step gcModelStep) {
			final Requirement hmmFile = new Requirement.FileRequirement(
																																	"Hidden Markov Model File Must Exist",
																																	proj.HMM_FILENAME.getValue());
			final Requirement pfbStepReq = new Requirement.StepRequirement(pfbStep);
			final Requirement pfbFileReq = new Requirement.FileRequirement(
																																		 "PFB File Must Exist",
																																		 proj.CUSTOM_PFB_FILENAME.getValue());
			final Requirement gcModelStepReq = new Requirement.StepRequirement(gcModelStep);
			final Requirement gcModelFileReq = new Requirement.FileRequirement(
																																				 "GCMODEL File Must Exist",
																																				 proj.GC_MODEL_FILENAME.getValue());
			final Requirement callingTypeReq = new Requirement.EnumRequirement(
																																				 CNVCaller.CNV_SCOPE_DESC,
																																				 CNVCaller.CALLING_SCOPE.AUTOSOMAL);
			final Requirement useCentroidsReq = new Requirement.OptionalBoolRequirement(
																																									"If calling chromosomal CNVs, use sex-specific centroids to recalculate LRR/BAF values?",
																																									true);
			final Requirement outputFileReq = new Requirement.OutputFileRequirement("Output filename.",
																																							"cnvs/genvisis.cnv") {
				@Override
				public boolean checkRequirement(String arg, Set<Step> stepSelections,
																				Map<Step, Map<Requirement, String>> variables) {
					return super.checkRequirement(proj.PROJECT_DIRECTORY.getValue() + arg, stepSelections,
																				variables);
				}
			};

			final RequirementSet reqSet = RequirementSetBuilder.and()
																												 .add(hmmFile)
																												 .add(RequirementSetBuilder.or()
																																									 .add(pfbStepReq)
																																									 .add(pfbFileReq))
																												 .add(RequirementSetBuilder.or()
																																									 .add(gcModelStepReq)
																																									 .add(gcModelFileReq))
																												 .add(callingTypeReq)
																												 .add(useCentroidsReq)
																												 .add(getNumThreadsReq())
																												 .add(outputFileReq);

			return register(new Step("Call CNVs", "", reqSet,
															 EnumSet.of(Requirement.Flag.MEMORY, Requirement.Flag.MULTITHREADED),
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

					CALLING_SCOPE scope = CALLING_SCOPE.valueOf(variables.get(this).get(callingTypeReq));

					(new File(ext.parseDirectoryOfFile(proj.PROJECT_DIRECTORY.getValue() + output))).mkdirs();

					String[] samples = proj.getSamples();
					boolean useCentroids = Boolean.valueOf(variables.get(this).get(useCentroidsReq));
					Centroids[] cents = new Centroids[] {null, null};
					if (useCentroids) {
						if (Files.exists(proj.SEX_CENTROIDS_FEMALE_FILENAME.getValue())
								&& Files.exists(proj.SEX_CENTROIDS_MALE_FILENAME.getValue())) {
							cents[0] = Centroids.load(proj.SEX_CENTROIDS_MALE_FILENAME.getValue());
							cents[1] = Centroids.load(proj.SEX_CENTROIDS_FEMALE_FILENAME.getValue());
						}
					}

					if (scope != CALLING_SCOPE.CHROMOSOMAL) {
						CNVCaller.callAutosomalCNVs(proj, output, samples, null, null, null,
																				CNVCaller.DEFAULT_MIN_SITES, CNVCaller.DEFAULT_MIN_CONF,
																				PFB_MANAGEMENT_TYPE.PENNCNV_DEFAULT, numThreads, 1);
						String file = proj.PROJECT_DIRECTORY.getValue() + output;
						if (Files.exists(file)) {
							proj.CNV_FILENAMES.addValue(file);
						}
					}
					if (scope != CALLING_SCOPE.AUTOSOMAL) {
						CNVCaller.callGenomeCnvs(proj, output, cents, null, CNVCaller.DEFAULT_MIN_SITES,
																		 CNVCaller.DEFAULT_MIN_CONF,
																		 PFB_MANAGEMENT_TYPE.PENNCNV_DEFAULT, numThreads, 1);
						String[] files = {
															proj.PROJECT_DIRECTORY.getValue() + output + "_23M.cnv",
															proj.PROJECT_DIRECTORY.getValue() + output + "_23F.cnv",
															proj.PROJECT_DIRECTORY.getValue() + output + "_24M.cnv"
						};
						for (String f : files) {
							if (Files.exists(f)) {
								proj.CNV_FILENAMES.addValue(f);
							}
						}
					}

					proj.saveProperties(new Property[] {proj.CNV_FILENAMES});
				}

				@Override
				public String getCommandLine(Project proj, Map<Step, Map<Requirement, String>> variables) {
					String kvCmd = "";

					String hmmP = proj.HMM_FILENAME.getValue();
					String hmmG = variables.get(this).get(hmmFile);
					if (hmmG != null && !hmmP.equals(hmmG)) {
						kvCmd += " HMM_FILENAME=" + hmmG;
					}
					String pfbP = proj.CUSTOM_PFB_FILENAME.getValue();
					String pfbG = variables.get(this).get(pfbFileReq);
					if (pfbG != null && !pfbP.equals(pfbG)) {
						kvCmd += " CUSTOM_PFB_FILENAME=" + pfbG;
					}
					String gcmP = proj.GC_MODEL_FILENAME.getValue();
					String gcmG = variables.get(this).get(gcModelFileReq);
					if (gcmG != null && !gcmP.equals(gcmG)) {
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

					boolean useCentroids = Boolean.valueOf(variables.get(this).get(useCentroidsReq));
					CALLING_SCOPE scope = CALLING_SCOPE.valueOf(variables.get(this).get(callingTypeReq));

					String autoCmd = cmd.append(Files.getRunString())
															.append(" cnv.hmm.CNVCaller proj=" + projPropFile)
															.append(" out=" + variables.get(this).get(outputFileReq)).append(" ")
															.append(PSF.Ext.NUM_THREADS_COMMAND).append(numThreads).toString();
					String genomeCmd = autoCmd + " -genome";
					if (!useCentroids) {
						genomeCmd += " -noCentroids";
					}

					switch (scope) {
						case AUTOSOMAL:
							return autoCmd;
						case CHROMOSOMAL:
							return genomeCmd;
						case BOTH:
							return autoCmd + "\n" + genomeCmd;
					}
					return autoCmd;
				}

				@Override
				public boolean checkIfOutputExists(Map<Step, Map<Requirement, String>> variables) {
					String output = variables.get(this).get(outputFileReq);
					return Files.exists(proj.PROJECT_DIRECTORY.getValue() + output);
				}

			});
		}

		private Step generatePCCorrectedProjectStep(final Step parseSamplesStep) {
			final Requirement parseSamplesStepReq = new Requirement.StepRequirement(parseSamplesStep);
			final Requirement numPCsReq = new Requirement.PosIntRequirement(
																																			"Number of principal components for correction.",
																																			MitoPipeline.DEFAULT_NUM_COMPONENTS);
			final Requirement outputBaseReq = new Requirement.OutputFileRequirement(
																																							"Output file path (relative to project directory) and baseName for principal components correction files",
																																							MitoPipeline.FILE_BASE) {
				@Override
				public boolean checkRequirement(String arg, Set<Step> stepSelections,
																				Map<Step, Map<Requirement, String>> variables) {
					String outputBase = proj.PROJECT_DIRECTORY.getValue() + arg;
					String finalReport = outputBase + PCA.FILE_EXTs[0];
					return super.checkRequirement(finalReport, stepSelections, variables);
				}
			};
			final Requirement callrateReq = new Requirement.DoubleRequirement(
																																				"Call-rate filter for determining high-quality markers",
																																				MitoPipeline.DEFAULT_MKR_CALLRATE_FILTER,
																																				0.0, 1.0);
			final Requirement recomputeLrrReq = new Requirement.OptionalBoolRequirement(
																																									"Re-compute Log-R Ratio values? (usually false if LRRs already exist)",
																																									false);
			final Requirement tempDirReq = new Requirement.OptionalFileRequirement(
																																						 "Temporary directory for intermediate files (which tend to be very large)",
																																						 "");
			final Requirement correctionStrategyReq = new Requirement.EnumRequirement("Correction Type",
																																								CORRECTION_TYPE.XY);
			final Requirement sexChromosomeStrategyReq = new Requirement.EnumRequirement(
																																									 "Sex Chromosome Strategy",
																																									 CHROMOSOME_X_STRATEGY.BIOLOGICAL);
			final Requirement setupCNVCalling = new Requirement.OptionalBoolRequirement(
																																									"Create script with steps to process corrected data and call CNVs?",
																																									false);

			final RequirementSet reqSet = RequirementSetBuilder.and()
																												 .add(parseSamplesStepReq)
																												 .add(numPCsReq)
																												 .add(outputBaseReq)
																												 .add(callrateReq)
																												 .add(recomputeLrrReq)
																												 .add(tempDirReq)
																												 .add(correctionStrategyReq)
																												 .add(sexChromosomeStrategyReq)
																												 .add(getNumThreadsReq())
																												 .add(setupCNVCalling);

			return register(new Step("Create PC-Corrected Project", "", reqSet,
															 EnumSet.of(Requirement.Flag.MEMORY, Requirement.Flag.RUNTIME,
																					Requirement.Flag.MULTITHREADED), priority()) {

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
					boolean cnvCalling = Boolean.parseBoolean(variables.get(this).get(setupCNVCalling));
					String retMsg = PRoCtOR.shadow(proj, tmpDir, outputBase, markerCallRateFilter,
																				 recomputeLRRPCs, type, strategy, numComponents,
																				 totalThreads, cnvCalling);
					if (retMsg != null && !"".equals(retMsg)) {
						throw new RuntimeException(retMsg);
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

					boolean cnvCalling = Boolean.parseBoolean(variables.get(this).get(setupCNVCalling));
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
					if (cnvCalling) {
						cmd.append(" -callCNVs");
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
			final Requirement parseSamplesStepReq = new Requirement.StepRequirement(parseSamplesStep);
			final RequirementSet reqSet = RequirementSetBuilder.and().add(parseSamplesStepReq);
			return register(new Step("Generate AB Lookup File", "", reqSet,
															 EnumSet.of(Requirement.Flag.RUNTIME), priority()) {

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
						throw new RuntimeException(
																			 "Failed to fill in missing alleles - please check log for more info.");
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

	private SortedSet<Step> generateSteps(boolean allowCorrectionStep) {
		StepBuilder sb = new StepBuilder();

		Step parseSamplesStep;
		Step markerBlastStep;
		if (proj.getArrayType() == ARRAY.AFFY_GW6 || proj.getArrayType() == ARRAY.AFFY_GW6_CN
				|| proj.getArrayType() == ARRAY.AFFY_AXIOM) {
			parseSamplesStep = sb.generateParseSamplesStep();
			markerBlastStep = sb.generateAffyMarkerBlastAnnotationStep(parseSamplesStep);
		} else {
			Step markerPositions = sb.generateIlluminaMarkerPositionsStep();
			parseSamplesStep = sb.generateParseSamplesStep(markerPositions);
			markerBlastStep = sb.generateIlluminaMarkerBlastAnnotationStep(parseSamplesStep);
		}
		Step createSampleDataStep = sb.generateCreateSampleDataStep(parseSamplesStep);
		Step transposeStep = sb.generateTransposeStep(parseSamplesStep);
		Step gcModelStep = sb.generateGCModelStep();
		Step sampleQCStep = sb.generateSampleQCStep(parseSamplesStep);
		sb.generateMarkerQCStep(parseSamplesStep);
		sb.generateSexChecksStep(parseSamplesStep, markerBlastStep, createSampleDataStep,
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
		sb.generateMitoCNEstimateStep(transposeStep);
		Step pfbStep = sb.generatePFBStep(parseSamplesStep);
		sb.generateSexCentroidsStep();
		sb.generateCNVStep(pfbStep, gcModelStep);
		if (allowCorrectionStep) {
			sb.generatePCCorrectedProjectStep(parseSamplesStep);
		}

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


	public static void setupImputation(String projectProperties) {
		Project proj = new Project(projectProperties);
		StepBuilder sb = (new GenvisisWorkflow(proj, null)).new StepBuilder();

		Step parseSamples = sb.generateParseSamplesStep();
		Step transpose = sb.generateTransposeStep(parseSamples);
		Step sampleData = sb.generateCreateSampleDataStep(parseSamples);
		Step blast = sb.generateAffyMarkerBlastAnnotationStep(parseSamples);
		Step sampleQc = sb.generateSampleQCStep(parseSamples);
		Step markerQc = sb.generateMarkerQCStep(parseSamples);
		Step sexChecks = sb.generateSexChecksStep(parseSamples, blast, sampleData, transpose,
																							sampleQc);
		Step exportPlink = sb.generatePlinkExportStep(parseSamples);
		Step gwasQc = sb.generateGwasQCStep(exportPlink);
		Step ancestry = sb.generateAncestryStep(gwasQc);
		Step faqcStep = sb.generateFurtherAnalysisQCStep(exportPlink, gwasQc, ancestry);

		Map<Step, Map<Requirement, String>> varMap = new HashMap<>();
		varMap.put(sampleQc, sampleQc.getDefaultRequirementValues());
		varMap.put(markerQc, markerQc.getDefaultRequirementValues());
		varMap.put(sexChecks, sexChecks.getDefaultRequirementValues());
		varMap.put(exportPlink, exportPlink.getDefaultRequirementValues());
		varMap.put(gwasQc, gwasQc.getDefaultRequirementValues());
		varMap.put(ancestry, ancestry.getDefaultRequirementValues());
		varMap.put(faqcStep, faqcStep.getDefaultRequirementValues());

		String s1 = sampleQc.getCommandLine(proj, varMap);
		String s2 = markerQc.getCommandLine(proj, varMap);
		String s3 = sexChecks.getCommandLine(proj, varMap);
		String s4 = exportPlink.getCommandLine(proj, varMap);
		String s5 = gwasQc.getCommandLine(proj, varMap);
		String s6 = ancestry.getCommandLine(proj, varMap);
		String s7 = faqcStep.getCommandLine(proj, varMap);

		String file = proj.PROJECT_DIRECTORY.getValue() + "ImputationPipeline.";
		String suggFile = file + ext.getTimestampForFilename() + ".pbs";
		String runFile = file + ext.getTimestampForFilename() + ".run";

		StringBuilder output = new StringBuilder("## Genvisis Project Pipeline - Stepwise Commands\n\n");

		addStepInfo(output, sampleQc, s1);
		addStepInfo(output, markerQc, s2);
		addStepInfo(output, sexChecks, s3);
		addStepInfo(output, exportPlink, s4);
		addStepInfo(output, gwasQc, s5);
		addStepInfo(output, ancestry, s6);
		addStepInfo(output, faqcStep, s7);

		Qsub.qsubDefaults(suggFile, output.toString());
		Files.write(output.toString(), runFile);
	}

	public static void setupCNVCalling(String projectProperties) {
		Project pcProj = new Project(projectProperties);
		StepBuilder sb = (new GenvisisWorkflow(pcProj, null)).new StepBuilder();
		Step transpose = sb.generateTransposeStep(null);
		// Create new sample data, run sex checks?
		Step gc = sb.generateGCModelStep();
		Step pfb = sb.generatePFBStep(null);
		Step cent = sb.generateSexCentroidsStep();
		Step cnv = sb.generateCNVStep(pfb, gc);
		Map<Step, Map<Requirement, String>> stepOpts = new HashMap<>();
		HashMap<Requirement, String> cnvOpts = new HashMap<>();
		List<Requirement> reqs = cnv.getRequirements().getFlatRequirementsList();
		for (Requirement req : reqs) {
			if (req.getDescription().equals(CNVCaller.CNV_SCOPE_DESC)) {
				cnvOpts.put(req, CNVCaller.CALLING_SCOPE.BOTH.toString());
			} else if (req.getDescription().equals(NUM_THREADS_DESC)) {
				cnvOpts.put(req, "" + (Runtime.getRuntime().availableProcessors() - 1));
			}
		}
		stepOpts.put(cnv, cnvOpts);

		String s1 = transpose.getCommandLine(pcProj, null);
		String s2 = gc.getCommandLine(pcProj, null);
		String s3 = pfb.getCommandLine(pcProj, null);
		String s4 = cent.getCommandLine(pcProj, null);
		String s5 = cnv.getCommandLine(pcProj, stepOpts);

		String file = pcProj.PROJECT_DIRECTORY.getValue() + "CNVCallingPipeline.";
		String suggFile = file + ext.getTimestampForFilename() + ".pbs";
		String runFile = file + ext.getTimestampForFilename() + ".run";

		StringBuilder output = new StringBuilder("## Genvisis Project Pipeline - Stepwise Commands\n\n");

		addStepInfo(output, transpose, s1);
		addStepInfo(output, gc, s2);
		addStepInfo(output, pfb, s3);
		addStepInfo(output, cent, s4);
		addStepInfo(output, cnv, s5);

		Qsub.qsubDefaults(suggFile, output.toString());
		Files.write(output.toString(), runFile);
	}

	public static void addStepInfo(StringBuilder output, Step step, String stepCmd) {
		output.append("## ").append(step.getName()).append("\n");
		output.append("echo \" start ").append(step.getName()).append(" at: \" `date`")
					.append("\n");
		output.append(stepCmd).append("\n");
		output.append("echo \" end ").append(step.getName()).append(" at: \" `date`")
					.append("\n");
		output.append("\n\n");
	}


}
