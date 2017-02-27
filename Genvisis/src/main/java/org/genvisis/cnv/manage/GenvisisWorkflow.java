package org.genvisis.cnv.manage;

import java.io.File;
import java.util.Collections;
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
import org.genvisis.cnv.manage.GenvisisWorkflowStep.BoolRequirement;
import org.genvisis.cnv.manage.GenvisisWorkflowStep.DoubleRequirement;
import org.genvisis.cnv.manage.GenvisisWorkflowStep.EnumRequirement;
import org.genvisis.cnv.manage.GenvisisWorkflowStep.FLAG;
import org.genvisis.cnv.manage.GenvisisWorkflowStep.FileRequirement;
import org.genvisis.cnv.manage.GenvisisWorkflowStep.IntRequirement;
import org.genvisis.cnv.manage.GenvisisWorkflowStep.OptionalBoolRequirement;
import org.genvisis.cnv.manage.GenvisisWorkflowStep.OptionalFileRequirement;
import org.genvisis.cnv.manage.GenvisisWorkflowStep.OutputFileRequirement;
import org.genvisis.cnv.manage.GenvisisWorkflowStep.Requirement;
import org.genvisis.cnv.manage.GenvisisWorkflowStep.RequirementInputType;
import org.genvisis.cnv.manage.GenvisisWorkflowStep.StepRequirement;
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

	private static final String PROJ_PROP_UPDATE_STR = " org.genvisis.cnv.filesys.Project proj=";
	final Project proj;
	private final SortedSet<GenvisisWorkflowStep> steps;
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
		numThreadsReq = new Requirement("Number of Threads to Use", RequirementInputType.NUMBER,
																		proj.NUM_THREADS.getValue()) {

			@Override
			public boolean checkRequirement(String arg, Set<GenvisisWorkflowStep> stepSelections,
																			Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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

	private SortedSet<GenvisisWorkflowStep> generateSteps() {
		SortedSet<GenvisisWorkflowStep> buildSteps = Sets.newTreeSet();
		double priority = 1.0;
		GenvisisWorkflowStep parseSamplesStep;
		if (proj.getArrayType() == ARRAY.AFFY_GW6 || proj.getArrayType() == ARRAY.AFFY_GW6_CN) {
			parseSamplesStep = generateAffyParseSamplesStep(priority++);
		} else {
			GenvisisWorkflowStep markerPositions = generateMarkerPositionsStep(priority++);
			buildSteps.add(markerPositions);
			parseSamplesStep = generateIlluminaParseSamplesStep(priority++, markerPositions);
		}
		GenvisisWorkflowStep createSampleDataStep = generateCreateSampleDataStep(priority++,
																																						 parseSamplesStep);
		GenvisisWorkflowStep transposeStep = generateTransposeStep(priority++, parseSamplesStep);
		GenvisisWorkflowStep gcModelStep = generateGCModelStep(priority++);
		GenvisisWorkflowStep sampleQCStep = generateSampleQCStep(priority++, parseSamplesStep);
		GenvisisWorkflowStep markerQCStep = generateMarkerQCStep(priority++, parseSamplesStep);
		GenvisisWorkflowStep sexChecksStep = generateSexChecksStep(priority++, parseSamplesStep,
																															 createSampleDataStep, transposeStep,
																															 sampleQCStep);
		GenvisisWorkflowStep abLookupStep = generateABLookupStep(priority++, parseSamplesStep);
		GenvisisWorkflowStep plinkExportStep = generatePlinkExportStep(priority++, parseSamplesStep);
		GenvisisWorkflowStep gwasQCStep = generateGwasQCStep(priority++, plinkExportStep);
		GenvisisWorkflowStep mosaicArmsStep = generateMosaicArmsStep(priority++, parseSamplesStep);
		GenvisisWorkflowStep annotateSampleDataStep = generateAnnotateSampleDataStep(priority++,
																																								 sampleQCStep,
																																								 createSampleDataStep,
																																								 gwasQCStep);
		GenvisisWorkflowStep createPCsStep = generateCreatePCsStep(priority++, transposeStep);
		GenvisisWorkflowStep pfbStep = generatePFBStep(priority++, parseSamplesStep);
		GenvisisWorkflowStep sexCentroidsStep = generateSexCentroidsStep(priority++, gcModelStep);
		GenvisisWorkflowStep cnvStep = generateCNVStep(priority++, pfbStep, gcModelStep);
		GenvisisWorkflowStep shadowStep = generateShadowStep(priority++, parseSamplesStep);

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

	private GenvisisWorkflowStep generateMarkerPositionsStep(final double priority) {
		final Requirement snpMapReq = new FileRequirement("An Illumina SNP_map file.",
																											proj.getLocationOfSNP_Map(false));
		final Requirement manifestReq = new FileRequirement("An Illumina Manifest file.",
																												proj.getLocationOfSNP_Map(false));
		return new GenvisisWorkflowStep("Create Marker Positions (if not already exists)", "",
																		new Requirement[][] {{snpMapReq}, {manifestReq}}, null, null,
																		priority) {

			@Override
			public void setNecessaryPreRunProperties(Project proj,
																							 Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
				// not needed for step
			}

			@Override
			public void run(Project proj, Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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
			public boolean checkIfOutputExists(Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
				return Files.exists(proj.MARKER_POSITION_FILENAME.getValue(false, false));
			}

			@Override
			public String getCommandLine(Project proj,
																	 Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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

	private GenvisisWorkflowStep generateIlluminaParseSamplesStep(double priority,
																																final GenvisisWorkflowStep markerPositionsStep) {
		final Requirement markerPositionsStepReq = new StepRequirement(markerPositionsStep);

		final Requirement markerPositionsReq = new FileRequirement("Parsed markerPositions file must already exist.",
																															 proj.MARKER_POSITION_FILENAME.getValue(false,
																																																			false));

		return new GenvisisWorkflowStep("Parse Illumina Sample Files", "",
																		new Requirement[][] {{markerPositionsStepReq,
																													markerPositionsReq},
																												 {numThreadsReq}},
																		new GenvisisWorkflowStep[] {markerPositionsStep},
																		new GenvisisWorkflowStep.FLAG[] {GenvisisWorkflowStep.FLAG.MEMORY,
																																		 GenvisisWorkflowStep.FLAG.RUNTIME,
																																		 GenvisisWorkflowStep.FLAG.MEMORY},
																		priority) {

			@Override
			public void setNecessaryPreRunProperties(Project proj,
																							 Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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
			public void run(Project proj, Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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
			public boolean checkIfOutputExists(Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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
																	 Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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

	private GenvisisWorkflowStep generateAffyParseSamplesStep(double priority) {

		final Requirement markerPositionsReq = new Requirement("markerPositions file must already exist.",
																													 RequirementInputType.FILE,
																													 proj.MARKER_POSITION_FILENAME.getValue(false,
																																																	false)) {

			@Override
			public boolean checkRequirement(String arg, Set<GenvisisWorkflowStep> stepSelections,
																			Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
				return Files.exists(arg);
			}
		};

		return new GenvisisWorkflowStep("Parse Sample Files", "",
																		new Requirement[][] {{markerPositionsReq},
																												 {numThreadsReq}},
																		null,
																		new GenvisisWorkflowStep.FLAG[] {GenvisisWorkflowStep.FLAG.MEMORY,
																																		 GenvisisWorkflowStep.FLAG.RUNTIME,
																																		 GenvisisWorkflowStep.FLAG.MEMORY},
																		priority) {

			@Override
			public void setNecessaryPreRunProperties(Project proj,
																							 Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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
			public void run(Project proj, Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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
			public boolean checkIfOutputExists(Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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
																	 Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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

	private GenvisisWorkflowStep generateCreateSampleDataStep(double priority,
																														final GenvisisWorkflowStep parseSamplesStep) {
		final Requirement parseSamplesStepReq = new StepRequirement(parseSamplesStep);

		final Requirement createMinimalSampleDataReq = new BoolRequirement("Create a minimal SampleData.txt file from sample files",
																																			 true);

		final String pedPreset = proj.PEDIGREE_FILENAME.getValue();

		final Requirement pedigreeReq = new Requirement("Either a Pedigree.dat file, or any file with a header containing all of the following elements (in any order):  \""
																										+ ArrayUtils.toStr(MitoPipeline.PED_INPUT, ", ")
																										+ "\"", RequirementInputType.FILE,
																										pedPreset) {

			@Override
			public boolean checkRequirement(String arg, Set<GenvisisWorkflowStep> stepSelections,
																			Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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
			public boolean checkRequirement(String arg, Set<GenvisisWorkflowStep> stepSelections,
																			Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
				return Files.exists(arg);
			}
		};
		return new GenvisisWorkflowStep("Create SampleData.txt File", "",
																		new Requirement[][] {{parseSamplesStepReq},
																												 {createMinimalSampleDataReq, pedigreeReq,
																													sampMapReq}},
																		new GenvisisWorkflowStep[] {parseSamplesStep},
																		new GenvisisWorkflowStep.FLAG[] {},
																		priority) {

			@Override
			public void setNecessaryPreRunProperties(Project proj,
																							 Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
				// Nothing to do
			}

			@Override
			public void run(Project proj, Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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
			public boolean checkIfOutputExists(Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
				return Files.exists(proj.SAMPLE_DATA_FILENAME.getValue(false, false));
			}

			@Override
			public String getCommandLine(Project proj,
																	 Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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

	private GenvisisWorkflowStep generateTransposeStep(double priority,
																										 final GenvisisWorkflowStep parseSamplesStep) {
		return new GenvisisWorkflowStep("Transpose Data into Marker-Dominant Files", "",
																		new Requirement[][] {{new StepRequirement(parseSamplesStep)}},
																		new GenvisisWorkflowStep[] {parseSamplesStep},
																		new FLAG[] {FLAG.MEMORY}, priority) {

			@Override
			public void setNecessaryPreRunProperties(Project proj,
																							 Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
				// Nothing to do here
			}

			@Override
			public void run(Project proj, Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
				proj.getLog().report("Transposing data");
				TransposeData.transposeData(proj, 2000000000, false); // compact if no LRR was provided
			}

			@Override
			public String getCommandLine(Project proj,
																	 Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
				String projPropFile = proj.getPropertyFilename();
				StringBuilder cmd = new StringBuilder();
				return cmd.append(Files.getRunString()).append(" cnv.manage.TransposeData -transpose proj="
																											 + projPropFile + " max=" + 2000000000)
									.toString();
			}

			@Override
			public boolean checkIfOutputExists(Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
				return Files.list(proj.MARKER_DATA_DIRECTORY.getValue(false, false),
													MarkerData.MARKER_DATA_FILE_EXTENSION, false).length > 0;
			}
		};
	}

	private GenvisisWorkflowStep generateGCModelStep(double priority) {
		final Requirement gcBaseReq = new FileRequirement("A GC Base file must exist.",
																											Resources.genome(proj.GENOME_BUILD_VERSION.getValue(),
																																			 proj.getLog())
																															 .getModelBase().get());
		final Requirement gcModelOutputReq = new OutputFileRequirement("GCModel output file must be specified.",
																																	 proj.GC_MODEL_FILENAME.getValue());

		return new GenvisisWorkflowStep("Compute GCMODEL File", "",
																		new Requirement[][] {{gcBaseReq},
																												 {gcModelOutputReq}},
																		null, new FLAG[] {}, priority) {

			@Override
			public void setNecessaryPreRunProperties(Project proj,
																							 Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
				String setGCOutputFile = proj.GC_MODEL_FILENAME.getValue();
				String gcOutputFile = variables.get(this).get(gcModelOutputReq);
				if (!ext.verifyDirFormat(setGCOutputFile).equals(gcOutputFile)) {
					proj.GC_MODEL_FILENAME.setValue(gcOutputFile);
				}
			}

			@Override
			public void run(Project proj, Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
				String gcBaseFile = variables.get(this).get(gcBaseReq);
				String gcOutputFile = variables.get(this).get(gcModelOutputReq);
				org.genvisis.cnv.analysis.PennCNV.gcModel(proj, gcBaseFile, gcOutputFile, 100);
			}

			@Override
			public String getCommandLine(Project proj,
																	 Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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
			public boolean checkIfOutputExists(Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
				String gcOutputFile = variables.get(this).get(gcModelOutputReq);
				return Files.exists(gcOutputFile);
			}

		};
	}

	private GenvisisWorkflowStep generateSampleQCStep(double priority,
																										final GenvisisWorkflowStep parseSamplesStep) {
		final Requirement parseSamplesStepReq = new StepRequirement(parseSamplesStep);

		return new GenvisisWorkflowStep("Run Sample QC Metrics", "",
																		new Requirement[][] {{parseSamplesStepReq},
																												 {numThreadsReq},},
																		new GenvisisWorkflowStep[] {parseSamplesStep},
																		new FLAG[] {FLAG.MULTITHREADED}, priority) {

			@Override
			public void setNecessaryPreRunProperties(Project proj,
																							 Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
				int numThreads = resolveThreads(variables.get(this).get(numThreadsReq));
				maybeSetProjNumThreads(numThreads);
			}

			@Override
			public void run(Project proj, Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
				proj.getLog().report("Running LrrSd");
				int numThreads = proj.NUM_THREADS.getValue();
				LrrSd.init(proj, null, null, numThreads, false);
			}

			@Override
			public String getCommandLine(Project proj,
																	 Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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
			public boolean checkIfOutputExists(Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
				return Files.exists(proj.SAMPLE_QC_FILENAME.getValue(false, false));
			}
		};
	}

	private GenvisisWorkflowStep generateMarkerQCStep(double priority,
																										final GenvisisWorkflowStep parseSamplesStep) {
		final Requirement parseSamplesStepReq = new StepRequirement(parseSamplesStep);

		final Requirement exportAllReq = new OptionalBoolRequirement("Export all markers in project.",
																																 true);

		final Requirement targetMarkersReq = new FileRequirement("A targetMarkers files listing the markers to QC.",
																														 proj.TARGET_MARKERS_FILENAMES.getValue()[0]);

		return new GenvisisWorkflowStep("Run Marker QC Metrics", "",
																		new Requirement[][] {{parseSamplesStepReq},
																												 {exportAllReq, targetMarkersReq},
																												 {numThreadsReq}},
																		new GenvisisWorkflowStep[] {parseSamplesStep},
																		new FLAG[] {FLAG.MULTITHREADED},
																		priority) {

			@Override
			public void setNecessaryPreRunProperties(Project proj,
																							 Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
				// Nothing to do here
			}

			@Override
			public void run(Project proj, Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
				boolean allMarkers = Boolean.parseBoolean(variables.get(this).get(exportAllReq));
				String tgtFile = allMarkers ? null : variables.get(this).get(targetMarkersReq);
				int numThreads = resolveThreads(variables.get(this).get(numThreadsReq));
				MarkerMetrics.fullQC(proj, null, tgtFile, true, numThreads);
			}

			@Override
			public String getCommandLine(Project proj,
																	 Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
				boolean allMarkers = Boolean.parseBoolean(variables.get(this).get(exportAllReq));
				String tgtFile = variables.get(this).get(targetMarkersReq);
				int numThreads = resolveThreads(variables.get(this).get(numThreadsReq));
				return Files.getRunString() + " cnv.qc.MarkerMetrics -fullQC" + " proj="
							 + proj.getPropertyFilename() + (allMarkers ? "" : " markers=" + tgtFile) + " "
							 + PSF.Ext.NUM_THREADS_COMMAND + numThreads;
			}

			@Override
			public boolean checkIfOutputExists(Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
				String markerMetricsFile = proj.MARKER_METRICS_FILENAME.getValue();
				return Files.exists(markerMetricsFile);
			}

		};
	}

	private GenvisisWorkflowStep generateSexChecksStep(double priority,
																										 final GenvisisWorkflowStep parseSamplesStep,
																										 final GenvisisWorkflowStep sampleDataStep,
																										 final GenvisisWorkflowStep transposeStep,
																										 final GenvisisWorkflowStep sampleQCStep) {
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

		return new GenvisisWorkflowStep("Run Sex Checks", "",
																		new Requirement[][] {{parseSamplesStepReq},
																												 {sampleDataStepReq},
																												 {transposeStepReq},
																												 {sampleQCStepReq},
																												 {addToSampleDataReq},
																												 {noCrossHybeReq,
																													oneHittersReq,
																													blastVCFReq}},
																		new GenvisisWorkflowStep[] {parseSamplesStep, sampleDataStep,
																																transposeStep, sampleQCStep},
																		new FLAG[] {},
																		priority) {

			@Override
			public void setNecessaryPreRunProperties(Project proj,
																							 Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
				// Nothing to do here
			}

			@Override
			public void run(Project proj, Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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
																	 Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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
			public boolean checkIfOutputExists(Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
				return Files.exists(proj.SEXCHECK_RESULTS_FILENAME.getValue());
			}

		};
	}

	private GenvisisWorkflowStep generateABLookupStep(double priority,
																										final GenvisisWorkflowStep parseSamplesStep) {
		final Requirement parseSamplesStepReq = new StepRequirement(parseSamplesStep);
		return new GenvisisWorkflowStep("Generate AB Lookup File", "",
																		new Requirement[][] {{parseSamplesStepReq}},
																		new GenvisisWorkflowStep[] {parseSamplesStep},
																		new FLAG[] {FLAG.RUNTIME}, priority) {

			@Override
			public void setNecessaryPreRunProperties(Project proj,
																							 Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
				// Nothing to do here
			}

			@Override
			public void run(Project proj, Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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
																	 Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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
			public boolean checkIfOutputExists(Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
				return Files.exists(proj.AB_LOOKUP_FILENAME.getValue(false, false));
			}
		};

	}

	private GenvisisWorkflowStep generatePlinkExportStep(double priority,
																											 final GenvisisWorkflowStep parseSamplesStep) {
		final Requirement parseSamplesStepReq = new StepRequirement(parseSamplesStep);
		final Requirement pedigreeRequirement = new FileRequirement("A pedigree.dat file must exist.",
																																proj.PEDIGREE_FILENAME.getValue(false,
																																																false));
		final Requirement createPedigreeRequirement = new BoolRequirement("Create a minimal pedigree.dat file [will pull information from SexChecks step results].",
																																			false);

		return new GenvisisWorkflowStep("Create PLINK Files", "",
																		new Requirement[][] {{parseSamplesStepReq},
																												 {pedigreeRequirement,
																													createPedigreeRequirement}},
																		new GenvisisWorkflowStep[] {parseSamplesStep},
																		new FLAG[] {FLAG.MEMORY}, priority) {

			@Override
			public void setNecessaryPreRunProperties(Project proj,
																							 Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
				if (!Boolean.parseBoolean(variables.get(this).get(createPedigreeRequirement))) {
					String projPedFile = proj.PEDIGREE_FILENAME.getValue(false, false);
					String pedFile = variables.get(this).get(pedigreeRequirement);
					if (!pedFile.equals(projPedFile)) {
						proj.PEDIGREE_FILENAME.setValue(pedFile);
					}
				}
			}

			@Override
			public void run(Project proj, Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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
																	 Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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
			public boolean checkIfOutputExists(Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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

	private GenvisisWorkflowStep generateGwasQCStep(double priority,
																									GenvisisWorkflowStep plinkExportStep) {
		// TODO: Move Ancestry to its own step
		final Requirement plinkExportStepReq = new StepRequirement(plinkExportStep);
		final Requirement genomeForRelatedsReq = new OptionalBoolRequirement("Keep genome info for unrelateds only",
																																				 false);
		final Requirement skipAncestryReq = new BoolRequirement("Skip ancestry checks", false);
		final Requirement putativeWhitesReq = new FileRequirement("File with FID/IID pairs of putative white samples",
																															"");
		final Requirement hapMapFoundersReq = new FileRequirement("PLINK root of HapMap founders",
																															Ancestry.DEFAULT_HAPMAP_PLINKROOT) {
			@Override
			public boolean checkRequirement(String arg, Set<GenvisisWorkflowStep> stepSelections,
																			Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
				String hapMapPlinkRoot = arg;
				int dotIndex = hapMapPlinkRoot.lastIndexOf('.');
				if (dotIndex > 0
						&& PSF.Plink.getPlinkBedBimFamSet("").contains(hapMapPlinkRoot.substring(dotIndex))) {
					hapMapPlinkRoot = hapMapPlinkRoot.substring(0, dotIndex);
				}
				return Files.checkAllFiles("", PSF.Plink.getPlinkBedBimFamSet(hapMapPlinkRoot), false, log);
			}
		};

		return new GenvisisWorkflowStep("Run GWAS QC", "",
																		new Requirement[][] {{plinkExportStepReq},
																												 {genomeForRelatedsReq},
																												 {skipAncestryReq,
																													putativeWhitesReq,
																													hapMapFoundersReq}},
																		new GenvisisWorkflowStep[] {plinkExportStep}, new FLAG[] {},
																		priority) {

			@Override
			public void setNecessaryPreRunProperties(Project proj,
																							 Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
				// not needed for step
			}

			@Override
			public void run(Project proj, Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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
																	 Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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
			public boolean checkIfOutputExists(Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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

	private GenvisisWorkflowStep generateMosaicArmsStep(double priority,
																											final GenvisisWorkflowStep parseSamplesStep) {
		final Requirement parseSamplesStepReq = new StepRequirement(parseSamplesStep);
		return new GenvisisWorkflowStep("Create Mosaic Arms File", "",
																		new Requirement[][] {{parseSamplesStepReq},
																												 {numThreadsReq}},
																		new GenvisisWorkflowStep[] {parseSamplesStep},
																		new FLAG[] {FLAG.MULTITHREADED}, priority) {

			@Override
			public void setNecessaryPreRunProperties(Project proj,
																							 Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {

				int numThreads = resolveThreads(variables.get(this).get(numThreadsReq));
				maybeSetProjNumThreads(numThreads);
			}

			@Override
			public void run(Project proj, Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
				Mosaicism.findOutliers(proj);
			}

			@Override
			public String getCommandLine(Project proj,
																	 Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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
			public boolean checkIfOutputExists(Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
				return Files.exists(proj.RESULTS_DIRECTORY.getValue(false, false) + "Mosaicism.xln");
			}
		};
	}

	private GenvisisWorkflowStep generateAnnotateSampleDataStep(double priority,
																															final GenvisisWorkflowStep sampleQCStep,
																															final GenvisisWorkflowStep createSampleDataStep,
																															final GenvisisWorkflowStep gwasQCStep) {
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
			public boolean checkRequirement(String arg, Set<GenvisisWorkflowStep> stepSelections,
																			Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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
		final Requirement numQReq = new IntRequirement("Number of Quantiles to Generate", 10, 1,
																									 Integer.MAX_VALUE);
		final Requirement replaceFIDIIDReq = new OptionalBoolRequirement("Replace FID and IID with data from Pedigree",
																																		 false);
		return new GenvisisWorkflowStep("Annotate Sample Data File", "",
																		new Requirement[][] {{sampleQCStepReq},
																												 {createSampleDataStepReq},
																												 {skipIDingDuplicatesReq, gwasQCStepReq},
																												 {notGcCorrectedLrrSdReq,
																													gcCorrectedLrrSdReq},
																												 {lrrSdThresholdReq},
																												 {callrateThresholdReq},
																												 {numQReq},
																												 {replaceFIDIIDReq}},
																		new GenvisisWorkflowStep[] {sampleQCStep,
																																createSampleDataStep,
																																gwasQCStep},
																		new FLAG[] {}, priority) {

			@Override
			public void setNecessaryPreRunProperties(Project proj,
																							 Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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
			public void run(Project proj, Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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
																	 Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {

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
			public boolean checkIfOutputExists(Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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

	private GenvisisWorkflowStep generateCreatePCsStep(double priority,
																										 GenvisisWorkflowStep transposeStep) {
		// FIXME http://genvisis.org/MitoPipeline/#illumina_marker_lists has illumina markers.. this
		// should be linked to, or
		// these steps split or something...
		final Requirement transposeStepReq = new StepRequirement(transposeStep);
		final Requirement medianMarkersReq = new FileRequirement("MedianMarkers file must exist.", "");
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
		final Requirement gcRegressionDistanceReq = new IntRequirement("Regression distance for the GC adjustment",
																																	 GcAdjustor.DEFAULT_REGRESSION_DISTANCE[0],
																																	 1, Integer.MAX_VALUE);
		final Requirement pcSelectionSamplesReq = new OptionalFileRequirement("A file listing a subset of samples (DNA ID) to use for determining optimal PC selection, typically a list of unrelated and single race samples. If a list is not provided, only samples passing sample qc thresholds will be used.",
																																					"");
		final Requirement externalBetaFileReq = new OptionalFileRequirement("An external beta file to optimize PC selection.",
																																				"");

		return new GenvisisWorkflowStep("Create Principal Components File and Mitochondrial Copy-Number Estimates File",
																		"",
																		new Requirement[][] {{transposeStepReq},
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
																		new GenvisisWorkflowStep[] {transposeStep},
																		new FLAG[] {FLAG.MULTITHREADED}, priority) {

			@Override
			public void setNecessaryPreRunProperties(Project proj,
																							 Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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
			public void run(Project proj, Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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
																	 Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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
			public boolean checkIfOutputExists(Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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

	private GenvisisWorkflowStep generatePFBStep(double priority,
																							 final GenvisisWorkflowStep parseSamplesStep) {
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
		return new GenvisisWorkflowStep("Compute Population BAF files", "",
																		new Requirement[][] {{parseSamplesStepReq,
																													sampleSubsetReq},
																												 {outputFileReq}},
																		new GenvisisWorkflowStep[] {parseSamplesStep},
																		new FLAG[] {}, priority) {

			@Override
			public void setNecessaryPreRunProperties(Project proj,
																							 Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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
			public void run(Project proj, Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
				org.genvisis.cnv.analysis.PennCNV.populationBAF(proj);
			}

			@Override
			public String getCommandLine(Project proj,
																	 Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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
			public boolean checkIfOutputExists(Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
				String subSampFile = variables.get(this).get(sampleSubsetReq);
				String pfbOutputFile = variables.get(this).get(outputFileReq);
				return Files.exists(pfbOutputFile) || Files.exists(ext.rootOf(subSampFile) + ".pfb");
			}
		};
	}

	private GenvisisWorkflowStep generateSexCentroidsStep(double priority,
																												GenvisisWorkflowStep gcModelStep) {
		final Requirement gcModelStepReq = new StepRequirement(gcModelStep);
		final Requirement gcModelFileReq = new FileRequirement("Full GC Model File.",
																													 proj.GC_MODEL_FILENAME.getValue());
		return new GenvisisWorkflowStep("Create Sex-Specific Centroids; Filter PFB and GCMODEL Files",
																		"",
																		new Requirement[][] {{gcModelStepReq,
																													gcModelFileReq},
																												 {numThreadsReq},},
																		new GenvisisWorkflowStep[] {gcModelStep},
																		new FLAG[] {FLAG.RUNTIME,
																								FLAG.MULTITHREADED,},
																		priority) {

			@Override
			public void setNecessaryPreRunProperties(Project proj,
																							 Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {

				int numThreads = resolveThreads(variables.get(this).get(numThreadsReq));
				maybeSetProjNumThreads(numThreads);
			}

			@Override
			public void run(Project proj, Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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
																	 Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {

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
			public boolean checkIfOutputExists(Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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

	private GenvisisWorkflowStep generateCNVStep(double priority, GenvisisWorkflowStep pfbStep,
																							 GenvisisWorkflowStep gcModelStep) {
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
			public boolean checkRequirement(String arg, Set<GenvisisWorkflowStep> stepSelections,
																			Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
				return super.checkRequirement(proj.PROJECT_DIRECTORY.getValue() + arg, stepSelections,
																			variables);
			}
		};

		return new GenvisisWorkflowStep("Call CNVs", "",
																		new Requirement[][] {{hmmFile},
																												 {pfbStepReq,
																													pfbFileReq},
																												 {gcModelStepReq,
																													gcModelFileReq},
																												 {numThreadsReq},
																												 {outputFileReq}},
																		new GenvisisWorkflowStep[] {pfbStep, gcModelStep},
																		new FLAG[] {FLAG.MEMORY, FLAG.MULTITHREADED,}, priority) {

			@Override
			public void setNecessaryPreRunProperties(Project proj,
																							 Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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
			public void run(Project proj, Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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
																	 Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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
			public boolean checkIfOutputExists(Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
				String output = variables.get(this).get(outputFileReq);
				return Files.exists(proj.PROJECT_DIRECTORY.getValue() + output);
			}
		};
	}

	private GenvisisWorkflowStep generateShadowStep(double priority,
																									final GenvisisWorkflowStep parseSamplesStep) {
		final Requirement parseSamplesStepReq = new StepRequirement(parseSamplesStep);
		final Requirement numPCsReq = new IntRequirement("Number of principal components for correction.",
																										 MitoPipeline.DEFAULT_NUM_COMPONENTS, 1,
																										 Integer.MAX_VALUE);
		final Requirement outputBaseReq = new OutputFileRequirement("Output file path (relative to project directory) and baseName for principal components correction files",
																																MitoPipeline.FILE_BASE) {
			@Override
			public boolean checkRequirement(String arg, Set<GenvisisWorkflowStep> stepSelections,
																			Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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

		return new GenvisisWorkflowStep("Create 'Shadow' Project", "",
																		new Requirement[][] {{parseSamplesStepReq},
																												 {numPCsReq},
																												 {outputBaseReq},
																												 {callrateReq},
																												 {recomputeLrrReq},
																												 {tempDirReq},
																												 {correctionStrategyReq},
																												 {sexChromosomeStrategyReq},
																												 {numThreadsReq},},
																		new GenvisisWorkflowStep[] {parseSamplesStep},
																		new FLAG[] {FLAG.MEMORY, FLAG.RUNTIME}, priority) {

			@Override
			public void setNecessaryPreRunProperties(Project proj,
																							 Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
				// not needed for step
			}

			@Override
			public void run(Project proj, Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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
																	 Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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
			public boolean checkIfOutputExists(Map<GenvisisWorkflowStep, Map<Requirement, String>> variables) {
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
