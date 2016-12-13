package org.genvisis.cnv.manage;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.genvisis.CLI;
import org.genvisis.cnv.Launch;
import org.genvisis.cnv.analysis.AnalysisFormats;
import org.genvisis.cnv.analysis.Mosaicism;
import org.genvisis.cnv.analysis.pca.PCA;
import org.genvisis.cnv.analysis.pca.PCAPrep;
import org.genvisis.cnv.filesys.ABLookup;
import org.genvisis.cnv.filesys.ABLookup.ABSource;
import org.genvisis.cnv.filesys.Centroids;
import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.Pedigree;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.gui.GenvisisWorkflowGUI;
import org.genvisis.cnv.hmm.CNVCaller;
import org.genvisis.cnv.hmm.CNVCaller.PFB_MANAGEMENT_TYPE;
import org.genvisis.cnv.manage.Resources.Resource;
import org.genvisis.cnv.prop.Property;
import org.genvisis.cnv.qc.GcAdjustor;
import org.genvisis.cnv.qc.LrrSd;
import org.genvisis.cnv.qc.MarkerBlastQC;
import org.genvisis.cnv.qc.MarkerMetrics;
import org.genvisis.cnv.qc.SampleQC;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.Array;
import org.genvisis.common.Elision;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.ext;
import org.genvisis.gwas.Ancestry;
import org.genvisis.gwas.PlinkMendelianChecker;
import org.genvisis.gwas.Qc;

import com.google.common.base.Joiner;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;

public class GenvisisWorkflow {

	private static final String PROJ_PROP_UPDATE_STR = " cnv.filesys.Project proj=";
	private static final String SAMPLE_STEP_REQ_MSG = "[Parse Sample Files] step must have been run already or must be selected and valid.";
	Project proj;
	Logger log;
	private final Launch launch;

	static final STEP S1I_CREATE_MKR_POS = new STEP("Create Marker Positions (if not already exists)",
																									"",
																									new String[][] {{"An Illumina SNP_map file."}},
																									new RequirementInputType[][] {{RequirementInputType.FILE}}) {

		@Override
		public void setNecessaryPreRunProperties(	Project proj, Map<STEP, List<String>> variables) {
			// not needed for step
		}

		@Override
		public void run(Project proj, Map<STEP, List<String>> variables) {
			proj.getLog().report("Generating marker positions file");
			String filename = variables.get(this).get(0);
			org.genvisis.cnv.manage.Markers.generateMarkerPositions(proj, filename);
		}

		@Override
		public boolean[][] checkRequirements(Project proj, Map<STEP, Boolean> stepSelections, Map<STEP, List<String>> variables) {
			return new boolean[][] {{Files.exists(variables.get(this).get(0))}};
		}

		@Override
		public Object[] getRequirementDefaults(Project proj) {
			String filename = proj.getLocationOfSNP_Map(false);
			if (filename == null) {
				filename = "";
			}
			return new Object[] {filename};
		}

		@Override
		public boolean checkIfOutputExists(Project proj, Map<STEP, List<String>> variables) {
			boolean mkrPosFile = Files.exists(proj.MARKER_POSITION_FILENAME.getValue(false, false));
			// boolean mkrSetFile =
			// Files.exists(proj.MARKERSET_FILENAME.getValue(false,
			// false));
			return mkrPosFile/* && mkrSetFile */;
		};

		@Override
		public String getCommandLine(Project proj, Map<STEP, List<String>> variables) {
			String projFile = proj.getPropertyFilename();
			String filename = variables.get(this).get(0);
			return Files.getRunString() + " cnv.manage.Markers proj=" + projFile + " snps=" + filename;
		}

	};

	static final STEP S2I_PARSE_SAMPLES = new STEP(	"Parse Illumina Sample Files", "",
																									new String[][] {{	"[Create Marker Positions] step must be selected and valid.",
																																		"Parsed markerPositions file must already exist."},
																																	{"Number of Threads to Use"}},
																									new RequirementInputType[][] {{	RequirementInputType.NONE,
																																									RequirementInputType.FILE},
																																								{RequirementInputType.NUMBER}},
																									S1I_CREATE_MKR_POS) {

		@Override
		public void setNecessaryPreRunProperties(	Project proj,
																							Map<STEP, List<String>> variables) {
			String projFile = proj.MARKER_POSITION_FILENAME.getValue(false, false);
			String mkrFile = variables.get(this).get(0);
			mkrFile = ext.verifyDirFormat(mkrFile);
			mkrFile = mkrFile.substring(0, mkrFile.length() - 1);
			if (!mkrFile.equals(projFile)) {
				proj.MARKER_POSITION_FILENAME.setValue(mkrFile);
			}
			int numThreads = checkIntArgOrNeg1(variables.get(this).get(1));
			if (numThreads <= 0) {
				numThreads = proj.NUM_THREADS.getValue();
			}
			if (numThreads != proj.NUM_THREADS.getValue()) {
				proj.NUM_THREADS.setValue(numThreads);
			}
		}

		@Override
		public void run(Project proj, Map<STEP, List<String>> variables) {
			int numThreads = proj.NUM_THREADS.getValue();
			proj.getLog().report("Parsing sample files");
			int retCode = org.genvisis.cnv.manage.SourceFileParser.createFiles(proj, numThreads);
			switch (retCode) {
				case 0:
					setFailed();
					failReasons.add("Operation failure, please check log for more information.");
					break;
				case 6:
					failReasons.add("ABLookup required but wasn't found.");
					break;
				case 1:
				default:
					break;
			}
		}

		@Override
		public boolean[][] checkRequirements(	Project proj, Map<STEP, Boolean> stepSelections,
																					Map<STEP, List<String>> variables) {
			int numThreads = checkIntArgOrNeg1(variables.get(this).get(1));
			return new boolean[][] {{stepSelections.get(S1I_CREATE_MKR_POS)
																	&& S1I_CREATE_MKR_POS.hasRequirements(proj, stepSelections,
																																			variables),
																Files.exists(ext.verifyDirFormat(variables.get(this).get(0))),},
															{numThreads > 0}};
		}

		@Override
		public Object[] getRequirementDefaults(Project proj) {
			return new Object[] {	proj.MARKER_POSITION_FILENAME.getValue(false, false),
														proj.NUM_THREADS.getValue()};
		}

		@Override
		public boolean checkIfOutputExists(Project proj, Map<STEP, List<String>> variables) {
			String sampleDirectory = proj.SAMPLE_DIRECTORY.getValue(false, false);
			boolean mkrSetFile = Files.exists(proj.MARKERSET_FILENAME.getValue(false, false));
			return mkrSetFile	&& Files.exists(sampleDirectory)
							&& Files.list(sampleDirectory, Sample.SAMPLE_FILE_EXTENSION, false).length > 0
							&& proj.getSampleList() != null && proj.getSampleList().getSamples().length > 0;
		}

		@Override
		public String getCommandLine(Project proj, Map<STEP, List<String>> variables) {
			String projPropFile = proj.getPropertyFilename();
			StringBuilder kvCmd =
													new StringBuilder(Files.getRunString())	.append(PROJ_PROP_UPDATE_STR)
																																	.append(projPropFile);
			StringBuilder kvPairs = new StringBuilder();
			String projFile = proj.MARKER_POSITION_FILENAME.getValue(false, false);
			String mkrFile = variables.get(this).get(0);
			mkrFile = ext.verifyDirFormat(mkrFile);
			mkrFile = mkrFile.substring(0, mkrFile.length() - 1);
			if (!mkrFile.equals(projFile)) {
				kvPairs.append(" MARKER_POSITION_FILENAME=").append(mkrFile);
			}
			int numThreads = checkIntArgOrNeg1(variables.get(this).get(1));
			if (numThreads <= 0) {
				numThreads = proj.NUM_THREADS.getValue();
			}
			if (numThreads != proj.NUM_THREADS.getValue()) {
				kvPairs.append(" ").append(proj.NUM_THREADS.getName()).append("=").append(numThreads);
			}
			StringBuilder command = new StringBuilder();
			if (kvPairs.length() != 0) {
				command.append(kvCmd).append(kvPairs).append("\n");
			}
			command	.append(Files.getRunString()).append(" cnv.manage.SourceFileParser proj=")
							.append(projPropFile).append(" ").append(PSF.Ext.NUM_THREADS_COMMAND).append(numThreads);
			return command.toString();
		}

	};

	static final STEP S2A_PARSE_SAMPLES =
																			new STEP(	"Parse Sample Files", "",
																								new String[][] {{"markerPositions file must already exist."},
																																{"Number of Threads to Use"}},
																								new RequirementInputType[][] {{RequirementInputType.FILE},
																																							{RequirementInputType.NUMBER}}) {

																				@Override
																				public void setNecessaryPreRunProperties(	Project proj,
																																									Map<STEP, List<String>> variables) {
																					String projFile =
																													proj.MARKER_POSITION_FILENAME.getValue(	false,
																																																	false);
																					String mkrFile = variables.get(this).get(0);
																					mkrFile = ext.verifyDirFormat(mkrFile);
																					mkrFile = mkrFile.substring(0, mkrFile.length() - 1);
																					if (!mkrFile.equals(projFile)) {
																						proj.MARKER_POSITION_FILENAME.setValue(mkrFile);
																					}
																					int numThreads = checkIntArgOrNeg1(variables.get(this).get(1));
																					if (numThreads <= 0) {
																						numThreads = proj.NUM_THREADS.getValue();
																					}
																					if (numThreads != proj.NUM_THREADS.getValue()) {
																						proj.NUM_THREADS.setValue(numThreads);
																					}
																				}

																				@Override
																				public void run(Project proj, Map<STEP, List<String>> variables) {
																					int numThreads = checkIntArgOrNeg1(variables.get(this).get(1));
																					if (numThreads <= 0) {
																						numThreads = proj.NUM_THREADS.getValue();
																					}
																					proj.getLog().report("Parsing sample files");
																					int retCode = org.genvisis.cnv.manage.SourceFileParser.createFiles(proj, numThreads);
																					switch (retCode) {
																						case 0:
																							setFailed();
																							failReasons.add("Operation failure, please check log for more information.");
																							break;
																						case 6:
																							failReasons.add("ABLookup required but wasn't found.");
																							break;
																						case 1:
																						default:
																							break;
																					}
																				}

																				@Override
																				public boolean[][] checkRequirements(	Project proj,
																																							Map<STEP, Boolean> stepSelections,
																																							Map<STEP, List<String>> variables) {

																					int numThreads = checkIntArgOrNeg1(variables.get(this).get(1));
																					return new boolean[][] {{Files.exists(ext.verifyDirFormat(variables.get(this).get(0))),},
																																	{numThreads > 0}};
																				}

																				@Override
																				public Object[] getRequirementDefaults(Project proj) {
																					return new Object[] {	proj.MARKER_POSITION_FILENAME.getValue(false,
																																																			false),
																																proj.NUM_THREADS.getValue()};
																				}

																				@Override
																				public boolean checkIfOutputExists(	Project proj,
																																						Map<STEP, List<String>> variables) {
																					String sampleDirectory = proj.SAMPLE_DIRECTORY.getValue(false,
																																																	false);
																					boolean mkrSetFile = Files.exists(proj.MARKERSET_FILENAME.getValue(	false,
																																																							false));
																					return mkrSetFile	&& Files.exists(sampleDirectory)
																									&& Files.list(sampleDirectory,
																																Sample.SAMPLE_FILE_EXTENSION, false).length > 0
																									&& proj.getSampleList() != null
																									&& proj.getSampleList().getSamples().length > 0;
																				}

																				@Override
																				public String getCommandLine(	Project proj,
																																			Map<STEP, List<String>> variables) {
																					String projPropFile = proj.getPropertyFilename();
																					StringBuilder kvCmd =
																															new StringBuilder(Files.getRunString())	.append(PROJ_PROP_UPDATE_STR)
																																																			.append(projPropFile);
																					StringBuilder kvPairs = new StringBuilder();
																					String projFile =
																													proj.MARKER_POSITION_FILENAME.getValue(	false,
																																																	false);
																					String mkrFile = variables.get(this).get(0);
																					mkrFile = ext.verifyDirFormat(mkrFile);
																					mkrFile = mkrFile.substring(0, mkrFile.length() - 1);
																					if (!mkrFile.equals(projFile)) {
																						kvPairs	.append(" MARKER_POSITION_FILENAME=")
																										.append(mkrFile);
																					}

																					int numThreads = checkIntArgOrNeg1(variables.get(this).get(1));
																					if (numThreads <= 0) {
																						numThreads = proj.NUM_THREADS.getValue();
																					}
																					if (numThreads != proj.NUM_THREADS.getValue()) {
																						kvPairs.append(" ").append(proj.NUM_THREADS.getName()).append("=").append(numThreads);
																					}
																					StringBuilder command = new StringBuilder();
																					if (kvPairs.length() != 0) {
																						command.append(kvCmd).append(kvPairs).append("\n");
																					}
																					command	.append(Files.getRunString())
																									.append(" cnv.manage.SourceFileParser proj=")
																									.append(projPropFile).append(" ")
																									.append(PSF.Ext.NUM_THREADS_COMMAND).append(numThreads);
																					return command.toString();
																				}

																			};

	static final STEP S3_CREATE_SAMPLEDATA = new STEP("Create SampleData.txt File", "",
																										new String[][] {{"[Parse Sample Files] step must have been run already or must be selected and valid"},
																																		{	"Create a minimal SampleData.txt file from sample files",
																																			"Either a Pedigree.dat file, or any file with a header containing all of the following elements (in any order):  \""
																																																																	+ Array.toStr(MitoPipeline.PED_INPUT,
																																																																							", ")
																																																																+ "\"",
																																			"A Sample_Map.csv file, with at least two columns having headers \""
																																																																					+ MitoPipeline.SAMPLEMAP_INPUT[1]
																																																																				+ "\" and \""
																																																																				+ MitoPipeline.SAMPLEMAP_INPUT[2]
																																																																				+ "\""}},
																										new RequirementInputType[][] {{RequirementInputType.NONE},
																																									{	RequirementInputType.BOOL,
																																										RequirementInputType.FILE,
																																										RequirementInputType.FILE}},
																										S2I_PARSE_SAMPLES, S2A_PARSE_SAMPLES) {

		@Override
		public void setNecessaryPreRunProperties(	Project proj,
																							Map<STEP, List<String>> variables) {
			// Nothing to do
		}

		@Override
		public void run(Project proj, Map<STEP, List<String>> variables) {
			Boolean minimal = Boolean.parseBoolean(variables.get(this).get(0));
			String pedFile = minimal ? null : variables.get(this).get(1);
			String sampleMapCsv = minimal ? null : variables.get(this).get(2);

			proj.getLog().report("Creating SampleData.txt");
			try {
				int retStat = SampleData.createSampleData(pedFile, sampleMapCsv, proj);
				if (retStat == -1) {
					setFailed();
					failReasons.add("SampleData already exists - please delete and try again.");
					return;
				}
			} catch (Elision e) {
				String msg = e.getMessage();
				setFailed();
				failReasons.add(msg);
				return;
			}
		}

		@Override
		public boolean[][] checkRequirements(	Project proj, Map<STEP, Boolean> stepSelections,
																					Map<STEP, List<String>> variables) {
			String sampDir = proj.SAMPLE_DIRECTORY.getValue();
			STEP parseStep = stepSelections.containsKey(S2I_PARSE_SAMPLES)	? S2I_PARSE_SAMPLES
																																			: S2A_PARSE_SAMPLES;
			boolean checkStepParseSamples = stepSelections.get(parseStep)
																			&& parseStep.hasRequirements(proj, stepSelections, variables);
			return new boolean[][] {{checkStepParseSamples
																|| (Files.exists(sampDir) && Files.list(sampDir, ".sampRAF",
																																				proj.JAR_STATUS.getValue()).length > 0)},
															{	Boolean.parseBoolean(variables.get(this).get(0)),
																Files.exists(variables.get(this).get(1)),
																Files.exists(variables.get(this).get(2))}};
		}

		@Override
		public Object[] getRequirementDefaults(Project proj) {
			String pedPreset = Files.exists(proj.PEDIGREE_FILENAME.getValue())	? proj.PEDIGREE_FILENAME.getValue()
																																					: "";
			String sampMapPreset = "".equals(pedPreset) ? getLocationOfSampleMap(proj) : null; // check
																																													// for
																																													// SampleMap
																																													// only
																																													// if
																																													// we
																																													// haven't
																																													// found
																																													// a
																																													// pedigree
																																													// file
			if (sampMapPreset == null) {
				sampMapPreset = "";
			}
			return new Object[] {true, pedPreset, sampMapPreset};
		}

		private String getLocationOfSampleMap(Project proj) {
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


		@Override
		public boolean checkIfOutputExists(Project proj, Map<STEP, List<String>> variables) {
			return Files.exists(proj.SAMPLE_DATA_FILENAME.getValue(false, false));
		}

		@Override
		public String getCommandLine(Project proj, Map<STEP, List<String>> variables) {
			String projPropFile = proj.getPropertyFilename();
			Boolean minimal = Boolean.parseBoolean(variables.get(this).get(0));
			String pedFile = minimal ? "" : variables.get(this).get(1);
			String sampleMapCsv = minimal ? "" : variables.get(this).get(2);
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

	static final STEP S4_TRANSPOSE_TO_MDF = new STEP(	"Transpose Data into Marker-Dominant Files", "",
																										new String[][] {{SAMPLE_STEP_REQ_MSG}},
																										new RequirementInputType[][] {{RequirementInputType.NONE}},
																										S2I_PARSE_SAMPLES, S2A_PARSE_SAMPLES) {

		@Override
		public void setNecessaryPreRunProperties(	Project proj,
																							Map<STEP, List<String>> variables) {
			// Nothing to do here
		}

		@Override
		public void run(Project proj, Map<STEP, List<String>> variables) {
			proj.getLog().report("Transposing data");
			TransposeData.transposeData(proj, 2000000000, false); // compact if no LRR was provided
		}

		@Override
		public boolean[][] checkRequirements(	Project proj, Map<STEP, Boolean> stepSelections,
																					Map<STEP, List<String>> variables) {
			String sampDir = proj.SAMPLE_DIRECTORY.getValue();
			STEP parseStep = stepSelections.containsKey(S2I_PARSE_SAMPLES)	? S2I_PARSE_SAMPLES
																																			: S2A_PARSE_SAMPLES;
			boolean checkStepParseSamples = stepSelections.get(parseStep)
																			&& parseStep.hasRequirements(proj, stepSelections, variables);
			return new boolean[][] {{checkStepParseSamples
																|| (Files.exists(sampDir)
																		&& Files.list(sampDir, ".sampRAF",
																									proj.JAR_STATUS.getValue()).length > 0),}};
		}

		@Override
		public Object[] getRequirementDefaults(Project proj) {
			return new Object[] {};
		}

		@Override
		public boolean checkIfOutputExists(Project proj, Map<STEP, List<String>> variables) {
			return Files.list(proj.MARKER_DATA_DIRECTORY.getValue(false, false),
												MarkerData.MARKER_DATA_FILE_EXTENSION, false).length > 0;
		}

		@Override
		public String getCommandLine(Project proj, Map<STEP, List<String>> variables) {
			String projPropFile = proj.getPropertyFilename();
			StringBuilder cmd = new StringBuilder();
			return cmd.append(Files.getRunString()).append(" cnv.manage.TransposeData -transpose proj="
																											+ projPropFile + " max=" + 2000000000)
								.toString();
		}

	};


	static final STEP S5_COMPUTE_GCMODEL = new STEP("Compute GCMODEL File", "",
																									new String[][] {{"A GC Base file must exist."},
																																	{"GCModel output file must be specified."}},
																									new RequirementInputType[][] {{RequirementInputType.FILE},
																																								{RequirementInputType.FILE}}) {

		@Override
		public void setNecessaryPreRunProperties(	Project proj,
																							Map<STEP, List<String>> variables) {
			String setGCOutputFile = proj.GC_MODEL_FILENAME.getValue();
			String gcOutputFile = variables.get(this).get(1);
			if (!ext.verifyDirFormat(setGCOutputFile).equals(gcOutputFile)) {
				proj.GC_MODEL_FILENAME.setValue(gcOutputFile);
			}
		}

		@Override
		public void run(Project proj, Map<STEP, List<String>> variables) {
			String gcBaseFile = variables.get(this).get(0);
			String gcOutputFile = variables.get(this).get(1);
			org.genvisis.cnv.analysis.PennCNV.gcModel(proj, gcBaseFile, gcOutputFile, 100);
		}

		@Override
		public boolean[][] checkRequirements(	Project proj, Map<STEP, Boolean> stepSelections,
																					Map<STEP, List<String>> variables) {
			String gcBaseFile = variables.get(this).get(0);
			String gcOutputFile = variables.get(this).get(1);
			return new boolean[][] {{Files.exists(gcBaseFile)}, {!Files.exists(gcOutputFile)}};
		}

		@Override
		public Object[] getRequirementDefaults(Project proj) {
			Resource gc5Base = Resources.genome(proj.GENOME_BUILD_VERSION.getValue(), proj.getLog())
																	.getModelBase();
			if (gc5Base.isAvailable(true)) {
				return new Object[] {gc5Base.get(), proj.GC_MODEL_FILENAME.getValue()};
			}
			return null;
		}

		@Override
		public boolean checkIfOutputExists(Project proj, Map<STEP, List<String>> variables) {
			String gcOutputFile = variables.get(this).get(1);
			boolean gcExists = Files.exists(gcOutputFile);
			return gcExists;
		}

		@Override
		public String getCommandLine(Project proj, Map<STEP, List<String>> variables) {
			String kvCmd = "";

			String setGCOutputFile = proj.GC_MODEL_FILENAME.getValue();
			String gcOutputFile = variables.get(this).get(1);
			if (!ext.verifyDirFormat(setGCOutputFile).equals(gcOutputFile)) {
				kvCmd += " GC_MODEL_FILENAME=" + gcOutputFile;
			}

			String projPropFile = proj.getPropertyFilename();
			StringBuilder cmd = new StringBuilder();
			if (kvCmd.length() > 0) {
				cmd	.append(Files.getRunString()).append(PROJ_PROP_UPDATE_STR + projPropFile)
						.append(kvCmd).append("\n");
			}
			String gcBaseFile = variables.get(this).get(0);
			return cmd.append(Files.getRunString())
								.append(" cnv.analysis.PennCNV proj="	+ proj.getPropertyFilename() + " log="
												+ proj.getLog().getFilename() + " gc5base=" + gcBaseFile)
								.toString();
		}

	};
	static final STEP S6_SAMPLE_QC = new STEP("Run Sample QC Metrics", "",
																						new String[][] {{SAMPLE_STEP_REQ_MSG},
																														{"Number of threads to use."},},
																						new RequirementInputType[][] {{RequirementInputType.NONE},
																																					{RequirementInputType.NUMBER}},
																						S2I_PARSE_SAMPLES, S2A_PARSE_SAMPLES) {

		@Override
		public void setNecessaryPreRunProperties(	Project proj,
																							Map<STEP, List<String>> variables) {
			int numThreads = checkIntArgOrNeg1(variables.get(this).get(0));
			if (numThreads <= 0) {
				numThreads = proj.NUM_THREADS.getValue();
			}
			if (numThreads != proj.NUM_THREADS.getValue()) {
				proj.NUM_THREADS.setValue(numThreads);
			}
		}

		@Override
		public void run(Project proj, Map<STEP, List<String>> variables) {
			proj.getLog().report("Running LrrSd");
			int numThreads = proj.NUM_THREADS.getValue();
			LrrSd.init(proj, null, null, numThreads, false);
		}

		@Override
		public boolean[][] checkRequirements(	Project proj, Map<STEP, Boolean> stepSelections,
																					Map<STEP, List<String>> variables) {

			int numThreads = checkIntArgOrNeg1(variables.get(this).get(0));
			
			String sampDir = proj.SAMPLE_DIRECTORY.getValue();
			STEP parseStep = stepSelections.containsKey(S2I_PARSE_SAMPLES)	? S2I_PARSE_SAMPLES
																																			: S2A_PARSE_SAMPLES;
			boolean checkStepParseSamples = stepSelections.get(parseStep)
																			&& parseStep.hasRequirements(proj, stepSelections, variables);
			return new boolean[][] {{checkStepParseSamples
																|| (Files.exists(sampDir) && Files.list(sampDir, ".sampRAF",
																																				proj.JAR_STATUS.getValue()).length > 0),},
															{numThreads > 0}};
		}

		@Override
		public Object[] getRequirementDefaults(Project proj) {
			return new Object[] {proj.NUM_THREADS.getValue()};
		}

		@Override
		public boolean checkIfOutputExists(Project proj, Map<STEP, List<String>> variables) {
			return Files.exists(proj.SAMPLE_QC_FILENAME.getValue(false, false));
		}

		@Override
		public String getCommandLine(Project proj, Map<STEP, List<String>> variables) {
			int numThreads = checkIntArgOrNeg1(variables.get(this).get(0));
			if (numThreads <= 0) {
				numThreads = proj.NUM_THREADS.getValue();
			}
			String projPropFile = proj.getPropertyFilename();
			StringBuilder cmd = new StringBuilder();
			cmd	.append(Files.getRunString()).append(" cnv.qc.LrrSd").append(" proj=").append(projPropFile)
					.append(" ").append(PSF.Ext.NUM_THREADS_COMMAND).append(numThreads).append(" projectMarkers=TRUE");
			return cmd.toString();
		}
	};

	static final STEP S7_MARKER_QC = new STEP("Run Marker QC Metrics", "",
																						new String[][] {{SAMPLE_STEP_REQ_MSG},
																														{	"Export all markers in project.",
																															"A targetMarkers files listing the markers to QC."},
																														{"Number of threads to use."}},
																						new RequirementInputType[][] {{RequirementInputType.NONE},
																																					{	RequirementInputType.NONE,
																																						RequirementInputType.FILE},
																																					{RequirementInputType.NUMBER}},
																						S2I_PARSE_SAMPLES) {

		@Override
		public void setNecessaryPreRunProperties(	Project proj,
																							Map<STEP, List<String>> variables) {
			// Nothing to do here
		}

		@Override
		public void run(Project proj, Map<STEP, List<String>> variables) {
			String tgtFile = variables.get(this).get(0);
			tgtFile = "".equals(tgtFile) ? null : tgtFile;
			int numThreads = checkIntArgOrNeg1(variables.get(this).get(1));
			if (numThreads <= 0) {
				numThreads = proj.NUM_THREADS.getValue();
			}
			MarkerMetrics.fullQC(proj, null, tgtFile, true, numThreads);
		}

		@Override
		public Object[] getRequirementDefaults(Project proj) {
			return new Object[] {	""/* proj.TARGET_MARKERS_FILENAMES.getValue()[0] */,
														proj.NUM_THREADS.getValue()};
		}

		@Override
		public boolean[][] checkRequirements(	Project proj, Map<STEP, Boolean> stepSelections,
																					Map<STEP, List<String>> variables) {
			String mkrPosFile = proj.MARKER_POSITION_FILENAME.getValue();
			String sampList = proj.SAMPLELIST_FILENAME.getValue();
			String tgtFile = variables.get(this).get(0);

			int numThreads = checkIntArgOrNeg1(variables.get(this).get(1));
			boolean step12 = Files.exists(mkrPosFile);
			boolean step2 = stepSelections.get(S2I_PARSE_SAMPLES)
											&& S2I_PARSE_SAMPLES.hasRequirements(proj, stepSelections, variables);
			boolean step22 = Files.exists(sampList);
			boolean step3 = Files.exists(tgtFile);
			return new boolean[][] {{step2 || step12 && step22}, {true, step3}, {numThreads > 0}};
		}

		@Override
		public boolean checkIfOutputExists(Project proj, Map<STEP, List<String>> variables) {
			String markerMetricsFile = proj.MARKER_METRICS_FILENAME.getValue();
			return Files.exists(markerMetricsFile);
		}

		@Override
		public String getCommandLine(Project proj, Map<STEP, List<String>> variables) {
			String tgtFile = variables.get(this).get(0);
			tgtFile = "".equals(tgtFile) ? null : tgtFile;
			int numThreads = checkIntArgOrNeg1(variables.get(this).get(1));
			if (numThreads <= 0) {
				numThreads = proj.NUM_THREADS.getValue();
			}
			return Files.getRunString()	+ " cnv.qc.MarkerMetrics -fullQC" + " proj="
							+ proj.getPropertyFilename() + (tgtFile == null ? "" : " markers=" + tgtFile)
							+ " " + PSF.Ext.NUM_THREADS_COMMAND + numThreads;
		}
	};

	static final STEP S8_SEX_CHECKS = new STEP(	"Run Sex Checks", "",
																							new String[][] {{SAMPLE_STEP_REQ_MSG},
																															{"[Create SampleData.txt File] step must have been run already or must be selected and valid."},
																															{"Add Estimated Sex to Sample Data"},
																															{	"Use only X and Y chromosome R values to identify sex discriminating markers",
																																"List of markers that do not cross hybridize",
																																"BLAST annotation VCF to generate list of markers that do not cross hybridize from"}},
																							new RequirementInputType[][] {{RequirementInputType.NONE},
																																						{RequirementInputType.NONE},
																																						{RequirementInputType.BOOL},
																																						{	RequirementInputType.BOOL,
																																							RequirementInputType.FILE,
																																							RequirementInputType.FILE}},
																							S2I_PARSE_SAMPLES, S2A_PARSE_SAMPLES,
																							S3_CREATE_SAMPLEDATA) {

		@Override
		public void setNecessaryPreRunProperties(	Project proj,
																							Map<STEP, List<String>> variables) {
			// Nothing to do here
		}

		@Override
		public void run(Project proj, Map<STEP, List<String>> variables) {
			proj.getLog().report("Running SexCheck");
			boolean addToSampleData = Boolean.parseBoolean(variables.get(this).get(0));
			String discriminatingMarkersFile;
			if (Boolean.parseBoolean(variables.get(this).get(1))) {
				discriminatingMarkersFile = null;
			} else {
				discriminatingMarkersFile = variables.get(this).get(2);
				if (!Files.exists(discriminatingMarkersFile)) {
					MarkerBlastQC.getOneHitWonders(	proj, variables.get(this).get(3),
																					discriminatingMarkersFile, 0.8, proj.getLog());
				}
			}
			org.genvisis.cnv.qc.SexChecks.sexCheck(proj, addToSampleData, discriminatingMarkersFile);
		}

		@Override
		public boolean[][] checkRequirements(	Project proj, Map<STEP, Boolean> stepSelections,
																					Map<STEP, List<String>> variables) {
			String sampDir = proj.SAMPLE_DIRECTORY.getValue();
			String sampDataFile = proj.SAMPLE_DATA_FILENAME.getValue(false, false);
			STEP parseStep = stepSelections.containsKey(S2I_PARSE_SAMPLES)	? S2I_PARSE_SAMPLES
																																			: S2A_PARSE_SAMPLES;
			boolean checkStepParseSamples = stepSelections.get(parseStep)
																			&& parseStep.hasRequirements(proj, stepSelections, variables);
			boolean useRValues = Boolean.parseBoolean(variables.get(this).get(1));
			String discriminatingMarkersFile = variables.get(this).get(2);
			String blastVCFFile = variables.get(this).get(3);
			return new boolean[][] {{checkStepParseSamples
																|| (Files.exists(sampDir) && Files.list(sampDir, ".sampRAF",
																																				proj.JAR_STATUS.getValue()).length > 0)},
															{stepSelections.get(S3_CREATE_SAMPLEDATA)
																&& S3_CREATE_SAMPLEDATA.hasRequirements(proj, stepSelections,
																																				variables)
																|| Files.exists(sampDataFile)},
															{true},
															{	useRValues, !useRValues && Files.exists(discriminatingMarkersFile),
																!useRValues && !Files.exists(discriminatingMarkersFile) && Files.exists(blastVCFFile)}};
		}

		@Override
		public Object[] getRequirementDefaults(Project proj) {
			return new Object[] {	true, false,
														MarkerBlastQC.defaultOneHitWondersFilename(proj.BLAST_ANNOTATION_FILENAME.getValue()),
														proj.BLAST_ANNOTATION_FILENAME.getValue()};
		}

		@Override
		public boolean checkIfOutputExists(Project proj, Map<STEP, List<String>> variables) {
			return Files.exists(proj.SEXCHECK_RESULTS_FILENAME.getValue());
		}

		@Override
		public String getCommandLine(Project proj, Map<STEP, List<String>> variables) {
			String projPropFile = proj.getPropertyFilename();
			StringBuilder cmd = new StringBuilder();
			boolean addToSampleData = Boolean.parseBoolean(variables.get(this).get(0));
			String discriminatingMarkersFile;
			if (Boolean.parseBoolean(variables.get(this).get(1))) {
				discriminatingMarkersFile = null;
			} else {
				discriminatingMarkersFile = variables.get(this).get(2);
				if (!Files.exists(discriminatingMarkersFile)) {
					cmd	.append(Files.getRunString()).append(" cnv.qc.MarkerBlastQC proj="	+ projPropFile
																										+ " blastVCF=" + variables.get(this).get(3))
							.append("\n");
				}
			}
			return cmd.append(Files.getRunString())
								.append(" cnv.qc.SexChecks -check proj=" + projPropFile).toString()
								+ (discriminatingMarkersFile == null	? ""
																										: " useMarkers=" + discriminatingMarkersFile)
							+ (addToSampleData ? "" : " -skipSampleData");
		}

	};

	static final STEP S9_GENERATE_ABLOOKUP = new STEP("Generate AB Lookup File", "",
																										new String[][] {{SAMPLE_STEP_REQ_MSG}},
																										new RequirementInputType[][] {{RequirementInputType.NONE}},
																										S2I_PARSE_SAMPLES, S2A_PARSE_SAMPLES) {

		@Override
		public void setNecessaryPreRunProperties(	Project proj,
																							Map<STEP, List<String>> variables) {
			// Nothing to do here
		}

		@Override
		public void run(Project proj, Map<STEP, List<String>> variables) {
			String filename = proj.PROJECT_DIRECTORY.getValue()	+ ext.addToRoot(ABLookup.DEFAULT_AB_FILE, "_parsed");
			ABLookup.parseABLookup(proj, ABSource.VCF, filename);

			if (ABLookup.fillInMissingAlleles(proj, filename, proj.getLocationOfSNP_Map(true), false)) {
				ABLookup.applyABLookupToFullSampleFiles(proj, filename);
			} else {
				setFailed();
			}
		}

		@Override
		public Object[] getRequirementDefaults(Project proj) {
			return new Object[] {};
		}

		@Override
		public boolean[][] checkRequirements(	Project proj, Map<STEP, Boolean> stepSelections,
																					Map<STEP, List<String>> variables) {
			String mkrPosFile = proj.MARKER_POSITION_FILENAME.getValue();
			String sampDir = proj.SAMPLE_DIRECTORY.getValue();
			boolean checkStepParseSamples = stepSelections.get(S2I_PARSE_SAMPLES)
																			&& S2I_PARSE_SAMPLES.hasRequirements(	proj, stepSelections,
																																						variables);
			return new boolean[][] {{checkStepParseSamples
																|| (Files.exists(mkrPosFile)	&& Files.exists(sampDir)
																		&& Files.list(sampDir, ".sampRAF",
																									proj.JAR_STATUS.getValue()).length > 0)},};
		}

		@Override
		public boolean checkIfOutputExists(Project proj, Map<STEP, List<String>> variables) {
			return Files.exists(proj.AB_LOOKUP_FILENAME.getValue(false, false));
		}

		@Override
		public String getCommandLine(Project proj, Map<STEP, List<String>> variables) {
			String filename = proj.PROJECT_DIRECTORY.getValue()	+ ext.addToRoot(ABLookup.DEFAULT_AB_FILE, "_parsed");
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
			
			List<String> commandProp = Lists.newArrayList(ImmutableList.of(Files.getRunString(), Project.class.getName(), CLI.formCmdLineArg(CLI.ARG_PROJ, projFile)));
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
	};

	static final STEP S10_RUN_PLINK =
																	new STEP(	"Create PLINK Files", "",
																						new String[][] {{SAMPLE_STEP_REQ_MSG},
																														// {},
																														{	"A pedigree.dat file is must exist.",
																															"Create a minimal pedigree.dat file [will pull information from SexChecks step results]."}},
																						new RequirementInputType[][] {{RequirementInputType.NONE},
																																					{	RequirementInputType.FILE,
																																						RequirementInputType.BOOL}},
																						S2I_PARSE_SAMPLES, S2A_PARSE_SAMPLES) {

																		@Override
																		public void setNecessaryPreRunProperties(	Project proj,
																																							Map<STEP, List<String>> variables) {
																			if (!Boolean.parseBoolean(variables.get(this).get(1))) {
																				String projPedFile = proj.PEDIGREE_FILENAME.getValue(	false,
																																															false);
																				String pedFile = variables.get(this).get(0);
																				if (!pedFile.equals(projPedFile)) {
																					proj.PEDIGREE_FILENAME.setValue(pedFile);
																				}
																			}
																		}

																		@Override
																		public void run(Project proj,
																										Map<STEP, List<String>> variables) {
																			if (Boolean.parseBoolean(variables.get(this).get(1))) {
																				proj.getLog().report("Creating Pedigree File");
																				Pedigree.build(proj, null, null, false);
																			}
																			if (!Files.exists(proj.PEDIGREE_FILENAME.getValue())) {
																				setFailed();
																				failReasons.add("Creation of Pedigree file in [Create/Run PLINK Files] step failed.");
																				return;
																			}

																			proj.getLog().report("Running PLINK");

																			boolean create = PlinkData.saveGenvisisToPlinkBedSet(	proj,
																																														"plink/plink",
																																														null,
																																														null,
																																														-1,
																																														true);
																			if (!create) {
																				setFailed();
																				failReasons.add("Creation of initial PLINK files failed.");
																			}
																			proj.PLINK_DIR_FILEROOTS.addValue(proj.PROJECT_DIRECTORY.getValue()
																																				+ "plink/plink");
																		}

																		@Override
																		public boolean[][] checkRequirements(	Project proj,
																																					Map<STEP, Boolean> stepSelections,
																																					Map<STEP, List<String>> variables) {
																			String sampDir = proj.SAMPLE_DIRECTORY.getValue();
																			String pedFile = variables.get(this).get(0);
																			STEP parseStep = stepSelections.containsKey(S2I_PARSE_SAMPLES)	? S2I_PARSE_SAMPLES
																																																			: S2A_PARSE_SAMPLES;
																			boolean checkStepParseSamples = stepSelections.get(parseStep)
																																			&& parseStep.hasRequirements(	proj,
																																																		stepSelections,
																																																		variables);
																			return new boolean[][] {{checkStepParseSamples
																																|| (Files.exists(sampDir)
																																		&& Files.list(sampDir,
																																									".sampRAF",
																																									proj.JAR_STATUS.getValue()).length > 0)},
																															{	Files.exists(pedFile),
																																Boolean.parseBoolean(variables.get(this)
																																															.get(1))}};
																		}

																		@Override
																		public Object[] getRequirementDefaults(Project proj) {
																			return new Object[] {	proj.PEDIGREE_FILENAME.getValue(false,
																																														false),
																														false};
																		}

																		@Override
																		public boolean checkIfOutputExists(	Project proj,
																																				Map<STEP, List<String>> variables) {
																			// String fileCheck1 =
																			// proj.PROJECT_DIRECTORY.getValue()+"gwas.map";
																			String fileCheck2 = proj.PROJECT_DIRECTORY.getValue()
																													+ "plink/plink.bed";
																			// String fileCheck3 =
																			// proj.PROJECT_DIRECTORY.getValue()+"genome/";
																			boolean pedCheck = Boolean.parseBoolean(variables.get(this)
																																									.get(1))	? Files.exists(proj.PEDIGREE_FILENAME.getValue())
																																														: true;
																			return /* Files.exists(fileCheck1) && */ Files.exists(fileCheck2)
																							/* && Files.exists(fileCheck3) */ /*
																																								 * && Files.list(
																																								 * fileCheck3,
																																								 * ".bed",
																																								 * false).length > 0
																																								 */ && pedCheck;
																		}

																		@Override
																		public String getCommandLine(	Project proj,
																																	Map<STEP, List<String>> variables) {
																			String kvCmd = "";

																			if (!Boolean.parseBoolean(variables.get(this).get(1))) {
																				String projPedFile = proj.PEDIGREE_FILENAME.getValue(	false,
																																															false);
																				String pedFile = variables.get(this).get(0);
																				if (!pedFile.equals(projPedFile)) {
																					kvCmd += " PEDIGREE_FILENAME=" + pedFile;
																				}
																			}

																			String projPropFile = proj.getPropertyFilename();
																			StringBuilder cmd = new StringBuilder();
																			if (kvCmd.length() > 0) {
																				cmd	.append(Files.getRunString())
																						.append(PROJ_PROP_UPDATE_STR)
																						.append(projPropFile).append(kvCmd).append("\n");
																			}
																			if (Boolean.parseBoolean(variables.get(this).get(1))) {
																				cmd	.append(Files.getRunString())
																						.append(" cnv.filesys.Pedigree proj=")
																						.append(projPropFile).append("\n");
																			}
																			cmd	.append(Files.getRunString())
																					.append(" cnv.manage.PlinkData -genvisisToBed plinkdata=plink/plink gcthreshold=-1 proj=")
																					.append(proj.getPropertyFilename());
																			return cmd.toString();
																		}
																	};

	static final STEP S11_GWAS_QC = new STEP(	"Run GWAS QC", "",
																						new String[][] {{"[Create/Run PLINK Files] step must have been run already or must be selected and valid."},
																														{"Keep genome info for unrelateds only?"},
																														{	"Skip ancestry checks",
																															"File with FID/IID pairs of putative white samples",
																															"PLINK root of HapMap founders"}},
																						new RequirementInputType[][] {{RequirementInputType.NONE},
																																					{RequirementInputType.BOOL},
																																					{	RequirementInputType.BOOL,
																																						RequirementInputType.FILE,
																																						RequirementInputType.FILE}},
																						S10_RUN_PLINK) {

		@Override
		public void setNecessaryPreRunProperties(	Project proj,
																							Map<STEP, List<String>> variables) {
			// not needed for step
		}

		@Override
		public void run(Project proj, Map<STEP, List<String>> variables) {
			String dir = getPlinkDir(proj, variables);
			boolean keepUnrelatedsOnly = Boolean.parseBoolean(variables.get(this).get(0));
			boolean skipAncestry = Boolean.parseBoolean(variables.get(this).get(1));
			String putativeWhites = variables.get(this).get(2);
			String hapMapPlinkRoot = variables.get(this).get(3);
			if (hapMapPlinkRoot.lastIndexOf('.') > 0) {
				String extension = hapMapPlinkRoot.substring(hapMapPlinkRoot.lastIndexOf('.'));
				if (".bed".equals(extension) || ".bim".equals(extension) || ".fam".equals(extension)) {
					hapMapPlinkRoot = hapMapPlinkRoot.substring(0, hapMapPlinkRoot.lastIndexOf('.'));
				}
			}
			Qc.fullGamut(dir, null, keepUnrelatedsOnly, proj.getLog());
			if (new File(dir + "genome/plink.genome").exists()) {
				proj.GENOME_CLUSTER_FILENAME.setValue(dir + "genome/plink.genome");
				proj.saveProperties();
			}
			if (!keepUnrelatedsOnly) {
				new PlinkMendelianChecker(proj).run();
			}
			if (!skipAncestry) {
				String ancestryDir = dir + Qc.ANCESTRY_DIR;
				Ancestry.runPipeline(	ancestryDir, putativeWhites, hapMapPlinkRoot, proj,
															new Logger(ancestryDir + "ancestry.log"));
			}
		}

		@Override
		public boolean[][] checkRequirements(	Project proj, Map<STEP, Boolean> stepSelections,
																					Map<STEP, List<String>> variables) {
			String dir = getPlinkDir(proj, variables);
			boolean skipAncestry = Boolean.parseBoolean(variables.get(this).get(1));
			boolean putativeWhitesExists = Files.exists(variables.get(this).get(2));
			String hapMapPlinkRoot = variables.get(this).get(3);
			if (hapMapPlinkRoot.lastIndexOf('.') > 0) {
				String extension = hapMapPlinkRoot.substring(hapMapPlinkRoot.lastIndexOf('.'));
				if (".bed".equals(extension) || ".bim".equals(extension) || ".fam".equals(extension)) {
					hapMapPlinkRoot = hapMapPlinkRoot.substring(0, hapMapPlinkRoot.lastIndexOf('.'));
				}
			}
			boolean hapMapPlinkRootExists = Files.exists(hapMapPlinkRoot + ".bed")
																				&& Files.exists(hapMapPlinkRoot + ".bim")
																			&& Files.exists(hapMapPlinkRoot + ".fam");
			boolean files = Files.exists(dir);
			return new boolean[][] {{(stepSelections.get(S10_RUN_PLINK)
																&& S10_RUN_PLINK.hasRequirements(proj, stepSelections, variables))
																|| files},
															{true}, {	skipAncestry, putativeWhitesExists && hapMapPlinkRootExists,
																				hapMapPlinkRootExists && putativeWhitesExists}};
		}

		@Override
		public Object[] getRequirementDefaults(Project proj) {
			return new Object[] {false, false, "", Ancestry.DEFAULT_HAPMAP_PLINKROOT};
		}

		@Override
		public boolean checkIfOutputExists(Project proj, Map<STEP, List<String>> variables) {
			String dir = getPlinkDir(proj, variables);
			boolean allExist = true;
			boolean skipAncestry = Boolean.parseBoolean(variables.get(this).get(1));
			folders: for (int i = 0; i < org.genvisis.gwas.Qc.FOLDERS_CREATED.length; i++) {
				for (int j = 0; j < org.genvisis.gwas.Qc.FILES_CREATED[i].length; j++) {
					if (!Files.exists(dir	+ org.genvisis.gwas.Qc.FOLDERS_CREATED[i]
														+ org.genvisis.gwas.Qc.FILES_CREATED[i][j])) {
						allExist = false;
						break folders;
					}
				}
			}
			if (!skipAncestry && (!Files.exists(dir + Qc.ANCESTRY_DIR + "freqsByRace.xln")
						|| !Files.exists(dir + Qc.ANCESTRY_DIR + "raceImputations.mds"))) {
				allExist = false;
			}
			return allExist;
		}

		@Override
		public String getCommandLine(Project proj, Map<STEP, List<String>> variables) {
			String dir = getPlinkDir(proj, variables);
			boolean keepUnrelatedsOnly = Boolean.parseBoolean(variables.get(this).get(0));
			boolean skipAncestry = Boolean.parseBoolean(variables.get(this).get(1));
			String putativeWhites = variables.get(this).get(2);
			String hapMapPlinkRoot = variables.get(this).get(3);
			if (hapMapPlinkRoot.lastIndexOf('.') > 0) {
				String extension = hapMapPlinkRoot.substring(hapMapPlinkRoot.lastIndexOf('.'));
				if (".bed".equals(extension) || ".bim".equals(extension) || ".fam".equals(extension)) {
					hapMapPlinkRoot = hapMapPlinkRoot.substring(0, hapMapPlinkRoot.lastIndexOf('.'));
				}
			}

			String command = Files.getRunString()	+ " gwas.Qc dir=" + dir
												+ " keepGenomeInfoForRelatedsOnly=" + keepUnrelatedsOnly;
			command += "\n";
			command += Files.getRunString() + " org.genvisis.cnv.filesys.Project proj=" + proj.getPropertyFilename();
			command += " " + proj.GENOME_CLUSTER_FILENAME.getName() + "=" + dir + "genome/plink.genome";
			if (!keepUnrelatedsOnly) {
				command += "\n";
				command += Files.getRunString() + " org.genvisis.gwas.PlinkMendelianChecker proj=" + proj.getPropertyFilename();
			}
			if (!skipAncestry) {
				String ancestryDir = dir + Qc.ANCESTRY_DIR;
				command += "\n";
				command += Files.getRunString() + " gwas.Ancestry -runPipeline dir=" + ancestryDir;
				command += " putativeWhites=" + putativeWhites;
				command += " proj=" + proj.getPropertyFilename();
				command += " hapMapPlinkRoot=" + hapMapPlinkRoot;
				command += " log=" + ancestryDir + "ancestry.log";
			}
			return command;
		}

		private String getPlinkDir(Project proj, Map<STEP, List<String>> variables) {
			// This directory is also used in S10_ANNOTATE_SAMPLE_DATA, any changes should be reflected
			// there
			String dir = "plink/";// variables.get(this).get(0);
			if (!dir.startsWith("/") && !dir.contains(":")) {
				dir = ext.verifyDirFormat(proj.PROJECT_DIRECTORY.getValue() + dir);
			}
			return dir;
		}

	};

	static final STEP S12_MOSAIC_ARMS = new STEP(	"Create Mosaic Arms File", "",
																								new String[][] {{SAMPLE_STEP_REQ_MSG},
																																{"Number of threads to use."}},
																								new RequirementInputType[][] {{RequirementInputType.NONE},
																																							{RequirementInputType.NUMBER}},
																								S2I_PARSE_SAMPLES) {

		@Override
		public void setNecessaryPreRunProperties(	Project proj,
																							Map<STEP, List<String>> variables) {

			int numThreads = checkIntArgOrNeg1(variables.get(this).get(0));
			if (numThreads <= 0) {
				numThreads = proj.NUM_THREADS.getValue();
			}
			if (numThreads != proj.NUM_THREADS.getValue()) {
				proj.NUM_THREADS.setValue(numThreads);
			}
		}

		@Override
		public void run(Project proj, Map<STEP, List<String>> variables) {
			Mosaicism.findOutliers(proj);
		}

		@Override
		public boolean[][] checkRequirements(	Project proj, Map<STEP, Boolean> stepSelections,
																					Map<STEP, List<String>> variables) {
			boolean checkStepParseSamples = stepSelections.get(S2I_PARSE_SAMPLES)
																			&& S2I_PARSE_SAMPLES.hasRequirements(	proj, stepSelections,
																																						variables);
			String mkrPosFile = proj.MARKERSET_FILENAME.getValue();
			boolean step11 = Files.exists(mkrPosFile);

			int numThreads = checkIntArgOrNeg1(variables.get(this).get(0));
			return new boolean[][] {{checkStepParseSamples || step11}, {numThreads > 0}};
		}

		@Override
		public boolean checkIfOutputExists(Project proj, Map<STEP, List<String>> variables) {
			boolean outputCheck = Files.exists(proj.RESULTS_DIRECTORY.getValue(false, false)
																					+ "Mosaicism.xln");
			return outputCheck;
		}

		@Override
		public Object[] getRequirementDefaults(Project proj) {
			return new Object[] {proj.NUM_THREADS.getValue()};
		}

		@Override
		public String getCommandLine(Project proj, Map<STEP, List<String>> variables) {
			String kvCmd = "";


			int numThreads = checkIntArgOrNeg1(variables.get(this).get(0));
			if (numThreads <= 0) {
				numThreads = proj.NUM_THREADS.getValue();
			}
			if (numThreads != proj.NUM_THREADS.getValue()) {
				kvCmd += " " + proj.NUM_THREADS.getName() + "=" + numThreads;
			}

			String projPropFile = proj.getPropertyFilename();
			StringBuilder cmd = new StringBuilder();
			if (kvCmd.length() > 0) {
				cmd	.append(Files.getRunString()).append(PROJ_PROP_UPDATE_STR + projPropFile)
						.append(kvCmd).append("\n");
			}
			return cmd.append(Files.getRunString())
								.append(" cnv.analysis.Mosaicism proj=" + proj.getPropertyFilename()).toString();
		}

	};
	static final STEP S13_ANNOTATE_SAMPLE_DATA = new STEP("Annotate Sample Data File", "",
																												new String[][] {{"["	+ S6_SAMPLE_QC.stepName
																																					+ "] step must have been run already or must be selected and valid."},
																																				{"["
																																						+ S3_CREATE_SAMPLEDATA.stepName
																																					+ "] step must have been run already or must be selected and valid."},
																																				{	"Skip identifying duplicates?",
																																					"[" + S11_GWAS_QC.stepName + "] step must have been run already or must be selected and valid."},
																																				{	"Do not use GC corrected LRR SD?",
																																					"GC Corrected LRR SD must exist in Sample QC File"},
																																				{"LRR SD Threshold"},
																																				{"Callrate Threshold"},
																																				{"Number of Quantiles to Generate"},
																																				{"Replace FID and IID with data from Pedigree?"}},
																												new RequirementInputType[][] {{RequirementInputType.NONE},
																																											{RequirementInputType.NONE},
																																											{	RequirementInputType.BOOL,
																																												RequirementInputType.NONE},
																																											{	RequirementInputType.BOOL,
																																												RequirementInputType.NONE},
																																											{RequirementInputType.NUMBER},
																																											{RequirementInputType.NUMBER},
																																											{RequirementInputType.NUMBER},
																																											{RequirementInputType.BOOL}},
																												S3_CREATE_SAMPLEDATA, S6_SAMPLE_QC) {

		@Override
		public void setNecessaryPreRunProperties(	Project proj,
																							Map<STEP, List<String>> variables) {
			double projLrrSdThreshold = proj.LRRSD_CUTOFF.getValue();
			double lrrSdThreshold = Double.parseDouble(variables.get(this).get(2));
			double projCallrateThreshold = proj.SAMPLE_CALLRATE_THRESHOLD.getValue();
			double callrateThreshold = Double.parseDouble(variables.get(this).get(3));

			if (projLrrSdThreshold != lrrSdThreshold) {
				proj.LRRSD_CUTOFF.setValue(lrrSdThreshold);
			}
			if (projCallrateThreshold != callrateThreshold) {
				proj.SAMPLE_CALLRATE_THRESHOLD.setValue(callrateThreshold);
			}
		}

		@Override
		public void run(Project proj, Map<STEP, List<String>> variables) {
			boolean checkDuplicates = !Boolean.parseBoolean(variables.get(this).get(0));
			String duplicatesSetFile = null;
			if (checkDuplicates) {
				String dir = "plink/";
				duplicatesSetFile = proj.PROJECT_DIRECTORY.getValue()	+ dir
														+ "/quality_control/genome/plink.genome_duplicatesSet.dat";
			}
			boolean gcCorrectedLrrSd = !Boolean.parseBoolean(variables.get(this).get(1));
			int numQ = Integer.parseInt(variables.get(this).get(4));
			boolean correctFidIids = Boolean.parseBoolean(variables.get(this).get(5));
			SampleQC.parseAndAddToSampleData(	proj, numQ, 0, false, gcCorrectedLrrSd, duplicatesSetFile,
																				correctFidIids);
		}

		@Override
		public boolean[][] checkRequirements(	Project proj, Map<STEP, Boolean> stepSelections,
																					Map<STEP, List<String>> variables) {
			boolean checkStepSampleQC = stepSelections.get(S6_SAMPLE_QC)
																	&& S6_SAMPLE_QC.hasRequirements(proj, stepSelections, variables);
			String sampleQCFile = proj.SAMPLE_QC_FILENAME.getValue();
			boolean sampleQCFileExists = Files.exists(sampleQCFile);
			String sampleDataFile = proj.SAMPLE_DATA_FILENAME.getValue();
			boolean checkDuplicates = !Boolean.parseBoolean(variables.get(this).get(0));
			String dir = "plink/";
			String duplicatesSetFile = proj.PROJECT_DIRECTORY.getValue()	+ dir
																	+ "/quality_control/genome/plink.genome_duplicatesSet.dat";
			boolean gcCorrectedLrrSd = !Boolean.parseBoolean(variables.get(this).get(1));
			boolean gcCorrectedLrrSdExists = false;
			if (sampleQCFileExists
					&& ext.indexOfStr("LRR_SD_Post_Correction",
														Files.getHeaderOfFile(sampleQCFile, proj.getLog())) != -1) {
				gcCorrectedLrrSdExists = true;
			}
			double lrrSdThreshold = checkDoubleArgOrNeg1(variables.get(this).get(2));
			double callrateThreshold = checkDoubleArgOrNeg1(variables.get(this).get(3));
			int numQ = checkIntArgOrNeg1(variables.get(this).get(4));
			return new boolean[][] {{checkStepSampleQC || sampleQCFileExists},
															{stepSelections.get(S3_CREATE_SAMPLEDATA)
																&& S3_CREATE_SAMPLEDATA.hasRequirements(proj, stepSelections,
																																				variables)
																|| Files.exists(sampleDataFile)},
															{	!checkDuplicates,
																stepSelections.get(S11_GWAS_QC)
																									&& S11_GWAS_QC.hasRequirements(	proj,
																																									stepSelections, variables)
																									|| Files.exists(duplicatesSetFile)},
															{!gcCorrectedLrrSd, gcCorrectedLrrSdExists},
															{lrrSdThreshold > proj.LRRSD_CUTOFF.getMinValue()
																&& lrrSdThreshold < proj.LRRSD_CUTOFF.getMaxValue()},
															{callrateThreshold > proj.SAMPLE_CALLRATE_THRESHOLD.getMinValue()
																&& callrateThreshold < proj.SAMPLE_CALLRATE_THRESHOLD.getMaxValue()},
															{numQ > 0}, {true}};
		}

		@Override
		public Object[] getRequirementDefaults(Project proj) {
			return new Object[] {	false, false, proj.LRRSD_CUTOFF.getValue(),
														proj.SAMPLE_CALLRATE_THRESHOLD.getValue(), 10, false};
		}

		@Override
		public boolean checkIfOutputExists(Project proj, Map<STEP, List<String>> variables) {
			String sampleDataFile = proj.SAMPLE_DATA_FILENAME.getValue();
			if (!Files.exists(sampleDataFile)) {
				return false;
			}
			boolean checkDuplicates = !Boolean.parseBoolean(variables.get(this).get(0));
			String[] header = Files.getHeaderOfFile(sampleDataFile, proj.getLog());
			if (checkDuplicates && ext.indexOfStr("DuplicateId", header, false, true) == -1) {
				return false;
			}
			String[] reqHdr = {
			                   
			};
			if (ext.indexOfStr("Class=Exclude", header, false, true) == -1) {
				return false;
			}
			if (ext.indexOfStr("ExcludeNote", header, false, true) == -1) {
				return false;
			}
			if (ext.indexOfStr("Use", header, false, true) == -1) {
				return false;
			}
			if (ext.indexOfStr("UseNote", header, false, true) == -1) {
				return false;
			}
			if (ext.indexOfStr("Use_cnv", header, false, true) == -1) {
				return false;
			}
			if (ext.indexOfStr("Use_cnvNote", header, false, true) == -1) {
				return false;
			}
			return true;
		}

		@Override
		public String getCommandLine(Project proj, Map<STEP, List<String>> variables) {

			double projLrrSdThreshold = proj.LRRSD_CUTOFF.getValue();
			double lrrSdThreshold = Double.parseDouble(variables.get(this).get(2));
			double projCallrateThreshold = proj.SAMPLE_CALLRATE_THRESHOLD.getValue();
			double callrateThreshold = Double.parseDouble(variables.get(this).get(3));

			String projPropFile = proj.getPropertyFilename();

			boolean checkDuplicates = !Boolean.parseBoolean(variables.get(this).get(0));
			String duplicatesSetFile = null;
			if (checkDuplicates) {
				String dir = "plink/";
				duplicatesSetFile = proj.PROJECT_DIRECTORY.getValue()	+ dir
														+ "/quality_control/genome/plink.genome_duplicatesSet.dat";
			}
			boolean gcCorrectedLrrSd = !Boolean.parseBoolean(variables.get(this).get(1));
			int numQ = Integer.parseInt(variables.get(this).get(4));
			boolean correctFidIids = Boolean.parseBoolean(variables.get(this).get(5));

			String kvCmd = "";

			if (projLrrSdThreshold != lrrSdThreshold) {
				kvCmd += " LRRSD_CUTOFF=" + lrrSdThreshold;
			}
			if (projCallrateThreshold != callrateThreshold) {
				kvCmd += " SAMPLE_CALLRATE_THRESHOLD=" + callrateThreshold;
			}

			StringBuilder cmd = new StringBuilder();
			if (kvCmd.length() > 0) {
				cmd	.append(Files.getRunString()).append(PROJ_PROP_UPDATE_STR + projPropFile)
						.append(kvCmd).append("\n");
			}
			cmd	.append(Files.getRunString())
					.append(" cnv.qc.SampleQC proj="	+ projPropFile + " numQ=" + numQ + " justQuantiles=false"
									+ " gcCorrectedLrrSd=" + gcCorrectedLrrSd + " duplicatesSetFile="
									+ duplicatesSetFile + " correctFidIids=" + correctFidIids);
			return cmd.toString();
		}

	};

	// FIXME http://genvisis.org/MitoPipeline/#illumina_marker_lists has illumina markers.. this should be linked to, or
	// these steps split or something...
	static final STEP S14_CREATE_PCS = new STEP("Create Principal Components File and Mitochondrial Copy-Number Estimates File",
																							"",
																							new String[][] {{"[Transpose Data into Marker-Dominant Files] step must have been run already or must be selected and valid."},
																															{"MedianMarkers file must exist."},
																															{"LRR SD threshold to filter samples."},
																															{"Call rate threshold to filter markers."},
																															{"Compute PCs with samples passing QC only?"},
																															{"Should impute mean value for NaN?"},
																															{"Should recompute Log-R ratio for PC markers?"},
																															{"Should recompute Log-R ratio for median markers?"},
																															{"Homozygous only?"},
																															{"Regression distance for the GC adjustment"},
																															{"Number of threads to use"},
																															{"A file listing a subset of samples (DNA ID) to use for determining optimal PC selection, typically a list of unrelated and single race samples. If a list is not provided, only samples passing sample qc thresholds will be used."},
																															{"An external beta file to optimize PC selection."},},
																							new RequirementInputType[][] {{RequirementInputType.NONE},
																																						{RequirementInputType.FILE},
																																						{RequirementInputType.NUMBER},
																																						{RequirementInputType.NUMBER},
																																						{RequirementInputType.BOOL},
																																						{RequirementInputType.BOOL},
																																						{RequirementInputType.BOOL},
																																						{RequirementInputType.BOOL},
																																						{RequirementInputType.BOOL},
																																						{RequirementInputType.NUMBER},
																																						{RequirementInputType.NUMBER},
																																						{RequirementInputType.FILE},
																																						{RequirementInputType.FILE},},
																							S4_TRANSPOSE_TO_MDF) {

		@Override
		public void setNecessaryPreRunProperties(	Project proj,
																							Map<STEP, List<String>> variables) {
			double sampleLRRSdFilter = Double.parseDouble(variables.get(this).get(1));
			if (sampleLRRSdFilter < 0) {
				switch (proj.ARRAY_TYPE.getValue()) {
					case AFFY_GW6:
					case AFFY_GW6_CN:
						proj.LRRSD_CUTOFF.setValue(0.35);
						proj.getLog()
								.reportTimeInfo("Setting "	+ proj.LRRSD_CUTOFF.getName()
																+ " to default 0.35 for array " + proj.ARRAY_TYPE.getValue());
						break;
					case ILLUMINA:
						proj.LRRSD_CUTOFF.setValue(0.30);
						proj.getLog()
								.reportTimeInfo("Setting "	+ proj.LRRSD_CUTOFF.getName()
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
		public void run(Project proj, Map<STEP, List<String>> variables) {
			String medianMarkers = variables.get(this).get(0);
			double markerCallRateFilter = Double.parseDouble(variables.get(this).get(2));
			boolean gcCorrect = Boolean.parseBoolean(variables.get(this).get(3));
			boolean imputeMeanForNaN = Boolean.parseBoolean(variables.get(this).get(4));
			boolean recomputeLRR_PCs = Boolean.parseBoolean(variables.get(this).get(5));
			boolean recomputeLRR_Median = Boolean.parseBoolean(variables.get(this).get(6));
			boolean homozygousOnly = Boolean.parseBoolean(variables.get(this).get(7));
			int bpGcModel = GcAdjustor.GcModel.DEFAULT_GC_MODEL_BIN_FASTA;// Integer.parseInt(variables.get(this).get(8));
			int regressionDistance = Integer.parseInt(variables.get(this).get(8));
			int numComponents = MitoPipeline.DEFAULT_NUM_COMPONENTS;// Integer.parseInt(variables.get(this).get(10));
			int numThreads = Integer.parseInt(variables.get(this).get(9));
			String outputBase = MitoPipeline.FILE_BASE;

			String betaOptFile = variables.get(this).get(10);
			String betaFile = variables.get(this).get(11);

			boolean markerQC = true;
			double[] pvalOpt = MitoPipeline.DEFAULT_PVAL_OPTS;
			String pedFile = null;
			String useFile = null;
			boolean sampLrr = true;
			boolean plot = false;
			int retCode = PCAPrep.prepPCA(proj, numThreads, outputBase, markerQC, markerCallRateFilter,
																		useFile, proj.getSampleList(), proj.getLog());
			if (retCode == 42) { // TODO remove magic number
				MitoPipeline.estimateMtDNACN(	proj, numThreads, medianMarkers, numComponents, outputBase,
																			homozygousOnly, markerCallRateFilter, betaOptFile, pedFile,
																			recomputeLRR_PCs, recomputeLRR_Median, sampLrr,
																			imputeMeanForNaN, gcCorrect, bpGcModel, regressionDistance,
																			proj.GENOME_BUILD_VERSION.getValue(), pvalOpt,
																			betaFile, plot, false, proj.getLog());
			} else {
				setFailed();
			}
		}

		@Override
		public boolean[][] checkRequirements(	Project proj, Map<STEP, Boolean> stepSelections,
																					Map<STEP, List<String>> variables) {
			String markerDir = proj.MARKER_DATA_DIRECTORY.getValue();
			String medianMkrs = variables.get(this).get(0).trim();
			String betaOptFile = variables.get(this).get(10).trim();
			String betaFile = variables.get(this).get(11).trim();
			double lrrThresh = checkDoubleArgOrNeg1(variables.get(this).get(1));
			double callrate = checkDoubleArgOrNeg1(variables.get(this).get(2));
			int regressionDistance = checkIntArgOrNeg1(variables.get(this).get(8));
			int numThreads = checkIntArgOrNeg1(variables.get(this).get(9));
			return new boolean[][] {{(stepSelections.get(S4_TRANSPOSE_TO_MDF)
																&& S4_TRANSPOSE_TO_MDF.hasRequirements(	proj, stepSelections,
																																				variables))
																|| Files.exists(markerDir)},
															{Files.exists(medianMkrs)}, {lrrThresh > 0}, {callrate > 0}, {true}, // TRUE
																																																		// or
																																																		// FALSE
																																																		// are
																																																		// both
																																																		// valid
																																																		// selections
															{true}, {true}, {true}, {true}, {regressionDistance > 0},
															{numThreads > 0},
															{"".equals(betaOptFile) || Files.exists(betaOptFile)},
															{"".equals(betaFile) || Files.exists(betaFile)},};
		}

		@Override
		public Object[] getRequirementDefaults(Project proj) {
			return new Object[] {	"", proj.LRRSD_CUTOFF.getValue(),
														MitoPipeline.DEFAULT_MKR_CALLRATE_FILTER, "true", "true", "true",
														"true", "true", GcAdjustor.DEFAULT_REGRESSION_DISTANCE[0],
														proj.NUM_THREADS.getValue(), "", "",};
		}

		@Override
		public boolean checkIfOutputExists(Project proj, Map<STEP, List<String>> variables) {
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

		@Override
		public String getCommandLine(Project proj, Map<STEP, List<String>> variables) {
			String medianMarkers = variables.get(this).get(0);
			double lrrSD = Double.parseDouble(variables.get(this).get(1));
			double markerCallRateFilter = Double.parseDouble(variables.get(this).get(2));
			boolean gcCorrect = Boolean.parseBoolean(variables.get(this).get(3));
			boolean imputeMeanForNaN = Boolean.parseBoolean(variables.get(this).get(4));
			boolean recomputeLRR_PCs = Boolean.parseBoolean(variables.get(this).get(5));
			boolean recomputeLRR_Median = Boolean.parseBoolean(variables.get(this).get(6));
			boolean homozygousOnly = Boolean.parseBoolean(variables.get(this).get(7));
			int bpGcModel = GcAdjustor.GcModel.DEFAULT_GC_MODEL_BIN_FASTA;
			int regressionDistance = Integer.parseInt(variables.get(this).get(8));
			int numComponents = MitoPipeline.DEFAULT_NUM_COMPONENTS;
			int numThreads = Integer.parseInt(variables.get(this).get(9));
			String outputBase = proj.PROJECT_DIRECTORY.getValue() + MitoPipeline.FILE_BASE;

			String betaOptFile = variables.get(this).get(10).trim();
			String betaFile = variables.get(this).get(11).trim();
			boolean sampLrr = true;


			String projPropFile = proj.getPropertyFilename();
			StringBuilder cmd = new StringBuilder();
			cmd	.append(Files.getRunString()).append(" org.genvisis.cnv.manage.MitoPipeline")
					.append(" proj=").append(projPropFile).append(" mitochondrialMarkers=")
					.append(medianMarkers).append(" numComponents=").append(numComponents)
					.append(" imputeMeanForNaN=").append(imputeMeanForNaN).append(" recomputeLRR_PCs=")
					.append(recomputeLRR_PCs).append(" recomputeLRR_Median=").append(recomputeLRR_Median)
					.append(" gcCorrect=").append(gcCorrect).append(" bpGcModel=").append(bpGcModel)
					.append(" LRRSD=").append(lrrSD).append(" markerCallRate=").append(markerCallRateFilter)
					.append(" regressionDistance=").append(regressionDistance).append(" sampLRR=")
					.append(sampLrr).append(" ").append(PSF.Ext.NUM_THREADS_COMMAND).append(numThreads).append(" log=")
					.append(proj.getLog().getFilename()).append(" output=").append(outputBase);
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

	};



	static final STEP S15_COMPUTE_PFB = new STEP(	"Compute Population BAF files", "",
																								new String[][] {{	SAMPLE_STEP_REQ_MSG,
																																	"A Sample subset file must exist."},
																																{"PFB (population BAF) output file must be specified."}},
																								new RequirementInputType[][] {{	RequirementInputType.NONE,
																																								RequirementInputType.FILE},
																																							{RequirementInputType.FILE}},
																								S2I_PARSE_SAMPLES, S2A_PARSE_SAMPLES) {

		@Override
		public void setNecessaryPreRunProperties(	Project proj,
																							Map<STEP, List<String>> variables) {
			String setSubSampFile = proj.SAMPLE_SUBSET_FILENAME.getValue();
			String subSampFile = variables.get(this).get(0);
			String setPFBFile = proj.CUSTOM_PFB_FILENAME.getValue();
			String pfbOutputFile = variables.get(this).get(1);

			if (!ext.verifyDirFormat(setSubSampFile).equals(subSampFile)) {
				proj.SAMPLE_SUBSET_FILENAME.setValue(subSampFile);
			}
			if (!ext.verifyDirFormat(setPFBFile).equals(pfbOutputFile)) {
				proj.CUSTOM_PFB_FILENAME.setValue(pfbOutputFile);
			}
		}

		@Override
		public void run(Project proj, Map<STEP, List<String>> variables) {
			org.genvisis.cnv.analysis.PennCNV.populationBAF(proj);
		}

		@Override
		public Object[] getRequirementDefaults(Project proj) {
			return new Object[] {	proj.SAMPLE_SUBSET_FILENAME.getValue(),
														Files.exists(proj.SAMPLE_SUBSET_FILENAME.getValue())	? ext.rootOf(proj.SAMPLE_SUBSET_FILENAME.getValue())
																																										+ ".pfb"
																																									: proj.CUSTOM_PFB_FILENAME.getValue()};
		}

		@Override
		public boolean[][] checkRequirements(	Project proj, Map<STEP, Boolean> stepSelections,
																					Map<STEP, List<String>> variables) {
			String sampListFile = proj.SAMPLELIST_FILENAME.getValue();
			String subSampFile = variables.get(this).get(0);
			String pfbOutputFile = variables.get(this).get(1);

			STEP parseStep = stepSelections.containsKey(S2I_PARSE_SAMPLES)	? S2I_PARSE_SAMPLES
																																			: S2A_PARSE_SAMPLES;
			boolean checkStepParseSamples = stepSelections.get(parseStep)
																			&& parseStep.hasRequirements(proj, stepSelections, variables);
			boolean step12 = Files.exists(sampListFile);
			boolean step13 = Files.exists(subSampFile);
			boolean step21 = !Files.exists(pfbOutputFile);

			return new boolean[][] {{checkStepParseSamples || step12, step13}, {step21}};
		}

		@Override
		public boolean checkIfOutputExists(Project proj, Map<STEP, List<String>> variables) {
			String subSampFile = variables.get(this).get(0);
			String pfbOutputFile = variables.get(this).get(1);
			boolean pfbExists = Files.exists(pfbOutputFile)
													|| Files.exists(ext.rootOf(subSampFile) + ".pfb");
			return pfbExists;
		}

		@Override
		public String getCommandLine(Project proj, Map<STEP, List<String>> variables) {
			String kvCmd = "";

			String setSubSampFile = proj.SAMPLE_SUBSET_FILENAME.getValue();
			String subSampFile = variables.get(this).get(0);
			String setPFBFile = proj.CUSTOM_PFB_FILENAME.getValue();
			String pfbOutputFile = variables.get(this).get(1);

			if (!ext.verifyDirFormat(setSubSampFile).equals(subSampFile)) {
				kvCmd += " SAMPLE_SUBSET_FILENAME=" + subSampFile;
			}
			if (!ext.verifyDirFormat(setPFBFile).equals(pfbOutputFile)) {
				kvCmd += " CUSTOM_PFB_FILENAME=" + pfbOutputFile;
			}

			String projPropFile = proj.getPropertyFilename();
			StringBuilder cmd = new StringBuilder();
			if (kvCmd.length() > 0) {
				cmd	.append(Files.getRunString()).append(PROJ_PROP_UPDATE_STR + projPropFile)
						.append(kvCmd).append("\n");
			}
			return cmd.append(Files.getRunString())
								.append(" cnv.analysis.PennCNV -pfb proj="	+ proj.getPropertyFilename() + " log="
												+ proj.getLog().getFilename())
								.toString();
		}

	};
	static final STEP S16_SEX_CENTROIDS_PFB_GCMODEL =
																									new STEP(	"Create Sex-Specific Centroids; Filter PFB and GCMODEL Files",
																														"",
																														new String[][] {{"["
																																								+ S5_COMPUTE_GCMODEL.stepName
																																							+ "] must be selected and valid.",
																																							"Full GC Model File."},
																																						{"Number of threads to use."},},
																														new RequirementInputType[][] {{	RequirementInputType.NONE,
																																														RequirementInputType.FILE},
																																													{RequirementInputType.NUMBER}},
																														S5_COMPUTE_GCMODEL) {
																										@Override
																										public void setNecessaryPreRunProperties(	Project proj,
																																															Map<STEP, List<String>> variables) {

																											int numThreads = checkIntArgOrNeg1(variables.get(this).get(1));
																											if (numThreads <= 0) {
																												numThreads = proj.NUM_THREADS.getValue();
																											}
																											if (numThreads != proj.NUM_THREADS.getValue()) {
																												proj.NUM_THREADS.setValue(numThreads);
																											}
																										}

																										@Override
																										public void run(Project proj,
																																		Map<STEP, List<String>> variables) {
																											String malePFB, femalePFB, centFilePathM,
																													centFilePathF, newGCFile;
																											// pennData =
																											// proj.getProperty(proj.PENNCNV_DATA_DIRECTORY);
																											// sexDir = pennData + "sexSpecific/";
																											String outputDir = proj.DATA_DIRECTORY.getValue();
																											newGCFile = outputDir + "sexSpecific.gcModel";
																											malePFB = outputDir + "males.pfb";
																											femalePFB = outputDir + "females.pfb";
																											centFilePathM = outputDir
																																			+ "sexSpecific_Male.cent";
																											centFilePathF = outputDir
																																			+ "sexSpecific_Female.cent";

																											int numThreads = checkIntArgOrNeg1(variables.get(this).get(1));
																											if (numThreads <= 0) {
																												numThreads = proj.NUM_THREADS.getValue();
																											}
																											String gcModelFile = variables.get(this)
																																										.get(0);
																											Centroids.computeSexSpecificCentroids(proj,
																																														AnalysisFormats.getChromosomalMarkersOnly(proj),
																																														new String[] {malePFB,
																																																					femalePFB},
																																														new String[] {centFilePathM,
																																																					centFilePathF},
																																														true,
																																														numThreads);

																											AnalysisFormats.filterSexSpecificGCModel(	proj,
																																																gcModelFile,
																																																newGCFile);
																										}

																										@Override
																										public boolean[][] checkRequirements(	Project proj,
																																													Map<STEP, Boolean> stepSelections,
																																													Map<STEP, List<String>> variables) {
																											boolean checkStepGCModel = stepSelections.get(S5_COMPUTE_GCMODEL)
																																									&& S5_COMPUTE_GCMODEL.hasRequirements(proj,
																																																												stepSelections, variables);

																											int numThreads = checkIntArgOrNeg1(variables.get(this).get(1));
																											String gcModelFile = variables.get(this)
																																										.get(0);
																											return new boolean[][] {{	checkStepGCModel,
																																								Files.exists(gcModelFile)},
																																							{numThreads > 0}};
																										}

																										@Override
																										public Object[] getRequirementDefaults(Project proj) {
																											return new Object[] {	proj.GC_MODEL_FILENAME.getValue(),
																																						proj.NUM_THREADS.getValue()};
																										}

																										@Override
																										public boolean checkIfOutputExists(	Project proj,
																																												Map<STEP, List<String>> variables) {
																											String malePFB, femalePFB, centFilePathM,
																													centFilePathF, newGCFile;
																											// pennData =
																											// proj.getProperty(proj.PENNCNV_DATA_DIRECTORY);
																											// sexDir = pennData + "sexSpecific/";
																											String outputDir = proj.DATA_DIRECTORY.getValue();
																											malePFB = outputDir + "males.pfb";
																											femalePFB = outputDir + "females.pfb";
																											centFilePathM = outputDir
																																			+ "sexSpecific_Male.cent";
																											centFilePathF = outputDir
																																			+ "sexSpecific_Female.cent";
																											newGCFile = outputDir + "sexSpecific.gcModel";
																											return Files.exists(malePFB)
																																&& Files.exists(femalePFB)
																															&& Files.exists(centFilePathM)
																															&& Files.exists(centFilePathF)
																															&& Files.exists(newGCFile);
																										}

																										@Override
																										public String getCommandLine(	Project proj,
																																									Map<STEP, List<String>> variables) {

																											int numThreads = checkIntArgOrNeg1(variables.get(this).get(1));
																											if (numThreads <= 0) {
																												numThreads = proj.NUM_THREADS.getValue();
																											}
																											String mainCmd = Files.getRunString()
																																					+ " cnv.filesys.Centroids proj="
																																				+ proj.getPropertyFilename()
																																				+ " -sexSpecific " + PSF.Ext.NUM_THREADS_COMMAND + numThreads;
																											String gcModelFile = variables.get(this)
																																										.get(0);
																											String gcCmd = Files.getRunString()
																																				+ " cnv.analysis.AnalysisFormats proj="
																																			+ proj.getPropertyFilename()
																																			+ " gcmodel=" + gcModelFile;
																											return mainCmd + "\n" + gcCmd;
																										}

																									};

	static final STEP S17_CNV_CALLING = new STEP(	"Call CNVs", "",
																								new String[][] {{"Hidden Markov Model File Must Exist"},
																																{	"[Compute Population BAF File] step must be selected and valid",
																																	"PFB File Must Exist"},
																																{	"[Compute GCMODEL File] step must be selected and valid",
																																	"GCMODEL File Must Exist"},
																																{"Number of threads To use."},
																																{"Output filename."}},
																								new RequirementInputType[][] {{RequirementInputType.FILE},
																																							{	RequirementInputType.NONE,
																																								RequirementInputType.FILE},
																																							{	RequirementInputType.NONE,
																																								RequirementInputType.FILE},
																																							{RequirementInputType.NUMBER},
																																							{RequirementInputType.FILE},},
																								S5_COMPUTE_GCMODEL, S15_COMPUTE_PFB) {
		@Override
		public void setNecessaryPreRunProperties(	Project proj,
																							Map<STEP, List<String>> variables) {
			String hmm_P = proj.HMM_FILENAME.getValue();
			String hmm_G = variables.get(this).get(0);
			if (!hmm_P.equals(hmm_G)) {
				proj.HMM_FILENAME.setValue(hmm_G);
			}
			String pfb_P = proj.CUSTOM_PFB_FILENAME.getValue();
			String pfb_G = variables.get(this).get(1);
			if (!pfb_P.equals(pfb_G)) {
				proj.CUSTOM_PFB_FILENAME.setValue(pfb_G);
			}
			String gcm_P = proj.GC_MODEL_FILENAME.getValue();
			String gcm_G = variables.get(this).get(2);
			if (!gcm_P.equals(gcm_G)) {
				proj.GC_MODEL_FILENAME.setValue(gcm_G);
			}
			int numThreads = checkIntArgOrNeg1(variables.get(this).get(3));
			if (numThreads <= 0) {
				numThreads = proj.NUM_THREADS.getValue();
			}
			if (numThreads != proj.NUM_THREADS.getValue()) {
				proj.NUM_THREADS.setValue(numThreads);
			}
		}

		@Override
		public void run(Project proj, Map<STEP, List<String>> variables) {
			int numThreads = checkIntArgOrNeg1(variables.get(this).get(3));
			if (numThreads <= 0) {
				numThreads = proj.NUM_THREADS.getValue();
			}
			if (numThreads != proj.NUM_THREADS.getValue()) {
				proj.NUM_THREADS.setValue(numThreads);
			}
			String output = variables.get(this).get(4); // gets PROJ_DIR prepended, so NOT ABSOLUTE
			(new File(ext.parseDirectoryOfFile(proj.PROJECT_DIRECTORY.getValue() + output))).mkdirs();
			CNVCaller.callAutosomalCNVs(proj, output, proj.getSamples(), null, null,
																	CNVCaller.DEFUALT_MIN_SITES, CNVCaller.DEFUALT_MIN_CONF,
																	PFB_MANAGEMENT_TYPE.PENNCNV_DEFAULT, numThreads, 1);// TODO, sex
																																											// specific
																																											// centroids,etc
			proj.CNV_FILENAMES.addValue(proj.PROJECT_DIRECTORY.getValue() + output);
			proj.saveProperties(new Property[] {proj.CNV_FILENAMES});
		}

		@Override
		public boolean[][] checkRequirements(	Project proj, Map<STEP, Boolean> stepSelections,
																					Map<STEP, List<String>> variables) {
			boolean checkHMM = Files.exists(variables.get(this).get(0));
			boolean checkPFB1 = stepSelections.get(S15_COMPUTE_PFB).booleanValue();
			boolean checkPFB2 = Files.exists(variables.get(this).get(1));
			boolean checkGC1 = stepSelections.get(S5_COMPUTE_GCMODEL).booleanValue();
			boolean checkGC2 = Files.exists(variables.get(this).get(2));
			int numThreads = checkIntArgOrNeg1(variables.get(this).get(3));
			return new boolean[][] {{checkHMM}, {checkPFB1, checkPFB2}, {checkGC1, checkGC2},
															{numThreads > 0}, {!Files.exists(variables.get(this).get(4))},};
		}

		@Override
		public Object[] getRequirementDefaults(Project proj) {
			return new Object[] {	proj.HMM_FILENAME.getValue(), proj.CUSTOM_PFB_FILENAME.getValue(),
														proj.GC_MODEL_FILENAME.getValue(), proj.NUM_THREADS.getValue(),
														"cnvs/genvisis.cnv"};
		}

		@Override
		public boolean checkIfOutputExists(Project proj, Map<STEP, List<String>> variables) {
			String output = variables.get(this).get(4);
			return Files.exists(proj.PROJECT_DIRECTORY.getValue() + output);
		}

		@Override
		public String getCommandLine(Project proj, Map<STEP, List<String>> variables) {
			String kvCmd = "";

			String hmm_P = proj.HMM_FILENAME.getValue();
			String hmm_G = variables.get(this).get(0);
			if (!hmm_P.equals(hmm_G)) {
				kvCmd += " HMM_FILENAME=" + hmm_G;
			}
			String pfb_P = proj.CUSTOM_PFB_FILENAME.getValue();
			String pfb_G = variables.get(this).get(1);
			if (!pfb_P.equals(pfb_G)) {
				kvCmd += " CUSTOM_PFB_FILENAME=" + pfb_G;
			}
			String gcm_P = proj.GC_MODEL_FILENAME.getValue();
			String gcm_G = variables.get(this).get(2);
			if (!gcm_P.equals(gcm_G)) {
				kvCmd += " GC_MODEL_FILENAME=" + gcm_G;
			}

			int numThreads = checkIntArgOrNeg1(variables.get(this).get(3));
			if (numThreads <= 0) {
				numThreads = proj.NUM_THREADS.getValue();
			}
			if (numThreads != proj.NUM_THREADS.getValue()) {
				kvCmd += " " + proj.NUM_THREADS.getName() + "=" + numThreads;
			}
			String projPropFile = proj.getPropertyFilename();
			StringBuilder cmd = new StringBuilder();
			if (kvCmd.length() > 0) {
				cmd	.append(Files.getRunString()).append(PROJ_PROP_UPDATE_STR + projPropFile)
						.append(kvCmd).append("\n");
			}
			return cmd.append(Files.getRunString()).append(" cnv.hmm.CNVCaller proj=" + projPropFile)
								.append(" out=" + variables.get(this).get(4)).append(" ").append(PSF.Ext.NUM_THREADS_COMMAND).append(numThreads)
								.toString();
		}

	};

	static final STEP S18_SHADOW_SAMPLES = new STEP("Create 'Shadow' Project", "",
																									new String[][] {{SAMPLE_STEP_REQ_MSG},
																																	{"Number of principal components for correction."},
																																	{"Output file path (relative to project directory) and baseName for principal components correction files"},
																																	{"Call-rate filter for determining high-quality markers"},
																																	{"Re-compute Log-R Ratio values? (usually false if LRRs already exist)"},
																																	{"Temporary directory for intermediate files (which tend to be very large)"},
																																	{"Number of threads to use."},},
																									new RequirementInputType[][] {{RequirementInputType.NONE},
																																								{RequirementInputType.NUMBER},
																																								{RequirementInputType.STRING},
																																								{RequirementInputType.NUMBER},
																																								{RequirementInputType.BOOL},
																																								{RequirementInputType.DIR},
																																								{RequirementInputType.NUMBER},},
																									S2I_PARSE_SAMPLES, S2A_PARSE_SAMPLES) {
		@Override
		public void setNecessaryPreRunProperties(	Project proj,
																							Map<STEP, List<String>> variables) {
			// not needed for step
		}

		@Override
		public void run(Project proj, Map<STEP, List<String>> variables) {
			int numComponents = Integer.parseInt(variables.get(this).get(0));
			String outputBase = variables.get(this).get(1);
			double markerCallRateFilter = Double.parseDouble(variables.get(this).get(2));
			boolean recomputeLRR_PCs = Boolean.parseBoolean(variables.get(this).get(3));
			String tmpDir = "".equals(variables.get(this).get(4).trim())	? null
																																		: variables.get(this).get(4);
			int totalThreads = Integer.parseInt(variables.get(this).get(5));
			PRoCtOR.shadow(	proj, tmpDir, outputBase, markerCallRateFilter, recomputeLRR_PCs,
											numComponents, totalThreads);
		}

		@Override
		public boolean[][] checkRequirements(	Project proj, Map<STEP, Boolean> stepSelections,
																					Map<STEP, List<String>> variables) {
			int numComponents = checkIntArgOrNeg1(variables.get(this).get(0));
			double markerCallRateFilter = checkDoubleArgOrNeg1(variables.get(this).get(2));
			String tmpDir = "".equals(variables.get(this).get(4).trim())	? null
																																		: variables.get(this).get(4);
			int totalThreads = checkIntArgOrNeg1(variables.get(this).get(5));
			String sampDir = proj.SAMPLE_DIRECTORY.getValue();
			STEP parseStep = stepSelections.containsKey(S2I_PARSE_SAMPLES)	? S2I_PARSE_SAMPLES
																																			: S2A_PARSE_SAMPLES;
			boolean checkStepParseSamples = stepSelections.get(parseStep)
																			&& parseStep.hasRequirements(proj, stepSelections, variables);
			return new boolean[][] {{checkStepParseSamples
																|| (Files.exists(sampDir) && Files.list(sampDir, ".sampRAF",
																																				proj.JAR_STATUS.getValue()).length > 0)},
															{numComponents > 0}, {true}, // TODO outputBase check?
															{markerCallRateFilter >= 0}, {true},
															{tmpDir == null || Files.exists(tmpDir)}, {totalThreads > 0},};
		}

		@Override
		public Object[] getRequirementDefaults(Project proj) {
			return new Object[] {	MitoPipeline.DEFAULT_NUM_COMPONENTS, // numComponents
														MitoPipeline.FILE_BASE, // output base
														MitoPipeline.DEFAULT_MKR_CALLRATE_FILTER, // call rate
														false, // recomputeLRR
														"", // tempDir
														Runtime.getRuntime().availableProcessors() // numThreads
			};
		}

		@Override
		public boolean checkIfOutputExists(Project proj, Map<STEP, List<String>> variables) {
			String outputBase = proj.PROJECT_DIRECTORY.getValue() + MitoPipeline.FILE_BASE;
			String finalReport = outputBase + PCA.FILE_EXTs[0];
			return Files.exists(finalReport);
		}

		@Override
		public String getCommandLine(Project proj, Map<STEP, List<String>> variables) {
			int numComponents = Integer.parseInt(variables.get(this).get(0));
			String outputBase = variables.get(this).get(1);
			double markerCallRateFilter = Double.parseDouble(variables.get(this).get(2));
			boolean recomputeLRR_PCs = Boolean.parseBoolean(variables.get(this).get(3));
			String tmpDir = "".equals(variables.get(this).get(4).trim())	? null
																																		: variables.get(this).get(4);
			int totalThreads = Integer.parseInt(variables.get(this).get(5));

			String projPropFile = proj.getPropertyFilename();
			StringBuilder cmd = new StringBuilder();
			cmd	.append(Files.getRunString()).append(" org.genvisis.cnv.manage.PRoCtOR").append(" proj=")
					.append(projPropFile).append(" numComponents=").append(numComponents)
					.append(" outputBase=").append(outputBase).append(" callrate=")
					.append(markerCallRateFilter).append(" recomputeLRR=").append(recomputeLRR_PCs)
					.append(" numThreads=").append(totalThreads);
			if (tmpDir != null) {
				cmd.append(" tmp=").append(tmpDir);
			}

			return cmd.toString();
		}
	};
	static final STEP S00_TEMPLATE = new STEP("", "", new String[][] {

	}, new RequirementInputType[][] {

	}) {
		@Override
		public void setNecessaryPreRunProperties(	Project proj, Map<STEP, List<String>> variables) {
			// INSERT CODE HERE
		}

		@Override
		public void run(Project proj, Map<STEP, List<String>> variables) {
			// INSERT CODE HERE
		}

		@Override
		public boolean[][] checkRequirements(	Project proj, Map<STEP, Boolean> stepSelections,
																					Map<STEP, List<String>> variables) {
			return new boolean[][] {};
		}

		@Override
		public Object[] getRequirementDefaults(Project proj) {
			return new Object[0];
		}

		@Override
		public boolean checkIfOutputExists(Project proj, Map<STEP, List<String>> variables) {
			return false;
		}

		@Override
		public String getCommandLine(Project proj, Map<STEP, List<String>> variables) {
			return null;
		}

	};

	public abstract static class STEP {
		public String stepName;
		public String stepDesc;
		public String[][] reqs;
		private boolean failed = false;
		protected ArrayList<String> failReasons = new ArrayList<String>();
		private final Set<STEP> relatedSteps;
		public RequirementInputType[][] reqTypes;

		public boolean getFailed() {
			return failed;
		}

		protected void setFailed() {
			failed = true;
		}

		public List<String> getFailureMessages() {
			return failReasons;
		}

		protected int checkIntArgOrNeg1(String val) {
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
		
		public abstract void setNecessaryPreRunProperties(Project proj, Map<STEP, List<String>> variables);

		public abstract void run(Project proj, Map<STEP, List<String>> variables);

		public void gracefulDeath(Project proj) {
			return;
		}

		public abstract boolean[][] checkRequirements(Project proj, Map<STEP, Boolean> stepSelections, Map<STEP, List<String>> variables);

		public boolean hasRequirements(	Project proj, Map<STEP, Boolean> stepSelections, Map<STEP, List<String>> variables) {
			if (stepSelections.get(this) == null || variables.get(this) == null) {
				return false;
			}
			int sum = 0;
			boolean[][] reqrs = checkRequirements(proj, stepSelections, variables);
			for (boolean[] req : reqrs) {
				sum += Array.booleanArraySum(req) > 0 ? 1 : 0;
			}
			return sum == reqs.length;
		}

		/**
		 * @return An array of requirements. At least one element of each subarray must be met to
		 *         satisfy the step pre-requisites - effectively this means elements of the first array
		 *         are AND'd together, while elements of the second array are OR'd.
		 */
		public String[][] getRequirements() {
			// TODO unify requirement names, AND/OR structure, input types and default values to avoid
			// maintaining these parallel arrays
			return reqs;
		}

		public RequirementInputType[][] getRequirementInputTypes() {
			return reqTypes;
		}

		public abstract boolean checkIfOutputExists(Project proj, Map<STEP, List<String>> variables);

		public void resetRun() {
			failed = false;
			failReasons.clear();
		}

		/**
		 * Get the default values for requirements in the order they're set [skipping NONE input types]
		 *
		 * @param proj
		 * @return
		 */
		public abstract Object[] getRequirementDefaults(Project proj);

		public abstract String getCommandLine(Project proj, Map<STEP, List<String>> variables);

		STEP(	String name, String desc, String[][] requirements, RequirementInputType[][] reqTypes,
					STEP... requiredSteps) {
			stepName = name;
			stepDesc = desc;
			reqs = requirements;
			this.reqTypes = reqTypes;
			final Set<STEP> steps = new HashSet<STEP>();
			steps.add(this);
			for (final STEP s : requiredSteps) {
				steps.add(s);
				steps.addAll(s.getRelatedSteps());
			}
			relatedSteps = Collections.unmodifiableSet(steps);
		}

		/**
		 * @return A {@link Collection} of the complete network of {@link STEP}s related to this
		 *         {@code STEP} - including this {@code STEP}, direct and transitive dependencies.
		 */
		public Collection<STEP> getRelatedSteps() {
			return relatedSteps;
		}


	}

	public enum RequirementInputType {
		NONE(), FILE(), DIR(), STRING(), NUMBER(), BOOL()
	}

	public GenvisisWorkflow(Project project, Launch launch) {
		proj = project;
		log = project == null ? new Logger() : project.getLog();
		this.launch = launch;
	}

	public void showDialogAndRun() {
		GenvisisWorkflowGUI gui;
		gui = new GenvisisWorkflowGUI(proj, launch);
		if (!gui.getCancelled()) {
			gui.setModal(true);
			gui.setVisible(true);

			if (gui.getCancelled()) {
				return;
			}
		}
	}

	private static STEP[] ILLUMINA_STEPS = {S1I_CREATE_MKR_POS, S2I_PARSE_SAMPLES,
																					S3_CREATE_SAMPLEDATA, S4_TRANSPOSE_TO_MDF,
																					S5_COMPUTE_GCMODEL, S6_SAMPLE_QC, S7_MARKER_QC,
																					S8_SEX_CHECKS, S9_GENERATE_ABLOOKUP, S10_RUN_PLINK,
																					S11_GWAS_QC, S12_MOSAIC_ARMS, S13_ANNOTATE_SAMPLE_DATA,
																					S14_CREATE_PCS, S15_COMPUTE_PFB,
																					// S14_CREATE_MT_CN_EST,
																					S16_SEX_CENTROIDS_PFB_GCMODEL, S17_CNV_CALLING,
																					S18_SHADOW_SAMPLES,};

	private static STEP[] AFFY_STEPS = {S2A_PARSE_SAMPLES, S3_CREATE_SAMPLEDATA, S4_TRANSPOSE_TO_MDF,
																			S5_COMPUTE_GCMODEL, S6_SAMPLE_QC, S7_MARKER_QC, S8_SEX_CHECKS,
																			S9_GENERATE_ABLOOKUP, S10_RUN_PLINK, S11_GWAS_QC,
																			S12_MOSAIC_ARMS, S13_ANNOTATE_SAMPLE_DATA, S14_CREATE_PCS,
																			S15_COMPUTE_PFB,
																			// S14_CREATE_MT_CN_EST,
																			S16_SEX_CENTROIDS_PFB_GCMODEL, S17_CNV_CALLING,
																			S18_SHADOW_SAMPLES,};

	public static STEP[] getStepsForProject(Project proj) {
		switch (proj.ARRAY_TYPE.getValue()) {
			case AFFY_GW6:
			case AFFY_GW6_CN:
				return AFFY_STEPS;
			case ILLUMINA:
			default:
				return ILLUMINA_STEPS;
		}
	}

}
