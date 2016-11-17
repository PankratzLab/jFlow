package org.genvisis.seq.analysis;

import java.io.File;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.concurrent.Callable;

import org.genvisis.common.Array;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.WorkerHive;
import org.genvisis.common.WorkerTrain.AbstractProducer;
import org.genvisis.common.ext;
import org.genvisis.seq.analysis.PlinkSeqUtils.PseqPhenoTypes;
import org.genvisis.seq.analysis.PlinkSeqUtils.PseqProject;
import org.genvisis.seq.analysis.PlinkSeqUtils.PseqProject.PROPERTIES;
import org.genvisis.seq.manage.VCFOps.VcfPopulation;

/**
 * Wrapper for the plinkSeq package (https://atgu.mgh.harvard.edu/plinkseq/index.shtml)
 *
 */
public class PlinkSeq implements Serializable {
	/**
	 *
	 */
	private static final long serialVersionUID = 1L;
	public static String PSEQ = "pseq";
	private static final String NEW_PROJECT = "new-project";
	private static final String LOAD_VCF = "load-vcf";
	private static final String LOAD_PHENO = "load-pheno";
	private static final String LOAD_LOC = "loc-load";
	private static final String VAR_DROP_ALL_SETS = "var-drop-all-sets";
	private static final String VAR_SET = "var-set";
	private static final String FILE = "--file";
	private static final String GROUP = "--group";

	private static final String V_ASSOC = "v-assoc";
	private static final String BURDEN = "assoc";

	private static final String PHENOTYPE = "--phenotype";
	private static final String PERM = "--perm";
	private static final String TESTS = "--tests";
	private static final String REQ_FILE = ".reg";
	/**
	 * Command for how the burden tests are collapsed in the locdb
	 */
	private static final String LOC_GROUP = "loc.group=";
	private static final String MAF = "maf=";

	private static final String MASK = "--mask";
	private static final String VAR_MASK = "var.req=";
	private static final String V_STATS = "v-stats";
	private static final String I_STATS = "i-stats";

	public enum ANALYSIS_TYPES {
															/**
															 * Single variant association test
															 */
															V_ASSOC,
															/**
															 * Burden analysis
															 */
															BURDEN,
															/**
															 * Summarize individuals in the project
															 */
															I_SUMMARY,
															/**
															 * Summarize variants in the project
															 */
															V_SUMMARY
	}

	/**
	 * See : https://atgu.mgh.harvard.edu/plinkseq/assoc.shtml
	 *
	 */
	public enum BURDEN_Tests {
														BURDEN
		// , UNIQ, VT, FW, CALPHA, SUMSTAT, FRQWGT
	}

	public enum LOAD_TYPES {
													VCF, PHENO, LOC_DB
	}

	// public static enum BUILDS {
	// HG_19, HG_18
	// }

	// private BUILDS build;
	private final Logger log;
	private final boolean fail, overwriteExisting, verbose;

	public PlinkSeq(boolean overwriteExisting, boolean verbose, Logger log) {
		this.overwriteExisting = overwriteExisting;
		this.verbose = verbose;
		this.log = log;
		fail = !verify();
	}

	private boolean dropAllVarSets(PseqProject pseqProject) {
		boolean progress = false;
		String[] command = new String[] {PSEQ, pseqProject.getProjectNameForPseq(), VAR_DROP_ALL_SETS};
		progress = CmdLine.runCommandWithFileChecks(command, "", null, null, true, false, false, log);
		return progress;
	}

	private boolean loadVarSet(PseqProject pseqProject, String varSetFile) {
		boolean progress = false;
		String[] command = new String[] {	PSEQ, pseqProject.getProjectNameForPseq(), VAR_SET, FILE,
																			varSetFile};
		progress = CmdLine.runCommandWithFileChecks(command, "", new String[] {varSetFile}, null, true,
																								false, false, log);
		return progress;
	}

	public boolean eraseAndLoadVarSet(PseqProject pseqProject, String varSetFile) {
		boolean progress = false;
		progress = dropAllVarSets(pseqProject);
		if (progress) {
			progress = loadVarSet(pseqProject, varSetFile);
		}
		return progress;

	}

	/**
	 * Initializes the project by creating directory structures and loading dbs (if neccesary)
	 */
	private boolean createProject(PseqProject pseqProject, boolean loadVCF, boolean loadReq) {
		boolean init = true;
		if (!fail) {
			init = createNewProject(pseqProject);
			if (init) {
				pseqProject.load();
				init = pseqProject.isLoaded();
				if (init && loadVCF) {
					init = loadData(pseqProject, LOAD_TYPES.VCF, null);
				}
				if (init) {
					init = loadData(pseqProject, LOAD_TYPES.PHENO, pseqProject.getPhenoFile());
				}
				if (init) {
					String[] reqFiles = Files.listFullPaths(pseqProject.getResourceDirectory(), REQ_FILE,
																									false);

					if (init && loadReq && reqFiles != null && reqFiles.length > 0) {
						log.reportTimeInfo("Found the following "	+ REQ_FILE
																+ " files to load into the loc db: ");
						log.reportTimeInfo(Array.toStr(reqFiles));
						for (String reqFile : reqFiles) {
							if (init) {
								init = loadData(pseqProject, LOAD_TYPES.LOC_DB, reqFile);
							}
						}
					}
				}
			}
		} else {
			init = false;
		}
		return init;
	}

	public boolean createNewProject(PseqProject pseqProject) {
		boolean init;
		String currentProject = pseqProject.getProjectDirectory() + pseqProject.getProjectName();
		String[] vcfs = getVCFCommands(pseqProject);
		String[] command = new String[] {	PSEQ,
																			pseqProject.getProjectDirectory()
																						+ pseqProject.getProjectName(),
																			NEW_PROJECT, pseqProject.getCommandFor(PROPERTIES.RESOURCES),
																			pseqProject.getResourceDirectory()};
		command = Array.concatAll(command, vcfs);
		String[] outputs = new String[] {currentProject + "." + PSEQ};
		String[] inputs = pseqProject.getVcfs();
		init = CmdLine.runCommandWithFileChecks(command, "", inputs, outputs, verbose,
																						overwriteExisting, false, log);
		return init;
	}

	public PlinkSeqWorker[] fullGamutAssoc(	PseqProject pseqProject, String[] locGroups,
																					String[] varMasks, int numPerm, String mac,
																					String outputRoot, int numThreads) {
		return fullGamutAssoc(pseqProject, locGroups, varMasks, numPerm, mac, outputRoot, true,
													numThreads);
	}

	public PlinkSeqWorker[] fullGamutAssoc(	PseqProject pseqProject, String[] locGroups,
																					String[] varMasks, int numPerm, String mac,
																					String outputRoot, boolean execute, int numThreads) {
		ArrayList<PlinkSeqWorker> workers = new ArrayList<PlinkSeqWorker>();
		PseqPhenoTypes[] pseqPhenoTypes = pseqProject.getPhenotypes();
		if (pseqPhenoTypes != null) {
			for (PseqPhenoTypes pseqPhenoType : pseqPhenoTypes) {
				if (varMasks == null) {
					workers.add(generateAWorker(pseqProject, ANALYSIS_TYPES.V_ASSOC, null, null, null,
																			pseqPhenoType.getName(), numPerm, "0", outputRoot,
																			overwriteExisting, log));
					for (String locGroup : locGroups) {
						workers.add(generateAWorker(pseqProject, ANALYSIS_TYPES.BURDEN, BURDEN_Tests.values(),
																				locGroup, null, pseqPhenoType.getName(), numPerm, mac,
																				outputRoot, overwriteExisting, log));
					}
				} else {
					for (String varMask : varMasks) {
						workers.add(generateAWorker(pseqProject, ANALYSIS_TYPES.V_ASSOC, null, null, varMask,
																				pseqPhenoType.getName(), numPerm, "0", outputRoot,
																				overwriteExisting, log));
						for (String locGroup : locGroups) {
							workers.add(generateAWorker(pseqProject, ANALYSIS_TYPES.BURDEN, BURDEN_Tests.values(),
																					locGroup, varMask, pseqPhenoType.getName(), numPerm, mac,
																					outputRoot, overwriteExisting, log));
						}
					}
				}
			}
		}
		workers.add(generateAWorker(pseqProject, ANALYSIS_TYPES.I_SUMMARY, null, null, null, null, 0,
																"0", outputRoot, overwriteExisting, log));
		workers.add(generateAWorker(pseqProject, ANALYSIS_TYPES.V_SUMMARY, null, null, null, null, 0,
																"0", outputRoot, overwriteExisting, log));
		if (execute) {
			WorkerHive<PlinkSeqWorker> assocHive = new WorkerHive<PlinkSeq.PlinkSeqWorker>(	numThreads, 10,
																																											log);
			assocHive.addCallables(workers.toArray(new PlinkSeqWorker[workers.size()]));
			assocHive.execute(true);
			ArrayList<PlinkSeqWorker> complete = assocHive.getResults();
			return complete.toArray(new PlinkSeqWorker[complete.size()]);
		} else {
			return workers.toArray(new PlinkSeqWorker[workers.size()]);
		}
	}

	private static String[] getVCFCommands(PseqProject pseqProject) {
		String[] vcfs = new String[pseqProject.getVcfs().length * 2];
		int index = 0;
		for (int i = 0; i < pseqProject.getVcfs().length; i++) {
			vcfs[index] = pseqProject.getCommandFor(PROPERTIES.VCF);
			index++;
			vcfs[index] = pseqProject.getVcfs()[i];
			index++;
		}
		return vcfs;
	}

	public boolean loadData(PseqProject pseqProject, LOAD_TYPES lTypes, String file) {
		boolean loaded = true;
		if (!fail) {
			String loadCommand = "";
			String[] inputs = null;
			boolean overideOverWrite = false;
			if (pseqProject.isLoaded()) {
				switch (lTypes) {
					case VCF:
						loadCommand = LOAD_VCF;
						// pseqProject.getLog().reportTimeError("JOHN TURN THIS OFF LATER");
						// return true;
						break;
					case PHENO:
						loadCommand = LOAD_PHENO;
						overideOverWrite = true;
						// pseqProject.getLog().reportTimeError("JOHN TURN THIS OFF LATER");
						// return true;
						break;
					case LOC_DB:
						loadCommand = LOAD_LOC;
						overideOverWrite = true;
						// pseqProject.getLog().reportTimeError("JOHN TURN THIS OFF LATER");
						// return true;
						break;
					default:
						log.reportError("Invalid load Type " + lTypes);
						break;
				}
				String[] command = new String[] {PSEQ, pseqProject.getProjectNameForPseq(), loadCommand};
				if (file != null) {
					command = Array.concatAll(command, new String[] {FILE, file});
					if (lTypes == LOAD_TYPES.LOC_DB) {
						command = Array.concatAll(command, new String[] {GROUP, ext.rootOf(file)});
					}
					inputs = new String[] {file};
				}
				CmdLine.runCommandWithFileChecks(	command, "", inputs, null, verbose, overideOverWrite,
																					false, log);
			} else {
				loaded = false;
				log.reportError("The project "	+ pseqProject.getProjectName() + " in directory "
														+ pseqProject.getProjectDirectory() + " has not been loaded");
			}
		}
		return loaded;
	}

	private static String[] addPermCommand(int numPerm, String[] command) {
		if (numPerm != 0) {
			return Array.concatAll(command, new String[] {PERM, numPerm + ""});
		} else {
			return command;
		}
	}

	private static String[] getPhenoCommand(String pheno) {
		return new String[] {PHENOTYPE, pheno};
	}

	/**
	 * @param pseqProject current project
	 * @param type the type of analysis to perform
	 * @param bTests if the type is {@link ANALYSIS_TYPES#BURDEN}, the burden tests to perform
	 * @param locGroups if the type is {@link ANALYSIS_TYPES#BURDEN} or {@link ANALYSIS_TYPES#V_ASSOC}
	 *        , the groupings to collapse to (KEGG, refseq,etc...)
	 * @param phenotype the phenotype for the test, must be contained the projects phenotype file...
	 * @param numPerm number of permutations for association testing
	 * @param outputRoot output root of the file names
	 * @param overwriteExisting
	 * @param log
	 * @return
	 */
	public static PlinkSeqWorker generateAWorker(	PseqProject pseqProject, ANALYSIS_TYPES type,
																								BURDEN_Tests[] bTests, String locGroups,
																								String varMask, String phenotype, int numPerm,
																								String maf, String outputRoot,
																								boolean overwriteExisting, Logger log) {
		String outputDirectory = ext.parseDirectoryOfFile(pseqProject.getFilename()) + "assoc/";
		new File(outputDirectory).mkdirs();
		String outputFile = outputDirectory + outputRoot + "." + type;
		String[] commandBase = new String[] {PSEQ, pseqProject.getProjectNameForPseq()};
		boolean validWorker = true;
		switch (type) {
			case V_ASSOC:
				if (phenotype == null) {
					log.reportError("Must provide a phentype for " + type);
					validWorker = false;
				} else {
					outputFile += "." + phenotype + ".txt";
					commandBase = Array.concatAll(commandBase, new String[] {V_ASSOC},
																				getPhenoCommand(phenotype));
					if (locGroups != null) {
						commandBase = Array.concatAll(commandBase, new String[] {MASK, LOC_GROUP + locGroups});
						ext.addToRoot(outputFile, locGroups);

					}
					// if (varMask != null) {
					// ext.addToRoot(outputFile, varMask);
					// commandBase = Array.concatAll(commandBase, new String[] { MASK, varMask });
					// }
				}
				commandBase = addPermCommand(numPerm, commandBase);
				break;
			case BURDEN:
				if (phenotype == null || locGroups == null) {
					log.reportError("Must provide a phentype and a group for " + type);
					validWorker = false;
				} else {
					outputFile += "."	+ phenotype + "." + ext.replaceWithLinuxSafeCharacters(locGroups, true)
												+ ".txt";
					commandBase = Array.concatAll(commandBase, new String[] {BURDEN},
																				getPhenoCommand(phenotype));
					commandBase = Array.concatAll(commandBase,
																				new String[] {MASK, LOC_GROUP + locGroups,
																											maf.equals("0") ? "" : MAF + "0:" + maf});
					if (varMask != null) {
						outputFile = ext.addToRoot(outputFile, "." + varMask);
						commandBase = Array.concatAll(commandBase, new String[] {VAR_MASK + varMask});
					}
					if (!maf.equals("0")) {
						outputFile = ext.addToRoot(outputFile, "maf_" + maf);
					}
				}
				if (bTests != null) {
					String[] tests = new String[bTests.length + 1]; // note that we skip FRQWGT here
					tests[0] = TESTS;
					int index = 1;
					for (BURDEN_Tests bTest : bTests) {
						// if (bTests[i] != BURDEN_Tests.FRQWGT) {
						tests[index] = (bTest + "").toLowerCase();
						index++;
						// }
					}
					commandBase = Array.concatAll(commandBase, tests);
				}
				commandBase = addPermCommand(numPerm, commandBase);

				break;
			case I_SUMMARY:
				commandBase = Array.concatAll(commandBase, new String[] {I_STATS});
				break;
			case V_SUMMARY:
				commandBase = Array.concatAll(commandBase, new String[] {V_STATS});
				break;
			default:
				break;
		}

		if (validWorker) {
			// we always pump the commands to a .bat since no --out option
			commandBase = Array.concatAll(commandBase, new String[] {PSF.Ext.CARROT, outputFile});
			String bat = ext.addToRoot(outputFile, ".bat");
			commandBase = CmdLine.prepareBatchForCommandLine(commandBase, bat, true, log);
			return new PlinkSeqWorker(type, commandBase, new String[] {bat}, new String[] {outputFile},
																true, overwriteExisting, false, log);
		} else {
			return null;
		}
	}

	/**
	 * Submits jobs for pseq analyses
	 *
	 */
	public static class PlinkSeqWorker implements Callable<PlinkSeqWorker> {
		private final ANALYSIS_TYPES type;
		private final String[] command;
		private final String[] inputFiles;
		private final String[] outputFiles;
		private boolean success;
		private final boolean verbose;
		private final boolean overWriteExisting;
		private final boolean skipReporting;
		private final Logger log;

		private PlinkSeqWorker(	ANALYSIS_TYPES type, String[] command, String[] inputFiles,
														String[] outputFiles, boolean verbose, boolean overWriteExisting,
														boolean skipReporting, Logger log) {
			super();
			this.type = type;
			this.command = command;
			this.inputFiles = inputFiles;
			this.outputFiles = outputFiles;
			success = false;
			this.verbose = verbose;
			this.overWriteExisting = overWriteExisting;
			this.skipReporting = skipReporting;
			this.log = log;
		}

		public ANALYSIS_TYPES getType() {
			return type;
		}

		public String[] getOutputFiles() {
			return outputFiles;
		}

		@Override
		public PlinkSeqWorker call() throws Exception {
			success = CmdLine.runCommandWithFileChecks(	command, "", inputFiles, outputFiles, verbose,
																									overWriteExisting, skipReporting, log);
			if (success) {
				switch (type) {
					case BURDEN:
						// System.out.println(outputFiles[0]);
						// PlinkSeqBurdenResults pSeqBurdenResults = new PlinkSeqBurdenResults(outputFiles[0],
						// log);
						// String burdenSummary = ext.parseDirectoryOfFile(inputFiles[0]) + "burden.summary";
						// PrintWriter writer = new PrintWriter(new FileWriter(burdenSummary, true));
						// HashSet<String> genesPassing = pSeqBurdenResults.getGenesPassing();
						// writer.print(Array.toStr(command, " ") + "\t" + pSeqBurdenResults.getNumTests() +
						// "\t" + pSeqBurdenResults.getBonferoniP() + "\t" + genesPassing.size());
						// for (String gene : genesPassing) {
						// writer.print("\t" + gene);
						// }
						// writer.println();
						// writer.close();
						break;
					case I_SUMMARY:
						break;
					case V_ASSOC:
						break;
					case V_SUMMARY:
						break;
					default:
						break;

				}
				isSuccess();

			}
			return this;
			// TODO Auto-generated method stub
		}

		private boolean isSuccess() {
			return success;
		}

	}

	private boolean verify() {
		boolean verified = true;
		if (!CmdLine.run(PSEQ, "")) {
			log.reportError("It is assumed that the program "	+ PSEQ
													+ " can be found on the system's path, please install before continuing");
			verified = false;
		}
		return verified;
	}

	public static PseqProject initialize(	PlinkSeq pSeq, String projName, String directory, String vcf,
																				VcfPopulation vpop, String resourceDirectory,
																				boolean loadVCF, boolean loadReq, Logger log) {
		String projectName = projName == null	? PlinkSeqUtils.PSEQ_PROJECT + ext.rootOf(vcf)
																					: PlinkSeqUtils.PSEQ_PROJECT + projName;
		String projectDirectory = directory == null	? ext.parseDirectoryOfFile(vcf) + projectName + "/"
																								: directory;
		new File(projectDirectory).mkdirs();
		String phenoFile = null;
		if (vpop != null) {
			phenoFile = projectDirectory + ext.rootOf(vpop.getFileName()) + ".pheno";
			vpop.generatePlinkSeqPheno(phenoFile);
		}
		return initialize(pSeq, projectName, projectDirectory, new String[] {vcf},
											vpop == null ? null : phenoFile, resourceDirectory, loadVCF, loadReq, log);
	}

	/**
	 * Initialize a pseq project
	 *
	 * @param pSeq valid {@link PlinkSeq} object
	 * @param projName name for the project
	 * @param vcf vcf file to use(full path)
	 * @param phenotypeFile phenotype file to use (full path)
	 * @param resourceDirectory directory containing the PlinkSeq databases ( from
	 *        https://atgu.mgh.harvard.edu/plinkseq/resources.shtml)
	 * @param log
	 * @return
	 */
	private static PseqProject initialize(PlinkSeq pSeq, String projName, String vcf,
																				String phenoFile, String resourceDirectory, boolean loadVCF,
																				boolean loadReq, Logger log) {
		String projectName = projName == null ? ext.rootOf(vcf) : PlinkSeqUtils.PSEQ_PROJECT + projName;
		String projectDirectory = ext.parseDirectoryOfFile(vcf) + projectName + "/";
		return initialize(pSeq, projectName, projectDirectory, new String[] {vcf}, phenoFile,
											resourceDirectory, loadVCF, loadReq, log);
	}

	public static PseqProject initialize(	PlinkSeq pSeq, String projectName, String projectDirectory,
																				String[] vcfs, String phenoFile, String resourceDirectory,
																				boolean loadVCF, boolean loadReq, Logger log) {
		new File(projectDirectory).mkdirs();
		String filename = projectDirectory + projectName + "." + PlinkSeq.PSEQ;
		PseqProject pseqProject = new PseqProject(filename, vcfs, phenoFile, resourceDirectory, log);
		pSeq.createProject(pseqProject, loadVCF, loadReq);
		return pseqProject;
	}

	public static void runPlinkSeq(	String projName, String vcf, String phenoFile,
																	String resourceDirectory, String outputRoot, String[] locGroups,
																	boolean overwriteExisting, boolean loadVCF, boolean loadReq,
																	String mac, int numThreads, Logger log) {
		PlinkSeq plinkSeq = new PlinkSeq(overwriteExisting, true, log);
		PlinkSeqUtils.PseqProject pseqProject = initialize(	plinkSeq, projName, vcf, phenoFile,
																												resourceDirectory, loadVCF, loadReq, log);
		plinkSeq.fullGamutAssoc(pseqProject, locGroups, null, -1, mac, outputRoot, numThreads);
	}

	public static class PlinkSeqProducer extends AbstractProducer<PlinkSeqWorker> {
		private final PlinkSeqWorker[] plinkSeqWorkers;
		private int index;
		// private Logger log;

		public PlinkSeqProducer(PlinkSeqWorker[] plinkSeqWorkers, Logger log) {
			super();
			this.plinkSeqWorkers = plinkSeqWorkers;
			// this.log = log;
			index = 0;
		}

		@Override
		public boolean hasNext() {
			// TODO Auto-generated method stub
			return index < plinkSeqWorkers.length;
		}

		@Override
		public Callable<PlinkSeqWorker> next() {
			PlinkSeqWorker tmp = plinkSeqWorkers[index];
			index++;
			return tmp;
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String vcf = null;
		String phenoFile = null;
		String resourceDirectory = null;
		String projName = "pseq_project";
		String outputRoot = "pseq_analysis";
		String[] locGroups = new String[] {"refseq"};
		int numThreads = 4;
		String mac = "0";
		boolean overwriteExisting = false;

		// TODO, use just a project file Name to load everything, loc.group args (refseq)
		String usage = "\n" + "seq.analysis.PlinkSeq requires 3 arguments\n";
		usage += "   (1) full path to a vcf file (i.e. vcf= (no default))\n" + "";
		usage += "   (2) full path to a plinkseq phenotype file (i.e. phenoFile= (no default))\n" + "";
		usage +=
					"   (3) full path to a plinkseq resource directory  (i.e. resourceDirectory= (no default))\n"
							+ "";
		usage += "   (4) overwrite existing files  (i.e. -overwriteExisting (default))\n" + "";
		usage += "   (5) analysis output root  (i.e. outputRoot=" + outputRoot + " (default))\n" + "";
		usage += "   (6) project name  (i.e. projName=" + projName + " (default))\n" + "";
		usage += "   (7) overwrite existing files  (i.e. -overwriteExisting (default))\n" + "";
		usage += "   (8) comma -delimted loc groups to test in the association  (i.e. locGroups="
							+ Array.toStr(locGroups, ",") + " (default))\n" + "";
		usage += "   (9) mac cutoff for burden tests  (i.e. mac=" + mac + " (default))\n" + "";
		usage += "   (10) " + PSF.Ext.getNumThreadsCommand(10, numThreads);

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("vcf=")) {
				vcf = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("phenoFile=")) {
				phenoFile = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("outputRoot=")) {
				outputRoot = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("overwriteExisting=")) {
				overwriteExisting = true;
				numArgs--;
			} else if (arg.startsWith("resourceDirectory=")) {
				resourceDirectory = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("locGroups=")) {
				locGroups = ext.parseStringArg(arg, "").split(",");
				numArgs--;
			} else if (arg.startsWith("projName=")) {
				projName = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
				numThreads = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("mac=")) {
				mac = ext.parseStringArg(arg, "");
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + arg);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			runPlinkSeq(projName, vcf, phenoFile, resourceDirectory, outputRoot, locGroups,
									overwriteExisting, true, true, mac, numThreads, new Logger());
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void test() {
		// System.out.println(CmdLine.getCmdLocation("java"));
		String testVCF = "/home/tsaim/shared/Project_Tsai_Project_021/vcf_filter/ExomeChip_AIMs/joint_genotypes.SNP.recal.INDEL.recal.eff.gatk.hg19_multianno.GATK_filter.EXOME_AIMS.vcf";
		String testPheno = "/panfs/roc/groups/14/tsaim/shared/Project_Tsai_Project_021/vcf_filter/ExomeChip_AIMs/pseqProj_PseqTest/testPheno.txt";
		String resourceDirectory = "/home/pankrat2/public/bin/pseqRef/hg19/";
		String projName = "PseqTest";
		Logger log = new Logger();
		PlinkSeq plinkSeq = new PlinkSeq(true, true, log);
		PseqProject pseqProject = initialize(	plinkSeq, projName, testVCF, testPheno, resourceDirectory,
																					true, true, log);
		plinkSeq.fullGamutAssoc(pseqProject, new String[] {"refseq"}, null, -1, "2", "test", 4);

	}

}
//
// private static PlinkSeqWorker generateWorker(PseqProject pseqProject, ANALYSIS_TYPES types,
// BURDEN_Tests[] bTests, String locGroups, String phenotype, int numPerm, String outputRoot,
// boolean overwriteExisting, Logger log) {
// String outputDirectory = ext.parseDirectoryOfFile(pseqProject.getFilename()) + "assoc/";
// new File(outputDirectory).mkdirs();
// String outputFileBase = outputDirectory + outputRoot + "." + types;
// String[] commandBase = new String[] { PSEQ, pseqProject.getProjectNameForPseq() };
// ArrayList<PlinkSeqWorker> pSeqWorkers = new ArrayList<PlinkSeqWorker>();
// switch (types) {
// case V_ASSOC:
//
// break;
// case BURDEN:
// if (bTests != null) {
// String[] tests = new String[bTests.length + 1];
// tests[0] = TESTS;
// for (int j = 0; j < bTests.length; j++) {
// tests[j + 1] = (bTests[j] + "").toLowerCase();
// }
// commandBase = Array.concatAll(commandBase, tests);
// }
// break;
// case I_SUMMARY:
// pSeqWorkers.add(getSummaryWorker(pseqProject, types, outputRoot, overwriteExisting, log,
// outputDirectory, I_STATS));
// break;
// case V_SUMMARY:
// break;
// default:
// break;
// }
//
// if (types == ANALYSIS_TYPES.V_ASSOC || types == ANALYSIS_TYPES.BURDEN) {
// for (int i = 0; i < pseqProject.getPhenotypes().length; i++) {
//
// String outputFile = outputDirectory + outputRoot + "." + types + "." +
// pseqProject.getPhenotypes()[i].getName() + ".txt";
// String[] command = new String[] { PSEQ, pseqProject.getProjectNameForPseq(), V_ASSOC, PHENOTYPE,
// pseqProject.getPhenotypes()[i].getName() };
//
// }
// }
//
// switch (types) {
// case V_ASSOC:
// if (bTests != null) {
// log.reportTimeError("Burden tests not valid for single-variant analysis type " +
// ANALYSIS_TYPES.V_ASSOC);
// }
// for (int i = 0; i < pseqProject.getPhenotypes().length; i++) {
// String outputFile = outputDirectory + outputRoot + "." + types + "." +
// pseqProject.getPhenotypes()[i].getName() + ".txt";
// String[] command = new String[] { PSEQ, pseqProject.getProjectNameForPseq(), V_ASSOC, PHENOTYPE,
// pseqProject.getPhenotypes()[i].getName() };
// if (numPerm != 0) {
// command = Array.concatAll(command, getnperms(numPerm));
// }
// command = Array.concatAll(command, new String[] { PSF.Ext.CARROT, outputFile });
// String bat = ext.addToRoot(outputFile, ".bat");
// command = CmdLine.prepareBatchForCommandLine(command, bat, true, log);
// pSeqWorkers.add(new PlinkSeqWorker(command, new String[] { bat }, new String[] { outputFile },
// true, overwriteExisting, false, log));
// }
//
// break;
// case BURDEN:
// if (locGroups == null || locGroups.length == 0) {
// log.reportTimeError("Must provide at least one collapsable set type for burden testing");
// } else {
// for (int i = 0; i < pseqProject.getPhenotypes().length; i++) {
// String[] command = new String[] { PSEQ, pseqProject.getProjectNameForPseq(), BURDEN, PHENOTYPE,
// pseqProject.getPhenotypes()[i].getName() };
//
// if (bTests != null) {
// String[] tests = new String[bTests.length + 1];
// tests[0] = TESTS;
// for (int j = 0; j < bTests.length; j++) {
// tests[j + 1] = (bTests[j] + "").toLowerCase();
// }
// command = Array.concatAll(command, tests);
// }
// for (int j = 0; j < locGroups.length; j++) {
// String outputFile = outputDirectory + outputRoot + "." + types + "." +
// pseqProject.getPhenotypes()[i].getName() + "." + locGroups[j] + ".txt";
//
// String[] finalCommand = Array.concatAll(command, new String[] { MASK, LOC_GROUP + locGroups[j],
// PSF.Ext.CARROT, outputFile });
// String bat = ext.addToRoot(outputFile, ".bat");
// finalCommand = CmdLine.prepareBatchForCommandLine(finalCommand, bat, true, log);
// pSeqWorkers.add(new PlinkSeqWorker(finalCommand, new String[] { bat }, new String[] { outputFile
// }, true, overwriteExisting, false, log));
// }
// }
// }
// break;
//
// case I_SUMMARY:
//
// break;
// case V_SUMMARY:
// pSeqWorkers.add(getSummaryWorker(pseqProject, types, outputRoot, overwriteExisting, log,
// outputDirectory, V_STATS));
//
// break;
// default:
// break;
// }
// return pSeqWorkers.toArray(new PlinkSeqWorker[pSeqWorkers.size()]);
// }
//
// private static PlinkSeqWorker getSummaryCommand(PseqProject pseqProject, ANALYSIS_TYPES types,
// String outputRoot, boolean overwriteExisting, Logger log, String outputDirectory, String type) {
//
// return new PlinkSeqWorker(finalCommandV, new String[] { batV }, new String[] { outputFileV },
// true, overwriteExisting, false, log);
// }
