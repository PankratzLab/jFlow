package seq.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.concurrent.Callable;

import seq.analysis.PlinkSeqUtils.PlinkSeqBurdenResults;
import seq.analysis.PlinkSeqUtils.PseqPhenoTypes;
import seq.analysis.PlinkSeqUtils.PseqProject;
import seq.analysis.PlinkSeqUtils.PseqProject.PROPERTIES;
import common.Array;
import common.CmdLine;
import common.Files;
import common.Logger;
import common.PSF;
import common.WorkerHive;
import common.ext;

/**
 * Wrapper for the plinkSeq package (https://atgu.mgh.harvard.edu/plinkseq/index.shtml)
 *
 */
public class PlinkSeq {
	public static String PSEQ = "pseq";
	private static final String NEW_PROJECT = "new-project";
	private static final String LOAD_VCF = "load-vcf";
	private static final String LOAD_PHENO = "load-pheno";
	private static final String LOAD_LOC = "loc-load";

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
	private static final String MAC = "mac=";

	private static final String MASK = "--mask";
	private static final String V_STATS = "v-stats";
	private static final String I_STATS = "i-stats";

	private enum ANALYSIS_TYPES {
		/**
		 * Single variant association test
		 */
		V_ASSOC, /**
		 * Burden analysis
		 */
		BURDEN, /**
		 * Summarize individuals in the project
		 */
		I_SUMMARY, /**
		 * Summarize variants in the project
		 */
		V_SUMMARY
	}

	/**
	 * See : https://atgu.mgh.harvard.edu/plinkseq/assoc.shtml
	 *
	 */
	private enum BURDEN_Tests {
		BURDEN, UNIQ, VT, FW, CALPHA, SUMSTAT
	}

	private enum LOAD_TYPES {
		VCF, PHENO, LOC_DB
	}

	// public static enum BUILDS {
	// HG_19, HG_18
	// }

	// private BUILDS build;
	private Logger log;
	private boolean fail, overwriteExisting, verbose;

	public PlinkSeq(boolean overwriteExisting, boolean verbose, Logger log) {
		this.overwriteExisting = overwriteExisting;
		this.verbose = verbose;
		this.log = log;
		this.fail = !verify();
	}

	/**
	 * Initializes the project by creating directory structures and loading dbs (if neccesary)
	 */
	private boolean createProject(PseqProject pseqProject) {
		boolean init = true;
		if (!fail) {
			String currentProject = pseqProject.getProjectDirectory() + pseqProject.getProjectName();
			String[] vcfs = getVCFCommands(pseqProject);
			String[] command = new String[] { PSEQ, pseqProject.getProjectDirectory() + pseqProject.getProjectName(), NEW_PROJECT, pseqProject.getCommandFor(PROPERTIES.RESOURCES), pseqProject.getResourceDirectory() };
			command = Array.concatAll(command, vcfs);
			String[] outputs = new String[] { currentProject + "." + PSEQ };
			String[] inputs = pseqProject.getVcfs();
			init = CmdLine.runCommandWithFileChecks(command, "", inputs, outputs, verbose, overwriteExisting, false, log);
			if (init) {
				pseqProject.load();
				init = pseqProject.isLoaded();
				if (init) {
					init = loadData(pseqProject, LOAD_TYPES.VCF, null);
				}
				if (init) {
					init = loadData(pseqProject, LOAD_TYPES.PHENO, pseqProject.getPhenoFile());
				}
				if (init) {
					String[] reqFiles = Files.listFullPaths(pseqProject.getResourceDirectory(), REQ_FILE, false);

					if (init && reqFiles != null && reqFiles.length > 0) {
						log.reportTimeInfo("Found the following " + REQ_FILE + " files to load into the loc db: ");
						log.reportTimeInfo(Array.toStr(reqFiles));
						for (int i = 0; i < reqFiles.length; i++) {
							if (init) {
								init = loadData(pseqProject, LOAD_TYPES.LOC_DB, reqFiles[i]);
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

	private void fullGamutAssoc(PseqProject pseqProject, String[] locGroups, int numPerm, String mac, String outputRoot, int numThreads) {
		WorkerHive<PlinkSeqWorker> assocHive = new WorkerHive<PlinkSeq.PlinkSeqWorker>(numThreads, 10, log);
		PseqPhenoTypes[] pseqPhenoTypes = pseqProject.getPhenotypes();
		for (int i = 0; i < pseqPhenoTypes.length; i++) {
			assocHive.addCallable(generateAWorker(pseqProject, ANALYSIS_TYPES.V_ASSOC, null, null, pseqPhenoTypes[i].getName(), numPerm, "0", outputRoot, overwriteExisting, log));
			for (int j = 0; j < locGroups.length; j++) {
				assocHive.addCallable(generateAWorker(pseqProject, ANALYSIS_TYPES.BURDEN, BURDEN_Tests.values(), locGroups[j], pseqPhenoTypes[i].getName(), numPerm, mac, outputRoot, overwriteExisting, log));
			}
		}
		assocHive.addCallable(generateAWorker(pseqProject, ANALYSIS_TYPES.I_SUMMARY, null, null, null, 0, "0", outputRoot, overwriteExisting, log));
		assocHive.addCallable(generateAWorker(pseqProject, ANALYSIS_TYPES.V_SUMMARY, null, null, null, 0, "0", outputRoot, overwriteExisting, log));
		assocHive.execute(true);
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

	private boolean loadData(PseqProject pseqProject, LOAD_TYPES lTypes, String file) {
		boolean loaded = true;
		if (!fail) {
			String loadCommand = "";
			String[] inputs = null;
			boolean overideOverWrite = false;
			if (pseqProject.isLoaded()) {
				switch (lTypes) {
				case VCF:
					loadCommand = LOAD_VCF;
					break;
				case PHENO:
					loadCommand = LOAD_PHENO;
					overideOverWrite = true;
					break;
				case LOC_DB:
					loadCommand = LOAD_LOC;
					overideOverWrite = true;
					break;
				default:
					log.reportTimeError("Invalid load Type");
					break;
				}
				String[] command = new String[] { PSEQ, pseqProject.getProjectNameForPseq(), loadCommand };
				if (file != null) {
					command = Array.concatAll(command, new String[] { FILE, file });
					if (lTypes == LOAD_TYPES.LOC_DB) {
						command = Array.concatAll(command, new String[] { GROUP, ext.rootOf(file) });
					}
					inputs = new String[] { file };
				}
				CmdLine.runCommandWithFileChecks(command, "", inputs, null, verbose, overideOverWrite, false, log);
			} else {
				loaded = false;
				log.reportTimeError("The project " + pseqProject.getProjectName() + " in directory " + pseqProject.getProjectDirectory() + " has not been loaded");
			}
		}
		return loaded;
	}

	private static String[] addPermCommand(int numPerm, String[] command) {
		if (numPerm != 0) {
			return Array.concatAll(command, new String[] { PERM, numPerm + "" });
		} else {
			return command;
		}
	}

	private static String[] getPhenoCommand(String pheno) {
		return new String[] { PHENOTYPE, pheno };
	}

	/**
	 * @param pseqProject
	 *            current project
	 * @param type
	 *            the type of analysis to perform
	 * @param bTests
	 *            if the type is {@link ANALYSIS_TYPES#BURDEN}, the burden tests to perform
	 * @param locGroups
	 *            if the type is {@link ANALYSIS_TYPES#BURDEN} or {@link ANALYSIS_TYPES#V_ASSOC} , the groupings to collapse to (KEGG, refseq,etc...)
	 * @param phenotype
	 *            the phenotype for the test, must be contained the projects phenotype file...
	 * @param numPerm
	 *            number of permutations for association testing
	 * @param outputRoot
	 *            output root of the file names
	 * @param overwriteExisting
	 * @param log
	 * @return
	 */
	private static PlinkSeqWorker generateAWorker(PseqProject pseqProject, ANALYSIS_TYPES type, BURDEN_Tests[] bTests, String locGroups, String phenotype, int numPerm, String mac, String outputRoot, boolean overwriteExisting, Logger log) {
		String outputDirectory = ext.parseDirectoryOfFile(pseqProject.getFilename()) + "assoc/";
		new File(outputDirectory).mkdirs();
		String outputFile = outputDirectory + outputRoot + "." + type;
		String[] commandBase = new String[] { PSEQ, pseqProject.getProjectNameForPseq() };
		boolean validWorker = true;
		switch (type) {
		case V_ASSOC:
			if (phenotype == null) {
				log.reportTimeError("Must provide a phentype for " + type);
				validWorker = false;
			} else {
				outputFile += "." + phenotype + ".txt";
				commandBase = Array.concatAll(commandBase, new String[] { V_ASSOC }, getPhenoCommand(phenotype));
				if (locGroups != null) {
					commandBase = Array.concatAll(commandBase, new String[] { MASK, LOC_GROUP + locGroups });
				}
			}
			commandBase = addPermCommand(numPerm, commandBase);
			break;
		case BURDEN:
			if (phenotype == null || locGroups == null) {
				log.reportTimeError("Must provide a phentype and a group for " + type);
				validWorker = false;
			} else {
				outputFile += "." + phenotype + "." + ext.replaceWithLinuxSafeCharacters(locGroups, true) + ".txt";
				commandBase = Array.concatAll(commandBase, new String[] { BURDEN }, getPhenoCommand(phenotype));
				commandBase = Array.concatAll(commandBase, new String[] { MASK, LOC_GROUP + locGroups, mac.equals("0") ? "": MAC + mac });
			}
			if (bTests != null) {
				String[] tests = new String[bTests.length + 1];
				tests[0] = TESTS;
				for (int i = 0; i < bTests.length; i++) {
					tests[i + 1] = (bTests[i] + "").toLowerCase();
				}
				commandBase = Array.concatAll(commandBase, tests);
			}
			commandBase = addPermCommand(numPerm, commandBase);

			break;
		case I_SUMMARY:
			commandBase = Array.concatAll(commandBase, new String[] { I_STATS });
			break;
		case V_SUMMARY:
			commandBase = Array.concatAll(commandBase, new String[] { V_STATS });
			break;
		default:
			break;
		}

		if (validWorker) {
			// we always pump the commands to a .bat since no --out option
			commandBase = Array.concatAll(commandBase, new String[] { PSF.Ext.CARROT, outputFile });
			String bat = ext.addToRoot(outputFile, ".bat");
			commandBase = CmdLine.prepareBatchForCommandLine(commandBase, bat, true, log);
			return new PlinkSeqWorker(type, commandBase, new String[] { bat }, new String[] { outputFile }, true, overwriteExisting, false, log);
		} else {
			return null;
		}
	}

	/**
	 * Submits jobs for pseq analyses
	 *
	 */
	private static class PlinkSeqWorker implements Callable<PlinkSeqWorker> {
		private ANALYSIS_TYPES type;
		private String[] command;
		private String[] inputFiles;
		private String[] outputFiles;
		private boolean success, verbose, overWriteExisting, skipReporting;
		private Logger log;

		private PlinkSeqWorker(ANALYSIS_TYPES type, String[] command, String[] inputFiles, String[] outputFiles, boolean verbose, boolean overWriteExisting, boolean skipReporting, Logger log) {
			super();
			this.type = type;
			this.command = command;
			this.inputFiles = inputFiles;
			this.outputFiles = outputFiles;
			this.success = false;
			this.verbose = verbose;
			this.overWriteExisting = overWriteExisting;
			this.skipReporting = skipReporting;
			this.log = log;
		}

		@Override
		public PlinkSeqWorker call() throws Exception {
			success = CmdLine.runCommandWithFileChecks(command, "", inputFiles, outputFiles, verbose, overWriteExisting, skipReporting, log);
			if (success) {
				switch (type) {
				case BURDEN:
					System.out.println(outputFiles[0]);
					PlinkSeqBurdenResults pSeqBurdenResults = new PlinkSeqBurdenResults(outputFiles[0], log);
					String burdenSummary = ext.parseDirectoryOfFile(inputFiles[0]) + "burden.summary";
					PrintWriter writer = new PrintWriter(new FileWriter(burdenSummary, true));
					HashSet<String> genesPassing = pSeqBurdenResults.getGenesPassing();
					writer.print(Array.toStr(command, " ") + "\t" + pSeqBurdenResults.getNumTests() + "\t" + pSeqBurdenResults.getBonferoniP() + "\t" + genesPassing.size());
					for (String gene : genesPassing) {
						writer.print("\t" + gene);
					}
					writer.println();
					writer.close();
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
			log.reportTimeError("It is assumed that the program " + PSEQ + " can be found on the system's path, please install before continuing");
			verified = false;
		}
		return verified;
	}

	/**
	 * Initialize a pseq project
	 * 
	 * @param pSeq
	 *            valid {@link PlinkSeq} object
	 * @param projName
	 *            name for the project
	 * @param vcf
	 *            vcf file to use(full path)
	 * @param phenotypeFile
	 *            phenotype file to use (full path)
	 * @param resourceDirectory
	 *            directory containing the PlinkSeq databases ( from https://atgu.mgh.harvard.edu/plinkseq/resources.shtml)
	 * @param log
	 * @return
	 */
	private static PseqProject initialize(PlinkSeq pSeq, String projName, String vcf, String phenoFile, String resourceDirectory, Logger log) {
		String projectName = projName == null ? ext.rootOf(vcf) : PlinkSeqUtils.PSEQ_PROJECT + projName;
		String projectDirectory = ext.parseDirectoryOfFile(vcf) + projectName + "/";
		return initialize(pSeq, projectName, projectDirectory, new String[] { vcf }, phenoFile, resourceDirectory, log);
	}

	private static PseqProject initialize(PlinkSeq pSeq, String projectName, String projectDirectory, String[] vcfs, String phenoFile, String resourceDirectory, Logger log) {
		new File(projectDirectory).mkdirs();
		String filename = projectDirectory + projectName + "." + PlinkSeq.PSEQ;
		PseqProject pseqProject = new PseqProject(filename, vcfs, phenoFile, resourceDirectory, log);
		pSeq.createProject(pseqProject);
		return pseqProject;
	}

	public static void runPlinkSeq(String projName, String vcf, String phenoFile, String resourceDirectory, String outputRoot, String[] locGroups, boolean overwriteExisting, String mac, int numThreads, Logger log) {
		PlinkSeq plinkSeq = new PlinkSeq(overwriteExisting, true, log);
		PlinkSeqUtils.PseqProject pseqProject = initialize(plinkSeq, projName, vcf, phenoFile, resourceDirectory, log);
		plinkSeq.fullGamutAssoc(pseqProject, locGroups, -1, mac, outputRoot, numThreads);
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String vcf = null;
		String phenoFile = null;
		String resourceDirectory = null;
		String projName = "pseq_project";
		String outputRoot = "pseq_analysis";
		String[] locGroups = new String[] { "refseq" };
		int numThreads = 4;
		String mac = "0";
		boolean overwriteExisting = false;

		// TODO, use just a project file Name to load everything, loc.group args (refseq)
		String usage = "\n" + "seq.analysis.PlinkSeq requires 3 arguments\n";
		usage += "   (1) full path to a vcf file (i.e. vcf= (no default))\n" + "";
		usage += "   (2) full path to a plinkseq phenotype file (i.e. phenoFile= (no default))\n" + "";
		usage += "   (3) full path to a plinkseq resource directory  (i.e. resourceDirectory= (no default))\n" + "";
		usage += "   (4) overwrite existing files  (i.e. -overwriteExisting (default))\n" + "";
		usage += "   (5) analysis output root  (i.e. outputRoot=" + outputRoot + " (default))\n" + "";
		usage += "   (6) project name  (i.e. projName=" + projName + " (default))\n" + "";
		usage += "   (7) overwrite existing files  (i.e. -overwriteExisting (default))\n" + "";
		usage += "   (8) comma -delimted loc groups to test in the association  (i.e. locGroups=" + Array.toStr(locGroups, ",") + " (default))\n" + "";
		usage += "   (9) mac cutoff for burden tests  (i.e. mac=" + mac + " (default))\n" + "";
		usage += "   (10) " + PSF.Ext.getNumThreadsCommand(10, numThreads);

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("vcf=")) {
				vcf = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("phenoFile=")) {
				phenoFile = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("outputRoot=")) {
				outputRoot = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("overwriteExisting=")) {
				overwriteExisting = true;
				numArgs--;
			} else if (args[i].startsWith("resourceDirectory=")) {
				resourceDirectory = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("locGroups=")) {
				locGroups = ext.parseStringArg(args[i], "").split(",");
				numArgs--;
			} else if (args[i].startsWith("projName=")) {
				projName = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
				numThreads = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("mac=")) {
				mac = ext.parseStringArg(args[i], "");
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			runPlinkSeq(projName, vcf, phenoFile, resourceDirectory, outputRoot, locGroups, overwriteExisting, mac, numThreads, new Logger());
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
		PseqProject pseqProject = initialize(plinkSeq, projName, testVCF, testPheno, resourceDirectory, log);
		plinkSeq.fullGamutAssoc(pseqProject, new String[] { "refseq" }, -1, "2", "test", 4);

	}

}
//
// private static PlinkSeqWorker generateWorker(PseqProject pseqProject, ANALYSIS_TYPES types, BURDEN_Tests[] bTests, String locGroups, String phenotype, int numPerm, String outputRoot, boolean overwriteExisting, Logger log) {
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
// pSeqWorkers.add(getSummaryWorker(pseqProject, types, outputRoot, overwriteExisting, log, outputDirectory, I_STATS));
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
// String outputFile = outputDirectory + outputRoot + "." + types + "." + pseqProject.getPhenotypes()[i].getName() + ".txt";
// String[] command = new String[] { PSEQ, pseqProject.getProjectNameForPseq(), V_ASSOC, PHENOTYPE, pseqProject.getPhenotypes()[i].getName() };
//
// }
// }
//
// switch (types) {
// case V_ASSOC:
// if (bTests != null) {
// log.reportTimeError("Burden tests not valid for single-variant analysis type " + ANALYSIS_TYPES.V_ASSOC);
// }
// for (int i = 0; i < pseqProject.getPhenotypes().length; i++) {
// String outputFile = outputDirectory + outputRoot + "." + types + "." + pseqProject.getPhenotypes()[i].getName() + ".txt";
// String[] command = new String[] { PSEQ, pseqProject.getProjectNameForPseq(), V_ASSOC, PHENOTYPE, pseqProject.getPhenotypes()[i].getName() };
// if (numPerm != 0) {
// command = Array.concatAll(command, getnperms(numPerm));
// }
// command = Array.concatAll(command, new String[] { PSF.Ext.CARROT, outputFile });
// String bat = ext.addToRoot(outputFile, ".bat");
// command = CmdLine.prepareBatchForCommandLine(command, bat, true, log);
// pSeqWorkers.add(new PlinkSeqWorker(command, new String[] { bat }, new String[] { outputFile }, true, overwriteExisting, false, log));
// }
//
// break;
// case BURDEN:
// if (locGroups == null || locGroups.length == 0) {
// log.reportTimeError("Must provide at least one collapsable set type for burden testing");
// } else {
// for (int i = 0; i < pseqProject.getPhenotypes().length; i++) {
// String[] command = new String[] { PSEQ, pseqProject.getProjectNameForPseq(), BURDEN, PHENOTYPE, pseqProject.getPhenotypes()[i].getName() };
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
// String outputFile = outputDirectory + outputRoot + "." + types + "." + pseqProject.getPhenotypes()[i].getName() + "." + locGroups[j] + ".txt";
//
// String[] finalCommand = Array.concatAll(command, new String[] { MASK, LOC_GROUP + locGroups[j], PSF.Ext.CARROT, outputFile });
// String bat = ext.addToRoot(outputFile, ".bat");
// finalCommand = CmdLine.prepareBatchForCommandLine(finalCommand, bat, true, log);
// pSeqWorkers.add(new PlinkSeqWorker(finalCommand, new String[] { bat }, new String[] { outputFile }, true, overwriteExisting, false, log));
// }
// }
// }
// break;
//
// case I_SUMMARY:
//
// break;
// case V_SUMMARY:
// pSeqWorkers.add(getSummaryWorker(pseqProject, types, outputRoot, overwriteExisting, log, outputDirectory, V_STATS));
//
// break;
// default:
// break;
// }
// return pSeqWorkers.toArray(new PlinkSeqWorker[pSeqWorkers.size()]);
// }
//
// private static PlinkSeqWorker getSummaryCommand(PseqProject pseqProject, ANALYSIS_TYPES types, String outputRoot, boolean overwriteExisting, Logger log, String outputDirectory, String type) {
//
// return new PlinkSeqWorker(finalCommandV, new String[] { batV }, new String[] { outputFileV }, true, overwriteExisting, false, log);
// }
