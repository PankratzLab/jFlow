package org.genvisis.gwas;

import java.util.ArrayList;
import java.util.concurrent.Callable;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.WorkerHive;
import org.genvisis.common.ext;

/**
 * Class for handling 23 and me files (https://www.cog-genomics.org/plink2/input#23file)<br>
 * Converts to plink format and merges all available files
 */
public class MeAnd23 {
	private static final String PROBLEMATIC_MARKERS = "problematic_markers.txt";

	public static void generatePlinkFiles(String directory, String outputRoot, int numThreads,
																				Logger log) {
		String[] allFiles = Files.toFullPaths(Files.list(directory, ".txt", false), directory);
		log.reportTimeInfo("Found " + allFiles.length + " potential 23andMe files with extension .txt");

		PlinkGenerator[] plinkGenerators = new PlinkGenerator[allFiles.length];
		for (int i = 0; i < plinkGenerators.length; i++) {
			plinkGenerators[i] = new PlinkGenerator(allFiles[i], log);
		}

		WorkerHive<PlinkGenerator> hive = new WorkerHive<PlinkGenerator>(numThreads, 10, log);
		hive.addCallables(plinkGenerators);
		hive.execute(true);
		plinkGenerators = hive.getResults().toArray(new PlinkGenerator[plinkGenerators.length]);

		log.reportTimeInfo("Generating list of files to merge for non-failing conversions");

		ArrayList<String> toMerge = new ArrayList<String>();
		for (PlinkGenerator plinkGenerator : plinkGenerators) {
			if (plinkGenerator.isCreated()) {
				toMerge.add(ArrayUtils.toStr(plinkGenerator.getPlinks()));
			}
		}
		log.reportTimeInfo("Found " + toMerge.size() + " file(s) to merge");

		String[] inputs = toMerge.toArray(new String[toMerge.size()]);
		String toMergeFile = directory + "23AndMeMerge.txt";
		Files.writeArray(inputs, toMergeFile);
		String mergedRoot = Files.isRelativePath(outputRoot) ? directory + outputRoot : outputRoot;
		mergeResults(toMergeFile, mergedRoot, log);

	}

	private static boolean mergeResults(String toMergeFile, String mergedRoot, Logger log) {
		log.reportTimeInfo("Attempting to merge files listed in "	+ toMergeFile + " to plink root "
												+ mergedRoot);
		boolean merged = true;
		boolean addedExclude = false;
		String[] command = new String[] {	"plink2", "--noweb", "--make-bed", "--merge-list", toMergeFile,
																			"--out", mergedRoot};
		String problematicMarkers = ext.parseDirectoryOfFile(toMergeFile) + PROBLEMATIC_MARKERS;
		if (Files.exists(problematicMarkers)) {
			command = addExclude(command, problematicMarkers);
			addedExclude = true;
		}
		merged = CmdLine.runCommandWithFileChecks(command, "", new String[] {toMergeFile},
																							PSF.Plink.getPlinkBedBimFam(mergedRoot), true, true,
																							false, log);
		if (!merged) {
			String[] missnps = Files.list(ext.parseDirectoryOfFile(mergedRoot), ".missnp", false);
			if (missnps.length > 0) {
				if (!addedExclude) {
					command = addExclude(command, problematicMarkers);
				}
				if (!Files.exists(problematicMarkers)) {
					Files.write("", problematicMarkers);
				}
				Files.copyFile(problematicMarkers, problematicMarkers + "2");
				String[] cats = ArrayUtils.concatAll(missnps, new String[] {problematicMarkers + "2"});
				Files.cat(cats, problematicMarkers, null, log);
				log.reportTimeInfo("Job failed, but since .missnp file(s) exist, we are going to try and just remove these snps (and any others listed in "
														+ problematicMarkers + ") ");
				merged = CmdLine.runCommandWithFileChecks(command, "", new String[] {toMergeFile},
																									PSF.Plink.getPlinkBedBimFam(mergedRoot), true,
																									true, false, log);
			}
		}
		if (!merged) {
			log.reportError("Could not merge files. If the issues was strand flips /.missnp related, try running again until all problematic markers have been found and removed");
		}
		return merged;

	}

	private static String[] addExclude(String[] command, String problematicMarkers) {
		command = ArrayUtils.concatAll(command, new String[] {"--exclude", problematicMarkers});
		return command;
	}

	/**
	 * Generates plink files for as single sample from 23 and me results
	 *
	 */
	private static class PlinkGenerator implements Callable<PlinkGenerator> {
		private final String input;
		private final String outRoot;
		private String name;
		private boolean created;
		private final Logger log;
		private final String[] plinks;

		public PlinkGenerator(String input, Logger log) {
			super();
			this.input = input;
			name = ext.rootOf(input);
			if (name.length() >= 38) {
				name = name.replaceAll("genome_Patient__Full_", "");
				name = name.replaceAll("genome_", "");
				name = name.replaceAll("Full_", "");
			}
			this.log = log;
			plinks = PSF.Plink.getPlinkBedBimFam(ext.rootOf(input, false));
			outRoot = ext.rootOf(input, false);
			created = false;
			getOutRoot();

		}

		public String getOutRoot() {
			return outRoot;
		}

		public boolean isCreated() {
			return created;
		}

		public String[] getPlinks() {
			return plinks;
		}

		@Override
		public PlinkGenerator call() throws Exception {
			String[] command = new String[] {	"plink2", "--make-bed", "--23file", input, name, name,
																				"--out", outRoot, "--noweb"};
			String problematicMarkers = ext.parseDirectoryOfFile(input) + PROBLEMATIC_MARKERS;
			if (Files.exists(problematicMarkers)) {
				command = addExclude(command, problematicMarkers);
			}
			log.reportTimeInfo(Thread.currentThread().getName()	+ ":" + name + " ( "
													+ ArrayUtils.toStr(command, " ") + " )");
			created = CmdLine.runCommandWithFileChecks(	command, "", new String[] {input}, plinks, true,
																									true, true, log);
			return this;
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String directory = null;
		String outputRoot = "MeAnd23_merge";
		int numThreads = 1;// be careful, not tested when excluding problematic markers
		Logger log;

		String usage = "\n" + "gwas.MeAnd23 requires 0-1 arguments\n";
		usage += "   (1) input directory with 23AndMe files  (i.e. directory= (no default))\n" + "";
		usage += "   (2) output Root (i.e. outputRoot=" + outputRoot + " (default))\n" + "";
		// usage += " (3) number of threads (i.e. " + PSF.Ext.NUM_THREADS_COMMAND + numThreads + "
		// (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("directory=")) {
				directory = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("outputRoot=")) {
				outputRoot = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith(PSF.Ext.NUM_THREADS_COMMAND)) {
				numThreads = ext.parseIntArg(arg);
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
			if (directory == null) {
				System.err.println("Error - must provide an input directory");
			} else {
				log = new Logger(directory + outputRoot + ".genvisis.log");
				generatePlinkFiles(directory, outputRoot, numThreads, log);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
