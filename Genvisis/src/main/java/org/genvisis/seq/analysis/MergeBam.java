package seq.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import common.Array;
import common.CmdLine;
import common.Files;
import common.Logger;
import common.ext;

//TODO when we bring in the htsjdk we can switch to merging bamfiles natively 
public class MergeBam {
	public static final String SAMTOOLS_COMMAND = "samtools=";
	private static final String SAMTOOLS_LOCATION = "samtools";

	private static final String MERGE = "merge";
	private static final String VIEW = "view";
	private static final String REHEADER = "reheader";

	private static final String BAM = ".bam";
	// private static final String R = "-r";// to add read groups while merging
	private static final String H = "-H";
	private static final String SMALL_H = "-h";
	private String samtoolsLocation;
	private boolean verbose, fail, overwriteExisting;
	private Logger log;

	public MergeBam(String samtoolsLocation, boolean overwriteExisting, boolean verbose, Logger log) {
		super();
		this.samtoolsLocation = samtoolsLocation;
		this.verbose = verbose;
		this.log = log;
		this.overwriteExisting = overwriteExisting;
		this.log = log;
		this.fail = !validSamtools();
	}

	public String getSamtoolsLocation() {
		return samtoolsLocation;
	}

	public boolean isFail() {
		return fail;
	}

	public void setFail(boolean fail) {
		this.fail = fail;
	}

	public boolean mergeSomeBams(String[] baseIds, String newBaseID, String[] inputBams, String outputBam, Logger altLog) {
		boolean progress = true;
		if (!fail) {
			if (inputBams.length >= 2) {
				ArrayList<String> tmpCommand = new ArrayList<String>();
				tmpCommand.add(SAMTOOLS_LOCATION);
				tmpCommand.add(MERGE);

				// tmpCommand.add(R);
				String outputHeader = ext.addToRoot(outputBam, ".header");

				tmpCommand.add(outputBam);
				for (int i = 0; i < inputBams.length; i++) {
					tmpCommand.add(inputBams[i]);
					// TODO
					// TODO
					// TODO

					getFullHeader(samtoolsLocation, inputBams[i], outputHeader, i == 0, altLog);// TODO, remove this when RGs are formatted properly

				}
				if (baseIds.length != 1) {
					String originalHeader = outputHeader;
					outputHeader = outputHeader + ".rh";
					String regex = "";
					String sed = buildSedCommand(baseIds, newBaseID, outputHeader, originalHeader, regex);
					log.report(ext.getTime() + " Info - running command " + sed + " to properly format header with new ID");
					String batFile = outputHeader + ".sed";
					Files.write(sed, batFile);
					Files.chmod(batFile);
					progress = CmdLine.runCommandWithFileChecks(new String[] { batFile }, "", new String[] { originalHeader }, new String[] { outputHeader }, verbose, true, false, altLog);
				}
				if (progress) {
					tmpCommand.add(SMALL_H);
					tmpCommand.add(outputHeader);
					String[] command = tmpCommand.toArray(new String[tmpCommand.size()]);
					progress = CmdLine.runCommandWithFileChecks(command, "", inputBams, new String[] { outputBam }, verbose, overwriteExisting, true, (altLog == null ? log : altLog));
				}
			} else {
				log.report(ext.getTime() + " Info - since there were less than two input bams, " + Array.toStr(inputBams) + " will not be merged");
				progress = true;
			}
		} else {
			log.reportError("Error - could not merge bam files...most likely samtools was not seen");
			progress = false;
		}
		return progress;
	}

	private String buildSedCommand(String[] baseIds, String newBaseID, String outputHeader, String originalHeader, String regex) {
		for (int i = 0; i < baseIds.length; i++) {
			if (i == 0) {
				regex = "SM:" + baseIds[i];
			} else {
				regex += "\\|SM:" + baseIds[i];
			}
		}
		String sed = "sed \"s/" + regex + "/SM:" + newBaseID + "/\" " + originalHeader + " > " + outputHeader;
		return sed;
	}

	public BamMerger mergeABam(String[] baseIds, String newBaseID, String[] inputBams, String outputDir, String mergeStage, Logger altLog) {
		altLog.report("trying to merge " + newBaseID + "\t" + Array.toStr(inputBams));
		BamMerger bamMerger = new BamMerger(newBaseID, outputDir, inputBams, altLog);
		bamMerger.parse(mergeStage);
		boolean progress = true;
		if (!fail && bamMerger.shouldMerge()) {
			progress = mergeSomeBams(baseIds, newBaseID, bamMerger.getInputBams(), bamMerger.getOutputBam(), altLog);
		}
		bamMerger.setFail(!progress);
		return bamMerger;
	}

	public ReHeader reHejaderBamFilePriorToMerge(String bamFile, String oldSM, String newSM, Logger log) {
		ReHeader reHeader = new ReHeader(bamFile);
		reHeader.parse();
		boolean progress = true;
		String[] output = new String[] { reHeader.getReHeaderBam() };
		log.report(ext.getTime() + " Info - sample name " + oldSM + " will be replaced with " + newSM + " in  the header of " + bamFile);
		String sed = "sed \"s/SM:" + oldSM + "/SM:" + newSM + "/\"";
		String[] command = new String[] { samtoolsLocation, VIEW, H, bamFile, "|", sed, "|", samtoolsLocation, REHEADER, "-", bamFile, ">", reHeader.getReHeaderBam() };
		String bat = ext.rootOf(bamFile, false) + ".rh.bat";
		Files.write(Array.toStr(command, " "), bat);
		Files.chmod(bat);
		progress = CmdLine.runCommandWithFileChecks(new String[] { bat }, "", new String[] { bamFile, bat }, output, verbose, overwriteExisting, false, log);
		reHeader.setFail(!progress);
		return reHeader;
	}

	private boolean validSamtools() {
		if (samtoolsLocation != null && !samtoolsLocation.equals("") && !samtoolsLocation.equals(SAMTOOLS_LOCATION)) {
			if (Files.exists(samtoolsLocation)) {
				if (verbose) {
					log.report(ext.getTime() + " Info - using samtools located at " + samtoolsLocation);
				}
				return true;
			} else {
				log.reportError("Error - invalid samtools location, could not find" + samtoolsLocation);
				return false;
			}
		} else {
			if (CmdLine.run(SAMTOOLS_LOCATION, "")) {
				if (verbose) {
					log.report(ext.getTime() + " Info - using samtools set by the system's path varaible");
				}
				samtoolsLocation = SAMTOOLS_LOCATION;
				return true;
			} else {
				log.reportError("Error - a path to samtools was not supplied and bwa was not detected on the system's path");
				return false;
			}
		}
	}

	public static class ReHeader {
		public static final String RE_HEADER = ".rh";
		private String inputBam, reHeaderBam;
		private boolean fail;

		public ReHeader(String inputBam) {
			super();
			this.inputBam = inputBam;
			this.fail = false;
		}

		public void parse() {
			this.reHeaderBam = ext.addToRoot(inputBam, RE_HEADER);
		}

		public boolean isFail() {
			return fail;
		}

		public void setFail(boolean fail) {
			this.fail = fail;
		}

		public String getReHeaderBam() {
			return reHeaderBam;
		}

	}

	public static class BamMerger {
		public static final String[] MERGE_STAGES = { "rrd_lane_" };
		public static final String MERGE = ".merge";

		private String baseId, outputDir, outputBam;
		private String[] inputBams;
		private boolean shouldMerge, fail;
		private Logger log;

		public BamMerger(String baseId, String outputDir, String[] inputBams, Logger log) {
			super();
			this.baseId = baseId;
			this.outputDir = outputDir;
			this.inputBams = inputBams;
			this.shouldMerge = inputBams.length > 1 && inputBams != null;
			this.fail = false;
			this.log = log;
		}

		public void parse(String mergeStage) {
			new File(outputDir).mkdirs();
			if (shouldMerge) {
				this.outputBam = outputDir + mergeStage + baseId + MERGE + BAM;
			} else {
				if (inputBams != null) {
					this.outputBam = inputBams[0];
				}
			}
		}

		public String getBaseId() {
			return baseId;
		}

		public void setFail(boolean fail) {
			this.fail = fail;
		}

		public boolean isFail() {
			return fail;
		}

		public String getOutputDir() {
			return outputDir;
		}

		public String getOutputBam() {
			return outputBam;
		}

		public String[] getInputBams() {
			return inputBams;
		}

		public boolean shouldMerge() {
			return shouldMerge;
		}

		public Logger getLog() {
			return log;
		}

	}

	private static String getFullHeader(String samtoolsLocation, String inputBam, String outputHeader, boolean full, Logger log) {
		String headerFile = ext.addToRoot(inputBam, ".header");
		String[] fullHeaderCommand = new String[] { samtoolsLocation, "view", H, inputBam, ">", headerFile };
		String bat = ext.addToRoot(inputBam, ".header.bat");
		Files.write(Array.toStr(fullHeaderCommand, " "), bat);
		Files.chmod(bat);
		CmdLine.run(bat, "");
		// CmdLine.runCommandWithFileChecks(fullHeaderCommand, "", new String[] { inputBam }, new String[] { headerFile }, true, true, true, log);
		if (full) {
			Files.copyFile(headerFile, outputHeader);
		} else {
			String readGroup = null;
			try {
				BufferedReader reader = Files.getAppropriateReader(headerFile);
				while (reader.ready()) {
					String line = reader.readLine().trim();
					if (line.startsWith("@RG")) {
						readGroup = line;
					}
				}
				reader.close();
				try {
					PrintWriter writer = new PrintWriter(new FileWriter(outputHeader, true));
					writer.println(readGroup);
					writer.close();
				} catch (Exception e) {
					log.reportError("Error writing to " + outputHeader);
					log.reportException(e);
				}
			} catch (FileNotFoundException fnfe) {
				log.reportError("Error: file \"" + headerFile + "\" not found in current directory");
			} catch (IOException ioe) {
				log.reportError("Error reading file \"" + headerFile + "\"");
			}
		}
		return outputHeader;
	}

//	private static void checkSampleName(String headerFile, String baseId, Logger log) {
//		String readGroup;
//		try {
//			BufferedReader reader = Files.getAppropriateReader(headerFile);
//			while (reader.ready()) {
//				String line = reader.readLine().trim();
//				if (line.startsWith("@RG")) {
//					readGroup = line;
//				}
//			}
//			reader.close();
//
//		} catch (FileNotFoundException fnfe) {
//			log.reportError("Error: file \"" + headerFile + "\" not found in current directory");
//		} catch (IOException ioe) {
//			log.reportError("Error reading file \"" + headerFile + "\"");
//		}
//	}

	// private static String replaceSMifNeeded(String baseID, String readGroupLine, Logger log) {
	// String[] curRGInfo = readGroupLine.split("\t");
	// //int indexOfSample = ext.indexOfStartsWith(, curRGInfo, false);
	// // if (indexOfSample >= 0) {
	// // // if(!)
	// // } else {
	// // log.reportError("Error -  could not find the sample \"SM:\" id in readgroup line " + readGroupLine + ", this should not happen");
	// // }
	//
	// }

}
