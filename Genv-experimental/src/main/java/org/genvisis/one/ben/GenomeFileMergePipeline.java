package org.genvisis.one.ben;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map.Entry;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.gwas.Qc;


public class GenomeFileMergePipeline {

	private static final String[] RELATEDS_COLUMNS = {"FID1", "IID1", "FID2", "IID2", "P(IBD=0)",
																										"P(IBD=1)", "P(IBD=2)", "PI_HAT", "Type"};
	private static final String[] GENOME_COLUMNS = {"FID1", "IID1", "FID2", "IID2", "Z0", "Z1", "Z2",
																									"PI_HAT"};

	private static final String SEARCHFOR_GENOME = "plink.genome";
	private static final String SEARCHFOR_RELATEDS_SUFF = "genome_relateds.xln";

	private final HashMap<String, String> idMap = new HashMap<String, String>();
	private final ArrayList<Project> plinkProjects = new ArrayList<Project>();
	private final ArrayList<Project> qcProjects = new ArrayList<Project>();
	private final ArrayList<Project> projects = new ArrayList<Project>();
	private final HashSet<String> dataLabelSet = new HashSet<String>();
	private final ArrayList<GenomeFileSet> files = new ArrayList<GenomeFileMergePipeline.GenomeFileSet>();

	private boolean runPlinkOrQCIfMissing = false;
	private String outputFile = "audit.xln";
	// private String prevAuditFile = null;

	private final Logger log = new Logger();

	private static class GenomeFileSet {
		public GenomeFileSet(String name2, String relFile, String genFile, Logger log) {
			name = name2;
			relatedsFile = relFile;
			genomeFile = genFile;
			relLineCount = Files.countLines(relatedsFile, 1);
			// this.genLineCount = Files.countLines(genomeFile, 1);
			log.report("Counted " + relLineCount + " lines in " + relatedsFile);
		}

		String name;
		String relatedsFile;
		String genomeFile;
		int relLineCount;
		// int genLineCount;
	}

	public void loadIDLookupFile(String file, boolean ignoreFirstLine) {
		String[][] ids = HashVec.loadFileToStringMatrix(file, ignoreFirstLine, new int[] {0, 1}, false);
		for (int i = 0; i < ids.length; i++) {
			if (idMap.containsKey(ids[i][0])) {
				log.report("Warning - duplicate FID found on line " + i + ": " + ids[i][0]);
			}
			if (idMap.containsKey(ids[i][1])) {
				log.report("Warning - duplicate IID found on line " + i + ": " + ids[i][1]);
			}
			idMap.put(ids[i][0], ids[i][1]);
			idMap.put(ids[i][1], ids[i][0]);
		}
	}

	public void setRunPlinkOrQCIfMissing(boolean run) {
		runPlinkOrQCIfMissing = run;
	}

	public void addProject(String propertiesFile) {
		Project proj = new Project(propertiesFile, false);
		String plinkDir = proj.PROJECT_DIRECTORY.getValue() + "plink/";
		String genDir = plinkDir + Qc.GENOME_DIR;
		String genFile = genDir = "plink.genome";
		String[] rel = Files.list(genDir, SEARCHFOR_RELATEDS_SUFF, false);
		boolean plink = false;
		boolean qc = false;
		String msg;
		if (!Files.exists(plinkDir)) {
			msg = "no \"plink/\" directory in project directory " + proj.PROJECT_DIRECTORY.getValue();
			if (runPlinkOrQCIfMissing) {
				log.reportTime("Warning - " + msg);
				log.reportTime("PLINK files will be created and QC'd for project "
												+ proj.PROJECT_NAME.getValue());
				plink = true;
			} else {
				log.reportError("Error - " + msg);
				return;
			}
		} else if (!Files.exists(genDir)) {
			msg = "no \"genome/\" directory in QC folders for PLINK data; looking for " + genDir;
			if (runPlinkOrQCIfMissing) {
				log.reportTime("Warning - " + msg);
				log.reportTime("PLINK files will be QC'd for project " + proj.PROJECT_NAME.getValue());
				qc = true;
			} else {
				log.reportError("Error - " + msg);
				return;
			}
		} else if (!Files.exists(genFile)) {
			msg = "no plink.genome file found for project " + proj.PROJECT_NAME.getValue();
			if (runPlinkOrQCIfMissing) {
				log.reportTime("Warning - " + msg);
				log.reportTime("PLINK files will be QC'd for project "	+ proj.PROJECT_NAME.getValue()
												+ "; however, the output folders were found but the plink.genome file was missing, so other errors may be present.");
				qc = true;
			} else {
				log.reportError("Error - " + msg);
				return;
			}
		} else if (rel == null || rel.length == 0 || "".equals(rel[0]) || !Files.exists(rel[0])) {
			msg = "no *.genome_relateds.xln file found for project " + proj.PROJECT_NAME.getValue();
			if (runPlinkOrQCIfMissing) {
				log.reportTime("Warning - " + msg);
				log.reportTime("PLINK files will be QC'd for project "	+ proj.PROJECT_NAME.getValue()
												+ "; however, the output folders were found but the *.genome_relateds.xln file was missing, so other errors may be present.");
				qc = true;
			} else {
				log.reportError("Error - " + msg);
				return;
			}
		}
		projects.add(proj);
		if (plink) {
			plinkProjects.add(proj);
		}
		if (qc) {
			qcProjects.add(proj);
		}
	}

	private boolean checkRelatedsFile(String relFile) {
		if (new File(relFile).length() == 0) {
			log.reportError("Error - genome relateds file "	+ relFile
													+ " is empty! Data from this file will not be included in the final output.");
			return false;
		}

		return true;
	}

	private boolean checkGenomeFile(String genFile) {
		BufferedReader reader;
		String line;
		int[] factors;

		if (new File(genFile).length() == 0) {
			log.reportError("Error - genome file "	+ genFile
													+ " is empty! Data from this file will not be included in the final output.");
			return false;
		}

		try {
			reader = Files.getAppropriateReader(genFile);
			line = reader.readLine();
			reader.close();
			factors = ext.indexFactors(GENOME_COLUMNS, line.trim().split("[\\s]+"), false, false);
			for (int i = 0; i < factors.length; i++) {
				if (factors[i] == -1) {
					log.reportError("Error - column "	+ GENOME_COLUMNS[i]
															+ " was missing from plink.genome file: " + genFile
															+ " ; data from this file will not be included in final output.");
					return false;
				}
			}
			return true;
		} catch (IOException e) {
			log.reportException(e);
			return false;
		}
	}

	public void addFiles(String name, String genomeDir) {
		String genDir = ext.verifyDirFormat(genomeDir);
		if (!Files.exists(genDir)) {
			log.reportError("Error - genome file directory \""	+ genomeDir
													+ "\" not found! Data set will be dropped: " + name);
			return;
		}
		if (!Files.exists(genDir + SEARCHFOR_GENOME)) {
			log.reportError("Error - genome file \""	+ (genDir + SEARCHFOR_GENOME)
													+ "\" not found! Data set will be dropped: " + name);
			return;
		}
		String[] rel = Files.list(genDir, SEARCHFOR_RELATEDS_SUFF, false);
		if (rel == null || rel.length == 0 || "".equals(rel[0]) || !Files.exists(genDir + rel[0])) {
			log.reportError("Error - genome relateds file matching filename pattern \"*"
														+ SEARCHFOR_RELATEDS_SUFF + "\" not found in directory " + genomeDir
													+ "! Data set will be dropped: " + name);
			return;
		}
		if (dataLabelSet.contains(name)) {
			log.reportError("Error - duplicate data set name: "	+ name
													+ ". Directory source will be dropped: " + genomeDir);
			return;
		}
		boolean genPass = checkGenomeFile(genDir + SEARCHFOR_GENOME);
		if (genPass) {
			boolean relPass = checkRelatedsFile(genDir + rel[0]);
			if (relPass) {
				dataLabelSet.add(name);
				files.add(new GenomeFileSet(name, genDir + rel[0], genDir + SEARCHFOR_GENOME, log));
			}
		}
	}

	public void addFiles(String name, String relatedsFile, String genomeFile) {
		if (!Files.exists(genomeFile)) {
			log.reportError("Error - genome file \""	+ genomeFile
													+ "\" not found! Data set will be dropped: " + name);
			return;
		}
		if (!Files.exists(relatedsFile)) {
			log.reportError("Error - genome relateds file \""	+ relatedsFile
													+ "\" not found! Data set will be dropped: " + name);
			return;
		}
		if (dataLabelSet.contains(name)) {
			log.reportError("Error - duplicate data set name: "	+ name
													+ ". File set will be dropped: " + relatedsFile + "  |  " + genomeFile);
			return;
		}
		boolean genPass = checkGenomeFile(genomeFile);
		if (genPass) {
			boolean relPass = checkRelatedsFile(relatedsFile);
			if (relPass) {
				dataLabelSet.add(name);
				files.add(new GenomeFileSet(name, relatedsFile, genomeFile, log));
			}
		}
	}

	// public void addPreviousRunFile(String prevFile) {
	// this.prevAuditFile = prevFile;
	// }

	public void setOutputFile(String outFile) {
		outputFile = outFile;
	}

	public void run() {
		BufferedReader reader;
		String line, outKey1, outKey2;
		String fid1, fid2, iid1, iid2, type;
		String ibd0, ibd1, ibd2, piHt;
		String[] outLine;
		int outLineCount;
		int projInd0, projInd1, projInd2, projInd3, typeInd;
		int[] factors;
		HashMap<String, String[]> outLineMap;
		PrintWriter writer;

		// for all projects in plinkProjects, run PLINK export, add all to qcProjects

		// for all projects in qcProjects, run gwas.Qc

		for (int p = 0; p < projects.size(); p++) {
			String name = projects.get(p).PROJECT_NAME.getValue();
			String dir = projects.get(p).PROJECT_DIRECTORY.getValue() + "plink/" + Qc.GENOME_DIR;
			addFiles(name, dir);
		}

		System.gc();


		// if not null, load previous audit file -- IGNORE FOR NOW
		// load all genome_relateds files
		// search all genome files for pairings in genome_relateds (and prev audit file if applicable)
		// files
		// write output

		BigInteger max = new BigInteger("0");
		for (GenomeFileSet f : files) {
			max.add(new BigInteger("" + f.relLineCount)); // only outputting data for relateds entries
		}
		outLineCount = 5 + (5 * files.size()); // fid1 iid1 fid2 iid2 + 5 columns per projects (ibd0,
																						// ibd1, ibd2, pi_hat, type) + diffFlag
		outLineMap = new HashMap<String, String[]>(max.intValue(), 0.95f);

		for (int p = 0; p < files.size(); p++) {
			log.reportTime("Loading relateds data from " + files.get(p).relatedsFile);
			projInd0 = 4 + p * 5 + 0;
			projInd1 = projInd0 + 1;
			projInd2 = projInd1 + 1;
			projInd3 = projInd2 + 1;
			typeInd = projInd3 + 1;

			int[] cols = ext.indexFactors(RELATEDS_COLUMNS,
																		Files.getHeaderOfFile(files.get(p).relatedsFile, log), false,
																		log, true, false);
			String[][] relatedsData = HashVec.loadFileToStringMatrix(	files.get(p).relatedsFile, true,
																																cols, false);

			for (String[] parts : relatedsData) {
				fid1 = parts[0];
				iid1 = parts[1];
				fid2 = parts[2];
				iid2 = parts[3];
				ibd0 = parts[4];
				ibd1 = parts[5];
				ibd2 = parts[6];
				piHt = parts[7];
				type = parts[8];

				if (idMap.containsKey(iid1) && !idMap.containsKey(fid1)) {
					fid1 = idMap.get(iid1);
				} else if (idMap.containsKey(fid1) && !idMap.containsKey(iid1)) {
					iid1 = idMap.get(fid1);
				}
				if (idMap.containsKey(iid2) && !idMap.containsKey(fid2)) {
					fid2 = idMap.get(iid2);
				} else if (idMap.containsKey(fid2) && !idMap.containsKey(iid2)) {
					iid2 = idMap.get(fid2);
				}

				outKey1 = fid1 + "|" + iid1 + "||" + fid2 + "|" + iid2;
				outKey2 = fid2 + "|" + iid2 + "||" + fid1 + "|" + iid1;

				if (outLineMap.containsKey(outKey1)) {
					outLine = outLineMap.get(outKey1);
				} else if (outLineMap.containsKey(outKey2)) {
					outLine = outLineMap.get(outKey2);
					// TODO will IBD/PIHAT need to be altered due to flipped ids?
				} else {
					outLine = ArrayUtils.stringArray(outLineCount, ".");
					outLine[0] = fid1;
					outLine[1] = iid1;
					outLine[2] = fid2;
					outLine[3] = iid2;
					outLineMap.put(outKey1, outLine);
				}

				outLine[projInd0] = ibd0;
				outLine[projInd1] = ibd1;
				outLine[projInd2] = ibd2;
				outLine[projInd3] = piHt;
				outLine[typeInd] = type;
			}

			System.gc();
		}

		String[] parts;
		for (int p = 0; p < files.size(); p++) {
			log.reportTime("Scanning genome file \""	+ files.get(p).genomeFile
											+ "\" for possible matches.");
			projInd0 = 4 + p * 5 + 0;
			projInd1 = projInd0 + 1;
			projInd2 = projInd1 + 1;
			projInd3 = projInd2 + 1;
			typeInd = projInd3 + 1;
			try {
				reader = Files.getAppropriateReader(files.get(p).genomeFile);
				line = reader.readLine();
				factors = ext.indexFactors(GENOME_COLUMNS, line.trim().split("[\\s]+"), false, false);

				while ((line = reader.readLine()) != null) {
					line = line.trim();
					if ("".equals(line)) {
						continue;
					}
					parts = ext.splitLine(line, "[\\s]+", log);

					fid1 = parts[factors[0]];
					iid1 = parts[factors[1]];
					fid2 = parts[factors[2]];
					iid2 = parts[factors[3]];
					ibd0 = parts[factors[4]];
					ibd1 = parts[factors[5]];
					ibd2 = parts[factors[6]];
					piHt = parts[factors[7]];

					if (idMap.containsKey(iid1) && !idMap.containsKey(fid1)) {
						fid1 = idMap.get(iid1);
					} else if (idMap.containsKey(fid1) && !idMap.containsKey(iid1)) {
						iid1 = idMap.get(fid1);
					}
					if (idMap.containsKey(iid2) && !idMap.containsKey(fid2)) {
						fid2 = idMap.get(iid2);
					} else if (idMap.containsKey(fid2) && !idMap.containsKey(iid2)) {
						iid2 = idMap.get(fid2);
					}

					outKey1 = fid1 + "|" + iid1 + "||" + fid2 + "|" + iid2;
					outKey2 = fid2 + "|" + iid2 + "||" + fid1 + "|" + iid1;

					if (outLineMap.containsKey(outKey1)) {
						outLine = outLineMap.get(outKey1);
					} else if (outLineMap.containsKey(outKey2)) {
						outLine = outLineMap.get(outKey2);
						// TODO will IBD/PIHAT need to be altered due to flipped ids?
					} else {
						continue; // do not add entries from genome file unless already existing in relateds
											// file
					}

					outLine[projInd0] = ibd0;
					outLine[projInd1] = ibd1;
					outLine[projInd2] = ibd2;
					outLine[projInd3] = piHt;
					if ("".equals(outLine[typeInd]) || ".".equals(outLine[typeInd])) {
						outLine[typeInd] = "Unrelated";
					}
				}

				reader.close();
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}

			System.gc();
		}

		// write outmap (options for exclusions? split pops, ea/aa, etc?)
		log.reportTime("Writing output to " + outputFile);

		writer = Files.getAppropriateWriter(outputFile);

		outLine = ArrayUtils.stringArray(outLineCount, "");
		for (int p = 0; p < files.size(); p++) {
			outLine[4 + (5 * p)] = files.get(p).name;
		}
		writer.println(ArrayUtils.toStr(outLine, "\t"));
		outLine = ArrayUtils.stringArray(outLineCount, "");
		outLine[0] = "FID1";
		outLine[1] = "IID1";
		outLine[2] = "FID2";
		outLine[3] = "IID2";
		for (int i = 0; i < files.size(); i++) {
			int fileInd = i * 5;
			outLine[4 + fileInd] = "P(IBD=0)"; // 4 9 14
			outLine[4 + fileInd + 1] = "P(IBD=1)"; // 5 10
			outLine[4 + fileInd + 2] = "P(IBD=2)"; // 6 11
			outLine[4 + fileInd + 3] = "PI_HAT"; // 7 12
			outLine[4 + fileInd + 4] = "Type"; // 8 13
		}
		outLine[outLine.length - 1] = "DIFF_FLAG";
		writer.println(ArrayUtils.toStr(outLine, "\t"));

		HashSet<String> types = new HashSet<String>();
		for (Entry<String, String[]> lines : outLineMap.entrySet()) {
			outLine = lines.getValue();
			types.clear();
			for (int p = 0; p < files.size(); p++) {
				typeInd = 4 + p * 5 + 4;
				if (!".".equals(outLine[typeInd].trim())) {
					types.add(outLine[typeInd]);
				}
			}
			outLine[outLine.length - 1] = (types.size() > 1) ? "1" : "0";
			writer.println(ArrayUtils.toStr(outLine, "\t"));
		}

		writer.flush();
		writer.close();

	}

	private static void loadDataFile(GenomeFileMergePipeline gfmp, String file) {
		String[] lines = HashVec.loadFileToStringArray(file, false, null, false);
		// first line is ids file
		gfmp.loadIDLookupFile(lines[0], false);
		// middle lines are {LABEL <tab> RELATEDS_FILE <tab> GENOME_FILE}
		for (int i = 1; i < lines.length - 1; i++) {
			String[] parts = lines[i].split("\t");
			gfmp.addFiles(parts[0], parts[1], parts[2]);
		}
		// last line is output file
		gfmp.setOutputFile(lines[lines.length - 1]);
	}

	public static void main(String[] args) {

		GenomeFileMergePipeline gfmp = new GenomeFileMergePipeline();
		gfmp.setRunPlinkOrQCIfMissing(false);
		loadDataFile(gfmp, "/scratch.global/cole0482/genomeFiles/data.txt");
		gfmp.run();

	}


}
