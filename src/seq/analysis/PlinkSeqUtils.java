package seq.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Properties;

import seq.analysis.PlinkSeq.BURDEN_Tests;
import common.Array;
import common.Files;
import common.HashVec;
import common.Logger;
import common.Positions;
import common.ext;
import filesys.GeneData;
import filesys.GeneTrack;

/**
 * Handles the project and phenotype system used by PlinkSeq
 *
 */
public class PlinkSeqUtils {

	public static final String PSEQ_PROJECT = "pseqProj_";

	/**
	 * Class for the phenotypes used in pseq, manages the file for later association testing
	 *
	 */
	public static class PseqPhenoTypes {
		private String name, type, missingValue, description;
		private Logger log;

		public PseqPhenoTypes(String[] phenoLine, Logger log) {
			this(phenoLine[0].replace("##", ""), phenoLine[1], phenoLine[2], phenoLine[3], log);
		}

		public PseqPhenoTypes(String name, String type, String missingValue, String description, Logger log) {
			super();
			this.name = name;
			this.type = type;
			this.missingValue = missingValue;
			this.description = description;
			this.log = log;
		}

		public String getName() {
			return name;
		}

		public String getType() {
			return type;
		}

		public String getMissingValue() {
			return missingValue;
		}

		public String getDescription() {
			return description;
		}

		public Logger getLog() {
			return log;
		}

		/**
		 * @param phenoFile
		 *            pseq formated phenotype file (See: https://atgu.mgh.harvard.edu/plinkseq/input.shtml )
		 * @param log
		 * @return
		 */
		public static PseqPhenoTypes[] loadPhenos(String phenoFile, Logger log) {
			ArrayList<PseqPhenoTypes> phenotypes = new ArrayList<PlinkSeqUtils.PseqPhenoTypes>();
			if (phenoFile == null) {
				log.reportTimeError("No phenotype file was provided, continuing carefully");
				return null;
			}
			try {
				BufferedReader reader = Files.getAppropriateReader(phenoFile);
				boolean done = false;
				while (reader.ready() && !done) {
					String line = reader.readLine().trim();
					if (line.startsWith("##")) {
						String[] phenoLine = line.split(",");
						if (phenoLine.length != 4) {
							log.reportTimeError("Phenotype lines starting with ## must have 4 comma-delimited entries");
							return null;
						} else {
							phenotypes.add(new PseqPhenoTypes(phenoLine, log));
						}
					} else {
						done = true;
					}

				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				log.reportError("Error: file \"" + phenoFile + "\" not found in current directory");
			} catch (IOException ioe) {
				log.reportError("Error reading file \"" + phenoFile + "\"");
			}
			if (phenotypes.size() > 0) {
				log.reportTimeInfo("Found " + phenotypes.size() + " phenotype(s) listed in " + phenoFile);
			} else {
				log.reportTimeError("Did not find any valid phenotype headers listed in " + phenoFile);
				log.reportTimeError("\tSee: 	https://atgu.mgh.harvard.edu/plinkseq/input.shtml ");

				return null;
			}
			return phenotypes.toArray(new PseqPhenoTypes[phenotypes.size()]);
		}

	}

	public static class PseqProject extends Properties {
		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;
		private Logger log;
		private boolean loaded;
		private String filename, projectName, projectDirectory, resourceDirectory, phenoFile;
		private String[] vcfs;
		private PseqPhenoTypes[] phenotypes;

		public static enum PROPERTIES {
			PROJN, OUTPUT, RESOURCES, METAMETA, VARDB, INDDB, LOCDB, REFDB, SEQDB, VCF
		}

		public PseqProject(String filename, String[] vcfs, String phenoFile, String resourceDirectory, Logger log) {
			super();
			this.filename = filename;
			this.projectName = ext.rootOf(filename);
			this.projectDirectory = ext.parseDirectoryOfFile(filename);
			this.resourceDirectory = resourceDirectory;
			this.loaded = false;
			this.vcfs = vcfs;
			this.log = log;
			this.phenoFile = phenoFile;
			this.phenotypes = PseqPhenoTypes.loadPhenos(phenoFile, log);
		}

		public boolean isLoaded() {
			return loaded;
		}

		public String getProjectNameForPseq() {
			return ext.rootOf(filename, false);
		}

		public String getPhenoFile() {
			return phenoFile;
		}

		public PseqPhenoTypes[] getPhenotypes() {
			return phenotypes;
		}

		/**
		 * Convert a property into the appropriate plink seq command
		 * 
		 * @param prop
		 * @return
		 */
		public String getCommandFor(PROPERTIES prop) {
			return ("--" + prop).toLowerCase();
		}

		public void load() {
			try {
				load(Files.getAppropriateReader(filename));
				loaded = true;
			} catch (FileNotFoundException e) {
				log.reportFileNotFound(filename);
				e.printStackTrace();
			} catch (IOException e) {
				log.reportIOException(filename);
				e.printStackTrace();
			}
		}

		public String getProperty(PROPERTIES property) {
			String prop = getProperty(property + "");
			if (prop == null) {
				log.reportTimeError("Invalid pseq property " + property);
			}
			return prop;
		}

		public void setProperty(PROPERTIES property, String value) {
			setProperty(property + "", value);
		}

		public String getFilename() {
			return filename;
		}

		public String getProjectName() {
			return projectName;
		}

		public String getProjectDirectory() {
			return projectDirectory;
		}

		public String getResourceDirectory() {
			return resourceDirectory;
		}

		public String[] getVcfs() {
			return vcfs;
		}
	}

	public static void generatePlinkSeqLoc(GeneTrack geneTrack, String locFileName, Logger log) {
		log.reportTimeInfo("Attempting to parse genetrack to" + locFileName);
		if (!Files.exists(locFileName)) {
			try {
				// new File(ext.parseDirectoryOfFile(locFileName))
				new File(ext.parseDirectoryOfFile(locFileName)).mkdirs();
				GeneData[][] gDatas = geneTrack.getGenes();
				PrintWriter writer = new PrintWriter(new FileWriter(locFileName));
				writer.println("#CHR\tPOS1\tPOS2\tID");
				for (int i = 0; i < gDatas.length; i++) {
					for (int j = 0; j < gDatas[i].length; j++) {
						writer.println(Positions.getChromosomeUCSC(gDatas[i][j].getChr(), true) + "\t" + gDatas[i][j].getStart() + "\t" + gDatas[i][j].getStop() + "\t" + gDatas[i][j].getGeneName());
					}
				}
				writer.close();
			} catch (Exception e) {
				log.reportError("Error writing to " + locFileName);
				log.reportException(e);
			}
		} else {
			log.reportFileExists(locFileName);
		}

	}

	public static class PlinkSeqBurdenSummary {
		private static final String[] HEADER = new String[] { "LOCUS", "POS", "ALIAS", "NVAR", "TEST", "P", "I", "DESC" };
		private String analysis;
		private String resultsFile;
		private Hashtable<String, Integer> locMap;
		private ArrayList<PlinkSeqLocSummary> locSummaries;
		private Logger log;

		public PlinkSeqBurdenSummary(String analysis, String resultsFile, Logger log) {
			super();
			this.analysis = analysis;
			this.resultsFile = resultsFile;
			this.locMap = new Hashtable<String, Integer>();
			this.locSummaries = new ArrayList<PlinkSeqUtils.PlinkSeqLocSummary>();
			this.log = log;
		}

		public void load() {
			if (Files.exists(resultsFile)) {
				String[] header = Files.getLineContaining(resultsFile, "\t", HEADER, log);
				if (header != null) {
					int[] indices = ext.indexFactors(HEADER, header, true, false);
					String curLoc = null;
					int index = 0;
					try {
						BufferedReader reader = Files.getAppropriateReader(resultsFile);
						while (reader.ready()) {
							String[] line = reader.readLine().trim().split("[\\s]+");
							if (line.length == HEADER.length) {
								try {
									String pos = line[indices[1]];
									String alias = line[indices[2]];
									int nvar = Integer.parseInt(line[indices[3]]);
									BURDEN_Tests type = BURDEN_Tests.valueOf(line[indices[4]]);
									double p = Double.parseDouble(line[indices[5]]);
									double i = Double.parseDouble(line[indices[6]]);
									if (curLoc == null || !curLoc.equals(line[indices[0]])) {
										curLoc = line[indices[0]];
										if (locMap.containsKey(curLoc)) {
											log.reportTimeError("Multiple entries for " + curLoc);
											return;
										} else {
											locSummaries.add(new PlinkSeqLocSummary(curLoc, pos, nvar, alias));
											locMap.put(curLoc, index);
											index++;
										}
									}
									locSummaries.get(locSummaries.size() - 1).addTest(curLoc, pos, nvar, type.toString(), p, i, alias, type, log);

								} catch (NumberFormatException nfe) {
									log.reportTimeInfo("Invalid number on line " + Array.toStr(line));
									log.reportException(nfe);
								}

							}
						}
						reader.close();

					} catch (FileNotFoundException fnfe) {
						log.reportError("Error: file \"" + resultsFile + "\" not found in current directory");
					} catch (IOException ioe) {
						log.reportError("Error reading file \"" + resultsFile + "\"");
					}
				} else {
				}
			} else {
			}
		}
	}

	public static class PlinkSeqLocSummary {
		private String locus;
		private String pos;
		private int NVAR;
		private String alias;
		private PlinkSeqTestSummary[] summaries;

		public PlinkSeqLocSummary(String locus, String pos, int nVAR, String alias) {
			super();
			this.locus = locus;
			this.pos = pos;
			this.NVAR = nVAR;
			this.alias = alias;
			this.summaries = new PlinkSeqTestSummary[BURDEN_Tests.values().length];
		}

		public void addTest(String alocus, String apos, int aNVAR, String test, double P, double I, String aAlias, BURDEN_Tests type, Logger log) {
			if (!locus.equals(alocus) || pos.equals(apos) || aNVAR == NVAR || !alias.equals(aAlias)) {
				log.reportTimeError("Mismatched tests being added...");
			} else {
				switch (type) {
				case BURDEN:
					if (summaries[0] != null) {
						log.reportTimeError("Multiple " + type + " tests for locus...");
					} else {
						summaries[0] = new PlinkSeqTestSummary(test, P, I, type);
					}
					break;
				case CALPHA:
					if (summaries[1] != null) {
						log.reportTimeError("Multiple " + type + " tests for locus...");
					} else {
						summaries[1] = new PlinkSeqTestSummary(test, P, I, type);
					}
					break;
				case FW:
					if (summaries[2] != null) {
						log.reportTimeError("Multiple " + type + " tests for locus...");
					} else {
						summaries[2] = new PlinkSeqTestSummary(test, P, I, type);
					}
					break;
				case SUMSTAT:
					if (summaries[3] != null) {
						log.reportTimeError("Multiple " + type + " tests for locus...");
					} else {
						summaries[3] = new PlinkSeqTestSummary(test, P, I, type);
					}
					break;
				case UNIQ:
					if (summaries[4] != null) {
						log.reportTimeError("Multiple " + type + " tests for locus...");
					} else {
						summaries[4] = new PlinkSeqTestSummary(test, P, I, type);
					}
					break;
				case VT:
					if (summaries[5] != null) {
						log.reportTimeError("Multiple " + type + " tests for locus...");
					} else {
						summaries[5] = new PlinkSeqTestSummary(test, P, I, type);
					}
					break;
				default:
					break;

				}
			}
		}
	}

	public static class PlinkSeqTestSummary {
		private String test;
		private double P;
		private double I;
		private BURDEN_Tests type;

		public PlinkSeqTestSummary(String test, double p, double i, BURDEN_Tests type) {
			super();
			this.test = test;
			P = p;
			I = i;
			this.type = type;
		}

	}

	public static class PlinkSeqBurdenResults {
		private static final String[] HEADER = new String[] { "LOCUS", "POS", "ALIAS", "NVAR", "TEST", "P", "I", "DESC" };
		private String resultsFile;
		private int numTests;
		private double bonferoniP;
		private HashSet<String> genesPassing;
		private Logger log;

		public PlinkSeqBurdenResults(String resultsFile, Logger log) {
			this.resultsFile = resultsFile;
			this.log = log;

			this.numTests = determineNumGenes();
			this.bonferoniP = determineBonferoni(numTests);
			this.genesPassing = determineGenesPassing();

		}

		public HashSet<String> getGenesPassing() {
			return genesPassing;
		}

		public int getNumTests() {
			return numTests;
		}

		public double getBonferoniP() {
			return bonferoniP;
		}

		private int determineNumGenes() {
			log.reportTimeInfo("loading results from " + resultsFile);
			if (Files.exists(resultsFile)) {
				String[] header = Files.getLineContaining(resultsFile, "\t", HEADER, log);
				if (header != null) {
					int locusIndex = ext.indexOfStr(HEADER[0], header);
					return HashVec.loadFileToStringArray(resultsFile, true, new int[] { locusIndex }, true).length;
				} else {
					return 0;
				}
			} else {
				return 0;
			}
		}

		private HashSet<String> determineGenesPassing() {
			HashSet<String> passingGenes = new HashSet<String>();

			if (Files.exists(resultsFile)) {
				String[] header = Files.getLineContaining(resultsFile, "\t", HEADER, log);
				if (header != null) {
					int locusIndex = ext.indexOfStr(HEADER[0], header);
					int pvalIndex = ext.indexOfStr(HEADER[5], header);
					try {
						BufferedReader reader = Files.getAppropriateReader(resultsFile);
						while (reader.ready()) {
							String[] line = reader.readLine().trim().split("[\\s]+");
							if (line.length == HEADER.length) {
								try {
									double pval = Double.parseDouble(line[pvalIndex]);
									if (pval < bonferoniP) {
										passingGenes.add(line[locusIndex]);
									}
								} catch (NumberFormatException nfe) {

								}
							}
						}
						reader.close();
						return passingGenes;

					} catch (FileNotFoundException fnfe) {
						log.reportError("Error: file \"" + resultsFile + "\" not found in current directory");
						return passingGenes;
					} catch (IOException ioe) {
						log.reportError("Error reading file \"" + resultsFile + "\"");
						return passingGenes;
					}
				} else {
					return passingGenes;
				}
			} else {
				return passingGenes;
			}
		}

		private static double determineBonferoni(int numTests) {

			return numTests > 0 ? (double) 0.05 / numTests : 0;
		}
	}

}
