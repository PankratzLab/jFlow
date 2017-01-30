package org.genvisis.seq.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Properties;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.common.ext;
import org.genvisis.filesys.GeneData;
import org.genvisis.seq.analysis.PlinkSeq.BURDEN_Tests;
import org.genvisis.seq.pathway.GenomeRegions;

/**
 * Handles the project and phenotype system used by PlinkSeq
 *
 */
public class PlinkSeqUtils {

	public static final String PSEQ_PROJECT = "pseqProj_";
	public static final String GENVISIS_GENE = "_GENVISIS_GENE";
	public static final String GENVISIS = "_GENVISIS";

	public static final String GENVISIS_PATHWAY = "_GENVISIS_PATHWAY";

	/**
	 * Class for the phenotypes used in pseq, manages the file for later association testing
	 *
	 */
	public static class PseqPhenoTypes implements Serializable {
		/**
		 *
		 */
		private static final long serialVersionUID = 1L;
		private final String name, type, missingValue, description;
		private final Logger log;

		public PseqPhenoTypes(String[] phenoLine, Logger log) {
			this(phenoLine[0].replace("##", ""), phenoLine[1], phenoLine[2], phenoLine[3], log);
		}

		public PseqPhenoTypes(String name, String type, String missingValue, String description,
													Logger log) {
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
		 * @param phenoFile pseq formated phenotype file (See:
		 *        https://atgu.mgh.harvard.edu/plinkseq/input.shtml )
		 * @param log
		 * @return
		 */
		public static PseqPhenoTypes[] loadPhenos(String phenoFile, Logger log) {
			ArrayList<PseqPhenoTypes> phenotypes = new ArrayList<PlinkSeqUtils.PseqPhenoTypes>();
			if (phenoFile == null) {
				log.reportError("No phenotype file was provided, continuing carefully");
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
							log.reportError("Phenotype lines starting with ## must have 4 comma-delimited entries");
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
				log.reportError("Did not find any valid phenotype headers listed in " + phenoFile);
				log.reportError("\tSee: 	https://atgu.mgh.harvard.edu/plinkseq/input.shtml ");

				return null;
			}
			return phenotypes.toArray(new PseqPhenoTypes[phenotypes.size()]);
		}

	}

	public static class PseqProject extends Properties implements Serializable {
		/**
		 *
		 */
		private static final long serialVersionUID = 1L;
		private final Logger log;
		private boolean loaded;
		private final String filename, projectName, projectDirectory, resourceDirectory, phenoFile;
		private final String[] vcfs;
		private final PseqPhenoTypes[] phenotypes;

		public static enum PROPERTIES {
																		PROJN, OUTPUT, RESOURCES, METAMETA, VARDB, INDDB, LOCDB, REFDB, SEQDB, VCF
		}

		public PseqProject(	String filename, String[] vcfs, String phenoFile, String resourceDirectory,
												Logger log) {
			super();
			this.filename = filename;
			projectName = ext.rootOf(filename);
			projectDirectory = ext.parseDirectoryOfFile(filename);
			this.resourceDirectory = resourceDirectory;
			loaded = false;
			this.vcfs = vcfs;
			this.log = log;
			this.phenoFile = phenoFile;
			phenotypes = PseqPhenoTypes.loadPhenos(phenoFile, log);
		}

		public Logger getLog() {
			return log;
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
				log.reportError("Invalid pseq property " + property);
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

	public static void generatePlinkSeqLoc(GenomeRegions gRegions, String locFileName, Logger log) {
		log.reportTimeInfo("Attempting to parse genome regions to" + locFileName);
		// if (!Files.exists(locFileName)) {
		try {
			// new File(ext.parseDirectoryOfFile(locFileName))
			new File(ext.parseDirectoryOfFile(locFileName)).mkdirs();
			GeneData[][] gDatas = gRegions.getGeneTrack().getGenes();
			PrintWriter writer = new PrintWriter(new FileWriter(locFileName));
			writer.println("#CHR\tPOS1\tPOS2\tID");
			for (GeneData[] gData : gDatas) {
				for (GeneData element : gData) {
					writer.println(Positions.getChromosomeUCSC(element.getChr(), true)	+ "\t"
													+ element.getStart() + "\t" + element.getStop() + "\t"
													+ element.getGeneName() + GENVISIS_GENE);
				}
			}
			// }
			// Pathway[] ways = gRegions.getPathways().getPathways();
			// for (int i = 0; i < ways.length; i++) {
			// GeneData[] pathGenes = ways[i].getLoci();
			// for (int j = 0; j < pathGenes.length; j++) {
			// writer.println(Positions.getChromosomeUCSC(pathGenes[j].getChr(), true) + "\t" +
			// pathGenes[j].getStart() + "\t" + pathGenes[j].getStop() + "\t" + ways[i].getPathwayName() +
			// GENVISIS_PATHWAY);
			// }
			// }
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + locFileName);
			log.reportException(e);
		}
		// } else {
		// log.reportFileExists(locFileName);
		// }

	}

	public static class PlinkSeqBurdenSummary implements Serializable {
		/**
		 *
		 */
		private static final long serialVersionUID = 1L;
		private static final String[] HEADER = new String[] {	"LOCUS", "POS", "ALIAS", "NVAR", "TEST",
																													"P", "I", "DESC"};
		public static final double[] I_THRESHOLDS = new double[] {0.5, 0.1, 0.01, 0.001};
		private final String analysis;
		private final String resultsFile;
		private final Hashtable<String, Integer> locMap;
		private final ArrayList<PlinkSeqLocSummary> locSummaries;
		private final int[] numTotalTests;
		private final int[][] numTestsAtI;
		private final Logger log;

		public PlinkSeqBurdenSummary(String analysis, String resultsFile, Logger log) {
			super();
			this.analysis = analysis;
			this.resultsFile = resultsFile;
			locMap = new Hashtable<String, Integer>();
			locSummaries = new ArrayList<PlinkSeqUtils.PlinkSeqLocSummary>();
			numTotalTests = new int[BURDEN_Tests.values().length];
			numTestsAtI = new int[BURDEN_Tests.values().length][I_THRESHOLDS.length];
			this.log = log;
		}

		public String getAnalysis() {
			return analysis;
		}

		public boolean hasSummaryFor(String locus) {
			return locMap.containsKey(locus);
		}

		public PlinkSeqLocSummary getPlinkSeqLocSummaryFor(String locus) {
			return locSummaries.get(locMap.get(locus));
		}

		public void correctPvalues() {
			// int bonfFullNumTests = locSummaries.size();
			for (int i = 0; i < BURDEN_Tests.values().length; i++) {
				for (int j = 0; j < locSummaries.size(); j++) {
					PlinkSeqTestSummary curSummary = locSummaries.get(j).getSummaries()[i];
					if (curSummary != null) {
						numTotalTests[i]++;
						for (int k = 0; k < I_THRESHOLDS.length; k++) {
							if (curSummary.getI() < I_THRESHOLDS[k]) {
								numTestsAtI[i][k]++;
							}
						}
					}
				}
				for (int j = 0; j < locSummaries.size(); j++) {
					PlinkSeqTestSummary curSummary = locSummaries.get(j).getSummaries()[i];
					if (curSummary != null) {
						curSummary.setBonfFull(numTotalTests[i]);
						curSummary.setbonfsI(numTestsAtI[i]);
					}
				}
				log.reportTimeInfo("Correction summary for " + analysis + " in " + resultsFile + ":");
				log.reportTimeInfo("Full number of tests performed: " + numTotalTests[i]);
				for (int j = 0; j < I_THRESHOLDS.length; j++) {
					log.reportTimeInfo("Ithreshold " + I_THRESHOLDS[j] + ": " + numTestsAtI[i][j]);
				}
			}

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
						reader.readLine();
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
									String desc = line[indices[7]];
									if (curLoc == null || !curLoc.equals(line[indices[0]])) {
										curLoc = line[indices[0]];
										if (locMap.containsKey(curLoc)) {
											log.reportError("Multiple entries for " + curLoc);
											return;
										} else {
											locSummaries.add(new PlinkSeqLocSummary(curLoc, pos, nvar, alias));
											locMap.put(curLoc, index);
											index++;
										}
									}
									locSummaries.get(locSummaries.size() - 1).addTest(curLoc, pos, nvar,
																																		type.toString(), p, i, alias,
																																		desc, type, log);

								} catch (NumberFormatException nfe) {
									log.reportTimeInfo("Invalid number on line " + ArrayUtils.toStr(line));
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

	public static class PlinkSeqLocSummary implements Serializable {
		/**
		 *
		 */
		private static final long serialVersionUID = 1L;
		private String locus;
		private String pos;
		private int NVAR;
		private String alias;
		private PlinkSeqTestSummary[] summaries;

		public PlinkSeqLocSummary(String locus, String pos, int nVAR, String alias) {
			super();
			this.locus = locus;
			this.pos = pos;
			NVAR = nVAR;
			this.alias = alias;
			summaries = new PlinkSeqTestSummary[BURDEN_Tests.values().length];
		}

		public String getLocus() {
			return locus;
		}

		public void setLocus(String locus) {
			this.locus = locus;
		}

		public String getPos() {
			return pos;
		}

		public void setPos(String pos) {
			this.pos = pos;
		}

		public int getNVAR() {
			return NVAR;
		}

		public void setNVAR(int nVAR) {
			NVAR = nVAR;
		}

		public String getAlias() {
			return alias;
		}

		public void setAlias(String alias) {
			this.alias = alias;
		}

		public PlinkSeqTestSummary[] getSummaries() {
			return summaries;
		}

		public void setSummaries(PlinkSeqTestSummary[] summaries) {
			this.summaries = summaries;
		}

		public void addTest(String alocus, String apos, int aNVAR, String test, double P, double I,
												String aAlias, String desc, BURDEN_Tests type, Logger log) {
			if (!locus.equals(alocus) || !pos.equals(apos) || aNVAR != NVAR || !alias.equals(aAlias)) {
				log.reportError("Mismatched tests being added...");
				log.reportError(locus + "->" + alocus);
				log.reportError(pos + "->" + apos);
				log.reportError(aNVAR + "->" + NVAR);
				log.reportError(alias + "->" + alias);

			} else {
				PlinkSeqTestSummary plinkSeqTestSummary = new PlinkSeqTestSummary(test, P, I, desc, type);
				// BURDEN, UNIQ, VT, FW, CALPHA, SUMSTAT,FRQWGT

				for (int i = 0; i < BURDEN_Tests.values().length; i++) {
					if (type == BURDEN_Tests.values()[i]) {
						if (summaries[i] != null) {
							log.reportError("Multiple " + type + " tests for locus...");
						} else {
							summaries[i] = plinkSeqTestSummary;
						}
					}
				}
			}
		}
	}

	public static class PlinkSeqTestSummary implements Serializable {
		/**
		 *
		 */
		private static final long serialVersionUID = 1L;
		public static final String[] SUMMARY = new String[] {"P_VAL", "I"};
		private String test;
		private double P;
		private double I;
		private double bonfFull;
		private double[] bonfsI;

		private String desc;
		private BURDEN_Tests type;

		public PlinkSeqTestSummary(String test, double p, double i, String desc, BURDEN_Tests type) {
			super();
			this.test = test;
			P = p;
			I = i;
			this.desc = desc;
			this.type = type;
		}

		public String getTest() {
			return test;
		}

		public void setBonfFull(int numTests) {
			if (numTests > 0) {
				bonfFull = Math.max(P * numTests, 1);
			} else {
				bonfFull = P;
			}
		}

		public void setbonfsI(int[] tests) {
			bonfsI = new double[tests.length];
			for (int i = 0; i < tests.length; i++) {
				if (tests[i] > 0) {
					bonfsI[i] = Math.max(P * tests[i], 1);
				} else {
					bonfsI[i] = P;
				}
			}
		}

		public double getBonfFull() {
			return bonfFull;
		}

		public double[] getBonfsI() {
			return bonfsI;
		}

		public void setTest(String test) {
			this.test = test;
		}

		public double getP() {
			return P;
		}

		public void setP(double p) {
			P = p;
		}

		public double getI() {
			return I;
		}

		public void setI(double i) {
			I = i;
		}

		public String getDesc() {
			return desc;
		}

		public void setDesc(String desc) {
			this.desc = desc;
		}

		public BURDEN_Tests getType() {
			return type;
		}

		public void setType(BURDEN_Tests type) {
			this.type = type;
		}

	}

	public static class PlinkSeqBurdenResults {
		private static final String[] HEADER = new String[] {	"LOCUS", "POS", "ALIAS", "NVAR", "TEST",
																													"P", "I", "DESC"};
		private final String resultsFile;
		private final int numTests;
		private final double bonferoniP;
		private final HashSet<String> genesPassing;
		private final Logger log;

		public PlinkSeqBurdenResults(String resultsFile, Logger log) {
			this.resultsFile = resultsFile;
			this.log = log;

			numTests = determineNumGenes();
			bonferoniP = determineBonferoni(numTests);
			genesPassing = determineGenesPassing();

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
					return HashVec.loadFileToStringArray(	resultsFile, true, new int[] {locusIndex},
																								true).length;
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

			return numTests > 0 ? 0.05 / numTests : 0;
		}
	}

}
