package seq.analysis;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Properties;

import common.Files;
import common.Logger;
import common.ext;

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

}
