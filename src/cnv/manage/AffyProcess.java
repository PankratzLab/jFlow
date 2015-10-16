package cnv.manage;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.concurrent.Callable;

import common.Array;
import common.Files;
import common.Logger;
import common.WorkerTrain;
import common.ext;
import common.WorkerTrain.Producer;
import cnv.filesys.Project;
import cnv.filesys.Project.ARRAY;
import cnv.filesys.Sample;

/**
 * @author lane0212
 * 
 */
public class AffyProcess {
	private static final String[] AFFY_CHP_HEADER = new String[] { "Probe Set ID", "Call Codes", "Forward Strand Base Calls", "Confidence", "Signal A", "Signal B" };
	private static final String[] AFFY_CN_CHP_HEADER = new String[] { "ProbeSet", "Log2Ratio" };

	private static final String CN_C5_PATTERN = "CN5.CNCHP";
	private static final String[][] REPLACES = new String[][] { { "-v2", "" }, { ".birdseed", "" } };
	private Project proj;
	private String[] chpFiles;
	private String[] cn5chpFiles;
	private String[] combinedOutputFiles;
	private String delimiter;
	private boolean valid;
	private Logger log;

	public AffyProcess(Project proj, String[] chpFiles, String delimiter, Logger log) {
		super();
		this.proj = proj;
		this.chpFiles = chpFiles;
		this.delimiter = delimiter;
		this.log = log;
	}

	public void combineFirst() {
		if (valid) {
			combineChp(proj, chpFiles[0], cn5chpFiles[0], combinedOutputFiles[0], delimiter, log);
		}
	}

	public void combineAll(int numthreads) {
		if (valid) {
			CombineProducer producer = new CombineProducer(proj, chpFiles, cn5chpFiles, combinedOutputFiles, delimiter, log);
			WorkerTrain<Boolean> train = new WorkerTrain<Boolean>(producer, numthreads, numthreads, log);
			while (train.hasNext()) {
				train.next();
			}
		}
	}

	public String[] getCombinedOutputFiles() {
		return combinedOutputFiles;
	}

	public void matchCn() {
		if (proj.getArrayType() != ARRAY.AFFY_GW6_CN) {
			log.reportTimeError("Array type must be " + ARRAY.AFFY_GW6_CN + ", halting");
			valid = false;

		} else if (chpFiles == null || chpFiles.length == 0) {
			log.reportTimeError("No files detected to parse");
			valid = false;
		} else if (!verifyAffy(chpFiles, "CHP", AFFY_CHP_HEADER, log)) {
			valid = false;
		} else {
			String dir = ext.parseDirectoryOfFile(chpFiles[0]);
			String[] allFiles = Files.list(dir, null, false);

			log.reportTimeInfo("Attempting to find matched " + CN_C5_PATTERN + " files for " + chpFiles.length + " chp format files in " + dir);
			this.cn5chpFiles = new String[chpFiles.length];
			this.combinedOutputFiles = new String[chpFiles.length];
			boolean foundAll = true;
			for (int i = 0; i < chpFiles.length; i++) {
				if (chpFiles[i].contains(CN_C5_PATTERN)) {
					log.reportTimeWarning(chpFiles[i] + " had " + CN_C5_PATTERN + " in the filename, this pattern is used to match chp files to CN5.CNCHP files and is likely an internal error ");
				}
				if (Files.getLineContaining(chpFiles[i], "\t", AFFY_CHP_HEADER, log) == null) {
					log.reportTimeError("Currently chp files must have at least the following columns " + Array.toStr(AFFY_CHP_HEADER));
					valid = false;
				} else {
					String toMatch = ext.replaceAllWith(ext.rootOf(chpFiles[i]), REPLACES);
					String match = null;
					String tmpExt = "." + proj.getArrayType() + ".tmp.gz";
					proj.SOURCE_FILENAME_EXTENSION.setValue(tmpExt);
					for (int j = 0; j < allFiles.length; j++) {
						if (allFiles[j].startsWith(toMatch) && allFiles[j].contains(CN_C5_PATTERN) && !allFiles[j].endsWith(tmpExt)) {
							match = allFiles[j];
							break;
						}
					}
					if (match == null) {
						foundAll = false;
						log.reportTimeError("Could not find matching " + CN_C5_PATTERN + " file for CHP file " + chpFiles[i]);
					} else {
						cn5chpFiles[i] = dir + match;
						combinedOutputFiles[i] = dir + toMatch + tmpExt;
					}
				}

			}
			valid = foundAll;
			if (valid) {
				valid = verifyAffy(cn5chpFiles, CN_C5_PATTERN, AFFY_CN_CHP_HEADER, log);
				if (valid) {
					log.reportTimeInfo("Matched all " + chpFiles.length + " CHP files to corresponding " + CN_C5_PATTERN + " files");
				}
			}
		}
	}

	private static class CombineProducer implements Producer<Boolean> {
		private Project proj;
		private String[] chpFiles;
		private String[] cn5chpFiles;
		private String[] combinedOutputFiles;
		private String delimiter;
		private Logger log;
		private int index;

		public CombineProducer(Project proj, String[] chpFiles, String[] cn5chpFiles, String[] combinedOutputFiles, String delimiter, Logger log) {
			super();
			this.proj = proj;
			this.chpFiles = chpFiles;
			this.cn5chpFiles = cn5chpFiles;
			this.combinedOutputFiles = combinedOutputFiles;
			this.delimiter = delimiter;
			this.log = log;
			this.index = 0;
		}

		@Override
		public boolean hasNext() {
			// TODO Auto-generated method stub
			return index > chpFiles.length;
		}

		@Override
		public Callable<Boolean> next() {
			final String currentChp = chpFiles[index];
			final String currentCn5chp = cn5chpFiles[index];
			final String currentOutput = combinedOutputFiles[index];
			final Project currentProj = proj;
			Callable<Boolean> callable = new Callable<Boolean>() {

				@Override
				public Boolean call() throws Exception {
					// TODO Auto-generated method stub
					return combineChp(currentProj, currentChp, currentCn5chp, currentOutput, delimiter, log);
				}
			};
			index++;

			return callable;
		}

		@Override
		public void remove() {

		}

		@Override
		public void shutdown() {

		}

	}

	private static boolean combineChp(final Project proj, final String chp, final String cnC5Chp, final String output, final String delimiter, final Logger log) {
		String[] headerChp = Files.getLineContaining(chp, delimiter, AFFY_CHP_HEADER, log);
		String[] headerCNChp = Files.getLineContaining(cnC5Chp, delimiter, AFFY_CN_CHP_HEADER, log);
		if (headerChp != null && headerCNChp != null) {
			try {
				PrintWriter writer = Files.getAppropriateWriter(output);

				BufferedReader reader = Files.getAppropriateReader(chp);
				String[] line;
				do {
					line = reader.readLine().trim().split("\t", -1);
				} while (reader.ready() && ext.indexFactors(SourceFileParser.SNP_HEADER_OPTIONS, line, false, true, false, false)[0] == -1);

				writer.println(Array.toStr(AFFY_CHP_HEADER));
				int[] indices = ext.indexFactors(AFFY_CHP_HEADER, line, true, false);
				int numSnps = 0;
				while (reader.ready()) {
					line = reader.readLine().trim().split(delimiter);
					writer.println(Array.toStr(Array.subArray(line, indices)));
					numSnps++;
				}

				reader.close();
				reader = Files.getAppropriateReader(cnC5Chp);
				do {
					line = reader.readLine().trim().split("\t", -1);
				} while (reader.ready() && ext.indexFactors(SourceFileParser.SNP_HEADER_OPTIONS, line, false, true, false, false)[0] == -1);

				indices = ext.indexFactors(AFFY_CN_CHP_HEADER, line, true, false);
				int numCN = 0;
				while (reader.ready()) {
					line = reader.readLine().trim().split(delimiter);

					if (proj.getArrayType().isCNOnly(line[indices[0]])) {
						numCN++;
						double log2Ratio = Double.NaN;
						try {
							log2Ratio = Double.parseDouble(line[indices[1]]);
							log2Ratio = Math.pow(log2Ratio, 2);
						} catch (NumberFormatException nfe) {
							log.reportTimeError("Could not parse Log2Ratio on line " + Array.toStr(line));
						}
						for (int i = 0; i < AFFY_CHP_HEADER.length; i++) {
							switch (i) {
							case 0:
								writer.print(line[indices[0]]);// SNP ID
								break;
							case 1:
								writer.print("\t" + Sample.ALT_NULLS[0]);// Genotype
								break;
							case 2:
								writer.print("\t" + Sample.ALT_NULLS[0]);// Genotype
								break;
							case 3:
								writer.print("\t" + 0);// GC
								break;
							case 4:
								writer.print("\t" + log2Ratio);//
								break;
							case 5:
								writer.print("\t" + log2Ratio);//
								break;
							default:
								log.reportTimeError("Internal error, mismatched CN reporting on index " + i);
								break;
							}
						}
						writer.println();
					}
				}
				reader.close();
				writer.close();
				log.reportTimeInfo("Combined " + numSnps + " SNP probesets and " + numCN + " CN probes to " + output);
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

			return true;
		} else {
			log.reportTimeError("Could not find appropriate headers");
			return false;
		}
	}

	public static boolean verifyAffy(String[] files, String type, String[] header, Logger log) {
		boolean valid = true;
		for (int i = 0; i < files.length; i++) {
			if (Files.getLineContaining(files[i], "\t", header, log) == null) {
				log.reportTimeError("Currently " + type + " files must have at least the following columns " + Array.toStr(header));
				log.reportTimeError(files[i] + " did not contain a line with  the required columns");

				return false;
			}
		}
		return valid;
	}
}
