package cnv.filesys;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;

import cnv.var.SampleData;
import common.Files;
import common.Logger;
import common.ext;

public class Pedigree {

	/**
	 * Typically corresponding to each DNA in the project and null if not
	 * represented;
	 */
	private PedigreeEntry[] pedigreeEntries;

	private Pedigree(PedigreeEntry[] pedigreeEntries, Logger log) {
		super();
		this.pedigreeEntries = pedigreeEntries;
	}

	public PedigreeEntry[] getPedigreeEntries() {
		return pedigreeEntries;
	}

	public static Pedigree loadPedigree(Project proj, String ped) {
		Pedigree pedigree = null;
		SampleData sampleData = proj.getSampleData(0, false);
		String[] samples = proj.getSamples();
		int numAdded = 0;
		

		try {
			BufferedReader reader = Files.getAppropriateReader(ped);
			String delim = ext.determineDelimiter(Files.getFirstNLinesOfFile(
					ped, 1, proj.getLog())[0]);
			int lineNum = 0;
			PedigreeEntry[] pedigreeEntries = new PedigreeEntry[samples.length];
			while (reader.ready()) {
				lineNum++;
				String[] line = reader.readLine().trim().split(delim);
				if (line.length < 7) {
					proj.getLog()
							.reportTimeError(
									"Detected that pedigree file "
											+ ped
											+ " exists, but starting at line "
											+ (lineNum)
											+ (line.length < 3 ? ""
													: " (individual " + line[0]
															+ "-" + line[1]
															+ ")")
											+ " there are only "
											+ line.length
											+ " columns in pedigree file '"
											+ proj.PEDIGREE_FILENAME.getValue()
											+ "'.\n"
											+ "  Pedigree files require 7 columns with no header: FID IID FA MO SEX PHENO DNA\n"
											+ "  where DNA is the sample name associated with the genotypic data (see the "
											+ proj.SAMPLE_DIRECTORY.getValue(
													false, true)
											+ " directory for examples)");
					reader.close();
					return pedigree;
				} else {
					String FID = line[0];

					String faFidIid = FID + line[2];
					String moFidIid = FID + line[3];
					String DNA = line[6];
					int iDNAIndex = getSampleIndex(DNA, sampleData, samples);
					int faDNAIndex = getSampleIndex(faFidIid, sampleData,
							samples);
					int moDNAIndex = getSampleIndex(moFidIid, sampleData,
							samples);

					PedigreeEntry pedigreeEntry = new PedigreeEntry(FID,
							line[1], line[2], line[3], line[4], line[4], DNA,
							iDNAIndex, faDNAIndex, moDNAIndex);
					if (iDNAIndex >= 0) {
						pedigreeEntries[iDNAIndex] = pedigreeEntry;
						numAdded++;
					}
				}

			}
		} catch (FileNotFoundException e) {
			proj.getLog().reportFileNotFound(ped);

		} catch (IOException e) {
			proj.getLog().reportException(e);
			e.printStackTrace();
		}
		proj.getLog()
				.reportTimeInfo("loaded " + numAdded + " pedigree entries");
		return pedigree;
	}

	private static int getSampleIndex(String sample, SampleData sampleData,
			String[] projectSamples) {
		int sampleIndex = -1;
		if (sampleData.lookup(sample) != null) {
			sampleIndex = ext.indexOfStr(sample, projectSamples);
		}
		return sampleIndex;
	}

	public static class PedigreeEntry {
		private String FID;
		private String IID;
		private String FA;
		private String MO;
		private String SEX;
		private String PHENO;
		private String iDNA;
		private int iDNAIndex;
		private int faDNAIndex;
		private int moDNAIndex;

		public PedigreeEntry(String fID, String iID, String fA, String mO,
				String sEX, String pHENO, String iDNA, int iDNAIndex,
				int faDNAIndex, int moDNAIndex) {
			super();
			this.FID = fID;
			this.IID = iID;
			this.FA = fA;
			this.MO = mO;
			this.SEX = sEX;
			this.PHENO = pHENO;
			this.iDNA = iDNA;
			this.iDNAIndex = iDNAIndex;
			this.faDNAIndex = faDNAIndex;
			this.moDNAIndex = moDNAIndex;

		}

		public String getFID() {
			return FID;
		}

		public String getIID() {
			return IID;
		}

		public String getFA() {
			return FA;
		}

		public String getMO() {
			return MO;
		}

		public String getSEX() {
			return SEX;
		}

		public String getPHENO() {
			return PHENO;
		}

		public String getiDNA() {
			return iDNA;
		}

		public int getiDNAIndex() {
			return iDNAIndex;
		}

		public int getFaDNAIndex() {
			return faDNAIndex;
		}

		public int getMoDNAIndex() {
			return moDNAIndex;
		}

	}

}
