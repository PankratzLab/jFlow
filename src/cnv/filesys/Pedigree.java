package cnv.filesys;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;

import cnv.qc.MendelErrors;
import cnv.qc.MendelErrors.MendelErrorCheck;
import cnv.var.SampleData;
import common.Files;
import common.Logger;
import common.ext;

public class Pedigree {

	/**
	 * Typically corresponding to each DNA in the project and null if not represented;
	 */
	private Project proj;
	private PedigreeEntry[] pedigreeEntries;
	private boolean projectOrder;

	private Pedigree(Project proj, PedigreeEntry[] pedigreeEntries, boolean projOrder, Logger log) {
		super();
		this.proj = proj;
		this.pedigreeEntries = pedigreeEntries;
		this.projectOrder = projOrder;
	}

	public PedigreeEntry[] getPedigreeEntries() {
		return pedigreeEntries;
	}

	/**
	 * @param markerData
	 *            evaluate mendel errors against this marker
	 * @param samplesToCheck
	 *            can be null, if not null -> only these samples will be used in the error checks
	 * @param sex
	 *            can be null, but used for chr23 checks on male samples
	 * @param clusterFilters
	 *            can be null
	 * @param gcThreshold
	 * @return
	 */
	public MendelErrorCheck[] checkMendelErrors(MarkerData markerData, boolean[] samplesToCheck, String[] sex, ClusterFilterCollection clusterFilters, float gcThreshold) {
		MendelErrorCheck[] mendelErrorChecks = new MendelErrorCheck[proj.getSamples().length];
		byte[] genotypes = markerData.getAbGenotypesAfterFilters(clusterFilters, markerData.getMarkerName(), gcThreshold);

		if (!projectOrder) {
			proj.getLog().reportTimeError("Pedigree file must be in project order, internal error");
			return null;
		} else {
			
			for (int i = 0; i < pedigreeEntries.length; i++) {
				if (pedigreeEntries[i] != null && (samplesToCheck == null || samplesToCheck[i]) && pedigreeEntries[i].getiDNAIndex() >= 0) {
					int sampleIndex = pedigreeEntries[i].getiDNAIndex();
					int faDNAIndex = pedigreeEntries[i].getFaDNAIndex();
					int moDNAIndex = pedigreeEntries[i].getMoDNAIndex();

					byte faGenotype = -1;
					if (faDNAIndex >= 0 && (samplesToCheck == null || samplesToCheck[faDNAIndex])) {
						faGenotype = genotypes[faDNAIndex];
					}
					byte moGenotype = -1;
					if (moDNAIndex >= 0 && (samplesToCheck == null || samplesToCheck[moDNAIndex])) {
						moGenotype = genotypes[moDNAIndex];
					}
					int sampleSex = -1;
					try {
						if (sex != null) {
							sampleSex = Integer.parseInt(sex[pedigreeEntries[i].getiDNAIndex()]);

						}
					} catch (NumberFormatException nfe) {

					}
					//System.out.println(faGenotype+"\t"+moGenotype);
					MendelErrors mendelErrors = new MendelErrors(markerData.getChr(), sampleSex, genotypes[sampleIndex], faGenotype, moGenotype);
					mendelErrorChecks[i] = mendelErrors.checkMendelError();
				} else {
					mendelErrorChecks[i] = new MendelErrors(markerData.getChr(), -1, (byte) -1, (byte) -1, (byte) -1).checkMendelError();
				}
			}
		}
		return mendelErrorChecks;
	}

	public static Pedigree loadPedigree(Project proj, String ped) {
		Pedigree pedigree = null;
		SampleData sampleData = proj.getSampleData(0, false);
		String[] samples = proj.getSamples();
		int numAdded = 0;
		int numTrio =0;
		int numChildOffspring =0;

		try {
			BufferedReader reader = Files.getAppropriateReader(ped);
			String delim = ext.determineDelimiter(Files.getFirstNLinesOfFile(ped, 1, proj.getLog())[0]);
			int lineNum = 0;
			PedigreeEntry[] pedigreeEntries = new PedigreeEntry[samples.length];
			while (reader.ready()) {
				lineNum++;
				String[] line = reader.readLine().trim().split(delim);
				if (line.length < 7) {
					proj.getLog().reportTimeError("Detected that pedigree file " + ped + " exists, but starting at line " + (lineNum) + (line.length < 3 ? "" : " (individual " + line[0] + "-" + line[1] + ")") + " there are only " + line.length + " columns in pedigree file '" + proj.PEDIGREE_FILENAME.getValue() + "'.\n" + "  Pedigree files require 7 columns with no header: FID IID FA MO SEX PHENO DNA\n" + "  where DNA is the sample name associated with the genotypic data (see the " + proj.SAMPLE_DIRECTORY.getValue(false, true) + " directory for examples)");
					reader.close();
					return pedigree;
				} else {
					String FID = line[0];

					String faFidIid = FID +"\t"+ line[2];
					String moFidIid = FID +"\t"+ line[3];
				
					String DNA = line[6];
					int iDNAIndex = getSampleIndex(DNA, sampleData, samples);
					int faDNAIndex = getSampleIndex(faFidIid, sampleData, samples);
					int moDNAIndex = getSampleIndex(moFidIid, sampleData, samples);

					PedigreeEntry pedigreeEntry = new PedigreeEntry(FID, line[1], line[2], line[3], line[4], line[4], DNA, iDNAIndex, faDNAIndex, moDNAIndex);
					if (iDNAIndex >= 0) {
						if(faDNAIndex>=0){
							numChildOffspring++;
						}if(moDNAIndex>=0){
							numChildOffspring++;
						}
						if(moDNAIndex>=0&&faDNAIndex>=0){
							numTrio++;
						}
						pedigreeEntries[iDNAIndex] = pedigreeEntry;
						numAdded++;
					}
				}
			}
			pedigree = new Pedigree(proj, pedigreeEntries, true, proj.getLog());
		} catch (FileNotFoundException e) {
			proj.getLog().reportFileNotFound(ped);

		} catch (IOException e) {
			proj.getLog().reportException(e);
			e.printStackTrace();
		}
		proj.getLog().reportTimeInfo("loaded " + numAdded + " pedigree entries; Parent Offspring pairs: " + numChildOffspring + "; Trios: " + numTrio);
		return pedigree;
	}

	private static int getSampleIndex(String sample, SampleData sampleData, String[] projectSamples) {
		int sampleIndex = -1;
		if (sampleData.lookup(sample) != null) {
			sampleIndex = ext.indexOfStr(sampleData.lookup(sample)[0], projectSamples);
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

		public PedigreeEntry(String fID, String iID, String fA, String mO, String sEX, String pHENO, String iDNA, int iDNAIndex, int faDNAIndex, int moDNAIndex) {
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
