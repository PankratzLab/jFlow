package org.genvisis.sra;

import org.genvisis.seq.NGSSample;
import org.genvisis.seq.SeqVariables.ASSAY_TYPE;
import org.genvisis.seq.SeqVariables.ASSEMBLY_NAME;
import org.genvisis.seq.SeqVariables.PLATFORM;

/**
 * Stores info about a sample from the SRA archive
 *
 */
public class SRASample extends NGSSample {

	private final String runS;
	private final String submittedSampleID;
	private String sraFile;

	/**
	 * @param runS the SRA ID
	 * @param submittedSampleID study ID
	 * @param aName
	 * @param aType
	 *
	 * @param platform sequencing platform
	 */
	public SRASample(	String runS, String submittedSampleID, ASSEMBLY_NAME aName, ASSAY_TYPE aType,
										PLATFORM platform) {
		super(aName, aType, platform);
		this.runS = runS;
		this.submittedSampleID = submittedSampleID;
		sraFile = null;
	}



	public void setSraFile(String sraFile) {
		this.sraFile = sraFile;
	}



	public String getSraFile() {
		return sraFile;
	}



	public String getRunS() {
		return runS;
	}

	public String getSubmittedSampleID() {
		return submittedSampleID;
	}

	@Override
	public String toString() {
		return "SRASample [runS="	+ runS + ", submittedSampleID=" + submittedSampleID + ", aName="
						+ getaName() + ", aType=" + getaType() + ", platform=" + getPlatform() + "]";
	}

}
