package org.genvisis.sra;

import org.genvisis.seq.SeqVariables.ASSAY_TYPE;
import org.genvisis.seq.SeqVariables.ASSEMBLY_NAME;
import org.genvisis.seq.SeqVariables.PLATFORM;

/**
 * Stores info about a sample from the SRA archive
 *
 */
public class SRASample {

	private String runS;
	private String submittedSampleID;
	private ASSEMBLY_NAME aName;
	private ASSAY_TYPE aType;
	private PLATFORM platform;

	/**
	 * @param runS
	 *            the SRA ID
	 * @param submittedSampleID
	 *            study ID
	 * @param aName
	 * @param aType
	 * 
	 * @param platform
	 *            sequencing platform
	 */
	public SRASample(String runS, String submittedSampleID, ASSEMBLY_NAME aName, ASSAY_TYPE aType, PLATFORM platform) {
		super();
		this.runS = runS;
		this.submittedSampleID = submittedSampleID;
		this.aName = aName;
		this.aType = aType;
		this.platform = platform;
	}

	public String getRunS() {
		return runS;
	}

	public String getSubmittedSampleID() {
		return submittedSampleID;
	}

	public ASSEMBLY_NAME getaName() {
		return aName;
	}

	public PLATFORM getPlatform() {
		return platform;
	}

	public ASSAY_TYPE getaType() {
		return aType;
	}

	@Override
	public String toString() {
		return "SRASample [runS=" + runS + ", submittedSampleID=" + submittedSampleID + ", aName=" + aName + ", aType="
				+ aType + ", platform=" + platform + "]";
	}

}
