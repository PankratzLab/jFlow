package org.genvisis.seq;

import org.genvisis.seq.SeqVariables.ASSAY_TYPE;
import org.genvisis.seq.SeqVariables.ASSEMBLY_NAME;
import org.genvisis.seq.SeqVariables.PLATFORM;

/**
 * Stores basic information about a Next Generation sequencing sample
 *
 */
public class NGSSample {

	private ASSEMBLY_NAME aName;
	private ASSAY_TYPE aType;
	private PLATFORM platform;

	/**
	 * @param aName
	 *            see {@link ASSEMBLY_NAME}
	 * @param aType
	 *            see {@link ASSAY_TYPE}
	 * @param platform
	 *            see {@link PLATFORM}
	 */
	public NGSSample(ASSEMBLY_NAME aName, ASSAY_TYPE aType, PLATFORM platform) {
		super();
		this.aName = aName;
		this.aType = aType;
		this.platform = platform;
	}

	public ASSEMBLY_NAME getaName() {
		return aName;
	}

	public ASSAY_TYPE getaType() {
		return aType;
	}

	public PLATFORM getPlatform() {
		return platform;
	}

}
