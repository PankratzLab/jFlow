package org.genvisis.seq.analysis.genage;

import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.seq.SeqVariables.ASSAY_TYPE;
import org.genvisis.seq.SeqVariables.ASSEMBLY_NAME;

/**
 * Going to be the pipeline of execution for a single input bam file
 *
 */
public class Pipeline {

	public static void pipeline(String inputBam, String referenceGenome, String captureBed, ASSAY_TYPE type,
			ASSEMBLY_NAME aName, Logger log) {
		if (!Files.exists(inputBam)) {
			throw new IllegalArgumentException("Bam file " + inputBam + " must exist");
		}

		if (!Files.exists(referenceGenome)) {
			throw new IllegalArgumentException("Reference Genome " + inputBam + " must exist");
		} else {
			log.reportTimeWarning("Assuming " + referenceGenome + " matches assembly type " + aName);
		}
		if (type == ASSAY_TYPE.WXS && (!Files.exists(captureBed))) {
			throw new IllegalArgumentException(captureBed + " must exist");
		}
		

	}

}
