package org.genvisis.sra;

/**
 * Common variables from the SRA
 *
 */
public class SRAVariables {

	private SRAVariables() {

	}

	/**
	 * Tracks which version of the genome a sample was aligned to.
	 *
	 */
	public enum ASSEMBLY_NAME {
		GRCH37, NCBI36, OTHER
	}

	/**
	 * Type of assay
	 *
	 */
	public enum ASSAY_TYPE {
		/**
		 * Whole exome sequencing
		 */
		WXS,
		/**
		 * Whole genome sequencing
		 */
		WGS
	}

	/**
	 * Sequencing platform
	 *
	 */
	public enum PLATFORM {
		ILLUMINA, ABI_SOLID
	}

}
