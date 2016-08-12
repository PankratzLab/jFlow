package org.genvisis.seq;

/**
 * Common seq variables
 *
 */
public class SeqVariables {

	private SeqVariables() {

	}

	/**
	 * Tracks which version of the genome a sample was aligned to.
	 * 
	 * 
	 * Used to determine mitochondrial sequences, and non-autosomal X and Ys,
	 * etc
	 *
	 * 
	 */
	public enum ASSEMBLY_NAME {

		GRCH37("MT", "X", "Y"), NCBI36("NA", "NA", "NA"), HG19("chrMT", "chrX", "chrY"), OTHER("NA", "NA", "NA");

		private String mitoContig;
		private String xContig;
		private String yContig;

		private ASSEMBLY_NAME(String mitoContig, String xContig, String yContig) {
			this.mitoContig = mitoContig;
			this.xContig = xContig;
			this.yContig = yContig;
		}

		/**
		 * @return mito contig string
		 */
		public String getMitoContig() {
			return mitoContig;
		}

		/**
		 * @return X contig string
		 */
		public String getxContig() {
			return xContig;
		}

		/**
		 * @return Y contig string
		 */
		public String getyContig() {
			return yContig;
		}

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
