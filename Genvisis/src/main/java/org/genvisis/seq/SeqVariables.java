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
   * Used to determine mitochondrial sequences, and non-autosomal X and Ys, etc
   *
   *
   */
  public enum ASSEMBLY_NAME {

                             GRCH37("MT", "X", "Y",
                                    false), NCBI36("NA", "NA", "NA",
                                                   false), HG19("chrMT", "chrX", "chrY",
                                                                true), OTHER("NA", "NA", "NA",
                                                                             false);

    private String mitoContig;
    private String xContig;
    private String yContig;
    private boolean addChr;

    private ASSEMBLY_NAME(String mitoContig, String xContig, String yContig, boolean addChr) {
      this.mitoContig = mitoContig;
      this.xContig = xContig;
      this.yContig = yContig;
      this.addChr = addChr;
    }

    /**
     * @return mito contig string
     */
    public String getMitoContig() {
      return mitoContig;
    }


    /**
     * @return whether to add "chr" for all contigs
     */
    public boolean addChr() {
      return addChr;
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
