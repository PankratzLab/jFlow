/**
 * 
 */
package org.genvisis.seq.cnv;

/**
 * Typcial getters for packages (SVtyper/Lumpy) using split/discordant reads for SV detection
 */
public interface PairedEndSVAnalysis {

  /**
   * @return the discordant bam file
   */
  public String getDiscordantBam();

  /**
   * @return the splitter bam
   */
  public String getSplitterBam();

  /**
   * @return the base bam
   */
  public String getBaseBam();

  /**
   * @return did it fail
   */
  public boolean isFail();
}
