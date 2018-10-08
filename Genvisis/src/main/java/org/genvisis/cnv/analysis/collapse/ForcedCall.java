/**
 * 
 */
package org.genvisis.cnv.analysis.collapse;

import org.genvisis.cnv.filesys.CNVariant;

/**
 * Specific type of {@link CNVariant} that was force called by{@link ForcedCalling}
 */
public class ForcedCall extends CNVariant {

  /**
   * 
   */
  private static final long serialVersionUID = 1L;

  /**
   * @param builder
   */
  public ForcedCall(CNVBuilder builder) {
    super(builder);
  }
}
