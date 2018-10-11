package org.genvisis.cnv.gwas.windows;

import org.genvisis.cnv.filesys.Project;
import org.pankratzlab.common.Logger;
import org.pankratzlab.utils.gwas.windows.MultiHitWindows;
import org.pankratzlab.utils.gwas.windows.PFileStructure;
import org.pankratzlab.utils.gwas.windows.WindowThreshold;

/**
 * Utility class for convenience methods that can operate on a project property
 */
public final class ProjectUtils {

  /**
   * Use the given {@link Project} to find logging and centromere splitting info.
   *
   * @see #generateHitWindows(String, PFileStructure, WindowThreshold, int, String, Logger)
   */
  public static void generateHitWindows(String outSuffix, PFileStructure pFile,
                                        WindowThreshold hitParams, Project proj) {
    MultiHitWindows.generateHitWindows(outSuffix, pFile, hitParams,
                                       proj.GENOME_BUILD_VERSION.getValue().getBuildInt(),
                                       proj.PLINK_DIR_FILEROOTS.getValueString() + "plink.bim",
                                       proj.getLog());
  }
}
