package org.genvisis.cnv.startup;

import java.util.List;

/**
 * Interface for classes that will be registered with {@link StartupValidation}, to perform a
 * specific check at startup.
 * <p>
 * NOTE: Currently all StartupChecks must be manually registered in {@link StartupValidation}. If a
 * more extensible model is desired, recommend using the Context in SciJava common
 * (http://www.scijava.org/).
 * </p>
 */
public interface StartupCheck {

  /**
   * Perform any necessary validation.
   *
   * @return A list of messages to report based on the check results.
   */
  List<String> check();

  /**
   * @return True if this startup check requires an internet connection
   */
  boolean requiresRemote();
}
