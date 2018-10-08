package org.genvisis.cnv.startup;

/**
 * Marker interface for responding to problems during {@link StartupValidation}.
 */
public interface StartupErrorHandler {
  // FIXME replace with lambda in Java 8

  /**
   * This method will be called if any problems happen during validation.
   */
  void handleWarnings(String warning);
}
