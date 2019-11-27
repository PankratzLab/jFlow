package org.genvisis.common.gui;

import org.pankratzlab.common.Logger;

/**
 * This helper class is a catch-all for any exceptions that would not otherwise be caught.
 */
public final class ExceptionHandler implements Thread.UncaughtExceptionHandler {

  public static final String X11_ERROR_MSG_FORE = "Error occurred with X11 forwarding - ";
  public static final String X11_ERROR_DISABLED = "it's likely that X11 forwarding is disabled; please check your SSH client settings and try again.";
  public static final String X11_ERROR_XMING_REC = "it's likely that X11 forwarding is enabled but you are missing an X11 forwarding server (we recommend Xming - http://sourceforge.net/projects/xming/)";

  public Logger log;

  public void setLog(Logger log) {
    this.log = log;
  }

  @Override
  public void uncaughtException(Thread t, Throwable e) {
    if (log != null) {
      log.reportError("Uncaught Exception in Thread {" + t.getName() + "}: " + e.getMessage());
      log.reportException(e);
    } else {
      System.err.println("Error - Uncaught Exception in Thread {" + t.getName() + "}:");
      e.printStackTrace();
    }
  }
}