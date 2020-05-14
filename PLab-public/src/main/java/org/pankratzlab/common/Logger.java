package org.pankratzlab.common;

import java.io.File;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.Serializable;

import javax.swing.JTextArea;

import com.google.common.base.Joiner;
import com.google.common.base.Strings;

public class Logger implements Serializable {

  private static final long serialVersionUID = 1L;
  public static final int LEVEL_ONE = 1;
  public static final int DEFAULT_LOG_LEVEL = 10;

  private final String filename;
  private final boolean logging;
  private int level;
  private JTextArea textArea;

  public Logger() {
    this(null, false);
  }

  public Logger(String filename) {
    this(filename, false);
  }

  public Logger(String filename, boolean append) {
    this(filename, append, DEFAULT_LOG_LEVEL);
  }

  public Logger(String filename, boolean append, int level) {
    this.filename = filename;
    logging = filename != null;
    if (logging && !append) {
      new File(filename).delete();
    }
    this.level = level;
    textArea = null;
  }

  public void linkTextArea(JTextArea text) {
    textArea = text;
  }

  public String getFilename() {
    return filename;
  }

  public int getLevel() {
    return level;
  }

  public void setLevel(int level) {
    this.level = level;
  }

  public void reportTime(String str) {
    report(ext.getTime() + "]\t" + str, true, true);
  }

  /**
   * @param str report this string with a time stamp and info message
   */
  public void reportTimeInfo(String str) {
    report(ext.getTime() + "]\t Info - " + str, true, true);
  }

  /**
   * @param str report this string with a time stamp and info message
   */
  public void reportTimeInfo(String str, int level) {
    report(ext.getTime() + "]\t Info - " + str, true, true, level);
  }

  /**
   * @param str report this string with a time stamp and info message
   * @return the time elapsed
   */
  public String reportTimeElapsed(long time) {
    String elapsed = ext.getTimeElapsed(time);
    reportTimeInfo("Time elapsed: " + elapsed);
    return elapsed;
  }

  /**
   * @param str report this string with a time stamp and info message
   * @return the time elapsed
   */
  public String reportTimeElapsed(String prepend, long time) {
    String elapsed = ext.getTimeElapsed(time);
    reportTimeInfo(prepend + elapsed);
    return elapsed;
  }

  /**
   * @param file report this file already exists with a time stamp and info message
   */
  public void reportFileExists(String file) {
    reportTimeInfo("File " + file + " already exists...");
  }

  /**
   * @param str report this string with a time stamp and warning message
   */
  public void reportTimeWarning(String str) {
    report(ext.getTime() + "]\t Warning - " + str, true, true);
  }

  private static String getVersion() {
    String v = "[Genvisis - version unknown] ";
    try {
      v = "[Genvisis " + LauncherManifest.loadGenvisisManifest().getVersion() + "] ";
    } catch (Exception e) {

    }
    return v;
  }

  /**
   * @param filename report that this file was not found with a time stamp and error message
   */
  public void reportFileNotFound(String filename) {
    String str = "file \"" + filename + "\" not found in current directory";
    reportError(str);
  }

  /**
   * @param filenames report the files that do not exist in this string array
   */
  public void reportFilesNotFound(String[] filenames) {
    if (filenames == null) {
      reportError("No file handles provided");
    } else {
      boolean haveMissing = false;
      for (int i = 0; i < filenames.length; i++) {
        if (!Files.exists(filenames[i])) {
          reportFileNotFound(filenames[i]);
          haveMissing = true;
        }
      }
      if (!haveMissing) {
        reportTimeInfo("All " + filenames.length + " files exist");
      }
    }
  }

  /**
   * @param filename report that this file had an IO err with time stamp and error message
   */
  public void reportIOException(String filename) {
    String str = "could not read file \"" + filename + "\"";
    reportError(str);
  }

  public void report(String str) {
    report(str, true, true);
  }

  public void report(String str, boolean line, boolean reportToScreen) {
    report(str, line, reportToScreen, 0);
  }

  public void report(String str, boolean line, boolean reportToScreen, int levelRequiredToReport) {
    PrintWriter writer;

    if (level >= levelRequiredToReport && reportToScreen) {
      if (line) {
        System.out.println(str);
      } else {
        System.out.print(str);
      }
    }

    if (level >= levelRequiredToReport && textArea != null) {
      textArea.setText(textArea.getText() + str + (line ? "\r\n" : ""));
      textArea.setCaretPosition(textArea.getDocument().getLength());
    }

    if (level >= levelRequiredToReport && logging) {
      try {
        writer = Files.openAppropriateWriter(filename, true);
        if (line) {
          writer.println(str);
        } else {
          writer.print(str);
        }
        writer.close();
      } catch (Exception e) {
        System.err.println("Error writing to " + filename);
        e.printStackTrace();
      }
    }
  }

  /**
   * @see #reportError(boolean, boolean, int, String...)
   */
  public void reportError(String err) {
    reportError(true, true, err);
  }

  /**
   * @see #reportError(boolean, boolean, int, String...)
   */
  public void reportError(String... err) {
    reportError(true, true, err);
  }

  /**
   * @see #reportError(boolean, boolean, int, String...)
   */
  public void reportError(String err, boolean line, boolean reportToScreen) {
    reportError(line, reportToScreen, err);
  }

  /**
   * @see #reportError(boolean, boolean, int, String...)
   */
  public void reportError(boolean line, boolean reportToScreen, String... err) {
    reportError(line, reportToScreen, 0, err);
  }

  /**
   * @see #reportError(boolean, boolean, int, String...)
   */
  public void reportError(String err, boolean line, boolean reportToScreen,
                          int levelRequiredToReport) {
    reportError(line, reportToScreen, levelRequiredToReport, err);
  }

  /**
   * @param line use {@link PrintStream#println()} to print error
   * @param reportToScreen report error to screen
   * @param levelRequiredToReport minimum logger level to report error
   * @param err lines of error message to print
   */
  public void reportError(boolean line, boolean reportToScreen, int levelRequiredToReport,
                          String... errLines) {
    PrintWriter writer;
    final String errorPrefix = "Error -";
    final boolean alreadyHasErrorPrefix = errLines.length > 0
                                          && errLines[0].regionMatches(true, 0, errorPrefix, 0, 5);
    String firstPrefix = getVersion() + ext.getTime() + (alreadyHasErrorPrefix ? "" : errorPrefix);
    String otherPrefix = Strings.repeat(" ", firstPrefix.length());
    String msg = firstPrefix + "\t"
                 + Joiner.on("\n" + otherPrefix + "\t").skipNulls().join(errLines);

    if (level >= levelRequiredToReport && reportToScreen) {
      if (line) {
        System.err.println(msg);
      } else {
        System.err.print(msg);
      }
    }

    if (level >= levelRequiredToReport && textArea != null) {
      textArea.setText(textArea.getText() + msg + (line ? "\r\n" : ""));
      textArea.setCaretPosition(textArea.getDocument().getLength());
    }

    if (level >= levelRequiredToReport && logging) {
      try {
        writer = Files.openAppropriateWriter(filename, true);
        if (line) {
          writer.println(msg);
        } else {
          writer.print(msg);
        }
        writer.close();
      } catch (Exception e) {
        System.err.println("Error writing to " + filename);
        e.printStackTrace();
      }
    }
  }

  public void reportException(Throwable e) {
    reportException(e, 0);
  }

  public void reportException(Throwable e, int levelRequiredToReport) {
    PrintWriter writer;
    String msg = getVersion() + ext.getTime() + "\t" + e.getMessage();

    System.err.println(msg);
    e.printStackTrace();
    if (level >= levelRequiredToReport && logging) {
      e.printStackTrace();
      try {
        writer = Files.openAppropriateWriter(filename, true);
        writer.println(msg);
        e.printStackTrace(writer);
        writer.close();
      } catch (Exception e2) {
        System.err.println("Error writing to " + filename);
        e2.printStackTrace();
      }
    }

    if (level >= 11) {
      System.exit(1);
    }
  }

  public void timestamp() {
    report(ext.getDate() + "\t" + ext.getTime());
  }

  public long memoryTotal() {
    long memory;

    report("Total heap size is: "
           + ext.prettyUpSize(memory = Runtime.getRuntime().totalMemory(), 1));

    return memory;
  }

  public long memoryFree() {
    long memory;

    report("Free heap size is: " + ext.prettyUpSize(memory = Runtime.getRuntime().freeMemory(), 1));

    return memory;
  }

  /**
   * @return percent free memory of the maximum available
   */
  public double memoryPercentTotalFree() {
    double used = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
    double percentTotalFree = 100 - (100 * used / Runtime.getRuntime().maxMemory());
    report("Percent free total heap size is: " + ext.formDeci(percentTotalFree, 1) + "%");

    return percentTotalFree;
  }

  public double memoryPercentFree() {
    double percentFree;

    percentFree = ((float) 100 * Runtime.getRuntime().freeMemory()
                   / Runtime.getRuntime().totalMemory());

    report("Percent free heap size is: " + ext.formDeci(percentFree, 1) + "%");

    return percentFree;
  }

  public long memoryUsed() {
    long memory;

    report("Used heap size is: " + ext.prettyUpSize(
                                                    memory = (Runtime.getRuntime().totalMemory()
                                                              - Runtime.getRuntime().freeMemory()),
                                                    1));

    return memory;
  }

  public long memoryMax() {
    long memory;

    report("Max available heap size is: "
           + ext.prettyUpSize(memory = Runtime.getRuntime().maxMemory(), 1));

    return memory;
  }

}
