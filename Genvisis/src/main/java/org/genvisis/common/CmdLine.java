package org.genvisis.common;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class CmdLine {

  public static class Builder {

    private boolean verbose = false;
    private boolean overWriteExistingOutput = false;
    private boolean skipReporting = false;
    private boolean treatEmptyAsMissing = false;
    private boolean ignoreIllegalStateExceptions = false;
    private PrintStream inOS = null;
    private PrintStream errOS = null;
    private Logger log;

    /**
     * @param log
     */
    public Builder(Logger log) {
      this.log = log;
    }

    /**
     * @return a new {@link CmdLine} with the current setup of this {@link Builder}
     */
    public CmdLine build() {
      return new CmdLine(this);
    }

    /**
     * Set verbosity to true
     * 
     * @return this {@link Builder}
     */
    public Builder verbose() {
      verbose = true;
      return this;
    }

    /**
     * Set overWriteExistingOutput to true, {@link Command}s will be run even if their
     * {@link Command#getExpectedOutputFiles()} already exist
     * 
     * @return this {@link Builder}
     */
    public Builder overWriteExistingOutput() {
      overWriteExistingOutput = true;
      return this;
    }

    /**
     * Set skipReporting to true, the output from {@link Command}s will not be logged
     * 
     * @return this {@link Builder}
     */
    public Builder skipReporting() {
      skipReporting = true;
      return this;
    }

    /**
     * Set treatEmptyAsMissing to true, empty files will be treated as missing when resolving
     * {@link Command#getNecessaryInputFiles()} and {@link Command#getExpectedOutputFiles()}
     * 
     * @return this {@link Builder}
     */
    public Builder treatEmptyAsMissing() {
      treatEmptyAsMissing = true;
      return this;
    }

    /**
     * Send the stin stream to inOS
     * 
     * @param inOS
     * @return this {@link Builder}
     */
    public Builder inOS(PrintStream inOS) {
      this.inOS = inOS;
      return this;
    }

    /**
     * Send the sterr stream to errOS
     * 
     * @param errOS
     * @return this {@link Builder}
     */
    public Builder errOS(PrintStream errOS) {
      this.errOS = errOS;
      return this;
    }

    /**
     * Set ignoreIllegalStateExceptions to true, {@link IllegalStateException}s will be caught and
     * not re-thrown
     * 
     * @return this {@link Builder}
     */
    public Builder ignoreIllegalStateExceptions() {
      ignoreIllegalStateExceptions = true;
      return this;
    }

    /**
     * Set verbosity
     * 
     * @param verbose new setting for verbosity
     * @return this {@link Builder}
     */
    private Builder setVerbose(boolean verbose) {
      this.verbose = verbose;
      return this;
    }

    /**
     * Set overWriteExistingOutput, when true {@link Command}s will be run even if their
     * {@link Command#getExpectedOutputFiles()} already exist
     * 
     * @param overWriteExistingOutput new setting for overWriteExistingOutput
     * @return this {@link Builder}
     */
    private Builder setOverWriteExistingOutput(boolean overWriteExistingOutput) {
      this.overWriteExistingOutput = overWriteExistingOutput;
      return this;
    }

    /**
     * Set skipReporting, when true the output from {@link Command}s will not be logged
     * 
     * @param skipReporting new value for skipReporting
     * @return this {@link Builder}
     */
    private Builder setSkipReporting(boolean skipReporting) {
      this.skipReporting = skipReporting;
      return this;
    }

    /**
     * Set treatEmptyAsMissing, when true empty files will be treated as missing when resolving
     * {@link Command#getNecessaryInputFiles()} and {@link Command#getExpectedOutputFiles()}
     * 
     * @param treatEmptyAsMissing new value for treatEmptyAsMissing
     * @return this {@link Builder}
     */
    private Builder setTreatEmptyAsMissing(boolean treatEmptyAsMissing) {
      this.treatEmptyAsMissing = treatEmptyAsMissing;
      return this;
    }

    /**
     * Set ignoreIllegalStateExceptions, when true {@link IllegalStateException}s will be caught and
     * not re-thrown
     * 
     * @param ignoreIllegalStateExceptions new value for ignoreIllegalStateExceptions
     * @return this {@link Builder}
     */
    private Builder setIgnoreIllegalStateExceptions(boolean ignoreIllegalStateExceptions) {
      this.ignoreIllegalStateExceptions = ignoreIllegalStateExceptions;
      return this;
    }

  }

  private final boolean verbose;
  private final boolean overWriteExistingOutput;
  private final boolean skipReporting;
  private final boolean treatEmptyAsMissing;
  private final PrintStream inOS;
  private final PrintStream errOS;
  private final boolean ignoreIllegalStateExceptions;
  private final Logger log;

  private CmdLine(Builder builder) {
    super();
    verbose = builder.verbose;
    overWriteExistingOutput = builder.overWriteExistingOutput;
    skipReporting = builder.skipReporting;
    treatEmptyAsMissing = builder.treatEmptyAsMissing;
    ignoreIllegalStateExceptions = builder.ignoreIllegalStateExceptions;
    inOS = builder.inOS;
    errOS = builder.errOS;
    log = builder.log;
  }

  /**
   * Run the specified {@link Command}
   * 
   * @param command The {@link Command} to run
   * @return true on success, false on failure
   */
  public boolean run(Command command) {
    boolean success = false;
    if (overWriteExistingOutput || command.getExpectedOutputFiles().isEmpty()
        || !Files.exists(command.getDir(), command.getExpectedOutputFiles(), treatEmptyAsMissing)) {
      if (command.getNecessaryInputFiles().isEmpty()
          || Files.exists(command.getDir(), command.getNecessaryInputFiles(),
                          treatEmptyAsMissing)) {
        if (verbose) {
          log.report(ext.getTime() + " Info - running command "
                     + ArrayUtils.toStr(command.getElements(), " "));
        }
        if (runCommand(command)) {
          if (!command.getExpectedOutputFiles().isEmpty()
              && !Files.exists(command.getDir(), command.getExpectedOutputFiles(),
                               treatEmptyAsMissing)) {
            log.reportError("Error - the command " + ArrayUtils.toStr(command.getElements(), " ")
                            + " appeared to run, but could not find all necessary output files in "
                            + command.getDir() + ":"
                            + IterableUtils.toStr(command.getExpectedOutputFiles(), "\n"));
          } else {
            if (verbose) {
              log.report(ext.getTime() + " Info - finished running command "
                         + ArrayUtils.toStr(command.getElements(), " "));
            }
            success = true;
          }
        } else {
          log.reportError("Error - the command " + ArrayUtils.toStr(command.getElements(), " ")
                          + " has failed");
        }
      } else {
        log.reportError("Error - could not find all necessary input files in " + command.getDir()
                        + ":\n" + IterableUtils.toStr(command.getNecessaryInputFiles(), "\n"));
      }
    } else {
      if (verbose) {
        log.report(ext.getTime()
                   + " Info - all of the expected output files exist and the overwrite option was not flagged, skipping:");
        log.report("COMMAND SKIPPED: " + ArrayUtils.toStr(command.getElements(), " "));
      }
      success = true;
    }
    return success;
  }

  private boolean runCommand(Command command) {

    Process proc;
    InputStream in, err;
    boolean finish;
    boolean noError;
    String charSet = "UTF-8";

    noError = true;
    String dir = command.getDir();
    if (dir.equals("")) {
      dir = "./";
    }
    String[] commandArray = command.getElements();
    for (String element : commandArray) {
      if (element.contains(">") || element.contains("|")) {
        if (Files.isWindows()) {
          log.reportError("FYI - the Runtime.exec command will likely not work, since it contains a pipe or redirect, write command to a file and exec that instead");
          break;
        } else if (!commandArray[0].startsWith("/bin/")) {
          log.reportTimeWarning("FYI - the command may not work as it contains a pipe or redirect and isn't prefaced by a shell invocation (such as \"/bin/bash\").  Try prefacing the command with \"/bin/bash -c\" and wrapping the command in quotes.");
          break;
        }
      }
    }

    if (Files.isWindows() && !(commandArray[0].equals("cmd") && commandArray[1].equals("/c"))) {
      log.reportTimeWarning("FYI - the command may not work as it is not prefixed with \"cmd /c\" - this prefix is required to run command-line programs on Windows systems from within Java.");
    }

    try {
      proc = Runtime.getRuntime().exec(commandArray, null, new File(dir));
      in = proc.getInputStream();
      err = proc.getErrorStream();
      finish = false;
      byte[] b;
      while (!finish) {
        try {
          while (in.available() > 0 || err.available() > 0) {
            while (in.available() > 0) {
              b = new byte[in.available()];
              in.read(b);
              if (inOS != null) {
                inOS.print(new String(b, charSet));
              }
              if (!skipReporting) {
                log.report(new String(b, charSet), false, true);
              }
              b = null;
            }
            while (err.available() > 0) {
              b = new byte[err.available()];
              err.read(b);
              if (errOS != null) {
                errOS.print(new String(b, charSet));
              }
              if (!skipReporting) {
                log.report(new String(b, charSet), false, true);
              }
              b = null;
            }
          }
          proc.exitValue();
          finish = true;
        } catch (IllegalThreadStateException e) {
          try {
            Thread.sleep(10);
          } catch (InterruptedException e2) {}
          if (ignoreIllegalStateExceptions) {
            finish = true;
            noError = false;
          }
        }
      }
    } catch (IOException ioe) {
      String message;

      message = ioe.getMessage();
      if (message.startsWith("Cannot run program ")) {
        message = message.substring(20);
        log.reportError("Error - The program \"" + message.substring(0, message.indexOf("\""))
                        + "\" is not installed or is not accessible "
                        + message.substring(message.indexOf("("), message.indexOf(")") + 1));
        noError = false;
      } else {
        log.reportException(ioe);
        noError = false;
      }
    }

    return noError;

  }

  /**
   * @param log
   * @return a basic {@link CmdLine} that doesn't override any of the default settings
   */
  public static CmdLine basic(Logger log) {
    return new Builder(log).build();
  }

  /**
   * A convenience for {@link Builder#Builder(Logger)}
   * 
   * @param log
   * @return a {@link Builder} to build a new {@link CmdLine}
   */
  public static Builder builder(Logger log) {
    return new Builder(log);
  }

  public static String getCmdLocation(String commmand) {
    Map<String, String> env = System.getenv();
    for (String envName : env.keySet()) {
      System.out.format("%s=%s%n", envName, env.get(envName));
    }

    if (env.containsKey(commmand)) {
      return env.get(commmand);
    } else {
      return null;
    }
  }

  /**
   * @param batFile where the the command will be written
   * @param verbose report the command written to the batFile
   * @param log
   * @param commands an array representing the command to run
   * @return String[] of the batFile
   * @deprecated Build a {@link Command} using {@link Command.Builder} including a call to
   *             {@link Command.Builder#batch(String, Logger)} to build and run a batch file
   */
  @Deprecated
  public static String[] prepareBatchForCommandLine(String batFile, boolean verbose, Logger log,
                                                    String... commands) {
    if (verbose) {
      log.report(ext.getTime() + " Info - running command " + ArrayUtils.toStr(commands, " ")
                 + "\nUsing file " + batFile);
    }
    Files.write(ArrayUtils.toStr(commands, " "), batFile);
    Files.chmod(batFile);
    return new String[] {batFile};
  }

  /**
   * @deprecated Build a {@link CmdLine} using {@link Builder} and a {@link Command} using
   *             {@link Command.Builder} and then call {@link #run(Command)} on the {@link CmdLine}
   *             instead. This reduces the possibility for mistakes with repeated args of the same
   *             type and makes for readable code
   */
  @Deprecated
  public static boolean run(String command, String dir) {
    return run(command, dir, null);
  }

  /**
   * @deprecated Build a {@link CmdLine} using {@link Builder} and a {@link Command} using
   *             {@link Command.Builder} and then call {@link #run(Command)} on the {@link CmdLine}
   *             instead. This reduces the possibility for mistakes with repeated args of the same
   *             type and makes for readable code
   */
  @Deprecated
  public static boolean run(String command, String dir, PrintStream os) {
    return run(command, dir, os, false);
  }

  /**
   * @deprecated Build a {@link CmdLine} using {@link Builder} and a {@link Command} using
   *             {@link Command.Builder} and then call {@link #run(Command)} on the {@link CmdLine}
   *             instead. This reduces the possibility for mistakes with repeated args of the same
   *             type and makes for readable code
   */
  @Deprecated
  public static boolean run(String command, String dir, PrintStream os,
                            boolean ignoreIllegalStateExceptions) {
    return run(command, dir, os, os, new Logger(), ignoreIllegalStateExceptions);
  }

  /**
   * @deprecated Build a {@link CmdLine} using {@link Builder} and a {@link Command} using
   *             {@link Command.Builder} and then call {@link #run(Command)} on the {@link CmdLine}
   *             instead. This reduces the possibility for mistakes with repeated args of the same
   *             type and makes for readable code
   */
  @Deprecated
  public static boolean run(Collection<String> commands, String dir, PrintStream inOs,
                            PrintStream errOS, Logger log, boolean ignoreIllegalStateExceptions) {
    return run(commands.toArray(new String[commands.size()]), dir, inOs, errOS, log,
               ignoreIllegalStateExceptions);
  }

  public static boolean run(String command, String dir, PrintStream inOs, PrintStream errOS,
                            Logger log, boolean ignoreIllegalStateExceptions) {
    // StringTokenizer st = new StringTokenizer(command, " \t\n\r\f");
    // String[] cmdarray = new String[st.countTokens()];
    // for (int i = 0; st.hasMoreTokens(); i++) {
    // cmdarray[i] = st.nextToken();
    // }
    String regex = "[\"\']([^\"\']*)[\"\']|(\\S+)";
    Matcher m = Pattern.compile(regex).matcher(command);
    ArrayList<String> cmdList = new ArrayList<String>();
    while (m.find()) {
      if (m.group(1) != null) {
        cmdList.add(m.group(1));
      } else {
        cmdList.add(m.group(2));
      }
    }
    String[] cmdarray = cmdList.toArray(new String[cmdList.size()]);
    return run(cmdarray, dir, inOs, errOS, log, ignoreIllegalStateExceptions);
  }

  /**
   * @deprecated Build a {@link CmdLine} using {@link Builder} and a {@link Command} using
   *             {@link Command.Builder} and then call {@link #run(Command)} on the {@link CmdLine}
   *             instead. This reduces the possibility for mistakes with repeated args of the same
   *             type and makes for readable code
   */
  @Deprecated
  public static boolean run(String[] commandArray, String dir, PrintStream inOS, PrintStream errOS,
                            Logger log, boolean ignoreIllegalStateExceptions) {
    return builder(log == null ? new Logger()
                               : log).inOS(inOS).errOS(errOS)
                                     .setIgnoreIllegalStateExceptions(ignoreIllegalStateExceptions)
                                     .build().run(Command.builder(commandArray).dir(dir).build());
  }

  /**
   * @param commandArray the String array of commands with associated commands grouped at a given
   *          index (String[] commands = new
   *          String[]{"myFavoriteCommand","input=one.txt","output=2.txt"}
   * @param dir directory to run the command in, and also directory to check for existing files
   * @param necessaryInputFiles will check these files for existence, will fail if they do not all
   *          exist
   * @param expectedOutputFiles will check these files for existence, will skip the command if all
   *          exist and overWriteExistingOutput is not flagged
   * @param verbose
   * @param overWriteExistingOutput
   * @param skipReporting if set to true, the cmdline reporting will be sent to the log, otherwise
   *          it will be ignored. Nice to report when error checking, but there is a high
   *          probability of reduced performance and strange output if always used
   * @param log
   * @return
   * @deprecated Build a {@link CmdLine} using {@link Builder} and a {@link Command} using
   *             {@link Command.Builder} and then call {@link #run(Command)} on the {@link CmdLine}
   *             instead. This reduces the possibility for mistakes with repeated args of the same
   *             type and makes for readable code
   */
  @Deprecated
  public static boolean runCommandWithFileChecks(String[] commandArray, String dir,
                                                 String[] necessaryInputFiles,
                                                 String[] expectedOutputFiles, boolean verbose,
                                                 boolean overWriteExistingOutput,
                                                 boolean skipReporting, Logger log) {
    return runCommandWithFileChecks(commandArray, dir, necessaryInputFiles, expectedOutputFiles,
                                    verbose, overWriteExistingOutput, skipReporting, true, log);
  }

  /**
   * @param commandArray the String array of commands with associated commands grouped at a given
   *          index (String[] commands = new
   *          String[]{"myFavoriteCommand","input=one.txt","output=2.txt"}
   * @param dir directory to run the command in, and also directory to check for existing files
   * @param necessaryInputFiles will check these files for existence, will fail if they do not all
   *          exist
   * @param expectedOutputFiles will check these files for existence, will skip the command if all
   *          exist and overWriteExistingOutput is not flagged
   * @param verbose
   * @param overWriteExistingOutput
   * @param skipReporting if set to true, the cmdline reporting will be sent to the log, otherwise
   *          it will be ignored. Nice to report when error checking, but there is a high
   *          probability of reduced performance and strange output if always used
   * @param treatEmptyAsMissing if true, 0 sized files will be set to missing. If 0 sized files are
   *          used as placeholders, set to false
   * @param log
   * @return
   * @deprecated Build a {@link CmdLine} using {@link Builder} and a {@link Command} using
   *             {@link Command.Builder} and then call {@link #run(Command)} on the {@link CmdLine}
   *             instead. This reduces the possibility for mistakes with repeated args of the same
   *             type and makes for readable code
   */
  @Deprecated
  public static boolean runCommandWithFileChecks(String[] commandArray, String dir,
                                                 String[] necessaryInputFiles,
                                                 String[] expectedOutputFiles, boolean verbose,
                                                 boolean overWriteExistingOutput,
                                                 boolean skipReporting, boolean treatEmptyAsMissing,
                                                 Logger log) {
    return builder(log).setVerbose(verbose).setOverWriteExistingOutput(overWriteExistingOutput)
                       .setSkipReporting(skipReporting).setTreatEmptyAsMissing(treatEmptyAsMissing)
                       .build()
                       .run(Command.builder(commandArray).necessaryInputFiles(necessaryInputFiles)
                                   .expectedOutputFiles(expectedOutputFiles).dir(dir).build());
  }

  /**
   * @param commandList the commands as a List of Strings, spaces will be inserted between each
   *          element
   * @param dir directory to run the command in, and also directory to check for existing files
   * @param necessaryInputFiles will check these files for existence, will fail if they do not all
   *          exist
   * @param expectedOutputFiles will check these files for existence, will skip the command if all
   *          exist and overWriteExistingOutput is not flagged
   * @param verbose
   * @param overWriteExistingOutput
   * @param treatEmptyAsMissing if true, 0 sized files will be set to missing. If 0 sized files are
   *          used as placeholders, set to false
   * @param log
   * @return true on success, false on failure
   * @deprecated Build a {@link CmdLine} using {@link Builder} and a {@link Command} using
   *             {@link Command.Builder} and then call {@link #run(Command)} on the {@link CmdLine}
   *             instead. This reduces the possibility for mistakes with repeated args of the same
   *             type and makes for readable code
   */
  @Deprecated
  public static boolean runCommandWithFileChecks(List<String> commandList, String dir,
                                                 Collection<String> necessaryInputFiles,
                                                 Collection<String> expectedOutputFiles,
                                                 boolean verbose, boolean overWriteExistingOutput,
                                                 boolean treatEmptyAsMissing, Logger log) {
    return builder(log).setVerbose(verbose).setOverWriteExistingOutput(overWriteExistingOutput)
                       .setTreatEmptyAsMissing(treatEmptyAsMissing).build()
                       .run(Command.builder(commandList).necessaryInputFiles(necessaryInputFiles)
                                   .expectedOutputFiles(expectedOutputFiles).dir(dir).build());
  }

  /**
   * @param commandList the commands as a List of Strings, spaces will be inserted between each
   *          element
   * @param dir directory to run the command in, and also directory to check for existing files
   * @param necessaryInputFiles will check these files for existence, will fail if they do not all
   *          exist
   * @param expectedOutputFiles will check these files for existence, will skip the command if all
   *          exist and overWriteExistingOutput is not flagged
   * @param verbose
   * @param overWriteExistingOutput
   * @param skipReporting if set to true, the cmdline reporting will be sent to the log, otherwise
   *          it will be ignored. Nice to report when error checking, but there is a high
   *          probability of reduced performance and strange output if always used
   * @param treatEmptyAsMissing if true, 0 sized files will be set to missing. If 0 sized files are
   *          used as placeholders, set to false
   * @param log
   * @return
   * @deprecated Build a {@link CmdLine} using {@link Builder} and a {@link Command} using
   *             {@link Command.Builder} and then call {@link #run(Command)} on the {@link CmdLine}
   *             instead. This reduces the possibility for mistakes with repeated args of the same
   *             type and makes for readable code
   */
  @Deprecated
  public static boolean runCommandWithFileChecks(List<String> commandList, String dir,
                                                 Collection<String> necessaryInputFiles,
                                                 Collection<String> expectedOutputFiles,
                                                 boolean verbose, boolean overWriteExistingOutput,
                                                 boolean skipReporting, boolean treatEmptyAsMissing,
                                                 Logger log) {
    return builder(log).setVerbose(verbose).setOverWriteExistingOutput(overWriteExistingOutput)
                       .setSkipReporting(skipReporting).setTreatEmptyAsMissing(treatEmptyAsMissing)
                       .build()
                       .run(Command.builder(commandList).necessaryInputFiles(necessaryInputFiles)
                                   .expectedOutputFiles(expectedOutputFiles).dir(dir).build());
  }

  /**
   * @deprecated Build a {@link CmdLine} using {@link Builder} and a {@link Command} using
   *             {@link Command.Builder} and then call {@link #run(Command)} on the {@link CmdLine}
   *             instead. This reduces the possibility for mistakes with repeated args of the same
   *             type and makes for readable code
   */
  @Deprecated
  public static boolean runDefaults(String command, String dir) {
    return run(command, dir, System.out, System.err, null, false);
  }

  /**
   * @deprecated Build a {@link CmdLine} using {@link Builder} and a {@link Command} using
   *             {@link Command.Builder} and then call {@link #run(Command)} on the {@link CmdLine}
   *             instead. This reduces the possibility for mistakes with repeated args of the same
   *             type and makes for readable code
   */
  @Deprecated
  public static boolean runDefaults(String command, String dir, Logger log) {
    return run(command, dir, System.out, System.err, log, false);
  }

  /**
   * @deprecated Build a {@link CmdLine} using {@link Builder} and a {@link Command} using
   *             {@link Command.Builder} and then call {@link #run(Command)} on the {@link CmdLine}
   *             instead. This reduces the possibility for mistakes with repeated args of the same
   *             type and makes for readable code
   */
  @Deprecated
  public static boolean runDefaults(Collection<String> command, String dir) {
    return run(command, dir, System.out, System.err, null, false);
  }

  /**
   * @deprecated Build a {@link CmdLine} using {@link Builder} and a {@link Command} using
   *             {@link Command.Builder} and then call {@link #run(Command)} on the {@link CmdLine}
   *             instead. This reduces the possibility for mistakes with repeated args of the same
   *             type and makes for readable code
   */
  @Deprecated
  public static boolean runDefaults(Collection<String> command, String dir, Logger log) {
    return run(command, dir, System.out, System.err, log, false);
  }

}
