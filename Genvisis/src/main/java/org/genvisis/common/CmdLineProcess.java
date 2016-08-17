package org.genvisis.common;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.Iterator;

public class CmdLineProcess implements Iterator<String> {
  /**
   * Builder for the run, defaults are set below
   *
   */
  public static class Builder {
    private String dir = null;
    private StandardInputProvider STIN = null;
    private INPUT_Mode inputMode = INPUT_Mode.NO_STIN;
    private OUTPUT_Mode outputMode = OUTPUT_Mode.NO_STOUT_CAPTURE;
    private ERR_Mode errorMode = ERR_Mode.STERR_CAPTURE_BY_LOG;
    private Logger log = new Logger();
    private boolean verbose = false;
    private int lineBufferSize = 10000;

    public CmdLineProcess build(String[] commandArray) {
      return new CmdLineProcess(this, commandArray);
    }

    public Builder dir(String dir) {
      this.dir = dir;
      return this;
    }

    public Builder errorMode(ERR_Mode errorMode) {
      this.errorMode = errorMode;
      return this;
    }

    public Builder inputMode(INPUT_Mode inputMode) {
      this.inputMode = inputMode;
      return this;
    }

    public Builder lineBufferSize(int lineBufferSize) {
      this.lineBufferSize = lineBufferSize;
      return this;
    }

    public Builder log(Logger log) {
      this.log = log;
      return this;
    }

    public Builder outputMode(OUTPUT_Mode outputMode) {
      this.outputMode = outputMode;
      return this;
    }

    public Builder STIN(StandardInputProvider STIN) {
      this.STIN = STIN;
      return this;
    }

    public Builder verbose(boolean verbose) {
      this.verbose = verbose;
      return this;
    }
  }

  public enum ERR_Mode {
    /**
     * Report errors to log
     */
    STERR_CAPTURE_BY_LOG, NO_STERR_CAPTURE
  }

  public enum INPUT_Mode {
    /**
     * Pass input to the command
     */
    STIN, NO_STIN
  }

  public enum INPUT_TYPE {
    /**
     * Standard in comes from a file
     */
    FILE,
    /**
     * Standard in comes from a string
     */
    STRING
  }

  class InputWriter extends Thread {

    public InputWriter() {}

    @Override
    public void run() {
      int added = 0;

      while (STIN.hasNext()) {
        String in = STIN.next();
        stIn.println(in);
        added++;
        if (verbose && added % 10000 == 0) {
          log.reportTimeInfo("Added " + added + " from STIN");
        }
        // if (added >= lineBufferSize) {
        // try {
        //
        // } catch (InterruptedException e) {
        // // TODO Auto-generated catch block
        // e.printStackTrace();
        // }
        // added = 0;
        // }
        stIn.flush();
      }
      stIn.close();
    }
  }
  public enum OUTPUT_Mode {
    /**
     * The output is captured and returned line by line
     */
    STOUT_CAPTURE_ITERATOR,
    /**
     * The output is not captured (inherited)
     */
    NO_STOUT_CAPTURE
  }
  /**
   * Interface that can be passed to the command line to provide input
   *
   */
  public static interface StandardInputProvider extends Iterator<String> {

  }

  private final String[] commandArray;
  private final String dir;
  private String iterOutputLine;
  private final StandardInputProvider STIN;
  private PrintWriter stIn;
  private BufferedReader stOut;
  private BufferedReader stErr;
  private final INPUT_Mode inputMode;
  private final OUTPUT_Mode outputMode;
  private final ERR_Mode errorMode;
  private final int lineBufferSize;
  private InputWriter inputWriter;

  private Process proc;

  private final boolean verbose, fail;

  private final Logger log;

  private CmdLineProcess(Builder builder, String[] commandArray) {
    this.commandArray = commandArray;
    dir = builder.dir;
    STIN = builder.STIN;
    inputMode = builder.inputMode;
    outputMode = builder.outputMode;
    errorMode = builder.errorMode;
    log = builder.log;
    verbose = builder.verbose;
    iterOutputLine = null;
    lineBufferSize = builder.lineBufferSize;
    fail = !initProcess();
  }

  private void flushErrorStream() {

    if (stErr != null) {
      String line;
      try {
        while (stErr.ready() && (line = stErr.readLine()) != null) {
          switch (errorMode) {
            case STERR_CAPTURE_BY_LOG:
              if (line != null) {
                // System.out.println(line);
                log.reportTimeError(line);
              }
              break;
            default:
              break;

          }
        }
      } catch (IOException e) {
        log.reportTimeError(
            "Could not read error stream for process " + Array.toStr(commandArray, " "));
        e.printStackTrace();
      }
    }
  }

  public InputWriter getInputWriter() {
    return inputWriter;
  }

  public int getLineBufferSize() {
    return lineBufferSize;
  }

  private void handleOutputStream() {
    switch (outputMode) {
      case NO_STOUT_CAPTURE:
        break;
      case STOUT_CAPTURE_ITERATOR:
        try {
          if (proc.getInputStream() != null && proc.getInputStream().available() > 0) {
            log.reportTimeError(
                "The output mode was set to iterate, but there is still available bytes.");
          }
        } catch (IOException e) {
          // TODO Auto-generated catch block
          e.printStackTrace();
        }
        break;
      default:
        break;

    }
  }

  @Override
  public boolean hasNext() {
    boolean next = true;
    flushErrorStream();

    if (!fail && stOut != null) {
      try {
        iterOutputLine = stOut.readLine();
        next = iterOutputLine != null;
      } catch (IOException e2) {
        next = false;
        // TODO Auto-generated catch block
        e2.printStackTrace();
      }
    } else {
      next = false;
    }

    return next;
  }

  /**
   * Sets up error streams,input streams, and output streams
   * 
   * @return process was initiated
   */
  private boolean initProcess() {
    boolean init = true;
    if (verbose) {
      log.reportTimeInfo("Attempting to run command " + Array.toStr(commandArray, " "));
    }
    ProcessBuilder probuilder = new ProcessBuilder(commandArray);

    if (dir != null && dir != "") {
      probuilder.directory(new File(dir));
    }
    setupStandardError(probuilder);

    try {
      proc = probuilder.start();
      switch (inputMode) {
        case NO_STIN:
          break;
        case STIN:

          if (STIN != null) {
            if (verbose) {
              log.reportTimeInfo(
                  "Passing input paramaters to command:  " + Array.toStr(commandArray));

            }
            stIn = new PrintWriter(proc.getOutputStream());
            writeInput();
            if (verbose) {
              log.reportTimeInfo(
                  "Finished passing input paramaters to command:  " + Array.toStr(commandArray));
            }

          } else {
            log.reportTimeError("A standard input source was not provided");
          }

          break;
        default:
          break;

      }

      if (outputMode == OUTPUT_Mode.STOUT_CAPTURE_ITERATOR) {
        if (verbose) {
          log.reportTimeInfo("Preparing standard out");
        }
        stOut = new BufferedReader(new InputStreamReader(proc.getInputStream(), "UTF-8"));

        if (verbose) {
          log.reportTimeInfo("Finished preparing standard out");
        }
      }
      if (errorMode == ERR_Mode.STERR_CAPTURE_BY_LOG) {
        stErr = new BufferedReader(new InputStreamReader(proc.getErrorStream(), "UTF-8"));
      }

    } catch (IOException e) {
      init = false;
      log.reportTimeError("Could not initialize process " + Array.toStr(commandArray));
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    return init;
  }

  @Override
  public String next() {
    if (!fail && stOut != null) {
      return iterOutputLine;

    } else {
      return null;
    }
  }

  @Override
  public void remove() {
    // TODO Auto-generated method stub

  }

  //

  private void setupStandardError(ProcessBuilder probuilder) {
    switch (errorMode) {
      case NO_STERR_CAPTURE:
        probuilder.redirectErrorStream(true);
        break;

      case STERR_CAPTURE_BY_LOG:
        break;
      default:
        log.reportTimeWarning("Invalid error capture mode, directing to standard output...");
        probuilder.redirectErrorStream(true);
        break;
    }
  }

  public boolean waitFor() {
    boolean error = false;
    try {
      // System.err.println("flush Error");

      flushErrorStream();
      // System.err.println("flush Out");

      handleOutputStream();
      // System.err.println("finished Out");

      proc.waitFor(); // wait for process to complete
      // System.err.println("finished Wait");

    } catch (InterruptedException e) {
      log.reportException(e);
      System.err.println(e); // "Can'tHappen"
      error = true;
    }
    return error;
  }

  private void writeInput() {
    new InputWriter().start();
  }
}

// try {
// if (proc.exitValue() == 0) {// if the process is done
// try {
// // System.out.println("has next");
//
// next = stOut.ready();// if there is anything to be read
// } catch (IOException e) {
// e.printStackTrace();
// }
// }
// } catch (IllegalThreadStateException e) {
// try {
// Thread.sleep(10);
// } catch (InterruptedException e1) {
// e1.printStackTrace();
// }
// }
// public static class StandardInputBasic implements StandardInputStreamReader {
//
// private String in;
// private INPUT_TYPE type;
// private Logger log;
//
// public StandardInputBasic(String in, INPUT_TYPE type, Logger log) {
// super();
// this.in = in;
// this.type = type;
// this.log = log;
// }
//
// @Override
// public BufferedReader getStandardInputStreamReader() {
// BufferedReader reader = null;
// switch (type) {
// case FILE:
// try {
// reader = Files.getAppropriateReader(in);
// } catch (FileNotFoundException e) {
// log.reportFileNotFound(in);
// e.printStackTrace();
// }
// break;
// case STRING:
// InputStream is = new ByteArrayInputStream(in.getBytes());
// reader = new BufferedReader(new InputStreamReader(is));
// break;
// default:
// reader = null;
// break;
//
// }
// return reader;
// }
// // }
//
// while (!next && !done) {
// try {
// next = stOut.ready();
// int test = 1 + proc.exitValue();
// done = true; //
// next = stOut.ready();
// } catch (IllegalThreadStateException illegalThreadStateException) {
// try {
// Thread.sleep(100);
// } catch (InterruptedException e) {
// // TODO Auto-generated catch block
// e.printStackTrace();
// }
// }
// }

// boolean done = false;
// if (verbose) {
// log.reportTimeInfo("Begining to flush iteration's error stream");
// }
// flushErrorStream();
// if (verbose) {
// log.reportTimeInfo("Finished flushing iterations error stream");
// }
// if (!fail && stOut != null) {
//
// iterOutputLine = stOut.readLine();
//
// try {
// try {
// proc.exitValue();
// next = proc.getInputStream().available() > 0;
//
// // System.out.println("next");
// } catch (IllegalThreadStateException e) {
// try {
// Thread.sleep(10);
// } catch (InterruptedException e1) {
// // TODO Auto-generated catch block
// e1.printStackTrace();
// }
//
// }
// } catch (IOException e1) {
// // TODO Auto-generated catch block
// e1.printStackTrace();
// }
//
// } else {
// next = false;
// }
// System.out.println("Passing this on " + next);
// try {
// boolean done = false;
// while (!done && !stOut.ready()) {// fill the buffer
// System.err.println("Here2");
//
// try {
// System.err.println("HI");
// proc.exitValue();
// done = true;
// // System.out.println("next");
// } catch (IllegalThreadStateException e) {
// try {
// Thread.sleep(100);
// } catch (InterruptedException e1) {
// // TODO Auto-generated catch block
// e1.printStackTrace();
// }
// }
//
// }
// System.err.println("Here3");
//
// String line = stOut.readLine();
// System.err.println("Here4");
//
// if (line == null) {
// log.reportTimeError("There was no standard out available for this process");
// }
// System.err.println("LINE " + line);
//
// return line;
// } catch (IOException e) {
// fail = true;
// log.reportTimeError("Could not read output stream for command " + Array.toStr(commandArray, "
// "));
// e.printStackTrace();
// e.printStackTrace();
// return null;
//
// }
