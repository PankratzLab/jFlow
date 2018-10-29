package org.pankratzlab.common.startup;

import java.io.IOException;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;

/**
 * Utility class for executing {@link StartupCheck} plugins and reporting the results.
 */
public final class StartupValidation {

  /*
   * Constants representing state of the validation check. Originally used enums but there seem to
   * be annoyances regarding synchronized locks and enums.
   */
  private static final int NOT_STARTED = -1;
  private static final int IN_PROGRESS = 0;
  private static final int DONE = 1;

  // Locking object for multi-thread syncrhonization. Originally tried using warningString, but
  // immutable objects seem to have issues obtaining object monitors on multiple threads.
  private static final Object LOCK = new Object();

  private static int status = -1;
  private static String warningString = "";
  private static List<StartupErrorHandler> handlers = new ArrayList<>();

  private StartupValidation() {
    // prevent instantiation of utility class
  }

  /**
   * Blocking method to wait for validation to complete and return the result of validation.
   *
   * @return True iff validation passed with no problems.
   */
  public static boolean passed() {
    if (status == IN_PROGRESS) {
      System.out.println("Waiting for startup validation to complete");

      while (status == IN_PROGRESS) {
        try {
          synchronized (LOCK) {
            // Wait for LOCK.notify()
            LOCK.wait();
          }
        } catch (InterruptedException e) {
          e.printStackTrace();
        }
      }
    }

    return warningString.isEmpty();
  }

  /**
   * Perform validation, calling {@link StartupCheck#doCheck()} for all registered instances and
   * recording the results.
   *
   * @param onError - {@link StartupErrorHandler#handleWarnings(String)} will be called if any
   *          problems occur.
   */
  public static void validate(StartupErrorHandler onError, Consumer<Boolean> preloadMethod,
                              StartupCheck... toCheck) {
    switch (status) {
      case DONE:
        report(onError);
        return;
      case NOT_STARTED:
        startValidation(preloadMethod, toCheck);
        break;
      default:
        break;
    }

    register(onError);
  }

  /**
   * Synchronize-locked helper method to start validating in a separate thread
   * 
   * @param preloadMethod
   */
  private static synchronized void startValidation(Consumer<Boolean> preloadMethod,
                                                   StartupCheck... toCheck) {
    if (status == NOT_STARTED) {
      status = IN_PROGRESS;
      new Thread(new Runnable() {

        @Override
        public void run() {
          doValidation(preloadMethod, toCheck);
        }

      }).start();
    }
  }

  /**
   * @return True if we have internet access
   */
  private static boolean netCheck() {
    boolean haveInternet = true;

    // Check for internet connection, which will determine which checks we perform
    try {
      URL url = new URL("http://genvisis.org/index.html");
      HttpURLConnection httpConn = (HttpURLConnection) url.openConnection();
      httpConn.setConnectTimeout(10 * 1000);
      httpConn.connect();
      httpConn.disconnect();
    } catch (IOException e) {
      haveInternet = false;
    }
    return haveInternet;
  }

  /**
   * Perform validation checks
   * 
   * @param preloadMethod
   */
  private static void doValidation(Consumer<Boolean> preloadMethod, StartupCheck... toCheck) {
    List<String> warnings = new ArrayList<>();
    boolean haveInternet = netCheck();

    // Run all checks
    for (StartupCheck check : toCheck) {
      if (haveInternet || !check.requiresRemote()) {
        List<String> output = check.check();
        if (!output.isEmpty()) {
          warnings.addAll(output);
        }
      }
    }

    // Create the warning message
    buildWarningString(warnings);

    // Change the status and report any problems
    // Needs to be synchronized with the register method to avoid threading problems.
    synchronized (handlers) {
      status = DONE;
      for (StartupErrorHandler handler : handlers) {
        report(handler);
      }
    }
    synchronized (LOCK) {
      // Wake up any threads waiting for #passed() checks.
      LOCK.notifyAll();
    }
    preloadMethod.accept(haveInternet);
  }

  /**
   * Helper method to register a {@link StartupErrorHandler} to run after validation is complete. If
   * validation is already complete, immediately calls {@link #report(StartupErrorHandler)}.
   */
  private static void register(StartupErrorHandler onError) {
    synchronized (handlers) {
      if (status != DONE) {
        handlers.add(onError);
      } else {
        report(onError);
      }
    }
  }

  /**
   * Helper method so individual {@link StartupErrorHandler}s don't have to check whether or not
   * validation had errors when calling {@link StartupErrorHandler#handleWarnings(String)}.
   */
  private static void report(StartupErrorHandler onError) {
    if (!passed()) {
      onError.handleWarnings(warningString);
    }
  }

  /**
   * Helper method to convert a list of errors to a formatted string.
   */
  private static void buildWarningString(List<String> errors) {
    StringBuilder sb = new StringBuilder();
    for (String msg : errors) {
      sb.append(msg);
      sb.append("\n");
    }
    warningString = sb.toString();
  }
}