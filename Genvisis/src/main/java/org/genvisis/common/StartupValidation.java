package org.genvisis.common;

import java.io.IOException;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;

import org.genvisis.cnv.manage.Resources.ResourceVersionCheck;
import org.genvisis.common.HttpUpdate.RemoteVersionCheck;

/**
 * Utility class for executing {@link StartupCheck} plugins and reporting the results.
 *
 */
public final class StartupValidation {

	private static boolean validated = false;
	private static String warningString = "";
	private static StartupCheck[] toCheck = new StartupCheck[] {new ResourceVersionCheck(),
																															new RemoteVersionCheck()};

	private StartupValidation() {
		// prevent instantiation of utility class
	}

	/**
	 * Perform validation, calling {@link StartupCheck#doCheck()} for all registered instances and
	 * recording the results.
	 *
	 * @return True if there were errors during validation.
	 */
	public static boolean validate() {
		if (!validated) {
			validated = true;
			List<String> warnings = new ArrayList<String>();
			boolean haveInternet = true;


			try {
				URL url = new URL("http://genvisis.org/index.html");
				HttpURLConnection httpConn = (HttpURLConnection) url.openConnection();
				httpConn.setConnectTimeout(10 * 1000);
				httpConn.connect();
				httpConn.disconnect();
			} catch (IOException e) {
				haveInternet = false;
			}

			for (StartupCheck check : toCheck) {
				if (haveInternet || !check.requiresRemote()) {
					List<String> output = check.check();
					if (!output.isEmpty()) {
						warnings.addAll(output);
					}
				}
			}
			buildWarningString(warnings);
		}
		return !warningString.isEmpty();
	}

	/**
	 * @return A formatted string containing all warnings that occurred during {@link #validate()}.
	 */
	public static String warnings() {
		return warningString;
	}

	/**
	 * Helper method to convert the array to an output string.
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
