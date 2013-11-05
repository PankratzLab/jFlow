package stats;

import common.Files;
import common.Logger;

public class Rscript {
	public static final String[] RSCRIPT_EXECS = {
//		"/soft/R/3.0.1/bin/Rscript", // MSI
		"/soft/R/2.15.1/bin/Rscript", // Itasca nodes can only see this
		"/share/apps/R-3.0.1/bin/Rscript", // psych
		"/share/apps/src/R-3.0.1/bin/Rscript", // alcatraz
	};
	
	public static final String[] R_EXECS = {
//		"/soft/R/3.0.1/bin/R", // MSI
		"/soft/R/2.15.1/bin/R", // Itasca nodes can only see this
		"/share/apps/R-3.0.1/bin/R", // psych
		"/share/apps/src/R-3.0.1/bin/R", // alcatraz
	};
	
	public static String getRscriptExecutable(Logger log) {
		if (System.getProperty("os.name").startsWith("Windows")) {
			return "Rscript";
		} else {
			for (int i = 0; i < RSCRIPT_EXECS.length; i++) {
				if (Files.exists(RSCRIPT_EXECS[i])) {
					return RSCRIPT_EXECS[i];
				}
			}
		}

		log.report("Warning - could not determine hostname, assuming R executbale is simply 'Rscript'");
		
		return "Rscript";
	}

	public static String getRExecutable(Logger log) {
		if (System.getProperty("os.name").startsWith("Windows")) {
			return "R";
		} else {
			for (int i = 0; i < R_EXECS.length; i++) {
				if (Files.exists(R_EXECS[i])) {
					return R_EXECS[i];
				}
			}
		}

		log.report("Warning - could not determine hostname, assuming R executbale is simply 'Rscript'");
		
		return "R";
	}
}
