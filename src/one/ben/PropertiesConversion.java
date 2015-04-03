package one.ben;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;

import common.Files;

public class PropertiesConversion {
	
	private static HashMap<String, String> replacements = new HashMap<String, String>() {{
		put("DISPLAY_QUANTILES", "DISPLAY_QUANTILES_YESNO");
		put("DISPLAY_ROTATED_QQ", "DISPLAY_ROTATED_QQ_YESNO");
		put("DISPLAY_STANDARD_QQ", "DISPLAY_STANDARD_QQ_YESNO");
		put("FID_ALIAS", "FID_ALIAS_LIST");
		put("GC_THRESHOLD", "GC_THRESHOLD_0d0_1d0_D");
		put("IID_ALIAS", "IID_ALIAS_LIST");
		put("INTENSITY_PC_NUM_COMPONENTS", "INTENSITY_PC_NUM_COMPONENTS_0_10000_I");
		put("JAR_STATUS", "JAR_STATUS_YESNO");
		put("LOG_LEVEL", "LOG_LEVEL_n1_12_I");
		put("LONG_FORMAT", "LONG_FORMAT_YESNO");
		put("LRRSD_CUTOFF", "LRRSD_CUTOFF_0d0_1d0_D");
		put("MAX_MARKERS_LOADED_PER_CYCLE", "MAX_MARKERS_LOADED_PER_CYCLE_1_10000_I");
		put("MAX_MEMORY_USED_TO_LOAD_MARKER_DATA", "MAX_MEMORY_USED_TO_LOAD_MARKER_DATA_8_65536_I");
		put("NUM_THREADS", "NUM_THREADS_1_99_I");
		put("PARSE_AT_AT_SYMBOL", "PARSE_AT_AT_SYMBOL_YESNO");
		put("QQ_MAX_NEG_LOG10_PVALUE", "QQ_MAX_NEG_LOG10_PVALUE_1_10000_I");
		put("SAMPLE_ALIAS", "SAMPLE_ALIAS_LIST");
		put("TWOD_LOADED_VARIABLES", "TWOD_LOADED_VARIABLES_LIST");
		put("WINDOW_AROUND_SNP_TO_OPEN_IN_TRAILER", "WINDOW_AROUND_SNP_TO_OPEN_IN_TRAILER_100_1000000_I");
	}};

	public static void run(String in, String out) throws IOException {
		BufferedReader reader = Files.getAppropriateReader(in);
		PrintWriter writer = Files.getAppropriateWriter(out);
		
		StringBuilder sb = new StringBuilder();
		
		String line = null;
		while ((line = reader.readLine()) != null) {
			sb.append(line).append(System.getProperty("line.separator"));
		}
		String finalWrite = sb.toString();
		for (java.util.Map.Entry<String, String> replacement : replacements.entrySet()) {
			finalWrite = finalWrite.replace(replacement.getKey(), replacement.getValue());
		}
		writer.println(finalWrite);
		writer.flush();
		writer.close();
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String filenameSearch = "example.properties";
		String filenameReplace = "example_replaced.properties";

		String usage = "\n" + 
						"one.ben.PropertiesConversion requires 2 arguments\n" + 
						"   (1) filename to read (i.e. fileIn=" + filenameSearch + " (default))\n" + 
						"   (2) filename to write (i.e. fileOut=" + filenameReplace + " (default))\n" + 
						"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("fileIn=")) {
				filenameSearch = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("fileOut=")) {
				filenameReplace = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			run(filenameSearch, filenameReplace);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
}
