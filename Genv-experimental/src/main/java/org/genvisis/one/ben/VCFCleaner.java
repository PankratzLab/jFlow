package org.genvisis.one.ben;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;

import org.genvisis.common.Files;
import org.genvisis.common.ext;

public class VCFCleaner {


	public static String fixVCF(String vcfFile) throws IOException {
		BufferedReader reader = Files.getAppropriateReader(vcfFile);
		String outFile = ext.rootOf(vcfFile, false) + "_fixed.vcf"
										 + (vcfFile.endsWith(".gz") ? ".gz" : "");
		PrintWriter writer = Files.getAppropriateWriter(outFile);
		String line = null;
		boolean pastChrome = false;
		StringBuilder sb;
		while ((line = reader.readLine()) != null) {
			if (pastChrome) {
				line = line.replace("\t\t\t", "\t").replaceAll("\t\t", "\t");
				int ind = line.indexOf(",");
				if (ind > 0) {
					sb = new StringBuilder(line.substring(0, ind));
					sb.append(line.substring(line.indexOf("\t", ind)));
					line = sb.toString();
				}
			} else if (line.startsWith("FORMAT")) {
				line = fixFormatLine(line);
			} else if (line.startsWith("#CHROM")) {
				pastChrome = true;
			}
			writer.println(line);
		}
		writer.flush();
		writer.close();
		reader.close();
		return outFile;
	}

	private static String fixFormatLine(String formatLine) {
		String[] parts = formatLine.split(",");
		StringBuilder sb = new StringBuilder("##");
		int typeInd = -1;
		for (int i = 0; i < parts.length; i++) {
			if (parts[i].startsWith("Type=")) {
				typeInd = i;
			} else if (parts[i].startsWith("Number=")) {
				sb.append(parts[i]).append(",").append(parts[typeInd]);
				if (i < parts.length - 1) {
					sb.append(",");
				}
			} else {
				sb.append(parts[i]);
				if (i < parts.length - 1) {
					sb.append(",");
				}
			}
		}
		return sb.toString();
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "chr1.vcf.gz";

		String usage = "\n" +
									 "org.genvisis.one.ben.VCFFixer requires 0-1 arguments\n" +
									 "   (1) filename (i.e. file=" + filename + " (default))\n" +
									 "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h")
					|| args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
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
			fixVCF(filename);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
