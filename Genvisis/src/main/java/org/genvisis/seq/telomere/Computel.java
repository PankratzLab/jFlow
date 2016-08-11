package org.genvisis.seq.telomere;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Scanner;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.Array;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.seq.manage.BamOps;

/**
 * Wrapper for the computel telomere length estimator ala
 * https://github.com/lilit-nersisyan/computel
 *
 *
 * Computel is a little annoying to set up, so this will reflect that...
 */
public class Computel {

	private static final String SAM_TO_FASTQ_LOC = "SamToFastq.jar";

	private Computel() {

	}

	private static boolean convertToFasta(String inputBam, String computelLoc, String r1, String r2, Logger log) {
		String[] inputs = new String[] { inputBam };
		String[] outputs = new String[] { r1, r2 };
		ArrayList<String> command = new ArrayList<String>();
		command.add("java");
		command.add("-jar");
		command.add(computelLoc + SAM_TO_FASTQ_LOC);
		command.add("I=");
		command.add(inputBam);
		command.add("F=");
		command.add(r1);
		command.add("F2=");
		command.add(r2);

		return CmdLine.runCommandWithFileChecks(Array.toStringArray(command), "", inputs, outputs, true, false, false,
				log);
	}

	private static void runComputel(String inputBam, String outputDir, String computelDirectory, Logger log) {
		String finalOutDirectory = outputDir + ext.rootOf(inputBam) + "/";
		new File(finalOutDirectory).mkdirs();
		try {
			copyDirectory(new File(computelDirectory), new File(finalOutDirectory));

		} catch (IOException e) {
			log.reportException(e);
		}

		try {
		
			String r1 = finalOutDirectory + "src/examples/analysis_reads1.fq";
			String r2 = finalOutDirectory + "src/examples/analysis_reads2.fq";
			boolean converted = convertToFasta(inputBam, finalOutDirectory, r1, r2, log);
			System.out.println(converted);
			if (converted) {
				int readLength = BamOps.estimateReadSize(inputBam, 100000, log);
				System.out.println(readLength);

				Scanner s = new Scanner(new File(finalOutDirectory + "src/examples/config_unix.txt"));
				
				String config = s.useDelimiter("\\Z").next();
				s.close();
				System.out.println(config);
				log.report(config);

			}
		} catch (FileNotFoundException e) {
			log.reportException(e);
		}

	}

	private static void copyDirectory(File sourceLocation, File targetLocation) throws IOException {

		if (sourceLocation.isDirectory()) {
			if (!targetLocation.exists()) {
				targetLocation.mkdir();
			}

			String[] children = sourceLocation.list();
			for (int i = 0; i < children.length; i++) {
				copyDirectory(new File(sourceLocation, children[i]), new File(targetLocation, children[i]));
			}
		} else {

			InputStream in = new FileInputStream(sourceLocation);
			OutputStream out = new FileOutputStream(targetLocation);

			// Copy the bits from instream to outstream
			byte[] buf = new byte[1024];
			int len;
			while ((len = in.read(buf)) > 0) {
				out.write(buf, 0, len);
			}
			in.close();
			out.close();
		}
	}

	public static void test(String bam, String computelLocation, String outDir) {
		runComputel(bam, outDir, computelLocation, new Logger(outDir + "computel.Log"));
	}

	public static void main(String[] args) {
		String targetBam = "/Volumes/Beta/data/aric_sra/bams/SRR1654226.bam";
		String computelLocation = "/Users/Kitty/git/computel/";
		String outDir = "/Volumes/Beta/data/aric_sra/test/computel/";
		new File(outDir).mkdirs();
		test(targetBam, computelLocation, outDir);
	}

	// bedtools bamtofastq -i aln.qsort.bam \
	// -fq aln.end1.fq \
	// -fq2 aln.end2.fq
}
