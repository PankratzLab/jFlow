package org.genvisis.seq.telomere;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import org.genvisis.common.Logger;
import org.genvisis.common.ext;

/**
 * Wrapper for the computel telomere length estimator ala
 * https://github.com/lilit-nersisyan/computel
 *
 *
 * Computel is a little annoying to set up, so this will reflect that...
 */
public class Computel {

	private Computel() {

	}

	private static void runComputel(String inputBam, String outputDir, String computelDirectory, Logger log) {
		String finalOutDirectory = outputDir + ext.rootOf(inputBam);
		new File(finalOutDirectory).mkdirs();
		try {
			copyDirectory(new File(computelDirectory), new File(finalOutDirectory));
			
		} catch (IOException e) {
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

	// bedtools bamtofastq -i aln.qsort.bam \
	// -fq aln.end1.fq \
	// -fq2 aln.end2.fq
}
