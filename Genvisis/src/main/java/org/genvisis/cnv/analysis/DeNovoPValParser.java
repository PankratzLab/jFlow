package org.genvisis.cnv.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import org.genvisis.CLI;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.ext;

/**
 * This utility takes the .trio.cnv output from {@link cnvTrio} and the .log outputs from PennCNV
 * denovo analysis (http://penncnv.openbioinformatics.org/en/latest/user-guide/denovo/) and outputs
 * a .scored.cnv that includes all output from .trio.cnv with a corresponding denovo entry, and
 * tacks on columns for number of markers with father origin, mother origin and the called p-value.
 */
public class DeNovoPValParser {
	public static void main(String... args) {
		final String trioOut = "trio";
		final String denovoDir = "denovo";

		CLI c = new CLI(PennCNVFamilies.class);
		c.addArg(trioOut, "A .trio.cnv file from running the cnvTrio QC metrics", true);
		c.addArg(denovoDir,
		         "A directory containing .gen and .log files from the de novo PennCNV output", true);

		c.parseWithExit(args);

		mergeResults(c.get(trioOut), c.get(denovoDir));
	}

	/**
	 * Read each line of the trio output. Use the FID, IID and CHR, BP1 and BP2 to match up .log files
	 * in the denovo directory. Add columns for the number of regions from each parent and p-value,
	 * and write the output to a new .score.cnv file.
	 */
	private static void mergeResults(String trioOut, String denovoDir) {
		try {
			BufferedReader reader = Files.getAppropriateReader(trioOut);
			PrintWriter writer = getWriter(trioOut);
			Map<String, File> logs = getLogs(denovoDir);

			String[] trioHeader = reader.readLine().split("\t");
			String[] outHeader = Arrays.copyOf(trioHeader, trioHeader.length + 3);
			outHeader[trioHeader.length] = "FA_ORIGIN";
			outHeader[trioHeader.length+1] = "MO_ORIGIN";
			outHeader[trioHeader.length+2] = "P-Val";
			writer.println(Array.toStr(outHeader, "\t"));
			int[] idxs = ext.indexFactors(new String[]{"FID", "IID", "CHR", "BP1", "BP2"}, trioHeader, true, false);
			String[] colsOfInterest = new String[]{"Paternal_origin(F)=", "Maternal_origin(M)=", "P-value="};

			while (reader.ready()) {
				String[] line = reader.readLine().split("\t");
				String key = line[idxs[0]];
				for (int i=1; i<idxs.length; i++) {
					key += "_" + line[idxs[i]];
				}
				key += ".log";
				File f = logs.get(key);
				if (f != null) {
					String[] lineOut = Arrays.copyOf(line, line.length + 3);
					Arrays.fill(lineOut, line.length, lineOut.length, ".");
					BufferedReader dnvReader = Files.getAppropriateReader(f.getAbsolutePath());
					boolean foundPval = false;
					while (dnvReader.ready() && !foundPval) {
						String dnvLine = dnvReader.readLine();
						if (dnvLine.contains(colsOfInterest[0])) {
							String[] dnvSplit = dnvLine.split(" ");
							int[] outIdxs = ext.indexFactors(colsOfInterest, dnvSplit, true, false);
							for (int i=0; i<outIdxs.length; i++) {
								lineOut[line.length + i] = dnvSplit[outIdxs[i] + 1];
							}
							foundPval = true;
						}
					}
					writer.println(Array.toStr(lineOut, "\t"));
					dnvReader.close();
				}
			}

			writer.close();
			reader.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	/**
	 * Create a map of file name (format: FID_IID_CHR_BP1_BP2.log) to file.
	 */
	private static Map<String, File> getLogs(String denovoDir) {
		Map<String, File> m = new HashMap<String, File>();
		File[] files = new File(denovoDir).listFiles(new FilenameFilter() {
			
			@Override
			public boolean accept(File dir, String name) {
				return name.endsWith(".log");
			}
		});

		for (File f : files) {
			m.put(f.getName(), f);
			
		}

		return m;
	}

	private static PrintWriter getWriter(String trioOut) throws IOException {
		File f = new File(trioOut);
		String name = f.getName();
		name = name.substring(0, name.indexOf('.'));
		String parent = f.getCanonicalFile().getParent();
		String outFile = (parent == null ? "" : parent + File.separator) + name + ".score.cnv";
		System.out.println("Writing unified denovo cnv output to: " + outFile);
		return Files.getAppropriateWriter(outFile);
	}
}
