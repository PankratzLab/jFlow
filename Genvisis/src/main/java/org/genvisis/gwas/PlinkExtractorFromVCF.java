package org.genvisis.gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;

import org.genvisis.common.Files;
import org.genvisis.common.Matrix;
import org.genvisis.common.Positions;
import org.genvisis.common.ext;

public class PlinkExtractorFromVCF {

	private static int[][] loadRegions(String rgns) {
		BufferedReader reader;
		String line;
		ArrayList<int[]> pos;
		int[] temp;
		int lineIndex;

		pos = new ArrayList<int[]>();
		line = null;
		lineIndex = 0;
		try {
			reader = Files.getAppropriateReader(rgns);
			while ((line = reader.readLine()) != null) {
				lineIndex++;
				temp = Positions.parseUCSClocation(line);
				if (temp == null || temp.length < 3 || temp[0] == -1 || temp[1] == -1 || temp[2] == -1) {
					System.err.println("Error - invalid UCSC region at line " + lineIndex + ": " + line);
					continue;
				}
				pos.add(temp);
			}
		} catch (IOException e) {
			e.printStackTrace();
			return new int[0][];
		}

		return pos.toArray(new int[pos.size()][]);
	}

	// assuming .vcf.gz paired with .vcf.gz.tbi files
	private static String[] get1000GFiles(final String srcDir, int[] chrs) {
		String[] files, results;
		String search;

		files = (new File(srcDir)).list(new FilenameFilter() {
			@Override
			public boolean accept(File dir, String name) {
				boolean hasChr = name.contains("chr");
				if (!hasChr) {
					return false;
				}
				boolean format = name.endsWith(".vcf.gz");
				if (!format) {
					return false;
				}
				boolean hasTbi = Files.exists(srcDir + name + ".tbi");
				return hasTbi;
			}
		});
		if (files.length == 0) {
			System.err.println("Error - no valid data files found!");
			return null;
		}
		System.out.println("Found " + files.length + " valid data files");

		results = new String[chrs.length];
		for (int i = 0; i < chrs.length; i++) {
			search = ".chr" + chrs[i] + ".";
			for (String file : files) {
				if (file.contains(search)) {
					results[i] = file;
					break;
				}
			}
			if (results[i] == null) {
				System.err.println("Error - unable to find 1000G data file for chromosome "	+ chrs[i]
														+ ".  Please ensure that data files names include \".chr#.\".");
			}
		}

		return results;
	}


	// plink2 --noweb --vcf <> --chr <> --from-bp <> --to-bp <> --make-bed --out
	public static void run(String dir, String rgns, String ids) {
		StringBuilder sb;
		String plinkRoot;

		String current = ext.verifyDirFormat((new File("./")).getAbsolutePath());
		(new File(current + "temp/")).mkdirs();

		int[][] ucsc = loadRegions(rgns);
		int[] chrs = Matrix.extractColumn(ucsc, 0);
		String[] dataFiles = get1000GFiles(dir, chrs);
		if (dataFiles == null) {
			return;
		}

		ArrayList<String> commands = new ArrayList<String>();
		ArrayList<String> roots = new ArrayList<String>();

		commands.add("cd " + current);
		for (int i = 0; i < ucsc.length; i++) {
			if (dataFiles[i] == null) {
				continue;
			}
			plinkRoot = "plink" + ucsc[i][0] + ucsc[i][1] + ucsc[i][2];
			sb = new StringBuilder("plink2 --noweb").append(" --vcf ").append(dir).append(dataFiles[i])
																							.append(" --chr ").append(ucsc[i][0])
																							.append(" --from-bp ").append(ucsc[i][1])
																							.append(" --to-bp ").append(ucsc[i][2]);
			if (ids != null) {
				sb.append(" --keep ").append(ids);
			}
			sb.append(" --make-bed --out ").append(current + "temp/" + plinkRoot);
			commands.add(sb.toString());
			roots.add("temp/" + plinkRoot);
		}

		String mergeRoot = roots.remove(0);
		Files.writeIterable(roots, current + "temp/mergeList.txt");

		String plinkMerge = "plink2 --noweb --bfile "	+ current + mergeRoot + " --merge-list " + current
												+ "temp/mergeList.txt --make-bed --out " + current + "merged";
		commands.add("");
		commands.add(plinkMerge);
		commands.add("");
		commands.add("rm -r " + current + "temp/");

		Files.writeIterable(commands, current + "run.sh");
		Files.chmod(current + "run.sh");
		Files.makeQsub(current + "run.sh", false, 1, 1, false, null, false);
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String dir = "";
		String rgnFile = "regions.txt";
		String idsFile = null;

		String usage = "\n"	+ "gwas.PlinkExtractorFromVCF requires 0-1 arguments\n"
										+ "   (1) Directory of 1000G VCF files (i.e. dir=" + dir + " (default))\n"
										+ "   (2) UCSC regions filename (i.e. regions=" + rgnFile + " (default))\n"
										+ "   (3) (Optional) IDs list filename (i.e. ids=" + idsFile + " (default))\n"
										+ "" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("dir=")) {
				dir = ext.verifyDirFormat(arg.split("=")[1]);
				numArgs--;
			} else if (arg.startsWith("regions=")) {
				rgnFile = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("ids=")) {
				idsFile = arg.split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + arg);
			}
		}
		if (args.length == 0 || numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}

		run(dir, rgnFile, idsFile);
	}


}
