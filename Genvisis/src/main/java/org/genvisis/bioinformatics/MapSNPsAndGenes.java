// -Xms2048M -Xmx2048M
package org.genvisis.bioinformatics;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Vector;

import org.genvisis.common.Aliases;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.PSF;
import org.genvisis.common.Positions;
import org.genvisis.common.ProgressMonitor;
import org.genvisis.common.Sort;
import org.genvisis.common.ext;
import org.genvisis.seq.VCFExport;
import org.genvisis.seq.analysis.ANNOVAR;
import org.genvisis.seq.analysis.GATK_Genotyper;
import org.genvisis.seq.analysis.SNPEFF;

public class MapSNPsAndGenes {
	// public static final boolean AUTOMATICALLY_ADD_ONE = true;
	// public static final boolean AUTOMATICALLY_ADD_ONE = false;
	public static final String GENES36_FILENAME = "genes36.xln";
	public static final String GENES37_FILENAME = "genes37.xln";

	public static final int DEFAULT_WIGGLE_ROOM = 15000;
	public static final int UCSC_WINDOW = DEFAULT_WIGGLE_ROOM;

	public static String getGeneDB(int build, Logger log) {
		String geneDB;

		geneDB = GENES37_FILENAME;
		switch (build) {
			case 36:
				geneDB = GENES36_FILENAME;
				break;
			case 37:
				geneDB = GENES37_FILENAME;
				break;
			default:
				log.reportError("Error - unknown build '"	+ build
												+ "'; using the default instead (b37/hg19)");
				break;
		}

		return Aliases.getPathToFileInReferenceDirectory(geneDB, true, log);
	}

	public static String getSNPVCF(int build, Logger log) {
		String snpVCF = null;
		boolean satisfied;

		satisfied = false;
		while (!satisfied) {
			switch (build) {
				case 36:
					throw new IllegalArgumentException("Error - VCF unavailable for build 36.");
				case 37:
					snpVCF = ParseSNPlocations.DEFAULT_B37_VCF_FILENAME;
					satisfied = true;
					break;
				default:
					log.reportError("Error - unknown build '"	+ build + "'; using the default instead (build "
													+ ParseSNPlocations.DEFAULT_BUILD + ")");
					build = ParseSNPlocations.DEFAULT_BUILD;
					break;
			}
		}

		return Aliases.getPathToFileInReferenceDirectory(snpVCF, true, log);
	}

	public static String getSNPDB(int build, Logger log) {
		String snpDB = null;
		boolean satisfied;

		satisfied = false;
		while (!satisfied) {
			switch (build) {
				case 36:
					snpDB = ParseSNPlocations.DEFAULT_B36_DB_FILENAME;
					satisfied = true;
					break;
				case 37:
					snpDB = ParseSNPlocations.DEFAULT_B37_DB_FILENAME;
					satisfied = true;
					break;
				default:
					log.reportError("Error - unknown build '"	+ build + "'; using the default instead (build "
													+ ParseSNPlocations.DEFAULT_BUILD + ")");
					build = ParseSNPlocations.DEFAULT_BUILD;
					break;
			}
		}

		return Aliases.getPathToFileInReferenceDirectory(snpDB, true, log);
	}

	public static String getUnmappedVCF(Logger log) {
		return Aliases.getPathToFileInReferenceDirectory(	ParseSNPlocations.DEFAULT_UNMAPPED_VCF_FILENAME,
																											true, log);
	}

	public static String getMergeVCF(Logger log) {
		return Aliases.getPathToFileInReferenceDirectory(	ParseSNPlocations.DEFAULT_MERGE_VCF_FILENAME,
																											true, log);
	}

	public static String getMergeDB(Logger log) {
		return Aliases.getPathToFileInReferenceDirectory(	ParseSNPlocations.DEFAULT_MERGE_FILENAME, true,
																											log);
	}

	public static String[] mapSNPsToGenesLoosely(	int[][] markerPositions, int wiggleRoom, int build,
																								Logger log) {
		return Matrix.extractColumn(mapSNPsToGenes(markerPositions, build, wiggleRoom, log), 1);
	}

	public static String[][] mapSNPsToGenes(int[][] markerPositions, int build, int wiggleRoom,
																					Logger log) {
		return mapSNPsToGenes(markerPositions, getGeneDB(build, log), wiggleRoom, log);
	}

	/**
	 * Map chromosome positions to genes.
	 *
	 * @param markerPositions a two dimensional array of marker positions. An example of the format is
	 *        0 1 0 chr pos 1 chr pos ...
	 * @param geneDB the database containing the start and stop positions of the genes
	 * @param wiggleRoom the number of basepairs on either side of the position to search for genes
	 * @param log the Logger
	 * @return a matrix of genes String[0][] is the nearest match String[1][] are all matches
	 *         String[2][] are all of the exact matches
	 */
	public static String[][] mapSNPsToGenes(int[][] markerPositions, String geneDB, int wiggleRoom,
																					Logger log) {
		BufferedReader reader;
		String[] line, genes, distances;
		int[] chr_start_stop;
		int dist;
		String[][] finalGenes;
		String[] geneNames;
		int[] dists;
		int[] order;

		genes = Array.stringArray(markerPositions.length, "");
		distances = Array.stringArray(markerPositions.length, "");
		try {
			reader = new BufferedReader(new FileReader(geneDB));
			reader.readLine();
			chr_start_stop = new int[3];
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				if (!line[3].equals(".")) {
					chr_start_stop[0] =
														line[2].equals("X")	? 23
																								: (line[2].equals("Y")	? 24
																																				: (line[2].equals("XY")	? 25
																																																: (line[2].equals("MT")	? 26
																																																												: (line[2].equals("Un")	? 27
																																																																								: Integer.parseInt(line[2])))));
					chr_start_stop[1] = Integer.parseInt(line[3]);
					chr_start_stop[2] = Integer.parseInt(line[4]);

					for (int j = 0; j < markerPositions.length; j++) {
						if (markerPositions[j][0] == chr_start_stop[0]
									&& markerPositions[j][1] >= chr_start_stop[1] - wiggleRoom
								&& markerPositions[j][1] <= chr_start_stop[2] + wiggleRoom) {
							genes[j] += (genes[j].equals("") ? "" : "|") + line[1];
							if (markerPositions[j][1] < chr_start_stop[1]) {
								dist = (line[5].equals("+") ? -1 : 1) * (chr_start_stop[1] - markerPositions[j][1]);
							} else if (markerPositions[j][1] > chr_start_stop[2]) {
								dist = (line[5].equals("+") ? 1 : -1) * (markerPositions[j][1] - chr_start_stop[2]);
							} else {
								dist = 0;
							}

							// distances[j] += (distances[j].equals("")?"":"|")+dist+","+line[5];
							distances[j] += (distances[j].equals("") ? "" : "|") + dist;
						}
					}
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + geneDB + "\" not found");
			System.exit(1);
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + geneDB + "\"");
			System.exit(2);
		}

		finalGenes = new String[markerPositions.length][3];
		for (int i = 0; i < markerPositions.length; i++) {
			if (genes[i].equals("")) {
				finalGenes[i][0] = ".";
				finalGenes[i][1] = ".";
				finalGenes[i][2] = ".";
			} else {
				geneNames = genes[i].split("\\|");
				dists = Array.toIntArray(distances[i].split("\\|"));
				for (int j = 0; j < geneNames.length; j++) {
					if (dists[j] < 0) {
						geneNames[j] += " " + ext.prettyUpDistance(-1 * dists[j], 1) + " upstream";
					}
					if (dists[j] > 0) {
						geneNames[j] += " " + ext.prettyUpDistance(dists[j], 1) + " downstream";
					}
				}
				for (int j = 0; j < dists.length; j++) {
					dists[j] = Math.abs(dists[j]);
				}
				order = Sort.getSortedIndices(dists);
				for (int j = 0; j < geneNames.length; j++) {
					if (dists[j] == 0) {
						finalGenes[i][0] = finalGenes[i][0] == null	? geneNames[j]
																												: finalGenes[i][0] + "|" + geneNames[j];
						finalGenes[i][1] = finalGenes[i][1] == null	? geneNames[j]
																												: finalGenes[i][1] + "|" + geneNames[j];
					}
				}
				if (finalGenes[i][1] == null) {
					finalGenes[i][1] = geneNames[order[0]];
				}
				finalGenes[i][2] = Array.toStr(Sort.getOrdered(geneNames, order), "|");
			}
		}

		return finalGenes;
	}


	public static void procSNPsToGenes(	String dir, String snps, int wiggleRoom, int build, Logger log,
																			boolean useVCF, boolean snpEff, boolean gatk,
																			String snpEffLoc, String annovarLoc, String swapFile,
																			boolean exportToXLN) {
		PrintWriter writer;
		String[] line;
		String[] data, markers;
		String[][] genes;
		int[][] markerPositions;

		ProgressMonitor monitor = new ProgressMonitor(null, log);
		if (useVCF) {
			System.out.println("Processing with VCF files...");
			ParseSNPlocations.parseSNPlocations(dir	+ snps, getSNPVCF(build, log), getUnmappedVCF(log),
																					getMergeVCF(log), log, monitor);
		} else {
			ParseSNPlocations.lowMemParse(dir + snps, getSNPDB(build, log), getMergeDB(log), true, log);
		}

		String output = null;
		if (gatk) {
			String input = prepInput(dir + snps);
			if (input == null) {
				// TODO error occurred!
			}
			output = GATK_Genotyper.annotateOnly(	input, "", "", PSF.Ext.DEFAULT_MEMORY_MB, snpEffLoc,
																						snpEffLoc, annovarLoc, SNPEFF.BUILDS[0], true, false,
																						log);
		} else if (snpEff) {
			output = SNPEffAnnotation.pipeline(ext.rootOf(dir + snps, false)	+ "_positions.xln",
																					SNPEffAnnotation.getDefaultConfigFile(), log);
		}
		if (output != null && !"".equals(output) && Files.exists(output)) {
			if (swapFile != null && !"".equals(swapFile) && Files.exists(swapFile)) {
				Files.replaceAll(output, ext.rootOf(output) + "_converted.vcf", swapFile, log);
				output = ext.rootOf(output) + "_converted.vcf";
			}
			if (exportToXLN) {
				try {
					VCFExport.exportToXLN(output);
				} catch (IOException e) {
					log.reportIOException(output);
					log.reportException(e);
				}
			}
		}

		data = Array.toStringArray(HashVec.loadFileToVec(ext.rootOf(dir + snps, false)
																											+ "_positions.xln", false, false, false));

		ArrayList<String> mkrList = new ArrayList<String>();
		ArrayList<int[]> posList = new ArrayList<int[]>();
		for (int i = 0; i < data.length - 1; i++) {
			line = data[i + 1].trim().split("[\\s]+");
			if (!ext.isMissingValue(line[1]) && !ext.isMissingValue(line[2])) {
				mkrList.add(line[0]);
				posList.add(new int[] {Positions.chromosomeNumber(line[1]), Integer.parseInt(line[2])});
			}
		}
		markers = mkrList.toArray(new String[0]);
		markerPositions = posList.toArray(new int[0][]);

		genes = mapSNPsToGenes(markerPositions, getGeneDB(build, log), wiggleRoom, log);

		try {
			System.out.println(dir + ext.rootOf(snps) + "_genes.xln");
			writer = new PrintWriter(new FileWriter(dir + ext.rootOf(snps) + "_genes.xln"));
			writer.println("SNP\tChr\tPosition\tGene(s)\t"	+ UCSC_WINDOW + "\tClosest"
											+ "\t<- dynamically linked basepair buffer in UCSC hyperlink");
			for (int i = 0; i < markers.length; i++) {
				// writer.println(markers[i]+"\t"+markerPositions[i][0]+"\t"+(AUTOMATICALLY_ADD_ONE?markerPositions[i][1]+1:markerPositions[i][1])+"\t"+(genes[i].equals("")?".":genes[i])+"\t=HYPERLINK(CONCATENATE(\"http://genome.ucsc.edu/cgi-bin/hgTracks?position=chr"+markerPositions[i][0]+":\","+markerPositions[i][1]+"-$E$1+1,\"-\","+markerPositions[i][1]+"+$E$1),CONCATENATE(\"chr"+markerPositions[i][0]+":\","+markerPositions[i][1]+"-$E$1+1,\"-\","+markerPositions[i][1]+"+$E$1))");
				writer.println(markers[i]	+ "\t" + markerPositions[i][0] + "\t" + markerPositions[i][1]
												+ "\t" + (genes[i][0] == null ? "." : genes[i][0]) + "\t"
												+ (genes[i][1] == null ? "." : genes[i][1]) + "\t"
												+ (genes[i][2] == null ? "." : genes[i][2])
												+ "\t=HYPERLINK(CONCATENATE(\"http://genome.ucsc.edu/cgi-bin/hgTracks?position=chr"
												+ markerPositions[i][0] + ":\"," + markerPositions[i][1] + "-$E$1+1,\"-\","
												+ markerPositions[i][1] + "+$E$1),CONCATENATE(\"chr" + markerPositions[i][0]
												+ ":\"," + markerPositions[i][1] + "-$E$1+1,\"-\"," + markerPositions[i][1]
												+ "+$E$1))");
			}
			writer.flush();
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing genes file");
			log.reportException(e);
		}
	}

	private static String prepInput(String origFile) {
		BufferedReader reader;
		PrintWriter writer;

		String input = ext.rootOf(origFile, false) + "_positions.xln";
		String output = ext.rootOf(origFile, false) + "_gatk.vcf";

		boolean complete = false;
		try {
			reader = Files.getAppropriateReader(input);
			writer = Files.getAppropriateWriter(output);

			reader.readLine();
			writer.println("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
			String line = null;
			while ((line = reader.readLine()) != null) {
				String[] parts = line.split("\t");
				writer.println(parts[1]	+ "\t" + parts[2] + "\t" + parts[0] + "\t"
												+ (parts.length >= 5 ? parts[3] : ".") + "\t"
												+ (parts.length >= 5 ? parts[4] : ".") + "\t.\t.\t.");
			}
			writer.flush();
			writer.close();
			reader.close();
			complete = true;
		} catch (IOException e) {
			e.printStackTrace();
		}
		return !complete ? null : output;
	}

	public static void fromParameters(String filename, Logger log) {
		Vector<String> params;

		String snpEffLoc = Aliases.getPathToFileInReferenceDirectory(	SNPEFF.SNP_EFF, true,
																																	new Logger());
		if (snpEffLoc == null) {
			snpEffLoc = Aliases.getPathToFileInReferenceDirectory("snpEff/"	+ SNPEFF.SNP_EFF, true,
																														new Logger());
		}
		if (snpEffLoc == null) {
			snpEffLoc = "";
		}
		String annovarLoc = Aliases.getPathToFileInReferenceDirectory(ANNOVAR.TABLE_ANNOVAR, true,
																																	new Logger());
		if (annovarLoc == null) {
			annovarLoc = Aliases.getPathToFileInReferenceDirectory("ANNOVAR/annovar/"
																																+ ANNOVAR.TABLE_ANNOVAR, true,
																															new Logger());
		}
		if (annovarLoc == null) {
			annovarLoc = "";
		}

		params = Files.parseControlFile(filename, "snps",
																		new String[] {"file=list.snps", "dir=", "win=15000", "build=37",
																									"vcf=true", "snpeff=false", "annotate=false",
																									"snpeffLoc=" + snpEffLoc,
																									"annovarLoc=" + annovarLoc,},
																		log);

		if (params != null) {
			params.add("log=" + log.getFilename());
			main(Array.toStringArray(params));
		}
	}


	public static void main(String[] args) {
		int numArgs = args.length;
		String dir = "";
		String filename = "list.snps";
		int wiggleRoom = 15000;
		String temp;
		byte build = 37;
		Logger log;
		boolean vcf = true;
		String snpEffLoc = Aliases.getPathToFileInReferenceDirectory("snpEff/"	+ SNPEFF.SNP_EFF, false,
																																	new Logger());
		if (snpEffLoc == null) {
			snpEffLoc = "";
		}
		String annovarLoc = Aliases.getPathToFileInReferenceDirectory("ANNOVAR/annovar/"
																																		+ ANNOVAR.TABLE_ANNOVAR, false,
																																	new Logger());
		if (annovarLoc == null) {
			annovarLoc = "";
		}
		boolean snpeff = false;
		boolean gatk = false;
		String swap = null;
		boolean xln = false;
		String logfile = null;

		String usage = "\n"	+ "bioinformatics.MapSNPsAndGenes requires 0-1 arguments\n"
										+ "   (1) directory (i.e. dir=" + dir + " (default))\n"
										+ "   (2) filename (i.e. file=" + filename + " (default))\n"
										+ "   (3) # bp up and down stream to count as an associated gene (i.e. win="
										+ wiggleRoom + " (default))\n"
										+ "   (4) build # of the NCBI gene map file (i.e. build=" + build
										+ " (default))\n"
										+ "   (5) should use vcf files instead of serialized database files (i.e. vcf=true (default))\n"
										+
										// " (6) create additional output file with annotations from SNPEFF (i.e.
										// snpeff=true (not the default; mutually-exclusive with 'gatk' option))\n" +
										"   (6) create additional output file with annotations from GATK, SNPEFF, and ANNOVAR (i.e. annotate=true (not the default; mutually-exclusive with 'snpeff' option))\n"
										+ "   (7) Location of SNPEFF program in filesystem; used if 'gatk' option set to TRUE (i.e. snpeffLoc="
										+ snpEffLoc + " (default))\n"
										+ "   (8) Location of ANNOVAR program in filesystem; used if 'gatk' option set to TRUE (i.e. annovarLoc="
										+ annovarLoc + " (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("dir=")) {
				dir = ext.parseStringArg(arg, "");
				numArgs--;
			} else if (arg.startsWith("file=")) {
				filename = ext.parseStringArg(arg, null);
				numArgs--;
			} else if (arg.startsWith("win=")) {
				wiggleRoom = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("build=")) {
				build = ext.parseByteArg(arg);
				numArgs--;
			} else if (arg.startsWith("vcf=")) {
				vcf = ext.parseBooleanArg(arg);
				numArgs--;
			} else if (arg.startsWith("snpeff=")) {
				snpeff = ext.parseBooleanArg(arg);
				numArgs--;
			} else if (arg.startsWith("annotate=")) {
				gatk = ext.parseBooleanArg(arg);
				numArgs--;
			} else if (arg.startsWith("snpeffLoc=")) {
				snpEffLoc = ext.parseStringArg(arg, null);
				numArgs--;
			} else if (arg.startsWith("annovarLoc=")) {
				annovarLoc = ext.parseStringArg(arg, null);
				numArgs--;
			} else if (arg.startsWith("swap=")) {
				swap = ext.parseStringArg(arg, null);
				numArgs--;
			} else if (arg.startsWith("xln=")) {
				xln = ext.parseBooleanArg(arg);
				numArgs--;
			} else if (arg.startsWith("log=")) {
				logfile = ext.parseStringArg(arg, null);
				numArgs--;
			} else {
				System.err.println("Error - don't know what to do with argument '" + arg + "'");
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			return;
		}

		try {
			if (filename.endsWith(".snps")) {
				log = new Logger(filename + ".log");
				temp = filename.substring(0, filename.length() - (".snps").length());
				if (filename.indexOf(".b36.") > 0) {
					build = 36;
				} else if (filename.indexOf(".b37.") > 0) {
					build = 37;
				} else {
					log.reportError("Warning - using the default build (b37/hg19) since the file '"	+ filename
													+ "' does not explicitly specify one using the convention \"filename.b37.20000.snps\"");
					build = 37;
				}
				if (temp.indexOf(".") > 0) {
					try {
						wiggleRoom = Integer.parseInt(temp.substring(temp.lastIndexOf(".") + 1));
					} catch (NumberFormatException nfe) {
						wiggleRoom = DEFAULT_WIGGLE_ROOM;
					}
				}
				procSNPsToGenes(dir, filename, wiggleRoom, build, log, vcf, snpeff, gatk, snpEffLoc,
												annovarLoc, swap, xln);
			} else {
				log = new Logger(logfile);
				procSNPsToGenes(dir, filename, wiggleRoom, build, log, vcf, snpeff, gatk, snpEffLoc,
												annovarLoc, swap, xln);
			}
		} catch (Exception e) {
			e.printStackTrace();
			if (filename.endsWith(".snps")) {
				ext.waitForResponse();
			}
		}
	}

}
