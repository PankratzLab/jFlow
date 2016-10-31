package org.genvisis.seq;

import java.io.BufferedReader;
import java.io.EOFException;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.bioinformatics.SeattleSeq;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.Positions;
import org.genvisis.common.Sort;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;
import org.genvisis.filesys.SegmentLists;

public class Polymutt {
	public static final String[] POLYMUTT_VCF_HEADER = {"#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
																											"FILTER", "INFO", "FORMAT", "2", "3", "1"};
	public static final String[] COMMON_VCF_LOCATIONS = {	"/lustre/polymutt/",
																												"/project/scratch/polymutt/",
																												"c:/scratch/polymutt/",
																												"D:/Logan/Polymutt/results/", "./"};
	public static final String[] GENOTYPE_ORDER = {	"A/A", "A/C", "A/G", "A/T", "C/C", "C/G", "C/T",
																									"G/G", "G/T", "T/T"};

	public static void batchAllGlf(String sourceDir) {
		String[] files;
		String root, dir, filename;
		Vector<String> v = new Vector<String>();

		if (!sourceDir.startsWith("/")) {
			System.err.println("Error - sourceDirectory needs to explicit");
			return;
		}

		dir = ext.pwd();
		new File(dir + "batches/").mkdirs();
		sourceDir = ext.verifyDirFormat(sourceDir);
		v = new Vector<String>();
		files = Files.list(sourceDir, ".bam", false);
		for (String file : files) {
			root = ext.rootOf(file);
			filename = root + "_batch";
			Files.write("~/bin/samtools-hybrid pileup "	+ sourceDir + root
									+ ".bam -g -f ~/bin/ref/hg19_canonical.fa > " + dir + root + ".bam.glf",
									dir + "batches/" + filename);
			Files.chmod(dir + "batches/" + filename);
			v.add(dir + "batches/" + filename);
		}
		// Files.qsubMultiple(v, null, dir+"batches/", "glfAll", 16, true, "sb", -1, 62000, 2);
		Files.qsubMultiple(v, null, dir + "batches/", "glfAll", 8, true, null, -1, 62000, 2);
	}

	private static void batchPolymutt(String triosFile) {
		String[][] iterations;
		String commands, scratch;
		Vector<String> v;

		scratch = "/lustre/polymutt/";
		scratch = "/project/scratch/polymutt/";

		v = new Vector<String>();
		Files.write("T\tGLF_Index", "same_for_all.dat");
		Files.writeArray(	new String[] {"fam	3	0	0	1	3", "fam	2	0	0	2	2", "fam	1	3	2	2	1"},
											"same_for_all.ped");
		iterations = HashVec.loadFileToStringMatrix(triosFile, false, new int[] {0, 1, 2}, false);
		for (String[] iteration : iterations) {
			Files.writeArray(	new String[] {"1 "	+ iteration[0] + ".bam.glf",
																			"2 " + iteration[1] + ".bam.glf",
																			"3 " + iteration[2] + ".bam.glf"},
												iteration[0] + ".gif");

			commands = "cd "	+ ext.pwd() + "\n"
									+ "~/bin/polymutt -p same_for_all.ped -d same_for_all.dat -g " + iteration[0]
									+ ".gif --nthreads 8  --out_vcf " + scratch + iteration[0]
									+ ".denovo.out.vcf --denovo --rate_denovo 1.5e-07\n" + "cd " + scratch + "\n"
									+ "more " + iteration[0]
									+ ".denovo.out.vcf | grep -v 'DQ=-0' | grep -v 'DQ=0.00' > " + iteration[0]
									+ ".vcf\n" + "gzip " + iteration[0] + ".denovo.out.vcf";

			Files.qsub("batches/run_" + iteration[0], commands, 22000, 48, 8);
			Files.chmod("batches/run_" + iteration[0]);
			v.add("qsub batches/run_" + iteration[0]);
		}
		v.insertElementAt("cd " + ext.pwd(), 0);
		Files.writeArray(Array.toStringArray(v), "master.allRuns");
		Files.chmod("master.allRuns");
	}

	public static void filterDenovo(String triosFile, String vcfDir, String annotationDir) {
		String[][] iterations;
		Logger log;
		Hashtable<String, String[]> parsedAnnotations;
		Vector<String> unknownAnnotations, finishedAnnotations;
		String filename, temp;

		log = new Logger(ext.parseDirectoryOfFile(triosFile) + "parseAllPolymuttInDir.log");

		if (!Files.exists(vcfDir)) {
			System.err.println("Error - vcf directory not found: " + vcfDir);
			return;
		}

		if (annotationDir == null) {
			log.reportError("annotationDir directory set to the cwd");
			annotationDir = "./";
		} else if (!Files.exists(annotationDir)) {
			log.reportError("Warning - annotationDir directory not found: " + annotationDir);
			log.reportError("          annotationDir directory reset to the cwd");
		}

		annotationDir = ext.verifyDirFormat(annotationDir);

		if (Files.isDirectory(triosFile)) {
			log.reportError("Error - defined trios \"file\" is actually a directory");
			return;
		} else if (!Files.exists(triosFile)) {
			log.reportError("Error - trios file not found: " + triosFile);
			return;
		}

		unknownAnnotations = new Vector<String>();
		finishedAnnotations = new Vector<String>();
		parsedAnnotations = SeattleSeq.loadAllBadAnnotationInDir(annotationDir, log);
		iterations = HashVec.loadFileToStringMatrix(triosFile, false, new int[] {0, 1, 2}, false);
		for (String[] iteration : iterations) {
			if (Files.exists(vcfDir + iteration[0] + ".vcf")) {
				parseDenovo(vcfDir	+ iteration[0] + ".vcf", parsedAnnotations, unknownAnnotations,
										finishedAnnotations, log);
			} else {
				System.err.println("Error - '" + vcfDir + iteration[0] + ".vcf' does not exist");
			}
		}
		if (unknownAnnotations.size() > 0) {
			filename = "SeattleSeq_" + ext.getTimestampForFilename() + ".input";
			Files.writeArray(Array.toStringArray(unknownAnnotations), annotationDir + filename);
		}
		temp = "Filename\tSampleRoot\tMarkerName\tChr\tPosition\tREF\tALT\tQUAL\tMapQual\tDenovoQual";
		for (int i = 1; i <= 3; i++) {
			temp += "\tGeno" + i;
		}
		temp += "\tFunction\tinDbSNP";
		for (int i = 1; i <= 3; i++) {
			temp += "\tReadDepth" + i;
		}
		temp += "\tReadDepthTotal";
		for (int i = 1; i <= 3; i++) {
			temp += "\tGenQual" + i;
		}
		for (int i = 1; i <= 3; i++) {
			temp += "\tPhred" + i;
		}
		for (int i = 1; i <= 3; i++) {
			temp += "\tNextBestGeno" + i;
		}

		finishedAnnotations.insertElementAt(temp, 0);
		Files.writeArray(	Array.toStringArray(finishedAnnotations),
											ext.parseDirectoryOfFile(triosFile) + "summary.xln");
	}

	private static void parseDenovo(String filename, Hashtable<String, String[]> parsedAnnotations,
																	Vector<String> unknownAnnotations,
																	Vector<String> finishedAnnotations, Logger log) {
		BufferedReader reader;
		String[] line, subline;
		String temp, trav;
		String markerName, refAllele, denovoAllele, chr, position;
		String genotype, parentalGenotypes, genotypes;
		boolean denovo;
		double dqScore, qualityScore, mapQuality;
		int[] readDepths; // 1 2 3 total
		double[] genotypeQualities;
		int[] phreds, phredOfCall;
		String[] nextBestGenotype;
		int indexOfPhred;

		try {
			reader = new BufferedReader(new FileReader(filename));
			do {
				temp = reader.readLine();
			} while (!temp.startsWith("#CHROM"));

			if (!ext.checkHeader(temp.split("[\\s]+"), POLYMUTT_VCF_HEADER, false)) {
				log.reportError("Problem with header for file: " + filename);
				reader.close();
				return;
			}
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				chr = line[0];
				position = line[1];
				markerName = chr + ":" + position;
				refAllele = line[3];
				genotypes = "";
				parentalGenotypes = "";
				denovoAllele = "";
				for (int i = 0; i < 3; i++) {
					genotype = ext.replaceAllWith(line[9 + i].split(":")[0], "/", "");
					if (i < 2) {
						parentalGenotypes += genotype;
					} else {
						for (int j = 0; j < 2; j++) {
							denovo = true;
							for (int k = 0; k < parentalGenotypes.length(); k++) {
								if (genotype.charAt(j) == parentalGenotypes.charAt(k)) {
									denovo = false;
								}
							}
							for (int k = 0; k < denovoAllele.length(); k++) {
								if (genotype.charAt(j) == denovoAllele.charAt(k)) {
									log.reportError("Warning - double denovo for "	+ markerName + " in "
																	+ ext.removeDirectoryInfo(filename));
									denovo = false;
								}
							}
							if (denovo) {
								denovoAllele += genotype.charAt(j);
							}
						}
					}
					genotypes += "\t" + genotype;
				}

				// double[] genotypeQualities, phredOfCall, nextBestPhred;
				// String[] phreds;

				qualityScore = Double.parseDouble(line[5]);


				subline = line[7].split(";");

				readDepths = new int[4];
				readDepths[3] = Integer.parseInt(subline[2].split("=")[1]);
				mapQuality = Double.parseDouble(subline[3].split("=")[1]);
				dqScore = Double.parseDouble(subline[4].split("=")[1]);

				if (!line[8].equals("GT:GQ:DP:PL")) {
					log.reportError("Error - column 9 is " + line[8] + " instead of GT:GQ:DP:PL");
				}

				phredOfCall = new int[3];
				nextBestGenotype = new String[3];
				genotypeQualities = new double[3];
				for (int i = 0; i < 3; i++) {
					subline = line[9 + i].split(":");
					phreds = Array.toIntArray(subline[3].split(","));
					indexOfPhred = ext.indexOfStr(subline[0], GENOTYPE_ORDER);
					phredOfCall[i] = phreds[indexOfPhred];
					phreds[indexOfPhred] = 999;
					nextBestGenotype[i] = GENOTYPE_ORDER[Sort.getSortedIndices(phreds)[0]];
					genotypeQualities[i] = Double.parseDouble(subline[1]);
					readDepths[i] = Integer.parseInt(subline[2]);
				}

				if (denovoAllele.length() > 1) {
					log.reportError("Warning - double heterozygote denovo for "	+ markerName + " in "
													+ ext.removeDirectoryInfo(filename));
				}

				for (int i = 0; i < denovoAllele.length(); i++) {
					trav = markerName + "_" + refAllele + "_" + denovoAllele.charAt(i);
					if (!parsedAnnotations.containsKey(trav)) {
						unknownAnnotations.add(trav	+ "\t" + chr + "\t" + position + "\t" + refAllele + "\t"
																		+ denovoAllele.charAt(i));
					} else {
						if (parsedAnnotations.get(trav).length > 0) {
							finishedAnnotations.add(ext.rootOf(filename)	+ "\t" + decrypt(ext.rootOf(filename))
																			+ "\t" + trav + "\t" + chr + "\t" + position + "\t"
																			+ refAllele + "\t" + denovoAllele.charAt(i) + "\t"
																			+ qualityScore + "\t" + mapQuality + "\t" + dqScore + ""
																			+ genotypes + "\t" + parsedAnnotations.get(trav)[0] + "\t"
																			+ Array.toStr(readDepths) + "\t"
																			+ Array.toStr(genotypeQualities) + "\t"
																			+ Array.toStr(phredOfCall) + "\t"
																			+ Array.toStr(nextBestGenotype));
						}
					}

				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + filename + "\" not found in current directory");
			return;
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + filename + "\"");
			return;
		}
	}

	public static String decrypt(String sample) {
		String str;

		str = sample;
		if (str.startsWith("dedup_") || str.startsWith("rrd_")) {
			str = str.substring(str.indexOf("_"));
		}
		if (str.indexOf("_L00") > 0) {
			str = str.substring(0, str.lastIndexOf("_"));
		}

		return str;
	}

	public static String[] getMostLikelyGenotypesFromPhredScores(String phredScores, Logger log) {
		Vector<String> mostLikelyGenotypes;
		String[] stringScores;
		int[] scores, order;
		int bestScore, count;

		stringScores = phredScores.split(",");
		scores = new int[stringScores.length];
		for (int i = 0; i < stringScores.length; i++) {
			try {
				scores[i] = Integer.parseInt(stringScores[i]);
			} catch (NumberFormatException nfe) {
				System.err.println("Error - failed to convert '"	+ stringScores[i]
														+ "' into an integer. Line: " + phredScores);
			}
		}

		if (scores.length != GENOTYPE_ORDER.length) {
			log.reportError("Error - mismatched number of phred scores. Found "	+ scores.length
											+ " when expecting " + GENOTYPE_ORDER.length + ". scores: " + phredScores);
		}

		mostLikelyGenotypes = new Vector<String>();
		order = Sort.getSortedIndices(scores);
		bestScore = scores[order[0]];
		count = 0;
		while (count < GENOTYPE_ORDER.length && scores[order[count]] == bestScore) {
			mostLikelyGenotypes.add(GENOTYPE_ORDER[order[count]]);
			count++;
		}
		if (count == GENOTYPE_ORDER.length) {
			if (Array.max(scores) > 0) {
				log.report("There was no likely genotype. Line: " + phredScores);
			}
		}

		return Array.toStringArray(mostLikelyGenotypes);
	}

	private static void findAllDenovo(String dir) {
		String[] filenames;

		filenames = Files.list(dir, ".vcf.gz", false);

		Files.qsub(	"findAll", ext.pwd(), -1,
								"java -cp ~/"							+ org.genvisis.common.PSF.Java.GENVISIS
																					+ " seq.Polymutt findDenovo=[%0]",
								Matrix.toMatrix(filenames), 1000, 12);
	}

	private static void findAllLocalDenovo(String dir) {
		String[] filenames;

		filenames = Files.list(dir, "vcf_denovo.vcf", false);

		for (String filename : filenames) {
			System.out.println(filename);
			findDenovo(dir + filename, 8);
		}
	}

	private static void findDenovo(String filename, int minReadDepth) {
		BufferedReader reader;
		PrintWriter[] writers;
		String[] line;
		String temp;
		String markerName, denovoAllele, chr, prevChr, position;
		String genotype, parentalGenotypes;
		String[] mostLikelyGenotypes;
		boolean denovo;
		Logger log;
		String outRoot;
		int[] readDepths;

		log = new Logger(ext.rootOf(filename, false) + ".log");

		try {
			reader = Files.getAppropriateReader(filename);
			new File(ext.parseDirectoryOfFile(filename) + "denovos/").mkdirs();
			outRoot = ext.parseDirectoryOfFile(filename) + "denovos/" + ext.rootOf(filename, true);
			writers = new PrintWriter[3];
			writers[0] = new PrintWriter(new FileWriter(outRoot + "_denovo.vcf"));
			writers[1] = new PrintWriter(new FileWriter(outRoot + "_ambiguousGenotypesInChild.vcf"));
			// writers[2] = new PrintWriter(new FileWriter(outRoot+"_lowCoverage.vcf"));
			do {
				temp = reader.readLine();
				writers[0].println(temp);
			} while (reader.ready() && !temp.startsWith("#CHROM"));

			if (!ext.checkHeader(temp.split("[\\s]+"), POLYMUTT_VCF_HEADER, false)) {
				log.reportError("Problem with header for file: " + filename);
				reader.close();
				Files.closeAll(writers);
				return;
			}
			prevChr = "";
			while (reader.ready()) {
				temp = reader.readLine();
				line = temp.trim().split("[\\s]+");
				chr = line[0];
				if (!chr.equals(prevChr)) {
					log.report(ext.getTime() + "\t" + chr);
					prevChr = chr;
				}

				position = line[1];
				markerName = chr + ":" + position;
				if (line.length < 11) {
					log.reportError("VCF file '" + filename + "' truncated at marker " + markerName);
					reader.close();
					Files.closeAll(writers);
					return;
				}
				parentalGenotypes = "";
				denovoAllele = "";
				readDepths = new int[3];
				for (int i = 0; i < 3; i++) {
					if (line[9 + i].split(":").length < 4
							|| line[9 + i].split(":")[3].split(",").length < 10) {
						log.reportError("VCF file '" + filename + "' truncated at marker " + markerName);
						reader.close();
						Files.closeAll(writers);
						return;
					}
					mostLikelyGenotypes =
															getMostLikelyGenotypesFromPhredScores(line[9 + i].split(":")[3], log);
					readDepths[i] = Integer.parseInt(line[9 + i].split(":")[2]);
					if (mostLikelyGenotypes.length == GENOTYPE_ORDER.length) {
						// writers[2].println(temp);
					}
					if (i == 2 && mostLikelyGenotypes.length > 1) {
						if (mostLikelyGenotypes.length == GENOTYPE_ORDER.length) {
							mostLikelyGenotypes = new String[] {line[9 + i].split(":")[0]};
						} else {
							// log.reportError("Warning - there are more than one 'most likely genotype' for child
							// in "+ext.rootOf(filename)+" for marker "+markerName+" so it will not be called
							// denovo. Line: "+line[9+i].split(":")[3]);
							writers[1].println(temp);
						}
					}
					for (String mostLikelyGenotype : mostLikelyGenotypes) {
						genotype = ext.replaceAllWith(mostLikelyGenotype, "/", "");

						if (i < 2) {
							parentalGenotypes += genotype;
						} else {
							for (int j = 0; j < 2; j++) {
								denovo = true;
								for (int k = 0; k < parentalGenotypes.length(); k++) {
									if (genotype.charAt(j) == parentalGenotypes.charAt(k)) {
										denovo = false;
									}
								}
								for (int k = 0; k < denovoAllele.length(); k++) {
									if (genotype.charAt(j) == denovoAllele.charAt(k)) {
										log.reportError("Warning - double denovo for "	+ markerName + " in "
																		+ ext.removeDirectoryInfo(filename));
										denovo = false;
									}
								}
								// denovo and reasonably high mapping quality and high genotype quality for kid and
								// at least 5 reads for all three individuals
								if (denovo	&& Double.parseDouble(line[7].split(";")[3].split("=")[1]) > 50
										&& Double.parseDouble(line[9 + 2].split(":")[1]) > 50
										&& Array.min(readDepths) >= minReadDepth) {
									writers[0].println(temp);
								}
							}
						}
					}
				}
			}
			reader.close();
			Files.closeAll(writers);
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + filename + "\" not found in current directory");
			return;
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + filename + "\"");
			return;
		}
	}

	public static void mergeAllPossible(String controlFile) {
		BufferedReader reader;
		String[] line;
		Vector<String> v;
		String pattern;
		String[] iterations;
		int[][] nCr;
		String command, fill;
		String dir;

		v = new Vector<String>();
		try {
			reader = Files.getAppropriateReader(controlFile);
			dir = ext.verifyDirFormat(reader.readLine());
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				pattern = line[0];
				iterations = Array.subArray(line, 1);

				for (String iteration : iterations) {
					fill = ext.replaceAllWith(pattern, "[%%]", iteration);
					if (!Files.exists(dir + fill)) {
						System.err.println("Error - '" + fill + "' does not exist");
					}
				}

				for (int r = 2; r <= iterations.length; r++) {
					nCr = Array.nCr_indices(iterations.length, r);

					for (int[] element : nCr) {
						fill = iterations[element[0]];
						for (int k = 1; k < element.length; k++) {
							fill += "_" + iterations[element[k]];
						}
						System.out.print("\t" + ext.replaceAllWith(pattern, "[%%]", fill));
						command = "samtools merge";
						command += " /scratch/polymutt/" + ext.replaceAllWith(pattern, "[%%]", fill);
						for (int element2 : element) {
							command += " " + ext.replaceAllWith(pattern, "[%%]", iterations[element2]);
						}
						v.add(command);
					}
				}
				System.out.println();
			}
			reader.close();

			// Files.qsub("mergeBams", dir, -1, "[%0]", Matrix.toMatrix(Array.toStringArray(v)), 1500, 4);
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + controlFile + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + controlFile + "\"");
			System.exit(2);
		}
	}

	public static void assessCoverage(String filename, String bedfile, String coverages) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String temp;
		int count;
		SegmentLists segmentLists;
		Segment[][] segments;
		Logger log;
		byte chr, prevChr;
		int[] readCounts, thresholds;
		int[][] runningCounts;
		int onTargetReads;
		boolean allQualify;

		log = new Logger("assessingCoverageFor_" + ext.rootOf(filename) + ".out");

		if (bedfile == null) {
			segments = null;
		} else if (Files.exists(bedfile + ".ser")) {
			log.report("Reading in preserialized " + bedfile + ".ser");
			segments = SegmentLists.load(bedfile + ".ser", false).getLists();
		} else {
			log.report("Importing " + bedfile);
			segmentLists = SegmentLists.parseSegmentList(bedfile, 0, 1, 2, true);
			segmentLists.serialize(bedfile + ".ser");
			segments = segmentLists.getLists();
		}

		try {
			thresholds = Array.toIntArray(coverages.split(","));
		} catch (Exception e) {
			log.reportError("Error - failed to parse coverage thresholds: " + coverages);
			log.reportException(e);
			return;
		}
		onTargetReads = 0;
		runningCounts = new int[thresholds.length][4];

		prevChr = -1;
		temp = "not even started";
		try {
			if (Files.exists(ext.rootOf(filename, false) + "_onTargetCoverage.vcf.gz")) {
				log.report("Reading directly from '"	+ ext.rootOf(filename, false)
										+ "_onTargetCoverage.vcf.gz" + "'");
				reader = Files.getAppropriateReader(ext.rootOf(filename, false)
																						+ "_onTargetCoverage.vcf.gz");
				writer = null;
			} else {
				reader = Files.getAppropriateReader(filename);
				writer = Files.getAppropriateWriter(ext.rootOf(filename, false)
																						+ "_onTargetCoverage.vcf.gz");
				do {
					temp = reader.readLine();
				} while (!temp.startsWith("#CHROM"));

				if (!ext.checkHeader(temp.split("[\\s]+"), POLYMUTT_VCF_HEADER, false)) {
					log.reportError("Problem with header for file: " + filename);
					reader.close();
					writer.close();
					return;
				}
			}
			count = 0;
			while (reader.ready()) {
				if (count % 10000000 == 0) {
					log.report(ext.getTimestampForFilename() + "\t" + count);
				}
				try {
					temp = reader.readLine();
				} catch (EOFException eofe) {
					log.report("Truncated at line: " + count);
					reader.close();
					continue;
				}
				line = temp.trim().split("[\\s]+");
				if (line.length < 5) {
					log.report("Truncated at line: " + count);
					continue;
				}
				chr = Positions.chromosomeNumber(line[0], log);
				if (chr != prevChr) {
					log.report(ext.getTimestampForFilename()	+ "\tStarting chromosome "
											+ Positions.CHR_CODES[chr]);
					prevChr = chr;
				}
				if (writer == null) {
					readCounts = new int[] {Integer.parseInt(line[2]), Integer.parseInt(line[3]),
																	Integer.parseInt(line[4])};
				} else {
					if (segments == null
							|| (segments[chr] != null
									&& Segment.overlapsAny(	new Segment(chr, Integer.parseInt(line[1]),
																											Integer.parseInt(line[1])),
																					segments[chr]))) {
						readCounts = new int[] {Integer.parseInt(line[9].split(":")[2]),
																		Integer.parseInt(line[10].split(":")[2]),
																		Integer.parseInt(line[11].split(":")[2])};
						writer.println(line[0] + "\t" + line[1] + "\t" + Array.toStr(readCounts));
					} else {
						readCounts = null;
					}
				}
				if (readCounts != null) {
					for (int i = 0; i < thresholds.length; i++) {
						allQualify = true;
						for (int j = 0; j < 3; j++) {
							if (readCounts[j] > thresholds[i]) {
								runningCounts[i][j]++;
							} else {
								allQualify = false;
							}
						}
						if (allQualify) {
							runningCounts[i][3]++;
						}
					}
					onTargetReads++;
				}
				count++;
			}
			reader.close();
			if (writer != null) {
				writer.close();
			}
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
			System.exit(1);
		} catch (Exception e) {
			System.err.println("Error reading file \"" + filename + "\" at line:");
			System.err.println(temp);
			e.printStackTrace();
		}

		for (int i = 0; i < thresholds.length; i++) {
			log.report("# reads >= " + thresholds[i], false, true);
			for (int j = 0; j < 4; j++) {
				log.report("\t"	+ ext.formPercent((double) runningCounts[i][j] / (double) onTargetReads, 1),
										false, true);
			}
			log.report("");
		}
	}

	// bedfile can be null
	public static void assessAllPossibleCoverage(	String controlFile, String coverages,
																								String bedfile) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		int[] thresholds;
		Logger log;
		int[][] nCr;

		Vector<String> filenames;
		Vector<int[]> iterations;
		String pattern;
		String[] reps;
		String fill;
		String dir;
		String[][][] results;
		boolean problem;
		Segment[][] bedList;

		if (bedfile != null) {
			bedList = SegmentLists.parseSegmentList(bedfile, 0, 1, 2, false).getLists();
		} else {
			bedList = null;
		}

		try {
			reader = Files.getAppropriateReader(controlFile);
			dir = ext.verifyDirFormat(reader.readLine());
			log = new Logger(dir + "assessingAllPossibleCoverageFor_" + ext.rootOf(controlFile) + ".out");

			try {
				thresholds = Array.toIntArray(coverages.split(","));
			} catch (Exception e) {
				log.reportError("Error - failed to parse coverage thresholds: " + coverages);
				log.reportException(e);
				return;
			}

			writer = new PrintWriter(new FileWriter(dir + ext.rootOf(controlFile) + "_coverage.xln"));
			writer.println("Trio\tCombination\t# lanes\tReadDepth threshold\tFather\tMother\tChild\tEntire Trio");
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				pattern = line[0];
				reps = Array.subArray(line, 1);

				problem = false;
				filenames = new Vector<String>();
				for (String rep : reps) {
					fill = ext.replaceAllWith(pattern, "[%%]", rep);
					if (Files.exists(dir + fill)) {
						filenames.add(dir + fill);
					} else {
						System.err.println("Error - '" + fill + "' does not exist");
						problem = true;
					}
				}
				if (problem) {
					reader.close();
					writer.close();
					return;
				}

				iterations = new Vector<int[]>();
				for (int r = 1; r <= reps.length; r++) {
					nCr = Array.nCr_indices(reps.length, r);
					for (int[] element : nCr) {
						iterations.add(element);
					}
				}
				nCr = Matrix.toMatrix(iterations);

				results = parseAllPossibleCoverages(Array.toStringArray(filenames), nCr, thresholds,
																						bedList, log);

				for (int i = 0; i < nCr.length; i++) {
					fill = reps[nCr[i][0]];
					for (int k = 1; k < nCr[i].length; k++) {
						fill += "_" + reps[nCr[i][k]];
					}
					for (int j = 0; j < thresholds.length; j++) {
						writer.println(pattern	+ "\t" + fill + "\t" + nCr[i].length + "\t>=" + thresholds[j]
														+ "\t" + Array.toStr(results[i][j]));
					}
				}

			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + controlFile + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + controlFile + "\"");
			System.exit(2);
		}
	}

	private static String[][][] parseAllPossibleCoverages(String[] filenames, int[][] nCr,
																												int[] thresholds, Segment[][] bedList,
																												Logger log) {
		BufferedReader[] readers;
		String[][] currents;
		String temp;
		byte chr;
		byte[] chrs;
		int minPosition;
		int[] positions;

		int[] proxyCounts;
		int[][] readCounts;
		int[][][] runningCounts;
		int onTargetReads;
		boolean allQualify, done;

		String[][][] results;

		onTargetReads = 1;
		runningCounts = new int[nCr.length][thresholds.length][4];

		chr = -1;
		temp = "not even started";
		readers = new BufferedReader[filenames.length];
		currents = new String[filenames.length][];
		chrs = new byte[filenames.length];
		positions = new int[filenames.length];
		readCounts = new int[filenames.length][];
		proxyCounts = new int[4];
		try {
			for (int i = 0; i < filenames.length; i++) {
				readers[i] = Files.getAppropriateReader(filenames[i]);
				temp = readers[i].readLine();
				currents[i] = temp.split("[\\s]+");
				chrs[i] = Positions.chromosomeNumber(currents[i][0]);
				positions[i] = Integer.parseInt(currents[i][1]);
			}
			done = false;
			while (!done) {
				if (onTargetReads % 10000000 == 0) {
					log.report(ext.getTimestampForFilename() + "\t" + onTargetReads);
				}

				if (chr != Array.min(chrs)) {
					chr = Array.min(chrs);
					log.report(ext.getTimestampForFilename()	+ "\tStarting chromosome "
											+ Positions.CHR_CODES[chr]);
				}

				minPosition = Integer.MAX_VALUE;
				for (int i = 0; i < filenames.length; i++) {
					if (chrs[i] == chr && positions[i] < minPosition) {
						minPosition = positions[i];
					}
				}

				for (int i = 0; i < filenames.length; i++) {
					if (chrs[i] == chr && positions[i] == minPosition) {
						readCounts[i] = new int[] {	Integer.parseInt(currents[i][2]),
																				Integer.parseInt(currents[i][3]),
																				Integer.parseInt(currents[i][4])};
						try {
							temp = readers[i].readLine();
							currents[i] = temp.split("[\\s]+");
							if (currents[i].length < 5) {
								log.report("Not enough columns for line: " + Array.toStr(currents[i]));
								done = true;
							}
						} catch (Exception e) {
							log.report("File '"	+ filenames[i] + "' is truncated. Last valid position was: " + chr
													+ ":" + minPosition);
							log.reportException(e);
							done = true;
						}
						chrs[i] = Positions.chromosomeNumber(currents[i][0]);
						positions[i] = Integer.parseInt(currents[i][1]);
						if (!readers[i].ready()) {
							log.report("First file to end was '"	+ filenames[i] + "'. Last valid position was: "
													+ chr + ":" + minPosition);
							done = true;
						}
					} else {
						readCounts[i] = new int[] {0, 0, 0};
					}
				}

				Segment seg = new Segment(chr, minPosition, minPosition);
				if (bedList == null
						|| bedList[seg.getChr()] != null && Segment.overlapsAny(seg, bedList[seg.getChr()])) {
					for (int i = 0; i < nCr.length; i++) {
						for (int k = 0; k < 3; k++) {
							proxyCounts[k] = 0;
							for (int j = 0; j < nCr[i].length; j++) {
								proxyCounts[k] += readCounts[nCr[i][j]][k];
							}
						}
						for (int j = 0; j < thresholds.length; j++) {
							allQualify = true;
							for (int k = 0; k < 3; k++) {
								if (proxyCounts[k] > thresholds[j]) {
									runningCounts[i][j][k]++;
								} else {
									allQualify = false;
								}
							}
							if (allQualify) {
								runningCounts[i][j][3]++;
							}
						}
					}
					onTargetReads++;
				}
			}
			for (BufferedReader reader : readers) {
				reader.close();
			}
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file not found ");
			log.reportException(fnfe);
		} catch (Exception e) {
			System.err.println("Error reading file at line:");
			System.err.println(temp);
			e.printStackTrace();
		}

		results = new String[nCr.length][thresholds.length][4];
		for (int i = 0; i < nCr.length; i++) {
			for (int j = 0; j < thresholds.length; j++) {
				for (int k = 0; k < 4; k++) {
					results[i][j][k] = ext.formPercent((double) runningCounts[i][j][k]
																							/ (double) onTargetReads, 1);
				}
			}
		}
		return results;
	}

	public static void batchAssess() {
		Vector<String> v = new Vector<String>();
		int count;
		String[] files;
		String dir;

		dir = ext.pwd();

		count = 0;
		files = Files.list(dir, ".out.vcf.gz", false);
		new File(dir + "chunks/").mkdirs();
		for (int i = 0; i < files.length; i++) {
			if (!Files.exists(dir + ext.rootOf(files[i], false) + "_onTargetCoverage.vcf.gz")) {
				System.out.println(files[i]);
				count++;
				Files.qsub(dir	+ "chunks/runAssess" + count + ".qsub",
										"cd "																		+ dir + "\nmodule load java\njava -cp ~/"
																														+ org.genvisis.common.PSF.Java.GENVISIS
																														+ " seq.Polymutt assess=" + files[i],
										3000, 6, 1);
				v.add(dir + "chunks/runAssess" + count + ".qsub");
			}
		}
		Files.qsubMultiple(v, null, dir + "chunks/", "chunkAssess", 6, false, null, -1, 12000, 6);
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String triosFile = "trios_list.txt";
		boolean batchGLF = false;
		boolean batchPolymutt = false;
		boolean filterDenovo = false;
		String bamDir = "./";
		String vcfDir = Files.firstDirectoryThatExists(COMMON_VCF_LOCATIONS, true, false, new Logger());
		String annotationDir = null;
		String findDenovo = null;
		boolean findAll = false;
		String controlFile = null;
		int minReadDepth = 5;
		String assess = null;
		// String bedfile = "~/bin/ref/slimS04380219_Regions.bed";
		String bedfile = "/home/spectorl/pankrat0/bin/ref/slimS04380219_Regions.bed";
		String coverages = "8,10,12,16,20";
		boolean batchAssess = false;
		String allAssess = null;

		allAssess = "D:/Logan/Polymutt/coverage/list.txt";
		bedfile = "D:/Logan/Cosmic/cosmic.bed";
		assessAllPossibleCoverage(allAssess, coverages, bedfile);
		System.exit(1);

		String usage = "\n"	+ "seq.Polymutt requires 0-1 arguments\n"
										+ "   (1) list of trio root names (i.e. trios=" + triosFile + " (default))\n"
										+ "   (2) source directory containing the .bam files (i.e. bamDir=" + bamDir
										+ " (default))\n"
										+ "   (3) directory containing VCF files to create or parse (i.e. vcfDir="
										+ vcfDir + " (default))\n"
										+ "   (4) directory with SeattleSeq annotation (i.e. annotationDir=annotation/ (not the default))\n"
										+ " AND:\n"
										+ "   (5) batch BAM->GLF generation (requires bamDir) (i.e. -batchGLF (not the default))\n"
										+ " OR:\n"
										+ "   (5) batch Polymutt runs (requires vcfDir to output to) (i.e. -batchPolymutt (not the default))\n"
										+ " OR:\n"
										+ "   (5) batch BAM->GLF generation (requires bamDir) (i.e. -batchGLF (not the default))\n"
										+ " OR:\n"
										+ "   (1) find de novo events using Phred scores instead of DQ value (i.e. findDenovo=trio612.vcf (not the default))\n"
										+ "   (2) minimum read depth across all members of trio (i.e. minReadDepth="
										+ minReadDepth + " (default))\n" + " OR:\n"
										+ "   (1) find all via Phred scores in directory (i.e. -findAll (not the default))\n"
										+ "   (2) directory containing VCF files to create or parse (i.e. vcfDir="
										+ vcfDir + " (default))\n" + " OR:\n"
										+ "   (1) merge all possible bam files with a given pattern (i.e. controlFile=example.dat (not the default))\n"
										+ "         Example line: pattern_N00[%%] 5 6 7 8\n" + " OR:\n"
										+ "   (1) assess percent coverage for a particular vcf file (i.e. assess=somefile.vcf.gz (not the default))\n"
										+ "   (2) use a bedfile to filter out off target regions (i.e. bed=" + bedfile
										+ " (default))\n"
										+ "   (3) summarize percent coverage at various read depth coverages (i.e. coverage="
										+ coverages + " (default))\n" + " OR:\n"
										+ "   (1) batch assess coverage for all .vcf.gz files in folder (i.e. -batchAssess (not the default))\n"
										+ " OR:\n"
										+ "   (1) assess all possible combinations of coverage (i.e. allAssess=control_list.txt (not the default))\n"
										+ "   (2) various read depth coverages (i.e. coverage=" + coverages
										+ " (default))\n"
										+ "   (3) (optional) use monified bedfile to focus on specific areas of the genes/exone (i.e. bed="
										+ bedfile + " (default))\n" + "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("bamDir=")) {
				bamDir = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("trios=")) {
				triosFile = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("-batchGLF")) {
				batchGLF = true;
				numArgs--;
			} else if (arg.startsWith("-batchPolymutt")) {
				batchPolymutt = true;
				numArgs--;
			} else if (arg.startsWith("-filterDenovo")) {
				filterDenovo = true;
				numArgs--;
			} else if (arg.startsWith("-findAll")) {
				findAll = true;
				numArgs--;
			} else if (arg.startsWith("vcfDir=")) {
				vcfDir = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("annotationDir=")) {
				annotationDir = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("findDenovo=")) {
				findDenovo = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("minReadDepth=")) {
				minReadDepth = ext.parseIntArg(arg);
				numArgs--;
			} else if (arg.startsWith("controlFile=")) {
				controlFile = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("assess=")) {
				assess = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("coverage=")) {
				coverages = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("bed=")) {
				bedfile = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("-batchAssess")) {
				batchAssess = true;
				numArgs--;
			} else if (arg.startsWith("allAssess=")) {
				allAssess = arg.split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + arg);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}

		// filterDenovo = true;
		// vcfDir = "D:/Logan/Polymutt/results/";
		// triosFile = "D:/Logan/Polymutt/results/trios_list.txt";
		// annotationDir = "D:/Logan/Polymutt/results/annotation/";
		//
		// findDenovo = "D:/Logan/Polymutt/results/dedup_13_04236.vcf";

		// vcfDir = "D:/Logan/Polymutt/results/denovos/";
		// findAllLocalDenovo(vcfDir);
		// System.exit(1);

		// bedfile = "D:/Logan/Polymutt/coverage/slimS04380219_Regions.bed";
		// assessCoverage("D:/Logan/Polymutt/coverage/rrd_13_04236_L005.denovo.out.vcf.gz", bedfile,
		// coverages);
		// assessCoverage("D:/Logan/Polymutt/coverage/rrd_13_04236_L006.denovo.out.vcf.gz", bedfile,
		// coverages);
		// assessCoverage("D:/Logan/Polymutt/coverage/rrd_13_04236_L007.denovo.out.vcf.gz", bedfile,
		// coverages);
		// assessCoverage("D:/Logan/Polymutt/coverage/rrd_13_04236_L008.denovo.out.vcf.gz", bedfile,
		// coverages);
		// System.exit(1);

		// mergeAllPossible("example.dat");
		// System.exit(1);

		// assessAllPossibleCoverage("D:/Logan/Polymutt/coverage/list.txt", coverages);
		// System.exit(1);

		try {
			if (batchGLF) {
				batchAllGlf(bamDir);
			} else if (batchPolymutt) {
				System.err.println("Error - don't run until fixing the sex issue, default is currently to female");
				System.exit(1);
				batchPolymutt(triosFile);
			} else if (findAll) {
				if (Files.isWindows()) {
					findAllLocalDenovo(vcfDir);
				} else {
					findAllDenovo(vcfDir);
				}
			} else if (findDenovo != null) {
				findDenovo(findDenovo, minReadDepth);
			} else if (filterDenovo) {
				filterDenovo(triosFile, vcfDir, annotationDir);
			} else if (controlFile != null) {
				mergeAllPossible(controlFile);
			} else if (assess != null) {
				assessCoverage(assess, bedfile, coverages);
			} else if (batchAssess) {
				batchAssess();
			} else if (allAssess != null) {
				assessAllPossibleCoverage(allAssess, coverages, bedfile);
			} else {
				System.out.println("no subroutine selected");
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
