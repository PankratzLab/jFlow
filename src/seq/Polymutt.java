package seq;

import java.io.*;
import java.util.*;

import bioinformatics.SeattleSeq;
import common.*;

public class Polymutt {
	public static final String[] POLYMUTT_VCF_HEADER = {"#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "2", "3", "1"};
	public static final String[] COMMON_VCF_LOCATIONS = {"/lustre/polymutt/", "/project/scratch/polymutt/", "c:/scratch/polymutt/", "D:/Logan/Polymutt/results/", "./"};
	public static final String[] GENOTYPE_ORDER = {"A/A", "A/C", "A/G", "A/T", "C/C", "C/G", "C/T", "G/G", "G/T", "T/T"};
	
	public static void batchAllGlf(String sourceDir) {
		String[] files;
		String root, dir, filename;
		Vector<String> v = new Vector<String>();
		
		if (!sourceDir.startsWith("/")) {
			System.err.println("Error - sourceDirectory needs to explicit");
			return;
		}
		
		dir = ext.pwd();
		new File(dir+"batches/").mkdirs();
		sourceDir = ext.verifyDirFormat(sourceDir);
		v = new Vector<String>();
		files = Files.list(sourceDir, ".bam", false);
		for (int i = 0; i < files.length; i++) {
			root = ext.rootOf(files[i]);
			filename = root+"_batch";
			Files.write("~/bin/samtools-hybrid pileup "+sourceDir+root+".bam -g -f ~/bin/ref/hg19_canonical.fa > "+dir+root+".bam.glf", dir+"batches/"+filename);
			Files.chmod(dir+"batches/"+filename);
			v.add(dir+"batches/"+filename);
		}
//		Files.qsubMultiple(v, null, dir+"batches/", "glfAll", 16, true, "sb", -1, 62000, 2);		
		Files.qsubMultiple(v, null, dir+"batches/", "glfAll", 8, true, null, -1, 62000, 2);		
	}

	private static void batchPolymutt(String triosFile) {
		String[][] iterations;
		String commands, scratch;
		Vector<String> v;
		
		scratch = "/lustre/polymutt/";
		scratch = "/project/scratch/polymutt/";
		
		v = new Vector<String>();
		Files.write("T\tGLF_Index", "same_for_all.dat");
		Files.writeList(new String[] {
				"fam	3	0	0	1	3",
				"fam	2	0	0	2	2",
				"fam	1	3	2	2	1"
		}, "same_for_all.ped");
		iterations = HashVec.loadFileToStringMatrix(triosFile, false, new int[] {0,1,2}, false);
		for (int i = 0; i < iterations.length; i++) {
			Files.writeList(new String[] {
					"1 "+iterations[i][0]+".bam.glf",
					"2 "+iterations[i][1]+".bam.glf",
					"3 "+iterations[i][2]+".bam.glf"
			}, iterations[i][0]+".gif");

			commands = "cd "+ext.pwd()+"\n"
					+"~/bin/polymutt -p same_for_all.ped -d same_for_all.dat -g "+iterations[i][0]+".gif --nthreads 8  --out_vcf "+scratch+iterations[i][0]+".denovo.out.vcf --denovo --rate_denovo 1.5e-07\n"
					+"cd "+scratch+"\n"
					+"more "+iterations[i][0]+".denovo.out.vcf | grep -v 'DQ=-0' | grep -v 'DQ=0.00' > "+iterations[i][0]+".vcf\n"
					+"gzip "+iterations[i][0]+".denovo.out.vcf";

			Files.qsub("batches/run_"+iterations[i][0], commands, 22000, 48, 8);
			Files.chmod("batches/run_"+iterations[i][0]);
			v.add("qsub batches/run_"+iterations[i][0]);
		}
		v.insertElementAt("cd "+ext.pwd(), 0);
		Files.writeList(Array.toStringArray(v), "master.allRuns");
		Files.chmod("master.allRuns");
	}

	public static void filterDenovo(String triosFile, String vcfDir, String annotationDir) {
		String[][] iterations;
		Logger log;
		Hashtable<String, String[]> parsedAnnotations;
		Vector<String> unknownAnnotations, finishedAnnotations;
		String filename, temp;

		log = new Logger(ext.parseDirectoryOfFile(triosFile)+"parseAllPolymuttInDir.log");
		
		if (!Files.exists(vcfDir)) {
			System.err.println("Error - vcf directory not found: "+vcfDir);
			return;
		}

		if (annotationDir == null) {
			log.reportError("annotationDir directory set to the cwd");
			annotationDir = "./";
		} else if (!Files.exists(annotationDir)) {
			log.reportError("Warning - annotationDir directory not found: "+annotationDir);
			log.reportError("          annotationDir directory reset to the cwd");
		}
		
		annotationDir = ext.verifyDirFormat(annotationDir);
		
		if (Files.isDirectory(triosFile)) {
			log.reportError("Error - defined trios \"file\" is actually a directory");
			return;
		} else if (!Files.exists(triosFile)) {
			log.reportError("Error - trios file not found: "+triosFile);
			return;
		}
		
		unknownAnnotations = new Vector<String>();
		finishedAnnotations = new Vector<String>();
		parsedAnnotations = SeattleSeq.loadAllAnnotationInDir(annotationDir, log);
		iterations = HashVec.loadFileToStringMatrix(triosFile, false, new int[] {0,1,2}, false);
		for (int i = 0; i < iterations.length; i++) {
			if (Files.exists(vcfDir+iterations[i][0]+".vcf")) {
				parseDenovo(vcfDir+iterations[i][0]+".vcf", parsedAnnotations, unknownAnnotations, finishedAnnotations, log);
			} else {
				System.err.println("Error - '"+vcfDir+iterations[i][0]+".vcf' does not exist");
			}
		}
		if (unknownAnnotations.size() > 0) {
			filename = "SeattleSeq_"+ext.getTimestampForFilename()+".input";
			Files.writeList(Array.toStringArray(unknownAnnotations), annotationDir+filename);
		}
		temp = "Filename\tSampleRoot\tMarkerName\tChr\tPosition\tREF\tALT\tQUAL\tMapQual\tDenovoQual";
		for (int i = 1; i <= 3; i++) {
			temp += "\tGeno"+i;
		}
		temp += "\tFunction\tinDbSNP";
		for (int i = 1; i <= 3; i++) {
			temp += "\tReadDepth"+i;
		}
		temp += "\tReadDepthTotal";
		for (int i = 1; i <= 3; i++) {
			temp += "\tGenQual"+i;
		}
		for (int i = 1; i <= 3; i++) {
			temp += "\tPhred"+i;
		}
		for (int i = 1; i <= 3; i++) {
			temp += "\tNextBestGeno"+i;
		}
		
		finishedAnnotations.insertElementAt(temp, 0);
		Files.writeList(Array.toStringArray(finishedAnnotations), ext.parseDirectoryOfFile(triosFile)+"summary.xln");
	}

	private static void parseDenovo(String filename, Hashtable<String, String[]> parsedAnnotations, Vector<String> unknownAnnotations, Vector<String> finishedAnnotations, Logger log) {
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
				log.reportError("Problem with header for file: "+filename);
				reader.close();
				return;
			}
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				chr = line[0];
				position = line[1];
				markerName = chr+":"+position;
				refAllele = line[3];
				genotypes = "";
				parentalGenotypes = "";
				denovoAllele = "";
				for (int i = 0; i < 3; i++) {
					genotype = ext.replaceAllWith(line[9+i].split(":")[0], "/","");
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
									log.reportError("Warning - double denovo for "+markerName+" in "+ext.removeDirectoryInfo(filename));
									denovo = false;
								}
							}
							if (denovo) {
								denovoAllele += genotype.charAt(j);
							}
						}
					}
					genotypes += "\t"+genotype;
				}

//				double[] genotypeQualities, phredOfCall, nextBestPhred;
//				String[] phreds;
				
				qualityScore = Double.parseDouble(line[5]);

				
				subline = line[7].split(";");

				readDepths = new int[4];
				readDepths[3] = Integer.parseInt(subline[2].split("=")[1]);
				mapQuality = Double.parseDouble(subline[3].split("=")[1]);
				dqScore = Double.parseDouble(subline[4].split("=")[1]);
				
				if (!line[8].equals("GT:GQ:DP:PL")) {
					log.reportError("Error - column 9 is "+line[8]+" instead of GT:GQ:DP:PL");
				}
				
				phredOfCall = new int[3];
				nextBestGenotype = new String[3];
				genotypeQualities = new double[3];
				for (int i = 0; i < 3; i++) {
					subline = line[9+i].split(":");
					phreds = Array.toIntArray(subline[3].split(","));
					indexOfPhred = ext.indexOfStr(subline[0], GENOTYPE_ORDER);
					phredOfCall[i] = phreds[indexOfPhred];
					phreds[indexOfPhred] = 999;
					nextBestGenotype[i] = GENOTYPE_ORDER[Sort.quicksort(phreds)[0]];
					genotypeQualities[i] = Double.parseDouble(subline[1]);
					readDepths[i] = Integer.parseInt(subline[2]);
				}
				
				if (denovoAllele.length() > 1) {
					log.reportError("Warning - double heterozygote denovo for "+markerName+" in "+ext.removeDirectoryInfo(filename));
				}
				
				for (int i = 0; i < denovoAllele.length(); i++) {
					trav = markerName+"_"+refAllele+"_"+denovoAllele.charAt(i);
					if (!parsedAnnotations.containsKey(trav)) {
						unknownAnnotations.add(trav+"\t"+chr+"\t"+position+"\t"+refAllele+"\t"+denovoAllele.charAt(i));
					} else {
						if (parsedAnnotations.get(trav).length > 0) {
							finishedAnnotations.add(ext.rootOf(filename)+"\t"+decrypt(ext.rootOf(filename))+"\t"+trav+"\t"+chr+"\t"+position+"\t"+refAllele+"\t"+denovoAllele.charAt(i)+"\t"+qualityScore+"\t"+mapQuality+"\t"+dqScore+""+genotypes+"\t"+parsedAnnotations.get(trav)[0]+"\t"+Array.toStr(readDepths)+"\t"+Array.toStr(genotypeQualities)+"\t"+Array.toStr(phredOfCall)+"\t"+Array.toStr(nextBestGenotype));
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
		for (int i = 0; i<stringScores.length; i++) {
			try {
				scores[i] = Integer.parseInt(stringScores[i]);
			} catch (NumberFormatException nfe) {
				System.err.println("Error - failed to convert '"+stringScores[i]+"' into an integer. Line: "+phredScores);
			}
		}
		
		if (scores.length != GENOTYPE_ORDER.length) {
			log.reportError("Error - mismatched number of phred scores. Found "+scores.length+" when expecting "+GENOTYPE_ORDER.length+". scores: "+phredScores);
		}
		
		mostLikelyGenotypes = new Vector<String>();
		order = Sort.quicksort(scores);
		bestScore = scores[order[0]];
		count = 0;
		while (count < GENOTYPE_ORDER.length && scores[order[count]] == bestScore) {
			mostLikelyGenotypes.add(GENOTYPE_ORDER[order[count]]);
			count++;
		}
		if (count == GENOTYPE_ORDER.length) {
			if (Array.max(scores) > 0) {
				log.report("There was no likely genotype. Line: "+phredScores);
			}
		}
		
		return Array.toStringArray(mostLikelyGenotypes);
	}

	private static void findAllDenovo(String dir) {
		String[] filenames;
		
		filenames = Files.list(dir, ".vcf.gz", false);
		
		Files.qsub("findAll", ext.pwd(), -1, "java -cp ~/park.jar seq.Polymutt findDenovo=[%0]", Matrix.toMatrix(filenames), 1000, 12);
	}
	
	private static void findAllLocalDenovo(String dir) {
		String[] filenames;
		
		filenames = Files.list(dir, "vcf_denovo.vcf", false);
		
		for (int i = 0; i < filenames.length; i++) {
			System.out.println(filenames[i]);
			findDenovo(dir+filenames[i], 8);
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
		
		log = new Logger(ext.rootOf(filename, false)+".log");
		
		try {
			reader = Files.getAppropriateReader(filename);
			new File(ext.parseDirectoryOfFile(filename)+"denovos/").mkdirs();
			outRoot = ext.parseDirectoryOfFile(filename)+"denovos/"+ext.rootOf(filename, true);
			writers = new PrintWriter[3];
			writers[0] = new PrintWriter(new FileWriter(outRoot+"_denovo.vcf"));
			writers[1] = new PrintWriter(new FileWriter(outRoot+"_ambiguousGenotypesInChild.vcf"));
//			writers[2] = new PrintWriter(new FileWriter(outRoot+"_lowCoverage.vcf"));
			do {
				temp = reader.readLine();
				writers[0].println(temp);
			} while (reader.ready() && !temp.startsWith("#CHROM"));
			
			if (!ext.checkHeader(temp.split("[\\s]+"), POLYMUTT_VCF_HEADER, false)) {
				log.reportError("Problem with header for file: "+filename);
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
					log.report(ext.getTime()+"\t"+chr);
					prevChr = chr;
				}
				
				position = line[1];
				markerName = chr+":"+position;
				if (line.length < 11) {
					log.reportError("VCF file '"+filename+"' truncated at marker "+markerName);
					reader.close();
					Files.closeAll(writers);
					return;
				}
				parentalGenotypes = "";
				denovoAllele = "";
				readDepths = new int[3];
				for (int i = 0; i < 3; i++) {
					if (line[9+i].split(":").length < 4 || line[9+i].split(":")[3].split(",").length < 10) {
						log.reportError("VCF file '"+filename+"' truncated at marker "+markerName);
						reader.close();
						Files.closeAll(writers);
						return;
					}
					mostLikelyGenotypes = getMostLikelyGenotypesFromPhredScores(line[9+i].split(":")[3], log);
					readDepths[i] = Integer.parseInt(line[9+i].split(":")[2]);
					if (mostLikelyGenotypes.length == GENOTYPE_ORDER.length) {
//						writers[2].println(temp);
					}
					if (i == 2 && mostLikelyGenotypes.length > 1) {
						if (mostLikelyGenotypes.length == GENOTYPE_ORDER.length) {
							mostLikelyGenotypes = new String[] {line[9+i].split(":")[0]};
//							writers[2].println(temp);
						} else {
//							log.reportError("Warning - there are more than one 'most likely genotype' for child in "+ext.rootOf(filename)+" for marker "+markerName+" so it will not be called denovo. Line: "+line[9+i].split(":")[3]);
							writers[1].println(temp);
						}
						// TODO check if this is possible then dump and consider
					}
					for (int g = 0; g < mostLikelyGenotypes.length; g++) {
						genotype = ext.replaceAllWith(mostLikelyGenotypes[g], "/","");
						
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
										log.reportError("Warning - double denovo for "+markerName+" in "+ext.removeDirectoryInfo(filename));
										denovo = false;
									}
								}
								//  denovo  and         reasonably high mapping quality        and         high genotype quality for kid                       and  at least 5 reads for all three individuals 
								if (denovo && Double.parseDouble(line[7].split(";")[3].split("=")[1]) > 50 && Double.parseDouble(line[9+2].split(":")[1]) > 50 && Array.min(readDepths) >= minReadDepth) {
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

				for (int i = 0; i < iterations.length; i++) {
					fill = ext.replaceAllWith(pattern, "[%%]", iterations[i]);
					if (!Files.exists(dir+fill)) {
						System.err.println("Error - '"+fill+"' does not exist");
					}
				}
				
				for (int r = 2; r <= iterations.length; r++) {
					nCr = Array.nCr_indices(iterations.length, r);
					
					for (int i = 0; i < nCr.length; i++) {
						fill = iterations[nCr[i][0]];
						for (int k = 1; k < nCr[i].length; k++) {
							fill += "_"+iterations[nCr[i][k]];
						}
						System.out.print("\t"+ext.replaceAllWith(pattern, "[%%]", fill));
						command = "samtools merge";
						command += " /scratch/polymutt/"+ext.replaceAllWith(pattern, "[%%]", fill);
						for (int k = 0; k < nCr[i].length; k++) {
							command += " "+ext.replaceAllWith(pattern, "[%%]", iterations[nCr[i][k]]);
						}
						v.add(command);
					}
				}
				System.out.println();
			}
			reader.close();

//			Files.qsub("mergeBams", dir, -1, "[%0]", Matrix.toMatrix(Array.toStringArray(v)), 1500, 4);
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + controlFile + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + controlFile + "\"");
			System.exit(2);
		}
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String triosFile = "trios_list.txt";
		boolean batchGLF = false;
		boolean batchPolymutt = false;
		boolean filterDenovo = false;
		String bamDir = "./";
		String vcfDir = Files.firstDirectoryThatExists(COMMON_VCF_LOCATIONS, true, false);
		String annotationDir = null;
		String findDenovo = null;
		boolean findAll = false;
		String controlFile = null;
		int minReadDepth = 5;
		
		String usage = "\n" + 
			"seq.Polymutt requires 0-1 arguments\n" + 
			"   (1) list of trio root names (i.e. trios=" + triosFile + " (default))\n" +
			"   (2) source directory containing the .bam files (i.e. bamDir=" + bamDir + " (default))\n" + 
			"   (3) directory containing VCF files to create or parse (i.e. vcfDir="+vcfDir+" (default))\n" +
			"   (4) directory with SeattleSeq annotation (i.e. annotationDir=annotation/ (not the default))\n" +
			" AND:\n" +
			"   (5) batch BAM->GLF generation (requires bamDir) (i.e. -batchGLF (not the default))\n" +
			" OR:\n" +
			"   (5) batch Polymutt runs (requires vcfDir to output to) (i.e. -batchPolymutt (not the default))\n" +
			" OR:\n" +
			"   (5) batch BAM->GLF generation (requires bamDir) (i.e. -batchGLF (not the default))\n" +
			" OR:\n" +
			"   (1) find de novo events using Phred scores instead of DQ value (i.e. findDenovo=trio612.vcf (not the default))\n" +
			"   (2) minimum read depth across all members of trio (i.e. minReadDepth="+minReadDepth+" (default))\n" +
			" OR:\n" +
			"   (1) find all via Phred scores in directory (i.e. -findAll (not the default))\n" +
			"   (2) directory containing VCF files to create or parse (i.e. vcfDir="+vcfDir+" (default))\n" +
			" OR:\n" +
			"   (1) merge all possible bam files with a given pattern (i.e. controlFile=example.dat (not the default))\n" +
			"         Example line: pattern_N00[%%] 5 6 7 8\n" +
			"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("bamDir=")) {
				bamDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("trios=")) {
				triosFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("-batchGLF")) {
				batchGLF = true;
				numArgs--;
			} else if (args[i].startsWith("-batchPolymutt")) {
				batchPolymutt = true;
				numArgs--;
			} else if (args[i].startsWith("-filterDenovo")) {
				filterDenovo = true;
				numArgs--;
			} else if (args[i].startsWith("-findAll")) {
				findAll = true;
				numArgs--;
			} else if (args[i].startsWith("vcfDir=")) {
				vcfDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("annotationDir=")) {
				annotationDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("findDenovo=")) {
				findDenovo = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("minReadDepth=")) {
				minReadDepth = ext.parseIntArg(args[i]);
				numArgs--;				
			} else if (args[i].startsWith("controlFile=")) {
				controlFile = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		
//		filterDenovo = true;
//		vcfDir = "D:/Logan/Polymutt/results/";
//		triosFile = "D:/Logan/Polymutt/results/trios_list.txt";
//		annotationDir = "D:/Logan/Polymutt/results/annotation/";
//		
//		findDenovo = "D:/Logan/Polymutt/results/dedup_13_04236.vcf";

		vcfDir = "D:/Logan/Polymutt/results/denovos/";
		findAllLocalDenovo(vcfDir);
		System.exit(1);
		
//		mergeAllPossible("example.dat");
//		System.exit(1);
		
		try {
			if (batchGLF) {
				batchAllGlf(bamDir);
			} else if (batchPolymutt) {
				System.err.println("Error - don't run until fixing the sex issue, default is currently to female");
//				System.exit(1);
				batchPolymutt(triosFile);
			} else if (findAll) {
				findAllDenovo(vcfDir);
			} else if (findDenovo != null) {
				findDenovo(findDenovo, minReadDepth);
			} else if (filterDenovo) {
				filterDenovo(triosFile, vcfDir, annotationDir);
			} else if (controlFile != null) {
				mergeAllPossible(controlFile);
			} else {
				System.out.println("no subroutine selected");
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
