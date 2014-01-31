package seq;

import java.io.*;
import java.util.*;

import bioinformatics.SeattleSeq;
import common.*;

public class Polymutt {
	public static final String[] POLYMUTT_VCF_HEADER = {"#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "2", "3", "1"};
	public static final String[] COMMON_VCF_LOCATIONS = {"/lustre/polymutt/", "/project/scratch/polymutt/", "c:/scratch/polymutt/", "./"};
	
	public static void batchAllGlf(String sourceDir) {
		String[] files;
		String root, dir, filename;
		Vector<String> v = new Vector<String>();
		
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
		Vector<String> unknownAnnotations;
		String filename;

		log = new Logger("parseAllPolymuttInDir.log");
		
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
		parsedAnnotations = SeattleSeq.loadAllAnnotationInDir(annotationDir, log);
		iterations = HashVec.loadFileToStringMatrix(triosFile, false, new int[] {0,1,2}, false);
		for (int i = 0; i < iterations.length; i++) {
			if (Files.exists(vcfDir+iterations[i][0]+".vcf")) {
				parseDenovo(vcfDir+iterations[i][0]+".vcf", parsedAnnotations, unknownAnnotations, log);
			} else {
				System.err.println("Error - '"+vcfDir+iterations[i][0]+".vcf' does not exist");
			}
		}
		if (unknownAnnotations.size() > 0) {
			filename = "SeattleSeq_"+ext.getTimestampForFilename()+".input";
			Files.writeList(Array.toStringArray(unknownAnnotations), annotationDir+filename);
		}
	}

	private static void parseDenovo(String filename, Hashtable<String, String[]> parsedAnnotations, Vector<String> unknownAnnotations, Logger log) {
		BufferedReader reader;
//		PrintWriter writer;
		String[] line;
		String temp; // , trav;
//		Hashtable<String, String> hash = new Hashtable<String, String>();
//		Vector<String> v = new Vector<String>();
//		int count;
//		long time;
		String markerName, refAllele, denovoAllele;
		String genotype, parentalGenotypes, genotypes;
		boolean denovo;
		double dqScore, qualityScore;
		
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
				markerName = line[0]+":"+line[1];
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
							if (denovo) {
								denovoAllele += genotype.charAt(j);
							}
						}
					}
					genotypes += "\t"+genotype;
				}

				qualityScore = Double.parseDouble(line[5]);
				
				line = line[7].split("=");
				
				dqScore = Double.parseDouble(line[line.length-1]);
				
				if (denovoAllele.length() == 0) {
					denovoAllele = ".";
				} else {
					markerName += "_"+refAllele+"_"+denovoAllele;
				}
				
				
				if (!parsedAnnotations.containsKey(markerName)) {
					unknownAnnotations.add(markerName+"\t"+refAllele+"\t"+qualityScore+"\t"+dqScore+"\t"+denovoAllele+"\t"+genotypes);
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

	public static void main(String[] args) {
		int numArgs = args.length;
		String triosFile = "trios_list.txt";
		boolean batchGLF = false;
		boolean batchPolymutt = false;
		boolean filterDenovo = false;
		String bamDir = "./";
		String vcfDir = Files.firstDirectoryThatExists(COMMON_VCF_LOCATIONS, true, false);
		String annotationDir = null;
		
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
			} else if (args[i].startsWith("vcfDir=")) {
				vcfDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("annotationDir=")) {
				annotationDir = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		
		filterDenovo = true;
		vcfDir = "C:/scratch/LA/Polymutt/";
		triosFile = "C:/scratch/LA/Polymutt/trios.txt";
		
		try {
			if (batchGLF) {
				batchAllGlf(bamDir);
			} else if (batchPolymutt) {
				System.err.println("Error - don't run until fixing the sex issue, default is currently to female");
				System.exit(1);
				batchPolymutt(triosFile);
			} else if (filterDenovo) {
				filterDenovo(triosFile, vcfDir, annotationDir);
			} else {
				System.out.println("no subroutine selected");
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
