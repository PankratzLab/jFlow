package seq;


import java.io.*;
import java.util.*;
import common.*;

public class Polymutt {
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

	public static void filterDenovo(String triosFile) {
//		String[][] iterations;
//		String commands, scratch;
//		Vector<String> v;
//		
//		scratch = "/lustre/polymutt/";
//		scratch = "/project/scratch/polymutt/";
//		
//		iterations = HashVec.loadFileToStringMatrix(triosFile, false, new int[] {0,1,2}, false);
//		for (int i = 0; i < iterations.length; i++) {
//			if (Files.exists(scratch+iterations[i][0]+".vcf")) {
//				
//			} else {
//				System.err.println("Error - '"+scratch+iterations[i][0]+".vcf' does not exist");
//			}
//		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
//		String logfile = null;
//		Logger log;
		String sourceDir = null;
		String triosFile = "trios_list.txt";
		boolean filterDenovo = false;

		String usage = "\n" + 
			"seq.Polymutt requires 0-1 arguments\n" + 
			"   (1) source directory containing the .bam files (i.e. sourceDir=" + sourceDir + " (default))\n" + 
			"   (2) list of trio root names (i.e. trios=" + triosFile + " (default))\n" +
			"   (3) filter the VCF results (i.e. -filterDenovo (not the default))\n" +
			"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("sourceDir=")) {
				sourceDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("trios=")) {
				triosFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("-filterDenovo")) {
				filterDenovo = true;
				numArgs--;
//			} else if (args[i].startsWith("log=")) {
//				logfile = args[i].split("=")[1];
//				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
//			log = new Logger(logfile);
			if (sourceDir != null) {
				batchAllGlf(sourceDir);
			} else if (filterDenovo) {
				filterDenovo(triosFile);
			} else {
				System.err.println("Error - don't run until fixing the sex issue, default is currently to female");
				System.exit(1);
				batchPolymutt(triosFile);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
