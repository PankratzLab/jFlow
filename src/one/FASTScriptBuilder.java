package one;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;

import common.Files;

public class FASTScriptBuilder {
	
	String FAST_LOC = "/home/pankarne/saonlib2/1000genomes/LnFXI/FAST.1.8.mc/FAST";
	String dir = "/home/pankarne/chandap/ARIC.whites.impute2/";
	String indivFile = "~/ordered9489.indiv";
	String traitFile = "~/ordered9489.trait";
	String filePattern = ".impute2.gz";
	String runDir = "/home/pankarne/saonlib2/1000genomes/LnFXI/FAST/FAST_pC/";
	int covarCount = 4;
	
	public FASTScriptBuilder(String FASTloc, String dataDir, String indivFile, String traitFile, String dataFileSuffix, String runDir, int covarCount) {
		this.FAST_LOC = FASTloc;
		this.dir = dataDir;
		this.indivFile = indivFile;
		this.traitFile = traitFile;
		this.filePattern = dataFileSuffix;
		this.runDir = runDir;
		this.covarCount = covarCount;
	}
	
	public void run() throws IOException {
		String[] dataFiles = (new File(dir)).list(new FilenameFilter() {
			@Override
			public boolean accept(File dir, String name) {
				return name.endsWith(filePattern);
			}
		});
		
		PrintWriter scriptInputWriter = new PrintWriter(new FileWriter(runDir + "input.txt"));
		
		for (int i = 0; i < dataFiles.length; i++) {
			String chr = dataFiles[i].substring(3, 5);
			if (chr.charAt(1) == '.') {
				chr = "" + chr.charAt(0);
			}
		
			StringBuilder fastString = new StringBuilder(FAST_LOC)
				.append(" --mode genotype --impute2-geno-file ")
				.append(dir)
				.append(dataFiles[i])
				.append(" --impute2-info-file ")
				.append(dir)
				.append(dataFiles[i].substring(0, dataFiles[i].length() - 3))
				.append("_info --indiv-file ")
				.append(indivFile)
				.append(" --trait-file ")
				.append(traitFile)
				.append(" --num-covariates ")
				.append(covarCount)
				.append(" --linear-snp ")
				.append(" --chr ")
				.append(chr)
				.append(" --out-file ")
				.append(runDir)
				.append("output/")
				.append(dataFiles[i].substring(0, dataFiles[i].length() - 3))
				.append(".out");
			
			scriptInputWriter.println(fastString.toString());
			
		}
		
		scriptInputWriter.flush();
		scriptInputWriter.close();
		
		int threads = 48;
		String command = "java -cp ~/park.jar one.ScriptExecutor file="+runDir+"input.txt token=null threads="+threads;
		
		Files.qsub(runDir + "master.qsub", command, 10000, 4, threads);
		
	}
	
	private static void concatResults(String resultsDir, String resultsFile) {
		String[] filenames = (new File(resultsDir)).list(new FilenameFilter() {
			@Override
			public boolean accept(File dir, String name) {
				return name.endsWith(".Linear.txt");
			}
		});
		
		PrintWriter writer = Files.getAppropriateWriter(resultsDir + resultsFile);
		boolean first = true;
		System.out.print("Concatenating results files: <");
		for (String str : filenames) {
			BufferedReader reader;
			try {
				reader = Files.getAppropriateReader(resultsDir + str);
				String line = reader.readLine();
				if (first) {
					writer.println(line);
					first = false;
				}
				while((line = reader.readLine()) != null) {
					writer.println(line);
				}
				reader.close();
				System.out.print("-");
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		System.out.println(">");
		writer.flush();
		writer.close();
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String fast = "~/FAST";
		String data = "~/data/";
		String indiv = "~/indiv.txt";
		String trait = "~/trait.txt";
		String suffix = ".impute2.gz";
		String run = "~/runFAST/";
		int covars = 0;
		String results = "~/FAST/output/";
		String out = "finalResults.txt";
		boolean concat = false;
		
		String usage = "one.FASTScriptBuilder requires 7 or 2 arguments\n" + 
					   "   (1) Full-path to FAST script (including /FAST) (i.e. fast=" + fast + " (default))\n" + 
					   "   (2) Full-path to data directory (i.e. data=" + data + " (default))\n" +
					   "   (3) Full-path to .indiv file (i.e. indiv=" + indiv + " (default))\n" +
					   "   (4) Full-path to .trait file (i.e. trait=" + trait + " (default))\n" +
					   "   (5) Suffix by which to identify data files in the data directory (i.e. suffix=" + suffix + " (default))\n" +
					   "   (6) Full-path to the directory in which you want to run these scripts (must include a folder named 'output') (i.e. rundir=" + run + " (default))\n" +
					   "   (7) Number of covariates in .trait file (i.e. covars=" + covars + " (default))\n" +
					   " OR " +
					   "   (1) Path to directory with results files (i.e. results=" + results + " (default))\n" +
					   "   (2) Desired name of concatenated result file (i.e. out=" + out + " (default))\n" +
					   "   (3) flag to indicate results processing is desired (i.e. -concat (not the default))\n" +
					   "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("fast=")) {
				fast = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("data=")) {
				data = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("indiv=")) {
				indiv = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("trait=")) {
				trait = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("suffix=")) {
				suffix = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("rundir=")) {
				run = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("results=")) {
				results = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("out=")) {
				out = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("covars=")) {
				covars = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("-concat")) {
				concat = true;
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
			if (concat) {
				concatResults(results, out);
			} else {
				new FASTScriptBuilder(fast, data, indiv, trait, suffix, run, covars).run();
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	
}