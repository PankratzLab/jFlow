package gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigInteger;
import java.util.Arrays;
import java.util.Comparator;

import parse.GenParser;
import common.Aliases;
import common.Array;
import common.Files;
import common.Logger;
import common.ext;

public class FAST {

	private static final String CHARGE_FORMAT = " 'SNP.id'=Markername 'Chr'=Chr 'Pos'=Pos $#<>=N 'Coded.Allele'=Effect_allele 'NonCoded.Allele'=Other_allele 'Coded.Af'=EAF 'Qual'=Imp_info 'Beta' 'Se'=SE 'Pval'=PValue";
	
	public static final String[] FORMATS = new String[]{CHARGE_FORMAT}; 
	
	public static void runParser(String FORMAT, String concattedResultsFile, String outFileName, int count) {
		String finalFormat = concattedResultsFile + " tab " + outFileName + FORMAT.replace("<>", count + "");
		String[] args = ext.removeQuotes(finalFormat).trim().split("[\\s]+");
		GenParser.parse(args, new Logger());
	}
	
	
	String FAST_LOC = "FAST";
	String dir = "/home/pankarne/chandap/ARIC.whites.impute2/";
	String indivFile = "~/ordered9489.indiv";
	String traitFile = "~/ordered9489.trait";
	String filePattern = ".impute2.gz";
	String runDir = "/home/pankarne/saonlib2/1000genomes/LnFXI/FAST/FAST_pC/";
	int covarCount = 4;
	
	public FAST(String FASTloc, String dataDir, String indivFile, String traitFile, String dataFileSuffix, String runDir, int covarCount) {
		this.FAST_LOC = FASTloc;
		this.dir = ext.verifyDirFormat(dataDir);
		this.indivFile = indivFile;
		this.traitFile = traitFile;
		this.filePattern = dataFileSuffix;
		this.runDir = ext.verifyDirFormat(runDir);
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
		
		int threads = 24;
		String command = "java -cp ~/park.jar one.ScriptExecutor file="+runDir+"input.txt token=took threads="+threads;
		
		Files.qsub(runDir + "master.qsub", command, 10000, 8, threads);
		
		(new File(runDir + "output/")).mkdirs();
	}
	
	private static void concatResults(String resultsDirectory, String resultsFile, double pvalThreshold, boolean writePValThresh, boolean runHitWindows) {
		String resultsDir = ext.verifyDirFormat(resultsDirectory);
		String[] filenames = (new File(resultsDir)).list(new FilenameFilter() {
			@Override
			public boolean accept(File dir, String name) {
				return name.endsWith(".Linear.txt");
			}
		});
		Arrays.sort(filenames, new Comparator<String>() {
			@Override
			public int compare(String o1, String o2) {
				String[] pts1 = o1.split("\\.");
				String[] pts2 = o2.split("\\.");
				
				Integer chr1 = Integer.valueOf(pts1[0].substring(3));
				Integer chr2 = Integer.valueOf(pts2[0].substring(3));
				
				int chrComp = chr1.compareTo(chr2);
				if (chrComp != 0) return chrComp;
				
				BigInteger pos1 = new BigInteger(pts1[1]);
				BigInteger pos2 = new BigInteger(pts2[1]);

				int posComp1 = pos1.compareTo(pos2);
				if (posComp1 != 0) return posComp1;
				
				BigInteger pos12 = new BigInteger(pts1[2]);
				BigInteger pos22 = new BigInteger(pts2[2]);
				
				return pos12.compareTo(pos22);
			}
		});
		
		PrintWriter writer = Files.getAppropriateWriter(resultsDir + resultsFile);
		int[] indices = null;
		
		PrintWriter pvalWriter = writePValThresh ? Files.getAppropriateWriter(resultsDir + "meetsPVal_" + ".out") : null;
		boolean first = true;
		System.out.print("Concatenating results files: <");
		for (String str : filenames) {
			BufferedReader reader;
			try {
				reader = Files.getAppropriateReader(resultsDir + str);
				String line = reader.readLine();
				if (first) {
					if (writePValThresh) {
						indices = ext.indexFactors(new String[][] {Aliases.PVALUES}, line.split("[\\s]+"), false, true, true, false);
					}
					writer.println(line);
					first = false;
				}
				while((line = reader.readLine()) != null) {
					if (writePValThresh && indices[0] != -1) {
						double pval = Double.parseDouble(line.split("[\\s]+")[indices[0]]);
						if (pval <= pvalThreshold) {
							pvalWriter.println(line);
						}
					}
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
		if (writePValThresh) {
			pvalWriter.flush();
			pvalWriter.close();
		}
		
		if (runHitWindows) {
			String[][] results = HitWindows.determine(resultsDir + resultsFile, 0.00000005f, 500000, 0.000005f, new String[0]);
			Files.writeMatrix(results, resultsDir + "hits.out", "\t");
		}
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
		
		int format = 0;
		boolean convert = false;
		int count = 0;
		
		double pval = 0.0001;
		boolean printPVals = false;
		
		boolean runHitWindows = false;
//		default p<5E-8, window p<5E-6
		//  0.00000005
		//  0.000005
		
		String usage = "gwas.FAST requires 7, 3, or 5 arguments\n" + 
					   "   (1) Full-path to FAST script (including /FAST) (i.e. fast=" + fast + " (default))\n" + 
					   "   (2) Full-path to data directory (i.e. data=" + data + " (default))\n" +
					   "   (3) Full-path to .indiv file (i.e. indiv=" + indiv + " (default))\n" +
					   "   (4) Full-path to .trait file (i.e. trait=" + trait + " (default))\n" +
					   "   (5) Suffix by which to identify data files in the data directory (i.e. suffix=" + suffix + " (default))\n" +
					   "   (6) Full-path to the directory in which you want to run these scripts (must include a folder named 'output') (i.e. rundir=" + run + " (default))\n" +
					   "   (7) Number of covariates in .trait file (i.e. covars=" + covars + " (default))\n" +
					   " OR \n" +
					   "   (1) Flag to indicate results processing is desired (i.e. -concat (not the default))\n" +
					   "   (2) Path to directory with results files (i.e. results=" + results + " (default))\n" +
					   "   (3) Desired name of concatenated result file (i.e. out=" + out + " (default))\n" +
					   " OR \n" +
					   "   (1) Flag to indicate format conversion processing is desired (i.e. -convert (not the default))\n" +
					   "   (2) Path to concatenated result files (i.e. results=" + results + " (default))\n" +
					   "   (3) Desired name of processed result file (i.e. out=" + out + " (default))\n" +
					   "   (4) Format flag: (i.e. format=" + format + " (default))\n" + 
					   "              0: CHARGE format \n" +
					   "   (5) Number of individuals in analysis (i.e. count=" + covars + " (not the default))\n" +
					   " OR \n" +
					   "   -concat and -convert can be combined:\n" +
					   "   (1) Both -concat and -convert flags\n" +
					   "   (2) Path to directory with results files (i.e. results=" + results + " (default))\n" +
					   "   (3) Desired name of processed result file (i.e. out=" + out + " (default))\n" +
					   "   (4) Format flag: (i.e. format=" + format + " (default))\n" + 
					   "           FORMATS:\n" + 
					   "               0: CHARGE format \n" +
					   "   (5) Number of individuals in analysis (i.e. count=" + count + " (not the default))\n" +
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
			} else if (args[i].startsWith("format=")) {
				format = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("count=")) {
				count = Integer.parseInt(args[i].split("=")[1]);
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
			} else if (args[i].startsWith("-convert")) {
				convert = true;
				numArgs--;
			} else if (args[i].startsWith("-writePVals")) {
				printPVals = true;
				numArgs--;
			} else if (args[i].startsWith("-hitWindows")) {
				runHitWindows = true;
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
			if (concat && convert) {
				String midOut = "concatenated.result";
				concatResults(results, midOut, pval, printPVals, runHitWindows);
				runParser(FORMATS[format], midOut, out, count);
			} else if (concat) {
				concatResults(results, out, pval, printPVals, runHitWindows);
			} else if (convert) {
				runParser(FORMATS[format], results, out, count);
			} else {
				new FAST(fast, data, indiv, trait, suffix, run, covars).run();
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
}
	