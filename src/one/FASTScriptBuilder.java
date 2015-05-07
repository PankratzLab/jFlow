package one;

import java.io.File;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;

public class FASTScriptBuilder {
	
	static String FAST_LOC = "/home/pankarne/saonlib2/1000genomes/LnFXI/FAST.1.8.mc/FAST";
	static String dir = "/home/pankarne/chandap/ARIC.whites.impute2/";
	static String indivFile = "~/ordered9489.indiv";
	static String traitFile = "~/ordered9489.trait";
	static String filePattern = ".impute2.gz";
	static String runDir = "/home/pankarne/saonlib2/1000genomes/LnFXI/FAST_pC/";
	static int covarCount = 4;
	
	public static void run() throws IOException {
		String[] gzFiles = (new File(dir)).list(new FilenameFilter() {
			@Override
			public boolean accept(File dir, String name) {
				return name.endsWith(filePattern);
			}
		});
		
		PrintWriter scriptInputWriter = new PrintWriter(new FileWriter(runDir + "input.txt"));
		
		for (int i = 0; i < gzFiles.length; i++) {
			String chr = gzFiles[i].substring(3, 5);
			if (chr.charAt(1) == '.') {
				chr = "" + chr.charAt(0);
			}
		
			StringBuilder fastString = new StringBuilder(FAST_LOC)
				.append(" --mode genotype --impute2-geno-file ")
				.append(dir)
				.append(gzFiles[i])
				.append(" --impute2-info-file ")
				.append(dir)
				.append(gzFiles[i].substring(0, gzFiles[i].length() - 3))
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
				.append(gzFiles[i].substring(0, gzFiles[i].length() - 3))
				.append(".out");
			
			scriptInputWriter.println(fastString.toString());
			
		}
		
		scriptInputWriter.flush();
		scriptInputWriter.close();
			
	}
	
	public static void main(String[] args) {
		try {
			new FASTScriptBuilder().run();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
}