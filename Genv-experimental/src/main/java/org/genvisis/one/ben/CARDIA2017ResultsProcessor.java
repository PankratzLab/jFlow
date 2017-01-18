package org.genvisis.one.ben;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.ext;
import org.genvisis.stats.ProbDist;

public class CARDIA2017ResultsProcessor {
	
	String[] IN_HDR = {
	                   "name", "A1", "A2", "Freq1", "MAF", "Quality", 
	                   "Rsq", // 6
	                   "n", 
	                   "Mean_predictor_allele", // 8 ----- EffAF??
	                   "beta_SNP_add", 
	                   "sebeta_SNP_add", // 10
	                   "beta_SNP_DPA", 
	                   "sebeta_SNP_DPA", // 12
	                   "cov_SNP_int_SNP_DPA", 
	                   "chi2_SNP" // 14
	};
	String[] OUT_HDR = {
	                    "chr", 
	                    "pos", 
	                    "marker_name", 
	                    "strand", 
	                    "base_allele", 
	                    "effect_allele", 
	                    "N", // IN_7
	                    "effect_allele_freq", 
	                    "imputation_type", 
	                    "Imputation_value", // IN_6
	                    "beta_main", // IN_9
	                    "se_main", // IN_10
	                    "beta_int", // IN_11
	                    "se_int", // IN_12
	                    "cov", // IN_13
	                    "chi_2df", // IN_14 
	                    "P_2df", // ProbDist.NormDist(IN_11/IN_12)
	                    "chi_P_2df" // ProbDist.chiDist(IN_14, 1)
  };
	String STRAND = "+";
	String TYPE = "1";
	String dir = "/scratch.global/cole0482/CARDIA_DATA/";
	static final String EA = "EA/";
	static final String AA = "AA/";
	String LOOKUP = "vcfLookup/";
	String[] MODEL_VARS = {"FEV1", "FVC", "FEV1", "FVC"};
	String[] NUTR_VARS = {"DHA", "DHA", "DPA", "DPA"};
	
	static final class SNP {
		String name;
		int chr;
		int pos;
		String base;
		String eff;
		
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + ((name == null) ? 0 : name.hashCode());
			return result;
		}
		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			SNP other = (SNP) obj;
			if (name == null) {
				if (other.name != null)
					return false;
			} else if (!name.equals(other.name))
				return false;
			return true;
		}
		
	}
	
	private String getAnc(String prefix) {
		if (AA.equals(prefix)) {
			return "AA";
		} else {
			return "Cauc";
		}
	}
	
	private void process(String prefix) throws IOException {
		BufferedReader reader;
		PrintWriter writer;
		PrintWriter missWriter;
		String dir = "/scratch.global/cole0482/CARDIA_2017/";
		String modelPref = "model_";
		int numModels = 4;
		
		HashMap<Integer, PrintWriter> writers = new HashMap<>();
		HashMap<Integer, PrintWriter> missWriters = new HashMap<>();
		
		for (int chr = 1; chr < 23; chr++) {
			HashMap<String, SNP> info = processChrFile(chr, prefix);
			for (int i = 0; i < numModels; i++) {
				String modelDir = dir + prefix + modelPref + (i + 1) + "/";
				String regFile = modelDir + "regression_chr" + chr + ".out_add.out.txt";
				String outFile = modelDir + "out/" + "CARDIA_" + MODEL_VARS[i] + "_" + getAnc(prefix) + "_" + NUTR_VARS[i] + "_" + (i + 1) + "_DDMMYYYY.txt";
				
				reader = Files.getAppropriateReader(regFile);
				writer = writers.get(i);
				if (writer == null) {
					writer = Files.getAppropriateWriter(outFile);
					writers.put(i, writer);
				}
				missWriter = missWriters.get(i);
				if (missWriter == null) {
					missWriter = Files.getAppropriateWriter(ext.rootOf(outFile, false) + "_missing.txt");
					missWriters.put(i, missWriter);
				}
				writer.println(Array.toStr(OUT_HDR));
				reader.readLine();
				String line;
				String[] parts;
				while ((line = reader.readLine()) != null) {
					parts = line.split("[\\s]+");
					StringBuilder sb = new StringBuilder();
					SNP snp = info.get(parts[0]);
					if (snp == null) {
						System.err.println("Error - missing SNP INFO for SNP: " + parts[0] + " in chr " + chr + " results file for EA.");
						continue;
					}
					if (snp.pos == -1) {
						missWriter.println(parts[0]);
					}
					sb.append(snp.chr).append("\t")
						.append(snp.pos).append("\t")
						.append(parts[0]).append("\t")
						.append(STRAND).append("\t")
						.append(snp.base).append("\t")
						.append(snp.eff).append("\t")
						.append(parts[7]).append("\t")
						.append(parts[8]).append("\t")
						.append(TYPE).append("\t")
						.append(parts[6]).append("\t")
						.append(parts[9]).append("\t")
						.append(parts[10]).append("\t")
						.append(parts[11]).append("\t")
						.append(parts[12]).append("\t")
						.append(parts[13]).append("\t")
						.append(parts[14]).append("\t")
						.append(ext.isMissingValue(parts[11]) || ext.isMissingValue(parts[12]) ? "NaN" : ProbDist.NormDist(Double.parseDouble(parts[11]) / Double.parseDouble(parts[12]))).append("\t")
						.append(ext.isMissingValue(parts[14]) ? "NaN" : ProbDist.ChiDist(Double.parseDouble(parts[14]), 1));
					writer.println(sb.toString());
				}
				reader.close();
			}
		}
		
		for (PrintWriter wr : writers.values()) {
			wr.flush();
			wr.close();
		}
		
	}
	
	private HashMap<String, SNP> processChrFile(int chr, String prefix) {
		String infoFile = dir + prefix + "chr" + chr + ".info"; 
		String posFile = dir + prefix + LOOKUP + "chr" + chr + "_positions.xln";
		
		HashMap<String, SNP> map = new HashMap<>();
		
		String[][] matr = HashVec.loadFileToStringMatrix(infoFile, true, new int[]{0, 1, 2, 3, 4}, false);
		for (String[] line : matr) {
			SNP snp = new SNP();
			snp.name = line[0];
			snp.base = line[1]; // assuming A1 is reference
			snp.eff = line[2]; // assuming A2 is predictor
			map.put(snp.name, snp);
		}
		
		matr = HashVec.loadFileToStringMatrix(posFile, true, new int[]{0, 1, 2}, false);
		for (String[] line : matr) {
			SNP snp = map.get(line[0]);
			if (snp == null) {
				System.err.println("Error - snp found in position file that wasn't found in the info file: " + line[0]);
				snp = new SNP();
				snp.name = line[0];
				snp.base = ".";
				snp.eff = ".";
			}
			if (".".equals(line[1])) {
				if (snp.name.contains(":")) {
					String[] pts = snp.name.split(":");
					try {
						snp.chr = Integer.parseInt(pts[0]);
					} catch (NumberFormatException e) {
						snp.chr = chr;
					}
					snp.pos = Integer.parseInt(pts[1]);
				} else {
					snp.chr = chr;
					snp.pos = -1;
				}
			} else {
				try {
					snp.chr = Integer.parseInt(line[1]);
				} catch (NumberFormatException e) {
					snp.chr = chr;
				}
				snp.pos = Integer.parseInt(line[2]);
			}
		}
		
		return map;
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;

		String usage = "FAILED";
		String pop = null; 
		
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h")
					|| args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("pop=")) {
				pop = args[i].endsWith("EA") ? EA : AA;
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0 || pop == null) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			new CARDIA2017ResultsProcessor().process(pop);
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	
}
