package org.genvisis.one.ben;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.gwas.parsing.AbstractColumnFilter;
import org.genvisis.gwas.parsing.AbstractFileColumn;
import org.genvisis.gwas.parsing.AliasedFileColumn;
import org.genvisis.gwas.parsing.ColumnFilter;
import org.genvisis.gwas.parsing.FileColumn;
import org.genvisis.gwas.parsing.FileParserFactory;
import org.genvisis.gwas.parsing.FixedValueColumn;
import org.genvisis.gwas.parsing.ParseFailureException;
import org.genvisis.gwas.parsing.StandardFileColumns;

public class CARDIA2Processor {


	static String[] IN_HDR = {"name", "A1", "A2", "Freq1", "MAF", "Quality", "Rsq", // 6
														"n", "Mean_predictor_allele", // 8 ----- EffAF??
														"beta_SNP_add", "sebeta_SNP_add", // 10
														"beta_SNP_DPA", "sebeta_SNP_DPA", // 12
														"cov_SNP_int_SNP_DPA", "chi2_SNP" // 14
	};

	static String[] OUT_HDR = {"chr", "pos", "marker_name", "strand", "base_allele", "effect_allele",
														 "N", // IN_7
														 "effect_allele_freq", "imputation_type", "Imputation_value", // IN_6
														 "beta_main", // IN_9
														 "se_main", // IN_10
														 "beta_int", // IN_11
														 "se_int", // IN_12
														 "cov", // IN_13
														 "chi_2df", // IN_14
														 "chi_P_2df",};
	static String STRAND = "+";
	static String TYPE = "1";
	static String posDirTempl_EA = "/home/pankarne/cole0482/probabelCARDIA/EA/info/cardia_impmetrics_chr#.out";
	static String posDirTempl_AA = "/scratch.global/cole0482/CARDIA/data/map/chr#_positions.xln";
	static String resultsTempl = "regression_chr#.out_add.out.txt";

	static int NUM_CHRS_EA = 22;
	static int NUM_CHRS_AA = 22;

	public static void processEA(String baseDir, String resultsDir) {
		Logger log = new Logger();
		String outputFile = baseDir + resultsDir + "combined.results";
		if (!Files.exists(outputFile)) {
			for (int i = 1; i < (NUM_CHRS_EA + 1); i++) {
				try {
					processChr(baseDir, resultsDir, null, i, log, true);
				} catch (IOException e) {
					log.reportError("ERROR PROCESSING EA CHR " + i);
					log.reportException(e);
					throw new RuntimeException();
				}
			}
		}

		String manPlot = Files.getRunString() + " org.genvisis.cnv.plots.ManhattanPlot file="
										 + outputFile + " -screenshot size=2000,1440 screen=" + baseDir + resultsDir
										 + "Manhattan_"
										 + resultsDir.replace('/', '_').substring(0, resultsDir.length() - 1) + ".jpeg";
		String manScr = baseDir + resultsDir + "createManhattanPlot.sh";
		if (!Files.exists(manScr)) {
			Files.write(manPlot, manScr);
			Files.chmod(manScr);
		}

		String qqPlot = Files.getRunString() + " org.genvisis.cnv.plots.QQPlot files=" + outputFile
										+ " outFile=" + baseDir + resultsDir + "QQ_"
										+ resultsDir.replace('/', '_').substring(0, resultsDir.length() - 1) + ".jpeg";

		String qqScr = baseDir + resultsDir + "createQQPlot.sh";
		if (!Files.exists(qqScr)) {
			Files.write(qqPlot, qqScr);
			Files.chmod(qqScr);
		}
	}

	public static void processAA(String baseDir, String resultsDir) {
		Logger log = new Logger();
		String dropsFile = null; // "/scratch.global/cole0482/CARDIA/drops.AA";
		String outputFile = baseDir + resultsDir + "combined.results";
		if (!Files.exists(outputFile)) {
			for (int i = 1; i < (NUM_CHRS_AA + 1); i++) {
				try {
					processChr(baseDir, resultsDir, dropsFile, i, log, false);
				} catch (IOException e) {
					log.reportError("ERROR PROCESSING AA CHR " + i);
					log.reportException(e);
					throw new RuntimeException();
				}
			}
		}

		String manPlot = Files.getRunString() + " org.genvisis.cnv.plots.ManhattanPlot file="
										 + outputFile + " -screenshot size=2000,1440 screen=" + baseDir + resultsDir
										 + "Manhattan_"
										 + resultsDir.replace('/', '_').substring(0, resultsDir.length() - 1) + ".jpeg";
		String manScr = baseDir + resultsDir + "createManhattanPlot.sh";
		if (!Files.exists(manScr)) {
			Files.write(manPlot, manScr);
			Files.chmod(manScr);
		}

		String qqPlot = Files.getRunString() + " org.genvisis.cnv.plots.QQPlot files=" + outputFile
										+ " outFile=" + baseDir + resultsDir + "QQ_"
										+ resultsDir.replace('/', '_').substring(0, resultsDir.length() - 1) + ".jpeg";

		String qqScr = baseDir + resultsDir + "createQQPlot.sh";
		if (!Files.exists(qqScr)) {
			Files.write(qqPlot, qqScr);
			Files.chmod(qqScr);
		}
	}

	private static FileColumn<String> getDoubleFormatterColumn(String name, String alias, int num) {
		return new AliasedFileColumn(name, alias) {
			@Override
			public String getValue(String[] line) throws ParseFailureException {
				String vStr = super.getValue(line);
				try {
					return ext.formDeci(Double.parseDouble(vStr), num);
				} catch (NumberFormatException e) {
					throw new ParseFailureException(e);
				}
			}
		};
	}

	private static void processChr(String baseDir, String resultsDir, String dropsFile, int chr,
																 Logger log, boolean EA) throws IOException {
		Set<String> drops = dropsFile == null
																					? new HashSet<>()
																					: HashVec.loadToHashSet(HashVec.loadFileToStringArray(dropsFile,
																																																false,
																																																new int[] {0},
																																																false));
		String map = (EA ? posDirTempl_EA : posDirTempl_AA).replace("#", Integer.toString(chr));
		log.reportTime("Loading map file: " + map);
		Map<String, String> posMap = HashVec.loadFileColumnToMap(map, 0, 2, EA, log);
		log.reportTime("Loaded " + posMap.size() + " positions.");
		Map<String, String> missMap = new HashMap<>();

		ChiSquaredDistribution csd = new ChiSquaredDistribution(2);

		FileColumn<String> chrC = new FixedValueColumn("chr", Integer.toString(chr));
		FileColumn<String> snp = new AliasedFileColumn("marker_name", "name");
		FileColumn<Integer> pos = new AbstractFileColumn<Integer>("pos") {
			@Override
			public Integer getValue(String[] line) throws ParseFailureException {
				String snpV = snp.getValue(line);
				String posV = posMap.get(snpV);
				if (posV == null) {
					if (snpV.contains(":")) {
						posV = snpV.split(":")[1];
					} else if (snpV.equals("rs55949118;rs60110500")) {
						posV = "24735753";
					} else if (snpV.equals("rs66965879;rs71696953")) {
						posV = "39223254";
					} else {
						log.reportError("Missing position for snp: " + snpV);
					}
				}
				if (posV.equals(".") && missMap.containsKey(snpV)) {
					posV = missMap.get(snpV);
				}
				return Integer.parseInt(posV);
			}

			@Override
			public void initialize(Map<String, Integer> headerMap) {
				// no-op
			}

			@Override
			public int hashCode() {
				final int prime = 31;
				int result = 1;
				result = prime * result + ((getName() == null) ? 0 : getName().hashCode());
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
				FileColumn<?> other = (FileColumn<?>) obj;
				if (getName() == null) {
					if (other.getName() != null)
						return false;
				} else if (!getName().equals(other.getName()))
					return false;
				return true;
			}
		};
		FileColumn<String> strand = new FixedValueColumn("strand", STRAND);
		FileColumn<String> a1 = new AliasedFileColumn("base_allele", "a1");
		FileColumn<String> a2 = new AliasedFileColumn("effect_allele", "a2");
		FileColumn<Integer> n = StandardFileColumns.n("N");
		FileColumn<String> maf = getDoubleFormatterColumn("effect_allele_freq", "Mean_predictor_allele",
																											4);
		FileColumn<String> type = new FixedValueColumn("inputation_type", TYPE);
		FileColumn<String> rsq = new AliasedFileColumn("Rsq", "Imputation_value");
		FileColumn<String> betaMain = getDoubleFormatterColumn("beta_main", "beta_SNP_add", 5);
		FileColumn<String> seMain = getDoubleFormatterColumn("se_main", "sebeta_SNP_add", 5);
		FileColumn<String> betaInt = getDoubleFormatterColumn("beta_int", "beta_SNP_DPA", 5);
		FileColumn<String> seInt = getDoubleFormatterColumn("se_int", "sebeta_SNP_DPA", 5);
		FileColumn<String> cov = new AliasedFileColumn("cov", "cov_SNP_int_SNP_DPA");

		FileColumn<String> pVals = new AliasedFileColumn("chi_2df\tchiP_2df", "chi2_SNP") {
			@Override
			public String getValue(String[] line) throws ParseFailureException {
				String chiStr = super.getValue(line);
				double chi = ext.isMissingValue(chiStr) ? Double.NaN : Double.parseDouble(chiStr);
				double pval = 1 - csd.cumulativeProbability(chi);
				String pvalStr = Double.isNaN(pval) ? "NA"
																						: pval > 0.001 ? ext.formDeci(pval, 4)
																													 : ext.formSciNot(pval, 4, false);
				return (Double.isNaN(chr) ? "NA" : chiStr) + "\t" + pvalStr;
			}
		};

		ColumnFilter dropFilter = new AbstractColumnFilter(snp) {
			@Override
			public boolean filter(Map<FileColumn<?>, String> values) {
				return !drops.contains(values.get(snp));
			}
		};

		String file = baseDir + resultsDir + resultsTempl.replace("#", Integer.toString(chr));
		log.report("Processing results file " + file);
		String outputFile = baseDir + resultsDir + "combined.results";

		FileParserFactory.setup(file, chrC, pos, snp, strand, a1, a2, n, maf, type, rsq, betaMain,
														seMain, betaInt, seInt, cov, pVals)
										 .parseFailValue("NA").filter(dropFilter).build()
										 .parseToFile(outputFile, "\t", chr == 1, chr != 1);
	}

	public static void main(String[] args) {

		boolean ea = false;
		boolean aa = false;

		for (String arg : args) {
			if (arg.startsWith("pops=")) {
				String[] pops = arg.substring("pops=".length()).split(",");
				for (String p : pops) {
					if (p.equalsIgnoreCase("ea")) {
						ea = true;
					}
					if (p.equalsIgnoreCase("aa")) {
						aa = true;
					}
				}
			}
		}

		String dir = "/scratch.global/cole0482/CARDIA/GxE/";
		String[] resultsDirsEA = {
															"EA/FEV1/DHA/",
															"EA/FEV1/DPA/",
															"EA/FEV1/EPA/",
															"EA/FEV1/Fish/",
															"EA/FVC/DHA/",
															"EA/FVC/DPA/",
															"EA/FVC/EPA/",
															"EA/FVC/Fish/",
		};
		String[] resultsDirsAA = {
															"AA/FEV1/DHA/",
															"AA/FEV1/DPA/",
															"AA/FEV1/EPA/",
															"AA/FEV1/Fish/",
															"AA/FVC/DHA/",
															"AA/FVC/DPA/",
															"AA/FVC/EPA/",
															"AA/FVC/Fish/",
		};
		if (aa) {
			for (String d : resultsDirsAA) {
				processAA(dir, d);
			}
			System.out.println("Done with AA!");
		}
		if (ea) {
			for (String d : resultsDirsEA) {
				processEA(dir, d);
			}
		}
	}

}
