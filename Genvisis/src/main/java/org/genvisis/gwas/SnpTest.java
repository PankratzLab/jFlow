package org.genvisis.gwas;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.genvisis.CLI;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.gwas.parsing.AliasedFileColumn;
import org.genvisis.gwas.parsing.DoubleWrapperColumn;
import org.genvisis.gwas.parsing.FileColumn;
import org.genvisis.gwas.parsing.FileParserFactory;
import org.genvisis.gwas.parsing.FixedValueColumn;
import org.genvisis.gwas.parsing.StandardFileColumns;
import org.genvisis.qsub.Qsub;


public class SnpTest {

	class SnpTestCommand {
		private String snpTestLocation;
		private String dataFile;
		private String sampleFile;
		private boolean isVCFData;
		private String vcfDataField = "GP";
		private String phenoName;
		private String[] covars = null;
		private String outputFile;

		public SnpTestCommand(String snpTestLocation, String dataFile, String sampleFile,
													boolean isVCFData, String vcfDataField, String phenoName,
													String[] covars, int chr, String outputFile) {
			this.snpTestLocation = snpTestLocation;
			this.dataFile = dataFile;
			this.sampleFile = sampleFile;
			this.isVCFData = isVCFData;
			this.vcfDataField = vcfDataField;
			this.phenoName = phenoName;
			this.covars = covars;
			this.outputFile = outputFile;
		}

		private String getCommand() {
			StringBuilder snpTestString = new StringBuilder(snpTestLocation)
																																			.append(" -data ")
																																			.append(dataFile)
																																			.append(" ")
																																			.append(sampleFile)
																																			.append(" -frequentist 1 ")
																																			.append(" -method expected ")
																																			.append(" -use_raw_phenotypes ");
			if (isVCFData) {
				snpTestString.append(" -genotype_field ").append(vcfDataField);
			}
			snpTestString.append(" -pheno ").append(phenoName);
			if (covars != null) {
				snpTestString.append(" -cov_names");
				for (String cov : covars) {
					snpTestString.append(" ").append(cov);
				}
			} else {
				snpTestString.append(" -cov_all");
			}
			snpTestString.append(" lower_sample_limit");
			snpTestString.append(" -o ").append(outputFile);

			return snpTestString.toString();
		}
	}

	private static FileColumn<?>[] getFileColumns(String sampleCount) {
		return new FileColumn<?>[] {
																StandardFileColumns.snp("SNP"),
																StandardFileColumns.chr("CHR"),
																StandardFileColumns.pos("POS"),
																new FixedValueColumn("STRAND", "FWD"),
																StandardFileColumns.a1("EFFECT_ALLELE"),
																StandardFileColumns.a2("OTHER_ALLELE"),
																new FixedValueColumn("N", sampleCount),
																new DoubleWrapperColumn(new AliasedFileColumn("EAF",
																																							new String[] {"all_maf"})),
																new DoubleWrapperColumn(new AliasedFileColumn("BETA",
																																							new String[] {"frequentist_add_beta_1"})),
																new DoubleWrapperColumn(new AliasedFileColumn("SE",
																																							new String[] {"frequentist_add_se_1"})),
																new DoubleWrapperColumn(new AliasedFileColumn("PVAL",
																																							new String[] {"frequentist_add_pvalue"}))
		};
	}

	private static void processResults(String dir) throws IOException {
		Logger log = new Logger();
		String[] sampleFiles = Files.list(dir, ".sample");
		final String sampleCount = sampleFiles.length > 1
																											? "N"
																											: Integer.toString(Files.countLines(dir
																																													+ sampleFiles[0],
																																													2));
		if (sampleCount.equals("N")) {
			log.reportTimeWarning("Could not determine sample file used.  Number of samples will be set to 'N' in results file.");
		}

		String outputTempl = "output_chr#.out";
		int fileInd = 0;
		for (int i = 1; i < 28; i++) {
			String outFile = outputTempl.replace("#", Integer.toString(i));
			if (Files.exists(dir + outFile)) {
				FileParserFactory.setup(dir + outFile, getFileColumns(sampleCount))
												 .skipPrefix("#").build()
												 .parseToFile(dir + "combined.results", "\t", fileInd == 0, fileInd != 0);
				fileInd++;
			}
		}
	}


	Logger log = new Logger();

	String snpTestExec;
	String dataDirectory;
	String dataFileTemplate;
	String dataFileExtension;
	String repl = "#";
	String sampleFile;
	String pheno;
	String[] covars = null;
	boolean isVCFOverride = false;
	String vcfField = "GP";

	public void setDataDirectory(String dataDirectory) {
		this.dataDirectory = ext.verifyDirFormat(dataDirectory);
	}

	public void setDataFileTemplate(String dataFileTemplate) {
		this.dataFileTemplate = dataFileTemplate;
	}

	public void setDataFileExtension(String dataFileExtension) {
		this.dataFileExtension = dataFileExtension;
	}

	public void setRepl(String repl) {
		this.repl = repl;
	}

	public void setSampleFile(String sampleFile) {
		this.sampleFile = sampleFile;
	}

	public void setPheno(String pheno) {
		this.pheno = pheno;
	}

	public void setCovars(String[] covars) {
		this.covars = covars;
	}

	private void checkNeeds() throws IllegalStateException {
		if (snpTestExec == null) {
			throw new IllegalStateException("SnpTest executable not specified; set using \""
																			+ ARG_SNPTEST + "\" argument and try again.");
		}
		if (!Files.exists(snpTestExec)) {
			throw new IllegalStateException("SnpTest executable not found in specified location: \""
																			+ snpTestExec + "\"");
		}
		if (dataDirectory == null) {
			throw new IllegalStateException("Data file directory not specified; set using \""
																			+ ARG_DATADIR + "\" argument and try again.");
		}
		if (!Files.exists(dataDirectory)) {
			throw new IllegalStateException("Specified data file directory not found: \"" + dataDirectory
																			+ "\"");
		}
		if (dataFileTemplate != null && !dataFileTemplate.contains(repl)) {
			throw new IllegalStateException(
																			"Data file template \""
																			+ dataFileTemplate
																			+ "\" does not contain the specified special character \" + repl + \".");
		}
		if (dataFileExtension != null) {
			if (new File(dataDirectory).list((d, f) -> {
				return f.endsWith(dataFileExtension);
			}).length == 0) {
				throw new IllegalStateException("No files found with extension \"" + dataFileExtension
																				+ "\" in specified data directory: \"" + dataDirectory
																				+ "\".");
			}
		}
		if (sampleFile == null) {
			throw new IllegalStateException("Sample file not specified; set using \"" + ARG_SAMPLE
																			+ "\" argument and try again.");
		}
		if (!Files.exists(sampleFile)) {
			throw new IllegalStateException("Specified sample file not found: " + sampleFile + "\"");
		}
		if (pheno == null) {
			throw new IllegalStateException("Phenotype variable not specified; set using \"" + ARG_PHENO
																			+ "\" argument and try again.");
		}
		String[] hdr = null;
		if (pheno.equals("")
				|| ext.indexOfStr(pheno, hdr = Files.getHeaderOfFile(sampleFile, log)) == -1) {
			throw new IllegalStateException("Invalid phenotype specified: \""
																			+ pheno
																			+ "\"."
																			+ (hdr == null ? "" : " Sample file header: "
																														+ ArrayUtils.toStr(hdr)));
		}
	}

	private Map<Integer, String> loadFiles() {
		Map<Integer, String> returnFiles = new HashMap<>();
		if (dataFileTemplate != null) {
			for (int i = 1; i < 27; i++) {
				String dF = dataDirectory + dataFileTemplate.replace(repl, Integer.toString(i));
				if (!Files.exists(dF)) {
					log.reportTimeWarning("No data file found for chr" + i + "; expected " + dF);
					continue;
				}
				returnFiles.put(i, dF);
			}
		} else if (dataFileExtension != null) {
			Pattern chrPattern = Pattern.compile(".*chr([0-2[XYM]]?[0-9[Y]]*).*");
			String[] files = new File(dataDirectory).list((d, f) -> {
				return f.endsWith(dataFileExtension);
			});
			for (String f : files) {
				Matcher m = chrPattern.matcher(f);
				if (m.matches()) {
					returnFiles.put(Integer.valueOf(m.group(1)), f);
				} else {
					log.reportError("Couldn't determined chr number from filename: \"" + f + "\".");
				}
			}
		}
		return returnFiles;
	}

	public void run() {
		checkNeeds();
		Map<Integer, String> dataFiles = loadFiles();
		List<String> cmds = new ArrayList<>();
		List<String> cmdOuts = new ArrayList<>();
		for (Entry<Integer, String> file : dataFiles.entrySet()) {
			String outFile = "./output_chr" + file.getKey() + ".out";
			cmds.add(new SnpTestCommand(snpTestExec, file.getValue(), sampleFile,
																	file.getValue().contains(".vcf") || isVCFOverride, vcfField,
																	pheno, covars, file.getKey(), outFile).getCommand());
			cmdOuts.add(outFile);
		}
		Files.writeIterable(cmds, "./input.txt");

		StringBuilder scriptExecCmd = new StringBuilder("cd ").append(ext.pwd())
																													.append("\n")
																													.append(Files.getRunString())
																													.append(" org.genvisis.one.ScriptExecutor token=finito threads="
																																	+ Runtime.getRuntime()
																																					 .availableProcessors());

		Qsub.qsubDefaults("./runSnpTest.qsub", scriptExecCmd.toString());

		String parseScr = "./parseSnpTest.sh";
		Files.write(Files.getRunString() + " org.genvisis.gwas.SnpTest dir=" + ext.pwd() + " -parse",
								parseScr);
		Files.chmod(parseScr);
	}

	private static final String ARG_SNPTEST = "snpTest";
	private static final String ARG_DATADIR = "dataDir";
	private static final String ARG_DATA = "data";
	private static final String ARG_DATAEXT = "dataExt";
	private static final String ARG_SAMPLE = "sample";
	private static final String ARG_PHENO = "pheno";
	private static final String ARG_COVARS = "covars";
	private static final String ARG_REPL = "repl";
	private static final String ARG_VCFOVERRIDE = "vcfOverride";
	private static final String ARG_VCFFIELD = "vcfField";

	private static final String DESC_SNPTEST = "SnpTest executable (full path included unless on path)";
	private static final String DESC_DATADIR = "Data file directory";
	private static final String DESC_DATA = "Data filename template (e.g. chr#.vcf.gz) - assumes chromosomally-separated data files.";
	private static final String DESC_DATAEXT = "Extension of data files in data directory.";
	private static final String DESC_SAMPLE = "Sample file (SnpTest-format)";
	private static final String DESC_PHENO = "Name of phenotype to test";
	private static final String DESC_COVARS = "Covariates to test, leave blank to use all available (uses the -cov_all flag)";
	private static final String DESC_REPL = "Special character in data file template into which chromosome number will be placed.";
	private static final String DESC_VCFOVERRIDE = "Override to specify if the data files are or are not VCF files";
	private static final String DESC_VCFFIELD = "Genotype field to use in VCF files (defaults to 'GP' / genotype probabilities)";

	public static void main(String[] args) throws IOException {
		CLI cli = new CLI(SnpTest.class);

		boolean parse = false;
		String parseDir = null;
		for (String arg : args) {
			if (arg.startsWith("-parse")) {
				parse = true;
			}
			if (arg.startsWith("dir=")) {
				parseDir = arg.split("=")[1];
			}
		}
		if ((parse && parseDir == null) || (!parse && parseDir != null)) {
			throw new RuntimeException(
																 "Invalid parsing options; both -parse flag and dir= argument must be set!");
		} else if (parse && parseDir != null) {
			processResults(parseDir);
			return;
		}

		cli.addArg(ARG_SNPTEST, DESC_SNPTEST);
		cli.addArg(ARG_DATADIR, DESC_DATADIR);
		cli.addArg(ARG_DATA, DESC_DATA);
		cli.addArg(ARG_DATAEXT, DESC_DATAEXT);
		cli.addArg(ARG_SAMPLE, DESC_SAMPLE);
		cli.addArg(ARG_PHENO, DESC_PHENO);

		cli.addArg(ARG_COVARS, DESC_COVARS, false);
		cli.addArg(ARG_REPL, DESC_REPL, false);
		cli.addArg(ARG_VCFOVERRIDE, DESC_VCFOVERRIDE, false);
		cli.addArg(ARG_VCFFIELD, DESC_VCFFIELD, false);

		cli.addGroup(ARG_DATA, ARG_DATAEXT);

		cli.parseWithExit(args);

		String snpTest = cli.get(ARG_SNPTEST);
		String dataDir = cli.get(ARG_DATADIR);
		String data = cli.has(ARG_DATA) ? cli.get(ARG_DATA) : null;
		String ext = cli.has(ARG_DATAEXT) ? cli.get(ARG_DATAEXT) : null;
		String samp = cli.get(ARG_SAMPLE);
		String phen = cli.get(ARG_PHENO);

		SnpTest run = new SnpTest();
		run.snpTestExec = snpTest;
		run.setDataDirectory(dataDir);
		if (data != null) {
			run.setDataFileTemplate(data);
		} else if (ext != null) {
			run.setDataFileExtension(ext);
		}
		run.setSampleFile(samp);
		run.setPheno(phen);

		if (cli.has(ARG_COVARS)) {
			run.setCovars(cli.get(ARG_COVARS).split(","));
		}
		if (cli.has(ARG_REPL)) {
			run.setRepl(cli.get(ARG_REPL));
		}
		if (cli.has(ARG_VCFOVERRIDE)) {
			run.isVCFOverride = Boolean.parseBoolean(cli.get(ARG_VCFOVERRIDE));
		}
		if (cli.has(ARG_VCFFIELD)) {
			run.vcfField = cli.get(ARG_VCFFIELD);
		}

		run.run();
	}

}
