package org.pankratzlab.internal.gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.CLI;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.parsing.StandardFileColumns;
import org.pankratzlab.common.qsub.Qsub;
import org.pankratzlab.fileparser.AliasedFileColumn;
import org.pankratzlab.fileparser.DoubleWrapperColumn;
import org.pankratzlab.fileparser.FileColumn;
import org.pankratzlab.fileparser.FileLink;
import org.pankratzlab.fileparser.FileParserFactory;
import org.pankratzlab.fileparser.FixedValueColumn;

import htsjdk.variant.vcf.VCFFileReader;

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
    private String excludeFile;
    private String includeFile;

    public SnpTestCommand(String snpTestLocation, String dataFile, String sampleFile,
                          boolean isVCFData, String vcfDataField, String phenoName, String[] covars,
                          int chr, String outputFile, String excludeFile, String includeFile) {
      this.snpTestLocation = snpTestLocation;
      this.dataFile = dataFile;
      this.sampleFile = sampleFile;
      this.isVCFData = isVCFData;
      this.vcfDataField = vcfDataField;
      this.phenoName = phenoName;
      this.covars = covars;
      this.outputFile = outputFile;
      this.excludeFile = excludeFile;
      this.includeFile = includeFile;
    }

    private String getCommand() {
      StringBuilder snpTestString = new StringBuilder(snpTestLocation).append(" -data ")
                                                                      .append(dataFile).append(" ")
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
        log.report("No covariates specified, using all available covariates.");
        snpTestString.append(" -cov_all");
      }
      snpTestString.append(" -lower_sample_limit 10");
      snpTestString.append(" -o ").append(outputFile);
      if (excludeFile != null) {
        snpTestString.append(" -exclude_samples ").append(excludeFile);
      }
      if (includeFile != null) {
        snpTestString.append(" -include_samples ").append(includeFile);
      }
      snpTestString.append(" -chunk 200 ");

      return snpTestString.toString();
    }
  }

  private static FileColumn<?>[] getFileColumns(String sampleCount) {
    return new FileColumn<?>[] {StandardFileColumns.snp("SNP"), StandardFileColumns.chr("CHR"),
                                StandardFileColumns.pos("POS"),
                                new FixedValueColumn("STRAND", "FWD"),
                                new AliasedFileColumn("EFFECT_ALLELE", "alleleB"),
                                new AliasedFileColumn("OTHER_ALLELE", "alleleA"),
                                new FixedValueColumn("N", sampleCount),
                                new DoubleWrapperColumn(new AliasedFileColumn("BETA",
                                                                              new String[] {"frequentist_add_beta_1"})),
                                new DoubleWrapperColumn(new AliasedFileColumn("SE",
                                                                              new String[] {"frequentist_add_se_1"})),
                                new DoubleWrapperColumn(new AliasedFileColumn("PVAL",
                                                                              new String[] {"frequentist_add_pvalue"})),};
  }

  private static List<FileColumn<?>> getOutputOrder(FileColumn<?>[] resultsColumns,
                                                    FileColumn<String> freqCol,
                                                    FileColumn<String> rsqCol) {
    List<FileColumn<?>> order = new ArrayList<>();
    order.add(resultsColumns[0]);
    order.add(resultsColumns[1]);
    order.add(resultsColumns[2]);
    order.add(resultsColumns[3]);
    order.add(resultsColumns[4]);
    order.add(resultsColumns[5]);
    order.add(resultsColumns[6]);
    order.add(freqCol);
    order.add(resultsColumns[7]);
    order.add(resultsColumns[8]);
    order.add(resultsColumns[9]);
    order.add(rsqCol);
    return order;
  }

  private static void processResults(String dir, String infoDirAndTempl) throws IOException {
    Logger log = new Logger();
    String[] sampleFiles = Files.list(dir, ".sample");
    final String sampleCount = sampleFiles.length > 1 ? "N"
                                                      : Integer.toString(Files.countLines(dir
                                                                                          + sampleFiles[0],
                                                                                          2));
    if (sampleCount.equals("N")) {
      log.reportTimeWarning("Could not determine sample file used.  Number of samples will be set to 'N' in results file.");
    }

    String outputTempl = "output_chr#.out.gz";
    int fileInd = 0;
    for (int i = 1; i < 28; i++) {
      String outFile = outputTempl.replace("#", Integer.toString(i));
      if (Files.exists(dir + outFile)) {
        FileColumn<?>[] resultsColumns = getFileColumns(sampleCount);
        AliasedFileColumn snpCol = StandardFileColumns.snp("SNP");
        FileColumn<String> freqCol = new AliasedFileColumn("EAF", "ALT_frq");
        FileColumn<String> rsqCol = new AliasedFileColumn("Rsq", "Rsq");

        List<FileColumn<?>> outputOrder = getOutputOrder(resultsColumns, freqCol, rsqCol);

        FileParserFactory.setup(dir + outFile, resultsColumns).skipPrefix("#")
                         .link(FileLink.setup(infoDirAndTempl.replace("#", Integer.toString(i)))
                                       .keys(snpCol).values(freqCol, rsqCol))
                         .build().parseToFile(dir + "combined.results", "\t", fileInd == 0,
                                              fileInd != 0, outputOrder);
        fileInd++;
      }
    }
  }

  Logger log = new Logger();

  String snpTestExec;
  String dataDirectory;
  String dataFileTemplate;
  String infoFileTemplate;
  String dataFileExtension;
  String repl = "#";
  String sampleFile;
  String excludesFile;
  String includesFile;
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

  public void setExcludes(String excludesFile) {
    this.excludesFile = excludesFile;
  }

  public void setIncludes(String includesFile) {
    this.includesFile = includesFile;
  }

  private void checkNeeds() throws IllegalStateException {
    if (snpTestExec == null) {
      throw new IllegalStateException("SnpTest executable not specified; set using \"" + ARG_SNPTEST
                                      + "\" argument and try again.");
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
      throw new IllegalStateException("Data file template \"" + dataFileTemplate
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
      throw new IllegalStateException("Invalid phenotype specified: \"" + pheno + "\"."
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

  public static void generateSampleSheet(String vcfPath, int numCovs, String covFile,
                                         String sexFile, String phenoFile) {
    VCFFileReader vcfFileReader = new VCFFileReader(new File(vcfPath), false);
    ArrayList<String> samples = vcfFileReader.getFileHeader().getSampleNamesInOrder();

    String outFile = "./snptest_sample_sheet.txt";
    String covHeader = "";
    String covHeaderTwo = "";
    int sampleSex = 0;
    String samplePheno;
    ArrayList<Double> currentCovs;
    HashMap<String, ArrayList<Double>> covLookup = null;
    HashMap<String, Integer> sexLookup = null;
    HashMap<String, String> phenoLookup = null;
    try {
      System.out.println("Building cov map");
      covLookup = buildCovMap(covFile);
      System.out.println("Building sex map");
      sexLookup = buildSexMap(sexFile);
      System.out.println("Building pheno map");
      phenoLookup = buildPhenoMap(phenoFile);
    } catch (IOException e) {
      e.printStackTrace();
    }
    for (int i = 1; i <= numCovs; i++) {
      covHeader += "cov" + i + "\t";
      covHeaderTwo += "C" + "\t";
    }

    String header = "ID_1" + "\t" + "ID_2" + "\t" + "missing" + "\t" + "sex" + "\t" + covHeader
                    + "phenotype" + "\t" + "race" + "\t";
    String headerTwo = "0" + "\t" + "0" + "\t" + "0" + "\t" + "D" + "\t" + covHeaderTwo + "B" + "\t"
                       + "D";

    PrintWriter writer = Files.getAppropriateWriter(outFile);
    System.out.println("Writing to file");
    writer.write(header + "\n" + headerTwo + "\n");
    if (covLookup != null && sexLookup != null && phenoLookup != null) {
      for (String sample : samples) {
        samplePheno = phenoLookup.get(sample);
        writer.write(sample + "\t" + sample + "\t" + "NA" + "\t");
        if (sexLookup.containsKey(sample)) {
          writer.write(sexLookup.get(sample) + "\t");
        } else {
          writer.write("NA" + "\t");
        }
        for (int i = 0; i < numCovs; i++) {
          currentCovs = covLookup.get(sample);
          if (currentCovs != null) {
            writer.write(currentCovs.get(i) + "\t");
          } else {
            writer.write("NA" + "\t");
          }
        }
        if (phenoLookup.containsKey(sample)) {
          writer.write(samplePheno + "\n");
        } else {
          writer.write("NA" + "\t");
        }
      }
    } else
      throw new IllegalStateException();
    vcfFileReader.close();
    writer.close();
  }

  private static HashMap<String, ArrayList<Double>> buildCovMap(String covFile) throws IOException {
    HashMap<String, ArrayList<Double>> sampleCovMap = new HashMap<String, ArrayList<Double>>();

    try (BufferedReader covReader = Files.getAppropriateReader(covFile)) {
      covReader.readLine();
      String line = covReader.readLine();
      String[] lineArray;
      while (line != null) {
        lineArray = line.trim().split("\t");
        ArrayList<Double> pcs = new ArrayList<Double>();
        for (int i = 2; i < lineArray.length; i++) {
          try {
            pcs.add(Double.parseDouble(lineArray[i]));
          } catch (NumberFormatException nfe) {
            // improper file pc format
            nfe.printStackTrace();
          }

        }
        sampleCovMap.put(lineArray[0], pcs);
        line = covReader.readLine();
      }
    }

    return sampleCovMap;
  }

  private static HashMap<String, Integer> buildSexMap(String sexChecksFile) throws IOException {
    HashMap<String, Integer> sampleSexMap = new HashMap<String, Integer>();

    try (BufferedReader sexReader = Files.getAppropriateReader(sexChecksFile)) {
      sexReader.readLine();
      String line;
      String[] lineArray;
      String fid = null;
      int sex = 0;
      while ((line = sexReader.readLine()) != null) {
        lineArray = line.trim().split("\t");
        try {
          sex = Integer.parseInt(lineArray[3]);
        } catch (NumberFormatException nfe) {
          nfe.printStackTrace();
        }
        fid = lineArray[1];
        sampleSexMap.put(fid, sex);
      }
    }
    return sampleSexMap;
  }

  private static HashMap<String, String> buildPhenoMap(String phenoFile) throws IOException {
    // phenoFile input : Fid\tIid\tpheno (binary e.g. case or control)
    HashMap<String, String> samplePhenoMap = new HashMap<String, String>();

    try (BufferedReader phenoReader = Files.getAppropriateReader(phenoFile)) {
      String line;
      String[] lineArray;
      String fid = null;
      String pheno = null;
      while ((line = phenoReader.readLine()) != null) {
        lineArray = line.trim().split("\t");
        pheno = lineArray[2];
        fid = lineArray[0];
        samplePhenoMap.put(fid, pheno);
      }
    }
    return samplePhenoMap;
  }

  public void run() {
    checkNeeds();
    Map<Integer, String> dataFiles = loadFiles();
    List<String> cmds = new ArrayList<>();
    List<String> cmdOuts = new ArrayList<>();
    for (Entry<Integer, String> file : dataFiles.entrySet()) {
      String outFile = "./output_chr" + file.getKey() + ".out";
      String excludes = this.excludesFile;
      String includes = this.includesFile;
      cmds.add(new SnpTestCommand(snpTestExec, file.getValue(), sampleFile,
                                  file.getValue().contains(".vcf")
                                                                            || isVCFOverride,
                                  vcfField, pheno, covars, file.getKey(), outFile, excludes,
                                  includes).getCommand());

      cmdOuts.add(outFile);
    }
    Files.writeIterable(cmds, "./input.txt");

    StringBuilder scriptExecCmd = new StringBuilder("cd ").append(ext.pwd()).append("\n")
                                                          .append(Files.getRunString())
                                                          .append(" org.genvisis.one.ScriptExecutor token=finito threads="
                                                                  + Runtime.getRuntime()
                                                                           .availableProcessors());

    Qsub.qsubDefaults("./runSnpTest.qsub", scriptExecCmd.toString());

    String parseScr = "./parseSnpTest.sh";
    Files.write(Files.getRunString() + " org.genvisis.gwas.SnpTest dir=" + ext.pwd() + " info="
                + dataDirectory + infoFileTemplate + " -parse", parseScr);
    Files.chmod(parseScr);
  }

  private static final String ARG_SNPTEST = "snpTest";
  private static final String ARG_DATADIR = "dataDir";
  private static final String ARG_DATA = "data";
  private static final String ARG_INFOTEMP = "infoTemplate";
  private static final String ARG_DATAEXT = "dataExt";
  private static final String ARG_SAMPLE = "sample";
  private static final String ARG_PHENO = "pheno";
  private static final String ARG_COVARS = "covars";
  private static final String ARG_REPL = "repl";
  private static final String ARG_VCFOVERRIDE = "vcfOverride";
  private static final String ARG_VCFFIELD = "vcfField";
  private static final String ARG_EXCLUDES = "exclude";
  private static final String ARG_INCLUDES = "include";
  private static final String ARG_COVS = "covs";
  private static final String ARG_SEX = "sex";
  private static final String ARG_PHENO_LIST = "phenoList";
  private static final String ARG_NUM_COV = "numCovs";
  private static final String ARG_VCF_PATH = "vcfPath";

  private static final String DESC_SNPTEST = "SnpTest executable (full path included unless on path)";
  private static final String DESC_DATADIR = "Data file directory";
  private static final String DESC_DATA = "Data filename template (e.g. chr#.vcf.gz) - assumes chromosomally-separated data files.";
  private static final String DESC_INFOTEMP = "Info filename template (e.g. chr#.info.gz) - assumes chromosomally-separated data files.";
  private static final String DESC_DATAEXT = "Extension of data files in data directory.";
  private static final String DESC_SAMPLE = "Sample file (SnpTest-format)";
  private static final String DESC_PHENO = "Name of phenotype to test";
  private static final String DESC_COVARS = "Covariates to test, leave blank to use all available (uses the -cov_all flag)";
  private static final String DESC_REPL = "Special character in data file template into which chromosome number will be placed.";
  private static final String DESC_VCFOVERRIDE = "Override to specify if the data files are or are not VCF files";
  private static final String DESC_VCFFIELD = "Genotype field to use in VCF files (defaults to 'GP' / genotype probabilities)";
  private static final String DESC_EXCLUDES = "File with a list of samples that should be excluded from analysis. These should match samples from the sample file.";
  private static final String DESC_INCLUDES = "File with a list of samples that should be included. Samples not in this file will be excluded.";
  private static final String DESC_COVS = "A covariates file to be used in generating the snptest sample sheet.";
  private static final String DESC_SEX = "A sex checks file to be used in generating the snptest sample sheet.";
  private static final String DESC_PHENO_LIST = "A file with a fid iid and phenotype to be used in generating the snptest sample sheet.";
  private static final String DESC_NUM_COV = "The number of covariates to be included in the snptest sample sheet";
  private static final String DESC_VCF_PATH = "The path to the vcf files. Used in -generate.";

  public static void main(String[] args) throws IOException {
    CLI cli = new CLI(SnpTest.class);

    boolean parse = false;
    String parseDir = null;
    String infoDirAndTemp = null;
    boolean generate = false;
    String sexFile = null;
    String phenoFile = null;
    String covFile = null;
    String vcfPath = null;
    int numCov = 0;
    for (String arg : args) {
      if (arg.startsWith("-parse")) {
        parse = true;
      }
      if (arg.startsWith("dir=")) {
        parseDir = arg.split("=")[1];
      }
      if (arg.startsWith("info=")) {
        infoDirAndTemp = arg.split("=")[1];
      }
      if (arg.startsWith("-generate")) {
        generate = true;
      }
      if (arg.startsWith("sex=")) {
        sexFile = arg.split("=")[1];
      }
      if (arg.startsWith("phenoList=")) {
        phenoFile = arg.split("=")[1];
      }
      if (arg.startsWith("covFile=")) {
        covFile = arg.split("=")[1];
      }
      if (arg.startsWith("vcfPath=")) {
        vcfPath = arg.split("=")[1];
      }
      if (arg.startsWith("numCovs=")) {
        numCov = Integer.parseInt(arg.split("=")[1]);
      }
    }
    boolean invalidParse = false;
    if (parse) {
      if (parseDir == null || infoDirAndTemp == null) {
        invalidParse = true;
      }
    } else if (generate) {

    } else {
      if (parseDir != null || infoDirAndTemp != null) {
        invalidParse = true;
      }
    }
    if (invalidParse) {
      throw new RuntimeException("Invalid parsing options; -parse flag, dir=, and info= arguments must be set!");
    } else if (parse && parseDir != null && infoDirAndTemp != null) {
      processResults(parseDir, infoDirAndTemp);
      return;
    }

    if (generate) {
      generateSampleSheet(vcfPath, numCov, covFile, sexFile, phenoFile);
      return;
    }

    cli.addArg(ARG_SNPTEST, DESC_SNPTEST);
    cli.addArg(ARG_DATADIR, DESC_DATADIR);
    cli.addArg(ARG_DATA, DESC_DATA);
    cli.addArg(ARG_INFOTEMP, DESC_INFOTEMP);
    cli.addArg(ARG_DATAEXT, DESC_DATAEXT);
    cli.addArg(ARG_SAMPLE, DESC_SAMPLE);
    cli.addArg(ARG_PHENO, DESC_PHENO);

    cli.addArg(ARG_COVARS, DESC_COVARS, false);
    cli.addArg(ARG_REPL, DESC_REPL, false);
    cli.addArg(ARG_VCFOVERRIDE, DESC_VCFOVERRIDE, false);
    cli.addArg(ARG_VCFFIELD, DESC_VCFFIELD, false);
    cli.addArg(ARG_EXCLUDES, DESC_EXCLUDES, false);
    cli.addArg(ARG_INCLUDES, DESC_INCLUDES, false);
    cli.addArg(ARG_COVS, DESC_COVS, false);
    cli.addArg(ARG_SEX, DESC_SEX, false);
    cli.addArg(ARG_PHENO_LIST, DESC_PHENO_LIST, false);
    cli.addArg(ARG_VCF_PATH, DESC_VCF_PATH, false);
    cli.addArg(ARG_NUM_COV, DESC_NUM_COV, false);

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
    if (cli.has(ARG_EXCLUDES)) {
      run.setExcludes(cli.get(ARG_EXCLUDES));
    }
    if (cli.has(ARG_INCLUDES)) {
      run.setIncludes(cli.get(ARG_INCLUDES));
    }

    run.run();

  }

}
