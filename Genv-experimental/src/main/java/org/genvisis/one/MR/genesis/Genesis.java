package org.genvisis.one.MR.genesis;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

import org.apache.commons.cli.ParseException;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.CLI;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.qsub.Qsub;
import org.pankratzlab.common.stats.Rscript;

public class Genesis {

  private static int batchChrs(String pheno, String libpath, String gds, int chr, String dir) {
    String batchDir = dir + "/batches/";
    if ((!Files.exists(batchDir)) || (!Files.isDirectory(batchDir))) {
      new File(batchDir).mkdirs();
    }
    String resultsDir = dir + "/results/";
    if ((!Files.exists(resultsDir)) || (!Files.isDirectory(resultsDir))) {
      new File(resultsDir).mkdirs();
    }
    List<String> commands = new ArrayList<String>();
    for (int i = 1; i <= 24; i++) {
      String rCode = ".libPaths(\"/home/pankrat2/mstimson/libs\")\r\n";
      rCode = rCode + "source(\"https://bioconductor.org/biocLite.R\")\r\n";
      rCode = rCode + "biocLite(\"GENESIS\")\r\n";
      rCode = rCode + "library(GENESIS)\r\n";
      rCode = rCode + "library(GWASTools)\r\n";
      rCode = rCode + "\r\n";
      rCode = rCode + "load(\"" + dir + "/" + pheno + "_nullmod.RData\")\r\n";
      rCode = rCode + "\r\n";
      rCode = rCode + "genoReader <- GdsGenotypeReader(filename=\""
              + ext.replaceAllWith(gds, "#", new StringBuilder(String.valueOf(chr)).toString())
              + "\")\r\n";
      rCode = rCode + "genos <- GenotypeData(genoReader)\r\n";
      rCode = rCode + "\r\n";
      rCode = rCode + "interval <- as.integer(nsnp(genos)/24)\r\n";
      rCode = rCode
              + "assoc <- assocTestMM(genoData=genos, nullMMobj=nullmod, test=\"Wald\", snp.include=c(getSnpID(genos, interval * "
              + (i - 1) + " + 1" + "):getSnpID(genos, "
              + (i == 24 ? "nsnp(genos)" : new StringBuilder("interval * ").append(i).toString())
              + ")))\r\n";
      rCode = rCode + "\r\n";
      rCode = rCode + "write.table(assoc, \"" + resultsDir + pheno + "_chr" + chr + "-" + i
              + ".csv\", sep=\",\")";
      try {
        if (!Files.exists(resultsDir + pheno + "_chr" + chr + "-" + i + ".csv")) {
          String filename = batchDir + pheno + "_chr" + chr + "-" + i + ".R";
          PrintWriter out = new PrintWriter(new FileOutputStream(filename));
          out.println(rCode);
          out.close();

          commands.add(Rscript.getRscriptExecutable(new Logger()) + " --no-save " + filename);
        }
      } catch (IOException e) {
        e.printStackTrace();
      }
    }
    Qsub.qsubExecutor(batchDir, commands, null, batchDir + "run_" + pheno + "_chr" + chr, 16, 24000,
                      96.0D);
    return commands.size();
  }

  private static void nullmod(String phenoFile, String pheno, String[] covars, String libpath,
                              String relMatrix, String dir) {
    String rCode = ".libPaths(\"" + libpath + "\")\r\n";
    rCode = rCode + "source(\"https://bioconductor.org/biocLite.R\")\r\n";
    rCode = rCode + "biocLite(\"GENESIS\")\r\n";
    rCode = rCode + "library(GENESIS)\r\n";
    rCode = rCode + "library(GWASTools)\r\n";
    rCode = rCode + "\r\n";
    rCode = rCode + "data <- read.table(\"" + dir + "/" + phenoFile
            + "\", sep=\"\\t\", header=T)\r\n";
    rCode = rCode + "row.names(data) <- data$scanID\r\n\r\n";
    rCode = rCode + "data <- data[order(data$scanID),]\r\n";
    rCode = rCode + "load(\"" + relMatrix + "\")\r\n";
    rCode = rCode + "\r\n";

    rCode = rCode
            + "mat <- mat[row.names(mat) %in% data$scanID, colnames(mat) %in% data$scanID]\r\n";
    rCode = rCode + "data <- data[data$scanID %in% colnames(mat),]\r\n\r\n";

    rCode = rCode + "scanAnnot <- ScanAnnotationDataFrame(data)\r\n\r\n";
    rCode = rCode + "nullmod <- fitNullMM(scanData = scanAnnot, outcome = \"" + pheno
            + "\", covars=c(\"" + ArrayUtils.toStr(covars, "\",\"")
            + "\"), covMatList = mat, family=gaussian)\r\n" + "\r\n";
    rCode = rCode + "save(nullmod, file=\"" + dir + "/" + pheno + "_nullmod.RData\")\r\n";
    try {
      PrintWriter out = new PrintWriter(new FileOutputStream(dir + "/" + pheno + "_nullmod.R"));
      out.println(rCode);
      out.close();

      Qsub.qsub(dir + "/run_" + pheno + "_nullmod.qsub",
                Rscript.getRscriptExecutable(new Logger()) + " --no-save " + dir + "/" + pheno
                                                         + "_nullmod.R",
                24000, 24.0D, 1);
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  private static void stitch(String pheno, String dir, String libpath, String gds) {
    String[] gdsHalves = gds.split("#");
    String rCode = ".libPaths(\"" + libpath + "\")\n"
                   + "source(\"https://bioconductor.org/biocLite.R\")\n" + "biocLite(\"GENESIS\")\n"
                   + "library(GENESIS)\n" + "library(GWASTools)\n" + "library(SeqVarTools)\n" + "\n"
                   + "setwd(\"" + dir + "/results/\")\n" + "for(i in 1:23) {\n"
                   + "  pattern = paste(\"" + pheno + "_chr\", i, \"-\\\\d+.csv\", sep=\"\")\n"
                   + "  files <- list.files(pattern=pattern)\n" + "  file <- read.csv(files[1])\n"
                   + "\n" + "  for (f in files[-1]) {\n" + "    df <- read.csv(f)\n"
                   + "    file <- rbind(file, df)\n" + "  }\n" + "\n"
                   + "  file <- file[order(file$snpID),]\n" + "\n" + "  gds <- seqOpen(paste(\""
                   + gdsHalves[0] + "\",i,\"" + gdsHalves[1] + "\",sep=\"\"))\n"
                   + "  seqData <- SeqVarData(gds)\n" + "\n"
                   + "  seqSetFilter(seqData, file$snpID)\n" + "\n"
                   + "  file$pos <- seqGetData(seqData, \"position\")\n"
                   + "  file$allele <- seqGetData(seqData, \"allele\")\n" + "\n"
                   + "  write.table(file, paste(\"" + pheno
                   + "_chr\",i,\".csv\",sep=\"\"), sep=\",\", row.names=F)\n"
                   + "  seqClose(seqData)\n" + "}";

    try {
      PrintWriter out = new PrintWriter(new FileOutputStream(dir + "/" + pheno + "_stitch.R"));
      out.println(rCode);
      out.close();

      Qsub.qsub(dir + "/run_" + pheno + "_stitch.qsub",
                Rscript.getRscriptExecutable(new Logger()) + " --no-save " + dir + "/" + pheno
                                                        + "_stitch.R",
                24000, 24.0D, 1);
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  public static void main(String[] args) throws ParseException {
    CLI c = new CLI("GENESIS");
    c.addArg("libpath", "Path to library folder where GENESIS is installed");
    c.addArg("pheno", "Column name of the phenotype to be analyzed");
    c.addArg("phenoFile", "Name of the phenotype file");
    c.addArg("relMatrix", "Name of the GENESIS .RData relationship matrix");
    c.addArg("covars", "A comma separated list of covariates to be included in the model");
    c.addArg("dir", "Path to the working directory");
    c.addArg("gds", "Path to the GENESIS .gds genotype matrix");

    c.parse(args);

    nullmod(c.get("phenoFile"), c.get("pheno"), c.get("covars").split(","), c.get("libpath"),
            c.get("relMatrix"), c.get("dir"));

    Vector<String> v = new Vector<String>();
    v.add("cd " + c.get("dir") + "/batches/");
    for (int i = 1; i <= 23; i++) {
      int num = batchChrs(c.get("pheno"), c.get("libpath"), c.get("gds"), i, c.get("dir"));
      if (num > 0) v.add("qsub run_" + c.get("pheno") + "_chr" + i + ".pbs");
    }
    Files.writeArray(ArrayUtils.toStringArray(v), c.get("dir") + "/scriptAll");
    Files.chmod("scriptAll");
    stitch(c.get("pheno"), c.get("dir"), c.get("libpath"), c.get("gds"));
  }
}
