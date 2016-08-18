package org.genvisis.gwas;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Date;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.IntVector;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.ext;
import org.genvisis.stats.Rscript;

public class SkatMetaPrimary {

  public static void batch(String cohort, String genos, String phenoFilename, String snpInfo,
                           int qsubMem, double qsubWalltime) {
    String phenoDir;
    String phenoRoot;
    String resultDir;
    String currentGeno;
    String currentSnpInfo;
    String rCode;
    String batchDir;
    PrintWriter out;
    Vector<String> v, consolidateVector;
    String filename, commands;
    String[][] iterations;
    boolean foundGenos;
    boolean foundSnpInfo;
    String consolidate;
    Vector<String> jobNamesWithAbsolutePaths;
    IntVector jobSizes;

    if (!new File(phenoFilename).exists()) {
      System.err.println("Error - File not found " + phenoFilename);
      return;
    }

    // create directory called dir+root+"/"; example: "c:/diffpath/pheno_F7_studyIDs/"
    phenoDir = ext.parseDirectoryOfFile(new File(phenoFilename).getAbsolutePath());
    phenoRoot = ext.rootOf(phenoFilename);
    if (!new File(phenoDir + phenoRoot + "/").exists()
        || !new File(phenoDir + phenoRoot + "/").isDirectory()) {
      new File(phenoDir + phenoRoot + "/").mkdirs();
    }

    resultDir = phenoDir + phenoRoot + "/results/";
    if (!new File(resultDir).exists() || !new File(resultDir).isDirectory()) {
      new File(resultDir).mkdirs();
    }

    batchDir = phenoDir + phenoRoot + "/batchFiles/";
    if (!new File(batchDir).exists() || !new File(batchDir).isDirectory()) {
      new File(batchDir).mkdirs();
    }

    try {
      v = new Vector<String>();
      consolidateVector = new Vector<String>();
      // generate batch files in dir+root+"/batchFiles/"; example:
      // "c:/diffpath/pheno_F7_studyIDs/batchFiles/chr1.R"
      foundGenos = false;
      foundSnpInfo = false;
      for (int i = 1; i <= 26; i++) {
        currentGeno = ext.replaceAllWith(genos, "#", i + ""); // Files.getAppropriateWriter();
        if (snpInfo.contains("#")) {
          currentSnpInfo = ext.replaceAllWith(snpInfo, "#", i + "");
        } else {
          currentSnpInfo = snpInfo;
        }
        if (new File(currentGeno).exists() && new File(currentSnpInfo).exists()) {
          foundGenos = true;
          foundSnpInfo = true;
          rCode = "library(\"skatMeta\")\n" + "setwd(\"" + resultDir + "\")\n" + "\n"
                  + (currentSnpInfo.toLowerCase()
                                   .endsWith(".rdata") ? "obj_name <- load(\"" + currentSnpInfo
                                                         + "\")\n" + "SNPInfo <- get(obj_name)\n"
                                                         + "rm(list=obj_name)\n" + "rm(obj_name)\n"
                                                       : "SNPInfo <- read.csv(\"" + currentSnpInfo
                                                         + "\", header=T, as.is=T)\n")
                  + "\n"

                  + (genos.toLowerCase().endsWith(".rdata")
                                                            ? "genoName <- load(\"" + currentGeno
                                                              + "\")\n" + "Z <- get(genoName)\n"
                                                              + "names <- colnames(Z)\n"
                                                              + "for (i in 1:ncol(Z)) {\n"
                                                              + "    names[i] <- paste(\"chr\", names[i], sep=\"\")\n"
                                                              + "}\n"
                                                            : "Z <- t(read.csv(\"" + currentGeno
                                                              + "\", header=T, as.is=T, row.names=1))\n"
                                                              + "names <- colnames(Z)\n"
                                                              + "for (i in 1:ncol(Z)) {\n"
                                                              + "    tmp <- names[i]\n"
                                                              + "    if (\"_\" == substr(tmp, start=nchar(tmp)-1, stop=nchar(tmp)-1)) {\n"
                                                              + "        names[i] = substr(tmp, start=1, stop=nchar(tmp)-2);\n"
                                                              + "    }\n" + "}\n")
                  + "colnames(Z) <- names\n" + "\n" + "pheno <- read.csv(\"" + phenoFilename
                  + "\", header=T, as.is=T, row.names=1)\n" + "xphen <- na.omit(pheno)\n"
                  + "merged <- merge(xphen, Z, by=\"row.names\")\n"
                  + "mPheno <- merged[,1:ncol(pheno)+1]\n" + "names <- colnames(pheno)\n"
                  + "if (length(names)>1) {\n" + "    formu <- paste(names[1], \"~\", names[2])\n"
                  + "    for (i in 3:length(names)) {\n"
                  + "        formu <- paste(formu, \"+\", names[i])\n" + "    }\n" + "} else {\n"
                  + "    len <- length(mPheno)\n" + "    mPheno <- c(mPheno, rep(1, len))\n"
                  + "    dim(mPheno) <- c(len, 2)\n" + "    mPheno <- as.data.frame(mPheno)\n"
                  + "    colnames(mPheno) <- c(names[1], \"dummy\")\n"
                  + "    formu <- paste(names[1], \"~ 1\")\n" + "}\n" + "\n"
                  + "offset <- 1+ncol(pheno)\n" + "mGeno <- merged[,1:ncol(Z)+offset]\n" + "\n"
                  + cohort + "_chr" + i
                  + " <- skatCohort(Z=mGeno, formula(formu), SNPInfo=SNPInfo, snpNames=\"SNP\", aggregateBy=\"SKATgene\", data=mPheno)\n"
                  // + "results <- singlesnpMeta(" + cohort + "_chr" + i + ", SNPInfo=SNPInfo,
                  // snpNames = \"Name\", cohortBetas = TRUE)\n"
                  + "results <- burdenMeta(" + cohort + "_chr" + i
                  + ", aggregateBy=\"SKATgene\", mafRange = c(0,0.05), SNPInfo=subset(SNPInfo, sc_nonsynSplice==TRUE), snpNames=\"SNP\", wts = 1)\n"
                  + "write.table(results, \"" + cohort + "_chr" + i
                  + "_beforeSave_results.csv\", sep=\",\", row.names = F)\n" + "save(" + cohort
                  + "_chr" + i + ", file=\"" + cohort + "_chr" + i
                  + ".RData\", compress=\"bzip2\")";


          // consolidate won't run if it's not added
          filename = batchDir + cohort + "_chr" + i + ".R";
          if (!Files.exists(resultDir + cohort + "_chr" + i + ".RData")) {
            out = new PrintWriter(new FileOutputStream(filename));
            out.println(rCode);
            out.close();

            // then run the R code, if it only takes a few minutes then on your machine
            // Runtime.getRuntime().exec("C:/Progra~1/R/R-3.0.1/bin/Rscript " + batchDir + "chr" + i
            // + ".R");
            v.add(filename);
          }
          consolidateVector.add(filename);

        }
      }
      if (!foundGenos) {
        System.err.println("Error - Files not found " + genos);
      }
      if (!foundSnpInfo) {
        System.err.println("Error - Files not found " + snpInfo);
      }

      iterations = Matrix.toMatrix(Array.toStringArray(v));
      System.out.println(iterations.length + "\tremaining to run for " + cohort);
      if (Files.isWindows()) {
        commands = "Rscript --no-save [%0]";
        Files.batchIt(batchDir + "run", "", 5, commands, iterations);
      } else {
        // commands = "/soft/R/3.0.1/bin/Rscript --no-save [%0]";
        commands = Rscript.getRscriptExecutable(new Logger()) + " --no-save [%0]";
        // Files.qsub("checkObject", dir, -1, commands, iterations, qsubMem, qsubWalltime);
        Files.qsub(batchDir + "run_" + cohort, batchDir, -1, commands, iterations, qsubMem,
                   qsubWalltime);
        if (iterations.length == 0) {
          new File(batchDir + "master.run_" + cohort).renameTo(new File(batchDir + "master.run_"
                                                                        + cohort + ".bak"));
        }
      }

      if (iterations.length > 0) {
        commands = "cd " + batchDir + "\n";
        if (iterations.length == 1) {
          commands += "./run_" + cohort + ".qsub\n";
        } else {
          for (int j = 0; j < iterations.length; j++) {
            commands += "./run_" + cohort + "_" + (j + 1) + ".qsub\n";
          }
        }
        Files.qsub(batchDir + "finishUpOnSB_" + cohort, commands, 60000,
                   (int) Math.ceil(iterations.length / 2.0), 16);
      } else {
        new File(batchDir + "finishUpOnSB_" + cohort).delete();
      }

      iterations = Matrix.toMatrix(Array.toStringArray(consolidateVector));
      v = new Vector<String>();
      jobNamesWithAbsolutePaths = new Vector<String>();
      jobSizes = new IntVector();
      v.add("setwd(\"" + resultDir + "\")");
      consolidate = cohort + "<- c(";
      for (int i = 0; i < iterations.length; i++) {
        v.add("load(\"" + ext.rootOf(iterations[i][0]) + ".RData\")");
        consolidate += (i == 0 ? "" : ", ") + ext.rootOf(iterations[i][0]);
        jobNamesWithAbsolutePaths.add(batchDir + "run_" + cohort + "_" + (i + 1) + ".qsub");
        jobSizes.add(Integer.MAX_VALUE
                     - (int) new File(ext.replaceAllWith(genos, "#", (i + 1) + "")).length());
      }
      consolidate += ")";
      v.add(consolidate);
      v.add("class(" + cohort + ") <- \"skatCohort\"");
      v.add("save(" + cohort + ", file=\"" + cohort + ".RData\", compress=\"bzip2\")");
      Files.writeList(Array.toStringArray(v), batchDir + "mergeRdataFiles.R");
      commands = Rscript.getRscriptExecutable(new Logger()) + " --no-save " + batchDir
                 + "mergeRdataFiles.R";
      Files.qsub(batchDir + "run_mergeRdataFiles_" + cohort, commands, qsubMem * 4, qsubWalltime,
                 1);
      Files.qsubMultiple(jobNamesWithAbsolutePaths, jobSizes, batchDir,
                         batchDir + "chunk_" + cohort, 8, true, null, -1, qsubMem, qsubWalltime);
      Files.qsubMultiple(jobNamesWithAbsolutePaths, jobSizes, batchDir,
                         batchDir + "chunkSB256_" + cohort, 16, true, "sb256", -1, qsubMem,
                         qsubWalltime);
      Files.qsubMultiple(jobNamesWithAbsolutePaths, jobSizes, batchDir,
                         batchDir + "chunkSB_" + cohort, 16, true, "sb", -1, qsubMem, qsubWalltime);

    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  public static void batchMany(String cohort, String genos, String phenosCommaDelimited,
                               String racesCommaDelimited, String snpInfo, int qsubMem,
                               double qsubWalltime) {
    String[] phenos, races;
    Vector<String> v;

    phenos = phenosCommaDelimited.split(",");
    races = racesCommaDelimited.split(",");

    v = new Vector<String>();
    for (String pheno : phenos) {
      for (String race : races) {
        try {
          batch(cohort + "_" + race + "_" + pheno, ext.replaceAllWith(genos, "[%race]", race),
                ext.pwd() + cohort + "_" + race + "_" + pheno + ".csv", snpInfo, qsubMem,
                qsubWalltime);
        } catch (Exception e) {
          System.err.println("Error - failed to script up " + pheno + "/" + race);
        }
        v.add("cd " + cohort + "_" + race + "_" + pheno + "/batchFiles/");
        v.add("./master.run_" + cohort + "_" + race + "_" + pheno);
        v.add("cd ../../");
        v.add("");
      }
    }
    Files.writeList(Array.toStringArray(v), "scriptAll");
    Files.chmod("scriptAll");

    v = new Vector<String>();
    for (String pheno : phenos) {
      for (String race : races) {
        v.add("cd " + cohort + "_" + race + "_" + pheno + "/batchFiles/");
        v.add("./master.chunkSB_" + cohort + "_" + race + "_" + pheno);
        v.add("cd ../../");
        v.add("");
      }
    }
    Files.writeList(Array.toStringArray(v), "scriptAllItasca");
    Files.chmod("scriptAllItasca");

    v = new Vector<String>();
    for (String pheno : phenos) {
      for (String race : races) {
        if (Files.exists(cohort + "_" + race + "_" + pheno + "/batchFiles/finishUpOnSB_" + cohort
                         + "_" + race + "_" + pheno)) {
          v.add("cd " + cohort + "_" + race + "_" + pheno + "/batchFiles/");
          v.add("qsub -q sb finishUpOnSB_" + cohort + "_" + race + "_" + pheno);
          v.add("cd ../../");
          v.add("");
        }
      }
    }
    Files.writeList(Array.toStringArray(v), "finishUpOnSB");
    Files.chmod("finishUpOnSB");

    v = new Vector<String>();
    for (String pheno : phenos) {
      for (String race : races) {
        v.add("cd " + cohort + "_" + race + "_" + pheno + "/batchFiles/");
        v.add(Rscript.getRscriptExecutable(new Logger()) + " --no-save mergeRdataFiles.R");
        v.add("mv ../results/" + cohort + "_" + race + "_" + pheno + ".RData ../../");
        v.add("cd ../../");
        v.add("");
      }
    }
    Files.writeList(Array.toStringArray(v), "mergeAll");
    Files.chmod("mergeAll");

    v = new Vector<String>();
    for (String pheno : phenos) {
      for (String race : races) {
        v.add("cd " + cohort + "_" + race + "_" + pheno + "/results/");
        v.add("tar -zcvf ../../" + cohort + "_" + race + "_" + pheno + "_"
              + ext.getDate(new Date(), "") + ".tar.gz " + cohort + "_" + race + "_" + pheno
              + "*.RData");
        v.add("cd ../../");
        v.add("");
      }
    }
    Files.writeList(Array.toStringArray(v), "packageUpAll");
    Files.chmod("packageUpAll");
  }

  public static void additionalModels(String cohort, String phenosCommaDelimited, String snpInfo,
                                      int qsubMem, double qsubWalltime) {
    String[] phenos;
    Vector<String> v;

    phenos = phenosCommaDelimited.split(",");

    v = new Vector<String>();
    for (String pheno : phenos) {
      try {
        batchAdditionals(pheno, cohort, snpInfo, qsubMem, qsubWalltime);
      } catch (Exception e) {
        System.err.println("Error - failed to script up " + pheno);
      }
      v.add("cd " + pheno + "/batchFiles/");
      v.add("./master.run_additionals");
      v.add("cd ../../");
      v.add("");
    }
    Files.writeList(Array.toStringArray(v), "addAll");
    Files.chmod("addAll");
  }

  private static void batchAdditionals(String phenoDirectory, String cohort, String snpInfo,
                                       int qsubMem, double qsubWalltime) {
    String phenoDir;
    // String phenoRoot;
    String resultDir;
    // String currentGeno;
    String currentSnpInfo;
    String rCode;
    String batchDir;
    PrintWriter out;
    Vector<String> v;
    String filename, commands;
    String[][] iterations;
    // boolean foundGenos;
    // boolean foundSnpInfo;
    // String consolidate;
    // Vector<String> jobNamesWithAbsolutePaths;
    // IntVector jobSizes;

    if (!new File(phenoDirectory).exists()) {
      System.err.println("Error - directory not found: " + phenoDirectory);
      return;
    }

    // create directory called dir+root+"/"; example: "c:/diffpath/pheno_F7_studyIDs/"
    phenoDir =
        ext.verifyDirFormat(ext.parseDirectoryOfFile(new File(phenoDirectory).getAbsolutePath()));
    resultDir = phenoDir + "/results/";
    batchDir = phenoDir + "/batchFiles/";

    try {
      v = new Vector<String>();
      // generate batch files in dir+root+"/batchFiles/"; example:
      // "c:/diffpath/pheno_F7_studyIDs/batchFiles/chr1.R"
      // foundGenos = false;
      // foundSnpInfo = false;
      for (int i = 1; i <= 26; i++) {
        if (snpInfo.contains("#")) {
          currentSnpInfo = ext.replaceAllWith(snpInfo, "#", i + "");
        } else {
          currentSnpInfo = snpInfo;
        }
        if (new File(resultDir + cohort + "_chr" + i + ".RData").exists()) {
          // foundGenos = true;
          // foundSnpInfo = true;
          rCode = "library(\"skatMeta\")\n" + "setwd(\"" + resultDir + "\")\n" + "\n"
                  + (currentSnpInfo.toLowerCase()
                                   .endsWith(".rdata") ? "obj_name <- load(\"" + currentSnpInfo
                                                         + "\")\n" + "SNPInfo <- get(obj_name)\n"
                                                         + "rm(list=obj_name)\n" + "rm(obj_name)\n"
                                                       : "SNPInfo <- read.csv(\"" + currentSnpInfo
                                                         + "\", header=T, as.is=T)\n")
                  + "\n"

                  + "cohortName <- load(\"" + resultDir + cohort + "_chr" + i + ".RData" + "\")\n"
                  + "cohort <- get(cohortName)\n"

                  // + "results <- singlesnpMeta(cohort, SNPInfo=SNPInfo, snpNames = \"Name\",
                  // cohortBetas = TRUE)\n"
                  // + "write.table(results, \"" + cohort + "_chr" + i + "_SingleSnp_results.csv\",
                  // sep=\",\", row.names = F)\n"
                  // + "results <- burdenMeta(" + cohort + "_chr" + i + ", aggregateBy=\"SKATgene\",
                  // mafRange = c(0,0.05), SNPInfo=subset(SNPInfo, sc_nonsynSplice==TRUE),
                  // snpNames=\"SNP\", wts = 1)\n"
                  // + "write.table(results, \"" + cohort + "_chr" + i + "burden_results.csv\",
                  // sep=\",\", row.names = F)\n"
                  // + "results <- singlesnpMeta(cohort, SNPInfo=SNPInfo, snpNames = \"Name\",
                  // cohortBetas = TRUE)\n"
                  // + "write.table(results, \"" + cohort + "_chr" + i + "_SingleSnp_results.csv\",
                  // sep=\",\", row.names = F)\n"

                  + "results <- singlesnpMeta(cohort, SNPInfo=SNPInfo_ExomeChipV5, snpNames = \"Name\", aggregateBy=\"SKATgene\", cohortBetas = TRUE)\n"
                  + "write.table( results, \"SingleSNP_chr\"+(i)+\".csv\", sep=\",\", row.names = F)\n"

                  + "results <- burdenMeta(cohort, SNPInfo=subset(SNPInfo_ExomeChipV5, sc_nonsynSplice==TRUE), snpNames = \"Name\", aggregateBy=\"SKATgene\", mafRange = c(0,0.01), wts = 1)\n"
                  + "write.table( results, \"T1Count_chr\"+(i)+\".csv\", sep=\",\", row.names = F)\n"

                  + "results <- burdenMeta(cohort, SNPInfo=subset(SNPInfo_ExomeChipV5, sc_nonsynSplice==TRUE), snpNames = \"Name\", aggregateBy=\"SKATgene\", mafRange = c(0,0.05), wts = 1)\n"
                  + "write.table( results, \"T5Count_chr\"+(i)+\".csv\", sep=\",\", row.names = F)\n"

                  + "results <- skatMeta(cohort, SNPInfo=subset(SNPInfo_ExomeChipV5, sc_nonsynSplice==TRUE), snpNames = \"Name\", aggregateBy=\"SKATgene\", mafRange = c(0,0.01), wts = 1)\n"
                  + "write.table( results, \"SKAT_T1_chr\"+(i)+\".csv\", sep=\",\", row.names = F)\n"

                  + "results <- skatMeta(cohort, SNPInfo=subset(SNPInfo_ExomeChipV5, sc_nonsynSplice==TRUE), snpNames = \"Name\", aggregateBy=\"SKATgene\", mafRange = c(0,0.05), wts = 1)\n"
                  + "write.table( results, \"SKAT_T5_chr\"+(i)+\".csv\", sep=\",\", row.names = F)\n"


          ;


          filename = batchDir + "additionals_chr" + i + ".R";
          out = new PrintWriter(new FileOutputStream(filename));
          out.println(rCode);
          out.close();

          v.add(filename);

        }
      }

      iterations = Matrix.toMatrix(Array.toStringArray(v));
      System.out.println(iterations.length + "\tremaining to run for " + cohort);
      if (Files.isWindows()) {
        commands = "Rscript --no-save [%0]";
        Files.batchIt(batchDir + "run", "", 5, commands, iterations);
      } else {
        commands = Rscript.getRscriptExecutable(new Logger()) + " --no-save [%0]";
        Files.qsub(batchDir + "run_additionals", batchDir, -1, commands, iterations, qsubMem,
                   qsubWalltime);
        if (iterations.length == 0) {
          new File(batchDir
                   + "master.run_additionals").renameTo(new File(batchDir
                                                                 + "master.run_additionals.bak"));
        }
      }

      // v = new Vector<String>();
      // jobNamesWithAbsolutePaths = new Vector<String>();
      // jobSizes = new IntVector();
      // v.add("setwd(\"" + resultDir + "\")");
      // consolidate = cohort + "<- c(";
      // for (int i = 0; i < iterations.length; i++) {
      // v.add("load(\"" + ext.rootOf(iterations[i][0])+ ".RData\")");
      // consolidate += (i==0?"":", ")+ext.rootOf(iterations[i][0]);
      // jobNamesWithAbsolutePaths.add(batchDir + "run_" + cohort+"_"+(i+1)+".qsub");
      // jobSizes.add((int)new File(ext.replaceAllWith(genos, "#", (i+1)+"")).length());
      // }
      // consolidate += ")";
      // v.add(consolidate);
      // v.add("class("+cohort + ") <- \"skatCohort\"");
      // v.add("save(" + cohort + ", file=\"" + cohort + ".RData\", compress=\"bzip2\")");
      // Files.writeList(Array.toStringArray(v), batchDir + "mergeRdataFiles.R");
      // commands = Rscript.getRscriptExecutable(new Logger())+" --no-save "+batchDir +
      // "mergeRdataFiles.R";
      // Files.qsub(batchDir + "run_mergeRdataFiles_" + cohort, commands, qsubMem*4, qsubWalltime);
      // Files.qsubMultiple(jobNamesWithAbsolutePaths, jobSizes, batchDir, batchDir+"chunk_"+cohort,
      // 8, true, 2000, 1);

    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  public static void renameMany(String cohort, String phenosCommaDelimited,
                                String racesCommaDelimited, String insert) {
    String[] phenos, races;

    phenos = phenosCommaDelimited.split(",");
    races = racesCommaDelimited.split(",");

    for (String pheno : phenos) {
      for (String race : races) {
        new File(cohort + "_" + race + "_"
                 + pheno + ".RData")
                                    .renameTo(new File(cohort + insert + "_" + race + "_" + pheno
                                                       + ".RData"));
      }
    }
  }

  public static void main(String[] args) {
    String cohort;
    String genos;
    String pheno;
    String snpInfo;
    String phenos = null;
    String races = "EA,AA";
    int qsubMem;
    double qsubWalltime;
    boolean additionals = false;
    String rename = null;

    cohort = "ARIC";
    genos = "D:/SkatMeta/genotypes_blacks_AA/AA_ARIC_noJHS_chr#t.csv";
    pheno = "D:/SkatMeta/results_hemostasis/pheno_F7_studyIDs.csv";
    snpInfo = "N:/statgen/skatMeta/fullExample/SNPInfo_HumanExome-12v1_rev5_justGeneSpliced.csv";
    qsubMem = 15000;
    qsubWalltime = 2;

    String usage = "\n" + "gwas.SkatMetaPrimary requires 4 arguments\n"
                   + "   (1) cohort name (i.e. cohort=" + cohort + " (not the default))\n"
                   + "   (2) genotype file name (i.e. geno=" + genos + " (not the default))\n"
                   + "   (3) phenotype file name (i.e. pheno=" + pheno + " (not the default))\n"
                   + "   (4) snpInfo file name (i.e. snpInfo=" + snpInfo + " (not the default))\n"
                   + "   (5) qsub memory in megabytes (i.e. qsubmem=" + qsubMem + " (default))\n"
                   + "   (6) qsub walltime in hours (i.e. qsubwalltime=" + qsubWalltime
                   + " (default))\n"
                   + "   (7) (optional) phenotypes to run (i.e. phenos=CRP,ICAM1,TNFA (not the default))\n"
                   + "   (8) (optional) races to run (i.e. races=" + races + " (default))\n"
                   + " ( the final optional parameters require a fixed phenotype format of [cohort name]_[race]_[phenotype]_.csv )\n"
                   + " ( as well as a genotype file with the format of regular_filename_[%race]_chr#.csv )\n"
                   + "   (9) (optional) run additional models for Rdata files (i.e. -additionals  (not the default))\n"
                   + " OR\n"
                   + "   (9) (optional) rename Rdata files with this insert (i.e. rename=42PCs (not the default))\n"
                   + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("cohort=")) {
        cohort = arg.split("=")[1];
      } else if (arg.startsWith("geno=")) {
        genos = arg.split("=")[1];
      } else if (arg.startsWith("pheno=")) {
        pheno = arg.split("=")[1];
      } else if (arg.startsWith("snpInfo=")) {
        snpInfo = arg.split("=")[1];
      } else if (arg.startsWith("qsubmem=")) {
        qsubMem = Integer.parseInt(arg.split("=")[1]);
      } else if (arg.startsWith("qsubwalltime=")) {
        qsubWalltime = Double.parseDouble(arg.split("=")[1]);
      } else if (arg.startsWith("phenos=")) {
        phenos = arg.split("=")[1];
      } else if (arg.startsWith("races=")) {
        races = arg.split("=")[1];
      } else if (arg.startsWith("-additionals")) {
        additionals = true;
      } else if (arg.startsWith("rename=")) {
        rename = arg.split("=")[1];
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }

    if (phenos != null) {
      if (rename != null) {
        renameMany(cohort, phenos, races, rename);
      } else if (additionals) {
        additionalModels(cohort, phenos, snpInfo, qsubMem, qsubWalltime);
      } else {
        batchMany(cohort, genos, phenos, races, snpInfo, qsubMem, qsubWalltime);
      }
    } else {
      batch(cohort, genos, pheno, snpInfo, qsubMem, qsubWalltime);
    }
  }

}
