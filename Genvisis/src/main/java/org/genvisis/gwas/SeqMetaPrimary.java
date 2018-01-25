package org.genvisis.gwas;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Vector;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.IntVector;
import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;
import org.genvisis.common.ext;
import org.genvisis.filesys.MetaAnalysisParams;
import org.genvisis.qsub.Qsub;
import org.genvisis.stats.Rscript;

public class SeqMetaPrimary {

  public static void batch(String cohort, String genos, String phenoFilename, String snpInfo,
                           int qsubMem, double qsubWalltime, String queue, boolean usePrep2,
                           String saveName, String conditionals, String resultDir) {
    String phenoDir;
    String phenoRoot;
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
    String name;
    String chrName;

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

    if (resultDir == null) {
      resultDir = phenoDir + phenoRoot + "/results/";
    } else {
      resultDir = new File(resultDir).getAbsolutePath() + "/";
    }
    if (!new File(resultDir).exists() || !new File(resultDir).isDirectory()) {
      new File(resultDir).mkdirs();
    }

    batchDir = phenoDir + phenoRoot + "/batchFiles/";
    if (!new File(batchDir).exists() || !new File(batchDir).isDirectory()) {
      new File(batchDir).mkdirs();
    }

    conditionals = conditionals == null ? null : new File(conditionals).getAbsolutePath();

    try {
      v = new Vector<>();
      consolidateVector = new Vector<>();
      // generate batch files in dir+root+"/batchFiles/"; example:
      // "c:/diffpath/pheno_F7_studyIDs/batchFiles/chr1.R"
      foundGenos = false;
      foundSnpInfo = false;
      for (int i = 1; i <= 26; i++) {
        if (i == 23) {
          chrName = "X";
        } else if (i == 24) {
          chrName = "Y";
        } else if (i == 25) {
          chrName = "XY";
        } else if (i == 26) {
          chrName = "MT";
        } else {
          chrName = "" + i;
        }
        name = (saveName == null ? cohort + "_chr" : saveName) + chrName;
        currentGeno = ext.replaceAllWith(genos, "#", i + ""); // Files.getAppropriateWriter();
        if (snpInfo.contains("#")) {
          currentSnpInfo = ext.replaceAllWith(snpInfo, "#", i + "");
        } else {
          currentSnpInfo = snpInfo;
        }

        if (new File(currentGeno).exists()) {
          foundGenos = true;
          if (new File(currentSnpInfo).exists()) {
            foundSnpInfo = true;
          } else {
            System.err.println("Error - Files not found " + currentSnpInfo
                               + "; this chromosome will be skipped");
            foundSnpInfo = false;
          }
        } else {
          System.err.println("Error - Files not found " + currentGeno
                             + "; this chromosome will be skipped");
          foundGenos = false;
        }

        if (foundGenos && foundSnpInfo) {
          rCode = rCode(currentSnpInfo, currentGeno, phenoFilename, conditionals, resultDir, i + "",
                        name, usePrep2);

          // consolidate won't run if it's not added
          filename = batchDir + name + ".R";
          if (!Files.exists(resultDir + name + ".RData")) {
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

      iterations = Matrix.toMatrix(ArrayUtils.toStringArray(v));
      System.out.println(iterations.length + "\tremaining to run for " + cohort);

      jobNamesWithAbsolutePaths = new Vector<String>();
      jobSizes = new IntVector();

      for (String f : v) {
        jobNamesWithAbsolutePaths.add(Rscript.getRscriptExecutable(new Logger()) + " --no-save "
                                      + f);
        jobSizes.add(qsubMem);
      }

      Qsub.qsubExecutor(batchDir, jobNamesWithAbsolutePaths, jobSizes, batchDir + "/batchChrs", 24,
                        62000, qsubWalltime * 2);

      if (iterations.length > 0) {
        commands = "cd " + batchDir + "\n";
        if (iterations.length == 1) {
          commands += "./run_" + cohort + ".qsub\n";
        } else {
          for (int j = 0; j < iterations.length; j++) {
            commands += "./run_" + cohort + "_" + (j + 1) + ".qsub\n";
          }
        }
        Qsub.qsub(batchDir + "finishWithHigherMem_" + cohort, commands, 60000,
                  (int) Math.ceil(iterations.length / 2.0), 1);
      } else {
        new File(batchDir + "finishWithHigherMem_" + cohort).delete();
      }

      name = saveName == null ? cohort : saveName;
      iterations = Matrix.toMatrix(ArrayUtils.toStringArray(consolidateVector));
      v = new Vector<>();
      jobNamesWithAbsolutePaths = new Vector<>();
      jobSizes = new IntVector();
      v.add("setwd(\"" + resultDir + "\")");
      consolidate = name + "<- c(";
      for (int i = 0; i < iterations.length; i++) {
        v.add("load(\"" + ext.rootOf(iterations[i][0]) + ".RData\")");
        consolidate += (i == 0 ? "" : ", ") + ext.rootOf(iterations[i][0]);
        jobNamesWithAbsolutePaths.add(batchDir + "run_" + name + "_" + (i + 1) + ".qsub");
        jobSizes.add(Integer.MAX_VALUE
                     - (int) new File(ext.replaceAllWith(genos, "#", (i + 1) + "")).length());
      }
      consolidate += ")";
      v.add(consolidate);
      v.add("class(" + name + ") <- \"skatCohort\"");
      v.add("save(" + name + ", file=\"" + name + ".RData\", compress=\"bzip2\")");
      Files.writeArray(ArrayUtils.toStringArray(v), batchDir + "mergeRdataFiles.R");
      commands = Rscript.getRscriptExecutable(new Logger()) + " --no-save " + batchDir
                 + "mergeRdataFiles.R";
      Qsub.qsub(batchDir + "run_mergeRdataFiles_" + name, commands, qsubMem * 4, qsubWalltime, 1);
      // Files.qsubMultiple( jobNamesWithAbsolutePaths, jobSizes, batchDir,
      // batchDir + "chunk_" + cohort, 8, true, null, -1, qsubMem, qsubWalltime);
      // Files.qsubMultiple( jobNamesWithAbsolutePaths, jobSizes, batchDir,
      // batchDir + "chunkSB256_" + cohort, 16, true, "sb256", -1, qsubMem,
      // qsubWalltime);
      // Files.qsubMultiple( jobNamesWithAbsolutePaths, jobSizes, batchDir,
      // batchDir + "chunkSB_" + cohort, 16, true, "sb", -1, qsubMem,
      // qsubWalltime);

    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  public static void batchMany(String cohort, String genos, String phenosCommaDelimited,
                               String racesCommaDelimited, String snpInfo, int qsubMem,
                               double qsubWalltime, String queue, boolean usePrep2, String saveName,
                               String conditionals, String resultDir) {
    String[] phenos, races;
    Vector<String> v;

    phenos = phenosCommaDelimited.split(",");
    races = racesCommaDelimited.split(",");

    v = new Vector<>();
    Vector<String> bash = new Vector<String>();
    bash.add("#!/bin/bash");
    boolean firstID = true;

    for (String pheno : phenos) {
      for (String race : races) {
        try {
          String sn = saveName;
          String r = resultDir;
          if (saveName != null) {
            SimpleDateFormat s = new SimpleDateFormat("ddMMMyyyy");
            Date now = new Date();
            String d = s.format(now).toUpperCase();
            sn = ext.replaceAllWith(saveName, "[%race]", (race.equals("AA") ? "AFA" : "EUR"));
            sn = ext.replaceAllWith(sn, "[%date]", d);
          }
          if (resultDir != null) {
            r = ext.replaceAllWith(r, "[%cohort]", cohort);
            r = ext.replaceAllWith(r, "[%pheno]", pheno);
            r = ext.replaceAllWith(r, "[%race]", race);
          }
          batch(cohort + "_" + race + "_" + pheno, ext.replaceAllWith(genos, "[%race]", race),
                ext.pwd() + cohort + "_" + race + "_" + pheno + ".csv", snpInfo, qsubMem,
                qsubWalltime, queue, usePrep2, sn, conditionals, r);
        } catch (Exception e) {
          System.err.println("Error - failed to script up " + pheno + "/" + race);
          e.printStackTrace();
        }
        v.add("cd " + cohort + "_" + race + "_" + pheno + "/batchFiles/");
        v.add("qsub batchChrs.pbs");
        v.add("cd ../../");
        v.add("");

        if (firstID) {
          bash.add("IDS=\"\"");
          firstID = false;
        } else {
          bash.add("IDS += \":\"");
        }
        bash.add("IDS += `qsub " + cohort + "_" + race + "_" + pheno + "/batchFiles/`");
      }
    }
    bash.add("echo $IDS");

    Files.writeArray(ArrayUtils.toStringArray(v), "scriptAll");
    Files.chmod("scriptAll");

    Files.writeArray(ArrayUtils.toStringArray(bash), "bashScriptAll");
    Files.chmod("bashScriptAll");

    // v = new Vector<String>();
    // for (String pheno : phenos) {
    // for (String race : races) {
    // v.add("cd " + cohort + "_" + race + "_" + pheno + "/batchFiles/");
    // v.add("./master.chunkSB_" + cohort + "_" + race + "_" + pheno);
    // v.add("cd ../../");
    // v.add("");
    // }
    // }
    // Files.writeArray(Array.toStringArray(v), "scriptAllItasca");
    // Files.chmod("scriptAllItasca");
    //
    v = new Vector<>();
    for (String pheno : phenos) {
      for (String race : races) {
        if (Files.exists(cohort + "_" + race + "_" + pheno + "/batchFiles/finishWithHigherMem_"
                         + cohort + "_" + race + "_" + pheno)) {
          v.add("cd " + cohort + "_" + race + "_" + pheno + "/batchFiles/");
          v.add("qsub finishWithHigherMem_" + cohort + "_" + race + "_" + pheno);
          v.add("cd ../../");
          v.add("");
        }
      }
    }
    Files.writeArray(ArrayUtils.toStringArray(v), "finishWithHigherMem");
    Files.chmod("finishWithHigherMem");

    v = new Vector<>();
    for (String pheno : phenos) {
      for (String race : races) {
        v.add("cd " + cohort + "_" + race + "_" + pheno + "/batchFiles/");
        v.add(Rscript.getRscriptExecutable(new Logger()) + " --no-save mergeRdataFiles.R");
        v.add("mv ../results/" + cohort + "_" + race + "_" + pheno + ".RData ../../");
        v.add("cd ../../");
        v.add("");
      }
    }
    Files.writeArray(ArrayUtils.toStringArray(v), "mergeAll");
    Files.chmod("mergeAll");

    v = new Vector<>();
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
    Files.writeArray(ArrayUtils.toStringArray(v), "packageUpAll");
    Files.chmod("packageUpAll");
  }

  public static String run(MetaAnalysisParams maps, String pheno, String[] chrs, String condFile,
                           String resultDir, boolean usePrep2, Logger log) {
    String rCode, currentGeno, currentSnpInfo;
    boolean runByChr = maps.runningByChr();
    String snpInfo = maps.getSnpInfoChrsDir()
                     + (runByChr ? "snpInfo_chr#.RData" : maps.getSnpInfoFilename());
    String genos = maps.getGenos();
    String[][] races = maps.getRacesWithFilenameAliases();
    String[] studies = maps.getStudies();
    String primaryDir = maps.getPrimaryDir();
    String rFiles = "";
    String saveName = "";

    if (studies == null || studies.length == 0) studies = new String[] {""};

    for (String study : studies) {
      for (String[] race : races) {
        String name;
        String[] phenos;
        int i = 0;
        do {
          name = (study.equals("") ? "" : study) + (pheno == null ? race[i] : race[i] + pheno);
          phenos = Files.list(primaryDir, name, ".csv", false, false);
          i++;
        } while (phenos.length == 0 && i < race.length);

        for (String phenoFilename : phenos) {
          for (String c : chrs) {
            currentGeno = ext.replaceAllWith(ext.replaceAllWith(genos, "#", c), "[%race]",
                                             race[0].equals("Whites") ? "EA" : "AA");
            currentSnpInfo = ext.replaceAllWith(snpInfo, "#", c);
            saveName = ext.rootOf(phenoFilename) + "_chr" + c;

            rCode = rCode(currentSnpInfo, currentGeno, phenoFilename, condFile, resultDir, c,
                          saveName, true);

            String filename = resultDir + saveName + ".R";
            try {
              PrintWriter out = Files.getAppropriateWriter(filename);

              out.write(rCode);
              rFiles += Rscript.getRscriptExecutable(log) + " --no-save " + filename + "\n";
              out.close();
            } catch (Exception e) {
              log.reportError("Problem writing to " + filename);
            }
          }
        }

      }
    }

    return rFiles;
  }

  private static String rCode(String snpInfo, String genos, String phenoFilename, String condFile,
                              String resultDir, String chr, String saveName, boolean usePrep2) {
    phenoFilename = new File(phenoFilename).getAbsolutePath();
    resultDir = new File(resultDir).getAbsolutePath() + "/";

    String rCode = "library(\"seqMeta\")\n" + "library(\"methods\")\n" + "setwd(\"" + resultDir
                   + "\")\n" + "\n" + (snpInfo.toLowerCase()
                                              .endsWith(".rdata") ? "obj_name <- load(\"" + snpInfo + "\")\n" + "SNPInfo <- get(obj_name)\n" + "rm(list=obj_name)\n" + "rm(obj_name)\n" : "SNPInfo <- read.csv(\"" + snpInfo + "\", header=T, as.is=T)\n")
                   + "\n"
                   + (genos.toLowerCase().endsWith(".rdata")
                      || genos.toLowerCase().endsWith(".rda")
                                                              ? "genoName <- load(\"" + genos
                                                                + "\")\n" + "Z <- get(genoName)\n"
                                                                + "if (any(grepl(\":\", rownames(Z)))){\n"
                                                                + "    Z <- t(Z)\n" + "    }\n"
                                                                + "percent_miss <- mean(colnames(Z) %in% SNPInfo[,\"SNP\"])\n"
                                                                + "if (percent_miss == 0) {\n"
                                                                + "    names <- colnames(Z)\n"
                                                                + "    for (i in 1:ncol(Z)) {\n"
                                                                + "        names[i] <- paste(\"chr\", names[i], sep=\"\")\n"
                                                                + "    }\n"
                                                                + "    colnames(Z) <- names\n"
                                                                + "}\n"
                                                              : "Z <- t(read.csv(\"" + genos
                                                                + "\", header=T, as.is=T, row.names=1))\n"
                                                                + "names <- colnames(Z)\n"
                                                                + "for (i in 1:ncol(Z)) {\n"
                                                                + "    tmp <- names[i]\n"
                                                                + "    if (\"_\" == substr(tmp, start=nchar(tmp)-1, stop=nchar(tmp)-1)) {\n"
                                                                + "        names[i] = substr(tmp, start=1, stop=nchar(tmp)-2);\n"
                                                                + "    }\n" + "}\n"
                                                                + "colnames(Z) <- names\n")
                   + "\n" + "pheno <- read.csv(\"" + phenoFilename
                   + "\", header=T, as.is=T, row.names=1, check.names=FALSE)\n"
                   + "xphen <- na.omit(pheno)\n" + "b <- nrow(unique(pheno[1])) \n"
                   + "if (b > 2) { \n" + "  family <- \"gaussian\" \n" + "} else {\n"
                   + "  family <- \"binomial\"\n" + "}\n"
                   + "merged <- merge(xphen, Z, by=\"row.names\")\n"
                   + "mPheno <- merged[,1:ncol(pheno)+1]\n" + "offset <- 1+ncol(pheno)\n"
                   + "mGeno <- merged[,1:ncol(Z)+offset]\n" + "\n" + "\n"
                   + "chrs <- c(\"Chr\", \"chr\", \"CHROM\")\n"
                   + "positions <- c(\"POS\", \"pos\", \"Pos\", \"start\")\n"
                   + "chrindex <- which(colnames(SNPInfo) %in% chrs, arr.ind=T)\n"
                   + "posindex <- which(colnames(SNPInfo) %in% positions, arr.ind=T)\n"
                   + "colnames(SNPInfo)[chrindex] <- \"chr\"\n"
                   + "colnames(SNPInfo)[posindex] <- \"pos\"\n"
                   + "conditions <- SNPInfo[0,c(\"SNP\", \"SKATgene\", \"chr\")]\n"
                   + conditionals(condFile, Integer.parseInt(chr)) + "\n"
                   + "names <- colnames(mPheno)\n"
                   + "coxy <- sum(names %in% c(\"time\", \"status\"))\n"
                   + "if (length(names)>1) {\n" + "    if (coxy == 2) {\n"
                   + "        formu <- paste(\"Surv(time,status)\", \"~\")\n"
                   + "        if (length(names)>2) {\n"
                   + "            formu <- paste(formu, names[3])\n"
                   + "            for (i in 4:length(names)) {\n"
                   + "                formu <- paste(formu, \"+\", names[i])\n" + "            }\n"
                   + "        } else {\n" + "                formu <- paste(formu, \"1\")\n"
                   + "        }\n" + "    } else {\n"
                   + "        formu <- paste(names[1], \"~\", names[2])\n"
                   + "        for (i in 3:length(names)) {\n"
                   + "            formu <- paste(formu, \"+\", names[i])\n" + "        }\n"
                   + "    }\n" + "} else {\n" + "    len <- length(mPheno)\n"
                   + "    mPheno <- c(mPheno, rep(1, len))\n" + "    dim(mPheno) <- c(len, 2)\n"
                   + "    mPheno <- as.data.frame(mPheno)\n"
                   + "    colnames(mPheno) <- c(names[1], \"dummy\")\n"
                   + "    formu <- paste(names[1], \"~ 1\")\n" + "}\n" + "\n" + "formu \n"
                   + "if (nrow(conditions) > 0 & any(colnames(mGeno) %in% conditions$SNP)) {\n"
                   + "  if(family == \"binomial\") { f<-binomial() } else { f<-gaussian() }\n"
                   + "  message(\"conditions detected; using prepCondScores\")\n" + "  " + saveName
                   + "<- prepCondScores(Z=mGeno, formula(formu), SNPInfo=SNPInfo, snpNames=\"SNP\", "
                   + "aggregateBy=\"SKATgene\", data=mPheno, adjustments=conditions, family=f)\n"
                   + "} else if (coxy == 2) {\n"
                   + "    message(\"time to event data detected; using a cox model\")\n" + "    "
                   + saveName
                   + " <- prepCox(Z=mGeno, formula(formu), SNPInfo=SNPInfo, snpNames=\"SNP\", aggregateBy=\"SKATgene\", data=mPheno)\n"
                   + "} else {\n" + "    message(\"using "
                   + (usePrep2 ? "prepScores2" : "traditional prepScores method") + "\")\n    "
                   + saveName + " <- prepScores" + (usePrep2 ? "2" : "")
                   + "(Z=mGeno, formula(formu), SNPInfo=SNPInfo, snpNames=\"SNP\", aggregateBy=\"SKATgene\", data=mPheno, family=family"
                   + (chr.equals("23") ? ", male=mPheno$SEX" : "") + ")\n" + "}\n"
                   + "results <- burdenMeta(" + saveName
                   + ", aggregateBy=\"SKATgene\", mafRange = c(0,0.05), SNPInfo=SNPInfo, snpNames=\"SNP\", wts = 1)\n"
                   + "write.table(results, \"" + resultDir + saveName
                   + "_beforeSave_results.csv\", sep=\",\", row.names = F)\n" + "save(" + saveName
                   + ", file=\"" + resultDir + saveName + ".RData\", compress=\"bzip2\")";
    return rCode;
  }

  public static void additionalModels(String cohort, String phenosCommaDelimited, String snpInfo,
                                      int qsubMem, double qsubWalltime) {
    String[] phenos;
    Vector<String> v;

    phenos = phenosCommaDelimited.split(",");

    v = new Vector<>();
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
    Files.writeArray(ArrayUtils.toStringArray(v), "addAll");
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
    phenoDir = ext.verifyDirFormat(ext.parseDirectoryOfFile(new File(phenoDirectory).getAbsolutePath()));
    resultDir = phenoDir + "/results/";
    batchDir = phenoDir + "/batchFiles/";

    try {
      v = new Vector<>();
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
          rCode = "library(\"seqMeta\")\n" + "library(\"methods\")\n" + "setwd(\"" + resultDir
                  + "\")\n" + "\n"
                  + (currentSnpInfo.toLowerCase().endsWith(".rdata")
                                                                     ? "obj_name <- load(\""
                                                                       + currentSnpInfo + "\")\n"
                                                                       + "SNPInfo <- get(obj_name)\n"
                                                                       + "rm(list=obj_name)\n"
                                                                       + "rm(obj_name)\n"
                                                                     : "SNPInfo <- read.csv(\""
                                                                       + currentSnpInfo
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

      iterations = Matrix.toMatrix(ArrayUtils.toStringArray(v));
      System.out.println(iterations.length + "\tremaining to run for " + cohort);
      if (Files.isWindows()) {
        commands = "Rscript --no-save [%0]";
        Files.batchIt(batchDir + "run", "", 5, commands, iterations);
      } else {
        commands = Rscript.getRscriptExecutable(new Logger()) + " --no-save [%0]";
        Qsub.qsub(batchDir + "run_additionals", batchDir, -1, commands, iterations, qsubMem,
                  qsubWalltime);
        if (iterations.length == 0) {
          new File(batchDir
                   + "master.run_additionals").renameTo(new File(batchDir + "master.run_additionals.bak"));
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
        new File(cohort + "_" + race + "_" + pheno + ".RData").renameTo(new File(cohort + insert
                                                                                 + "_" + race + "_"
                                                                                 + pheno
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
    String queue = null;
    boolean usePrep2 = true;
    String saveName = null;
    String markers = null;
    String conditionals = null;

    cohort = "ARIC";
    genos = "D:/SkatMeta/genotypes_blacks_AA/AA_ARIC_noJHS_chr#t.csv";
    pheno = "D:/SkatMeta/results_hemostasis/pheno_F7_studyIDs.csv";
    snpInfo = "N:/statgen/skatMeta/fullExample/SNPInfo_HumanExome-12v1_rev5_justGeneSpliced.csv";
    qsubMem = 15000;
    qsubWalltime = 2;

    String usage = "\n" + "gwas.SeqMetaPrimary requires 4 arguments\n"
                   + "   (1) cohort name (i.e. cohort=" + cohort + " (not the default))\n"
                   + "   (2) genotype file name (i.e. geno=" + genos + " (not the default))\n"
                   + "   (3) phenotype file name (i.e. pheno=" + pheno + " (not the default))\n"
                   + "   (4) snpInfo file name (i.e. snpInfo=" + snpInfo + " (not the default))\n"
                   + "   (5) qsub memory in megabytes (i.e. qsubmem=" + qsubMem + " (default))\n"
                   + "   (6) qsub walltime in hours (i.e. qsubwalltime=" + qsubWalltime
                   + " (default))\n" + "   (7) queue to use (i.e. queue=" + queue + " (default))\n"
                   + "   (8) (optional) phenotypes to run (i.e. phenos=CRP,ICAM1,TNFA (not the default))\n"
                   + "   (9) (optional) races to run (i.e. races=" + races + " (default))\n"
                   + " ( the final optional parameters require a fixed phenotype format of [cohort name]_[race]_[phenotype]_.csv )\n"
                   + " ( as well as a genotype file with the format of regular_filename_[%race]_chr#.csv )\n"
                   + "   (10) (optional) run additional models for Rdata files (i.e. -additionals  (not the default))\n"
                   + " OR\n"
                   + "   (10) (optional) rename Rdata files with this insert (i.e. rename=42PCs (not the default))\n"
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
      } else if (arg.startsWith("queue=")) {
        queue = ext.parseStringArg(arg, null);
      } else if (arg.startsWith("-prep2")) {
        usePrep2 = false;
      } else if (arg.startsWith("savename=")) {
        saveName = ext.parseStringArg(arg);
      } else if (arg.startsWith("markers=")) {
        markers = ext.parseStringArg(arg);
      } else if (arg.startsWith("conditionals=")) {
        conditionals = ext.parseStringArg(arg);
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }

    if (markers != null) {
      dumpMarkers(phenos, cohort, genos, markers, races);

    } else if (phenos != null) {
      if (rename != null) {
        renameMany(cohort, phenos, races, rename);
      } else if (additionals) {
        additionalModels(cohort, phenos, snpInfo, qsubMem, qsubWalltime);
      } else {
        batchMany(cohort, genos, phenos, races, snpInfo, qsubMem, qsubWalltime, queue, usePrep2,
                  saveName, conditionals, null);
      }
    } else {
      batch(cohort, genos, pheno, snpInfo, qsubMem, qsubWalltime, queue, usePrep2, saveName,
            conditionals, null);
    }
  }

  private static String dumpPhenosWithCovars(String phenoFilename, String cohort, String genos,
                                             String markers, String race) {
    String phenoRoot, outfile, loop, filename;
    String[] g;
    boolean byChr;

    // TODO Auto-generated method stub
    genos = ext.replaceAllWith(genos, "[%race]", race);

    phenoRoot = ext.rootOf(phenoFilename);

    g = genos.contains("#") ? null : genos.split("#");

    byChr = g != null;

    outfile = phenoRoot + "_covars.csv";
    String genoLoadCmd = "geno <- "
                         + (genos.endsWith(".csv") ? "t(read.csv(name, header=T, as.is=T, row.names=1))"
                                                   : "load(name)");

    if (byChr) {
      String genoFile = "paste(" + g[0] + ", c, " + g[1] + ", sep=\"\")";

      loop = "for(c in chrs) {\n" + "  m <- markers$Marker[markers$Chr == c]\n" + "  name <- "
             + genoFile + "\n" + "  " + genoLoadCmd + "\n"
             + "  results <- geno[,colnames(geno) %in% m]\n"
             + "  final <- merge(xphen, results, by=\"row.names\")" + "  out <- paste(\"" + outfile
             + "_chr\", c, \".csv\", sep=\"\")\n"
             + "  write.table(final, out, sep=\",\", row.names=T, col.names=T)\n" + "}";
    } else {
      loop = "name <- " + genos + "\n " + genoLoadCmd + "\n"
             + "results <- geno[,colnames(geno) %in% markers$Marker]\n"
             + "final <- merge(xphen, results, by=\"row.names\")" + "write.table(final, " + outfile
             + ".csv, sep=\",\", row.names=T, col.names=T)\n";
    }

    String rCode = "setwd(\"" + ext.pwd() + "\")\n" + "markers <- read.table(\"" + markers
                   + "\", col.names=c(\"Markers\", \"Chr\")) \n" + "chrs <- unique(markers$Chr)\n"
                   + "names <- markers$Markers\n" + "pheno <- read.csv(\"" + phenoFilename
                   + "\", header=T, as.is=T, row.names=1)\n" + "xphen <- na.omit(pheno)\n" + loop;

    try {
      filename = phenoRoot + "_dumpCovars.R";
      PrintWriter out = new PrintWriter(new FileOutputStream(filename));
      out.println(rCode);
      out.close();
      return filename;

    } catch (Exception e) {
      e.printStackTrace();
      System.exit(1);
    }

    return null;
  }

  private static void dumpMarkers(String phenosCommaDelimited, String cohort, String genos,
                                  String markers, String racesCommaDelimited) {
    String[] phenos, races;
    String phenoFileName;
    Vector<String> v;

    phenos = phenosCommaDelimited == null ? null : phenosCommaDelimited.split(",");
    races = racesCommaDelimited.split(",");

    v = new Vector<String>();

    for (String r : races) {
      v.add(dumpByRace(genos, markers, r));
      if (phenos != null) {
        for (String p : phenos) {
          phenoFileName = ext.pwd() + cohort + "_" + r + "_" + p + ".csv";
          v.add(dumpPhenosWithCovars(phenoFileName, cohort, genos, markers, r));
        }
      }

      Vector<String> jobNamesWithAbsolutePaths = new Vector<String>();
      IntVector jobSizes = new IntVector();

      for (String f : v) {
        jobNamesWithAbsolutePaths.add(Rscript.getRscriptExecutable(new Logger()) + " --no-save "
                                      + f);
        jobSizes.add(8000);
      }

      Qsub.qsubExecutor(ext.pwd(), jobNamesWithAbsolutePaths, jobSizes, "run_dumpMarkers", 1, 62000,
                        2);
    }
  }

  private static String dumpByRace(String genos, String markers, String race) {
    String[] g;
    String rCode, genoFile, dir, outfile, loop, filename;
    boolean byChr = false;

    dir = ext.pwd() + "markersOfInterest/";

    if (!new File(dir).exists() || !new File(dir).isDirectory()) {
      new File(dir).mkdirs();
    }

    genos = ext.replaceAllWith(genos, "[%race]", race);
    g = genos.contains("#") ? null : genos.split("#");

    byChr = g != null;

    outfile = race + "_markersOfInterest";
    String genoLoadCmd = "geno <- "
                         + (genos.endsWith(".csv") ? "t(read.csv(name, header=T, as.is=T, row.names=1))"
                                                   : "load(name)");

    if (byChr) {
      genoFile = "paste(" + g[0] + ", c, " + g[1] + ", sep=\"\")";

      loop = "for(c in chrs) {\n" + "  m <- markers$Marker[markers$Chr == c]\n" + "  name <- "
             + genoFile + "\n" + "  " + genoLoadCmd + "\n"
             + "  results <- geno[rownames(geno) %in% m,]\n" + "  out <- paste(\"" + outfile
             + "_chr\", c, \".csv\", sep=\"\")\n"
             + "  write.table(results, out, sep=\",\", row.names=T, col.names=T)\n" + "}";
    } else {
      loop = "name <-" + genos + "\n" + genoLoadCmd + "\n"
             + "results <- geno[, colnames(geno) %in% markers$Marker]\n" + "write.table(results, "
             + outfile + ".csv, sep=\",\", row.names=T, col.names=T)\n";
    }

    rCode = "setwd(\"" + dir + "\")" + "markers <- read.table(\"" + markers
            + "\", col.names=c(\"Markers\", \"Chr\")) \n" + "chrs <- unique(markers$Chr)\n"
            + "names <- markers$Markers\n" + loop;

    try {
      filename = race + "_dumpMarkers.R";
      PrintWriter out = new PrintWriter(new FileOutputStream(filename));
      out.println(rCode);
      out.close();
      return filename;

    } catch (Exception e) {
      e.printStackTrace();
      System.exit(1);
    }

    return null;
  }

  private static String conditionals(String condFile, int chr) {
    String rdata = "";

    if (condFile != null) {
      rdata = "condMarkers <- read.table(\"" + condFile + "\", header=T)\n" + "\n"
              + "snpsonchr <- SNPInfo[SNPInfo$chr == " + chr + ",]\n"
              + "for(i in 1:nrow(condMarkers)) {\n" + "  if (condMarkers$chr[i] == " + chr + "){\n"
              + "    cond <- snpsonchr[(snpsonchr$pos >= condMarkers$start[i] & "
              + "snpsonchr$pos <= condMarkers$end[i]), c(\"SNP\", \"SKATgene\")]\n"
              + "    cond$SNP <- rep(condMarkers$SNP[i], nrow(cond))\n" + "    cond<-unique(cond)\n"
              + "    conditions <- rbind(conditions, cond)\n"
              + "    if (any(colnames(mPheno) %in% condMarkers$SNP[i])) {\n"
              + "       col <- which(colnames(mPheno) %in% condMarkers$SNP[i])\n"
              + "       mGeno$NEW <- mPheno[,col]\n"
              + "       colnames(mGeno)[colnames(mGeno) == \"NEW\"] <-  colnames(mPheno)[col]\n"
              + "       if (col == ncol(mPheno)) {\n" + "           mPheno <- mPheno[, 1:(col-1)]\n"
              + "       } else { \n"
              + "           mPheno <- mPheno[, c(1:(col-1), (col+1):ncol(mPheno))]\n" + "       }\n"
              + "    }\n" + "  } else if (any(colnames(mPheno) %in% condMarkers$SNP[i])) {\n"
              + "    col <- which(colnames(mPheno) %in% condMarkers$SNP[i])\n"
              + "    if (col == ncol(mPheno)) {\n" + "        mPheno <- mPheno[, 1:(col-1)]\n"
              + "    } else { \n"
              + "        mPheno <- mPheno[, c(1:(col-1), (col+1):ncol(mPheno))]\n" + "    }\n"
              + "  }\n" + "}\n";
    }

    return rdata;
  }

}
