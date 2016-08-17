package org.genvisis.one;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FilenameFilter;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Vector;

import org.genvisis.common.Elision;
import org.genvisis.common.Files;
import org.genvisis.common.IntVector;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class CALiCo_SOL {

  public static void generateQsubMultiple(Vector<String> jobNamesWithAbsolutePaths, String batchDir,
      String modelName, int qsubMemInMBs, double qsubWalltimeInHours, Logger log) {
    IntVector jobSizes;
    Vector<String> GenoFileNames;
    Vector<Integer> fileSizes;
    String line;
    File file;
    int index;
    int genoFileSize;

    jobSizes = new IntVector();
    GenoFileNames = new Vector<String>();
    fileSizes = new Vector<Integer>();
    for (int i = 0; i < jobNamesWithAbsolutePaths.size(); i++) {
      line = jobNamesWithAbsolutePaths.elementAt(i);
      line = line.substring(line.indexOf(" --GENO " + 8));
      line = line.substring(0, line.indexOf(""));
      index = GenoFileNames.indexOf(line);
      if (index >= 0) {
        ;
        jobSizes.add(fileSizes.elementAt(index));
      } else {
        file = new File(line);
        if (file.exists()) {
          genoFileSize = Integer.MAX_VALUE - (int) file.length();
          jobSizes.add(genoFileSize);
          GenoFileNames.add(line);
          fileSizes.add(genoFileSize);
        } else {
          log.reportError("Error - file not found: " + file.getAbsolutePath());
          System.exit(0);
        }
      }
    }

    Files.qsubMultiple(jobNamesWithAbsolutePaths, jobSizes, batchDir,
        batchDir + "chunk_" + modelName, 8, true, null, -1, qsubMemInMBs, qsubWalltimeInHours);
    Files.qsubMultiple(jobNamesWithAbsolutePaths, jobSizes, batchDir,
        batchDir + "chunkSB256_" + modelName, 16, true, "sb256", -1, qsubMemInMBs,
        qsubWalltimeInHours);
    Files.qsubMultiple(jobNamesWithAbsolutePaths, jobSizes, batchDir,
        batchDir + "chunkSB_" + modelName, 16, true, "sb", -1, qsubMemInMBs, qsubWalltimeInHours);
  }

  public static void generateQsubs(String commands, String root, String[][] iterations,
      Logger log) {
    Files.qsub(root, commands, iterations);
  }

  // public static void generateQsubMultiple(Vector<String> jobNamesWithAbsolutePaths, String
  // batchDir, String modelName, int qsubMemInMBs, double qsubWalltimeInHours, Logger log) {
  // Files.qsubMultiple(jobNamesWithAbsolutePaths, null, batchDir, batchDir+"chunk_"+modelName, 8,
  // true, null, -1, qsubMemInMBs, qsubWalltimeInHours);
  // Files.qsubMultiple(jobNamesWithAbsolutePaths, null, batchDir, batchDir+"chunkSB256_"+modelName,
  // 16, true, "sb256", -1, qsubMemInMBs, qsubWalltimeInHours);
  // Files.qsubMultiple(jobNamesWithAbsolutePaths, null, batchDir, batchDir+"chunkSB_"+modelName,
  // 16, true, "sb", -1, qsubMemInMBs, qsubWalltimeInHours);
  // }

  public static Vector<String> generateScripts(String commandFullPath, String[] methods,
      String phenoDir, String phenoNameExtOrFullPath, String genoDir, String genoNameExtOrFullPath,
      String probMatrixDir, String probMatrixNameExtOrFullPath, String outDir, String batchDir,
      Logger log) {
    PrintWriter writer;
    File[] genoFiles = null;
    File[] phenoFiles = null;
    File[] probMatrixFiles = null;
    String[] genoLabels;
    String[] phenoLabels;
    String[] probMatrixLabels;
    String resultFullPath;
    Vector<String> script;

    if (log == null) {
      log = new Logger();
    }

    try {
      genoFiles = getListOfFilesInADirectory(genoDir, genoNameExtOrFullPath, log);
      phenoFiles = getListOfFilesInADirectory(phenoDir, phenoNameExtOrFullPath, log);
      probMatrixFiles = getListOfFilesInADirectory(probMatrixDir, probMatrixNameExtOrFullPath, log);
    } catch (Elision e) {
      e.printStackTrace();
      System.exit(1);
    }
    genoLabels = getLabels(genoFiles, log);
    phenoLabels = getLabels(phenoFiles, log);
    probMatrixLabels = getLabels(probMatrixFiles, log);

    if (methods == null) {
      methods = new String[] {"LS"};
    }

    script = new Vector<String>();
    for (int i = 0; (phenoFiles != null) && (i < phenoFiles.length); i++) {
      for (int j = 0; j < genoFiles.length; j++) {
        for (String method : methods) {
          for (int l = 0; l < probMatrixFiles.length; l++) {
            resultFullPath = outDir + "Pheno_" + phenoLabels[i] + "_Geno_" + genoLabels[j]
                + "_Method_" + method + "_PMatrix_" + probMatrixLabels[l];
            script.add(commandFullPath + " --METHOD " + method + " --PHENO "
                + phenoFiles[i].getAbsolutePath() + " --GENO " + genoFiles[j].getAbsolutePath()
                + " --PROBMATRIX " + probMatrixFiles[l].getAbsolutePath() + " --OUT "
                + resultFullPath + ".out --LOG " + resultFullPath
                + ".log --CAF 0.01 --TRIM 0.02 --TRIM2 0.0004");
          }
        }
      }
    }

    if (batchDir != null) {
      if (!batchDir.endsWith("/")) {
        batchDir += "/";
      }
      try {
        writer = new PrintWriter(new FileOutputStream(batchDir + "script.txt"));
        for (int i = 0; i < script.size(); i++) {
          writer.println(script.elementAt(i));
        }
        writer.close();
      } catch (FileNotFoundException e) {
        e.printStackTrace();
      }
    }
    return script;
  }

  /**
   * Given the file names, find the minimum length segment of each of them that distinguish it.
   * 
   * @param fileFullPaths
   * @param log
   * @return
   */
  public static String[] getLabels(File[] fileFullPaths, Logger log) {
    String[] probMatrixLabels;
    int len;
    int start;
    int startMin;
    int startMax;
    int strLength;

    if (fileFullPaths.length == 1) {
      probMatrixLabels = new String[] {ext.rootOf(fileFullPaths[0].getName())};
    } else {
      startMin = 0;
      probMatrixLabels = new String[fileFullPaths.length];
      for (int i = 0; i < probMatrixLabels.length; i++) {
        // probMatrixLabels[i] = ext.rootOf(probMatrixFullPathes[i].getName()).toLowerCase();
        probMatrixLabels[i] = ext.rootOf(fileFullPaths[i].getName());
        if (probMatrixLabels[i].length() > startMin) {
          startMin = probMatrixLabels[i].length();
        }
      }

      startMax = 0;
      strLength = 1;
      for (int i = 0; i < (probMatrixLabels.length - 1); i++) {
        start = probMatrixLabels[i].length();
        for (int j = (i + 1); j < probMatrixLabels.length; j++) {
          len = Math.min(probMatrixLabels[i].length(), probMatrixLabels[j].length());
          for (int k = 0; k < len; k++) {
            if (!(probMatrixLabels[i].substring(k, k + strLength)
                .equals(probMatrixLabels[j].substring(k, k + strLength))) && k < start) {
              start = k;
              break;
            }
          }
        }
        if (startMin > start) {
          startMin = start;
        }
        if (startMax < start) {
          startMax = start;
        }
      }
      startMax++;

      for (int i = 0; i < probMatrixLabels.length; i++) {
        probMatrixLabels[i] = probMatrixLabels[i].substring(startMin,
            Math.min(startMax, probMatrixLabels[i].length()));
      }
    }

    return probMatrixLabels;
  }

  // public static String[] getLabels(File[] probMatrixFullPathes, Logger log) {
  // String[] probMatrixLabels;
  // String probMatrixName1;
  // String probMatrixName2;
  // int len;
  // int strLength;
  // boolean isFound;
  //
  // if (probMatrixFullPathes.length == 1) {
  // probMatrixLabels = new String[] {ext.rootOf(probMatrixFullPathes[0].getName())};
  // } else {
  // strLength = 1;
  // probMatrixLabels = new String[probMatrixFullPathes.length];
  // probMatrixName1 = ext.rootOf(probMatrixFullPathes[0].getName()).toLowerCase();
  // for (int i = 1; i < probMatrixLabels.length; i++) {
  // probMatrixName2 = ext.rootOf(probMatrixFullPathes[i].getName()).toLowerCase();
  // len = Math.min(probMatrixName1.length(), probMatrixName2.length());
  // for (int j = 0; j < len; j++) {
  // isFound = false;
  // if(!(probMatrixName2.substring(j, j + strLength).equals(probMatrixName1.substring(j, j +
  // strLength)))) {
  // if(i == 1) {
  // probMatrixLabels[0] = probMatrixName1.charAt(j) + "";
  // if (probMatrixLabels.length == 2) {
  // probMatrixLabels[1] = probMatrixName2.charAt(j) + "";
  // }
  // } else {
  // for (int k = 0; k < (i-1); k++) {
  // if(probMatrixLabels[k].equals(probMatrixName1.charAt(j) + "")) {
  // isFound = true;
  // break;
  // }
  // }
  // if (!isFound) {
  // probMatrixLabels[i-1] = probMatrixName1.charAt(j) + "";
  // if ((i+1) == probMatrixLabels.length) {
  // probMatrixLabels[i] = probMatrixName2.charAt(j) + "";
  // }
  // break;
  // } else {
  // strLength ++;
  // i = 1;
  // }
  // }
  // }
  // }
  // probMatrixName1 = probMatrixName2;
  // }
  // }
  //
  // return probMatrixLabels;
  // }

  /**
   * List all the files in a directory with the file name extension specified. Alternatively, you
   * can leave the directory as null, and specify a full path in filenameOrNameExtToList, then the
   * method will check whether the file specified by exists.
   * 
   * @param directory The directory to list files in.
   * @param filenameOrNameExtToList There ways to use this parameter: 1) leave as null, then you
   *        will get list of all the files; 2) specify some file name extension, then you will get
   *        list of the files with the file name extension; 3) specify the full path of a file while
   *        leave directory as null, then you will know whether the file exists.
   * @param log
   * @return
   * @throws Elision
   */
  public static File[] getListOfFilesInADirectory(String directory,
      final String filenameOrNameExtToList, Logger log) throws Elision {
    File[] phenoFiles = null;

    if (log == null) {
      log = new Logger();
    }

    if (directory != null) {
      if (!directory.endsWith("/")) {
        directory += "/";
      }
      phenoFiles = new File(directory).listFiles(new FilenameFilter() {
        @Override
        public boolean accept(File file, String filename) {
          if (filenameOrNameExtToList == null) {
            return true;
          } else if (filename.toLowerCase().endsWith(filenameOrNameExtToList.toLowerCase())) {
            return true;
          } else {
            return false;
          }
        }
      });
    } else if (filenameOrNameExtToList == null) {
      throw new Elision("Error - neither directory nor file name is specified.");
    } else {
      phenoFiles = new File[] {new File(filenameOrNameExtToList)};
      if (!phenoFiles[0].exists()) {
        throw new Elision("Error - file not found: " + filenameOrNameExtToList);
      }
    }

    return phenoFiles;
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String commandFullPath;
    String[] methods;
    String phenoDir;
    String phenoNameExtOrFullPath;
    String genoDir;
    String genoNameExtOrFullPath;
    String probMatrixDir;
    String probMatrixNameExtOrFullPath;
    String outDir;
    String batchDir;
    int qsubMemInMBs;
    double qsubWalltimeInHours;
    Logger log;
    boolean isExisted;
    boolean isSucceeded;
    Vector<String> jobNamesWithAbsolutePaths;

    commandFullPath = "/home/pankrat2/shared/bin/SOLReg";
    methods = new String[] {"LS", "logistic"};
    // phenoDir = "/home/pankrat2/shared/CALiCo_SOL/T2DM/xln_files/";
    phenoDir = "D:/CALiCo_SOL/T2DM/xln_files/";
    phenoNameExtOrFullPath = ".dat";
    // genoDir = "/home/pankrat2/shared/bin/installs/sol/solregfiles_131031/geno/";
    genoDir = "D:/CALiCo_SOL/solregfiles_131031/solregfiles_131031/geno/";
    genoNameExtOrFullPath = ".txt";
    // probMatrixDir = "/home/pankrat2/shared/bin/installs/sol/solregfiles_131031/sampling_prob/";
    probMatrixDir = "D:/CALiCo_SOL/sampling_prob/";
    probMatrixNameExtOrFullPath = ".txt";
    outDir = "/home/pankrat2/shared/CALiCo_SOL/T2DM/results/";
    // batchDir = "/home/pankrat2/shared/CALiCo_SOL/T2DM/batches/";
    batchDir = "D:/CALiCo_SOL/";
    qsubMemInMBs = 15000;
    qsubWalltimeInHours = 24;

    String usage = "\n" + "one.CALiCo_SOL requires 10 - 12 arguments\n"
        + "   (1) full path of the command (i.e. command=" + commandFullPath + " (not default))\n"
        + "   (2) methods of analysis (i.e. methods=" + methods[0] + ";" + methods[1]
        + " (add additional ones separated by ',' or ';' with no space) (default))\n"
        + "   (3) directory of the pheno files (i.e. phenodir=" + phenoDir + " (not default))\n"
        + "   (4) name extension of the pheno files, or full path of a specific pheno file (i.e. phenoext="
        + phenoNameExtOrFullPath + " (not default))\n"
        + "   (5) directory of the pheno files (i.e. genodir=" + genoDir + " (not default))\n"
        + "   (6) name extension of the pheno files, or full path of a specific geno file (i.e. genoext="
        + genoNameExtOrFullPath + " (not default))\n"
        + "   (7) directory of the probability matrix files (i.e. probmaxtrixdir=" + probMatrixDir
        + " (not default))\n"
        + "   (8) name extension of the probability matrix files, or full path of a specific probability matrix file (i.e. probmatrixext="
        + probMatrixNameExtOrFullPath + " (not default))\n"
        + "   (9) directory of output files (i.e. outdir=" + outDir + " (not default))\n"
        + "   (10) directory of batch files (i.e. batchdir=" + batchDir + " (not default))\n"
        + "   (11) (optional) qsub memory size (i.e. qsubmem=" + qsubMemInMBs
        + " (default; in megabytes))\n" + "   (12) (optional) qsub walltime (i.e. qsubwalltime="
        + qsubWalltimeInHours + " (default; in hours))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("command=")) {
        commandFullPath = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("methods=")) {
        methods = arg.split("=")[1].split("[,;]");
        numArgs--;
      } else if (arg.startsWith("phenodir=")) {
        phenoDir = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("phenoext=")) {
        phenoNameExtOrFullPath = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("genodir=")) {
        genoDir = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("genoext=")) {
        genoNameExtOrFullPath = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("probmatrixdir=")) {
        probMatrixDir = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("probmatrixext=")) {
        probMatrixNameExtOrFullPath = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("outdir=")) {
        outDir = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("batchdir=")) {
        batchDir = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("qsubmem=")) {
        qsubMemInMBs = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      } else if (arg.startsWith("qsubwalltime=")) {
        qsubWalltimeInHours = Double.parseDouble(arg.split("=")[1]);
        numArgs--;
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }

    isExisted = false;
    isSucceeded = true;
    try {
      if (!new File(outDir).exists()) {
        if (!new File(outDir).mkdir()) {
          log = new Logger();
          isSucceeded = false;
        } else {
          log = new Logger(outDir + "Genvisis_CALiCo_SOL_"
              + (new SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date())) + ".log");
          isExisted = true;
        }
      } else {
        log = new Logger(outDir + "Genvisis_CALiCo_SOL_"
            + (new SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date())) + ".log");
      }

      log.report("Genvisis (R) 2014. \nCALiCo_SOL analysis "
          + (new SimpleDateFormat("MM/dd/yyyy HH:mm:ss").format(new Date()))
          + "\n-Pheno files directory: " + phenoDir + "\n-Pheno file name extension or full path: "
          + phenoNameExtOrFullPath + "\n-Geno files directory: " + genoDir
          + "\n-Geno file name extension or full path: " + genoNameExtOrFullPath
          + "\n-Probability matrix files directory: " + probMatrixDir
          + "\n-Probabliity matrix file name extension or full path: " + probMatrixNameExtOrFullPath
          + "\n-Output files directory: " + outDir + "\n-Batch files directory: " + batchDir);

      if (isExisted) {
        log.reportError("Warning --- Directory " + outDir
            + " already exists. Existing files might be reused or overwritten.");
      } else if (!isSucceeded) {
        log.reportError("Warning --- Cannot create the directory " + outDir);
      } else {
        log.report("Creating result directory " + outDir);
      }

      jobNamesWithAbsolutePaths = generateScripts(commandFullPath, methods, phenoDir,
          phenoNameExtOrFullPath, genoDir, genoNameExtOrFullPath, probMatrixDir,
          probMatrixNameExtOrFullPath, outDir, batchDir, log);
      generateQsubMultiple(jobNamesWithAbsolutePaths, batchDir, "CALiCo_SOL_ver6", qsubMemInMBs,
          qsubWalltimeInHours, log);
      log.report("\nCalico_SOL analysis is finished.");

    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
