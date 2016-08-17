package org.genvisis.cnv.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.Array;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Positions;
import org.genvisis.common.StringVector;
import org.genvisis.common.ext;
import org.genvisis.filesys.CNVariant;

public class DeNovoCNV {

  // private static Project proj;

  public static void batch(String pedigreeOfTrio) {
    String command;
    String[][] iterations;

    // command = "perl C:/penncnv/detect_cnv.pl -joint -hmm C:/penncnv/lib/hh550.hmm -pfb
    // C:/penncnv/lib/hh550.hg18.pfb -gcmodel C:/penncnv/lib/hh550.hg18.gcmodel [%1] [%2] [%3] -out
    // [%1].jointcnv";
    command =
        "perl ../bin/penncnv/detect_cnv.pl -joint -hmm ../bin/penncnv/lib/hhall.hmm -pfb ../bin/custom.pfb -gcmodel ../bin/custom.gcmodel [%1] [%2] [%0] -out ../results/[%0].jointcnv -log ../results/[%0].log";

    iterations = HashVec.loadFileToStringMatrix(pedigreeOfTrio, true, new int[] {4, 5, 6}, false);

    org.genvisis.common.Files.qsub("denovo", "/share/bulk/gedi/pankr018/denovo/penn_data", 65,
                                   command, iterations, 2500, 2);
  }

  /*
   * 12/6/2013: Cannot remember what this method does. First, the same function is available in
   * another method, parsePennCnvResults(...). Second, looks like it is incomplete, - there is no
   * any writer in this method.
   * 
   * @param pennCnvResultDir
   * 
   * @param cnvFileExtensions
   * 
   * @param log
   */
  public static void detectDenovoCnv(String pennCnvResultDir, String cnvFileExtensions,
                                     Logger log) {
    File dir;
    File[] listOfFiles;
    BufferedReader reader;
    String[] line;
    Vector<String[]> cnv;
    int j, k;
    boolean found;

    if (log == null) {
      log = new Logger();
    }
    dir = new File(pennCnvResultDir);
    listOfFiles = dir.listFiles();

    for (File listOfFile : listOfFiles) {
      if (listOfFile.isFile() && listOfFile.toString().contains(cnvFileExtensions)) {
        cnv = new Vector<String[]>();
        try {
          reader = new BufferedReader(new FileReader(listOfFile.getPath()));
          while (reader.ready()) {
            line = reader.readLine().trim().split("[\\s]+");
            if (!line[0].contains("NOTICE")) {
              cnv.add(line);
            } else {
              break;
            }
          }
          reader.close();
        } catch (FileNotFoundException fnfe) {
          log.reportError("Error: file not found " + listOfFile.getPath());
        } catch (IOException ioe) {
          log.reportError("Error reading " + listOfFile.getPath());
        }
        j = 0;
        while (j < cnv.size()) {
          if (cnv.elementAt(j)[7].contains("offspring")) {
            break;
          }
          j++;
        }
        k = j;
        while (j < cnv.size()) {
          found = false;
          for (int l = 0; l < k; l++) {
            if (cnv.elementAt(j)[0].split(":")[0].equals(cnv.elementAt(l)[0].split(":")[0])
                && ((Integer.parseInt(cnv.elementAt(j)[0].split(":")[1].split("-")[0]) >= Integer.parseInt(cnv.elementAt(l)[0].split(":")[1].split("-")[0])
                     && Integer.parseInt(cnv.elementAt(j)[0].split(":")[1].split("-")[0]) <= Integer.parseInt(cnv.elementAt(l)[0].split(":")[1].split("-")[1]))
                    || (Integer.parseInt(cnv.elementAt(j)[0].split(":")[1].split("-")[1]) >= Integer.parseInt(cnv.elementAt(l)[0].split(":")[1].split("-")[0])
                        && Integer.parseInt(cnv.elementAt(j)[0].split(":")[1].split("-")[1]) <= Integer.parseInt(cnv.elementAt(l)[0].split(":")[1].split("-")[1]))
                    || (Integer.parseInt(cnv.elementAt(l)[0].split(":")[1].split("-")[0]) >= Integer.parseInt(cnv.elementAt(j)[0].split(":")[1].split("-")[0])
                        && Integer.parseInt(cnv.elementAt(l)[0].split(":")[1].split("-")[0]) <= Integer.parseInt(cnv.elementAt(j)[0].split(":")[1].split("-")[1]))
                    || (Integer.parseInt(cnv.elementAt(l)[0].split(":")[1].split("-")[1]) >= Integer.parseInt(cnv.elementAt(j)[0].split(":")[1].split("-")[0])
                        && Integer.parseInt(cnv.elementAt(l)[0].split(":")[1].split("-")[1]) <= Integer.parseInt(cnv.elementAt(j)[0].split(":")[1].split("-")[1])))) {
              found = true;
              break;
            }
          }
          if (!found) {
            log.report(listOfFile.getName() + "\t" + cnv.elementAt(j)[0]);
          }
          j++;
        }
      } else if (listOfFile.isDirectory()) {
        log.report("Directory " + listOfFile.getName());
      }
    }
    log.report(ext.getTime() + "\tgenerated De Novo CNV output for " + pennCnvResultDir + "*"
               + cnvFileExtensions);
  }

  public static void generateFamFileForTrios(Project proj, String triosPedigreeFileFullPath,
                                             boolean pedigreeHasHeader) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    SampleData sampleData;

    sampleData = proj.getSampleData(2, false);
    try {
      reader = new BufferedReader(new FileReader(triosPedigreeFileFullPath));
      if (pedigreeHasHeader) {
        reader.readLine();
      }
      writer =
          new PrintWriter(new FileWriter(proj.DATA_DIRECTORY.getValue(false, true) + "trios.fam"));
      while (reader.ready()) {
        line = reader.readLine().split("\t");
        writer.println(line[0] + "\t" + line[1] + "\t" + line[2] + "\t" + line[3] + "\t"
                       + sampleData.getSexForIndividual(line[4]) + "\t1");
      }
      writer.close();
      reader.close();

    } catch (FileNotFoundException fnfe) {
      return;
    } catch (IOException ioe) {
      return;
    }
  }

  public static void generateListOfTrios(String input_pedigreeFullPath,
                                         String output_listOfTrioSetsFullPath, Logger log) {
    BufferedReader reader;
    String[] line;
    PrintWriter writer;
    StringVector fId, iId, faId, moId, dna;

    // log.report("Generating list of trio sets.");
    // trav = new HashVec();
    try {
      // reader = Files.getReader(new FileReader(filename), proj.getJarStatus(), true, false);
      reader = new BufferedReader(new FileReader(input_pedigreeFullPath));
      // reader.readLine();
      // TODO detect the column number of FID IID FAID MOID
      fId = new StringVector();
      iId = new StringVector();
      faId = new StringVector();
      moId = new StringVector();
      dna = new StringVector();
      while (reader.ready()) {
        line = reader.readLine().trim().split("\t", -1);
        // if (line[6]!=null && !line[6].equals(".")) {
        if (line[4] != null && !line[4].equals(".")) {
          fId.add(line[0]);
          iId.add(line[0] + "\t" + line[1]);
          faId.add(line[2]);
          moId.add(line[3]);
          // faId.add(line[0]+"\t"+line[2]);
          // moId.add(line[0]+"\t"+line[3]);
          dna.add(line[6]);
          // dna.add(line[4]);
        }
      }
      reader.close();

      // writer = new PrintWriter(new FileWriter(ext.parseDirectoryOfFile(input_pedigreeFullPath) +
      // "PedigreeOfTrios.txt"));
      writer = new PrintWriter(new FileWriter(output_listOfTrioSetsFullPath));
      writer.println("fId\tiId\tfaId\tmoId\tiDna\tfaDna\tmoDna");
      for (int i = 0; i < iId.size(); i++) {
        if (iId.contains(fId.elementAt(i) + "\t" + faId.elementAt(i))
            && iId.contains(fId.elementAt(i) + "\t" + moId.elementAt(i))) {
          writer.println(iId.elementAt(i) + "\t" + faId.elementAt(i) + "\t" + moId.elementAt(i)
                         + "\t" + dna.elementAt(i) + "\t"
                         + dna.elementAt(iId.indexOf(fId.elementAt(i) + "\t" + faId.elementAt(i)))
                         + "\t"
                         + dna.elementAt(iId.indexOf(fId.elementAt(i) + "\t" + moId.elementAt(i))));
        }
      }
      writer.flush();
      writer.close();
    } catch (FileNotFoundException fnfe) {
      log.reportError("Error: file \"" + input_pedigreeFullPath
                      + "\" not found in current directory");
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + input_pedigreeFullPath + "\"");
    }
    log.report(ext.getTime() + "\tgenerated list of trio sets " + output_listOfTrioSetsFullPath);
  }

  public static void main(String... args) {
    String projPropertyFileFullPath;
    Project proj;
    String gcBaseFileFullPath;
    String gcModelFullPath;
    int numWindowUnits;
    String pedigreeFullPath;
    String pennDataDir;
    String pennOutputDir;
    String pennOutputFilenameRoot;
    String pennBinDir;
    String trioPedigreeFullPath;
    boolean parseOrNot;
    Logger log;

    // projPropertyFileFullPath = "C:/workspace/Genvisis/projects/OSv2.properties";
    projPropertyFileFullPath = "/home/pankrat2/zxu/projects/gedi_gwas.properties";
    proj = new Project(projPropertyFileFullPath, false);
    pedigreeFullPath = proj.DATA_DIRECTORY.getValue(false, true) + "pedigree.dat";
    trioPedigreeFullPath = ext.parseDirectoryOfFile(pedigreeFullPath) + "pedigreeTrios.dat";
    // gcBaseFileFullPath = "D:/PennCNV_Related/GcBase/gc5Base_hg19.txt";
    gcBaseFileFullPath = "/home/pankrat2/zxu/PennCNV_Related/GcBase/gc5Base_hg19.txt";
    numWindowUnits = 0;
    pennOutputFilenameRoot = "gedi_all.rawcnv";
    pennBinDir = "C:/penncnv/";
    pennDataDir = proj.PENNCNV_DATA_DIRECTORY.getValue(false, true);
    pennOutputDir = proj.PENNCNV_RESULTS_DIRECTORY.getValue(false, true);
    parseOrNot = false;

    String usage =
        "\n" + "cnv.analysis.DeNovoCNV analysis normally contains the following steps, split into several rounds, within each of which all the steps will be finished by just calling the program once."
                   + "1st round:" + "   (1) Search the pedigree file for list of trio sets;"
                   + "   (2) From .sampRAF data generate 'custom.pbf';"
                   + "   (3) From GC Base standarded file, generate custom GC model 'custom.gcmodel';"
                   + "   (4) For Linux/Unix users: generate scripts to run PennCNV for trio CNV and joint CNV detection;"
                   + "       For Windows users: call PennCNV to run trio CNV and joint CNV detection;"
                   + "   (5) For Linux/Unix users: the programm ends here, and you have to mannually run PennCNV using the scripts either through interactive Bash mode or backend batch mode;"
                   + "       For Windows users: wait till the PennCNV finishes;" + "" + "2nd round:"
                   + "   (6) Parse the PennCNV results for De Novo CNV;" + ""
                   + "Pre-requisition of cnv.analysis.DeNovoCNV:"
                   + "   (1) Entries of PENNCNV_DATA_DIRECTORY and PENNCNV_RESULTS_DIRECTORY in the Genvisis project '.properties' file;"
                   + "   (2) Genvisis .sampRAF data, if 'custom.pbf' and 'custom.gcmodel' have not been generated;"
                   + "" + "cnv.analysis.DeNovoCNV requires 3 arguments\n"
                   + "   (1) The project's property file's name (i.e. proj="
                   + projPropertyFileFullPath + " (not the default));\n"
                   + "   (2) The output pedigree file name (i.e. pedigree=" + pedigreeFullPath
                   + " (not the default));\n"
                   + "   (3) Directory of PennCNV software (i.e. pennBinDir=" + pennBinDir
                   + " (not the default));\n"
                   + "   (4) The standard GC Base file name  (i.e. gcbase=" + gcBaseFileFullPath
                   + " (not the default));\n"
                   + "   (5) Number of window units for GC model (i.e. nwin=" + numWindowUnits
                   + " (not the default));\n"
                   + "   (6) To parse the PennCNV trio CNV and joint CNV results (i.e. parse="
                   + parseOrNot
                   + " (default)), after you submitted the batch jobs and get results.\n"
                   + "   (7) To generate list of De Novo CNV (i.e. parse=" + parseOrNot
                   + " (default)), after you submitted the batch jobs and get results.\n" + ""
                   + "Note:"
                   + "   (1) Population BAF is to be stored at the project's root directory with the name 'custom.pfb' (not changeable);\n"
                   + "   (2) Custom GC Model is to be stored at the project's root directory with the name 'custom.gcmodel' (not changeable);\n"
                   + "   (3) The PennCNV data directory is specified in the project's properties file (i.e. "
                   + pennDataDir + " (not the default));\n"
                   + "   (4) The PennCNV trios output directory is specified in the project's properties file (i.e. "
                   + pennOutputDir + " (not the default)).\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("proj=")) {
        projPropertyFileFullPath = arg.split("=")[1];
      } else if (arg.startsWith("pennBinDir=")) {
        pennBinDir = arg.split("=")[1];
      } else if (arg.startsWith("nwin=")) {
        numWindowUnits = Integer.parseInt(arg.split("=")[1]);
      } else if (arg.startsWith("gcbase=")) {
        gcBaseFileFullPath = arg.split("=")[1];
      } else if (arg.startsWith("pedigree=")) {
        pedigreeFullPath = arg.split("=")[1];
      } else if (arg.startsWith("penndatadir=")) {
        pennDataDir = arg.split("=")[1];
      } else if (arg.startsWith("pennoutdir=")) {
        pennOutputDir = arg.split("=")[1];
      } else if (arg.startsWith("outfile=")) {
        pennOutputFilenameRoot = arg.split("=")[1];
      } else if (arg.startsWith("parse=")) {
        parseOrNot = Boolean.parseBoolean(arg.split("=")[1]);
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }

    proj = new Project(projPropertyFileFullPath, false);
    log = proj.getLog();
    log.report(ext.getTime() + "\tstarting DeNovoCNV.java");
    pedigreeFullPath = proj.DATA_DIRECTORY.getValue(false, true) + "pedigree.dat";
    trioPedigreeFullPath = ext.parseDirectoryOfFile(pedigreeFullPath) + "pedigreeTrios.dat";
    gcModelFullPath = proj.PROJECT_DIRECTORY.getValue() + "custom.gcmodel";
    pennDataDir = proj.PENNCNV_DATA_DIRECTORY.getValue(false, true);
    pennOutputDir = proj.PENNCNV_RESULTS_DIRECTORY.getValue(false, true);


    if (!new File(trioPedigreeFullPath).exists()) {
      generateListOfTrios(pedigreeFullPath, trioPedigreeFullPath, log);
      verifyDataAvailabilityForListOfTrioSets(proj, trioPedigreeFullPath);
      generateFamFileForTrios(proj, trioPedigreeFullPath, true);
    } else {
      log.report("Skipped generating list of trio sets, - found " + trioPedigreeFullPath);
    }

    if (!new File(proj.PROJECT_DIRECTORY.getValue() + "custom.pfb").exists()) {
      PennCNV.populationBAF(proj);
    } else {
      log.report("Skipped generating custom population BAF, - found " + trioPedigreeFullPath);
    }

    if (!new File(gcModelFullPath).exists()) {
      PennCNV.gcModel(proj, gcBaseFileFullPath, gcModelFullPath, numWindowUnits);
    } else {
      log.report("Skipped generating custom GC Model, - found " + trioPedigreeFullPath);
    }

    if (!parseOrNot) {
      runPennCnv(proj, trioPedigreeFullPath, pennDataDir, gcModelFullPath, pennOutputDir,
                 pennOutputFilenameRoot, pennBinDir, -1, 30);
    } else {
      // batch(pedigreeFullPath);
      parsePennCnvResult(proj, pennOutputDir, trioPedigreeFullPath, ".jointcnv");
      parsePennCnvResult(proj, pennOutputDir, trioPedigreeFullPath, ".triocnv");

      // generateDenovoList(proj, proj.getDir(Project.DATA_DIRECTORY) + "denovo_joint.cnv", true,
      // log);

      // detectDenovoCnv(usage, ".triocnv", log);
      // detectDenovoCnv(usage, ".jointcnv", log);
    }

    log.report(ext.getTime() + "\tfinished DeNovoCNV.java.");
  }

  public static void parsePennCnvResult(Project proj, String pennCnvResultDir,
                                        String pennCnvResultFileNameExt) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    String temp, currentOffSpring;
    Vector<String> warnings = new Vector<String>();
    int[] position;
    String score;
    SampleData sampleData;
    String famIndPair;
    Hashtable<String, String> offspringCnv;
    String[] ids;
    String[] filenames;
    boolean offspringExisted;
    Logger log = proj.getLog();

    sampleData = proj.getSampleData(2, false);
    if (pennCnvResultFileNameExt.startsWith(".")) {
      pennCnvResultFileNameExt = pennCnvResultFileNameExt.substring(1);
    }
    try {
      filenames = Files.list(pennCnvResultDir, pennCnvResultFileNameExt, false);
      writer =
          new PrintWriter(new FileWriter(proj.DATA_DIRECTORY.getValue(false, true) + "denovo_"
                                         + pennCnvResultFileNameExt.replace("cnv", "") + ".cnv"));
      writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER));
      offspringCnv = new Hashtable<String, String>();
      for (String filename : filenames) {
        // currentOffspring = ext.rootOf(filenames[i]);
        offspringExisted = false;
        reader = new BufferedReader(new FileReader(pennCnvResultDir + filename));
        while (reader.ready()) {
          temp = reader.readLine();
          if (!temp.startsWith("NOTICE:")) {
            temp = PennCNV.translateDerivedSamples(temp, log);
            line = temp.trim().split("[\\s]+");
            if (line[7].equalsIgnoreCase("offspring")) {
              if (line[8].split("=")[1].startsWith("33")) {
                position = Positions.parseUCSClocation(line[0]);
                currentOffSpring = line[4];
                currentOffSpring =
                    currentOffSpring.substring(currentOffSpring.lastIndexOf("/") + 1);
                ids = sampleData.lookup(currentOffSpring);
                if (ids == null) {
                  // if (! offspringCnv.containsKey(trav)) {
                  // if (log != null) {
                  // log.reportError("Error - '" + trav + "' was not found in " +
                  // proj.getFilename(proj.SAMPLE_DATA_FILENAME));
                  // } else {
                  // System.out.println("Error - '" + trav + "' was not found in " +
                  // proj.getFilename(proj.SAMPLE_DATA_FILENAME));
                  // }
                  // }
                  famIndPair = currentOffSpring + "\t" + currentOffSpring;
                } else {
                  famIndPair = ids[1];
                }
                if (!offspringExisted) {
                  offspringCnv.put(currentOffSpring, famIndPair);
                  offspringExisted = true;
                }

                if (line.length < 8 || !line[7].startsWith("conf=")) {
                  score = "127";
                  if (!warnings.contains(currentOffSpring) && warnings.size() < 10) {
                    if (log != null) {
                      log.reportError("Warning - no conf estimates for " + currentOffSpring);
                    } else {
                      System.out.println("Warning - no conf estimates for " + currentOffSpring);
                    }
                    warnings.add(currentOffSpring);
                  }
                } else {
                  score = ext.formDeci(Double.parseDouble(line[7].substring(5)), 4, true);
                }
                writer.println(famIndPair + "\t" + position[0] + "\t" + position[1] + "\t"
                               + position[2] + "\t" + line[3].substring(line[3].indexOf("=") + 1)
                               + "\t" + score + "\t" + line[1].substring(7));
              }
            }
          }
        }
        reader.close();
      }
      writer.close();

      // FilterCalls.stdFilters(dir, ext.rootOf(filename)+".cnv", MAKE_UCSC_TRACKS);

    } catch (FileNotFoundException fnfe) {
      return;
    } catch (IOException ioe) {
      return;
    }
  }

  // /**
  // * Does not work. Still work in progress.
  // * @param proj
  // * @param cnvFileFullPath
  // * @param cnvFileHasHeader
  // * @param log
  // */
  // public static void generateDenovoList(Project proj, String cnvFileFullPath, boolean
  // cnvFileHasHeader) {
  // BufferedReader reader;
  // PrintWriter writer;
  // String[] line;
  // String temp, currentOffSpring;
  // int[] position;
  // String score;
  // SampleData sampleData;
  // String famIndPair;
  // Hashtable<String, String> offspringCnv;
  // String[] ids;
  // String[] filenames;
  // boolean offspringExisted;
  // Logger log = proj.getLog();
  // sampleData = proj.getSampleData(2, false);
  // try {
  // reader = new BufferedReader(new FileReader(cnvFileFullPath));
  // if (cnvFileHasHeader) {
  // reader.readLine();
  // }
  // writer = new PrintWriter(new FileWriter(ext.parseDirectoryOfFile(cnvFileFullPath) +
  // ext.rootOf(cnvFileFullPath) + "_list.txt"));
  // writer.write("DNA\tPosition\tComments");
  // while (reader.ready()) {
  // line = reader.readLine().split("\t");
  //// if(offspringCnv.containsKey(line[4])) {
  //// writer.println(offspringCnv.get(line[4]) + "\t" + sampleData.lookup(line[5])[2] + "\t" +
  // sampleData.lookup(line[6])[2] + "\t" + sampleData.getSexForIndividual(line[4]) + "\t1");
  //// writer.println(offspringCnv.get(line[4]) + "\t" + sampleData.lookup(line[5])[2] + "\t" +
  // sampleData.lookup(line[6])[2] + "\t" + sampleData.getSexForIndividual(line[4]) + "\t1");
  //// writer.println(offspringCnv.get(line[4]) + "\t" + sampleData.lookup(line[5])[2] + "\t" +
  // sampleData.lookup(line[6])[2] + "\t" + sampleData.getSexForIndividual(line[4]) + "\t1");
  //// }
  // }
  // writer.close();
  // reader.close();
  //
  // } catch (FileNotFoundException fnfe) {
  // return;
  // } catch (IOException ioe) {
  // return;
  // }
  // }

  public static void parsePennCnvResult(Project proj, String pennCnvResultDir,
                                        String triosPedigreeFullPath,
                                        String pennCnvResultFileNameExt) {
    BufferedReader reader;
    PrintWriter writer1, writer2;
    String[] line;
    String temp, currentOffspring, currentFather = null, currentMother = null;
    // Vector<String> warnings = new Vector<String>(); // TODO check to see if -conf is an option
    // for trio/joint calling
    int[] position;
    String score;
    SampleData sampleData;
    String famIndPair = null;
    Hashtable<String, String[]> triosPedigree;
    String[] ids;
    String[] filenames;
    String[] currentParents;
    Logger log = proj.getLog();

    sampleData = proj.getSampleData(2, false);
    if (pennCnvResultFileNameExt.startsWith(".")) {
      pennCnvResultFileNameExt = pennCnvResultFileNameExt.substring(1);
    }
    try {
      triosPedigree = new Hashtable<String, String[]>();
      reader = new BufferedReader(new FileReader(triosPedigreeFullPath));
      reader.readLine();
      while (reader.ready()) {
        line = reader.readLine().split("\t");
        triosPedigree.put(line[4], new String[] {line[5], line[6]});
      }

      filenames = Files.list(pennCnvResultDir, pennCnvResultFileNameExt, false);
      writer1 =
          new PrintWriter(new FileWriter(proj.DATA_DIRECTORY.getValue(false, true) + "denovo_"
                                         + pennCnvResultFileNameExt.replace("cnv", "") + ".cnv"));
      writer1.println(Array.toStr(CNVariant.PLINK_CNV_HEADER));
      writer2 = new PrintWriter(new FileWriter(proj.DATA_DIRECTORY.getValue(false, true) + "denovo_"
                                               + pennCnvResultFileNameExt.replace("cnv", "")
                                               + "_list.txt"));
      writer2.println("DNA\tPosition\tComments");
      for (String filename : filenames) {
        currentOffspring = null;
        reader = new BufferedReader(new FileReader(pennCnvResultDir + filename));
        while (reader.ready()) {
          temp = reader.readLine();
          if (!temp.startsWith("NOTICE:")) {
            temp = PennCNV.translateDerivedSamples(temp, log);
            line = temp.trim().split("[\\s]+");
            if (line[7].equalsIgnoreCase("offspring")) {
              if (line[8].split("=")[1].startsWith("33")) {
                position = Positions.parseUCSClocation(line[0]);
                if (currentOffspring == null) {
                  currentOffspring = line[4];
                  currentOffspring =
                      currentOffspring.substring(currentOffspring.lastIndexOf("/") + 1);
                  currentParents = triosPedigree.get(currentOffspring);
                  currentFather = currentParents[0];
                  currentMother = currentParents[1];
                  ids = sampleData.lookup(currentOffspring);
                  if (ids == null) {
                    // if (! offspringCnv.containsKey(trav)) {
                    // if (log != null) {
                    // log.reportError("Error - '" + trav + "' was not found in " +
                    // proj.getFilename(proj.SAMPLE_DATA_FILENAME));
                    // } else {
                    // System.out.println("Error - '" + trav + "' was not found in " +
                    // proj.getFilename(proj.SAMPLE_DATA_FILENAME));
                    // }
                    // }
                    famIndPair = currentOffspring + "\t" + currentOffspring;
                  } else {
                    famIndPair = ids[1];
                  }

                }
                if (line.length < 8 || !line[7].startsWith("conf=")) {
                  score = "127";
                  // if (!warnings.contains(trav) && warnings.size() < 10) {
                  // if (log != null) {
                  // log.reportError("Warning - no conf estimates for " + trav);
                  // } else {
                  // System.out.println("Warning - no conf estimates for " + trav);
                  // }
                  // warnings.add(trav);
                  // }
                } else {
                  score = ext.formDeci(Double.parseDouble(line[7].substring(5)), 4, true);
                }
                writer1.println(famIndPair + "\t" + position[0] + "\t" + position[1] + "\t"
                                + position[2] + "\t" + line[3].substring(line[3].indexOf("=") + 1)
                                + "\t" + score + "\t" + line[1].substring(7));
                writer2.println(currentOffspring + "\t" + line[0] + "\toffspring");
                writer2.println(currentFather + "\t" + line[0] + "\tfather");
                writer2.println(currentMother + "\t" + line[0] + "\tmother");
              }
            }
          }
        }
        reader.close();
      }
      writer1.close();
      writer2.close();
      // FilterCalls.stdFilters(dir, ext.rootOf(filename)+".cnv", MAKE_UCSC_TRACKS);
    } catch (FileNotFoundException fnfe) {
      return;
    } catch (IOException ioe) {
      return;
    }

    log.report(ext.getTime() + "\tparsed PennCNV results for " + pennCnvResultDir + "*"
               + pennCnvResultFileNameExt);
  }

  /**
   * Call PennCNV to run trios analysis
   * 
   * @param proj
   * @param listOfTrioSetsFullPath
   * @param pennDataDir
   * @param pennOutDir
   * @param pennOutFileNameRoot
   * @param pennBinDir
   */
  public static void runPennCnv(Project proj, String listOfTrioSetsFullPath, String pennDataDir,
                                String gcModelFullPath, String pennOutDir,
                                String pennOutFileNameRoot, String pennBinDir, int numQsubFiles,
                                int numCommandsPerQsubFile) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    String commands;
    Vector<String[]> iterationsVec;
    String[][] iterations;
    boolean gzip;
    Vector<String> jobNamesWithAbsolutePaths;
    String[] files;
    int counter;
    boolean found;
    boolean jarStatus;
    String pennDataFileExtension;
    String commandPrefix;
    String commandSuffix;
    Logger log = proj.getLog();

    if (!(new File(pennOutDir)).exists()) {
      new File(pennOutDir).mkdir();
    }
    // PennCNV.populationBAF(proj);
    // PennCNV.gcModel(proj, "D:/PennCNV_Related/GcModel/gc5Base_hg18.txt", "penn_output/gcmodel",
    // 0);

    // gzip = proj.getBoolean(proj.PENNCNV_GZIP_YESNO);
    gzip = proj.PENNCNV_GZIP_YESNO.getValue();
    if (gzip) {
      if (!(new File(pennOutDir + "lists/")).exists()) {
        new File(pennOutDir + "lists/").mkdir();
      }
      pennDataFileExtension = ".gz";
      commandPrefix = "`gunzip -c ";
      commandSuffix = ".gz`";
    } else {
      pennDataFileExtension = "";
      commandPrefix = "";
      commandSuffix = "";
    }
    jarStatus = proj.JAR_STATUS.getValue();

    iterationsVec = new Vector<String[]>();
    counter = 0;
    try {
      reader = new BufferedReader(new FileReader(listOfTrioSetsFullPath));
      reader.readLine();
      // ext.checkHeader(line, expected, kill)
      while (reader.ready()) {
        counter++;
        found = true;
        line = reader.readLine().trim().split("\t", -1);

        for (int i = 4; i <= 6; i++) {
          if (!(new File(pennDataDir + line[i] + pennDataFileExtension)).exists()) {
            if (Files.exists(proj.SAMPLE_DIRECTORY.getValue(false, true) + line[i]
                             + Sample.SAMPLE_FILE_EXTENSION, jarStatus)) {
              org.genvisis.cnv.analysis.AnalysisFormats.penncnv(proj, new String[] {line[i]}, null,
                                                                null, 1); // TODO How to generate
                                                                          // .gz format?
            } else {
              log.reportError("warning - skipped the following trio set due to no data avaiable (iDNA: "
                              + line[4] + "\tFaDNA: " + line[5] + "\tMoDNA:" + line[6] + ")");
              found = false;
              break;
            }
          }
        }

        if (found) {
          iterationsVec.add(new String[] {line[5], line[6], line[4]});

          writer =
              new PrintWriter(new FileOutputStream(pennOutDir + "lists/" + line[4] + "_solo.lst"));
          writer.println(commandPrefix + pennDataDir + line[5] + commandSuffix + "\n"
                         + commandPrefix + pennDataDir + line[6] + commandSuffix + "\n"
                         + commandPrefix + pennDataDir + line[4] + commandSuffix);
          writer.close();

          writer =
              new PrintWriter(new FileOutputStream(pennOutDir + "lists/" + line[4] + "_trio.lst"));
          writer.println(commandPrefix + pennDataDir + line[5] + commandSuffix + "\t"
                         + commandPrefix + pennDataDir + line[6] + commandSuffix + "\t"
                         + commandPrefix + pennDataDir + line[4] + commandSuffix);
          writer.close();
        }
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      log.reportError("Error: file \"" + listOfTrioSetsFullPath
                      + "\" not found in current directory");
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + listOfTrioSetsFullPath + "\"");
    }
    log.report("Total number of trio sets: " + counter + "\t(" + listOfTrioSetsFullPath + ")");

    iterations = new String[iterationsVec.size()][3];
    for (int i = 0; i < iterations.length; i++) {
      iterations[i] = iterationsVec.elementAt(i);
    }

    if (numQsubFiles >= 0) {
      if (numCommandsPerQsubFile >= 0) {
        log.reportError("Warning: cannot specify numQsubFiles and numCommandsPerQsubFile at the same time. Program continues by ignoring numCommandsPerQsubFile.");
      }
    } else if (numCommandsPerQsubFile >= 0) {
      numQsubFiles = (int) Math.ceil((double) iterations.length * 3 / numCommandsPerQsubFile);
    } else {
      log.reportError("Error: must specify either numQsubFiles or numCommandsPerQsubFile. Program exits.");
      return;
    }

    if (Files.isWindows()) {
      log.report(ext.getTime() + "\tstarting PennCNV for CNV detection.");
      for (String[] iteration : iterations) {
        CmdLine.run("perl " + pennBinDir + "detect_cnv.pl -test -hmm " + pennBinDir
                    + "lib/hh550.hmm -pfb " + proj.PROJECT_DIRECTORY.getValue()
                    + "custom.pfb -gcmodel " + gcModelFullPath + " " + pennDataDir + iteration[0]
                    + " " + pennDataDir + iteration[1] + " " + pennDataDir + iteration[2] + " -out "
                    + pennOutDir + iteration[2] + ".rawcnv", pennOutDir);
        CmdLine.run("perl " + pennBinDir + "detect_cnv.pl -trio -hmm " + pennBinDir
                    + "lib/hh550.hmm -pfb " + proj.PROJECT_DIRECTORY.getValue()
                    + "custom.pfb -gcmodel " + gcModelFullPath + " -cnv " + pennOutDir
                    + iteration[0] + ".rawcnv " + pennDataDir + iteration[0] + " " + pennDataDir
                    + iteration[1] + " " + pennDataDir + iteration[2] + " -out " + pennOutDir
                    + iteration[2] + ".triocnv", pennOutDir);
        CmdLine.run("perl " + pennBinDir + "detect_cnv.pl -joint -hmm " + pennBinDir
                    + "lib/hh550.hmm -pfb " + proj.PROJECT_DIRECTORY.getValue()
                    + "custom.pfb -gcmodel " + gcModelFullPath + " " + pennDataDir + iteration[0]
                    + " " + pennDataDir + iteration[1] + " " + pennDataDir + iteration[2] + " -out "
                    + pennOutDir + iteration[2] + ".jointcnv", pennOutDir);
      }

      // (new File(cnvOutputDir+"resultAll.jointcnv")).renameTo(new File(cnvOutputDir +
      // "resultAll_tmp.jointcnv"));
      // common.Files.cat(new String[] {cnvOutputDir+line[4]+".jointcnv",
      // cnvOutputDir+"resultAll_tmp.jointcnv"}, cnvOutputDir+"resultAll.jointcnv", null, null);
      // (new File(cnvOutputDir + line[4] + ".jointcnv")).delete();
      // if (!(new File("penn_output/resultAll_tmp.jointcnv")).delete()) {
      // System.err.println("Error deleting the temporaryfile 'resultAll_tmp.jointcnv'.");
      // }

      log.report(ext.getTime() + "\tfinished running PennCNV. Output location: " + pennOutDir);
    } else {
      // Runtime.getRuntime().exec("/home/pankrat2/shared/penncnv/detect_cnv.pl -joint -hmm
      // /home/pankrat2/shared/penncnv/lib/hh550.hmm -pfb " + proj.getProjectDir() + "custom.pfb
      // -gcmodel " + gcBaseFileFullPath + " " + pennDataDir + line[4] + " " + pennDataDir + line[5]
      // + " " + pennDataDir + line[6] + " -out " + line[4] + ".jointcnv");
      // Runtime.getRuntime().exec("/home/pankrat2/shared/penncnv/detect_cnv.pl -trio -hmm
      // /home/pankrat2/shared/penncnv/lib/hh550.hmm -pfb " + proj.getProjectDir() + "custom.pfb
      // -cnv " + pennOutDir + pennOutFileNameRoot + " " + pennDataDir + line[4] + " " + pennDataDir
      // + line[5] + " " + pennDataDir + line[6] + " -out " + line[4] + ".triocnv");

      commands = pennBinDir + "detect_cnv.pl -test -hmm " + pennBinDir + "lib/hh550.hmm -pfb "
                 + proj.PROJECT_DIRECTORY.getValue() + "custom.pfb -gcmodel " + gcModelFullPath
                 + " -list " + pennOutDir + "lists/[%2]_solo.lst -out " + pennOutDir
                 + "[%2].rawcnv\n" + pennBinDir + "detect_cnv.pl -trio -hmm " + pennBinDir
                 + "lib/hh550.hmm -pfb " + proj.PROJECT_DIRECTORY.getValue()
                 + "custom.pfb -gcmodel " + gcModelFullPath + " -cnv " + pennDataDir
                 + "[%2].rawcnv -list " + pennOutDir + "/lists/[%2]_trio.lst -out " + pennOutDir
                 + "[%2].triocnv\n" + pennBinDir + "detect_cnv.pl -joint -hmm " + pennBinDir
                 + "lib/hh550.hmm -pfb " + proj.PROJECT_DIRECTORY.getValue()
                 + "custom.pfb -gcmodel " + gcModelFullPath + " -list " + pennOutDir
                 + "lists/[%2]_trio.lst -out " + pennOutDir + "[%2].jointcnv\n";

      if (!(new File(pennOutDir + "scripts/")).exists()) {
        new File(pennOutDir + "scripts/").mkdir();
      }
      // Files.qsub("runPennCNV", proj.getProjectDir() + "penn_scripts/", numQsubFiles, commands,
      // iterations, 2000, 24);
      // Files.qsubMultiple("chunkCNV", Array.stringArraySequence(numQsubFiles, proj.getProjectDir()
      // + "penn_scripts/runPennCNV_", ".qsub"), 8, -1, 22000, 24);
      Files.qsub(pennOutDir + "scripts/runPennCNV", pennOutDir + "scripts/", numQsubFiles, commands,
                 iterations, 2000, 24);
      // Files.qsubMultiple("chunkCNV", Array.stringArraySequence(numQsubFiles, pennOutDir +
      // "scripts/runPennCNV_", ".qsub"), 2000, 24);

      files = Files.list(pennOutDir + "scripts/", "runPennCNV", ".qsub", false, false);
      jobNamesWithAbsolutePaths = new Vector<String>(files.length);
      for (String file : files) {
        jobNamesWithAbsolutePaths.add(pennOutDir + "scripts/" + file);
      }
      Files.qsubMultiple(jobNamesWithAbsolutePaths, null, pennOutDir + "scripts/", "chunkCNV", 8,
                         true, null, -1, 22000, 24);

      log.report(ext.getTime()
                 + "\tfinished generating the batch scripts to run PennCNV for trio CNV and joint CNV detection.\n                 output location: "
                 + pennOutDir
                 + "scripts/\nYour next step after this program: submit/run the corresponding scripts.");
    }
  }

  public static void verifyDataAvailabilityForListOfTrioSets(Project proj,
                                                             String originalListOfTrioSetsFullPath) {
    BufferedReader reader;
    String line1;
    String[] line;
    // String sampDataReady;
    String noSampData;
    PrintWriter writer;
    String[] sampleList;
    boolean[] found;
    String currentId;
    int counter;
    String originalDirAndRoot;
    String originalExt;
    String filename;
    Logger log = proj.getLog();

    originalDirAndRoot = ext.parseDirectoryOfFile(originalListOfTrioSetsFullPath)
                         + ext.rootOf(originalListOfTrioSetsFullPath);
    counter = originalListOfTrioSetsFullPath.lastIndexOf(".");
    if (counter == -1 || counter < originalListOfTrioSetsFullPath.lastIndexOf("/")) {
      originalExt = "";
    } else {
      originalExt = originalListOfTrioSetsFullPath.substring(counter);
    }
    sampleList = proj.getSampleList().getSamples();
    counter = 0;
    // sampDataReady = "";
    noSampData = "";
    try {
      // reader = Files.getReader(new FileReader(filename), proj.getJarStatus(), true, false);
      reader = new BufferedReader(new FileReader(originalListOfTrioSetsFullPath));
      reader.readLine();
      filename = originalDirAndRoot + "_Tmp" + originalExt;
      writer = new PrintWriter(new FileWriter(filename));
      writer.println("fId\tiId\tfaId\tmoId\tiDna\tfaDna\tmoDna");
      while (reader.ready()) {
        line1 = reader.readLine();
        line = line1.split("\t", -1);
        found = new boolean[3];
        for (int i = 0; i < 3; i++) {
          currentId = line[4 + i];
          for (String element : sampleList) {
            if (element.equals(currentId)) {
              found[i] = true;
              break;
            }
          }
        }
        if (found[0] && found[1] && found[2]) {
          writer.println(line1);
        } else {
          counter++;
          noSampData += (line1 + "\n");
        }
      }
      reader.close();
      writer.flush();
      writer.close();

      if (counter > 0) {
        new File(originalListOfTrioSetsFullPath).renameTo(new File(originalListOfTrioSetsFullPath
                                                                   + "_wDataAndNoData"
                                                                   + originalExt));
        new File(filename).renameTo(new File(originalListOfTrioSetsFullPath));
        filename = originalDirAndRoot + "_NoDataOnly" + originalExt;
        writer = new PrintWriter(new FileWriter(filename));
        writer.println("fId\tiId\tfaId\tmoId\tiDna\tfaDna\tmoDna");
        writer.print(noSampData);
        writer.flush();
        writer.close();
        log.report("There are " + counter
                   + " set(s) of trio(s) are removed from the list due to no corresponding "
                   + Sample.SAMPLE_FILE_EXTENSION + " file(s).\nCheck the following for detail: "
                   + filename);
      } else {
        new File(filename).delete();
      }
    } catch (FileNotFoundException fnfe) {
      log.reportError("Error: file \"" + originalListOfTrioSetsFullPath
                      + "\" not found in current directory");
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + originalListOfTrioSetsFullPath + "\"");
    }
    log.report(ext.getTime() + "\tverified data availability for the list of trio sets "
               + originalListOfTrioSetsFullPath);
  }

}
