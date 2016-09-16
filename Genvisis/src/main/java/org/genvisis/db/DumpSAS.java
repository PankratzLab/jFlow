package org.genvisis.db;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.Array;
import org.genvisis.common.CmdLine;
import org.genvisis.common.CountHash;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.Sort;
import org.genvisis.common.ext;
import org.genvisis.parse.LookupTable;

public class DumpSAS {
  // public static final String DEFAULT_SAS_EXECUTABLE =
  // "/share/apps/src/sas/SASFoundation/9.2/bin/sas_en";
  public static final String DEFAULT_SAS_EXECUTABLE = "sas";
  public static final String[] CONTENTS_HEADER =
                                               {"LIBNAME", "MEMNAME", "MEMLABEL", "TYPEMEM", "NAME",
                                                "TYPE", "LENGTH", "VARNUM", "LABEL", "FORMAT",
                                                "FORMATL", "FORMATD", "INFORMAT", "INFORML",
                                                "INFORMD", "JUST", "NPOS", "NOBS", "ENGINE",
                                                "CRDATE", "MODATE", "DELOBS", "IDXUSAGE", "MEMTYPE",
                                                "IDXCOUNT", "PROTECT", "FLAGS", "COMPRESS", "REUSE",
                                                "SORTED", "SORTEDBY", "CHARSET", "COLLATE",
                                                "NODUPKEY", "NODUPREC", "ENCRYPT", "POINTOBS",
                                                "GENMAX", "GENNUM", "GENNEXT"};
  public static final String[] ALL_CONTENTS_HEADER = {"MEMNAME", "NAME", "LABEL", "TYPE", "LENGTH",
                                                      "VARNUM"};

  public static void dump(String dir, String sasExecutable, String password) {
    PrintWriter writer;
    String[] files;

    dir = ext.verifyDirFormat(dir);
    files = Files.list(dir, ".sas7bdat", false);
    dir = dir.substring(0, dir.length() - 1);

    System.out.println("Found " + files.length + " with the extension .sas7bdat to convert");
    try {
      writer = new PrintWriter(new FileWriter("dumpAll.sas"));

      writer.println("options nofmterr;");
      // writer.println("options PS=1000000;"); // if printing to .lst file and you don't want page
      // breaks

      writer.println("");
      writer.println("libname a '" + dir + "';");
      writer.println("");
      for (String file : files) {
        System.out.println(file);
        if (!Files.isWindows()) {
          new File(file).renameTo(new File(file.toLowerCase()));
        }
        writer.println("DATA " + ext.rootOf(file) + "; SET a." + ext.rootOf(file)
                       + (password != null ? " (pw=\"" + password + "\")" : "") + "; RUN;");

        writer.println("PROC EXPORT DATA=" + ext.rootOf(file));
        writer.println("OUTFILE=\"" + ext.rootOf(file) + ".xln\"");
        writer.println("DBMS=TAB REPLACE;");
        writer.println("RUN;");
        writer.println("");

        writer.println("PROC CONTENTS DATA=" + ext.rootOf(file) + " out=contents;");
        writer.println("RUN;");
        writer.println("");

        writer.println("PROC EXPORT DATA=contents");
        writer.println("OUTFILE=\"" + ext.rootOf(file) + "_contents.xln\"");
        writer.println("DBMS=TAB REPLACE;");
        writer.println("RUN;");

        // writer.println("PROC EXPORT DATA=a."+ext.rootOf(files[i]));
        // writer.println("OUTFILE=\""+ext.rootOf(files[i])+".xln\"");
        // writer.println("DBMS=TAB REPLACE;");
        // writer.println("RUN;");
        // writer.println("");
        //
        // writer.println("PROC CONTENTS DATA=a."+ext.rootOf(files[i])+" out=a.contents;");
        // writer.println("RUN;");
        // writer.println("");
        //
        // writer.println("PROC EXPORT DATA=a.contents");
        // writer.println("OUTFILE=\""+ext.rootOf(files[i])+"_contents.xln\"");
        // writer.println("DBMS=TAB REPLACE;");
        // writer.println("RUN;");
        writer.println("");
      }

      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + "dumpAll.sas");
      e.printStackTrace();
    }

    Files.write(sasExecutable + " -noterminal " + "dumpAll.sas", "batchForDump");
    Files.chmod("batchForDump");

    CmdLine.run("./batchForDump", dir);
    new File("contents.sas7bdat").delete();
  }

  public static void tallyAll(String dir, String filename, String tableFormat) {
    BufferedReader reader;
    PrintWriter iWriter, writer;
    String[] line;
    Hashtable<String, String> typedWhites, typedBlacks, descriptions, validBlackFemales,
        validBlackMales, validWhiteFemales, validWhiteMales;
    Hashtable<String, int[]> cells;
    Hashtable<String, Vector<String>> hash, tableOrder;
    Vector<String> v = new Vector<String>();
    String[] dbs, header, vars;
    Logger log;
    int[] indices;
    int[][] counts;
    int[] cellCounts;
    String[] keys;
    CountHash ch;

    log = new Logger(dir + ext.rootOf(filename) + ".log");
    hash = HashVec.loadFileToHashVec(dir + filename, 0, new int[] {1}, "\t", false, false);
    descriptions = HashVec.loadFileToHashString(dir + filename, new int[] {0, 1}, new int[] {2},
                                                false, "\t", false, false, false);
    dbs = HashVec.getKeys(hash);
    cells = new Hashtable<String, int[]>();
    validBlackFemales = new Hashtable<String, String>();
    validBlackMales = new Hashtable<String, String>();
    validWhiteFemales = new Hashtable<String, String>();
    validWhiteMales = new Hashtable<String, String>();
    ch = new CountHash();
    try {
      writer = new PrintWriter(new FileWriter(dir + ext.rootOf(filename) + ".out"));
      for (String db : dbs) {
        writer.println(db);
        header = Files.getHeaderOfFile(dir + db.toLowerCase() + ".xln", "\t", log);
        v = hash.get(db);
        vars = Array.toStringArray(v);
        v.insertElementAt("CAReID", 0);
        // v.insertElementAt("ID", 0);
        indices = ext.indexFactors(Array.toStringArray(v), header, false, true);
        writer.println("\t" + Array.toStr(Array.subArray(header, indices)));
        writer.print("\t");
        for (int j = 1; j < indices.length; j++) {
          writer.print("\t" + descriptions.get(db + "\t" + vars[j - 1]));
        }
        writer.println();
        try {
          iWriter = new PrintWriter(new FileWriter(dir + db + "_subset.xln"));
          try {
            reader = new BufferedReader(new FileReader(dir + db.toLowerCase() + ".xln"));
            if (Files.exists(dir + db.substring(0, db.indexOf("_")) + "_Whites" + ".fam", false)) {
              typedWhites = HashVec.loadFileToHashString(dir + db.substring(0, db.indexOf("_"))
                                                         + "_Whites" + ".fam", 1, new int[] {4},
                                                         "\t", false);
            } else {
              typedWhites = new Hashtable<String, String>();
            }
            if (Files.exists(dir + db.substring(0, db.indexOf("_")) + "_Blacks" + ".fam", false)) {
              typedBlacks = HashVec.loadFileToHashString(dir + db.substring(0, db.indexOf("_"))
                                                         + "_Blacks" + ".fam", 1, new int[] {4},
                                                         "\t", false);
            } else {
              typedBlacks = new Hashtable<String, String>();
            }
            counts = new int[vars.length + 1][3];
            while (reader.ready()) {
              line = reader.readLine().split("\t", -1);
              line = Array.subArray(line, indices);
              iWriter.println(Array.toStr(line));
              for (int j = 0; j < indices.length; j++) {
                if (!line[j].equals("")) {
                  counts[j][0]++;
                  if (typedWhites.containsKey(line[0])) {
                    counts[j][1]++;
                    if (typedWhites.get(line[0]).equals("2")) {
                      validWhiteFemales.put(line[0], "");
                    } else if (typedWhites.get(line[0]).equals("1")) {
                      validWhiteMales.put(line[0], "");
                    } else {
                      System.err.println(line[0] + "\t" + typedWhites.get(line[0]) + "\t" + db
                                         + "\twhite");
                      ch.add(typedWhites.get(line[0]) + "\t" + db + "\twhite");
                    }
                  }
                  if (typedBlacks.containsKey(line[0])) {
                    counts[j][2]++;
                    if (typedBlacks.get(line[0]).equals("2")) {
                      validBlackFemales.put(line[0], "");
                    } else if (typedBlacks.get(line[0]).equals("1")) {
                      validBlackMales.put(line[0], "");
                    } else {
                      System.err.println(line[0] + "\t" + typedBlacks.get(line[0]) + "\t" + db
                                         + "\tblack");
                      ch.add(typedWhites.get(line[0]) + "\t" + db + "\tblack");
                    }
                  }
                }
              }
            }
            writer.print("N");
            for (int j = 0; j < indices.length; j++) {
              writer.print("\t" + counts[j][0]);
            }
            writer.println();
            writer.print("N+genoWhites");
            for (int j = 0; j < indices.length; j++) {
              writer.print("\t" + counts[j][1]);
            }
            writer.println();
            writer.print("N+genoBlacks");
            for (int j = 0; j < indices.length; j++) {
              writer.print("\t" + counts[j][2]);
            }
            writer.println();

            for (int j = 1; j < indices.length; j++) {
              cellCounts = new int[] {counts[0][0], counts[0][1], counts[0][2], counts[j][0],
                                      counts[j][1], counts[j][2]};
              cells.put(db + "\t" + v.elementAt(j), cellCounts);
            }

            reader.close();
          } catch (FileNotFoundException fnfe) {
            System.err.println("Error: file \"" + dir + db + ".xln"
                               + "\" not found in current directory");
            System.exit(1);
          } catch (IOException ioe) {
            System.err.println("Error reading file \"" + dir + db + ".xln" + "\"");
            System.exit(2);
          }
          iWriter.close();
        } catch (Exception e) {
          System.err.println("Error writing to " + dir + db + "_subset.xln");
          e.printStackTrace();
        }
        writer.println();
        writer.println();
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + dir + ext.rootOf(filename) + ".out");
      e.printStackTrace();
    }

    System.out.println("Total valid black females: " + validBlackFemales.size());
    System.out.println("Total valid black males: " + validBlackMales.size());
    System.out.println("Total valid white females: " + validWhiteFemales.size());
    System.out.println("Total valid white males: " + validWhiteMales.size());
    keys = ch.getValues();
    for (String key : keys) {
      System.out.println(key + " " + ch.getCount(key));
    }


    try {
      writer = new PrintWriter(new FileWriter(dir + ext.rootOf(filename) + "_tables.xln"));
      tableOrder = HashVec.loadFileToHashVec(dir + tableFormat, new int[] {0}, new int[] {1, 2},
                                             "\t", false, false);
      descriptions = HashVec.loadFileToHashString(dir + tableFormat, new int[] {1, 2},
                                                  new int[] {3}, false, "\t", false, false, false);
      keys = HashVec.getKeys(tableOrder);
      for (String key : keys) {
        v = tableOrder.get(key);
        counts = new int[v.size()][4];
        writer.print(key);
        for (int j = 0; j < counts.length; j++) {
          writer.print("\t" + descriptions.get(v.elementAt(j)));
          counts[j] = cells.get(v.elementAt(j));
        }
        writer.println();
        for (int k = 0; k < 8; k++) {
          switch (k) {
            case 0:
              writer.print("Total in cohort");
              break;
            case 1:
              writer.print("N Whites with IBC");
              break;
            case 2:
              writer.print("N Blacks with IBC");
              break;
            case 3:
              writer.print("Total with Pheno");
              break;
            case 4:
              writer.print("N Whites with both Pheno and IBC");
              break;
            case 5:
              writer.print("N Blacks with both Pheno and IBC");
              break;
            case 6:
              writer.print("N Whites need to phenotype");
              break;
            case 7:
              writer.print("N Blacks need to phenotype");
              break;
          }
          for (int[] count : counts) {
            if (k >= 6) {
              writer.print("\t" + (count[k - 5] - count[k - 2]));
            } else {
              writer.print("\t" + count[k]);
            }
          }
          writer.println();
        }
        writer.println();

      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + dir + ext.rootOf(filename) + "_tables.xln");
      e.printStackTrace();
    }
  }

  public static void procAllContents(String dir) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    Hashtable<String, String> hash;
    String[] files;

    files = new File(dir).list(new FilenameFilter() {
      @Override
      public boolean accept(File file, String filename) {
        return filename.endsWith("_contents.xln") && !filename.equals("all_contents.xln");
      }
    });

    try {
      writer = new PrintWriter(new FileWriter(dir + "all_contents.xln"));
      writer.println(Array.toStr(ALL_CONTENTS_HEADER));
      for (String file : files) {
        System.out.println("Parsing " + file);
        try {
          reader = new BufferedReader(new FileReader(dir + file));
          ext.checkHeader(reader.readLine().trim().split("\t", -1), CONTENTS_HEADER, false);
          hash = new Hashtable<String, String>();
          while (reader.ready()) {
            line = ext.removeQuotesFromExcelToken(reader.readLine(), new String[][] {{"\t", ""}})
                      .split("\t", -1);
            hash.put(line[7], line[1] + "\t" + line[4] + "\t" + line[8] + "\t" + line[5] + "\t"
                              + line[6] + "\t" + line[7]);
          }
          reader.close();
          for (int j = 1; j <= hash.size(); j++) {
            writer.println(hash.get(j + ""));
          }
        } catch (FileNotFoundException fnfe) {
          System.err.println("Error: file \"" + dir + file + "\" not found in current directory");
          System.exit(1);
        } catch (IOException ioe) {
          System.err.println("Error reading file \"" + dir + file + "\"");
          System.exit(2);
        }
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + dir + "all_contents.xln");
      e.printStackTrace();
    }
  }

  public static void catAllAlls(String dir, String manifest) {
    BufferedReader reader;
    PrintWriter writer;
    String[][] mani;

    mani = HashVec.loadFileToStringMatrix(dir + manifest, false, new int[] {0, 1}, false);
    try {
      writer = new PrintWriter(new FileWriter(dir + "master_contents.xln"));
      writer.println("Group\t" + Array.toStr(ALL_CONTENTS_HEADER));
      for (String[] element : mani) {
        try {
          reader = new BufferedReader(new FileReader(dir + element[1] + "all_contents.xln"));
          ext.checkHeader(reader.readLine().trim().split("[\\s]+"), ALL_CONTENTS_HEADER, true);
          if (!reader.ready()) {
            System.err.println("'" + element[1] + "all_contents.xln" + "' contianed no elements");
          }
          while (reader.ready()) {
            writer.println(element[0] + "\t" + reader.readLine());
          }
          reader.close();
        } catch (FileNotFoundException fnfe) {
          System.err.println("Error: file \"" + dir + element[1] + "all_contents.xln"
                             + "\" not found in current directory");
          System.exit(1);
        } catch (IOException ioe) {
          System.err.println("Error reading file \"" + dir + element[1] + "all_contents.xln"
                             + "\"");
          System.exit(2);
        }
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + dir + "master_contents.xln");
      e.printStackTrace();
    }
  }

  public static void mergeFromParameters(String filename, Logger log) {
    PrintWriter writer;
    String[] line;
    Hashtable<String, Hashtable<String, String>> hashes;
    Hashtable<String, String> hash, forcedIDs, masterIDs;
    Vector<String> params, files, ids;
    boolean error, foundAnID;
    String file;
    String[] keys, header;
    int[] indices;
    String crffile, outfile, idfile;
    String rootDir, idName;

    params = Files.parseControlFile(filename, "sas",
                                    new String[] {"ID CARDIA_ID otherPossibleNamesForTheSameID ids=id_superset.dat crf=parse.crf out=output.xln rootDir=N:/cardia/phenotypes/",
                                                  "#SUBDIRECTORY (blank tab if all in root directory)\tFILENAME\tVARIABLE\tVARIABLE DESCRIPTION(optional)\tForce_which_ID_to_use(optional)",
                                                  "Y15\tF1F05\tF05CTB1B\t10ML LAVENDER TB, 2ML PLASMA-ISOPROS\tCARDIA_ID",
                                                  "Y15\tF1ISOP\tFLFISOP\tISOPROSTANE PG/ML",
                                                  "Y20_CORE\tG3F05\tG05YTB1A\tRED STRIPE LABEL, 10ML LAVENDER TB, 1.6ML PLASMA-ISOPROSTANE"},
                                    log);
    if (params != null) {
      line = params.remove(0).trim().split("[\\s]+");
      crffile = ext.rootOf(filename) + "_parse.crf";
      outfile = ext.rootOf(filename) + "_out.xln";
      idfile = ext.rootOf(filename) + "_ids.dat";
      ids = new Vector<String>();
      rootDir = ext.parseDirectoryOfFile(filename);
      for (String element : line) {
        if (element.startsWith("crf=")) {
          crffile = ext.parseStringArg(element, ext.rootOf(filename) + "_parse.crf");
        } else if (element.startsWith("out=")) {
          outfile = ext.parseStringArg(element, ext.rootOf(filename) + "_out.xln");
        } else if (element.startsWith("ids=")) {
          idfile = ext.parseStringArg(element, ext.rootOf(filename) + "_ids.dat");
        } else if (element.startsWith("rootDir=")) {
          rootDir = ext.parseStringArg(element, "./");
        } else {
          ids.add(element);
        }
      }
      error = false;
      forcedIDs = new Hashtable<String, String>();
      hashes = new Hashtable<String, Hashtable<String, String>>();
      files = new Vector<String>();
      for (int i = 0; i < params.size(); i++) {
        line = params.elementAt(i).trim().split("\t");
        if (line[0].equals("")) {
          line[0] = ".";
        }
        line[0] = ext.verifyDirFormat(line[0]);
        if (!new File(rootDir + line[0]).exists() || !new File(rootDir + line[0]).isDirectory()) {
          log.reportError("Error - could not find subdirectory '" + line[0] + "' within '" + rootDir
                          + "'");
          error = true;
        }
        file = line[0] + line[1] + ".xln";
        if (!new File(rootDir + file).exists()) {
          file = line[0] + "dump/" + line[1] + ".xln";
          if (!new File(rootDir + file).exists()) {
            log.reportError("Error - could not find file '" + line[1] + ".xln" + "' in '" + rootDir
                            + line[0] + "' or '" + rootDir + line[0] + "dump/'");
            error = true;
          }
        }
        HashVec.addIfAbsent(file, files);
        HashVec.addToHashHash(hashes, file, line[2], line.length > 3 ? line[3] : "");
        if (line.length > 4) {
          forcedIDs.put(file, line[4]);
        }
      }
      try {
        writer = new PrintWriter(new FileWriter(crffile));
        writer.println("lookup");
        writer.println(idfile + " head=IID out=" + outfile);
        masterIDs = new Hashtable<String, String>();
        for (int i = 0; i < files.size(); i++) {
          file = files.elementAt(i);
          header = Files.getHeaderOfFile(rootDir + file, log);
          writer.print(rootDir + file);
          if (forcedIDs.containsKey(file)) {
            idName = forcedIDs.get(file);
            writer.print(" '" + idName + "'");
            if (idName.equals("")) {
              log.reportError("Error - ID was set to an empty TAB, there must be an extra or a trailing tab for file '"
                              + file + "'");
              error = true;
            } else if (ext.indexOfStr(idName, header) == -1) {
              log.reportError("\nError - Since there were more than 4 columns for file '" + file
                              + "', the algorithm assumes that you are specifying a specific ID to use (in this case '"
                              + idName
                              + "'). However, there was no such column header in the file.");
              error = true;
            }
            indices = ext.indexFactors(Array.toStringArray(ids), header, true, log, false, false);
          } else {
            foundAnID = false;
            // TODO If case insensitivity is necessary, make sure all downstream calls can and do
            // ignore case (including GenParser)
            indices = ext.indexFactors(Array.toStringArray(ids), header, true, log, false, false);
            for (int j = 0; j < ids.size(); j++) {
              if (indices[j] != -1) {
                if (foundAnID) {
                  log.reportError("More than one \"valid\" ID available for file '" + file
                                  + "'; only the first will be used");
                } else {
                  writer.print(" '" + ids.elementAt(j) + "'");
                  foundAnID = true;
                  keys = Array.toStringArray(HashVec.loadFileToVec(rootDir + file, true,
                                                                   new int[] {indices[j]}, false,
                                                                   false, false, "\t"));
                  for (String key : keys) {
                    masterIDs.put(key, "");
                  }
                }
              }
            }
            if (!foundAnID) {
              log.reportError("Error - could not find a valid ID to match on within '" + file
                              + "'");
              error = true;
            }
          }
          writer.print(" tab");

          hash = hashes.get(file);
          keys = HashVec.getKeys(hash);
          indices = Sort.putInOrder(ext.indexFactors(keys, header, true, false));
          for (int indice : indices) {
            writer.print(" '" + header[indice] + "'=" + header[indice]
                         + (hash.get(header[indice]).equals("") ? ""
                                                                : "_"
                                                                  + ext.replaceAllWith(hash.get(header[indice]),
                                                                                       " ", "_")));
          }
          writer.println();
        }
        writer.close();
        Files.writeList(HashVec.getKeys(masterIDs, true, false), idfile);
        if (error) {
          log.reportError("\nFailed to generate all bits of " + crffile
                          + "; the hits algorithm was not attempted");
          return;
        }
        LookupTable.fromParameters(crffile, new Logger(ext.rootOf(crffile) + ".log"));
      } catch (Exception e) {
        System.err.println("Error writing to " + crffile);
        e.printStackTrace();
      }
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String dir = "./";
    String sasExecutable = DEFAULT_SAS_EXECUTABLE;
    String password = null;

    String usage = "\n" + "db.DumpSAS requires 0-1 arguments\n"
                   + "   (1) directory containing .sas7bdat files to dump (i.e. dir=" + dir
                   + " (default))\n" + "   (2) path to SAS executable (i.e. executable="
                   + sasExecutable + " (default))\n"
                   + "   (3) (optional) password, if one is required for datasets (i.e. password=whatever (not the default))\n"
                   + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("dir=")) {
        dir = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("executable=")) {
        sasExecutable = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("password=")) {
        password = arg.split("=")[1];
        numArgs--;
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }

    try {
      if (Files.isWindows()) {
        // tallyAll("D:/xCARe_transfer/sas/all2/", "tallyho.txt", "tableFormat.txt");
        // procAllContents("N:/cardia/phenotypes/Y0/dump/");
        // procAllContents("N:/cardia/phenotypes/Y2/dump/");
        // procAllContents("N:/cardia/phenotypes/Y5/dump/");
        // procAllContents("N:/cardia/phenotypes/Y7/dump/");
        // procAllContents("N:/cardia/phenotypes/Y10/dump/");
        // procAllContents("N:/cardia/phenotypes/Y15/dump/");
        // procAllContents("N:/cardia/phenotypes/Y20/CORE/dump/");
        // procAllContents("N:/cardia/phenotypes/Y20/CFS/dump/");
        // procAllContents("N:/cardia/phenotypes/Year 25 Dec 1 2011/dump/");
        // catAllAlls("N:/cardia/phenotypes/", "manifest.txt");
        // tallyAll("D:/Myron/CARDIA/ICAM1/", "tallyho.txt", "tableFormat.txt");
        // procAllContents("D:/PAGE/SOL/00src/dumps/");
        // procAllContents("D:/PAGE/SOL/00src/dumps/test/");
        mergeFromParameters("parseBloodCellTraits.crf", new Logger("parseBloodCellTraits_.log"));
      } else {
        dump(dir, sasExecutable, password);
        procAllContents("./");
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

}
