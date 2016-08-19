// File order: DNA FamID IndID

package org.genvisis.park;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.ext;

public class CheckIDsAgainstDNAs {
  // public static final String DNA_YEAR_BUG_SOURCE = "C:\\Documents and Settings\\npankrat\\My
  // Documents\\tWork\\Global PD files\\year_bug_dnas.txt";
  public static final String DNA_YEAR_BUG_SOURCE =
                                                 "C:\\Documents and Settings\\npankrat\\My Documents\\1_CRFdb\\Global PD files\\year_bug_dnas.txt";

  private final Hashtable<String, Vector<Vector<String>>> hashMatches;
  private final Hashtable<String, String> hashDNAs;
  private final Hashtable<String, String> yearBugs;
  private boolean yearIssue;

  public CheckIDsAgainstDNAs() {
    BufferedReader reader;
    String[] line;
    String trav;
    Vector<Vector<String>> idMatches;

    hashMatches = new Hashtable<String, Vector<Vector<String>>>();
    hashDNAs = new Hashtable<String, String>();
    yearBugs = new Hashtable<String, String>();
    yearIssue = false;

    loadNinfoType1(1);
    loadNinfoType1(5);

    try {
      reader = tools.getNinfoReader(3, false);
      while (reader.ready()) {
        line = reader.readLine().split("[\\s]+");
        if (line[0].equals("SWAPPED")) {
          trav = line[1] + "-" + line[2];
          if (!hashMatches.containsKey(trav)) {
            System.err.println("Error - '" + trav
                               + "' is listed as a swap in ninfo3 but was not found in ninfo1");
          } else {
            idMatches = hashMatches.get(trav);
            if (line[4].equals("resent")) {
              idMatches.elementAt(3).add(line[3]);
            } else {
              idMatches.elementAt(0).remove(line[3]);
              idMatches.elementAt(1).add(line[3]);
              if (!line[4].equals("destroy")) {
                idMatches.elementAt(0).add(line[4]);
                idMatches.elementAt(2).add(line[4]);
              }
            }
          }
        }
      }
      reader.close();
    } catch (IOException ioe) {
      System.err.println("Error parsing ninfo3 file");
      System.exit(3);
    }

    try {
      reader = new BufferedReader(new FileReader(DNA_YEAR_BUG_SOURCE));
      ext.checkHeader(reader.readLine().trim().split("[\\s]+"),
                      new String[] {"UniqueID", "FamID", "IndID", "WrongYear", "CorrectYear"},
                      true);
      while (reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");
        yearBugs.put(line[1] + ":" + line[2] + ":" + line[3], line[4]);
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + DNA_YEAR_BUG_SOURCE
                         + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + DNA_YEAR_BUG_SOURCE + "\"");
      System.exit(2);
    }
  }

  public void loadNinfoType1(int ninfo) {
    BufferedReader reader;
    String[] line;
    String trav;
    Vector<Vector<String>> idMatches;

    try {
      reader = tools.getNinfoReader(ninfo, false);
      while (reader.ready()) {
        line = reader.readLine().split("[\\s]+");
        trav = line[0] + "-" + line[1];
        idMatches = HashVec.newVecVecString(4);
        line = line[4].split(":");
        for (int i = 0; i < line.length; i++) {
          if (!line[i].equals("")) { // recently added to address
            // weird new database issue with
            // blank DNAs
            // idMatches.elementAt(0).add(line[i]); // all DNAs are
            // added to the first Vector
            idMatches.elementAt(0).add(line[i].toUpperCase()); // all
            // DNAs
            // are
            // added
            // to
            // the
            // first
            // Vector
          }
          hashDNAs.put(line[i], trav);
        }
        hashMatches.put(trav, idMatches);
      }
      reader.close();
    } catch (IOException ioe) {
      System.err.println("Error parsing ninfo1 file");
      System.exit(3);
    }
  }

  public String checkPair(String famInd, String dna, boolean verbose) {
    String[] id;

    if (famInd.startsWith("ND") || dna.startsWith("ND")) {
      // return famID+"\t"+famID+"\t"+famID;
      return "";
    }
    id = org.genvisis.park.tools.getFamID(famInd);
    if (dna.contains("AD") || famInd.startsWith("71")) {
      return id[0] + "\t" + id[1] + "\t" + dna;
    }

    return checkPair(id[0], id[1], dna, verbose);
  }

  public String checkUsingUniqueID(String uniqueID, String dna, boolean verbose) {
    String[] id;

    if (uniqueID.startsWith("ND") || dna.startsWith("ND")) {
      // return famID+"\t"+famID+"\t"+famID;
      return "";
    }
    id = new String[] {uniqueID.substring(0, 5), Integer.parseInt(uniqueID.substring(5)) + ""};
    // if (dna.contains("AD") || uniqueID.startsWith("71")) {
    if (dna.contains("AD")) {
      return id[0] + "\t" + id[1] + "\t" + dna;
    }

    return checkPair(id[0], id[1], dna, verbose);
  }

  public String checkPair(String famID, String indID, String dna, boolean verbose) {
    Vector<Vector<String>> idMatches;
    String trav = famID + "-" + indID;

    if (famID.startsWith("ND") || dna.startsWith("ND")) {
      // return famID+"\t"+famID+"\t"+famID;
      return "";
    }

    if (!hashMatches.containsKey(trav)) {
      System.err.print(trav + " appears to have been typed, but was not found in ninfo1 file");
      if (hashDNAs.containsKey(dna)) {
        System.err.println("; its DNA '" + dna + "' maps to individual " + hashDNAs.get(dna));
        return dna + "\t" + famID + "\t" + indID + ":" + dna + "\t"
               + hashDNAs.get(dna).replace('-', '\t');
      } else {
        System.err.println("; its DNA '" + dna + "' similarly doesn't match to squat");
        return dna + "\t" + famID + "\t" + indID + ":???????????????????????????????";
      }
    } else if (dna.equals("")) {
      idMatches = hashMatches.get(trav);
      if (idMatches.elementAt(0).size() == 1) {
        if (verbose) {
          System.err.println(trav + " has no DNA number, suggesting replacement");
        }
        // return
        // dna+"\t"+famID+"\t"+indID+":"+idMatches.elementAt(0).elementAt(0)+"\t"+famID+"\t"+indID;
        return "";
      } else {
        if (verbose) {
          System.err.println("Having trouble filling in the blank for " + trav
                             + " -- possibilities include:");
        }
        // temp = dna+"\t"+famID+"\t"+indID;
        for (int i = 0; i < idMatches.elementAt(0).size(); i++) {
          if (verbose) {
            System.err.print(" " + idMatches.elementAt(0).elementAt(i));
          }
          // temp += ":"+idMatches.elementAt(0).elementAt(i)+"\t"+famID+"\t"+indID;
        }
        // return temp;
        return "";
      }
    } else {
      idMatches = hashMatches.get(trav);
      if (idMatches.elementAt(0).contains(dna)) {
        if (idMatches.elementAt(3).contains(dna)) {
          if (verbose) {
            System.out.println("Make sure " + dna + " (for " + trav
                               + ") is data from the resent vial; there were problems with data before Jan 1st, 2007");
          }
        } else if (idMatches.elementAt(2).contains(dna)) {
          if (verbose) {
            System.out.println(dna + " was correctly swapped for " + trav + "");
          }
        }

      } else {
        try {
          if (yearBugs.containsKey(famID + ":" + indID + ":" + dna)) {
            if (!yearIssue) {
              yearIssue = true;
              System.err.println("Error - note this set of samples has that weird year issue!!!!");
            }
            return "yearbug\t" + yearBugs.get(famID + ":" + indID + ":" + dna) + "\t" + famID + "\t"
                   + indID;
          }
          // if (dna.substring(4,6).equals("PD")) {
          // int year = Integer.parseInt(dna.substring(0,4));
          // int sample = Integer.parseInt(dna.substring(6));
          // if (dna.equals("2003PD1006") || (year == 2003 && sample
          // >= 1179 || sample <= 1203)) {
          // if (!yearIssue) {
          // yearIssue = true;
          // System.err.println("Error - note this set of samples has
          // that weird year issue!!!!");
          // }
          // return
          // dna+"\t"+famID+"\t"+indID+":"+idMatches.elementAt(0).elementAt(0)+"\t"+famID+"\t"+indID;
          // return "";
          // }
          // }
        } catch (Exception e) {
        }

        System.err.print(dna + " was typed as " + trav);
        if (idMatches.elementAt(1).contains(dna)) {
          System.err.println(", but this relationship was found to be flawed");
          return dna + "\t" + famID + "\t" + indID + ":???????????????????????????????";
        } else {
          if (idMatches.elementAt(0).size() == 1) {
            System.err.println("; did you mean " + idMatches.elementAt(0).elementAt(0)
                               + (hashDNAs.containsKey(dna) ? " or " + hashDNAs.get(dna) : ""));
            return dna + "\t" + famID + "\t" + indID + ":" + idMatches.elementAt(0).elementAt(0)
                   + "\t" + famID + "\t" + indID
                   + (hashDNAs.containsKey(dna) ? ":" + dna + "\t"
                                                  + hashDNAs.get(dna).replace('-', '\t')
                                                : "");
          } else {
            System.err.print("; is this a typo? Valid choices include:");
            // temp = dna+"\t"+famID+"\t"+indID;
            for (int i = 0; i < idMatches.elementAt(0).size(); i++) {
              System.err.print(" " + idMatches.elementAt(0).elementAt(i));
              // temp += ":"+idMatches.elementAt(0).elementAt(i)+"\t"+famID+"\t"+indID;
            }
            System.err.println((hashDNAs.containsKey(dna) ? " or perhaps " + hashDNAs.get(dna)
                                                          : ""));
            return (hashDNAs.containsKey(dna) ? ":" + dna + "\t" + famID + "\t" + indID : "");
          }
        }
      }
    }
    return "";
  }

  public static void checkChromes(boolean verbose) throws IOException {
    BufferedReader reader = null;
    PrintWriter writer = null;
    String[] line;
    String suggestion;
    int count = 0;
    CheckIDsAgainstDNAs idCheck = new CheckIDsAgainstDNAs();

    while (count <= 23 && !new File("chromosome" + (++count) + ".dat").exists()) {
      ;
    }
    if (count > 23) {
      System.err.println("Error no valid chromosome#.dat files in current directory");
      System.exit(1);
    }

    writer = new PrintWriter(new FileWriter("fixit_suggestions.dat"));
    reader = new BufferedReader(new FileReader("chromosome" + count + ".dat"));
    reader.readLine();
    reader.readLine();
    while (reader.ready()) {
      line = reader.readLine().split("[\\s]+");
      suggestion = idCheck.checkPair(line[1], line[2], line[0], verbose);
      if (suggestion.length() > 0) {
        writer.println(suggestion);
      }
    }
    reader.close();
    writer.close();

  }

  public static void checkFile(String filename) {
    BufferedReader reader = null;
    String[] line, split;
    CheckIDsAgainstDNAs idCheck;
    boolean labFormat = false;

    try {
      idCheck = new CheckIDsAgainstDNAs();
      reader = new BufferedReader(new FileReader("check this.txt"));
      while (reader.ready()) {
        line = reader.readLine().split("[\\s]+");
        if (line.length == 2 || line[1].indexOf("-") > 0) {
          labFormat = true;
          split = tools.getFamID(line[1]);
          line = new String[] {line[0], split[0], split[1]};
        } else if (labFormat) {
          System.err.println("Error - some of these IDs have the lab format (i.e. 70345-1) and others do not");
        }
        idCheck.checkPair(line[1], line[2], line[0], true);
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + "check this.txt" + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + "check this.txt" + "\"");
      System.exit(2);
    }

  }

  public static void fixit(String fixit) throws IOException {
    BufferedReader reader;
    PrintWriter writer;
    String bakfile, skipped = "";
    Hashtable<String, String[]> hash = new Hashtable<String, String[]>();
    String[] line, trav;
    String temp;

    try {
      reader = new BufferedReader(new FileReader(fixit));
      while (reader.ready()) {
        temp = reader.readLine();
        line = temp.split(":", -1);
        if (line.length == 1 || line[0].split("\t", -1).length != 3
            || line[1].split("\t", -1).length != 3) {
          System.err.println("Error - malformed fixit syntax, expecting 3 tab delimited Strings on each side of a single semicolon(:)");
          System.err.println("        '" + temp + "'");
        } else if (line.length > 2) {
          System.err.println("Error - malformed fixit syntax, requires disambiguation:");
          System.err.println("        '" + temp + "'");
        } else {
          trav = line[0].split("\t", -1);
          hash.put(trav[0] + ":" + trav[1] + ":" + trav[2], line[1].split("\t", -1));
        }
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error - could not find " + fixit + " in current directory");
      System.exit(2);
    } catch (IOException ioe) {
      System.err.println("Error parsing " + fixit + "");
      System.exit(3);
    }

    for (int i = 1; i <= 23; i++) {
      try {
        bakfile = Files.getBakFilename("chromosome" + i, "checkIDsAgaintsDNAs");
        new File("chromosome" + i + ".dat").renameTo(new File(bakfile));
        reader = new BufferedReader(new FileReader(bakfile));
        writer = new PrintWriter(new FileWriter("chromosome" + i + ".dat"));
        writer.println(reader.readLine());
        writer.println(reader.readLine());
        while (reader.ready()) {
          line = reader.readLine().split("[\\s]+");
          temp = line[0] + ":" + line[1] + ":" + line[2];
          if (hash.containsKey(temp)) {
            trav = hash.get(temp);
            line[0] = trav[0];
            line[1] = trav[1];
            line[2] = trav[2];
          }
          for (int j = 0; j < line.length; j++) {
            writer.print((j == 0 ? "" : "\t") + line[j]);
          }
          writer.println();
        }
        reader.close();
        writer.close();
      } catch (FileNotFoundException fnfe) {
        skipped += " " + i;
      } catch (IOException ioe) {
        System.err.println("Error parsing " + "chromosome" + i + ".dat" + "");
        System.exit(3);
      }
    }
    System.err.println("Skipped chromosomes:" + skipped);
  }

  public static void main(String[] args) throws IOException {
    int numArgs = args.length;
    boolean verbose = false;
    String filename = "check this.txt";
    String fixit = "fix this.txt";

    String usage = "\n" + "park.init.CheckIDsAgainstDNAs requires 0-1 arguments\n"
                   + "     (1) verbose (i.e. verbose=" + verbose + " (default))\n"
                   + "     (2) name of file to check (i.e. file=" + filename + " (default))\n"
                   + "     (3) name of file containing corrections  (i.e. fixit=" + fixit
                   + " (default))\n" + "";

    for (int i = 0; i < args.length; i++) {
      if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h")
          || args[i].equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (args[i].startsWith("verbose=")) {
        verbose = !args[i].split("=")[1].toLowerCase().equals("false");
        numArgs--;
      } else if (args[i].startsWith("file=")) {
        filename = args[i].split("=")[1];
        numArgs--;
      } else if (args[i].startsWith("fixit=")) {
        fixit = args[i].split("=")[1];
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      if (new File(fixit).exists()) {
        fixit(fixit);
      } else if (new File(filename).exists()) {
        checkFile(filename);
      } else {
        CheckIDsAgainstDNAs.checkChromes(verbose);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
