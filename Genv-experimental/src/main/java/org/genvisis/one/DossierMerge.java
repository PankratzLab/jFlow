package org.genvisis.one;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Internat;
import org.genvisis.common.Logger;
import org.genvisis.common.RealTextFormatWriter;
import org.genvisis.common.ext;

public class DossierMerge {
  private static void getCitations(String dir, boolean banned) {
    PrintWriter writer;
    String[] line;
    String[] pmids;
    String pmid, doi;
    String[] results;
    String scopusCitationCount, googleScholarCitationCount, numGS_results, title, temp;
    Logger log = new Logger();
    int[] indices;
    int index;
    ArrayList<String> authors;

    pmids =
        HashVec.loadFileToStringArray(dir + "personal bibliography with DOI from Endnote.utf8.txt",
                                      false, false, new int[] {0, 1}, false, false, "\t");

    String filename = "citations" + ext.getTimestampForFilename() + ".dat";
    try {
      writer = Files.getAppropriateWriter(dir + filename);
      writer.println("PMID\tScopus\tGoogle\tUnique Google Scholar result\tdoi");
      for (int i = 0; i < pmids.length; i++) {
        line = pmids[i].split("\t", -1);
        pmid = line[0];
        doi = line[1].equals("") ? null : line[1];

        System.out.println((i + 1) + " of " + pmids.length + "\tQuerying PubMed for PMID:" + pmid
                           + "  doi:" + doi);
        try {
          Thread.sleep(1000);
        } catch (InterruptedException ie) {
        }

        results = Internat.getPage("http://www.ncbi.nlm.nih.gov/pubmed/" + pmid);
        // Files.writeList(results, dir+"pubmed_"+pmids[i]+".out");

        doi = title = null;
        authors = new ArrayList<String>();
        try {
          for (int j = 0; (doi == null || title == null || authors.size() == 0)
                          && j < results.length; j++) {
            if (doi == null && results[j].contains("doi: ")) {
              doi = results[j].substring(results[j].indexOf("doi: ") + 5);
              doi = doi.split("[\\s]+")[0];
              if (doi.contains("</div>")) {
                doi = doi.substring(0, doi.indexOf("</div>"));
              }
              if (doi.endsWith("\"")) {
                doi = doi.substring(0, doi.length() - 1);
              }
              if (doi.endsWith(".")) {
                doi = doi.substring(0, doi.length() - 1);
              }
              System.out.println("   found missing doi: " + doi);
            }
            if (title == null && results[j].contains("<title>")) {
              title = results[j].substring(results[j].indexOf("<title>") + 7);
              index = title.indexOf(".  - PubMed - NCBI");
              if (index == -1) {
                System.out.println("Error - PubMed format appears to have changed, no longer getting title between <title> and '.  - PubMed - NCBI'");
              } else {
                title = title.substring(0, index);
              }
            }
            if (results[j].contains("%5BAuthor%5D")) {
              temp = results[j];

              while (temp.contains("%5BAuthor%5D")) {
                temp = temp.substring(temp.indexOf("%5BAuthor%5D"));
                temp = temp.substring(temp.indexOf(">") + 1);
                authors.add(temp.substring(0, temp.indexOf("<")));
              }
              // System.err.println(" parsed "+authors.size()+" authors");
            }
          }
        } catch (Exception e) {
          System.err.println("    ran into trouble parsing PubMed entry for " + pmid + " / " + doi);
          e.printStackTrace();
        }

        if (doi == null) {
          System.err.println("	no doi on PubMed");
        }

        if (doi == null) {
          System.out.println((i + 1) + " of " + pmids.length + "\tGetting Scopus result for " + pmid
                             + " using PMID");
          results =
              Internat.getPage("http://api.elsevier.com/content/search/index:SCOPUS?query=PMID("
                               + pmid + ")&apiKey=c7af0f4beab764ecf68568961c2a21ea");
        } else {
          System.out.println((i + 1) + " of " + pmids.length + "\tGetting Scopus result for " + pmid
                             + " using doi: " + doi);
          results =
              Internat.getPage("http://api.elsevier.com/content/search/index:SCOPUS?query=DOI("
                               + doi + ")&apiKey=c7af0f4beab764ecf68568961c2a21ea");
        }
        // Files.writeList(results, dir+"scopus.out");
        //

        scopusCitationCount = results[0];
        try {
          scopusCitationCount =
              results[0].substring(results[0].indexOf("\"citedby-count\":") + 16 + 1);
          scopusCitationCount = scopusCitationCount.substring(0, scopusCitationCount.indexOf("\""));
        } catch (Exception e) {
          System.err.println("Error - failed to find scopus citation count for " + pmid);
          e.printStackTrace();
        }

        doi = results[0];
        try {
          doi = results[0].substring(results[0].indexOf("\"prism:doi\":") + 12 + 1);
          doi = doi.substring(0, doi.indexOf("\""));
        } catch (Exception e) {
          System.err.println("Error - failed to find doi from scopus for " + pmid);
          e.printStackTrace();
        }

        // doi = "10.1016/S0140-6736(05)17828-3";


        // override

        if (doi != null && !doi.equals("ults")) {
          // if you get black listed by Google, then just reverse the commenting of these two lines
          if (banned) {
            results =
                new String[] {"https://scholar.google.com/scholar?q=http%3A%2F%2Fdx.doi.org%2F"
                              + doi};
          } else {
            results =
                Internat.getMash("https://scholar.google.com/scholar?q=http%3A%2F%2Fdx.doi.org%2F"
                                 + doi, new Logger());
          }
        } else {
          if (banned) {
            results = new String[] {"https://scholar.google.com/scholar?q=\""
                                    + ext.replaceAllWith(title, " ", "+") + "\""};
          } else {
            System.out.println("Getting Google Scholar result for " + pmid
                               + " at https://scholar.google.com/scholar?q=\""
                               + ext.replaceAllWith(title, " ", "+") + "\"");
            results = Internat.getMash("https://scholar.google.com/scholar?q=\""
                                       + ext.replaceAllWith(title, " ", "+") + "\"", new Logger());
          }
        }

        googleScholarCitationCount = results[0];

        if (banned) {
          numGS_results = "1";
        } else {
          numGS_results = "missing";
          try {
            numGS_results = results[0].substring(0, results[0].indexOf(" result"));
            indices = ext.indicesWithinString(">", numGS_results);
            numGS_results = numGS_results.substring(indices[indices.length - 1] + 1);
          } catch (Exception e) {
            System.err.println("Error - failed to find the number of Google Scholar results from the query for "
                               + pmid);
            e.printStackTrace();
          }

          googleScholarCitationCount = "missing";
          try {
            googleScholarCitationCount = results[0].substring(results[0].indexOf("Cited by ") + 9);
            googleScholarCitationCount =
                googleScholarCitationCount.substring(0, googleScholarCitationCount.indexOf("</a>"));
          } catch (Exception e) {
            System.err.println("Error - failed to find the number of Google Scholar results from the query for "
                               + pmids[i]);
            e.printStackTrace();
          }
        }

        if (numGS_results.equals("Showing the best")) {
          numGS_results = "1";
        }

        if (googleScholarCitationCount.equals("missing") || !numGS_results.equals("1")) {
          System.err.println("Error - failed to properly parse google scholar results for " + pmid);
          Files.writeList(results, dir + "scholar_" + pmid + ".out");
        }

        // googleScholarCitationCount = numGS_results = "1";
        writer.println(pmid + "\t" + scopusCitationCount + "\t" + googleScholarCitationCount + "\t"
                       + (numGS_results.equals("1") ? "TRUE" : "FALSE (" + numGS_results + ")")
                       + "\t" + doi + "\t" + title);
        writer.flush();
      }
      writer.close();
    } catch (Exception e) {
      log.reportError("Error writing to " + dir + filename);
      log.reportException(e);
    }

  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String dir = "D:/umn/promotion/citations/";
    boolean rtf = true;
    boolean getCitations = false;
    boolean banned = true;

    String usage = "\n" + "one.AuthorCorral requires 2+ arguments\n"
                   + "   (1) working directory (i.e. dir=" + dir + " (default))\n" + ""
                   + "   (2) output in RTF format (used for superscripts) (i.e. rtf=" + rtf
                   + " (default))\n";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("dir=")) {
        dir = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("banned=")) {
        banned = ext.parseBooleanArg(arg);
        numArgs--;
      } else if (arg.startsWith("rtf=")) {
        rtf = Boolean.valueOf(arg.split("=")[1]);
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
      if (getCitations) {
        getCitations(dir, banned);
      } else {
        run(dir, rtf);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public static void run(String dir, boolean rtfOutput) throws IOException {
    RealTextFormatWriter writer;
    String[] keys, bits, line;
    Hashtable<String, String> citationHash, timesCitedHash, impactHash, rolesHash,
        abbreviationsHash;
    HashSet<String> authorsToBoldHash, jointFirstAuthorsHash;
    String pmid, citation, journal, timesCited, impactFactor, role;
    Logger log;
    boolean jointFirstAuthor;
    int count;

    writer = new RealTextFormatWriter(dir + "formattedDossier" + (rtfOutput ? ".rtf" : ".out"),
                                      rtfOutput);
    log = new Logger(dir + "formattedDossier.log");

    // citationHash = HashVec.loadFileToHashString(dir+"pubmed_result.utf8.csv", new int[] {9}, new
    // int[] {2, 0, 3}, true, "\t", true, false, false);
    citationHash =
        HashVec.loadFileToHashString(dir + "personal bibliography with DOI from Endnote.utf8.txt",
                                     new int[] {0}, new int[] {2, 3, 4}, false, "\t", true, false,
                                     false);
    abbreviationsHash =
        HashVec.loadFileToHashString(dir + "Abbreviations.txt", new int[] {0}, new int[] {1}, false,
                                     null, false, false, false);

    impactHash = HashVec.loadFileToHashString(dir + "impactFactors.dat", new int[] {0},
                                              new int[] {2, 1}, false, " in ", false, false, false);
    timesCitedHash =
        HashVec.loadFileToHashString(dir + "citations.dat", new int[] {0}, new int[] {1, 2}, false,
                                     "\t", false, false, false);
    rolesHash = HashVec.loadFileToHashString(dir + "roles.dat", new int[] {0}, new int[] {1}, false,
                                             null, false, false, false);


    authorsToBoldHash = HashVec.loadToHashSet(HashVec.loadFileToStringArray(dir + "wordsToBold.dat",
                                                                            false, null, false));
    jointFirstAuthorsHash = HashVec.loadFileToHashSet(dir + "jointFirstAuthors.dat", false);

    new File(dir + "journalsInDatabase.xln").delete();


    count = 0;
    keys = HashVec.getKeys(citationHash, true, true);
    for (int i = 0; i < keys.length; i++) {
      pmid = keys[keys.length - i - 1];
      citation = citationHash.remove(pmid);
      impactFactor = impactHash.remove(pmid);
      timesCited = timesCitedHash.remove(pmid);
      role = rolesHash.remove(pmid);
      jointFirstAuthor = jointFirstAuthorsHash.contains(pmid);

      if (role == null || !role.equals("None")) {
        count++;
        System.out.println(pmid);
        bits = citation.split("\t", -1);

        line = bits[0].substring(0, bits[0].length() - 1).split(",");
        for (int j = 0; j < line.length; j++) {
          String author = line[j].trim();
          if (authorsToBoldHash.contains(author)) {
            line[j] = "<b>" + line[j] + "</>";
          }
        }
        if (jointFirstAuthor) {
          line[0] += "<super>�</>";
          line[1] += "<super>�</>";
        }
        bits[0] = Array.toStr(line, ",") + ".";

        if (bits[2].contains(" doi:")) {
          bits[2] = bits[2].substring(0, bits[2].indexOf(" doi:"))
                    + (bits[2].contains("[Epub ahead of print]") ? " [Epub ahead of print]" : "");
        }

        journal = bits[2].substring(0, bits[2].indexOf("."));
        if (abbreviationsHash.containsKey(journal)) {
          ext.appendToFile(journal + "\t" + abbreviationsHash.get(journal),
                           dir + "journalsInDatabase.xln");
          journal = abbreviationsHash.get(journal);
        } else if (!abbreviationsHash.contains(journal)) {
          ext.appendToFile(journal + "\t.", dir + "journalsInDaabase.xln");
          log.reportError("Error journal '" + journal
                          + "' was not found as an abbreviation or as a full journal name");
        }
        bits[2] = "<u>" + journal + "</>" + bits[2].substring(bits[2].indexOf("."));

        writer.println(count + "." + (rtfOutput ? "\\tab " : "\t") + Array.toStr(bits, " "));

        if (jointFirstAuthor) {
          writer.println((rtfOutput ? "\\bullet \\tab " : "\t") + "<super>�</>Joint first author");
        }

        if (impactFactor != null) {
          writer.println((rtfOutput ? "\\bullet \\tab " : "\t") + "Journal Impact Factor: "
                         + impactFactor);
        } else {
          writer.println((rtfOutput ? "\\bullet \\tab " : "\t")
                         + "Journal Impact Factor: Not available");
          log.reportError("Error - no impact factor available for " + pmid + " (" + citation + ")");
        }

        if (timesCited != null) {
          line = timesCited.split("[\\s]+");
          writer.println((rtfOutput ? "\\bullet \\tab " : "\t") + "Times Cited: " + line[0]
                         + " (Scopus) / " + line[1] + " (Google Scholar)");
        } else {
          writer.println((rtfOutput ? "\\bullet \\tab " : "\t") + "Times Cited: Not available");
          log.reportError("Error - no times cited available for " + pmid + " (" + citation + ")");
        }

        if (role != null) {
          writer.println((rtfOutput ? "\\bullet \\tab " : "\t") + "Role: " + role);
        } else {
          writer.println((rtfOutput ? "\\bullet \\tab " : "\t") + "Role: ");
          log.reportError("Error - no role available for " + pmid + " (" + citation + ")");
        }

        writer.newParagraph();
      }
    }

    writer.close();
  }
}
