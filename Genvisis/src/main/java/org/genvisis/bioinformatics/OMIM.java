package org.genvisis.bioinformatics;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;

/**
 * @author lane0212
 *
 *         Query da OMIM db if ya want
 *
 */
public class OMIM {

  public static class OMIMGene {
    private final String[] geneSymbols;
    private final String CytogeneticLocation;
    private final String geneStatus;
    private final String title;
    private final String MIMNumber;
    private final String method;
    private final String disorders;
    private final String status;
    private final String references;
    private final String[] all;

    public OMIMGene(String[] all, String geneSymbol[], String cytogeneticLocation,
                    String geneStatus, String title, String mIMNumber, String method,
                    String disorders, String references, String diseaseEvidence) {
      super();
      this.all = all;
      geneSymbols = geneSymbol;
      CytogeneticLocation = cytogeneticLocation;
      this.geneStatus = geneStatus;
      this.title = title;
      MIMNumber = mIMNumber;
      this.method = method;
      this.disorders = disorders;
      this.references = references;
      status = diseaseEvidence;
    }

    public String[] getAll() {
      return all;
    }

    public String getCytogeneticLocation() {
      return CytogeneticLocation;
    }

    public String getDisorders() {
      return disorders;
    }

    public String getGeneStatus() {
      return geneStatus;
    }

    public String[] getGeneSymbols() {
      return geneSymbols;
    }

    public String getMethod() {
      return method;
    }

    public String getMIMNumber() {
      return MIMNumber;
    }

    public String getReferences() {
      return references;
    }

    public String getStatus() {
      return status;
    }

    public String getTitle() {
      return title;
    }

  }

  private static final String GENE_MAP_FILE = "geneMap2.txt";
  private static final String OMIM_TEXT = "omim.txt.Z";

  private static Hashtable<String, ArrayList<OMIMGene>> loadGeneOmimMap(String filename,
                                                                        Logger log) {
    Hashtable<String, ArrayList<OMIMGene>> gHashtable =
        new Hashtable<String, ArrayList<OMIMGene>>();
    Hashtable<String, String> status = new Hashtable<String, String>();
    status.put("C", "confirmed - observed in at least two laboratories or in several families");
    status.put("P", "based on evidence from one laboratory or one family");
    status.put("I", "inconsistent - results of different laboratories disagree");
    status.put("L",
               "limbo - evidence not as strong as that provisional, but included for heuristic reasons");

    try {
      BufferedReader reader = Files.getAppropriateReader(filename);
      while (reader.ready()) {
        String[] line = reader.readLine().trim().split("\\|");
        String[] tmp = new String[15];
        Arrays.fill(tmp, "NA");
        System.arraycopy(line, 0, tmp, 0, line.length - 1);
        if (tmp.length != 15) {
          throw new IllegalArgumentException(filename + " must have 15 entries per line");
        }
        String[] genes = Array.unique(line[5].split(","));
        for (int i = 0; i < genes.length; i++) {
          if (!gHashtable.containsKey(genes[i])) {
            gHashtable.put(genes[i], new ArrayList<OMIM.OMIMGene>());
          }
          OMIMGene mGene =
              new OMIMGene(line, genes, tmp[4], tmp[6], tmp[7], tmp[8], tmp[9], tmp[11], tmp[13],
                           status.containsKey(tmp[6]) ? status.get(tmp[6]) : tmp[6]);
          gHashtable.get(genes[i]).add(mGene);
        }
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      log.reportError("Error: file \"" + filename + "\" not found in current directory");
      return null;
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + filename + "\"");
      return null;
    }
    return gHashtable;
  }

  public static void main(String[] args) {
    String dir = "C:/bin/ref/OMIM/";
    new OMIM(dir, new Logger());
  }

  private final Hashtable<String, ArrayList<OMIMGene>> gHashtable;

  private final String omimDir;

  private final Logger log;

  // private static class OMIMDisorder {
  // private String name;
  // private String[] genes;
  // private String[] geneMims;
  // private String cytoLoc;
  //
  // }

  public OMIM(String omimDir, Logger log) {
    super();
    this.omimDir = omimDir;
    this.log = log;
    verify();
    gHashtable = loadGeneOmimMap(omimDir + GENE_MAP_FILE, log);
    log.reportTimeInfo("Loaded omim db " + omimDir + GENE_MAP_FILE);
  }

  public ArrayList<OMIMGene> getOmimGene(String gene) {
    if (gHashtable.containsKey(gene)) {
      return gHashtable.get(gene);
    } else {
      ArrayList<OMIMGene> blank = new ArrayList<OMIM.OMIMGene>();
      blank.add(new OMIMGene(new String[] {"NA"}, new String[] {"NA"}, "NA", "NA", "NA", "NA", "NA",
                             "NA", "NA", "NA"));
      return blank;
    }
  }

  private boolean verify() {
    if (!Files.exists(omimDir + GENE_MAP_FILE)) {
      log.reportFileNotFound(omimDir + GENE_MAP_FILE);
      throw new IllegalArgumentException();
    }
    if (!Files.exists(omimDir + OMIM_TEXT)) {
      log.reportFileNotFound(omimDir + OMIM_TEXT);
      throw new IllegalArgumentException();
    }
    return true;
  }
}
