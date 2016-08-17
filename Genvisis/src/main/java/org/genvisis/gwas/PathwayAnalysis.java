package org.genvisis.gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Date;
import java.util.Hashtable;
import java.util.Vector;

import org.genvisis.bioinformatics.MapGenesToSNPs;
import org.genvisis.common.Array;
import org.genvisis.common.CmdLine;
import org.genvisis.common.CountVector;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;
import org.genvisis.filesys.SerialFloatArray;
import org.genvisis.filesys.SerialStringArray;
import org.genvisis.filesys.SnpMarkerSet;

class KEGGpathway {
  public Hashtable<String, String> pathwayLookup;
  public Hashtable<String, Vector<String>> pathways, genes, snpsToGenes, genesToSNPs;

  public KEGGpathway(String ko_file, boolean gene_ids_instead_of_gene_names) {
    BufferedReader reader;
    String[] line;
    String temp;
    Vector<String> v;
    String pathway, pathwayName;
    String gene, entry;

    pathwayLookup = new Hashtable<String, String>();
    pathways = new Hashtable<String, Vector<String>>();
    genes = new Hashtable<String, Vector<String>>();
    v = new Vector<String>();
    entry = "the very beginning";
    try {
      reader = new BufferedReader(new FileReader(ko_file));
      temp = reader.readLine();
      while (reader.ready()) {
        if (temp.startsWith("PATHWAY")) {
          temp = temp.substring(("PATHWAY").length());
          if (v.size() > 0) {
            System.err.println("Error - did not reset after ENTRY " + entry + " starting at: "
                               + temp);
          }
          v = new Vector<String>();
          while (temp.startsWith(" ")) {
            temp = temp.trim();
            pathway = temp.split("[\\s]+")[0];
            pathwayName = temp.substring(pathway.length()).trim();
            if (!pathwayLookup.containsKey(pathway)) {
              pathwayLookup.put(pathway, pathwayName);
            } else if (!pathwayLookup.get(pathway).equals(pathwayName)) {
              System.err.println("Error - mismatched ID (" + pathway + ") with pathway name (was '"
                                 + pathwayLookup.get(pathway) + "', now '" + pathwayName + "')");
            }
            v.add(pathway);
            temp = reader.readLine();
          }
        }
        if (temp.startsWith("ENTRY")) {
          v = new Vector<String>();
          entry = temp.trim().split("[\\s]+")[1];
          // System.err.println(entry);
        }
        if (temp.startsWith("GENES")) {
          temp = temp.substring(("GENES").length());
          while (temp.startsWith(" ")) {
            temp = temp.trim();
            if (temp.startsWith("HSA:")) {
              line = temp.split("[\\s]+");
              for (int i = 1; i < line.length; i++) {
                if (line[i].indexOf("(") == -1 || line[i].indexOf(")") == -1) {
                  gene = line[i];
                  // System.err.println("Error - no gene name for '"+line[i]+"' in entry "+entry+"
                  // for pathways: "+ext.listWithCommas(Array.toStringArray(v)));
                } else {
                  if (gene_ids_instead_of_gene_names) {
                    gene = line[i].substring(0, line[i].indexOf("("));
                  } else {
                    gene = line[i].substring(line[i].indexOf("(") + 1, line[i].indexOf(")"));
                  }
                }
                for (int j = 0; j < v.size(); j++) {
                  HashVec.addToHashVec(genes, gene, v.elementAt(j), true); // count both ways and
                                                                           // compare
                  HashVec.addToHashVec(pathways, v.elementAt(j), gene, true); // count both ways and
                                                                              // compare
                }
              }
            }
            temp = reader.readLine();
          }
        }
        temp = reader.readLine();
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + ko_file + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + ko_file + "\"");
      System.exit(2);
    }
  }

  public KEGGpathway(String ko2_file, String ko2_names_file) {
    BufferedReader reader;
    String[] line;

    pathwayLookup = HashVec.loadFileToHashString(ko2_names_file, false);
    pathways = new Hashtable<String, Vector<String>>();
    genes = new Hashtable<String, Vector<String>>();
    try {
      reader = new BufferedReader(new FileReader(ko2_file));
      while (reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");
        for (int i = 1; i < line.length; i++) {
          HashVec.addToHashVec(genes, line[0], "ko" + line[i], true); // count both ways and compare
          HashVec.addToHashVec(pathways, "ko" + line[i], line[0], true); // count both ways and
                                                                         // compare
        }
      }
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + ko2_file + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + ko2_file + "\"");
      System.exit(2);
    }
  }

  public void matchMarkersToGenes(String mapfile, String geneLocFile, int offset) {
    String[] line;
    Hashtable<String, Vector<String>> genes;
    Vector<String> v;
    String[][] genesAndSegs;
    Segment seg;
    byte maxChr;
    String[][] geneIDs;
    Segment[][] segs;
    SnpMarkerSet markerSet;
    String[] markerNames;
    byte[] chrs;
    int[] positions;
    CountVector noChr;
    String[] noChrValues;
    int[] noChrCounts;

    genesAndSegs = HashVec.loadFileToStringMatrix(geneLocFile, true, new int[] {0, 1}, false);
    genes = new Hashtable<String, Vector<String>>();
    maxChr = -1;
    for (String[] genesAndSeg : genesAndSegs) {
      seg = new Segment(genesAndSeg[1]);
      if (seg.getChr() > maxChr) {
        maxChr = seg.getChr();
      }
      HashVec.addToHashVec(genes, seg.getChr() + "", genesAndSeg[0] + "\t" + genesAndSeg[1], false);
    }
    geneIDs = new String[maxChr + 1][];
    segs = new Segment[maxChr + 1][];
    for (int i = 0; i < segs.length; i++) {
      v = genes.get(i + "");
      if (v == null) {
        geneIDs[i] = null;
        segs[i] = null;
      } else {
        geneIDs[i] = new String[v.size()];
        segs[i] = new Segment[v.size()];
        for (int j = 0; j < segs[i].length; j++) {
          line = v.elementAt(j).split("[\\s]+");
          geneIDs[i][j] = line[0];
          segs[i][j] = new Segment(line[1]);
        }
      }
    }

    snpsToGenes = new Hashtable<String, Vector<String>>();
    genesToSNPs = new Hashtable<String, Vector<String>>();

    markerSet = new SnpMarkerSet(mapfile);
    markerNames = markerSet.getMarkerNames();
    chrs = markerSet.getChrs();
    positions = markerSet.getPositions();
    noChr = new CountVector();
    for (int i = 0; i < markerNames.length; i++) {
      seg = new Segment(chrs[i], positions[i] - offset, positions[i] + offset);
      if (segs.length <= chrs[i] || segs[chrs[i]] == null) {
        noChr.add(chrs[i] + "");
      } else {
        for (int j = 0; j < segs[chrs[i]].length; j++) {
          try {
            if (seg.overlaps(segs[chrs[i]][j])) {
              HashVec.addToHashVec(snpsToGenes, markerNames[i], geneIDs[chrs[i]][j], true);
              HashVec.addToHashVec(genesToSNPs, geneIDs[chrs[i]][j], markerNames[i], true);
            }
          } catch (Exception e) {
            System.err.println("Error - problem trying to overlap " + seg.getUCSClocation()
                               + " and " + geneIDs[chrs[i]][j]);
          }
        }
      }
    }
    noChrValues = noChr.getValues();
    noChrCounts = noChr.getCounts();
    for (int i = 0; i < noChrValues.length; i++) {
      System.err.println("Error - no genes on chromosome " + noChrValues[i] + " (affects "
                         + noChrCounts[i] + " snps)");
    }
  }
}


public class PathwayAnalysis {
  // public static final String ROOT_DIRECTORY = "/home/npankrat/"; // galileo
  // public static final String ROOT_DIRECTORY = "/export/home/npankrat/"; // alcatraz
  public static final String ROOT_DIRECTORY = "/state/partition1/npankrat/"; // indiviudal nodes

  // public static final String KEGG_KO = ROOT_DIRECTORY+"NCBI/KEGG/ko.txt";
  // public static final String KEGG_KO2 = ROOT_DIRECTORY+"NCBI/KEGG/altFiles/hsa_gene_map.tab";
  // public static final String KEGG_KO2_TITLES = ROOT_DIRECTORY+"NCBI/KEGG/altFiles/map_title.tab";
  // public static final String KEGG_GROUPINGS = ROOT_DIRECTORY+"NCBI/KEGG/raw_groupings.txt";
  // public static final String KEGG_GENE_POSITIONS = ROOT_DIRECTORY+"NCBI/KEGG/geneLocs.xln";
  public static final String KEGG_KO = ROOT_DIRECTORY + "NCBI/KEGG/ko.txt";
  public static final String KEGG_KO2 = ROOT_DIRECTORY + "NCBI/KEGG/altFiles/hsa_gene_map.tab";
  public static final String KEGG_KO2_TITLES = ROOT_DIRECTORY + "NCBI/KEGG/altFiles/map_title.tab";
  public static final String KEGG_GROUPINGS = ROOT_DIRECTORY + "NCBI/KEGG/raw_groupings.txt";
  public static final String KEGG_GENE_POSITIONS = ROOT_DIRECTORY + "NCBI/KEGG/geneLocs.xln";
  public static final int DEFAULT_BUFFER = 15000;
  public static final int DEFAULT_P_THRESHOLD = 1;
  public static final int DEFAULT_NRSS_INDEX_THRESHOLD = 1;


  public static double calculateScore(float[] pvals, boolean[] included, double p_threshold,
                                      double nrss_index_threshold) {
    double score = 0;

    for (int i = 0; i < included.length; i++) {
      if (included[i] && !Float.isNaN(pvals[i])) {
        if (pvals[i] < p_threshold) {
          score += -1 * Math.log10(pvals[i]);
        }
      }
    }

    return score;
  }

  public static void checkConcentrationFromGeneList(String filename, String groupings_file) {
    BufferedReader reader;
    PrintWriter writer;
    String temp;
    String pathway;
    String dir, main, header, subheader;
    KEGGpathway kegg;
    Hashtable<String, Vector<String>> hash;
    Vector<String> v;
    int countGene, countTimesGene;

    // kegg = new KEGGpathway(keggFile, true);
    kegg = new KEGGpathway(KEGG_KO2, KEGG_KO2_TITLES);
    dir = ext.verifyDirFormat(ext.parseDirectoryOfFile(filename));
    hash = HashVec.loadFileToHashVec(filename, 0, new int[] {0}, "", false, false);

    Files.writeList(Array.toStringArray(kegg.pathways.get("ko05010")), dir + "AlzGenes.xln");

    try {
      reader = new BufferedReader(new FileReader(groupings_file));
      writer = new PrintWriter(new FileWriter(dir + "groupings2.xln"));
      writer.println("Main Heading\tHeader\tPathway\tPathwayName\tSubHeader\t# genes\t# genes Covered\t# times genes covered");
      temp = reader.readLine();
      while (reader.ready()) {
        if (temp.startsWith("<h4>")) {
          main = temp.substring(temp.indexOf("<h4>") + 4, temp.indexOf("</h4>"));
          writer.println(main);
        } else if (temp.startsWith("<b>")) {
          header = temp.substring(temp.indexOf("<b>") + 3, temp.indexOf("</b>"));
          writer.println("\t" + header);
        } else if (temp.contains("/kegg/pathway/")) {
          pathway = "ko" + temp.substring(temp.indexOf(".html") - 5, temp.indexOf(".html"));
          subheader = temp.substring(temp.indexOf(".html") + 7, temp.indexOf("</a>"));
          countGene = countTimesGene = 0;
          if (kegg.pathways.containsKey(pathway)) {
            v = kegg.pathways.get(pathway);
            for (int i = 0; i < v.size(); i++) {
              if (hash.containsKey(v.elementAt(i))) {
                countGene++;
                countTimesGene += hash.get(v.elementAt(i)).size();
              }
            }
          }
          writer.print("\t" + "\t" + pathway + "\t" + kegg.pathwayLookup.get(pathway) + "\t"
                       + subheader + "\t"
                       + (kegg.pathways.containsKey(pathway) ? kegg.pathways.get(pathway).size()
                                                             : 0));
          writer.println("\t" + countGene + "\t" + countTimesGene);
        }
        temp = reader.readLine();
      }
      reader.close();
      writer.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + groupings_file + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + groupings_file + "\"");
      System.exit(2);
    }
  }

  public static void compress(String filename, String[] markerNames, boolean deleteOriginal) {
    BufferedReader reader;
    String[] line;
    int count;
    int[] indices;
    float[] results;
    String root;

    try {
      reader = new BufferedReader(new FileReader(filename));
      indices = ext.indexFactors(new String[] {"SNP", "TEST", "P"},
                                 reader.readLine().trim().split("[\\s]+"), false, new Logger(),
                                 false, false);
      results = new float[markerNames.length];
      root = ext.rootOf(filename);
      if (root.endsWith(".assoc")) {
        root = ext.rootOf(root);
      }
      count = 0;
      while (reader.ready()) {
        line = reader.readLine().trim().split("[\\s]+");
        if (indices[1] == -1 || line[indices[1]].equals("ADD")) {
          if (!line[indices[0]].equals(markerNames[count])) {
            System.err.println("Error - marker mismatch at marker " + (count + 1) + " (expecting "
                               + markerNames[count] + ", found " + line[indices[0]] + ")");
            reader.close();
            return;
          }
          results[count] =
              line[indices[2]].equals("NA") ? Float.NaN : Float.parseFloat(line[indices[2]]);
          count++;
        }
      }
      new SerialFloatArray(results).serialize(root + ".results");
      reader.close();
      if (deleteOriginal && new File(root + ".results").exists()) {
        new File(root + ".dat").delete();
        new File(ext.rootOf(root) + "_cov" + root.substring(root.indexOf(".")) + ".dat").delete();
        new File(root + ".log").delete();
        new File(root + ".hh").delete();
        new File(filename).delete();
      }
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + filename + "\" not found in current directory");
      return;
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + filename + "\"");
      return;
    }
  }

  public static void compressAll(String pheno, String covars) {
    String[] markerNames;
    String filename;
    boolean quant, done;
    int count;

    quant = Array.determineType(pheno, 2, true) == 1;
    markerNames =
        SerialStringArray.load(ext.rootOf(pheno, false) + ".markerNames", false).getArray();

    count = 0;
    done = false;
    while (!done) {
      count++;
      if (covars == null) {
        filename = ext.rootOf(pheno, false) + "." + count + "." + (quant ? "q" : "") + "assoc";
      } else {
        filename =
            ext.rootOf(pheno, false) + "." + count + ".assoc." + (quant ? "linear" : "logistic");
      }
      if (new File(filename).exists()) {
        compress(filename, markerNames, true);
      } else {
        done = true;
      }

    }
  }

  public static void evaluateAllSets(String pheno, String mapfile, int offset, String results,
                                     String repDir, String groupings_file, double p_thresh,
                                     double nrss_thresh) {
    BufferedReader reader;
    PrintWriter writer;
    String temp;
    String pathway;
    String main, header, subheader;
    KEGGpathway kegg;
    Vector<String> genesInPathway, snpsInGene, snpsInPathway;
    int countGene;
    String pvalue;

    kegg = new KEGGpathway(KEGG_KO2, KEGG_KO2_TITLES);
    System.out.print("Mapping snps to genes...");
    kegg.matchMarkersToGenes(mapfile, KEGG_GENE_POSITIONS, offset);
    System.out.println("done");

    try {
      reader = new BufferedReader(new FileReader(groupings_file));
      writer = new PrintWriter(new FileWriter("AllPathways.xln"));
      writer.println("Main Heading\tHeader\tPathway\tPathwayName\tSubHeader\t# genes\t# genes Covered\t# snps in pathway\tEMP1");
      temp = reader.readLine();
      while (reader.ready()) {
        if (temp.startsWith("<h4>")) {
          main = temp.substring(temp.indexOf("<h4>") + 4, temp.indexOf("</h4>"));
          writer.println(main);
        } else if (temp.startsWith("<b>")) {
          header = temp.substring(temp.indexOf("<b>") + 3, temp.indexOf("</b>"));
          writer.println("\t" + header);
        } else if (temp.contains("/kegg/pathway/")) {
          pathway = "ko" + temp.substring(temp.indexOf(".html") - 5, temp.indexOf(".html"));
          subheader = temp.substring(temp.indexOf(".html") + 7, temp.indexOf("</a>"));
          countGene = 0;
          snpsInPathway = new Vector<String>();
          if (kegg.pathways.containsKey(pathway)) {
            genesInPathway = kegg.pathways.get(pathway);
            for (int i = 0; i < genesInPathway.size(); i++) {
              if (kegg.genesToSNPs.containsKey(genesInPathway.elementAt(i))) {
                countGene++;
                snpsInGene = kegg.genesToSNPs.get(genesInPathway.elementAt(i));
                for (int j = 0; j < snpsInGene.size(); j++) {
                  HashVec.addIfAbsent(snpsInGene.elementAt(j), snpsInPathway);
                }
              }
            }
          }

          System.out.println(subheader + "\t" + snpsInPathway.size());
          pvalue = evaluateGeneset(pheno, Array.toStringArray(snpsInPathway), results, repDir,
                                   p_thresh, nrss_thresh);

          writer.print("\t" + "\t" + pathway + "\t" + kegg.pathwayLookup.get(pathway) + "\t"
                       + subheader + "\t"
                       + (kegg.pathways.containsKey(pathway) ? kegg.pathways.get(pathway).size()
                                                             : 0));
          writer.println("\t" + countGene + "\t" + snpsInPathway.size() + "\t" + pvalue);
          writer.flush();
        }
        temp = reader.readLine();
      }
      reader.close();
      writer.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + groupings_file + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + groupings_file + "\"");
      System.exit(2);
    }
  }

  public static String evaluateGeneset(String pheno, String[] markerList, String results,
                                       String repDir, double p_thresh, double nrss_thresh) {
    int count, hits;
    long time;
    float[] pvals;
    double score, travScore;
    String[] markerNames;
    boolean[] included;
    String filename;
    Hashtable<String, String> markerSetHash;
    boolean done;
    double p;
    String root;

    if (markerList.length == 0) {
      System.out.println("score = 0; permutation aborted");
      return "1.00";
    }

    markerSetHash = new Hashtable<String, String>();
    for (String element : markerList) {
      markerSetHash.put(element, "");
    }


    time = new Date().getTime();
    // markerNames = HashVec.loadFileToStringArray(results, true, new int[] {1}, false);
    markerNames = HashVec.loadFileToStringArray(repDir + "plink.bim", false, new int[] {1}, false);

    if (!new File(repDir + ext.rootOf(pheno, false) + ".markerNames").exists()) {
      System.err.println("Error - need file '" + repDir + ext.rootOf(pheno, false)
                         + ".markerNames' to verify marker order");
      return "error";
    }
    if (Array.equals(markerNames,
                     SerialStringArray.load(repDir + ext.rootOf(pheno, false) + ".markerNames",
                                            false)
                                      .getArray(),
                     false)) {
      System.out.println("Marker order matches with replicates");
    } else {
      System.err.println("Error - marker order does not match between results file and replicate files");
      return "error";
    }

    included = new boolean[markerNames.length];
    for (int i = 0; i < markerNames.length; i++) {
      included[i] = markerSetHash.containsKey(markerNames[i]);
    }

    root = ext.rootOf(results, false);
    if (root.endsWith(".assoc")) {
      root = ext.rootOf(root, false);
    }

    if (!new File(root + ".results").exists()) {
      compress(results, markerNames, false);
    }
    pvals = SerialFloatArray.load(root + ".results", false).getArray();
    score = calculateScore(pvals, included, p_thresh, nrss_thresh);
    // System.out.println("score = "+score);

    hits = 0;
    count = 0;
    done = false;
    while (!done) {
      count++;
      filename = repDir + ext.rootOf(pheno) + "." + count + ".results";
      if (count % 1000 == 0) {
        p = ((double) hits + 1) / ((double) count + 1);
        System.out.println(hits + "/" + count + "=" + ext.prettyP(p));
        if (hits > 100) {
          done = true;
        }
      }
      if (new File(filename).exists()) {
        pvals = SerialFloatArray.load(filename, false).getArray();
        travScore = calculateScore(pvals, included, p_thresh, nrss_thresh);
        if (travScore >= score) {
          hits++;
        }
      } else {
        done = true;
      }
    }
    System.out.println("Finished in " + ext.getTimeElapsed(time));
    System.out.println("score = " + score + "; final permutation: " + hits + "/" + count + "="
                       + ext.prettyP(((double) hits + 1) / ((double) count + 1)));

    return ext.prettyP(((double) hits + 1) / ((double) count + 1));
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String pheno = "pheno.dat";
    String covars = null;
    String results = "plink.qassoc";
    String repDir = "./";
    String geneset = "";
    boolean permAndRun = false;
    boolean compress = false;
    boolean parseKEGG = false;
    // String list = "";
    String ids = "";
    String gene_info = ROOT_DIRECTORY + "NCBI/NCBI/seq_gene.md";
    boolean parseGeneInfo = false;
    boolean allSets = false;
    boolean linkCovToPheno = true;
    boolean sex = false;
    String mapfile = "plink.bim";
    int offset = 15000;
    double p_threshold = DEFAULT_P_THRESHOLD;
    double nrss_index_threshold = DEFAULT_NRSS_INDEX_THRESHOLD;

    String usage = "\n" + "gwas.PathwayAnalysis requires 0-1 arguments\n"
                   + "   (1) phenotype filename (i.e. pheno=" + pheno + " (default))\n"
                   + "   (2) covariates filename (i.e. covars=" + covars
                   + " (optional) (default))\n" + "   (2b) include sex as a covariate (i.e. sex="
                   + sex + " (optional) (default))\n"
                   + "   (3) link covariates to pheno, not geno (i.e. link=" + linkCovToPheno
                   + " (default))\n"
                   + "   (3) permute pheno and run until plug is pulled (i.e. -permAndRun (not the default))\n"
                   + "   (4) compress by keeping only pvalues (i.e. -compress (not the default))\n"
                   + "   (4b) compress all existing result files (i.e. -compress (same flag; not the default))\n"
                   + "   (5) gene-set with which to evaluate (i.e. evaluate=markerList.txt (not the default))\n"
                   + "   (6) result set to compare against (i.e. results=" + results
                   + " (default))\n"
                   + "   (7) name of directory containing replicates (i.e. repDir=" + repDir
                   + " (default))\n"
                   + "   (8) parse gene location info (i.e. -parseGeneInfo (not the default))\n" +
                   // " (9) evaluate set of snps (i.e. set=list.txt (not the default))\n"+
                   "   (10) evaluate all KEGG sets (i.e. -allSets (not the default))\n"
                   + "   (11) mapfile to map to KEGG sets (i.e. map=" + mapfile + " (default))\n"
                   + "   (12) max offset from SNP position to gene start/stop (i.e. offset="
                   + offset + " (default))\n"
                   + "   (13) p-value threshold for inclusion in score (i.e. p_thresh="
                   + p_threshold + " (default))\n"
                   + "   (14) index threhsold for NRS score (i.e. nrss_thresh="
                   + nrss_index_threshold + " (default; 1=nrss_not_used))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
        // } else if (args[i].startsWith("list=")) {
        // list = args[i].split("=")[1];
        // numArgs--;
      } else if (arg.startsWith("pheno=")) {
        pheno = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("covars=")) {
        covars = ext.parseStringArg(arg, null);
        numArgs--;
      } else if (arg.startsWith("sex=")) {
        sex = ext.parseBooleanArg(arg);
        numArgs--;
      } else if (arg.startsWith("link=")) {
        linkCovToPheno = ext.parseBooleanArg(arg);
        numArgs--;
      } else if (arg.startsWith("results=")) {
        results = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("evaluate=")) {
        geneset = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("repDir=")) {
        repDir = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("-permAndRun")) {
        permAndRun = true;
        numArgs--;
      } else if (arg.startsWith("-compress")) {
        compress = true;
        numArgs--;
      } else if (arg.startsWith("-parseGeneInfo")) {
        parseGeneInfo = true;
        numArgs--;
      } else if (arg.startsWith("-allSets")) {
        allSets = true;
        numArgs--;
      } else if (arg.startsWith("mapfile=")) {
        mapfile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("offset=")) {
        offset = ext.parseIntArg(arg);
        numArgs--;
      } else if (arg.startsWith("p_thresh=")) {
        p_threshold = ext.parseDoubleArg(arg);
        numArgs--;
      } else if (arg.startsWith("nrss_thresh=")) {
        nrss_index_threshold = ext.parseDoubleArg(arg);
        numArgs--;
      } else {
        System.err.println("Error - what to do with " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    // geneset = "subset.txt";
    // results = "actual.qassoc";
    // ids = "D:\\umn\\Myron\\MethylationGrant\\AlzGeneWebsite.txt";

    // parseGeneInfo = true;

    // compress = true;
    // permAndRun = true;
    // covars = "covars.dat";
    // linkCovToPheno = true;

    // allSets = true;

    // parseSNPsInPathway("ko04740", "D:\\LOAD\\Pathway\\plink.bim", 25000);

    try {
      if (!ids.equals("")) {
        checkConcentrationFromGeneList(ids, KEGG_GROUPINGS);
        // } else if (!list.equals("")) {
        // runList(list);
      } else if (parseGeneInfo) {
        parseGeneLocs(gene_info, KEGG_GENE_POSITIONS);
      } else if (parseKEGG) {
        parseKEGG(KEGG_KO, KEGG_GROUPINGS);
      } else if (permAndRun) {
        permAndRun(pheno, covars, sex, linkCovToPheno, compress);
      } else if (compress) {
        compressAll(pheno, covars);
      } else if (allSets) {
        evaluateAllSets(pheno, mapfile, offset, results, repDir, KEGG_GROUPINGS, p_threshold,
                        nrss_index_threshold);
      } else if (!geneset.equals("")) {
        evaluateGeneset(pheno, HashVec.loadFileToStringArray(geneset, false, new int[] {0}, false),
                        results, repDir, p_threshold, nrss_index_threshold);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public static void parseGeneLocs(String filename, String output) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;

    try {
      reader = new BufferedReader(new FileReader(filename));
      writer = new PrintWriter(new FileWriter(output));
      reader.readLine();
      writer.println("GeneID\tlocation");
      while (reader.ready()) {
        line = reader.readLine().trim().split("\t", -1);
        if (line[11].equals("GENE")) {
          writer.println(line[10].substring(7) + "\tchr" + line[1] + ":" + line[2] + "-" + line[3]);
        }
      }
      writer.close();
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + filename + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + filename + "\"");
      System.exit(2);
    }
  }

  // public static void runList(String list) {
  //// BufferedReader reader;
  // PrintWriter writer;
  //// String[] line;
  //// String temp, trav;
  //// Hashtable<String,String> hash = new Hashtable<String,String>();
  //// Vector<String> v = new Vector<String>();
  //// int count;
  //// long time;
  // String[] sets;
  //
  // sets = HashVec.loadFileToStringArray(list, false, new int[] {0}, false);
  // try {
  // writer = new PrintWriter(new FileWriter(ext.rootOf(list)+"_results.xln"));
  // writer.println("Pathway\ticam_whites\tpsel_whites");
  // for (int i = 0; i<sets.length; i++) {
  // writer.println(sets[i]+"\t"+
  // evaluateGeneset("icam_whites/pheno.dat",
  // HashVec.loadFileToStringArray("pathways/"+sets[i]+".snps", false, new int[] {0}, false),
  // "icam_whites/actual.qassoc", "icam_whites/") +"\t"+
  // evaluateGeneset("psel_whites/pheno.dat",
  // HashVec.loadFileToStringArray("pathways/"+sets[i]+".snps", false, new int[] {0}, false),
  // "psel_whites/actual.qassoc", "psel_whites/")
  // );
  // writer.flush();
  // }
  // writer.close();
  // } catch (Exception e) {
  // System.err.println("Error writing to "+ext.rootOf(list)+"_results.xln");
  // e.printStackTrace();
  // }
  // }

  public static void parseKEGG(String ko_file, String groupings_file) {
    BufferedReader reader;
    PrintWriter writer;
    String[] keys;
    String temp;
    String pathway;
    String dir, main, header, subheader;
    KEGGpathway kegg;

    kegg = new KEGGpathway(ko_file, false);
    dir = ext.verifyDirFormat(ext.parseDirectoryOfFile(ko_file));

    new File(dir + "pathways/").mkdirs();
    keys = HashVec.getKeys(kegg.pathways);
    for (String key : keys) {
      Files.writeList(Array.toStringArray(kegg.pathways.get(key)),
                      dir + "pathways/" + key + ".list");
      try {
        writer = new PrintWriter(new FileWriter(dir + "pathways/" + key + ".crf"));
        writer.println("genes");
        writer.println(dir + "pathways/" + key + ".list 0 out=" + dir + "pathways/" + key
                       + ".snps");
        writer.println(dir + "plink.bim 1 0 3");
        writer.println("buffer=" + DEFAULT_BUFFER);
        writer.close();
      } catch (Exception e) {
        System.err.println("Error writing to " + dir + "pathways/" + key + ".crf");
        e.printStackTrace();
      }
      MapGenesToSNPs.filter(dir + "pathways/" + key + ".crf", new Logger());
    }

    try {
      writer = new PrintWriter(new FileWriter(dir + "genes.xln"));
      keys = HashVec.getKeys(kegg.genes);
      for (String key : keys) {
        writer.println(key + "\t" + kegg.genes.get(key).size() + "\t"
                       + Array.toStr(Array.toStringArray(kegg.genes.get(key))));
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing to " + dir + "genes.xln");
      e.printStackTrace();
    }

    try {
      reader = new BufferedReader(new FileReader(groupings_file));
      writer = new PrintWriter(new FileWriter(dir + "groupings.xln"));
      writer.println("Main Heading\tHeader\tPathway\tPathwayName\tSubHeader\t# genes");
      temp = reader.readLine();
      while (reader.ready()) {
        if (temp.startsWith("<h4>")) {
          main = temp.substring(temp.indexOf("<h4>") + 4, temp.indexOf("</h4>"));
          writer.println(main);
        } else if (temp.startsWith("<b>")) {
          header = temp.substring(temp.indexOf("<b>") + 3, temp.indexOf("</b>"));
          writer.println("\t" + header);
        } else if (temp.contains("/kegg/pathway/")) {
          pathway = "ko" + temp.substring(temp.indexOf(".html") - 5, temp.indexOf(".html"));
          subheader = temp.substring(temp.indexOf(".html") + 7, temp.indexOf("</a>"));
          writer.println("\t" + "\t" + pathway + "\t" + kegg.pathwayLookup.get(pathway) + "\t"
                         + subheader + "\t"
                         + (kegg.pathways.containsKey(pathway) ? kegg.pathways.get(pathway).size()
                                                               : 0));
        }
        temp = reader.readLine();
      }
      reader.close();
      writer.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + groupings_file + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + groupings_file + "\"");
      System.exit(2);
    }
  }

  public static void parseSNPsInPathway(String pathway, String mapfile, int offset) {
    Vector<String> snpsInPathway, genesInPathway, snpsInGene;
    KEGGpathway kegg;

    kegg = new KEGGpathway(KEGG_KO2, KEGG_KO2_TITLES);
    System.out.print("Mapping snps to genes...");
    kegg.matchMarkersToGenes(mapfile, KEGG_GENE_POSITIONS, offset);
    System.out.println("done");

    snpsInPathway = new Vector<String>();
    if (kegg.pathways.containsKey(pathway)) {
      genesInPathway = kegg.pathways.get(pathway);
      for (int i = 0; i < genesInPathway.size(); i++) {
        if (kegg.genesToSNPs.containsKey(genesInPathway.elementAt(i))) {
          snpsInGene = kegg.genesToSNPs.get(genesInPathway.elementAt(i));
          for (int j = 0; j < snpsInGene.size(); j++) {
            HashVec.addIfAbsent(snpsInGene.elementAt(j), snpsInPathway);
          }
        }
      }
      Files.writeList(Array.toStringArray(snpsInPathway),
                      ext.parseDirectoryOfFile(mapfile) + pathway + "_" + ((int) (offset / 1000.0))
                                                          + "K.dat");
    } else {
      System.err.println("Error - '" + pathway + "' is not a valid pathway");
    }
  }

  public static void permAndRun(String pheno, String covars, boolean sex, boolean linkCovToPheno,
                                boolean compress) {
    String[] markerNames, filenames;
    String root, exten, results;
    boolean quant;
    long time;
    int rep;

    quant = Array.determineType(pheno, 2, true) == 1;

    time = new Date().getTime();
    Files.writeList(new String[0], "plug");
    root = ext.rootOf(pheno, false);
    exten = pheno.substring(pheno.lastIndexOf("."));
    if (new File(root + ".markerNames").exists()) {
      markerNames = SerialStringArray.load(root + ".markerNames", false).getArray();
    } else {
      markerNames = HashVec.loadFileToStringArray("plink.bim", false, new int[] {1}, false);
      new SerialStringArray(markerNames).serialize(root + ".markerNames");
    }

    rep = 0;
    while (new File("plug").exists()) {
      if (new File("wait").exists()) {
        try {
          Thread.sleep(5000);
        } catch (Exception e) {
        }
      } else {
        Files.writeList(new String[0], "wait");
        rep = Files.findNextRep(new String[] {root + ".#" + exten, root + ".#.results"}, -1, rep);
        filenames = permutePheno(pheno, linkCovToPheno ? covars : null, rep);
        new File("wait").delete();

        if (covars == null && !sex) {
          CmdLine.run(ROOT_DIRECTORY + "bin/plink --bfile plink --noweb --assoc --pheno "
                      + filenames[0] + " --out " + root + "." + rep, "./");
          results = root + "." + rep + "." + (quant ? "q" : "") + "assoc";
        } else {
          CmdLine.run(ROOT_DIRECTORY + "bin/plink --bfile plink --noweb --"
                      + (quant ? "linear" : "logistic") + " --pheno " + filenames[0] + " --covar "
                      + (linkCovToPheno ? filenames[1] : covars) + (sex ? " --sex" : "") + " --out "
                      + root + "." + rep, "./");
          results = root + "." + rep + ".assoc." + (quant ? "linear" : "logistic");
        }
        if (!new File(results).exists()) {
          System.err.println("Fumble at rep " + rep);
          try {
            Thread.sleep(20000);
          } catch (Exception e) {
          }
          for (String filename : filenames) {
            if (new File(filename).exists()) {
              new File(filename).delete();
            }
          }
          // Files.writeList(new String[0], "error at "+filenames[0]+" unfortunately.out");
          // return;
        } else if (compress) {
          compress(results, markerNames, true);
        }
      }
    }
    System.out.println("Ran for " + ext.getTimeElapsed(time));
    System.out.println("Produced " + rep + " replicates");
  }

  public static String[] permutePheno(String pheno, String covars, int rep) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line;
    String[] trait;
    int[] keys, indices;
    String[] filenames;


    filenames = new String[2];
    filenames[0] = ext.rootOf(pheno, false) + "." + rep + pheno.substring(pheno.lastIndexOf("."));
    filenames[1] =
        ext.rootOf(pheno, false) + "_cov." + rep + pheno.substring(pheno.lastIndexOf("."));
    trait = HashVec.loadFileToStringArray(pheno, true, new int[] {2}, false);
    line = Files.getHeaderOfFile(pheno, "[\\s]+", new Logger());
    if (!line[0].equals("FID") || !line[1].equals("IID")) {
      System.err.println("Error - need to use FID/IID in both pheno and covar files");
      Files.writeList(new String[] {"Error - need to use FID/IID in both pheno and covar files"},
                      "NEED_TO_USE_FID-IID_IN_BOTH_PHENO_AND_COVAR_FILES.txt");
      new File("wait").delete();
      System.exit(1);
    }
    keys = Array.random(trait.length);
    try {
      reader = new BufferedReader(new FileReader(pheno));
      writer = new PrintWriter(new FileWriter(filenames[0]));
      writer.println(reader.readLine());
      for (int i = 0; i < trait.length; i++) {
        line = reader.readLine().trim().split("[\\s]+");
        writer.println(line[0] + "\t" + line[1] + "\t" + trait[keys[i]]);
      }
      writer.close();
      reader.close();
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + pheno + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + pheno + "\"");
      System.exit(2);
    }

    if (covars != null) {
      if (!new File("pheno_covar_checks_out").exists()) {
        if (!Array.equals(HashVec.loadFileToStringArray(pheno, true, new int[] {0, 1}, false),
                          HashVec.loadFileToStringArray(covars, true, new int[] {0, 1}, false),
                          false)) {
          System.err.println("Error - the covars files needs to match the ids in the pheno file line for line");
          Files.writeList(new String[] {"Error - the covars files needs to match the ids in the pheno file line for line"},
                          "NEED_TO_SYNC_PHENO_AND_COVAR_FILES.txt");
          new File("wait").delete();
          System.exit(1);
        }
        Files.writeList(new String[] {}, "pheno_covar_checks_out");
      }

      line = Files.getHeaderOfFile(covars, "[\\s]+", new Logger());
      indices = new int[line.length - 2];
      for (int i = 0; i < indices.length; i++) {
        indices[i] = 2 + i;
      }
      trait = HashVec.loadFileToStringArray(covars, true, indices, false);

      try {
        reader = new BufferedReader(new FileReader(covars));
        writer = new PrintWriter(new FileWriter(filenames[1]));
        writer.println(reader.readLine());
        for (int i = 0; i < trait.length; i++) {
          line = reader.readLine().trim().split("[\\s]+");
          writer.println(line[0] + "\t" + line[1] + "\t" + trait[keys[i]]);
        }
        writer.close();
        reader.close();
      } catch (FileNotFoundException fnfe) {
        System.err.println("Error: file \"" + covars + "\" not found in current directory");
        System.exit(1);
      } catch (IOException ioe) {
        System.err.println("Error reading file \"" + covars + "\"");
        System.exit(2);
      }
    }

    return filenames;
  }
}
