package org.pankratzlab.internal.gwas;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Hashtable;
import org.apache.commons.cli.ParseException;
import org.genvisis.cnv.plots.ManhattanPlot;
import org.genvisis.cnv.plots.QQPlot;
import org.pankratzlab.common.Aliases;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.CLI;
import org.pankratzlab.common.CmdLine;
import org.pankratzlab.common.Elision;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.Matrix;
import org.pankratzlab.common.ext;

public class ParseMXCResults {

  private static String usage = "";

  private static String[] getKeys(String[][] matrix, String colname, Logger log) {
    String[] headers = matrix[0];
    int index = -1;
    for (int i = 0; i < headers.length; i++)
      if (colname.equals(headers[i])) index = i;

    if (index == -1) {
      log.report("Could not find column: " + colname);
      System.exit(0);
    }

    return Matrix.extractColumn(matrix, index);
  }

  private static String[][] loadGenes(String filename, boolean omitHeader, Logger log) {

    String[] line = Files.getHeaderOfFile(filename, log);
    int[] cols = ext.indexFactors(new String[][] {new String[] {"Assembly_name", "Gene",
                                                                "SKATgene"},
                                                  Aliases.POSITIONS_START, Aliases.POSITIONS_STOP,
                                                  Aliases.CHRS},
                                  line, true, false, false, false, log);

    return HashVec.loadFileToStringMatrix(filename, omitHeader, cols);
  }

  private static String[][] removeNAs(String[][] matrix, int[] cols) {
    boolean[] rowsToKeep = new boolean[matrix.length];
    if (cols == null || cols.length == 0) return matrix;

    for (int c : cols) {
      for (int i = 0; i < matrix.length; i++) {
        rowsToKeep[i] = !(matrix[i][c].equals("NA") || matrix[i][c].equals("."));
      }
    }

    matrix = Matrix.subset(matrix, rowsToKeep);
    return matrix;
  }

  private static void addMetalHits(String posfile, String mxcfile, String metalfile,
                                   String genesFile, Logger log) {
    String outDir = ext.parseDirectoryOfFile(mxcfile) + "results_parsed/";
    if (!new File(outDir).exists()) {
      new File(outDir).mkdir();
    }

    // read in the metal gwas results file
    String[][] metal = null;

    String[] line = Files.getHeaderOfFile(metalfile, log);

    int[] valueIndices = ext.indexFactors(new String[][] {Aliases.MARKER_NAMES, Aliases.PVALUES},
                                          line, false, true, false);

    metal = HashVec.loadFileToStringMatrix(metalfile, true, valueIndices);

    float bf = (float) 0.05 / metal.length;

    String[][] positions = HashVec.loadFileToStringMatrix(posfile, false, null);

    String[] keys = Matrix.extractColumn(metal, 0);
    String[][] results = null;

    // combine snps, pvals, and positions;
    try {
      results = Files.combineInMemory(keys, positions, "NA", true, true, log);
    } catch (Elision e) {
      log.reportError(e.getMessage());
    }

    // metal should now be of the form MarkerName, P-value, Chr, Position
    metal = ArrayUtils.append(metal, results);
    metal = removeNAs(metal, new int[] {2, 3});
    Arrays.sort(metal, new Comparator<String[]>() {

      @Override
      public int compare(final String[] entry1, final String[] entry2) {
        final String chr1 = entry1[2];
        final String chr2 = entry2[2];
        if (chr1.equals(chr2)) {
          final int pos1 = Integer.parseInt(entry1[3]);
          final int pos2 = Integer.parseInt(entry2[3]);
          return pos1 - pos2;
        }
        return chr1.compareTo(chr2);
      }
    });

    String[][] mxc = HashVec.loadFileToStringMatrix(mxcfile, false, null);

    // get position ranges for genes
    String[][] genePositions = loadGenes(genesFile, true, log);
    genePositions = removeNAs(genePositions, new int[] {1});

    Arrays.sort(genePositions, new Comparator<String[]>() {

      @Override
      public int compare(final String[] entry1, final String[] entry2) {
        final String chr1 = entry1[3];
        final String chr2 = entry2[3];
        if (chr1.equals(chr2)) {
          final int pos1 = Integer.parseInt(entry1[1]);
          final int pos2 = Integer.parseInt(entry2[1]);
          return pos1 - pos2;
        }
        return chr1.compareTo(chr2);
      }
    });

    // maps gene to related snps
    String[] mxcGenes = getKeys(mxc, "gene_name", log);
    int startPos;
    int endPos;
    int snpPos;
    double pval;
    String chr, snpChr;
    String[] snp, g;
    int start = 1;
    String[][] sig = new String[genePositions.length][];

    for (int gene = 0; gene < genePositions.length; gene++) {
      g = genePositions[gene];
      startPos = Integer.parseInt(g[1]) < 250000 ? 0 : Integer.parseInt(g[1]) - 250000;
      endPos = Integer.parseInt(g[2]) + 250000;
      chr = g[3];

      sig[gene] = new String[] {g[0], "0", "0", g[1], g[3]};

      while (start < metal.length && !metal[start][2].equals(chr)) {
        start++;
      }

      for (int j = start; j < metal.length; j++) {
        snp = metal[j];
        pval = Double.parseDouble(snp[1]);
        snpChr = snp[2];
        snpPos = Integer.parseInt(snp[3]);

        if (!snpChr.equals(chr)) break;

        // genes are sorted by start position, so if we're not there yet, we can move it up
        if (snpPos < startPos) {
          start = j;
          continue;
        }

        if (snpPos > endPos) break;

        if (pval < bf * 100) {
          int[] numSig = new int[] {Integer.parseInt(sig[gene][1]), Integer.parseInt(sig[gene][2])};
          if (pval < bf) sig[gene][1] = numSig[0] + 1 + "";

          sig[gene][2] = numSig[1] + 1 + "";
        }
      }
    }

    // append num of sig snps to mxc file
    try {
      results = Files.combineInMemory(mxcGenes, sig, "NA", true, true, log);
      results[0] = new String[] {"numSig", "numSug", "pos", "chr"};
    } catch (Elision e) {
      log.reportError(e.getMessage());
    }

    mxc = ArrayUtils.append(mxc, results);
    // write to mxc file
    Files.writeMatrix(mxc, outDir + ext.removeDirectoryInfo(mxcfile), ",");

    try {
      BufferedWriter writer = new BufferedWriter(new FileWriter(outDir + "no_sig_markers_"
                                                                + ext.removeDirectoryInfo(mxcfile)));
      writer.write(ArrayUtils.toStr(mxc[0], ",") + "\n");
      for (String[] s : mxc) {
        if (s[s.length - 4].equals("0")) writer.write(ArrayUtils.toStr(s, ",") + "\n");
      }
      writer.close();
    } catch (IOException e) {
      log.reportError("Unable to write to " + "no_sig_markers_" + ext.removeDirectoryInfo(mxcfile));
      e.printStackTrace();
    }

    mplot(outDir + "no_sig_markers_" + ext.removeDirectoryInfo(mxcfile),
          outDir + ext.rootOf(mxcfile) + ".png", log);
    qqplot(mxcfile, outDir + ext.rootOf(mxcfile) + ".png", log);

    // generateRScript(outDir + ext.rootOf(mxcfile) + "_plot.R");
  }

  private static void run(String data, String db, String posmap, String covar, String out,
                          boolean overwrite, String mxcFolder, Logger log) {
    String gwas_folder = ext.parseDirectoryOfFile(new File(data).getAbsolutePath());

    db = new File(db).getAbsolutePath();
    covar = new File(covar).getAbsolutePath();

    out = new File(out).getAbsolutePath();

    String py = "MetaXcan.py";

    String delimiter = Files.determineDelimiter(data, log);
    String[] header = Files.getHeaderOfFile(data, delimiter, log);
    int[] index = ext.indexFactors(new String[][] {Aliases.MARKER_NAMES, Aliases.REF_ALLELES,
                                                   Aliases.ALT_ALLELES, Aliases.EFFECTS,
                                                   Aliases.STD_ERRS, Aliases.PVALUES},
                                   header, false, true, false);

    if (index[0] == -1) {
      log.reportError("Unable to find SNP column. Aborting.");
      System.exit(1);
    } else if (index[1] == -1 || index[2] == -1) {
      log.reportError("Unable to find allele columns. Aborting.");
      System.exit(1);
    } else if ((index[5] == -1 && index[4] == -1) || index[3] == -1) {
      log.reportError("Unable to find necessary pval or effect columns. Beta column and either P-value or SE columns are required. Aborting.");
      System.exit(1);
    }

    String markerColumn = header[index[0]];
    String a1 = header[index[1]];
    String a2 = header[index[2]];
    String effect = " --beta_column " + header[index[3]]
                    + (index[5] == -1 ? " --se_column " + header[index[4]]
                                      : " --pvalue_column " + header[index[5]]);

    // build mxc command
    String command = "./" + py + " --model_db_path " + db + " --covariance " + covar
                     + " --gwas_folder " + gwas_folder + " --gwas_file_pattern "
                     + new File(data).getName() + " --output_file " + out
                     + " --effect_allele_column " + a1 + " --non_effect_allele_column " + a2
                     + " --snp_column " + markerColumn + effect + (overwrite ? " --overwrite" : "");

    log.report(command);

    // run MetaXcan on the given inputs
    boolean runSuccess = CmdLine.run(command, ext.parseDirectoryOfFile(mxcFolder), null, null,
                                     new Logger(new File("").getAbsolutePath() + "MetaXcan.log"),
                                     false);

    if (!runSuccess || !new File(out).exists()) {
      log.reportError("Encountered a problem running MetaXcan with the given inputs.");
      System.exit(1);
    }
  }

  private static void qqplot(String filename, String out, Logger log) {
    String[] header = Files.getHeaderOfFile(filename, log);
    int[] cols = ext.indexFactors(new String[] {"pvalue"}, header, false);
    String[] data = HashVec.loadFileToStringArray(filename, true, cols, false);
    data = ArrayUtils.removeMissingValues(data);

    try {
      QQPlot qq = new QQPlot(new String[] {ext.rootOf(filename)},
                             new double[][] {ArrayUtils.toDoubleArray(data)},
                             new boolean[][] {ArrayUtils.booleanArray(data.length, true)}, true,
                             false, false, Float.MAX_VALUE, log);
      qq.screenCap(ext.addToRoot(out, "_qq"));
    } catch (Exception e) {

    }
  }

  private static void mplot(String filename, String out, Logger log) {

    try {
      ManhattanPlot mp = new ManhattanPlot(null);
      boolean load = mp.loadFileAuto(filename);
      if (!load) {
        log.reportError("Unable to load " + filename + " for manhattan plot.");
        return;
      }
      while (!mp.isDataLoaded()) {
        Thread.sleep(200);
      }

      mp.getManPan().setSize(800, 400);
      mp.screenshot(ext.addToRoot(out, "_manhattan"));
    } catch (InterruptedException e) {
      log.reportError("Problem creating manhattan plot from " + filename);
      e.printStackTrace();
    }

  }

  private static void combine(String pattern, Logger log) {
    String[] files = Files.list(ext.parseDirectoryOfFile(pattern), ext.rootOf(pattern), ".csv",
                                true, true);

    String[][] data, genes = null;
    String[] header, keys;
    int[] cols;
    String[][] combined = null;
    String label;
    for (String s : files) {
      label = s.split(pattern)[1];
      header = Files.getHeaderOfFile(s, log);

      cols = ext.indexFactors(new String[] {"gene_name", "pvalue", "effect_size", "numSig",
                                            "numSug"},
                              header, false);

      data = HashVec.loadFileToStringMatrix(s, false, cols);
      data[0] = new String[] {"Gene", "p_" + label, "beta_" + label, "sig_" + label,
                              "sugg_" + label};

      if (genes == null || data.length > genes.length) {
        cols = ext.indexFactors(new String[][] {new String[] {"gene_name"}, Aliases.CHRS,
                                                Aliases.POSITIONS},
                                header, false, true, false);

        genes = HashVec.loadFileToStringMatrix(s, false, cols);
        genes[0] = new String[] {"Gene", "Chr", "pos"};

        // we want our keys to be the longest set of genes so far
        // so if we have an existing list, swap it with data before combining
        if (combined != null) {
          String[][] temp = combined;
          combined = data;
          data = temp;
        }
      }

      if (combined == null) {
        combined = data;
      } else {
        try {
          keys = getKeys(combined, "Gene", log);
          combined = ArrayUtils.append(combined,
                                       Files.combineInMemory(keys, data, ".", true, true, log));
        } catch (Elision e) {
          log.reportError("Something went wrong while combining data from file " + s);
          e.printStackTrace();
        }
      }
    }

    try {
      keys = getKeys(genes, "Gene", log);
      combined = ArrayUtils.append(genes,
                                   Files.combineInMemory(keys, combined, ".", true, true, log));
    } catch (Elision e) {
      log.reportError("Unable to map combined files to genes.");
      e.printStackTrace();
    }

    Matrix.writeToFile(combined, pattern + "combined.txt");
  }

  private static void index(String indicesFile, String mxcFile, int range, Logger log) {
    // read in indices from which to pull data
    String[] header = Files.getHeaderOfFile(indicesFile, log);
    int[] cols = ext.indexFactors(new String[][] {Aliases.GENE_UNITS, Aliases.POSITIONS,
                                                  Aliases.CHRS},
                                  header, false, true, false);
    if (ArrayUtils.min(cols) == -1) {
      log.reportError("Expected " + indicesFile
                      + " to be of the form: Gene Label, Position, Chromosome");
      System.exit(1);
    }

    String[][] indices = HashVec.loadFileToStringMatrix(indicesFile, true, cols);
    indices = removeNAs(indices, new int[] {1, 2});
    Arrays.sort(indices, new Comparator<String[]>() {

      @Override
      public int compare(final String[] entry1, final String[] entry2) {
        final int chr1 = Integer.parseInt(entry1[2]);
        final int chr2 = Integer.parseInt(entry2[2]);
        if (chr1 == chr2) {
          final int pos1 = Integer.parseInt(entry1[1]);
          final int pos2 = Integer.parseInt(entry2[1]);
          return pos1 - pos2;
        }
        return chr1 - chr2;
      }
    });

    // read in mxc file and sort it by chr/pos
    header = Files.getHeaderOfFile(mxcFile, log);
    cols = ext.indexFactors(new String[][] {new String[] {"gene_name"}, Aliases.POSITIONS,
                                            Aliases.CHRS, Aliases.PVALUES},
                            header, false, true, false);

    if (ArrayUtils.min(cols) == -1) {
      log.reportError("Invalid header for " + mxcFile);
      System.exit(1);
    }

    String[][] mxc = HashVec.loadFileToStringMatrix(mxcFile, true, cols);
    mxc = removeNAs(mxc, new int[] {1, 2});
    Arrays.sort(mxc, new Comparator<String[]>() {

      @Override
      public int compare(final String[] entry1, final String[] entry2) {
        final int chr1 = Integer.parseInt(entry1[2]);
        final int chr2 = Integer.parseInt(entry2[2]);
        if (chr1 == chr2) {
          final int pos1 = Integer.parseInt(entry1[1]);
          final int pos2 = Integer.parseInt(entry2[1]);
          return pos1 - pos2;
        }
        return chr1 - chr2;
      }
    });

    Hashtable<String, ArrayList<String[]>> regions = new Hashtable<>();
    int start = 0;
    // for each index, generate a list of mxc genes in the region
    // calculate the lambda for this region
    // calculate the bf threshold
    for (String[] i : indices) {
      String label = i[0];
      int startPos = Integer.max(Integer.parseInt(i[1]) - range, 0);
      int endPos = Integer.parseInt(i[1]) + range;
      int chr = Integer.parseInt(i[2]);

      for (int j = start; j < mxc.length; j++) {
        // check if this gene is in range
        int c = Integer.parseInt(mxc[j][2]);
        int pos = Integer.parseInt(mxc[j][1]);
        if (c < chr) continue;
        else if (c > chr || pos > endPos) break;

        if (pos < startPos) continue;

        // if this is the first one, add it and update the start pos
        if (!regions.containsKey(label)) {
          regions.put(label, new ArrayList<String[]>());
          start = j;
        }

        regions.get(label).add(mxc[j]);
      }
    }

    try {
      BufferedWriter out = new BufferedWriter(new FileWriter(ext.rootOf(indicesFile)
                                                             + "_regions.txt"));

      out.write("Label\tGene\tPos\tChr\tPvalue\tLambda\tSig Threshold\n");
      for (String key : regions.keySet()) {
        ArrayList<String[]> g = regions.get(key);
        double[] pvals = new double[g.size()];
        for (int j = 0; j < g.size(); j++) {
          pvals[j] = Double.parseDouble(g.get(j)[3]);
        }

        double lambda = ArrayUtils.lambda(pvals);
        double bf = 0.05 / pvals.length;
        for (String[] s : g) {
          out.write(key + "\t" + ArrayUtils.toStr(s) + "\t" + lambda + "\t" + bf + "\n");
        }
      }
      out.close();
    } catch (IOException e) {
      log.report("Error writing to " + ext.rootOf(indicesFile) + "_regions.txt");
    }
  }

  public static void main(String[] args) throws IOException {
    if (args.length == 0) {
      System.out.println(usage);
      System.exit(1);
    }
    Logger log = new Logger();

    CLI cli = new CLI("MetaXcan");

    cli.addArgWithDefault("data",
                          "File or folder of summary data to be analyzed (data=Metal_results.tbl (default))",
                          "Metal_results.tbl");
    cli.addArgWithDefault("posmap", "Map file containing chr/pos for included SNPs",
                          "/panfs/roc/groups/5/pankrat2/mstimson/parkinsons/data/1000G_PD.map");
    cli.addArgWithDefault("db",
                          "MetaXcan weights database (db=DGN-HapMap-2015/DGN-WB_0.5.db (default)",
                          "DGN-HapMap-2015/DGN-WB_0.5.db");
    cli.addArgWithDefault("covar", "path to MetaXcan covariate file",
                          "/panfs/roc/groups/5/pankrat2/shared/bin/MetaXcan/mxc_presets/DGN-HapMap-2015/DGN-WB.txt.gz");
    cli.addArgWithDefault("genes", "Path to build containing gene positions",
                          "/home/pankrat2/public/bin/NCBI/genes37.xln");
    cli.addArgWithDefault("mxc", "path to folder containing MetaXcan executable",
                          "/panfs/roc/groups/5/pankrat2/shared/bin/MetaXcan/software/");
    cli.addArgWithDefault("freq", "Allelic frequency file", "freq.tbl");
    cli.addArgWithDefault("ref", "Reference file for allele verification (optional)", "1000G.xln");
    cli.addArgWithDefault("indexFile",
                          "File to be used with -index of the form Gene Label / Pos / Chr",
                          "index.txt");
    cli.addArgWithDefault("range", "Size of range to search around an index with -index", 250000);
    cli.addArg("pattern", "File name pattern to be used with -combine");
    cli.addFlag("overwrite", "Overwrite existing MetaXcan output");
    cli.addFlag("verify", "Verify allele order and strand before running MetaXcan");
    cli.addFlag("combine", "Combine multiple results files");
    cli.addFlag("index", "Create an index file based on provided regions");
    cli.addArg("out", "Path to output file");

    try {
      cli.parse(args);
    } catch (ParseException e) {
      log.reportError("Problem parsing command line arguments");
      e.printStackTrace();
    }

    String data = cli.get("data");
    String posmap = cli.get("posmap");
    String out = cli.get("out");
    String db = cli.get("db");
    String covar = cli.get("covar");
    String mxcFolder = cli.get("mxc");
    String genesFile = cli.get("genes");
    String freqFile = cli.get("freq");
    String refFile = cli.get("ref");
    String indexFile = cli.get("indexFile");
    int range = cli.getI("range");
    String pattern = cli.get("pattern");
    boolean overwrite = cli.has("overwrite");
    boolean verify = cli.has("verify");
    boolean combine = cli.has("combine");
    boolean index = cli.has("index");

    if (combine && pattern != null) {
      combine(pattern, new Logger(ext.parseDirectoryOfFile(pattern) + "combine.log"));
    } else if (index) {
      index(indexFile, data, range, new Logger("index.log"));
    } else {
      if (verify) {
        AlleleVerification.verifyAlleles(data, refFile, freqFile, posmap, false, log);
        data = ext.rootOf(data, false) + "_allele_verified.txt";
      }

      run(data, db, posmap, covar, out, overwrite, mxcFolder, log);
      // take the output mxc file and find the number of hits for each gene range
      addMetalHits(posmap, out, data, genesFile,
                   new Logger(ext.parseDirectoryOfFile(data) + "parseMXC.log"));
    }

  }
}
