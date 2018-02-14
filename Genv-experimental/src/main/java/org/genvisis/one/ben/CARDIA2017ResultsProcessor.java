package org.genvisis.one.ben;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.PSF;
import org.genvisis.common.ext;

public class CARDIA2017ResultsProcessor {

  String[] IN_HDR = {"name", "A1", "A2", "Freq1", "MAF", "Quality", "Rsq", // 6
                     "n", "Mean_predictor_allele", // 8 ----- EffAF??
                     "beta_SNP_add", "sebeta_SNP_add", // 10
                     "beta_SNP_DPA", "sebeta_SNP_DPA", // 12
                     "cov_SNP_int_SNP_DPA", "chi2_SNP" // 14
  };
  String[] OUT_HDR = {"chr", "pos", "marker_name", "strand", "base_allele", "effect_allele", "N", // IN_7
                      "effect_allele_freq", "imputation_type", "Imputation_value", // IN_6
                      "beta_main", // IN_9
                      "se_main", // IN_10
                      "beta_int", // IN_11
                      "se_int", // IN_12
                      "cov", // IN_13
                      "chi_2df", // IN_14
                      "chi_P_2df",};
  String STRAND = "+";
  String TYPE = "1";
  String dir = "/scratch.global/cole0482/CARDIA/data/";
  static final String EA = "EA/";
  static final String AA = "AA/";
  String[][] MODEL_NUTR = {{"FEV1", "DHA"}, {"FEV1", "DPA"}, {"FEV1", "EPA"}, {"FEV1", "Fish"},
                           {"FVC", "DHA"}, {"FVC", "DPA"}, {"FVC", "EPA"}, {"FVC", "Fish"}};

  int IMPUTATION_QUALITY_IN_INDEX = 5;

  static final class SNP {

    String name;
    int chr;
    int pos;
    String base;
    String eff;

    @Override
    public int hashCode() {
      final int prime = 31;
      int result = 1;
      result = prime * result + ((name == null) ? 0 : name.hashCode());
      return result;
    }

    @Override
    public boolean equals(Object obj) {
      if (this == obj) return true;
      if (obj == null) return false;
      if (getClass() != obj.getClass()) return false;
      SNP other = (SNP) obj;
      if (name == null) {
        if (other.name != null) return false;
      } else if (!name.equals(other.name)) return false;
      return true;
    }

  }

  private String getAnc(String prefix) {
    if (AA.equals(prefix)) {
      return "AA";
    } else {
      return "Cauc";
    }
  }

  private void process(String dir, String prefix) throws IOException {
    BufferedReader reader;
    PrintWriter writer;
    PrintWriter missWriter;
    String modelPref = "model_";

    int numChrs = EA.equals(prefix) ? 23 : 24;

    String form = "ddMMMyy";
    SimpleDateFormat sdf = new SimpleDateFormat(form);
    String date = sdf.format(new Date()).toUpperCase();

    String skipFile = this.dir + prefix + "removeSnpsLowQual.txt";
    HashSet<String> ignore = Files.exists(skipFile) ? HashVec.loadFileToHashSet(skipFile, false)
                                                    : new HashSet<String>();

    HashMap<Integer, PrintWriter> writers = new HashMap<>();
    HashMap<Integer, PrintWriter> missWriters = new HashMap<>();
    ChiSquaredDistribution csd = new ChiSquaredDistribution(2);

    for (int chr = 1; chr < numChrs; chr++) {
      HashMap<String, SNP> info = processChrFile(chr, prefix);
      for (int i = 0; i < MODEL_NUTR.length; i++) {
        String modelDir = dir + prefix + modelPref + (i + 1) + "/";
        String regFile = modelDir + "regression_chr" + chr + ".out_add.out.txt";
        String outFile = modelDir + "out/" + "CARDIA_" + MODEL_NUTR[i][0] + "_" + getAnc(prefix)
                         + "_" + MODEL_NUTR[i][1] + "_" + (i + 1) + "_" + date + ".txt";

        reader = Files.getAppropriateReader(regFile);
        writer = writers.get(i);
        if (writer == null) {
          writer = Files.getAppropriateWriter(outFile);
          writers.put(i, writer);
        }
        missWriter = missWriters.get(i);
        if (missWriter == null) {
          missWriter = Files.getAppropriateWriter(ext.rootOf(outFile, false) + "_missing.txt");
          missWriters.put(i, missWriter);
        }
        if (chr == 1) {
          writer.println(ArrayUtils.toStr(OUT_HDR));
        }
        reader.readLine();
        String line;
        String[] parts;
        while ((line = reader.readLine()) != null) {
          parts = line.split(PSF.Regex.GREEDY_WHITESPACE);
          StringBuilder sb = new StringBuilder();
          SNP snp = info.get(parts[0]);
          if (snp == null) {
            System.err.println("Error - missing SNP INFO for SNP: " + parts[0] + " in chr " + chr
                               + " results file for EA.");
            continue;
          }
          if (ignore.contains(parts[0])) {
            continue;
          }
          if (snp.pos == -1 || snp.chr == -1) {
            missWriter.println(parts[0]);
          }

          double chi = ext.isMissingValue(parts[14]) ? Double.NaN : Double.parseDouble(parts[14]);
          double pval = 1 - csd.cumulativeProbability(chi);
          String pvalStr = Double.isNaN(pval) ? "NA"
                                              : pval > 0.001 ? ext.formDeci(pval, 4)
                                                             : ext.formSciNot(pval, 4, false);
          String eaf = ext.isMissingValue(parts[8]) ? "NA"
                                                    : ext.formDeci(Double.parseDouble(parts[8]), 4);

          sb.append(snp.chr).append("\t").append(snp.pos).append("\t").append(parts[0]).append("\t")
            .append(STRAND).append("\t");
          if (EA.equals(prefix)) {
            sb.append(snp.base.charAt(0)).append("\t").append(snp.eff.charAt(0)).append("\t");
          } else {
            sb.append(snp.eff.charAt(0)).append("\t").append(snp.base.charAt(0)).append("\t");
          }
          sb.append(parts[7]).append("\t").append(eaf).append("\t").append(TYPE).append("\t")
            .append(parts[IMPUTATION_QUALITY_IN_INDEX]).append("\t")
            .append(ext.isMissingValue(parts[9]) ? "NA"
                                                 : ext.formDeci(Double.parseDouble(parts[9]), 5))
            .append("\t")
            .append(ext.isMissingValue(parts[10]) ? "NA"
                                                  : ext.formDeci(Double.parseDouble(parts[10]), 5))
            .append("\t")
            .append(ext.isMissingValue(parts[11]) ? "NA"
                                                  : ext.formDeci(Double.parseDouble(parts[11]), 5))
            .append("\t")
            .append(ext.isMissingValue(parts[12]) ? "NA"
                                                  : ext.formDeci(Double.parseDouble(parts[12]), 5))
            .append("\t").append(parts[13]).append("\t").append(parts[14]).append("\t")
            .append(pvalStr);

          writer.println(sb.toString());
        }
        reader.close();
      }
    }

    for (PrintWriter wr : writers.values()) {
      wr.flush();
      wr.close();
    }
    for (PrintWriter wr : missWriters.values()) {
      wr.flush();
      wr.close();
    }

  }

  String mapDirAA = "/scratch.global/cole0482/CARDIA_DATA/AA/src/map_1000G_mach2dat_qtl/";
  String mapDirEA = dir + EA + "map/";

  private HashMap<String, SNP> processChrFile(int chr, String prefix) {
    String infoFile = dir + prefix + "chr" + chr + ".info";
    boolean isEA = EA.equals(prefix);
    String mapFile = (isEA ? mapDirEA : mapDirAA) + "chr" + chr + ".map";

    HashMap<String, SNP> mapMap = new HashMap<>();
    HashMap<String, SNP> map = new HashMap<>();

    String[][] matr = HashVec.loadFileToStringMatrix(mapFile, true, new int[] {0, 1, 2, 3, 4});
    for (String[] line : matr) {
      SNP snp = new SNP();
      snp.name = isEA ? line[0] : line[2];
      snp.chr = "X".equals(isEA ? line[1] : line[0]) ? 23
                                                     : Integer.parseInt(isEA ? line[1] : line[0]);
      snp.pos = Integer.parseInt(isEA ? line[2] : line[1]);
      snp.base = line[3];
      snp.eff = line[4];
      mapMap.put(snp.name, snp);
    }

    matr = HashVec.loadFileToStringMatrix(infoFile, true, new int[] {0, 1, 2, 3, 4});
    for (String[] line : matr) {
      SNP snp = mapMap.get(line[0]);
      if (snp == null) {
        System.err.println("Error - snp found in .info file that wasn't found in the map file: "
                           + line[0]);
        snp = new SNP();
        snp.name = line[0];
        snp.chr = chr;
        snp.pos = -1;
      }
      snp.base = line[1];
      snp.eff = line[2];
      map.put(line[0], snp);
    }

    return map;
  }

  public static void combineChrXDose(String dir) throws IOException {
    // String dir = "/scratch.global/cole0482/CARDIA_DATA/AA/";
    String femD = "chrX.female.dose";
    String malD = "chrX.male.dose";
    String outD = "chrX.dose";
    // String ids = "ids.txt";
    String ids = "AA.ids";

    PrintWriter writerD;
    BufferedReader readerF;
    BufferedReader readerM;

    String[] idsArr = HashVec.loadFileToStringArray(dir + ids, false, null, false);
    HashMap<String, Integer> idIndexLookup = new HashMap<>();
    for (int i = 0; i < idsArr.length; i++) {
      idIndexLookup.put(idsArr[i], i);
    }
    boolean[] found = ArrayUtils.booleanArray(idsArr.length, false);
    ArrayList<String> outLines = new ArrayList<>((int) (idsArr.length * 1.8));
    for (int i = 0; i < idsArr.length; i++) {
      outLines.add(null);
    }

    readerF = Files.getAppropriateReader(dir + femD);
    // no header
    String line;
    while ((line = readerF.readLine()) != null) {
      String id = line.substring(0, line.indexOf(' ')).replace("->", "\t");
      Integer ind = idIndexLookup.get(id);
      if (ind == null) {
        System.err.println("Error - ID found in female dosage data that wasn't present in the IDs file: "
                           + id);
      } else {
        outLines.set(ind, line);
        found[ind] = true;
      }
    }
    readerF.close();

    readerM = Files.getAppropriateReader(dir + malD);
    // no header
    while ((line = readerM.readLine()) != null) {
      String id = line.substring(0, line.indexOf(' ')).replace("->", "\t");
      Integer ind = idIndexLookup.get(id);
      if (ind == null) {
        System.err.println("Error - ID found in male dosage data that wasn't present in the IDs file: "
                           + id);
      } else {
        outLines.set(ind, line);
        found[ind] = true;
      }
    }
    readerM.close();

    writerD = Files.getAppropriateWriter(dir + outD);
    for (int i = 0; i < idsArr.length; i++) {
      if (!found[i]) {
        System.err.println("Error - missing data for ID " + idsArr[i] + " on line " + i);
      } else {
        writerD.println(outLines.get(i));
      }
    }
    writerD.flush();
    writerD.close();

  }

  public static void combineChrXInfo(String dir) throws IOException {
    // String dir = "/scratch.global/cole0482/CARDIA_DATA/AA/";

    String femI = "chrX.female.info";
    String malI = "chrX.male.info";
    String outI = "chrX.info";

    PrintWriter writerI;

    BufferedReader readerMI = Files.getAppropriateReader(dir + malI);
    BufferedReader readerFI = Files.getAppropriateReader(dir + femI);
    writerI = Files.getAppropriateWriter(dir + outI);

    String lineM, lineF;
    readerFI.readLine();
    writerI.println(readerMI.readLine());
    while ((lineM = readerMI.readLine()) != null && (lineF = readerFI.readLine()) != null) {
      String[] partsM = lineM.split(PSF.Regex.GREEDY_WHITESPACE);
      String[] partsF = lineF.split(PSF.Regex.GREEDY_WHITESPACE);
      StringBuilder newLine = new StringBuilder();
      if (!partsM[0].equals(partsF[0])) {
        System.err.println("Error - mismatched markers: " + partsM[0] + " / " + partsF[0]);
      }

      String mkr = partsF[0];
      String a1M = partsM[1];
      String a1F = partsF[1];
      String a2M = partsM[2];
      String a2F = partsF[2];
      String freq1MStr = partsM[3];
      String freq1FStr = partsF[3];
      String mafMStr = partsM[4];
      String mafFStr = partsF[4];
      String nMStr = partsM[5];
      String nFStr = partsF[5];
      String rSqMStr = partsM[6];
      String rSqFStr = partsF[6];

      Double freq1M = Double.parseDouble(freq1MStr);
      Double freq1F = Double.parseDouble(freq1FStr);
      Double mafM = Double.parseDouble(mafMStr);
      Double mafF = Double.parseDouble(mafFStr);
      int nM = Integer.parseInt(nMStr);
      int nF = Integer.parseInt(nFStr);
      Double rSqM = "-".equals(rSqMStr) ? Double.NaN : Double.parseDouble(rSqMStr);
      Double rSqF = "-".equals(rSqFStr) ? Double.NaN : Double.parseDouble(rSqFStr);

      newLine.append(mkr).append("\t");
      if (a1M.equals(a1F)) {
        newLine.append(a1M).append("\t");
        newLine.append(a2M).append("\t");
        newLine.append(combine(freq1M, freq1F, nM, nF)).append("\t");
        newLine.append(combine(mafM, mafF, nM, nF)).append("\t");
        newLine.append(nM + nF).append("\t");
        newLine.append(combine(rSqM, rSqF, nM, nF));
      } else if (a1M.equals(a2F)) {
        newLine.append(a1F).append("\t");
        newLine.append(a2F).append("\t");
        newLine.append(combine(1 - freq1M, freq1F, nM, nF)).append("\t");
        newLine.append(combine(1 - mafM, mafF, nM, nF)).append("\t");
        newLine.append(nM + nF).append("\t");
        newLine.append(combine(rSqM, rSqF, nM, nF));
      } else {
        System.err.println("Error - mismatched alleles for marker " + mkr + " Male: " + a1M + "/"
                           + a2M + " | Female: " + a1F + "/" + a2F);
      }
      writerI.println(newLine.toString());
    }
    writerI.flush();
    writerI.close();
    readerMI.close();
    readerFI.close();
  }

  private static String combine(double m, double f, int mN, int fN) {
    if (Double.isNaN(m)) {
      if (Double.isNaN(f)) {
        return "-";
      } else {
        return ext.formDeci(f, 4);
      }
    }
    if (Double.isNaN(f)) {
      if (Double.isNaN(m)) {
        return "-";
      } else {
        return ext.formDeci(m, 4);
      }
    }
    return ext.formDeci(((m * mN) + (f * fN)) / (mN + fN), 4);
  }

  public static void main(String[] args) {
    int numArgs = args.length;

    String usage = "FAILED";
    String dir = "./";
    String pop = null;

    for (int i = 0; i < args.length; i++) {
      if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h")
          || args[i].equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (args[i].startsWith("pop=")) {
        pop = args[i].endsWith("EA") ? EA : AA;
        numArgs--;
      } else if (args[i].startsWith("dir=")) {
        dir = args[i].split("=")[1];
        numArgs--;
      } else {
        System.err.println("Error - invalid argument: " + args[i]);
      }
    }
    if (numArgs != 0 || pop == null) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      new CARDIA2017ResultsProcessor().process(dir, pop);

    } catch (Exception e) {
      e.printStackTrace();
    }
  }

}
