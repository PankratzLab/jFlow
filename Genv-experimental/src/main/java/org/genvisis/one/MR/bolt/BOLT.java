package org.genvisis.one.MR.bolt;

import java.io.File;
import java.util.ArrayList;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.qsub.Qsub;

public class BOLT {

  // FIXME: Assumes all covariates are quantitative
  private static String boltCommand(String phenoFile, String fam, String out, String[] covars) {
    String covarCommand = "";
    for (int i = 3; i < covars.length; i++) {
      covarCommand += "--qCovarCol=" + covars[i] + " \\\n";
    }
    String command = "cd /panfs/roc/groups/5/pankrat2/shared/ukbb/BOLT-LMM_v2.3/ \n" + "./bolt \\\n"
                     + "--bed=/panfs/roc/scratch/bb/all/EGAD00010001226/001/ukb_cal_chr{1:22}_v2.bed \\\n"
                     + "--bim=/panfs/roc/scratch/bb/all/EGAD00010001226/001/ukb_snp_chr{1:22}_v2.bim \\\n"
                     + "--fam=" + fam + " \\\n"
                     + "--remove=/panfs/roc/groups/5/pankrat2/shared/ukbb/bolt_results/bolt.in_plink_but_not_imputed.FID_IID.976.txt \\\n"
                     + "--phenoFile=" + phenoFile + " \\\n" + "--phenoCol=" + covars[2] + " \\\n"
                     + "--covarFile=" + phenoFile + " \\\n" + "--covarMaxLevels=30 \\\n"
                     + covarCommand
                     + "--LDscoresFile=/panfs/roc/groups/5/pankrat2/shared/ukbb/BOLT-LMM_v2.3/tables/LDSCORE.1000G_EUR.tab.gz \\\n"
                     + "--geneticMapFile=/panfs/roc/groups/5/pankrat2/shared/ukbb/BOLT-LMM_v2.3/tables/genetic_map_hg19_withX.txt.gz \\\n"
                     + "--lmmForceNonInf \\\n" + "--numThreads=8 \\\n" + "--statsFile=" + out
                     + ".stats.gz \\\n"
                     + "--bgenFile=/panfs/roc/scratch/bb/all/EGAD00010001225/001/ukb_imp_chr{1:22}_v2.bgen \\\n"
                     + "--bgenMinMAF=1e-3 \\\n" + "--bgenMinINFO=0.3 \\\n"
                     + "--sampleFile=/panfs/roc/scratch/bb/all/EGAD00010001225/001/ukb1773_imp_chr1_v2_s487398.sample \\\n"
                     + "--statsFileBgenSnps=" + out + ".bgen.stats.gz \\\n" + "--verboseStats";

    return command;
  }

  // let's just extract the phenotype data in here, shall we?
  private static String extractPhenos(String filename, String[][] fields, String caseCode) {
    // phenotype array is of the form "field id", "common name"

    String[] header = Files.getHeaderOfFile(filename, new Logger());
    ArrayList<String> finalHeader = new ArrayList<String>();
    ArrayList<String> cols = new ArrayList<String>();

    // bolt expects the first two columns to be FID and IID
    // grab ukb's 'eid' column for this purpose, and ensure it's first
    int eid = ext.indexOfStr("eid", header, false, false);
    cols.add("" + eid);
    cols.add("" + eid);
    finalHeader.add("FID");
    finalHeader.add("IID");

    for (int i = 0; i < fields.length; i++) {
      String field = fields[i][0];
      String name = fields[i][1];
      String f = field;

      // add characters to ensure we're not matching columns too aggressively
      // assumes the csv format output by ukb-extract
      if (!f.contains("-")) f += "-";
      if (!f.startsWith("\"")) f = "\"" + field;

      int[] c = ext.indicesOfStr(f, header, false, false);

      // add all matching columns to the column list, and the corresponding common name to the
      // header list
      for (int j : c) {
        cols.add("" + j);
        finalHeader.add(ext.removeAndSimplifyQuotes(ext.replaceAllWith(header[j], field, name),
                                                    null));
      }
    }

    String[][] extract = HashVec.loadFileToStringMatrix(filename, false,
                                                        ArrayUtils.toIntArray(cols), ",", 1000,
                                                        true, true);

    boolean[] keep = null;
    // set case and control status for the phenotype
    if (caseCode != null) {
      keep = ArrayUtils.booleanArray(extract[0].length, true);
      // while we're in a pheno column, check for the case code
      for (int j = 2; j < extract[0].length && extract[0][j].contains(fields[0][0]); j++) {
        // we only want one phenotype column so if this is past the first, dump it
        if (j > 2) keep[j] = false;
        for (int i = 1; i < extract.length; i++) {
          if (extract[i][j] != null && extract[i][j].equals(caseCode)) extract[i][2] = "2";
          else if (extract[i][j] == null || !extract[i][j].equals("2")) extract[i][2] = "1";
        }
      }
    }

    for (int i = 0; i < extract[0].length; i++) {
      extract[0][i] = finalHeader.get(i);
    }

    Files.writeMatrix(extract, finalHeader.get(2) + "_extracted.tab", "\t", "NA", keep);
    return finalHeader.get(2) + "_extracted.tab";
  }

  public static void main(String[] args) {
    // read in args
    String out = null; // output file root
    String fam = "/panfs/roc/scratch/bb/chr21/ukb1773_l2r_chrY_v2_s488374_col6corrected.fam";
    String fieldsFile = "/panfs/roc/groups/5/pankrat2/shared/ukbb/phenotypes/height_fields.txt";
    String phenoFile = "/panfs/roc/groups/5/pankrat2/shared/ukbb/phenotypes/ukb9520.csv.gz";
    String caseCode = null; // if the phenotype is a code that must be converted to case/control,
                           // this is the code for a case

    for (String arg : args) {
      if (arg.startsWith("out=")) out = ext.parseStringArg(arg);
      else if (arg.startsWith("fam=")) fam = ext.parseStringArg(arg);
      else if (arg.startsWith("fieldFile=")) fieldsFile = ext.parseStringArg(arg);
      else if (arg.startsWith("phenoFile=")) phenoFile = ext.parseStringArg(arg);
      else if (arg.startsWith("caseCode=")) caseCode = ext.parseStringArg(arg);
      else System.err.println("Unknown argument: " + arg);
    }

    String[][] fields = HashVec.loadFileToStringMatrix(fieldsFile, false, new int[] {0, 1}, "\t", 0,
                                                       false);

    String extracted = extractPhenos(phenoFile, fields, caseCode);
    String[] header = Files.getHeaderOfFile(extracted, null);
    String pheno = header[2];
    if (out == null) {
      out = new File("").getAbsolutePath() + "/" + pheno;
    }
    String command = boltCommand(new File(pheno + "_extracted.tab").getAbsolutePath(), fam, out,
                                 header);

    // write bolt command to a file
    Qsub.qsub("bolt.qsub", command, 64000, 24, 24);
    // qsub it
    CmdLine.run("qsub -q pankratz bolt.qsub", new File("").getAbsolutePath());
  }
}
