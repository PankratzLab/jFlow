package org.genvisis.cnv.manage;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.List;
import java.util.concurrent.Callable;
import org.genvisis.cnv.filesys.ClusterFilterCollection;
import org.genvisis.cnv.filesys.MarkerSetInfo;
import org.genvisis.cnv.filesys.Pedigree;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.ExtProjectDataParser.ProjectDataParserBuilder;
import org.genvisis.cnv.qc.MarkerMetrics;
import org.genvisis.seq.manage.VCFOps;
import org.genvisis.seq.manage.VCOps;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.CmdLine;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.WorkerTrain;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.WorkerTrain.AbstractProducer;
import org.pankratzlab.gwas.RelationAncestryQc;
import org.pankratzlab.phenoprep.PhenoPrep;
import org.pankratzlab.shared.stats.StatsCrossTabs;
import org.pankratzlab.shared.stats.Rscript.COLUMNS_MULTIPLOT;
import org.pankratzlab.shared.stats.Rscript.PLOT_DEVICE;
import org.pankratzlab.shared.stats.Rscript.RScatter;
import org.pankratzlab.shared.stats.Rscript.RScatters;
import org.pankratzlab.shared.stats.Rscript.Restrictions;
import org.pankratzlab.shared.stats.Rscript.SCATTER_TYPE;
import org.pankratzlab.shared.stats.StatsCrossTabs.STAT_TYPE;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

/**
 * Run a gwas of original (Array) genotypes for mitochondrial cn for estimates from particular
 * pc-estimates/pedigree/covars <br>
 * Currently hard-codes the recommended analysis (PC15,PC150 for natural order PCs and PC15, PCLast
 * for Stepwise)<br>
 * Requires Plink2
 */
public class MitoGWAS {

  private static final String SUB_DIR = "_eval/typed/";

  /**
   * @param proj
   * @param pedFile ped file with dna, full path
   * @param pcFile pc file ( must have valid eval data structure), full path
   * @param covarFile covariates to use, full path
   * @param outputDir where the analysis will happen, full path
   */
  public static void analyze(Project proj, String pedFile, String pcFile, String covarFile,
                             String outputDir, boolean qc, boolean perm, int numthreads) {
    Pedigree ped = new Pedigree(proj, pedFile);

    String fullOut = outputDir + ext.rootOf(proj.getPropertyFilename()) + "/";
    new File(fullOut).mkdirs();

    String root = fullOut + ext.rootOf(proj.getPropertyFilename());
    String newPed = root + ".pedigree.dat";

    ped.writeToFile(newPed, true);
    proj.PEDIGREE_FILENAME.setValue(newPed);

    String plinkPed = root + ".ped";
    String plinkMap = root + ".map";
    if (!Files.exists(plinkMap) || !Files.exists(plinkPed)) {
      MarkerSetInfo markerSet = proj.getMarkerSet();
      String[] markerNames = markerSet.getMarkerNames();
      byte[] chr = markerSet.getChrs();
      ArrayList<String> markersToAnalyze = new ArrayList<>();
      int numSkipped = 0;
      for (int i = 0; i < markerNames.length; i++) {
        if (proj.getArrayType().isCNOnly(markerNames[i]) || chr[i] < 1) {
          numSkipped++;
        } else {
          markersToAnalyze.add(markerNames[i]);
        }
      }
      proj.getLog().reportTimeInfo(numSkipped + " copy number or chr0 only probes  were removed, "
                                   + markersToAnalyze.size() + " remaining");
      String exportList = root + "_markers.txt";
      Files.writeArray(ArrayUtils.toStringArray(markersToAnalyze), exportList);
      String blankCluster = root + ".blankCluster.ser";
      proj.getLog().reportTimeWarning("Using blank cluster filter file to ensure AB lookup");
      new ClusterFilterCollection().serialize(blankCluster);
      proj.GC_THRESHOLD.setValue(0.01);
      proj.getLog().reportTimeWarning("Setting gc threshold to 0.01");
      PlinkData.saveGenvisisToPlinkPedSet(proj, root, "", blankCluster, exportList, numthreads);
    } else {
      proj.getLog().reportTimeInfo(plinkMap + " and " + plinkPed + "exist, skipping");
    }
    proj.MARKER_METRICS_FILENAME.setValue(root + "markerMetricsFull.txt");
    if (!Files.exists(proj.MARKER_METRICS_FILENAME.getValue())) {
      MarkerMetrics.fullQC(proj, null, null, false, numthreads);
    }

    double maf = (double) 5 / ped.getDnas().length;
    double geno = 0.1;
    double mind = 0.1;
    proj.getLog()
        .reportTimeInfo("using maf of " + maf + " (MAF= " + maf + "  of " + ped.getDnas().length);

    ArrayList<String> plinkConverCommand = new ArrayList<>();
    plinkConverCommand.add("plink2");
    plinkConverCommand.add("--file");
    plinkConverCommand.add(root);
    plinkConverCommand.add("--make-bed");
    plinkConverCommand.add("--out");
    root = root + "_maf_" + ext.roundToSignificantFigures(maf, 5) + "_geno_" + geno + "_mind_"
           + mind;
    plinkConverCommand.add(root);
    plinkConverCommand.add("--maf");
    plinkConverCommand.add(maf + "");
    plinkConverCommand.add("--geno");
    plinkConverCommand.add(geno + "");
    plinkConverCommand.add("--mind");
    plinkConverCommand.add("" + mind);

    String[] in = new String[] {plinkPed, plinkMap};
    String fam = root + ".fam";
    String bim = root + ".bim";
    String bed = root + ".bed";

    Files.copyFileUsingFileChannels(fam, fullOut + "plink.fam", proj.getLog());
    if (!Files.exists(fullOut + "plink.bed")) {
      Files.copyFileUsingFileChannels(bim, fullOut + "plink.bim", proj.getLog());
      Files.copyFileUsingFileChannels(bed, fullOut + "plink.bed", proj.getLog());
    }
    if (qc) {
      RelationAncestryQc.fullGamut(fullOut, "plink", true, proj.getLog());
    }

    String[] out = new String[] {fam, bim, bed};
    CmdLine.runCommandWithFileChecks(ArrayUtils.toStringArray(plinkConverCommand), "", in, out,
                                     true, false, false, proj.getLog());
    ArrayList<PlinkAssoc> plinkCommands = new ArrayList<>();
    String subDir = ext.rootOf(pcFile, false) + SUB_DIR;
    String natty = subDir + "WITHOUT_BUILDERS_NATURAL_WITHOUT_INDEPS_finalSummary.estimates.txt.gz";
    String[] nattyTitles = new String[] {"PC15", "PC150"};
    plinkCommands.addAll(generateCommands(proj, root, natty, nattyTitles, covarFile, ped, perm));

    String stepWise = subDir
                      + "WITHOUT_BUILDERS_STEPWISE_RANK_R2_WITHOUT_INDEPS_finalSummary.estimates.txt.gz";
    String[] header = Files.getHeaderOfFile(stepWise, proj.getLog());
    String[] stepWiseTitles = new String[] {"PC15", header[header.length - 1]};

    plinkCommands.addAll(generateCommands(proj, root, stepWise, stepWiseTitles, covarFile, ped,
                                          perm));

    PlinkAssocProducer producer = new PlinkAssocProducer(plinkCommands, proj.getLog());
    try (WorkerTrain<Boolean> train = new WorkerTrain<>(producer, numthreads, 2, proj.getLog())) {
      while (train.hasNext()) {
        train.next();
      }
    }

    summarize(proj, root, plinkCommands);
  }

  private static void summarize(Project proj, String root, List<PlinkAssoc> plinkCommands) {
    // MarkerMetrics.fullQC(proj, null, null, false, 24);
    String[] qcpass = HashVec.loadFileToStringArray(ext.parseDirectoryOfFile(root)
                                                    + "ldPruning/plink.prune.out", false,
                                                    new int[] {0}, true);
    qcpass = ArrayUtils.concatAll(qcpass,
                                  HashVec.loadFileToStringArray(ext.parseDirectoryOfFile(root)
                                                                + "ldPruning/plink.prune.in", false,
                                                                new int[] {0}, true));
    Hashtable<String, String> qcHash = new Hashtable<>();
    for (String qcpas : qcpass) {
      qcHash.put(qcpas, qcpas);
    }
    proj.getLog().reportTimeInfo("Computing correlation matrix of results");
    double[][] emp1s = new double[plinkCommands.size()][];
    String[] empTitles = new String[plinkCommands.size()];
    String outQQ = root + "summaryPvals";

    String pvalDB = outQQ + ".txt";
    String pvalQQ = outQQ + ".qq.txt";
    // /String plot = outQQ + ".jpeg";
    String outTabs = root + "stabs.correl.txt";
    for (int j = 0; j < plinkCommands.size(); j++) {
      String results = plinkCommands.get(j).getOutputs()[0];
      empTitles[j] = ext.removeDirectoryInfo(results);
    }
    empTitles = ArrayUtils.untag(empTitles, true, true);
    String[] pvalFiles = new String[plinkCommands.size()];
    for (int i = 0; i < pvalFiles.length; i++) {
      pvalFiles[i] = ext.parseDirectoryOfFile(pvalDB) + empTitles[i] + ".pvalsQQ.txt";
    }
    MarkerSetInfo markerSet = proj.getMarkerSet();

    String[] empLogP = ArrayUtils.tagOn(empTitles, "p_", null);
    // if (!Files.exists(pvalDB) || !Files.exists(pvalQQ)) {
    for (int j = 0; j < plinkCommands.size(); j++) {
      if (!Files.exists(pvalFiles[j])) {
        String results = plinkCommands.get(j).getOutputs()[0];
        proj.getLog().reportTimeInfo("Loading " + results);
        ProjectDataParserBuilder builderPermResults = new ProjectDataParserBuilder();
        builderPermResults.treatAllNumeric(false);
        builderPermResults.sampleBased(false);
        builderPermResults.hasHeader(true);
        builderPermResults.dataKeyColumnName("SNP");
        builderPermResults.requireAll(false);
        builderPermResults.separator(PSF.Regex.GREEDY_WHITESPACE);
        builderPermResults.setInvalidNumericToNaN(true);
        builderPermResults.firstEntryOnly(true);
        builderPermResults.numericDataTitles(new String[] {"P"});
        ExtProjectDataParser permParser;
        try {
          permParser = builderPermResults.build(proj, results);
          permParser.determineIndicesFromTitles();
          permParser.loadData();
          emp1s[j] = permParser.getNumericDataForTitle("P");

          Files.writeArray(ArrayUtils.toStringArray(emp1s[j]), pvalFiles[j]);
        } catch (FileNotFoundException e) {
          e.printStackTrace();
          return;
        }
      } else {
        emp1s[j] = ArrayUtils.toDoubleArray(HashVec.loadFileToStringArray(pvalFiles[j], false,
                                                                          new int[] {0}, false));
      }
      if (emp1s[j].length != markerSet.getMarkerNames().length) {
        throw new IllegalArgumentException("Invalid loading of pvalues");
      }

    }
    if (!Files.exists(outTabs)) {
      StatsCrossTabs sTabs = new StatsCrossTabs(emp1s, null, null, empTitles,
                                                STAT_TYPE.SPEARMAN_CORREL, true, proj.getLog());
      sTabs.computeTable();
      sTabs.dumpTables(outTabs);
    }
    boolean[] valids = ArrayUtils.booleanArray(emp1s[0].length, true);
    VCFFileReader annoReader = null;
    if (Files.exists(proj.BLAST_ANNOTATION_FILENAME.getValue())) {
      annoReader = new VCFFileReader(new File(proj.BLAST_ANNOTATION_FILENAME.getValue()), true);
    }
    String[] annotations = annoReader == null ? new String[] {}
                                              : VCFOps.getAnnotationKeys(proj.BLAST_ANNOTATION_FILENAME.getValue(),
                                                                         proj.getLog())[0];
    ExtProjectDataParser qcParser = MarkerMetrics.developParser(proj,
                                                                proj.MARKER_METRICS_FILENAME.getValue());
    int[] qcindices = qcParser.getTypedFileParser().getNumericColumns()[0];
    CloseableIterator<VariantContext> vcIter = annoReader == null ? null : annoReader.iterator();
    try {

      PrintWriter writer = Files.openAppropriateWriter(pvalDB);
      PrintWriter writerSig = Files.openAppropriateWriter(pvalDB + ".sig");
      PrintWriter writerSigQc = Files.openAppropriateWriter(pvalDB + ".sig.qc");

      String header = "SNP\tCHR\tBP\tREF\tALT\tPlinkPruneOut\t" + ArrayUtils.toStr(empTitles)
                      + "\tMinPval\tMaxPval\tMaxDiffPval\tStdPval\t"
                      + ArrayUtils.toStr(MarkerMetrics.FULL_QC_BASE_HEADER) + "\t"
                      + ArrayUtils.toStr(annotations);
      writer.println(header);
      writerSig.println(header);
      writerSigQc.println(header);
      int numInvalid = 0;
      for (int i = 0; i < emp1s[0].length; i++) {
        if (i % 10000 == 0) {
          proj.getLog().reportTimeInfo("Summarized " + i + " markers");
        }
        if (vcIter != null && !vcIter.hasNext()) {
          writer.close();
          writerSig.close();
          writerSigQc.close();
          annoReader.close();
          vcIter.close();
          throw new IllegalArgumentException("Invalid marker annotation "
                                             + proj.BLAST_ANNOTATION_FILENAME.getValue());
        }
        VariantContext markVC = vcIter == null ? null : vcIter.next();
        if (vcIter != null && !markVC.getID().equals(markerSet.getMarkerNames()[i])) {
          writer.close();
          writerSig.close();
          writerSigQc.close();
          annoReader.close();
          vcIter.close();
          throw new IllegalArgumentException("Mismatched file/project order for "
                                             + proj.BLAST_ANNOTATION_FILENAME.getValue());
        }

        boolean valid = true;
        ArrayList<String> line = new ArrayList<>(100);
        double[] pvals = new double[emp1s.length];

        for (int j = 0; j < emp1s.length; j++) {
          if (valid) {
            valid = Double.isFinite(emp1s[j][i]) && markerSet.getChrs()[i] > 0;
            if (!valid) {
              numInvalid++;
              valids[i] = false;
            }
            if (j == 0) {
              line.add(markerSet.getMarkerNames()[i]);
              line.add(markerSet.getChrs()[i] + "");
              line.add(markerSet.getPositions()[i] + "");
              line.add(markVC == null ? "NA" : markVC.getReference().getDisplayString() + "");
              line.add(markVC == null ? "NA" : markVC.getAlternateAlleles().toString() + "");
              if (qcHash.containsKey(markerSet.getMarkerNames()[i])) {
                line.add(0 + "");
              } else {
                line.add(1 + "");

              }
            }
            line.add(emp1s[j][i] + "");
            pvals[j] = emp1s[j][i];
          }
        }
        pvals = ArrayUtils.removeNaN(pvals);
        double max = ArrayUtils.max(pvals);
        double min = ArrayUtils.min(pvals);
        double maxDiff = max - min;
        double std = ArrayUtils.stdev(pvals);
        line.add(min + "");
        line.add(max + "");
        line.add(maxDiff + "");
        line.add(std + "");
        for (int qcindice : qcindices) {
          line.add(qcParser.getNumericData()[qcindice][i] + "");
        }
        if (valid) {
          String data = ArrayUtils.toStr(ArrayUtils.toStringArray(line))
                        + (markVC == null ? ""
                                          : "\t"
                                            + ArrayUtils.toStr(VCOps.getAnnotationsFor(annotations,
                                                                                       markVC,
                                                                                       ".")));
          writer.println(data);
          if (min < 0.05) {
            writerSig.println(data);
            if (qcHash.containsKey(markerSet.getMarkerNames()[i])) {
              writerSigQc.println(data);
            }
          }
        }
      }
      if (numInvalid > 0) {
        proj.getLog()
            .reportTimeWarning(numInvalid
                               + " markers had an NaN or chr 0 in one of the gwas runs ( or was not included in the analysis), removed");
      }
      writer.close();
      writerSig.close();
      writerSigQc.close();
    } catch (Exception e) {
      proj.getLog().reportError("Error writing to " + pvalDB);
      proj.getLog().reportException(e);
    }
    if (annoReader != null) {
      annoReader.close();
    }

    try {

      for (int i = 0; i < emp1s.length; i++) {
        emp1s[i] = ArrayUtils.subArray(emp1s[i], valids);
        Arrays.sort(emp1s[i]);
      }
      PrintWriter writer = Files.openAppropriateWriter(pvalQQ);
      writer.println("RANK\t" + ArrayUtils.toStr(empLogP));
      for (int i = 0; i < emp1s[0].length; i++) {
        StringBuilder builder = new StringBuilder();
        double ranke = i + 1;
        double rankP = -1 * Math.log10(ranke / emp1s[0].length);
        builder.append(rankP);
        for (int j = 0; j < emp1s.length; j++) {
          double plog = -1 * Math.log10(emp1s[j][i]);
          builder.append("\t" + plog);
        }
        writer.println(builder.toString());
      }
      writer.close();
    } catch (Exception e) {
      proj.getLog().reportError("Error writing to " + pvalQQ);
      proj.getLog().reportException(e);
    }
    // }
    proj.QQ_FILENAMES.setValue(pvalFiles);
    proj.saveProperties();
    ArrayList<RScatter> rs = new ArrayList<>();

    String out = ext.addToRoot(pvalDB, ".qc");

    String[][] qcGroups = new String[][] {{"meanTheta_AA", "meanTheta_AB", "meanTheta_BB"},
                                          {"diffTheta_AB-AA", "diffTheta_BB-AB"},
                                          {"sdTheta_AA", "sdTheta_AB", "sdTheta_BB"},
                                          {"meanR_AA", "meanR_AB", "meanR_BB"}, {"LRR_SD"},
                                          {"MARKER_GC_CONTENT"}, {"LRR_SEX_z"}};

    boolean overwrite = true;
    String gwasQCFlag = "Markers passing gwas qc";
    for (int j = 0; j < qcGroups.length; j++) {
      String currentOutMIN = ext.addToRoot(out, "METRICS_MIN" + j);
      RScatter rScatterMIN = new RScatter(pvalDB, currentOutMIN + ".rscript",
                                          ext.removeDirectoryInfo(currentOutMIN),
                                          currentOutMIN + ".jpeg", "MinPval", qcGroups[j], null,
                                          SCATTER_TYPE.POINT, proj.getLog());
      rScatterMIN.setOverWriteExisting(overwrite);
      rScatterMIN.setrScriptLoc("/panfs/roc/itascasoft/R/3.2.1/bin/Rscript");
      rScatterMIN.setxLabel("Min pval");
      rScatterMIN.execute();
      rs.add(rScatterMIN);

      String currentOutPruneMIN = ext.addToRoot(out, "GWAS_QC_METRICS_GROUP_MIN" + j);
      RScatter rScatterPruneMIN = new RScatter(pvalDB, currentOutPruneMIN + ".rscript",
                                               ext.removeDirectoryInfo(currentOutPruneMIN),
                                               currentOutPruneMIN + ".jpeg", "MinPval", qcGroups[j],
                                               null, SCATTER_TYPE.POINT, proj.getLog());

      rScatterPruneMIN.setOverWriteExisting(overwrite);
      Restrictions restrictionsMIN = new Restrictions(new String[] {"PlinkPruneOut"},
                                                      new double[] {0}, new String[] {"!="}, null);
      rScatterPruneMIN.setRestrictions(new Restrictions[] {restrictionsMIN});
      rScatterPruneMIN.setxLabel(gwasQCFlag + "Min pval");
      rScatterPruneMIN.setrScriptLoc("/panfs/roc/itascasoft/R/3.2.1/bin/Rscript");
      rScatterPruneMIN.execute();
      rs.add(rScatterPruneMIN);

      String currentOutMAX = ext.addToRoot(out, "METRICS_MAX" + j);
      RScatter rScatterMAX = new RScatter(pvalDB, currentOutMAX + ".rscript",
                                          ext.removeDirectoryInfo(currentOutMAX),
                                          currentOutMAX + ".jpeg", "MaxPval", qcGroups[j], null,
                                          SCATTER_TYPE.POINT, proj.getLog());
      rScatterMAX.setOverWriteExisting(overwrite);
      rScatterMAX.setrScriptLoc("/panfs/roc/itascasoft/R/3.2.1/bin/Rscript");
      rScatterMAX.setxLabel("Max pval");
      rScatterMAX.execute();
      rs.add(rScatterMAX);

      String currentOutPruneMAX = ext.addToRoot(out, "GWAS_QC_METRICS_GROUP_MAX" + j);
      RScatter rScatterPruneMAX = new RScatter(pvalDB, currentOutPruneMAX + ".rscript",
                                               ext.removeDirectoryInfo(currentOutPruneMAX),
                                               currentOutPruneMAX + ".jpeg", "MaxPval", qcGroups[j],
                                               null, SCATTER_TYPE.POINT, proj.getLog());

      rScatterPruneMAX.setOverWriteExisting(overwrite);
      Restrictions restrictionsMAX = new Restrictions(new String[] {"PlinkPruneOut"},
                                                      new double[] {0}, new String[] {"!="}, null);
      rScatterPruneMAX.setRestrictions(new Restrictions[] {restrictionsMAX});
      rScatterPruneMAX.setxLabel(gwasQCFlag + "Max pvals");
      rScatterPruneMAX.setrScriptLoc("/panfs/roc/itascasoft/R/3.2.1/bin/Rscript");
      rScatterPruneMAX.execute();
      rs.add(rScatterPruneMAX);

      String currentOutSD = ext.addToRoot(out, "METRICS_SD" + j);
      RScatter rScatterSD = new RScatter(pvalDB, currentOutSD + ".rscript",
                                         ext.removeDirectoryInfo(currentOutSD),
                                         currentOutSD + ".jpeg", "StdPval", qcGroups[j], null,
                                         SCATTER_TYPE.POINT, proj.getLog());
      rScatterSD.setOverWriteExisting(overwrite);
      rScatterSD.setrScriptLoc("/panfs/roc/itascasoft/R/3.2.1/bin/Rscript");
      rScatterSD.setxLabel("SD pvals");
      rScatterSD.execute();
      rs.add(rScatterSD);

      String currentOutPruneSD = ext.addToRoot(out, "GWAS_QC_METRICS_GROUP_SD" + j);
      RScatter rScatterPruneSD = new RScatter(pvalDB, currentOutPruneSD + ".rscript",
                                              ext.removeDirectoryInfo(currentOutPruneSD),
                                              currentOutPruneSD + ".jpeg", "StdPval", qcGroups[j],
                                              null, SCATTER_TYPE.POINT, proj.getLog());

      rScatterPruneSD.setOverWriteExisting(overwrite);
      Restrictions restrictionsSD = new Restrictions(new String[] {"PlinkPruneOut"},
                                                     new double[] {0}, new String[] {"!="}, null);
      rScatterPruneSD.setRestrictions(new Restrictions[] {restrictionsSD});
      rScatterPruneSD.setxLabel(gwasQCFlag + "SD pvals");
      rScatterPruneSD.setrScriptLoc("/panfs/roc/itascasoft/R/3.2.1/bin/Rscript");
      rScatterPruneSD.execute();
      rs.add(rScatterPruneSD);

      String currentOut = ext.addToRoot(out, "METRICS_" + j);
      RScatter rScatter = new RScatter(pvalDB, currentOut + ".rscript",
                                       ext.removeDirectoryInfo(currentOut), currentOut + ".jpeg",
                                       "MaxDiffPval", qcGroups[j], null, SCATTER_TYPE.POINT,
                                       proj.getLog());
      rScatter.setOverWriteExisting(overwrite);
      rScatter.setrScriptLoc("/panfs/roc/itascasoft/R/3.2.1/bin/Rscript");
      rScatter.setxLabel("Maximum Pval difference");
      rScatter.execute();
      rs.add(rScatter);

      String currentOutPrune = ext.addToRoot(out, "GWAS_QC_METRICS_GROUP_" + j);
      RScatter rScatterPrune = new RScatter(pvalDB, currentOutPrune + ".rscript",
                                            ext.removeDirectoryInfo(currentOutPrune),
                                            currentOutPrune + ".jpeg", "MaxDiffPval", qcGroups[j],
                                            null, SCATTER_TYPE.POINT, proj.getLog());
      rScatterPrune.setOverWriteExisting(overwrite);
      Restrictions restrictions = new Restrictions(new String[] {"PlinkPruneOut"}, new double[] {0},
                                                   new String[] {"!="}, null);
      rScatterPrune.setRestrictions(new Restrictions[] {restrictions});
      rScatterPrune.setxLabel(gwasQCFlag + "Maximum Pval difference");
      rScatterPrune.setrScriptLoc("/panfs/roc/itascasoft/R/3.2.1/bin/Rscript");
      rScatterPrune.execute();
      rs.add(rScatterPrune);

    }

    RScatters rsScatters = new RScatters(rs.toArray(new RScatter[rs.size()]), out + "finalRscript",
                                         out + "final.pdf", COLUMNS_MULTIPLOT.COLUMNS_MULTIPLOT_1,
                                         PLOT_DEVICE.PDF, proj.getLog());
    rsScatters.execute();
    // String[] script = generateManhatQQScript(pvalDB, empTitles, plot);
    // String rscript = plot + ".rscript";
    // Files.writeList(script, rscript);
    // CmdLine.runCommandWithFileChecks(new String[] { "/panfs/roc/itascasoft/R/3.1.1/bin/Rscript",
    // rscript }, "", null, new String[] { plot }, true, true, false, proj.getLog());
  }

  // private static String[] generateMadnhatQQScript(String db, String[] pvalColumns, String output)
  // {
  // ArrayList<String> command = new ArrayList<String>();
  // ArrayList<String> order = new ArrayList<String>();
  //
  // command.add("library(qqman)");
  //
  // String main = "data";
  // command.add(main + "=read.table(\"" + db + "\", header=TRUE)");
  // //
  // for (int i = 0; i < pvalColumns.length; i++) {
  // command.add("jpeg(file=\"" + ext.parseDirectoryOfFile(output) + pvalColumns[i] + ".jpeg" +
  // "\",height=2000,width=2500)");
  // // ,onefile = TRUE
  // command.add("op <- par(mfrow=c(1,2))");
  //
  // order.add(pvalColumns[i]);
  // String man = pvalColumns[i] + "man =manhattan(" + main + ", p =\"" + pvalColumns[i] + "\")";
  // command.add(man);
  // command.add(pvalColumns[i] + "q =qq(" + main + "$" + pvalColumns[i] + ")");
  // command.add(pvalColumns[i] + "man");
  // command.add("title(main = \"" + pvalColumns[i] + "\")");
  // command.add(pvalColumns[i] + "q");
  // command.add("title(main = \"" + pvalColumns[i] + "\")");
  // command.add("par(op)");
  // command.add("dev.off()");
  // }
  //
  // Files.writeList(Array.toStringArray(order), output + ".order.txt");
  // return Array.toStringArray(command);
  // }

  private static ArrayList<PlinkAssoc> generateCommands(Project proj, String root,
                                                        String mtPhenoFile, String[] titles,
                                                        String covarFile, Pedigree ped,
                                                        boolean perm) {
    ArrayList<PlinkAssoc> plinkCommands = new ArrayList<>();

    ProjectDataParserBuilder builderCurrent = new ProjectDataParserBuilder();
    builderCurrent.numericDataTitles(titles);
    builderCurrent.sampleBased(true);
    builderCurrent.dataKeyColumnName("DNA");
    builderCurrent.treatAllNumeric(false);
    builderCurrent.requireAll(true);
    ExtProjectDataParser parser;
    try {
      parser = builderCurrent.build(proj, mtPhenoFile);
      parser.determineIndicesFromTitles();
      parser.loadData();

      ProjectDataParserBuilder builderCovar = new ProjectDataParserBuilder();
      builderCovar.sampleBased(true);
      builderCovar.hasHeader(true);
      builderCovar.dataKeyColumnName("DNA");
      builderCovar.requireAll(false);
      builderCovar.numericDataTitles(ArrayUtils.subArray(Files.getHeaderOfFile(covarFile,
                                                                               proj.getLog()),
                                                         1));
      ExtProjectDataParser covarParser = builderCovar.build(proj, covarFile);
      covarParser.determineIndicesFromTitles();
      covarParser.loadData();
      boolean[] hasVarianceWithinPed = ArrayUtils.booleanArray(covarParser.getNumericDataTitles().length,
                                                               false);
      for (int i = 0; i < hasVarianceWithinPed.length; i++) {
        Hashtable<String, String> varHash = new Hashtable<>();
        for (int j = 0; j < ped.getDnas().length; j++) {
          int sampIndex = ext.indexOfStr(ped.getDnas()[j], proj.getSamples());
          String val = covarParser.getNumericDataForTitle(covarParser.getNumericDataTitles()[i])[sampIndex]
                       + "";
          varHash.put(val, val);
          if (varHash.size() > 1) {
            hasVarianceWithinPed[i] = true;
            break;
          }
        }
        if (varHash.size() < 2) {
          proj.getLog().reportTimeWarning(covarParser.getNumericDataTitles()[i]
                                          + " had no variance in ped samples, removing");
        }
      }

      for (String title : titles) {
        String outCurrent = root + ext.rootOf(mtPhenoFile) + "_" + title + ".txt";
        for (int j = 0; j < covarParser.getNumericDataTitles().length; j++) {
          if (hasVarianceWithinPed[j]) {
            outCurrent = ext.addToRoot(outCurrent, "_" + covarParser.getNumericDataTitles()[j]);
          }
        }
        try {
          PrintWriter writer = Files.openAppropriateWriter(outCurrent);
          writer.print("FID\tIID\t" + title);
          for (int j = 0; j < covarParser.getNumericDataTitles().length; j++) {
            if (hasVarianceWithinPed[j]) {
              writer.print("\t" + covarParser.getNumericDataTitles()[j]);
            }
          }
          writer.println();
          for (int j = 0; j < ped.getDnas().length; j++) {
            int sampIndex = ext.indexOfStr(ped.getDnas()[j], proj.getSamples());
            writer.print(ped.getFID(j) + "\t" + ped.getIID(j) + "\t"
                         + parser.getNumericDataForTitle(title)[sampIndex]);
            for (int k = 0; k < covarParser.getNumericDataTitles().length; k++) {
              if (hasVarianceWithinPed[k]) {
                writer.print("\t"
                             + covarParser.getNumericDataForTitle(covarParser.getNumericDataTitles()[k])[sampIndex]);
              }
            }
            writer.println();
          }
          writer.close();

          String covarTitles = ArrayUtils.toStr(ArrayUtils.subArray(covarParser.getNumericDataTitles(),
                                                                    hasVarianceWithinPed),
                                                ",");
          String outPrepReg = ext.addToRoot(outCurrent, ".prepped");

          PlinkAssoc regCommand = prepareAssoc(root, false, outCurrent, title, covarTitles,
                                               outPrepReg, perm, proj.getLog());
          plinkCommands.add(regCommand);
          String outPrepInv = ext.addToRoot(outCurrent, ".prepped.inv");
          PlinkAssoc invCommand = prepareAssoc(root, true, outCurrent, title, covarTitles,
                                               outPrepInv, perm, proj.getLog());
          plinkCommands.add(invCommand);

        } catch (Exception e) {
          proj.getLog().reportError("Error writing to " + outCurrent);
          proj.getLog().reportException(e);
        }
      }

    } catch (FileNotFoundException e) {
      e.printStackTrace();
    }
    return plinkCommands;

  }

  private static class PlinkAssocProducer extends AbstractProducer<Boolean> {

    private final List<PlinkAssoc> assocs;
    private int index;

    public PlinkAssocProducer(List<PlinkAssoc> assocs, Logger log) {
      super();
      this.assocs = assocs;
      index = 0;

    }

    @Override
    public boolean hasNext() {
      return index < assocs.size();
    }

    @Override
    public Callable<Boolean> next() {
      PlinkAssoc toReturn = assocs.get(index);
      index++;
      return toReturn;
    }
  }

  private static class PlinkAssoc implements Callable<Boolean> {

    private final List<String> command;
    private final String[] inputs;
    private final String[] outputs;
    private final Logger log;

    public PlinkAssoc(List<String> command, String[] inputs, String[] outputs, Logger log) {
      super();
      this.command = command;
      this.inputs = inputs;
      this.outputs = outputs;
      this.log = log;
    }

    public String[] getOutputs() {
      return outputs;
    }

    @Override
    public Boolean call() throws Exception {
      return CmdLine.runCommandWithFileChecks(ArrayUtils.toStringArray(command), "", inputs,
                                              outputs, true, false, false, log);
    }
  }

  private static PlinkAssoc prepareAssoc(String root, boolean inverse, String inputDb, String pheno,
                                         String covars, String output, boolean perm, Logger log) {
    String processed = ext.addToRoot(output, "_pheno");
    if (!Files.exists(processed)) {
      PhenoPrep.parse("", inputDb, "IID", pheno, null, 3.0, false, false, true, false, inverse,
                      covars, root + ".fam", true, true, false, false, true, false, null, output,
                      true, false, false, false, false, null, false, log);
    }

    ArrayList<String> plink = new ArrayList<>();
    plink.add("plink2");

    plink.add("--linear");
    if (perm) {
      plink.add("perm");
    }

    plink.add("--bfile");
    plink.add(root);
    plink.add("--covar-name");
    plink.add(covars);
    plink.add("--covar");
    plink.add(inputDb);
    plink.add("--pheno");
    plink.add(processed);
    plink.add("--out");
    plink.add(ext.rootOf(output, false));
    plink.add("--threads");
    plink.add(3 + "");

    String[] inputs = new String[] {processed, inputDb};
    String[] outputs = new String[] {ext.rootOf(output, false) + ".assoc.linear"};
    if (perm) {
      outputs = ArrayUtils.concatAll(outputs, new String[] {ext.rootOf(output, false)
                                                            + ".assoc.linear.perm"});
    }
    PlinkAssoc assoc = new PlinkAssoc(plink, inputs, outputs, log);
    return assoc;
  }

  // public static void test() {
  // Project proj = new Project("/home/pankrat2/lanej/projects/gedi_gwas.properties", false);
  // String ped = "/panfs/roc/groups/5/pankrat2/shared/gedi_gwas/mitoAnalyze/gedi_gwas.ped";
  // String covar = "/panfs/roc/groups/5/pankrat2/shared/gedi_gwas/mitoAnalyze/gedi_gwas.covar";
  // int numthreads = 12;
  // String outputDir = "/scratch.global/lanej/mitoAnalyze/";
  // String pcFile = proj.PROJECT_DIRECTORY.getValue() +
  // "gedi_gwasALL_1000PCs_OHW_40_ws15_gc_corrected_recomp.PCs.extrapolated.txt";
  // // analyze(proj, ped, pcFile, covar, outputDir, numthreads);
  // }

  public static void main(String[] args) {

    int numArgs = args.length;
    String filename = null;
    String ped = null;
    String outputDir = null;
    int numthreads = 24;
    String pcFile = null;
    String covar = null;
    boolean qc = true;
    boolean perm = false;

    String usage = "\n" + "one.JL.MitoAnalyze requires 0-1 arguments\n";
    usage += "   (1) filename (i.e. proj= (no default))\n" + "";
    usage += "   (2) ped (i.e. ped= (no default))\n" + "";
    usage += "   (3) output directory (i.e. out= (no default))\n" + "";
    usage += "   (4) numthreads (i.e. numthreads=" + numthreads + " ( default))\n" + "";
    usage += "   (5) PC file (i.e. pc= (no default))\n" + "";
    usage += "   (6) covarFile (i.e. cov= (no default))\n" + "";
    usage += "   (7) qc  (i.e. qc=" + qc + " (no default))\n" + "";
    usage += "   (8) permute  (i.e. perm=" + perm + " (no default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      }

      else if (arg.startsWith("proj=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("ped=")) {
        ped = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("out=")) {
        outputDir = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("pc=")) {
        pcFile = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("cov=")) {
        covar = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("qc=")) {
        qc = ext.parseBooleanArg(arg);
        numArgs--;
      } else if (arg.startsWith("perm=")) {
        perm = ext.parseBooleanArg(arg);
        numArgs--;
      } else if (arg.startsWith("numthreads=")) {
        numthreads = ext.parseIntArg(arg);
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
      Project proj = new Project(filename);
      analyze(proj, ped, pcFile, covar, outputDir, qc, perm, numthreads);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
