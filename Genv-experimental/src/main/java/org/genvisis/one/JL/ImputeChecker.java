package org.genvisis.one.JL;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Hashtable;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.stats.Rscript.COLUMNS_MULTIPLOT;
import org.genvisis.stats.Rscript.PLOT_DEVICE;
import org.genvisis.stats.Rscript.RScatter;
import org.genvisis.stats.Rscript.RScatters;
import org.genvisis.stats.Rscript.Restrictions;
import org.genvisis.stats.Rscript.SCATTER_TYPE;

public class ImputeChecker {

  public static void main(String[] args) {
    String dir = "C:/data/misc/loganImmunoChip/";
    String immunoChip = dir + "5-0130 9K Immune-Inflammation Panel R1_5.txt";
    String aricImpute = dir + "all_aric_1000G_variants.info";
    run(dir, immunoChip, aricImpute);
  }

  public static void run(String dir, String immunoChip, String aricImpute) {
    Logger log = new Logger(dir + "log.log");
    int idIndex = ext.indexOfStr("External Id", Files.getHeaderOfFile(immunoChip, log));
    String[][] immunoRs = HashVec.loadFileToStringMatrix(immunoChip, false, null, false);
    Hashtable<String, Integer> index = new Hashtable<String, Integer>();
    for (int i = 0; i < immunoRs.length; i++) {
      index.put(immunoRs[i][idIndex], i);
    }
    String[] headerAric = Files.getHeaderOfFile(aricImpute, "[\\s]+", log);

    log.reportTimeInfo("Found " + immunoRs.length + " ids from " + immunoChip);
    log.reportTimeInfo("Found " + headerAric.length + " headers from " + aricImpute);

    String out = dir + ext.rootOf(immunoChip) + ".imputationSummary";

    if (!Files.exists(out)) {
      boolean[] foundMask = Array.booleanArray(immunoRs.length, false);
      int snpIndex = ext.indexOfStr("rs_id", headerAric);
      try {
        PrintWriter writer = new PrintWriter(new FileWriter(out));
        writer.println(
            Array.toStr(immunoRs[0]) + "\t" + Array.toStr(headerAric) + "\tFound\tDirectMatch");
        BufferedReader reader = Files.getAppropriateReader(aricImpute);
        int totalFound = 0;
        while (reader.ready()) {
          String[] line = reader.readLine().trim().split("[\\s]+");
          String snp = line[snpIndex];
          if (index.containsKey(snp)) {
            writer.println(Array.toStr(immunoRs[index.get(snp)]) + "\t" + Array.toStr(line)
                + "\t1\t" + (line[0] == snp ? 1 : 0));
            foundMask[index.get(snp)] = true;
            totalFound++;
            log.reportTimeInfo("Found " + snp + " (" + totalFound + "of " + immunoRs.length + ")");

          }
        }

        String[] blank = Array.stringArray(headerAric.length, "NA");
        for (int i = 0; i < foundMask.length; i++) {
          if (!foundMask[i]) {
            writer.println(Array.toStr(immunoRs[i]) + "\t" + Array.toStr(blank) + "\t0\t0");
          }
        }

        writer.close();
        reader.close();
      } catch (FileNotFoundException fnfe) {
        log.reportError("Error: file \"" + aricImpute + "\" not found in current directory");
        return;
      } catch (IOException ioe) {
        log.reportError("Error reading file \"" + aricImpute + "\"");
        return;
      }
    }

    String rootR = dir + "r";
    ArrayList<RScatter> rsArrayList = new ArrayList<RScatter>();
    String[] aricType = HashVec.loadFileToStringArray(out, true,
        new int[] {ext.indexOfStr("type", Files.getHeaderOfFile(out, log))}, false);
    String[] types = new String[] {"exp_freq_a1", "info", "certainty"};

    Restrictions[] restrictionTN = new Restrictions[] {
        new Restrictions(new String[] {"type"}, new double[] {0}, new String[] {"=="}, null)};
    String title = " DirectMatch=" + Array.countIf(aricType, "2") + " Impute="
        + Array.countIf(aricType, "0") + " No match = " + Array.countIf(aricType, "NA");

    for (String type : types) {
      String outtype = rootR + type;
      RScatter rsScatter = new RScatter(out, outtype + "rscript", ext.rootOf(outtype),
          outtype + ".jpg", "type", new String[] {type}, SCATTER_TYPE.HIST, log);
      rsScatter.setRestrictions(restrictionTN);
      rsScatter.setTitle(title);
      rsScatter.setyLabel("Number of imputed Markers");
      rsScatter.setxLabel(type + " of imputted markers");
      rsScatter.setOverWriteExisting(true);
      // rsScatter.execute();
      // rsArrayList.add(rsScatter);

      String outtypeBox = rootR + type + "_box";
      RScatter rsScatterBox = new RScatter(out, outtypeBox + "rscript", ext.rootOf(outtypeBox),
          outtypeBox + ".jpg", "type", new String[] {type}, SCATTER_TYPE.BOX, log);
      rsScatterBox.setRestrictions(restrictionTN);
      rsScatterBox.setTitle(title);
      rsScatterBox.setyLabel(type);
      rsScatter.setxLabel(type + " of imputed markers");

      // rsScatterBox.setxLabel(types[i]);
      rsScatterBox.setOverWriteExisting(true);
      // rsScatterBox.execute();
      // rsArrayList.add(rsScatterBox);

    }

    String outtypeFreqCer = rootR + "certfreq";
    RScatter rsScatterBoxFreq = new RScatter(out, outtypeFreqCer + "rscript",
        ext.rootOf(outtypeFreqCer), outtypeFreqCer + ".jpg", "exp_freq_a1",
        new String[] {"certainty"}, SCATTER_TYPE.POINT, log);
    rsScatterBoxFreq.setRestrictions(restrictionTN);
    rsScatterBoxFreq.setTitle(title);
    rsScatterBoxFreq.setxLabel("exp_freq_a1");
    rsScatterBoxFreq.setyLabel("certainty");
    rsScatterBoxFreq.setOverWriteExisting(true);
    rsScatterBoxFreq.execute();
    rsScatterBoxFreq.setyRange(new double[] {0, 1});
    rsArrayList.add(rsScatterBoxFreq);

    String outtypeFreqOmfp = rootR + "infofreq";
    RScatter rsScatterBoxINF = new RScatter(out, outtypeFreqOmfp + "rscript",
        ext.rootOf(outtypeFreqOmfp), outtypeFreqOmfp + ".jpg", "exp_freq_a1", new String[] {"info"},
        SCATTER_TYPE.POINT, log);
    rsScatterBoxINF.setRestrictions(restrictionTN);
    rsScatterBoxINF.setTitle(title);
    rsScatterBoxINF.setxLabel("exp_freq_a1");
    rsScatterBoxINF.setyLabel("info");
    rsScatterBoxINF.setOverWriteExisting(true);
    rsScatterBoxINF.execute();
    rsScatterBoxINF.setyRange(new double[] {0, 1});

    rsArrayList.add(rsScatterBoxINF);

    String outtypeFreqOmfpce = rootR + "freqCertinfo";
    RScatter rsScatterBoxINFd = new RScatter(out, outtypeFreqOmfpce + "rscript",
        ext.rootOf(outtypeFreqOmfpce), outtypeFreqOmfpce + ".jpg", "info",
        new String[] {"certainty"}, "exp_freq_a1", SCATTER_TYPE.POINT, log);
    rsScatterBoxINFd.setRestrictions(restrictionTN);
    rsScatterBoxINFd.setTitle(title);
    rsScatterBoxINFd.setxLabel("info");
    rsScatterBoxINFd.setyLabel("certainty");
    rsScatterBoxINFd.setOverWriteExisting(true);
    rsScatterBoxINFd.setyRange(new double[] {0, 1});
    rsScatterBoxINFd.setxRange(new double[] {0, 1});

    rsScatterBoxINFd.execute();
    rsArrayList.add(rsScatterBoxINFd);

    String finalR = dir + "final";

    RScatters rsScatters =
        new RScatters(rsArrayList.toArray(new RScatter[rsArrayList.size()]), finalR + ".rscript",
            finalR + ".pdf", COLUMNS_MULTIPLOT.COLUMNS_MULTIPLOT_1, PLOT_DEVICE.PDF, log);
    rsScatters.execute();
  }
}
