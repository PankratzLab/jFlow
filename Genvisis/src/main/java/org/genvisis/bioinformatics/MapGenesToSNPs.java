package org.genvisis.bioinformatics;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.List;
import java.util.Vector;
import org.genvisis.common.Aliases;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.Positions;
import org.genvisis.common.ext;
import org.genvisis.filesys.GeneSet;
import org.genvisis.filesys.GeneTrack;

public class MapGenesToSNPs {

  public static final int DEFAULT_BUFFER = 15000;

  public static void filter(String filename, Logger log) {
    BufferedReader reader;
    PrintWriter writer;
    String[] line, genes;
    String genesFile, lookupFile, outputFilename;
    Hashtable<String, String> locOverwrite;
    Vector<String> snps;
    int col, snpCol, chrCol, posCol, locCol;
    boolean ignoreFirstLine, commaDelimited, passes;
    GeneTrack track;
    int[][] locs;
    int buffer;
    List<String> paramV;

    buffer = -1;

    paramV = Files.parseControlFile(filename, "genes",
                                    new String[] {"candidates.txt header , 0 out=file.out",
                                                  "plink.bim 1 0 3", "buffer=" + DEFAULT_BUFFER},
                                    log);
    if (paramV == null) {
      return;
    }
    track = GeneTrack.load(Aliases.getPathToFileInReferenceDirectory(GeneSet.REFSEQ_TRACK, true,
                                                                     log));

    line = paramV.remove(0).trim().split(PSF.Regex.GREEDY_WHITESPACE);
    genesFile = line[0];
    col = -1;
    locCol = -1;
    commaDelimited = false;
    ignoreFirstLine = false;
    outputFilename = ext.rootOf(genesFile) + "_described.xln";
    for (int j = 1; j < line.length; j++) {
      if (line[j].equals("header")) {
        ignoreFirstLine = true;
      } else if (line[j].equals(",")) {
        commaDelimited = true;
      } else if (line[j].startsWith("out=")) {
        outputFilename = line[j].split("=")[1];
      } else if (col == -1) {
        col = Integer.parseInt(line[j]);
      } else if (locCol == -1) {
        locCol = Integer.parseInt(line[j]);
      } else {
        System.err.println("Error - what am I supposed to do with '" + line[j] + "'?");
      }
    }
    if (col == -1) {
      col = 0;
    }
    if (locCol == -1) {
      locCol = col + 1;
    }
    log.report("Loading genes from '" + genesFile + "'");
    genes = HashVec.loadFileToStringArray(genesFile, ignoreFirstLine, new int[] {col}, true, false,
                                          commaDelimited ? "," : PSF.Regex.GREEDY_WHITESPACE);
    locOverwrite = HashVec.loadFileToHashString(genesFile, new int[] {col}, new int[] {locCol},
                                                false, "\t", false, true);

    locs = new int[genes.length][];
    log.report("Found " + genes.length + " genes to interrogate");
    for (int i = 0; i < genes.length; i++) {
      if (locOverwrite.get(genes[i]).equals(".")) {
        locs[i] = track.lookupPosition(genes[i]);
      } else {
        locs[i] = Positions.parseUCSClocation(locOverwrite.get(genes[i]));
      }
      log.report(genes[i] + "\t" + Positions.getUCSCformat(locs[i]));
    }
    log.report("Looking up SNPs within these genes...");

    line = paramV.remove(0).trim().split(PSF.Regex.GREEDY_WHITESPACE);
    lookupFile = line[0];
    snpCol = Integer.parseInt(line[1]);
    chrCol = Positions.chromosomeNumber(line[2]);
    posCol = Integer.parseInt(line[3]);

    for (int i = 0; i < paramV.size(); i++) {
      if (paramV.get(i).startsWith("buffer=")) {
        buffer = Integer.parseInt(paramV.get(i).split("=")[1]);
        paramV.remove(i);
        i--;
      }
    }

    if (buffer == -1) {
      System.out.println("No basepair buffer defined; using default (" + DEFAULT_BUFFER + " bp)");
      buffer = DEFAULT_BUFFER;
    }

    snps = new Vector<String>(1000);
    try {
      reader = new BufferedReader(new FileReader(lookupFile));
      writer = Files.openAppropriateWriter(ext.rootOf(genesFile) + "_SNPs.xln");
      writer.println(reader.readLine());
      while (reader.ready()) {
        line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
        passes = false;
        for (int i = 0; i < locs.length && !passes; i++) {
          if (Positions.chromosomeNumber(line[chrCol]) == locs[i][0]
              && Integer.parseInt(line[posCol]) >= locs[i][1] - buffer
              && Integer.parseInt(line[posCol]) <= locs[i][2] + buffer) {
            passes = true;
          }
        }
        if (passes) {
          snps.add(line[snpCol]);
        }
      }
      reader.close();
      writer.close();
    } catch (FileNotFoundException fnfe) {
      log.reportError("Error: file \"" + lookupFile + "\" not found in current directory");
      log.reportException(fnfe);
      return;
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + lookupFile + "\"");
      log.reportException(ioe);
      return;
    }
    log.report("Found " + snps.size() + " SNPs in these genes");

    Files.combine(ArrayUtils.toStringArray(snps), ArrayUtils.toStringArray(paramV), "SNP",
                  outputFilename, log, true);
  }
}
