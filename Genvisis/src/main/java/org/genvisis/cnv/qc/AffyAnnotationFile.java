package org.genvisis.cnv.qc;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.genvisis.cnv.filesys.AllelePair;
import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.GenomicPosition;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.filesys.Positions;

import htsjdk.variant.variantcontext.Allele;

public class AffyAnnotationFile {

  public static enum IMPORT_SCHEME {
    PROBESET_ID("Probe Set ID"), RS_OR_AFFY_ID("Affy SNP ID", "dbSNP RS ID");

    private IMPORT_SCHEME(String... headerFactors) {
      this.headerFactors = ArrayUtils.appendToArray(new String[] {"Chromosome", "Physical Position",
                                                                  "Allele A", "Allele B",
                                                                  "Ref Allele", "Flank"},
                                                    headerFactors);
      this.chr = 0;
      this.pos = 1;
      this.a = 2;
      this.b = 3;
      this.ref = 4;
      this.flank = 5;
      this.id1 = 6;
      this.id2 = 6 + (headerFactors.length > 1 ? 1 : 0);
    }

    public String[] headerFactors;
    public int id1;
    public int id2;
    public int chr;
    public int pos;
    public int a;
    public int b;
    public int ref;
    public int flank;
  }

  private final String annotFile;
  private final IMPORT_SCHEME scheme;

  Logger log = new Logger();

  /**
   * 
   * @param annotFile
   * @param log
   */
  public AffyAnnotationFile(String annotFile, IMPORT_SCHEME scheme, Logger log) {
    this.annotFile = annotFile;
    this.scheme = scheme;
  }

  public List<Marker> load() throws IOException {
    int[] hdrInds = null;
    List<Marker> markers = new ArrayList<>();
    String line = null;
    boolean foundHeader = false;
    String[] parts;
    String mkrName;
    int lines = 0;
    BufferedReader reader = Files.getAppropriateReader(annotFile);
    while ((line = reader.readLine()) != null) {
      lines++;
      if (line.charAt(0) == '#') continue;
      if (line.charAt(0) == '%') continue;
      if (!foundHeader) {
        foundHeader = true; // first nonComment line
        hdrInds = ext.indexFactors(scheme.headerFactors,
                                   ext.splitCommasIntelligently(line, true, log), false);

        continue;
      }
      parts = ext.splitCommasIntelligently(line, true, log);
      switch (scheme) {
        default:
        case PROBESET_ID:
          mkrName = parts[hdrInds[scheme.id1]];
          break;
        case RS_OR_AFFY_ID:
          mkrName = isMissing(parts[hdrInds[scheme.id2]]) ? parts[hdrInds[scheme.id1]]
                                                          : parts[hdrInds[scheme.id2]];
          break;
      }

      GenomicPosition gPos;
      gPos = new GenomicPosition(!isMissing(parts[hdrInds[scheme.chr]]) ? Positions.chromosomeNumber(parts[hdrInds[scheme.chr]],
                                                                                                     log)
                                                                        : (byte) 0,
                                 !isMissing(parts[hdrInds[scheme.pos]]) ? Integer.parseInt(parts[hdrInds[scheme.pos]])
                                                                        : 0);

      String a = parts[hdrInds[scheme.a]];
      String b = parts[hdrInds[scheme.b]];
      boolean aMiss = isMissing(a);
      boolean bMiss = isMissing(b);
      boolean aRef = parts[hdrInds[scheme.a]].equals(parts[hdrInds[scheme.ref]]);

      if (aMiss && bMiss) {
        // skip
        a = "N";
        b = "N";
      } else if (aMiss) {
        a = Character.toString(parts[hdrInds[scheme.flank]].charAt(parts[hdrInds[scheme.flank]].indexOf('[')
                                                                   - 1));
        b = a + b;
      } else if (bMiss) {
        throw new IllegalStateException("B allele is an unexpected indel for marker " + mkrName);
      }

      Marker m = new Marker(mkrName, gPos,
                            AllelePair.of(Allele.create(a, aRef), Allele.create(b, !aRef)));

      markers.add(m);
    }
    if (lines == 0 && annotFile.endsWith(".zip")) {
      throw new RuntimeException("Coudn't read the marker annotation file at " + annotFile
                                 + ".  The zipped annotation file contains two files: the annotation file itself and a README file.  "
                                 + "Please unzip these files, set the annotation file argument to point to the unzipped annotation file, and rerun the pipeline.");
    }
    reader.close();
    return markers;
  }

  private boolean isMissing(String value) {
    return value.equals("---") || value.equals("\"---\"") || value.equals("-")
           || value.equals("\"-\"");
  }

}
