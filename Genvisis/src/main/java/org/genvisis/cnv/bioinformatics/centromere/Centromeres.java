/**
 * 
 */
package org.genvisis.cnv.bioinformatics.centromere;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import org.genvisis.cnv.manage.Resources.GENOME_BUILD;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Logger;
import org.pankratzlab.shared.filesys.Positions;
import org.pankratzlab.shared.filesys.Segment;

/**
 * Methods for parsing centromere info
 */
public class Centromeres {

  private Centromeres() {
    //
  }

  /**
   * @param build
   * @param log
   * @return array of midpoints for each centromere
   */
  public static int[] getCentromereMidPoints(GENOME_BUILD build, Logger log) {
    List<Centromere> cents = loadCentromeres(build, log);
    int[] mids = new int[cents.size()];
    for (int i = 0; i < mids.length; i++) {
      mids[i] = cents.get(i).getMidPoint();
    }
    return mids;
  }

  // TODO, could really have a chromosomal arm class for the below method/parsing arms
  /**
   * Helper methods to get chromosome length for centromeres, which can be useful for creating
   * chromosomal arm segments
   * 
   * @param build
   * @param log
   * @return array of
   */
  public static int[] getChrLength(GENOME_BUILD build, Logger log) {
    List<Centromere> cents = loadCentromeres(build, log);
    int[] lengths = new int[cents.size()];
    for (int i = 0; i < lengths.length; i++) {
      switch (build) {
        case HG18:
          lengths[i] = Positions.CHROMOSOME_LENGTHS_B36_HG18[cents.get(i).getChr()];
          break;
        case HG19:
          lengths[i] = Positions.CHROMOSOME_LENGTHS_B37_HG19[cents.get(i).getChr()];
          break;
        default:
          throw new IllegalArgumentException("Invalid genome build");

      }
    }
    return lengths;
  }

  static class Centromere extends Segment {

    /**
     * 
     */
    private static final long serialVersionUID = 1L;
    private Segment p;
    private Segment q;
    private GENOME_BUILD build;
    private int chrLength;

    private Centromere(Segment p, Segment q, GENOME_BUILD build) {
      super(p.getChr(), p.getStart(), q.getStop());
      this.p = p;
      this.q = q;
      this.build = build;

      if (p.getChr() != q.getChr()) {
        throw new IllegalArgumentException("Invalid chr");
      }
      if (p.getStop() != q.getStart()) {
        throw new IllegalArgumentException("Invalid start stop");
      }
    }

    private int getMidPoint() {
      return p.getStop();
    }

    Segment getP() {
      return p;
    }

    Segment getQ() {
      return q;
    }

  }

  private static InputStream loadCents(GENOME_BUILD build) {
    return Centromeres.class.getResourceAsStream(build.getBuild() + ".cent.txt");
  }

  /**
   * Load a list of {@link Centromere} for a given {@link GENOME_BUILD}
   * 
   * @param build
   * @param log
   * @return
   */
  static List<Centromere> loadCentromeres(GENOME_BUILD build, Logger log) {

    List<Centromere> centromeres = new ArrayList<>();
    BufferedReader br;
    try {
      br = new BufferedReader(new InputStreamReader(loadCents(build)));
      String sCurrentLine;

      Segment p = null;
      Segment q = null;
      while ((sCurrentLine = br.readLine()) != null) {
        if (!sCurrentLine.startsWith("#")) {
          String[] info = sCurrentLine.trim().split("\t");
          if (p == null) {
            if (!info[3].startsWith("p")) {
              throw new IllegalArgumentException("Invalid p arm for " + ArrayUtils.toStr(info));
            }
            p = fromLine(info);
          } else {
            if (!info[3].startsWith("q")) {
              throw new IllegalArgumentException("Invalid q arm for " + ArrayUtils.toStr(info));
            }
            q = fromLine(info);
          }
          if (p != null && q != null) {
            centromeres.add(new Centromere(p, q, build));
            p = null;
            q = null;
          }

        }
      }
      br.close();
    } catch (IOException e) {
      log.reportException(e);
    }
    centromeres.add(new Centromere(new Segment((byte) 25, -1, -1), new Segment((byte) 25, -1, -1),
                                   build));
    centromeres.add(new Centromere(new Segment((byte) 26, -1, -1), new Segment((byte) 26, -1, -1),
                                   build));
    centromeres.add(new Centromere(new Segment((byte) 0, -1, -1), new Segment((byte) 0, -1, -1),
                                   build));
    Collections.sort(centromeres);
    return centromeres;

  }

  private static Segment fromLine(String[] info) {
    return new Segment(info[0], Integer.parseInt(info[1]), Integer.parseInt(info[2]));
  }

}
