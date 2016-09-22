package org.genvisis.common;

import java.net.URLEncoder;

import org.genvisis.filesys.Segment;
import org.genvisis.filesys.SnpMarkerSet;

public class Positions {
  public static final String[] CHR_CODES = {"Un", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                                            "11", "12", "13", "14", "15", "16", "17", "18", "19",
                                            "20", "21", "22", "X", "Y", "XY", "MT"};
  public static final int[] CENTROMERE_MIDPOINTS = {-1, 124300000, 93300000, 91700000, 50700000,
                                                    47700000, 60500000, 59100000, 45200000,
                                                    51800000, 40300000, 52900000, 35400000,
                                                    16000000, 15600000, 17000000, 38200000,
                                                    22700000, 16100000, 28500000, 27100000,
                                                    12300000, 11800000, 60000000, 11500000, -1, -1};
  // public static final int[] CHROMOSOME_LENGTHS_B36_HG18 = { -1, 247200000, 242800000, 199500000,
  // 191300000, 180900000, 170900000, 158900000, 146300000, 140300000, 135400000, 134500000,
  // 132300000, 114200000, 106400000, 100400000, 88900000, 78700000, 76200000, 63900000, 62500000,
  // 47000000, 49600000, 154600000, 57800000, 154900000, 20000, -1 };
  // public static final int[] CHROMOSOME_LENGTHS_B37_HG19 = { -1, 247200000, 242800000, 199500000,
  // 191300000, 180900000, 170900000, 158900000, 146300000, 140300000, 135400000, 134500000,
  // 132300000, 114200000, 106400000, 100400000, 88900000, 78700000, 76200000, 63900000, 62500000,
  // 47000000, 49600000, 154600000, 57800000, 154900000, 20000, -1 };
  // public static final int[] CHROMOSOME_LENGTHS_B38_HG38 = { -1, 249000000, 242200000, 199500000,
  // 191300000, 180900000, 170900000, 158900000, 146300000, 140300000, 135400000, 134500000,
  // 132300000, 114200000, 106400000, 100400000, 88900000, 78700000, 76200000, 63900000, 62500000,
  // 47000000, 49600000, 154600000, 57800000, 154900000, 20000, -1 };
  public static final int[] CHROMOSOME_LENGTHS_B36_HG18 =
                                                        {-1, 247300000, 243000000, 199600000,
                                                         191300000, 180900000, 170900000, 158900000,
                                                         146300000, 140300000, 135400000, 134500000,
                                                         132400000, 114200000, 106400000, 100400000,
                                                         88900000, 78800000, 76200000, 63900000,
                                                         62500000, 47000000, 49700000, 155000000,
                                                         57800000, 155000000, 20000, -1};
  public static final int[] CHROMOSOME_LENGTHS_B37_HG19 =
                                                        {-1, 249300000, 243200000, 198100000,
                                                         191200000, 181000000, 171200000, 159200000,
                                                         146400000, 141300000, 135600000, 135100000,
                                                         133900000, 115200000, 107400000, 102600000,
                                                         90400000, 81200000, 78100000, 59200000,
                                                         63100000, 48200000, 51400000, 155300000,
                                                         59400000, 155300000, 20000, -1};
  public static final int[] CHROMOSOME_LENGTHS_B38_HG38 =
                                                        {-1, 249000000, 242200000, 198300000,
                                                         190300000, 181600000, 170900000, 159400000,
                                                         145200000, 138400000, 133800000, 135100000,
                                                         133300000, 114400000, 107100000, 102000000,
                                                         90400000, 83300000, 80400000, 58700000,
                                                         64500000, 46800000, 50900000, 156100000,
                                                         57300000, 156100000, 20000, -1};
  public static final int[] CHROMOSOME_LENGTHS_MAX = {-1, 249300000, 243200000, 199600000,
                                                      191300000, 181600000, 171200000, 159400000,
                                                      146400000, 141300000, 135600000, 135100000,
                                                      133900000, 115200000, 107400000, 102600000,
                                                      90400000, 83300000, 80400000, 63900000,
                                                      64500000, 48200000, 51400000, 156100000,
                                                      59400000, 156100000, 20000, -1};

  // public static final String[] CENTROMERE_MIDPOINT_SEGS = {"chr1:125,000,000", "chr2:93,000,000",
  // "chr3:91,000,000", "chr4:50,000,000", "chr5:48,000,000", "chr6:60,000,000", "chr7:59,500,000",
  // "chr8:45,000,000", "chr9:53,000,000", "chr10:40,000,000", "chr11:54,000,000",
  // "chr12:35,000,000", "chr13:15,600,000", "chr14:16,000,000", "chr15:16,000,000",
  // "chr16:37,500,000", "chr17:22,700,000", "chr18:16,500,000", "chr19:28,300,000",
  // "chr20:27,000,000", "chr21:11,500,000", "chr22:13,000,000", "chrX:61,000,000",
  // "chrY:11,800,000"};
  // public static final String[] CENTROMERE_MIDPOINT_SEGS_HG18 = {"chr1:124,300,000",
  // "chr2:93,300,000", "chr3:91,700,000", "chr4:50,700,000", "chr5:47,700,000", "chr6:60,500,000",
  // "chr7:59,100,000", "chr8:45,200,000", "chr9:51,800,000", "chr10:40,300,000",
  // "chr11:52,900,000", "chr12:35,400,000", "chr13:16,000,000", "chr14:15,600,000",
  // "chr15:17,000,000", "chr16:38,200,000", "chr17:22,700,000", "chr18:16,100,000",
  // "chr19:28,500,000", "chr20:27,100,000", "chr21:12,300,000", "chr22:11,800,000",
  // "chrX:59,500,000", "chrY:11,800,000"};

  public static final int[][] CENTROMERE_BOUNDARIES_FROM_SNPS_B36_HG18 = { // determined using the
                                                                           // Omni map
                                                                           {0, 0},
                                                                           {121186705, 141477078},
                                                                           {91682828, 94691119},
                                                                           {90576572, 94994003},
                                                                           {49353840, 52355341},
                                                                           {46440100, 49442477},
                                                                           {58887738, 61940508},
                                                                           {58023003, 61060840},
                                                                           {43956838, 46958273},
                                                                           {46992793, 65231255},
                                                                           {39194226, 41677980},
                                                                           {51447829, 54468566},
                                                                           {34747647, 36144166},
                                                                           {10384996, 17922259},
                                                                           {0, 18070422},
                                                                           {0, 18275409},
                                                                           {35141900, 44943958},
                                                                           {22699477, 22707804},
                                                                           {15399965, 16770020},
                                                                           {24421125, 32424385},
                                                                           {26267039, 28039018},
                                                                           {10208131, 13260157},
                                                                           {0, 14430353},
                                                                           {58588152, 61607037},
                                                                           {11799053, 12116660},
                                                                           {0, 0}, {0, 0}};

  public static final int[][] CENTROMERE_BOUNDARIES_FROM_SNPS_B37_HG19 = { // determined using the
                                                                           // Omni 2.5 map
                                                                           {0, 0},
                                                                           {121357107, 142603938},
                                                                           {92304211, 95350864},
                                                                           {90501225, 93512116},
                                                                           {49628667, 52682617},
                                                                           {46404562, 49447784},
                                                                           {58775254, 61880512},
                                                                           {58042660, 61055273},
                                                                           {43822051, 46842124},
                                                                           {46386250, 65629772},
                                                                           {39154535, 42386853},
                                                                           {51591253, 54707899},
                                                                           {34853011, 37857837},
                                                                           {0, 19058717},
                                                                           {0, 19000480},
                                                                           {0, 20005287},
                                                                           {35285609, 46416238},
                                                                           {22252439, 25268096},
                                                                           {15400035, 18511035},
                                                                           {24599528, 27735604},
                                                                           {26316921, 29423716},
                                                                           {11186430, 14338959},
                                                                           {0, 16114244},
                                                                           {58563705, 61685579},
                                                                           {24522209, 27009547},
                                                                           {0, 0}, {0, 0}};

  public static int[] parseUCSClocation(String str) {
    int chr, start, stop;

    // if (!str.startsWith("chr") || !str.contains(":") ||
    // !str.contains("-")) {
    if (!str.startsWith("chr")) {
      System.err.println("Error - '" + str + "' is not a proper UCSC position");
    }

    try {
      str = ext.replaceAllWith(str, ",", "");
      if (str.contains(":-") || str.contains("--")) {
        chr = Positions.chromosomeNumber(str.substring(3, str.indexOf(":")));
        start = -1;
        stop = -1;
      } else if (str.contains(":")) {
        chr = Positions.chromosomeNumber(str.substring(3, str.indexOf(":")));
        if (str.contains("-")) {
          try {
            start = Integer.parseInt(str.substring(str.indexOf(":") + 1, str.indexOf("-")));
          } catch (NumberFormatException e) {
            start = -1;
          }
          try {
            stop = Integer.parseInt(str.substring(str.indexOf("-") + 1));
          } catch (NumberFormatException e) {
            stop = -1;
          }
        } else {
          start = stop = Integer.parseInt(str.substring(str.indexOf(":") + 1));
        }
      } else if (str.endsWith("p")) {
        chr = Positions.chromosomeNumber(str.substring(3, str.length() - 1));
        start = 0;
        stop = CENTROMERE_MIDPOINTS[chr];
      } else if (str.endsWith("q")) {
        chr = Positions.chromosomeNumber(str.substring(3, str.length() - 1));
        start = CENTROMERE_MIDPOINTS[chr];
        stop = CHROMOSOME_LENGTHS_B36_HG18[chr];
      } else {
        chr = Positions.chromosomeNumber(str.substring(3));
        start = -1;
        stop = -1;
      }
    } catch (Exception e) {
      System.err.println("Error - '" + str + "' is not a proper UCSC position");
      return null;
    }

    return new int[] {chr, start, stop};
  }

  public static String getUCSCformat(String[] chr_pos) {
    int[] pos = new int[chr_pos.length];

    pos[0] = chromosomeNumber(chr_pos[0]);
    for (int i = 1; i < chr_pos.length; i++) {
      pos[i] = Integer.parseInt(chr_pos[i]);
    }

    return getUCSCformat(pos);
  }

  public static String getUCSCformat(int[] pos) {
    if (pos.length < 1 && pos.length > 3) {
      System.err.println("Error - could not make a valid UCSC position from '"
                         + Array.toStr(pos, "/") + "' (need 1-3 integers)");
      return null;
    }
    return "chr" + getChromosomeUCSC(pos[0], false) + (pos.length > 1
                                                                      ? ":" + pos[1]
                                                                        + (pos.length == 3 ? "-"
                                                                                             + pos[2]
                                                                                           : "")
                                                                      : "");
  }

  public static String getChromosomeUCSC(int chr, boolean addChr) {
    return getChromosomeUCSC(chr, addChr, false);
  }

  /**
   * This function is necessary due to UCSC representing the following chromosomes as follows
   * <p>
   * 23 -> "X"
   * <p>
   * 24 -> "Y"
   * <p>
   * 25 -> "X" (Pseudo-autosomal: in relation to X chromosome)
   * <p>
   * 26 - > "M"
   * <p>
   * Note that chromosome 0 will remain 0, and will not be found in the UCSC Genome database
   *
   * @param chr integer representing the chromosome to parse
   * @param addChr add "chr" to the string returned
   * @param pseudoAutosmalSetToX in some internal cases, it is useful to keep the "chrXY"
   *        designation
   * @return the original chromosome if it is less than 23. Else, the parsed version, or "" if
   *         unable to parse;
   */

  public static String getChromosomeUCSC(int chr, boolean addChr, boolean pseudoAutosmalSetToX) {
    String chrUCSC = "";
    if (chr < 23) {
      chrUCSC = chr + "";
    } else if (chr == 23) {
      chrUCSC = "X";
    } else if (chr == 24) {
      chrUCSC = "Y";
    } else if (chr == 25) {
      chrUCSC = pseudoAutosmalSetToX ? "X" : "XY";
    } else if (chr == 26) {
      chrUCSC = "M";
    }
    return addChr ? "chr" + chrUCSC : chrUCSC;
  }

  public static String getUCSClink(int[] pos) {
    return getUCSClink(pos, "hg18");
  }

  public static String getUCSClink(int[] pos, String hgBuild) {
    if (pos.length != 3) {
      System.err.println("Error - not a valid UCSC position");
      return null;
    }
    return "http://genome.ucsc.edu/cgi-bin/hgTracks?db=" + hgBuild + "&position="
           + getUCSCformat(pos);
  }

  public static String getUCSClinkInExcel(int[] pos) {
    return getUCSClinkInExcel(pos, "hg18");
  }

  public static String getUCSClinkInExcel(int[] pos, String hgBuild) {
    if (pos.length != 3) {
      System.err.println("Error - not a valid UCSC position");
      return null;
    }
    return "=HYPERLINK(\"" + getUCSClink(pos, hgBuild) + "\", \"" + getUCSCformat(pos) + "\")";
  }

  public static String getUCSCUploadLink(int[] pos, String filename) {
    // In case the encoding fails, have a default URL we can send them to
    String gatewayURL = "http://genome.ucsc.edu/cgi-bin/hgGateway";

    try {
      if (pos.length != 3) {
        System.err.println("Error - not a valid UCSC position");
        return null;
      }

      // URL encode the UCSC string, then plug it into the UCSC upload link
      String encodedUCSC = URLEncoder.encode(getUCSCformat(pos), "UTF-8");
      String encodedFilename = URLEncoder.encode(filename, "UTF-8");
      String uploadURL =
                       "http://genome.ucsc.edu/cgi-bin/hgCustom?hgHubConnect.destUrl=..%2Fcgi-bin%2FhgTracks&clade=mammal&org=Human&db=hg18&"
                         + encodedUCSC + "&hgt.customFile=" + encodedFilename
                         + "&hgt.positionInput=enter+position%2C+gene+symbol+or+search+terms&hgt.suggestTrack=knownGene&hgsid=366816683&pix=1499";
      return uploadURL;
    } catch (Exception e) {
      return gatewayURL;
    }
  }

  public static byte chromosomeNumber(String chromosome) {
    return chromosomeNumber(chromosome, new Logger());
  }

  public static byte chromosomeNumber(String chromosome, Logger log) {
    return chromosomeNumber(chromosome, true, log);
  }

  public static byte chromosomeNumber(String chromosome, boolean verbose, Logger log) {
    byte chr = -1;

    if (chromosome.startsWith("chr")) {
      chromosome = chromosome.substring(3);
    }
    if (chromosome.endsWith("p") || chromosome.endsWith("q")) {
      chromosome = chromosome.substring(0, chromosome.length() - 1);
    }

    if (chromosome.equals("XY") || chromosome.equals("PAR")) {
      chr = 25;
    } else if (chromosome.equalsIgnoreCase("X")) {
      chr = 23;
    } else if (chromosome.equalsIgnoreCase("Y")) {
      chr = 24;
    } else if (chromosome.equalsIgnoreCase("MT") || chromosome.equalsIgnoreCase("M")
               || chromosome.equalsIgnoreCase("Mito")) {
      chr = 26;
    } else if (chromosome.equalsIgnoreCase("un") || chromosome.equalsIgnoreCase("multi")
               || chromosome.equalsIgnoreCase("altonly") || chromosome.equalsIgnoreCase("noton")
               || chromosome.contains("random")) {
      chr = 0;
    } else {
      try {
        chr = Byte.parseByte(chromosome);
        if (chr < 0 || chr > 26) {
          log.reportError("Error - chromosome number '" + chromosome
                          + "' out of range for homo sapiens", true, true, 10);
        }
      } catch (NumberFormatException nfe) {
        if (verbose) {
          log.reportError("Error - '" + chromosome + "' is an invalid chromosome", true, false, 10);
        }
      }
    }

    return chr;
  }

  public static String chromosomeNumberInverse(int chr) {
    return Positions.CHR_CODES[chr];
  }

  // public static Segment[] loadCentromereMidpoints(String[] midpointPositions) {
  // Segment[] midpoints = new Segment[midpointPositions.length];
  //
  // for (int i = 0; i<midpoints.length; i++) {
  // midpoints[i] = new Segment(midpointPositions[i]);
  // }
  //
  // return midpoints;
  // }
  //
  public static Segment[] computeCentromereMidpoints(int[][] centromereBoundaries) {
    Segment[] midpoints;

    midpoints = new Segment[centromereBoundaries.length];
    for (int i = 0; i < midpoints.length; i++) {
      midpoints[i] = new Segment("chr" + i + ":"
                                 + (int) Math.floor((centromereBoundaries[i][1]
                                                     - centromereBoundaries[i][0])
                                                    / 2.0 + centromereBoundaries[i][0]));
    }

    return midpoints;
  }

  public static int[][] determineCentromereBoundariesFromMarkerSet(byte[] chrs, int[] positions,
                                                                   int build, Logger log) {
    byte chr;
    int[] midpointEstimates;
    Segment[] midpointSegments;
    int[][] centromereBoundaries;
    int markerPosition;

    if (build == 36) {
      centromereBoundaries = CENTROMERE_BOUNDARIES_FROM_SNPS_B36_HG18;
    } else if (build == 37) {
      centromereBoundaries = CENTROMERE_BOUNDARIES_FROM_SNPS_B37_HG19;
    } else {
      log.reportError("Error - build needs to be either 36 or 37 at this point in time, not '"
                      + build + "'");
      return null;
    }

    midpointEstimates = new int[27];
    midpointSegments = computeCentromereMidpoints(centromereBoundaries);
    for (chr = 0; chr < 27; chr++) {
      centromereBoundaries[chr] = new int[] {-1, Integer.MAX_VALUE};
      midpointEstimates[chr] = midpointSegments[chr].getStart();
    }

    for (int i = 0; i < positions.length; i++) {
      chr = chrs[i];
      markerPosition = positions[i];
      if (markerPosition < midpointEstimates[chr]
          && markerPosition > centromereBoundaries[chr][0]) {
        centromereBoundaries[chr][0] = markerPosition;
      }
      if (markerPosition > midpointEstimates[chr]
          && markerPosition < centromereBoundaries[chr][1]) {
        centromereBoundaries[chr][1] = markerPosition;
      }
    }

    for (chr = 0; chr < 27; chr++) {
      if (centromereBoundaries[chr][1] == Integer.MAX_VALUE) {
        centromereBoundaries[chr] = new int[] {0, 0};
      }
    }

    return centromereBoundaries;
  }

  public static int[][] determineCentromereBoundariesFromMarkerSet(String markerSetFilename,
                                                                   int build, Logger log) {
    byte chr;
    byte[] chrs;
    int[] positions;
    int[] midpointEstimates;
    Segment[] midpointSegments;
    int[][] centromereBoundaries;
    int markerPosition;
    SnpMarkerSet markerSet;

    if (build == 36) {
      centromereBoundaries = CENTROMERE_BOUNDARIES_FROM_SNPS_B36_HG18;
    } else if (build == 37) {
      centromereBoundaries = CENTROMERE_BOUNDARIES_FROM_SNPS_B37_HG19;
    } else {
      log.reportError("Error - build needs to be either 36 or 37 at this point in time, not '"
                      + build + "'");
      return null;
    }

    if (markerSetFilename == null || !Files.exists(markerSetFilename)) {
      log.reportError("Warning - since "
                      + (markerSetFilename == null ? "markerSetFilename was set to null"
                                                   : " the file '" + markerSetFilename
                                                     + "' could not be found")
                      + "; then the default centromere boundaries for build " + build
                      + " will be used");
    } else {
      markerSet = new SnpMarkerSet(markerSetFilename, false, log);
      chrs = markerSet.getChrs();
      positions = markerSet.getPositions();

      midpointEstimates = new int[27];
      midpointSegments = computeCentromereMidpoints(centromereBoundaries);
      for (chr = 0; chr < 27; chr++) {
        centromereBoundaries[chr] = new int[] {-1, Integer.MAX_VALUE};
        midpointEstimates[chr] = midpointSegments[chr].getStart();
      }

      for (int i = 0; i < positions.length; i++) {
        chr = chrs[i];
        markerPosition = positions[i];
        if (markerPosition < midpointEstimates[chr]
            && markerPosition > centromereBoundaries[chr][0]) {
          centromereBoundaries[chr][0] = markerPosition;
        }
        if (markerPosition > midpointEstimates[chr]
            && markerPosition < centromereBoundaries[chr][1]) {
          centromereBoundaries[chr][1] = markerPosition;
        }
      }

      for (chr = 0; chr < 27; chr++) {
        if (centromereBoundaries[chr][1] == Integer.MAX_VALUE) {
          centromereBoundaries[chr] = new int[] {0, 0};
        }
      }
    }

    return centromereBoundaries;
  }

  public static int[][] determineCentromereBoundariesFromMarkerSet(SnpMarkerSet markerSet,
                                                                   int build, Logger log) {
    byte chr;
    byte[] chrs;
    int[] positions;
    int[] midpointEstimates;
    Segment[] midpointSegments;
    int[][] centromereBoundaries;
    int markerPosition;

    if (build == 36) {
      centromereBoundaries = CENTROMERE_BOUNDARIES_FROM_SNPS_B36_HG18;
    } else if (build == 37) {
      centromereBoundaries = CENTROMERE_BOUNDARIES_FROM_SNPS_B37_HG19;
    } else {
      log.reportError("Error - build needs to be either 36 or 37 at this point in time, not '"
                      + build + "'");
      return null;
    }

    chrs = markerSet.getChrs();
    positions = markerSet.getPositions();

    midpointEstimates = new int[27];
    midpointSegments = computeCentromereMidpoints(centromereBoundaries);
    for (chr = 0; chr < 27; chr++) {
      centromereBoundaries[chr] = new int[] {-1, Integer.MAX_VALUE};
      midpointEstimates[chr] = midpointSegments[chr].getStart();
    }

    for (int i = 0; i < positions.length; i++) {
      chr = chrs[i];
      markerPosition = positions[i];
      if (markerPosition < midpointEstimates[chr]
          && markerPosition > centromereBoundaries[chr][0]) {
        centromereBoundaries[chr][0] = markerPosition;
      }
      if (markerPosition > midpointEstimates[chr]
          && markerPosition < centromereBoundaries[chr][1]) {
        centromereBoundaries[chr][1] = markerPosition;
      }
    }

    for (chr = 0; chr < 27; chr++) {
      if (centromereBoundaries[chr][1] == Integer.MAX_VALUE) {
        centromereBoundaries[chr] = new int[] {0, 0};
      }
    }

    return centromereBoundaries;
  }

}
