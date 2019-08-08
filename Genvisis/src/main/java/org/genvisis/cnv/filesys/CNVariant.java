package org.genvisis.cnv.filesys;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Vector;
import java.util.stream.Collectors;

import org.genvisis.cnv.annotator.AnnotatedSegment;
import org.genvisis.cnv.hmm.PennHmm;
import org.genvisis.seq.manage.BEDFileReader;
import org.pankratzlab.common.ByteVector;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.IntVector;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.PSF;
import org.pankratzlab.common.SerializedFiles;
import org.pankratzlab.common.StringVector;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.filesys.LocusSet;
import org.pankratzlab.common.filesys.Positions;

public class CNVariant extends AnnotatedSegment {

  @Override
  public int hashCode() {
    final int prime = 31;
    int result = super.hashCode();
    result = prime * result + cn;
    result = prime * result + ((familyID == null) ? 0 : familyID.hashCode());
    result = prime * result + ((individualID == null) ? 0 : individualID.hashCode());
    result = prime * result + numMarkers;
    long temp;
    temp = Double.doubleToLongBits(score);
    result = prime * result + (int) (temp ^ (temp >>> 32));
    result = prime * result + source;
    return result;
  }

  @Override
  public boolean equals(Object obj) {
    if (this == obj) {
      return true;
    }
    if (!super.equals(obj)) {
      return false;
    }
    if (getClass() != obj.getClass()) {
      return false;
    }
    CNVariant other = (CNVariant) obj;
    if (cn != other.cn) {
      return false;
    }
    if (familyID == null) {
      if (other.familyID != null) {
        return false;
      }
    } else if (!familyID.equals(other.familyID)) {
      return false;
    }
    if (individualID == null) {
      if (other.individualID != null) {
        return false;
      }
    } else if (!individualID.equals(other.individualID)) {
      return false;
    }
    if (numMarkers != other.numMarkers) {
      return false;
    }
    if (Double.doubleToLongBits(score) != Double.doubleToLongBits(other.score)) {
      return false;
    }
    if (source != other.source) {
      return false;
    }
    return true;
  }

  public static final long serialVersionUID = 1L;
  public static final String[] PLINK_CNV_HEADER = {"FID", "IID", "CHR", "BP1", "BP2", "TYPE",
                                                   "SCORE", "SITES"};

  private String familyID;
  private String individualID;
  private int cn;
  private double score;
  protected int numMarkers;
  private int source;

  public CNVariant(String[] plinkLine) {
    this(plinkLine, 0);
  }

  public CNVariant(String[] plinkLine, int sourceIndex) {
    super(Positions.chromosomeNumber(plinkLine[2]), Integer.parseInt(plinkLine[3]),
          Integer.parseInt(plinkLine[4]));
    familyID = plinkLine[0];
    individualID = plinkLine[1];
    cn = Integer.parseInt(plinkLine[5]);
    score = Double.parseDouble(plinkLine[6]);
    numMarkers = Integer.parseInt(plinkLine[7]);
    source = sourceIndex;
  }

  public CNVariant(String familyID, String individualID, byte chr, int start, int stop, int cn,
                   double score, int numMarkers, int source) {
    super(chr, start, stop);
    this.familyID = familyID;
    this.individualID = individualID;
    this.cn = cn;
    this.score = score;
    this.numMarkers = numMarkers;
    this.source = source;
  }

  public CNVariant(CNVariant cnv) {
    super(cnv.getChr(), cnv.getStart(), cnv.getStop());
    familyID = cnv.familyID;
    individualID = cnv.individualID;
    cn = cnv.cn;
    score = cnv.score;
    numMarkers = cnv.numMarkers;
    source = cnv.source;
  }

  public CNVariant(String ucscLocation) {
    super(ucscLocation);
  }

  public void setFamilyID(String famID) {
    familyID = famID;
  }

  public String getFamilyID() {
    return familyID;
  }

  public void setIndividualID(String indID) {
    individualID = indID;
  }

  public String getIndividualID() {
    return individualID;
  }

  public int getSource() { // generally not entered, since loading from a serialized file
    return source;
  }

  public void setSource(int sourceIndex) {
    source = sourceIndex;
  }

  public int getCN() {
    return cn;
  }

  public double getScore() {
    return score;
  }

  public int getNumMarkers() {
    return numMarkers;
  }

  public boolean equalsIncludingIndividual(CNVariant cnv) {
    return cnv.familyID.equals(familyID) && cnv.individualID.equals(individualID)
           && cnv.getChr() == getChr() && cnv.getStart() == getStart() && cnv.getStop() == getStop()
           && cnv.cn == cn;
  }

  public boolean overlapsLocAndIndividual(CNVariant cnv) {
    return familyID.equals(cnv.familyID) && individualID.equals(cnv.individualID)
           && amountOfOverlapInBasepairs(cnv) > 0;
  }

  public boolean overlapsLocAndIndividualSignificantly(CNVariant cnv) {
    boolean overlapsLocAndIndividualSignificantly = familyID.equals(cnv.familyID)
                                                    && individualID.equals(cnv.individualID)
                                                    && significantOverlap(cnv);
    return overlapsLocAndIndividualSignificantly;
  }

  public String toPlinkFormat() {
    return familyID + "\t" + individualID + "\t" + getChr() + "\t" + getStart() + "\t" + getStop()
           + "\t" + cn + "\t" + ext.formDeci(score, 5) + "\t" + numMarkers;
  }

  /**
   * Assumes that {@link CNVariant#individualID} is actually the DNA id
   *
   * @return region list ready for trailer, with the score as the comment
   */
  public String[] toTrailerFormat() {
    String[] tTrail = new String[3];
    tTrail[0] = individualID;
    tTrail[1] = getUCSClocation();
    tTrail[2] = score + "";
    return tTrail;
  }

  public String getFingerprint() {
    return familyID + "_" + individualID + "_" + getChr() + "_" + getStart() + "_" + getStop() + "_"
           + cn + "_" + numMarkers;
  }

  public static CNVariant[] toCNVariantArray(List<CNVariant> v) {
    return v.toArray(new CNVariant[v.size()]);
  }

  public static CNVariant[] sortCNVsInPlace(CNVariant[] array) {
    Arrays.sort(array);

    return array;
  }

  public static CNVariant[] loadPlinkFile(String filename) {
    return CNVariant.sortCNVsInPlace(CNVariant.toCNVariantArray(loadPlinkFile(filename, null,
                                                                              true)));
  }

  public static CNVariant[] loadPlinkFile(String filename, boolean includeLOH) {
    return CNVariant.sortCNVsInPlace(CNVariant.toCNVariantArray(loadPlinkFile(filename, null,
                                                                              includeLOH)));
  }

  public static List<CNVariant> loadPlinkFile(String filename, Hashtable<String, String> sampleHash,
                                              boolean includeLOH) {
    BufferedReader reader;
    List<CNVariant> v = null;
    String[] line;

    v = new ArrayList<>();
    try {
      reader = Files.getReader(filename, true, true);

      reader.mark(1000);
      line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
      if (!line[2].toLowerCase().equals("chr") && Positions.chromosomeNumber(line[2]) != -1) {
        reader.reset();
      }
      while (reader.ready()) {
        line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
        if (sampleHash == null || sampleHash.containsKey(line[0] + "\t" + line[1])) {
          CNVariant var = new CNVariant(line);
          if (!includeLOH && var.getCN() == PennHmm.LOH_FLAG) {
            var = null;
            continue;
          }
          v.add(new CNVariant(line));
        }
      }
      reader.close();

      return v;
    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + filename + "\" not found in current directory");
      fnfe.printStackTrace();
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + filename + "\"");
      ioe.printStackTrace();
    }

    return null;
  }

  public static List<CNVariant> loadFromBED(String bedFile) {
    return loadFromBED(bedFile, "_");
  }

  public static List<CNVariant> loadFromBED(String bedFile, String filenameDelim) {
    String id = ext.removeDirectoryInfo(bedFile).split(filenameDelim)[0]; // ID is defined as the
                                                                          // first token
                                                                          // in the filename
    try (BEDFileReader reader = new BEDFileReader(bedFile, false);) {

      return reader.iterator().stream().map(bf -> {
        return new CNVariant(id, id, Positions.chromosomeNumber(bf.getContig()), bf.getStart(),
                             bf.getEnd(), 2, -1, -1, 0);
      }).collect(Collectors.toList());
    }
  }

  // TODO what was this created for? Incomplete
  public static void mergeCNVs(String filename) {
    BufferedReader reader;
    String[] line;
    // PrintWriter writer;
    // Hashtable<String, Vector<CNVariant>> hash;
    // hash = new Hashtable<String, Vector<CNVariant>>();
    StringVector markerNames;
    ByteVector chrs;
    IntVector positions;

    // hash.add(fID+"\t"+iId, Vector<CNVariant>);
    try {
      reader = new BufferedReader(new FileReader(filename));

      markerNames = new StringVector();
      chrs = new ByteVector();
      positions = new IntVector();

      reader.readLine();
      while (reader.ready()) {
        line = reader.readLine().trim().split(PSF.Regex.GREEDY_WHITESPACE);
        markerNames.add(line[0]);
        chrs.add(Positions.chromosomeNumber(line[1]));
        positions.add(Integer.parseInt(line[2]));
      }
      reader.close();
      /*
       * new MarkerSet(markerNames, chrs, positions).serialize(filename+".ser"); writer =
       * Files.openAppropriateWriter("????"); for (int i = 0; i<temp.length; i++) {
       * writer.println(temp[i].toPlinkFormat()); } writer.close();
       */

    } catch (FileNotFoundException fnfe) {
      System.err.println("Error: file \"" + filename + "\" not found in current directory");
      System.exit(1);
    } catch (IOException ioe) {
      System.err.println("Error reading file \"" + filename + "\"");
      System.exit(2);
    }
  }

  public static void findConsensus(String file1, String file2) {
    PrintWriter writer;
    CNVariant[] list1, list2, consensus;
    Vector<CNVariant> v = new Vector<>();

    list1 = loadPlinkFile(file1);
    list2 = loadPlinkFile(file2);

    for (CNVariant element : list1) {
      for (CNVariant element2 : list2) {
        if (element.overlapsLocAndIndividual(element2)) {
          v.add(element);
        }
      }
    }
    consensus = sortCNVsInPlace(CNVariant.toCNVariantArray(v));

    try {
      writer = Files.openAppropriateWriter(ext.rootOf(file1) + "_" + ext.rootOf(file2)
                                           + "_consensus.cnv");
      for (CNVariant consensu : consensus) {
        writer.println(consensu.toPlinkFormat());
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing consensus");
      e.printStackTrace();
    }
  }

  public enum CONSENSUS_TYPE {
    /**
     * Check copy number state,
     */
    CN_AWARE,
    /**
     * Do not check copy number state
     */
    NOT_CN_AWARE;
  }

  public enum OVERLAP_TYPE {
    /**
     * Must have significant overlap
     */
    OVERLAP_LOC_AND_INDIVIDUAL_SIGNIFICANTLY,
    /**
     * Any overlap
     */
    OVERLAP_LOC_AND_INDIVIDUAL;
  }

  public static MatchResults findSignificantConsensus(String file1, String file2,
                                                      CONSENSUS_TYPE cType, OVERLAP_TYPE oType) {
    String output = ext.rootOf(file1) + "_" + ext.rootOf(file2) + "_signif_consensus.cnv";
    return findSignificantConsensus(file1, file2, output, cType, oType);
  }

  public static MatchResults findSignificantConsensus(String file1, String file2, String output,
                                                      CONSENSUS_TYPE cType, OVERLAP_TYPE oType) {
    PrintWriter writer;
    CNVariant[] list1, list2;
    HashSet<CNVariant> matched1 = new HashSet<>();
    HashSet<CNVariant> matched2 = new HashSet<>();

    HashSet<CNVariant> unmatched1 = new HashSet<>();
    HashSet<CNVariant> unmatched2 = new HashSet<>();
    Vector<String> outputLines = new Vector<>();

    list1 = loadPlinkFile(file1);
    list2 = loadPlinkFile(file2);

    for (CNVariant element : list1) {
      for (CNVariant element2 : list2) {
        boolean match = false;
        switch (oType) {
          case OVERLAP_LOC_AND_INDIVIDUAL:
            match = element.overlapsLocAndIndividual(element2);
            break;
          case OVERLAP_LOC_AND_INDIVIDUAL_SIGNIFICANTLY:
            match = element.overlapsLocAndIndividualSignificantly(element2);
            break;
          default:
            break;
        }
        if (match) {
          switch (cType) {
            case CN_AWARE:
              match = element.getCN() == element2.getCN();
              break;
            case NOT_CN_AWARE:
              break;
            default:
              break;

          }
        }
        if (match) {
          matched1.add(element);
          matched2.add(element2);
        }
      }
    }

    for (CNVariant element : list1) {
      if (matched1.contains(element)) {
        outputLines.add(element.toPlinkFormat() + "\t3");
      } else {
        unmatched1.add(element);
      }
    }

    for (CNVariant element : list2) {
      if (matched2.contains(element)) {
        outputLines.add(element.toPlinkFormat() + "\t3");
      } else {
        unmatched2.add(element);
      }
    }

    Iterator<CNVariant> unmatchedIterator = unmatched1.iterator();
    while (unmatchedIterator.hasNext()) {
      outputLines.add(unmatchedIterator.next().toPlinkFormat() + "\t1");
    }
    unmatchedIterator = unmatched2.iterator();
    while (unmatchedIterator.hasNext()) {
      outputLines.add(unmatchedIterator.next().toPlinkFormat() + "\t2");
    }

    try {
      writer = Files.openAppropriateWriter(output);
      for (int i = 0; i < outputLines.size(); i++) {
        writer.println(outputLines.get(i));
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing significant-consensus file");
      e.printStackTrace();
    }
    return new MatchResults(file1, file2, matched1, matched2, unmatched1, unmatched2);

  }

  public static class MatchResults implements Serializable {

    /**
     *
     */
    private static final long serialVersionUID = 1L;
    private final String file1;
    private final String file2;
    private final HashSet<CNVariant> matched1;
    private final HashSet<CNVariant> matched2;

    private final HashSet<CNVariant> unmatched1;
    private final HashSet<CNVariant> unmatched2;

    public MatchResults(String file1, String file2, HashSet<CNVariant> matched1,
                        HashSet<CNVariant> matched2, HashSet<CNVariant> unmatched1,
                        HashSet<CNVariant> unmatched2) {
      super();
      this.file1 = file1;
      this.file2 = file2;
      this.matched1 = matched1;
      this.matched2 = matched2;
      this.unmatched1 = unmatched1;
      this.unmatched2 = unmatched2;
    }

    public String getFile1() {
      return file1;
    }

    public String getFile2() {
      return file2;
    }

    public HashSet<CNVariant> getMatched1() {
      return matched1;
    }

    public HashSet<CNVariant> getMatched2() {
      return matched2;
    }

    public HashSet<CNVariant> getUnmatched1() {
      return unmatched1;
    }

    public HashSet<CNVariant> getUnmatched2() {
      return unmatched2;
    }

    public void writeSerial(String fileName) {
      SerializedFiles.writeSerial(this, fileName, true);
    }

    public static MatchResults readSerial(String fileName, Logger log) {
      return (MatchResults) SerializedFiles.readSerial(fileName, log, false, true);
    }

  }

  public static void findSignificantConsensus(String file1, String file2, boolean checkLarger) {
    PrintWriter writer;
    CNVariant[] list1, list2;
    HashSet<CNVariant> matched1 = new HashSet<>();
    HashSet<CNVariant> matched2 = new HashSet<>();

    HashSet<CNVariant> unmatched1 = new HashSet<>();
    HashSet<CNVariant> unmatched2 = new HashSet<>();
    Vector<String> outputLines = new Vector<>();

    list1 = loadPlinkFile(file1);
    list2 = loadPlinkFile(file2);

    for (CNVariant element : list1) {
      for (CNVariant element2 : list2) {

        if (element.familyID.equals(element2.familyID)
            && element.individualID.equals(element2.individualID)
            && element.significantOverlap(element2, checkLarger)) {
          matched1.add(element);
          matched2.add(element2);
        }
      }
    }

    for (CNVariant element : list1) {
      if (matched1.contains(element)) {
        outputLines.add(element.toPlinkFormat() + "\t3");
      } else {
        unmatched1.add(element);
      }
    }

    for (int i = 0; i < list2.length; i++) {
      if (matched2.contains(list2[i])) {
        outputLines.add(list2[i].toPlinkFormat() + "\t3");
      } else {
        unmatched2.add(list1[i]);
      }
    }

    Iterator<CNVariant> unmatchedIterator = unmatched1.iterator();
    while (unmatchedIterator.hasNext()) {
      outputLines.add(unmatchedIterator.next().toPlinkFormat() + "\t1");
    }
    unmatchedIterator = unmatched2.iterator();
    while (unmatchedIterator.hasNext()) {
      outputLines.add(unmatchedIterator.next().toPlinkFormat() + "\t2");
    }

    try {
      writer = Files.openAppropriateWriter(ext.rootOf(file1) + "_" + ext.rootOf(file2)
                                           + "_signif_consensus.cnv");
      for (int i = 0; i < outputLines.size(); i++) {
        writer.println(outputLines.get(i));
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing significant-consensus file");
      e.printStackTrace();
    }

  }

  public static CNVariant[] putInOrder(CNVariant[] array, int[] order) {
    CNVariant[] newArray;

    newArray = new CNVariant[array.length];
    for (int i = 0; i < order.length; i++) {
      newArray[i] = array[order[i]];
    }

    return newArray;
  }

  /**
   * Copies and sorts the given array. For in-place sorting use {@link #sortInPlaceByQuality}
   *
   * @param array sort by the quality scores
   * @return Sorted array
   */
  public static <T extends CNVariant> T[] sortByQuality(final T[] array, final boolean reverse) {
    T[] toSort = Arrays.copyOf(array, array.length);
    return sortInPlaceByQuality(toSort, reverse);
  }

  public static <T extends CNVariant> T[] sortInPlaceByQuality(final T[] array,
                                                               final boolean reverse) {
    Arrays.sort(array, new Comparator<T>() {

      @Override
      public int compare(T o1, T o2) {
        T t1 = reverse ? o2 : o1;
        T t2 = reverse ? o1 : o2;
        return Double.compare(t1.getScore(), t2.getScore());
      }

    });
    return array;
  }

  @Override
  public CNVariant clone() {
    return new CNVariant(familyID, individualID, getChr(), getStart(), getStop(), cn, score,
                         numMarkers, source);
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename1 = "CNVariant1.dat";
    String filename2 = "CNVariant2.dat";
    boolean signif = false;

    String usage = "\n" + "cnv.var.CNVariant requires 2+ arguments\n"
                   + "   (1) first CNV filename (i.e. file1=" + filename1 + " (not the default))\n"
                   + "   (2) second CNV filename (i.e. file2=" + filename2 + " (not the default))\n"
                   + "   (3) (Optional) restrict to significant overlap (i.e. -sig (not the default))\n"
                   + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("file1=")) {
        filename1 = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("file2=")) {
        filename1 = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("-sig")) {
        signif = true;
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
      if (signif) {
        findSignificantConsensus(filename1, filename2, CONSENSUS_TYPE.NOT_CN_AWARE,
                                 OVERLAP_TYPE.OVERLAP_LOC_AND_INDIVIDUAL_SIGNIFICANTLY);
      } else {
        findConsensus(filename1, filename2);
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  // public static void main(String[] args) {
  // String file1 = "C:\\Documents and Settings\\npankrat\\My
  // Documents\\CNV\\penncnv\\again_noGenderProblems\\conf_0kb_5SNP_10.0.cnv";
  // String file2 = "C:\\Documents and Settings\\npankrat\\My
  // Documents\\CNV\\quantisnp\\noGenderProblems\\conf_0kb_5SNP_10.0.cnv";
  // // String file1 = "Penn_conf_0kb_5SNP_10.0.cnv";
  // // String file2 = "Quanti_conf_0kb_5SNP_10.0.cnv";
  //
  // try {
  // findConsensus(file1, file2);
  // } catch (Exception e) {
  // e.printStackTrace();
  // }
  // }

  // Return the Individual ID when getting as a string
  @Override
  public String toString() {
    return getIndividualID();
  }

  @Override
  public String toAnalysisString() {
    return toPlinkFormat();
  }

  @Override
  public String[] getHeader() {
    return PLINK_CNV_HEADER;
  }

  @SuppressWarnings("unchecked")
  public static LocusSet<CNVariant> readSerialCnvSet(String filename, Logger log) {
    return ((LocusSet<CNVariant>) SerializedFiles.readSerial(filename, log, false, true));
  }

  public static LocusSet<CNVariant> loadLocSet(String cnvFile, Logger log) {
    CNVariant[] cnvs = CNVariant.loadPlinkFile(cnvFile);
    LocusSet<CNVariant> cLocusSet = new LocusSet<>(cnvs, true, log);
    return cLocusSet;

  }

  public static HashSet<String> getUniqueInds(String cnvFile, Logger log) {
    LocusSet<CNVariant> cLocusSet = CNVariant.loadLocSet(cnvFile, log);

    return getUniqueInds(cLocusSet, log);
  }

  public static <T extends CNVariant> HashSet<String> getUniqueInds(LocusSet<T> cLocusSet,
                                                                    Logger log) {
    HashSet<String> uniqueInds = new HashSet<>();
    for (int i = 0; i < cLocusSet.getLoci().length; i++) {
      CNVariant tmp = cLocusSet.getLoci()[i];
      uniqueInds.add(tmp.getFamilyID() + "\t" + tmp.getIndividualID());
    }
    return uniqueInds;
  }

  /**
   * Will break the {@link LocusSet} into a map between FID/IID and {@link CNVariant}s belonging to
   * that sample
   * 
   * @param set
   * @param log
   * @return
   */
  public static <T extends CNVariant> Map<String, LocusSet<T>> breakIntoInds(LocusSet<T> set,
                                                                             Logger log) {
    HashMap<String, ArrayList<T>> cnvSplits = new HashMap<>();
    for (int i = 0; i < set.getLoci().length; i++) {
      T tmp = set.getLoci()[i];
      String key = tmp.getFamilyID() + "\t" + tmp.getIndividualID();
      if (!cnvSplits.containsKey(key)) {
        cnvSplits.put(key, new ArrayList<T>());
      }
      cnvSplits.get(key).add(tmp);
    }
    HashMap<String, LocusSet<T>> setSplit = new HashMap<>();

    for (String fidIid : cnvSplits.keySet()) {
      ArrayList<T> indCNVs = cnvSplits.get(fidIid);
      if (!indCNVs.isEmpty()) {
        LocusSet<T> indSet = new LocusSet<T>(indCNVs, true, log) {

          /**
           * 
           */
          private static final long serialVersionUID = 1L;

        };
        setSplit.put(fidIid, indSet);
      }
    }
    return setSplit;
  }

  /**
   * Useful when the terms of the cnv are not completely known, such as in CNVCalling
   */
  public static class CNVBuilder {

    private byte chr = -1;
    private int start = -1;
    private int stop = -1;
    private String familyID = null;
    private String individualID = null;
    private int cn = -22;
    private double score = -22;
    private int numMarkers = -1;
    private int source = -1;
    private int startIndex = -1;// for tracking external indices, such as positions in project
    private int stopIndex = -1;// for tracking external indices, such as positions in project

    public CNVBuilder() {

    }

    public CNVBuilder(CNVariant cnVariant) {
      chr = cnVariant.getChr();
      start = cnVariant.getStart();
      stop = cnVariant.getStop();
      familyID = cnVariant.familyID;
      individualID = cnVariant.individualID;
      cn = cnVariant.cn;
      score = cnVariant.score;
      numMarkers = cnVariant.numMarkers;
      source = cnVariant.source;
    }

    public CNVBuilder familyID(String familyID) {
      this.familyID = familyID;
      return this;
    }

    public CNVBuilder individualID(String individualID) {
      this.individualID = individualID;
      return this;
    }

    public CNVBuilder start(int start) {
      this.start = start;
      return this;
    }

    public CNVBuilder stop(int stop) {
      this.stop = stop;
      return this;
    }

    public CNVBuilder startIndex(int startIndex) {
      this.startIndex = startIndex;
      return this;
    }

    public int getStartIndex() {
      return startIndex;
    }

    public int getStopIndex() {
      return stopIndex;
    }

    public CNVBuilder stopIndex(int stopIndex) {
      this.stopIndex = stopIndex;
      return this;
    }

    public CNVBuilder chr(byte chr) {
      this.chr = chr;
      return this;
    }

    public CNVBuilder cn(int cn) {
      this.cn = cn;
      return this;
    }

    public CNVBuilder score(double score) {
      this.score = score;
      return this;
    }

    public CNVBuilder numMarkers(int numMarkers) {
      this.numMarkers = numMarkers;
      return this;
    }

    public CNVBuilder source(int source) {
      this.source = source;
      return this;
    }

    public int getStart() {
      return start;
    }

    public int getStop() {
      return stop;
    }

    public int getCn() {
      return cn;
    }

    public int getNumMarkers() {
      return numMarkers;
    }

    public CNVariant build() {
      return new CNVariant(this);
    }
  }

  protected CNVariant(CNVBuilder builder) {
    super(builder.chr, builder.start, builder.stop);
    familyID = builder.familyID;
    individualID = builder.individualID;
    cn = builder.cn;
    score = builder.score;
    numMarkers = builder.numMarkers;
    source = builder.source;
  }

}
