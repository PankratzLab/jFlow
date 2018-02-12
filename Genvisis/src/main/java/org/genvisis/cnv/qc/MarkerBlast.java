package org.genvisis.cnv.qc;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import org.genvisis.cnv.annotation.markers.AnalysisParams;
import org.genvisis.cnv.annotation.markers.AnnotationData;
import org.genvisis.cnv.annotation.markers.AnnotationFileWriter;
import org.genvisis.cnv.annotation.markers.BlastAnnotationTypes.PROBE_TAG;
import org.genvisis.cnv.annotation.markers.BlastAnnotationTypes.TOP_BOT;
import org.genvisis.cnv.annotation.markers.BlastAnnotationWriter;
import org.genvisis.cnv.annotation.markers.BlastParams;
import org.genvisis.cnv.annotation.markers.LocusAnnotation;
import org.genvisis.cnv.annotation.markers.LocusAnnotation.Builder;
import org.genvisis.cnv.annotation.markers.MarkerGCAnnotation;
import org.genvisis.cnv.filesys.MarkerDetailSet;
import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.genvisis.cnv.manage.ExtProjectDataParser;
import org.genvisis.cnv.manage.Markers;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.GenomicPosition;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.Positions;
import org.genvisis.common.WorkerHive;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.analysis.Blast;
import org.genvisis.seq.analysis.Blast.BlastWorker;
import org.genvisis.seq.analysis.Blast.FastaEntry;
import org.genvisis.seq.manage.ReferenceGenome;
import org.genvisis.seq.manage.StrandOps;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.Allele;

/**
 * As starts handling more than just annotatating with Blast results, we should maybe re-name this
 * class
 */
public abstract class MarkerBlast {

  protected static final String ARG_BLAST_WORD_SIZE = "blastWordSize";
  protected static final String DESC_BLAST_WORD_SIZE = "word size for initial match in blast db";

  protected static final String ARG_REPORT_WORD_SIZE = "reportWordSize";
  protected static final String DESC_REPORT_WORD_SIZE = "minimum number of base pairs with an exact match to report";

  protected static final String ARG_MAX_ALIGNMENTS = "maxAlignmentsReported";
  protected static final String DESC_MAX_ALIGNMENTS = "the maximum number of alignments to summarize";
  protected static final int DEFAULT_MAX_ALIGNMENTS = 20;

  protected static final String ARG_MARKER_POSITIONS_OVERRIDE = "overrideMarkerPositions";
  protected static final String DESC_MARKER_POSITIONS_OVERRIDE = "marker positions file to override MarkerSet positions with (e.g. for a different genome build)";
  protected static final String EXAMPLE_MARKER_POSITIONS_OVERRID = "markerPositions.txt";

  protected static final String FLAG_SKIP_GC_ANNOTATION = "noGCAnno";
  protected static final String DESC_SKIP_GC_ANNOTATION = "skip annotating the summary file with probe and regional GC content";
  protected static final boolean DEFAULT_ANNOTATE_GC_CONTENT = true;

  protected static final String FLAG_SKIP_BLAST = "noBlast";
  protected static final String DESC_SKIP_BLAST = "skip the blast analysis and simply add the probe annotations";
  protected static final boolean DEFAULT_DO_BLAST = true;

  protected static final boolean DEFAULT_REPORT_TO_TEMPORARY_FILE = true;

  public enum FILE_SEQUENCE_TYPE {
    // /**
    // * like HumanOmni2-5-8-v1-2-A-Strand-Report-FDT.txt
    // */
    // STRAND_REPORT,
    /**
     * Like HumanExome-12-v1-0-B.csv
     */
    MANIFEST_FILE,
    /**
     * Like GenomeWideSNP_6.probe_tab
     */
    AFFY_ANNOT;
  }

  protected final Project proj;
  protected final int blastWordSize;
  protected final int reportWordSize;
  protected final int maxAlignmentsReported;
  protected final boolean reportToTmp;
  protected final boolean annotateGCContent;
  protected final boolean doBlast;
  protected final int numThreads;
  protected final Logger log;
  private final Map<String, GenomicPosition> overrideMarkerPositions;
  private MarkerDetailSet naiveMarkerSet = null;

  /**
   * @param proj
   * @param blastWordSize
   * @param reportWordSize
   * @param maxAlignmentsReported
   * @param reportToTmp
   * @param annotateGCContent
   * @param doBlast
   * @param numThreads
   */
  protected MarkerBlast(Project proj, int blastWordSize, int reportWordSize,
                        int maxAlignmentsReported, boolean reportToTmp, boolean annotateGCContent,
                        boolean doBlast, int numThreads) {
    super();
    this.proj = new Project(proj.getPropertyFilename()) {

      @Override
      public synchronized MarkerDetailSet getMarkerSet() {
        if (naiveMarkerSet == null) {
          naiveMarkerSet = loadNaiveMarkerSet();
        }
        return naiveMarkerSet;
      };
    };
    this.blastWordSize = blastWordSize;
    this.reportWordSize = reportWordSize;
    this.maxAlignmentsReported = maxAlignmentsReported;
    this.reportToTmp = reportToTmp;
    this.annotateGCContent = annotateGCContent;
    this.doBlast = doBlast;
    this.numThreads = numThreads;
    this.log = proj.getLog();
    this.overrideMarkerPositions = Maps.newHashMap();
  }

  protected abstract String getSourceString();

  protected abstract String getNameBase();

  @SuppressWarnings("deprecation")
  private MarkerDetailSet loadNaiveMarkerSet() {
    // This method intentionally accesses proj.MARKERSET_FILENAME to circumvent proj.getMarkerSet()
    // and return the naive MarkerSet, with positions not based on any previous BLAST annotation
    String markerSetFile = proj.MARKERSET_FILENAME.getValue();
    if (!Files.exists(markerSetFile)) {
      throw new IllegalStateException(markerSetFile + " does not exist. Cannot run "
                                      + this.getClass().getName() + " without a "
                                      + MarkerSet.class.getName());
    }
    MarkerDetailSet naiveMarkerSet = new MarkerDetailSet(MarkerSet.load(markerSetFile));
    if (overrideMarkerPositions.isEmpty()) {
      return naiveMarkerSet;
    }
    List<MarkerDetailSet.Marker> markers = Lists.newArrayListWithCapacity(naiveMarkerSet.getMarkers()
                                                                                        .size());

    for (MarkerDetailSet.Marker naiveMarker : naiveMarkerSet.getMarkers()) {
      String name = naiveMarker.getName();
      MarkerDetailSet.Marker marker;
      if (overrideMarkerPositions.containsKey(name)) {
        marker = new MarkerDetailSet.Marker(name, overrideMarkerPositions.get(name));
      } else {
        marker = naiveMarker;
      }
      markers.add(marker);
    }
    // This sort makes the fingerprint not match that of the underlying MarkerSet but is necessary
    // in order to generate a sorted VCF. This could be fixed downstream or the fingerprint should
    // be replaced by a better metric of matching to a project
    Collections.sort(markers);
    return new MarkerDetailSet(markers);

  }

  public void clearMarkerPositionOverrides() {
    naiveMarkerSet = null;
    overrideMarkerPositions.clear();
  }

  public void overrideMarkerPositions(String markerPositionsFile) {
    naiveMarkerSet = null;
    Map<String, String> markerStringMap = Markers.loadFileToHashString(markerPositionsFile, log);
    for (Map.Entry<String, String> markerEntry : markerStringMap.entrySet()) {
      String name = markerEntry.getKey();
      String[] chrPos = markerEntry.getValue().split(PSF.Regex.GREEDY_WHITESPACE);
      byte chr = Positions.chromosomeNumber(chrPos[0], log);
      int position = Integer.parseInt(chrPos[1]);
      if (overrideMarkerPositions.put(name, new GenomicPosition(chr, position)) != null) {
        log.reportTimeWarning("Marker " + name + " has been overriden more than once");
      }
    }

  }

  /**
   * {@link MarkerBlast#blastEm(Project, String, FILE_SEQUENCE_TYPE, int, int, int, int, boolean, boolean, boolean)}
   */
  public MarkerBlastResult blastEm() {
    String fastaDb = proj.getReferenceGenomeFASTAFilename();
    if (!Files.exists(fastaDb) && doBlast) {
      log.reportError("Unable to find or obtain reference genome");
      return null;
    } else {
      double evalueCutoff = 10000;
      BlastParams blastParams = new BlastParams(getSourceString(), fastaDb, maxAlignmentsReported,
                                                reportWordSize, blastWordSize,
                                                ext.getTimestampForFilename(), evalueCutoff,
                                                proj.getMarkerSet().getFingerprint(), "", log);
      Blast blast = new Blast(fastaDb, blastWordSize, reportWordSize, log, true, true);
      blast.setEvalue(evalueCutoff);// we rely on the wordSize instead
      String dir = proj.PROJECT_DIRECTORY.getValue() + "Blasts/";
      new File(dir).mkdirs();

      String root = dir + getNameBase() + "." + ext.rootOf(proj.getReferenceGenomeFASTAFilename())
                    + ".blasted.ws." + blastWordSize + ".rep." + reportWordSize;
      String output = root + ".blasted";

      String[] tmps = new String[numThreads];
      for (int i = 0; i < numThreads; i++) {
        tmps[i] = root + ".tmp" + i;
        if (!Files.exists(tmps[i])) {
          tmps[i] = tmps[i] + ".gz";
        }
        // log.reportTimeInfo("REmember gz");
      }

      if (doBlast && !Files.exists("", tmps)) {
        MarkerFastaEntry[] fastaEntries = getMarkerFastaEntries(blastParams, false);
        List<MarkerFastaEntry[]> splits = ArrayUtils.splitUpArray(fastaEntries, numThreads, log);

        ArrayList<BlastWorker> workers = new ArrayList<Blast.BlastWorker>();
        if (fastaEntries != null && fastaEntries.length > 0) {
          for (int i = 0; i < splits.size(); i++) {
            if (!Files.exists(tmps[i])) {
              workers.add(new BlastWorker(blast, splits.get(i), reportToTmp ? tmps[i] : null));
            } else {
              log.reportTimeWarning("Skipping index " + i + ", " + tmps[i] + " exists");
            }
          }
        }

        if (workers.size() > 0) {
          WorkerHive<Blast.BlastResultsSummary[]> hive = new WorkerHive<Blast.BlastResultsSummary[]>(numThreads,
                                                                                                     10,
                                                                                                     log);
          hive.addCallables(workers.toArray(new BlastWorker[workers.size()]));
          hive.setReportEvery(1);
          hive.execute(true);
          hive.getResults();
        }
      } else {
        for (String tmp : tmps) {
          if (doBlast) {
            log.reportFileExists(tmp);
          }
        }
      }
      log.reportTimeWarning("Assuming that all sequences have length "
                            + proj.getArrayType().getProbeLength() + " based on array type "
                            + proj.getArrayType());
      MarkerBlastResult result = new MarkerBlastResult(output, blastWordSize, reportWordSize,
                                                       proj.getArrayType().getProbeLength());
      result.setTmpFiles(tmps);

      final String blastAnnotationFile = proj.BLAST_ANNOTATION_FILENAME.getValue();
      if (Files.exists(blastAnnotationFile)) {
        String blastAnnotationDir = ext.parseDirectoryOfFile(blastAnnotationFile);
        Files.backup(ext.removeDirectoryInfo(blastAnnotationFile), blastAnnotationDir,
                     blastAnnotationDir, true);
      }
      MarkerFastaEntry[] entries = getMarkerFastaEntries(blastParams, true);

      log.reportTimeInfo("Summarizing blast results to " + blastAnnotationFile);
      BlastAnnotationWriter blastAnnotationWriter = new BlastAnnotationWriter(proj,
                                                                              new AnalysisParams[] {blastParams},
                                                                              blastAnnotationFile,
                                                                              doBlast ? tmps
                                                                                      : new String[] {},
                                                                              entries,
                                                                              reportWordSize,
                                                                              proj.getArrayType()
                                                                                  .getProbeLength(),
                                                                              proj.getArrayType()
                                                                                  .getProbeLength(),
                                                                              maxAlignmentsReported);

      if (proj.getArrayType() != ARRAY.ILLUMINA) {
        log.reportTimeWarning("Did not detect array type " + ARRAY.ILLUMINA
                              + " , probe sequence annotation may not reflect the true design since multiple designs may be reported");
      }
      blastAnnotationWriter.summarizeResultFiles(false);
      blastAnnotationWriter.close();
      if (annotateGCContent) {
        annotateGCContent();
      }

      // } else {
      // log.reportFileExists(proj.BLAST_ANNOTATION_FILENAME.getValue());
      // }
      return result;
    }
  }

  public static int getDefaultWordSize(Project proj) {
    return getDefaultWordSize(proj.getArrayType());
  }

  public static int getDefaultWordSize(ARRAY arrayType) {
    switch (arrayType) {
      case AFFY_GW6:
      case AFFY_GW6_CN:
        return 15;
      case ILLUMINA:
        return 25;
      default:
        throw new IllegalArgumentException("Invalid array type " + arrayType);
    }
  }

  /**
   * Adds the gc content annotation to the summary file, does some repeat loading but oh well.
   */
  protected void annotateGCContent() {
    MarkerFastaEntry[] fastaEntries = getMarkerFastaEntries(null, false);
    MarkerDetailSet markerSet = proj.getMarkerSet();
    String[] markerNames = markerSet.getMarkerNames();
    Map<String, Integer> indices = markerSet.getMarkerIndices();
    // ReferenceGenome referenceGenome = new ReferenceGenome(fastaDb, log);
    LocusAnnotation[] gcAnnotations = new LocusAnnotation[markerNames.length];
    for (MarkerFastaEntry fastaEntrie : fastaEntries) {
      String marker = fastaEntrie.getName();
      PROBE_TAG tag = PROBE_TAG.parseMarkerTag(marker, log);
      marker = marker.substring(0, marker.length() - tag.getTag().length());
      int index = indices.get(marker);
      double gcContent = fastaEntrie.getGCMinusInterrogationPosition();
      Builder builder = new Builder();
      AnnotationData annotationData = MarkerGCAnnotation.getGCAnnotationDatas();
      annotationData.setData(gcContent + "");
      builder.annotations(new AnnotationData[] {annotationData});
      MarkerGCAnnotation markerGCAnnotation = new MarkerGCAnnotation(builder, marker,
                                                                     fastaEntrie.getMarkerSegment());
      gcAnnotations[index] = markerGCAnnotation;
    }
    AnnotationFileWriter writer = new AnnotationFileWriter(proj, null,
                                                           new AnnotationData[] {MarkerGCAnnotation.getGCAnnotationDatas()},
                                                           proj.BLAST_ANNOTATION_FILENAME.getValue(),
                                                           false) {};
    for (int i = 0; i < gcAnnotations.length; i++) {
      if (gcAnnotations[i] != null) {// currently some markers may not be represented, such as affy
        // CN markers
        writer.write(gcAnnotations[i], false, true);

      } else {
        double gcContent = Double.NaN;
        Builder builder = new Builder();
        AnnotationData annotationData = MarkerGCAnnotation.getGCAnnotationDatas();
        annotationData.setData(gcContent + "");
        builder.annotations(new AnnotationData[] {annotationData});
        MarkerGCAnnotation blank = new MarkerGCAnnotation(builder, markerNames[i],
                                                          new Segment(markerSet.getChrs()[i],
                                                                      markerSet.getPositions()[i],
                                                                      markerSet.getPositions()[i]));
        writer.write(blank, false, true);
      }
    }
    writer.close();

  }

  public static class MarkerBlastResult {

    private final String output;
    private String[] tmpFiles;
    private final int blastWordSize;
    private final int reportWordSize;
    private final int sequenceSize;

    public MarkerBlastResult(String output, int blastWordSize, int reportWordSize,
                             int sequenceSize) {
      super();
      this.output = output;
      this.blastWordSize = blastWordSize;
      this.reportWordSize = reportWordSize;
      this.sequenceSize = sequenceSize;
    }

    public int getSequenceSize() {
      return sequenceSize;
    }

    public void setTmpFiles(String[] tmpFiles) {
      this.tmpFiles = tmpFiles;
    }

    public String[] getTmpFiles() {
      return tmpFiles;
    }

    public int getBlastWordSize() {
      return blastWordSize;
    }

    public int getReportWordSize() {
      return reportWordSize;
    }

    public String getOutput() {
      return output;
    }
  }

  protected static String getRefStrand(ExtProjectDataParser parser) {
    String strand = null;
    if (parser.hasStringDataForTitle("RefStrand")) {
      strand = "RefStrand";
    } else if (parser.hasStringDataForTitle("GenomicStrand")) {
      strand = "GenomicStrand";
    }
    return strand;
  }

  /**
   * {@link MarkerBlast#getMarkerFastaEntries(BlastParams, boolean)}
   * 
   * @param blastParams TODO
   */
  protected abstract MarkerFastaEntry[] getMarkerFastaEntries(BlastParams blastParams,
                                                              boolean alleleLookup);

  protected static String parseStrandFromIlmnID(String[] fullId, TOP_BOT topBotRef) {
    String refStrand = null;
    String tb = fullId[fullId.length - 3];
    String fr = fullId[fullId.length - 2];
    if ((!tb.equals("B") && !tb.equals("T") && !tb.equals("M"))
        || (!fr.equals("F") && !fr.equals("R"))) {
      throw new IllegalStateException("Invalid IlmnID parsing asumption, splitting with _ for "
                                      + ArrayUtils.toStr(fullId, "_") + "\t" + tb + "\t" + fr
                                      + "\nPlease update TOP/BOT strand designation");
    }

    switch (topBotRef) {
      case BOT:
      case MINUS:
        if (((tb.equals("B") || tb.equals("M")) && fr.equals("F"))
            || tb.equals("T") && fr.equals("R")) {
          refStrand = "+";
        } else if (((tb.equals("B") || tb.equals("M")) && fr.equals("R"))
                   || tb.equals("T") && fr.equals("F")) {
          refStrand = "-";
        }

        break;
      case PLUS:
      case TOP:
        if ((tb.equals("B") && fr.equals("R")) || tb.equals("T") && fr.equals("F")) {
          refStrand = "+";
        } else if ((tb.equals("B") && fr.equals("F")) || tb.equals("T") && fr.equals("R")) {
          refStrand = "-";
        }
        break;
      case NA:
        break;
      default:
        break;

    }
    if (refStrand == null) {
      throw new IllegalStateException("Could not parse ref TOP/BOT  " + topBotRef + " using IlmnID"
                                      + ArrayUtils.toStr(fullId, "_") + " to a strand");
    }
    return refStrand;
  }

  /**
   * @author lane0212 Parse alleles to ref/alt/A/B with strand flip for alts if needed
   */
  static class AlleleParser {

    /**
     * 
     */
    private final Segment loc;
    private final String aS;
    private final String bS;
    private final String pseqB;
    private final String markerName;
    private final ReferenceGenome referenceGenome;
    private Allele A;
    private Allele B;
    private Allele ref;
    private Allele[] alts;

    protected AlleleParser(String markerName, Segment loc, String aS, String bS, String pseqB,
                           ReferenceGenome referenceGenome) {
      super();
      this.markerName = markerName;
      this.loc = loc;
      this.aS = aS;
      this.pseqB = pseqB;
      this.bS = bS;
      this.referenceGenome = referenceGenome;
    }

    public Allele getA() {
      return A;
    }

    public Allele getB() {
      return B;
    }

    public Allele getRef() {
      return ref;
    }

    public Allele[] getAlts() {
      return alts;
    }

    private boolean isIndel() {
      if (aS.equals("I") || bS.equals("I")) {
        return true;
      }
      if (aS.equals("D") || bS.equals("D")) {
        return true;
      }
      return false;
    }

    private void parseIndel(ARRAY array, Strand strand, String sourceSeq) {
      if (sourceSeq == null) {
        throw new IllegalArgumentException("Must have source seq to get indel");
      } else if (!pseqB.equals("")) {
        throw new IllegalArgumentException("Assumptions violated for single probe indel ");

      } else if (!(aS.equals("I") && bS.equals("D")) && !(aS.equals("D") && bS.equals("I"))) {
        throw new IllegalArgumentException("Assumptions violated for insertion deletion codes " + aS
                                           + " and " + bS);

      } else {

        int start = sourceSeq.indexOf("[") + 1;
        int stop = sourceSeq.indexOf("]");

        String indel = sourceSeq.substring(start, stop);
        if (!indel.startsWith("-/") && !indel.startsWith("+/")) {
          throw new IllegalArgumentException("Indel source seq must start with - or + for "
                                             + indel);
        }
        indel = indel.substring(2);
        int len = indel.length() - 1;
        if (loc.getChr() > 0) {

          Segment newQuery = new Segment(loc.getChr(), loc.getStart() - 1, loc.getStop() + len);// need
                                                                                                // to
                                                                                                // grab
                                                                                                // preceeding
                                                                                                // ref
                                                                                                // bp
          String[] tmp = referenceGenome == null ? ArrayUtils.stringArray(indel.length() + 1, "N")
                                                 : referenceGenome.getSequenceFor(newQuery);
          if (tmp.length - 1 != indel.length()) {
            throw new IllegalStateException("Invalid reference indel query; REF -> "
                                            + ArrayUtils.toStr(tmp) + "\t indel -> " + indel);
          }
          boolean insertionisRef = true;
          for (int i = 0; i < indel.length(); i++) {
            if (!tmp[i + 1].equalsIgnoreCase(indel.charAt(i) + "")
                && !tmp[i + 1].equalsIgnoreCase(StrandOps.flipIfNeeded(indel.charAt(i) + "", strand,
                                                                       false))) {
              insertionisRef = false;
              break;
            }
          }
          ref = Allele.create(ArrayUtils.toStr(tmp, ""), true);
          if (insertionisRef) {
            if (aS.equals("I")) {

              String tmpA = ArrayUtils.toStr(tmp, "");
              String tmpB = tmp[0];
              A = Allele.create(StrandOps.flipsIfNeeded(tmpA, strand, false), false);
              B = Allele.create(StrandOps.flipsIfNeeded(tmpB, strand, false), false);
              alts = new Allele[1];
              alts[0] = Allele.create(tmp[0], false);

            } else {
              String tmpA = tmp[0];
              String tmpB = ArrayUtils.toStr(tmp, "");
              A = Allele.create(StrandOps.flipsIfNeeded(tmpA, strand, false), false);
              B = Allele.create(StrandOps.flipsIfNeeded(tmpB, strand, false), false);
              alts = new Allele[1];
              alts[0] = Allele.create(tmp[0], false);
            }
          } else {
            if (aS.equals("I")) {
              String tmpA = StrandOps.flipIfNeeded(tmp[0], strand, false) + indel;
              String tmpB = StrandOps.flipIfNeeded(tmp[0], strand, false);
              A = Allele.create(tmpA, false);
              B = Allele.create(tmpB, false);
              alts = new Allele[2];
              alts[0] = Allele.create(tmp[0], false);
              alts[1] = Allele.create(tmp[0] + StrandOps.flipsIfNeeded(indel, strand, false),
                                      false);
            } else {
              String tmpA = StrandOps.flipIfNeeded(tmp[0], strand, false);
              String tmpB = StrandOps.flipIfNeeded(tmp[0], strand, false) + indel;
              A = Allele.create(tmpA, false);
              B = Allele.create(tmpB, false);
              alts = new Allele[2];
              alts[0] = Allele.create(tmp[0], false);
              alts[1] = Allele.create(tmp[0] + StrandOps.flipsIfNeeded(indel, strand, false),
                                      false);
            }
          }
        } else {

          ref = Allele.create("NN", true);
          if (aS.equals("I")) {
            String tmpA = "N" + indel;
            String tmpB = "N";
            A = Allele.create(tmpA, false);
            B = Allele.create(tmpB, false);
            alts = new Allele[2];
            alts[0] = Allele.create("N" + StrandOps.flipsIfNeeded(indel, strand, false), false);
            alts[1] = Allele.create("N", false);
          } else {
            String tmpA = "N";
            String tmpB = "N" + indel;
            A = Allele.create(tmpA, false);
            B = Allele.create(tmpB, false);
            alts = new Allele[2];
            alts[0] = Allele.create("N" + StrandOps.flipsIfNeeded(indel, strand, false), false);
            alts[1] = Allele.create("N", false);
          }
        }
      }
    }

    protected void parse(ARRAY array, Strand strand, String sourceSeq) {
      if (loc.getChr() > 0) {
        if (!isIndel()) {
          String[] tmp;
          if (referenceGenome == null || loc.getStart() > referenceGenome.getContigLength(loc)) {
            tmp = new String[] {"N"};
          } else {
            tmp = referenceGenome.getSequenceFor(loc);
          }
          if (tmp.length != 1) {// don't think we need multiple
            throw new IllegalArgumentException("base query must be length one ("
                                               + loc.getUCSClocation() + " returned "
                                               + ArrayUtils.toStr(tmp));
          } else if (array.isCNOnly(markerName)) {

            ref = Allele.create(tmp[0], true);
            alts = new Allele[1];
            alts[0] = getCNPAllele(array);// Symbolic allele
            A = Allele.create(alts[0], false);
            B = Allele.create(alts[0], false);

          } else {
            String flipA = StrandOps.flipIfNeeded(aS, strand, false);
            String flipB = StrandOps.flipIfNeeded(bS, strand, false);

            ref = Allele.create(tmp[0], true);
            A = Allele.create(aS, false);
            B = Allele.create(bS, false);
            if (Allele.create(flipA, false).equals(ref, true)) {// A is ref
              alts = new Allele[1];
              alts[0] = Allele.create(flipB);
            } else if (Allele.create(flipB, false).equals(ref, true)) {// B is ref
              alts = new Allele[1];
              alts[0] = Allele.create(flipA);
            } else {// neither are ref
              alts = new Allele[2];
              alts[0] = Allele.create(flipA, false);
              alts[1] = Allele.create(flipB, false);
            }
          }
        } else {
          parseIndel(array, strand, sourceSeq);
        }

      } else {
        if (array.isCNOnly(markerName)) {

          ref = Allele.create("N", true);
          alts = new Allele[1];
          alts[0] = getCNPAllele(array);// Symbolic allele
          A = Allele.create(alts[0], false);
          B = Allele.create(alts[0], false);

        } else if (isIndel()) {
          parseIndel(array, strand, sourceSeq);
        } else {
          ref = Allele.create("N", true);
          alts = new Allele[2];
          alts[0] = Allele.create(StrandOps.flipIfNeeded(aS, strand, false), false);
          alts[1] = Allele.create(StrandOps.flipIfNeeded(bS, strand, false), false);
          A = Allele.create(aS, false);
          B = Allele.create(bS, false);
        }
      }
    }

    private Allele getCNPAllele(ARRAY array) {
      return Allele.create("<" + array + "_CNP>", false);
    }
  }

  public static class MarkerFastaEntry extends FastaEntry {

    private final int interrogationPosition;// Illumina = sequence length, affy = somewhere in the
                                            // middle
    private final Segment markerSegment;
    private final String seqB;
    private final Strand strand;
    private final TOP_BOT topBotProbe;
    private final TOP_BOT topBotRef;
    private final Allele A;
    private final Allele B;
    private final Allele ref;
    private final Allele[] alts;

    // IlmnStrand
    public MarkerFastaEntry(String name, String seqA, String seqB, Strand strand,
                            int interrogationPosition, Segment markerSegment, TOP_BOT topBotProbe,
                            TOP_BOT topBotRef, Allele A, Allele B, Allele ref, Allele[] alts) {
      super(name, seqA);
      this.seqB = seqB;
      this.strand = strand;
      this.interrogationPosition = interrogationPosition;
      this.markerSegment = markerSegment;
      this.topBotProbe = topBotProbe;
      this.topBotRef = topBotRef;
      this.A = A;
      this.B = B;
      this.ref = ref;
      this.alts = alts;
    }

    public String getSeqB() {
      return seqB;
    }

    public Allele getRef() {
      return ref;
    }

    public Allele[] getAlts() {
      return alts;
    }

    public Strand getStrand() {
      return strand;
    }

    public Allele getA() {
      return A;
    }

    public Allele getB() {
      return B;
    }

    public Segment getMarkerSegment() {

      return markerSegment;
    }

    public int getInterrogationPosition() {
      return interrogationPosition;
    }

    public TOP_BOT getTopBotProbe() {
      return topBotProbe;
    }

    public TOP_BOT getTopBotRef() {
      return topBotRef;
    }

    public double getGCMinusInterrogationPosition() {
      int GorC = 0;
      int total = 0;
      for (int i = 0; i < sequence.length(); i++) {
        if (i != interrogationPosition) {
          total++;
          if ((sequence.charAt(i) + "").equalsIgnoreCase("G")
              || (sequence.charAt(i) + "").equalsIgnoreCase("C")) {
            GorC++;
          }
        }
      }
      return (double) GorC / total;
    }

  }

  protected static Project prepareDummyProject(String csv, FILE_SEQUENCE_TYPE type) {
    if (Files.exists(csv)) {
      String dir = ext.parseDirectoryOfFile(csv);
      Project proj = new Project();
      proj.PROJECT_DIRECTORY.setValue(dir);
      Logger log = proj.getLog();
      switch (type) {
        case AFFY_ANNOT:
          System.err.println("Not implemented for AFFY, yet");
          return null;
        // proj.ARRAY_TYPE.setValue(ARRAY.AFFY_GW6);
        case MANIFEST_FILE:
          proj.ARRAY_TYPE.setValue(ARRAY.ILLUMINA);
          log.reportTimeWarning("Extracting marker positions from " + csv);
          String[] markerNames = extractMarkerPositionsFromManifest(csv, proj.getArrayType(), type,
                                                                    proj.MARKER_POSITION_FILENAME.getValue(),
                                                                    ",", log);
          Markers.orderMarkers(markerNames, proj.MARKER_POSITION_FILENAME.getValue(),
                               proj.MARKERSET_FILENAME.getValue(true, true), log);

          break;
        default:
          break;
      }
      return proj;
    } else {
      System.err.println("could not find " + csv);
      return null;
    }
  }

  /**
   * @param csv a .csv manifest file <br>
   *          Example Head; <br>
   *          Illumina, Inc. <br>
   *          [Heading] <br>
   *          Descriptor File Name,HumanOmni1-Quad_v1-0-Multi_H.bpm <br>
   *          Assay Format,Infinium HD Super <br>
   *          Date Manufactured,5/2/2011 <br>
   *          Loci Count ,1134514 <br>
   *          [Assay] <br>
   *          IlmnID,Name,IlmnStrand,SNP,AddressA_ID,AlleleA_ProbeSeq,AddressB_ID,AlleleB_ProbeSeq,
   *          GenomeBuild,Chr,MapInfo,Ploidy,Species,Source,SourceVersion,SourceStrand,SourceSeq,
   *          TopGenomicSeq,BeadSetID,Exp_Clusters,Intensity_Only,RefStrand <br>
   * @param array must be {@link ARRAY#ILLUMINA}
   * @param type must be {@link FILE_SEQUENCE_TYPE#MANIFEST_FILE}
   * @param output the output file, typically a projects marker positions
   * @param delimiter , typically ","
   * @param log
   * @return a String[] of marker names found in the manifest file
   */
  public static String[] extractMarkerPositionsFromManifest(String csv, ARRAY array,
                                                            FILE_SEQUENCE_TYPE type, String output,
                                                            String delimiter, Logger log) {

    if (type != FILE_SEQUENCE_TYPE.MANIFEST_FILE || array != ARRAY.ILLUMINA) {
      throw new IllegalArgumentException("This method should only be used in preparing marker positions for an "
                                         + ARRAY.ILLUMINA + " array using a "
                                         + FILE_SEQUENCE_TYPE.MANIFEST_FILE);
    }
    String[] required = new String[] {"Name", "Chr", "MapInfo"};
    ArrayList<String> markerNames = new ArrayList<String>();
    try {
      BufferedReader reader = Files.getAppropriateReader(csv);
      boolean start = false;
      int[] extract = new int[required.length];
      PrintWriter writer = Files.openAppropriateWriter(output);
      while (reader.ready()) {
        String[] line = reader.readLine().trim().split(delimiter);

        if (!start
            && ArrayUtils.countIf(ext.indexFactors(required, line, true, log, false), -1) == 0) {
          start = true;
          extract = ext.indexFactors(required, line, true, log, false);
          writer.println("Name\tChr\tPosition");
        } else if (start) {
          String[] lineMP = null;
          try {
            lineMP = new String[required.length];
            lineMP[0] = line[extract[0]];
            lineMP[1] = line[extract[1]];
            lineMP[2] = line[extract[2]];
            if (lineMP[2].equals("0")) {
              lineMP[1] = "0";
            }
          } catch (ArrayIndexOutOfBoundsException aOfBoundsException) {
            log.reportTimeWarning("Skipping line " + ArrayUtils.toStr(line));
            lineMP = null;
          }
          if (lineMP != null && lineMP[0] != null) {
            markerNames.add(lineMP[0]);
            writer.println(ArrayUtils.toStr(lineMP));
          }
        }
      }

      reader.close();
      writer.close();
      if (!start) {
        throw new IllegalStateException("Could not find required header subset "
                                        + ArrayUtils.toStr(required, ",") + " in " + csv);
      }

    } catch (FileNotFoundException fnfe) {
      log.reportError("Error: file \"" + csv + "\" not found in current directory");
      return null;
    } catch (IOException ioe) {
      log.reportError("Error reading file \"" + csv + "\"");
      return null;
    }
    return ArrayUtils.toStringArray(markerNames);
  }

  // public static void test() {
  // String file = "/home/pankrat2/shared/aric_exome_chip/Blast/HumanExome-12-v1-0-B.csv";
  // String fastaDb = "/home/pankrat2/public/bin/ref/hg19_canonical.fa";
  // blastEm(new Project("/home/pankrat2/lanej/projects/aric_exome.properties", false), file,
  // FILE_SEQUENCE_TYPE.MANIFEST_FILE, fastaDb, 2);
  // }

  /**
   * @deprecated Use {@link IlluminaMarkerBlast#main(String...)} or
   *             {@link AffyMarkerBlast#main(String...)} instead
   */
  @Deprecated
  public static void main(String[] args) {
    System.err.println("This class has been abstracted with implemention split out by array.");
    System.err.println("Use " + IlluminaMarkerBlast.class.getName() + " or "
                       + AffyMarkerBlast.class.getName() + " instead.");
  }

}
