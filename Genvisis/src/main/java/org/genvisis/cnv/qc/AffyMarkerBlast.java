package org.genvisis.cnv.qc;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import org.genvisis.CLI;
import org.genvisis.CLI.Arg;
import org.genvisis.cnv.annotation.markers.BlastAnnotationTypes.PROBE_TAG;
import org.genvisis.cnv.annotation.markers.BlastAnnotationTypes.TOP_BOT;
import org.genvisis.cnv.annotation.markers.BlastParams;
import org.genvisis.cnv.filesys.MarkerSetInfo;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.genvisis.cnv.manage.ExtProjectDataParser;
import org.genvisis.cnv.manage.ExtProjectDataParser.ProjectDataParserBuilder;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.ReferenceGenome;
import com.google.common.collect.ImmutableMap;
import htsjdk.tribble.annotation.Strand;

public class AffyMarkerBlast extends MarkerBlast {

  private static final String PROBE_MARKER_NAME = "PROBESET_ID";
  private static final String PROBE_SEQUENCE = "PROBE_SEQUENCE";

  private static final String ANNOT_MARKER_NAME = "Probe Set ID";
  private static final String ANNOT_ALLELE_A = "Allele A";
  private static final String ANNOT_ALLELE_B = "Allele B";
  private static final String ANNOT_STRAND = "Strand";

  private static final ImmutableMap<String, Strand> STRAND_ANNOTS = ImmutableMap.of("+",
                                                                                    Strand.POSITIVE,
                                                                                    "-",
                                                                                    Strand.NEGATIVE);

  public static final String ARG_PROBE_FILE = "probeFile";
  public static final String DESC_PROBE_FILE = "an Affymetrix Probe Set file";
  public static final String EXAMPLE_PROBE_FILE = "GenomeWideSNP_6.probe_tab";

  public static final String ARG_ANNOT_FILE = "annotFile";
  public static final String DESC_ANNOT_FILE = "an Affymetrix Annotation file";
  public static final String EXAMPLE_ANNOT_FILE = "GenomeWideSNP_6.na35.annot.csv";

  private static final int EXAMPLE_WORD_SIZE = MarkerBlast.getDefaultWordSize(ARRAY.AFFY_GW6);

  private final String probeFile;
  private final String annotFile;

  public AffyMarkerBlast(Project proj, int numThreads, String probeFile, String annotFile) {
    this(proj, getDefaultWordSize(proj), getDefaultWordSize(proj), DEFAULT_MAX_ALIGNMENTS,
         DEFAULT_REPORT_TO_TEMPORARY_FILE, DEFAULT_ANNOTATE_GC_CONTENT, DEFAULT_DO_BLAST,
         numThreads, probeFile, annotFile);
  }

  /**
   * @param proj
   * @param blastWordSize
   * @param reportWordSize
   * @param maxAlignmentsReported
   * @param reportToTmp
   * @param annotateGCContent
   * @param doBlast
   * @param numThreads
   * @param probeFile
   * @param annotFile
   */
  public AffyMarkerBlast(Project proj, int blastWordSize, int reportWordSize,
                         int maxAlignmentsReported, boolean reportToTmp, boolean annotateGCContent,
                         boolean doBlast, int numThreads, String probeFile, String annotFile) {
    super(proj, blastWordSize, reportWordSize, maxAlignmentsReported, reportToTmp,
          annotateGCContent, doBlast, numThreads);
    if (proj.getArrayType() != ARRAY.AFFY_GW6 && proj.getArrayType() != ARRAY.AFFY_GW6_CN) {
      log.reportError("Array type was set to " + proj.getArrayType() + " and this file is for "
                      + ARRAY.AFFY_GW6 + " or " + ARRAY.AFFY_GW6_CN);
      throw new IllegalArgumentException("Cannot generate " + AffyMarkerBlast.class.getName()
                                         + " for non-Affy array");
    }
    this.probeFile = probeFile;
    this.annotFile = annotFile;
  }

  @Override
  protected String getSourceString() {
    return probeFile + "," + annotFile;
  }

  @Override
  protected String getNameBase() {
    return ext.rootOf(probeFile) + "_" + ext.rootOf(annotFile);
  }

  @Override
  protected MarkerFastaEntry[] getMarkerFastaEntries(BlastParams params, boolean alleleLookup) {
    int seqLength = ARRAY.AFFY_GW6.getProbeLength();
    ExtProjectDataParser probeFileParser;
    try {
      probeFileParser = probeFileParser().build(proj, probeFile);
    } catch (FileNotFoundException e) {
      log.reportFileNotFound(probeFile);
      e.printStackTrace();
      return null;
    }
    ExtProjectDataParser annotFileParser;
    try {
      annotFileParser = annotFileParser().build(proj, annotFile);
    } catch (FileNotFoundException e) {
      log.reportFileNotFound(annotFile);
      e.printStackTrace();
      return null;
    }

    probeFileParser.determineIndicesFromTitles();
    probeFileParser.loadData();
    annotFileParser.determineIndicesFromTitles();
    annotFileParser.loadData();
    ArrayList<MarkerFastaEntry> entries = new ArrayList<MarkerFastaEntry>(ArrayUtils.booleanArraySum(probeFileParser.getDataPresent()));
    MarkerSetInfo markerSet = proj.getMarkerSet();
    ReferenceGenome referenceGenome = Files.exists(proj.getReferenceGenomeFASTAFilename()) ? new ReferenceGenome(proj.getReferenceGenomeFASTAFilename(),
                                                                                                                 log)
                                                                                           : null;
    if (referenceGenome == null) {
      log.reportTimeWarning("A reference genome was not found");
    }
    String refStrandTitle = getRefStrand(probeFileParser);
    if (refStrandTitle == null && params != null) {
      params.setNotes("Warning, the positive and negative strands of the probe design are actually forward and reverse designations, due to parsing IlmnIDs ");
    }
    byte[] chrs = markerSet.getChrs();
    int[] positions = markerSet.getPositions();
    for (int i = 0; i < probeFileParser.getDataPresent().length; i++) {
      if (probeFileParser.getDataPresent()[i]) {
        if ((i + 1) % 10000 == 0) {
          if (alleleLookup && referenceGenome != null) {
            log.reportTimeInfo("Loaded " + (i + 1) + " reference alleles from "
                               + referenceGenome.getReferenceFasta());
          } else {
            log.reportTimeInfo("Loaded " + (i + 1) + " probes from " + probeFile);
          }
        }
        String markerName = probeFileParser.getDataToLoad()[i];
        Segment markerSegment = new Segment(chrs[i], positions[i], positions[i]);
        String[] tmpSeq = ArrayUtils.unique(probeFileParser.getStringDataForTitle(PROBE_SEQUENCE)[i].split("\t"));

        if (tmpSeq.length != 2) {
          log.reportError("Marker " + markerName + " did not have 2 unique probe designs");
          log.reportError("found the following " + markerName + "\t" + ArrayUtils.toStr(tmpSeq));
          return null;
        } else {
          if (tmpSeq[0].length() != seqLength || tmpSeq[1].length() != seqLength) {
            log.reportError("Sequence " + tmpSeq[0] + " or " + tmpSeq[1] + "  did not have length "
                            + proj.ARRAY_TYPE.getValue().getProbeLength());
            return null;
          }
          int interrogationPosition = -1;
          for (int j = 0; j < tmpSeq[0].length(); j++) {
            if (tmpSeq[0].charAt(j) != tmpSeq[1].charAt(j)) {
              if (interrogationPosition != -1) {
                log.reportError("Multiple interrogation position for " + markerName);
                return null;
              }
              interrogationPosition = j;
            }
          }
          Strand strand = STRAND_ANNOTS.get(annotFileParser.getStringDataForTitle(ANNOT_STRAND)[i]);
          if (strand == null) {
            log.reportError("Invalid Annot Strand for " + markerName + ": "
                            + annotFileParser.getStringDataForTitle(ANNOT_STRAND)[i]
                            + " , marker will not be annotated");
            continue;
          }
          // probeFileParser and annotFileParser should both be indexed by marker indices
          String aS = annotFileParser.getStringDataForTitle(ANNOT_ALLELE_A)[i];
          String bS = annotFileParser.getStringDataForTitle(ANNOT_ALLELE_B)[i];
          AlleleParser alleleParser = new AlleleParser(markerName, markerSegment, aS, bS, tmpSeq[1],
                                                       referenceGenome);
          if (alleleLookup) {
            alleleParser.parse(proj.getArrayType(), strand, null);
          }
          entries.add(new MarkerFastaEntry(markerName + PROBE_TAG.A.getTag(), tmpSeq[0], tmpSeq[1],
                                           strand, interrogationPosition, markerSegment, TOP_BOT.NA,
                                           TOP_BOT.NA, alleleParser.getA(), alleleParser.getB(),
                                           alleleParser.getRef(), alleleParser.getAlts()));

          // TODO: Figure out if this should be tmpSeq[0] for A (or B?)
          entries.add(new MarkerFastaEntry(markerName + PROBE_TAG.B.getTag(), tmpSeq[1], tmpSeq[1],
                                           strand, interrogationPosition, markerSegment, TOP_BOT.NA,
                                           TOP_BOT.NA, alleleParser.getA(), alleleParser.getB(),
                                           alleleParser.getRef(), alleleParser.getAlts()));
        }
      }
    }

    log.reportTimeInfo("Found " + entries.size() + " marker sequences");
    return entries.toArray(new MarkerFastaEntry[entries.size()]);
  }

  protected ProjectDataParserBuilder probeFileParser() {
    ExtProjectDataParser.ProjectDataParserBuilder builder = new ExtProjectDataParser.ProjectDataParserBuilder();
    builder.separator("\t");
    builder.dataKeyColumnName(PROBE_MARKER_NAME);
    builder.stringDataTitles(PROBE_SEQUENCE);
    builder.concatMultipleStringEntries(true);
    builder.sampleBased(false);
    builder.treatAllNumeric(false);
    builder.requireAll(false);
    builder.verbose(false);
    return builder;
  }

  protected ProjectDataParserBuilder annotFileParser() {
    ExtProjectDataParser.ProjectDataParserBuilder builder = new ExtProjectDataParser.ProjectDataParserBuilder();
    builder.separator(",");
    builder.dataKeyColumnName(ANNOT_MARKER_NAME);
    builder.headerFlags(ANNOT_MARKER_NAME, ANNOT_ALLELE_A, ANNOT_ALLELE_B, ANNOT_STRAND);
    builder.stringDataTitles(ANNOT_ALLELE_A, ANNOT_ALLELE_B, ANNOT_STRAND);
    builder.sampleBased(false);
    builder.treatAllNumeric(false);
    builder.requireAll(false);
    builder.verbose(false);
    return builder;
  }

  public static void main(String... args) {

    CLI c = new CLI(IlluminaMarkerBlast.class);

    c.addArg(CLI.ARG_PROJ, CLI.DESC_PROJ, CLI.EXAMPLE_PROJ, true, Arg.FILE);
    c.addArg(ARG_ANNOT_FILE, DESC_ANNOT_FILE, EXAMPLE_ANNOT_FILE, true, Arg.FILE);
    c.addArg(ARG_PROBE_FILE, DESC_PROBE_FILE, EXAMPLE_PROBE_FILE, true, Arg.FILE);
    c.addArgWithDefault(CLI.ARG_THREADS, CLI.DESC_THREADS, CLI.EXAMPLE_THREADS);
    c.addArgWithDefault(ARG_BLAST_WORD_SIZE, DESC_BLAST_WORD_SIZE, EXAMPLE_WORD_SIZE);
    c.addArgWithDefault(ARG_REPORT_WORD_SIZE, DESC_REPORT_WORD_SIZE, EXAMPLE_WORD_SIZE);
    c.addArgWithDefault(ARG_MAX_ALIGNMENTS, DESC_MAX_ALIGNMENTS, DEFAULT_MAX_ALIGNMENTS);
    c.addArg(ARG_MARKER_POSITIONS_OVERRIDE, DESC_MARKER_POSITIONS_OVERRIDE,
             EXAMPLE_MARKER_POSITIONS_OVERRID, CLI.Arg.FILE);
    c.addFlag(FLAG_SKIP_GC_ANNOTATION, DESC_SKIP_GC_ANNOTATION);
    c.addFlag(FLAG_SKIP_BLAST, DESC_SKIP_BLAST);

    c.parseWithExit(args);

    Project proj = new Project(c.get(CLI.ARG_PROJ));
    String probeFile = c.get(ARG_PROBE_FILE);
    String annotFile = c.get(ARG_ANNOT_FILE);
    int blastWordSize = c.getI(ARG_BLAST_WORD_SIZE);
    int reportWordSize = c.getI(ARG_REPORT_WORD_SIZE);
    int maxAlignmentsReported = c.getI(ARG_MAX_ALIGNMENTS);
    int numThreads = c.getI(CLI.ARG_THREADS);
    String markerPositionsOverride = c.get(ARG_MARKER_POSITIONS_OVERRIDE);
    boolean annotateGCContent = !c.has(FLAG_SKIP_GC_ANNOTATION);
    boolean doBlast = !c.has(FLAG_SKIP_BLAST);

    try {
      MarkerBlast markerBlast = new AffyMarkerBlast(proj, blastWordSize, reportWordSize,
                                                    maxAlignmentsReported,
                                                    MarkerBlast.DEFAULT_REPORT_TO_TEMPORARY_FILE,
                                                    annotateGCContent, doBlast, numThreads,
                                                    probeFile, annotFile);
      if (markerPositionsOverride != null) {
        markerBlast.overrideMarkerPositions(markerPositionsOverride);
      }
      markerBlast.blastEm();
    } catch (Exception e) {
      e.printStackTrace();
    }

  }

}
