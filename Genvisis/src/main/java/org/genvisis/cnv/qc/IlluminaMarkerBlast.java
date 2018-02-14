package org.genvisis.cnv.qc;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import org.genvisis.CLI;
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
import htsjdk.tribble.annotation.Strand;

public class IlluminaMarkerBlast extends MarkerBlast {

  public static final String ARG_MANIFEST = "fileSeq";
  public static final String DESC_MANIFEST = "full path to an Illumina manifest file";
  public static final String EXAMPLE_MANIFEST = "HumanExome-12-v1-0-B.csv";
  public static final int EXAMPLE_WORD_SIZE = MarkerBlast.getDefaultWordSize(ARRAY.ILLUMINA);

  private final String manifestFile;

  public IlluminaMarkerBlast(Project proj, int numThreads, String manifestFile) {
    this(proj, getDefaultWordSize(proj), getDefaultWordSize(proj), DEFAULT_MAX_ALIGNMENTS,
         DEFAULT_REPORT_TO_TEMPORARY_FILE, DEFAULT_ANNOTATE_GC_CONTENT, DEFAULT_DO_BLAST,
         numThreads, manifestFile);

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
   * @param manifestFile
   */
  public IlluminaMarkerBlast(Project proj, int blastWordSize, int reportWordSize,
                             int maxAlignmentsReported, boolean reportToTmp,
                             boolean annotateGCContent, boolean doBlast, int numThreads,
                             String manifestFile) {
    super(proj, blastWordSize, reportWordSize, maxAlignmentsReported, reportToTmp,
          annotateGCContent, doBlast, numThreads);
    this.manifestFile = manifestFile;
  }

  @Override
  protected String getSourceString() {
    return manifestFile;
  }

  @Override
  protected String getNameBase() {
    return ext.rootOf(manifestFile, true);
  }

  @Override
  protected MarkerFastaEntry[] getMarkerFastaEntries(BlastParams params, boolean alleleLookup) {
    ExtProjectDataParser.ProjectDataParserBuilder builder = formatParser();
    MarkerFastaEntry[] fastaEntries = null;
    try {
      int seqLength = ARRAY.ILLUMINA.getProbeLength();

      ExtProjectDataParser parser = builder.build(proj, manifestFile);
      parser.determineIndicesFromTitles();
      parser.loadData();
      ArrayList<MarkerFastaEntry> entries = new ArrayList<MarkerFastaEntry>(ArrayUtils.booleanArraySum(parser.getDataPresent()));
      MarkerSetInfo markerSet = proj.getMarkerSet();
      ReferenceGenome referenceGenome = Files.exists(proj.getReferenceGenomeFASTAFilename()) ? new ReferenceGenome(proj.getReferenceGenomeFASTAFilename(),
                                                                                                                   log)
                                                                                             : null;
      if (referenceGenome == null) {
        log.reportTimeWarning("A reference genome was not found");
      }
      String refStrandTitle = getRefStrand(parser);
      if (refStrandTitle == null && params != null) {
        params.setNotes("Warning, the positive and negative strands of the probe design are actually forward and reverse designations, due to parsing IlmnIDs ");
      }
      for (int i = 0; i < parser.getDataPresent().length; i++) {
        if (parser.getDataPresent()[i]) {
          if (alleleLookup) {
            if ((i + 1) % 10000 == 0 && referenceGenome != null) {
              log.reportTimeInfo("Loaded " + (i + 1) + " reference alleles from "
                                 + referenceGenome.getReferenceFasta());
            }
          }
          String seqA = null;
          String seqB = null;
          String markerName = parser.getDataToLoad()[i];
          Segment markerSegment = new Segment(markerSet.getChrs()[i], markerSet.getPositions()[i],
                                              markerSet.getPositions()[i]);
          seqA = parser.getStringDataForTitle("AlleleA_ProbeSeq")[i];
          seqB = parser.getStringDataForTitle("AlleleB_ProbeSeq")[i];

          String refStrand = null;
          TOP_BOT topBotProbe = TOP_BOT.valueOf(parser.getStringDataForTitle("IlmnStrand")[i].toUpperCase());
          TOP_BOT topBotRef = TOP_BOT.valueOf(parser.getStringDataForTitle("SourceStrand")[i].toUpperCase());
          if (refStrandTitle != null) {
            refStrand = parser.getStringDataForTitle(refStrandTitle)[i];
          } else {
            refStrand = parseStrandFromIlmnID(parser.getStringDataForTitle("IlmnID")[i].split("_"),
                                              topBotProbe);
          }
          Strand strand = null;
          if (refStrand.equals("+")) {
            strand = Strand.POSITIVE;
          } else if (refStrand.equals("-")) {
            strand = Strand.NEGATIVE;
          } else {
            log.reportError("Invalid RefStrand " + refStrand);
            return null;
          }
          if (seqA.length() != seqLength) {
            log.reportError("Sequence " + seqA + " did not have length "
                            + proj.ARRAY_TYPE.getValue().getProbeLength() + " " + markerName);
            return null;
          }

          String[] snp = parser.getStringDataForTitle("SNP")[i].replaceAll("\\[", "")
                                                               .replaceAll("\\]", "").split("/");
          AlleleParser alleleParser = new AlleleParser(markerName, markerSegment, snp[0], snp[1],
                                                       seqB, referenceGenome);
          if (alleleLookup) {
            alleleParser.parse(proj.getArrayType(), strand,
                               parser.getStringDataForTitle("SourceSeq")[i]);
          }
          if (seqB.equals("")) {// not ambiguous
            entries.add(new MarkerFastaEntry(markerName + PROBE_TAG.BOTH.getTag(), seqA, seqB,
                                             strand, seqLength, markerSegment, topBotProbe,
                                             topBotRef, alleleParser.getA(), alleleParser.getB(),
                                             alleleParser.getRef(), alleleParser.getAlts()));

          } else {// interrogationPosition is the last bp for Ambiguous SNPS
            entries.add(new MarkerFastaEntry(markerName + PROBE_TAG.A.getTag(), seqA, seqB, strand,
                                             seqLength - 1, markerSegment, topBotProbe, topBotRef,
                                             alleleParser.getA(), alleleParser.getB(),
                                             alleleParser.getRef(), alleleParser.getAlts()));
            entries.add(new MarkerFastaEntry(markerName + PROBE_TAG.B.getTag(), seqB, seqB, strand,
                                             seqLength - 1, markerSegment, topBotProbe, topBotRef,
                                             alleleParser.getA(), alleleParser.getB(),
                                             alleleParser.getRef(), alleleParser.getAlts()));
            if (seqB.length() != seqLength) {
              log.reportError("Sequence " + seqB + " did not have length "
                              + proj.ARRAY_TYPE.getValue().getProbeLength() + " " + markerName);
              return null;
            }
          }

        }
      }
      fastaEntries = entries.toArray(new MarkerFastaEntry[entries.size()]);
    } catch (FileNotFoundException e) {
      log.reportFileNotFound(manifestFile);
      e.printStackTrace();
    }
    log.reportTimeInfo("Found " + fastaEntries.length + " marker sequences");
    return fastaEntries;
  }

  private ProjectDataParserBuilder formatParser() {
    ExtProjectDataParser.ProjectDataParserBuilder builder = new ExtProjectDataParser.ProjectDataParserBuilder();
    if (proj.getArrayType() != ARRAY.ILLUMINA) {
      log.reportError("Array type was set to " + proj.getArrayType() + " and this file is for "
                      + ARRAY.ILLUMINA);
      return null;
    }
    builder.separator(",");
    builder.dataKeyColumnName("Name");
    String[] headerFlags = new String[] {"Name", "AlleleA_ProbeSeq"};
    builder.headerFlags(headerFlags);
    String ref = null;
    String[] header = Files.getLineContaining(manifestFile, ",", headerFlags, log);
    if (header != null) {
      if (ext.indexOfStr("RefStrand", header) >= 0) {
        ref = "RefStrand";
      } else if (ext.indexOfStr("GenomicStrand", header) >= 0) {
        ref = "GenomicStrand";
      }

    } else {
      log.reportError("Header of " + manifestFile + " not found");
      return null;
    }
    String[] dataTitles = new String[] {"AlleleA_ProbeSeq", "AlleleB_ProbeSeq", "SNP", "IlmnStrand",
                                        "SourceStrand", "SourceSeq", "IlmnID"};
    if (ref != null) {
      dataTitles = ArrayUtils.concatAll(dataTitles, new String[] {ref});
    } else {
      log.reportTimeWarning(manifestFile
                            + " did not have a column for RefStrand, will determine strand from IlmnID but these will be inacurate, forward assigned to +");
      // log.reportTimeError("Actully we are stopping, get the strands from
      // http://www.well.ox.ac.uk/~wrayner/strand/");
      // return null;
    }
    builder.stringDataTitles(dataTitles);

    builder.sampleBased(false);
    builder.treatAllNumeric(false);
    builder.requireAll(false);
    builder.verbose(false);
    return builder;
  }

  public static void main(String... args) {
    CLI c = new CLI(IlluminaMarkerBlast.class);

    c.addArg(CLI.ARG_PROJ, CLI.DESC_PROJ, CLI.EXAMPLE_PROJ, true, CLI.Arg.FILE);
    c.addArg(ARG_MANIFEST, DESC_MANIFEST, EXAMPLE_MANIFEST, true, CLI.Arg.FILE);
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
    String fileSeq = c.get(ARG_MANIFEST);
    int blastWordSize = c.getI(ARG_BLAST_WORD_SIZE);
    int reportWordSize = c.getI(ARG_REPORT_WORD_SIZE);
    int maxAlignmentsReported = c.getI(ARG_MAX_ALIGNMENTS);
    int numThreads = c.getI(CLI.ARG_THREADS);
    String markerPositionsOverride = c.get(ARG_MARKER_POSITIONS_OVERRIDE);
    boolean annotateGCContent = !c.has(FLAG_SKIP_GC_ANNOTATION);
    boolean doBlast = !c.has(FLAG_SKIP_BLAST);

    try {
      MarkerBlast markerBlast = new IlluminaMarkerBlast(proj, blastWordSize, reportWordSize,
                                                        maxAlignmentsReported,
                                                        MarkerBlast.DEFAULT_REPORT_TO_TEMPORARY_FILE,
                                                        annotateGCContent, doBlast, numThreads,
                                                        fileSeq);
      if (markerPositionsOverride != null) {
        markerBlast.overrideMarkerPositions(markerPositionsOverride);
      }
      markerBlast.blastEm();
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

}
