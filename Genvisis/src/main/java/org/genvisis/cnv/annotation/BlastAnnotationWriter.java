package org.genvisis.cnv.annotation;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Hashtable;

import org.genvisis.cnv.annotation.BlastAnnotationTypes.BLAST_ANNOTATION_TYPES;
import org.genvisis.cnv.annotation.BlastAnnotationTypes.BlastAnnotation;
import org.genvisis.cnv.annotation.BlastAnnotationTypes.PROBE_TAG;
import org.genvisis.cnv.annotation.LocusAnnotation.Builder;
import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.genvisis.cnv.qc.MarkerBlast.MarkerFastaEntry;
import org.genvisis.common.Array;
import org.genvisis.common.ArraySpecialList.ArrayBlastAnnotationList;
import org.genvisis.common.Files;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.analysis.Blast;
import org.genvisis.seq.analysis.Blast.BlastResults;
import org.genvisis.seq.manage.CigarOps;
import org.genvisis.stats.Histogram.DynamicHistogram;

import htsjdk.samtools.Cigar;

/**
 * Handles
 */
public class BlastAnnotationWriter extends AnnotationFileWriter {

  /**
   * @param proj
   * @return initialized blast summaries for all markers
   */
  private static LocusAnnotation[] initializeSummaries(Project proj, int minAlignmentLength) {
    MarkerSet markerSet = proj.getMarkerSet();
    // TODO ablookup for annovar etc...
    byte[] chrs = markerSet.getChrs();
    int[] pos = markerSet.getPositions();
    String[] markerNames = proj.getMarkerNames();
    LocusAnnotation[] anDatas = new LocusAnnotation[markerNames.length];
    for (int i = 0; i < anDatas.length; i++) {
      Builder builder = new Builder();

      DynamicHistogram histogram =
          new DynamicHistogram(minAlignmentLength, proj.getArrayType().getProbeLength(), 0);
      MarkerBlastHistogramAnnotation markerBlastHistogramAnnotation =
          new MarkerBlastHistogramAnnotation(MarkerBlastHistogramAnnotation.DEFAULT_NAME,
                                             MarkerBlastHistogramAnnotation.DEFAULT_DESCRIPTION, histogram);
      MarkerEvalueHistogramAnnotation markerEvalueHistogramAnnotation =
          new MarkerEvalueHistogramAnnotation(MarkerEvalueHistogramAnnotation.DEFAULT_NAME,
                                              MarkerEvalueHistogramAnnotation.DEFAULT_DESCRIPTION);

      builder.annotations(Array.concatAll(BlastAnnotationTypes.getAnnotationDatas(),
                                          new AnnotationData[] {markerBlastHistogramAnnotation},
                                          new AnnotationData[] {markerEvalueHistogramAnnotation},
                                          new AnnotationData[] {MarkerSeqAnnotation.getDefault()}));
      Segment markerSeg = new Segment(chrs[i], pos[i], pos[i]);
      anDatas[i] = builder.build(markerNames[i], markerSeg);
    }
    return anDatas;
  }

  private static void summarize(Project proj, int maxAlignmentsReported, LocusAnnotation[] anDatas,
                                MarkerBlastHistogramAnnotation[] mHistogramAnnotations,
                                MarkerEvalueHistogramAnnotation[] mevHistogramAnnotations,
                                MarkerFastaEntry[] markerFastaEntries,
                                ArrayBlastAnnotationList[][] intLists) {
    Hashtable<String, Integer> markerSeqIndices = new Hashtable<String, Integer>();
    if (markerFastaEntries != null) {
      for (int i = 0; i < markerFastaEntries.length; i++) {
        PROBE_TAG tag = PROBE_TAG.parseMarkerTag(markerFastaEntries[i].getName(), proj.getLog());
        if (tag != PROBE_TAG.B) {
          String marker =
              markerFastaEntries[i].getName().substring(0, markerFastaEntries[i].getName().length()
                                                           - tag.getTag().length());
          markerSeqIndices.put(marker, i);
        }
      }
    }
    for (int i = 0; i < anDatas.length; i++) {
      if (i % 100000 == 0) {
        proj.getLog().reportTimeInfo("Summarized " + i + " markers ");
        proj.getLog().reportTimeInfo(" Free memory " + proj.getLog().memoryPercentFree());
      }
      for (int j = 0; j < BLAST_ANNOTATION_TYPES.values().length; j++) {
        BlastAnnotation[] tmp = intLists[i][j].toArray(new BlastAnnotation[intLists[i][j].size()]);
        if (tmp.length > 0) {
          Cigar[] cigarstmp = new Cigar[tmp.length];
          for (int k = 0; k < tmp.length; k++) {
            cigarstmp[k] = tmp[k].getCigar();
          }
          int[] order = CigarOps.sortByRefMatchLength(cigarstmp);

          BlastAnnotation[] toReport =
              new BlastAnnotation[Math.min(tmp.length, maxAlignmentsReported)];
          for (int k = 0; k < toReport.length; k++) {
            toReport[k] = tmp[order[k]];
          }

          anDatas[i].getAnnotations()[j].setData(Array.toStr(BlastAnnotation.toAnnotationString(toReport),
                                                             BLAST_ANNOTATION_TYPES.values()[j].getSep()));
        } else {
          anDatas[i].getAnnotations()[j].setData(BLAST_ANNOTATION_TYPES.values()[j].getDefaultValue());
        }
      }
      if (mHistogramAnnotations[i] != null) {
        mHistogramAnnotations[i].setDataToHistogram(true);
        anDatas[i].addAnnotation(mHistogramAnnotations[i]);
      }
      if (mevHistogramAnnotations[i] != null) {
        mevHistogramAnnotations[i].setDataToExactHistogram();
        anDatas[i].addAnnotation(mevHistogramAnnotations[i]);
      }

      if (markerSeqIndices.containsKey(anDatas[i].getLocusName())) {
        int seqIndex = markerSeqIndices.get(anDatas[i].getLocusName());
        MarkerFastaEntry mfe = markerFastaEntries[seqIndex];
        MarkerSeqAnnotation markerSeqAnnotation = new MarkerSeqAnnotation();
        markerSeqAnnotation.setDesignData(mfe.getSequence(), mfe.getSeqB(),
                                          mfe.getInterrogationPosition(), mfe.getStrand(),
                                          mfe.getTopBotProbe(), mfe.getTopBotRef(), mfe.getA(),
                                          mfe.getB());
        anDatas[i].addAnnotation(markerSeqAnnotation);
        anDatas[i].setRef(mfe.getRef());
        anDatas[i].setAlts(mfe.getAlts());
      }
    }
  }

  /**
   * @return the {@link LocusAnnotation } array for the markers in the project. This could also be
   *         passed to another {@link AnnotationFileWriter} instead of writing here
   */
  private static LocusAnnotation[] summarizeResultFile(Project proj, String[] blastResultFile,
                                                       MarkerFastaEntry[] markerFastaEntries,
                                                       int minAlignmentLength, int maxGaps,
                                                       int maxMismatches,
                                                       int maxAlignmentsReported) {
    int seqLength = proj.getArrayType().getProbeLength();
    LocusAnnotation[] anDatas = initializeSummaries(proj, minAlignmentLength);
    MarkerBlastHistogramAnnotation[] mHistogramAnnotations =
        new MarkerBlastHistogramAnnotation[proj.getMarkerNames().length];
    MarkerEvalueHistogramAnnotation[] mevHistogramAnnotations =
        new MarkerEvalueHistogramAnnotation[proj.getMarkerNames().length];
    ArrayBlastAnnotationList[][] intLists =
        new ArrayBlastAnnotationList[anDatas.length][BLAST_ANNOTATION_TYPES.values().length];
    for (int i = 0; i < intLists.length; i++) {
      for (int j = 0; j < intLists[i].length; j++) {
        intLists[i][j] = new ArrayBlastAnnotationList(maxAlignmentsReported);
      }
    }

    Hashtable<String, Integer> markerIndices = proj.getMarkerIndices();

    for (String element : blastResultFile) {
      proj.getLog().reportTimeInfo("Loading results from " + element);
      try {
        BufferedReader reader = Files.getAppropriateReader(element);
        int numEntries = 0;
        while (reader.ready()) {
          String[] line = reader.readLine().trim().split("\t");
          if ((line.length == Blast.BLAST_HEADER.length - 1
               || line.length == Blast.BLAST_HEADER.length)
              && !line[0].startsWith(Blast.BLAST_HEADER[0])) {
            numEntries++;
            if (numEntries % 1000000 == 0) {
              proj.getLog().reportTimeInfo("Processed " + numEntries + " blast results");
              proj.getLog().memoryPercentTotalFree();
              // reader.close();
              // break;
            }
            BlastResults blastResults = new BlastResults(line, proj.getLog());

            if (blastResults.getAlignmentLength() >= minAlignmentLength
                && blastResults.getGapOpens() <= maxGaps
                && blastResults.getMismatches() <= maxMismatches) {
              String marker = blastResults.getQueryID();
              PROBE_TAG tag = null;
              if (proj.getArrayType() == ARRAY.AFFY_GW6
                  || proj.getArrayType() == ARRAY.AFFY_GW6_CN) {
                if (marker.endsWith(PROBE_TAG.A.getTag())
                    || marker.endsWith(PROBE_TAG.B.getTag())) {
                  tag = PROBE_TAG.parseMarkerTag(marker, proj.getLog());
                  marker = marker.substring(0, marker.length() - tag.getTag().length());

                } else {
                  proj.getLog()
                      .reportTimeError("Query id did not end with " + PROBE_TAG.A.getTag() + "or "
                                       + PROBE_TAG.B.getTag()
                                       + " which is required for an AFFY array");
                }
              } else if (proj.getArrayType() == ARRAY.ILLUMINA) {
                try {
                  tag = PROBE_TAG.parseMarkerTag(marker, proj.getLog());
                  marker = marker.substring(0, marker.length() - tag.getTag().length());

                } catch (IllegalArgumentException ile) {
                  proj.getLog()
                      .reportTimeError("Query id (" + marker
                                       + ") did not end one of the following which is required for an ILLUMINA array");
                  for (int i = 0; i < PROBE_TAG.values().length; i++) {
                    proj.getLog().report(PROBE_TAG.values()[i].getTag());
                  }
                }
              }

              int markerIndex = markerIndices.get(marker);

              Segment markerSeg = anDatas[markerIndex].getSeg().getBufferedSegment(1);
              if (mHistogramAnnotations[markerIndex] == null) {
                DynamicHistogram histogram =
                    new DynamicHistogram(minAlignmentLength, proj.getArrayType().getProbeLength(),
                                         0);
                mHistogramAnnotations[markerIndex] =
                    new MarkerBlastHistogramAnnotation(MarkerBlastHistogramAnnotation.DEFAULT_NAME,
                                                       MarkerBlastHistogramAnnotation.DEFAULT_DESCRIPTION,
                                                       histogram);
              }
              if (mevHistogramAnnotations[markerIndex] == null) {
                mevHistogramAnnotations[markerIndex] =
                    new MarkerEvalueHistogramAnnotation(MarkerEvalueHistogramAnnotation.DEFAULT_NAME,
                                                        MarkerEvalueHistogramAnnotation.DEFAULT_DESCRIPTION);
              }
              mevHistogramAnnotations[markerIndex].addExactHistogram(blastResults.getEvalue());
              mHistogramAnnotations[markerIndex].getDynamicHistogram()
                                                .addDataPointToHistogram(blastResults.getAlignmentLength());
              for (int i = 0; i < BLAST_ANNOTATION_TYPES.values().length; i++) {
                if (BlastAnnotationTypes.shouldBeAnnotatedAs(proj, blastResults,
                                                             BLAST_ANNOTATION_TYPES.values()[i],
                                                             markerSeg, proj.getLog())) {
                  BlastAnnotation blastAnnotation =
                      new BlastAnnotation(CigarOps.convertBtopToCigar(blastResults, seqLength,
                                                                      proj.getLog()),
                                          blastResults.getSegment(), blastResults.determineStrand(),
                                          tag, blastResults.getEvalue());
                  if (intLists[markerIndex][i].size() >= maxAlignmentsReported) {// we now start
                                                                                 // bounding by
                                                                                 // e-value
                    double eval = blastResults.getEvalue();
                    if (!Double.isNaN(eval) && eval < intLists[markerIndex][i].getMaxEval()) {// otherwise
                                                                                              // we
                                                                                              // skip
                                                                                              // it
                      // System.out.println("Add\t" + eval + "\t" + blastAnnotation.getCigar());
                      // System.out.println("OlD\t" + intLists[markerIndex][i].getMaxEval() + "\t" +
                      // intLists[markerIndex][i].get(intLists[markerIndex][i].getMaxEvalIndex()).getCigar());
                      intLists[markerIndex][i].remove(intLists[markerIndex][i].getMaxEvalIndex());
                      intLists[markerIndex][i].add(blastAnnotation);
                      intLists[markerIndex][i].update();
                    }
                  } else {
                    intLists[markerIndex][i].add(blastAnnotation);
                  }

                }
              }
            }
          }
        }
        reader.close();
      } catch (FileNotFoundException fnfe) {
        proj.getLog().reportError("Error: file \"" + element + "\" not found in current directory");
        return null;
      } catch (IOException ioe) {
        proj.getLog().reportError("Error reading file \"" + element + "\"");
        return null;
      }
    }
    summarize(proj, maxAlignmentsReported, anDatas, mHistogramAnnotations, mevHistogramAnnotations,
              markerFastaEntries, intLists);
    return anDatas;
  }

  private final Project proj;
  private final String[] blastResultFiles;
  private MarkerFastaEntry[] markerFastaEntries;
  private final int minAlignmentLength;

  private final int maxGaps;

  private final int maxMismatches;

  private final int maxAlignmentsReported;

  public BlastAnnotationWriter(Project proj, AnalysisParams[] analysisParams, String outputFile,
                               String[] blastResultFiles, int minAlignmentLength, int maxGaps,
                               int maxMismatches, int maxAlignmentsReported) {
    super(proj, analysisParams,
          Array.concatAll(BlastAnnotationTypes.getBaseAnnotations(),
                          new Annotation[] {MarkerBlastHistogramAnnotation.getDefaultBlastAnnotation()},
                          new Annotation[] {MarkerSeqAnnotation.getDefault()},
                          new Annotation[] {MarkerEvalueHistogramAnnotation.getDefaultBlastAnnotation()}),
          outputFile, false);
    this.proj = proj;
    this.blastResultFiles = blastResultFiles;
    this.minAlignmentLength = minAlignmentLength;
    this.maxGaps = maxGaps;
    this.maxMismatches = maxMismatches;
    this.maxAlignmentsReported = maxAlignmentsReported;
  }

  public void setMarkerFastaEntries(MarkerFastaEntry[] markerFastaEntries) {
    this.markerFastaEntries = markerFastaEntries;
  }

  public void summarizeResultFiles(boolean skipDefaultValue) {
    proj.getLog().reportTimeInfo(" Free memory " + proj.getLog().memoryPercentFree());
    proj.getLog()
        .reportTimeInfo("Max number of alignments reported set to " + maxAlignmentsReported);
    LocusAnnotation[] annotations =
        summarizeResultFile(proj, blastResultFiles, markerFastaEntries, minAlignmentLength, maxGaps,
                            maxMismatches, maxAlignmentsReported);
    for (int i = 0; i < annotations.length; i++) {
      if ((i + 1) % 100000 == 0) {
        proj.getLog().reportTimeInfo("Written " + i + " markers ");
        proj.getLog().reportTimeInfo(" Free memory " + proj.getLog().memoryPercentFree());

      }

      write(annotations[i], skipDefaultValue, false);
    }
  }

}
