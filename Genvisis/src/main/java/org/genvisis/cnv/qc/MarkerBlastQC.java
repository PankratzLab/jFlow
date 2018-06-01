package org.genvisis.cnv.qc;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.stream.Collectors;
import org.genvisis.cnv.annotation.markers.AnnotationFileLoader.QUERY_TYPE;
import org.genvisis.cnv.annotation.markers.AnnotationParser;
import org.genvisis.cnv.annotation.markers.BlastAnnotationTypes.BLAST_ANNOTATION_TYPES;
import org.genvisis.cnv.annotation.markers.BlastAnnotationTypes.BlastAnnotation;
import org.genvisis.cnv.annotation.markers.BlastAnnotationTypes.PROBE_TAG;
import org.genvisis.cnv.annotation.markers.MarkerAnnotationLoader;
import org.genvisis.cnv.annotation.markers.MarkerBlastAnnotation;
import org.genvisis.cnv.annotation.markers.MarkerGCAnnotation;
import org.genvisis.cnv.filesys.MarkerDetailSet;
import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.genvisis.cnv.filesys.MarkerSetInfo;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import com.google.common.collect.Lists;

public class MarkerBlastQC {

  private MarkerBlastQC() {
    throw new IllegalAccessError("Utility Class");
  }

  public static final double DEFAULT_CROSS_HYBE_THRESHOLD = 0.8;

  public static String defaultOneHitWondersFilename(String blastVCF) {
    return ext.rootRootOf(blastVCF) + "_OneHitWonders.txt";
  }

  public static void getOneHitWonders(Project proj, String blastVCF, String outFile,
                                      double crossHybePercent, Logger log) {
    if (outFile == null) {
      outFile = defaultOneHitWondersFilename(blastVCF);
    }
    List<String> oneHitters = getOneHitWonders(proj, blastVCF, crossHybePercent, log);
    log.report("Writing results to " + outFile);
    Files.writeIterable(oneHitters, outFile);
  }

  static class QCResults {

    public Set<String> getMatches() {
      return matches;
    }

    public Set<String> getMissing() {
      return missing;
    }

    public Set<String> getBad() {
      return bad;
    }

    public Set<String> getAmbig() {
      return ambig;
    }

    public Set<String> getPerfect() {
      return perfect;
    }

    private final Set<String> perfect = new HashSet<>();
    private final Set<String> matches = new HashSet<>();
    private final Set<String> missing = new HashSet<>();
    private final Set<String> bad = new HashSet<>();
    private final Set<String> ambig = new HashSet<>();

  }

  public static QCResults getSingleHitMarkers(Project proj, String blastVCF) {
    if (blastVCF == null) {
      blastVCF = proj.BLAST_ANNOTATION_FILENAME.getValue();
    }
    MarkerDetailSet markerSet = proj.getMarkerSet();
    Map<String, Marker> markerMap = markerSet.getMarkerNameMap();
    List<Marker> markers = markerSet.getMarkers();
    List<String> markerNames = markers.stream().map(Marker::getName).collect(Collectors.toList());
    MarkerAnnotationLoader markerAnnotationLoader = new MarkerAnnotationLoader(null,
                                                                               proj.BLAST_ANNOTATION_FILENAME.getValue(),
                                                                               proj.getMarkerSet(),
                                                                               true, proj.getLog());
    markerAnnotationLoader.setReportEvery(500000);

    Map<String, MarkerGCAnnotation> gcAnnotations = MarkerGCAnnotation.initForMarkers(proj,
                                                                                      markerNames);
    Map<String, MarkerBlastAnnotation> blastResults = MarkerBlastAnnotation.initForMarkers(markerNames);
    List<Map<String, ? extends AnnotationParser>> parsers = Lists.newArrayList();
    parsers.add(gcAnnotations);
    parsers.add(blastResults);

    markerAnnotationLoader.fillAnnotations(null, parsers, QUERY_TYPE.ONE_TO_ONE);

    QCResults result = new QCResults();
    for (int i = 0; i < markers.size(); i++) {
      String m = markers.get(i).getName();
      MarkerBlastAnnotation annot = blastResults.get(m);
      if (annot == null) {
        result.getMissing().add(m);
        continue;
      }
      List<BlastAnnotation> perfectMatches = annot.getAnnotationsFor(BLAST_ANNOTATION_TYPES.PERFECT_MATCH,
                                                                     proj.getLog());
      if (perfectMatches.size() == 1) {
        result.getPerfect().add(m);
      } else if (!perfectMatches.isEmpty()) {
        BlastAnnotation best = MarkerDetailSet.closestChrMatch(markerMap.get(m).getChr(),
                                                               markerMap.get(m).getPosition(),
                                                               perfectMatches);
        if (best != null) {
          result.getMatches().add(m);
        } else {
          result.getAmbig().add(m);
        }
      } else {
        result.getBad().add(m);
      }
    }

    return result;
  }

  public static List<String> getOneHitWonders(Project proj, String blastVCF,
                                              double crossHybePercent, Logger log) {
    if (blastVCF == null) {
      blastVCF = proj.BLAST_ANNOTATION_FILENAME.getValue();
    }
    MarkerSetInfo markerSet = proj.getMarkerSet();
    String[] markerNames = markerSet.getMarkerNames();
    MarkerAnnotationLoader markerAnnotationLoader = new MarkerAnnotationLoader(null,
                                                                               proj.BLAST_ANNOTATION_FILENAME.getValue(),
                                                                               proj.getMarkerSet(),
                                                                               true, log);
    markerAnnotationLoader.setReportEvery(500000);
    Map<String, MarkerGCAnnotation> gcAnnotations = MarkerGCAnnotation.initForMarkers(proj,
                                                                                      markerNames);
    Map<String, MarkerBlastAnnotation> blastResults = MarkerBlastAnnotation.initForMarkers(markerNames);
    List<Map<String, ? extends AnnotationParser>> parsers = Lists.newArrayList();
    parsers.add(gcAnnotations);
    parsers.add(blastResults);
    ArrayList<String> oneHitters = new ArrayList<>();

    // TODO: Update to more efficient Annotation Loading for entire file
    markerAnnotationLoader.fillAnnotations(null, parsers, QUERY_TYPE.ONE_TO_ONE);

    for (Entry<String, MarkerBlastAnnotation> blastResult : blastResults.entrySet()) {
      MarkerBlastAnnotation current = blastResult.getValue();
      String markerName = blastResult.getKey();
      List<BlastAnnotation> perfectMatches = current.getAnnotationsFor(BLAST_ANNOTATION_TYPES.PERFECT_MATCH,
                                                                       log);
      if (perfectMatches.size() == 1) {
        int[] alignmentHistogram = current.getAlignmentHistogram(proj);
        int sub = (int) Math.round(crossHybePercent * alignmentHistogram.length);
        int numHits = ArrayUtils.sum(ArrayUtils.subArray(alignmentHistogram, sub));
        if (numHits == 1) {
          // Perfect match is the only hit within crossHybePercent
          oneHitters.add(markerName);
        } else if (numHits == 2) {
          PROBE_TAG perfectTag = perfectMatches.get(0).getTag();
          if (perfectTag != PROBE_TAG.BOTH) {
            for (BlastAnnotation annotation : current.getAnnotationsFor(BLAST_ANNOTATION_TYPES.ON_T_ALIGNMENTS_NON_PERFECT,
                                                                        log)) {
              if (annotation.getTag() != perfectTag
                  && annotation.getRefLoc().getSize() == proj.getArrayType().getProbeLength() - 1) {
                // Second match is an on target read to the second probe, missing only the target
                // base
                oneHitters.add(markerName);
                break;
              }
            }
          }
        }
      }
    }

    log.reportTime(oneHitters.size() + " one hit wonder markers identified out of "
                   + markerNames.length + " total markers");
    return oneHitters;
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = null;
    String blastVCF = null;
    double crossHybridizationThreshold = DEFAULT_CROSS_HYBE_THRESHOLD;
    String usage = "\n" + "cnv.qc.MarkerBlast requires 3 arguments\n";
    usage += "   (1) Project file name (i.e. proj=" + filename + " (default))\n" + "";
    usage += "   (2) Blast.vcf filename  (i.e. blastVCF=" + "./blast.vcf.gz "
             + " (default based on project properties))\n" + "";
    usage += "   (3) Cross hybridization threshold  (i.e. crossHybridizationThreshold="
             + crossHybridizationThreshold + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("proj=")) {
        filename = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith("blastVCF=")) {
        blastVCF = ext.parseStringArg(arg, "");
        numArgs--;
      } else if (arg.startsWith("crossHybridizationThreshold=")) {
        crossHybridizationThreshold = ext.parseDoubleArg(arg);
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
      Project proj = null;

      if (filename == null) {
        System.err.println(usage);
        System.exit(1);
      } else {
        proj = new Project(filename);
      }
      if (proj == null) {
        System.err.println("Invalid project");
        System.exit(1);
      }
      if (blastVCF == null) {
        blastVCF = proj.BLAST_ANNOTATION_FILENAME.getValue();
      }
      getOneHitWonders(proj, blastVCF, defaultOneHitWondersFilename(blastVCF),
                       crossHybridizationThreshold, proj.getLog());
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
