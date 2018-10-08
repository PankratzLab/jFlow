/**
 * 
 */
package org.genvisis.one.JL.cnv;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import org.genvisis.cnv.LocusSet;
import org.genvisis.cnv.LocusSet.TO_STRING_TYPE;
import org.genvisis.cnv.analysis.BeastScore;
import org.genvisis.cnv.annotation.segments.GDIAnnotator;
import org.genvisis.cnv.annotation.segments.GeneAnnotator;
import org.genvisis.cnv.annotation.segments.SegmentAnnotationKeys;
import org.genvisis.cnv.annotation.segments.SegmentAnotation;
import org.genvisis.cnv.annotation.segments.WESMappabilityAnnotator;
import org.genvisis.cnv.filesys.MarkerSet.PreparedMarkerSet;
import org.genvisis.cnv.filesys.CNVariant;
import org.genvisis.cnv.filesys.CNVariantAnnotated;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Project.ARRAY;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.seq.manage.BamImport.NGS_MARKER_TYPE;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.CLI;
import org.pankratzlab.common.CLI.Arg;
import org.pankratzlab.shared.filesys.Segment;
import com.google.common.primitives.Ints;

/**
 * Methods to annotate CNVariants with features of interest (Genes overlapped, mappability, etc)
 */
public class CNVAnnotation {

  /**
   * 
   */
  private static final String MAX_PROBLEMATIC_OVERLAP = "maxProblematicOverlap";
  /**
   * 
   */
  private static final String PROBLEM_REGIONS = "problemRegions";
  /**
   * 
   */
  private static final String FREQ_CNV = "freqCnvFiles";
  /**
   * 
   */
  private static final String CNV_FILE = "cnvs";

  /**
   * Annotate a given cnv file
   * 
   * @param proj {@link Project}
   * @param cnvFile file containing {@link CNVariant} to annotate
   * @param freqCnvFile use this file containing {@link CNVariant} to assign frequencies
   * @param problemRegionsFile bed file of problematic regions
   * @param outDir output directory
   * @param maxProblematicOverlap flag cnvs that have problematic regions overlapping more than this
   *          amount
   * @param threads number of threads
   */
  public static void annotate(Project proj, String cnvFile, String[] freqCnvFiles,
                              String problemRegionsFile, String outDir,
                              double maxProblematicOverlap, int threads) {
    new File(outDir).mkdirs();

    LocusSet<CNVariant> all = CNVariant.loadLocSet(cnvFile, proj.getLog());
    Map<String, LocusSet<CNVariant>> freqMap = new HashMap<>();
    for (String freqCnvFile : freqCnvFiles) {
      LocusSet<CNVariant> freqs = CNVariant.loadLocSet(freqCnvFile, proj.getLog());
      freqMap.put(freqCnvFile, freqs);
      proj.getLog().reportTimeInfo("Loaded " + freqs.getLoci().length
                                   + " cnvs for frequency annotation from " + freqCnvFile);
    }

    Map<String, LocusSet<CNVariant>> set = CNVariant.breakIntoInds(all, proj.getLog());
    proj.getLog().reportTimeInfo("Loaded " + all.getLoci().length + " cnvs across " + set.size()
                                 + " individuals to annotate");
    GeneAnnotator geneAnnotator = GeneAnnotator.getDefaultAnnotator(proj.GENOME_BUILD_VERSION.getValue(),
                                                                    proj.getLog());
    GDIAnnotator gdiAnnotator = GDIAnnotator.getDefaultGDIAnnotator(proj.getLog());
    WESMappabilityAnnotator wesMappabilityAnnotator = WESMappabilityAnnotator.getDefaultAnnotator(geneAnnotator,
                                                                                                  proj.GENOME_BUILD_VERSION.getValue(),
                                                                                                  proj.getLog());
    LocusSet<Segment> probs = new LocusSet<>(Segment.loadUCSCregions(problemRegionsFile, 0, false,
                                                                     proj.getLog()),
                                             true, proj.getLog());
    PreparedMarkerSet mds = PreparedMarkerSet.getPreparedMarkerSet(proj.getMarkerSet());
    ArrayList<CNVariantAnnotated> annotateds = new ArrayList<>();

    for (String sample : proj.getSamples()) {
      String key = proj.getSampleData(false).lookupFIDIID(sample);

      proj.getLog().reportTimeInfo("Annotating for sample " + sample);
      if (set.containsKey(key)) {
        CNVariant[] loci = set.get(key).getLoci();
        proj.getCNMarkersMask();
        int[][] cnvIndices = new int[loci.length][];
        Map<String, Integer> track = proj.getMarkerIndices();
        for (int i = 0; i < loci.length; i++) {
          String[] namesIn = mds.getMarkersIn(loci[i], mds.getIndicesByChr());
          ArrayList<Integer> nonVariant = new ArrayList<>();

          for (String element : namesIn) {
            NGS_MARKER_TYPE type = NGS_MARKER_TYPE.getType(element);
            if (proj.ARRAY_TYPE.getValue() != ARRAY.NGS || type != NGS_MARKER_TYPE.VARIANT_SITE) {
              nonVariant.add(track.get(element));
            }
          }
          cnvIndices[i] = Ints.toArray(nonVariant);
        }
        Sample samp = proj.getFullSampleFromRandomAccessFile(sample);
        BeastScore beastScore = null;
        if (samp != null) {
          beastScore = new BeastScore(samp.getLRRs(), mds.getIndicesByChr(), cnvIndices,
                                      proj.getLog());
          beastScore.computeBeastScores();
        }

        for (int i = 0; i < cnvIndices.length; i++) {
          CNVariant cnv = loci[i];
          SegmentAnotation segmentAnotation = new SegmentAnotation();
          if (beastScore != null) {
            segmentAnotation.getAttributes()
                            .put(SegmentAnnotationKeys.BEAST.toString(),
                                 new ArrayList<>(Arrays.asList(Float.toString(beastScore.getBeastHeights()[i]))));
          }

          wesMappabilityAnnotator.annotate(cnv, segmentAnotation);
          geneAnnotator.annotate(cnv, segmentAnotation);
          gdiAnnotator.annotate(cnv, segmentAnotation);
          Segment[] probSegs = probs.getOverLappingLoci(cnv);

          if (probSegs != null && probSegs.length > 0) {
            for (Segment problem : probSegs) {
              double bpOlap = (double) cnv.amountOfOverlapInBasepairs(problem) / cnv.getSize();
              if (bpOlap > maxProblematicOverlap) {
                segmentAnotation.getAttributes()
                                .put(SegmentAnnotationKeys.PROBLEM_REGION.toString(),
                                     new ArrayList<>(Arrays.asList("TRUE")));
                break;
              }
            }
          }
          CNVariantAnnotated cnVariantAnnotated = new CNVariantAnnotated(cnv, segmentAnotation);
          annotateds.add(cnVariantAnnotated);
        }
      }
    }
    wesMappabilityAnnotator.close();
    new LocusSet<>(annotateds, true, proj.getLog()).writeRegions(outDir + ext.rootOf(cnvFile)
                                                                 + "BEAST.annotated.cnvs",
                                                                 TO_STRING_TYPE.REGULAR, true,
                                                                 proj.getLog());

  }

  /**
   * @param args
   */
  public static void main(String[] args) {
    CLI c = new CLI(CNVAnnotation.class);
    c.addArg(CLI.ARG_PROJ, CLI.DESC_PROJ, true, Arg.FILE);
    c.addArg(CNV_FILE, "cnvFile to annotate", true, Arg.FILE);
    c.addArg(FREQ_CNV,
             "comma-delimited list cnvFiles with external/internal cnvs to use for frequency annotation, can be the same file as "
                       + CNV_FILE,
             true, Arg.FILE);
    c.addArg(PROBLEM_REGIONS, "list of problematic regions", true, Arg.FILE);
    c.addArg(CLI.ARG_OUTDIR, CLI.DESC_OUTDIR, true, Arg.FILE);

    c.addArgWithDefault(CLI.ARG_THREADS, CLI.DESC_THREADS, "3");
    c.addArgWithDefault(MAX_PROBLEMATIC_OVERLAP,
                        "if a cnv is covered by a problematic regions more than this, flag it",
                        "0.5");

    c.parseWithExit(args);
    annotate(new Project(c.get(CLI.ARG_PROJ)), c.get(CNV_FILE), c.get(FREQ_CNV).split(","),
             c.get(PROBLEM_REGIONS), c.get(CLI.ARG_OUTDIR), c.getD(MAX_PROBLEMATIC_OVERLAP),
             c.getI(CLI.ARG_THREADS));
  }

}
