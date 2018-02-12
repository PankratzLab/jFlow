/**
 * 
 */
package org.genvisis.cnv.analysis.pod;

import java.awt.Color;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.genvisis.CLI;
import org.genvisis.cnv.analysis.pod.InformativeBAF.BAF_STRATEGY;
import org.genvisis.cnv.analysis.pod.InformativeBAF.InformativeResult;
import org.genvisis.cnv.analysis.pod.PODAnalysis.PODResults;
import org.genvisis.cnv.analysis.pod.PODAnalysis.PODResults.Builder;
import org.genvisis.cnv.analysis.pod.PODAnalysis.SEARCH_SPACE_TYPE;
import org.genvisis.cnv.analysis.pod.PODGenotype.Genotype;
import org.genvisis.cnv.analysis.pod.PODGenotype.POD;
import org.genvisis.cnv.filesys.MarkerDetailSet;
import org.genvisis.cnv.filesys.Pedigree;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.plots.ColorExt.ColorItem;
import org.genvisis.cnv.plots.ColorExt.ColorManager;
import org.genvisis.cnv.var.MosaicRegion;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.filesys.CNVariant;
import org.genvisis.filesys.LocusSet;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.BEDFileReader;

/**
 * Method to annotate segments of the genome for Parent of Origin
 */
public class PODAnnotator {

  enum TYPE {
    /**
     * .bed formatted file
     */
    BED,
    /**
     * {@link CNVariant} file, as loaded by {@link CNVariant#loadLocSet(String, Logger)}
     */
    CNVARIANT,
    /**
     * {@link MosaicRegion} file, as loaded by {@link MosaicRegion#loadMosLocSet(String, Logger)}
     */
    MOSAIC_REGION;
  }

  static class PODRegion<T extends Segment> extends Segment {

    /**
     * 
     */
    private static final long serialVersionUID = 1L;
    private T t;
    private PODResults podResults;

    PODRegion(T t, PODResults podResults) {
      super(t.getUCSClocation());
      this.t = t;
      this.podResults = podResults;
    }

    @Override
    public String toAnalysisString() {
      return t.toAnalysisString() + "\t" + podResults.toString();
    }

    @Override
    public String[] getHeader() {
      return ArrayUtils.concatAll(t.getHeader(), PODResults.getHeader());
    }
  }

  private static <T extends CNVariant> LocusSet<T> filter(LocusSet<T> set, int minNumMarkers,
                                                          Logger log) {
    List<T> filtered = new ArrayList<>();
    for (T t : set.getLoci()) {
      if (t.getNumMarkers() > minNumMarkers) {
        filtered.add(t);
      }
    }
    return filtered.isEmpty() ? null : new LocusSet<>(filtered, true, log);
  }

  private static <T extends Segment> List<PODRegion<T>> annotate(Project proj, String offDNA,
                                                                 String moDNA, String faDNA,
                                                                 LocusSet<T> set,
                                                                 MarkerDetailSet markerDetailSet,
                                                                 SEARCH_SPACE_TYPE sType) {

    List<PODResults> results = PODAnalysis.analyze(proj, markerDetailSet,
                                                   Arrays.asList(set.getLoci()), offDNA, moDNA,
                                                   faDNA, sType, InformativeBAF.CHEBYSHEV);

    List<PODRegion<T>> regionResults = new ArrayList<>();
    for (int i = 0; i < set.getLoci().length; i++) {
      if (!set.getLoci()[i].getUCSClocation().equals(results.get(i).getSearchSpaceString())) {
        throw new IllegalStateException("Mismatched segments returned: expected "
                                        + set.getLoci()[i].getUCSClocation() + " and got "
                                        + results.get(i).getSearchSpaceString());
      } else {
        regionResults.add(new PODRegion<T>(set.getLoci()[i], results.get(i)));
      }
    }
    return regionResults;
  }

  private static class DNATrio {

    private String offDNA;
    private String moDNA;
    private String faDNA;

    private DNATrio(String offDNA, String moDNA, String faDNA) {
      super();
      this.offDNA = offDNA;
      this.moDNA = moDNA;
      this.faDNA = faDNA;
    }

    private boolean isValid() {
      return offDNA != null && (moDNA != null || faDNA != null);
    }

    private boolean isFullTrio() {
      return offDNA != null && moDNA != null && faDNA != null;
    }
  }

  private static DNATrio getTrioFor(Project proj, String dna, Pedigree pedigree) {
    String offDNA = null;
    String moDNA = null;
    String faDNA = null;
    int offIndex = ext.indexOfStr(dna, pedigree.getDnas());
    if (offIndex >= 0) {
      offDNA = dna;
      int moIndex = pedigree.getMoDNAIndex(offIndex);
      int faIndex = pedigree.getFaDNAIndex(offIndex);
      if (moIndex >= 0) {
        moDNA = proj.getSamples()[moIndex];
      }
      if (faIndex >= 0) {
        faDNA = proj.getSamples()[faIndex];
      }
    }
    return new DNATrio(offDNA, moDNA, faDNA);
  }

  private static List<DNATrio> getDNATrios(Project proj, Pedigree pedigree) {
    List<DNATrio> trios = new ArrayList<>();
    int numFull = 0;
    for (String samp : proj.getSamples()) {
      DNATrio trio = getTrioFor(proj, samp, pedigree);
      if (trio.isValid()) {
        trios.add(trio);
        if (trio.isFullTrio()) {
          numFull++;
        }
      }
    }
    proj.getLog().reportTimeInfo("Found " + trios.size() + " parent offspring pairs, " + numFull
                                 + " full trios");
    return trios;
  }

  private static <T extends Segment> void annotateBedTrios(Project proj, List<DNATrio> trios,
                                                           LocusSet<T> set,
                                                           MarkerDetailSet markerDetailSet,
                                                           String output, SEARCH_SPACE_TYPE sType) {
    PrintWriter writer = Files.getAppropriateWriter(output);
    boolean headerWritten = false;
    for (DNATrio trio : trios) {
      List<PODRegion<T>> results = annotate(proj, trio.offDNA, trio.moDNA, trio.faDNA, set,
                                            markerDetailSet, sType);
      for (PODRegion<T> pRegion : results) {
        if (!headerWritten) {
          headerWritten = true;
          writer.println(ArrayUtils.toStr(pRegion.getHeader()));
        }
        writer.println(pRegion.toAnalysisString());
      }
    }
    writer.close();
  }

  private static <T extends CNVariant> void annotateCNVarTrios(Project proj, List<DNATrio> trios,
                                                               LocusSet<T> initialSet,
                                                               MarkerDetailSet markerDetailSet,
                                                               String output, TYPE type,

                                                               int minNumMarkers) {

    LocusSet<T> filteredSet = filter(initialSet, minNumMarkers, proj.getLog());
    if (filteredSet != null) {
      PrintWriter writer = Files.getAppropriateWriter(output);
      boolean headerWritten = false;
      proj.getLog()
          .reportTimeInfo("As a defense against duplicate FIDs/IIDs, DNA is assumed to be listed in both columns");
      Map<String, LocusSet<T>> map = CNVariant.breakIntoInds(filteredSet, proj.getLog());
      for (DNATrio trio : trios) {
        if (map.containsKey(trio.offDNA + "\t" + trio.offDNA)) {
          List<PODRegion<T>> results = annotate(proj, trio.offDNA, trio.moDNA, trio.faDNA,
                                                map.get(trio.offDNA + "\t" + trio.offDNA),
                                                markerDetailSet, SEARCH_SPACE_TYPE.INDIVIDUAL);
          for (PODRegion<T> pRegion : results) {
            if (!headerWritten) {
              headerWritten = true;
              writer.println(ArrayUtils.toStr(pRegion.getHeader()));
            }
            writer.println(pRegion.toAnalysisString());
          }
        } else {
          proj.getLog().reportTimeInfo("No " + type + " regions found for DNA " + trio.offDNA);
        }
      }
      writer.close();
    } else {
      proj.getLog()
          .reportTimeWarning("No regions remain  - applied numMarker filter " + minNumMarkers);
    }
  }

  private static void annotateFile(Project proj, Pedigree pedigree, String segFile, TYPE type,
                                   int minNumMarkers, SEARCH_SPACE_TYPE sType) {
    List<DNATrio> trios = getDNATrios(proj, pedigree);
    String output = ext.addToRoot(segFile, "." + type + ".POD.minMarkers" + minNumMarkers);
    proj.getLog().reportTimeInfo("Reporting to " + output);
    MarkerDetailSet markerDetailSet = proj.getMarkerSet();
    switch (type) {
      case BED:
        BEDFileReader reader = new BEDFileReader(segFile, false);
        annotateBedTrios(proj, trios, reader.loadAll(proj.getLog()), markerDetailSet, output,
                         sType);
        reader.close();
        break;
      case CNVARIANT:
        if (sType != SEARCH_SPACE_TYPE.INDIVIDUAL) {
          throw new IllegalArgumentException("Invalid search space type " + sType + " for " + type);
        }
        annotateCNVarTrios(proj, trios, CNVariant.loadLocSet(segFile, proj.getLog()),
                           markerDetailSet, output, type, minNumMarkers);
        break;
      case MOSAIC_REGION:
        if (sType != SEARCH_SPACE_TYPE.INDIVIDUAL) {
          throw new IllegalArgumentException("Invalid search space type " + sType + " for " + type);
        }
        annotateCNVarTrios(proj, trios, MosaicRegion.loadMosLocSet(segFile, proj.getLog()),
                           markerDetailSet, output, type, minNumMarkers);
        break;
      default:
        throw new IllegalArgumentException("Invalid type " + type);

    }
  }

  /**
   * Get a {@link POD} based {@link ColorManager} for this sample, if possible (requires full trio)
   * 
   * @param proj
   * @param offDNA
   * @param pedigree
   * @return
   */
  public static ColorManager<String> getPODColors(Project proj, String offDNA, String pedigree) {
    Map<String, Integer> map = proj.getMarkerIndices();
    String[] names = proj.getMarkerNames();
    int[] projectIndicesToPOD = new int[map.keySet().size()];
    for (int i = 0; i < names.length; i++) {
      projectIndicesToPOD[i] = i;
    }
    proj.getLog().reportTimeInfo("Deriving informative marker status");
    InformativeResult informativeResult = getPODStatus(proj, offDNA, new Pedigree(proj, pedigree),
                                                       projectIndicesToPOD);
    proj.getLog().reportTimeInfo("Assigning color scheme");

    Hashtable<String, ColorItem<String>> manager = new Hashtable<>();
    Hashtable<String, String> lookup = new Hashtable<>();

    for (POD pod : POD.values()) {
      switch (pod) {
        case MATERNAL:
          manager.put(pod.toString(), new ColorItem<String>(pod.toString(), Color.BLUE));
          break;
        case NONE:
          manager.put(pod.toString(), new ColorItem<String>(pod.toString(), Color.LIGHT_GRAY));
          break;
        case PATERNAL:
          manager.put(pod.toString(), new ColorItem<String>(pod.toString(), Color.MAGENTA));
          break;
        default:
          break;
      }
    }
    HashMap<String, Integer> infoMap = new HashMap<>();
    int index = 0;
    for (Integer i : informativeResult.getInformatives()) {
      infoMap.put(names[i], index);
      index++;
    }
    for (String marker : map.keySet()) {

      if (infoMap.containsKey(marker)) {
        lookup.put(marker, informativeResult.getPods().get(infoMap.get(marker)).toString());
      } else {
        lookup.put(marker, POD.NONE.toString());
      }
    }
    proj.getLog().reportTimeInfo("finished assigning color scheme");
    ColorManager<String> colorManager = new ColorManager<String>(lookup, manager) {

      /**
       *
       */
      private static final long serialVersionUID = 1L;
    };
    colorManager.setColorBAF(true);
    return colorManager;
  }

  private static InformativeResult getPODStatus(Project proj, String offDNA, Pedigree pedigree,
                                                int[] projectIndicesToPOD) {

    proj.getLog().reportTimeInfo("Sample: " + offDNA);
    DNATrio dnaTrio = getTrioFor(proj, offDNA, pedigree);
    if (dnaTrio.isFullTrio()) {
      Sample off = proj.getFullSampleFromRandomAccessFile(dnaTrio.offDNA);
      Sample mo = proj.getFullSampleFromRandomAccessFile(dnaTrio.moDNA);
      Sample fa = proj.getFullSampleFromRandomAccessFile(dnaTrio.faDNA);
      InformativeResult informativeResultTmp = InformativeBAF.getInformativeIndices(ArrayUtils.toDoubleArray(off.getBAFs()),
                                                                                    off.getAB_Genotypes(),
                                                                                    BAF_STRATEGY.HET_ONLY,
                                                                                    InformativeBAF.CHEBYSHEV);
      InformativeResult informativeResult = PODAnalysis.bafGenotypePOD(projectIndicesToPOD,
                                                                       new Builder(), off, mo, fa,
                                                                       informativeResultTmp);
      byte[] genos = off.getAB_Genotypes();
      for (Integer i : informativeResult.getInformatives()) {
        if (Genotype.fromByte(genos[i]) != Genotype.AB) {
          throw new IllegalArgumentException("Invalid baf outlier test for non het genotype");
        }
      }
      return informativeResult;
    } else {
      proj.getLog()
          .reportTimeWarning("Can not determine POD without valid parents and offspring in pedigree file");
      return new InformativeResult(new ArrayList<>(), new ArrayList<>(), new NormalDistribution());

    }
  }

  /**
   * @param args
   */
  public static void main(String[] args) {
    CLI c = new CLI(PODAnnotator.class);
    c.addArg(CLI.ARG_PROJ, CLI.DESC_PROJ);
    c.addArg("pedigreeFile", "pedigree file to use");
    c.addArg("inputFile", "region based input file (bed, cnv, mos)");
    c.addArgWithDefault("fileType",
                        "Type of input file using, options are "
                                    + ArrayUtils.toStr(TYPE.values(), " ,"),
                        TYPE.MOSAIC_REGION.toString());

    c.addArgWithDefault("minMarkers",
                        "For " + TYPE.CNVARIANT + " and " + TYPE.MOSAIC_REGION
                                      + " only variants with more than this number of markers will be analyzed",
                        "30");

    c.addArgWithDefault("searchType",
                        "If using a " + TYPE.BED + " file, how to search/collapse the segment "
                                      + ArrayUtils.toStr(SEARCH_SPACE_TYPE.values(), " ,"),
                        SEARCH_SPACE_TYPE.INDIVIDUAL.toString());

    c.parseWithExit(args);
    Project proj = new Project(c.get(CLI.ARG_PROJ));
    Pedigree pedigree = new Pedigree(proj, c.get("pedigreeFile"));

    annotateFile(proj, pedigree, c.get("inputFile"), TYPE.valueOf(c.get("fileType")),
                 c.getI("minMarkers"), SEARCH_SPACE_TYPE.valueOf(c.get("searchType")));

  }
}
