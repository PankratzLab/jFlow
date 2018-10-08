package org.genvisis.cnv.imputation;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.apache.poi.util.IOUtils;
import org.genvisis.cnv.filesys.MarkerDetailSet;
import org.genvisis.cnv.filesys.MarkerDetailSet.Marker;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.seq.manage.StrandOps;
import org.pankratzlab.common.Aliases;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.GenomicPosition;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ProgressMonitor;
import org.pankratzlab.common.ext;
import org.pankratzlab.shared.filesys.Positions;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.Allele;

public class ImputationPrep {

  private static final String[][] REF_COLS;

  static {
    String[] markerNameAliases = Arrays.copyOf(Aliases.MARKER_NAMES,
                                               Aliases.MARKER_NAMES.length + 1);
    markerNameAliases[markerNameAliases.length - 1] = "ID";

    REF_COLS = new String[][] {Aliases.CHRS, Aliases.POSITIONS, markerNameAliases,
                               Aliases.REF_ALLELES, Aliases.ALT_ALLELES};
  }
  private static final Map<Character, Character> PALINDROMIC_PAIRS = ImmutableMap.of('A', 'T', 'T',
                                                                                     'A', 'G', 'C',
                                                                                     'C', 'G');

  private final Project proj;
  private final String referenceFile;
  private final Logger log;

  private final Set<Marker> matchingMarkers;

  static class ReferencePosition {

    private byte chr;
    private int position;
    private String id;
    private Allele ref;
    private Allele alt;

    /**
     * @param chr
     * @param position
     * @param id
     * @param ref
     * @param alt
     */
    public ReferencePosition(byte chr, int position, String id, Allele ref, Allele alt) {
      super();
      this.chr = chr;
      this.position = position;
      this.id = id;
      this.ref = ref;
      this.alt = alt;
    }

    public byte getChr() {
      return chr;
    }

    public int getPosition() {
      return position;
    }

    public String getId() {
      return id;
    }

    public Allele getRef() {
      return ref;
    }

    public Allele getAlt() {
      return alt;
    }

  }

  /**
   * @param proj Project to run for. Must have a {@link Project#BLAST_ANNOTATION_FILENAME} or an
   *          {@link MarkerDetailSet} that otherwise includes Ref/Alt alleles
   * @param referenceFile The reference panel / site list file (preferably gzipped). Must include
   *          chromosome, position, ref allele, and alt allele
   */
  public ImputationPrep(Project proj, String referenceFile) {
    super();
    this.proj = proj;
    this.referenceFile = referenceFile;
    this.log = proj.getLog();
    matchingMarkers = filterMarkers();
  }

  public Set<Marker> getMatchingMarkers() {
    return matchingMarkers;
  }

  private Map<Byte, Map<Integer, Set<ReferencePosition>>> readRefFile() {
    Set<GenomicPosition> markerSetPositions = proj.getMarkerSet().getGenomicPositionMap().keySet();
    log.report("Parsing reference panel file");
    long time = System.currentTimeMillis();
    Map<Byte, Map<Integer, Set<ReferencePosition>>> referencePositionsBuild = Maps.newHashMap();
    String taskName = "parseRefFile";
    proj.getProgressMonitor().beginDeterminateTask(taskName, "Parsing reference panel file",
                                                   Files.countLines(referenceFile, 0),
                                                   ProgressMonitor.DISPLAY_MODE.GUI_AND_CONSOLE);
    BufferedReader reader = null;
    try {
      reader = Files.getAppropriateReader(referenceFile);
      String header = null;
      do {
        header = reader.readLine();
        if (header == null) {
          log.reportError("Reference file is empty");
          throw new IllegalArgumentException("Reference file is empty");
        }
      } while (header.startsWith("##"));
      String delim = ext.determineDelimiter(header);
      int[] cols = ext.indexFactors(REF_COLS, header.split(delim), false, true, false);
      for (int index : cols) {
        if (index < 0) {
          throw new IllegalArgumentException("Invalid Reference File header");
        }
      }
      String temp;
      while ((temp = reader.readLine()) != null) {
        String[] refLine = ArrayUtils.subArray(temp.split(delim), cols);
        proj.getProgressMonitor().updateTask(taskName);
        byte chr = Positions.chromosomeNumber(refLine[0], log);
        int position;
        try {
          position = Integer.parseInt(refLine[1]);
        } catch (NumberFormatException e) {
          throw new IllegalStateException("Imputation reference file (" + referenceFile
                                          + ") contains a non-integer position: " + refLine[1]);
        }
        if (markerSetPositions.contains(new GenomicPosition(chr, position))) {
          String id = refLine[2];
          Allele ref = Allele.create(refLine[3], true);
          Allele alt = Allele.create(refLine[4], false);
          Map<Integer, Set<ReferencePosition>> posMap = referencePositionsBuild.get(chr);
          if (posMap == null) {
            posMap = Maps.newHashMap();
            referencePositionsBuild.put(chr, posMap);
          }
          Set<ReferencePosition> refPosSet = posMap.get(position);
          if (refPosSet == null) {
            refPosSet = Sets.newHashSet();
            posMap.put(position, refPosSet);
          }
          refPosSet.add(new ReferencePosition(chr, position, id, ref, alt));
        }
      }
    } catch (IOException ioe) {
      log.reportIOException(referenceFile);
    } finally {
      IOUtils.closeQuietly(reader);
    }
    proj.getProgressMonitor().endTask(taskName);
    log.report("Finished parsing Reference File in " + ext.getTimeElapsed(time));
    return referencePositionsBuild;
  }

  private Set<Marker> filterMarkers() {
    MarkerDetailSet markerSet = proj.getMarkerSet();
    List<Marker> markers = markerSet.markersAsList();
    Map<Byte, Map<Integer, Set<ReferencePosition>>> referencePositions = readRefFile();
    ImmutableSet.Builder<Marker> matchingMarkersBuilder = ImmutableSet.builder();
    int mismatchPos = 0;
    int alleleFlips = 0;
    int strandFlips = 0;
    int strandAlelleFlips = 0;
    int mismatchAlleles = 0;

    String taskName = "Matching Project Markers to Reference Panel";
    proj.getProgressMonitor().beginDeterminateTask(taskName, taskName, markers.size(),
                                                   ProgressMonitor.DISPLAY_MODE.GUI_AND_CONSOLE);
    List<Marker> noMatches = Lists.newArrayList();
    for (Marker marker : markers) {
      proj.getProgressMonitor().updateTask(taskName);
      Map<Integer, Set<ReferencePosition>> chrMap = referencePositions.get(marker.getChr());
      if (chrMap == null) {
        log.reportError("Warning - Chr " + marker.getChr() + " missing from reference file");
        referencePositions.put(marker.getChr(), new HashMap<Integer, Set<ReferencePosition>>());
        chrMap = referencePositions.get(marker.getChr());
      }
      Set<ReferencePosition> refMatches = chrMap.get(marker.getPosition());
      if (refMatches == null) {
        mismatchPos++;
        noMatches.add(marker);
      } else {
        // Copy the matches in order to pick off bad matches without removing from actual structure
        refMatches = Sets.newHashSet(refMatches);
        Allele ref = marker.getRef();
        Allele alt = marker.getAlt();
        // Find the best matched ReferencePosition from all of the matches
        if (matches(ref, alt, refMatches)) {
          matchingMarkersBuilder.add(marker);
          // Don't accept anything less than a proper match but count for reporting purposes
        } else if (matches(alt, ref, refMatches)) {
          alleleFlips++;
        } else if (alt != null && ref != null) {
          // TODO Try reversing as well (for Indels on wrong Strand)
          Allele refFlip = Allele.create(StrandOps.flipsIfNeeded(ref.getBaseString(),
                                                                 Strand.NEGATIVE, true));
          Allele altFlip = Allele.create(StrandOps.flipsIfNeeded(alt.getBaseString(),
                                                                 Strand.NEGATIVE, true));
          if (matches(refFlip, altFlip, refMatches)) {
            strandFlips++;
          } else if (matches(altFlip, refFlip, refMatches)) {
            strandAlelleFlips++;
          } else {
            mismatchAlleles++;
          }
        } else {
          mismatchAlleles++;
        }
      }
    }

    proj.getProgressMonitor().endTask(taskName);

    Set<Marker> matchingMarkersSet = matchingMarkersBuilder.build();

    log.report(mismatchPos + " positions did not match to a position in the reference set");
    log.report(mismatchAlleles + " positions had mismatched alleles and were excluded");
    log.report(alleleFlips + " positions would have matched after flipping the alleles");
    log.report(strandFlips + " positions would have matched after flipping the strand");
    log.report(strandAlelleFlips
               + " positions would have matched after flipping the alleles and strand");
    log.report("A total of " + (markers.size() - matchingMarkersSet.size())
               + " positions were excluded");
    log.report("");
    log.report(matchingMarkersSet.size()
               + " positions matched the reference set and will be used for imputation");
    return matchingMarkersSet;
  }

  private static boolean matches(Allele a, Allele b, Collection<ReferencePosition> refMatches) {
    if (a != null && b != null) {
      for (ReferencePosition refPos : refMatches) {
        if (a.basesMatch(refPos.getRef()) && b.basesMatch(refPos.getAlt())) {
          return true;
        }
      }
    }
    return false;
  }

  public static void validateRefFile(String referenceFile, Logger log) throws IOException,
                                                                       IllegalArgumentException {
    BufferedReader reader;
    reader = Files.getAppropriateReader(referenceFile);
    String header = reader.readLine();
    if (header == null) {
      log.reportError("Reference file is empty");
      throw new IllegalArgumentException("Reference file is empty");
    }
    String delim = ext.determineDelimiter(header);
    int[] cols = ext.indexFactors(REF_COLS, header.split(delim), false, true, false);
    for (int index : cols) {
      if (index < 0) {
        throw new IllegalArgumentException("Invalid Reference File header");
      }
    }
    reader.close();
  }

}
