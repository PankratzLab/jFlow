package org.genvisis.cnv.filesys;

import java.io.PrintWriter;
import java.io.Serializable;
import java.lang.ref.Reference;
import java.lang.ref.SoftReference;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.NavigableMap;
import java.util.NavigableSet;
import java.util.Set;
import java.util.SortedSet;
import java.util.stream.Stream;
import org.genvisis.cnv.annotation.markers.AnnotationFileLoader.QUERY_TYPE;
import org.genvisis.cnv.annotation.markers.AnnotationParser;
import org.genvisis.cnv.annotation.markers.BlastAnnotationTypes.BLAST_ANNOTATION_TYPES;
import org.genvisis.cnv.annotation.markers.BlastAnnotationTypes.BlastAnnotation;
import org.genvisis.cnv.annotation.markers.MarkerAnnotationLoader;
import org.genvisis.cnv.annotation.markers.MarkerBlastAnnotation;
import org.genvisis.cnv.annotation.markers.MarkerSeqAnnotation;
import org.genvisis.cnv.manage.Resources.GENOME_BUILD;
import org.genvisis.cnv.manage.TextExport;
import org.genvisis.seq.manage.ReferenceGenome;
import org.genvisis.seq.manage.StrandOps;
import org.pankratzlab.common.AllelePair;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.GenomicPosition;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.SerializedFiles;
import org.pankratzlab.core.CLI;
import org.pankratzlab.shared.filesys.Segment;
import com.google.common.base.Joiner;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.ImmutableSortedMap;
import com.google.common.collect.ImmutableSortedSet;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Ordering;
import com.google.common.collect.Sets;
import com.google.common.primitives.Bytes;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Floats;
import com.google.common.primitives.Ints;
import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.Allele;

public class MarkerDetailSet implements MarkerSetInfo, Serializable, TextExport {

  public static class Marker implements Serializable, Comparable<Marker> {

    private static final long serialVersionUID = 4L;

    private final String name;
    private final GenomicPosition genomicPosition;
    private final AllelePair allelePair;

    /**
     * @param name marker name
     * @param genomicPosition marker's position in the genome
     * @param allelePair marker's pair of alleles
     */
    public Marker(String name, GenomicPosition genomicPosition, AllelePair allelePair) {
      super();
      this.name = name;
      this.genomicPosition = genomicPosition;
      this.allelePair = allelePair;
    }

    public Marker(String name, GenomicPosition genomicPosition) {
      this(name, genomicPosition, AllelePair.of('A', 'B'));
    }

    public String getName() {
      return name;
    }

    public byte getChr() {
      return genomicPosition.getChr();
    }

    public int getPosition() {
      return genomicPosition.getPosition();
    }

    public GenomicPosition getGenomicPosition() {
      return genomicPosition;
    }

    public char getA() {
      return getAB()[0];
    }

    public char getB() {
      return getAB()[1];
    }

    public char[] getAB() {
      return new char[] {allelePair.getA(), allelePair.getB()};
    }

    public Allele getAlleleA() {
      return allelePair.getAlleleA();
    }

    public Allele getAlleleB() {
      return allelePair.getAlleleB();
    }

    public AllelePair.RefAllele getRefAllele() {
      return allelePair.getRefAllele();
    }

    public Allele getRef() {
      return allelePair.getRef();
    }

    public Allele getAlt() {
      return allelePair.getAlt();
    }

    @Override
    public int hashCode() {
      final int prime = 31;
      int result = 1;
      result = prime * result + ((allelePair == null) ? 0 : allelePair.hashCode());
      result = prime * result + ((genomicPosition == null) ? 0 : genomicPosition.hashCode());
      result = prime * result + ((name == null) ? 0 : name.hashCode());
      return result;
    }

    @Override
    public boolean equals(Object obj) {
      if (this == obj) return true;
      if (obj == null) return false;
      if (!(obj instanceof Marker)) return false;
      Marker other = (Marker) obj;
      if (allelePair == null) {
        if (other.allelePair != null) return false;
      } else if (!allelePair.equals(other.allelePair)) return false;
      if (genomicPosition == null) {
        if (other.genomicPosition != null) return false;
      } else if (!genomicPosition.equals(other.genomicPosition)) return false;
      if (name == null) {
        if (other.name != null) return false;
      } else if (!name.equals(other.name)) return false;
      return true;
    }

    @Override
    public int compareTo(Marker o) {
      int cmp = genomicPosition.compareTo(o.genomicPosition);
      if (cmp != 0) return cmp;
      cmp = name.compareTo(o.name);
      if (cmp != 0) return cmp;
      cmp = allelePair.compareTo(o.allelePair);
      if (cmp != 0) return cmp;
      cmp = name.compareTo(o.name);
      return cmp;
    }
  }

  public static final long serialVersionUID = 10L;

  public static final List<String> MARKER_POSITIONS_ISSUES_HEADER = Lists.newArrayList("Marker",
                                                                                       "DefinedChr",
                                                                                       "DefinedPos",
                                                                                       "BestMatchChr",
                                                                                       "BestMatchPos");

  private final ImmutableSet<Marker> markersSet;
  private final ImmutableList<Marker> markersList;
  private final int hashCode;
  private final long markerSetFingerprint;

  private transient Reference<Map<String, Marker>> markerNameMapRef = null;
  private transient Reference<NavigableMap<Byte, NavigableSet<Marker>>> chrMapRef = null;
  private transient Reference<NavigableMap<GenomicPosition, NavigableSet<Marker>>> genomicPositionMapRef = null;
  private transient Reference<Map<Marker, Integer>> markerIndexMapRef = null;
  private transient Reference<String[]> markerNameArrayRef = null;

  private transient Reference<byte[]> chrArrayRef = null;
  private transient Reference<int[]> positionArrayRef = null;

  public MarkerDetailSet(MarkerSetInfo markerSet) {
    this(markerSet.getMarkerNames(), markerSet.getChrs(), markerSet.getPositions(),
         markerSet.getFingerprint());
  }

  public MarkerDetailSet(MarkerSetInfo markerSet, char[][] abAlleles) {
    this(markerSet.getMarkerNames(), markerSet.getChrs(), markerSet.getPositions(), abAlleles,
         markerSet.getFingerprint());
  }

  public MarkerDetailSet(String[] markerNames, byte[] chrs, int[] positions) {
    this(markerNames, chrs, positions, null);
  }

  public MarkerDetailSet(String[] markerNames, byte[] chrs, int[] positions,
                         long markerSetFingerprint) {
    this(markerNames, chrs, positions, null, markerSetFingerprint);
  }

  @SuppressWarnings("deprecation")
  public MarkerDetailSet(String[] markerNames, byte[] chrs, int[] positions, char[][] abAlleles) {
    this(markerNames, chrs, positions, abAlleles, MarkerSet.fingerprint(markerNames));
  }

  public MarkerDetailSet(String[] markerNames, byte[] chrs, int[] positions, char[][] abAlleles,
                         long markerSetFingerprint) {
    super();
    int numMarkers = markerNames.length;
    if (numMarkers != chrs.length || numMarkers != positions.length
        || (abAlleles != null && numMarkers != abAlleles.length)) {
      throw new IllegalArgumentException(this.getClass().getName()
                                         + " cannot be constructed with mismatched list of Markers and positions or AB Alleles");
    }
    ImmutableSet.Builder<Marker> markersBuilder = ImmutableSet.builder();
    for (int i = 0; i < markerNames.length; i++) {
      String name = markerNames[i];
      GenomicPosition genomicPosition = new GenomicPosition(chrs[i], positions[i]);
      Marker marker;
      if (abAlleles != null) {
        marker = new Marker(name, genomicPosition, AllelePair.of(abAlleles[i][0], abAlleles[i][1]));
      } else {
        marker = new Marker(name, genomicPosition);
      }
      markersBuilder.add(marker);
    }
    markersSet = markersBuilder.build();
    if (markersSet.size() != markerNames.length) throw new IllegalStateException("Duplicate Markers cannot exist in MarkerDetailSet");
    markersList = markersSet.asList();
    this.markerSetFingerprint = markerSetFingerprint;
    hashCode = generateHashCode();
  }

  @SuppressWarnings("deprecation")
  public MarkerDetailSet(Collection<Marker> markers) {
    this.markersSet = ImmutableSet.copyOf(markers);
    if (markersSet.size() != Iterables.size(markers)) throw new IllegalStateException("Duplicate Markers cannot exist in MarkerDetailSet");
    this.markersList = markersSet.asList();
    this.markerSetFingerprint = MarkerSet.fingerprint(getMarkerNames());
    this.hashCode = generateHashCode();
  }

  public static BlastAnnotation closestChrMatch(byte naiveChr, int naivePosition,
                                                Iterable<BlastAnnotation> matches) {
    // If more than one match, try choosing one with chr that matches naiveChr (chr
    // won't be affected by build) and closest position to naivePosition
    BlastAnnotation bestMatch = null;
    for (BlastAnnotation annotation : matches) {
      if (annotation.getRefLoc().getChr() == naiveChr) {
        if (bestMatch == null) {
          bestMatch = annotation;
        } else {
          if (Math.abs(naivePosition
                       - annotation.getRefLoc().getStart()) < Math.abs(
                                                                       naivePosition
                                                                       - bestMatch.getRefLoc()
                                                                                  .getStart())) {
            bestMatch = annotation;
          }
        }
      }
    }
    return bestMatch;
  }

  private int generateHashCode() {
    final int prime = 31;
    int result = 1;
    result = prime * result + ((markersList == null) ? 0 : markersList.hashCode());
    return result;
  }

  private Map<String, Marker> generateMarkerNameMap() {
    ImmutableMap.Builder<String, Marker> markerNameMapBuilder = ImmutableMap.builder();
    for (Marker marker : markersList) {
      markerNameMapBuilder.put(marker.getName(), marker);
    }
    return markerNameMapBuilder.build();
  }

  private NavigableMap<Byte, NavigableSet<Marker>> generateChrMap() {
    Map<Byte, ImmutableSortedSet.Builder<Marker>> chrBuilders = Maps.newHashMap();
    for (Marker marker : markersList) {
      Byte chr = marker.getChr();
      if (!chrBuilders.containsKey(chr)) {
        chrBuilders.put(chr, ImmutableSortedSet.naturalOrder());
      }
      chrBuilders.get(chr).add(marker);
    }
    ImmutableSortedMap.Builder<Byte, NavigableSet<Marker>> chrMapBuilder = ImmutableSortedMap.naturalOrder();
    chrBuilders.forEach((k, v) -> chrMapBuilder.put(k, v.build()));
    return chrMapBuilder.build();
  }

  private NavigableMap<GenomicPosition, NavigableSet<Marker>> generateGenomicPositionMap() {
    Map<GenomicPosition, ImmutableSortedSet.Builder<Marker>> gpBuilders = Maps.newHashMapWithExpectedSize(markersList.size());
    for (Marker marker : markersList) {
      GenomicPosition pos = marker.getGenomicPosition();
      if (!gpBuilders.containsKey(pos)) {
        gpBuilders.put(pos, ImmutableSortedSet.naturalOrder());
      }
      gpBuilders.get(pos).add(marker);
    }

    ImmutableSortedMap.Builder<GenomicPosition, NavigableSet<Marker>> positionMapBuilder = ImmutableSortedMap.naturalOrder();
    gpBuilders.forEach((k, v) -> positionMapBuilder.put(k, v.build()));
    return positionMapBuilder.build();
  }

  private Map<Marker, Integer> generateMarkerIndexMap() {
    ImmutableMap.Builder<Marker, Integer> indexMapBuilder = ImmutableMap.builder();
    for (int i = 0; i < markersList.size(); i++) {
      indexMapBuilder.put(markersList.get(i), i);
    }
    return indexMapBuilder.build();
  }

  private String[] generateMarkerNameArray() {
    String[] markerNames = new String[markersList.size()];
    for (int i = 0; i < markersList.size(); i++) {
      markerNames[i] = markersList.get(i).getName();
    }
    return markerNames;
  }

  private byte[] generateChrArray() {
    byte[] chrs = new byte[markersList.size()];
    for (int i = 0; i < markersList.size(); i++) {
      chrs[i] = markersList.get(i).getChr();
    }
    return chrs;
  }

  private int[] generatePositionsArray() {
    int[] chrs = new int[markersList.size()];
    for (int i = 0; i < markersList.size(); i++) {
      chrs[i] = markersList.get(i).getPosition();
    }
    return chrs;
  }

  public static class BlastParser {

    private static final Joiner TAB_JOINER = Joiner.on('\t');

    private final Project proj;
    private final MarkerSetInfo naiveMarkerSet;
    private final String blastAnnotation;
    private final ReferenceGenome referenceGenome;
    private final Logger log;

    private int missingPositionCount;
    private int ambiguousPositionCount;
    private List<String> missingSeqMkrs;

    /**
     * @param proj Project to parse MarkerDetailSet from
     * @param naiveMarkerSet a MarkerSetInfo with the chr and pos supplied by the array in Project
     *          order
     * @param blastAnnotation blast.vcf filename
     * @param log
     */
    public BlastParser(Project proj, MarkerSetInfo naiveMarkerSet, String blastAnnotation,
                       Logger log) {
      super();
      this.proj = proj;
      this.naiveMarkerSet = naiveMarkerSet;
      this.blastAnnotation = blastAnnotation;
      this.referenceGenome = proj.getReferenceGenome();
      this.log = log;
    }

    private static int determinePositionOffset(MarkerBlastAnnotation markerBlastAnnotation) {
      int interrogationPosition = markerBlastAnnotation.getMarkerSeqAnnotation()
                                                       .getInterrogationPosition();
      int positionOffset;
      String seq = markerBlastAnnotation.getMarkerSeqAnnotation().getSequence();
      if (seq != null) {
        int probeLength = markerBlastAnnotation.getMarkerSeqAnnotation().getSequence().length();
        positionOffset = interrogationPosition - probeLength + 1;
      } else {

        positionOffset = 0;
      }
      return positionOffset;
    }

    /**
     * Currently, this method will simply return the passed {@code naivePosition} but still goes
     * through the process of identifying a best match in order to annotate issues. When issue #238
     * is resolved, this method will be able to return the best match chr and position, allowing
     * projects to dynamically use chr and position that do not match the "project order" of markers
     * 
     * @param markerBlastAnnotation
     * @param naivePosition
     * @param issuesWriter
     * @return {@link GenomicPosition} parsed from annotation
     */
    private GenomicPosition parseGenomicPosition(MarkerBlastAnnotation markerBlastAnnotation,
                                                 GenomicPosition naivePosition,
                                                 PrintWriter issuesWriter) {
      boolean ambiguousPosition = false;
      List<BlastAnnotation> perfectMatches = markerBlastAnnotation.getAnnotationsFor(BLAST_ANNOTATION_TYPES.PERFECT_MATCH,
                                                                                     log);
      BlastAnnotation useThis = null;
      BlastAnnotation bestMatch = null;
      if (perfectMatches.size() == 1) {
        useThis = bestMatch = perfectMatches.get(0);
      } else if (!perfectMatches.isEmpty()) {
        // TODO implement liftOver here to check which match is actually in naiveMarkerSet
        ambiguousPosition = true;
        bestMatch = closestChrMatch(naivePosition.getChr(), naivePosition.getPosition(),
                                    perfectMatches);
        if (bestMatch == null) {
          log.reportTimeWarning("Arbitrarily selecting first of " + perfectMatches.size()
                                + " perfect matches for " + markerBlastAnnotation.getMarkerName()
                                + ", none of which are on the same chromosome as indicated by the array");
          bestMatch = perfectMatches.get(0);
        } else {
          useThis = bestMatch;
        }
      } else {
        // Otherwise, choose lowest e-value match from non-perfect matches
        List<BlastAnnotation> onTargetMatches = markerBlastAnnotation.getAnnotationsFor(BLAST_ANNOTATION_TYPES.ON_T_ALIGNMENTS_NON_PERFECT,
                                                                                        log);
        List<BlastAnnotation> offTargetMatches = markerBlastAnnotation.getAnnotationsFor(BLAST_ANNOTATION_TYPES.OFF_T_ALIGNMENTS,
                                                                                         log);
        Set<BlastAnnotation> nonPerfectMatches = Sets.newHashSet();
        nonPerfectMatches.addAll(onTargetMatches);
        nonPerfectMatches.addAll(offTargetMatches);
        List<BlastAnnotation> bestMatches = Lists.newArrayList();
        double bestMatchEVal = Double.MAX_VALUE;
        for (BlastAnnotation annotation : nonPerfectMatches) {
          double eVal = annotation.geteValue();
          if (eVal < bestMatchEVal) {
            bestMatches = Lists.newArrayList(annotation);
            bestMatchEVal = eVal;
          } else if (eVal == bestMatchEVal) {
            bestMatches.add(annotation);
          }
        }
        if (bestMatches.size() == 1) {
          bestMatch = bestMatches.get(0);
        } else if (!bestMatches.isEmpty()) {
          ambiguousPosition = true;
          bestMatch = closestChrMatch(naivePosition.getChr(), naivePosition.getPosition(),
                                      bestMatches);
          if (bestMatch == null) {
            log.reportTimeWarning("Arbitrarily selecting first of " + bestMatches.size()
                                  + " best matches with identical E-values for "
                                  + markerBlastAnnotation.getMarkerName()
                                  + ", none of which are on the same chromosome as indicated by the array");
            bestMatch = bestMatches.get(0);
          }
        }
      }
      byte bestMatchChr;
      int bestMatchPosition;
      if (bestMatch == null) {
        missingPositionCount++;
        bestMatchChr = 0;
        bestMatchPosition = 0;
      } else {
        Segment seg = bestMatch.getRefLoc();
        bestMatchChr = seg.getChr();
        int positionOffset = determinePositionOffset(markerBlastAnnotation);
        bestMatchPosition = bestMatch.getEffectiveInterrogationPosition(positionOffset, log);
      }
      if (ambiguousPosition) {
        ambiguousPositionCount++;
      }
      if (bestMatchPosition != naivePosition.getPosition()
          || bestMatchChr != naivePosition.getChr()) {
        issuesWriter.println(TAB_JOINER.join(markerBlastAnnotation.getMarkerName(),
                                             naivePosition.getChr(), naivePosition.getPosition(),
                                             bestMatchChr, bestMatchPosition));
      }
      return naivePosition;

    }

    private Allele determineReferenceAllele(GenomicPosition genomicPosition) {
      String[] refSequence = referenceGenome.getSequenceFor(new Segment(genomicPosition), true);
      if (refSequence == null) return null;
      return Allele.create(Iterables.getOnlyElement(Arrays.asList(refSequence)), true);
    }

    private AllelePair parseAllelePair(MarkerBlastAnnotation markerBlastAnnotation,
                                       GenomicPosition genomicPosition) {
      Allele a;
      Allele b;

      if (markerBlastAnnotation.getMarkerSeqAnnotation().getSequence() == null) {
        missingSeqMkrs.add(markerBlastAnnotation.getMarkerName());
        a = Allele.NO_CALL;
        b = Allele.NO_CALL;
      } else {
        MarkerSeqAnnotation markerSeqAnnotation = markerBlastAnnotation.getMarkerSeqAnnotation();
        Strand strand = markerSeqAnnotation.getStrand();
        a = StrandOps.flipIfNeeded(markerSeqAnnotation.getA(), strand);
        b = StrandOps.flipIfNeeded(markerSeqAnnotation.getB(), strand);
        Allele ref = determineReferenceAllele(genomicPosition);
        // ensure that only ref Allele is set as reference status
        if (ref != null && ref.basesMatch(a)) {
          a = ref;
          b = Allele.create(b, true);
        } else if (ref != null && ref.basesMatch(b)) {
          a = Allele.create(a, true);
          b = ref;
        } else {
          a = Allele.create(a, true);
          b = Allele.create(b, true);
        }
      }

      return AllelePair.of(a, b);
    }

    /**
     * @return a MarkerDetailSet generated using this BLAST annotation and the order of the
     *         naiveMarkerSet
     */
    public MarkerDetailSet parse() {

      final String[] markerNames = naiveMarkerSet.getMarkerNames();
      List<Marker> markers = Lists.newArrayListWithCapacity(markerNames.length);
      if (!Files.exists(blastAnnotation)) {
        log.reportTimeWarning("Could not find " + blastAnnotation
                              + ", cannot generate BLAST Annotation based Marker Set");
        return null;
      }
      Map<String, MarkerBlastAnnotation> masterMarkerList = MarkerBlastAnnotation.initForMarkers(markerNames);
      MarkerAnnotationLoader annotationLoader = new MarkerAnnotationLoader(null, blastAnnotation,
                                                                           new MarkerDetailSet(naiveMarkerSet),
                                                                           true, log);
      annotationLoader.setReportEvery(10000);
      List<Map<String, ? extends AnnotationParser>> parsers = Lists.newArrayList();
      parsers.add(masterMarkerList);
      annotationLoader.fillAnnotations(null, parsers, QUERY_TYPE.ONE_TO_ONE);

      missingPositionCount = 0;
      ambiguousPositionCount = 0;
      missingSeqMkrs = new ArrayList<>();
      try (PrintWriter issuesWriter = Files.getAppropriateWriter(proj.PROJECT_DIRECTORY.getValue()
                                                                 + "MarkerPositionBLASTIssues.txt")) {
        issuesWriter.println(TAB_JOINER.join(MARKER_POSITIONS_ISSUES_HEADER));
        for (int i = 0; i < markerNames.length; i++) {
          MarkerBlastAnnotation markerBlastAnnotation = masterMarkerList.get(markerNames[i]);
          GenomicPosition genomicPosition = parseGenomicPosition(markerBlastAnnotation,
                                                                 new GenomicPosition(naiveMarkerSet.getChrs()[i],
                                                                                     naiveMarkerSet.getPositions()[i]),
                                                                 issuesWriter);
          AllelePair allelePair = parseAllelePair(markerBlastAnnotation, genomicPosition);
          markers.add(new Marker(markerNames[i], genomicPosition, allelePair));
        }
      }
      if (ambiguousPositionCount > 0) {
        log.reportError("Warning - there " + (ambiguousPositionCount > 1 ? "were " : "was ")
                        + ambiguousPositionCount + " marker"
                        + (ambiguousPositionCount > 1 ? "s" : "")
                        + " with ambiguous BLAST matches");
      }
      if (missingPositionCount > 0) {
        log.reportError("Warning - there " + (missingPositionCount > 1 ? "were " : "was ")
                        + missingPositionCount + " marker" + (missingPositionCount > 1 ? "s" : "")
                        + " missing a BLAST result and association position");
      }
      if (missingSeqMkrs.size() > 0) {
        log.reportError("Warning - there " + (missingSeqMkrs.size() > 1 ? "were " : "was ")
                        + missingSeqMkrs.size() + " marker"
                        + (missingSeqMkrs.size() > 1 ? "s " : "")
                        + " missing a BLAST probe sequence.");
      }

      MarkerDetailSet markerDetailSet = new MarkerDetailSet(markers);

      if (markerDetailSet.getFingerprint() != naiveMarkerSet.getFingerprint()) {
        throw new IllegalStateException(markerDetailSet.getClass().getName() + " fingerprint ("
                                        + markerDetailSet.getFingerprint()
                                        + ") does not match naive "
                                        + naiveMarkerSet.getClass().getName() + " fingerprint ("
                                        + naiveMarkerSet.getFingerprint() + ")");
      }

      return markerDetailSet;

    }

  }

  public static MarkerDetailSet load(String filename) {
    return (MarkerDetailSet) SerializedFiles.readSerial(filename, false);
  }

  public void serialize(String filename) {
    SerializedFiles.writeSerial(this, filename);
  }

  /**
   * @return a List of {@link Marker}s in "project" order
   */
  public List<Marker> markersAsList() {
    return markersList;
  }

  /**
   * @return a Set of {@link Marker}s in "project" order
   */
  public Set<Marker> markersAsSet() {
    return markersSet;
  }

  public Map<String, Marker> getMarkerNameMap() {
    Map<String, Marker> markerNameMap = markerNameMapRef == null ? null : markerNameMapRef.get();
    if (markerNameMap == null) {
      markerNameMap = generateMarkerNameMap();
      markerNameMapRef = new SoftReference<>(markerNameMap);
    }
    return markerNameMap;
  }

  /**
   * @return an unmodifiable {@link NavigableMap} from chromosome to a {@link NavigableSet} of
   *         {@link Marker}s on that chromosome
   */
  public NavigableMap<Byte, NavigableSet<Marker>> getChrMap() {
    NavigableMap<Byte, NavigableSet<Marker>> chrMap = chrMapRef == null ? null : chrMapRef.get();
    if (chrMap == null) {
      chrMap = generateChrMap();
      chrMapRef = new SoftReference<>(chrMap);
    }
    return chrMap;
  }

  /**
   * @return an unmodifiable {@link NavigableMap} from {@link GenomicPosition} to a
   *         {@link NavigableSet} of {@link Marker}s
   */
  public NavigableMap<GenomicPosition, NavigableSet<Marker>> getGenomicPositionMap() {
    NavigableMap<GenomicPosition, NavigableSet<Marker>> genomicPositionMap = genomicPositionMapRef == null ? null
                                                                                                           : genomicPositionMapRef.get();
    if (genomicPositionMap == null) {
      genomicPositionMap = generateGenomicPositionMap();
      genomicPositionMapRef = new SoftReference<>(genomicPositionMap);
    }
    return genomicPositionMap;
  }

  public Iterable<Marker> getMarkersSortedByChrPos() {
    return Iterables.concat(getGenomicPositionMap().values());
  }

  public Map<Marker, Integer> getMarkerIndexMap() {
    Map<Marker, Integer> markerIndexMap = markerIndexMapRef == null ? null
                                                                    : markerIndexMapRef.get();
    if (markerIndexMap == null) {
      markerIndexMap = generateMarkerIndexMap();
      markerIndexMapRef = new SoftReference<>(markerIndexMap);
    }
    return markerIndexMap;
  }

  /**
   * @return an {@link Iterable} view of all autosomal markers
   */
  public Iterable<Marker> getAutosomalMarkers() {
    return Iterables.concat(getChrMap().subMap((byte) 1, (byte) 23).values());
  }

  /**
   * @param genomeBuild
   * @return a {@link Stream} over all markers on the X chromosome outside of the Psuedoautosomal
   *         Regions
   */
  public Stream<Marker> getNonPARXMarkers(GENOME_BUILD genomeBuild) {
    return viewMarkersInSeg((byte) 23, genomeBuild.getPAR1End() + 1,
                            genomeBuild.getPAR2Start() - 1);
  }

  @Override
  public void exportToText(Project proj, String filename) {
    PrintWriter writer;

    try {
      writer = Files.getAppropriateWriter(filename);
      for (Marker marker : markersList) {
        writer.println(Joiner.on('\t').join(marker.getName(), marker.getChr(), marker.getPosition(),
                                            marker.getA(), marker.getB()));
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing " + filename);
      e.printStackTrace();
    }
  }

  /**
   * Unless a new collection is explicitly needed, use {@link #viewMarkersInSeg(Segment)} for a view
   * of {@link Marker}s in seg
   * 
   * @param seg Segment
   * @return {@link ImmutableSortedSet} of Markers in seg
   */
  public ImmutableSortedSet<Marker> getMarkersInSeg(Segment seg) {
    return viewMarkersInSeg(seg).collect(ImmutableSortedSet.toImmutableSortedSet(Ordering.natural()));
  }

  /**
   * Convenience method for iterating over all {@link Marker}s in a particular {@link Segment}
   *
   * @see #viewMarkersInSeg(GenomicPosition, GenomicPosition)
   */
  public Stream<Marker> viewMarkersInSeg(Segment seg) {
    return viewMarkersInSeg(seg.getChr(), seg.getStart(), seg.getStop());
  }

  /**
   * @see #viewMarkersInSeg(GenomicPosition, GenomicPosition)
   */
  public Stream<Marker> viewMarkersInSeg(byte chr, int start, int stop) {
    return viewMarkersInSeg(new GenomicPosition(chr, start), new GenomicPosition(chr, stop));
  }

  /**
   * @param start First position of the desired region
   * @param stop Last position of the desired region
   * @return an Iterable that iterates over the {@link Marker}s in the specified region
   */
  public Stream<Marker> viewMarkersInSeg(GenomicPosition start, GenomicPosition stop) {
    NavigableMap<GenomicPosition, NavigableSet<Marker>> positionMap = getGenomicPositionMap();
    return positionMap.subMap(start, true, stop, true).values().stream()
                      .flatMap(Collection::stream);
  }

  /**
   * @param markerMask a boolean[] in project order indicating whether to include a given Marker in
   *          the returned {@link Set}
   * @return a mutable {@link Set} of all {@link Marker}s for which the project index in markerMask
   *         is true, the iteration order of the returned {@link Set} is in project order
   */
  public Set<Marker> includeProjectOrderMask(boolean[] markerMask) {
    if (markerMask.length != markersList.size()) throw new IllegalArgumentException("Project order marker mask does not match size of project marker List");
    Set<Marker> includedMarkers = Sets.newLinkedHashSet();
    for (int i = 0; i < markerMask.length; i++) {
      if (markerMask[i]) includedMarkers.add(markersList.get(i));
    }
    return includedMarkers;
  }

  /**
   * @param markerMask a boolean[] in project order indicating whether to exclude a given Marker in
   *          the returned {@link Set}
   * @return a mutable {@link Set} of all {@link Marker}s for which the project index in markerMask
   *         is false, the iteration order of the returned {@link Set} is in project order
   */
  public Set<Marker> excludeProjectOrderMask(boolean[] markerMask) {
    if (markerMask.length != markersList.size()) throw new IllegalArgumentException("Project order marker mask does not match size of project marker List");
    Set<Marker> excludedMarkers = Sets.newLinkedHashSet();
    for (int i = 0; i < markerMask.length; i++) {
      if (!markerMask[i]) excludedMarkers.add(markersList.get(i));
    }
    return excludedMarkers;
  }

  /**
   * @param dataInProjOrder {@link List} of data of the same length as {@link #markersAsList()}
   * @return a mutable {@link Map} from {@link Marker} to the data at the project index of the
   *         {@link Marker} in dataInProjOrder, the iteration order of the returned {@link Map} is
   *         in project order
   */
  public <T> Map<Marker, T> mapProjectOrderData(List<T> dataInProjOrder) {
    if (dataInProjOrder.size() != markersList.size()) throw new IllegalArgumentException("Project order data does not match size of project marker List");
    Map<Marker, T> markerDataMap = Maps.newLinkedHashMapWithExpectedSize(markersList.size());
    for (int i = 0; i < markersList.size(); i++) {
      markerDataMap.put(markersList.get(i), dataInProjOrder.get(i));
    }
    return markerDataMap;
  }

  /**
   * @see #mapProjectOrderData(List)
   */
  public <T> Map<Marker, T> mapProjectOrderData(T[] dataInProjOrder) {
    return mapProjectOrderData(Arrays.asList(dataInProjOrder));
  }

  /**
   * @see #mapProjectOrderData(List)
   */
  public Map<Marker, Double> mapProjectOrderData(double[] dataInProjOrder) {
    return mapProjectOrderData(Doubles.asList(dataInProjOrder));
  }

  /**
   * @see #mapProjectOrderData(List)
   */
  public Map<Marker, Float> mapProjectOrderData(float[] dataInProjOrder) {
    return mapProjectOrderData(Floats.asList(dataInProjOrder));
  }

  /**
   * @see #mapProjectOrderData(List)
   */
  public Map<Marker, Byte> mapProjectOrderData(byte[] dataInProjOrder) {
    return mapProjectOrderData(Bytes.asList(dataInProjOrder));
  }

  /**
   * @see #mapProjectOrderData(List)
   */
  public Map<Marker, Integer> mapProjectOrderData(int[] dataInProjOrder) {
    return mapProjectOrderData(Ints.asList(dataInProjOrder));
  }

  @Override
  public int hashCode() {
    return hashCode;
  }

  @Override
  public boolean equals(Object obj) {
    if (this == obj) return true;
    if (obj == null) return false;
    if (getClass() != obj.getClass()) return false;
    MarkerDetailSet other = (MarkerDetailSet) obj;
    if (markersList == null) {
      if (other.markersList != null) return false;
    } else if (!markersList.equals(other.markersList)) return false;
    return true;
  }

  @Override
  public long getFingerprint() {
    return markerSetFingerprint;
  }

  @SuppressWarnings("deprecation")
  @Override
  public boolean checkFingerprint(Sample samp) {
    return MarkerSet.checkFingerprints(this, samp);
  }

  @Override
  @Deprecated
  public int[] getIndicesOfMarkersIn(Segment seg, int[][] indicesByChr, Logger log) {
    Map<Marker, Integer> markerIndexMap = getMarkerIndexMap();
    Stream<Marker> segMarkers = viewMarkersInSeg(seg);
    return segMarkers.mapToInt(markerIndexMap::get).toArray();
  }

  @Override
  @Deprecated
  public String[] getMarkersIn(Segment seg, int[][] indicesByChr) {
    Stream<Marker> segMarkers = viewMarkersInSeg(seg);
    return segMarkers.map(Marker::getName).toArray(String[]::new);
  }

  @Override
  @Deprecated
  /**
   * @deprecated use the {@link MarkerDetailSet.Marker}-based methods instead
   * @return a String[] of marker names in project order WARNING: This method will potentially
   *         return the same array multiple times, which could have been modified by a previous
   *         caller, use {@link #clearArrayRefs()} to prevent this
   */
  public String[] getMarkerNames() {
    // FIXME Maintaining a reference to a mutable returned array is dangerous. Fix anywhere that is
    // relying on this method not building a fresh array every time it is called and
    // remove/reimplement
    String[] markerNames = markerNameArrayRef == null ? null : markerNameArrayRef.get();
    if (markerNames == null) {
      markerNames = generateMarkerNameArray();
      markerNameArrayRef = new SoftReference<>(markerNames);
    }
    return markerNames;
  }

  @Override
  @Deprecated
  public byte[] getChrs() {
    // FIXME Maintaining a reference to a mutable returned array is dangerous. Fix anywhere that is
    // relying on this method not building a fresh array every time it is called and
    // remove/reimplement
    byte[] chrs = chrArrayRef == null ? null : chrArrayRef.get();
    if (chrs == null) {
      chrs = generateChrArray();
      chrArrayRef = new SoftReference<>(chrs);
    }
    return chrs;
  }

  @Override
  @Deprecated
  public int[] getPositions() {
    // FIXME Maintaining a reference to a mutable returned array is dangerous. Fix anywhere that is
    // relying on this method not building a fresh array every time it is called and
    // remove/reimplement
    int[] positions = positionArrayRef == null ? null : positionArrayRef.get();
    if (positions == null) {
      positions = generatePositionsArray();
      positionArrayRef = new SoftReference<>(positions);
    }
    return positions;
  }

  /**
   * Clears the references to the arrays returned by {@link #getMarkerNames()}, {@link #getChrs()},
   * and {@link #getPositions()} These are inherently dangerous as the returned arrays could be
   * modified by their caller and then returned again later. This method can be used to cut-off that
   * risk
   */
  public void clearArrayRefs() {
    markerNameArrayRef = null;
    chrArrayRef = null;
    positionArrayRef = null;
  }

  @Override
  @Deprecated
  public int[][] getPositionsByChr() {
    NavigableMap<Byte, NavigableSet<Marker>> chrMap = getChrMap();
    int[][] positionsByChr = new int[MarkerSet.CHR_INDICES][];

    for (byte chr = 0; chr < MarkerSet.CHR_INDICES; chr++) {
      SortedSet<Marker> markers = chrMap.get(chr);
      if (markers == null) markers = ImmutableSortedSet.of();
      positionsByChr[chr] = new int[markers.size()];
      int i = 0;
      for (Marker marker : markers) {
        positionsByChr[chr][i] = marker.getPosition();
        i++;
      }
    }
    return positionsByChr;
  }

  @Override
  @Deprecated
  public int[][] getIndicesByChr() {
    NavigableMap<Byte, NavigableSet<Marker>> chrMap = getChrMap();
    Map<Marker, Integer> markerIndexMap = getMarkerIndexMap();
    int[][] indicesByChr = new int[MarkerSet.CHR_INDICES][0];

    for (byte chr : chrMap.keySet()) {
      if (chr >= 0 && chr < MarkerSet.CHR_INDICES) {
        SortedSet<Marker> markers = chrMap.get(chr);
        indicesByChr[chr] = new int[markers.size()];
        int i = 0;
        for (Marker marker : markers) {
          indicesByChr[chr][i] = markerIndexMap.get(marker);
          i++;
        }
      }
    }
    return indicesByChr;
  }

  @Deprecated
  public char[][] getABAlleles() {
    char[][] abAlleles = new char[markersList.size()][2];
    for (int i = 0; i < markersList.size(); i++) {
      abAlleles[i] = markersList.get(i).getAB();
    }
    return abAlleles;
  }

  @Deprecated
  public Map<String, Integer> getMarkerIndices() {
    Map<Marker, Integer> markerIndices = getMarkerIndexMap();
    Map<String, Integer> indices = Maps.newHashMapWithExpectedSize(markerIndices.size());
    for (Map.Entry<Marker, Integer> entry : markerIndices.entrySet()) {
      indices.put(entry.getKey().getName(), entry.getValue());
    }
    return indices;
  }

  public static void main(String... args) {
    CLI cli = new CLI(MarkerDetailSet.class);
    cli.addArg(CLI.ARG_PROJ, CLI.DESC_PROJ, CLI.EXAMPLE_PROJ, true);
    cli.parseWithExit(args);

    Project proj = new Project(cli.get(CLI.ARG_PROJ));
    proj.getMarkerSet();
  }
}
