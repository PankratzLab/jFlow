package org.genvisis.cnv.qc;

import java.util.Collection;
import java.util.Map.Entry;
import java.util.Set;
import org.genvisis.cnv.filesys.ClusterFilterCollection;
import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.common.Logger;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.ImmutableSetMultimap;
import com.google.common.collect.Multimap;

public class DuplicateConcordance {

  private final ImmutableSetMultimap<Integer, Integer> discordantCalls;
  private final Logger log;

  public DuplicateConcordance(byte[] genotypes, Collection<Set<Integer>> duplicateIndices,
                              Logger log) {
    super();
    this.log = log;
    discordantCalls = calculateDiscordantCalls(genotypes, duplicateIndices);
  }

  /**
   * @param marker marker to check for concordance
   * @param duplicateIndices Sets of 2+ duplicate indices
   * @param log
   */
  public DuplicateConcordance(MarkerData marker, Collection<Set<Integer>> duplicateIndices,
                              Logger log) {
    super();
    this.log = log;
    byte[] genos = marker.getAbGenotypes();
    if (genos == null) {
      genos = marker.getForwardGenotypes();
      if (genos == null) {
        throw new IllegalArgumentException("No AB or forward genotypes found for marker: "
                                           + marker.getMarkerName());
      }
    }
    discordantCalls = calculateDiscordantCalls(genos, duplicateIndices);
  }

  /**
   * @param marker marker to check for concordance
   * @param duplicateIndices Sets of 2+ duplicate indices
   * @param clusterFilters clusterFilter to apply
   * @param log
   */
  public DuplicateConcordance(MarkerData marker, Collection<Set<Integer>> duplicateIndices,
                              ClusterFilterCollection clusterFilters, Logger log) {
    super();
    this.log = log;
    byte[] genos = marker.getAbGenotypesAfterFilters(clusterFilters, log);
    if (genos == null) {
      throw new IllegalArgumentException("No AB genotypes found for marker: "
                                         + marker.getMarkerName());
    }
    discordantCalls = calculateDiscordantCalls(genos, duplicateIndices);
  }

  /**
   * @return a {@link Multimap} from each discordant sample's index to the indices of the concordant
   *         duplicate(s). The genotype call with the most number of duplicates is considered
   *         concordant. When there is a simple pair of duplicates or equal genotype counts, the
   *         concordant genotype is arbitrary. Missing genotypes are not included. TODO:
   *         Intelligently pick the concordant genotype when equal (e.g. by highest frequency)
   */
  private static ImmutableSetMultimap<Integer, Integer> calculateDiscordantCalls(byte[] genotypes,
                                                                                 Collection<Set<Integer>> duplicateIndices) {
    ImmutableSetMultimap.Builder<Integer, Integer> discordantCalls = ImmutableSetMultimap.builder();
    for (Set<Integer> dupeSet : duplicateIndices) {
      HashMultimap<Byte, Integer> genoCallIndices = HashMultimap.create();
      for (int dupeIndex : dupeSet) {
        byte geno = genotypes[dupeIndex];
        if (geno != -1) {
          genoCallIndices.put(geno, dupeIndex);
        }
      }
      if (genoCallIndices.keySet().size() > 1) {
        int maxCount = -1;
        byte maxGeno = -1;
        for (Entry<Byte, Collection<Integer>> genoEntry : genoCallIndices.asMap().entrySet()) {
          int count = genoEntry.getValue().size();
          if (count > maxCount) {
            maxCount = count;
            maxGeno = genoEntry.getKey();
          }
        }
        Set<Integer> concordantIndices = genoCallIndices.removeAll(maxGeno);
        for (int discordantCallIndex : genoCallIndices.values()) {
          discordantCalls.putAll(discordantCallIndex, concordantIndices);
        }
      }
    }
    return discordantCalls.build();
  }

  public int calculateDiscordanceCount() {
    return discordantCalls.keySet().size();
  }

}
