package org.genvisis.cnv.qc;

import java.util.Collections;
import java.util.Map;
import java.util.Set;
import org.apache.commons.math3.stat.inference.ChiSquareTest;
import org.genvisis.cnv.filesys.ClusterFilterCollection;
import org.genvisis.cnv.filesys.MarkerData;
import org.pankratzlab.common.Logger;
import org.pankratzlab.shared.stats.ContingencyTable;
import org.pankratzlab.shared.stats.ProbDist;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.HashMultiset;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableSortedSet;
import com.google.common.collect.Multimaps;
import com.google.common.collect.Ordering;
import com.google.common.collect.Sets;

public class BatchEffects {

  public static class BatchEffect implements Comparable<BatchEffect> {

    private final double pValue;
    private final String batch;
    private final TestType testType;

    /**
     * @param pValue
     * @param batch
     * @param testType
     */
    private BatchEffect(double pValue, String batch, TestType testType) {
      super();
      this.pValue = pValue;
      this.batch = batch;
      this.testType = testType;
    }

    public double getpValue() {
      return pValue;
    }

    public String getBatch() {
      return batch;
    }

    public TestType getTestType() {
      return testType;
    }

    @Override
    public int hashCode() {
      final int prime = 31;
      int result = 1;
      result = prime * result + ((batch == null) ? 0 : batch.hashCode());
      long temp;
      temp = Double.doubleToLongBits(pValue);
      result = prime * result + (int) (temp ^ (temp >>> 32));
      result = prime * result + ((testType == null) ? 0 : testType.hashCode());
      return result;
    }

    @Override
    public boolean equals(Object obj) {
      if (this == obj) return true;
      if (obj == null) return false;
      if (getClass() != obj.getClass()) return false;
      BatchEffect other = (BatchEffect) obj;
      if (batch == null) {
        if (other.batch != null) return false;
      } else if (!batch.equals(other.batch)) return false;
      if (Double.doubleToLongBits(pValue) != Double.doubleToLongBits(other.pValue)) return false;
      if (testType != other.testType) return false;
      return true;
    }

    @Override
    public int compareTo(BatchEffect o) {
      int cmp = Double.compare(pValue, o.pValue);
      if (cmp != 0) return cmp;
      cmp = batch.compareTo(o.batch);
      if (cmp != 0) return cmp;
      return testType.compareTo(o.testType);
    }

  }

  public enum TestType {
    ALLELIC, MISSINGNESS;
  }

  private static final ChiSquareTest CHI_SQUARE = new ChiSquareTest();

  private final ImmutableMap<String, ImmutableSortedSet<BatchEffect>> notableEffects;

  public BatchEffects(byte[] genotypes, Map<Integer, String> indexBatches, double pValueThreshold,
                      Logger log) {
    super();
    notableEffects = calculateNotableEffects(genotypes, indexBatches, pValueThreshold);
  }

  public BatchEffects(MarkerData marker, ClusterFilterCollection clusterFilters,
                      Map<Integer, String> indexBatches, double pValueThreshold, Logger log) {
    super();
    byte[] genos = marker.getAbGenotypesAfterFilters(clusterFilters, log);
    if (genos == null) {
      throw new IllegalArgumentException("No AB genotypes found for marker: "
                                         + marker.getMarkerName());
    }
    notableEffects = calculateNotableEffects(genos, indexBatches, pValueThreshold);

  }

  public BatchEffects(MarkerData marker, Map<Integer, String> indexBatches, double pValueThreshold,
                      Logger log) {
    super();
    byte[] genos = marker.getAbGenotypes();
    if (genos == null) {
      genos = marker.getForwardGenotypes();
      if (genos == null) {
        throw new IllegalArgumentException("No AB or forward genotypes found for marker: "
                                           + marker.getMarkerName());
      }
    }
    notableEffects = calculateNotableEffects(marker.getAbGenotypes(), indexBatches,
                                             pValueThreshold);

  }

  /**
   * @return ImmutableMap from batch name to ImmutableSortedSet of {@link BatchEffect}s, with keys
   *         sorted by the minimum effect p-value
   */
  public ImmutableMap<String, ImmutableSortedSet<BatchEffect>> getNotableEffects() {
    return notableEffects;
  }

  public BatchEffect getMostNotableBatchEffect() {
    if (hasBatchEffects()) {
      return notableEffects.values().iterator().next().first();
    } else {
      return null;
    }
  }

  public boolean hasBatchEffects() {
    return !notableEffects.isEmpty();
  }

  private static ImmutableMap<String, ImmutableSortedSet<BatchEffect>> calculateNotableEffects(byte[] genos,
                                                                                               Map<Integer, String> indexBatches,
                                                                                               double pValueThreshold) {
    HashMultimap<String, Integer> batchIndices = HashMultimap.create();
    Multimaps.invertFrom(Multimaps.forMap(indexBatches), batchIndices);

    Set<String> batches = Sets.newHashSet();
    HashMultiset<String> missCounts = HashMultiset.create(batchIndices.keySet().size());
    HashMultiset<String> genoCounts = HashMultiset.create(batchIndices.keySet().size());
    HashMultiset<String> bAlleleCounts = HashMultiset.create(batchIndices.keySet().size());

    for (int i = 0; i < genos.length; i++) {
      String batch = indexBatches.get(i);
      if (batch != null) {
        batches.add(batch);
        byte geno = genos[i];
        if (0 <= geno && geno <= 2) {
          bAlleleCounts.add(batch, geno);
          genoCounts.add(batch);
        } else if (geno == -1) {
          missCounts.add(batch);
        } else {
          throw new IllegalArgumentException("Genotypes should be -1, 0, 1, or 2, found: " + geno);
        }
      }
    }
    ImmutableMap.Builder<String, ImmutableSortedSet<BatchEffect>> batchEffectsBuilder = ImmutableMap.builder();
    batchEffectsBuilder.orderEntriesByValue(Ordering.natural().onResultOf(Collections::min));
    int totalMisses = missCounts.size();
    int totalGenos = genoCounts.size();

    int totalBAlleles = bAlleleCounts.size();
    int totalAlleles = totalGenos * 2;

    for (String batch : batches) {
      int batchMisses = missCounts.count(batch);
      int batchGenos = genoCounts.count(batch);

      ImmutableSortedSet.Builder<BatchEffect> batchNotableEffectsBuilder = ImmutableSortedSet.naturalOrder();

      double missP = calcMissChiSquareP(batchMisses, batchGenos, totalMisses, totalGenos);
      if (missP <= pValueThreshold) {
        batchNotableEffectsBuilder.add(new BatchEffect(missP, batch, TestType.MISSINGNESS));
      }

      int batchBAlleles = bAlleleCounts.count(batch);
      int batchAlleles = batchGenos * 2;

      double allelicP = calcAllelicChiSquareP(batchBAlleles, batchAlleles, totalBAlleles,
                                              totalAlleles);
      if (allelicP <= pValueThreshold) {
        batchNotableEffectsBuilder.add(new BatchEffect(allelicP, batch, TestType.ALLELIC));
      }

      ImmutableSortedSet<BatchEffect> batchNotableEffects = batchNotableEffectsBuilder.build();
      if (!batchNotableEffects.isEmpty()) {
        batchEffectsBuilder.put(batch, batchNotableEffects);
      }

    }

    return batchEffectsBuilder.build();

  }

  private static double calcMissChiSquareP(int batchMisses, int batchGenos, int totalMisses,
                                           int totalGenos) {
    int otherMisses = totalMisses - batchMisses;
    int otherGenos = totalGenos - batchGenos;
    double[][] counts = new double[][] {{batchMisses, batchGenos}, {otherMisses, otherGenos}};
    return ProbDist.ChiDist(ContingencyTable.ChiSquare(counts, false, false), 1);
  }

  private static double calcAllelicChiSquareP(int batchBAlleles, int batchAlleles,
                                              int totalBAlleles, int totalAlleles) {
    int batchAAlleles = batchAlleles - batchBAlleles;
    int otherBAlleles = totalBAlleles - batchBAlleles;
    int otherAAlelles = totalAlleles - (batchAlleles + otherBAlleles);
    double[][] counts = new double[][] {{batchBAlleles, batchAAlleles},
                                        {otherBAlleles, otherAAlelles}};
    return ProbDist.ChiDist(ContingencyTable.ChiSquare(counts, false, false), 1);
  }

}
