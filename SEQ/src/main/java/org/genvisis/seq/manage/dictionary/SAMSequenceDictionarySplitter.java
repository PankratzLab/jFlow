/**
 * 
 */
package org.genvisis.seq.manage.dictionary;

import java.util.ArrayList;
import java.util.List;
import org.genvisis.cnv.filesys.Positions;
import org.pankratzlab.common.Logger;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

/**
 * Utility to split a {@link SAMSequenceDictionary} into approximately equal base pair bins
 */
public class SAMSequenceDictionarySplitter {

  private SAMSequenceDictionarySplitter() {

  }

  /**
   * Stores the chromosome as a string, which helps to query atypical contigs that would be parsed
   * to -1 by {@link Positions#chromosomeNumber(String)} if needed
   */
  public static class Query {

    private final String chrom;
    private final int start;
    private final int end;

    /**
     * @param chrom
     * @param start
     * @param end
     */
    private Query(String chrom, int start, int end) {
      super();
      this.chrom = chrom;
      this.start = start;
      this.end = end;
    }

    /**
     * @return the chrom
     */
    public String getChrom() {
      return chrom;
    }

    /**
     * @return the start
     */
    public int getStart() {
      return start;
    }

    /**
     * @return the end
     */
    public int getEnd() {
      return end;
    }

  }

  private static List<Integer> nSplit(int size, long aproxBasePairsPerPiece) {
    int pieces = (int) Math.ceil((double) size / aproxBasePairsPerPiece);

    List<Integer> bpSplits = new ArrayList<>();
    int stop = size - (size / pieces);
    for (int i = 1; i < stop; i += size / pieces) {
      bpSplits.add(i);
    }
    bpSplits.add(size + 1);

    return bpSplits;
  }

  /**
   * @param samSequenceDictionary {@link SAMSequenceDictionary} to break into chunks
   * @param aproximateNumberOfPieces the dictionary will be broken into ~ this many chunks
   * @param log
   * @return a list of {@link Query}
   */
  public static List<Query> parseToQueries(SAMSequenceDictionary samSequenceDictionary,
                                           int aproximateNumberOfPieces, Logger log) {
    List<Query> queries = new ArrayList<>();

    long totalLength = samSequenceDictionary.getReferenceLength();

    long aproxBasePairsPerPiece = totalLength / aproximateNumberOfPieces;

    log.reportTime("Total length of genome = " + totalLength + " over "
                   + samSequenceDictionary.getSequences().size() + " contigs");

    log.reportTime("Targeting ~ " + aproxBasePairsPerPiece + " base pairs per chunk");

    for (SAMSequenceRecord record : samSequenceDictionary.getSequences()) {
      List<Integer> splits = nSplit(record.getSequenceLength(), aproxBasePairsPerPiece);
      for (int j = 0; j < splits.size(); j++) {

        if (j > 0) {
          Query q = new Query(record.getSequenceName(), splits.get(j - 1), splits.get(j));
          queries.add(q);
        }
      }
      if (splits.get(splits.size() - 1) < record.getSequenceLength()) {

        throw new IllegalArgumentException("not totally covered, total ="
                                           + record.getSequenceLength());
      }
    }

    log.reportTimeInfo("Split to " + queries.size() + " total sections");

    return queries;

  }

}
