/**
 *
 */
package org.genvisis.seq.manage.mtdna;

import java.io.File;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Bundled;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;

/**
 * @author Kitty
 */
public class RCRS {

  private final ReferenceSequence rcrsRef;
  private final String[] bases;
  public static final String RESOURCE = "NC_012920.1.fasta";
  private static final String RESOURCE_INDEX = RESOURCE + ".fai";
  private static final String RESOURCE_DICT = "NC_012920.1.dict";

  private RCRS(ReferenceSequence rcrsRef, Logger log) {
    this.rcrsRef = rcrsRef;
    bases = ArrayUtils.decodeByteArray(rcrsRef.getBases(), log);
  }

  /**
   * @param log
   * @return the {@link RCRS} mitochondriome
   */
  public static RCRS getRCRS(Logger log) {
    return new RCRS(loadRCRS(), log);
  }

  public ReferenceSequence getRcrsRef() {
    return rcrsRef;
  }

  public String[] getBases() {
    return bases;
  }

  /**
   * @param outputDir where the ref will be written
   * @param log
   */
  public static void writeRef(String outputDir, Logger log) {
    new File(outputDir).mkdirs();
    for (String toWrite : new String[] {RESOURCE, RESOURCE_INDEX, RESOURCE_DICT}) {
      File f = Bundled.getFile(toWrite);
      Files.copyFileUsingFileChannels(f, new File(outputDir + toWrite), log);
    }
  }

  /**
   * @return the {@link ReferenceSequence} for the rcrs mitochondriome, from within the jar
   */
  private static ReferenceSequence loadRCRS() {
    File f = Bundled.getFile("NC_012920.1.fasta");
    FastaSequenceFile rcrs = new FastaSequenceFile(f, true);
    ReferenceSequence rcrsRef = rcrs.nextSequence();
    rcrs.close();
    return rcrsRef;
  }
}
