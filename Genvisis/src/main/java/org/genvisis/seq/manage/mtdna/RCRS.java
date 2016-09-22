/**
 * 
 */
package org.genvisis.seq.manage.mtdna;

import java.io.File;

import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;

/**
 * @author Kitty
 *
 */
public class RCRS {
  private ReferenceSequence rcrsRef;
  private String[] bases;
  public static final String RESOURCE = "NC_012920.1.fasta";
  private static final String RESOURCE_INDEX = RESOURCE + ".fai";
  private static final String RESOURCE_DICT = "NC_012920.1.dict";

  private RCRS(ReferenceSequence rcrsRef, Logger log) {
    this.rcrsRef = rcrsRef;
    this.bases = Array.decodeByteArray(rcrsRef.getBases(), log);
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
    File f = new File(RCRS.class.getResource(RESOURCE).getFile());
    Files.copyFileUsingFileChannels(f, new File(outputDir + RESOURCE), log);
    File fr = new File(RCRS.class.getResource(RESOURCE_INDEX).getFile());
    Files.copyFileUsingFileChannels(fr, new File(outputDir + RESOURCE_INDEX), log);
    File frd = new File(RCRS.class.getResource(RESOURCE_DICT).getFile());
    Files.copyFileUsingFileChannels(fr, new File(outputDir + RESOURCE_DICT), log);
  }

  /**
   * @return the {@link ReferenceSequence} for the rcrs mitochondriome, from within the jar
   */
  private static ReferenceSequence loadRCRS() {
    File f = new File(RCRS.class.getResource("NC_012920.1.fasta").getFile());
    FastaSequenceFile rcrs = new FastaSequenceFile(f, true);
    ReferenceSequence rcrsRef = rcrs.nextSequence();
    rcrs.close();
    return rcrsRef;
  }
}
