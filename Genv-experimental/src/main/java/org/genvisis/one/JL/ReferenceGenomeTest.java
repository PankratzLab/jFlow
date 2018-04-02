package org.genvisis.one.JL;

import java.util.ArrayList;
import org.genvisis.cnv.filesys.MarkerSetInfo;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.ReferenceGenome;

public class ReferenceGenomeTest {

  public static void testRef(Project proj) {

    // basic
    ReferenceGenome referenceGenome = new ReferenceGenome(proj.getReferenceGenomeFASTAFilename(),
                                                          proj.getLog());
    // String[] test1 = referenceGenome.getSequenceFor(new Segment((byte) 26, 1, 50));

    // other testing
    String[] test = referenceGenome.getSequenceFor(new Segment((byte) 26, 1, 50));
    String[] bah = Files.getFirstNLinesOfFile(proj.getReferenceGenomeFASTAFilename(), 2,
                                              proj.getLog());
    System.out.println(ArrayUtils.toStr(test));
    System.out.println(ArrayUtils.toStr(bah));

    for (int i = 0; i < test.length; i++) {
      System.out.println(test[i] + "\t" + bah[1].charAt(i));// these should be equal
    }
    System.out.println(ArrayUtils.toStr(test));
    System.out.println(ArrayUtils.toStr(bah));
    ArrayList<Segment> tseg = new ArrayList<>();

    for (int i = 1; i < 23; i++) {
      tseg.add(new Segment((byte) i,
                           referenceGenome.getIndexedFastaSequenceFile().getSequenceDictionary()
                                          .getSequence(i).getSequenceLength()
                                     - 10000,
                           referenceGenome.getIndexedFastaSequenceFile().getSequenceDictionary()
                                          .getSequence(i).getSequenceLength() - 9800));
    }

    int off = tseg.size();
    for (int i = 1; i < 23; i++) {
      Segment segment1 = tseg.get(i);
      Segment segment2 = tseg.get(off - i);

      tseg.add(segment1);
      tseg.add(segment2);
    }
    for (int i = 1; i < 23; i++) {
      tseg.add(new Segment((byte) i, 1, 100));
    }

    long time = System.currentTimeMillis();
    proj.getLog().reportTimeInfo("Starting da query");
    for (int i = 0; i < tseg.size(); i++) {
      long disTime = System.currentTimeMillis();
      System.out.println("On: " + tseg.get(i).getUCSClocation());
      String[] bases = referenceGenome.getSequenceFor(tseg.get(i));
      if (ArrayUtils.unique(bases).length > 2) {
        System.out.println(tseg.get(i).getUCSClocation() + "\t" + ArrayUtils.toStr(bases, ""));
      }
      proj.getLog().reportTimeElapsed(disTime);
    }
    proj.getLog().reportTimeElapsed(time);

    Segment[] markerSegs = new Segment[proj.getMarkerNames().length];
    MarkerSetInfo markerSet = proj.getMarkerSet();
    for (int i = 0; i < markerSegs.length; i++) {
      markerSegs[i] = new Segment(markerSet.getChrs()[i], markerSet.getPositions()[i] - 50,
                                  markerSet.getPositions()[i] + 50);
    }
    proj.getLog().reportTimeInfo("Going for every marker in the project n=" + markerSegs.length);
    time = System.currentTimeMillis();

    for (int i = 0; i < markerSegs.length; i++) {
      if (i % 1000 == 0) {
        proj.getLog().reportTimeInfo(i + "");
      }
      referenceGenome.getSequenceFor(markerSegs[i]);
    }
    proj.getLog().reportTimeElapsed(time);

  }

  public static void testRef2(Project proj) {

    // basic
    ReferenceGenome referenceGenome = new ReferenceGenome(proj.getReferenceGenomeFASTAFilename(),
                                                          proj.getLog());
    Segment[] markerSegs = new Segment[proj.getMarkerNames().length];
    MarkerSetInfo markerSet = proj.getMarkerSet();
    for (int i = 0; i < markerSegs.length; i++) {
      markerSegs[i] = new Segment(markerSet.getChrs()[i], markerSet.getPositions()[i] - 50,
                                  markerSet.getPositions()[i] + 50);
    }
    proj.getLog().reportTimeInfo("Going for every marker in the project n=" + markerSegs.length);
    long time = System.currentTimeMillis();
    for (int i = 0; i < markerSegs.length; i++) {
      if (i % 1000 == 0) {
        proj.getLog().reportTimeInfo(i + "");
      }
      referenceGenome.getSequenceFor(markerSegs[i]);
    }
    proj.getLog().reportTimeElapsed(time);
    proj.getLog().reportTimeInfo(" n=" + markerSegs.length);

  }

  public static void main(String[] args) {
    Project proj = new Project(args[0]);
    testRef2(proj);
  }
}
