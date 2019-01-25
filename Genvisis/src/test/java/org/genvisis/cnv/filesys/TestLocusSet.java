package org.genvisis.cnv.filesys;

import org.junit.Test;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.filesys.LocusSet;
import org.pankratzlab.common.filesys.Segment;
import junit.framework.Assert;

/**
 * Unit tests for the {@link LocusSet} class
 */

public class TestLocusSet {

  /**
   * @see LocusSet#hasNoOverlap()
   */
  @Test
  public void testHasNoOverlap() {
    LocusSet<Segment> testSet = new LocusSet<>(new Segment[] {new Segment((byte) 1, 10, 1000),
                                                              new Segment((byte) 1, 15, 100)},
                                               true, new Logger());
    Assert.assertFalse("Overlap exists, expected hasNoOverlap() to return false",
                       testSet.hasNoOverlap());

    testSet = new LocusSet<>(new Segment[] {new Segment((byte) 1, 10, 1000),
                                            new Segment((byte) 2, 15, 100)},
                             true, new Logger());
    Assert.assertTrue("Overlap does not exist, expected hasNoOverlap() to return true",
                      testSet.hasNoOverlap());
  }

  /**
   * @see LocusSet#hasOverlap()
   */
  @Test
  public void testHasOverlap() {
    LocusSet<Segment> testSet = new LocusSet<>(new Segment[] {new Segment((byte) 1, 10, 1000),
                                                              new Segment((byte) 1, 15, 100)},
                                               true, new Logger());
    Assert.assertTrue("Overlap exists, expected hasOverlap() to return true", testSet.hasOverlap());

    testSet = new LocusSet<>(new Segment[] {new Segment((byte) 1, 10, 1000),
                                            new Segment((byte) 2, 15, 100)},
                             true, new Logger());
    Assert.assertFalse("Overlap does not exist, expected hasOverlap() to return false",
                       testSet.hasOverlap());
  }

}
