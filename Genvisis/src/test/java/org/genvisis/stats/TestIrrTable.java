package org.genvisis.stats;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.qc.SampleQC;
import org.pankratzlab.shared.stats.IrrTable;
import org.pankratzlab.shared.stats.Quantiles;

/**
 * Manual test for {@link IrrTable}
 */
public class TestIrrTable {

  public static void test() {
    Project proj = new Project(null);
    SampleQC sampleQC = SampleQC.loadSampleQC(proj);
    Quantiles[] quantiles = Quantiles.qetQuantilesFor(100, sampleQC.getQcMatrix(),
                                                      sampleQC.getQctitles(), proj.getLog());
    IrrTable rIrrTable = new IrrTable(2, proj.getSamples().length, true, proj.getLog());
    rIrrTable.addRatings(0, quantiles[1].getQuantileMembershipAsRoundedInt());
    rIrrTable.addRatings(1, quantiles[1].getQuantileMembershipAsRoundedInt());
    rIrrTable.parseAgreement();
    System.out.println(rIrrTable.getPercentAgreementFor(10));
    System.out.println(rIrrTable.getCohensKappa());
  }

  public static void main(String[] args) {
    test();
  }

}
