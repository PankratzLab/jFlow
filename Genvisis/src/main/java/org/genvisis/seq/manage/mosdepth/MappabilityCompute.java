package org.genvisis.seq.manage.mosdepth;

import java.util.Arrays;
import java.util.List;
import org.genvisis.common.Logger;
import org.genvisis.filesys.LocusSet;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.BEDFileReader;
import org.genvisis.seq.manage.VCOps;
import org.genvisis.seq.manage.BEDFileReader.BEDFeatureSeg;
import htsjdk.variant.variantcontext.VariantContext;

public class MappabilityCompute {

  public static final int DEFAULT_BP_BUFFER = 150;

  private String file;
  private Logger log = new Logger();
  LocusSet<BEDFeatureSeg> segSet;

  private MappabilityCompute(String file) {
    this.file = file;
    if (file.equals("") || file.trim().equals("")) {
      throw new RuntimeException("File cannot be blank.");
    }
  }

  private void load() {
    BEDFileReader reader = new BEDFileReader(file, false);
    segSet = reader.loadAll(log);
    reader.close();
  }

  public double getAverageMap(VariantContext vc) {
    return getAverageMap(vc, DEFAULT_BP_BUFFER);
  }

  public double getAverageMap(VariantContext vc, int bpBufferAroundVariant) {
    List<BEDFeatureSeg> overlappingMapEntries = Arrays.asList(segSet.getOverLappingLoci(VCOps.getSegment(vc)
                                                                                             .getBufferedSegment(bpBufferAroundVariant)));
    //    BEDFeatureSeg[] overlappingMapEntries = segSet.getOverLappingLoci(VCOps.getSegment(vc)
    //                                                                           .getBufferedSegment(bpBufferAroundVariant));
    switch (overlappingMapEntries.size()) {
      case 0:
        return 0;
      case 1:
        return overlappingMapEntries.get(0).getBedFeature().getScore();
      default:
        return getAverageMap(overlappingMapEntries,
                             (VCOps.getSegment(vc).getBufferedSegment(bpBufferAroundVariant)));
    }
  }

  private static double getAverageMap(List<BEDFeatureSeg> overlappingMapEntries, Segment query) {
    int totalBp = 0;
    double sum = 0;
    for (BEDFeatureSeg bFeatureSeg : overlappingMapEntries) {
      int intersection = bFeatureSeg.getIntersection(query, new Logger()).getSize();
      totalBp += intersection;
      sum += intersection * bFeatureSeg.getBedFeature().getScore();
    }
    return sum / totalBp;
  }

  public static MappabilityCompute open(String file) {
    MappabilityCompute mc = new MappabilityCompute(file);
    mc.load();
    return mc;
  }

}
