package org.genvisis.cnv.annotation;

import org.genvisis.common.Array;
import org.genvisis.stats.Histogram.DynamicHistogram;

import htsjdk.variant.vcf.VCFHeaderLineType;

public abstract class HistogramAnnotation extends AnnotationData {

  private final DynamicHistogram dynamicHistogram;

  public HistogramAnnotation(String name, String description, DynamicHistogram dynamicHistogram) {
    super(VCFHeaderLineType.String, null, 1, name, description, DEFUALT_VALUE, DEFUALT_VALUE);
    this.dynamicHistogram = dynamicHistogram;
  }

  public DynamicHistogram getDynamicHistogram() {
    return dynamicHistogram;
  }

  public void setDataToHistogram(boolean truncate) {
    if (truncate) {
      String truncatedHistogram = "";
      int index = 0;
      while (index < dynamicHistogram.getCounts().length
          && dynamicHistogram.getCounts()[index] == 0) {
        index++;
      }
      for (int i = index; i < dynamicHistogram.getCounts().length; i++) {
        truncatedHistogram +=
            (i == index ? "" : DEFUALT_DELIMITER) + dynamicHistogram.getCounts()[i];
      }
      setData(truncatedHistogram);
    } else {
      setData(Array.toStr(dynamicHistogram.getCounts(), DEFUALT_DELIMITER));
    }
  }

}
