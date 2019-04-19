package org.genvisis.seq.manage;

import java.util.ArrayList;
import java.util.List;

import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.bed.SimpleBEDFeature;

public class AnnotatedBEDFeature extends SimpleBEDFeature implements BEDFeature {

  public AnnotatedBEDFeature(String chr, int start, int end) {
    super(start, end, chr);
    this.annots = new ArrayList<>();
  }

  List<String> annots;

  public void addAnnotation(String a) {
    this.annots.add(a);
  }

  public int count() {
    return this.annots.size();
  }

  public String getAnnotation(int index) {
    return this.annots.get(index);
  }

}
