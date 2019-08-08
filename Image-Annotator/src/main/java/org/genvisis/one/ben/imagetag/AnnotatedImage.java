package org.genvisis.one.ben.imagetag;

import java.util.ArrayList;

public class AnnotatedImage {

  public static class Annotation {

    final String annotation;
    final char mnemonic;

    public Annotation(String ann, String mnem) {
      this.annotation = ann;
      this.mnemonic = mnem.charAt(0);
    }

    public Annotation(String ann) {
      this.annotation = ann;
      this.mnemonic = Annotator.determineMnemonic(ann);
    }

    @Override
    public int hashCode() {
      final int prime = 31;
      int result = 1;
      result = prime * result + ((annotation == null) ? 0 : annotation.hashCode());
      return result;
    }

    @Override
    public boolean equals(Object obj) {
      if (this == obj) return true;
      if (obj == null) return false;
      if (getClass() != obj.getClass()) return false;
      Annotation other = (Annotation) obj;
      if (annotation == null) {
        if (other.annotation != null) return false;
      } else if (!annotation.equals(other.annotation)) return false;
      return true;
    }
  }

  private String name;
  private String imageFile;
  private ArrayList<AnnotatedImage.Annotation> annots;
  private final boolean isRoot;
  private boolean missing = false;

  private static final String DELIM = "|";

  public String exportToString() {
    StringBuilder sb = new StringBuilder();
    sb.append(name).append(DELIM).append(imageFile == null ? "" : imageFile);
    for (AnnotatedImage.Annotation a : annots) {
      sb.append(DELIM).append(a.annotation);
    }
    return sb.toString();
  }

  public AnnotatedImage(String ident, boolean isRoot) {
    this.name = ident;
    this.isRoot = isRoot;
    annots = new ArrayList<>();
  }

  public void setMissing(boolean miss) {
    missing = miss;
  }

  @Override
  public String toString() {
    return name + (missing ? " - MISSING" : "");
  }

  public boolean isRoot() {
    return isRoot;
  }

  public String getName() {
    return name;
  }

  public void setName(String name) {
    this.name = name;
  }

  public void setImageFile(String image) {
    this.imageFile = image;
  }

  public String getImageFile() {
    return this.imageFile;
  }

  public ArrayList<AnnotatedImage.Annotation> getAnnotations() {
    return this.annots;
  }
}
