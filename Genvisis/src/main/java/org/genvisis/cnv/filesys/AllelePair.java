package org.genvisis.cnv.filesys;

import java.io.Serializable;
import javax.annotation.Nullable;
import htsjdk.variant.variantcontext.Allele;

public class AllelePair implements Serializable, Comparable<AllelePair> {

  public enum RefAllele {
    A, B
  }

  private static final long serialVersionUID = 1L;

  private final Allele alleleA;
  private final Allele alleleB;
  private final char charA;
  private final char charB;
  private final AllelePair.RefAllele refAllele;

  private AllelePair(Allele a, Allele b) {
    this.alleleA = a;
    this.alleleB = b;
    if (a.isReference()) {
      if (b.isNonReference()) {
        refAllele = RefAllele.A;
      } else throw new IllegalArgumentException("Cannot construct " + this.getClass().getName()
                                                + " with two reference alleles");
    } else if (b.isReference()) {
      refAllele = RefAllele.B;
    } else {
      refAllele = null;
    }
    char[] ab = ABLookup.parseABFromAlleles(a, b);
    this.charA = ab[0];
    this.charB = ab[1];
  }

  private AllelePair(char a, char b, @Nullable AllelePair.RefAllele refAllele) {
    this.alleleA = parseAllele(a, RefAllele.A.equals(refAllele));
    this.alleleB = parseAllele(b, RefAllele.B.equals(refAllele));
    this.charA = a;
    this.charB = b;
    this.refAllele = refAllele;
  }

  private static Allele parseAllele(char charAllele, boolean reference) {
    try {
      return Allele.create(String.valueOf(charAllele), reference);
    } catch (IllegalArgumentException iae) {
      return Allele.NO_CALL;
    }
  }

  /**
   * @param a Allele for A (reference status is used to set ref allele)
   * @param b Allele for B (reference status is used to set ref allele)
   * @return {@link AllelePair} of a and b
   */
  public static AllelePair of(Allele a, Allele b) {
    return new AllelePair(a, b);
  }

  /**
   * @param a char for allele A
   * @param b char for allele B
   * @param refAllele identifies which allele is ref
   * @return {@link AllelePair} of a and b with refAllele
   */
  public static AllelePair of(char a, char b, AllelePair.RefAllele refAllele) {
    return new AllelePair(a, b, refAllele);
  }

  /**
   * @param a char for allele A
   * @param b char for allele B
   * @return {@link AllelePair} of a and b with no ref allele
   */
  public static AllelePair of(char a, char b) {
    return new AllelePair(a, b, null);
  }

  public Allele getAlleleA() {
    return alleleA;
  }

  public Allele getAlleleB() {
    return alleleB;
  }

  public char getA() {
    return charA;
  }

  public char getB() {
    return charB;
  }

  public AllelePair.RefAllele getRefAllele() {
    return refAllele;
  }

  public Allele getRef() {
    if (refAllele == null) {
      return null;
    }
    switch (refAllele) {
      case A:
        return alleleA;
      case B:
        return alleleB;
      default:
        throw new IllegalStateException("Undefined refAllele instance:" + refAllele.toString());

    }
  }

  public Allele getAlt() {
    if (refAllele == null) {
      return null;
    }
    switch (refAllele) {
      case A:
        return alleleB;
      case B:
        return alleleA;
      default:
        throw new IllegalStateException("Undefined refAllele instance:" + refAllele.toString());
    }
  }

  @Override
  public int hashCode() {
    final int prime = 31;
    int result = 1;
    result = prime * result + ((alleleA == null) ? 0 : alleleA.hashCode());
    result = prime * result + ((alleleB == null) ? 0 : alleleB.hashCode());
    result = prime * result + charA;
    result = prime * result + charB;
    result = prime * result + ((refAllele == null) ? 0 : refAllele.hashCode());
    return result;
  }

  @Override
  public boolean equals(Object obj) {
    if (this == obj) return true;
    if (obj == null) return false;
    if (!(obj instanceof AllelePair)) return false;
    AllelePair other = (AllelePair) obj;
    if (alleleA == null) {
      if (other.alleleA != null) return false;
    } else if (!alleleA.equals(other.alleleA)) return false;
    if (alleleB == null) {
      if (other.alleleB != null) return false;
    } else if (!alleleB.equals(other.alleleB)) return false;
    if (charA != other.charA) return false;
    if (charB != other.charB) return false;
    if (refAllele != other.refAllele) return false;
    return true;
  }

  @Override
  public int compareTo(AllelePair o) {
    int cmp = alleleA.compareTo(o.alleleA);
    if (cmp != 0) return cmp;
    cmp = alleleB.compareTo(o.alleleB);
    if (cmp != 0) return cmp;
    if (refAllele != null && o.refAllele != null) {
      cmp = refAllele.compareTo(o.refAllele);
      return cmp;
    } else {
      if (refAllele == null && o.refAllele != null) return -1;
      if (refAllele != null && o.refAllele == null) return 1;
      return 0;
    }
  }
}
