package org.genvisis.seq;

/**
 * Supported genome reference builds
 */
public enum GenomeBuild {
  HG38("hg38", 38, 2781479, 155701383),
  HG19("hg19", 37, 2699520, 154931044),
  HG18("hg18", 36, 2709521, 154584237);

  final String build;
  private final int buildInt;
  // Matched to the PAR definitions from plink: https://www.cog-genomics.org/plink/1.9/data#split_x
  private final int par1End;
  private final int par2Start;

  private GenomeBuild(String build, int buildInt, int par1_end, int par2_start) {
    this.build = build;
    this.buildInt = buildInt;
    this.par1End = par1_end;
    this.par2Start = par2_start;
  }

  public String getBuild() {
    return build;
  }

  public int getBuildInt() {
    return buildInt;
  }

  /**
   * @return the last position of x chromosome PAR 1
   */
  public int getPAR1End() {
    return par1End;
  }

  /**
   * @return the first position of x chromosome PAR 2
   */
  public int getPAR2Start() {
    return par2Start;
  }

}