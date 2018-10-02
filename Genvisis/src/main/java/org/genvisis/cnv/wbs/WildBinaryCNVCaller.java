/**
 * 
 */
package org.genvisis.cnv.wbs;

import org.genvisis.cnv.filesys.Project;
import org.pankratzlab.core.CLI;

/**
 * Class that utilizes Wild Binary Segmentation (WBS) to call CNVs
 * http://stats.lse.ac.uk/fryzlewicz/wbs/wbs.pdf Why WBS instead of Circular Binary Segmentation
 * (CBS)? It addresses similar issues as CBS (through an initial random interval scan - hence
 * "Wild") and is _much_ easier to implement. https://github.com/cran/wbs/ Compare to DNACopy
 * fortran +R code
 * https://github.com/Bioconductor-mirror/DNAcopy/tree/cf15bb37dbc08f607004f37cbeb4b70d1eece39b/src
 * and
 * https://github.com/Bioconductor-mirror/DNAcopy/tree/cf15bb37dbc08f607004f37cbeb4b70d1eece39b/R
 * Or, from Illumina Canvas
 * https://github.com/Illumina/canvas/tree/c9c1141927b4b541b8d48641aa20111aa1481a82/Src/Canvas/
 * CanvasPartition/CanvasPartition.net45
 */
public class WildBinaryCNVCaller {

  private static void callProject(Project proj) {

  }

  public static void main(String[] args) {
    CLI c = new CLI(WildBinaryCNVCaller.class);

    c.addArg(CLI.ARG_PROJ, CLI.DESC_PROJ, true);
    c.parseWithExit(args);
    // wake
    // callProject(new Project(filename, jar));

  }

}
