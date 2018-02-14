package org.genvisis.one.spencer.wiCell;

import java.util.List;
import com.googlecode.charts4j.collect.Lists;

public class ExtractHapMapSnps {

  public static void main(String[] args) {
    String baseDir = "/scratch.global/topmed/58807/topmed-dcc/exchange/phs001211_TOPMed_WGS_ARIC/Combined_Study_Data/Genotypes/freeze.5b/minDP10/";
    String hapMapSites = "/home/pankrat2/shared/resources/HapMap/unambiguousHapMapFounders.sites";
    String geneSTARSamples = baseDir + "HapMapOverlap/GeneSTARSamples.txt";

    List<String> cmds = Lists.newArrayList();
    for (int i = 1; i < 24; i++) {
      if (i == 12) {
        for (int j = 0; j < 15; j++)
          System.out.println("wait");
      }
      String chr = i == 23 ? "X" : String.valueOf(i);
      String inFile = baseDir + "freeze.5b.chr" + chr + ".pass_and_fail.gtonly.minDP10.bcf";
      String outFile = baseDir + "HapMapOverlap/freeze.5b.chr" + chr
                       + ".pass.gtonly.minDP10.GeneSTAR.HapMapOverlap.bcf.gz";
      System.out.println("bcftools view -f PASS " + " -i '%ID=@" + hapMapSites + "' -S "
                         + geneSTARSamples + " -O b -o " + outFile + " " + inFile + " &");
    }
    for (int i = 0; i < 15; i++)
      System.out.println("wait");
  }

}
