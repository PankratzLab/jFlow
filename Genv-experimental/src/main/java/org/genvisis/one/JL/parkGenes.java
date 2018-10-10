package org.genvisis.one.JL;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import org.genvisis.cnv.LocusSet;
import org.genvisis.cnv.Resources;
import org.genvisis.cnv.Resources.GENOME_BUILD;
import org.genvisis.cnv.filesys.CNVariant;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.shared.filesys.GeneData;
import org.pankratzlab.shared.filesys.GeneTrack;
import org.pankratzlab.shared.filesys.Segment;

public class parkGenes {

  public static void main(String[] args) {
    GeneTrack geneTrack = GeneTrack.load(Resources.genome(GENOME_BUILD.HG19, new Logger())
                                                  .getGTrack().get());

    LocusSet<CNVariant> set = CNVariant.loadLocSet("D:/data/NGRC/cnvs/decentCalls_centromeresBroken.cnv",
                                                   new Logger());
    ArrayList<CNVariant> found = new ArrayList<>();
    ArrayList<String> genes = new ArrayList<>();
    genes.add("HTRA2");
    genes.add("SNCA");
    genes.add("PINK1");
    genes.add("LRRK2");
    genes.add("ATP13A2");
    genes.add("PARK7");
    genes.add("PARK2");
    ArrayList<GeneData> genLocs = new ArrayList<>();
    for (String gene : genes) {
      if (geneTrack.lookupAllGeneData(gene).length == 0) {
        throw new IllegalArgumentException(gene);
      }
      GeneData[] d = geneTrack.lookupAllGeneData(gene);
      for (GeneData element : d) {
        genLocs.add(element);
      }
    }
    Segment[] segs = genLocs.toArray(new GeneData[genLocs.size()]);
    Arrays.sort(segs);
    for (int i = 0; i < set.getLoci().length; i++) {
      if (Segment.overlapsAny(set.getLoci()[i], segs)) {
        found.add(set.getLoci()[i]);
      }
    }
    System.out.println(found.size());

    String out = "D:/data/NGRC/cnvs/pdGenes.cnv";
    try {
      PrintWriter writer = Files.openAppropriateWriter(out);
      writer.println(ArrayUtils.toStr(CNVariant.PLINK_CNV_HEADER) + "\tGENE");
      for (CNVariant cnv : found) {
        for (GeneData gene : genLocs) {
          if (cnv.overlaps(gene)) {
            writer.println(cnv.toPlinkFormat() + "\t" + gene.getGeneName());
            System.out.println("HDF");
          }
        }
      }
      writer.close();
    } catch (Exception e) {

    }
  }
}
