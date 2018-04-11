package org.genvisis.one.JL;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import org.genvisis.cnv.manage.Resources;
import org.genvisis.cnv.manage.Resources.GENOME_BUILD;
import org.genvisis.common.Logger;
import org.genvisis.filesys.CNVariant;
import org.genvisis.filesys.GeneData;
import org.genvisis.filesys.GeneTrack;
import org.genvisis.filesys.LocusSet;
import org.genvisis.seq.manage.VCFOps.VcfPopulation;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;

public class SanityNICHD {

  //  private static List<String> geneNames = (List<String>) Arrays.asList(new String[] {"AIP", "MEN1",
  //                                                                                     "PRKAR1A",
  //                                                                                     "CDKN1B",
  //                                                                                     "CDKN2C",
  //                                                                                     "GPR101",
  //                                                                                     "USP8",
  //                                                                                     "CABLES1",
  //                                                                                     "TSC2", "BAI1",
  //                                                                                     "RASD1",
  //                                                                                     "DICER1"});
  private static List<String> geneNames = (List<String>) Arrays.asList(new String[] {"CDKN1B"});

  public static void main(String[] args) {
    String cnvFile = "/Volumes/Beta2/NGS/Cushings/cnvs/ExomeDepthAll.all.cnvs";

    VcfPopulation vpop = VcfPopulation.load("/Volumes/Beta2/NGS/Cushings/cnvs/clinicalAssoc/CUSHING_FREQ_ALL.vpop",
                                            POPULATION_TYPE.ANY, new Logger());

    String gtrackFile = Resources.genome(GENOME_BUILD.HG19, new Logger()).getGTrack().get();

    LocusSet<GeneData> geneSet = GeneTrack.load(gtrackFile).convertToLocusSet(new Logger());
    //    chr12:12870302-12875305
    //    chr12:12870203-12875316
    HashSet<String> notControls = new HashSet<>();
    notControls.add("EXCLUDE");
    notControls.add("EPP");
    ArrayList<CNVariant> cushings = new ArrayList<>();
    ArrayList<CNVariant> controls = new ArrayList<>();

    LocusSet<CNVariant> cnvs = CNVariant.loadLocSet(cnvFile, new Logger());
    System.out.println(cnvs.getLoci().length);
    GeneTrack g = GeneTrack.load(gtrackFile);
    for (String gene : geneNames) {

      GeneData[] genes = g.lookupAllGeneData(gene);

      for (GeneData gd : genes) {
        System.out.println(gd.getUCSClocation());
        CNVariant[] olaps = cnvs.getOverLappingLoci(gd);
        if (olaps != null) {
          for (CNVariant o : olaps) {
            System.out.println(o.toPlinkFormat());
          }
        }

      }

    }

  }

}
