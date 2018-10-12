package org.genvisis.one.JL.anotia;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Map;
import org.genvisis.cnv.LocusSet;
import org.genvisis.cnv.LocusSet.TO_STRING_TYPE;
import org.genvisis.cnv.filesys.CNVariant;
import org.genvisis.seq.cnv.CNVExtraInfo;
import org.genvisis.seq.cnv.CNVExtraInfo.EXTRA_INFO_TYPE;
import org.genvisis.seq.cnv.SeqCNVariant;
import org.genvisis.seq.manage.VCFOps.VcfPopulation;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.RETRIEVE_TYPE;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;
import org.pankratzlab.shared.filesys.GeneData;
import org.pankratzlab.shared.filesys.GeneTrack;
import org.pankratzlab.shared.filesys.Positions;

public class AnotiaCNVs {

  public static void main(String[] args) {
    String[] cnvFiles = new String[] {"/Volumes/Beta/data/ANOTIA/CNVs/CUSHINGS_FP_EXOME_DEPTH/results/ExomeDepthAll.all.cnvs",
                                      "/Volumes/Beta/data/ANOTIA/CNVs/CUSHINGS_FP_EXOME_DEPTH/results/ExomeDepthAll.all.noCNVR.cnvs"};

    for (String cnvFile : cnvFiles) {
      String vpopFile = "/Volumes/Beta/data/ANOTIA/CNVs/CUSHINGS_FP_EXOME_DEPTH/CUSHING_FP.vpop";
      String outDir = "/Volumes/Beta/data/ANOTIA/CNVs/CUSHINGS_FP_EXOME_DEPTH/cnvResults/";
      Logger log = new Logger(outDir + "log.log");
      VcfPopulation vpop = VcfPopulation.load(vpopFile, POPULATION_TYPE.ANY, log);
      Map<String, LocusSet<CNVariant>> set = CNVariant.breakIntoInds(CNVariant.loadLocSet(cnvFile,
                                                                                          log),
                                                                     log);
      ArrayList<CNVariant> anotia = new ArrayList<>();
      ArrayList<CNVariant> controls = new ArrayList<>();

      for (String ind : set.keySet()) {

        if (vpop.getPopulationForInd(ind.split("\t")[0], RETRIEVE_TYPE.SUPER) != null
            && vpop.getPopulationForInd(ind.split("\t")[0], RETRIEVE_TYPE.SUPER).length > 0) {
          String pop = vpop.getPopulationForInd(ind.split("\t")[0], RETRIEVE_TYPE.SUPER)[0];
          if (pop.equals("ANOTIA")) {
            Collections.addAll(anotia, set.get(ind).getLoci());

          } else if (!pop.equals("EXCLUDE")) {
            Collections.addAll(controls, set.get(ind).getLoci());
          }
        } else {
          log.reportTimeWarning("Sample " + ind + " not in vpop");
        }
      }

      LocusSet<CNVariant> controlSet = new LocusSet<>(controls, true, log);
      LocusSet<CNVariant> caseSet = new LocusSet<>(anotia, true, log);
      String outCNVsAll = outDir + ext.rootOf(cnvFile) + "AllCase.cnv";
      caseSet.writeRegions(outCNVsAll, TO_STRING_TYPE.REGULAR, true, log);

      String outCNVsAllControl = outDir + ext.rootOf(cnvFile) + "AllControl.cnv";
      controlSet.writeRegions(outCNVsAllControl, TO_STRING_TYPE.REGULAR, true, log);

      ArrayList<CNVariant> filteredAnotia = new ArrayList<>();

      for (CNVariant anotiaCnv : anotia) {
        if (anotiaCnv.getScore() > 5) {
          CNVariant[] overlaps = controlSet.getOverLappingLoci(anotiaCnv);
          if (overlaps == null || overlaps.length == 0) {
            filteredAnotia.add(anotiaCnv);
          } else {
            boolean add = true;
            for (CNVariant cnv : overlaps) {
              double bpOlap = (double) anotiaCnv.amountOfOverlapInBasepairs(cnv)
                              / anotiaCnv.getSize();
              // System.out.println(bpOlap);
              if (bpOlap > .5) {
                add = false;
              }
            }
            if (add) {
              filteredAnotia.add(anotiaCnv);
            }
          }
        }
      }
      String gdi = "/Volumes/Beta/data/ANOTIA/CNVs/CUSHINGS_FP_EXOME_DEPTH/GDI_percentile.txt";
      Hashtable<String, String> gdiHash = HashVec.loadFileToHashString(gdi, new int[] {0},
                                                                       new int[] {4}, false, "\t",
                                                                       true, true);
      LocusSet<GeneData> geneSet = GeneTrack.load("/Users/Kitty/workspace.other/Genvisis/Genvisis/resources/Genome/hg19/RefSeq_hg19.gtrack")
                                            .convertToLocusSet(log);
      ArrayList<SeqCNVariant> seqCNVariantsFiltered = new ArrayList<>();
      HashMap<String, HashSet<String>> counts = new HashMap<>();

      for (CNVariant filteredCNV : filteredAnotia) {
        GeneData[] genes = geneSet.getOverLappingLoci(filteredCNV);
        AnotiaEI geneInfo = new AnotiaEI(EXTRA_INFO_TYPE.EXOME_DEPTH, "GENE", "");
        AnotiaEI gdiInfo = new AnotiaEI(EXTRA_INFO_TYPE.EXOME_DEPTH, "GDI", "");
        AnotiaEI ucscLink = new AnotiaEI(EXTRA_INFO_TYPE.EXOME_DEPTH, "UCSC", "");
        int[] pos = new int[] {filteredCNV.getChr(), filteredCNV.getStart(), filteredCNV.getStop()};
        String link = Positions.getUCSClink(pos, "hg19");
        ucscLink.setdExtra("[" + filteredCNV.getUCSClocation() + "](" + link + ")");
        if (genes != null) {
          for (GeneData gene : genes) {
            if (!counts.containsKey(gene.getGeneName())) {
              counts.put(gene.getGeneName(), new HashSet<>());
            }
            geneInfo.setdExtra((geneInfo.getdExtra().equals("") ? gene.getGeneName()
                                                                : geneInfo.getdExtra() + ";"
                                                                  + gene.getGeneName()));
            counts.get(gene.getGeneName()).add(filteredCNV.getIndividualID());
            if (gdiHash.containsKey(gene.getGeneName())) {
              gdiInfo.setdExtra(gdiInfo.getdExtra().equals("") ? gdiHash.get(gene.getGeneName())
                                                               : gdiInfo.getdExtra() + ";"
                                                                 + gdiHash.get(gene.getGeneName()));
            } else {
              gdiInfo.setdExtra(gdiInfo.getdExtra().equals("") ? "NA"
                                                               : gdiInfo.getdExtra() + ";NA");

            }
            System.out.println(filteredCNV.getIndividualID() + "\t" + gene.getGeneName() + "\t"
                               + gdiHash.get(gene.getGeneName()));
          }
        }
        SeqCNVariant seqCNVariant = new SeqCNVariant(filteredCNV,
                                                     new AnotiaEI[] {geneInfo, gdiInfo, ucscLink});
        seqCNVariantsFiltered.add(seqCNVariant);
      }
      String outCounts = outDir + ext.rootOf(cnvFile) + "rareGeneCounts.txt";
      StringBuilder output = new StringBuilder("Gene\tNumAnotiaSamples\tSamples");
      for (String gene : counts.keySet()) {
        output.append("\n" + gene + "\t" + counts.get(gene).size() + "\t" + counts.get(gene));
      }
      Files.write(output.toString(), outCounts);

      LocusSet<SeqCNVariant> filtered = new LocusSet<>(seqCNVariantsFiltered, true, log);

      String outCNVs = outDir + ext.rootOf(cnvFile) + "rareGeneCounts.cnv";
      filtered.writeRegions(outCNVs, TO_STRING_TYPE.REGULAR, true, log);

    }
  }

  private static class AnotiaEI extends CNVExtraInfo {

    private static final long serialVersionUID = 1L;

    public AnotiaEI(EXTRA_INFO_TYPE type, String title, String extra) {
      super(EXTRA_INFO_TYPE.EXOME_DEPTH);
      dExtra = extra;
      sExtra = title;
    }
  }

}
