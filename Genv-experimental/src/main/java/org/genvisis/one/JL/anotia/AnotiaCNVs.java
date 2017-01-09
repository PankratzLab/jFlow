package org.genvisis.one.JL.anotia;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;

import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.filesys.CNVariant;
import org.genvisis.filesys.GeneData;
import org.genvisis.filesys.GeneTrack;
import org.genvisis.filesys.LocusSet;
import org.genvisis.seq.cnv.CNVExtraInfo.EXTRA_INFO_TYPE;
import org.genvisis.seq.cnv.SeqCNVariant;
import org.genvisis.seq.manage.VCFOps.VcfPopulation;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.RETRIEVE_TYPE;

public class AnotiaCNVs {

	public static void main(String[] args) {
		String[] cnvFiles = new String[] {
				"/Volumes/Beta/data/ANOTIA/CNVs/CUSHINGS_FP_EXOME_DEPTH/results/ExomeDepthAll.all.cnvs",
				"/Volumes/Beta/data/ANOTIA/CNVs/CUSHINGS_FP_EXOME_DEPTH/results/ExomeDepthAll.all.noCNVR.cnvs" };

		for (String cnvFile : cnvFiles) {
			String vpopFile = "/Volumes/Beta/data/ANOTIA/CNVs/CUSHINGS_FP_EXOME_DEPTH/CUSHING_FP.vpop";
			String outDir = "/Volumes/Beta/data/ANOTIA/CNVs/CUSHINGS_FP_EXOME_DEPTH/cnvResults/";
			Logger log = new Logger(outDir + "log.log");
			VcfPopulation vpop = VcfPopulation.load(vpopFile, POPULATION_TYPE.ANY, log);
			Hashtable<String, LocusSet<CNVariant>> set = CNVariant.breakIntoInds(CNVariant.loadLocSet(cnvFile, log),
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

			ArrayList<CNVariant> filteredAnotia = new ArrayList<>();

			for (CNVariant anotiaCnv : anotia) {
				CNVariant[] overlaps = controlSet.getOverLappingLoci(anotiaCnv);
				if (overlaps == null || overlaps.length == 0) {
					filteredAnotia.add(anotiaCnv);
				} else {
					boolean add = true;
					for (CNVariant cnv : overlaps) {
						double bpOlap = (double) anotiaCnv.amountOfOverlapInBasepairs(cnv) / anotiaCnv.getSize();
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
			String gdi = "/Volumes/Beta/data/ANOTIA/CNVs/CUSHINGS_FP_EXOME_DEPTH/GDI_percentile.txt";
			Hashtable<String, String> gdiHash = HashVec.loadFileToHashString(gdi, new int[] { 0 }, new int[] { 4 },
					false, "\t", true, false, true);
			LocusSet<GeneData> geneSet = GeneTrack
					.load("/Users/Kitty/workspace.other/Genvisis/Genvisis/resources/Genome/hg19/RefSeq_hg19.gtrack",
							false)
					.convertToLocusSet(log);
			ArrayList<SeqCNVariant> seqCNVariantsFiltered = new ArrayList<>();

			for (CNVariant filteredCNV : filteredAnotia) {
				GeneData[] genes = geneSet.getOverLappingLoci(filteredCNV);
//				EXTRA_INFO_TYPE
				if (genes != null) {
					for (GeneData gene : genes) {
						System.out.println(filteredCNV.getIndividualID()+"\t"+gene.getGeneName() + "\t" + gdiHash.get(gene.getGeneName()));
					}
				}
			}
			LocusSet<CNVariant> filtered = new LocusSet<>(filteredAnotia, true, log);

		}
	}

}
