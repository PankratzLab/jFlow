package org.genvisis.one.JL.cushing;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Map;

import org.genvisis.cnv.annotation.segments.GDIAnnotator;
import org.genvisis.cnv.annotation.segments.GeneAnnotator;
import org.genvisis.cnv.annotation.segments.SegmentAnnotationKeys;
import org.genvisis.cnv.annotation.segments.SegmentAnotation;
import org.genvisis.cnv.annotation.segments.WESMappabilityAnnotator;
import org.genvisis.cnv.manage.Resources;
import org.genvisis.cnv.manage.Resources.GENOME_BUILD;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.filesys.CNVariant;
import org.genvisis.filesys.CNVariantAnnotated;
import org.genvisis.filesys.GeneData;
import org.genvisis.filesys.GeneTrack;
import org.genvisis.filesys.LocusSet;
import org.genvisis.filesys.CNVariantAnnotated.TallyResult;
import org.genvisis.filesys.LocusSet.TO_STRING_TYPE;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.ReferenceGenome;
import org.genvisis.seq.manage.VCFOps.VcfPopulation;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import org.genvisis.seq.manage.VCFOps.VcfPopulation.RETRIEVE_TYPE;

public class CNVScan {

	public static void main(String[] args) {

		String cnvFile = "/Volumes/Beta/data/Cushings/cnvs/ExomeDepthAll.all.noCNVR.cnvs";
		// String cnvFile =
		// "/Volumes/Beta/data/Cushings/cnvs/ExomeDepthAll.all.cnvs";

		String g1000Cnvs = "/Volumes/Beta/data/Cushings/cnvs/GRCh37_hg19_variants_2015-07-23.txt.cnv";
		String problemRegions = "/Volumes/Beta/data/Cushings/cnvs/problematicRegions_hg19.dat";

		String outDir = "/Volumes/Beta/data/Cushings/cnvs/filtered/";

		Logger log = new Logger(outDir + "log.log");
		LocusSet<Segment> probs = new LocusSet<>(Segment.loadUCSCregions(problemRegions, 0, false, log), true, log);

		VcfPopulation vpop = VcfPopulation.load("/Volumes/Beta/data/Cushings/cnvs/CUSHING_FREQ_V2.vpop",
				POPULATION_TYPE.ANY, log);
		vpop.report();
		HashSet<String> lqs = HashVec
				.loadFileToHashSet("/Volumes/Beta/data/Cushings/cnvs/CUSHING_FREQ_V2.lowerQualitySamples.txt", false);
		Hashtable<String, LocusSet<CNVariant>> set = CNVariant.breakIntoInds(
				LocusSet.combine(CNVariant.loadLocSet(cnvFile, log), CNVariant.loadLocSet(g1000Cnvs, log), true, log),
				log);
		// Hashtable<String, LocusSet<CNVariant>> set =
		// CNVariant.breakIntoInds(CNVariant.loadLocSet(cnvFile, log),log);
		HashSet<String> notControls = new HashSet<>();
		notControls.add("EXCLUDE");
		notControls.add("EPP");

		ArrayList<CNVariant> cushings = new ArrayList<>();
		ArrayList<CNVariant> controls = new ArrayList<>();

		splitIntoCaseControl(log, vpop, set, cushings, controls, "CUSHING_FREQ_V2", notControls, lqs);

		LocusSet<CNVariant> controlSet = new LocusSet<>(controls, true, log);

		String outCNVsAllControl = outDir + ext.rootOf(cnvFile) + "AllControl.cnv";
		controlSet.writeRegions(outCNVsAllControl, TO_STRING_TYPE.REGULAR, true, log);

		LocusSet<CNVariant> caseSet = new LocusSet<>(cushings, true, log);
		String outCNVsAll = outDir + ext.rootOf(cnvFile) + "AllCase.cnv";
		caseSet.writeRegions(outCNVsAll, TO_STRING_TYPE.REGULAR, true, log);
		double minQual = 10;
		double maxOverlap = 0.5;
		double minMap = 0.75;
		ArrayList<CNVariant> filteredCase = filterOutControls(cushings, controlSet, minQual, maxOverlap);

		GeneAnnotator geneAnnotator = GeneAnnotator.getDefaultAnnotator(GENOME_BUILD.HG19, log);
		GDIAnnotator gdiAnnotator = GDIAnnotator.getDefaultGDIAnnotator(log);
		WESMappabilityAnnotator wesMappabilityAnnotator = WESMappabilityAnnotator.getDefaultAnnotator(geneAnnotator,
				GENOME_BUILD.HG19, log);

		ArrayList<CNVariantAnnotated> annotateds = new ArrayList<>();

		for (CNVariant cnv : filteredCase) {
			SegmentAnotation segmentAnotation = new SegmentAnotation();

			wesMappabilityAnnotator.annotate(cnv, segmentAnotation);
			double mapScore = Double.parseDouble(
					segmentAnotation.getAttributes().get(SegmentAnnotationKeys.MAPPABILITY.toString()).get(0));
			// System.out.println(segmentAnotation.getAttributes().get(SegmentAnnotationKeys.MAPPABILITY.toString()));
			if (!Double.isNaN(mapScore) && mapScore > minMap) {
				geneAnnotator.annotate(cnv, segmentAnotation);
				gdiAnnotator.annotate(cnv, segmentAnotation);
				Segment[] probSegs = probs.getOverLappingLoci(cnv);
				boolean use = true;
				if (probSegs != null && probSegs.length > 0) {
					for (Segment problem : probSegs) {
						double bpOlap = (double) cnv.amountOfOverlapInBasepairs(problem) / cnv.getSize();
						if (bpOlap > maxOverlap) {
							use = false;
						}
					}
				}
				if (use) {
					CNVariantAnnotated cnVariantAnnotated = new CNVariantAnnotated(cnv, segmentAnotation);
					annotateds.add(cnVariantAnnotated);
				}
			}

		}
		wesMappabilityAnnotator.close();
		new LocusSet<>(annotateds, true, log).writeRegions(outDir + ext.rootOf(cnvFile) + ".annotated.cnvs",
				TO_STRING_TYPE.REGULAR, true, log);

		String outTally = outDir + ext.rootOf(cnvFile) + ".annotated.tally." + SegmentAnnotationKeys.GENE.toString();
		Map<String, TallyResult> geneTally = CNVariantAnnotated.tallyAnnotation(annotateds, SegmentAnnotationKeys.GENE);
		StringBuilder builderTally = new StringBuilder(SegmentAnnotationKeys.GENE.toString()
				+ "\tCOUNT_ALL\tCOUNT_DEL\tCOUNT_DUP\tSAME_CN\tGDI_PERCENTILE\tGDI_RAW\tGDI_PHRED\tSAMPLE\tLOC\tSCORE\tDISTANCE_TO_END_OF_CHR");
		ReferenceGenome genome = new ReferenceGenome(GENOME_BUILD.HG19, log);
		for (String key : geneTally.keySet()) {
			String gdi = gdiAnnotator.getGdiLookup().containsKey(key) ? gdiAnnotator.getGdiLookup().get(key)
					: (SegmentAnnotationKeys.GENE.getMissingValue() + "\t"
							+ SegmentAnnotationKeys.GENE.getMissingValue() + "\t"
							+ SegmentAnnotationKeys.GENE.getMissingValue());
			TallyResult tallyResult = geneTally.get(key);
			builderTally.append("\n" + key + "\t" + tallyResult.getAllCN().size() + "\t" + tallyResult.getDelCN().size()
					+ "\t" + tallyResult.getDupCN().size() + "\t"
					+ (tallyResult.getAllCN().size() == tallyResult.getDupCN().size()
							|| tallyResult.getAllCN().size() == tallyResult.getDelCN().size())
					+ "\t" + Array.toStr(gdi.replaceAll("%", "").split(";")) + "\t");
			int len = -1;
			if (tallyResult.getAllLocs().size() == 1) {
				for (CNVariant seg : tallyResult.getAllLocs()) {
					builderTally
							.append(":" + seg.getIndividualID() + "\t" + seg.getUCSClocation() + "\t" + seg.getScore());
					len = Math.min(seg.getStart(), genome.getContigLength(seg) - seg.getStop());
				}
			}
			builderTally.append("\t" + len);
		}
		Files.write(builderTally.toString(), outTally);

	}

	private static void splitIntoCaseControl(Logger log, VcfPopulation vpop, Hashtable<String, LocusSet<CNVariant>> set,
			ArrayList<CNVariant> cushings, ArrayList<CNVariant> controls, String caseDef, HashSet<String> notControls,
			HashSet<String> lq) {
		for (String ind : set.keySet()) {
			if (!lq.contains(ind.split("\t")[0])) {
				if (vpop.getPopulationForInd(ind.split("\t")[0], RETRIEVE_TYPE.SUPER) != null
						&& vpop.getPopulationForInd(ind.split("\t")[0], RETRIEVE_TYPE.SUPER).length > 0) {
					String pop = vpop.getPopulationForInd(ind.split("\t")[0], RETRIEVE_TYPE.SUPER)[0];
					if (pop.equals(caseDef)) {
						Collections.addAll(cushings, set.get(ind).getLoci());

					} else if (!notControls.contains(pop)) {
						Collections.addAll(controls, set.get(ind).getLoci());
					}
				} else {
					log.reportTimeWarning("Sample " + ind + " not in vpop");
				}
			} else {
				log.reportTime("removing " + ind);
			}
		}

	}

	private static ArrayList<CNVariant> filterOutControls(ArrayList<CNVariant> caseSet, LocusSet<CNVariant> controlSet,
			double minScore, double maxOverlap) {
		ArrayList<CNVariant> filtered = new ArrayList<>();
		for (CNVariant caseCNV : caseSet) {
			if (caseCNV.getScore() > minScore) {
				CNVariant[] overlaps = controlSet.getOverLappingLoci(caseCNV);
				if (overlaps == null || overlaps.length == 0) {
					filtered.add(caseCNV);
				} else {
					boolean add = true;
					for (CNVariant cnv : overlaps) {
						double bpOlap = (double) caseCNV.amountOfOverlapInBasepairs(cnv) / caseCNV.getSize();
						if (bpOlap > maxOverlap) {
							add = false;
						}
					}
					if (add) {
						filtered.add(caseCNV);
					}
				}
			}
		}
		return filtered;
	}

}
