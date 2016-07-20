package org.genvisis.jlDev;

import java.util.ArrayList;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.HashVec;
import org.genvisis.common.ext;
import org.genvisis.filesys.CNVariant;
import org.genvisis.filesys.LocusSet;
import org.genvisis.filesys.LocusSet.TO_STRING_TYPE;

public class ARIC {

	public static void main(String[] args) {
		Project proj = new Project("/home/pankrat2/lanej/projects/shadowAric_gw6.properties", false);
		String phenoDef = "/home/pankrat2/shared/aric_gw6/shadow_ARICGenvisis_CEL_FULL/vte/vte.txt";
		String cnvs = "/home/pankrat2/shared/aric_gw6/shadow_ARICGenvisis_CEL_FULL/vte/penncnv.cnv";
		String mapFile = "/panfs/roc/groups/5/pankrat2/public/bin/ref/mapability/hg18/wgEncodeCrgMapabilityAlign100mer.bedgraph";
		String geneTrackFile = "/home/pankrat2/public/bin/NCBI/RefSeq_hg18.gtrack";

		String[] vte = HashVec.loadFileToStringArray(phenoDef, false, new int[] { 0 }, true);

		LocusSet<CNVariant> cnvsSet = CNVariant.loadLocSet(cnvs, proj.getLog());
		ArrayList<CNVariant> caseCNVs = new ArrayList<CNVariant>();
		ArrayList<CNVariant> controlCNVs = new ArrayList<CNVariant>();
		SampleData sampleData = proj.getSampleData(0, false);

		String cases = ext.addToRoot(cnvs, ".vte");
		String controls = ext.addToRoot(cnvs, ".controls");
		int numExcluded =0;
		for (int i = 0; i < cnvsSet.getLoci().length; i++) {
			String dna = cnvsSet.getLoci()[i].getFamilyID();
			if (!sampleData.individualShouldBeExcluded(dna)) {
				if (ext.indexOfStr(dna, vte) >= 0) {
					caseCNVs.add(cnvsSet.getLoci()[i]);
				} else {
					controlCNVs.add(cnvsSet.getLoci()[i]);
				}
			}else{
				numExcluded++;
			}

		}
		proj.getLog().reportTimeInfo("excluded " + numExcluded);

		LocusSet<CNVariant> caseSet = new LocusSet<CNVariant>(caseCNVs.toArray(new CNVariant[caseCNVs.size()]), true, proj.getLog()) {

			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

		};
		caseSet.writeRegions(cases, TO_STRING_TYPE.REGULAR, true, proj.getLog());
		LocusSet<CNVariant> controlSet = new LocusSet<CNVariant>(controlCNVs.toArray(new CNVariant[controlCNVs.size()]), true, proj.getLog()) {

			/**
			 * 
			 */
			private static final long serialVersionUID = 1L;

		};
		controlSet.writeRegions(controls, TO_STRING_TYPE.REGULAR, true, proj.getLog());

		cnvsSet.getStrictSegmentSet().writeRegions(cnvs + "regions.txt", TO_STRING_TYPE.REGULAR, false, proj.getLog());
		CushingCnvs.filter(proj, mapFile, cases, new String[] { controls }, geneTrackFile, cnvs + "regions.txt", proj.getLog());

	}

}