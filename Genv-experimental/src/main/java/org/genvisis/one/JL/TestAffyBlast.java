package org.genvisis.one.JL;

import org.genvisis.cnv.filesys.ABLookup;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.qc.AffyMarkerBlast;

public class TestAffyBlast {

	public static void main(String[] args) {
		Project proj = new Project("/home/pankrat2/lanej/projects/Aric_gw6.properties", false);
		String probeFile = "/home/pankrat2/public/bin/Arrays/affy_gw6/GenomeWideSNP_6.probe_tab";
		String annotFile = "/home/pankrat2/public/bin/Arrays/affy_gw6/GenomeWideSNP_6.na35.annot.csv";
		proj.BLAST_ANNOTATION_FILENAME.setValue(proj.DATA_DIRECTORY.getValue() + "JLblastTest.vcf.gz");

		new AffyMarkerBlast(proj, -1, -1, 1, false, false, false, 1, probeFile,
												annotFile).blastEm();
		new ABLookup().parseFromAnnotationVCF(proj);



	}

}
