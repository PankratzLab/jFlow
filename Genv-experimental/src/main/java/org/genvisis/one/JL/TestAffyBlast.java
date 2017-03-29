package org.genvisis.one.JL;

import org.genvisis.cnv.filesys.ABLookup;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.qc.AffyMarkerBlast;

public class TestAffyBlast {

	public static void main(String[] args) {
		Project proj = new Project("/home/pankrat2/lane0212/projects/Aric_gw6.properties", false);
		String probeFile = "/home/pankrat2/public/bin/Arrays/affy_gw6/GenomeWideSNP_6.probe_tab";
		String annotFile = "/home/pankrat2/public/bin/Arrays/affy_gw6/GenomeWideSNP_6.na35.annot.csv";
		proj.BLAST_ANNOTATION_FILENAME.setValue(proj.DATA_DIRECTORY.getValue() + "STKblastTest.vcf.gz");

		new AffyMarkerBlast(proj, 20, 20, AffyMarkerBlast.DEFAULT_MAX_ALIGNMENTS_REPORTED,
												AffyMarkerBlast.DEFAULT_REPORT_TO_TEMPORARY_FILE,
												AffyMarkerBlast.DEFAULT_ANNOTATE_GC_CONTENT,
												AffyMarkerBlast.DEFAULT_DO_BLAST, proj.NUM_THREADS.getValue(), probeFile,
												annotFile).blastEm();
		new ABLookup().parseFromAnnotationVCF(proj);
	}

}
