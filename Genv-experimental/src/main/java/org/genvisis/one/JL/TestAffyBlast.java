package org.genvisis.one.JL;

import org.genvisis.cnv.filesys.ABLookup;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.qc.MarkerBlast;
import org.genvisis.cnv.qc.MarkerBlast.FILE_SEQUENCE_TYPE;

public class TestAffyBlast {
	
	public static void main(String[] args) {
		Project proj= new Project("/home/pankrat2/lanej/projects/Aric_gw6.properties",false);
		String seqFile = "/home/pankrat2/public/bin/Arrays/affy_gw6/GenomeWideSNP_6.probe_tab";
		proj.BLAST_ANNOTATION_FILENAME.setValue(proj.DATA_DIRECTORY.getValue()+"JLblastTest.vcf.gz");
		
		
		MarkerBlast.blastEm(proj, seqFile, FILE_SEQUENCE_TYPE.AFFY_ANNOT, -1, -1, 1, 1, false, false, false);
		new ABLookup().parseFromAnnotationVCF(proj);
		
		
		
	}

}
