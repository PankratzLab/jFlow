package org.genvisis.one.JL;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.hmm.CNVCaller;
import org.genvisis.cnv.hmm.CNVCaller.PFB_MANAGEMENT_TYPE;

/**
 * Testing to correct https://github.com/PankratzLab/Genvisis/issues/88
 *
 */
public class TestCNVCaller {

	public static void main(String[] args) {
		Project proj = new Project(	"/Users/Kitty/workspace.other/Genvisis/Genvisis/projects/gedi_gwas.properties",
																false);

		CNVCaller.callAutosomalCNVs(proj, "test.cnvs", new String[] {"4874928039_R01C01"}, null, null,
																CNVCaller.DEFUALT_MIN_SITES, CNVCaller.DEFUALT_MIN_CONF,
																PFB_MANAGEMENT_TYPE.PENNCNV_DEFAULT, 2, 4);
	}

}
