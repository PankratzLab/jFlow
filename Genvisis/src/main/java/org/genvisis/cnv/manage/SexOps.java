/**
 * 
 */
package org.genvisis.cnv.manage;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.qc.SexChecks;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.Array;

/**
 * Helper methods for loading/parsing sample sex from sampleData
 *
 */
public class SexOps {


	/**
	 *
	 *
	 */
	public enum SEX_LOAD_TYPE {
															/**
															 * 0,1,2
															 */
															MAPPED_SEX,
															/**
															 * Number of X chromosomes, facilitates regressing out triple X, etc
															 */
															NUM_X_SEX;
	}

	/**
	 * @param proj
	 * @return the 0,1,2 sample sex as found in sample data for all samples in the project
	 */
	public static int[] getSampleSex(Project proj, SEX_LOAD_TYPE sType) {
		String[] samples = proj.getSamples();
		SampleData sampleData = proj.getSampleData(0, false);
		int[] sex = new int[samples.length];
		for (int i = 0; i < samples.length; i++) {
			sex[i] = sampleData.getSexForIndividual(samples[i]);
			if (sex[i] < 0) {
				sex[i] = 0;
			}
			switch (sType) {
				case MAPPED_SEX:
					sex[i] = SexChecks.getMappedSex(sex[i] + "");
					break;
				case NUM_X_SEX:
					sex[i] = SexChecks.getNumXChrSex(sex[i] + "");
					break;

				default:
					throw new IllegalArgumentException("Invalid load type " + sType);

			}
		}
		proj.getLog().reportTimeInfo(Array.toStr(Array.unique(Array.toStringArray(sex))));


		return sex;
	}

}
