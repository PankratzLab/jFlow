package org.genvisis.cnv.qc;

import java.util.Hashtable;

import org.genvisis.CLI;
import org.genvisis.cnv.Launch;
import org.genvisis.cnv.analysis.pca.PCMatrix;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

/**
 * Utility for applying estimated sex results from {@link SexChecks}.
 */
//TODO merge with SexChecks?
public class EstimatedSex {

	public static final String ARGS_PROJ = "proj";
	public static final String ARGS_USE_PED = "ped";


	/**
	 * Convert estimated sex values in the given {@link Project}'s sample data to definitive sex
	 * values.
	 */
	public static void apply(Project proj) {
		apply(proj, false);
	}

	/**
	 * As {@link #apply(Project)}. If {@code usePed} is true and a pedigree file exists for this
	 * project, any unknown estimated sex values will use the known value in the pedigree.
	 */
	public static void apply(Project proj, boolean usePed) {
		final String sampleFile = proj.SAMPLE_DATA_FILENAME.getValue();
		final Logger log = proj.getLog();

		String[] sampleHeader = Files.getHeaderOfFile(sampleFile, log);

		int sampleCol = ext.indexOfStr("DNA", sampleHeader, false, true);
		int estSexCol = ext.indexOfStr("CLASS=" + SexChecks.EST_SEX_HEADER, sampleHeader, false, true);

		if (estSexCol == -1) {
			log.reportError(EstimatedSex.class
			                + " - no estimated sex class found, can not apply estimated sex. Was SexChecks run?");
			return;
		}

		// Get estimated sex values for all samples
		Hashtable<String, String> sexMap = HashVec.loadFileToHashString(sampleFile, sampleCol,
		                                                                   new int[] {estSexCol}, "\t",
		                                                                   true, false);


		Hashtable<String, String> pedigreeMap = null;

		// Load pedigree file
		if (usePed) {
			final String pedFile = proj.PEDIGREE_FILENAME.getValue();
			if (Files.exists(pedFile)) {
				log.report("Loading Pedigree file, assuming standard pedigree.dat file format (FID, IID, FA, MO, SEX, PHENO, DNA)");
				pedigreeMap = HashVec.loadFileToHashString(pedFile, 6, new int[]{4}, "\t", false, false);
			}
		}

		// Map estimated sex to sex
		for (String k : sexMap.keySet()) {
			int v = SexChecks.getMappedSex(sexMap.get(k));
			if (v == 0 && pedigreeMap.containsKey(k)) {
				v = Integer.parseInt(pedigreeMap.get(k));
			}
			sexMap.put(k, Integer.toString(v));
		}

		SampleData sampleData = proj.getSampleData(SampleData.BASIC_CLASSES.length, false);
		sampleData.addData(sexMap, "DNA", new String[]{"CLASS=Sex"}, "0", null, log);
	}

	public static void main(String... args) {
		CLI c = new CLI(PCMatrix.class);
		c.addArgWithDefault(ARGS_PROJ, "project properties filename",
		                    Launch.getDefaultDebugProjectFile(false));
		c.addFlag(ARGS_USE_PED, "use pedigree for unknown sexes", false);

		c.parseWithExit(args);

		Project proj = new Project(c.get(ARGS_PROJ), false);
		apply(proj, c.has(ARGS_USE_PED));
	}
}
