package org.genvisis.one.JL.assoc;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.ext;
import org.genvisis.qsub.Qsub;
import org.genvisis.stats.Rscript;

/**
 * Placeholder for threading out genome wide assoc scan
 *
 */
public class AssocTests {

	public static void main(String[] args) {
		String dir = "/scratch.global/lanej/burdenQA/full/";
		String inputsFile = "/scratch.global/lanej/burdenQA/full/SPSS_inputs.txt";
		String rscript = "/scratch.global/lanej/burdenQA/full/msiBurden.R";

		String[] inputs = Files.listFullPaths(dir, "gene_counts.txt");
		String rFile = ArrayUtils.toStr(HashVec.loadFileToStringArray(rscript, false, null, false), "\n");

		for (String input : inputs) {
			String tmp = rFile;
			tmp = tmp.replaceAll("inputsFile = REPLACE",
					"inputsFile = read.delim(\"" + inputsFile + "\",stringsAsFactors = FALSE)");
			tmp = tmp.replaceAll("BurdenCountAllAllFilesAll = REPLACE",
					"BurdenCountAllAllFilesAll =" + Rscript.generateRVector(new String[] { input }, true));

			String rscriptTmp = ext.addToRoot(input, ".rscript");
			String qsub = ext.addToRoot(input, ".pbs");
			Files.write(tmp, rscriptTmp);
			Qsub.qsub(qsub, "Rscript " + rscriptTmp, 15000, 124, 1);
		}
	}

	// dir with .counts.txt
	// threads
	// pheno.txt
	// out dir

}
