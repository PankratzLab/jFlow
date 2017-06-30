package org.genvisis.gwas.mxc;

import java.io.File;

import org.genvisis.common.CmdLine;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class MetaXcan {
	// MetaXcan python call:
	// ./MetaXcan.py \
	// --model_db_path ../../mxc_presets/DGN-HapMap-2015/DGN-WB_0.5.db \
	// --covariance data/covariance.DGN-WB_0.5.txt.gz \
	// --gwas_folder ../../data \
	// --gwas_file_pattern ".*tbl" \
	// --beta_column Effect \
	// --pvalue_column "P-value" \
	// --output_file ../../mxc_results/zscore_run.csv \
	// --effect_allele_column Allele2 \
	// --non_effect_allele_column Allele1 \
	// --snp_column MarkerName \
	// --beta_sign_column Direction \
	// --overwrite


	public static void main(String[] args) {
		String usage = "Arguments \n (1) File of summary data to be analyzed (file=Metal_results.tbl (default))"
									 + "\n (2) Location of MetaXcan script (mxc=MetaXcan/software (default))"
									 + "\n (3) Output folder (eg out=results (default))"
									 + "\n (4) Covariance table for MetaXcan (covar=covariance.DGN-WB_0.5.txt.gz (default))"
									 + "\n (5) MetaXcan weights database (db=DGN-HapMap-2015/DGN-WB_0.5.db (default))"
									 + "\n (6) Beta column name (beta_column=Effect (default))"
									 + "\n (7) snp column name (snp_column=MarkerName (default))"
									 + "\n (8) Effect Allele column (effect_allele=Allele1 (default))"
									 + "\n (9) Non-effect Allele column (non_effect_allele=Allele2 (default))"
									 + "\n (10) Standard Error column (se=se (default))"
									 + "\n (11) Overwrite existing MetaXcan output (-overwrite)";


		String data = "Metal_results.tbl";
		String mxc_folder = "MetaXcan/software";
		String out = "results/mxc.csv";
		String db = "DGN-HapMap-2015/DGN-WB_0.5.db";
		String covar = "covariance.DGN-WB_0.5.txt.gz";
		String beta_column = "Effect";
		String snp_column = "MarkerName";
		String a1 = "Allele1";
		String a2 = "Allele2";
		String se = "se";
		String posmap = "/panfs/roc/groups/5/pankrat2/mstimson/parkinsons/data/1000G_PD.map";


		boolean overwrite = false;

		for (String arg : args) {
			if (arg.startsWith("data=")) {
				data = ext.parseStringArg(arg);
			} else if (arg.startsWith("mxc=")) {
				mxc_folder = ext.parseStringArg(arg);
			} else if (arg.startsWith("out=")) {
				out = ext.parseStringArg(arg);
			} else if (arg.startsWith("db=")) {
				db = ext.parseStringArg(arg);
			} else if (arg.startsWith("covar=")) {
				covar = ext.parseStringArg(arg);
			} else if (arg.startsWith("beta_column=")) {
				beta_column = ext.parseStringArg(arg);
			} else if (arg.startsWith("snp_column=")) {
				snp_column = ext.parseStringArg(arg);
			} else if (arg.startsWith("effect_allele=")) {
				a1 = ext.parseStringArg(arg);
			} else if (arg.startsWith("non_effect_allele=")) {
				a2 = ext.parseStringArg(arg);
			} else if (arg.startsWith("-overwrite")) {
				overwrite = true;
			} else if (arg.startsWith("se=")) {
				se = ext.parseStringArg(arg);
			} else {
				System.err.println(usage);
				System.exit(0);
			}
		}

		String gwas_folder = ext.parseDirectoryOfFile(new File(data).getAbsolutePath());

		db = new File(db).getAbsolutePath();
		covar = new File(covar).getAbsolutePath();

		out = new File(out).getAbsolutePath();

		String py = "MetaXcan.py";

		// build mxc command
		String command = "./" + py
										 + " --model_db_path " + db + " --covariance " + covar + " --gwas_folder "
										 + gwas_folder + " --gwas_file_pattern " + new File(data).getName()
										 + " --beta_column "
										 + beta_column + " --output_file " + out + " --effect_allele_column " + a1
										 + " --non_effect_allele_column " + a2 + " --snp_column " + snp_column
										 + " --se_column " + se
										 + (overwrite ? " --overwrite" : "");

		System.out.println(command);

		// run MetaXcan on the given inputs
		boolean runSuccess = CmdLine.run(command, ext.parseDirectoryOfFile(mxc_folder), null, null,
																		 new Logger("MetaXcan.log"), false);

		if (!runSuccess || !new File(out).exists()) {
			System.out.println("Error running MetaXcan with the given inputs.");
			System.exit(0);
		}


		// take the output mxc file and find the number of hits for each gene range
		ParseMXCResults.addMetalHits(posmap, out, covar, data,
																 new Logger("parseMXC.log"));


	}
}
