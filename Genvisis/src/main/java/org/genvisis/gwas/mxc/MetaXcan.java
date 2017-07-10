package org.genvisis.gwas.mxc;

import java.io.File;
import java.util.ArrayList;

import org.genvisis.cnv.Launch;
import org.genvisis.common.Files;
import org.genvisis.common.ext;
import org.genvisis.qsub.Qsub;

public class MetaXcan {
	private static void runMany(String dataFolder, String extension, String mxc_folder,
															String refFile, String freqFile, String db, String covar,
															String posmap, String out, boolean overwrite, boolean verify) {
		String[] files = Files.list(dataFolder, extension, false);
		int numProcs = Math.min(24, files.length);

		ArrayList<String> commands = new ArrayList<String>();
		for (String s : files) {
			String c = "\njava -jar " + Launch.getJarLocation()
								 + " gwas.mxc.ParseMXCResults data=" + dataFolder + s + " mxc="
								 + mxc_folder + " ref=" + refFile + " freq=" + freqFile + " db=" + db + " covar="
								 + covar + " posmap=" + posmap + " out="
								 + ext.addToRoot(out, "_" + Files.removeExtention(s))
								 + (overwrite ? " -overwrite" : "")
								 + (verify ? " -verify" : "");
			commands.add(c);
		}

		Files.writeIterable(commands, "batchMXC.chain");

		String script = "cd " + new File("").getAbsolutePath() + "\n"
										+ "java -jar ~/" + org.genvisis.common.PSF.Java.GENVISIS
										+ " one.ScriptExecutor file=batchMXC.chain threads=" + numProcs + "\n"
										+ "java -jar ~/" + org.genvisis.common.PSF.Java.GENVISIS
										+ " gwas.mxc.ParseMXCResults -combine pattern=" + ext.rootOf(out, false) + "_";

		Qsub.qsub("batchMXC.pbs", script, 32000,
							8, numProcs);
	}

	public static void main(String[] args) {
		String usage = "Arguments: \n (1) File or folder of summary data to be analyzed (data=Metal_results.tbl (default))"
									 + "\n (2) Location of MetaXcan script (mxc=MetaXcan/software (default))"
									 + "\n (3) Output folder (eg out=results (default))"
									 + "\n (4) Covariance table for MetaXcan (covar=covariance.DGN-WB_0.5.txt.gz (default))"
									 + "\n (5) MetaXcan weights database (db=DGN-HapMap-2015/DGN-WB_0.5.db (default))"
									 + "\n (6) Overwrite existing MetaXcan output (-overwrite)"
									 + "\n (7) Verify allele order and strand before running MetaXcan (-verify)"
									 + "\n\t Reference file for expected allele order (ref=1000G.xln (default))"
									 + "\n\t File containing allele frequencies (freq=freq.tbl (default))";

		String posmap = "/panfs/roc/groups/5/pankrat2/mstimson/parkinsons/data/1000G_PD.map";
		String freqFile = "freq.tbl";
		String refFile = "1000G.xln";
		String data = "Metal_results.tbl";
		String mxc_folder = "MetaXcan/software/";
		String out = "results/mxc.csv";
		String db = "DGN-HapMap-2015/DGN-WB_0.5.db";
		String covar = "covariance.DGN-WB_0.5.txt.gz";
		String extension = ".tbl";

		boolean verify = false;
		boolean overwrite = false;


		for (String arg : args) {
			if (arg.startsWith("data=")) {
				data = ext.parseStringArg(arg);
			} else if (arg.startsWith("ext=")) {
				extension = ext.parseStringArg(arg);
			} else if (arg.startsWith("mxc=")) {
				mxc_folder = ext.parseStringArg(arg);
			} else if (arg.startsWith("out=")) {
				out = ext.parseStringArg(arg);
			} else if (arg.startsWith("db=")) {
				db = ext.parseStringArg(arg);
			} else if (arg.startsWith("covar=")) {
				covar = ext.parseStringArg(arg);
			} else if (arg.startsWith("-overwrite")) {
				overwrite = true;
			} else if (arg.startsWith("-verify")) {
				verify = true;
			} else if (arg.startsWith("freq=")) {
				freqFile = ext.parseStringArg(arg);
			} else if (arg.startsWith("ref=")) {
				refFile = ext.parseStringArg(arg);
			} else if (arg.startsWith("posmap=")) {
				posmap = ext.parseStringArg(arg);
			} else {
				System.err.println(usage);
				System.exit(0);
			}
		}

		mxc_folder = new File(mxc_folder).getAbsolutePath();
		refFile = new File(refFile).getAbsolutePath();
		freqFile = new File(freqFile).getAbsolutePath();
		covar = new File(covar).getAbsolutePath();
		posmap = new File(posmap).getAbsolutePath();
		data = new File(data).getAbsolutePath();

		if (!mxc_folder.endsWith("/"))
			mxc_folder += "/";

		if (Files.isDirectory(data)) {
			if (!data.endsWith("/"))
				data += "/";
			runMany(data, extension, mxc_folder, refFile, freqFile, db, covar, posmap, out, overwrite,
							verify);
		} else {
			String command = "cd " + ext.parseDirectoryOfFile(data) + "\njava -jar ~/"
											 + org.genvisis.common.PSF.Java.GENVISIS + " gwas.mxc.ParseMXCResults data="
											 + data + " mxc=" + mxc_folder + " ref=" + refFile + " freq=" + freqFile
											 + " db=" + db + " covar=" + covar + " posmap=" + posmap + " out=" + out
											 + (overwrite ? " -overwrite" : "") + (verify ? " -verify" : "");
			Qsub.qsub("metaXcan.qsub", command, 32000, 2, 1);
		}

	}
}
