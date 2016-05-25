package one.JL;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.concurrent.Callable;

import common.Array;
import common.CmdLine;
import common.Files;
import common.HashVec;
import common.Logger;
import common.WorkerHive;
import common.ext;
import cnv.analysis.FilterCalls;
import filesys.CNVariant;
import filesys.LocusSet;

public class PlinkCNV {

	private static void run(String cnvFile, String sampFile) {
		final String dir = ext.parseDirectoryOfFile(cnvFile);
		String out = dir + "decentCalls_centromeresBroken.cnv";
		final Logger log = new Logger(dir + "cnv.log");
		if (!Files.exists(out)) {
			FilterCalls.fromParameters(dir + "filterCNVs.crf", new Logger());
		}
		String sampCNVs = ext.addToRoot(out, ".sampsDesired");
		HashSet<String> sampSet = HashVec.loadFileToHashSet(sampFile, false);
		if (!Files.exists(sampCNVs)) {
			LocusSet<CNVariant> cnvs = CNVariant.loadLocSet(out, log);

			try {
				int numTotal = 0;
				int numSamps = 0;
				PrintWriter writer = new PrintWriter(new FileWriter(sampCNVs));
				writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER));
				for (int i = 0; i < cnvs.getLoci().length; i++) {
					numTotal++;
					if (sampSet.contains(cnvs.getLoci()[i].getIndividualID())) {
						writer.println(cnvs.getLoci()[i].toPlinkFormat());
						numSamps++;
					}
				}
				writer.close();
				log.reportTimeInfo(numTotal + " - > " + numSamps);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		String[][] phenos = HashVec.loadFileToStringMatrix(dir + "pheno.dat", false, null, false);
		for (int i = 2; i < phenos[0].length; i++) {
			final String pheno = phenos[0][i];
			final String opDir = dir + phenos[0][i] + "/";
			new File(opDir).mkdirs();
			String phenoCNV = opDir + pheno + ".cnv";
			if (!Files.exists(phenoCNV)) {
				Files.copyFile(sampCNVs, phenoCNV);
			}
			ArrayList<String> fam = new ArrayList<String>();
			for (int j = 1; j < phenos.length; j++) {
				if (sampSet.contains(phenos[j][0])) {
					fam.add(phenos[j][0] + "\t" + phenos[j][0] + "\t" + 0 + "\t" + 0 + "\t" + phenos[j][1] + "\t" + phenos[j][i]);
				}
			}
			Files.writeArrayList(fam, opDir + pheno + ".fam");

			CmdLine.run(dir + "plink --cnv-list " + pheno + ".cnv --cnv-make-map --out " + pheno, opDir);

			Callable<Boolean> c1 = new Callable<Boolean>() {

				@Override
				public Boolean call() throws Exception {
					ArrayList<String> cmd = new ArrayList<String>();
					cmd.add(dir + "plink");
					cmd.add("--cfile");
					cmd.add(pheno);
					cmd.add("--cnv-indiv-perm");
					cmd.add("--mperm");
					cmd.add("10000");
					cmd.add("--out");
					cmd.add(pheno);
					return CmdLine.runCommandWithFileChecks(Array.toStringArray(cmd), opDir, null, new String[] { opDir + "VTPEDX.cnv.summary.mperm" }, true, false, false, log);
				}
			};
			Callable<Boolean> c2 = new Callable<Boolean>() {

				@Override
				public Boolean call() throws Exception {
					ArrayList<String> cmd = new ArrayList<String>();
					cmd.add(dir + "plink");
					cmd.add("--cfile");
					cmd.add(pheno);
					cmd.add("--mperm");
					cmd.add("10000");
					cmd.add("--out");
					cmd.add(pheno + "_position");
					return CmdLine.runCommandWithFileChecks(Array.toStringArray(cmd), opDir, null, new String[] { opDir + "VTPEDX.cnv.summary.mperm" }, true, false, false, log);

				}
			};
			Callable<Boolean> c3 = new Callable<Boolean>() {

				@Override
				public Boolean call() throws Exception {
					ArrayList<String> cmd = new ArrayList<String>();
					cmd.add(dir + "plink");
					cmd.add("--cfile");
					cmd.add(pheno);
					cmd.add("--mperm");
					cmd.add("10000");

					cmd.add("--cnv-test-window");
					cmd.add("200");
					cmd.add("--out");
					cmd.add(pheno + "_window");
					return CmdLine.runCommandWithFileChecks(Array.toStringArray(cmd), opDir, null, null, true, false, false, log);
				}
			};
			WorkerHive<Boolean> hive = new WorkerHive<Boolean>(3, 10, log);
			hive.addCallable(c1);
			hive.addCallable(c2);
			hive.addCallable(c3);
			hive.execute(true);
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String dir = "C:/data/ARIC/shadowCNVs/";
		String cnvFile = dir + "penncnv.cnv";
		String sampFile = dir + "whites.txt";

		String usage = "\n" +
				"one.JL.ARICCNV requires 0-1 arguments\n" +
				"   (1) cnvs (i.e. cnvs=" + cnvFile + " (default))\n" +
				"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				cnvFile = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			run(cnvFile, sampFile);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
