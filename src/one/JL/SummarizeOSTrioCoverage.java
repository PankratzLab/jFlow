package one.JL;

import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Set;

import javax.jms.IllegalStateException;

import seq.manage.VCFOps.VcfPopulation;
import seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import common.Files;
import common.Logger;
import common.ext;

public class SummarizeOSTrioCoverage {
	private static final int[] covTargets = new int[] { 10, 20, 30, 40 };

	// /home/spectorl/lanej0/OS_seq/bamQC
	public static void summarize(String vpopFile, String indir, String outDir,
			Logger log) throws IllegalStateException {
		VcfPopulation vpop = VcfPopulation.load(vpopFile, POPULATION_TYPE.ANY,
				log);
		vpop.report();
		String[] files = Files.listFullPaths(indir,
				".sorted.dedup.realigned.recal.txt", false);
		log.reportTimeInfo("Found " + files.length + " files to summarize");
		Hashtable<String, String> map = new Hashtable<String, String>();
		for (int i = 0; i < files.length; i++) {
			String key = ext.removeDirectoryInfo(files[i]).split("_")[0];
			map.put(key, files[i]);
			System.out.println(key);

		}
		ArrayList<String> famsNotFound = new ArrayList<String>();
		for (String fam : vpop.getSubPop().keySet()) {
			// mo,
			// fa
			Set<String> famInds = vpop.getSubPop().get(fam);
			String off = null;
			String offFile = null;
			String mo = null;
			String moFile = null;
			String fa = null;
			String faFile = null;
			boolean useFam = false;
			for (String ind : famInds) {
				ind = ind.replaceAll(".variant.*", "");
				if (ind.endsWith("C")) {
					if (off != null) {
						throw new IllegalStateException("Too many kids");
					}
					off = ind;
					if (map.containsKey(off)) {
						offFile = map.get(off);
						useFam = true;
					} else {
						log.reportTimeInfo("No file for ind " + ind);
					}
				}

				else if (ind.endsWith("D")) {
					if (fa != null) {
						throw new IllegalStateException("Too many pops");
					}
					fa = ind;
					if (map.containsKey(fa)) {
						faFile = map.get(fa);
					}
				} else if (ind.endsWith("M")) {
					if (mo != null) {
						throw new IllegalStateException("Too many mom");
					}
					mo = ind;
					if (map.containsKey(mo)) {
						moFile = map.get(mo);
					}
				} else {
					throw new IllegalStateException("Unknown");

				}
			}
			if (!useFam) {
				famsNotFound.add(fam);
			} else {

				String[] tags = new String[] { off, mo, fa };
				String coverageTrio = outDir + fam + ".covSummary.txt";
				int col = ext.indexOfStr("averageCoverage",
						Files.getHeaderOfFile(offFile, log));
				if (!Files.exists(coverageTrio)) {
					Files.paste(new String[] { offFile, moFile, faFile },
							coverageTrio, new int[] { col }, 0, tags, null, log);
				}

			}

		}

	}

	public static void main(String[] args) {
		String vpop = "/panfs/roc/groups/12/spectorl/lanej0/OS_seq/trioCoverage/dn.vpop";
		String inDir = "/panfs/roc/groups/12/spectorl/lanej0/OS_seq/bamQC/";
		String outDir = "/panfs/roc/groups/12/spectorl/lanej0/OS_seq/trioCoverage/";
		Logger log = new Logger(outDir + "tc.log");
		try {
			summarize(vpop, inDir, outDir, log);
		} catch (IllegalStateException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
