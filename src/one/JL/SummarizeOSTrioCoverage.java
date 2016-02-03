package one.JL;

import java.util.Hashtable;

import seq.manage.VCFOps.VcfPopulation;
import seq.manage.VCFOps.VcfPopulation.POPULATION_TYPE;
import common.Files;
import common.Logger;

public class SummarizeOSTrioCoverage {
	// /home/spectorl/lanej0/OS_seq/bamQC
	public static void summarize(String vpopFile, String indir, String outDir, Logger log) {
		VcfPopulation vpop = VcfPopulation.load(vpopFile, POPULATION_TYPE.ANY, log);

		String[] files = Files.listFullPaths(indir, ".sorted.dedup.realigned.recal.txt", false);
		Hashtable<String, String> map = new Hashtable<String, String>();
		for (int i = 0; i < files.length; i++) {
			map.put(files[i].split("_")[0], files[i]);
		}

	}

	public static void main(String[] args) {
		// String vpop
	}

}
