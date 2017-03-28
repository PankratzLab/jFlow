package org.genvisis.one.JL;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;

import org.genvisis.cnv.plots.QQPlot;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.gwas.MatchSamples;

import com.google.common.primitives.Doubles;

public class ControlSelection {
	private static int SAMPLE_SIZE = 550;
	// private static int SAMPLE_SIZE = 200;

	public static void main(String[] args) throws IOException {
		// plink2 --bfile ../plink --keep keeps.txt --make-bed --out ./plink
		// plink2 --bfile ../plink --extract range
		// /Volumes/Beta/data/controlSelection/plink/meregedCaptureRegions.txt
		// --keep cushingPheno.txt --make-bed --out plink
		// ./plink
		// QQPlot qqPlot1 = QQPlot.loadPvals(
		// new String[] {
		// "/Volumes/Beta/data/controlSelection/withCushing/antimatch/CUSHING_AND_OTHERS_v_ARIC/CUSHING_AND_OTHERS_v_ARICantioptimal.p.txt",
		// "/Volumes/Beta/data/controlSelection/withCushing/randomMatch/CUSHING_AND_OTHERS_v_ARIC/CUSHING_AND_OTHERS_v_ARICrandom.p.txt",
		// "/Volumes/Beta/data/controlSelection/withCushing/match/CUSHING_AND_OTHERS_v_ARIC/CUSHING_AND_OTHERS_v_ARICoptimal.p.txt"
		// },
		// "Hi", false, true, false, -1, false, Float.MAX_VALUE, new Logger());
		// qqPlot1.screenCap("/Volumes/Beta/data/controlSelection/withCushing/"
		// + "test.png");
		// System.exit(1);
		// //whites
		String dir = "/Volumes/Beta/data/controlSelection/plink/";
		String plink1 = "/Volumes/Beta/data/controlSelection/plinkARIC/quality_control/further_analysis_QC/";
		String plink2 = "/Volumes/Beta/data/controlSelection/plinkCushing/quality_control/further_analysis_QC/";
		String plinkNew = dir + "plinkMerge/";
		String root = "/Volumes/Beta/data/controlSelection/ARIC_Only_Special/";
		ArrayList<String> lambdaReport = new ArrayList<>();
		lambdaReport.add("Comparison\tmethod\tnumSamples\tlambda\tphenoFile");
		for (int i = 500; i < SAMPLE_SIZE; i += 50) {

			// String phenoFile =
			// "/Volumes/Beta/data/controlSelection/plink/centerPheno2.txt";

			String[] phenoFiles = Files.listFullPaths("/Volumes/Beta/data/controlSelection/centerPhenos/", ".txt",
					false);

			for (String phenoFile : phenoFiles) {
				String pheno = ext.rootOf(phenoFile);
				String rootOut = root + pheno + "_select" + i + "/";
				new File(plinkNew).mkdirs();
				new File(rootOut).mkdirs();
				String start = pheno.charAt(0) + "";
				String mds = "/Volumes/Beta/data/controlSelection/plinkARIC/quality_control/ancestry/unrelated/quality_control/genome/mds20.mds";
				String keepSamples = "/Volumes/Beta/data/controlSelection/plinkARIC/quality_control/genome/plink.genome_keep.dat";

				String plinkQc = "/Volumes/Beta/data/controlSelection/plinkARIC/quality_control/further_analysis_QC/plink_QCd";

				//
				Hashtable<String, String> have = HashVec.loadFileToHashString(mds, new int[] { 0 },
						new int[] { 0, 1, 2, 3, 4, 5, 6, 7 }, false, "\t", false, false, true);
				Hashtable<String, String> keeps = HashVec.loadFileToHashString(keepSamples, new int[] { 0 },
						new int[] { 0, 1, 2, 3, 4, 5, 6, 7 }, false, "\t", false, false, true);

				Hashtable<String, String> hash = HashVec.loadFileToHashString(phenoFile, new int[] { 0 },
						new int[] { 0, 1, 2, 3, 4, 5, 6, 7 }, false, "\t", false, false, true);

				HashMap<String, ArrayList<String>> centers = new HashMap<>();
				for (String sample : hash.keySet()) {
					String center = hash.get(sample).split("\t")[7];
					if (!center.equals("-1")) {
						if (!centers.containsKey(center)) {
							centers.put(center, new ArrayList<>());
						}
						if (have.containsKey(sample) && keeps.containsKey(sample)) {
							centers.get(center).add(hash.get(sample));
						}
					}
				}

				int min = i + 1;
				for (String cent : centers.keySet()) {
					if (centers.get(cent).size() < min) {
						System.out.println(min);
						min = centers.get(cent).size();
					}
				}
				ArrayList<String> qqFiles = new ArrayList<>();
				HashSet<String> comps = new HashSet<>();
				for (String cent1 : centers.keySet()) {
					if (cent1.startsWith(start)) {
						comps.add(cent1);

						ArrayList<String> vAllBarn = new ArrayList<>();
						for (String cent2 : centers.keySet()) {
							if (!cent1.equals(cent2)) {
								vAllBarn.addAll(centers.get(cent2));
							}
						}
						runit(lambdaReport, i, rootOut, mds, plinkQc, hash, centers, min, qqFiles, cent1,
								cent1 + "_v_ALL", vAllBarn, pheno);
						if (i < 1100) {
							for (String cent2 : centers.keySet()) {
								String centComp = cent1 + "_v_" + cent2;

								if (!cent1.equals(cent2) && !comps.contains(centComp) && !cent2.startsWith(start)) {
									comps.add(centComp);
									comps.add(cent2 + "_v_" + cent1);

									ArrayList<String> barns = new ArrayList<>();
									barns.addAll(centers.get(cent2));
									runit(lambdaReport, i, rootOut, mds, plinkQc, hash, centers, min, qqFiles, cent1,
											centComp, barns, pheno);
								}

							}
						}

					}
				}

				System.out.println(ArrayUtils.toStr(qqFiles));
				QQPlot qqPlot = QQPlot.loadPvals(ArrayUtils.toStringArray(qqFiles), "Hi", false, true, false, -1, false,
						Float.MAX_VALUE, new Logger());
				qqPlot.screenCap(rootOut + "test_pheno" + pheno + "_" + i + "samples.png");
				Files.copyFileUsingFileChannels(rootOut + "test_pheno" + pheno + "_" + i + "samples.png",
						root + "test_pheno" + pheno + "_" + i + "samples.png", new Logger());
			}
		}
		Files.writeIterable(lambdaReport, root + "lambdaReport.txt");

	}

	private static void runit(ArrayList<String> lambdaReport, int i, String rootOut, String mds, String plinkQc,
			Hashtable<String, String> hash, HashMap<String, ArrayList<String>> centers, int min,
			ArrayList<String> qqFiles, String cent1, String centComp, ArrayList<String> barns, String pheno) {
		String matchDir = rootOut + "match/" + centComp + "/";

		String[] run = new String[] { "C1", "C2" };
		double[] w = new double[] { 1, 1 };
		String antiOptimalDir = rootOut + "antimatch/" + centComp + "/";
		new File(antiOptimalDir).mkdirs();
		Files.copyFileUsingFileChannels(mds, antiOptimalDir + "mds20.mds", new Logger());

		Files.writeIterable(centers.get(cent1).subList(0, min - 1), antiOptimalDir + centComp + ".anchors.txt");
		Files.writeIterable(barns, antiOptimalDir + centComp + ".barns.txt");
		System.out.println("RUNNING match1");
		String matchFileanti = MatchSamples.matchMaker(antiOptimalDir, centComp + ".anchors.txt",
				centComp + ".barns.txt", ext.removeDirectoryInfo("mds20.mds"), run, w, false);
		matchFileanti = MatchSamples.normalizeDistances(antiOptimalDir, matchFileanti, 0, 100);
		System.out.println("RUNNING match3");

		String pairsAnti = antiOptimalDir + MatchSamples.matchPairs(antiOptimalDir, matchFileanti, false, true);
		System.out.println("RUNNING match4");

		String[] casesAnti = HashVec.loadFileToStringArray(pairsAnti, true, new int[] { 0 }, true);
		String[] controlsAnti = HashVec.loadFileToStringArray(pairsAnti, true, new int[] { 1 }, true);
		Files.copyFileUsingFileChannels(mds, antiOptimalDir + "mds20.mds", new Logger());
		run(plinkQc, mds, hash, qqFiles, centComp, antiOptimalDir, "antioptimal", casesAnti, controlsAnti, lambdaReport,
				i, pheno);

		new File(matchDir).mkdirs();
		Files.copyFileUsingFileChannels(mds, matchDir + "mds20.mds", new Logger());

		Files.writeIterable(centers.get(cent1).subList(0, min - 1), matchDir + centComp + ".anchors.txt");
		Files.writeIterable(barns, matchDir + centComp + ".barns.txt");
		System.out.println("RUNNING match1");
		String matchFile = MatchSamples.matchMaker(matchDir, centComp + ".anchors.txt", centComp + ".barns.txt",
				ext.removeDirectoryInfo("mds20.mds"), run, w, false);
		matchFile = MatchSamples.normalizeDistances(matchDir, matchFile, 0, 100);
		System.out.println("RUNNING match3");

		String pairs = matchDir + MatchSamples.matchPairs(matchDir, matchFile, true);
		System.out.println("RUNNING match4");

		System.out.println(matchFile);
		String[] cases = HashVec.loadFileToStringArray(pairs, true, new int[] { 0 }, true);
		String[] controls = HashVec.loadFileToStringArray(pairs, true, new int[] { 1 }, true);

		run(plinkQc, mds, hash, qqFiles, centComp, matchDir, "optimal", cases, controls, lambdaReport, i, pheno);

		String matchDirRandom = rootOut + "randomMatch/" + centComp + "/";
		new File(matchDirRandom).mkdirs();
		Files.copyFileUsingFileChannels(mds, matchDirRandom + "mds20.mds", new Logger());
		Collections.shuffle(barns);
		List<String> sub = barns.subList(0, Math.min(min, cases.length) - 1);
		Collections.sort(sub);
		run(plinkQc, mds, hash, qqFiles, centComp, matchDirRandom, "random", cases, ArrayUtils.toStringArray(sub),
				lambdaReport, i, pheno);
	}

	private static void run(String plinkQc, String mds, Hashtable<String, String> hash, ArrayList<String> qqFiles,
			String cent, String matchDir, String tag, String[] cases, String[] controls, ArrayList<String> lambdaReport,
			int i, String pheno) {

		ArrayList<String> newFam = new ArrayList<>();
		HashSet<String> cas = new HashSet<>();
		HashSet<String> cont = new HashSet<>();

		for (String sample : cases) {
			String[] c = hash.get(sample.split("\t")[0]).split("\t");
			c[5] = "2";
			newFam.add(ArrayUtils.toStr(ArrayUtils.subArray(c, 0, 6)));
			cas.add(sample.split("\t")[0]);
		}
		for (String sample : controls) {
			// System.out.println(sample);
			String[] c = hash.get(sample.split("\t")[0]).split("\t");
			c[5] = "1";
			newFam.add(ArrayUtils.toStr(ArrayUtils.subArray(c, 0, 6)));
			cont.add(sample.split("\t")[0]);
		}
		Files.writeIterable(newFam, matchDir + "keeps.txt");
		double maf = 40.0 / (controls.length + cases.length);
		// plink2 --bfile ../plink --keep keeps.txt --make-bed --out ./plink

		if (!Files.exists(matchDir + "plink.bed")) {
			ArrayList<String> command = new ArrayList<>();
			command.add("plink2");
			command.add("--bfile");
			command.add(plinkQc);
			command.add("--keep");
			command.add(matchDir + "keeps.txt");
			command.add("--make-bed");

			command.add("--out");

			command.add(matchDir + "plink");
			CmdLine.run(command, matchDir, null, null, new Logger(), false);
		}
		String[] ins = HashVec.loadFileToStringArray(matchDir + "plink.fam", false, null, false);
		// System.out.println(ArrayUtils.toStr(ins));
		newFam = new ArrayList<>();
		HashSet<String> have = new HashSet<>();
		for (String in : ins) {
			have.add(in);
		}
		for (String sample : ins) {
			// System.out.println(sample.split("\\s+")[0]);
			if (cas.contains(sample.split("\\s+")[0])) {
				String[] c = hash.get(sample.split("\\s+")[0]).split("\t");
				c[5] = "2";
				newFam.add(ArrayUtils.toStr(ArrayUtils.subArray(c, 0, 6)));
			}
			if (cont.contains(sample.split("\\s+")[0])) {
				String[] c = hash.get(sample.split("\\s+")[0]).split("\t");
				c[5] = "1";
				newFam.add(ArrayUtils.toStr(ArrayUtils.subArray(c, 0, 6)));
			}

		}

		// for (String sample : controls) {
		// // System.out.println(sample);
		// if (have.contains(sample)) {
		// String[] c = hash.get(sample.split("\t")[0]).split("\t");
		// c[5] = "1";
		// newFam.add(ArrayUtils.toStr(ArrayUtils.subArray(c, 0, 6)));
		// }
		// }
		Collections.sort(newFam);
		Files.writeIterable(newFam, matchDir + "plink.fam");

		if (!Files.exists(matchDir + "plink.assoc.logistic.tabs")) {

			ArrayList<String> plink = new ArrayList<String>();
			plink.add("plink2");

			plink.add("--logistic");
			// plink.add("mperm=100000");

			// plink.add("--adjust");
			// plink.add("mperm=100000");
			// plink.add("--all-pheno");
			// if (perm) {
			// plink.add("perm");
			// }
			plink.add("--maf");
			plink.add(maf + "");
			plink.add("--bfile");
			plink.add(matchDir + "plink");
			plink.add("--geno");
			plink.add(".01");
			plink.add("--mind");
			plink.add(".05");
			// plink.add("--covar-name");
			// plink.add(covars);
			// plink.add("--covar");
			// plink.add(inputDb);
			// plink.add("--pheno");
			// plink.add(processed);
			plink.add("--out");
			plink.add(matchDir + "plink");
			plink.add("--threads");
			plink.add(4 + "");
			plink.add("--covar");
			plink.add(mds);
			plink.add("--covar-name");
			// plink.add("C1,C2,C3,C4,C5,C6,C7,C8,C9,C10");
			plink.add("C1,C2");
			System.out.println(maf);
			CmdLine.run(plink, matchDir, null, null, new Logger(), false);
			// System.exit(1);

		}
		String t1 = "cat " + matchDir + "plink.assoc.logistic |grep \"ADD\\|SNP\"| tr -s ' ' '\t'|cut -f 2- > "
				+ matchDir + "plink.assoc.logistic.tabs";
		String t2 = "cut -f9 " + matchDir + "plink.assoc.logistic.tabs > " + matchDir + cent + tag + ".p.txt";
		org.genvisis.common.Files.write(t1 + "\n" + t2, matchDir + "tab.sh");
		org.genvisis.common.Files.chmod(matchDir + "tab.sh");
		CmdLine.run(matchDir + "tab.sh", matchDir);
		qqFiles.add(matchDir + cent + tag + ".p.txt");
		String[] ps = HashVec.loadFileToStringArray(matchDir + cent + tag + ".p.txt", true, new int[] { 0 }, false);
		ArrayList<Double> p = new ArrayList<>();
		for (String a : ps) {
			try {
				p.add(Double.parseDouble(a));
			} catch (NumberFormatException nfe) {

			}
		}
		lambdaReport.add(cent + "\t" + tag + "\t" + i + "\t" + ArrayUtils.lambda(Doubles.toArray(p)) + "\t" + pheno);
		// ArrayUtils.lambda(pvals[i])
	}

}

// FileUtils.copyDirectoryToDirectory(new
// File(ext.parseDirectoryOfFile(plink1)), new
// File(plinkNew+"ARIC/"));
// FileUtils.copyDirectoryToDirectory(new
// File(ext.parseDirectoryOfFile(plink2)), new
// File(plinkNew+"CUSHING/"));

// cat
// /Volumes/Beta/data/controlSelection/plinkCushing/quality_control/further_analysis_QC/plink_QCd.bim|tr
// -s ' ' '\t'|cut -f
// 2|sort>/Volumes/Beta/data/controlSelection/plink/cushingMarkers.txt
// cat
// /Volumes/Beta/data/controlSelection/plinkARIC/quality_control/further_analysis_QC/plink_QCd.bim|
// tr -s ' ' '\t'|cut -f 2|sort
// >/Volumes/Beta/data/controlSelection/plink/ARICMarkers.txt
// comm -12 cushingMarkers.txt ARICMarkers.txt >intersectMarkers.txt
// Qc.fullGamut("/Volumes/Beta/data/controlSelection/plinkCushing/",
// "plink", false, new Logger());
// MergeDatasets.checkForHomogeneity(null, new String[] { plink2,
// plink1
// }, plinkNew, "ALL(NP)", new Logger());

// System.exit(1);
// cat /Volumes/Beta/data/controlSelection/plink/fake/exclude50.txt
// /Volumes/Beta/data/controlSelection/plink/plinkMerge/lackOfHomogeneity.dat
// |sort |uniq > exclude.txt
// plink2 --bfile ../plink --keep centerPheno.txt --extract
// intersectMarkers.txt --exclude exclude.txt --make-bed --out
// ./plink
//
// plink2 --bfile ../plink --keep centerPheno.txt
// Qc.fullGamut(dir, "plink", false, new Logger());
// Qc.fullGamut("/Volumes/Beta/data/controlSelection/plinkARIC/quality_control/ancestry/unrelated/",
// "plink",
// false, new Logger());
// //
//
// System.exit(1);
// rs67076383