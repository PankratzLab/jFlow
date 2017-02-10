package org.genvisis.gwas;

import java.io.File;
import java.util.Hashtable;
import java.util.Map;
import java.util.Set;
import java.util.Vector;

import org.genvisis.cnv.analysis.pca.PCImputeRace;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.filesys.SnpMarkerSet;

import com.google.common.collect.Sets;



public class Ancestry {

	public static final String DEFAULT_HAPMAP_PLINKROOT =
																											"/home/pankrat2/shared/bin/HapMap/unambiguousHapMapFounders";
	public static final String RACE_IMPUTATIONAS_FILENAME = "raceImputations.mds";
	public static final String RACE_FREQS_FILENAME = "freqsByRace.xln";

	public static void runPipeline(String dir, String putativeWhitesFile, Project proj, Logger log) {
		runPipeline(dir, putativeWhitesFile, null, proj, log);
	}

	public static void runPipeline(	String dir, String putativeWhitesFile, String hapMapPlinkRoot,
																	Project proj, Logger log) {
		if (hapMapPlinkRoot == null) {
			hapMapPlinkRoot = DEFAULT_HAPMAP_PLINKROOT;
		}
		if (!Files.exists(dir + "homogeneity/" + MergeDatasets.CHI_SQUARE_DROPS_FILENAME)
				&& Files.list(dir + "homogeneity/", ".Rout", false).length == 0) {
			log.report("Running homogeneity checks...");
			checkHomogeneity(dir, putativeWhitesFile, dir + "plink", hapMapPlinkRoot, log);
		}
		String homogeneityDrops = parseHomogeneity(dir, log);
		mergeHapMap(dir, dir + "plink", hapMapPlinkRoot, homogeneityDrops, log);
		runEigenstrat(dir);
		imputeRace(dir, proj);
	}

	public static void checkHomogeneity(String dir, String putativeWhitesFile,
																			String projectPlinkRoot, String hapMapPlinkRoot, Logger log) {
		String homoDir = dir + "homogeneity/";
		String homoProjDir = homoDir + ext.removeDirectoryInfo(projectPlinkRoot) + "/";
		String homoHapMapDir = homoDir + ext.removeDirectoryInfo(hapMapPlinkRoot) + "/";
		new File(homoProjDir).mkdirs();
		new File(homoHapMapDir).mkdirs();
		CmdLine.runDefaults("plink2 --bfile "	+ projectPlinkRoot + " --keep " + putativeWhitesFile
												+ " --hardy", homoProjDir, log);
		CmdLine.runDefaults("plink2 --bfile "	+ hapMapPlinkRoot + " --keep "
												+ ext.parseDirectoryOfFile(hapMapPlinkRoot) + "CEUFounders.txt --hardy",
												homoHapMapDir, log);

		MergeDatasets.checkForHomogeneity(homoDir, null, null, "UNAFF", log);
	}

	public static String parseHomogeneity(String dir, Logger log) {
		int rOuts = Files.list(dir + "homogeneity/", ".Rout", false).length;
		if (rOuts == 0) {
			log.report("No Fisher's Exact results found, using Chi Square to choose homogeneous markers");
			return dir + "homogeneity/" + MergeDatasets.CHI_SQUARE_DROPS_FILENAME;
		}
		MergeDatasets.parseHomo(dir + "homogeneity/");
		return dir + "homogeneity/" + MergeDatasets.FISHER_OR_CHI_SQUARE_DROPS_FILENAME;
	}

	public static void mergeHapMap(	String dir, String projectPlinkRoot, String hapMapPlinkRoot,
																	String dropMarkersFile, Logger log) {
		if (!Files.exists(dir + "unrelateds.txt")) {
			log.reportError("Error - need a file called unrelateds.txt with FID and IID pairs before we can proceed");
			return;
		}
		if (!Files.exists(dir + "plink.bim_unambiguous.txt")) {
			log.report(ext.getTime() + "]\tGenerating list of unambiguous SNPs");
			new SnpMarkerSet(dir + "plink.bim", true, log).listUnambiguousMarkers(dir
																																							+ "plink.bim_unambiguous.txt",
																																						dropMarkersFile, true);
		}

		if (!Files.exists(dir + "unambiguousHapMap.bed")) {
			log.report(ext.getTime() + "]\tExtracting unambiguous SNPs for HapMap founders");
			CmdLine.runDefaults("plink2 --bfile "	+ hapMapPlinkRoot
													+ " --extract plink.bim_unambiguous.txt --make-bed --out unambiguousHapMap --noweb",
													dir, log);
		}

		if (!Files.exists(dir + "overlap.txt")) {
			log.report(ext.getTime() + "]\tGenerating list of overlapping SNPs");
			Files.writeIterable(HashVec.loadFileToVec(dir	+ "unambiguousHapMap.bim", false, new int[] {1},
																								false, false),
													dir + "overlap.txt");
		}

		if (Files.exists(dir + "combo.missnp")) {
			new File(dir + "combo.missnp").delete();
		}

		if (!Files.exists(dir + "unambiguous.bed")) {
			log.report(ext.getTime() + "]\tExtracting overlapping SNPs for study samples");
			CmdLine.runDefaults("plink2 --bfile plink --extract overlap.txt --make-bed --out unambiguous --noweb",
													dir, log);
		}

		if (!Files.exists(dir + "combo.missnp")) {
			log.report(ext.getTime() + "]\tMerging study data and HapMap data for overlapping SNPs");
			CmdLine.runDefaults("plink --bfile unambiguous --bmerge unambiguousHapMap.bed unambiguousHapMap.bim unambiguousHapMap.fam --out combo --noweb",
													dir, log);
		}

		if (Files.exists(dir + "combo.missnp")) {
			if (Files.exists(dir + "combo.1.missnp")) {
				new File(dir + "combo.1.missnp").delete();
			}
			log.report(ext.getTime() + "]\tChecking for flipped alleles");
			new File(dir + "combo.missnp").renameTo(new File(dir + "combo.1.missnp"));
			CmdLine.runDefaults("plink2 --bfile unambiguous --flip combo.1.missnp --make-bed --out unambiguousFlipped --noweb",
													dir, log);
			CmdLine.runDefaults("plink --bfile unambiguousFlipped --bmerge unambiguousHapMap.bed unambiguousHapMap.bim unambiguousHapMap.fam --make-bed --out combo --noweb",
													dir, log);
			if (Files.exists(dir + "combo.missnp")) {
				if (Files.exists(dir + "combo.2.missnp")) {
					new File(dir + "combo.2.missnp").delete();
				}
				log.report(ext.getTime() + "]\tDropping SNPs that cannot be resolved by flipping alleles");
				new File(dir + "combo.missnp").renameTo(new File(dir + "combo.2.missnp"));
				CmdLine.runDefaults("plink2 --bfile unambiguousFlipped --exclude combo.2.missnp --make-bed --out unambiguousFlippedDropped --noweb",
														dir, log);
				CmdLine.runDefaults("plink2 --bfile unambiguousHapMap --exclude combo.2.missnp --make-bed --out unambiguousDroppedHapMap --noweb",
														dir, log);
				CmdLine.runDefaults("plink --bfile unambiguousFlippedDropped --bmerge unambiguousDroppedHapMap.bed unambiguousDroppedHapMap.bim unambiguousDroppedHapMap.fam --make-bed --out combo --noweb",
														dir, log);
			}
		}

		if (!Files.exists(dir + "finalSNPs.txt")) {
			log.report(ext.getTime() + "]\tWriting final list of SNPs to use");
			Files.writeIterable(HashVec.loadFileToVec(dir	+ "combo.bim", false, new int[] {1}, false,
																								false),
													dir + "finalSNPs.txt");
		}
	}

	public static void runEigenstrat(String dir) {
		Logger log;

		dir = ext.verifyDirFormat(dir);
		log = new Logger(dir + "ancestry.log");

		String unrelatedsDir = dir + "unrelateds/";

		if (!Files.exists(unrelatedsDir + "unrelateds.txt")) {
			log.report(ext.getTime() + "]\tGenerating combined unrelateds.txt");
			new File(unrelatedsDir).mkdir();
			Vector<String> unrelateds = HashVec.loadFileToVec(dir	+ "unambiguousHapMap.fam", false,
																												new int[] {0, 1}, false, false);
			unrelateds.addAll(HashVec.loadFileToVec(dir + "unrelateds.txt", false, false, false));
			Files.writeIterable(unrelateds, unrelatedsDir + "unrelateds.txt");
		}

		if (!Files.exists(unrelatedsDir + "plink.bed")) {
			log.report(ext.getTime() + "]\tGenerating PLINK files based on combined unrelateds.txt");
			CmdLine.runDefaults("plink2 --bfile ../combo --keep unrelateds.txt --make-bed --noweb",
													unrelatedsDir, log);
		}

		if (!Files.exists(unrelatedsDir + "master")) {
			log.report(ext.getTime() + "]\tCreating Eigenstrat");
			CmdLine.runDefaults(Files.getRunString()	+ " gwas.Eigenstrat source=plink -create",
													unrelatedsDir, log);
		}

		if (!Files.exists(unrelatedsDir + "plink.pca.evec")) {
			log.report(ext.getTime() + "]\tRunning master");
			CmdLine.runDefaults("./master", unrelatedsDir, log);
			CmdLine.runDefaults(Files.getRunString()	+ " gwas.Eigenstrat convert=plink.pca.evec",
													unrelatedsDir, log);
		}

		if (!Files.exists(dir + "convertf.par")) {
			log.report(ext.getTime() + "]\tRunning convertf");
			CmdLine.runDefaults(Files.getRunString() + " gwas.Eigenstrat source=combo -create", dir, log);
			CmdLine.runDefaults("convertf -p convertf.par", dir, log);
		}

		if (!Files.exists(dir + "combo_fancy_postnormed_eigens.xln")) {
			CmdLine.runDefaults("plink2 --bfile unrelateds/plink --freq --out unrelateds/plink --noweb",
													dir, log);
			CmdLine.runDefaults(Files.getRunString()
														+ " gwas.Eigenstrat source=unrelateds/plink target=combo -parse -eigenFormat",
													dir, log);
		}
	}

	public static void imputeRace(String dir, Project proj) {
		if (!Files.exists(dir + RACE_IMPUTATIONAS_FILENAME)) {
			String[][] pcResults = HashVec.loadFileToStringMatrix(dir
																															+ "combo_fancy_postnormed_eigens.xln",
																														true, new int[] {0, 1, 2, 3}, false);

			String sd = proj.SAMPLE_DATA_FILENAME.getValue();
			String[] sdHeader = Files.getHeaderOfFile(sd, proj.getLog());
			int fidIndex = ext.indexOfStr("FID", sdHeader, false, true);
			int iidIndex = ext.indexOfStr("IID", sdHeader, false, true);
			int hapIndex = ext.indexOfStr("Class=HapMap;1=CEU;2=YRI;3=CHB;4=JPT", sdHeader, false, true);
			if (hapIndex < 0) {
				proj.getLog().reportError("Cannot impute: no HapMap column in sample data.");
				proj.getLog()
						.reportError("Please create a column with header: \"Class=HapMap;1=CEU;2=YRI;3=CHB;4=JPT\"");
				return;
			}
			Hashtable<String, Hashtable<String, String>> hapmaps = HashVec.loadFileToHashHash(sd,
																																												fidIndex,
																																												iidIndex,
																																												hapIndex,
																																												true);

			Set<PCImputeRace.Sample> samples = Sets.newHashSet();
			Set<PCImputeRace.Sample> europeans = Sets.newHashSet();
			Set<PCImputeRace.Sample> africans = Sets.newHashSet();
			Set<PCImputeRace.Sample> asians = Sets.newHashSet();

			for (int i = 0; i < pcResults.length; i++) {


				Map<String, String> iidTable = hapmaps.get(pcResults[i][0]);
				if (iidTable == null || !iidTable.containsKey(pcResults[i][1])) {
					continue;
				}
				String fid = pcResults[i][0];
				String iid = pcResults[i][1];
				double pc1;
				double pc2;
				try {
					pc1 = Double.parseDouble(pcResults[i][2]);
				} catch (NumberFormatException nfe) {
					pc1 = Double.NaN;
				}
				try {
					pc2 = Double.parseDouble(pcResults[i][3]);
				} catch (NumberFormatException nfe) {
					pc2 = Double.NaN;
				}

				PCImputeRace.Sample sample = new PCImputeRace.Sample(fid, iid, pc1, pc2);
				samples.add(sample);
				try {
					int race = Integer.parseInt(iidTable.get(pcResults[i][1]));
					switch (race) {
						case 1:
							europeans.add(sample);
							break;
						case 2:
							africans.add(sample);
							break;
						case 3:
						case 4:
							asians.add(sample);
							break;
						default:
							break;
					}
				} catch (NumberFormatException nfe) {
					// For non hap-map samples, don't add to a race set
				}
			}

			PCImputeRace pcir =
												new PCImputeRace(proj, samples, europeans, africans, asians, proj.getLog());
			pcir.correctPCsToRace(dir + RACE_IMPUTATIONAS_FILENAME);
		} else {
			proj.getLog().reportTimeWarning("Skipping imputation - output already exists: "
																			+ (dir + RACE_IMPUTATIONAS_FILENAME));
		}

		if (!Files.exists(dir + RACE_FREQS_FILENAME)) {
			PCImputeRace.freqsByRace(dir	+ RACE_IMPUTATIONAS_FILENAME, dir + "plink",
																dir + RACE_FREQS_FILENAME, proj.getLog());
		}
		{
			proj.getLog().reportTimeWarning("Skipping race freq calculation - output already exists: "
																			+ (dir + RACE_FREQS_FILENAME));
		}

	}


	public static void main(String[] args) {

		int numArgs = args.length;
		String dir = "./";
		String hapMapPlinkRoot = DEFAULT_HAPMAP_PLINKROOT;
		String putativeWhites = null;
		boolean checkHomo = false;
		boolean run = false;
		boolean runPipeline = false;
		boolean imputeRace = false;
		Project proj = null;
		String logfile = null;
		Logger log;

		String usage = "\n"	+ "gwas.Ancestry requires 3+ arguments\n"
										+ "   (1) Run directory with plink.* files and unrelateds.txt (i.e. dir=" + dir
										+ " (default))\n"
										+ "   (2) PLINK root of Unambiguous HapMap Founders (i.e. hapMapPlinkRoot="
										+ hapMapPlinkRoot + " (default))\n" + "   (3) Logfile (i.e. log="
										+ "ancestry.log" + " (default))\n" + "  AND\n"
										+ "   (4) Run full pipeline (i.e. -runPipeline (not the default, requires arguments for each step))\n"
										+ "  OR\n"
										+ "   (5) Check Homogeneity using Chi-Square (Generates PBS script to run Fisher's exact, if desired) (i.e. -checkHomo (not the default))\n"
										+ "   (6) File of FID/IID pairs of putative whites to use for finding homogenous markers by comparison to CEU (i.e. putativeWhites=whites.txt (not the default))\n"
										+ "  OR\n"
										+ "   (7) Parse homogeneity checks and run Eigenstrat (i.e. -run (not the default))\n"
										+ "  OR\n" + "   (8) Impute race (i.e. -imputeRace (not the default))\n"
										+ "   (9) Project properties file (i.e. proj=example.properties (not the default))\n"
										+ "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("dir=")) {
				dir = ext.parseStringArg(arg, "./");
				numArgs--;
			} else if (arg.startsWith("log=")) {
				logfile = ext.parseStringArg(arg, null);
				numArgs--;
			} else if (arg.startsWith("hapMapPlinkRoot=")) {
				hapMapPlinkRoot = ext.parseStringArg(arg, null);
				numArgs--;
			} else if (arg.startsWith("-runPipeline")) {
				runPipeline = true;
				numArgs--;
			} else if (arg.startsWith("-checkHomo")) {
				checkHomo = true;
				numArgs--;
			} else if (arg.startsWith("putativeWhites=")) {
				putativeWhites = ext.parseStringArg(arg, null);
				numArgs--;
			} else if (arg.startsWith("-run")) {
				run = true;
				numArgs--;
			} else if (arg.startsWith("-imputeRace")) {
				imputeRace = true;
				numArgs--;
			} else if (arg.startsWith("proj")) {
				proj = new Project(ext.parseStringArg(arg, "./"), false);
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + arg);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		if (logfile == null) {
			if (proj == null) {
				log = new Logger(dir + "ancestry.log");
			} else {
				log = proj.getLog();
			}
		} else {
			log = new Logger(logfile);
		}
		try {
			if (runPipeline) {
				runPipeline(dir, putativeWhites, hapMapPlinkRoot, proj, log);
			} else if (checkHomo) {
				checkHomogeneity(dir, putativeWhites, dir + "plink", hapMapPlinkRoot, log);
			} else if (run) {
				String homogeneityDrops = parseHomogeneity(dir, log);
				mergeHapMap(dir, dir + "plink", hapMapPlinkRoot, homogeneityDrops, log);
				runEigenstrat(dir);
			} else if (imputeRace) {
				if (proj == null) {
					System.err.println("Project required for race imputation");
					System.err.println(usage);
					System.exit(1);
				}
				imputeRace(dir, proj);
			} else {
				System.err.println(usage);
				System.exit(1);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

}
