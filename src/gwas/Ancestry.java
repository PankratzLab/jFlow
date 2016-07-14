package gwas;

import java.io.File;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Vector;

import cnv.analysis.pca.PCImputeRace;
import cnv.filesys.Project;
import common.Array;
import common.CmdLine;
import common.Files;
import common.HashVec;
import common.Logger;
import common.ext;
import filesys.SnpMarkerSet;

public class Ancestry {
	
	public static void checkHomogeneity(String dir, String projectPlinkRoot, String hapMapPlinkRoot, Logger log) {
		String homoDir = dir + "homogeneity/";
		String homoProjDir = homoDir + ext.removeDirectoryInfo(projectPlinkRoot) + "/";
		String homoHapMapDir = homoDir + ext.removeDirectoryInfo(hapMapPlinkRoot) + "/";
		new File(homoProjDir).mkdirs();
		new File(homoHapMapDir).mkdirs();
		CmdLine.runDefaults("plink2 --bfile " + projectPlinkRoot + " --keep " + dir + "whites.txt --hardy", homoProjDir, log);
		CmdLine.runDefaults("plink2 --bfile " + hapMapPlinkRoot + " --keep " + ext.parseDirectoryOfFile(hapMapPlinkRoot) + "CEUFounders.txt --hardy", homoHapMapDir, log);
		
		MergeDatasets.checkForHomogeneity(homoDir, null, null, "UNAFF", log);
	}
	
	public static void parseHomogeneity(String dir) {
		MergeDatasets.parseHomo(dir + "homogeneity/");
	}
	
	public static void mergeHapMap(String dir, String projectPlinkRoot, String hapMapPlinkRoot, String dropMarkersFile, Logger log) {
		if (!Files.exists(dir+"unrelateds.txt")) {
			log.reportError("Error - need a file called unrelateds.txt with FID and IID pairs before we can proceed");
			return;
		}
		if (!Files.exists(dir+"plink.bim_unambiguous.txt")) {
			log.report(ext.getTime() + "]\tGenerating list of unambiguous SNPs");
			new SnpMarkerSet(dir + "plink.bim", true, log).listUnambiguousMarkers(dir+"plink.bim_unambiguous.txt", dropMarkersFile, true);
		}
		
		if (!Files.exists(dir+"unambiguousHapMap.bed")) {
			log.report(ext.getTime() + "]\tExtracting unambiguous SNPs for HapMap founders");
			CmdLine.runDefaults("plink2 --bfile " + hapMapPlinkRoot + " --extract plink.bim_unambiguous.txt --make-bed --out unambiguousHapMap --noweb", dir, log);
		}
		
		if (!Files.exists(dir+"overlap.txt")) {
			log.report(ext.getTime() + "]\tGenerating list of overlapping SNPs");
			Files.writeVector(HashVec.loadFileToVec(dir+"unambiguousHapMap.bim", false, new int[] {1}, false, false), dir+"overlap.txt");
		}
		
		if (Files.exists(dir+"combo.missnp")) {
			new File(dir+"combo.missnp").delete();
		}
		
		if (!Files.exists(dir+"unambiguous.bed")) {
			log.report(ext.getTime() + "]\tExtracting overlapping SNPs for study samples");
			CmdLine.runDefaults("plink2 --bfile plink --extract overlap.txt --make-bed --out unambiguous --noweb", dir, log);
		}
		
		if (!Files.exists(dir+"combo.missnp")) {
			log.report(ext.getTime() + "]\tMerging study data and HapMap data for overlapping SNPs");
			CmdLine.runDefaults("plink --bfile unambiguous --bmerge unambiguousHapMap.bed unambiguousHapMap.bim unambiguousHapMap.fam --out combo --noweb", dir, log);
		}
		
		if (Files.exists(dir+"combo.missnp")) {
			if (Files.exists(dir+"combo.1.missnp")) {
				new File(dir+"combo.1.missnp").delete();
			}
			log.report(ext.getTime() + "]\tChecking for flipped alleles");
			new File(dir+"combo.missnp").renameTo(new File(dir+"combo.1.missnp"));
			CmdLine.runDefaults("plink2 --bfile unambiguous --flip combo.1.missnp --make-bed --out unambiguousFlipped --noweb", dir, log);
			CmdLine.runDefaults("plink --bfile unambiguousFlipped --bmerge unambiguousHapMap.bed unambiguousHapMap.bim unambiguousHapMap.fam --make-bed --out combo --noweb", dir, log);
			if (Files.exists(dir+"combo.missnp")) {
				if (Files.exists(dir+"combo.2.missnp")) {
					new File(dir+"combo.2.missnp").delete();
				}
				log.report(ext.getTime() + "]\tDropping SNPs that cannot be resolved by flipping alleles");
				new File(dir+"combo.missnp").renameTo(new File(dir+"combo.2.missnp"));
				CmdLine.runDefaults("plink2 --bfile unambiguousFlipped --exclude combo.2.missnp --make-bed --out unambiguousFlippedDropped --noweb", dir, log);
				CmdLine.runDefaults("plink2 --bfile unambiguousHapMap --exclude combo.2.missnp --make-bed --out unambiguousDroppedHapMap --noweb", dir, log);
				CmdLine.runDefaults("plink --bfile unambiguousFlippedDropped --bmerge unambiguousDroppedHapMap.bed unambiguousDroppedHapMap.bim unambiguousDroppedHapMap.fam --make-bed --out combo --noweb", dir, log);
			}
		}
		
		if (!Files.exists(dir+"finalSNPs.txt")) {
			log.report(ext.getTime() + "]\tWriting final list of SNPs to use");
			Files.writeVector(HashVec.loadFileToVec(dir+"combo.bim", false, new int[] {1}, false, false), dir+"finalSNPs.txt");
		}
	}
	
	public static void runEigenstrat(String dir) {
		Logger log;
		
		dir = ext.verifyDirFormat(dir);
		log = new Logger(dir+"ancestry.log");

		String unrelatedsDir = dir + "unrelateds/";
		
		if (!Files.exists(unrelatedsDir + "unrelateds.txt")) {
			log.report(ext.getTime() + "]\tGenerating combined unrelateds.txt");
			new File(unrelatedsDir).mkdir();
			Vector<String> unrelateds = HashVec.loadFileToVec(dir+"unambiguousHapMap.fam", false, new int[] {0, 1}, false, false);
			unrelateds.addAll(HashVec.loadFileToVec(dir + "unrelateds.txt", false, false, false));
			Files.writeVector(unrelateds, unrelatedsDir + "unrelateds.txt");
		}
		
		if (!Files.exists(unrelatedsDir + "plink.bed")) {
			log.report(ext.getTime() + "]\tGenerating PLINK files based on combined unrelateds.txt");
			CmdLine.runDefaults("plink2 --bfile ../combo --keep unrelateds.txt --make-bed --noweb", unrelatedsDir, log);
		}
		
		if (!Files.exists(unrelatedsDir + "master")) {
			log.report(ext.getTime() + "]\tCreating Eigenstrat");
			CmdLine.runDefaults("jcp gwas.Eigenstrat source=plink -create", unrelatedsDir, log);
		}
		
		if (!Files.exists(unrelatedsDir + "plink.pca.evec")) {
			log.report(ext.getTime() + "]\tRunning master");
			CmdLine.runDefaults("./master", unrelatedsDir, log);
			CmdLine.runDefaults("jcp gwas.Eigenstrat convert=plink.pca.evec", unrelatedsDir, log);
		}
		
		if (!Files.exists(dir+"convertf.par")) {
			log.report(ext.getTime() + "]\tRunning convertf");
			CmdLine.runDefaults("jcp gwas.Eigenstrat source=combo -create", dir, log);
			CmdLine.runDefaults("convertf -p convertf.par", dir, log);
		}
		
		if (!Files.exists(dir+"combo_fancy_postnormed_eigens.xln")) {
			CmdLine.runDefaults("plink2 --bfile unrelateds/plink --freq --out unrelateds/plink --noweb", dir, log);
			CmdLine.runDefaults("jcp gwas.Eigenstrat source=unrelateds/plink target=combo -parse -eigenFormat", dir, log);
		}
	}
	
	public static void imputeRace(String dir, Project proj) {
		String[][] pcResults = HashVec.loadFileToStringMatrix(dir+"combo_fancy_postnormed_eigens.xln", true, new int[] {0, 1, 2, 3}, false);
		
//		SampleData sampleData = proj.getSampleData(0, false);
//		int hapMapClass = ext.indexOfStr("HapMap", sampleData.getClasses());
		String sd = proj.SAMPLE_DATA_FILENAME.getValue();
		String[] sdHeader = Files.getHeaderOfFile(sd, proj.getLog());
		Hashtable<String, Hashtable<String, String>> hapmaps = HashVec.loadFileToHashHash(sd, ext.indexOfStr("FID", sdHeader, false, true), ext.indexOfStr("IID", sdHeader, false, true), ext.indexOfStr("Class=HapMap;1=CEU;2=YRI;3=CHB;4=JPT", sdHeader, false, true), true);
		
		String[] fidiids = new String[pcResults.length];
		double[] pc1 = new double[pcResults.length];
		double[] pc2 = new double[pcResults.length];
		ArrayList<Integer> europeans = new ArrayList<Integer>();
		ArrayList<Integer> africans = new ArrayList<Integer>();
		ArrayList<Integer> asians = new ArrayList<Integer>();
		for (int i = 0; i < pcResults.length; i++) {
			fidiids[i] = pcResults[i][0] + "\t" + pcResults[i][1];
			try {
				pc1[i] = Double.parseDouble(pcResults[i][2]);
			} catch (NumberFormatException nfe) {
				pc1[i] = Double.NaN;
			}
			try {
				pc2[i] = Double.parseDouble(pcResults[i][3]);
			} catch (NumberFormatException nfe) {
				pc2[i] = Double.NaN;
			}
			
			try {
				int race = Integer.parseInt(hapmaps.get(pcResults[i][0]).get(pcResults[i][1]));
				switch(race) {
					case 1:
						europeans.add(i);
						break;
					case 2:
						africans.add(i);
						break;
					case 3:
						asians.add(i);
					case 4:
						asians.add(i);
						break;
					default:
						break;
				}
			} catch (NumberFormatException nfe) { }
			
			
		}
		PCImputeRace pcir = new PCImputeRace(proj, fidiids, pc1, pc2, Array.toIntArray(europeans), Array.toIntArray(africans), Array.toIntArray(asians), proj.getLog());
		pcir.correctPCsToRace(dir+"raceImputations.mds");
		
		PCImputeRace.freqsByRace(dir+"raceImputations.mds", dir + "plink", dir + "freqsByRace.xln", proj.getLog());
		
	}
	

	public static void main(String[] args) {

		int numArgs = args.length;
		String dir = "./";
		String hapMapPlinkRoot = "/home/pankrat2/shared/bin/HapMap/unambiguousHapMapFounders";
		boolean checkHomo = false;
		boolean run = false;
		boolean imputeRace = false;
		Project proj = null;
		String logfile = null;
		Logger log;
		
		String usage = "\n" +
				"gwas.Ancestry requires 3 arguments\n" +
				"   (1) Run directory with plink.* files and unrelateds.txt (i.e. dir=" + dir + " (default))\n" + 
				"   (2) PLINK root of Unambiguous HapMap Founders (i.e. hapMapPlinkRoot=" + hapMapPlinkRoot + " (default))\n" +
				"   (3) Logfile (i.e. logfile=" + "ancestry.log" + " (default))" +
				"  AND" +
				"   (4) Generate PBS script to check Homogeneity (i.e. -checkHomo (not the default))" +
				"  OR" +
				"   (5) Parse homogeneity checks and run Eigenstrat (i.e. -run (not the default))" +
				"  OR" +
				"   (6) Impute race (i.e. -imputeRace (not the default))" +
				"   (7) Project properties file (i.e. proj=example.properties (not the default))\n" + 
				"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("dir=")) {
				dir = ext.parseStringArg(args[i], "./");
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("hapMapPlinkRoot=")) {
				hapMapPlinkRoot = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("-checkHomo")) {
				checkHomo = true;
				numArgs--;
			} else if (args[i].startsWith("-run")) {
				run = true;
				numArgs--;
			} else if (args[i].startsWith("-imputeRace")) {
				imputeRace = true;
				numArgs--;
			} else if (args[i].startsWith("proj")) {
				proj = new Project(ext.parseStringArg(args[i], "./"), false);
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		if (logfile == null) {
			log = new Logger(dir+"ancestry.log");
		} else {
			log = new Logger(logfile);
		}
		try {
			if (checkHomo) {
				checkHomogeneity(dir, dir + "plink", hapMapPlinkRoot, log);
			} else if (run) {
				parseHomogeneity(dir);
				mergeHapMap(dir, dir + "plink", hapMapPlinkRoot, dir + "homogeneity/FisherOrChiSquareDrops.dat", log);
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
