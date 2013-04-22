// still need to code useResiduals for conPlink
// still need to code useResiduals for conDosge
// allow for an annotation file in the root to be passed
//     if equalsignorecase(chr) and position are present, use that instead
//     can contain Discovery if Replication or Replication if Discovery as well as Meta-results

package gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;

import stats.ProbDist;

import common.Array;
import common.CmdLine;
import common.CountHash;
import common.Files;
import common.HashVec;
import common.Logger;
import common.Matrix;
import common.Sort;
import common.Zip;
import common.ext;

import filesys.DosageData;
import filesys.SerialHash;
import filesys.SnpMarkerSet;

public class Conditional {
//	public static void addCountsAsCovariate(String baseDir, String originalCovariatesFile, String newCovariatesFile, String allPossibleSNPs, Logger log) {
//		addCountsAsCovariate(baseDir, "../plink", originalCovariatesFile, newCovariatesFile, allPossibleSNPs, log);
//	}
	
	public static void addCountsAsCovariate(String baseDir, String originalCovariatesFile, String newCovariatesFile, String allPossibleSNPs, Logger log) {
		String subDir;
		if (allPossibleSNPs == null) {
			subDir = "allMarkers/";
		} else {
			subDir = ext.rootOf(allPossibleSNPs, false)+"/";
		}
		addCountsAsCovariate(baseDir, subDir, "../plink", originalCovariatesFile, newCovariatesFile, allPossibleSNPs, log);
	}

	public static void addCountsAsCovariate(String baseDir, String subDir, String plinkFileDirPlusRoot, String originalCovariatesFile, String newCovariatesFile, String allPossibleSNPs, Logger log) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		Hashtable<String,String> hash;
		Vector<String> v;
		String[] markerNames;
		SnpMarkerSet markerSet;
		int[] indices;
		int numFixedCovariates;
		String dir;
		int countMissingCovariate;
		
		dir = baseDir + subDir;
		if (!new File(dir).exists()) {
			new File(dir).mkdirs();
		}

		if (!new File(plinkFileDirPlusRoot+".bed").exists() || !new File(plinkFileDirPlusRoot+".bim").exists() || !new File(plinkFileDirPlusRoot+".fam").exists()) {
			log.reportError("Error - a full complement of binary PLINK files are required to perform this algorithm");
			return;
		}
		if (allPossibleSNPs == null) {
			log.report("All markers in the dataset will be added to the covars file, since no subset was pre-defined");
			CmdLine.run("plink --bfile "+plinkFileDirPlusRoot+" --make-bed", dir);
		} else {
			markerNames = HashVec.loadFileToStringArray(baseDir + allPossibleSNPs, false, new int[] {0}, true);
			log.report("Curently using file '" + allPossibleSNPs + "', which has "+markerNames.length+" independent elements in it");
//			CmdLine.run("plink --bfile "+plinkFiles+" --extract ../"+(baseDir.equals("")?"":"../")+ext.removeDirectoryInfo(allPossibleSNPs)+" --make-bed", dir);
//			CmdLine.run("C:/plink/plink --bfile "+plinkFiles+" --extract ../"+(baseDir.equals("")?"":"../")+ext.removeDirectoryInfo(allPossibleSNPs)+" --make-bed", dir);
			CmdLine.run("C:/plink/plink --bfile " + plinkFileDirPlusRoot + " --extract ../" + ext.removeDirectoryInfo(allPossibleSNPs) + " --make-bed", dir);
		}
		
//		CmdLine.run("plink --bfile plink --recode", dir);
//		CmdLine.run("plink --bfile plink --freq", dir);
		CmdLine.run("C:/plink/plink --bfile plink --recode", dir);
		CmdLine.run("C:/plink/plink --bfile plink --freq", dir);
		markerSet = new SnpMarkerSet(dir+"plink.map");
		markerNames = markerSet.getMarkerNames();
		log.report("Found "+markerNames.length+" of these markers in the PLINK dataset");
		
		try {
			reader = new BufferedReader(new FileReader(dir+"plink.bim"));
			writer = new PrintWriter(new FileWriter(dir+"plink.crf"));
			writer.println("plink");
			writer.println("plink.ped");
			writer.println("plink.map");
			writer.println("plink.frq");
			writer.println("counts.xln maskModel=true maskSex=true maskAffectionStatus=true");
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				writer.println(line[1]+"\tADD");
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + dir+"plink.bim" + "\" not found in current directory");
			return;
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + dir+"plink.bim" + "\"");
			return;
		}
		CreateDatabaseFromPlink.createDatabase(dir, dir+"plink.crf", log);
    	if (!Array.equals(markerNames, HashVec.loadFileToStringArray(dir+"plink.frq", true, new int[] {1}, false), true)) {
    		log.reportError("Error - the freq file does not match the map file");
    		return;
    	}
    	
		v = new Vector<String>();
    	if (originalCovariatesFile == null) {
    		Files.copyFile(dir + "counts.xln", dir + newCovariatesFile);
    	} else {
    		line = Files.getHeaderOfFile(baseDir+originalCovariatesFile, "[\\s]+", log);
    		numFixedCovariates = line.length-2;
    		indices = new int[numFixedCovariates];
    		for (int i = 2; i < line.length; i++) {
    			indices[i-2] = i;
    			v.add(line[i]);
			}
    		hash = HashVec.loadFileToHashString(baseDir + originalCovariatesFile, new int[] {0,1}, indices, false, "\t", false, false, false);
    		try {
				reader = new BufferedReader(new FileReader(dir+"counts.xln"));
				writer = new PrintWriter(new FileWriter(dir + newCovariatesFile));
				countMissingCovariate = 0;
				while (reader.ready()) {
					line = reader.readLine().trim().split("[\\s]+");
					if (hash.containsKey(line[0]+"\t"+line[1])) {
						writer.println(Array.toStr(line)+"\t"+hash.get(line[0]+"\t"+line[1]));
					} else {
						if (countMissingCovariate < 5) {
							log.reportError("Error - no covariate information for individual "+line[0]+"-"+line[1], true, true, 8);
						} else if (countMissingCovariate == 5) {
							log.reportError("...");
						}
						countMissingCovariate++;
					}
				}
				reader.close();
				writer.close();
			} catch (FileNotFoundException fnfe) {
				System.err.println("Error: file \"" + "counts.xln" + "\" not found in current directory");
				return;
			} catch (IOException ioe) {
				System.err.println("Error reading file \"" + "counts.xln" + "\"");
				return;
			}
    		if (countMissingCovariate > 5) {
    			log.reportError("There were "+countMissingCovariate+" indiviudals that were in the genotype file that were not in the covariate file");
    		}
    	}
	}


//	public static void addCountsAsCovariate(String baseDir, String plinkFiles, String originalCovariatesFile, String newCovariatesFile, String allPossibleSNPs, Logger log) {
//		BufferedReader reader;
//		PrintWriter writer;
//		String[] line;
//		Hashtable<String,String> hash;
//		Vector<String> v;
//		String[] markerNames;
//		SnpMarkerSet markerSet;
//		int[] indices;
//		int numFixedCovariates;
//		String dir;
//		
//		if (allPossibleSNPs == null) {
//			dir = baseDir+"allMarkers/";
//		} else {
//			dir = baseDir+ext.rootOf(allPossibleSNPs, false)+"/";
//		}
//		new File(dir).mkdirs();
//		
//		if (!new File(plinkFiles+".bed").exists() || !new File(plinkFiles+".bim").exists() || !new File(plinkFiles+".fam").exists()) {
//			log.reportError("Error - a full complement of binary PLINK files are required to perform this algorithm");
//			return;
//		}
//		if (allPossibleSNPs == null) {
//			log.report("All markers in the dataset will be added to the covars file, since no subset was pre-defined");
//			CmdLine.run("plink --bfile "+plinkFiles+" --make-bed", dir);
//		} else {
//			markerNames = HashVec.loadFileToStringArray(baseDir+allPossibleSNPs, false, new int[] {0}, true);
//			log.report("Curently using file '"+allPossibleSNPs+"', which has "+markerNames.length+" independent elements in it");
////			CmdLine.run("plink --bfile "+plinkFiles+" --extract ../"+(baseDir.equals("")?"":"../")+ext.removeDirectoryInfo(allPossibleSNPs)+" --make-bed", dir);
////			CmdLine.run("C:/plink/plink --bfile "+plinkFiles+" --extract ../"+(baseDir.equals("")?"":"../")+ext.removeDirectoryInfo(allPossibleSNPs)+" --make-bed", dir);
//			CmdLine.run("C:/plink/plink --bfile "+plinkFiles+" --extract ../"+ext.removeDirectoryInfo(allPossibleSNPs)+" --make-bed", dir);
//		}
//		
////		CmdLine.run("plink --bfile plink --recode", dir);
////		CmdLine.run("plink --bfile plink --freq", dir);
//		CmdLine.run("C:/plink/plink --bfile plink --recode", dir);
//		CmdLine.run("C:/plink/plink --bfile plink --freq", dir);
//		markerSet = new SnpMarkerSet(dir+"plink.map");
//		markerNames = markerSet.getMarkerNames();
//		log.report("Found "+markerNames.length+" of these markers in the PLINK dataset");
//		
//		try {
//			reader = new BufferedReader(new FileReader(dir+"plink.bim"));
//			writer = new PrintWriter(new FileWriter(dir+"plink.crf"));
//			writer.println("plink");
//			writer.println("plink.ped");
//			writer.println("plink.map");
//			writer.println("plink.frq");
//			writer.println("counts.xln maskModel=true");
//			while (reader.ready()) {
//				line = reader.readLine().trim().split("[\\s]+");
//				writer.println(line[1]+"\tADD");
//			}
//			reader.close();
//			writer.close();
//		} catch (FileNotFoundException fnfe) {
//			log.reportError("Error: file \"" + dir+"plink.bim" + "\" not found in current directory");
//			return;
//		} catch (IOException ioe) {
//			log.reportError("Error reading file \"" + dir+"plink.bim" + "\"");
//			return;
//		}
//		CreateDatabaseFromPlink.createDatabase(dir, dir+"plink.crf", log);
//    	if (!Array.equals(markerNames, HashVec.loadFileToStringArray(dir+"plink.frq", true, new int[] {1}, false), true)) {
//    		log.reportError("Error - the freq file does not match the map file");
//    		return;
//    	}
//    	
//		v = new Vector<String>();
//    	if (originalCovariatesFile == null) {
//    		Files.copyFile(dir+"counts.xln", dir+newCovariatesFile);
//    	} else {
//    		line = Files.getHeaderOfFile(baseDir+originalCovariatesFile, "[\\s]+", log);
//    		numFixedCovariates = line.length-2;
//    		indices = new int[numFixedCovariates];
//    		for (int i = 2; i < line.length; i++) {
//    			indices[i-2] = i;
//    			v.add(line[i]);
//			}
//    		hash = HashVec.loadFileToHashString(baseDir+originalCovariatesFile, new int[] {0,1}, indices, false, "\t", false, false, false);
//    		try {
//				reader = new BufferedReader(new FileReader(dir+"counts.xln"));
//				writer = new PrintWriter(new FileWriter(dir+newCovariatesFile));
//				while (reader.ready()) {
//					line = reader.readLine().trim().split("[\\s]+");
//					if (hash.containsKey(line[0]+"\t"+line[1])) {
//						writer.println(Array.toStr(line)+"\t"+hash.get(line[0]+"\t"+line[1]));
//					} else {
//						log.reportError("Error - no covariate information for individual "+line[0]+"-"+line[1], true, true, 8);
//					}
//				}
//				reader.close();
//				writer.close();
//			} catch (FileNotFoundException fnfe) {
//				System.err.println("Error: file \"" + "counts.xln" + "\" not found in current directory");
//				return;
//			} catch (IOException ioe) {
//				System.err.println("Error reading file \"" + "counts.xln" + "\"");
//				return;
//			}
//    	}
//	}

	public static void addDosageAsCovariate(String dir, String infoFile, String doseFile, String originalPhenotypeFile, boolean justIID, String newPhenotypeFile, String[] snps, boolean allowMissingConditional, boolean makeResiduals, Logger log) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		Hashtable<String,String> hash;
		Hashtable<String,Hashtable<String,String>> hashes;
		Vector<String> v;
		int count;
		String[] markerNames;
		SnpMarkerSet markerSet;
		boolean error;
		boolean machFormat;
		String fid, iid;
		
		machFormat = infoFile.endsWith(".mlinfo");
		
		markerSet = new SnpMarkerSet(dir+infoFile, false, new Logger());		
		markerNames = markerSet.getMarkerNames();
		log.report("Analyzing "+dir+doseFile+","+infoFile+", which contain "+markerNames.length+" markers");
		
		v = new Vector<String>();
		for (int i = 0; i < snps.length; i++) {
			if (!new File(dir+snps[i]+".ser").exists()) {
				v.add(snps[i]+"\t1");
			}
		}
		Files.writeList(Array.toStringArray(v), dir+"condis.txt");
		hashes = new Hashtable<String, Hashtable<String,String>>();
		if (v.size() > 0) {
			if (Mach.extractSpecificMarkers(dir, dir+"condis.txt", dir+doseFile, dir+infoFile, false, log)) { // Array.toStringArray(v)
				markerNames = new SnpMarkerSet(dir+"condis"+(machFormat?".mlinfo":".info"), false, new Logger()).getMarkerNames();
				for (int i = 0; i < markerNames.length; i++) {
					hashes.put(markerNames[i], new Hashtable<String, String>());
				}
				try {
					reader = new BufferedReader(new FileReader(dir+"condis"+(machFormat?".mldose":".dose")));
					while (reader.ready()) {
						line = reader.readLine().trim().split("[\\s]+");
						for (int i = 0; i < markerNames.length; i++) {
							fid = line[0].substring(0, line[0].indexOf("->"));
							iid = line[0].substring(line[0].indexOf("->")+2);
							if (justIID) {
								hashes.get(markerNames[i]).put(iid, line[i+2]);
							} else {
								hashes.get(markerNames[i]).put(fid+"\t"+iid, line[i+2]);
							}
						}
					}
					reader.close();
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \"" + dir+"condis."+(machFormat?".mldose":".dose") + "\" not found in current directory");
					return;
				} catch (IOException ioe) {
					System.err.println("Error reading file \"" + dir+"condis."+(machFormat?".mldose":".dose") + "\"");
					return;
				}
				for (int i = 0; i < markerNames.length; i++) {
					SerialHash.createSerializedStringHash(dir+markerNames[i]+".ser", hashes.get(markerNames[i]));
				}
			} else {
//				log.reportError("Missing a marker");
			}
		}
		
		error = false;
		for (int i = 0; i < snps.length; i++) {
			if (!hashes.containsKey(snps[i])) {
				if (new File(dir+snps[i]+".ser").exists()) {
					log.report("Loading pre-serialized "+dir+snps[i]+".ser", true, false);
					hashes.put(snps[i], SerialHash.loadSerializedStringHash(dir+snps[i]+".ser"));
				} else {
					if (allowMissingConditional) {
						log.reportError("Warning - '"+snps[i]+"' not available for '"+dir+"'; this marker will not be included as a covariate for this dataset");
					} else {
						log.reportError("Error - indicated SNP is not present in pre-serialized form (i.e. '"+dir+snps[i]+".ser') and is also not present in the info file ('"+infoFile+"')");
						error = true;
					}
				}
			}
		}
		if (error) {
//			return;
			System.exit(1);
		}

		try {
			reader = new BufferedReader(new FileReader(dir+originalPhenotypeFile));
			writer = new PrintWriter(new FileWriter(dir+newPhenotypeFile));
			line = reader.readLine().trim().split("[\\s]+");
			writer.print(Array.toStr(line));
			for (int i = 0; i < snps.length; i++) {
				hash = hashes.get(snps[i]);
				if (hash != null) {
					writer.print("\t"+snps[i]);
				}
			}
			writer.println();
			count = 0;
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				writer.print(Array.toStr(line));
				for (int i = 0; i < snps.length; i++) {
					hash = hashes.get(snps[i]);
					if (hash != null) {
						if (justIID && hash.containsKey(line[0])) {
							writer.print("\t"+hash.get(line[0]));
						} else if (hash.containsKey(line[0]+"\t"+line[1])) {
							writer.print("\t"+hash.get(line[0]+"\t"+line[1]));
						} else {
							writer.print("\t.");
							error = true;
						}
					}
				}
				writer.println();
				if (error) {
					count++;
				}
			}
			
			if (count > 0) {
				log.reportError("Warning - there were "+count+" indiviudals(s) present in the phenotype file that were absent in the genotype file");
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + dir+originalPhenotypeFile + "\" not found in current directory");
			return;
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + dir+originalPhenotypeFile + "\"");
			return;
		}
	}
	
	// probabel flag is not current implemented, neither is family data
	public static void run(String[] dirs, String filenames, String[] snps, String annotationFile, boolean useResiduals, boolean run, boolean probabel, boolean iterate, boolean genomiccontrol, String pheno, String phenoMissingValue, String outfile, String covariates, boolean sexAsCovariate, boolean allowMissingConditional, String doseFile, String infoFile) {
		Hashtable<String,String> allMarkers, freqs, annotation;
		Hashtable<String,Hashtable<String,String>> allResults;
		Vector<String> results, v, covars, metaInputs;
		DosageData dosageData;
		PrintWriter writer;
		int runningTally, countSig;
		String[][] pvals, baseCovars;
		String[] files, markerNames, chrPositions, annotationHeader;
		boolean done, header;
		Logger log;
		int count, minIndex, chrIndex, posIndex;
		double minP;
		byte[] chrs;
		int[] positions;
		boolean[] annotationHeaderMask;
		Vector<String> models;
		String filename;
		SnpMarkerSet markerSet;
		boolean[] dosageSets, familySets;
		String[] keys;
		int[] order;
		String[] line;
		boolean justIID;

		log = new Logger("Conditionals.log");
		annotation = new Hashtable<String, String>();
		annotationHeader = new String[0];
		if (new File(annotationFile).exists()) {
			annotationHeader = Files.getHeaderOfFile(annotationFile, "[\\s]+", log);
			annotation = HashVec.loadFileToHashString(annotationFile, new int[] {0}, Array.subArray(Array.intArray(annotationHeader.length), 1), false, "\t", true, false, false);
			annotationHeader = Array.subArray(annotationHeader, 1);
		} else if (annotationFile != null) {
			log.reportError("Error - annotation file '"+annotationFile+"' not found; no annotation will occur");
		}
		if (dirs == null) {
			dirs = new String[] {""};
		} else {
			for (int i = 0; i < dirs.length; i++) {
				dirs[i] = ext.verifyDirFormat(dirs[i]);
			}
		}
		dosageSets = new boolean[dirs.length];
		familySets = new boolean[dirs.length];
		
		boolean preferGenotyped = false; ////////////////////////////////////////////
		for (int k = 0; k < dirs.length; k++) {
			if (preferGenotyped && new File(dirs[k]+"plink.bed").exists() && new File(dirs[k]+"plink.bim").exists() && new File(dirs[k]+"plink.fam").exists()) {
				dosageSets[k] = false;
			} else if (pheno != null && new File(dirs[k]+pheno).exists()) {
				dosageSets[k] = true;
			} else if (new File(dirs[k]+"plink.bed").exists() && new File(dirs[k]+"plink.bim").exists() && new File(dirs[k]+"plink.fam").exists()) {
				dosageSets[k] = false;
			} else {
				System.err.println("Error - not enough info for directory '"+dirs[k]+"'");
			}
			familySets[k] = new File(dirs[k]+ext.rootOf(doseFile)+".fhsR").exists();
		}
		
		files = filenames.trim().split(",");
		if (snps == null || snps.length == 0 || snps[0].equals("null") || snps[0].equals("")) {
			snps = new String[0];
		}
		models = new Vector<String>();
		for (int i = 0; i < files.length; i++) {
			if (files[i] != null && files[i].equals("null")) {
				files[i] = null;
				allMarkers = new Hashtable<String, String>();
			} else {
				allMarkers = HashVec.loadFileToHashNull(files[i], false);
			}
			baseCovars = new String[dirs.length][];
			for (int k = 0; k < dirs.length; k++) {
				if (!dosageSets[k]) {
					filename = "conCovars"+(files[i]==null?"":"_"+ext.rootOf(files[i], false))+".dat";
					log.report("Generating allele counts for "+dirs[k]);
					if (covariates != null && !new File(dirs[k]+covariates).exists()) {
						log.reportError("Error: file "+dirs[k]+covariates+" not found, not including for this dataset");
						addCountsAsCovariate(dirs[k], null, filename, files[i], log);
					} else {
						addCountsAsCovariate(dirs[k], covariates, filename, files[i], log);
						if (covariates != null) {
							baseCovars[k] = Array.subArray(Files.getHeaderOfFile(dirs[k]+covariates, "[\\s]+", log), 2);
						}
					}
				}
			}
			done = false;
			count = 1;
			results = new Vector<String>();
			runningTally = -1;
			v = Array.toStringVector(snps);
			freqs = new Hashtable<String, String>();
			allResults = new Hashtable<String, Hashtable<String,String>>();
			while (!done) {
				if (iterate) {
					log.report("Running"+(files[i]==null?"":" "+ext.rootOf(files[i]))+" iteration #"+count);
				}
				metaInputs = new Vector<String>();
				for (int k = 0; k < dirs.length; k++) {
					if (dosageSets[k]) {
						if (useResiduals) {
							System.err.println("Error - useResiduals is not yet coded for dosage analyses");
							return;
						}
						if (outfile == null || (dirs != null && dirs.length > 1) || i>0 || files[i] != null) {
							outfile = ext.rootOf(pheno, false)+(files[i]==null?"":"_"+ext.rootOf(files[i], false))+"_iteration"+count+".dat";
//						} else if (count > 1) {
//							outfile = ext.rootOf(outfile, false)+count+outfile.substring(outfile.lastIndexOf(".")+1);
						}
						if (new File(dirs[k]+pheno).exists()) {
							line = Files.getHeaderOfFile(dirs[k]+pheno, "[\\s]+", log);
							if (line[0].equals("FID")) {
								justIID = false;
							} else if (line[0].equals("ID") || line[0].equals("IID")) {
								justIID = true;
							} else {
								log.reportError("Could not determine if phenotype file links up with a single ID variable or a FID/ID pair; assuming pair");
								justIID = false;
							}
							addDosageAsCovariate(dirs[k], infoFile, doseFile, pheno, justIID, outfile, Array.toStringArray(v), allowMissingConditional, useResiduals, log);
						} else {
							log.reportError("Error: file "+dirs[k]+pheno+" not found");
							done = true;
						}
						if (run || iterate) {
							if (new File(dirs[k]+doseFile+".ser").exists()) {
								dosageData = DosageData.load(dirs[k]+doseFile+".ser", false);
							} else {
								dosageData = new DosageData(dirs[k]+doseFile, dirs[k]+"list.txt", dirs[k]+infoFile, true, log);
								if (new File(ext.rootOf(dirs[k]+doseFile, false)+".map").exists()) {
									dosageData.getMarkerSet().setPositions(new SnpMarkerSet(ext.rootOf(dirs[k]+doseFile, false)+".map", false, new Logger()));
								} else if (new File(ext.rootOf(dirs[k]+doseFile, false)+".bim").exists()) {
									dosageData.getMarkerSet().setPositions(new SnpMarkerSet(ext.rootOf(dirs[k]+doseFile, false)+".bim", false, new Logger()));
								} else {
									log.reportError("Warning - no map file found in "+dirs[k]+"; delete '"+dirs[k]+doseFile+".ser' if map is added; otherwise use annotation file for positions");
								}
								dosageData.serialize(dirs[k]+doseFile+".ser");
							}

							markerNames = dosageData.getMarkerSet().getMarkerNames();
							chrPositions = dosageData.getMarkerSet().getChrAndPositions();
							for (int j = 0; j < markerNames.length; j++) {
								if (files[i] == null || allMarkers.containsKey(markerNames[j])) {
									allMarkers.put(markerNames[j], chrPositions[j]);
								}
							}
							if (new File(dirs[k]+pheno).exists()) {
//								if (familySets[k]) {
									
									
									
									
									
									
									
									
									
//								} else {
									dosageData.analyze(dirs[k]+outfile, phenoMissingValue, files[i], false, log);
									metaInputs.add(dirs[k]+ext.rootOf(outfile)+".se.metal");
//								}
							}
						}
					} else {
						if (useResiduals) {
							System.err.println("Error - useResiduals is not yet coded for PLINK analyses");
							return;
						}

						outfile = (files[i]==null?"":ext.rootOf(files[i], false)+"_")+"iteration"+count+".dat";
						if (run || iterate) {
							covars = Array.toStringVector(baseCovars[k]);
							for (int j = 0; j < v.size(); j++) {
								covars.add(v.elementAt(j));
							}
							
							// need to edit this to run quantitative and custom traits
							CmdLine.run("plink --bfile plink --logistic"+(sexAsCovariate?" --sex":"")+" --ci 0.95 --out iteration"+count+(covars.size()==0?"":" --covar conCovars"+(files[i]==null?"":"_"+ext.rootOf(files[i], false))+".dat --covar-name "+Array.toStr(Array.toStringArray(covars), ",")), dirs[k]+ext.rootOf(files[i]));
							Metal.convertPlinkResults(dirs[k]+ext.rootOf(files[i], false)+"/", "iteration"+count+".assoc.logistic", "ADD", "logistic", "plink.frq", true, true, dirs[k]+ext.rootOf(outfile)+".se.metal", true);
							metaInputs.add(dirs[k]+ext.rootOf(outfile)+".se.metal");

							markerSet = new SnpMarkerSet(dirs[k]+"plink.bim", false, new Logger());
							markerNames = markerSet.getMarkerNames();
							chrPositions = markerSet.getChrAndPositions();
							for (int j = 0; j < markerNames.length; j++) {
								if (files[i] == null || allMarkers.containsKey(markerNames[j])) {
									allMarkers.put(markerNames[j], chrPositions[j]);
								}
							}
						}
					}
				}
				if ((run || iterate) && dirs.length > 1) {
					try {
						writer = new PrintWriter(new FileWriter(ext.rootOf(outfile, false)+".Nweighted.batch"));
//						writer.println("MARKER MARKER");
//						writer.println("ALLELE REF OTHER");
//						writer.println("WEIGHT N");
//						writer.println("EFFECT DIR");
//						writer.println("PVALUE PVALUE");
						writer.println("MARKER MarkerName");
						writer.println("ALLELE Allele1 Allele2");
						writer.println("WEIGHT Weight");
						writer.println("EFFECT Direction");
						writer.println("PVALUE P-value");
						writer.println("");
						for (int j = 0; j < dirs.length; j++) {
//							writer.println("PROCESS "+dirs[j]+ext.rootOf(outfile, false)+".se.metal");
							writer.println("PROCESS "+metaInputs.elementAt(j));
						}
						writer.println("OUTFILE "+ext.rootOf(outfile, false)+".Nweighted .out");
						writer.println("ANALYZE");
						writer.println("");
						writer.println("QUIT");
						writer.close();
						CmdLine.run("metal < "+ext.rootOf(outfile, false)+".Nweighted.batch", "./");

						writer = new PrintWriter(new FileWriter(ext.rootOf(outfile, false)+".InvVar.batch"));
//						writer.println("MARKER MARKER");
//						writer.println("ALLELE REF OTHER");
//						writer.println("EFFECT beta");
//						writer.println("STDERR SE");
//						writer.println("SCHEME STDERR");
						writer.println("MARKER MarkerName");
						writer.println("ALLELE Allele1 Allele2");
						writer.println("EFFECT Effect");
						writer.println("STDERR StdErr");
						writer.println("SCHEME STDERR");
						if (genomiccontrol) {
							writer.println("GENOMICCONTROL ON");
						}
						writer.println("");
						for (int j = 0; j < dirs.length; j++) {
							if (genomiccontrol) {
								if (new File(dirs[j]+ext.rootOf(pheno, false)+"_GENOMICCONTROL.txt").exists()) {
									writer.println("GENOMICCONTROL "+HashVec.loadFileToStringArray(dirs[j]+ext.rootOf(pheno, false)+"_GENOMICCONTROL.txt", false, new int[] {0}, false)[0]);
								} else {
									log.reportError("Error - missing "+dirs[j]+ext.rootOf(pheno, false)+"_GENOMICCONTROL.txt; assuming a value of 1.00");
									writer.println("GENOMICCONTROL 1");
								}
							}
//							writer.println("PROCESS "+dirs[j]+ext.rootOf(outfile, false)+".se.metal");
							writer.println("PROCESS "+metaInputs.elementAt(j));
						}
						writer.println("OUTFILE "+ext.rootOf(outfile, false)+".InvVar .out");
						writer.println("ANALYZE");
						writer.println("");
						writer.println("QUIT");
						writer.close();
						CmdLine.run("metal < "+ext.rootOf(outfile, false)+".InvVar.batch", "./");
					} catch (Exception e) {
						System.err.println("Error writing to " + ext.rootOf(outfile, false)+".Nweighted.batch");
						e.printStackTrace();
					}
				}				
				if (iterate) {
					if (dirs.length > 1) {
						filename = ext.rootOf(outfile, false)+".InvVar1.out";
						pvals = HashVec.loadFileToStringMatrix(filename, true, new int[] {0,5,1,2,3,4}, false);
					} else {
						filename = dirs[0]+ext.rootOf(outfile, false)+".se.metal";
						pvals = HashVec.loadFileToStringMatrix(filename, true, new int[] {0,5,1,2,6,7}, false);
					}
					countSig = 0;
					minP = 1;
					minIndex = -1;
					for (int j = 0; j < pvals.length; j++) {
						if (Double.parseDouble(pvals[j][1]) < 0.05 && Math.abs(Double.parseDouble(pvals[j][4]))<5) {
							countSig++;
							if (Double.parseDouble(pvals[j][1]) < minP) {
								minP = Double.parseDouble(pvals[j][1]);
								minIndex = j;						
							}
						}
						HashVec.addToHashHash(allResults, count+"", pvals[j][0], pvals[j][1]);
					}
					if (runningTally == -1 ) {
						runningTally = allMarkers.size();
					}
					if (results.size() > 0) {
						results.add(results.remove(results.size()-1)+"\t"+(runningTally-countSig));
					}
					runningTally = countSig;

					models.add(filename+"\t"+(files[i]==null?"null":files[i])+"\t"+(minIndex == -1?"minimum":pvals[minIndex][0])+"\t"+(v.size()==0?"null":Array.toStr(Array.toStringArray(v),",")));
					
					if (minIndex == -1) {
						done = true;
					} else {
						v.add(pvals[minIndex][0]);
						results.add(pvals[minIndex][0]+"\t"+allMarkers.get(pvals[minIndex][0])+"\t"+pvals[minIndex][2]+"/"+pvals[minIndex][3]+"\t"+(freqs.containsKey(pvals[minIndex][0])?freqs.get(pvals[minIndex][0]):".")+"\t"+ext.formDeci(Math.exp(Double.parseDouble(pvals[minIndex][4])), 2, true)+" ("+ext.formDeci(Math.exp(Double.parseDouble(pvals[minIndex][4])-Double.parseDouble(pvals[minIndex][5])), 2, true)+"-"+ext.formDeci(Math.exp(Double.parseDouble(pvals[minIndex][4])+Double.parseDouble(pvals[minIndex][5])), 2, true)+")"+"\t"+pvals[minIndex][1]);
					}
					count++;
				} else {
					done = true;
				}
				
				if (count > 1) {	// number of hard-coded iterations
					done = true;
				}
			}

			if (iterate) {
				filename = ext.rootOf(outfile, false).substring(0, ext.rootOf(outfile, false).lastIndexOf("_iteration"));
				try {
					writer = new PrintWriter(new FileWriter(filename+"_table.xln"));
					writer.println("Marker\tChr\tPosition\tAlleles (Ref/Other)\tMAF\tOR (95% CI)\tpval\t# SNPs tagged");
					for (int j = 0; j < results.size(); j++) {
						writer.println(results.elementAt(j));
					}
					writer.close();
				} catch (Exception e) {
					System.err.println("Error writing to " + filename+"_table.xln");
					e.printStackTrace();
				}
				try {
					writer = new PrintWriter(new FileWriter(filename+"_fullTable.xln"));
					keys = HashVec.getKeys(allMarkers, false, false);
					chrs = new byte[keys.length];
					positions = new int[keys.length];
					chrIndex = ext.indexOfStr("chr", annotationHeader, false, true);
					posIndex = ext.indexOfStr("position", annotationHeader, false, true);
					annotationHeaderMask = Array.booleanArray(annotationHeader.length, true);
					if (chrIndex >= 0 && posIndex >= 0) {
						annotationHeaderMask[chrIndex] = false;
						annotationHeaderMask[posIndex] = false;
					}
					for (int j = 0; j < keys.length; j++) {
						if (chrIndex >= 0 && posIndex >= 0 && annotation.containsKey(keys[j])) {
							line = annotation.get(keys[j]).split("[\\s]+");
							chrs[j] = Byte.parseByte(line[chrIndex]);
							positions[j] = Integer.parseInt(line[posIndex]);
							allMarkers.put(keys[j], line[chrIndex]+"\t"+line[posIndex]);
						} else {
							if (chrIndex >= 0 && posIndex >= 0) {
								log.reportError("Error - annotation file does not contain information for "+keys[j]);
							}
							line = allMarkers.get(keys[j]).split("[\\s]+");
							if (line[0].equals("")) {
								allMarkers.put(keys[j], ".\t.");
								log.reportError("Error - no position information for '"+keys[j]+"'");
								chrs[j] = 0;
								positions[j] = -1;
							} else {
								chrs[j] = Byte.parseByte(line[0]);
								positions[j] = Integer.parseInt(line[1]);
							}
						}
					}
					order = Sort.orderTwoLayers(chrs, positions);
					writer.print("Marker\tChr\tPosition"+(annotationHeader.length>0?"\t"+Array.toStr(annotationHeader, annotationHeaderMask, "\t", "."):""));
					for (int j = 1; j < count; j++) {
						writer.print("\titeration"+1);
					}
					writer.println();
					for (int j = 0; j < keys.length; j++) {
						writer.print(keys[order[j]]+"\t"+allMarkers.get(keys[order[j]]));
						if (annotationHeader.length>0) {
							writer.print("\t"+Array.toStr(annotation.containsKey(keys[order[j]])?annotation.get(keys[order[j]]).split("[\\s]+"):Array.stringArray(annotationHeader.length, "."), annotationHeaderMask, "\t", "."));
						}
						for (int k = 1; k < count && allResults.containsKey(k+""); k++) {
							try {
							if (allResults.get(k+"").containsKey(keys[order[j]])) {
								writer.print("\t"+allResults.get(k+"").get(keys[order[j]]));
							} else {
								writer.print("\t.");
							}
							} catch (Exception e) {
								e.printStackTrace();
							}
						}
						writer.println();
					}
					writer.close();
				} catch (Exception e) {
					System.err.println("Error writing to " + filename+"_fullTable.xln");
					e.printStackTrace();
				}
				try {
					header = !new File(ext.parseDirectoryOfFile(outfile)+"All_tables.xln").exists();
					writer = new PrintWriter(new FileWriter(ext.parseDirectoryOfFile(outfile)+"All_tables.xln", true));
					if (header) {
						writer.println("Region\tNumMarkers\tMarker\tChr\tPosition\tAlleles (Ref/Other)\tMAF\tOR (95% CI)\tpval\t# SNPs tagged");
					}
					for (int j = 0; j < results.size(); j++) {
						writer.println((j==0?(files[i]==null?"allMarkers":ext.rootOf(files[i], false)):"")+"\t"+(j==0?allMarkers.size():"")+"\t"+results.elementAt(j));
					}
					writer.close();
				} catch (Exception e) {
					System.err.println("Error writing to " + ext.rootOf(outfile, false)+"_table.xln");
					e.printStackTrace();
				}
			}
		}
		
		if (iterate) {
			models.insertElementAt("#"+ext.getDate()+" "+ext.getTime(), 0);
			Files.writeList(Array.toStringArray(models), ext.parseDirectoryOfFile(outfile)+"All_models.dat");
		}
	}

	public static void runModels(String[] dirs, String fileWithModels, String annotationFile, boolean useResiduals, boolean genomiccontrol, String pheno, String phenoMissingValue, String outfile, String covariates, boolean sexAsCovariate, boolean allowMissingConditional, String dose, String info) {
		PrintWriter writer, w2;
		String[][] pvals, models;
		Logger log;
		int minIndex, index;
		double minP;
		String metaFile, snpFile, target;
		String[] snps;
		String filename;
		
		boolean reportMinimums = false;
		
		log = new Logger(ext.rootOf(fileWithModels, false)+".log");
		models = HashVec.loadFileToStringMatrix(fileWithModels, true, new int[] {0,1,2,3}, false);
		try {
			w2 = new PrintWriter(new FileWriter(ext.rootOf(fileWithModels, false)+"_runs.xln"));
			w2.println("Region\tNumMarkers\tMarker\tbeta\tpval\tMarker\tmeta-beta\tmeta-pval");
			for (int i = 0; i < models.length; i++) {
				metaFile = models[i][0];
				snpFile = models[i][1];
				target = models[i][2];
				snps = models[i][3].trim().split(",");

				// not using probabel
				run(dirs, snpFile, snps, annotationFile, useResiduals, true, false, false, genomiccontrol, pheno, phenoMissingValue, outfile, covariates, sexAsCovariate, allowMissingConditional, dose, info);
				
				filename = (pheno==null?"":ext.rootOf(pheno, false)+"_")+(snpFile.equals("null")?"":ext.rootOf(snpFile, false)+"_")+"iteration"+1;
				if (dirs.length > 1) {
					filename = filename+".InvVar1.out";
					pvals = HashVec.loadFileToStringMatrix(filename, true, new int[] {0,5,1,2,3,4}, false);
				} else {
					filename = dirs[0]+filename+".se.metal";
					pvals = HashVec.loadFileToStringMatrix(filename, true, new int[] {0,5,1,2,6,7}, false);
				}
				minP = 1;
				minIndex = -1;
				index = -1;
				for (int j = 0; j < pvals.length; j++) {
					if (pvals[j][0].equals(target)) {
						index = j;
					}
					if (Math.abs(Double.parseDouble(pvals[j][4]))<5 && Double.parseDouble(pvals[j][1]) < minP) {
						minP = Double.parseDouble(pvals[j][1]);
						minIndex = j;						
					}
				}
				if (target.equals("minimum")) {
					index = minIndex;
				}
				if (reportMinimums || !target.equals("minimum")) {
					w2.print((pheno==null?"plink":pheno)+(snpFile.equals("null")?"":"_"+snpFile)+"\t"+pvals.length+"\t"+(target.equals("minimum")?"min=":"")+(index==-1?".\t.\t.":pvals[index][0]+"\t"+pvals[index][4]+"\t"+pvals[index][1]));
				}

				try {
					//need to either pass or detect what kinds of header it is (or you could freaking standardize it!!)
					writer = new PrintWriter(new FileWriter(ext.rootOf(filename, true)+"_"+ext.rootOf(metaFile, true)+".InvVar.batch"));
					writer.println("MARKER MarkerName");
					writer.println("ALLELE Allele1 Allele2");
					writer.println("EFFECT Effect");
					writer.println("STDERR StdErr");
					writer.println("SCHEME STDERR");
//					writer.println("MARKER MARKER");
//					writer.println("ALLELE REF OTHER");
//					writer.println("EFFECT beta");
//					writer.println("STDERR SE");
//					writer.println("SCHEME STDERR");
					writer.println("");
					writer.println("PROCESS "+metaFile);
					writer.println("");
					writer.println("PROCESS "+filename);
					writer.println("");
					writer.println("OUTFILE "+ext.rootOf(filename, true)+"_"+ext.rootOf(metaFile, true)+".InvVar .out");
					writer.println("ANALYZE");
					writer.println("");
					writer.println("QUIT");
					writer.close();
					CmdLine.run("metal < "+ext.rootOf(filename, true)+"_"+ext.rootOf(metaFile, true)+".InvVar.batch", "./");
				} catch (Exception e) {
					log.reportError("Error writing to " + ext.rootOf(filename, true)+"_"+ext.rootOf(metaFile, true)+".Nweighted.batch");
					log.reportException(e);
				}
				
				pvals = HashVec.loadFileToStringMatrix(ext.rootOf(filename, true)+"_"+ext.rootOf(metaFile, true)+".InvVar1.out", true, new int[] {0,5,1,2,3,4}, false);
				minP = 1;
				minIndex = -1;
				index = -1;
				for (int j = 0; j < pvals.length; j++) {
					if (pvals[j][0].equals(target)) {
						index = j;
					}
					if (Math.abs(Double.parseDouble(pvals[j][4]))<5 && Double.parseDouble(pvals[j][1]) < minP) {
						minP = Double.parseDouble(pvals[j][1]);
						minIndex = j;						
					}
				}
				if (target.equals("minimum")) {
					index = minIndex;
				}
				if (reportMinimums || !target.equals("minimum")) {
					w2.println("\t"+(index==-1?".\t.\t.":(target.equals("minimum")?"min=":"")+pvals[index][0]+"\t"+pvals[index][4]+"\t"+pvals[index][1]));
				}
			}
			
			w2.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + ext.rootOf(fileWithModels, false)+"_runs.xln");
			log.reportException(e);
		}
	}
	
	public static void runAllRegions(String indexSNPs, String pheno, int columnIndexContainingAllMarkersForTheModel) {
		PrintWriter writer;
		String[] markerNames, addlMarkers;
		String markers;
		
		if (!Files.exists(pheno, false)) {
			System.err.println("Error - pheno file '"+pheno+"' does not exist");
			return;
		}
		
		markerNames = HashVec.loadFileToStringArray(indexSNPs, false, new int[] {0}, false);
		if (columnIndexContainingAllMarkersForTheModel >= 0) {
			addlMarkers = HashVec.loadFileToStringArray(indexSNPs, false, new int[] {columnIndexContainingAllMarkersForTheModel}, false);
		} else {
			addlMarkers = null;
		}
		try {
			writer = new PrintWriter(new FileWriter("runAllRegions"));
			for (int i = 0; i < markerNames.length; i++) {
				if (Files.exists(markerNames[i], false)) {
					writer.println("echo \""+(i+1)+" of "+markerNames.length+": "+markerNames[i]+"\"");
					writer.println("cp list.txt "+markerNames[i]+"/");
					writer.println("cd "+markerNames[i]);
					markers = addlMarkers == null ? markerNames[i] : addlMarkers[i];
					writer.println("jcp gwas.Conditional snps="+ext.replaceAllWith(markers, "_", ":")+" dose=data.dose info=data.info pheno=../"+pheno+" out=conPheno.dat");
					writer.println("awk '{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7}' data.info > data.pinfo");
					Files.qsub(markerNames[i]+"/"+markerNames[i]+".qsub", Probabel.EXECS[0]+" -p conPheno.dat -d data.dose -i data.pinfo -c 1 -o "+markers);
					writer.println("qsub "+markerNames[i]+".qsub");
					writer.println("cd ..");
					writer.println();
				} else {
					System.err.println("Error - no directory for '"+markerNames[i]+"'; need to run splitByRegion first?");
				}
			}
			writer.close();
			Files.chmod("runAllRegions");
			System.out.println("run ./runAllRegions to proceed with analyses");
		} catch (Exception e) {
			System.err.println("Error writing to " + "runAllRegions");
			e.printStackTrace();
		}
	}
		
	public static void checkAllRegions(String indexSNPs, int columnIndexContainingAllMarkersForTheModel, String markerList, int markerNameIndex) {
		BufferedReader reader;
		PrintWriter writer;
		String[] markerNames, markers;
		CountHash ch;
		int exp, obs;
		String filename;
		
		ch = new CountHash();
		try {
			reader = new BufferedReader(new FileReader(markerList));
			while (reader.ready()) {
				ch.add(ext.replaceAllWith(reader.readLine().trim().split("[\\s]+")[markerNameIndex], ":", "_"));
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + markerList + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + markerList + "\"");
			System.exit(2);
		}
		
		markerNames = HashVec.loadFileToStringArray(indexSNPs, false, new int[] {0}, false);
		if (columnIndexContainingAllMarkersForTheModel == -1) {
			markers = HashVec.loadFileToStringArray(indexSNPs, false, new int[] {0}, false);
		} else {
			markers = HashVec.loadFileToStringArray(indexSNPs, false, new int[] {columnIndexContainingAllMarkersForTheModel}, false);
		}
		try {
			writer = new PrintWriter(new FileWriter("checkAllRegions.out"));
			writer.println("MarkerName\texpected\tobserved");
			for (int i = 0; i < markerNames.length; i++) {
				exp = ch.getCount(markerNames[i]);
				writer.print(markerNames[i]+"\t"+exp+"\t");
				filename = markerNames[i]+"/"+markers[i]+"_add.out.txt";
				if (Files.exists(filename, false)) {
					obs = Files.countLines(filename, true);
					writer.println(obs);
					if (exp != obs) {
						System.err.println("Error - mismatch for "+markerNames[i]+" expecting "+exp+", found "+obs);
					}					
				} else {
					System.err.println("Error - no results for '"+markerNames[i]+"'; could not find file "+filename);
					writer.println(".");
				}
			}
			writer.close();
			Files.chmod("runAllRegions");
		} catch (Exception e) {
			System.err.println("Error writing to " + "checkAllRegions.out");
			e.printStackTrace();
		}
	}
	
	public static void parseResults(String indexSNPs, int columnIndexContainingAllMarkersForTheModel, String name) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String[] regionNames, markers;
		
		regionNames = HashVec.loadFileToStringArray(indexSNPs, false, new int[] {0}, false);
		if (columnIndexContainingAllMarkersForTheModel == -1) {
			markers = HashVec.loadFileToStringArray(indexSNPs, false, new int[] {0}, false);
		} else {
			markers = HashVec.loadFileToStringArray(indexSNPs, false, new int[] {columnIndexContainingAllMarkersForTheModel}, false);
		}
		new File("results/").mkdirs();
		for (int i = 0; i < regionNames.length; i++) {
			try {
				reader = new BufferedReader(new FileReader(regionNames[i]+"/"+markers[i]+"_add.out.txt"));
				writer = new PrintWriter(new FileWriter("results/"+name+"_"+regionNames[i]+".txt"));
				line = reader.readLine().trim().split("[\\s]+");
				ext.checkHeader(line, Probabel.LOGIST_OUTPUT_HEADER, true);
				writer.println("MarkerName\tEffect_allele\tReference_allele\tBETA\tSE\tP");
				while (reader.ready()) {
					line = reader.readLine().trim().split("[\\s]+");
    		        writer.println(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[10]+"\t"+line[11]+"\t"+(line[12].equals("nan")||line[12].equals("-inf")?".":ProbDist.ChiDist(Math.max(Double.parseDouble(line[12]), 0), 1)));
				}
				writer.close();
				reader.close();
				Zip.gzip("results/"+name+"_"+regionNames[i]+".txt");
			} catch (FileNotFoundException fnfe) {
				System.err.println("Error: file \"" + regionNames[i]+"/"+markers[i]+"_add.out.txt" + "\" not found in current directory");
				System.exit(1);
			} catch (IOException ioe) {
				System.err.println("Error reading file \"" + regionNames[i]+"/"+markers[i]+"_add.out.txt" + "\"");
				System.exit(2);
			}
		}
	}
	
	public static void metaAllRegions(String controlFile, Logger log) {
		PrintWriter writer, w2, w3;
		PrintStream ps;
		String[] line;
		String trav;
        Vector<String> paramV;
        String regionListFilename;
        String[][] dirsAndPatterns;
        String[][] markersAndChrs;
        int minNumValids;
        int[] indices;
        String[][] results;
		double minValue;
		String minMarker, minInfo;
		int index;
		String study, filename;
		String markerName, effectAllele, refAllele, beta, se, pval, weight, freq, gc, snps;
		boolean minmaxfreq, averagefreq, verbose, schemeIsStdErr;

		paramV = Files.parseControlFile(controlFile, "meta", new String[] {
				"# Can set regionList file to null in order to just do a standard meta-analysis, just set first column (the directories) to .",
				"regionListWithAnyExtraInfoForFilenamePatterns.dat minNum=7",
				"progeni_genepd	progeni_genepd_[%1].txt", 
				"hihg	hihg_[%1].txt", 
				"ngrc	ngrc_[%1].txt",
				"Germany	GER_chr[%2]_[%1]",
				"Netherlands	NL_chr[%2]_[%1]",
				"NIA	NIA_chr#_[%1]_[%3]",
				"23andMe	TTAM_[%1]_[%3].txt",
				"",
				"# Variable names will be set as follows (alternate scheme is SAMPLESIZE):",
				"MARKER MarkerName",
				"ALLELE Effect_allele Reference_allele",
				"EFFECT BETA",
				"STDERR SE",
				"PVAL P-value",
				"WEIGHT N",
				"SCHEME STDERR",
				"GENOMICCONTROL OFF",
				"",
				"# Additional advanced options available",
				"FREQ EFFECT_ALLELE_FREQ",
				"MINMAXFREQ ON",
				"AVERAGEFREQ ON",
				"ADDFILTER SNP IN (rs10830963,rs563694)",
				"VERBOSE ON",
				}, log);
		
		if (paramV != null) {
			line = paramV.elementAt(0).split("[\\s]+");
			regionListFilename = line[0];
			minNumValids = -9;
			for (int i = 1; i < line.length; i++) {
				if (line[i].toLowerCase().startsWith("minnum=")) {
					minNumValids = ext.parseIntArg(line[i]);
				} else {
					System.err.println("Error - do not know what to do with parameter '"+line[i]+"'");
				}
			}
			
			schemeIsStdErr = true;
			minmaxfreq = false;
			averagefreq = false;
			verbose = false;
			markerName = effectAllele = refAllele = beta = se = pval = weight = freq = snps = null;
			gc = "OFF";


			for (int i = paramV.size()-1; i > 0; i--) {
				trav = paramV.elementAt(i).trim();
				if (trav.startsWith("MARKER ")) {
					markerName = paramV.remove(i).split("[\\s]+")[1];
				} else if (trav.startsWith("ALLELE ")) {
					line = paramV.remove(i).trim().split("[\\s]+");
					effectAllele = line[1];
					refAllele = line[2];
				} else if (trav.startsWith("EFFECT ")) {
					beta = paramV.remove(i).split("[\\s]+")[1];
				} else if (trav.startsWith("STDERR ")) {
					se = paramV.remove(i).split("[\\s]+")[1];
				} else if (trav.startsWith("PVAL ")) {
					pval = paramV.remove(i).split("[\\s]+")[1];
				} else if (trav.startsWith("WEIGHT ")) {
					weight = paramV.remove(i).split("[\\s]+")[1];
				} else if (trav.startsWith("SCHEME ")) {
					schemeIsStdErr = paramV.remove(i).split("[\\s]+")[1].equalsIgnoreCase("STDERR");
				} else if (trav.startsWith("GENOMICCONTROL ")) {
					gc = paramV.remove(i).split("[\\s]+")[1];
				} else if (trav.startsWith("FREQ ")) {
					freq = paramV.remove(i).split("[\\s]+")[1];
				} else if (trav.startsWith("ADDFILTER ")) {
					snps = trav;
				} else if (trav.startsWith("MINMAXFREQ ")) {
					minmaxfreq = paramV.remove(i).split("[\\s]+")[1].equalsIgnoreCase("ON");
				} else if (trav.startsWith("AVERAGEFREQ ")) {
					averagefreq = paramV.remove(i).split("[\\s]+")[1].equalsIgnoreCase("ON");
				} else if (trav.startsWith("VERBOSE ")) {
					verbose = paramV.remove(i).split("[\\s]+")[1].equalsIgnoreCase("ON");
				} else if (trav.split("[\\s]+").length != 2) {
					System.err.println("Error - processing '"+paramV.elementAt(i)+"'");
				}
			}			
			
			dirsAndPatterns = new String[paramV.size()-1][2];
			for (int i = 0; i < dirsAndPatterns.length; i++) {
				line = paramV.elementAt(1+i).trim().split("[\\s]+");
				if (line.length != 2) {
					System.err.println("Error - processing '"+paramV.elementAt(1+i)+"'");
					return;
				}
				line[0] = ext.verifyDirFormat(line[0]);
				dirsAndPatterns[i] = line;
			}
		
			if (regionListFilename.equalsIgnoreCase("null")) {
				markersAndChrs = new String[][] {{"Metal_SE"}};
			} else {
				markersAndChrs = HashVec.loadFileToStringMatrix(regionListFilename, false, Array.intArray(Files.getHeaderOfFile(regionListFilename, log).length), "\t", false, 100, false);
			}
			System.out.println("minNumValid="+minNumValids);

			new File("metas/").mkdir();
			try {
				ps = new PrintStream(new File(ext.rootOf(controlFile)+".log"));
				w2 = new PrintWriter(new FileWriter("minPvals.out"));
				w2.print("IndexSnp\tMinConditionalSNP\tMinConditionalPval\tAllele1\tAllele2\tConditionalBeta\tConditionalStdErr"); 
				for (int i = 0; i < dirsAndPatterns.length; i++) {
					study = dirsAndPatterns[i][0].substring(0, dirsAndPatterns[i][0].length()-1);
					w2.print("\t"+study+"_beta\t"+study+"_SE\t"+study+"_pval");
				}
				w2.println();
				w3 = new PrintWriter(new FileWriter("allPvals.out"));
				w3.println("IndexSnp(s)\tMinConditionalSNP\tAllele1\tAllele2\tConditionalBeta\tConditionalStdErr\tMinConditionalPval\tDirection"); 
				for (int i = 0; i < markersAndChrs.length; i++) {
					System.out.println("Processing "+(i+1)+" of "+markersAndChrs.length+" "+markersAndChrs[i][0]);
					if (!Files.exists("metas/"+markersAndChrs[i][0]+"_SE1.tbl", false)) {
						try {
							writer = new PrintWriter(new FileWriter("metas/"+markersAndChrs[i][0]+".metal"));
							writer.println("MARKER "+markerName);
							writer.println("ALLELE "+effectAllele+" "+refAllele);
							if (schemeIsStdErr) {
								writer.println("EFFECT "+beta);
								writer.println("STDERR "+se);
								writer.println("SCHEME STDERR");
							} else {
								writer.println("PVAL "+pval);
								writer.println("WEIGHT "+weight);
								writer.println("SCHEME SAMPLESIZE");
							}
							writer.println("GENOMICCONTROL "+gc);
							if (freq != null) {
								writer.println("FREQ "+freq);
							}
							if (minmaxfreq) {
								writer.println("MINMAXFREQ ON");
							}
							if (averagefreq) {
								writer.println("AVERAGEFREQ ON");
							}
							if (snps != null) {
								writer.println(snps);
							}
							if (verbose) {
								writer.println("VERBOSE ON");
							}
							writer.println();
							for (int j = 0; j < dirsAndPatterns.length; j++) {
								filename = dirsAndPatterns[j][1];
								for (int k = 0; k < markersAndChrs[i].length; k++) {
									filename = ext.replaceAllWith(filename, "[%"+(k+1)+"]", markersAndChrs[i][k]);
								}								
								writer.println("PROCESS "+dirsAndPatterns[j][0]+filename);
							}
							writer.println();
							writer.println("OUTFILE metas/"+markersAndChrs[i][0]+"_SE .tbl");
							writer.println("ANALYZE HETEROGENEITY");
							writer.close();
						} catch (Exception e) {
							System.err.println("Error writing to " + "metas/"+markersAndChrs[i][0]+".metal");
							e.printStackTrace();
						}
		
						CmdLine.run("metal < "+"metas/"+markersAndChrs[i][0]+".metal", "./", ps);
					}
					if (!Files.exists("metas/"+markersAndChrs[i][0]+"_SE1.tbl", false)) {
						log.reportError("Error - failed to run meta-analysis for "+markersAndChrs[i][0]);
					} else {
						indices = ext.indexFactors(new String[] {"MarkerName", "Direction", "P-value", "Allele1", "Allele2", "Effect", "StdErr"}, Files.getHeaderOfFile("metas/"+markersAndChrs[i][0]+"_SE1.tbl", "\t", log), false, true);
						results = HashVec.loadFileToStringMatrix("metas/"+markersAndChrs[i][0]+"_SE1.tbl", true, indices, false);
						
						minValue = 1;
						minMarker = null;
						minInfo = null;
						String indexSNPs = markersAndChrs[i][0];
						for (int j = 2; j < markersAndChrs[i].length; j++) {
							indexSNPs += ";"+markersAndChrs[i][j];
						}
						for (int j = 0; j < results.length; j++) {
							if (minNumValids < 0 || ext.replaceAllWith(results[j][1], "?", "").length() >= minNumValids) {
								if (Double.parseDouble(results[j][2]) < minValue) {
									minValue = Double.parseDouble(results[j][2]);
									minMarker = results[j][0];
									minInfo = results[j][3]+"\t"+results[j][4]+"\t"+results[j][5]+"\t"+results[j][6];
								} else if (Double.parseDouble(results[j][2]) == minValue) {
									minMarker += ";"+results[j][0];
								}
								w3.println(indexSNPs+"\t"+results[j][0]+"\t"+results[j][3]+"\t"+results[j][4]+"\t"+results[j][5]+"\t"+results[j][6]+"\t"+results[j][2]+"\t"+results[j][1]);
							}
						}
						if (minMarker == null) {
							log.reportError("Error - no result for "+markersAndChrs[i][0]+" had valid data for "+minNumValids+" or more studies");
						} else {
							if (minMarker.indexOf(";") > 0) {
								minMarker = minMarker.substring(0, minMarker.indexOf(";"));
								log.report("There was a tie for the minimum marker for the "+markersAndChrs[i][0]+" meta-analysis ("+ext.listWithCommas(minMarker.split(";"))+")");
							}
							minMarker = ext.replaceAllWith(minMarker, "_", ":");
						}
						w2.print(ext.replaceAllWith(markersAndChrs[i][0], "_", ":"));
						if (minMarker == null) {
							w2.print("\t.\t.\t.\t.");
						} else {
							w2.print("\t"+minMarker+"\t"+minValue+"\t"+minInfo);
							for (int j = 0; j < dirsAndPatterns.length; j++) {
								filename = dirsAndPatterns[j][0]+dirsAndPatterns[j][1];
								for (int k = 0; k < markersAndChrs[i].length; k++) {
									filename = ext.replaceAllWith(filename, "[%"+(k+1)+"]", markersAndChrs[i][k]);
								}								
								if (Files.exists(filename, false)) {
									indices = ext.indexFactors(new String[] {"MarkerName", "BETA", "SE", "P"}, Files.getHeaderOfFile(filename, "\t", log), false, true);
									results = HashVec.loadFileToStringMatrix(filename, true, indices, false);
									index = ext.indexOfStr(minMarker, Matrix.extractColumn(results, 0));
									if (index == -1) {
										System.err.println("Error - couldn't find '"+minMarker+"' in "+dirsAndPatterns[j][0]+ext.replaceAllWith(dirsAndPatterns[j][1], new String[][] {{"[name]", markersAndChrs[i][0]}, {"#", markersAndChrs[i][1]}}));
										w2.print("\t.\t.\t.");
									} else {
										w2.print("\t"+Array.toStr(Array.subArray(results[index], 1)));
									}
								} else {
									System.err.println("Error - file not found for analysis "+markersAndChrs[i][0]+" for "+dirsAndPatterns[j][0]);
									w2.print("\tmissing\tmissing\tmissing");
								}
							}
						}
						w2.println();
					}
				}			
				w2.close();
				w3.close();
			} catch (Exception e) {
				System.err.println("Error writing to " + "minPvals.out");
				e.printStackTrace();
			}
		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filenames = "null";
		String pheno = null;
		String phenoMissingValue = "NA";
		boolean sexAsCovariate = false;
		String covariates = null;
		boolean useResiduals = false;
		String[] snps = null;
		String dose = null;
		String info = null;
		String[] dirs = null;
		boolean run = false;
		boolean iterate = false;
		String outfile = null;
		boolean genomiccontrol = false;
		String annotation = "annotation.xln";
		String models = null;
		boolean allowMissingConditional = false;
		boolean probabel = false;
		String indexMarkers = null;
		String markerList = null;
		int markerNameIndex = -1;
		String parse = null;
		String metaAll = null;
		int columnIndexContainingAllMarkersForTheModel = -1;

//		filenames = "GBA.txt,GAK.txt,BST1.txt,SNCA.txt,chr5.txt,HLA.txt,chr12.txt,chr16.txt,MAPT.txt,RIT2.txt";
////		filenames = "SNCA.txt";
////		filenames = "RIT2.txt";
//		
////		dirs = new String[] {"CIDR/"};
////		dirs = new String[] {"CIDR/", "Fung/", "Sing550/", "LEAPS/", "Miami/", "NGRC/"};
//		dirs = new String[] {"CIDR/", "Fung/", "Sing550/", "Miami/", "NGRC/"};
////		mldose = "CIDR/chr4/MACH_step2_500_chr4.mldose";
////		mlinfo = "CIDR/chr4/MACH_step2_500_chr4.mlinfo";
//		mldose = "Final693.mldose";
//		mlinfo = "Final693.mlinfo";
//		pheno = "pheno.dat";
////		snps = new String[] {"rs356165"};
//		snps = new String[] {};
////		snps = new String[] {"rs6854087"};
//		run = true;
//		iterate = true;
//		genomiccontrol = true;
		
////		filenames = "GBA.txt,GAK.txt,BST1.txt,SNCA.txt,chr5.txt,HLA.txt,chr12.txt,chr16.txt,MAPT.txt,RIT2.txt";
////		filenames = "SNCA.txt";
////		filenames = "HLA.txt";
//		dirs = new String[] {"Replication/"};
//		sexAsCovariate = true;
//		covariates = "PC1_covar.dat";
//		run = true;
////		iterate = true;
////		snps = new String[] {"rs356220"};
////		snps = new String[] {"rs9275184"};
//		models = "final_combined_models.dat";
		
////		filenames = "HLA.txt";
//		dirs = new String[] {"CIDR/", "Fung/", "Sing550/", "Miami/", "NGRC/"};
////		snps = new String[] {"rs9275184"};
////		models = "replication_models.dat";
//		mldose = "Final693.mldose";
//		mlinfo = "Final693.mlinfo";
//		pheno = "pheno.dat";
//		run = true;
//		models = "final_combined_models.dat";
		
//		filenames = "GBA.txt,GAK.txt,BST1.txt,SNCA.txt,chr5.txt,HLA.txt,chr12.txt,chr16.txt,MAPT.txt,RIT2.txt";
////		filenames = "HLA.txt";
//		dirs = new String[] {"CIDR/", "Fung/", "Sing550/", "Miami/", "NGRC/", "Replication/"};
//		sexAsCovariate = true;
//		covariates = "PC1_covar.dat";
//		mldose = "Final693.mldose";
//		mlinfo = "Final693.mlinfo";
//		pheno = "pheno.dat";
//		iterate = true;
//		genomiccontrol = true;
		
//		dirs = new String[] {"CIDR/", "Fung/", "Sing550/", "Miami/", "NGRC/"};
//		sexAsCovariate=false;
//		covariates = null;
//		dose = "critHLA.dose";
//		info = "critHLA.minfo";
//		pheno = "pheno.dat";
//		run = true;
//		genomiccontrol = true;
//		iterate = true;
//		models = "models.dat";

////		dirs = new String[] {"whites/ARIC/", "whites/CARDIA/", "whites/CFS/", "whites/CHS/", "whites/FHS/", "whites/MESA/"};
////		dirs = new String[] {"blacks/CARDIA/", "blacks/CFS/", "blacks/CHS/", "blacks/MESA/"};
////		dirs = new String[] {"whites/ARIC/"};
//		dirs = new String[] {"whites/ARIC/", "whites/CARDIA/", "whites/CFS/", "whites/CHS/", "whites/FHS/", "whites/MESA/", "blacks/CARDIA/", "blacks/CHS/", "blacks/CFS/", "blacks/MESA/", "asians/MESA/", "hispanics/MESA/"};
////		dirs = new String[] {"whites/MESA/", "blacks/MESA/", "asians/MESA/", "hispanics/MESA/"};
//		//"blacks/CHS/", 
//		sexAsCovariate=false;
////		covariates = "plink_pheno.dat";
//		dose = "abo_icam.mldose";
//		info = "abo_icam.pinfo";
//		pheno = "plink_pheno.dat";
////		phenoMissingValue = null;
//		phenoMissingValue = "-9";
//		run = true;
//		genomiccontrol = false;
//		iterate = true;
////		models = "models.dat";
////		snps = new String[] {"rs5498"};
////		snps = new String[] {"rs5498", "rs1799969"};
//		snps = new String[] {"rs5498", "rs1799969", "rs651007"};
////		allowMissingConditional = false;
//		allowMissingConditional = true;

//		dirs = new String[] {"CIDR/", "Fung/", "Sing550/", "Miami/", "NGRC/"};
//		sexAsCovariate=false;
//		covariates = null;
//		dose = "hla.dose";
//		info = "hla.minfo";
//		pheno = "pheno.dat";
//		run = true;
//		genomiccontrol = true;
//		snps = new String[] {"rs2395163"};
//		iterate = false;

		
		String usage = "\n" + 
		"gwas.Conditional requires 0-1 arguments\n" + 
		"  Parameters for all types of runs:\n" + 
		"   (1) filenames of SNP subsets to analyze (i.e. files=" + filenames + " (default; set to null for all markers; use a comma to differentiate separate runs))\n" + 
		"   (2) generate residuals first instead of including everything in the genetic model (i.e. residuals="+useResiduals+" (default))\n" + 
		"   (3) perform only run or first iteration with specific SNPs (i.e. snps=rs6830724,rs1560489,rs754750 (default is snps=null))\n" + 
		"   (4) run the model as well (i.e. -run (not the default))\n" +
		"   (5) run dynamic iterative models instead (i.e. -iterate (not the default))\n" + 
		"   (6) marker annotation, especially chr and position fields (i.e. annotation="+annotation+" (default))\n" + 
		"   (7) directories to perform analyses and to meta-analyze (i.e. dirs=CIDR,Miami,NGRC (not the default; leave as null for a single dataset))\n" + 
		"   (8) use genomic control in meta-analysis (i.e. gc="+genomiccontrol+" (default; requires a '[pheno]_GENOMICCONTROL.txt' file in each directory that includes nothing but the lambda value))\n" + 
		"   (9) return a p-value for specific marker given a particular model (i.e. return=rs356220 (not the default))\n" + 
		"  (10) run specific models, e.g. to replicate previous conditionals (i.e. models=All_previous_models.dat (not the default))\n" + 
		"  Also, for PLINK data:\n" +
		"   (a) name of original covariates file (i.e. covariates=" + covariates + " (default))\n" +
		"   (b) name of new covariates file (i.e. out=conCovars.dat (default))\n" +
		"   (c) include sex as a covariate (i.e. sex=" + sexAsCovariate + " (default))\n" +
		"  Or for dosage data:\n" +
		"   (a) dosage file of consolidated imputed data (i.e. dose=chr4.mldose (not the default))\n" +
		"   (b) info file of consolidated imputed data (i.e. info=chr4.mlinfo (not the default))\n" +
		"   (c) phenotype file with any covariates (i.e. pheno=pheno.dat (not the default; FID IID PHENO COVAR1 COVAR2 ...))\n" +
		"   (d) name of new phenotype file (i.e. out=conPheno.dat (default))\n" +
		"   (e) if ProbABEL is installed, use that instead (i.e. -probabel (not the default))\n" +
		" OR\n" + 
		"   (1) run conditionals for all index markers in file (i.e. indexMarkers=" + indexMarkers + " (not the default; should be output file from splitByRegion))\n" + 
		"   (2) name of phenotype file (i.e. pheno="+pheno+" (default))\n" +
		"   (3) (optional) index of column containing all SNPs that should be included in the model (i.e. addlMarkers=1 (not the default))\n" +
		"       if this index is provided, then the original index SNP in the first column needs to be included\n" +
		"       the first column may then become the name of the directory (which is probably still just the name of the original SNP)\n" +
		" OR\n" + 
		"   (1) check to make sure that all conditionals completed with the expected number of markers (i.e. indexMarkers=" + indexMarkers + " (not the default; should be output file from splitByRegion))\n" + 
		"   (2) marker list that was originally extracted (i.e. markerList=extractInfo.dat (not the default))\n"+
		"   (3) index of the region name within the marker list file (i.e. markerNameIndex=[index of region name in markerList file] (not the default))\n"+
		"   (4) (optional) index of column containing all SNPs that should have been included in the model (i.e. addlMarkers=1 (not the default))\n" +
		" OR\n" + 
		"   (1) parse results converting chi square to p-value (i.e. indexMarkers=" + indexMarkers + " (not the default; should be output file from splitByRegion))\n" + 
		"   (2) name of study (i.e. parse=progeni_genepd (not the default))\n" + 
		"   (3) (optional) index of column containing all SNPs that should have been included in the model (i.e. addlMarkers=1 (not the default))\n" +
		" OR\n" + 
		"   (1) meta-analyze all regions using control file (i.e. metaAll=C:/mega/iteration1/control.dat (not the default))\n" + 
		"       contains keyword meta," +
		"		then splitByRegion filename," +
		"		then all subdirectories in first column and then a filename template where [name]=markerName and #=chr#\n" + 
		"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				return;
			} else if (args[i].startsWith("files=")) {
				filenames = ext.parseStringArg(args[i], "null");
				numArgs--;
			} else if (args[i].startsWith("residuals=")) {
				useResiduals = ext.parseBooleanArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("snps=")) {
				snps = args[i].split("=")[1].split(",");
//				if (snps[0].equals("null") || snps[0].equals("")) {
//					snps = null;
//				}
				numArgs--;
			} else if (args[i].startsWith("-run")) {
				run = true;
				numArgs--;
			} else if (args[i].startsWith("-iterate")) {
				iterate = true;
				numArgs--;
			} else if (args[i].startsWith("annotation=")) {
				annotation = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("covariates=")) {
				covariates = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("sex=")) {
				sexAsCovariate = ext.parseBooleanArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("pheno=")) {
				pheno = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("missingPheno=")) {
				phenoMissingValue = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("-allowmissingConditional")) {
				allowMissingConditional = true;
				numArgs--;
			} else if (args[i].startsWith("out=")) {
				outfile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("dose=")) {
				dose = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("info=")) {
				info = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("dirs=")) {
				dirs = args[i].split("=")[1].split(",");
				if (dirs.length == 1) {
					System.err.println("Warning - need to define multiple directories before a meta-analysis is performed");
				}
				numArgs--;
			} else if (args[i].startsWith("gc=")) {
				genomiccontrol = ext.parseBooleanArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("models=")) {
				models = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("indexMarkers=")) {
				indexMarkers = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("markerList=")) {
				markerList = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("addlMarkers=")) {
				columnIndexContainingAllMarkersForTheModel = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("markerNameIndex=")) {
				markerNameIndex = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("parse=")) {
				parse = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("-probabel")) {
				probabel = true;
				numArgs--;
			} else if (args[i].startsWith("metaAll=")) {
				metaAll = ext.parseStringArg(args[i], null);
				numArgs--;
			} else {
				System.err.println("Error - don't know what to do with argument: "+args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			return;
		}
		try {
//			dirs = new String[] {"rs17621489/"};
//			dose = "data.dose";
//			info = "data.pinfo";
//			pheno = "pheno.dat";
//			snps = new String[] {"rs17621489"};
//			outfile = "conPheno.dat";
//			metaAll = "meta.crf";
			if (markerNameIndex > 0) {
				checkAllRegions(indexMarkers, columnIndexContainingAllMarkersForTheModel, markerList, markerNameIndex);
			} else if (parse != null) {
				parseResults(indexMarkers, columnIndexContainingAllMarkersForTheModel, parse);
			} else if (indexMarkers != null) {
				runAllRegions(indexMarkers, pheno, columnIndexContainingAllMarkersForTheModel);
			} else if (metaAll != null) {
				metaAllRegions(metaAll, new Logger("metaAll.log"));
			} else if (models != null) {	// not convinced this works with the output from run(with iterate)
				runModels(dirs, models, annotation, useResiduals, genomiccontrol, pheno, phenoMissingValue, outfile, covariates, sexAsCovariate, allowMissingConditional, dose, info);
			} else {
				run(dirs, filenames, snps, annotation, useResiduals, run, probabel, iterate, genomiccontrol, pheno, phenoMissingValue, outfile, covariates, sexAsCovariate, allowMissingConditional, dose, info);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
