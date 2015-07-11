package one;

import gwas.Conditional;
import gwas.Metal;
import gwas.ResultsPackager;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Vector;

import bioinformatics.MapSNPsAndGenes;
import parse.GenParser;
import stats.RegressionModel;
import common.Aliases;
import common.Array;
import common.CmdLine;
import common.Files;
import common.HashVec;
import common.Logger;
import common.Sort;
import common.Unique;
import common.ext;

public class CALiCo {
	public static final String[] SEX_SYNONYMS = {"sex", "male", "female", "gender"};
	public static final String[] SAMPLE_ID_SYNONYMS = {"ID", "IID"};
	public static final String CONDITIONALS_TXT_FILE = "conditionals.txt";

	public static byte[] getChrs(String genoFileDirPlusRoot, String[] markers, Logger log) {
		BufferedReader reader;
		String[] line;
		byte[] chrs;
		
		chrs = new byte[markers.length];
		for (int i=0; i<chrs.length; i++) {
			chrs[i] = (byte) -1;
		}
		try {
			reader = new BufferedReader(new FileReader(genoFileDirPlusRoot + ".bim"));
			while (reader.ready()) {
				line = reader.readLine().split("[\\s]+");
				for (int i=0; i<markers.length; i++) {
					if (markers[i].equals(line[1])) {
						if (chrs[i] == (byte) -1) {
							chrs[i] = Byte.parseByte(line[0]);
						} else {
							log.reportError(markers[i] + " has appeared more than once in " + genoFileDirPlusRoot + ".bim");
						}
					}
				}
			}
			reader.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return chrs;
	}




	public static void parseAllFilesInDirectory(String genoPlinkDirPlusRoot, String phenoCovarDir, String resultDir, String scratchDir, String plinkCommand, Logger log) {
		String[] models = null;
		String[] args;
		String[] uniqueMarkers;
		String[] fileParameters;
		byte[] chrs;
		String[] line;
//		Logger log;
		BufferedReader reader;
		String trav;
		Hashtable<String,Double> minimumPvalueHash;
		Hashtable<String,int[]> markerPositionHash;
		int[][] markerPositions;
		String[] genes;
		PrintWriter writer;
		int[] idVariable;
		boolean sexAsCovariate;
		Hashtable<String, String[]> conditionals = null;
		String root, subDir;
		String[] markersConditional;
		String condGenoPlinkDirPlusRoot;
		Logger[] logs;

		args = new String[] {null, null, "'MarkerName'", "'Chr'", "'Position'", "'P-value'", "!'P-value'<0.001", "tab", "replace=."};
		new File(resultDir + "topHits.xln").delete();
		if (Files.exists(phenoCovarDir + CONDITIONALS_TXT_FILE)) {
			//TODO delete all the old topHists.xln

			log.report("Performing conditional corvariate analysis based on file: " + phenoCovarDir + CONDITIONALS_TXT_FILE);
			conditionals = new Hashtable<String, String[]>();
			try {
				reader = new BufferedReader(new FileReader(phenoCovarDir + CONDITIONALS_TXT_FILE));
				while (reader.ready()) {
					line = reader.readLine().split("[\\s]+");
					root = ext.rootOf(line[0]);
					conditionals.put(root + "_" + ext.replaceWithLinuxSafeCharacters(line[1], false), new String[] {root, line[1]});
				}
				reader.close();
			} catch (FileNotFoundException e1) {
				e1.printStackTrace();
			} catch (IOException e2) {
				e2.printStackTrace();
			}

			models = conditionals.keySet().toArray(new String[0]);
			logs = new Logger[models.length];
			line = new String[models.length];
			for (int i = 0; i < line.length; i ++) {
				line[i] = conditionals.get(models[i])[1];
			}
			chrs = getChrs(genoPlinkDirPlusRoot, line, log);

			idVariable = new int[1];
			for (int i=0; i<models.length; i++) {
				root = conditionals.get(models[i])[0];
				subDir = resultDir + models[i] + "/";
				if (!new File(resultDir + models[i] + ".out").exists()) {
					new File(subDir).mkdir();
					logs[i] = new Logger(subDir+"Genvisis_CALiCo_" + (new SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date())) + ".log");;
					
					args[0] = resultDir + models[i] + ".out";
					args[1] = "out=" + subDir + models[i] + "_hits.txt";
	
					markersConditional = conditionals.get(models[i])[1].split(",");
					if (markersConditional.length==1) {
						if (!new File(scratchDir + "plink_chr" + chrs[i] + ".bim").exists() || !new File(scratchDir + "plink_chr" + chrs[i] + ".bim").exists()) {
							CmdLine.run("plink --bfile " + genoPlinkDirPlusRoot + " --chr " + chrs[i] + " --make-bed --out plink_chr"+chrs[i], scratchDir);
							log.report("Generated new set of genotype files: " + scratchDir + "plink_chr" + chrs[i] + ".");
						} else {
							log.reportError("Warning --- Found existing set of genotype files: " + scratchDir + "plink_chr" + chrs[i] + ", and will use them for the analysis followed.");
						}
						condGenoPlinkDirPlusRoot = scratchDir + "plink_chr" + chrs[i];
					} else {
						condGenoPlinkDirPlusRoot = genoPlinkDirPlusRoot;
					}
					Files.writeList(markersConditional, resultDir + models[i] + ".txt");

					idVariable[0] = ext.indexOfStr(Files.getHeaderOfFile(phenoCovarDir + root + ".xln", null)[0], SAMPLE_ID_SYNONYMS);
					sexAsCovariate = parsePhenotypes(condGenoPlinkDirPlusRoot, phenoCovarDir + root + ".xln", SAMPLE_ID_SYNONYMS[idVariable[0]], subDir, logs[i]);
					Conditional.addCountsAsCovariate(resultDir, models[i] + "/", condGenoPlinkDirPlusRoot, models[i] + "/" + root + "_covars.dat", models[i] + "_covars.dat", models[i] + ".txt", logs[i]);
					new File(resultDir + models[i] + ".txt").delete();
					new File(resultDir  + models[i] + "/" + root + "_covars.dat").delete();
	
					runAndParseResults(condGenoPlinkDirPlusRoot, resultDir + models[i] + "/" + root, resultDir + models[i] + "/" + models[i] + "_covars.dat", sexAsCovariate, scratchDir, resultDir + models[i] + "/", models[i], true, plinkCommand, logs[i]);
//					new File(resultDir + models[i] + "/" + models[i] + ".out").renameTo(new File(resultDir + models[i] + ".out"));
	
					GenParser.parse(args, logs[i]);

				} else {
					log.reportError("\nWarning --- Found '" + resultDir + models[i] + ".out' already exists. Unless you rename or delete the file or directory,"
									+ "the program does not regenerate this file and the '" + resultDir + models[i] + "/" + models[i] + "_hits.txt' file, "
									+ "and will generate the following files based on the latter: '" + resultDir + "cat_hits.txt', '" + resultDir + "all_hits.txt' and '" + resultDir + "topHits.xln'");
				}
			}

			args = new String[models.length];
			for (int i=0; i<models.length; i++) {
				args[i] = resultDir + models[i] + "/" + models[i] + "_hits.txt";
			}

		} else {
			models = Files.list(phenoCovarDir, ".xln", false);
			idVariable = new int[models.length];
			for (int i=0; i<models.length; i++) {
				idVariable[i] = ext.indexOfStr(Files.getHeaderOfFile(phenoCovarDir + models[i], null)[0], SAMPLE_ID_SYNONYMS);
				if (idVariable[i] < 0) {
					System.out.println("Removing " + models[i] + ".");
					Array.removeFromArray(models, i);
					i--;
				} else if(idVariable[i] == 1) {
					System.out.println("");
				} else {
					
				}
			}
			
			for (int i=0; i<models.length; i++) {
				root = ext.rootOf(models[i]);
				args[0] = resultDir + root + ".out";
				args[1] = "out=" + resultDir + root + "_hits.txt";
				if (!Files.exists(resultDir + root + ".out")) {
					parsePhenotypes(genoPlinkDirPlusRoot, phenoCovarDir + models[i], SAMPLE_ID_SYNONYMS[idVariable[i]], scratchDir, log);
					sexAsCovariate = parsePhenotypes(genoPlinkDirPlusRoot, phenoCovarDir + models[i], SAMPLE_ID_SYNONYMS[idVariable[i]], resultDir, log);
//					runAndParseResults(genoFileDirPlusRoot, resultDir, ext.rootOf(files[i]), sexAsCovariate, scratchDir, log);
					runAndParseResults(genoPlinkDirPlusRoot, resultDir + root, resultDir + root + "_covars.dat", sexAsCovariate, scratchDir, resultDir, root, false, plinkCommand, log);
				} else {
					log.reportError("\nWarning --- Found '" + resultDir + models[i] + ".out' already exists. Unless you rename or delete the file or directory,"
							+ "the program does not regenerate this file and the '" + resultDir + models[i] + "_hits.txt' file, "
							+ "and will generate the following files based on the latter: '" + resultDir + "cat_hits.txt', '" + resultDir + "all_hits.txt' and '" + resultDir + "topHits.xln'");
				}

				GenParser.parse(args, log);
			}

			args = new String[models.length];
			for (int i=0; i<models.length; i++) {
				args[i] = resultDir + ext.rootOf(models[i]) + "_hits.txt";
			}
		}

//		combine lists, only those that are unique
		Files.cat(args, resultDir+"cat_hits.txt", Array.intArray(models.length, 1), null);
		
		minimumPvalueHash = new Hashtable<String,Double>();
		markerPositionHash = new Hashtable<String,int[]>(); // key=markerName, values=new int[] {chr,position,position}
		//populate both hashtables
        try {
			reader = new BufferedReader(new FileReader(resultDir+"cat_hits.txt"));
			while((trav = reader.readLine()) != null) {
				line = trav.trim().split("\\t");
				if (!markerPositionHash.containsKey(line[0])) {
					minimumPvalueHash.put(line[0], Double.parseDouble(line[3]));
					markerPositionHash.put(line[0], new int[] {Integer.parseInt(line[1]), Integer.parseInt(line[2]), Integer.parseInt(line[2])});
				} else if (minimumPvalueHash.get(line[0]) > Double.parseDouble(line[3])) {
					minimumPvalueHash.put(line[0], Double.parseDouble(line[3]));
				}
			}
			reader.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		uniqueMarkers = HashVec.getKeys(minimumPvalueHash, false, false);
		
		
		// get map information
		markerPositions = new int[uniqueMarkers.length][3];
		for (int i=0; i<markerPositions.length; i++) {
			markerPositions[i] = markerPositionHash.get(uniqueMarkers[i]);
		}
		genes = MapSNPsAndGenes.mapSNPsToGenes(markerPositions, 0, log);
		
		double[] minPvalues = new double[uniqueMarkers.length];
        try {
			writer = new PrintWriter(new FileWriter(resultDir+"minimumPvalues.txt"));
			writer.println("MarkerName\tminPval\tGenes");
			for (int i = 0; i < minPvalues.length; i++) {
				minPvalues[i] = minimumPvalueHash.get(uniqueMarkers[i]);
				writer.println(uniqueMarkers[i]+"\t"+minPvalues[i]+"\t"+(genes[i].equals("")?".":genes[i]));
			}
			writer.close();
			uniqueMarkers = Sort.putInOrder(Sort.quicksort(minPvalues), uniqueMarkers);
		} catch (IOException e) {
			e.printStackTrace();
		}

		Files.writeList(uniqueMarkers, resultDir+"all_hits.txt");
		
		fileParameters = new String[models.length+3];
		for (int i=0; i<models.length; i++) {
			root = ext.rootOf(models[i]);
			if (i==0) {
//				if (conditionals == null || conditionals.size() == 0) {
					fileParameters[0] = resultDir + root + ".out 'MarkerName' 'Chr' 'Position' 'Effect_allele' 'Reference_allele'";
//				} else {
//					fileParameters[0] = resultDir + root + "/" + root + ".out 'MarkerName' 'Chr' 'Position' 'Effect_allele' 'Reference_allele'";
//				}
				fileParameters[1] = resultDir + "minimumPvalues.txt 0 2=Gene(s) tab";
			}
//			if (conditionals == null || conditionals.size() == 0) {
				fileParameters[i+2] = resultDir + root + ".out 'MarkerName' 'N' 'Effect_allele_frequency' 'BETA'=beta_" + root + " 'P-value'=pval_" + root;
//			} else {
//				fileParameters[i+2] = resultDir + root + "/" + root + ".out 'MarkerName' 'N' 'Effect_allele_frequency' 'BETA'=beta_" + root + " 'P-value'=pval_" + root;
//			}
		}
		fileParameters[models.length+2] = resultDir+"minimumPvalues.txt 0 1=minPval tab";
		
		Files.combine(uniqueMarkers, fileParameters, null, "MarkerName", ".", resultDir + "topHits.xln", log, true, true, false);
	}


//	public static void parseAllFilesInDirectory_2(String genoFiles, String phenoCovarDir, String resultDir, String scratchDir, Logger log) {
//		String[] files;
//		String[] args;
//		String[] uniqueMarkers;
//		String[] fileParameters;
//		byte[] chrs;
//		int[] positions;
//		String[] line;
////		Logger log;
//		BufferedReader reader;
//		String trav;
//		Hashtable<String,Double> minimumPvalueHash;
//		Hashtable<String,int[]> markerPositionHash;
//		int[][] markerPositions;
//		String[] genes;
//		PrintWriter writer;
//		int[] idVariable;
//		boolean sexAsCovariate;
//		
////		log = new Logger(phenoCovarDir+"CALiCo.log");
//		new File(phenoCovarDir+"topHits.xln").delete();
//
//		files = Files.list(phenoCovarDir, ".xln", false);
//		idVariable = new int[files.length];
//		for (int i=0; i<files.length; i++) {
//			idVariable[i] = ext.indexOfStr(Files.getHeaderOfFile(phenoCovarDir + files[i], null)[0], SAMPLE_ID_SYNONYMS);
//			if (idVariable[i] < 0) {
//				System.out.println("Removing " + files[i] + ".");
//				Array.removeFromArray(files, i);
//				i--;
//			} else if(idVariable[i] == 1) {
//				System.out.println("");
//			} else {
//				
//			}
//		}
//		
//		args = new String[] {null, null, "'MarkerName'", "'Chr'", "'Position'", "'P-value'", "!'P-value'<0.001", "tab", "replace=."};
//		for (int i=0; i<files.length; i++) {
//			if (!Files.exists(phenoCovarDir+ext.rootOf(files[i])+".out")) {
//				sexAsCovariate = parsePhenotypes(genoFiles, phenoCovarDir+files[i], SAMPLE_ID_SYNONYMS[idVariable[i]], phenoCovarDir, log);
//				runAndParseResults(resultDir, ext.rootOf(files[i]), sexAsCovariate, genoFiles, scratchDir, log);
//			}
//			args[0] = phenoCovarDir+ext.rootOf(files[i])+".out";
//			args[1] = "out="+phenoCovarDir+ext.rootOf(files[i])+"_hits.txt";
//			GenParser.parse(args, log);
//		}
//		
////		combine lists, only those that are unique
//		args = new String[files.length];
//		for (int i=0; i<files.length; i++) {
//			args[i] = phenoCovarDir+ext.rootOf(files[i])+"_hits.txt";
//		}
//		Files.cat(args, phenoCovarDir+"cat_hits.txt", Array.intArray(files.length, 1), null);
//		
//		minimumPvalueHash = new Hashtable<String,Double>();
//		markerPositionHash = new Hashtable<String,int[]>(); // key=markerName, values=new int[] {chr,position,position}
//		//populate both hashtables
//        try {
//			reader = new BufferedReader(new FileReader(phenoCovarDir+"cat_hits.txt"));
//			while((trav = reader.readLine()) != null) {
//				line = trav.trim().split("\\t");
//				if (line[1].equals("7") && line[2].equals("44201262")) {
//					System.out.println();
//				}
//				if (!markerPositionHash.containsKey(line[0])) {
//					minimumPvalueHash.put(line[0], Double.parseDouble(line[3]));
//					markerPositionHash.put(line[0], new int[] {Integer.parseInt(line[1]), Integer.parseInt(line[2]), Integer.parseInt(line[2])});
//				} else if (minimumPvalueHash.get(line[0]) > Double.parseDouble(line[3])) {
//					minimumPvalueHash.put(line[0], Double.parseDouble(line[3]));
//				}
//			}
//			reader.close();
//		} catch (FileNotFoundException e) {
//			e.printStackTrace();
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//		uniqueMarkers = HashVec.getKeys(minimumPvalueHash, false, false);
//		
//		
//		// get map information
//		markerPositions = new int[uniqueMarkers.length][3];
//		for (int i=0; i<markerPositions.length; i++) {
//			markerPositions[i] = markerPositionHash.get(uniqueMarkers[i]);
//		}
//		genes = MapSNPsAndGenes.mapSNPsToGenes(markerPositions, 0, log);
//		
//		double[] minPvalues = new double[uniqueMarkers.length];
//        try {
//			writer = new PrintWriter(new FileWriter(phenoCovarDir+"minimumPvalues.txt"));
//			writer.println("MarkerName\tminPval\tGenes");
//			for (int i = 0; i < minPvalues.length; i++) {
//				minPvalues[i] = minimumPvalueHash.get(uniqueMarkers[i]);
//				writer.println(uniqueMarkers[i]+"\t"+minPvalues[i]+"\t"+(genes[i].equals("")?".":genes[i]));
//			}
//			writer.close();
//			uniqueMarkers = Sort.putInOrder(Sort.quicksort(minPvalues), uniqueMarkers);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
//
//
////		uniqueMarkers = HashVec.loadFileToStringArray(dir+"cat_hits.txt", false, new int[] {0,1,2}, true);
////		chrs = new byte[uniqueMarkers.length];
////		positions = new int[uniqueMarkers.length];
////		for (int i = 0; i < uniqueMarkers.length; i++) {
////			line = uniqueMarkers[i].split("\t");
////			chrs[i] = Positions.chromosomeNumber(line[1]);
////			positions[i] = Integer.parseInt(line[2]);
////			uniqueMarkers[i] = line[0];
////		}
////		uniqueMarkers = Sort.putInOrder(uniqueMarkers, Sort.orderTwoLayers(chrs, positions));
//		Files.writeList(uniqueMarkers, phenoCovarDir+"all_hits.txt");
//		
//		fileParameters = new String[files.length+3];
//		for (int i=0; i<files.length; i++) {
//			String root = ext.rootOf(files[i]);
//			if (i==0) {
//				fileParameters[0] = phenoCovarDir+root+".out 'MarkerName' 'Chr' 'Position' 'Effect_allele' 'Reference_allele'";
//				fileParameters[1] = phenoCovarDir+"minimumPvalues.txt 0 2=Gene(s) tab";
//			}
//			fileParameters[i+2] = phenoCovarDir+root+".out 'MarkerName' 'N' 'Effect_allele_frequency' 'BETA'=beta_"+root+" 'P-value'=pval_"+root;
//		}
//		fileParameters[files.length+2] = phenoCovarDir+"minimumPvalues.txt 0 1=minPval tab";
//		
//		Files.combine(uniqueMarkers, fileParameters, null, "MarkerName", phenoCovarDir+"topHits.xln", log, true, true, false, false);
//
//	}








//	public static void parseAllFilesInDirectory_1(String dir, String locationOfGenotypeFiles, String scratchDir) {
//		String[] files;
//		String[] args;
//		String[] uniqueMarkers;
//		String[] fileParameters;
//		byte[] chrs;
//		int[] positions;
//		String[] line;
//		Logger log;
//		
//		log = new Logger(dir+"CALiCo.log");
//
//		files = Files.list(dir, ".xln", false);
//
//		args = new String[] {null, null, "'MarkerName'", "'Chr'", "'Position'", "'P-value'", "!'P-value'<0.001", "tab", "replace=."};
//		for (int i=0; i<files.length; i++) {
//			if (!Files.exists(dir+ext.rootOf(files[i])+".out")) {
//				parsePhenotypes(dir, dir+files[i], SAMPLE_ID_SYNONYMS[0], locationOfGenotypeFiles, scratchDir);
//			}
//			args[0] = dir+ext.rootOf(files[i])+".out";
//			args[1] = "out="+dir+ext.rootOf(files[i])+"_hits.txt";
//			GenParser.parse(args, log);
//		}
//		
////		combine lists, only those that are unique
//		args = new String[files.length];
//		for (int i=0; i<files.length; i++) {
//			args[i] = dir+ext.rootOf(files[i])+"_hits.txt";
//		}
//		Files.cat(args, dir+"cat_hits.txt", Array.intArray(files.length, 1), null);
//		
//		uniqueMarkers = HashVec.loadFileToStringArray(dir+"cat_hits.txt", false, new int[] {0,1,2}, true);
//		chrs = new byte[uniqueMarkers.length];
//		positions = new int[uniqueMarkers.length];
//		for (int i = 0; i < uniqueMarkers.length; i++) {
//			line = uniqueMarkers[i].split("\t");
//			chrs[i] = Positions.chromosomeNumber(line[1]);
//			positions[i] = Integer.parseInt(line[2]);
//			uniqueMarkers[i] = line[0];
//		}
//		uniqueMarkers = Sort.putInOrder(uniqueMarkers, Sort.orderTwoLayers(chrs, positions));
//		Files.writeList(uniqueMarkers, dir+"all_hits.txt");
//		
//		fileParameters = new String[files.length];
//		for (int i=0; i<files.length; i++) {
//			String root = ext.rootOf(files[i]);
//			if (i==0) {
//				fileParameters[i] = dir+root+".out 'MarkerName' 'Chr' 'Position' 'Effect_allele' 'Reference_allele' 'N' 'Effect_allele_frequency' 'BETA'=beta_"+root+" 'P-value'=pval_"+root;
//			} else {
//				fileParameters[i] = dir+root+".out 'MarkerName' 'N' 'Effect_allele_frequency' 'BETA'=beta_"+root+" 'P-value'=pval_"+root;
//			}
//		}
//		
//		Files.combine(uniqueMarkers, fileParameters, null, "MarkerName", dir+"topHits.xln", log, true, true, true, false);
//
//	}









//	public static void parsePhenotypes(String genoFiles, String phenoCovarFilename, String idVarName, String scratchDir, Logger log) {
//		BufferedReader reader;
//		PrintWriter writer;
////		Logger log;
//		String[] header;
//		String temp;
//		String[] args;
//		Vector<String> v;
//		String root;
//		String[] phenos;
//		String[] line;
//		boolean isBinary;
//
//		root = ext.rootOf(phenoCovarFilename);
////		analysisDir = ext.parseDirectoryOfFile(filename);
////		log = new Logger(resultDir + root + ".log");
//		temp = Files.getFirstNLinesOfFile(phenoCovarFilename, 1, log)[0];
//		header = temp.trim().split("[\\s]+");
//
//		args = new String[] {phenoCovarFilename, "out="+scratchDir+root+"_pheno.dat", "'"+idVarName+"'=FID", "'"+idVarName+"'=IID", "'"+header[1]+"'", "tab", "replace=."};
//		GenParser.parse(args, log);
//		
//		// load phenotype to memory
//		phenos = HashVec.loadFileToStringArray(scratchDir+root+"_pheno.dat", true, new int[] {2}, false);
//		isBinary = RegressionModel.isBinaryTrait(phenos, log);
//		if (isBinary) {
//			if (Array.max(Array.toIntArray(phenos)) == 1 && Array.min(Array.toIntArray(phenos)) == 0) {
//				new File(scratchDir+root+"_pheno.dat").renameTo(new File(scratchDir+root+"_pheno.dat.temp"));
//				
//				try {
//					reader = new BufferedReader(new FileReader(scratchDir+root+"_pheno.dat.temp"));
//					writer = new PrintWriter(new FileWriter(scratchDir+root+"_pheno.dat"));
//					//header
//					writer.println(reader.readLine());
//					while (reader.ready()) {
//						// 1 -> 2
//						// 0 -> 1
//						// anything else becomes missing; if value is not among ext.MISSING_VALUES, then report error to log				
//						line = reader.readLine().split("\t");
//						if (line[2].equals("1")) {
//							line[2] = "2";
//						} else if (line[2].equals("0")) {
//							line[2] = "1";
//						} else {
//							for (int i=0; i<ext.MISSING_VALUES.length; i++) {
//								if (ext.MISSING_VALUES[i].contains(line[2])) {
//									log.reportError("Unrecognized pheno type: "+line[2]);
//								}
//							}
//							line[2] = "";
//						}
//						writer.println(line[0]+"\t"+line[1]+"\t"+line[2]);
//					}
//					writer.close();
//					reader.close();
//				} catch (FileNotFoundException fnfe) {
//					System.err.println("Error: file \"" + scratchDir+root+"_pheno.dat.temp" + "\" not found in current directory");
//					System.exit(1);
//				} catch (IOException ioe) {
//					System.err.println("Error reading file \"" + scratchDir+root+"_pheno.dat.temp" + "\"");
//					System.exit(2);
//				}
//			}
//		}
//	}
//
//	public static boolean parseCovars(String genoFiles, String phenoCovarFilename, String idVarName, String scratchDir, Logger log) {
////		Logger log;
//		String[] header;
//		String temp;
//		String[] args;
//		boolean sexAsCovariate;
//		Hashtable<String, String> hashPheno, hashCovariates;
//		String[] keys, covars;
////		String analysisDir;
//		ArrayList<Integer> keysToRemove;
//		Vector<String> v;
//		String root;
//		String[] line;
//
//		root = ext.rootOf(phenoCovarFilename);
////		analysisDir = ext.parseDirectoryOfFile(filename);
////		log = new Logger(resultDir + root + ".log");
//		temp = Files.getFirstNLinesOfFile(phenoCovarFilename, 1, log)[0];
//		header = temp.trim().split("[\\s]+");
//
//		args = new String[] {phenoCovarFilename, "out="+scratchDir+root+"_pheno.dat", "'"+idVarName+"'=FID", "'"+idVarName+"'=IID", "'"+header[1]+"'", "tab", "replace=."};
//		GenParser.parse(args, log);
//
//		sexAsCovariate = false;
//		v = new Vector<String>();
//		v.add(phenoCovarFilename);
//		v.add("out="+scratchDir+root+"_covars.dat");
//		v.add("'ID'=FID");
//		v.add("'ID'=IID");
//		for (int i=2; i<header.length; i++) {
//			if (ext.indexOfStr(header[i], SEX_SYNONYMS, false, true, log, false) >= 0) {
//				sexAsCovariate = true;
//			} else {
//				v.add ("'"+header[i]+"'");
//			}
//		}
//		v.add("tab");
//		v.add("replace=.");
//		args = Array.toStringArray(v);
//		GenParser.parse(args, log);
//
////		determine which samples have complete data for pheotype AND all covariates
//		keysToRemove = new ArrayList<Integer>();
//		hashPheno = HashVec.loadFileToHashString(scratchDir+root+"_pheno.dat", new int[] {0,1}, new int[] {2}, false, "\t", true, false, false);
//		hashCovariates = HashVec.loadFileToHashString(scratchDir+root+"_covars.dat", new int[] {0,1}, Array.subArray(Array.intArray(args.length-6), 2), false, "\t", true, false, false);
//		keys = HashVec.getKeys(hashPheno);
//		for (int i = 0; i < keys.length; i++) {
//			if (ext.isMissingValue(hashPheno.get(keys[i]))) {
//				keysToRemove.add(i);
//			} else if (hashCovariates.containsKey(keys[i])) {
//				covars = hashCovariates.get(keys[i]).trim().split("[\\s]+");
//				for (int j = 0; j < covars.length; j++) {
//					if (ext.isMissingValue(covars[j])) {
//						keysToRemove.add(i);
//					}
//				}
//			} else {
//				keysToRemove.add(i);
//			}
//		}
//		for (int i=0; i<keysToRemove.size(); i++) {
//			hashPheno.remove(keys[keysToRemove.get(i)]);
//		}
//		Files.writeList(HashVec.getKeys(hashPheno), scratchDir+root+"_used.dat");
//		log.report("There were "+hashPheno.size()+" samples with complete phenotypic data");
//
//		
////		write batch (run1.bat) using locationOfGenotypeFiles and scratchDir
////		add sex if necessary
//		return sexAsCovariate;
//	}
	









//	public static boolean parsePhenotypes(String locationOfGenotypeFiles, String phenoCovarFilename, String idVarName, String resultDir, Logger log) {
	public static boolean parsePhenotypes(String genoFiles, String phenoCovarFilename, String idVarName, String outputDir, Logger log) {
		BufferedReader reader;
		PrintWriter writer;
//		Logger log;
		String[] header;
		String temp;
		String[] args;
		boolean sexAsCovariate;
		Hashtable<String, String> hashPheno, hashCovariates;
		String[] keys, covars;
//		String analysisDir;
		ArrayList<Integer> keysToRemove;
		Vector<String> v;
		String root;
		String[] phenos;
		String[] line;
		boolean isBinary;

		root = ext.rootOf(phenoCovarFilename);
//		analysisDir = ext.parseDirectoryOfFile(filename);
//		log = new Logger(resultDir + root + ".log");
		temp = Files.getFirstNLinesOfFile(phenoCovarFilename, 1, log)[0];
		header = temp.trim().split("[\\s]+");

		args = new String[] {phenoCovarFilename, "out=" + outputDir + root + "_pheno.dat", "'" + idVarName + "'=FID", "'" + idVarName + "'=IID", "'" + header[1] + "'", "tab", "replace=."};
		GenParser.parse(args, log);
		
		// load phenotype to memory
		phenos = HashVec.loadFileToStringArray(outputDir+root+"_pheno.dat", true, new int[] {2}, false);
		isBinary = RegressionModel.isBinaryTrait(phenos, log);
		if (isBinary) {
			if (Array.max(Array.toIntArray(phenos)) == 1 && Array.min(Array.toIntArray(phenos)) == 0) {
				new File(outputDir + root + "_pheno.dat").renameTo(new File(outputDir + root + "_pheno.dat.temp"));
				
				try {
					reader = new BufferedReader(new FileReader(outputDir + root + "_pheno.dat.temp"));
					writer = new PrintWriter(new FileWriter(outputDir + root + "_pheno.dat"));
					//header
					writer.println(reader.readLine());
					while (reader.ready()) {
						// 1 -> 2
						// 0 -> 1
						// anything else becomes missing; if value is not among ext.MISSING_VALUES, then report error to log				
						line = reader.readLine().split("\t");
						if (line[2].equals("1")) {
							line[2] = "2";
						} else if (line[2].equals("0")) {
							line[2] = "1";
						} else {
							for (int i=0; i<ext.MISSING_VALUES.length; i++) {
								if (ext.MISSING_VALUES[i].contains(line[2])) {
									log.reportError("Unrecognized pheno type: "+line[2]);
								}
							}
							line[2] = "";
						}
						writer.println(line[0]+"\t"+line[1]+"\t"+line[2]);
					}
					writer.close();
					reader.close();
				} catch (FileNotFoundException fnfe) {
					System.err.println("Error: file \"" + outputDir + root + "_pheno.dat.temp" + "\" not found in current directory");
					System.exit(1);
				} catch (IOException ioe) {
					System.err.println("Error reading file \"" + outputDir + root + "_pheno.dat.temp" + "\"");
					System.exit(2);
				}
			}
		}
		
//		leave out phenotype
//		leave out sex if it exists (make a flag to use sex)
//		at end, add tab and replace

//		sexAsCovariate = false;
//		args = new String[header.length+3];
//		args[0] = filename;
//		args[1] = "out="+dir+"covars1.dat";
//		args[2] = "'ID'=FID";
//		args[3] = "'ID'=IID";
//		for (int i=3; i<header.length; i++) {
//			args[i+1] = "'"+header[i]+"'";
//		}
//		args[args.length-2] = "tab";
//		args[args.length-1] = "replace=.";
		
		sexAsCovariate = false;
		v = new Vector<String>();
		v.add(phenoCovarFilename);
		v.add("out=" + outputDir + root + "_covars.dat");
		v.add("'IID'=FID");
		v.add("'IID'=IID");
		for (int i=2; i<header.length; i++) {
			if (ext.indexOfStr(header[i], SEX_SYNONYMS, false, true, log, false) >= 0) {
				sexAsCovariate = true;
			} else {
				v.add ("'"+header[i]+"'");
			}
		}
		v.add("tab");
		v.add("replace=.");
		args = Array.toStringArray(v);
		GenParser.parse(args, log);
		
//		determine which samples have complete data for pheotype AND all covariates
		keysToRemove = new ArrayList<Integer>();
		hashPheno = HashVec.loadFileToHashString(outputDir + root + "_pheno.dat", new int[] {0,1}, new int[] {2}, false, "\t", true, false, false);
		hashCovariates = HashVec.loadFileToHashString(outputDir + root + "_covars.dat", new int[] {0,1}, Array.subArray(Array.intArray(args.length-6), 2), false, "\t", true, false, false);
		keys = HashVec.getKeys(hashPheno);
		for (int i = 0; i < keys.length; i++) {
			if (ext.isMissingValue(hashPheno.get(keys[i]))) {
				keysToRemove.add(i);
			} else if (hashCovariates.containsKey(keys[i])) {
				covars = hashCovariates.get(keys[i]).trim().split("[\\s]+");
				for (int j = 0; j < covars.length; j++) {
					if (ext.isMissingValue(covars[j])) {
						keysToRemove.add(i);
					}
				}
			} else {
				keysToRemove.add(i);
			}
		}
		for (int i=0; i<keysToRemove.size(); i++) {
			hashPheno.remove(keys[keysToRemove.get(i)]);
		}
		Files.writeList(HashVec.getKeys(hashPheno), outputDir + root + "_used.dat");
		log.report("There were " + hashPheno.size() + " samples with complete phenotypic data");

		
//		write batch (run1.bat) using locationOfGenotypeFiles and scratchDir
//		add sex if necessary
		
		return sexAsCovariate;
	}









//	public static void runAndParseResults(String genoFiles, String resultDir, String root, boolean sexAsCovariate, String scratchDir, Logger log) {
//		String[] commands;
//		
//		if (scratchDir == null) {
//			scratchDir = resultDir;
//		}
//		
//		commands = new String[] {//"Path = %Path%; C:\plink",
//				  "C:/PLINK/plink --bfile N:/statgen/CALiCo/filteredGenotypes/plink --pheno " + scratchDir + root + "_pheno.dat --covar " + scratchDir + root + "_covars.dat" + (sexAsCovariate?" --sex":"") + " --logistic --ci 0.95 --out " + scratchDir + root,
//				  "C:/PLINK/plink --bfile N:/statgen/CALiCo/filteredGenotypes/plink --keep " + scratchDir + root + "_used.dat --freq --out " + scratchDir + root + "_freq"
//				  };
//		Files.writeList(commands, scratchDir + "run_" + root + ".bat");
//		
//		//run batch
//		CmdLine.run(scratchDir + "run_" + root + ".bat", "");
//		
//		//parse results
//		if (new File(scratchDir + root + ".assoc.logistic").exists()) {
//			ResultsPackager.parseStdFormatFromPlink("", scratchDir + root + ".assoc.logistic", "ADD", genoFiles+".bim", scratchDir + root + "_freq.frq", null, resultDir + root + ".out", log);
//		} else {
//			ResultsPackager.parseStdFormatFromPlink("", scratchDir + root + ".assoc.linear", "ADD", genoFiles+".bim", scratchDir + root + "_freq.frq", null, resultDir + root + ".out", log);
//		}
//	}

	public static void runAndParseResults(String genoFileDirPlusRoot, String phenoFileDirPlusRoot, String covarFileDirPlusName, boolean sexAsCovariate, String scratchDir, String resultDir, String outFileRoot, boolean isConditional, String plinkCommand, Logger log) {
		String[] commands;
		String plinkResultFile;

		if(scratchDir==null) {
			scratchDir = resultDir;
		}

		if (plinkCommand == null) {
			plinkCommand = "plink";
		}

		commands = new String[] {//"Path = %Path%; C:\plink",
					plinkCommand + " --bfile " + genoFileDirPlusRoot + " --pheno " + phenoFileDirPlusRoot + "_pheno.dat --covar " + covarFileDirPlusName + (sexAsCovariate?" --sex":"") + " --logistic --ci 0.95 --out " + scratchDir + outFileRoot,
					plinkCommand + " --bfile " + genoFileDirPlusRoot + " --keep " + phenoFileDirPlusRoot + "_used.dat --freq --out " + scratchDir + outFileRoot + "_freq"
				  };
		Files.writeList(commands, resultDir + "run_" + outFileRoot + ".bat");
		Files.chmod(resultDir + "run_" + outFileRoot + ".bat");
		
		//run batch
		CmdLine.run(resultDir + "run_" + outFileRoot + ".bat", "");
		
		//parse results
//		if (new File(scratchDir + outFileRoot + ".assoc.logistic").exists()) {
//			ResultsPackager.parseStdFormatFromPlink("", scratchDir + outFileRoot + ".assoc.logistic", "ADD", genoFileDirPlusRoot+".bim", scratchDir + outFileRoot + "_freq.frq", null, 1, resultDir + outFileRoot + ".out", log);
//			ResultsPackager.parseStdFormatFromPlink("", scratchDir + outFileRoot + ".assoc.logistic", outFileRoot.replace("_", ":"), genoFileDirPlusRoot+".bim", scratchDir + outFileRoot + "_freq.frq", null, 1, resultDir + outFileRoot + "_cov.out", log);
//		} else {
//			ResultsPackager.parseStdFormatFromPlink("", scratchDir + outFileRoot + ".assoc.linear", "ADD", genoFileDirPlusRoot+".bim", scratchDir + outFileRoot + "_freq.frq", null, 1, resultDir + outFileRoot + ".out", log);
//			ResultsPackager.parseStdFormatFromPlink("", scratchDir + outFileRoot + ".assoc.linear", outFileRoot.replace("_", ":"), genoFileDirPlusRoot+".bim", scratchDir + outFileRoot + "_freq.frq", null, 1, resultDir + outFileRoot + "_cov.out", log);
//		}

		plinkResultFile = scratchDir + outFileRoot + ".assoc.linear";
		if (! new File(plinkResultFile).exists()) {
			plinkResultFile = scratchDir + outFileRoot + ".assoc.logistic";
		}

		if (isConditional) {
			resultDir = resultDir.substring(0, resultDir.lastIndexOf("/", resultDir.length()-2)+1); 
		}
		ResultsPackager.parseStdFormatFromPlink("", plinkResultFile, "ADD", genoFileDirPlusRoot+".bim", scratchDir + outFileRoot + "_freq.frq", null, 1, resultDir + outFileRoot + ".out", log);
		ResultsPackager.parseStdFormatFromPlink("", plinkResultFile, outFileRoot.substring(outFileRoot.indexOf("_") + 1).replace("_", ":"), genoFileDirPlusRoot+".bim", scratchDir + outFileRoot + "_freq.frq", null, 1, resultDir + outFileRoot + "_cov.out", log);
	}

	public static void parseResultsForPageGlucoseInsulinPaperFormat(String modelListFileName, String resultFilesDir, String outFileName) {
		String filename, model, snp, outLine, delimiterForOutputFile;
		String[] model1, modelLine, snpLine, indexLine, line;
		Logger log;
		Hashtable<String, String[]> modelList, resultsSnp, resultsIndex;
		Enumeration<String> models, snps;
		Vector<String> outFile;
		int numSamples;
		
		log = new Logger();

		if (outFileName.endsWith(".csv")) {
			delimiterForOutputFile = ",";
		} else {
			delimiterForOutputFile = "\t";
		}

		outFile = new Vector<String>();
		outFile.add("model" + delimiterForOutputFile + "SNP.original" + delimiterForOutputFile + "rsID" + delimiterForOutputFile + "CHR.build36" + delimiterForOutputFile + "CHR.build37" + delimiterForOutputFile + "bp.build36" + delimiterForOutputFile + "bp.build37" + delimiterForOutputFile + "analysis.locus" + delimiterForOutputFile  + "analysis.race" + delimiterForOutputFile+ "index.rsID" + delimiterForOutputFile + "snp.coded" + delimiterForOutputFile + "snp.noncoded" + delimiterForOutputFile + "snp.BETA" + delimiterForOutputFile + "snp.SE" + delimiterForOutputFile + "snp.STAT" + delimiterForOutputFile + "snp.P" + delimiterForOutputFile + "snp.L95" + delimiterForOutputFile + "snp.U95" + delimiterForOutputFile + "snp.NMISS" + delimiterForOutputFile + "snp.CAF" + delimiterForOutputFile + "snp.CAC" + delimiterForOutputFile + "index.BETA" + delimiterForOutputFile + "index.SE" + delimiterForOutputFile + "index.STAT" + delimiterForOutputFile + "index.P" + delimiterForOutputFile + "index.L95" + delimiterForOutputFile + "index.U95" + delimiterForOutputFile + "index.NMISS" + delimiterForOutputFile + "index.CAF" + delimiterForOutputFile + "index.CAC" + delimiterForOutputFile + "index.coded" + delimiterForOutputFile + "index.noncoded");
		modelList = SkatMeta.loadFile(modelListFileName, null, new String[] {"Name", "rsID"}, new String[] {"Locus", "Chr", "bp36.start", "bp36.stop"}, null, true, log);
		models = modelList.keys();
		while (models.hasMoreElements()) {
			model = models.nextElement();
			model1 = model.replaceAll(" ", "").split("\t");
//			if (!model1[1].startsWith("rs")) {
//				filename = model1[0] + "_rs" + model1[1];
//			} else {
				filename = model1[0] + "_" + ext.replaceWithLinuxSafeCharacters(model1[1], false);
//			}
			modelLine = modelList.get(model);
			resultsSnp = SkatMeta.loadFile(resultFilesDir + filename + ".out", null, new String[] {"MarkerName", "Chr", "Position"}, new String[] {"Effect_allele", "Reference_allele", "Effect_allele_frequency", "N", "BETA", "SE", "P-value"}, new String[] {"Chr ==" + modelLine[1], "Position >= " + modelLine[2], "Position <= " + modelLine[3]}, log);
			resultsIndex = SkatMeta.loadFile(resultFilesDir + filename + "_cov.out", null, new String[] {"MarkerName", "Chr", "Position"}, new String[] {"Effect_allele", "Reference_allele", "Effect_allele_frequency", "N", "BETA", "SE", "P-value"}, new String[] {"Chr ==" + modelLine[1], "Position >= " + modelLine[2], "Position <= " + modelLine[3]}, log);
			numSamples = Files.countLines(resultFilesDir + filename + "/" + model1[0] + "_pheno.dat", true);
			snps = resultsSnp.keys();
			while (snps.hasMoreElements()) {
				snp = snps.nextElement();
				snpLine = resultsSnp.get(snp);
				indexLine = resultsIndex.get(snp);
				line = snp.split("\t");
				outLine = model1[0]			//model
						+ delimiterForOutputFile				//SNP.original
						+ delimiterForOutputFile + line[0]	//rsID
						+ delimiterForOutputFile + line[1]	//CHR.build36
						+ delimiterForOutputFile				//CHR.build37
						+ delimiterForOutputFile + line[2]	//bp.build36
						+ delimiterForOutputFile				//bp.build37
						+ delimiterForOutputFile + modelLine[0]	//analysis.locus
						+ delimiterForOutputFile + "AA"			//analysis.race
						+ delimiterForOutputFile + model1[1]	//index.rsID
						+ delimiterForOutputFile + snpLine[0]	//snp.coded
						+ delimiterForOutputFile + snpLine[1]	//snp.noncoded
						+ delimiterForOutputFile + snpLine[4]	//snp.BETA
						+ delimiterForOutputFile + snpLine[5]	//snp.SE
						+ delimiterForOutputFile 				//snp.Stat
						+ delimiterForOutputFile + snpLine[6]	//snp.P
						+ ((snpLine[4].equals("NA") || snpLine[4].equals("NA"))?
								(delimiterForOutputFile + delimiterForOutputFile) :
								(delimiterForOutputFile + (Double.parseDouble(snpLine[4]) - 1.95 * Double.parseDouble(snpLine[5]))
								+ delimiterForOutputFile + (Double.parseDouble(snpLine[4]) + 1.95 * Double.parseDouble(snpLine[5]))))	//snp.L95	snp.U95
						+ delimiterForOutputFile + (numSamples - Integer.parseInt(snpLine[3]))	//snp.NMISS
						+ delimiterForOutputFile + snpLine[2]	//snp.CAF
						+ delimiterForOutputFile + (snpLine[2].equals("NA")? "" : (Double.parseDouble(snpLine[2]) * Integer.parseInt(snpLine[3]) * 2))	//snp.CAC
						+ delimiterForOutputFile + (indexLine == null || indexLine.length < 5 ? "NA" : indexLine[4])	//index.BETA
						+ delimiterForOutputFile + (indexLine == null || indexLine.length < 5 ? "NA" : indexLine[5])	//index.SE
						+ delimiterForOutputFile	//index.Stat
						+ delimiterForOutputFile + (indexLine == null || indexLine.length < 5 ? "NA" : indexLine[6])	//index.P
						+ ((indexLine == null || indexLine.length < 6 || indexLine[4].equals("NA")) ?
								(delimiterForOutputFile + delimiterForOutputFile) :
								(delimiterForOutputFile + (Double.parseDouble(indexLine[4]) - 1.95 * Double.parseDouble(indexLine[5]))
								+ delimiterForOutputFile + (Double.parseDouble(indexLine[4]) + 1.95 * Double.parseDouble(indexLine[5]))))	//index.L95	index.U95
						+ delimiterForOutputFile + (indexLine == null || indexLine.length < 5 ? "NA" : (numSamples - Integer.parseInt(indexLine[3])))	//index.NMISS
						+ delimiterForOutputFile + (indexLine == null || indexLine.length < 5 ? "NA" : indexLine[2]	//index.CAF
						+ delimiterForOutputFile + (indexLine == null || indexLine.length < 4 || indexLine[2].equals("NA")? "" : (Double.parseDouble(indexLine[2]) * Integer.parseInt(indexLine[3]) * 2)))	//index.CAC
						+ delimiterForOutputFile + (indexLine == null || indexLine.length < 5 ? "NA" : indexLine[0])	//index.coded
						+ delimiterForOutputFile + (indexLine == null || indexLine.length < 5 ? "NA" : indexLine[1])	//index.noncoded
						;
				outFile.add(outLine);
			}
		}
		Files.writeList(outFile.toArray(new String[0]), outFileName);
		log.report("Output is ready at: " + outFileName + "\nPlease use other software to sort by model ID, index.rsID, and positions.");
	}

	public static void metaAnalyzeSOL(String dir, String filePattern, String root, String mapFile) {
		String filename;
		Vector<String> v;
		int count;

		String[] sites = new String[] {"b", "c", "m", "s"};
		
		for (int i = 0; i < sites.length; i++) {
			if (Files.exists(dir+root+sites[i]+".out")) {
				System.out.println(root+sites[i]+".out already exists");
			} else {
				System.out.println("Creating "+root+sites[i]+".out");
				count = 1;
				v = new Vector<String>();
				filename = ext.replaceAllWith(filePattern, new String[][] {{"#", count+""}, {"%%", sites[i]}});
				while (Files.exists(dir+filename)) {
					v.add(dir+filename);
					count++;
					filename = ext.replaceAllWith(filePattern, new String[][] {{"#", count+""}, {"%%", sites[i]}});
				}
				Files.cat(Array.toStringArray(v), dir+root+sites[i]+".out", new int[0], new Logger());
			}
		}
		
		metaAnalyzeSOL(dir, root, mapFile);
	}
	
	public static void metaAnalyzeSOL(String dir, String root, String mapFile) {
		String filename, freqFile, outfile;
		Logger log;
		String[] inputFiles;
		Vector<String> qqFiles, lowCallrateMarkerFiles;
		
		log = new Logger();
		freqFile = null;
		dir = ext.verifyDirFormat(dir);
		
		String[] sites = new String[] {"b", "c", "m", "s"};
		
		qqFiles = new Vector<String>();
		lowCallrateMarkerFiles = new Vector<String>();
		inputFiles = new String[0];
		for (int i = 0; i < sites.length; i++) {
			filename = root+sites[i]+".out";
			outfile = root+sites[i]+"_parsed.out";
			if (Files.exists(dir+outfile)) {
				log.report(outfile+" already exists in "+dir);
			} else {
				log.report("Parsing "+filename);
				ResultsPackager.parseSOLformat(dir, filename, mapFile, freqFile, null, 1.0, 0.95, outfile, log);
			}
			qqFiles.add(outfile+",9="+ext.rootOf(outfile));
			lowCallrateMarkerFiles.add(dir+outfile+"lowCallRateMarkers.out");
			inputFiles = Array.addStrToArray(outfile, inputFiles);
		}
		Unique.proc(Array.toStringArray(lowCallrateMarkerFiles), null, null, dir+"lowCallrateMarkers.dat", null, true);
		if (Files.exists(dir+root+"_InvVar1.out")) {
			log.report(root+"_InvVar1.out already exists in "+dir);
		} else {
			log.report("Running inverse variance weighted meta-analysis...");
			Metal.metaAnalyze(dir, inputFiles, Aliases.MARKER_NAMES, root+"_InvVar", Metal.SE_ANALYSIS, null, true, log);
		}
		qqFiles.add(root+"_InvVar1.out"+",5="+ext.rootOf(root+"_InvVar"));
		if (Files.exists(dir+root+"_NWeighted1.out")) {
			log.report(root+"_NWeighted1.out already exists in "+dir);
		} else {
			log.report("Running sample size weighted meta-analysis...");
			Metal.metaAnalyze(dir, inputFiles, Aliases.MARKER_NAMES, root+"_NWeighted", Metal.PVAL_ANALYSIS, null, true, log);
		}
		qqFiles.add(root+"_NWeighted1.out"+",7="+ext.rootOf(root+"_NWeighted"));
		
		Files.write("java -cp C:/home/npankrat/vis.jar cnv.plots.QQPlot files=\""+Array.toStr(Array.toStringArray(qqFiles), ";")+"\" maxToPlot=10", dir+"plotQQs.bat");
		dir = dir.substring(dir.substring(0,dir.length()-1).lastIndexOf("/")+1, dir.length());
		log.report("QQ_FILENAMES="+dir+Array.toStr(Array.toStringArray(qqFiles), ";"+dir));
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String geno = "N:/statgen/CALiCo/filteredGenotypes/plink";
//		String genos = "N:/statgen/CALiCo_ARIC/filteredGenotypes/plink";
//		String phenoCovarFilename = "N:/statgen/CALiCo/BMI/Males.xln";
		String phenoCovarFilename = null;
		String phenoCovarDir = "N:/statgen/CALiCo/Age_at_menopause/xln/";
//		String phenoCovarDir = "N:/statgen/CALiCo_ARIC/LungFunction/xln/";
//		String scratchDir = phenoCovarDir.substring(0, phenoCovarDir.substring(0, phenoCovarDir.length()-2).lastIndexOf("/")) + "/scratches/";
		String scratchDir = "D:/scratch/";
//		String resultDir = phenoCovarDir.substring(0, phenoCovarDir.substring(0, phenoCovarDir.length()-2).lastIndexOf("/")) + "/results/";
		String resultDir = phenoCovarDir + (Files.exists(phenoCovarDir + CONDITIONALS_TXT_FILE)?"conditionals/":"results/");
		String plinkCommand;
		String modelListFileNameWhenParsingForPageGlucoseInsulinConditionalPaper = "D:/CALiCo_conditional/instructions/RequestedConditionalAnalyses.txt";
		String outputFileNameWhenParsingForPageGlucoseInsulinConditionalPaper = "D:/CALiCo_conditional/parsedResults.csv";
		Logger log;
		boolean exists = true;
		
//		metaAnalyzeSOL("D:/data/SOL/", "MODEL3wPCs_", "N:/statgen/CALICo_SOL/SOL-2013-04-05_Metabochip-mappingfile.txt");
//		metaAnalyzeSOL("D:/data/SOL/", "MODEL4wPCs_", "N:/statgen/CALICo_SOL/SOL-2013-04-05_Metabochip-mappingfile.txt");
//		System.exit(1);

//		metaAnalyzeSOL("D:/data/SOL/models12/", "MODEL1", "N:/statgen/CALICo_SOL/SOL-2013-04-05_Metabochip-mappingfile.txt");
//		metaAnalyzeSOL("D:/data/SOL/models12/", "MODEL2", "N:/statgen/CALICo_SOL/SOL-2013-04-05_Metabochip-mappingfile.txt");
		
//		metaAnalyzeSOL("D:/data/SOL/GlucoseInsulin/", "MODEL1", "N:/statgen/CALICo_SOL/SOL-2013-04-05_Metabochip-mappingfile.txt");
//		metaAnalyzeSOL("D:/data/SOL/GlucoseInsulin/", "MODEL2", "N:/statgen/CALICo_SOL/SOL-2013-04-05_Metabochip-mappingfile.txt");
//		metaAnalyzeSOL("D:/data/SOL/GlucoseInsulin/", "MODEL3", "N:/statgen/CALICo_SOL/SOL-2013-04-05_Metabochip-mappingfile.txt");
//		metaAnalyzeSOL("D:/data/SOL/GlucoseInsulin/", "MODEL4", "N:/statgen/CALICo_SOL/SOL-2013-04-05_Metabochip-mappingfile.txt");

//		metaAnalyzeSOL("N:/statgen/CALICo_SOL/T2DM/results_ver6/", "T2DM_MODEL1_SolGeno#_%%_ver6.out", "T2DM", "N:/statgen/CALICo_SOL/SOL-2013-04-05_Metabochip-mappingfile.txt");
//		System.exit(1);

		plinkCommand = "C:/plink/plink";

		String usage = "\n" +
		"one.CALiCo requires 0-1 arguments\n" +
		"  Note: to run conditional analysis, please save the file \"conditionals.txt\" in the directory of pheno and covariates.\n" +
		"   (1) plink command (i.e. plinkcommand=" + plinkCommand + " (default))\n" +
		"   (2) location of genotype files (i.e. geno=" + geno + " (default))\n" +
		"   (3) phenotype filename (i.e. phenocovarfile=" + phenoCovarFilename + " (default))\n" +
		" 	OR:\n" +
		"   (3) directory to perform for all *.xln files (i.e. phenocovardir=C:/test/ (not the default))\n" +
		"   (4) (optional) scratch directory (i.e. scratchdir=" + scratchDir + " (default; use null for directory of phenotype file))\n" +
		"   (5) (optional) results directory (i.e. resultdir=" + resultDir + " (default; use null for directory of phenotype file))\n" +
		"\n" +
		" To parse results per PAGE Glucose Insulin Paper's template:\n" +
		"   (1) full path to the model list file (i.e. pageglucoseinsulinmodellist=" + modelListFileNameWhenParsingForPageGlucoseInsulinConditionalPaper  + " (default))\n" +
		"   (2) results directory (i.e. resultdir=" + resultDir + " (default))\n" +
		"   (3) full path to the output file of this parsing (i.e. pageglucoseinsulinoutput=" + outputFileNameWhenParsingForPageGlucoseInsulinConditionalPaper + " (default))\n" +
		"";

		phenoCovarDir = null;
		phenoCovarFilename = null;
		scratchDir = null;
		resultDir = null;
		outputFileNameWhenParsingForPageGlucoseInsulinConditionalPaper = null;

		/**
		 * The following is for parsing the Results for Page Glucose Insulin Paper's Format.
		 */
//		parseResultsForPageGlucoseInsulinPaperFormat("D:/CALiCo_conditional/instructions/RequestedConditionalAnalyses.txt", "D:/CALiCo_conditional/results/", "D:/CALiCo_conditional/parsedResults.csv");
//		System.exit(0);
//		modelListFileNameWhenParsingForPageGlucoseInsulinConditionalPaper = "D:/CALiCo_conditional/instructions/RequestedConditionalAnalyses.txt";
//		resultDir = "D:/CALiCo_conditional/results/";
//		outputFileNameWhenParsingForPageGlucoseInsulinConditionalPaper = "D:/CALiCo_conditional/parsedResults.csv";

		/**
		 * The following 3 parameters are all what needed for conditional analysis, in addition to the "conditionals.txt" file mentioned in the help.
		 */
		geno = "D:/CALiCo_conditional/plink";
		phenoCovarDir = "D:/CALiCo_conditional/";
		plinkCommand = null;

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("plinkcommand=")) {
				plinkCommand = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("phenocovarfile=")) {
				phenoCovarFilename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("geno=")) {
				geno = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("scratchdir=")) {
				scratchDir = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("phenocovardir=")) {
				phenoCovarDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("resultdir=")) {
				resultDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("pageglucoseinsulinmodellist=")) {
				modelListFileNameWhenParsingForPageGlucoseInsulinConditionalPaper = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("pageglucoseinsulinoutput=")) {
				outputFileNameWhenParsingForPageGlucoseInsulinConditionalPaper = args[i].split("=")[1];
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
			if (resultDir == null) {
				resultDir = phenoCovarDir + "results/";
			}
			if(!new File(resultDir).exists()) {
				new File(resultDir).mkdir();
				exists = false;
			}
			
			if (scratchDir == null) {
				scratchDir = phenoCovarDir + "scratch/";
			}
			if(!new File(scratchDir).exists()) {
				new File(scratchDir).mkdir();
			}

			log = new Logger(resultDir + "Genvisis_CALiCo_" + (new SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date())) + ".log");
			log.report("Genvisis (R) 2013. \nCalico analysis "
						+ (new SimpleDateFormat("MM/dd/yyyy HH:mm:ss").format(new Date()))
						+ "\n\n-Geno data directory and prefix: " + geno
						+ "\n-Pheno and covariate data directory: " + phenoCovarFilename
						+ "\n-Scratch directory: " + scratchDir
						+ "\n-Result directory: " + resultDir);
			if (exists) {
				log.reportError("Warning --- Directory " + resultDir + " already exists. Existing files might be reused or overwritten.");
			} else {
				log.report("Creating result directory " + resultDir);
			}

			if(!new File(scratchDir).exists()) {
				new File(scratchDir).mkdir();
				log.report("Creating scratch directory " + scratchDir);
			} else {
				log.reportError("Warning --- Scratch directory " + scratchDir + " already exists. Existing files might be reused or overwritten.");
			}
			
			if (plinkCommand != null) {
//				parseAllFilesInDirectory(dir, genos, scratchDir);
//				if (Files.exists(phenoCovarDir + CONDITIONALS_TXT_FILE)) {
//					runConditional(genos, phenoCovarDir, resultDir, scratchDir, log);
//				} else {
					parseAllFilesInDirectory(geno, phenoCovarDir, resultDir, scratchDir, plinkCommand, log);
//				}

			} else if (phenoCovarFilename != null) {
				parsePhenotypes(geno, phenoCovarFilename, SAMPLE_ID_SYNONYMS[1], resultDir, log);

			} else {
			    System.out.println("Parse Results");
				parseResultsForPageGlucoseInsulinPaperFormat(modelListFileNameWhenParsingForPageGlucoseInsulinConditionalPaper, resultDir, outputFileNameWhenParsingForPageGlucoseInsulinConditionalPaper);
			}

			log.report("\nCalico analysis is finished.");

		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
