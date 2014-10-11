package link;

import java.io.*;
import java.util.*;

import parse.GenParser;
import stats.RegressionModel;
import common.*;
import filesys.FamilyStructure;

public class Heritability {
	public static final String DEFAULT_MERLIN_EXEC = "merlin";
	public static final String DEFAULT_SOLAR_EXEC = "solar";

	private static String computeWithMerlin(String dir, String pedfile, String pheno, String covars, String prefix, String merlinExec, Logger log) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		Hashtable<String, String> phenoHash, covarHash, seen, referenced;
		String[] keys, covarHeader;
		String temp;
		String estimate;
		
		estimate = ".";
		dir = ext.verifyDirFormat(dir);
		
		line = Files.getHeaderOfFile(dir+pheno, log);
		phenoHash = HashVec.loadFileToHashString(dir+pheno, new int[] {0,1}, new int[] {2}, pheno.endsWith(".csv"), "\t", true, false, false);
		covarHeader = null;
		if (covars != null) {
//			System.err.println("Error - covars not implemented for Merlin yet, will eventually construct a residual and run that instead; though it looks like it may work in the .dat file");
			covarHeader = Files.getHeaderOfFile(dir+covars, log);
			if (!covarHeader[0].equals("FID") || !covarHeader[1].equals("IID")) {
				System.err.println("Error - warning covariates file expected to start with FID and IID; may not merge properly without");
			}
			covarHash = HashVec.loadFileToHashString(dir+covars, new int[] {0,1}, Array.subArray(Array.intArray(covarHeader.length), 2), covars.endsWith(".csv"), "\t", true, false, false);
			covarHeader = Array.subArray(covarHeader, 2);
		} else {
			covarHash = new Hashtable<String, String>();
			covarHeader = new String[0];
		}
		
		seen = new Hashtable<String, String>();
		referenced = new Hashtable<String, String>();
		try {
			reader = new BufferedReader(new FileReader(pedfile));
			reader.mark(10000);
			line = reader.readLine().trim().split("[\\s]+");
			if (FamilyStructure.likelyPedHeader(line)) {
				if (log.getLevel() > 8) {
					log.report("Header row detected in the pedigree file");
				}
			} else {
				if (log.getLevel() > 8) {
					log.report("No header row detected in the pedigree file");
				}
				reader.reset();
			}
			writer = new PrintWriter(new FileWriter(dir+"re_chrom01.pre"));
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				writer.println(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]
						+"\t"+(phenoHash.containsKey(line[0]+"\t"+line[1])?phenoHash.get(line[0]+"\t"+line[1]):"x")
						+"\t"+(covarHash.containsKey(line[0]+"\t"+line[1])?covarHash.get(line[0]+"\t"+line[1]):Array.toStr(Array.stringArray(covarHeader.length, "x")))+
						"\t1\t1");
				seen.put(line[0]+"\t"+line[1], "");
				if (!line[2].equals("0")) {
					referenced.put(line[0]+"\t"+line[2], "1");
				}
				if (!line[3].equals("0")) {
					referenced.put(line[0]+"\t"+line[3], "2");
				}
			}
			reader.close();
			
			keys = HashVec.getKeys(referenced);
			for (int i = 0; i < keys.length; i++) {
				if (!seen.containsKey(keys[i])) {
					writer.println(keys[i]+"\t0\t0\t"+referenced.get(keys[i])+"\t"+(phenoHash.containsKey(keys[i])?phenoHash.get(keys[i]):"x")+"\t0\t0");
				}
			}
			
			writer.close();
			new LinkageMap(1, new String[] {"rs007"}, 2, new double[] {1}, false, false).createFile();
			Merlin.createMerlinFiles(1, dir+"chr01", Merlin.QUANTITATIVE_TRAIT, covarHeader);
			new File("map01.dat").delete();
			try {
				CmdLine.run(merlinExec+" -d chr01.dat -p re_chrom01.pre -m chr01.map --vc"+(covarHeader.length>0?" --useCovariates":"")+" --bits 1000 --megabytes 8000 --tabulate --step 5 --markerNames --information --prefix "+prefix, dir, new PrintStream(new File(dir+prefix+"_merlin.log")));
			} catch (IOException ioe) {
				log.reportError("Error - merlin failed to run/not available");
			}

			try {
				reader = new BufferedReader(new FileReader(dir+prefix+"_merlin.log"));
				while (reader.ready()) {
					temp = reader.readLine();
					if (temp.contains("h2 = ")) {
						estimate = temp.substring(temp.indexOf("h2 = ")+5, temp.length()-1);
						log.report("Merlin estimate: "+estimate);
					}
				}
				reader.close();
			} catch (FileNotFoundException fnfe) {
				log.reportError("Error - merlin failed to run");
				fnfe.printStackTrace();
			} catch (IOException ioe) {
				log.reportError("Error reading file \"" + dir+prefix+".log" + "\"; merlin probably ran into trouble");
			}
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + pedfile + "\" not found in current directory");
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + pedfile + "\"");
		}
		
		return estimate;
	}

	private static void computeWithSolar(String dir, String pedfile, String pheno, String covars, String prefix, String solarExec, Logger log) {
		BufferedReader reader;
		String temp;
		String[] ids, params, line, covarsHeader;
		String trait;

		dir = ext.verifyDirFormat(dir);
		
		GenParser.parse(new String[] {dir+pedfile, "out="+dir+"/"+prefix+"_fam.csv", "0=FAMID", "1=ID", "2=FA", "3=MO", "4=SEX"}, log);
		ids = HashVec.loadFileToStringArray(dir+pedfile, false, new int[] {1}, false);
		if (Array.unique(ids).length != ids.length) {
			System.err.println("Error - lookup for solar heritabilty screen currently requires unique IDs (found "+(ids.length-Array.unique(ids).length)+" that were not unique); aborting solar run");
			return;
		}
		
		params = new String[covars==null?2:3];
		params[0] = dir+pedfile+" 1 0=FAMID 1=ID skip=0";
		params[1] = dir+pheno+" 1 2";
		trait = Files.getHeaderOfFile(dir+pheno, log)[2];
		if (covars != null) {
			covarsHeader = Files.getHeaderOfFile(dir+covars, log);
			line = Array.toStringArray(Array.intArray(covarsHeader.length));
			line = Array.removeFromArray(line, 0);
			params[2] = dir+covars+" "+Array.toStr(line, " ");
		} else {
			covarsHeader = null;
		}
		Files.combine(ids, params, null, null, "", dir+"/"+prefix+"_ptypes.csv", log, true, true, true);
		Solar.cleanPhenotypeFile(dir+"/"+prefix+"_ptypes.csv");

//		Files.write("echo -e \"load ped "+prefix+"_fam.csv\\nautomodel "+prefix+"_ptypes.csv "+trait+"\\npolygenic -screen\\nquit\\n\" | "+solarExec+" > "+prefix+"_solar.log", dir+"/batch");
		Files.write("echo -e \"load ped "+prefix+"_fam.csv\\nload phenotype "+prefix+"_ptypes.csv\\ntrait "+trait+"\\n"+(covarsHeader != null?"covariates "+Array.toStr(Array.subArray(covarsHeader, 2), " ")+"\\n":"")+"polygenic -screen\\nquit\\n\" | "+solarExec+" > "+prefix+"_solar.log", dir+"/batch");
		if (!System.getProperty("os.name").startsWith("Windows")) {
			CmdLine.run("chmod +x batch", dir);
			CmdLine.run("./batch", dir);
		}
		
		try {
			reader = new BufferedReader(new FileReader(dir+trait+"/polygenic.out"));
			while (reader.ready()) {
				temp = reader.readLine();
				if (temp.contains("H2r is ")) {
					log.report("Solar  estimate: "+temp.substring(temp.indexOf("H2r is ")+7));
					
					temp = reader.readLine();
					while (!temp.contains("Loglikelihoods")) {
						log.report(temp);
						temp = reader.readLine();
					}
					
					while (reader.ready()) {
						temp = reader.readLine();
						if (temp.contains("Kurtosis")) {
							log.report(temp);
						}
					}
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error - solar failed to run");
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + trait+"/polygenic.out\"; solar probably ran into trouble");
		}
	}

	public static void getHeritabilitiesOfAllPhenosInADir(String dir, String pedigreeFile, String covars, Logger log) {
		String[] phenos, results;
		String merlinExec, solarExec;
		
		merlinExec = DEFAULT_MERLIN_EXEC;
		solarExec = DEFAULT_SOLAR_EXEC;
		phenos = Files.list(dir, ".txt", false);
		results = new String[phenos.length];
		for (int i = 0; i < phenos.length; i++) {
			results[i] = computeWithMerlin(dir, pedigreeFile, phenos[i], covars, ext.rootOf(phenos[i]), merlinExec, log);
		}
		Files.writeList(results, dir + "heritabilities.txt");
	}

	public static void fromParameters(String filename, Logger log) {
		PrintWriter writer, summary;
		String[] line;
		Hashtable<String,String> famIdHash;
		Vector<String> params;
		Vector<String[]> models;
		String pedigreeFile, dbFile;
		String merlinExec, solarExec;
		String temp;
		String[][] data;
		String root;
		String dbDelimiter;
		String[] dbHeader;
		boolean[] use;
		int numNotInPed;
		int[] indices;
		CountHash counter;
		String dir;
		String merlinEstimate;
		String[] solarEstimate;
		int numOfAllSamples;
		int numOfFamiliesSizedTwoOrAbove;
		int numOfFamiliesSizedOne;
		
		dir = ext.parseDirectoryOfFile(filename);
		
		params = Files.parseControlFile(filename, "heritability", new String[] {
				"ped=plink.fam",
				"db=input.xln",
				"# name of model (will be used to make a subdirectory) and dependent variable, followed by all covariates to be included",
				"model1 pta age sex PC1 PC2",
				"model2 icam1 age sex center",
				"",
				"# optional specifications of program locations; simply correct and un-comment",
				"# merlin_exec=/bin/merlin",
				"# solar_exec=/bin/solar"}, log);
		
		if (params != null) {
			pedigreeFile = null;
			dbFile = null;
			models = new Vector<String[]>();
			merlinExec = DEFAULT_MERLIN_EXEC;
			solarExec = DEFAULT_SOLAR_EXEC;
    		for (int i = 0; i < params.size(); i++) {
    			temp = params.elementAt(i);

    			if (temp.startsWith("ped=")) {
    				pedigreeFile = ext.parseStringArg(temp, "failed_to_define_pedigree_file");
    			} else if (temp.startsWith("db=")) {
    				dbFile = ext.parseStringArg(temp, "failed_to_define_db_file");
    			} else if (temp.startsWith("merlin_exec=")) {
    				merlinExec = ext.parseStringArg(temp, DEFAULT_MERLIN_EXEC);
    			} else if (temp.startsWith("solar_exec=")) {
    				solarExec = ext.parseStringArg(temp, DEFAULT_SOLAR_EXEC);
    			} else {
    				models.add(temp.trim().split("[\\s]+"));
    			}
			}
    		
    		if (!Files.exists(pedigreeFile)) {
    			log.report("Error - could not find pedigree file '"+pedigreeFile+"'; aborting");
    			return;
    		}
    		
    		if (!Files.exists(dbFile)) {
    			log.report("Error - could not find database file '"+dbFile+"'; aborting");
    			return;
    		}
    		
    		famIdHash = HashVec.loadFileToHashString(pedigreeFile, new int[] {1}, new int[] {0}, Files.determineDelimiter(pedigreeFile, log).equals(","), null, false, false, false);
    		
    		log.setLevel(8);

    		try {
				summary = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+"_summary.xln"));
				summary.println("Model\tMerlin_est.\tSolar_est.\tSolar_p\tn_Samples\tn_Families\tn_Families_size>1\tAverage_size_families_siez>1\tn_Families_size=1");
	    		for (int i = 0; i < models.size(); i++) {
	    			line = models.elementAt(i);
	    			if (line.length < 2) {
	    				log.reportError("Error - need at least two arguments (model->directory and the model, even if there are no more covariates), found:");
	    				log.reportError(Array.toStr(line));
	    			} else {
	    				root = line[0];
	    				root = ext.replaceWithLinuxSafeCharacters(root, false);
	    				new File(root).mkdirs();

	    				dbDelimiter = Files.determineDelimiter(dbFile, log);
	    				dbHeader = Files.getHeaderOfFile(dbFile, dbDelimiter, log);
	    				line[0] = dbHeader[0];
	    				indices = ext.indexFactors(line, dbHeader, false, log, true, false);
	    				if (Array.min(indices) == -1) {
	    					summary.close();
	    					return;
	    				}
	    				data = HashVec.loadFileToStringMatrix(dbFile, true, indices, dbDelimiter, false, 1000, false);
	    				use = RegressionModel.getRowsWithCompleteData(null, Matrix.prune(data, null, Array.subArray(Array.intArray(line.length), 1), log), log);
	    				
	    				numNotInPed = 0;
						counter = new CountHash();
	    				try {
							writer = new PrintWriter(new FileWriter(root+"/pheno.dat"));
							writer.println("FID\tIID\t"+line[1]);
							for (int j = 0; j < use.length; j++) {
								if (use[j]) {
									if (famIdHash.containsKey(data[j][0])) {
										writer.println(famIdHash.get(data[j][0])+"\t"+data[j][0]+"\t"+data[j][1]);
										counter.add(famIdHash.get(data[j][0]));
									} else {
										if (numNotInPed == 0) {
											log.reportError("Warning - the following samples were not found in the pedigree file:");
										}
										if (numNotInPed == 10) {
											log.reportError("...");
										} else if (numNotInPed < 10) {
											log.reportError(data[j][0]);
										}
										numNotInPed++;
										use[j] = false;
									}
								}
							}
							writer.close();
						} catch (Exception e) {
							System.err.println("Error writing to " + root+"/pheno.dat");
							e.printStackTrace();
						}
//	    				counter.getCounts();

	    				if (line.length > 2) {
		    				try {
								writer = new PrintWriter(new FileWriter(root+"/covars.dat"));
								writer.println("FID\tIID\t"+Array.toStr(Array.subArray(line, 2)));
								for (int j = 0; j < use.length; j++) {
									if (use[j]) {
										writer.println(famIdHash.get(data[j][0])+"\t"+data[j][0]+"\t"+Array.toStr(Array.subArray(data[j], 2)));
									}
								}
								writer.close();
							} catch (Exception e) {
								System.err.println("Error writing to " + root+"/pheno.dat");
								e.printStackTrace();
							}
		    				
							if (numNotInPed > 0) {
								log.reportError("There were "+numNotInPed+" sample(s) with valid and complete phenotype data that were not in the pedigree file");
							}
	    				}
						
						log.report("Heritability for "+root);
						log.report(line[1]+" ~ "+Array.toStr(Array.subArray(line, 2), " "));
						merlinEstimate = computeWithMerlin(dir+root, pedigreeFile, "pheno.dat", line.length > 2?"covars.dat":null, root, merlinExec, log);
//						solarEstimate = computeWithSolar(dir+root, pedigreeFile, "pheno.dat", line.length > 2?"covars.dat":null, root, solarExec, log);
						numOfAllSamples = counter.getTotalCount();
						numOfFamiliesSizedOne = counter.getSizeOfCountEquals(1);
						numOfFamiliesSizedTwoOrAbove = counter.getSizeOfCountGreaterThan(2);
						log.report("Number of samples: " + numOfAllSamples + "\nNumber of families: " + counter.getSize() + "\nNumber of families of size>=2: " + numOfFamiliesSizedTwoOrAbove + "\nAverage size of families of size>=2: " + ext.formDeci((numOfAllSamples - numOfFamiliesSizedOne) / (float) numOfFamiliesSizedTwoOrAbove, 3) + "\nNumber of families of size=1: " + numOfFamiliesSizedOne);
//						summary.println(root + "\t" + merlinEstimate + "\t" + solarEstimate[0] + "\t" + solarEstimate[1] + "\t" + numOfAllSamples + "\t" + counter.getSize() + "\t" + numOfFamiliesSizedTwoOrAbove + "\t" + String.format("%.3", ((float) (numOfAllSamples - numOfFamiliesSizedOne)) / numOfFamiliesSizedTwoOrAbove) + "\t" + numOfFamiliesSizedOne);
						summary.println(root + "\t" + merlinEstimate + "\t\t\t" + numOfAllSamples + "\t" + counter.getSize() + "\t" + numOfFamiliesSizedTwoOrAbove + "\t" + ext.formDeci((float) (numOfAllSamples - numOfFamiliesSizedOne) / numOfFamiliesSizedTwoOrAbove, 3) + "\t" + numOfFamiliesSizedOne);
						log.report("");
	    			}
				}				
				summary.close();
			} catch (Exception e) {
				System.err.println("Error writing to " + ext.rootOf(filename, false)+"_summary.xln");
				e.printStackTrace();
			}
    		
    		
		}
	}
	
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String phenoDir = "N:/statgen/Mitochondrial_CN_heritability/phenos/";
		String pheno = "pheno.dat";
		String covars = "covars.dat";
		String pedfile = "plink.fam";
//		String method = "merlin";
//		String method = "solar";
		String method = "both";
		Logger log;
		String logfile = null;
		String controlFile = null;
		String prefix = "vc-chr01";

//		fromParameters("N:/statgen/BOSS/phenotypes/PhenoPrep/word_recogition/heritability_noCovars.crf", new Logger());
//		System.exit(0);

		String usage = "\n" + 
		"link.Heritability requires 0-1 arguments\n" + 
		"   (1) pedigree filename (i.e. ped=" + pedfile + " (default))\n" + 
		"   (2) phenotype filename (i.e. pheno=" + pheno + " (default))\n" + 
		"   (3) covariates filename (i.e. covars=" + covars + " (default; set to null if not needed))\n" + 
		"   (4) method: merlin/solar/both (i.e. method=" + method + " (default))\n" + 
		"   (5) (optional) prefix for the files (i.e. prefix=" + prefix + " (default))\n" + 
		" OR:\n" + 
		"   (1) control file (i.e. crf=pheno.crf (not the default))\n" + 
		"";

		controlFile = null;
		method = "";
		phenoDir = null;

		phenoDir = "N:/statgen/Mitochondrial_CN_heritability/phenos/";
		pedfile = "N:/statgen/Mitochondrial_CN_heritability/pedigrees/pedigree.dat";
		covars = null;

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("ped=")) {
				pedfile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("pheno=")) {
				pheno = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("covars=")) {
				covars = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("method=")) {
				method = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("prefix=")) {
				method = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("crf=")) {
				controlFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("phenodir=")) {
				phenoDir = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		
//		if (args.length == 0) {
//			controlFile = "heritabilityTest.crf";
//			controlFile = "N:/statgen/BOSS/phenotypes/PhenoPrep/MMSE/heritability.crf";
//		}
		try {
			log = new Logger(logfile);
			if (controlFile != null) {
				fromParameters(controlFile, new Logger(ext.rootOf(controlFile, false)+".log"));
			} else {
				if (phenoDir != null) {
					getHeritabilitiesOfAllPhenosInADir(phenoDir, pedfile, covars, log);
				} else {
					if (method.equalsIgnoreCase("merlin") || method.equalsIgnoreCase("both")) {
						computeWithMerlin(pedfile, pheno, covars, prefix, "", DEFAULT_MERLIN_EXEC, log);
					}
					if (method.equalsIgnoreCase("solar") || method.equalsIgnoreCase("both")) {
						computeWithSolar(pedfile, pheno, covars, prefix, "", DEFAULT_SOLAR_EXEC, log);
					}
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
