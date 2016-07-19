package org.genvisis.gwas;

import java.io.*;
import java.util.*;

import org.genvisis.common.*;

public class GWAF {
	public static final String[] HEADER_IMPUTED = {"phen", "snp", "N", "AF", "h2q", "beta", "se", "pval"};
	public static final String[] HEADER_GENOTYPED = {"phen", "snp", "n0", "n1", "n2", "h2q", "beta", "se", "chisq", "df", "model", "pval"};
	public static final String[] HEADER_SUMMARY = {"MarkerName", "beta", "StdErr", "Pvalue", "N", "freqA1", "MAF"};
//	public static final String DELIMITER = "\t";
	public static final String DELIMITER = ",";

	public static final int DEFAULT_VERSION_USED = 1;

	public static void parse(String outfileTemplate, int startAt, boolean imputedNotGenotype, String outfile) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, header, expectedHeader;
		int count, trav, last;
		Logger log;
		boolean first;
		IntVector iv, headerIndices;
		double freq, maf;
		int numPer, lineOffset, n;
		
		if (outfile == null) {
			outfile = ext.rootOf(ext.replaceAllWithSafer(outfileTemplate, "#", ""))+(DELIMITER.equals(",")?".csv":".xln");
		}
		
		expectedHeader = imputedNotGenotype?HEADER_IMPUTED:HEADER_GENOTYPED;
		
		numPer = -1;
		first = true;
		log = new Logger();
		iv = new IntVector();
		try {
			writer = new PrintWriter(new FileWriter(outfile));
			trav = last = startAt;
			while (trav < last+20) {
				if (new File(ext.replaceAllWith(outfileTemplate, "#", trav+"")).exists()) {
					if (iv.size() > 0) {
						log.reportError("Skipping over file"+(iv.size()>1?"s":"")+" '"+ext.replaceAllWith(outfileTemplate, "#", iv.elementAt(0)+"")+"'"+(iv.size()>1?" through '"+ext.replaceAllWith(outfileTemplate, "#", iv.elementAt(iv.size()-1)+"")+"'":""));
						iv.clear();
					}
					try {
						reader = new BufferedReader(new FileReader(ext.replaceAllWith(outfileTemplate, "#", trav+"")));
						headerIndices = new IntVector();
						count = 0;
						while (reader.ready()) {
							line = reader.readLine().trim().split(",");
							if (Array.equals(line, expectedHeader, false)) {
								headerIndices.add(count);
							}
							count++;
						}
						reader.close();
						lineOffset = 0;
						if (headerIndices.size() == 0) {
							log.reportError("Error - file '"+ext.replaceAllWith(outfileTemplate, "#", trav+"")+"' does not have a proper header");
						} else if (headerIndices.size() > 1) {
							log.reportError("Error - file '"+ext.replaceAllWith(outfileTemplate, "#", trav+"")+"' has muliple headers (a total of "+headerIndices.size()+"); using only the final set of data");
							lineOffset = headerIndices.elementAt(headerIndices.size()-1);
						} else if (headerIndices.elementAt(0) != 0) {
							log.reportError("Error - file '"+ext.replaceAllWith(outfileTemplate, "#", trav+"")+"' has data before the expected header; using only the data after the final header");
							lineOffset = headerIndices.elementAt(0);
						}
						
						reader = new BufferedReader(new FileReader(ext.replaceAllWith(outfileTemplate, "#", trav+"")));
						for (int i = 0; i < lineOffset; i++) {
							reader.readLine();
						}
						header = reader.readLine().trim().split(",");
						
						ext.checkHeader(header, expectedHeader, true);
						if (first) {
							writer.println(Array.toStr(HEADER_SUMMARY, DELIMITER));
							first = false;
						}
						count = 0;
						while (reader.ready()) {
							count++;
							line = reader.readLine().trim().split(",", -1);
							if (line.length < expectedHeader.length) {
								System.err.println("Error - truncated file at marker number "+count+(line.length>1?" ("+line[1]+")":"")+" in file '"+ext.replaceAllWith(outfileTemplate, "#", trav+"")+"'");
							}
							if (imputedNotGenotype) {
								n = Integer.parseInt(line[2]);
								freq = Double.parseDouble(line[3]);
							} else {
								n = Integer.parseInt(line[2])+Integer.parseInt(line[3])+Integer.parseInt(line[4]);
								freq = ( Double.parseDouble(line[4])*2+Double.parseDouble(line[3])*1 )/(double)(n*2);
								
							}
							maf = freq>0.50?1-freq:freq;

							if (imputedNotGenotype) {
								writer.println(line[1]+DELIMITER+line[5]+DELIMITER+line[6]+DELIMITER+line[7]+DELIMITER+n+DELIMITER+freq+DELIMITER+maf);
							} else {
								writer.println(line[1]+DELIMITER+line[6]+DELIMITER+line[7]+DELIMITER+line[11]+DELIMITER+n+DELIMITER+freq+DELIMITER+maf);
							}
						}
						if (numPer == -1) {
							numPer = count;
						} else if (count != numPer) {
							log.report("File '"+ext.replaceAllWith(outfileTemplate, "#", trav+"")+"' had "+count+" markers, whereas previous files had "+numPer);
						}						
						reader.close();
					} catch (FileNotFoundException fnfe) {
						System.err.println("Error: file \"" + ext.replaceAllWith(outfileTemplate, "#", trav+"")+ "\" not found in current directory");
						System.exit(1);
					} catch (IOException ioe) {
						System.err.println("Error reading file \"" + ext.replaceAllWith(outfileTemplate, "#", trav+"")+ "\"");
						System.exit(2);
					}
					
					last = trav;
				} else {
					iv.add(trav);
				}
				if (trav < 0) {
					trav = Integer.MAX_VALUE;
				} else {
					trav++;
				}
			}
			if (first) {
				log.reportError("Error - no files were found matching the template '"+outfileTemplate+"'");
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + ext.rootOf(outfileTemplate)+".xln");
			e.printStackTrace();
		}
	}

	public static final String[] IBC_OUTPUT_FORMAT1 = {"CHR", "POS", "SNP", "STRAND (Illumina)", "STRAND (HapMap)", "N", "EFFECT_ALLELE1", "NON_EFFECT_ALLELE", "EA_FREQ", "BETA", "SE", "P_VAL"};
	public static final String[] TRADITIONAL_OUTPUT_FORMAT = {"Chr", "Position", "MarkerName", "Strand", "HapMapStrand", "N", "Effect_allele", "Reference_allele", "Freq1", "BETA", "SE", "P-value"};

	public static void parseToMetal(String resultsFile, String infoFile, String markersToReport, String outfile, boolean includeN, boolean includeFreq, boolean includeMAF) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String temp, trav;
		Hashtable<String, String> alleleHash;
		HashSet<String> markerHash;
		String delimiter;
		
		String dir = "";
		
		
		if (outfile == null) {
			outfile = ext.rootOf(resultsFile, false)+"_out.csv";
		}
		
		if (markersToReport != null) {
			markerHash = HashVec.loadFileToHashSet(dir+markersToReport, false);
		} else {
			markerHash = null;
		}
		
		alleleHash = HashVec.loadFileToHashString(dir+infoFile, new int[] {0}, new int[] {1,2}, false, "\t", true, false, false);
		
		try {
			reader = Files.getAppropriateReader(dir+resultsFile);
			writer = Files.getAppropriateWriter(dir+outfile);
			temp = reader.readLine().trim();
			delimiter = ext.determineDelimiter(temp);
			line = temp.split(delimiter);
			ext.checkHeader(line, GWAF.HEADER_SUMMARY, true);
//			writer.println(Array.toStr(IBC_OUTPUT_FORMAT));
			writer.print("MarkerName\tAllele1\tAllele2\tEffect\tStdErr\tP-value");
			if (includeN) {
				writer.print("\tN");
			}
			if (includeFreq) {
				writer.print("\tFreqA1");
			}
			if (includeMAF) {
				writer.print("\tMAF");
			}
			writer.println();
			while (reader.ready()) {
				line = reader.readLine().trim().split(delimiter);
				trav = line[0];
				if ((markerHash == null || markerHash.contains(trav)) && !line[3].equals("")) {
					writer.print(trav);
					if (alleleHash.containsKey(trav)) {
						writer.print("\t"+alleleHash.get(trav));
					} else {
						System.err.println("Error - no alleles from info file for "+trav);
						writer.print("\t.\t.");
					}
//					public static final String[] HEADER_SUMMARY = {"MarkerName", "beta", "StdErr", "Pvalue", "N", "freqA1", "MAF"};
					writer.print("\t"+line[1]+"\t"+line[2]+"\t"+line[3]);
					if (includeN) {
						writer.print("\t"+line[4]);
					}
					if (includeFreq) {
						writer.print("\t"+line[5]);
					}
					if (includeMAF) {
						writer.print("\t"+line[6]);
					}					
					writer.println();
				}
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + resultsFile + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + resultsFile + "\"");
			System.exit(2);
		}
	}
	
	
	public static void parseToMetal(String resultsFile, String mapFile, String originalFrqFile, String markersToReport, String outfile) {
//		BufferedReader reader;
//		PrintWriter writer;
//		String[] line;
//		String temp, trav;
//		Hashtable<String, String> markerHash, mapHash, originalFreqHash; // , customFreqHash;
//		String delimiter;
//		
//		String dir = "";
//		
//		
//		if (outfile == null) {
//			outfile = ext.rootOf(resultsFile, false)+"_out.csv";
//		}
//		
//		if (markersToReport != null) {
//			markerHash = HashVec.loadFileToHashNull(dir+markersToReport, false);
//		} else {
//			markerHash = null;
//		}
//		
//		mapHash = HashVec.loadFileToHashString(dir+mapFile, new int[] {1}, new int[] {0,3}, false, "\t", false, false, false);
//		originalFreqHash = HashVec.loadFileToHashString(dir+originalFrqFile, new int[] {1}, new int[] {2,3}, false, "\t", false, false, false); // add 4 if you want global frequency
////		originalFreqHash = HashVec.loadFileToHashString(originalFrqFile, new int[] {1}, customFrqFile==null?new int[] {2,3,4}:new int[] {2,3}, false, "\t", false, false, false); // add 4 if you want global frequency instead of custom Freq
//
//		
////		if (customFrqFile != null) {
////			System.err.println("Error - use of custom freq file is not currently implemented; still only computes from file");
////			customFreqHash = HashVec.loadFileToHashString(originalFrqFile, new int[] {1}, new int[] {2,3,4}, false, "\t", false, false, false); // add 4 if you want global frequency instead of custom Freq
////		} else {
////			customFreqHash = null;
////		}
//		
//		try {
//			reader = Files.getAppropriateReader(dir+resultsFile);
//			writer = Files.getAppropriateWriter(dir+outfile);
//			temp = reader.readLine().trim();
//			delimiter = ext.determineDelimiter(temp);
//			line = temp.split(delimiter);
//			ext.checkHeader(line, GWAF.HEADER_SUMMARY, true);
////			writer.println(Array.toStr(IBC_OUTPUT_FORMAT));
//			writer.println(Array.toStr(TRADITIONAL_OUTPUT_FORMAT));
//			while (reader.ready()) {
//				line = reader.readLine().trim().split(delimiter);
//				trav = line[0];
//				if ((markerHash == null || markerHash.containsKey(trav)) && !line[3].equals("")) {
//					if (mapHash.containsKey(trav)) {
//						writer.print(mapHash.get(trav));
//					} else if (mapHash.containsKey(ext.replaceAllWith(trav, ".", "-"))) {
//						writer.print(mapHash.get(ext.replaceAllWith(trav, ".", "-")));
//					} else {
//						System.err.println("Error - no map position for "+trav);
//						writer.print(".\t.");
//					}
//					writer.print("\t"+trav+"\tNA\t+\t"+line[4]);
//					if (originalFreqHash.containsKey(trav)) {
//						writer.print("\t"+originalFreqHash.get(trav));
//					} else if (originalFreqHash.containsKey(ext.replaceAllWith(trav, ".", "-"))) {
//						writer.print("\t"+originalFreqHash.get(ext.replaceAllWith(trav, ".", "-")));
//					} else {
//						System.err.println("Error - no alleles from original .frq file for "+trav);
//						writer.print("\t.\t.");
//					}
//					writer.println("\t"+line[5]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]);
//				}
//			}
//			reader.close();
//			writer.close();
//		} catch (FileNotFoundException fnfe) {
//			System.err.println("Error: file \"" + resultsFile + "\" not found in current directory");
//			System.exit(1);
//		} catch (IOException ioe) {
//			System.err.println("Error reading file \"" + resultsFile + "\"");
//			System.exit(2);
//		}
	}
	
	
	// use numBatches=1 if you want to run in another directory on Windows 
	public static void batch(String dir, String phenoFile, String pheno, String[] covars, String model, String geneticDataTemplate, int startAt, boolean imputedNotGenotype, String pedfile, String outfileTemplate, String rootTemplate, String[] nodesToUse, int numBatches, int versionOfGWAF) {
		PrintWriter writer;
		String[] list;
		int count; //, step;
		Vector<String> v;
		
		if (rootTemplate == null) {
			rootTemplate = "batches/"+pheno+"_file#";
		} else if (rootTemplate.endsWith(".qsub")) {
			rootTemplate = rootTemplate.substring(0, rootTemplate.lastIndexOf("."));
		}
		
		if (covars.length == 1 && covars[0].equalsIgnoreCase("null")) {
			covars = null;
		}
		
		new File("batches/").mkdirs();
		count = startAt;		
		try {
			if (!Files.exists(dir+"kmat.Rfile", false)) {
				writer = new PrintWriter(new FileWriter(dir+"createKmat.R"));
				if (versionOfGWAF == 2) {
					writer.println("library(kinship2)");
				} else {
					writer.println("library(kinship)");
				}
				writer.println("library(GWAF)");
				writer.println("cdata <- read.csv(\""+pedfile+"\", header=T)");
				writer.println("kmat <- makekinship(cdata$famid, cdata$id, cdata$fa, cdata$mo)");
				writer.println("kmat<-kmat*2");
				writer.println("save(kmat, file=\"kmat.Rfile\")");
				writer.close();
			}

			while (new File(ext.replaceAllWith(dir+geneticDataTemplate, "#", count+"")).exists() && (geneticDataTemplate.contains("#") || count == startAt)) {
				writer = new PrintWriter(new FileWriter(dir+ext.insertNumbers(rootTemplate, count)+".R"));
				if (versionOfGWAF == 2) {
					writer.println("library(kinship2)");
				} else {
					writer.println("library(kinship)");
				}
				writer.println("library(GWAF)");
				
				if (imputedNotGenotype) {
					writer.println("lme"+(versionOfGWAF==2?"pack":"")+".batch.imputed(\""+phenoFile+"\", \""+ext.replaceAllWith(geneticDataTemplate, "#", count+"")+"\", \""+pedfile+"\", \""+pheno+"\", \"kmat.Rfile\", covars="+(covars==null?"NULL":"c(\""+Array.toStr(covars, "\",\"")+"\")")+", \""+ext.replaceAllWith(outfileTemplate, "#", count+"")+"\", col.names=T, sep.ped=\",\", sep.phe=\",\", sep.gen=\",\")");
				} else {
					writer.println("lme"+(versionOfGWAF==2?"pack":"")+".batch(\""+phenoFile+"\", \""+ext.replaceAllWith(geneticDataTemplate, "#", count+"")+"\", \""+pedfile+"\", \""+pheno+"\", \"kmat.Rfile\", model=\""+model+"\", covars="+(covars==null?"NULL":"c(\""+Array.toStr(covars, "\",\"")+"\")")+", \""+ext.replaceAllWith(outfileTemplate, "#", count+"")+"\", col.names=T, sep.ped=\",\", sep.phe=\",\", sep.gen=\",\")");
				}

				writer.close();
				count++;
			}
			count--;
			System.out.println("last file seen was '"+ext.replaceAllWith(dir+geneticDataTemplate, "#", count+""));
			if (numBatches < 1) {
				if (nodesToUse == null) {
					v = Array.toStringVector(Files.qsub(dir, rootTemplate, startAt, count, "R --no-save < "+rootTemplate+".R > "+rootTemplate+".log", null, 5000, 12, null));
				} else {
					v = new Vector<String>();
//					step = (int)Math.ceil((double)((count-startAt)+1)/(double)nodesToUse.length);
//					for (int i = 0; i < nodesToUse.length; i++) {
//						list = Files.qsub("", null, i*step+startAt, i==nodesToUse.length-1?count:((i+1)*step+startAt-1), "R --no-save < batches/"+pheno+"_gwaf#.R > batches/"+pheno+"_file#.log", "batches/"+pheno+"_file", null, -1, nodesToUse[i]);
					for (int i = startAt; i <= count; i++) {
						list = Files.qsub(dir, rootTemplate, i, i, "R --no-save < "+rootTemplate+".R > "+rootTemplate+".log", null, 5000, 12, nodesToUse[i%nodesToUse.length]);
						for (int j = 0; j < list.length; j++) {
							v.add(list[j]);
						}
					}
				}
				if (!Files.exists(dir+"kmat.Rfile", false)) {
					v.insertElementAt("R --no-save < createKmat.R > createKmat.log", 0);
				}
				Files.writeList(Array.toStringArray(v), dir+"master."+pheno);
				Files.chmod(dir+"master."+pheno);
			} else {
				if (Files.isWindows()) {
					Files.write("R --no-save < createKmat.R > createKmat.log", dir+"createKmat.bat");
					Files.batchIt((numBatches==1?dir:"")+"run."+pheno, -1, startAt, count, numBatches, "R --no-save 0< "+rootTemplate+".R 1> "+rootTemplate+".log 2> "+rootTemplate+".err.log");
				} else {
					Files.batchIt((numBatches==1?dir:"")+"run."+pheno, -1, startAt, count, numBatches, "R --no-save < "+rootTemplate+".R > "+rootTemplate+".log");
				}
			}
		} catch (Exception e) {
			System.err.println("Error queuing up file " + count);
			e.printStackTrace();
		}
	}
	
	public static void splitFiles(String split, int blockSize, String geneticDataTemplate) {
		BufferedReader reader;
		PrintWriter writer;
		PrintWriter[] writers;
		String[] line;
		int numOutFiles;
		String delimiter;
		String dir;
		
		delimiter = Files.determineDelimiter(split, new Logger());
		dir = ext.parseDirectoryOfFile(geneticDataTemplate);
		new File(dir).mkdirs();
		
		try {
			reader = Files.getAppropriateReader(split);
			reader.mark(1000000);
			line = reader.readLine().trim().split(delimiter);
			numOutFiles = (int)Math.ceil((double)line.length/(double)blockSize);
			writers = new PrintWriter[numOutFiles];
			for (int i = 0; i < numOutFiles; i++) {
				writers[i] = Files.getAppropriateWriter(ext.replaceAllWith(geneticDataTemplate, "#", i+""));
			}
			reader.reset();
			while (reader.ready()) {
				line = reader.readLine().trim().split(delimiter, -1);
				writer = null;
				for (int i = 1; i < line.length; i++) {
					if (i % blockSize == 1) {
						writer = writers[i/blockSize];
						writer.print(line[0]);
					}
					writer.print(delimiter+line[i]);
				}
				for (int i = 0; i < numOutFiles; i++) {
					writers[i].println();
				}
			}
			reader.close();
			for (int i = 0; i < numOutFiles; i++) {
				writers[i].close();
			}
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + split + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + split + "\"");
			System.exit(2);
		}		
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String phenoFile = "pheno.csv";
		String pheno = "pheno";
		String[] covars = null;
		String geneticDataTemplate = "gwaf/file#.fhsR.gz";
		String pedfile = "pedfile.csv";
		String outfileTemplate = "results#.csv";
		String qsubRoot = null;
		int startAt = 0;
		String[] nodesToUse = null;
		int numBatches = -1;
		boolean parse = false;
		boolean imputedNotGenotype = true;
		String outfile = null;
		String split = null;
		int blockSize=5000;
		String model = "a";
		int versionOfGWAF = DEFAULT_VERSION_USED;

		String usage = "\n" + 
		"gwas.GWAF requires 0-1 arguments\n" + 
		"   (1) name of file with phenotypes (i.e. phenoFile=" + phenoFile + " (default; currently needs to be comma-delimited))\n" + 
		"   (2) name of phenotype to run (i.e. pheno=" + pheno + " (default))\n" + 
		"   (3) name of covariates to include (i.e. covars=" + covars + " (default; null leads to none; comma-delimited))\n" + 
		"   (4) genetic model to use: lowercase a/d/r/g (i.e. model=" + model + " (default))\n" + 
		"   (5) format of genotype filenames (i.e. genoPrimer=" + geneticDataTemplate + " (default; can be zipped))\n" + 
		"   (6) number to start looking for file pattern (i.e. startAt=" + startAt + " (default))\n" +
		"   (7) data is imputed not genotyped (i.e. imputed=" + imputedNotGenotype + " (default))\n" +
		"   (8) name of pedigree file (i.e. pedfile=" + pedfile + " (default))\n" + 
		"   (9) template for output files (i.e. outfileTemplate=" + outfileTemplate + " (default))\n" +
		"   (10) template for qsub files (i.e. qsubRoot=[trait]_file#.qsub (default))\n" +
		"   (11) nodes to use (i.e. nodesToUse=" + nodesToUse + " (default; qsubs only; full names, comma-delimited))\n" + 
		"   (12) number of batches to create (i.e. numBatches=" + numBatches + " (default; anything less than 1 leads to qsubs being made))\n" +
		"   (13) version of gwafToUse: 1 or 2 (i.e. gwafVersion=" + versionOfGWAF + " (default))\n" +
		" OR\n" +
		"   (1) parse results (i.e. -parseResults (not the default))\n" +
		"   (2) template for output files (i.e. outfileTemplate=" + outfileTemplate + " (default))\n" + 
		"   (3) number to start looking for file pattern (i.e. startAt=" + startAt + " (default))\n" + 
		"   (4) data is imputed not genotyped (i.e. imputed=" + imputedNotGenotype + " (default))\n" +
		"   (5) name of parsed output file (i.e. out=" + outfile + " (default))\n" +
		" OR\n" +
		"   (1) split file into subfiles (i.e. split=gwaf.csv (not the default))\n" +
		"   (2) size of chunks (i.e. size="+blockSize+" (default))\n" + 
		"   (3) template for output files (i.e. genoPrimer="+geneticDataTemplate+" (default))\n" + 
		"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("phenoFile=")) {
				phenoFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("pheno")) {
				pheno = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("covars=")) {
				covars = args[i].split("=")[1].split(",");
				numArgs--;
			} else if (args[i].startsWith("model=")) {
				model = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("genoPrimer=")) {
				geneticDataTemplate = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("startAt=")) {
				startAt = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("imputed=")) {
				imputedNotGenotype = ext.parseBooleanArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("pedfile=")) {
				pedfile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("outfileTemplate=")) {
				outfileTemplate = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("qsubRoot=")) {
				qsubRoot = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("nodesToUse=")) {
				nodesToUse = args[i].split("=")[1].split(",");
				numArgs--;
			} else if (args[i].startsWith("numBatches=")) {
				numBatches = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("-parseResults")) {
				parse = true;
				numArgs--;
			} else if (args[i].startsWith("out=")) {
				outfile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("split=")) {
				split = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("size=")) {
				blockSize = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("gwafVersion=")) {
				versionOfGWAF = ext.parseIntArg(args[i]);
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
//			outfileTemplate = "FHS_whites/inverseNormalizedPhenoWithConditionals_withCondi_results#.csv";
//			parse = true;
			
//			outfileTemplate = "pta_results#.csv";
//			parse = true;
//			imputedNotGenotype = false;

			if (parse) {
				parse(outfileTemplate, startAt, imputedNotGenotype, outfile);
			} else if (split != null) {
				splitFiles(split, blockSize, geneticDataTemplate);
			} else {
				batch("", phenoFile, pheno, covars, model, geneticDataTemplate, startAt, imputedNotGenotype, pedfile, outfileTemplate, qsubRoot, nodesToUse, numBatches, versionOfGWAF);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
