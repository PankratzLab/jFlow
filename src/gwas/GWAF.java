package gwas;

import java.io.*;
import java.util.*;

import common.*;

public class GWAF {
	public static final String[] HEADER_IMPUTED = {"phen", "snp", "N", "AF", "h2q", "beta", "se", "pval"};
	public static final String[] HEADER_GENOTYPED = {"cut and paste actual header"};
//	public static final String DELIMITER = "\t";
	public static final String DELIMITER = ",";
	
	public static void parse(String outfileTemplate, int startAt, boolean imputedNotGenotype, String outfile) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, header, expectedHeader;
		int count, trav, last;
		Logger log;
		boolean first;
		IntVector iv, headerIndices;
		double freq, maf;
		int numPer, lineOffset;
		
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
							writer.println("MarkerName"+DELIMITER+"beta"+DELIMITER+"StdErr"+DELIMITER+"Pvalue"+DELIMITER+"N"+DELIMITER+"freqA1"+DELIMITER+"MAF");
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
								freq = Double.parseDouble(line[3]);
							} else {
								freq = Double.parseDouble(line[3]); // update
							}
							maf = freq>0.50?1-freq:freq;
							if (imputedNotGenotype) {
								writer.println(line[1]+DELIMITER+line[5]+DELIMITER+line[6]+DELIMITER+line[7]+DELIMITER+line[2]+DELIMITER+freq+DELIMITER+maf);
							} else {
								writer.println(line[1]+DELIMITER+line[5]+DELIMITER+line[6]+DELIMITER+line[7]+DELIMITER+line[2]+DELIMITER+freq+DELIMITER+maf);
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
				trav++;
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
	
	public static void batch(String phenoFile, String pheno, String[] covars, String geneticDataTemplate, int startAt, boolean imputedNotGenotype, String pedfile, String outfileTemplate, String qsubRoot, String[] nodesToUse, int numBatches) {
		PrintWriter writer;
		String[] list;
		int count; //, step;
		Vector<String> v;
		
		if (qsubRoot == null) {
			qsubRoot = "batches/"+pheno+"_file#";
		} else if (qsubRoot.endsWith(".qsub")) {
			qsubRoot = qsubRoot.substring(0, qsubRoot.lastIndexOf("."));
		}		
		
		new File("batches/").mkdirs();
		count = startAt;		
		try {
			while (new File(ext.replaceAllWith(geneticDataTemplate, "#", count+"")).exists()) {
				writer = new PrintWriter(new FileWriter(ext.insertNumbers(qsubRoot, count)+".R"));
				writer.println("library(kinship)");
				writer.println("library(GWAF)");
				if (imputedNotGenotype) {
					writer.println("lme.batch.imputed(\""+phenoFile+"\", \""+ext.replaceAllWith(geneticDataTemplate, "#", count+"")+"\", \""+pedfile+"\", \""+pheno+"\", \"kmat.Rfile\", covars="+(covars==null?"NULL":"c(\""+Array.toStr(covars, "\",\"")+"\")")+", \""+ext.replaceAllWith(outfileTemplate, "#", count+"")+"\", col.names=T, sep.ped=\",\", sep.phe=\",\", sep.gen=\",\")");
				} else {
					writer.println("lme.batch(\""+phenoFile+"\", \""+ext.replaceAllWith(geneticDataTemplate, "#", count+"")+"\", \""+pedfile+"\", \""+pheno+"\", \"kmat.Rfile\", model=\"a\", covars="+(covars==null?"NULL":"c(\""+Array.toStr(covars, "\",\"")+"\")")+", \""+ext.replaceAllWith(outfileTemplate, "#", count+"")+"\", col.names=T, sep.ped=\",\", sep.phe=\",\", sep.gen=\",\")");
				}

				writer.close();
				count++;
			}
			count--;
			System.out.println("last file seen was '"+"gwaf/file"+count+".fhsR.gz'");
			if (numBatches < 1) {
				if (nodesToUse == null) {
					v = Array.toStringVector(Files.qsub("", qsubRoot, startAt, count, "R --no-save < "+qsubRoot+".R > "+qsubRoot+".log", null, -1, null));
				} else {
					v = new Vector<String>();
//					step = (int)Math.ceil((double)((count-startAt)+1)/(double)nodesToUse.length);
//					for (int i = 0; i < nodesToUse.length; i++) {
//						list = Files.qsub("", null, i*step+startAt, i==nodesToUse.length-1?count:((i+1)*step+startAt-1), "R --no-save < batches/"+pheno+"_gwaf#.R > batches/"+pheno+"_file#.log", "batches/"+pheno+"_file", null, -1, nodesToUse[i]);
					for (int i = startAt; i <= count; i++) {
						list = Files.qsub("", qsubRoot, i, i, "R --no-save < "+qsubRoot+".R > "+qsubRoot+".log", null, -1, nodesToUse[i%nodesToUse.length]);
						for (int j = 0; j < list.length; j++) {
							v.add(list[j]);
						}
					}
				}
				Files.writeList(Array.toStringArray(v), "master."+pheno);
				Files.chmod("master."+pheno);
			} else {
				Files.batchIt("run."+pheno, -1, startAt, count, numBatches, "R --no-save < batches/"+pheno+"_gwaf#.R > batches/"+pheno+"_file#.log");
			}
		} catch (Exception e) {
			System.err.println("Error queuing up file " + count);
			e.printStackTrace();
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

		String usage = "\n" + 
		"gwas.GWAF requires 0-1 arguments\n" + 
		"   (1) name of file with phenotypes (i.e. phenoFile=" + phenoFile + " (default; currently needs to be comma-delimited))\n" + 
		"   (2) name of phenotype to run (i.e. pheno=" + pheno + " (default))\n" + 
		"   (3) name of covariates to include (i.e. covars=" + covars + " (default; null leads to none; comma-delimited))\n" + 
		"   (4) format of genotype filenames (i.e. genoPrimer=" + geneticDataTemplate + " (default; can be zipped))\n" + 
		"   (5) number to start looking for file pattern (i.e. startAt=" + startAt + " (default))\n" +
		"   (6) data is imputed not genotyped (i.e. imputed=" + imputedNotGenotype + " (default))\n" +
		"   (7) name of pedigree file (i.e. pedfile=" + pedfile + " (default))\n" + 
		"   (8) template for output files (i.e. outfileTemplate=" + outfileTemplate + " (default))\n" +
		"   (9) template for qsub files (i.e. qsubRoot=[trait]_file#.qsub (default))\n" +
		"   (10) nodes to use (i.e. nodesToUse=" + nodesToUse + " (default; qsubs only; full names, comma-delimited))\n" + 
		"   (11) number of batches to create (i.e. numBatches=" + numBatches + " (default; anything less than 1 leads to qsubs being made))\n" +
		" OR\n" +
		"   (1) parse results (i.e. -parseResults (not the default))\n" +
		"   (2) template for output files (i.e. outfileTemplate=" + outfileTemplate + " (default))\n" + 
		"   (3) number to start looking for file pattern (i.e. startAt=" + startAt + " (default))\n" + 
		"   (4) data is imputed not genotyped (i.e. imputed=" + imputedNotGenotype + " (default))\n" +
		"   (5) name of parsed output file (i.e. out=" + outfile + " (default))\n" +
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
			
			if (parse) {
				parse(outfileTemplate, startAt, imputedNotGenotype, outfile);
			} else {
				batch(phenoFile, pheno, covars, geneticDataTemplate, startAt, imputedNotGenotype, pedfile, outfileTemplate, qsubRoot, nodesToUse, numBatches);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
