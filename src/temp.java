import java.io.*;
//import java.util.*;
import common.*;
import parse.GenParser;
import stats.LogisticRegression;
import stats.Nonparametric;
import stats.Ttest;

public class temp {
	public static void parseAll(String filename) {
		BufferedReader reader;
		PrintWriter writer;
		String temp;
		String root;
		String annotation;
		
		try {
			reader = new BufferedReader(new FileReader(filename));
			writer = new PrintWriter(new FileWriter(ext.rootOf(filename)+"_parsed.out"));
			root = "";
			annotation = "nada";
			while (reader.ready()) {
				temp = reader.readLine();
				if (temp.startsWith("rootLink=")) {
					root = temp.substring(9);
				} else if (temp.startsWith("<tr>")) {
					annotation = "";
					reader.readLine();
					reader.readLine();
					temp = reader.readLine().trim();
					temp = temp.substring(4, temp.length()-5);
					annotation += temp;
					reader.readLine();
					temp = reader.readLine().trim();
					temp = temp.substring(4, temp.length()-5);
					if (temp.length() > 0) {
						annotation += " ("+temp+")";
					}
				} else if (temp.indexOf(".mp4'>") > 0) {
					temp = temp.substring(0, temp.indexOf(".mp4'>")+4);
					temp = temp.substring(temp.lastIndexOf("href='")+6);
					writer.println(ext.link(root, temp)+"\t"+ext.replaceWithLinuxSafeCharacters(annotation+" "+ext.removeDirectoryInfo(temp), false));
				}
				
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			System.exit(2);
		}
	}

	public static void downloadAll(String filename, String dir) {
		String[][] files;
		
		files = HashVec.loadFileToStringMatrix(filename, false, new int[] {0,1}, "\t", false, 1000, false);
		
		for (int i = 0; i < files.length; i++) {
			Internat.downloadFile(files[i][0], dir+files[i][1]);
		}

	}
	
	public static void script() {
		for (int chr = 1; chr <= 22; chr++) {
			System.out.println("tar -zxvf chr"+chr+".tar.gz");
			System.out.println("cp list.txt chr"+chr+"/");
			System.out.println("cd chr"+chr+"/");
			System.out.println("gzip chr"+chr+".dose");
			System.out.println("gzip chr"+chr+".info");
			System.out.println("cd ..");
		}
	}

	public static void runPennCnv(String infile) {
		String[][] trios;
		String commands;
		
		trios = HashVec.loadFileToStringMatrix(infile, true, new int[] {4,5,6}, "\t", false, 100, false);
		for (int i = 0; i < trios.length; i++) {
			boolean fine = true;
			for (int j = 0; j < 3; j++) {
//				if (!Files.exists("D:/BOSS/TriosSamples/penn_data/"+trios[i][j], false)) {
				if (!Files.exists("penn_data/"+trios[i][j], false)) {
//				if (!Files.exists(trios[i][j], false)) {
					trios[i][j] +="*";
					fine = false;
				}
			}
			if (!fine) {
				System.out.println(i+"\t"+Array.toStr(trios[i]));
			}
		}
//		commands = "/home/npankrat/bin/detect_cnv.pl -joint -hmm /home/npankrat/bin/lib/hh550.hmm -pfb ../custom.pfb -gcmodel ../custom.gcmodel [%1] [%2] [%0] -out ../joint_results/[%0].jointcnv 2> ../joint_results/[%0].log";
//		Files.qsub("runJoint", commands, trios);

		commands = "/home/npankrat/bin/detect_cnv.pl -trio -hmm /home/npankrat/bin/lib/hh550.hmm -pfb ../custom.pfb -gcmodel ../custom.gcmodel -cnv penncnv.rawcnv [%1] [%2] [%0] -out ../trio_results/[%0].triocnv 2> ../trio_results/[%0].log";
		Files.qsub("runTrio", commands, trios);
	}
	
	public static void moveSamplesToDifferentFolder() {
		String[] samples;
		String dir;
		
		dir = "D:/data/GEDI/";
		samples = HashVec.loadFileToStringArray(dir+"samplesThatArePartOfTrios.txt", false, new int[] {0}, false);
		
		new File(dir+"penn_data/trios/").mkdirs();
		for (int i = 0; i < samples.length; i++) {
			new File(dir+"penn_data/"+samples[i]+".gz").renameTo(new File(dir+"penn_data/trios/"+samples[i]+".gz"));
		}
		
		String[] list, commands;
		
		list = HashVec.loadFileToStringArray("C:/GEDI_exome/00src/SentrixBarcode_list.txt", false, new int[] {0}, false);
		commands = new String[list.length];
		for (int i = 0; i < commands.length; i++) {
			commands[i] = "move "+list[i]+".zip batch1";
		}
		Files.writeList(commands, "C:/GEDI_exome/00src/move.out");
		
		
	}
	
	public static void schaid(String infile, String timeToMatch) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		
		try {
			reader = new BufferedReader(new FileReader(infile));
			writer = new PrintWriter(new FileWriter(infile+".out"));
//			line = reader.readLine().trim().split("[\\s]+");
//			for (int i = 1; i < line.length; i++) {
//				writer.print((i==1?"":"\t")+line[i]);
//			}
//			writer.print("\n");
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				if (line.length > 1) {
					for (int i = 1; i < line.length; i++) {
						writer.print((i==1?"":"\t")+line[i]);
					}
					writer.print("\n");
				}
			}
			reader.close();
			writer.close();
			new File(infile+".out").setLastModified(new File(timeToMatch).lastModified());
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + infile + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + infile + "\"");
			System.exit(2);
		}
	}

	
	public static void parseAllFilesInDorectory(String dir) {
//		search for .xln extension, run all
		
//		determine top hits for all models
//		args = model1.out out=hits1.txt 'MarkerName' !'P-value'<0.001
		
//		combine lists, find unique
//		Files.cat(originalFiles, finalFile, skips, log)
//		uniqueValues = Array.unique(array), write to file for our records
		
//		parse unique top hits from results files
//		Files.combine(keys, fileParameters, headers, unit, outputFilename, log, ignoreCase, finalHeader, hideIndex, outIsCommaDelimited)
//		Files.combine(uniqueValues, fileParameters, null, "MarkerName", outputFilename, log, true, true, true, false)

//		String[] fileParameters = new String[] {"model1.out 'MarkerName' 'Chr' 'Position' 'Effect_allele' 'Reference_allele' 'N' 'Effect_allele_frequency' 'BETA'=beta_M1 'P-value'=pval_M1",
//												"model2.out 'MarkerName' 'N' 'Effect_allele_frequency' 'BETA'=beta_M2 'P-value'=pval_M2",
//												"model3.out 'MarkerName' 'N' 'Effect_allele_frequency' 'BETA'=beta_M3 'P-value'=pval_M3",
//												"model4.out 'MarkerName' 'N' 'Effect_allele_frequency' 'BETA'=beta_M4 'P-value'=pval_M4"
//				};
	}
	
	public static void parseIt(String filename, String locationOfGenotypeFiles) {
		Logger log;
		String[] header;
		String temp;
		String[] args;
		
		log = new Logger(ext.rootOf(filename, false)+".log");	
		temp = Files.getFirstNLinesOfFile(filename, 1, log)[0];
		header = temp.trim().split("[\\s]+");
		
		args = new String[] {filename, "out=pheno1.dat", "'ID'=FID", "'ID'=IID", "'"+header[1]+"'", "tab", "replace=."};
		GenParser.parse(args, log);
		
		
		args = new String[header.length+4];
//		add filename, out=covars1.dat, FID, IID
//		leave out phenotype
//		leave out sex if it exists (make a flag to use sex)
//		at end, add tab and replace		
		
		GenParser.parse(args, log);
		
//		write batch (run1.bat) using locationOfGenotypeFiles
//		add sex if necessary
		
//		run batch
//		CmdLine.run("run1.bat", dir)
		
//		parse results
//		ResultsPackager.parseStdFormatFromPlink(dir, "D:/scratch/model1.assoc.linear", "ADD", "../filteredGenotypes/plink.bim", "D:/scratch/model1_freq.frq", null, outfile, log);
		
		
	}
	
	
	
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "source.txt";

		String usage = "\n"+
		".temp requires 0-1 arguments\n"+
		"   (1) filename (i.e. file=" + filename + " (default))\n"+
		"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
//			Internat.downloadFile("http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2708794/pdf/JPATH175000054.pdf", "D:/SH3GL2.pdf");
			
//			parseAll(filename);
//			downloadAll(ext.rootOf(filename)+"_parsed.out", "C:/Ezgi/");
//			downloadAll("catchup.txt", "C:/Ezgi/");

//			script();
//			runPennCnv("TriosForDenovoCnv.txt");
//			schaid("D:/tWork/LOAD/Schaid_request/pruned/load_fancy_eigens.xln");
//			schaid("D:/tWork/LOAD/Schaid_request/pruned/load_norm_fancy_eigens.xln");
//			schaid("D:/tWork/LOAD/Schaid_request/all_snps/all_load_fancy_eigens.xln");
//			schaid("D:/tWork/LOAD/Schaid_request/all_snps/all_load_fancy_postnormed_eigens.xln");
//			schaid("D:/tWork/LOAD/Schaid_request/pruned_snps/load_fancy_eigens.xln", "D:/tWork/LOAD/Schaid_request/pruned_snps/all_load_fancy_eigens.out");
//			schaid("D:/tWork/LOAD/Schaid_request/pruned_snps/load_fancy_postnormed_eigens.xln", "D:/tWork/LOAD/Schaid_request/pruned_snps/all_load_normalized_fancy_eigens.out");
			
//			moveSamplesToDifferentFolder();
			
			testEffectOfLNontype1error(100, 1.0);
			testEffectOfLNontype1error(100, 1.0);

			testEffectOfLNontype1error(1000, 1.0);
			testEffectOfLNontype1error(1000, 1.0);
			
			testEffectOfLNontype1error(100, 1.3);
			testEffectOfLNontype1error(100, 1.3);
			testEffectOfLNontype1error(100, 1.3);
			testEffectOfLNontype1error(100, 1.3);
			testEffectOfLNontype1error(100, 1.3);
			testEffectOfLNontype1error(100, 1.3);
			testEffectOfLNontype1error(100, 0.78);
			testEffectOfLNontype1error(100, 0.78);
			testEffectOfLNontype1error(100, 0.78);
			testEffectOfLNontype1error(100, 0.78);
			testEffectOfLNontype1error(100, 0.78);
			testEffectOfLNontype1error(100, 0.78);
			
			
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private static void testEffectOfLNontype1error(int sizePerRep, double multiplier) {
		int[] counts, tCounts, nCounts;
		double[][] betas;
		int reps;
		int[] deps;
		double[] indeps, indepsTransformed;
		LogisticRegression lr;
		Logger log;
		
		log = new Logger();
		
		System.out.print("x"+multiplier+" n="+sizePerRep+" ");
		
		reps = 100000;
//		sizePerRep = 1000;
		counts = new int[] {0, 0};
		tCounts = new int[] {0, 0};
		nCounts = new int[] {0, 0};
		betas = new double[2][reps];
		deps = new int[sizePerRep];
		indeps = new double[sizePerRep];
		indepsTransformed = new double[sizePerRep];
		for (int i = 0; i < reps; i++) {
			if (i % (reps/10) == 0) {
				System.out.print(".");
			}
			for (int j = 0; j < sizePerRep; j++) {
				deps[j] = Math.random()<0.2?1:0;
				indeps[j] = Math.random()*3;
				if (deps[j]>0.5) {
					indeps[j] = indeps[j]*multiplier;
				}
				indepsTransformed[j] = Math.exp(indeps[j]);
			}
//			System.out.println(Array.mean(indeps)+"\t"+Array.min(indeps)+"\t"+Array.max(indeps));
//			System.out.println(Array.mean(indepsTransformed)+"\t"+Array.min(indepsTransformed)+"\t"+Array.max(indepsTransformed));
			lr = new LogisticRegression(deps, indeps);
			betas[0][i] = lr.getBetas()[1];
			if (lr.getSigs()[1] < 0.05) {
				counts[0]++;
			}
			if (new Ttest(deps, indeps).getPvalue() < 0.05) {
				tCounts[0]++;
			}
			if (Nonparametric.runWilcoxonRankSumTest(deps, indeps, log) < 0.05) {
				nCounts[0]++;
			}
			lr = new LogisticRegression(deps, indepsTransformed);
			betas[1][i] = lr.getBetas()[1];
			if (lr.getSigs()[1] < 0.05) {
				counts[1]++;
			}
			if (new Ttest(deps, indepsTransformed).getPvalue() < 0.05) {
				tCounts[1]++;
			}
			if (Nonparametric.runWilcoxonRankSumTest(deps, indepsTransformed, log) < 0.05) {
				nCounts[1]++;
			}
		}
		System.out.println();
		System.out.println(counts[0]+"\t"+(double)counts[0]/reps); //+"\tmean beta: "+Array.mean(betas[0]));
		System.out.println(counts[1]+"\t"+(double)counts[1]/reps); //+"\tmean beta: "+Array.mean(betas[1]));
		System.out.println();
		System.out.println(tCounts[0]+"\t"+(double)tCounts[0]/reps);
		System.out.println(tCounts[1]+"\t"+(double)tCounts[1]/reps);
		System.out.println();
		System.out.println(nCounts[0]+"\t"+(double)nCounts[0]/reps);
		System.out.println(nCounts[1]+"\t"+(double)nCounts[1]/reps);
		
	}
}
