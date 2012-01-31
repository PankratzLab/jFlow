// -Xms1024M -Xmx1024M
package gwas;

import java.io.*;
import java.util.*;

import common.*;
import bioinformatics.Sequence;

public class Metal {
	public static final String[][] CONVERSION_REQS = { {"SNP", "Marker", "Name", "name"}, {"A1", "Allele"}, {"N", "NMISS"}, {"BETA", "ODDS", "OR"}, {"P", "pval", "p-val", "p-value"}};
	public static final String[][] SE_REQS = { {"SNP", "Marker", "Name", "MarkerName", "name"}, {"A1", "Allele1", "REF"}, {"A2", "Allele2", "OTHER"}, {"Effect", "beta", "beta_SNP_add"}, {"SE", "StdErr", "sebeta_SNP_add"}};
	public static final String[][] N_REQS = { {"SNP", "Marker", "Name", "MarkerName", "name"}, {"A1", "Allele1", "REF"}, {"A2", "Allele2", "OTHER"}, {"Effect", "Direction", "DIR"}, {"P", "pval", "p-val", "p-value"}, {"N", "NMISS"}};
	
	public static final String TEST = "TEST";
	public static final String[] VALID_ALLELES = {"A", "C", "G", "T", "I", "D"};
	public static final String[] NULL_ALLELES = {".", "-", "N", "NA"};
	public static final int STRAND_CONFIG_SAME = 1;
	public static final int STRAND_CONFIG_SAME_FLIPPED = 2;
	public static final int STRAND_CONFIG_OPPOSITE = 3;
	public static final int STRAND_CONFIG_OPPOSITE_FLIPPED = 4;
	public static final int STRAND_CONFIG_DIFFERENT_ALLELES = 5;
	public static final int STRAND_CONFIG_BOTH_NULL = 6;
	public static final int STRAND_CONFIG_SPECIAL_CASE = 7;
	public static final int[] FLIP = {1, 0};
	public static final String[] SUFFIXES = {"_A1", "_A2", "_freq", "_N", "_Rsq", "_effN"}; // this is the final output header, _effN is always computed

	public static final String[] TEXT = {"NADA", "STRAND_CONFIG_SAME", "STRAND_CONFIG_SAME_FLIPPED", "STRAND_CONFIG_OPPOSITE", "STRAND_CONFIG_OPPOSITE_FLIPPED", "STRAND_CONFIG_DIFFERENT_ALLELES", "STRAND_CONFIG_BOTH_NULL", "STRAND_CONFIG_SPECIAL_CASE"};
	
	public static void convert(String dir, String results, String test, String method, String freq, boolean useSE, boolean useMetalNomenclature, String outfile, boolean suppressWarning) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, header;
		Hashtable<String,String> hash = new Hashtable<String,String>();
		boolean logistic;
		int[] indices;
		int testIndex, seIndex;
		String modelAllele;
		String[] alleles, refAlleles;

		if (method.equalsIgnoreCase("logistic")) {
			logistic = true;
		} else if (method.equalsIgnoreCase("linear")) {
			logistic = false;
		} else {
			System.err.println("Error - '"+method+"' is not a valid method, looking for logistic/linear");
			logistic = false;
		}

		hash = HashVec.loadFileToHashString(dir+freq, 1, new int[] {2, 3}, "\t", true);

		try {
			reader = new BufferedReader(new FileReader(dir+results));
			writer = new PrintWriter(new FileWriter(outfile==null?dir+results+(useSE?".se":"")+".metal":outfile));
//			header = reader.readLine().trim().split("\t");
			header = reader.readLine().trim().split("[\\s]+");
			indices = ext.indexFactors(CONVERSION_REQS, header, false, true, true, true);
			if (useSE) {
				seIndex = ext.indexFactors(new String[][] {{"SE"}}, header, false, false, true, true)[0];
			} else {
				seIndex = -1;
			}
			testIndex = ext.indexOfStr(TEST, header);
			if (useMetalNomenclature) {
				writer.println("MarkerName\tAllele1\tAllele2\tWeight\tDirection\tP-value"+(useSE?"\tEffect\tStdErr":""));
			} else {
				writer.println("MARKER\tREF\tOTHER\tN\tDIR\tPVALUE"+(useSE?"\tbeta\tSE":""));
			}
			while (reader.ready()) {
//				line = reader.readLine().trim().split("\t");
				line = reader.readLine().trim().split("[\\s]+");
				if ((testIndex==-1||test==null||line[testIndex].equals(test))&&!line[indices[3]].equals("NA")) {
					if (hash.containsKey(line[indices[0]])) {
						modelAllele = line[indices[1]];
						refAlleles = hash.get(line[indices[0]]).split("\t");
						if (modelAllele.equals(refAlleles[0])) {
							alleles = new String[] {refAlleles[0], refAlleles[1]};
						} else if (modelAllele.equals(refAlleles[1])) {
							alleles = new String[] {refAlleles[1], refAlleles[0]};
						} else if (modelAllele.equals(Sequence.flip(refAlleles[0]))) {
							alleles = new String[] {modelAllele, Sequence.flip(refAlleles[1])};
						} else if (modelAllele.equals(Sequence.flip(refAlleles[1]))) {
							alleles = new String[] {modelAllele, Sequence.flip(refAlleles[0])};
						} else {
							System.err.println("Error - could not match up the model allele for "+line[indices[0]]+" ("+modelAllele+") with the alleles in the freq file ("+refAlleles[0]+"/"+refAlleles[1]+")");
							alleles = new String[] {"???", "???"};
						}
						writer.print(line[indices[0]]+"\t"+alleles[0]+"\t"+alleles[1]);
					} else {
						System.err.println("Error - no allele frequency for marker '"+line[indices[0]]+"'");
						writer.print(line[indices[0]]+"\t"+line[indices[1]]+"\tU");
					}
					writer.println("\t"+line[indices[2]]+"\t"+(line[indices[3]].equals("NA")||line[indices[3]].equals(".")?"?":(logistic?(Double.parseDouble(line[indices[3]])>1?"+":"-"):(Double.parseDouble(line[indices[3]])>0?"+":"-")))+"\t"+(line[indices[4]].equals(".")?"NA":line[indices[4]])+(useSE?"\t"+(logistic?Math.log(Double.parseDouble(line[indices[3]])):line[indices[3]])+"\t"+line[seIndex]:""));
				}
			}
			writer.close();
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+dir+results+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+dir+results+"\"");
			System.exit(2);
		}

		if (new File(dir+"ref."+results+".metal").exists()) {
			makeBatch(dir, results+".metal", "ref."+results+".metal");
		} else if (!suppressWarning){
			System.err.println("METAL batch file was not made because "+"ref."+results+".metal did not exist in "+dir);
		}
	}

	public static void makeBatch(String dir, String file1, String file2) {
		makeBatch(dir, file1, file2, "WEIGHT N", "WEIGHT N");
	}

	public static void makeBatch(String dir, String file1, String file2, int n1, int n2) {
		makeBatch(dir, file1, file2, "DEFAULTWEIGHT "+n1, "DEFAULTWEIGHT "+n2);
	}

	public static void makeBatch(String dir, String file1, String file2, String weight1, String weight2) {
		PrintWriter writer;

		try {
			writer = new PrintWriter(new FileWriter(dir+ext.rootRootOf(file1)+"_"+ext.rootRootOf(file2)+"_metal.batch"));
			writer.println("metal << EOT");
			writer.println("");
			writer.println("MARKER MARKER");
			writer.println("ALLELE REF OTHER");
			writer.println("EFFECT DIR");
			writer.println("PVALUE PVALUE");
			writer.println(weight1);
			writer.println("PROCESS "+file1);
			writer.println(weight2);
			writer.println("PROCESS "+file2);
			writer.println("OUTFILE "+ext.rootRootOf(file1)+"_"+ext.rootRootOf(file2)+".metal .out");
			writer.println("ANALYZE");
			writer.println("");
			writer.println("QUIT");
			writer.println("");
			writer.println("EOT");
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing metal batch file");
		}
	}
	
	public static void metaAnalyze(String dir, String[] filenames, String outputFile, boolean se, double[] defaultWeights, Logger log) {
		PrintWriter writer;
		String mappings;
		String[][] reqs;
		String[] header, trav, prev;
		int[] indices;

		try {
			writer = new PrintWriter(new FileWriter(dir+ext.rootOf(outputFile, false)+"_metal"+(se?"_se":"")+"_input.txt"));
			if (se) {
				mappings = "SCHEME STDERR\n";
				reqs = SE_REQS;
				if (defaultWeights != null) {
					log.reportError("Warning - included weights will be ignored, since the scheme was set to StdErr");
				}
			} else {
				mappings = "\n";
				reqs = N_REQS;
				if (defaultWeights.length != filenames.length) {
					log.reportError("Error - the number of weights included ("+defaultWeights.length+") did not match the number of filenames included ("+filenames.length+")");
					return;
				}
			}
			prev = Array.stringArray(reqs.length, "");

			for (int i = 0; i < filenames.length; i++) {
				if (!new File(dir+filenames[i]).exists()) {
					log.reportError("Error - file '"+filenames[i]+"' does not exist in directory '"+dir+"'");
					return;
				}
				header = Files.getHeaderOfFile(dir+filenames[i], "\t", log);
				indices = ext.indexFactors(reqs, header, false, true, true, log, false);
				if (!se && defaultWeights != null) {
					indices[5] = 0;
				}
				if (Array.min(indices) == -1) {
					log.reportError("Error parsing '"+dir+filenames[i]+"'");
					return;
				}
				trav = Array.subArray(header, indices);
				if (!trav[0].equals(prev[0])) {
					mappings += "MARKER "+trav[0]+"\n";
				}
				if (!trav[1].equals(prev[1]) || !trav[2].equals(prev[2])) {
					mappings += "ALLELE "+trav[1]+" "+trav[2]+"\n";
				}
				if (se) {
					if (!trav[3].equals(prev[3])) {
						mappings += "EFFECT "+trav[3]+"\n";
					}
					if (!trav[4].equals(prev[4])) {
						mappings += "STDERR "+trav[4]+"\n";
					}
				} else {
					if (!trav[3].equals(prev[3])) {
						mappings += "EFFECT "+trav[3]+"\n";
					}
					if (!trav[4].equals(prev[4])) {
						mappings += "PVALUE "+trav[4]+"\n";
					}
					if (defaultWeights != null) {
						mappings += "DEFAULTWEIGHT "+defaultWeights[i]+"\n";
					} else if (!trav[5].equals(prev[5])) {
						mappings += "WEIGHT "+trav[5]+"\n";
					}
				}
				
				if (mappings.length() > 3) {
					writer.println(mappings);
					mappings = "\n";
				}
				
				writer.println("PROCESS "+filenames[i]);
				prev = trav;				
			}
			writer.println();
			writer.println("OUTFILE "+outputFile+" .out");
			writer.println("ANALYZE");
			writer.println("");
			writer.println("QUIT");
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing metal batch file");
		}
		
		if (!CmdLine.run("metal < "+ext.rootOf(outputFile, false)+"_metal"+(se?"_se":"")+"_input.txt", dir, System.out, true)) {
			log.report("metal < "+ext.rootOf(outputFile, false)+"_metal"+(se?"_se":"")+"_input.txt");
		}
	}

	public static void checkStrand(String filename) {
		BufferedReader reader;
        PrintWriter writer;
        String trav;
        String[] line;
        Hashtable<String,Vector<String>> hash = new Hashtable<String,Vector<String>>();
        Vector<String> v;
        String[] header, roots, values;
        int[][] indices;
        int[] counts, order;
        int count;
        long time;
        
        try {
	        reader = new BufferedReader(new FileReader(filename));
	        writer = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+"_consensus.xln"));
	        header = reader.readLine().trim().split("[\\s]+");
	        v = new Vector<String>();
	        for (int i = 1; i<header.length; i++) {
	        	HashVec.addIfAbsent(header[i].substring(0, header[i].indexOf("_")), v);
            }
	        roots = Array.toStringArray(v);
	        indices = new int[roots.length][2];
	        for (int i = 0; i<roots.length; i++) {
	        	indices[i][0] = ext.indexOfStr(roots[i]+"_A", header);
	        	if (indices[i][0] == -1) {
	        		System.err.println("Error - expecting, but did not find, a column header called "+roots[i]+"_A");
	        	}
	        	indices[i][1] = ext.indexOfStr(roots[i]+"_B", header);
	        	if (indices[i][1] == -1) {
	        		System.err.println("Error - expecting, but did not find, a column header called "+roots[i]+"_B");
	        	}
            }
	        count = 0;
            time = new Date().getTime();
	        while (reader.ready()) {
	        	count++;
	        	if (count%1000 == 0) {
	        		System.out.println(count+"\t"+(new Date().getTime()-time)/(double)count);
	        	}
		        hash.clear();
	        	line = reader.readLine().trim().split("[\\s]+");
	        	
	        	for (int i = 0; i<roots.length; i++) {
	        		trav = line[indices[i][0]];
	        		if (ext.indexOfStr(trav, VALID_ALLELES) < ext.indexOfStr(line[indices[i][1]], VALID_ALLELES)) {
	        			trav += line[indices[i][1]];
	        		} else {
	        			trav = line[indices[i][1]]+trav;
	        		}
	        		HashVec.addToHashVec(hash, trav, roots[i], false);
                }
	        	values = HashVec.getKeys(hash);
	        	if (values.length == 1) {
		        	writer.println(line[0]+"\t"+values[0]);
	        	} else {
	        		counts = new int[values.length];
	        		for (int i = 0; i<values.length; i++) {
	        			counts[i] = hash.get(values[i]).size();
                    }
	        		order = Sort.quicksort(counts, Sort.DESCENDING);
		        	writer.print(line[0]+"\t"+values[order[0]]);
		        	for (int i = 1; i<values.length; i++) {
		        		writer.print("\t"+Array.toStr(Array.toStringArray(hash.get(values[order[i]])), ",")+"\t"+values[order[i]]);
                    }
		        	writer.println();
	        	}
	        }
	        reader.close();
            writer.close();
        } catch (FileNotFoundException fnfe) {
	        System.err.println("Error: file \""+filename+"\" not found in current directory");
	        System.exit(1);
        } catch (IOException ioe) {
	        System.err.println("Error reading file \""+filename+"\"");
	        System.exit(2);
        }
	}
	
	public static void splitCompFile(String filename, int numSplits) {
        Files.splitFile(filename, numSplits, 1, 1, ext.rootOf(filename, false), ".dat", false);
        Files.batchIt("strand", -1, 1, numSplits, 6, "jcp gwas.Metal check="+ext.rootOf(filename, false)+"#.dat");
	}

	public static void sortCompResults(String filename) {
		BufferedReader reader;
        PrintWriter writer, writer2;
        String[] line;
        
        try {
	        reader = new BufferedReader(new FileReader(filename));
	        writer = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+"_diffs.xln"));
	        writer2 = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+"_ambiguous.xln"));
	        while (reader.ready()) {
	        	line = reader.readLine().trim().split("[\\s]+");
	        	if (line.length > 2) {
	        		writer.println(Array.toStr(line));
	        	} else if (line[1].equals("AT") || line[1].equals("CG")) {
	        		writer2.println(Array.toStr(line));
	        	}
	        }
	        reader.close();
            writer.close();
            writer2.close();
        } catch (FileNotFoundException fnfe) {
	        System.err.println("Error: file \""+filename+"\" not found in current directory");
	        System.exit(1);
        } catch (IOException ioe) {
	        System.err.println("Error reading file \""+filename+"\"");
	        System.exit(2);
        }
	}
	
	public static void generateUniformsFromParamters(String filename, Logger log) {
		PrintWriter writer;
		Vector<String> params;
		String[] line;
		String newControlFile, indices;
		
		params = Files.parseControlFile(filename, "uniform", new String[] {"newControl.crf", "1 4=[pre]_A1 5=[pre]_A2 6=[pre]_freq 9=[pre]_N 7=[pre]_Rsq", "plink.EV/ARIC_Whites.results\taric_whites", "plink.EV/CARDIA_Whites.results\tcardia_whites", "R.EV/CHS_Whites.results\tchs_whites", "#generates control paramters such that each file that is listed uses the same indices and names them according to the designated prefix"}, log);
		if (params != null) {
			newControlFile = params.elementAt(0).trim();
			indices = params.elementAt(1).trim();
			try {
				writer = new PrintWriter(new FileWriter(newControlFile));
				writer.println("hits");
				writer.println("TBD");
				for (int j = 2; j < params.size(); j++) {
					line = params.elementAt(j).trim().split("[\\s]+");
					if (line.length > 2) {
						log.reportError("Error - requires exactly one or two paramters in lines 2 on: [filename and] prefix");
					} else if (line.length == 1) {
						writer.println(ext.replaceAllWith(indices, "[pre]", line[0]));
					} else {
						writer.println(line[0]+" "+ext.replaceAllWith(indices, "[pre]", line[1]));
					}
				}
				writer.close();
			} catch (Exception e) {
				System.err.println("Error writing to " + newControlFile);
				e.printStackTrace();
			}
		}
	}
	
	public static void generateInputFile(String filename, Logger log) {
		Vector<String> params;
		String[] files;
		
		params = Files.parseControlFile(filename, "metal", new String[] {"Discovery.se.metal", "Replication.se.metal", "outfile_root", "", "#Alternatively", "", "Discovery.assoc.logistic", "Replication.assoc.logistic", "outfile_root"}, log);

		if (params != null) {
			files = new String[] {"Discovery.se.metal", "Replication.se.metal", "something"};
			log.report(Array.toStr(Array.toStringArray(params)));
			for (int i = 0; i < files.length; i++) {
				if (params.size() >= i+1) {
					files[i] = params.elementAt(i).trim().split("[\\s]+")[0];
				}
			}
			Files.writeList(new String[] {"java -cp /home/npankrat/park.jar gwas.Metal test=ADD results="+files[0]+" method=logistic freq=plink.frq -se -metal out="+files[0]+".se.metal", "java -cp /home/npankrat/park.jar gwas.Metal test=ADD results="+files[1]+" method=logistic freq=plink.frq -se -metal out="+files[1]+".se.metal"}, ext.rootOf(filename)+"_convert.bat");
			Files.writeList(new String[] {"metal < "+ext.rootOf(filename)+"_metal_Nweighted.txt", "metal < "+ext.rootOf(filename)+"_metal_InvVar.txt"}, ext.rootOf(filename, false)+".bat");
//			Files.writeList(new String[] {"MARKER MARKER", "ALLELE REF OTHER", "WEIGHT N", "EFFECT DIR", "PVALUE PVALUE", "", "PROCESS "+files[0], "PROCESS "+files[1], "", "OUTFILE "+files[2]+".Nweighted .out", "ANALYZE", "", "QUIT"}, ext.rootOf(filename)+"_metal_Nweighted.txt");
//			Files.writeList(new String[] {"MARKER MARKER", "ALLELE REF OTHER", "EFFECT beta", "STDERR SE", "SCHEME STDERR", "GENOMICCONTROL OFF", "", "PROCESS "+files[0], "PROCESS "+files[1], "", "OUTFILE "+files[2]+".InvVar .out", "ANALYZE", "", "QUIT"}, ext.rootOf(filename)+"_metal_InvVar.txt");
			Files.writeList(new String[] {"MARKER MarkerName", "ALLELE Allele1 Allele2", "WEIGHT Weight", "EFFECT Direction", "PVALUE P-value", "", "PROCESS "+files[0], "PROCESS "+files[1], "", "OUTFILE "+files[2]+".Nweighted .out", "ANALYZE", "", "QUIT"}, ext.rootOf(filename)+"_metal_Nweighted.txt");
			Files.writeList(new String[] {"MARKER MarkerName", "ALLELE Allele1 Allele2", "EFFECT Effect", "STDERR StdErr", "SCHEME STDERR", "GENOMICCONTROL OFF", "", "PROCESS "+files[0], "PROCESS "+files[1], "", "OUTFILE "+files[2]+".InvVar .out", "ANALYZE", "", "QUIT"}, ext.rootOf(filename)+"_metal_InvVar.txt");
		}
	}
	
	public static void calculateWeightedAlleleFrequencyFromParamters(String filename, Logger log) {
		Vector<String> params;
		String[] line;
		String freqsFile;
		double rsqThresh, mafThresh;
		
		params = Files.parseControlFile(filename, "freq", new String[] {"freqs.xln Rsq>0.30 MAF>0.05", "", "#Looking for columns with the same prefix and the following suffixes (only Rsq is optional, which defaults to values of 1.0 -- i.e. genotyped/not imputed): _A1 _A2 _freq _N _Rsq"}, log);
		if (params != null) {
			rsqThresh = -1;
			mafThresh = -1;
			line = params.elementAt(0).trim().split("[\\s]+");
			freqsFile = line[0];
			for (int j = 1; j < line.length; j++) {
    			if (line[j].startsWith("Rsq>")) {
    				rsqThresh = Double.parseDouble(line[j].split(">")[1]);
    			}
    			if (line[j].startsWith("MAF>")) {
    				mafThresh = Double.parseDouble(line[j].split(">")[1]);
    			}
			}
    		calculateWeightedAlleleFrequency(freqsFile, rsqThresh, mafThresh, log);
		}
	}
	
	public static void calculateWeightedAlleleFrequency(String filename, double thresholdForRsq, double thresholdForMAFfile, Logger log) {
		BufferedReader reader;
        PrintWriter writer, w2, w3;
        String[] line;
        Vector<String> v;
        String[] header, roots;
        int[][] indices;
        int count;
        long time;
        double sumAlleles, sumEffN, effN, maf;
        String[] refAlleles;
        int[] countFlipped;
        
        try {
	        reader = new BufferedReader(new FileReader(filename));
	        writer = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+"_freq.xln"));
	        w2 = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+"_exceedingThreshold.dat"));
	        w3 = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+"_notExceedingThreshold.dat"));
	        writer.println("SNP\tA1\tA2\tfreqA1\teffN\tMAF");
	        header = reader.readLine().trim().split("[\\s]+");
	        v = new Vector<String>();
	        for (int i = 1; i<header.length; i++) {
	        	HashVec.addIfAbsent(header[i].substring(0, header[i].lastIndexOf("_")), v);
            }
	        roots = Array.toStringArray(v);
	        indices = new int[roots.length][5];
	        for (int i = 0; i<roots.length; i++) {
	        	for (int j = 0; j<SUFFIXES.length-1; j++) {
		        	indices[i][j] = ext.indexOfStr(roots[i]+SUFFIXES[j], header);
		        	if (j==4 && indices[i][j] == -1) {
		        		log.reportError("Warning - did not find an Rsq column header called "+roots[i]+SUFFIXES[j]+"; assuming a weight of 1 (i.e. genotyped/not imputed)");
		        	} else if (indices[i][j] == -1) {
		        		log.reportError("Error - expecting, but did not find, a column header called "+roots[i]+SUFFIXES[j]);
		        	}
                }
            }
	        count = 0;
            time = new Date().getTime();
            countFlipped = new int[roots.length];
	        while (reader.ready()) {
	        	count++;
	        	if (count%100000 == 0) {
	        		System.out.println(count+"\t"+ext.formDeci((new Date().getTime()-time)/(double)count, 4, true));
	        	}
//	        	line = reader.readLine().trim().split("[\\s]+");
	        	line = ext.replaceAllWith(reader.readLine(), ".-000", ".000") .trim().split("[\\s]+");
	        	
	        	sumAlleles = sumEffN = 0;
	        	refAlleles = new String[2];
	        	for (int i = 0; i<roots.length; i++) {
	        		// includes all instances where the freq is not missing and either the Rsq column is not present or the Rsq value is valid 
	        		if ((!line[indices[i][2]].equals(".") && !line[indices[i][2]].equals("NA")) && (indices[i][4] == -1 || (!line[indices[i][4]].equals(".") && !line[indices[i][4]].equals("NA") && Double.parseDouble(line[indices[i][4]]) > thresholdForRsq))) {
		        		effN = (indices[i][4] == -1?1:Math.min(Double.parseDouble(line[indices[i][4]]), 1))*Double.parseDouble(line[indices[i][3]]);
		        		if (effN < 0) {
		        			effN = 0;
		        		}
		        		switch (determineStrandConfig(new String[] {line[indices[i][0]].toUpperCase(), line[indices[i][1]].toUpperCase()}, refAlleles)) {
		        		case STRAND_CONFIG_SAME_FLIPPED:
		        			countFlipped[i]++;
		        		case STRAND_CONFIG_SAME:
		                    sumAlleles += Double.parseDouble(line[indices[i][2]])*effN;
		                    sumEffN += effN;
		                    break;
		        		case STRAND_CONFIG_OPPOSITE_FLIPPED:
		        			countFlipped[i]++;
		        		case STRAND_CONFIG_OPPOSITE:
		                    sumAlleles += (1-Double.parseDouble(line[indices[i][2]]))*effN;
		                    sumEffN += effN;
		                    break;
		        		case STRAND_CONFIG_BOTH_NULL:
		        			break;
		        		case STRAND_CONFIG_DIFFERENT_ALLELES:
		        			log.reportError("Error - "+roots[i]+" has different alleles ("+line[indices[i][0]]+"/"+line[indices[i][1]]+") than the rest ("+refAlleles[0]+"/"+refAlleles[1]+")");
		        			break;
		        		case STRAND_CONFIG_SPECIAL_CASE:
		        			log.reportError("Warning - Special case starting with "+roots[i]+": alleles ("+line[indices[i][0]]+"/"+line[indices[i][1]]+") where previous had only ("+refAlleles[0]+"/"+refAlleles[1]+")");
		                    break;
	                    default:
	                    	log.reportError("Error - unknown determineStrandConfig return code");
	                    	break;
	                    }
	        		}
                }
	        	maf = sumAlleles/sumEffN < 0.50?sumAlleles/sumEffN:1-sumAlleles/sumEffN;
	        	writer.println(line[0]+"\t"+(refAlleles[0]==null?".":refAlleles[0])+"\t"+(refAlleles[1]==null?".":refAlleles[1])+"\t"+(sumEffN > 0?ext.formDeci(sumAlleles/sumEffN, 4, true)+"\t"+ext.formDeci(sumEffN, 1, true)+"\t"+ext.formDeci(maf, 4, true):".\t0\t."));
	        	if (maf > thresholdForMAFfile) {
	        		w2.println(line[0]);
	        	} else {
	        		w3.println(line[0]);
	        	}
	        }
	        reader.close();
            writer.close();
            w2.close();
            w3.close();
            log.report("Number of SNPs flipped in each study:");
            for (int i = 0; i<roots.length; i++) {
            	log.report(roots[i]+": "+countFlipped[i]);
            }
        } catch (FileNotFoundException fnfe) {
	        log.reportError("Error: file \""+filename+"\" not found in current directory");
	        log.reportException(fnfe);
	        System.exit(1);
        } catch (IOException ioe) {
	        log.reportError("Error reading file \""+filename+"\"");
	        log.reportException(ioe);
	        System.exit(2);
        }
	}
	
	public static void compareAlleleFrequencies(String filename, double diffThreshold) {
		BufferedReader reader;
        PrintWriter writer;
        PrintWriter[] writers;
        String[] line;
        Vector<String> v;
        String[] header, roots;
        int[][] indices;
        int count;
        long time;
        String[] refAlleles;
        Logger log;
        int countFlipped;
        double[] freqs;
        boolean flipped, ambiguous;
        int numNull;
        
        log = new Logger(ext.rootOf(filename, false)+"_comp.log");
        try {
	        reader = new BufferedReader(new FileReader(filename));
	        writers = new PrintWriter[3];
	        writers[0] = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+"_comp.xln"));
	        writers[1] = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+"_compFlipped.xln"));
	        writers[2] = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+"_compAmbiguous.xln"));
	        header = reader.readLine().trim().split("[\\s]+");
	        v = new Vector<String>();
	        for (int i = 1; i<header.length; i++) {
	        	HashVec.addIfAbsent(header[i].substring(0, header[i].indexOf("_")), v);
            }
	        if (v.size() != 2) {
	        	System.err.println("Error - expecting two and only two roots to compare");
	        	return;
	        }
	        roots = Array.toStringArray(v);
	        indices = new int[roots.length][6];
	        for (int i = 0; i<roots.length; i++) {
	        	for (int j = 0; j<SUFFIXES.length; j++) {
		        	indices[i][j] = ext.indexOfStr(roots[i]+SUFFIXES[j], header);
		        	if (j<3 && indices[i][j] == -1) {
		        		log.reportError("Error - expecting, but did not find, a column header called "+roots[i]+SUFFIXES[j]);
		        	}
                }
            }
	        count = 0;
            time = new Date().getTime();
            countFlipped = 0;
            for (int w = 0; w<3; w++) {
            	writer = writers[w];
                writer.print("SNP\tref_A1\tref_A2");
                for (int i = 0; i<roots.length; i++) {
        			writer.print("\t"+roots[i]+SUFFIXES[2]);
                }
    			writer.print("\tflipped\tambiguous\tdiff\tdiff>"+diffThreshold+"_wFlip\tdiff>"+diffThreshold+"_woFlip");
                for (int i = 0; i<roots.length; i++) {
                	for (int j = 3; j<6; j++) {
                		if (indices[i][j] >= 0) {
                			writer.print("\t"+roots[i]+SUFFIXES[j]);
                		}
                    }
                }
                writer.println();
            }
	        while (reader.ready()) {
	        	count++;
	        	if (count%100000 == 0) {
	        		System.out.println(count+"\t"+ext.formDeci((new Date().getTime()-time)/(double)count, 4, true));
	        	}
	        	line = reader.readLine().trim().split("[\\s]+");
	        	
	        	refAlleles = new String[2];
	        	freqs = new double[] {Double.NaN, Double.NaN};
	        	flipped = false;
	        	for (int i = 0; i<roots.length; i++) {
	        		if (!line[indices[i][2]].equals(".") && !line[indices[i][2]].equals("NA")) {
		        		switch (determineStrandConfig(new String[] {line[indices[i][0]].toUpperCase(), line[indices[i][1]].toUpperCase()}, refAlleles)) {
		        		case STRAND_CONFIG_SAME_FLIPPED:
		        			countFlipped++;
		        			flipped = true;
		        		case STRAND_CONFIG_SAME:
		        			freqs[i] = Double.parseDouble(line[indices[i][2]]);
		                    break;
		        		case STRAND_CONFIG_OPPOSITE_FLIPPED:
		        			countFlipped++;
		        			flipped = true;
		        		case STRAND_CONFIG_OPPOSITE:
		        			freqs[i] = 1-Double.parseDouble(line[indices[i][2]]);
		                    break;
		        		case STRAND_CONFIG_BOTH_NULL:
		        			break;
		        		case STRAND_CONFIG_DIFFERENT_ALLELES:
		        			log.reportError("Error - "+roots[i]+" has different alleles ("+line[indices[i][0]]+"/"+line[indices[i][1]]+") than the rest ("+refAlleles[0]+"/"+refAlleles[1]+")");
		        			break;
		        		case STRAND_CONFIG_SPECIAL_CASE:
		        			log.reportError("Warning - Special case starting with "+roots[i]+": alleles ("+line[indices[i][0]]+"/"+line[indices[i][1]]+") where previous had only ("+refAlleles[0]+"/"+refAlleles[1]+")");
		                    break;
	                    default:
	                    	log.reportError("Error - unknown determineStrandConfig return code");
	                    	break;
	                    }
	        		}
                }
	            numNull = 0;
	            for (int i = 0; i<roots.length; i++) {
	            	if ((freqs[i]+"").equals("NaN")) {
	            		numNull++;
	            	}
	            }
				ambiguous = numNull == 0 && refAlleles[0].equals(Sequence.flip(refAlleles[1])); 
	        	for (int w = 0; w<3; w++) {
	        		writer = writers[w];
	        		if (w == 0 || (w == 1 && flipped) || (w == 2 && ambiguous)) {
			            writer.print(line[0]+"\t"+refAlleles[0]+"\t"+refAlleles[1]);
			            for (int i = 0; i<roots.length; i++) {
		            		writer.print("\t"+((freqs[i]+"").equals("NaN")?".":ext.formDeci(freqs[i], 4, true)));
			            }
						if (numNull > 0) {
							writer.print("\t.\t.\t.\t.\t.");
						} else {
							writer.print("\t"+(flipped?1:0));
							writer.print("\t"+(ambiguous?1:0));
							writer.print("\t"+ext.formDeci(Math.abs(freqs[0]-freqs[1]), 4, true));
							writer.print("\t"+(Math.abs(freqs[0]-freqs[1])>diffThreshold?1:0));
							writer.print("\t"+(Math.abs(freqs[0]-(flipped?1-freqs[1]:freqs[1]))>diffThreshold?1:0));
						}
			            for (int i = 0; i<roots.length; i++) {
			            	for (int j = 3; j<6; j++) {
			            		if (indices[i][j] >= 0) {
			            			writer.print("\t"+line[indices[i][j]]);
			            		}
			                }
			            }
			            writer.println();
	        		}
                }
	        }
	        reader.close();
	        for (int w = 0; w<3; w++) {
	        	writers[w].close();
            }
            log.report("Number of SNPs flipped: "+countFlipped);
        } catch (FileNotFoundException fnfe) {
	        log.reportError("Error: file \""+filename+"\" not found in current directory");
	        log.reportException(fnfe);
	        System.exit(1);
        } catch (IOException ioe) {
	        log.reportError("Error reading file \""+filename+"\"");
	        log.reportException(ioe);
	        System.exit(2);
        }
	}

	public static int determineStrandConfig(String[] alleles, String[] referenceAlleles) {
		String[] flipped;
		boolean[] nullChecks;
		int index;
		
		if (ext.indexOfStr(alleles[0], VALID_ALLELES)>=0 && ext.indexOfStr(alleles[1], VALID_ALLELES)>=0) {
			if (referenceAlleles[0] == null) {
				referenceAlleles[0] = alleles[0];
				referenceAlleles[1] = alleles[1];
				return STRAND_CONFIG_SAME;
			} else if (referenceAlleles[1] == null) {
				if (alleles[0].equals(referenceAlleles[0])) {
					referenceAlleles[1] = alleles[1];
//					return STRAND_CONFIG_SAME;
					return STRAND_CONFIG_SPECIAL_CASE;
				} else if (alleles[1].equals(referenceAlleles[0])) {
					referenceAlleles[1] = alleles[0];
//					return STRAND_CONFIG_OPPOSITE;
					return STRAND_CONFIG_SPECIAL_CASE;
				} else {
					flipped = new String[] {Sequence.flip(alleles[0]), Sequence.flip(alleles[1])};
					if (flipped[0].equals(referenceAlleles[0])) {
						referenceAlleles[1] = flipped[1];
//						return STRAND_CONFIG_SAME_FLIPPED;
						return STRAND_CONFIG_SPECIAL_CASE;
					} else if (flipped[1].equals(referenceAlleles[0])) {
						referenceAlleles[1] = flipped[0];
//						return STRAND_CONFIG_OPPOSITE_FLIPPED;
						return STRAND_CONFIG_SPECIAL_CASE;
					} else {
						return STRAND_CONFIG_DIFFERENT_ALLELES;
					}
				}
			} else {
				if (alleles[0].equals(referenceAlleles[0]) && alleles[1].equals(referenceAlleles[1])) {
					return STRAND_CONFIG_SAME;
				} else if (alleles[0].equals(referenceAlleles[1]) && alleles[1].equals(referenceAlleles[0])) {
					return STRAND_CONFIG_OPPOSITE;
				} else {
					flipped = new String[] {Sequence.flip(alleles[0]), Sequence.flip(alleles[1])};
					if (flipped[0].equals(referenceAlleles[0]) && flipped[1].equals(referenceAlleles[1])) {
						return STRAND_CONFIG_SAME_FLIPPED;
					} else if (flipped[0].equals(referenceAlleles[1]) && flipped[1].equals(referenceAlleles[0])) {
						return STRAND_CONFIG_OPPOSITE_FLIPPED;
					} else {
						return STRAND_CONFIG_DIFFERENT_ALLELES;
					}
				}
			}
		} else {
			nullChecks = new boolean[] {false, false};
			for (int i = 0; i<nullChecks.length; i++) {
				if (ext.indexOfStr(alleles[i], NULL_ALLELES) >= 0) {
					nullChecks[i] = true;
				} else if (ext.indexOfStr(alleles[i], VALID_ALLELES) == -1) {
					return STRAND_CONFIG_SPECIAL_CASE;
				}				
            }
			if (Array.booleanArraySum(nullChecks) == 1) {
				index = nullChecks[0]?1:0;
				if (referenceAlleles[0] == null) {
					referenceAlleles[0] = alleles[index];
					return index == 0?STRAND_CONFIG_SAME:STRAND_CONFIG_OPPOSITE;
				} else if (referenceAlleles[1] == null) {
					if (alleles[index].equals(referenceAlleles[0])) {
						return index == 0?STRAND_CONFIG_SAME:STRAND_CONFIG_OPPOSITE;
					} else {
						flipped = new String[] {Sequence.flip(alleles[index])};
						if (flipped[0].equals(referenceAlleles[0])) {
							return index == 0?STRAND_CONFIG_SAME_FLIPPED:STRAND_CONFIG_OPPOSITE_FLIPPED;
						} else {
							return STRAND_CONFIG_DIFFERENT_ALLELES;
						}
					}
				} else {
					if (alleles[index].equals(referenceAlleles[0])) {
						return index == 0?STRAND_CONFIG_SAME:STRAND_CONFIG_OPPOSITE;
					} else if (alleles[index].equals(referenceAlleles[1])) {
						return index == 1?STRAND_CONFIG_SAME:STRAND_CONFIG_OPPOSITE;
					} else {
						flipped = new String[] {Sequence.flip(alleles[index])};
						if (flipped[0].equals(referenceAlleles[0])) {
							return index == 0?STRAND_CONFIG_SAME_FLIPPED:STRAND_CONFIG_OPPOSITE_FLIPPED;
						} else if (flipped[0].equals(referenceAlleles[1])) {
							return index == 1?STRAND_CONFIG_SAME_FLIPPED:STRAND_CONFIG_OPPOSITE_FLIPPED;
						} else {
							return STRAND_CONFIG_DIFFERENT_ALLELES;
						}
					}
				}
			} else if (Array.booleanArraySum(nullChecks) == 2) {
				return STRAND_CONFIG_BOTH_NULL;
			} else {
				return STRAND_CONFIG_DIFFERENT_ALLELES;
			}
		}
	}
		
	public static void main(String[] args) {
		int numArgs = args.length;
		String results = "additive.assoc.logistic";
		String test = "ADD";
		String method = "logistic";
		String freq = "freq.frq";
//		String freq = "plink.frq";
//		String dir = "C:/Documents and Settings/npankrat/My Documents/gwas/udall/";
//		String dir = "C:/Documents and Settings/npankrat/My Documents/tWork/Consortium/analyses/Sing550/";
//		String dir = "C:/Documents and Settings/npankrat/My Documents/tWork/Consortium/analyses/";
		String dir = "";
		String file1 = "CIDR.metal";
		int weight1 = 1783;
		String file2 = "Fung.metal";
		int weight2 = 1000;
		boolean batch = false;
		boolean useSE = false;
		String strandFile = "";
		String fileToSplit = "";
		String calcFreq = "";
		String compFreq = "";
		int numSplits = 6;
		String fileToSort = "";
		double rsqThreshold = 0.3;
		double mafThreshold = 0.05;
		double diffThreshold = 0.20;
		String outfile = null;
		boolean useMetalNomenclature = true;
		String analyze = null;
		String output = "META_SE";
		boolean se = true;
		
//		strandFile = "C:/Documents and Settings/npankrat/My Documents/UMN/Folson/VTE_meta_analysis/finalAnalysis/16 assessing strand/snplist1_described3.xln";
//		fileToSort = "C:/Documents and Settings/npankrat/My Documents/UMN/Folson/VTE_meta_analysis/finalAnalysis/16 assessing strand/all_consensus.xln";
//		calcFreq = "C:/Documents and Settings/npankrat/My Documents/UMN/Folson/VTE_meta_analysis/finalAnalysis/16 assessing strand/compFreqs/snplist_described2.xln";
//		calcFreq = "C:/Documents and Settings/npankrat/My Documents/UMN/Folson/VTE_meta_analysis/finalAnalysis/16 assessing strand/compFreqs/snplist_described2_noARIC_noCHS.xln";
//		calcFreq = "C:/Documents and Settings/npankrat/My Documents/UMN/Folson/VTE_meta_analysis/finalAnalysis/16 assessing strand/compFreqs/snplist_described2_noARIC_noCHS.xln";
//		compFreq = "C:/Documents and Settings/npankrat/My Documents/UMN/Folson/VTE_meta_analysis/finalAnalysis/16 assessing strand/compFreqs/ARIC_comp.xln";
//		compFreq = "C:/Documents and Settings/npankrat/My Documents/UMN/Folson/VTE_meta_analysis/finalAnalysis/16 assessing strand/compFreqs/CHS_comp.xln";
//		calcFreq = "C:/Documents and Settings/npankrat/My Documents/UMN/Folson/VTE_meta_analysis/finalAnalysis/11 X chromosome/00src/snplist_described.xln";
//		calcFreq = "C:/Documents and Settings/npankrat/My Documents/tWork/Consortium/analysisOfImputation/alleleFreq/snplist_described.xln";
//		calcFreq = "D:/tWork/Consortium/Megas/alleleFreqInput.dat";
//		calcFreq = "D:/Myron/CARe/ICAM1/LOG_ICAM1/IBC/whites.xln";
//		calcFreq = "D:/mega/Consortium/Megas/alleleFreqInput.dat";

		String usage = "\n"+
		"gwas.Metal requires 0-1 arguments\n"+
		"   (0) directory (i.e. dir="+dir+" (default))\n"+
		"   (1) plink results file (i.e. results="+results+" (default))\n"+
		"   (2) TEST to parse (i.e. test="+test+" (default))\n"+
		"   (3) method (logistic/linear) (i.e. method="+method+" (default))\n"+
		"   (4) PLINK --freq output (i.e. freq="+freq+" (default))\n"+
		"   (5) include standard error (i.e. -se (not the default))\n"+
		"   (6) name of output file (i.e. out=[results_file].metal (default))\n"+
		"   (7) column headers in output use Metal nomenclature (i.e. metalNom="+useMetalNomenclature+" (default))\n"+
		" OR\n"+
		"   (1) create batch file (i.e. -batch (not the default))\n"+
		"   (2) file #1 (i.e. file1="+file1+" (default))\n"+
		"   (3) weight for file #1 (i.e. weight1="+weight1+" (default))\n"+
		"   (4) file #2 (i.e. file2="+file2+" (default))\n"+
		"   (5) weight for file #2 (i.e. weight2="+weight2+" (default))\n"+
		" OR\n"+
		"   (1) check strand (i.e. check=filename.dat (not the default))\n"+
		" OR\n"+
		"   (1) split file to be checked for strand (i.e. split=filename.dat (not the default))\n"+
		"   (2) number of files to split it into (i.e. numSplits="+numSplits+" (default))\n"+
		" OR\n"+
		"   (1) sort results of strand check (i.e. sort=filename.dat (not the default))\n"+
		" OR\n"+
		"   (1) calculate weighted allele frequency (i.e. calcFreq=filename.dat (not the default))\n"+
		"   (2) Rsq threshold (i.e. rsq="+rsqThreshold+" (default))\n"+
		"   (3) MAF threshold for separate file listing those markers greater than threshold (i.e. maf>"+mafThreshold+" (default))\n"+
		" OR\n"+
		"   (1) compare allele frequencies (i.e. compFreq=filename.dat (not the default))\n"+
		"   (2) allele frequency difference threshold (i.e. diff="+diffThreshold+" (default))\n"+
		" OR\n"+
		"   (1) meta-analyze a list of files (i.e. analyze=file1.metal,file2.metal,file3.metal (not the default))\n"+
		"   (2) name of output root (i.e. out="+output+" (default))\n"+
		"   (3) stderr scheme instead of sample-size weighting (i.e. se="+se+" (default))\n"+
		"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("results=")) {
				results = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("dir=")) {
				dir = ext.parseStringArg(args[i], "");
				numArgs--;
			} else if (args[i].startsWith("test=")) {
				test = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("method=")) {
				method = ext.parseStringArg(args[i], "logistic");
				numArgs--;
			} else if (args[i].startsWith("freq=")) {
				freq = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("-batch")) {
				batch = true;
				numArgs--;
			} else if (args[i].startsWith("file1=")) {
				file1 = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("weight1=")) {
				weight1 = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("file2=")) {
				file2 = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("weight2=")) {
				weight2 = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].toLowerCase().startsWith("-se")) {
				useSE = true;
				numArgs--;
			} else if (args[i].toLowerCase().startsWith("metalnom=")) {
				useMetalNomenclature = ext.parseBooleanArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("out=")) {
				outfile = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("check=")) {
				strandFile = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("rsq=")) {
				rsqThreshold = Double.parseDouble(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("maf>")) {
				mafThreshold = Double.parseDouble(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("split=")) {
				fileToSplit = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("numSplits=")) {
				numSplits = Integer.parseInt(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("sort=")) {
				fileToSort = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("calcFreq=")) {
				calcFreq = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("compFreq=")) {
				compFreq = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("diff=")) {
				diffThreshold = Double.parseDouble(args[i].split("=")[1]);
				numArgs--;
			} else if (args[i].startsWith("analyze=")) {
				analyze = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("se=")) {
				se = ext.parseBooleanArg(args[i]);
				numArgs--;
			} else {
				System.err.println("Error - don't know what to do with argument: '"+args[i]+"'");
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		
		try {
			if (batch) {
				makeBatch(dir, file1, file2, weight1, weight2);
			} else if (!fileToSplit.equals("")) {
				splitCompFile(fileToSplit, numSplits);
			} else if (!fileToSort.equals("")) {
				sortCompResults(fileToSort);
			} else if (!strandFile.equals("")) {
				checkStrand(strandFile);
			} else if (!calcFreq.equals("")) {
				calculateWeightedAlleleFrequency(calcFreq, rsqThreshold, mafThreshold, new Logger(ext.rootOf(calcFreq, false)+"_freq.log"));
			} else if (!compFreq.equals("")) {
				compareAlleleFrequencies(compFreq, diffThreshold);
			} else if (analyze != null) {
				metaAnalyze("./", analyze.split(","), output, se, null, new Logger());
			} else {
				convert(dir, results, test, method, freq, useSE, useMetalNomenclature, outfile, false);
			}
//			convert(dir, "IUhitsAllAdd.txt", null, "logistic", freq);
//			convert(dir, "IUhitsAllDom.txt", null, "logistic", freq);
//			convert(dir, "IUhitsAllRec.txt", null, "logistic", freq);
//			convert(dir, "UdallhitsAllAdd.txt", null, "logistic", freq);
//			convert(dir, "UdallhitsAllDom.txt", null, "logistic", freq);
//			convert(dir, "UdallhitsAllRec.txt", null, "logistic", freq);
//			makeBatch(dir, "IUhitsAllAdd.txt.metal", "UdallhitsAllAdd.txt.metal", 1783, 1000);
//			makeBatch(dir, "IUhitsAllDom.txt.metal", "UdallhitsAllDom.txt.metal", 1783, 1000);
//			makeBatch(dir, "CIDR.metal", "Fung.metal", 1783, 1000);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
