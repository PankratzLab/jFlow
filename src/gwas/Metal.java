// -Xms1024M -Xmx1024M -Xmx12G
package gwas;

import java.io.*;
import java.util.*;

import common.*;
import filesys.Hits;
import bioinformatics.MapSNPsAndGenes;
import bioinformatics.Sequence;

public class Metal {
	public static final String[][] CONVERSION_REQS = { {"SNP", "Marker", "Name", "name"}, {"A1", "Allele"}, {"N", "NMISS"}, {"BETA", "ODDS", "OR"}, {"P", "pval", "p-val", "p-value"}};

	public static final int SE_ANALYSIS = 0;
	public static final int WEIGHTED_SE_ANALYSIS = 1;
	public static final int PVAL_ANALYSIS = 2;
	
	public static final String[][][] REQS = {
		{ Aliases.EFFECTS, Aliases.STD_ERRS },
		{ {"wbeta"}, {"wse"} },
		{ {"beta", "Direction", "Effect", "DIR"}, Aliases.PVALUES, Aliases.NS },
	};

	public static final String[][] FREQS = {
		Aliases.ALLELE_FREQS,
		{"sampleAA", "fAllele11"},
		{"sampleAR", "fAllele12"},		
		{"sampleRR", "fAllele22"},		
	};
	
	public static final String TEST = "TEST";
	public static final String[] VALID_ALLELES = {"A", "C", "G", "T", "I", "D"};
	public static final String[] NULL_ALLELES = {".", "-", "N", "NA", "0"};
	public static final int STRAND_CONFIG_SAME_ORDER_SAME_STRAND = 1;
	public static final int STRAND_CONFIG_SAME_ORDER_FLIPPED_STRAND = 2;
	public static final int STRAND_CONFIG_OPPOSITE_ORDER_SAME_STRAND = 3;
	public static final int STRAND_CONFIG_OPPOSITE_ORDER_FLIPPED_STRAND = 4;
	public static final int STRAND_CONFIG_DIFFERENT_ALLELES = 5;
	public static final int STRAND_CONFIG_BOTH_NULL = 6;
	public static final int STRAND_CONFIG_SPECIAL_CASE = 7;
	public static final int[] FLIP = {1, 0};
	public static final String[] SUFFIXES = {"_A1", "_A2", "_freq", "_N", "_Rsq", "_effN"}; // this is the final output header, _effN is always computed

	public static final String[] TEXT = {"NADA", "STRAND_CONFIG_SAME", "STRAND_CONFIG_SAME_FLIPPED", "STRAND_CONFIG_OPPOSITE", "STRAND_CONFIG_OPPOSITE_FLIPPED", "STRAND_CONFIG_DIFFERENT_ALLELES", "STRAND_CONFIG_BOTH_NULL", "STRAND_CONFIG_SPECIAL_CASE"};

	public static void convertPlinkResults(String dir, String results, String test, String method, String freq, boolean useSE, boolean useMetalNomenclature, String outfile, boolean suppressWarning) {
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
//			System.err.println("METAL batch file was not made because "+"ref."+results+".metal did not exist in "+dir);
		}
	}
	
	public static void reformatResults(String filename, String[] unitOfAnlaysis, int defaultN, String output, Logger log) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		String delimiterIn, delimiterOut;
		String[] header;
		int[] markerIndices, alleleIndices, freqIndices;
		int[][] reqIndices;
		double freq;
		Hashtable<String,String> hash;
		int nSNPsIndex;
		
		delimiterIn = Files.determineDelimiter(filename, log);
		delimiterOut = Files.suggestDelimiter(output, log);
		
		try {
			reader = new BufferedReader(new FileReader(filename));
			writer = new PrintWriter(new FileWriter(output));

			hash = new Hashtable<String, String>();
			header = reader.readLine().trim().split(delimiterIn);
			
			markerIndices = ext.indexFactors(new String[][] {unitOfAnlaysis}, header, false, true, false, log, false);
			if (markerIndices[0] == -1) {
				log.reportError("Error - no unit of analysis in file "+filename+" ("+Array.toStr(unitOfAnlaysis, "/")+")");
				reader.close();
				writer.close();
				return;
			}

			alleleIndices = ext.indexFactors(Aliases.ALLELES, header, false, true, false, log, false);
			if (Array.min(alleleIndices) == -1 && Array.sum(alleleIndices) != -2) {
				log.reportError("Error - found only one allele in file "+filename+" (need "+Array.toStr(Aliases.ALLELES[0], "/")+" AND "+Array.toStr(Aliases.ALLELES[1], "/")+"); zeroing out and assuming uniform across studies...");
				alleleIndices[0] = -1;
				alleleIndices[1] = -1;
			}
			nSNPsIndex = ext.indexFactors(new String[][] {{"NALLELES", "num_variants"}}, header, false, true, false, log, false)[0];
			
			writer.print(unitOfAnlaysis[0]+delimiterOut+Aliases.ALLELES[0][0]+delimiterOut+Aliases.ALLELES[1][0]);
			reqIndices = new int[REQS.length][];
			for (int i = 0; i < REQS.length; i++) {
				reqIndices[i] = ext.indexFactors(REQS[i], header, false, true, false, log, false);
				if (i == PVAL_ANALYSIS && reqIndices[i][2] == -1 && defaultN >= 0) {
					reqIndices[i][2] = Integer.MAX_VALUE;
				}
				if (i == PVAL_ANALYSIS && ext.indexOfStr("mbpval", header) >= 0) {
					reqIndices[i][1] = ext.indexOfStr("mbpval", header);
				}
				if (i == SE_ANALYSIS && ext.indexOfStr("mbse", header) >= 0) {
					reqIndices[i][1] = ext.indexOfStr("mbse", header);
				}
				if (Array.min(reqIndices[i]) == -1) {
					reqIndices[i] = null;
				} else {
					for (int j = 0; j < reqIndices[i].length; j++) {
						if (hash.containsKey(REQS[i][j][0])) {
							reqIndices[i][j] = -1;
						} else {
							writer.print(delimiterOut+REQS[i][j][0]);
							hash.put(REQS[i][j][0], "");
						}
					}
				}
			}
			freqIndices = ext.indexFactors(FREQS, header, false, true, false, log, false);
			if (freqIndices[0] == -1 && Array.min(Array.subArray(freqIndices, 1)) == -1) {
				freqIndices = null;
			} else {
				writer.print(delimiterOut+FREQS[0][0]);
			}
			writer.println();
			while (reader.ready()) {
				line = reader.readLine().split(delimiterIn);
				if (nSNPsIndex == -1 || (!ext.isMissingValue(line[nSNPsIndex]) && Integer.parseInt(line[nSNPsIndex]) >= 5)) {
					writer.print(line[markerIndices[0]]);
					if (alleleIndices[1] == -1) {
						writer.print(delimiterOut+"G"+delimiterOut+"A");
					} else {
						writer.print(delimiterOut+line[alleleIndices[0]]+delimiterOut+line[alleleIndices[1]]);
					}
					for (int i = 0; i < reqIndices.length; i++) {
						if (reqIndices[i] != null) {
							for (int j = 0; j < reqIndices[i].length; j++) {
								if (reqIndices[i][j] != -1) {
									if (reqIndices[i][j] == Integer.MAX_VALUE) {
										writer.print(delimiterOut+defaultN);
									} else {
										writer.print(delimiterOut+line[reqIndices[i][j]]);
									}
								}
							}
						}
					}
					if (freqIndices != null) {
						if (freqIndices[0] != -1) {
							writer.print(delimiterOut+line[freqIndices[0]]);
						} else {
							try {
								freq = (Double.parseDouble(line[freqIndices[1]])*2+Double.parseDouble(line[freqIndices[2]])) / (Double.parseDouble(line[freqIndices[1]])*2+Double.parseDouble(line[freqIndices[2]])*2+Double.parseDouble(line[freqIndices[3]])*2);
							} catch (Exception e) {
								System.err.println("Error - invalid allele count ("+Array.subArray(line, Array.subArray(freqIndices, 1), "/")+") for marker "+line[markerIndices[0]]+" in file "+filename);
								freq = -999;
							}
							writer.print(delimiterOut+freq);
						}
					}
					writer.println();
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
			writer.println("ANALYZE HETEROGENEITY");
			writer.println("");
			writer.println("QUIT");
			writer.println("");
			writer.println("EOT");
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing metal batch file");
		}
	}
	
	public static void metaAnalyze(String dir, String[] filenames, String outputFile, boolean se, Logger log) {
		metaAnalyze(dir, filenames, Aliases.MARKER_NAMES, outputFile, SE_ANALYSIS, null, Array.doubleArray(filenames.length, -9), log);
	}
	
	public static void metaAnalyze(String dir, String[] filenames, String[] unitOfAnalysis, String outputFile, int analysisType, double[] defaultWeights, double[] gcValues/*boolean gcControlOn*/, Logger log) {
		PrintWriter writer;
		Vector<String> mappings;
		String[] header, travAlleles, prevAlleles, travReqs, prevReqs, travFreq, prevFreq;
		String travMarker, prevMarker;
		int[] indices;
		String filename;
		boolean travCommaDelimited, prevCommaDelimited;
		String batchFile;
		
//		filename = ext.rootOf(outputFile, false)+"_metal"+(anlaysisType == SE_ANALYSIS?"_se":(anlaysisType == WEIGHTED_SE_ANALYSIS?"_wse":""))+"_input.txt";
		filename = ext.rootOf(outputFile, false)+"_input.txt";

		try {
			writer = new PrintWriter(new FileWriter(dir+filename));
			prevMarker = "";
			prevAlleles = new String[] {"", ""};
			prevReqs = Array.stringArray(REQS[analysisType].length, "");
			prevCommaDelimited = false;
			travFreq = prevFreq = new String[] {"none"};
			
			if (gcValues == null || gcValues.length == 0) {
			   writer.println("GENOMICCONTROL OFF");
//			   writer.println("GENOMICCONTROL "+(gcControlOn?"ON":"OFF"));
			}
			if (analysisType == SE_ANALYSIS || analysisType == WEIGHTED_SE_ANALYSIS) {
				writer.println("SCHEME STDERR");
				if (defaultWeights != null) {
					log.reportError("Warning - included weights will be ignored, since the scheme was set to StdErr");
				}
			} else if (analysisType == PVAL_ANALYSIS) {
				writer.println("SCHEME SAMPLESIZE");
				if (defaultWeights != null && defaultWeights.length != filenames.length) {
					log.reportError("Error - the number of weights included ("+defaultWeights.length+") did not match the number of filenames included ("+filenames.length+")");
					writer.close();
					return;
				}
			} else {
				log.reportError("Error - unknown analysis type selected");
				writer.close();
				return;
			}

			for (int i = 0; i < filenames.length; i++) {
				mappings = new Vector<String>();

				if (!new File(dir+filenames[i]).exists()) {
					log.reportError("Error - file '"+filenames[i]+"' does not exist in directory '"+dir+"'");
					writer.close();
					return;
				}
				travCommaDelimited = filenames[i].endsWith(".csv");
				header = Files.getHeaderOfFile(dir+filenames[i], log);

				indices = ext.indexFactors(new String[][] {unitOfAnalysis}, header, false, true, true, log, false);
				if (indices[0] == -1) {
					log.reportError("Error parsing '"+dir+filenames[i]+"'");
					writer.close();
					return;
				}
				travMarker = header[indices[0]];
				if (!travMarker.equals(prevMarker)) {
					mappings.add("MARKER "+travMarker);
				}

				indices = ext.indexFactors(Aliases.ALLELES, header, false, true, true, log, false);
				if (Array.min(indices) == -1) {
					log.reportError("Error parsing '"+dir+filenames[i]+"'");
					writer.close();
					return;
				}
				travAlleles = Array.subArray(header, indices);
				if (!travAlleles[0].equals(prevAlleles[0]) || !travAlleles[1].equals(prevAlleles[1])) {
					mappings.add("ALLELE "+travAlleles[0]+" "+travAlleles[1]);
				}
				
				indices = ext.indexFactors(REQS[analysisType], header, false, true, true, log, false);
				if (analysisType == PVAL_ANALYSIS && defaultWeights != null) {
					indices[2] = 0;
				}
				if (Array.min(indices) == -1) {
					log.reportError("Error parsing '"+dir+filenames[i]+"'");
					writer.close();
					return;
				}
				travReqs = Array.subArray(header, indices);
				if (analysisType == SE_ANALYSIS || analysisType == WEIGHTED_SE_ANALYSIS) {
					if (!travReqs[0].equals(prevReqs[0])) {
						mappings.add("EFFECT "+travReqs[0]);
					}
					if (!travReqs[1].equals(prevReqs[1])) {
						mappings.add("STDERR "+travReqs[1]);
					}
				} else if (analysisType == PVAL_ANALYSIS) {
					if (!travReqs[0].equals(prevReqs[0])) {
						mappings.add("EFFECT "+travReqs[0]);
					}
					if (!travReqs[1].equals(prevReqs[1])) {
						mappings.add("PVALUE "+travReqs[1]);
					}
					if (defaultWeights != null) {
						mappings.add("DEFAULTWEIGHT "+defaultWeights[i]);
					} else if (!travReqs[2].equals(prevReqs[2])) {
						mappings.add("WEIGHT "+travReqs[2]);
					}

					indices = ext.indexFactors(FREQS, header, false, true, false, log, false);
					if (indices[0] != -1) {
						travFreq = Array.subArray(header, indices, "");
						if (!travFreq[0].equals(prevFreq[0])) {
							mappings.add("FREQLABEL "+travFreq[0]);
						}
						if (prevFreq[0].equals("none")) {
							mappings.add("AVERAGEFREQ ON");
						}
					}
				}
				
				if (travCommaDelimited != prevCommaDelimited) {
					if (travCommaDelimited) {
						mappings.add("SEPARATOR COMMA");
					} else {
						mappings.add("SEPARATOR WHITESPACE");
					}
				}
				
				if (mappings.size() > 0) {
					writer.println();
					for (int j = 0; j < mappings.size(); j++) {
						writer.println(mappings.elementAt(j));
					}
					writer.println();
				}
				
				if (gcValues[i] == -1) {
				    writer.println("GENOMICCONTROL OFF");
				    log.report("Setting GENOMICCONTROL to {OFF} for file {" + filenames[i] + "}");
				} else if (gcValues[i] == -9) {
				    writer.println("GENOMICCONTROL ON");
				    log.report("Setting GENOMICCONTROL to {ON} for file {" + filenames[i] + "}");
				} else {
				    writer.println("GENOMICCONTROL " + gcValues[i]);
				    log.report("Setting GENOMICCONTROL to {" + gcValues[i] + "} for file {" + filenames[i] + "}");
				}
				writer.println("PROCESS "+filenames[i]);
				prevMarker = travMarker;				
				prevAlleles = travAlleles;				
				prevReqs = travReqs;
				prevCommaDelimited = travCommaDelimited;
				prevFreq = travFreq;
			}
			writer.println();
			writer.println("OUTFILE "+outputFile+" .out");
			writer.println("ANALYZE HETEROGENEITY");
			writer.println("");
			writer.println("QUIT");
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing metal batch file");
			e.printStackTrace();
		}
		
		batchFile = "run_"+ext.rootOf(filename, false)+".bat";
		Files.write("metal < "+filename+" > "+ext.rootOf(filename, false)+".log", dir+batchFile);
		Files.chmod(dir+batchFile, false);
		if (!CmdLine.run("./" + batchFile, dir, System.out, false)) {
			log.report("metal < "+filename);
		}

//		int count = 0;
//		while (!Files.exists(outputFile+"1.out")) {
//			try {
//				Thread.sleep(500);
//			} catch (InterruptedException ie) {
//			}
//			count++;
//			if (count == 20) {
//				log.report("Still waiting for "+outputFile+"1.out to appear");
//			}
//		}
//		
//		long prevSize = -2;
//		long travSize = -1;
//		count = 0;
//		while (count < 4) {
//			travSize = new File(outputFile+"1.out").length();
//			if (travSize == prevSize) {
//				count++;
//			}
//			try {
//				Thread.sleep(500);
//			} catch (InterruptedException ie) {
//			}
//			prevSize = travSize;
//		}
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
				writer.println("lookup");
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
	
	public static void fromParameters(String filename, Logger log) {
		Vector<String> params;
		String[] inputFiles, tempFiles;
		String outputFile;
		Hits hits;
		String[] fileParameters;
		String[] header;
		int[] indices;
		BufferedReader reader;
		String[] line;
		Hashtable<String, int[]> markerPositionHash;
//		Vector<String> markers;
//		Vector<String> chrs;
//		Vector<String> positions;
//		Vector<String> geneNames;
		int[][] markerPositions;
		String[][] genes;
		PrintWriter writer;
		byte build;
		int[] chrPosition, trav;
		String[] hitList;
		int countMissing;
//		boolean gcControlOn;
		double thresholdForHits;
		double[] gcValues;
		int countMismatches;
		
		params = Files.parseControlFile(filename, "metal", new String[] {"outfile_root", "build=37",/* "genomic_control=TRUE",*/ "hits_p<=0.001", "file1.metal", "file2.txt", "file3.assoc.logistic"}, log);

		thresholdForHits = 0.001;
//		gcControlOn = true;
		build = -1;
		if (params != null) {
			outputFile = params.remove(0);
			for (int i = 0; i < params.size(); i++) {
				if (params.elementAt(i).startsWith("build=")) {
					build = ext.parseByteArg(params.elementAt(i));
					params.remove(i);
				}
//				if (params.elementAt(i).startsWith("genomic_control=")) {
//					gcControlOn = ext.parseBooleanArg(params.elementAt(i));
//					params.remove(i);
//				}
				if (params.elementAt(i).startsWith("hits_p<=")) {
					thresholdForHits = ext.parseDoubleArg(params.elementAt(i));
					params.remove(i);
				}
			}
			if (build == -1) {
				log.reportError("Warning - build was not specified, assuming build 37 (aka hg19)");
				build = 37;
			}
			tempFiles = Array.toStringArray(params);
			gcValues = Array.doubleArray(tempFiles.length, -9); // default to ON
			inputFiles = new String[tempFiles.length];//Array.toStringArray(params);
			for (int i = 0; i < tempFiles.length; i++) {
			    String[] parts = tempFiles[i].split("\t");
			    if (parts.length == 2) {
			        gcValues[i] = Double.parseDouble(parts[1]);
			    }
			    inputFiles[i] = parts[0];
			}
			
//			backupDir = "./backup";
//			Files.backup(outputFile+"_InvVar", null, backupDir);
//			Files.backup(outputFile+"_InvVar", null, backupDir);

//            String dir = ext.verifyDirFormat((new File("./")).getAbsolutePath());
			if (!Files.exists(outputFile+"_InvVar1.out")) {
	            log.report("Running inverse variance weighted meta-analysis...");
			    metaAnalyze("./", inputFiles, Aliases.MARKER_NAMES, outputFile+"_InvVar", SE_ANALYSIS, null, gcValues, log);
			} else {
			    log.report("Found inverse variance weighted meta-analysis results - skipping.");
			}
//			metaAnalyze(dir, inputFiles, Aliases.MARKER_NAMES, outputFile+"_InvVar", SE_ANALYSIS, null, gcControlOn, log);
			if (!Files.exists(outputFile+"_NWeighted1.out")) {
    			log.report("Running sample size weighted meta-analysis...");
    			metaAnalyze("./", inputFiles, Aliases.MARKER_NAMES, outputFile+"_NWeighted", PVAL_ANALYSIS, null, gcValues, log);
			} else {
			    log.report("Found sample size weighted meta-analysis results - skipping.");
			}
//			metaAnalyze(dir, inputFiles, Aliases.MARKER_NAMES, outputFile+"_NWeighted", PVAL_ANALYSIS, null, gcControlOn, log);
			
//			check to see if file exists, report error otherwise
//			Files.backup(filename, sourceDir, backupDir);
//			if (!Files.isFileReady(outputFile+"_InvVar.out", 5000, 2000) ||	!Files.isFileReady(outputFile+"_NWeighted.out", 5000, 2000)) {
//				System.err.println("Error - metaAnalysis gets hang or does not start.");
//			}
			
//			sort results for both, determine minimum-pvalue, report only those minP<0.001
			if (!Files.exists("hits.txt")) {
    			hits = new Hits();
    			hits.incorporateFromFile(outputFile+"_InvVar1.out", thresholdForHits, log);
//    			hits.incorporateFromFile(outputFile+"_NWeighted1.out", thresholdForHits, log);
    			hits.writeHits("hits.txt");
			} 
			
//			summarize data, Results/packager
//			markername, chr, postition, geneName
//			From InvVar results: Allele1, Allele2, Effect=Beta,StdErr,P-value,Direction,
//			From NWeighted results: Weight, P-value
//			From each individual input file: Allele1, Allele2, AlleleFreq, N, Beta, SE, P-value, Impute_info/RSq

//			markers = new Vector<String>();
//			chrs = new Vector<String>();
//			positions = new Vector<String>();
//			geneNames = new Vector<String>();
			markerPositionHash = new Hashtable<String, int[]>();

			countMismatches = 0;
			fileParameters = new String[4 + inputFiles.length];
			fileParameters[1] = "hits.txt 0 1=minPval skip=0";
			fileParameters[2] = outputFile + "_InvVar1.out 0 'Allele1' 'Allele2' 'Effect'=Beta 'StdErr' 'P-value' 'Direction'";
			fileParameters[3] = outputFile + "_NWeighted1.out 0 'Weight' 'P-value'";
			for (int i=0; i<inputFiles.length; i++) {
				header = Files.getHeaderOfFile(inputFiles[i], log);
				indices = ext.indexFactors(new String[][] {Aliases.ALLELES[0], Aliases.ALLELES[1], Aliases.ALLELE_FREQS, Aliases.NS, Aliases.EFFECTS, Aliases.STD_ERRS, Aliases.PVALUES, Aliases.IMPUTATION_EFFICIENCY}, header, true, false, true, true, log, false);
				fileParameters[i + 4] = inputFiles[i]+" 0";
				for (int j = 0; j < indices.length; j++) {
					if (indices[j] != -1) {
						fileParameters[i + 4] += " '"+header[indices[j]]+"'";
					}
				}

				indices = ext.indexFactors(new String[][] {Aliases.CHRS, Aliases.POSITIONS}, header, true, false, true, true, log, false);
				if (indices[0] != -1 && indices[0] != -1) {
					try {
						reader = Files.getAppropriateReader(inputFiles[i]);//new BufferedReader(new FileReader(inputFiles[i]));
						String hdr = reader.readLine(); // header
						String delim = ext.determineDelimiter(hdr);
						while (reader.ready()) {
							line = reader.readLine().split(delim);
							if (!ext.isMissingValue(line[indices[0]]) && !ext.isMissingValue(line[indices[1]])) {
								trav = new int[] {Integer.parseInt(line[indices[0]]), Integer.parseInt(line[indices[1]])};
								if (markerPositionHash.containsKey(line[0])) {
									chrPosition = markerPositionHash.get(line[0]);
									if (trav[0] != chrPosition[0] || trav[1] != chrPosition[1]) {
										if (countMismatches < 42) {
											log.reportError("Error - mismatched positions for marker "+line[0]+" ("+Array.toStr(trav, ":")+" versus "+Array.toStr(chrPosition, ":")+")");
										} else if (countMismatches == 42) {
											log.reportError("...");
										}
										countMismatches++;
									}
								}
								markerPositionHash.put(line[0], trav);
							}
						}
						reader.close();
					} catch (IOException e) {
						e.printStackTrace();
					}
				}
			}
			log.report("There were a total of "+countMismatches+" mismatches on position");
			
			
			hitList = HashVec.loadFileToStringArray("hits.txt", false, new int[] {0}, false);
			try {
				writer = new PrintWriter(new FileWriter("genes.txt"));
				writer.println("MarkerName\tChr\tPosition\tGene(s)");
				
				if (markerPositionHash.size()>0) {
					countMissing = 0;
					markerPositions = new int[hitList.length][];
					for (int i=0; i<hitList.length; i++) {
						if (markerPositionHash.containsKey(hitList[i])) {
							markerPositions[i] = markerPositionHash.get(hitList[i]);
						} else {
							markerPositions[i] = new int[] {0,i};
							countMissing++;
							if (countMissing < 10) {
								log.reportError("Warning - no map position available for "+hitList[i]);
							}
						}
					}
					if (countMissing >= 10) {
						log.reportError("No map positions for a total of "+countMissing+" markers");
					}
					System.out.println(ext.getTime()+"\tstarting genes");
					genes = MapSNPsAndGenes.mapSNPsToGenes(markerPositions, build, 50000, log);				
					System.out.println(ext.getTime()+"\tfinished genes");
					
					for (int i=0; i<genes.length; i++) {
						writer.println(hitList[i]+"\t"+Array.toStr(markerPositionHash.get(hitList[i]), "\t")+"\t"+genes[i][1]);
					}
				}
				writer.close();
			} catch (IOException e) {
				e.printStackTrace();
			}

			fileParameters[0] = "genes.txt" + " 0=MarkerName 1=Chr 2=Position 3=Gene(s)";
			Files.combine(hitList, fileParameters, null, "MarkerName", ".", "topHits.xln", log, true, true, false);

			String[][] results = HitWindows.determine("topHits.xln", 0.00000005f, 500000, 0.000005f, new String[0]);
			try {
                results = includeExtraInfoFromTopHits(results);
            } catch (IOException e) {
                System.err.println("Error - exception occured while incorporating topHits.xln data into topHitWindows.xln");
                e.printStackTrace();
            }
            Files.writeMatrix(results, "topHitWindows.out", "\t");
			
		}
	}
	
	private static String[][] includeExtraInfoFromTopHits(String[][] hwResults) throws IOException {
	    String[][] newResults = new String[hwResults.length][];
	    
	    int[] hwInd = ext.indexFactors(new String[][]{Aliases.MARKER_NAMES}, hwResults[0], false, true, false, false);
	    int mkrInd = hwInd[0];
	    
	    HashSet<String> hwMkrs = new HashSet<String>();
        for (int i = 1; i < hwResults.length; i++) {
	        hwMkrs.add(hwResults[i][mkrInd]);
	    }
        
	    HashMap<String, String[]> topHitsParts = new HashMap<String, String[]>();
	    
	    BufferedReader reader = Files.getAppropriateReader("topHits.xln");
	    String headerLine = reader.readLine().trim();
	    String delim = ext.determineDelimiter(headerLine);
	    String[] line = headerLine.split(delim);
	    int[] topInd = ext.indexFactors(new String[][]{Aliases.MARKER_NAMES, Aliases.CHRS, Aliases.POSITIONS}, line, false, true, false, false);
	    int topMkrInd = topInd[0];
	    String readLine = null;
	    while ((readLine = reader.readLine()) != null) {
	        line = readLine.split(delim);
	        if (hwMkrs.contains(line[topMkrInd])) {
	            topHitsParts.put(line[topMkrInd], line);
	        }
	    }
	    reader.close();
	    
	    String[] hdr = hwResults[0];
	    String[] pts = headerLine.split(delim);
	    String[] newHdr = new String[pts.length - 3];
	    int cnt = 0;
	    for (int i = 0; i < pts.length; i++) {
            if (i == topInd[0] || i == topInd[1] || i == topInd[2]) {
                continue;
            }
            newHdr[cnt] = pts[i];            
            cnt++;
	    }
	    newResults[0] = Array.concatAll(hdr, newHdr);
	    
	    for (int i = 1; i < hwResults.length; i++) {
	        String hwMkr = hwResults[i][mkrInd];
	        String[] topData = topHitsParts.get(hwMkr);
	        String[] newTopData = new String[topData.length - 3];
	        int cntInd = 0;
	        for (int j = 0; j < topData.length; j++) {
	            if (j == topInd[0] || j == topInd[1] || j == topInd[2]) {
	                continue;
	            }
	            newTopData[cntInd] = topData[j];
	            cntInd++;
	        }
	        
	        String[] newStringArray = Array.concatAll(hwResults[i], newTopData);
	        newResults[i] = newStringArray;
	    }
	    
	    return newResults;
	}
	
	
	public static void createUnionOfMap(String[] filenames, String mapOut) {

		
//		order = Sort.orderTwoLayers(chrs, positions);
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
			Files.writeList(new String[] {"MARKER MarkerName", "ALLELE Allele1 Allele2", "WEIGHT Weight", "EFFECT Direction", "PVALUE P-value", "", "PROCESS "+files[0], "PROCESS "+files[1], "", "OUTFILE "+files[2]+".Nweighted .out", "ANALYZE HETEROGENEITY", "", "QUIT"}, ext.rootOf(filename)+"_metal_Nweighted.txt");
			Files.writeList(new String[] {"MARKER MarkerName", "ALLELE Allele1 Allele2", "EFFECT Effect", "STDERR StdErr", "SCHEME STDERR", "GENOMICCONTROL OFF", "", "PROCESS "+files[0], "PROCESS "+files[1], "", "OUTFILE "+files[2]+".InvVar .out", "ANALYZE HETEROGENEITY", "", "QUIT"}, ext.rootOf(filename)+"_metal_InvVar.txt");
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
        String previousRef;
        
        try {
	        reader = new BufferedReader(new FileReader(filename));
	        writer = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+"_freq.xln"));
	        if (thresholdForMAFfile < 0) {
	        	w2 = null;
	        	w3 = null;
	        } else {
		        w2 = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+"_exceedingThreshold.dat"));
		        w3 = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+"_notExceedingThreshold.dat"));
	        }
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
			        	
			        	previousRef = refAlleles[0]+"/"+refAlleles[1];
		        		switch (determineStrandConfig(new String[] {line[indices[i][0]].toUpperCase(), line[indices[i][1]].toUpperCase()}, refAlleles)) {
		        		case STRAND_CONFIG_SAME_ORDER_FLIPPED_STRAND:
		        			countFlipped[i]++;
		        		case STRAND_CONFIG_SAME_ORDER_SAME_STRAND:
		                    sumAlleles += Double.parseDouble(line[indices[i][2]])*effN;
		                    sumEffN += effN;
		                    break;
		        		case STRAND_CONFIG_OPPOSITE_ORDER_FLIPPED_STRAND:
		        			countFlipped[i]++;
		        		case STRAND_CONFIG_OPPOSITE_ORDER_SAME_STRAND:
		                    sumAlleles += (1-Double.parseDouble(line[indices[i][2]]))*effN;
		                    sumEffN += effN;
		                    break;
		        		case STRAND_CONFIG_BOTH_NULL:
		        			break;
		        		case STRAND_CONFIG_DIFFERENT_ALLELES:
		        			log.reportError("Error - "+roots[i]+" has different alleles for "+line[0]+" ("+line[indices[i][0]]+"/"+line[indices[i][1]]+") than had been seen previously ("+previousRef+")");
		        			break;
		        		case STRAND_CONFIG_SPECIAL_CASE:
		        			log.reportError("Warning - Special case starting with "+roots[i]+": alleles for "+line[0]+" were "+line[indices[i][0]]+"/"+line[indices[i][1]]+" where only "+previousRef+" had been seem previously");
		                    break;
	                    default:
	                    	log.reportError("Error - unknown determineStrandConfig return code");
	                    	break;
	                    }
	        		}
                }
	        	maf = sumAlleles/sumEffN < 0.50?sumAlleles/sumEffN:1-sumAlleles/sumEffN;
	        	writer.println(line[0]+"\t"+(refAlleles[0]==null?".":refAlleles[0])+"\t"+(refAlleles[1]==null?".":refAlleles[1])+"\t"+(sumEffN > 0?ext.formDeci(sumAlleles/sumEffN, 4, true)+"\t"+ext.formDeci(sumEffN, 1, true)+"\t"+ext.formDeci(maf, 4, true):".\t0\t."));
	        	if (thresholdForMAFfile > 0) {
		        	if (maf > thresholdForMAFfile) {
		        		w2.println(line[0]);
		        	} else {
		        		w3.println(line[0]);
		        	}
	        	}
	        }
	        reader.close();
            writer.close();
	        if (thresholdForMAFfile > 0) {
	            w2.close();
	            w3.close();
	        }
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
		        		case STRAND_CONFIG_SAME_ORDER_FLIPPED_STRAND:
		        			countFlipped++;
		        			flipped = true;
		        		case STRAND_CONFIG_SAME_ORDER_SAME_STRAND:
		        			freqs[i] = Double.parseDouble(line[indices[i][2]]);
		                    break;
		        		case STRAND_CONFIG_OPPOSITE_ORDER_FLIPPED_STRAND:
		        			countFlipped++;
		        			flipped = true;
		        		case STRAND_CONFIG_OPPOSITE_ORDER_SAME_STRAND:
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

	// confusing terminology, here flipped means opposite strand, and opposite means flipped allele
	public static int determineStrandConfig(String[] alleles, String[] referenceAlleles) {
		String[] flipped;
		boolean[] nullChecks;
		int index;
		
		if (ext.indexOfStr(alleles[0], VALID_ALLELES)>=0 && ext.indexOfStr(alleles[1], VALID_ALLELES)>=0) {
			if (referenceAlleles[0] == null) {
				referenceAlleles[0] = alleles[0];
				referenceAlleles[1] = alleles[1];
				return STRAND_CONFIG_SAME_ORDER_SAME_STRAND;
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
					return STRAND_CONFIG_SAME_ORDER_SAME_STRAND;
				} else if (alleles[0].equals(referenceAlleles[1]) && alleles[1].equals(referenceAlleles[0])) {
					return STRAND_CONFIG_OPPOSITE_ORDER_SAME_STRAND;
				} else {
					flipped = new String[] {Sequence.flip(alleles[0]), Sequence.flip(alleles[1])};
					if (flipped[0].equals(referenceAlleles[0]) && flipped[1].equals(referenceAlleles[1])) {
						return STRAND_CONFIG_SAME_ORDER_FLIPPED_STRAND;
					} else if (flipped[0].equals(referenceAlleles[1]) && flipped[1].equals(referenceAlleles[0])) {
						return STRAND_CONFIG_OPPOSITE_ORDER_FLIPPED_STRAND;
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
					return index == 0?STRAND_CONFIG_SAME_ORDER_SAME_STRAND:STRAND_CONFIG_OPPOSITE_ORDER_SAME_STRAND;
				} else if (referenceAlleles[1] == null) {
					if (alleles[index].equals(referenceAlleles[0])) {
						return index == 0?STRAND_CONFIG_SAME_ORDER_SAME_STRAND:STRAND_CONFIG_OPPOSITE_ORDER_SAME_STRAND;
					} else {
						flipped = new String[] {Sequence.flip(alleles[index])};
						if (flipped[0].equals(referenceAlleles[0])) {
							return index == 0?STRAND_CONFIG_SAME_ORDER_FLIPPED_STRAND:STRAND_CONFIG_OPPOSITE_ORDER_FLIPPED_STRAND;
						} else {
							return STRAND_CONFIG_DIFFERENT_ALLELES;
						}
					}
				} else {
					if (alleles[index].equals(referenceAlleles[0])) {
						return index == 0?STRAND_CONFIG_SAME_ORDER_SAME_STRAND:STRAND_CONFIG_OPPOSITE_ORDER_SAME_STRAND;
					} else if (alleles[index].equals(referenceAlleles[1])) {
						return index == 1?STRAND_CONFIG_SAME_ORDER_SAME_STRAND:STRAND_CONFIG_OPPOSITE_ORDER_SAME_STRAND;
					} else {
						flipped = new String[] {Sequence.flip(alleles[index])};
						if (flipped[0].equals(referenceAlleles[0])) {
							return index == 0?STRAND_CONFIG_SAME_ORDER_FLIPPED_STRAND:STRAND_CONFIG_OPPOSITE_ORDER_FLIPPED_STRAND;
						} else if (flipped[0].equals(referenceAlleles[1])) {
							return index == 1?STRAND_CONFIG_SAME_ORDER_FLIPPED_STRAND:STRAND_CONFIG_OPPOSITE_ORDER_FLIPPED_STRAND;
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

	// this routine currently assumes that the strand is constant across files
	public static void compareResults(String[] filenames, Logger log) {
		String trav1, trav2;
		Hashtable<String, String> hash1, hash2;
		int count;
		long time;
		int[][][] indices;
		String[] keys;
		int[][] agreements;
		int[] matches;
		Hashtable<String,String> missHash1, missHash2;
		
		if (filenames.length != 2) {
			log.reportError("Error - must provide two and only two filenames for comparison");
			return;
		}

		indices = new int[filenames.length][2][];
		for (int i = 0; i < filenames.length; i++) {
			indices[i][0] = ext.indexFactors(new String[][] {{"MarkerName"}}, Files.getHeaderOfFile(filenames[i], log), false, true, true, log, true);
			indices[i][1] = ext.indexFactors(new String[][] {{"Direction"}}, Files.getHeaderOfFile(filenames[i], log), false, true, true, log, true);
		}

		time = new Date().getTime();
		log.report(ext.getTime()+"\tLoading "+filenames[0], false, true);
		hash1 = HashVec.loadFileToHashString(filenames[0], indices[0][0], indices[0][1], false, "\t", true, false, false);
		log.report("Finished loading "+hash1.size()+" markers in " + ext.getTimeElapsed(time));
		
		time = new Date().getTime();
		log.report(ext.getTime()+"\tLoading "+filenames[1], false, true);
		hash2 = HashVec.loadFileToHashString(filenames[1], indices[1][0], indices[1][1], false, "\t", true, false, false);
		log.report("Finished loading "+hash2.size()+" markers in " + ext.getTimeElapsed(time));
		
		log.report(ext.getTime()+"\tParsing markers");
		keys = HashVec.getKeys(hash1, false, false);
		log.report(ext.getTime()+"\tThe following keys from "+filenames[0]+" were not found in "+filenames[1]+":");
		count = 0;
		for (int i = 0; i < keys.length; i++) {
			if (!hash2.containsKey(keys[i])) {
				log.report(keys[i], true, false);
				hash1.remove(keys[i]);
				count++;
			}
		}		
		if (count == 0) {
			log.report("[none]");
		} else {
			log.report(count+" total");
		}
		log.report("");

		log.report(ext.getTime()+"\tParsing markers");
		keys = HashVec.getKeys(hash2, false, false);
		log.report(ext.getTime()+"\tThe following keys from "+filenames[1]+" were not found in "+filenames[0]+":");
		count = 0;
		for (int i = 0; i < keys.length; i++) {
			if (!hash1.containsKey(keys[i])) {
				log.report(keys[i]);
				hash2.remove(keys[i]);
				count++;
			}
		}		
		if (count == 0) {
			log.report("[none]");
		} else {
			log.report(count+" total");
		}
		log.report("");
		
		log.report(ext.getTime()+"\tParsing union of marker sets");
		keys = HashVec.getKeys(hash1, false, false);
		log.report(ext.getTime()+"\tDetermining differential counts for a combined "+keys.length+" markers");
		count = 0;
		agreements = null;
//		missings1 = new Vector<Hashtable<String,String>>();
//		missings2 = new Vector<Hashtable<String,String>>();
		for (int i = 0; i < keys.length; i++) {
			trav1 = hash1.get(keys[i]);
			trav2 = hash2.get(keys[i]);
			if (agreements == null) {
				agreements = new int[trav1.length()][trav2.length()];
//				for (int j = 0; j < trav1.length(); j++) {
//					missings1.add(new Hashtable<String, String>());
//				}
//				for (int k = 0; k < trav2.length(); k++) {
//					missings2.add(new Hashtable<String, String>());
//				}
			}
			for (int j = 0; j < trav1.length(); j++) {
				for (int k = 0; k < trav2.length(); k++) {
					agreements[j][k] += trav1.charAt(j) == trav2.charAt(k)?1:0;

//					if (trav2.charAt(k) == '?') {
//						missings2.elementAt(k).put(keys[i], "");
//					}
				}
//				if (trav1.charAt(j) == '?') {
//					missings1.elementAt(j).put(keys[i], "");
//				}
			}
		}
		
		matches = new int[agreements.length];
		log.report("F1_study\tF2_study\tnumAgree\tnumExclusiveToF1\tnumExclusiveToF2");
		for (int i = 0; i < agreements.length; i++) {
			matches[i] = Sort.quicksort(agreements[i], Sort.DESCENDING)[0];
			missHash1 = new Hashtable<String, String>();
			missHash2 = new Hashtable<String, String>();
			for (int j = 0; j < keys.length; j++) {
				trav1 = hash1.get(keys[j]);
				trav2 = hash2.get(keys[j]);

				if (trav1.charAt(i) != trav2.charAt(matches[i])) {
					if (trav2.charAt(matches[i]) == '?') {
						missHash2.put(keys[j], "");
					} else if (trav1.charAt(i) == '?') {
						missHash1.put(keys[j], "");
					} else {
						System.err.println("Mismatch: '"+trav1.charAt(i)+"' and '"+trav2.charAt(matches[i])+"' for marker "+keys[j]);
					}
				}
			}

			log.report("Study"+(i+1)+"\tStudy"+(matches[i]+1)+"\t"+agreements[i][matches[i]]+"\t"+missHash2.size()+"\t"+missHash1.size());
			if (missHash2.size() > 0) {
				Files.writeList(HashVec.getKeys(missHash2), "Study"+(i+1)+"_exclusiveTo.dat");
			}
			if (missHash1.size() > 0) {
				Files.writeList(HashVec.getKeys(missHash1), "Study"+(i+1)+"_missingFrom.dat");
			}
			
//			missHash1 = missings1.elementAt(i);
//			missHash2 = missings2.elementAt(matches[i]);
//			keys = HashVec.getKeys(missHash1);
//			for (int j = 0; j < keys.length; j++) {
//				if (missHash2.containsKey(keys[j])) {
//					missHash1.remove(keys[j]);
//					missHash2.remove(keys[j]);
//				}
//			}
//			log.report("Study"+(i+1)+"\tStudy"+(matches[i]+1)+"\t"+agreements[i][matches[i]]+"\t"+missHash2.size()+"\t"+missHash1.size());
//			Files.writeList(HashVec.getKeys(missHash2), "Study"+(i+1)+"_exclusiveTo.dat");
//			Files.writeList(HashVec.getKeys(missHash1), "Study"+(i+1)+"_missingFrom.dat");
		}

		log.report("");
		log.report(ext.getTime()+"\tdone!");
	}

	public static void countMissing(String filename, Logger log) {
		BufferedReader reader;
		PrintWriter writer, w2;
		String[] line;
		String trav;
		int present, absent;
		
		try {
			reader = Files.getAppropriateReader(filename);
			writer = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+"_counts.out"));
			w2 = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+"_lowCount.out"));
			writer.println("MarkerName\tnumPresent\tnumAbsent");
			line = reader.readLine().trim().split("[\\s]+");
			ext.checkHeader(line, new String[] {"MarkerName", "Allele1", "Allele2", "Effect", "StdErr", "P-value", "Direction"}, new int[] {0,1,2,3,4,5,6}, false, log, true);
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				trav = line[6];
				present = absent = 0;
				for (int i = 0; i < trav.length(); i++) {
					if (trav.charAt(i) == '+' || trav.charAt(i) == '-' || trav.charAt(i) == '0') {
						present++;
					} else if (trav.charAt(i) == '?') {
						absent++;
					} else {
						log.report("Error - don't know what to do with symbol: "+trav.charAt(i));
					}
				}
				writer.println(line[0]+"\t"+present+"\t"+absent);
				if (present < 3) {
					w2.println(line[0]+"\t"+present+"\t"+absent);
				}
			}
			reader.close();
			writer.close();
			w2.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + filename + "\" not found in current directory");
			return;
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + filename + "\"");
			return;
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
		String compResults = null;
		String countMissing = null;
		String[] unitOfAnalyses = Aliases.MARKER_NAMES;
//		boolean gcControlOn = true;
		double[] gcControl = null;
		
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

//		compResults = "META_ANALYSIS_beta_se_Final1.tbl,AllResults.txt";
//		countMissing = "D:/mega/Discovery_allResults_23AndMe_just10K.tbl";
		
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
		"   (2) unit of anlaysis (i.e. unit=gene (default=snp))\n"+
		"   (3) name of output root (i.e. out="+output+" (default))\n"+
		"   (4) stderr scheme instead of sample-size weighting (i.e. se="+se+" (default))\n"+
//		"   (5) use genomic control (i.e. gcControlOn="+gcControlOn+" (default))\n"+
		"   (5) use genomic control (i.e. gcControl=-9,-9,-9 (default, one value per file, -1 for OFF, -9 for ON, other values used as given))\n"+
		" OR\n"+
		"   (1) compare two sets of meta-analysis results (i.e. compResults=fileA1.tbl,fileB1.tbl (not the default))\n"+
		" OR\n"+
		"   (1) count the number of studies missing each variant (i.e. countMissing=fileA1.tbl (not the default))\n"+
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
			} else if (args[i].startsWith("gcControl=")) {
			    String[] ctrl = ext.parseStringArg(args[i], null).split(",");
			    gcControl = new double[ctrl.length];
			    for (int j = 0; j < ctrl.length; j++) {
			        gcControl[j] = Double.parseDouble(ctrl[j]);
			    }
			    numArgs--;
			} else if (args[i].startsWith("unit=")) {
				unitOfAnalyses = ext.parseStringArg(args[i], null).equalsIgnoreCase("gene")?Aliases.GENE_UNITS:Aliases.MARKER_NAMES;
				numArgs--;
			} else if (args[i].startsWith("se=")) {
				se = ext.parseBooleanArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("compResults=")) {
				compResults = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("countMissing=")) {
				countMissing = ext.parseStringArg(args[i], null);
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
			    String[] files = analyze.split(",");
			    if (gcControl == null) {
			        gcControl = Array.doubleArray(files.length, -9);
			    } else if (files.length != gcControl.length) {
			        System.err.println("ERROR - list of GC correction values must equal the length of the files list.  Values are {-1=OFF, -9=ON, other values used as given}");
			        System.exit(1);
			    }
				metaAnalyze("./", files, unitOfAnalyses, output, SE_ANALYSIS, null, gcControl, new Logger());
			} else if (compResults != null) {
				compareResults(compResults.split(","), new Logger("compareMetalResults.log"));
			} else if (countMissing != null) {
				countMissing(countMissing, new Logger("countMissing.log"));
			} else {
				convertPlinkResults(dir, results, test, method, freq, useSE, useMetalNomenclature, outfile, false);
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
