package bioinformatics;

import java.io.*;
import java.util.*;
import common.*;

public class SnpNexusParser {
	public static final String[] MARKS = {
		"Genomic Coordinates and External Links",
		"Hapmap CEU population",
		"Hapmap YRI population",
		"Hapmap JPT population",
		"Hapmap CHB population",
		"Consequences on Refseq",
		"Consequences on Ensembl",
		"Consequences on VEGA",
		"Consequences on AceView",
		"Consequences on UCSC",
		"Conserved Transcription Factor Binding Sites (TFBS)",
		"First-Exon and Promoter Prediction (FirstEF)",
		"miRBASE 11.0",
		"TargetScan miRNA Regulatory Sites",
		"microRNAs (miRNA Registry) / snoRNAs and scaRNAs (snoRNA-LBME-DB)",
		"Genetic Association Studies of Complex Diseases and Disorders (GAD)"
		};

	public static final String[][] HEADERS = {
		{"SNP_name", "Contig", "contigStart", "contigEnd", "Chromosome", "chromStart", "chromEnd", "dbSNP", "HapMap"},
		{"Name", "Genotype1", "Count", "Frequency", "Genotype2", "Count", "Frequency", "Genotype3", "Count", "Frequency", "Allele1", "Count", "Frequency", "Allele2", "Count", "Frequency"},
		{"Name", "Genotype1", "Count", "Frequency", "Genotype2", "Count", "Frequency", "Genotype3", "Count", "Frequency", "Allele1", "Count", "Frequency", "Allele2", "Count", "Frequency"},
		{"Name", "Genotype1", "Count", "Frequency", "Genotype2", "Count", "Frequency", "Genotype3", "Count", "Frequency", "Allele1", "Count", "Frequency", "Allele2", "Count", "Frequency"},
		{"Name", "Genotype1", "Count", "Frequency", "Genotype2", "Count", "Frequency", "Genotype3", "Count", "Frequency", "Allele1", "Count", "Frequency", "Allele2", "Count", "Frequency"},
		{"SNP_name", "Refseq_gene", "Refseq_transcript", "EntrezGene_ID", "SNP_Predicted_function", "cdna_position", "cds_position", "aa_position", "peptide_shift", "distanceTOsplice"},
		{"SNP_name", "Ensembl_gene", "Ensembl_transcript", "Symbol", "SNP_Predicted_function", "cdna_position", "cds_position", "aa_position", "peptide_shift", "distanceTOsplice"},
		{"SNP_name", "VEGA_gene", "VEGA_transcript", "Symbol", "SNP_Predicted_function", "cdna_position", "cds_position", "aa_position", "peptide_shift", "distanceTOsplice"},
		{"SNP_name", "AceView_gene", "SNP_Predicted_function", "cdna_position", "cds_position", "aa_position", "peptide_shift", "distanceTOsplice"},
		{"SNP_name", "UCSC_known_gene", "EntrezGene_ID", "SNP_Predicted_function", "cdna_position", "cds_position", "aa_position", "peptide_shift", "distanceTOsplice"},
		{"SNP_name", "TFBS_id", "chromStart", "chromEnd", "TFBS_Accession", "TFBS_Species", "TFBS_name", "SwissProt_Accession"},
		{"SNP_name", "chrom", "chromStart", "chromEnd", "Name", "Score", "Strand"},
		{"SNP_name", "chrom", "chromStart", "chromEnd", "Name", "Score", "Strand", "Type"},
		{"SNP_name", "chrom", "chromStart", "chromEnd", "Name", "Score", "Strand"},
		{"SNP_name", "chrom", "chromStart", "chromEnd", "Name", "Score", "Strand", "Type"},
		{"GAD_ID", "Association(Y/N)", "Phenotype", "Disease Class", "Gene", "Reference", "Pubmed", "Population", "EntrezGene_ID"}
		};
	
	public static final boolean[] HEADER_ON_SAME_LINE = {
		false,
		true,
		true,
		true,
		true,
		false,
		false,
		false,
		false,
		false,
		false,
		false,
		false,
		false,
		false,
		false,
	};

	public static final int[] KEY_TYPE = {
		0,
		1,
		1,
		1,
		1,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		2,
	};
	
	public static void parse(String dir, String input, String output) {
		BufferedReader reader;
        PrintWriter writer;
        String[] line;
        Hashtable<String,Hashtable<String,String[]>> hashes = new Hashtable<String,Hashtable<String,String[]>>();
        Hashtable<String,String[]> hash;
        Vector<String> v = new Vector<String>();
        int index, option;
        String[] data;
        double d;
        String[] keys;
        String snp, rsNum;
        
        for (int i = 0; i<MARKS.length; i++) {
        	hashes.put(MARKS[i], new Hashtable<String,String[]>());
        }
        try {
	        reader = new BufferedReader(new FileReader(dir+input));
	        option = -1;
	        while (reader.ready()) {
	        	line = reader.readLine().split("\\t", -1);
	        	index = ext.indexOfStartsWith(line[0], MARKS, true);
	        	if (index >= 0) {
	        		option = index;
	        		if (HEADER_ON_SAME_LINE[option]) {
	        			line = Array.subArray(line, 1);
	        		} else {
	        			line = reader.readLine().split("\\t", -1);
	        		}
	        		ext.checkHeader(line, HEADERS[option], true);
	        		System.out.println("Header checks out for "+MARKS[option]);
        			line = reader.readLine().split("\\t", -1);
	        	}
	        	hash = hashes.get(MARKS[option]);
	        	switch (option) {
	        	case -1:
	        		System.err.println("Error - could not find a valid option");
	        		System.exit(1);
	        		break;
	        	case 0:
	        		for (int i = 0; i<2; i++) {
	        			v = new Vector<String>();
		        		if (hash.containsKey(line[0])) {
		        			data = hash.get(line[0]);
		        			for (int j = 0; j<data.length; j++) {
		        				v.add(data[j]);
                            }
		        		}
		        		if (!line[7+i].equals("")) {
		        			HashVec.addIfAbsent(line[7+i], v);
		        		}
		        		hash.put(line[0], Array.toStringArray(v));
                    }
	        		break;
	        	case 1:
	        	case 2:
	        	case 3:
	        	case 4:
	        		d = Double.parseDouble(line[12]);
	        		hash.put(line[0], new String[] {d > 0.50?line[15]:line[12]});
	        		break;
	        	default:
	        		break;
	        	}
	        }
	        reader.close();
        } catch (FileNotFoundException fnfe) {
	        System.err.println("Error: file \""+dir+input+"\" not found in current directory");
	        System.exit(1);
        } catch (IOException ioe) {
	        System.err.println("Error reading file \""+dir+input+"\"");
	        System.exit(2);
        }
        
        try {
	        writer = new PrintWriter(new FileWriter(dir+output));
	        writer.print("SNP"+"\t"+"rsNumber");
        	for (int j = 1; j<MARKS.length; j++) {
    	        hash = hashes.get(MARKS[j]);
        		if (hash.size() > 0) {
        			writer.print("\t"+MARKS[j]);
        		}
        	}
        	writer.println();
	        hash = hashes.get(MARKS[0]);
	        keys = HashVec.getKeys(hash);
	        for (int i = 0; i<keys.length; i++) {
	        	snp = keys[i];
	        	line = hashes.get(MARKS[0]).get(snp);
	        	for (int j = 0; j<line.length; j++) {
	        		rsNum = line[j];
		        	writer.print(snp+"\t"+rsNum);		        	
		        	for (int k = 1; k<MARKS.length; k++) {
		    	        hash = hashes.get(MARKS[k]);
		        		if (hash.size() > 0) {
		        			switch (KEY_TYPE[k]) {
		        			case 0:
		        				if (hash.containsKey(snp)) {
		        					writer.print("\t"+hash.get(snp)[0]);
		        				} else {
		        					writer.print("\t.");
		        				}
		        				break;
		        			case 1:
		        				if (hash.containsKey(rsNum)) {
		        					writer.print("\t"+hash.get(rsNum)[0]);
		        				} else {
		        					writer.print("\t.");
		        				}
		        				break;
		        			default:
		        				System.err.println("Error - no lookup key defined for '"+MARKS[k]+"'");
		        				break;
		        			}
		        		}
	                }
		        	writer.println();
                }
            }
	        writer.close();
        } catch (Exception e) {
	        System.err.println("Error writing to "+dir+output);
	        e.printStackTrace();
        }
	}
	
	public static void main(String[] args) {
	    int numArgs = args.length;
	    String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\tWork\\SequencingProjectWithCIDR\\TestInDeepTrio\\";
//	    String filename = "results180.txt";
	    String filename = "resultsAll.txt";

	    
	    String usage = "\n"+
	    "bioinformatics.SnpNexusParser requires 0-1 arguments\n"+
	    "   (1) filename (i.e. file="+filename+" (default))\n"+
	    "";

	    for (int i = 0; i<args.length; i++) {
		    if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
			    System.err.println(usage);
			    System.exit(1);
		    } else if (args[i].startsWith("file=")) {
			    filename = args[i].split("=")[1];
			    numArgs--;
		    }
	    }
	    if (numArgs!=0) {
		    System.err.println(usage);
		    System.exit(1);
	    }
	    try {
		    parse(dir, filename, filename+"_proc.xln");
	    } catch (Exception e) {
		    e.printStackTrace();
	    }
    }
}
