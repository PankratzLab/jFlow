package org.genvisis.bioinformatics;

import java.io.*;
import java.util.*;

import org.genvisis.common.*;

import com.google.common.primitives.Chars;

public class SequenceVariants {
	public static final String[] ALLELE_SPECIFIC_FREQS = {"chr", "pos", "ref", "control_A_freq", "case_A_freq", "control_C_freq", "case_C_freq", "control_G_freq", "case_G_freq", "control_T_freq", "case_T_freq"};
	public static final double MIN_FREQ = 0.01;
	public static final double FREQ_BUFFER = 0.001;
	
	private static void parse(String filename, String annotationFile, double penetrance) {
		BufferedReader reader;
		PrintWriter writer, writer2;
		String[] line;
		String trav;
		Hashtable<String, Vector<String>> hash, dbsnpHash;
		double[][] freqs;
		int[] order;
		int refIndex;
		double altCaseFreq, altControlFreq, refCaseFreq, refControlFreq, mafCase, mafControl;
		
		if (annotationFile != null && Files.exists(annotationFile)) {
			hash = HashVec.loadFileToHashVec(annotationFile, new int[] {1,2,3,4}, new int[] {8}, "\t", true, true);
			dbsnpHash = HashVec.loadFileToHashVec(annotationFile, new int[] {1,2,3,4}, new int[] {0}, "\t", true, true);
		} else {
			hash = new Hashtable<String, Vector<String>>();
			dbsnpHash = new Hashtable<String, Vector<String>>();
		}
		
		try {
			reader = new BufferedReader(new FileReader(filename));
			writer = new PrintWriter(new FileWriter(filename+"_parsed.xln"));
			writer2 = new PrintWriter(new FileWriter(filename+"_SeattleSeq.input"));
			writer.println("Chr\tPosition\tRef\tAlt\tRefIsMostCommonAllele\tRefCasesFreq\tRefControlsFreq\tAltCasesFreq\tAltControlsFreq\tCallrateCases\tCallrateControls\tAdjCaseMAF\tAdjControlMAF\tOR\tEstFreq");
			ext.checkHeader(reader.readLine().trim().split("[\\s]+"), ALLELE_SPECIFIC_FREQS, true);
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				freqs = new double[3][4];
				freqs[0][0] = Double.parseDouble(line[4]); // case_A_freq
				freqs[1][0] = Double.parseDouble(line[3]); // control_A_freq
				freqs[2][0] = Array.sum(freqs[0]);
				freqs[0][1] = Double.parseDouble(line[6]); // case_C_freq
				freqs[1][1] = Double.parseDouble(line[5]); // control_C_freq
				freqs[2][1] = Array.sum(freqs[0]);
				freqs[0][2] = Double.parseDouble(line[8]); // case_G_freq
				freqs[1][2] = Double.parseDouble(line[7]); // control_G_freq
				freqs[2][2] = Array.sum(freqs[0]);
				freqs[0][3] = Double.parseDouble(line[10]); // case_T_freq
				freqs[1][3] = Double.parseDouble(line[9]); // control_T_freq
				freqs[2][3] = Array.sum(freqs[0]);
				
				refIndex = Chars.indexOf(Sequence.ALLELES, line[2].charAt(0));
				order = Sort.quicksort(freqs[2]);
				for (int i = 0; i < order.length; i++) {
					if (order[i] != refIndex && (freqs[0][order[i]] > MIN_FREQ || freqs[1][order[i]] > MIN_FREQ)) {
						refCaseFreq = freqs[0][refIndex];
						refControlFreq = freqs[1][refIndex];
						altCaseFreq = freqs[0][order[i]];
						altControlFreq = freqs[1][order[i]];
						
						mafCase = altCaseFreq/(refCaseFreq+altCaseFreq);
						mafControl = altControlFreq/(refControlFreq+altControlFreq);
						if (mafControl > 0.5) {
							mafCase = 1-mafCase;
							mafControl = 1-mafControl;
						}						
						
						trav = line[0].substring(3)+"\t"+line[1]+"\t"+line[2]+"\t"+Sequence.ALLELES[order[i]];
						writer.println(trav+"\t"+(order[0]==refIndex?"1":"0")+
								"\t"+refCaseFreq+"\t"+refControlFreq+"\t"+altCaseFreq+"\t"+altControlFreq+
								"\t"+(refCaseFreq+altCaseFreq)+"\t"+(refControlFreq+altControlFreq)+"\t"+mafCase+"\t"+mafControl+
								"\t"+ext.formDeci((mafCase+FREQ_BUFFER)/(mafControl+FREQ_BUFFER),2)+"\t"+ext.prettyP(mafCase*penetrance+mafControl*(1-penetrance))+
								"\t"+(dbsnpHash.containsKey(trav)?dbsnpHash.get(trav):".")+"\t"+(hash.containsKey(trav)?hash.get(trav):"."));

						writer2.println(trav);
					}
					
				}

			}
			writer.close();
			writer2.close();
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + filename + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + filename + "\"");
			System.exit(2);
		}
		
		
		
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "SequenceVariants.dat";

		String usage = "\n" + "bioinformatics.SequenceVariants requires 0-1 arguments\n" + "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = args[i].split("=")[1];
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
			parse("D:/Myron/Indian_Diabetes/SequencingPilot/alleleSpecificFrequencies.txt", "D:/Myron/Indian_Diabetes/SequencingPilot/SeattleSeqAnnotation.txt", 0.13);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
