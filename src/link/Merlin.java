package link;

import java.io.*;
//import java.util.*;
import common.*;

public class Merlin {
	public static final int CHR_START = 1;
	public static final int CHR_STOP  = 22;
//	public static final int CHR_START = 16;
//	public static final int CHR_STOP  = 16;
	public static final int UNKNOWN_TRAIT = 0;
	public static final int BINARY_TRAIT = 1;
	public static final int QUANTITATIVE_TRAIT = 2;
	
	public static void prepAll() {
        for (int i = CHR_START; i<=CHR_STOP; i++) {
        	createMerlinFiles(i, UNKNOWN_TRAIT);
        }
	}

	public static void batchAllBinary(String qsub) {
//		PrintWriter writer;
//		String chrome;
//		
//		try {
//	        writer = new PrintWriter(new FileWriter("MerlinAll.bat"));
//	        for (int i = CHR_START; i<=CHR_STOP; i++) {
//	        	chrome = ext.chrome(i);
////		        writer.println("java -cp \"C:\\Documents and Settings\\npankrat\\My Documents\\jProjects\\park.jar\" -Xms1024M -Xmx1024M  link.Merlin chr="+i);
//		        writer.println("jcp link.Merlin chr="+i);
//		        writer.println("merlin -d chr"+chrome+".dat -p re_chrom"+chrome+".pre -m chr"+chrome+".map --npl --tabulate --step 5 --markerNames --information --ibd --prefix merlin-chr"+chrome+"");
//		        writer.println();
//            }
//	        writer.close();
//        } catch (Exception e) {
//	        System.err.println("Error writing to "+"MerlinAll.bat");
//	        e.printStackTrace();
//        }
		
		String commands = "echo \"Starting chr# at...\"\n"
			+"date\n"
			+"jcp link.Merlin chr=#\n"
			+"merlin -d chr##.dat -p re_chrom##.pre -m chr##.map --npl --tabulate --step 5 --markerNames --information --ibd --prefix merlin-chr##\n"
			+"echo \"Finished chr# at...\"\n"
			+"date\n";
		Files.batchIt("batch", 0, CHR_START, CHR_STOP, 4, commands);
	}

	public static void batchAllVC(String qsub, boolean blade) {
		String commands;
		boolean win = System.getProperty("os.name").startsWith("Windows");

		if (qsub != null) {
			if (blade) {
				commands = ""+
				"cd "+qsub+"\n"+
				"java -cp /home/bc2/pankratz/park.jar link.Merlin chr=#\n"+
				"merlin -d chr##.dat -p re_chrom##.pre -m chr##.map --vc --tabulate --step 5 --markerNames --information --prefix vc-chr## > vc-chr##.log";
				Files.qsubBlade(qsub+"_vc", "java", CHR_START, CHR_STOP, commands, 1, 10);
			} else {
				commands = ""+
				"cd "+qsub+"\n"+
				Files.JCP+"link.Merlin chr=#\n"+
				"/share/apps/bin/merlin -d chr##.dat -p re_chrom##.pre -m chr##.map --vc --tabulate --step 5 --markerNames --information --prefix vc-chr## > vc-chr##.log";
				Files.qsub(qsub+"_vc", CHR_START, CHR_STOP, commands);
			}
		} else {
			commands = "echo \"Starting chr# at...\"\n"
				+(win?"date /t\ntime /t\n":"date\n")
				+"java -cp "+(win?"C:":"")+"/home/npankrat/park.jar link.Merlin chr=#\n"
				+"merlin -d chr##.dat -p re_chrom##.pre -m chr##.map --vc --tabulate --step 5 --markerNames --information --prefix vc-chr## > vc-chr##.log\n"
				+"echo \"Finished chr# at...\"\n"
				+(win?"date /t\ntime /t\n":"date\n");
			Files.batchIt("vc", 0, 8, CHR_STOP, 5, commands);
		}
	}

	public static void batchAllQuant(double[] quant, String qsub, boolean blade) {
		String commands;
		boolean win = System.getProperty("os.name").startsWith("Windows");
	
		if (qsub != null) {
			if (blade) {
				commands = "cd "+qsub+"\n"+
				"java -cp /home/bc2/pankratz/park.jar link.Merlin chr=#\n"+
				"merlin-regress -d chr##.dat -p re_chrom##.pre -m chr##.map --mean "+quant[0]+" --var "+quant[1]+" --her "+quant[2]+" --tabulate --step 5 --prefix regress-chr##";
				Files.qsubBlade(qsub+"_regress", null, CHR_START, CHR_STOP, commands, 1, 48);
			} else {
				commands = "cd "+qsub+"\n"+
				Files.JCP+"link.Merlin chr=#\n"+
				"/share/apps/bin/merlin-regress -d chr##.dat -p re_chrom##.pre -m chr##.map --mean "+quant[0]+" --var "+quant[1]+" --her "+quant[2]+" --tabulate --step 5 --prefix regress-chr##";
				Files.qsub(qsub+"_regress", CHR_START, CHR_STOP, commands);
			}
		} else {
			commands = "echo \"Starting chr# at...\"\n"
				+(win?"date /t\ntime /t\n":"date\n")
				+"java -cp "+(win?"C:":"")+"/home/npankrat/park.jar link.Merlin chr=#\n"
				+"merlin-regress -d chr##.dat -p re_chrom##.pre -m chr##.map --mean "+quant[0]+" --var "+quant[1]+" --her "+quant[2]+" --tabulate --step 5 --information --ibd --prefix regress-chr##\n"
				+"echo \"Finished chr# at...\"\n"
				+(win?"date /t\ntime /t\n":"date\n");
			Files.batchIt("regress", 0, CHR_START, CHR_STOP, 4, commands);
		}
	}

	public static void createMerlinFiles(int chr, int trait) {
		createMerlinFiles(chr, "chr"+ext.chrome(chr), trait);
	}
	
	public static void createMerlinFiles(int chr, String root, int trait) {
		PrintWriter writer;
		LinkageMap map;
		String[] markerNames;
		double[] positions;
		double[][] alleleFreqs;
		
		map = new LinkageMap(chr);
		markerNames = map.getMarkerNames();
		positions = map.getCumulativePositions(false);
		alleleFreqs = map.getAlleleFreqs();
		
		try {
			writer = new PrintWriter(new FileWriter(root+".dat"));
			switch (trait) {
			case UNKNOWN_TRAIT:
				writer.println(checkTrait(chr)?"A affection_status":"T trait");
				break;
			case BINARY_TRAIT:
				writer.println("A affection_status");
				break;
			case QUANTITATIVE_TRAIT:
				writer.println("T trait");
				break;
			}
			for (int i = 0; i<markerNames.length; i++) {
				writer.println("M "+markerNames[i]);
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing "+root+".dat");
			e.printStackTrace();
		}

		try {
			writer = new PrintWriter(new FileWriter(root+".map"));
			writer.println("CHROMOSOME\tMARKER\tPOSITION");
			for (int i = 0; i<markerNames.length; i++) {
				writer.println(chr+"\t"+markerNames[i]+"\t"+ext.formDeci(positions[i], 10, false));
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing "+root+".map");
			e.printStackTrace();
		}
		
		try {
			writer = new PrintWriter(new FileWriter(root+".freq"));
			for (int i = 0; i<markerNames.length; i++) {
				writer.println("M "+markerNames[i]);
				writer.println("F "+Array.toStr(alleleFreqs[i], 6, 6, " "));
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing "+root+".freq");
			e.printStackTrace();
		}
	}
	
	public static boolean checkTrait(int chr) {
		BufferedReader reader;
		boolean binary;
		String[] line;
		int conf;

		binary = true;
		conf = 0;
		try {
	        reader = new BufferedReader(new FileReader("re_chrom"+ext.chrome(chr)+".pre"));
	        while (binary && conf < 10 && reader.ready()) {
	        	line = reader.readLine().trim().split("[\\s]+");
	        	if (line[5].equals("0")) {
	        		
	        	} else if (line[5].equals("1") || line[5].equals("2")) {
	        		conf++;
	        	} else {
	        		binary = false;
	        	}
	        }
	        reader.close();
        } catch (FileNotFoundException fnfe) {
	        System.err.println("Error: file \""+"re_chrom"+ext.chrome(chr)+".pre"+"\" not found in current directory");
	        System.exit(1);
        } catch (IOException ioe) {
	        System.err.println("Error reading file \""+"re_chrom"+ext.chrome(chr)+".pre"+"\"");
	        System.exit(2);
        }

        return binary;
	}

	public static void main(String[] args) {
	    int numArgs = args.length;
	    int chr = -1;
	    boolean batch = false;
	    String qsub = null;
	    boolean vc = false;
	    boolean prep = false;
	    int trait = UNKNOWN_TRAIT;
	    double[] quant = null;
	    boolean blade = false;
	    
//	    vc = true;
//	    qsub = "PCA1";
	    chr=1;
	    
	    // PCA1
//	    qsub = "PCA1";
//	    quant = new double[] {0.070835749, 3.30957236, 0.4079};

	    // PCA1 boxcox
//	    qsub = "PCA1_boxcox";
//	    quant = new double[] {3.209175862, 0.541843525, 0.4079};
	    
	    // PCA2
//	    quant = new double[] {0.010311256, 1.196156432, 0.2450};

	    // PCA3
//	    quant = new double[] {0.016727564, 0.587978393, 0.4974};
	    
	    String usage = "\n"+
	    "link.Merlin requires 0-1 arguments\n"+
	    "   (1) generate files for specific chromosome (i.e. chr=2 (not the default))\n"+
	    "   (2) binary trait (i.e. -binary)\n"+
	    "   (3) quantitative trait (i.e. -quant)\n"+
	    "  OR\n"+
	    "   (1) generate files for all chromsomes (i.e. -prep (not the default))\n"+
	    "  OR\n"+
	    "   (1) batch create and analyze a binary trait (i.e. -batch (not the default))\n"+
	    "  OR\n"+
	    "   (1) batch create and analyze a quantitative trait using vc (i.e. -vc (not the default))\n"+
	    "   (2) make .qsub files for this directory instead of batch files (i.e. qsub=dir (not the default; don't use slash))\n"+
	    "   (3) use blade header for qsub instead of alc (i.e. -blade (not the default))\n"+
	    "  OR\n"+
	    "   (1) batch create and analyze a quantitative trait using h-e (i.e. quant=mean,variance,heritability)\n"+
	    "";

	    for (int i = 0; i<args.length; i++) {
		    if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
			    System.err.println(usage);
			    System.exit(1);
		    } else if (args[i].startsWith("chr=")) {
			    chr = Integer.parseInt(args[i].split("=")[1]);
			    numArgs--;
		    } else if (args[i].startsWith("-batch")) {
			    batch = true;
		    	trait = BINARY_TRAIT;
			    numArgs--;
		    } else if (args[i].startsWith("-prep")) {
			    prep = true;
			    numArgs--;
		    } else if (args[i].startsWith("qsub=")) {
			    qsub = args[i].split("=")[1];
			    numArgs--;
		    } else if (args[i].startsWith("-binary")) {
		    	trait = BINARY_TRAIT;
			    numArgs--;
		    } else if (args[i].startsWith("-vc")) {
			    vc = true;
			    numArgs--;
		    } else if (args[i].startsWith("-quant")) {
		    	trait = QUANTITATIVE_TRAIT;
			    numArgs--;
		    } else if (args[i].startsWith("quant=")) {
			    quant = Array.toDoubleArray(args[i].split("=")[1].split(","));
				if (quant.length == 3) {
				    numArgs--;
				} else {
					System.err.println("Error - batchAllQuant requires three values");
				}
		    } else if (args[i].startsWith("-blade")) {
			    blade = true;
			    numArgs--;
		    } else {
			    System.err.println("Error - don't know what to do with argument: "+args[i]);
		    }
	    }
	    if (numArgs!=0) {
		    System.err.println(usage);
		    System.exit(1);
	    }

	    try {
	    	if (prep) {
	    		prepAll();
	    	} else if (batch) {
	    		batchAllBinary(qsub);
	    	} else if (vc) {
	    		batchAllVC(qsub, blade);
	    	} else if (quant != null) {
	    		batchAllQuant(quant, qsub, blade);
	    	} else if (chr > 0) {
	    		createMerlinFiles(chr, trait);
	    	}
	    } catch (Exception e) {
		    e.printStackTrace();
	    }
    }
}
