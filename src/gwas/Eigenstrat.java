// at one point, I'd thought this required pre-normalization, as is necessary when comparison PCAs from SPSS; but that was not the case with Eigenstrat
// turned out there was just a bug in the Eigensoft package itself when markers were listed out of order (it sorts them, but the weights get off)
// to avoid this problem, make sure to perform a --make-bed after any merging (i.e. with HapMap data)
// after fancy MAF weighting, the final score is divided by the square root of the sample size

// probably want to start using a Logger at some point

package gwas;

import java.io.*;
import java.util.*;

import common.*;
import filesys.*;

public class Eigenstrat {
	public static boolean DEFAULT_FANCY_WEIGHTING = true;
	
	public static void generateFiles(String sourceRoot) {
        PrintWriter writer;
        FamilyStructure struct;
        byte[] affs;
        
        if (sourceRoot.equals("")) {
        	System.err.println("Error - no root of plink filenames provided");
        	System.exit(1);
        }
        
        struct = new FamilyStructure(sourceRoot+".fam");
        affs = struct.getAffections();
        for (int i = 0; i<affs.length; i++) {
        	if (affs[i] == -9) {
        		affs[i] = 1;
        	}
        }
        struct.setAffections(affs);
        struct.writeToFile(sourceRoot+".fam", false);
        
		try {
	        writer = new PrintWriter(new FileWriter("convertf.par"));
	        writer.println("genotypename:    "+sourceRoot+".bed");
	        writer.println("snpname:         "+sourceRoot+".bim");
	        writer.println("indivname:       "+sourceRoot+".fam");
	        writer.println("outputformat:    EIGENSTRAT");
	        writer.println("genotypeoutname: "+sourceRoot+".eigenstratgeno");
	        writer.println("snpoutname:      "+sourceRoot+".snp");
	        writer.println("indivoutname:    "+sourceRoot+".ind");
	        writer.println("familynames:     NO");
	        writer.close();
        } catch (Exception e) {
	        System.err.println("Error writing to "+"convertf.par");
	        e.printStackTrace();
        }

        try {
	        writer = new PrintWriter(new FileWriter("smartpca.par"));
	        writer.println("genotypename:     "+sourceRoot+".eigenstratgeno");
	        writer.println("snpname:          "+sourceRoot+".snp");
	        writer.println("indivname:        "+sourceRoot+".ind");
	        writer.println("evecoutname:      "+sourceRoot+".pca.evec");
	        writer.println("evaloutname:      "+sourceRoot+".eval");
	        writer.println("numoutevec:       20");
	        writer.println("numoutlieriter:   0");
	        writer.println("snpweightoutname: "+sourceRoot+".weights.out");
	        writer.close();
        } catch (Exception e) {
	        System.err.println("Error writing to "+"smartpca.par");
	        e.printStackTrace();
        }
        
        try {
	        writer = new PrintWriter(new FileWriter("master"));
	        writer.println("convertf -p convertf.par");
	        writer.println("smartpca -p smartpca.par");
	        writer.close();
	        Files.chmod("master");
        } catch (Exception e) {
	        System.err.println("Error writing to "+"master");
	        e.printStackTrace();
        }
        System.out.println("./master");
        System.out.println("    will run...");
        System.out.println("convertf -p convertf.par");
        System.out.println("smartpca -p smartpca.par");
	}

	public static void parse(String sourceRoot, String targetRoot, boolean eigenFormat, boolean fancy_weighting, boolean suppress) {
        BufferedReader reader;
        PrintWriter writer;
        String[] line;
        double[][] weights;
        String[] markerNames;
        double[] scores;
        int numEigens;
        double alleleCount;
        char[][] alleles;
        Hashtable<String, String> strandHash, mafHash, weightsHash;
        String str, maf;
        double[] mafs, fancyweights;
        Vector<String> ids;
        Vector<double[]> allScores;
        double[][] scoreMatrix = null;
        double[] data, array, means, stdevs;
        int count;
        double sampleSizeCorrection;
        
        System.out.println("Loading allele frequencies to use for missing genotypes...");
        mafHash = HashVec.loadFileToHashString(sourceRoot+".frq", 1, new int[] {2,3,4}, " ", false);
        System.out.println("Indexing alleles...");
        strandHash = HashVec.loadFileToHashString(sourceRoot+".snp", 0, new int[] {4,5}, " ", false);
        
        System.out.print("Parsing weights...");
        numEigens = Files.getHeaderOfFile(sourceRoot+".weights.out", "[\\s]+", new Logger(null)).length - 3;
        weightsHash = HashVec.loadFileToHashString(sourceRoot+".weights.out", 0, Array.subArray(Array.intArray(numEigens+3), 3, numEigens+3), "\t", false);
        System.out.println("found weights for "+weightsHash.size()+" markers");

        if (!new File(targetRoot+(eigenFormat?".bim":".map")).exists()) {
        	System.err.println("Error - could not find a '"+targetRoot+(eigenFormat?".bim":".map")+"' to pair up with '"+targetRoot+(eigenFormat?".eigenstratgeno":".ped")+"'; problem if markers are in a different order ("+targetRoot+".snp does not necessarily reflect the actual order in "+targetRoot+".eigenstratgeno)");
        	System.exit(1);
        }
    	markerNames = new SnpMarkerSet(targetRoot+(eigenFormat?".bim":".map"), false, new Logger(null)).getMarkerNames();
        weights = new double[markerNames.length][];
        alleles = new char[markerNames.length][2];
        mafs = new double[markerNames.length];
        fancyweights = new double[markerNames.length];

        count = 0;
        for (int i = 0; i < markerNames.length; i++) {
        	if (weightsHash.containsKey(markerNames[i])) {
        		weights[i] = Array.toDoubleArray(weightsHash.remove(markerNames[i]).split("[\\s]+"));

        		str = strandHash.get(markerNames[i]);
	        	if (str == null) {
	        		System.err.println("Error - '"+markerNames[i]+"' has a weight but is not present in the source .snp file");
	        		System.exit(1);
	        	}
	        	line = str.split("[\\s]+");
	        	alleles[i][0] = line[0].charAt(0);
	        	alleles[i][1] = line[1].charAt(0);
	        	
	        	if (fancy_weighting) {
	        		maf = mafHash.get(markerNames[i]);
	        		line = maf.split("[\\s]+");
		        	mafs[i] = Double.parseDouble(line[2]);
	        		if (!maf.startsWith(str)) {
	        			if (str.equals(line[1]+" "+line[0])) {
	        				mafs[i] = 1 - mafs[i];
	        			} else {
	        				System.err.println("Error - different alleles for "+markerNames[i]+" in "+sourceRoot+".weights.out and "+targetRoot+".frq");
	        			}
	        		}
	        		if (mafs[i] <0.00001 || mafs[i] > 0.99999) {
	        			mafs[i] = 0.00001;
	        		}	        		
		        	fancyweights[i] = Math.sqrt(mafs[i]*(1-mafs[i]));
	        	} else {
	        		mafs[i] = 0;
	        		fancyweights[i] = 1;
	        	}
        	} else {
        		weights[i] = null;
        		count++;
        	}
		}
        if (weightsHash.size() > 0) {
        	System.err.println("Error - the following markers were missing from the target files (i.e. "+targetRoot+(eigenFormat?".bim":".map")+"):");
        	line = HashVec.getKeys(weightsHash);
        	for (int j = 0; j < line.length; j++) {
        		System.err.println("        "+line[j]);
			}
        	System.exit(1);
        }
        if (count > 0) {
        	System.err.println("FYI - there were "+count+" markers in the target map ("+targetRoot+(eigenFormat?".bim":".map")+") that were not included in the eigenstrat PCA");
        }

        System.out.println("Computing eigenvalues...");
        ids = new Vector<String>();
        allScores = new Vector<double[]>();
        try {
	        if (eigenFormat) {
    	        reader = new BufferedReader(new FileReader(targetRoot+".eigenstratgeno"));
    	        if (new File(targetRoot+".fam").exists()) {
    	        	System.out.println("Lately, we've had trouble with the eigen .ind file not having both FID and IID; using plink-style .fam instead");
    	        	ids = HashVec.loadFileToVec(targetRoot+".fam", false, new int[] {0,1}, false, false);
        	        if (!Array.equals(Matrix.extractColumn(Array.splitStrings(Array.toStringArray(ids), true), 1), HashVec.loadFileToStringArray(targetRoot+".ind", false, new int[] {0}, false), false)) {
        	        	System.err.println("But unfortunately, "+targetRoot+".fam"+" doesn't match up with "+targetRoot+".ind; duplicating IID instead for now");
            	        ids = HashVec.loadFileToVec(targetRoot+".ind", false, new int[] {0, 0}, false, false);
        	        }
    	        }
	        	scoreMatrix = new double[ids.size()][numEigens];
	        	for (int i = 0; i<markerNames.length; i++) {
    	        	str = reader.readLine();
    	        	data = new double[ids.size()];
    	        	for (int j = 0; j<ids.size(); j++) {
    	        		if (str.charAt(j) == '9') {
    	        			data[j] = mafs[i]*2;
//	        				data[j] = Double.NaN;
    	        		} else {
	    	        		data[j] = Double.parseDouble(str.charAt(j)+"");
    	        		}
    	        	}
    	        	
//    	        	mafs[i] = Array.mean(data, true)/2;
//    	        	for (int j = 0; j<data.length; j++) {
//    	        		if (Double.isNaN(data[j])) {
//	    	        		data[j] = mafs[i]*2;
//    	        		}
//    	        	}
//	        		if (mafs[i] <0.00001 || mafs[i] > 0.99999) {
//	        			mafs[i] = 0.00001;
//	        		}	        		
//		        	fancyweights[i] = Math.sqrt(mafs[i]*(1-mafs[i]));
    	        	
    	        	for (int j = 0; j<ids.size(); j++) {
	        			for (int k = 0; k<numEigens; k++) {
	        				scoreMatrix[j][k] += weights[i][k]*(data[j]-mafs[i])/fancyweights[i];
	                    }
	        			
                    }
    	        }
    	        reader.close();
        	} else {
        		reader = new BufferedReader(new FileReader(targetRoot+".ped"));
    	        while (reader.ready()) {
    	        	str = reader.readLine();
    	        	line = str.trim().split("[\\s]+");
    	        	ids.add(line[0]+"\t"+line[1]);

    	        	scores = new double[numEigens];
    	        	for (int i = 0; i<markerNames.length; i++) {
    	        		alleleCount = computeAlleleCount(new char[] {line[6+i*2+0].charAt(0), line[6+i*2+1].charAt(0)}, alleles[i]);
    	        		if (weights[i] == null) {
    	        			// skip
    	        		} else if (alleleCount == -9) {
    	        			System.err.println("Error - alleles don't match up for "+markerNames[i]);
    	        		} else {
    	        			if (alleleCount == -1) {
    	        				alleleCount = mafs[i]*2;
    	        			}
    	        			for (int j = 0; j<numEigens; j++) {
    		        			scores[j] += weights[i][j]*(alleleCount-mafs[i])/fancyweights[i];
    	                    }
    	        		}
                    }
    	        	allScores.add(scores);
    	        }
    	        reader.close();
    	        scoreMatrix = Matrix.toDoubleArrays(allScores);        
        	}	        
	        
        } catch (FileNotFoundException fnfe) {
	        System.err.println("Error: file \""+targetRoot+".ped"+"\" not found in current directory");
	        System.exit(1);
        } catch (IOException ioe) {
	        System.err.println("Error reading file \""+targetRoot+".ped"+"\"");
	        System.exit(2);
        }
        
        System.out.println("Writing raw scores to file...");
        try {
	        writer = new PrintWriter(new FileWriter(targetRoot+(fancy_weighting?"_fancy":"")+"_eigens.xln"));
	        writer.println("FID\tIID\t"+Array.toStr(Array.stringArraySequence(numEigens, "C")));
	        for (int i = 0; i<ids.size(); i++) {
	        	writer.print(ids.elementAt(i));
        		for (int j = 0; j<numEigens; j++) {
        			writer.print("\t"+ext.formDeci(scoreMatrix[i][j], 4, true));
                }
        		writer.println();
            }
        writer.close();
        } catch (Exception e) {
	        System.err.println("Error writing to "+targetRoot+"_eigens.xln");
	        e.printStackTrace();
        }

        System.out.println("Writing normalized scores to file...");
        means = new double[numEigens];
        stdevs = new double[numEigens];
        for (int i = 0; i<numEigens; i++) {
        	array = Matrix.extractColumn(scoreMatrix, i);
        	means[i] = Array.mean(array);
        	stdevs[i] = Array.stdev(array);
        }
        try {
	        writer = new PrintWriter(new FileWriter(targetRoot+(fancy_weighting?"_fancy":"")+"_postnormed_eigens.xln"));
	        writer.println("FID\tIID\t"+Array.toStr(Array.stringArraySequence(numEigens, "C")));
	        sampleSizeCorrection = Math.sqrt(ids.size());
	        for (int i = 0; i<ids.size(); i++) {
	        	writer.print(ids.elementAt(i));
        		for (int j = 0; j<numEigens; j++) {
        			writer.print("\t"+ext.formDeci((scoreMatrix[i][j]-means[j])/stdevs[j]/sampleSizeCorrection, 4, true));
                }
        		writer.println();
            }
	        writer.close();
        } catch (Exception e) {
	        System.err.println("Error writing to "+targetRoot+"_eigens.xln");
	        e.printStackTrace();
        }
	}
	
	public static double computeAlleleCount(char[] observed, char[] possible) {
		int count;
		
		if (observed[0] == '0') {
			return -1;
		} else {
			count = 0;
			for (int i = 0; i<2; i++) {
				if (observed[i] == possible[0]) {
					count++;
				} else if (observed[i] != possible[1]) {
					return -9;
				}
            }
			return count;
		}
	}
	
	public static void convertToMDS(String filename) {
		FamilyStructure struct;
		BufferedReader reader;
        PrintWriter writer;
        String[] line;
        int count;
        String dir;
        String[][] ids;
        
        dir = ext.parseDirectoryOfFile(filename);        
        struct = new FamilyStructure(dir+ext.rootOf(ext.rootOf(filename))+".fam");
        ids = struct.getIds();
        try {
	        reader = new BufferedReader(new FileReader(filename));
	        writer = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+".mds"));
	        line = reader.readLine().trim().split("[\\s]+");
	        if (!line[0].equals("#eigvals:")) {
	        	System.err.println("Error - I do not think this file contains what you think it contains");
	        	System.exit(1);
	        }
	        writer.print("FID\tIID");
	        for (int i = 1; i<line.length; i++) {
	        	writer.print("\tPC"+i+"_"+line[i]);
            }
	        writer.println();
	        count = 0;
	        while (reader.ready()) {
	        	line = reader.readLine().trim().split("[\\s]+");
	        	if (!line[0].equals(ids[count][1])) {
	        		System.err.println("Error - mismatched indiviudals on line "+(count+1)+" (found "+line[0]+"; expecting "+ids[count][1]+"), can't interpolate Family IDs");
	        		System.exit(1);
	        	}
	        	writer.println(ids[count][0]+"\t"+ids[count][1]+"\t"+Array.toStr(Array.subArray(line, 1)));
	        	count++;
	        }
            writer.close();
	        reader.close();
        } catch (FileNotFoundException fnfe) {
	        System.err.println("Error: file \""+filename+"\" not found in current directory");
	        System.exit(1);
        } catch (IOException ioe) {
	        System.err.println("Error reading file \""+filename+"\"");
	        System.exit(2);
        }
	}
	
	public static void main(String[] args) {
	    int numArgs = args.length;
	    String sourceRoot = "combo";
	    String targetRoot = "plink";
	    boolean create = false;
	    boolean parse = false;
	    boolean eigenFormat = false;
	    String convert = "";
	    boolean fancy_weighting = DEFAULT_FANCY_WEIGHTING;
	    boolean suppress = false;

	    String usage = "\n"+
	    "gwas.Eigenstrat requires 0-1 arguments\n"+
	    "   (1) root name of the source one per family plink files (i.e. source="+sourceRoot+" (default))\n"+
	    "   (2) root name of the target all indiviudals plink files (i.e. target="+targetRoot+" (default))\n"+
	    "   (3) create files for smartpca run (i.e. -create (not the default))\n"+
	    "   (4) calculate eigen values for target population (i.e. -parse (not the default))\n"+
	    "   (5) target files are actually eigensoft format, not plink (i.e. -eigenFormat (not the default))\n"+
	    "   (6) perform fancy weighting when parsing (i.e. fancy="+fancy_weighting+" (default))\n"+
	    "   (7) .evec file to convert to .mds (i.e. convert=plink.pca.evec (not the default))\n"+
	    "   (8) suppress warning/kill regarding marker order (i.e. -suppress (not the default))\n"+
	    "";

	    for (int i = 0; i<args.length; i++) {
		    if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
			    System.err.println(usage);
			    System.exit(1);
		    } else if (args[i].startsWith("source=")) {
			    sourceRoot = args[i].split("=", -1)[1];
			    numArgs--;
		    } else if (args[i].startsWith("target=")) {
			    targetRoot = args[i].split("=", -1)[1];
			    numArgs--;
		    } else if (args[i].startsWith("-create")) {
			    create = true;
			    numArgs--;
		    } else if (args[i].equals("-parse")) {
			    parse = true;
			    numArgs--;
		    } else if (args[i].startsWith("-eigenFormat")) {
		    	eigenFormat = true;
			    numArgs--;
		    } else if (args[i].startsWith("fancy=")) {
		    	fancy_weighting = ext.parseBooleanArg(args[i]);
			    numArgs--;
		    } else if (args[i].startsWith("convert=")) {
			    convert = args[i].split("=", -1)[1];
			    numArgs--;
		    } else if (args[i].startsWith("-suppress")) {
		    	suppress = true;
			    numArgs--;
		    }
	    }
	    if (numArgs!=0) {
		    System.err.println(usage);
		    System.exit(1);
	    }

	    parse = true;
	    
	    
	    try {
	    	if (create) {
	    		generateFiles(sourceRoot);
	    	} else if (!convert.equals("")) {
	    		convertToMDS(convert);
	    	} else if (parse) {
	    		parse(sourceRoot, targetRoot, eigenFormat, fancy_weighting, suppress);
	    	}
	    } catch (Exception e) {
		    e.printStackTrace();
	    }
    }
}
