package cnv.filesys;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Hashtable;

import common.Array;
import common.Logger;
import common.ext;

public class ABLookup {
	public static final String DEFAULT_AB_FILE = "AB_lookup.dat";

	private String[] markerNames;
	private char[][] lookup;
	
	public ABLookup() {
		markerNames = new String[0];
		lookup = new char[0][0];
	}
	
	public ABLookup(String[] markerNames, String filename, boolean verbose, boolean allOrNothing) {
		Hashtable<String,char[]> lookupHash;
		int count;
		
		count = 0;
		this.markerNames = markerNames;
        lookup = new char[markerNames.length][];
        if (filename.toLowerCase().endsWith(".csv")) {
        	lookupHash = generateABLookupHashFromCSV(filename);
        } else {
        	lookupHash = generateABLookupHash(filename);
        }
        for (int i = 0; i<markerNames.length; i++) {
//        	System.out.println(markerNames[i]);
//        	System.out.println(lookupHash.containsKey(markerNames[i]));
//        	System.out.println(lookupHash.get(markerNames[i]));
    		lookup[i] = lookupHash.get(markerNames[i]);
    		if (lookup[i] == null) {
    			if (verbose) {
        			System.err.println("Error - no AB value for marker '"+markerNames[i]+"'");
    			}
    			count++;
    		}
        }
        if (count > 0) {
        	System.err.println("Warning - there "+(count>1?"were ":"was only ")+count+" marker"+(count>1?"s":"")+" without an AB value");
        }
        if (allOrNothing && count > 0) {
			lookup = null;
        }
	}

	public void parseFromGenotypeClusterCenters(Project proj) {
        PrintWriter writer, writer2;
        String trav;
        String[] samples;
        String[][] genotypesSeen;
        double[][] meanThetasForGenotypes;
        int[][] countsForGenotypes;
        Sample fsamp;
        byte[] genotypes;
        float[] thetas;
        char[] alleles;
        double travD;
        int order;
        MarkerSet markerSet;
        
        samples = proj.getSamples();
        markerSet = proj.getMarkerSet();
        markerNames = markerSet.getMarkerNames();
        
        genotypesSeen = new String[markerNames.length][3];
        meanThetasForGenotypes = new double[markerNames.length][3];
        countsForGenotypes = new int[markerNames.length][3];
        for (int i = 0; i<samples.length; i++) {
        	System.out.println((i+1)+" of "+samples.length);
        	fsamp = proj.getFullSampleFromRandomAccessFile(samples[i]);
        	thetas = fsamp.getThetas();
        	genotypes = fsamp.getForwardGenotypes();
        	
        	for (int j = 0; j<markerNames.length; j++) {
        		if (genotypes[j] > 0) {
        			trav = Sample.ALLELE_PAIRS[genotypes[j]];
        			if (trav.charAt(0) != trav.charAt(1)) {
        				if (genotypesSeen[j][1] == null || genotypesSeen[j][1].equals(trav)) {
        					genotypesSeen[j][1] = trav;
        			        meanThetasForGenotypes[j][1] += thetas[j];
        			        countsForGenotypes[j][1]++;
        				} else {
        					System.err.println("Error - different heterozygote ("+trav+") than the on seen previously ("+genotypesSeen[j][1]+") for marker "+markerNames[j]+" and sample "+samples[i]);
        				}
        			} else {
        				if (genotypesSeen[j][0] == null || genotypesSeen[j][0].equals(trav)) {
        					genotypesSeen[j][0] = trav;
        			        meanThetasForGenotypes[j][0] += thetas[j];
        			        countsForGenotypes[j][0]++;
        				} else if (genotypesSeen[j][2] == null || genotypesSeen[j][2].equals(trav)) {
            				genotypesSeen[j][2] = trav;
        			        meanThetasForGenotypes[j][2] += thetas[j];
        			        countsForGenotypes[j][2]++;
        				} else {
        					System.err.println("Error - different genotype ("+trav+") than anything seen before ("+Array.toStr(genotypesSeen[j], "/")+") for marker "+markerNames[j]+" and sample "+samples[i]);
        				}
        			}
        		}
            }
        }
        
        try {
	        writer = new PrintWriter(new FileWriter(proj.getProjectDir()+"AB_breakdown.xln"));
	        writer.println("Marker\tG11\tG12\tG22\tG11 counts\tG12 counts\tG22 counts\tMean Theta G11\tMean Theta G12\tMean Theta G22\torder\tA allele\tB allele");
	        writer2 = new PrintWriter(new FileWriter(proj.getProjectDir()+"poss_"+DEFAULT_AB_FILE));
	        writer2.println("Marker\tAlleleA\tAlleleB");
//	        int countMissing;
	        lookup = new char[markerNames.length][];
	        for (int i = 0; i<markerNames.length; i++) {
//	        	countMissing = 0;
	        	for (int j = 0; j<3; j++) {
	        		meanThetasForGenotypes[i][j] /= countsForGenotypes[i][j];
//	        		if (countsForGenotypes[i][j] == 0) {
//	        			countMissing++;
//	        		}
                }
	        	writer.print(markerNames[i]+"\t"+genotypesSeen[i][0]+"\t"+genotypesSeen[i][1]+"\t"+genotypesSeen[i][2]+"\t"+countsForGenotypes[i][0]+"\t"+countsForGenotypes[i][1]+"\t"+countsForGenotypes[i][2]+"\t"+meanThetasForGenotypes[i][0]+"\t"+meanThetasForGenotypes[i][1]+"\t"+meanThetasForGenotypes[i][2]);
	        	
    			alleles = new char[] {'N', 'N'};
	        	for (int j = 0; j<3; j++) {
	        		for (int k = 0; k<2; k++) {
        				if (countsForGenotypes[i][j] > 0) {
	        				if (alleles[0] == 'N') {
	        					alleles[0] = genotypesSeen[i][j].charAt(k);
	        				} else if (alleles[0] == genotypesSeen[i][j].charAt(k)) {
	        				} else if (alleles[1] == 'N') {
	        					alleles[1] = genotypesSeen[i][j].charAt(k);
	        				} else if (alleles[1] == genotypesSeen[i][j].charAt(k)) {
	        				} else {
	        					System.err.println("Error - too many alleles for "+markerNames[i]+" first "+alleles[0]+" and "+alleles[1]+" and now "+genotypesSeen[i][j].charAt(k));
	        				}
        				}
                    }
                }
	        	order = 0;
	        	travD = -1;
	        	for (int j = 0; j<3; j++) {
    				if (countsForGenotypes[i][j] > 0) {
    					if (travD == -1) {
    						travD = meanThetasForGenotypes[i][j];
    					} else if (meanThetasForGenotypes[i][j] > travD) {
    						if (order == 0) {
    							order = 1;
    						} else if (order != 1) {
    							order = -1;
    						}
    					} else if (meanThetasForGenotypes[i][j] < travD) {
    						if (order == 0) {
    							order = 2;
    						} else if (order != 2) {
    							order = -1;
    						}
    					} else {
    						System.err.println("Error - equal meanThetas?");
    					}
    				}
	        	}
	        	if (order == 0 && meanThetasForGenotypes[i][0] > 0.50) {
	        		order = 2;
	        	}
	        	if (order == -1 && countsForGenotypes[i][0] > 0 && countsForGenotypes[i][2] > 0) {
	        		order = meanThetasForGenotypes[i][0]<meanThetasForGenotypes[i][2]?-1:-2;
	        	}
	        	if (order==2||order==-2) {
	        		lookup[i] = new char[] {alleles[1], alleles[0]};
	        	} else {
	        		lookup[i] = alleles;
	        	}
	        	writer.println("\t"+order+"\t"+lookup[i][0]+"\t"+lookup[i][1]);
	        	writer2.println(markerNames[i]+"\t"+lookup[i][0]+"\t"+lookup[i][1]);
            }
	        writer.close();
	        writer2.close();
        } catch (Exception e) {
	        System.err.println("Error writing to either "+"AB_breakdown.xln"+" or "+"poss_"+DEFAULT_AB_FILE);
	        e.printStackTrace();
        }
		
	}
	
	public void parseFromOriginalGenotypes(Project proj) {
        String fullGenotype;
        String[] samples;
        String[][] genotypesSeen;
        int[][] countsForGenotypes;
        Sample fsamp;
        byte[] abGenotypes, fullGenotypes;
        MarkerSet markerSet;
        
        samples = proj.getSamples();
        markerSet = proj.getMarkerSet();
        markerNames = markerSet.getMarkerNames();
        
        genotypesSeen = new String[markerNames.length][3];
        countsForGenotypes = new int[markerNames.length][3];
        for (int i = 0; i<samples.length; i++) {
        	System.out.println((i+1)+" of "+samples.length);
        	fsamp = proj.getFullSampleFromRandomAccessFile(samples[i]);
        	abGenotypes = fsamp.getAB_Genotypes();
        	fullGenotypes = fsamp.getForwardGenotypes();
        	
        	for (int j = 0; j<markerNames.length; j++) {
        		if (abGenotypes[j] >= 0) {
        			fullGenotype = Sample.ALLELE_PAIRS[fullGenotypes[j]];
        			if (genotypesSeen[j][abGenotypes[j]] == null) {
        				genotypesSeen[j][abGenotypes[j]] = fullGenotype;
        			} else if (!fullGenotype.equals(genotypesSeen[j][abGenotypes[j]])) {
        				System.err.println("Error - different heterozygote ("+fullGenotype+") than the on seen previously ("+genotypesSeen[j][1]+") for marker "+markerNames[j]+" and sample "+samples[i]);
        			}
    				countsForGenotypes[j][abGenotypes[j]]++;
        		}
            }
        }
        
        lookup = new char[markerNames.length][2];
        for (int i = 0; i<markerNames.length; i++) {
			lookup[i] = new char[] {'N', 'N'};
			if (countsForGenotypes[i][1] > 0 && genotypesSeen[i][1].charAt(0) == genotypesSeen[i][1].charAt(1)) {
				System.err.println("Error - impossible heterozygote ("+genotypesSeen[i][1]+") for marker "+markerNames[i]);
			}
			if (countsForGenotypes[i][0] > 0) {
				if (genotypesSeen[i][0].charAt(0) != genotypesSeen[i][0].charAt(1)) {
					System.err.println("Error - impossible A/A homozygote ("+genotypesSeen[i][0]+") for marker "+markerNames[i]);
				} else {
					lookup[i][0] = genotypesSeen[i][0].charAt(0);
				}
			}
			if (countsForGenotypes[i][2] > 0) {
				if (genotypesSeen[i][2].charAt(0) != genotypesSeen[i][2].charAt(1)) {
					System.err.println("Error - impossible B/B homozygote ("+genotypesSeen[i][2]+") for marker "+markerNames[i]);
				} else {
					lookup[i][1] = genotypesSeen[i][2].charAt(0);
				}
			}
			if (countsForGenotypes[i][1] > 0 && countsForGenotypes[i][0]+countsForGenotypes[i][2]==0) {
				System.err.println("Error - only heterozygotes present for marker "+markerNames[i]+"; cannot determine which allele is the A allele versus the B allele");
			} else if (countsForGenotypes[i][1] > 0 && countsForGenotypes[i][0] == 0 ) {
				if (genotypesSeen[i][1].charAt(0) == lookup[i][1]) {
					lookup[i][0] = genotypesSeen[i][1].charAt(1);
				} else if (genotypesSeen[i][1].charAt(1) == lookup[i][1]) {
					lookup[i][0] = genotypesSeen[i][1].charAt(0);
				} else {
					System.err.println("Error - heterozygote for marker "+markerNames[i]+" ("+genotypesSeen[i][1]+") does not match up with the B allele ("+lookup[i][1]+")");
				} 
			} else if (countsForGenotypes[i][1] > 0 && countsForGenotypes[i][2] == 0 ) {
				if (genotypesSeen[i][1].charAt(0) == lookup[i][0]) {
					lookup[i][1] = genotypesSeen[i][1].charAt(1);
				} else if (genotypesSeen[i][1].charAt(1) == lookup[i][0]) {
					lookup[i][1] = genotypesSeen[i][1].charAt(0);
				} else {
					System.err.println("Error - heterozygote for marker "+markerNames[i]+" ("+genotypesSeen[i][1]+") does not match up with the A allele ("+lookup[i][1]+")");
				} 
			}
        }
	}

	public static Hashtable<String,char[]> generateABLookupHash(String filename) {
		BufferedReader reader;
		String[] line;
		Hashtable<String,char[]> hash;
		
        hash = new Hashtable<String,char[]>();
		try {
            reader = new BufferedReader(new FileReader(filename));
            while (reader.ready()) {
            	line = reader.readLine().trim().split("[\\s]+");
            	hash.put(line[0], new char[] {line[1].charAt(0), line[2].charAt(0)});
            }
            reader.close();
            return hash;
        } catch (FileNotFoundException fnfe) {
        	return null;
        } catch (IOException ioe) {
            System.err.println("Error reading file \""+DEFAULT_AB_FILE+"\"");
            return null;
        }
	}

	public static Hashtable<String,char[]> generateABLookupHashFromCSV(String filename) {
		BufferedReader reader;
		String[] line;
		Hashtable<String,char[]> hash;
		int[] indices;
		
        hash = new Hashtable<String,char[]>();
		try {
            reader = new BufferedReader(new FileReader(filename));
            line = reader.readLine().trim().split(",");
            indices = ext.indexFactors(new String[][] {{"Name", "MarkerName"}, {"SNP"}}, line, false, true, true, new Logger(), true);
            while (reader.ready()) {
            	line = reader.readLine().trim().split(",");
            	hash.put(line[indices[0]], new char[] {line[indices[1]].charAt(1), line[indices[1]].charAt(3)});
            }
            reader.close();
            return hash;
        } catch (FileNotFoundException fnfe) {
        	return null;
        } catch (IOException ioe) {
            System.err.println("Error reading file \""+DEFAULT_AB_FILE+"\"");
            return null;
        }
	}

	public char[][] getLookup() {
        return lookup;
	}
	
	public void writeToFile(String outfile) {
		PrintWriter writer;
		
		try {
	        writer = new PrintWriter(new FileWriter(outfile));
	        writer.println("Marker\tA\tB");
	        for (int i = 0; i<markerNames.length; i++) {
	        	writer.println(markerNames[i]+"\t"+lookup[i][0]+"\t"+lookup[i][1]);
            }
	        writer.close();
        } catch (Exception e) {
	        System.err.println("Error writing to "+outfile);
	        e.printStackTrace();
        }
	}
	

	public static void applyABLookupToFullSampleFiles(Project proj) {
		ABLookup abLookup;
        Sample fsamp;
        String[] samples;
        String genotype;
        byte[] forwardGenotypes, abGenotypes;
        
		abLookup = new ABLookup(proj.getMarkerNames(), proj.getFilename(Project.AB_LOOKUP_FILENAME), false, false);
        samples = proj.getSamples();
        for (int i=0; i<samples.length; i++) {
        	System.out.println((i+1)+" of "+samples.length);
        	fsamp = proj.getFullSampleFromRandomAccessFile(samples[i]);
        	forwardGenotypes = fsamp.getForwardGenotypes();
        	if (forwardGenotypes.length != abLookup.markerNames.length) {
        		System.err.println("Error - mismatched array lengths for forwardGenotypes and abLookup markerNames");
        	}
        	abGenotypes = new byte[forwardGenotypes.length];
        	for (int j=0; j<abLookup.markerNames.length; j++) {
        		if (forwardGenotypes[j] == 0) {
        			abGenotypes[j] = -1;
        		} else {
//	        		genotype = FullSample.ALLELE_PAIRS[forwardGenotypes[j]];
	        		genotype = Sample.ALLELE_PAIRS[forwardGenotypes[j]];
	        		for (int k = 0; k < 2; k++) {
	        			if (genotype.charAt(k) == abLookup.lookup[j][1]) {
	        				abGenotypes[j]++;
	        			} else if (genotype.charAt(k) != abLookup.lookup[j][0]) {
	        				System.err.println("Error - mismatched alleles found "+genotype.charAt(k)+" for marker "+abLookup.markerNames[j]+" for sample "+samples[i]+" (expecting "+abLookup.lookup[j][0]+"/"+abLookup.lookup[j][1]+")");
	        			}
					}
        		}
        	}
        	fsamp.setAB_Genotypes(abGenotypes);
        	fsamp.saveToRandomAccessFile(proj.getDir(Project.SAMPLE_DIRECTORY)+samples[i]+Sample.SAMPLE_DATA_FILE_EXTENSION);
        }
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		Project proj;
		String filename = Project.DEFAULT_PROJECT;
		boolean parseFromOriginalGenotypes = false;
		boolean parseFromGenotypeClusterCenters = false;
		String outfile = DEFAULT_AB_FILE;
		ABLookup abLookup;
		String mapFile = null;
		boolean applyAB = false;

		String usage = "\n" +
				"cnv.filesys.ABLookup requires 0-1 arguments\n" +
				"   (1) project file (i.e. proj="+filename+" (default))\n"+
				"   (2) name of output file (i.e. out="+outfile+" (default))\n" + 
				"  AND\n" + 
				"   (3) parse ABLookup from centroids (i.e. -parseFromGenotypeClusterCenters (not the default))\n" + 
				"  OR\n" + 
				"   (3) parse ABLookup from existing original genotypes (i.e. -parseFromOriginalGenotypes (not the default))\n" + 
				"  OR\n" + 
				"   (3) parse ABLookup from [Illumina] SNP Table (i.e. mapFile=SNP_Map.csv (not the default))\n" + 
				"  OR\n" + 
				"   (3) apply the project's AB lookup to all Sample files in project (i.e. -applyAB (not the default))\n" + 
				"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("out=")) {
				outfile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].equalsIgnoreCase("-parseFromOriginalGenotypes")) {
				parseFromOriginalGenotypes = true;
				numArgs--;
			} else if (args[i].equalsIgnoreCase("-parseFromGenotypeClusterCenters")) {
				parseFromGenotypeClusterCenters = true;
				numArgs--;
			} else if (args[i].startsWith("mapFile=")) {
				mapFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].equalsIgnoreCase("-applyAB")) {
				applyAB = true;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		
		try {
			proj = new Project(filename, false);
			if (mapFile != null) {
				abLookup = new ABLookup(proj.getMarkerNames(), mapFile, true, true);
				abLookup.writeToFile(proj.getProjectDir()+outfile);
			} else if (parseFromOriginalGenotypes) {
				abLookup = new ABLookup();
				abLookup.parseFromOriginalGenotypes(proj);
				abLookup.writeToFile(proj.getProjectDir()+outfile);
			} else if (parseFromGenotypeClusterCenters) {
				abLookup = new ABLookup();
				abLookup.parseFromGenotypeClusterCenters(proj);
				abLookup.writeToFile(proj.getProjectDir()+outfile);
			} else if (applyAB) {
				applyABLookupToFullSampleFiles(proj);
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
