package org.genvisis.bioinformatics;

import java.io.*;
import java.util.*;

import org.genvisis.common.*;
import org.genvisis.filesys.Segment;

public class FilterSNPsByRegion {

	public static final int DEFAULT_BUFFER = 0;

	public static void fromParameters(String filename, Logger log) {
		BufferedReader reader;
        PrintWriter writer;
        String[] line, regionNames;
        Segment[] regions;
        Segment mark;
        String regionsFile, lookupFile, outputFilename;
        Vector<String> snps, snpsWithRegionNumbers;
        int snpCol, chrCol, startCol, stopCol, posCol;
        boolean passes;
        int buffer;
        Vector<String> paramV;
        boolean header, ucsc, regionNumber;
        int regionNameIndex;
        
        buffer = -1;
        
		paramV = Files.parseControlFile(filename, "filterSNPs", new String[] {"regions.txt header 0 1 2 ucsc out=file.out regionNumber regionName=3", "plink.bim 1 0 3", "additonalAnnotation 0 5 6", "buffer="+DEFAULT_BUFFER}, log);
		if (paramV == null) {
			return;
		}

    	line = paramV.remove(0).trim().split("[\\s]+");
    	regionsFile = line[0];
    	chrCol = startCol = stopCol = -1;
    	header = false;
    	ucsc = false;
    	regionNumber = false;
    	regionNameIndex = -1;
    	outputFilename = ext.rootOf(regionsFile)+"_filtered.xln";
    	for (int j = 1; j<line.length; j++) {
    		if (line[j].equals("ucsc")) {
        		ucsc = true;
    		} else if (line[j].equals("header")) {
    			header = true;
    		} else if (line[j].equals("regionNumber")) {
    			regionNumber = true;
    		} else if (line[j].startsWith("regionName=")) {
    			regionNameIndex = ext.parseIntArg(line[j]);
    		} else if (line[j].startsWith("out=")) {
    			outputFilename = line[j].split("=")[1];
    		} else if (chrCol == -1) {
    			chrCol = Integer.parseInt(line[j]);
    		} else if (startCol == -1) {
    			startCol = Integer.parseInt(line[j]);
    		} else if (stopCol == -1) {
    			stopCol = Integer.parseInt(line[j]);
    		} else {
    			System.err.println("Error - what am I supposed to do with '"+line[j]+"'?");
    		}
        }
    	if (chrCol == -1) {
    		System.err.println("Warning - assuming regions are defined in columns 0, 1 and 2");
    		chrCol = 0;
    		startCol = 1;
    		stopCol = 2;
    	}
    	if (!ucsc && startCol == -1) {
    		System.err.println("Warning - since only one column was specified, assuming you meant to include ucsc flag but didn't");
    		ucsc = true;
    	}
    	log.report("Loading regions from '"+regionsFile+"'");
    	if (ucsc) {
    		regions = Segment.loadUCSCregions(regionsFile, chrCol, header, log);
    	} else {
    		regions = Segment.loadRegions(regionsFile, chrCol, startCol, stopCol, header);
    	}
    	if (regionNameIndex >= 0) {
    		regionNames = HashVec.loadFileToStringArray(regionsFile, header, new int[] {regionNameIndex}, false);
    	} else {
    		regionNames = null;
    	}

    	log.report("Looking up SNPs within these regions...");
    	line = paramV.remove(0).trim().split("[\\s]+");
    	lookupFile = line[0];
    	snpCol = Integer.parseInt(line[1]);
    	chrCol = Positions.chromosomeNumber(line[2]);
    	posCol = Integer.parseInt(line[3]);

    	for (int i = 0; i < paramV.size(); i++) {
        	if (paramV.elementAt(i).startsWith("buffer=")) {
        		buffer = Integer.parseInt(paramV.elementAt(i).split("=")[1]);
        		paramV.removeElementAt(i);
        		i--;
        	}
		}
        
        if (buffer == -1) {
        	System.out.println("No basepair buffer defined; using default ("+DEFAULT_BUFFER+" bp)");
        	buffer = DEFAULT_BUFFER;
        }
    	
    	snps = new Vector<String>(1000);
    	snpsWithRegionNumbers = new Vector<String>(1000);
    	snpsWithRegionNumbers.add("MarkerName\tRegionNumber"+(regionNameIndex>=0?"\tRegionName":""));
    	try {
            reader = new BufferedReader(new FileReader(lookupFile));
            writer = new PrintWriter(new FileWriter(outputFilename));
            while (reader.ready()) {
            	line = reader.readLine().trim().split("[\\s]+");
            	passes = false;
            	mark = new Segment(Positions.chromosomeNumber(line[chrCol]), Integer.parseInt(line[posCol]), Integer.parseInt(line[posCol]));
            	for (int i = 0; i<regions.length && !passes; i++) {
            		if (mark.overlaps(regions[i])) {
        				passes = true;
                		snpsWithRegionNumbers.add(line[snpCol]+"\t"+(i+1)+(regionNameIndex>=0?"\t"+regionNames[i]:""));
        			}
                }
            	if (passes) {
            		snps.add(line[snpCol]);
            	}
            }
            reader.close();
            writer.close();
        } catch (FileNotFoundException fnfe) {
        	log.reportError("Error: file \""+lookupFile+"\" not found in current directory");
        	log.reportException(fnfe);
            return;
        } catch (IOException ioe) {
        	log.reportError("Error reading file \""+lookupFile+"\"");
        	log.reportException(ioe);
            return;
        }
    	log.report("Found "+snps.size()+" SNPs within these regions");
    	
    	if (regionNumber) {
    		Files.writeList(Array.toStringArray(snpsWithRegionNumbers), "snp_region_matchup.dat");
    		paramV.add("snp_region_matchup.dat 0 1=RegionNumber"+(regionNameIndex>=0?"\t2=RegionName":""));
    	}
    	
    	Files.combine(Array.toStringArray(snps), Array.toStringArray(paramV), "MarkerName", outputFilename, log, true);
	}
}
