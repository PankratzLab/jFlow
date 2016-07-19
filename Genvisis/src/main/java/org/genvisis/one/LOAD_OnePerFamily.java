package org.genvisis.one;

import java.io.*;
import java.util.*;

import org.genvisis.common.*;

public class LOAD_OnePerFamily {
	public static void pick(String dir, String pedigreeFile, String phenotypeFile) {
        PrintWriter writer;
        String[] fams;
        String trav;
        Hashtable<String,String> dxHash, aooHash;
        Hashtable<String,Vector<String>> pedHash;
        Vector<String> v;
        int pick;
        String dx; 
        int dxLevel, minDxLevel, minAOO;
        
        dxHash = HashVec.loadFileToHashString(dir+phenotypeFile, "IndID", new String[] {"Dx"}, "\t");
        aooHash = HashVec.loadFileToHashString(dir+phenotypeFile, "IndID", new String[] {"AgeDem"}, "\t");
        
        pedHash = HashVec.loadFileToHashVec(dir+pedigreeFile, 0, new int[] {1}, "\t", false, false);
        try {
        	fams = HashVec.getKeys(pedHash);
	        writer = new PrintWriter(new FileWriter(dir+"picks.xln"));
	        for (int i = 0; i<fams.length; i++) {
	        	v = pedHash.get(fams[i]);
	        	pick = -1;
        		minAOO = 999;
        		minDxLevel = 99;
	        	for (int j = 0; j<v.size(); j++) {
	        		trav = v.elementAt(j);
	        		dx = dxHash.get(trav);
	        		if (dx != null) {
	        			dxLevel = Integer.parseInt(dx);
	        			if (dxLevel == 1 || dxLevel == 2 || dxLevel == 3 || dxLevel == 7) {
	        				if (dxLevel == minDxLevel && getAOO(aooHash.get(trav)) < minAOO) {
	        					pick = j;
	        					minAOO = getAOO(aooHash.get(trav));
	        				} else if (dxLevel < minDxLevel) {
	        					pick = j;
	        					minDxLevel = dxLevel;
	        					minAOO = getAOO(aooHash.get(trav));
	        				}
	        			} else if (dxLevel == 9 || (dxLevel > 3 && dxLevel < 7)) {
	        				// do nothing
	        			} else {
	        				System.err.println("Error - invalid dxLevel ("+dx+")");
	        			}
	        		}
                }
	        	writer.println(fams[i]+"\t"+(pick==-1?".\t.\t.":v.elementAt(pick)+"\t"+minDxLevel+"\t"+minAOO));
            }
	        writer.close();
        } catch (Exception e) {
	        System.err.println("Error writing to "+dir+"picks.xln");
	        e.printStackTrace();
        }
	}
	
	public static int getAOO(String aoo) {
		return aoo == null || aoo.equals(".")?999:Integer.parseInt(aoo);
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\LOAD\\";
		String pedigree = "finalPedigree_geno.pre-trimmed.pre";
		String phenotypeFile = "Phenotype.dat";
	
		String usage = "\n"+
		"one.LOAD_OnePerFamily requires 0-1 arguments\n"+
		"   (1) pedigree file (i.e. ped="+phenotypeFile+" (default))\n"+
		"   (2) phenotype file (i.e. pheno="+phenotypeFile+" (default))\n"+
		"";
	
		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("ped=")) {
				pedigree = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("pheno=")) {
				phenotypeFile = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			pick(dir, pedigree, phenotypeFile);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
