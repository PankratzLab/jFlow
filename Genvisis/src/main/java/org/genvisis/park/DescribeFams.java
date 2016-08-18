// incorporate the info from demographics at some point
package org.genvisis.park;

import java.io.*;
import java.util.*;

import org.genvisis.common.*;

public class DescribeFams {
	public static void describe() {
		BufferedReader reader;
        PrintWriter writer;
        String[] line;
        Hashtable<String,Vector<String>> famInds = new Hashtable<String,Vector<String>>();
//        Hashtable<String,Vector<String>> sibship = new Hashtable<String,Vector<String>>();
        
        Hashtable<String,String> dxs = new Hashtable<String,String>();
        Vector<String> v;
        String[] fams, sibships, vpdSibships;
//        String famid, uniqueid;
        int vpd, nvpd, rpd, noev, nrpd;
        String dx;
		String[][] majorSibships, majorVPDsibships;
		int sibIndex, vpdSibIndex;
		boolean dominantTransmission, sfh;
        
//        try {
//	        reader = tools.getNinfoReader(1);
//	        reader.readLine();
//	        while (reader.ready()) {
//	        	line = reader.readLine().trim().split("[\\s]+");
//	        	famid = line[0];
//	        	uniqueid = tools.getUniqueID(line[0], line[1]);
//	        	HashVec.addToHashVec(famInds, famid, uniqueid, true);
//	        }
//	        reader.close();
//        } catch (IOException ioe) {
//	        System.err.println("Error reading file ninfo1");
//	        System.exit(2);
//        }

		dxs = tools.getBestPDdx();
        
        try {
	        reader = tools.getNinfoReader(2);
	        reader.readLine();
	        while (reader.ready()) {
	        	line = reader.readLine().trim().split("[\\s]+");
//	        	famid = line[0];
//	        	v = famInds.get(famid);
//	        	uniqueid = tools.getUniqueID(line[0], line[1]);
//	        	if (v != null && v.contains(uniqueid)) {
//		        	HashVec.addToHashVec(sibship, famid, line[4]+"\t"+line[5], true);
//	        	}

	        	HashVec.addToHashVec(famInds, line[0], tools.getUniqueID(line[0], line[1]), true);
	        }
	        reader.close();
        } catch (IOException ioe) {
	        System.err.println("Error reading file ninfo1");
	        System.exit(2);
        }
        
		majorSibships = tools.getCombonentsOfTheMajorSibships(true, false);
		sibships = Matrix.extractColumn(majorSibships, 0);
		majorVPDsibships = tools.getCombonentsOfTheMajorSibships(true, true);
		vpdSibships = Matrix.extractColumn(majorVPDsibships, 0);
        try {
	        writer = new PrintWriter(new FileWriter(tools.CRF_DIR+"FamilyStats.xln"));
	        fams = HashVec.getKeys(famInds);
	        writer.println("FamID\t# Affected\t# Unaffected\tVPD\tNVPD\tRPD\tNOEV\tNRPD\tSizeOfLargestSibship\tSizeOfLargestVPDSibship\tDominantTransmission\tSFH");
	        for (int i = 0; i<fams.length; i++) {
	        	vpd = nvpd = rpd = noev = nrpd = 0;
	        	v = famInds.get(fams[i]);
	        	if (fams[i].equals("70012")) {
	        		System.out.println("hola");
	        	}
	        	for (int j = 0; j<v.size(); j++) {
	        		dx = dxs.get(v.elementAt(j));
	        		if (dx == null) {
//	        			System.out.println(v.elementAt(j));
	        		} else if (dx.equals("VPD")) {
	        			vpd++;
	        		} else if (dx.equals("NVPD")) {
	        			nvpd++;
	        		} else if (dx.equals("RPD")) {
	        			rpd++;
	        		} else if (dx.equals("NOEV")) {
	        			noev++;
	        		} else if (dx.equals("NRPD")) {
	        			nrpd++;
	        		} else {
	        			System.err.println("Error - don't know what to do with dx: "+dx);
	        		}
                }
	        	
	        	sibIndex = ext.indexOfStr(fams[i], sibships);
	        	vpdSibIndex = ext.indexOfStr(fams[i], vpdSibships);
//	        	dominantTransmission = sfh = false;
	        	dominantTransmission = (majorSibships[sibIndex].length > 1 
	        			&& (tools.isAffected(dxs.get(fams[i]+"\t"+majorSibships[sibIndex][1])).equals("1")
	        					|| tools.isAffected(dxs.get(fams[i]+"\t"+majorSibships[sibIndex][2])).equals("1"))
	        			|| majorVPDsibships[sibIndex].length > 1 
	        			&& (tools.isAffected(dxs.get(fams[i]+"\t"+majorVPDsibships[sibIndex][1])).equals("1")
	        					|| tools.isAffected(dxs.get(fams[i]+"\t"+majorVPDsibships[sibIndex][2])).equals("1")));
	        	sfh = majorVPDsibships[sibIndex].length - 3 >= 2 && (dominantTransmission || vpd+nvpd+rpd >= 4);
	        	
	        	if (vpd+nvpd+noev>0) {
		        	writer.println(fams[i]+"\t"+(vpd+nvpd+rpd)+"\t"+(noev+nvpd)+"\t"+vpd+"\t"+nvpd+"\t"+rpd+"\t"+noev+"\t"+nrpd
		        			+"\t"+(majorSibships[sibIndex].length > 1?majorSibships[sibIndex].length-3:0)
		        			+"\t"+(majorVPDsibships[vpdSibIndex].length > 1?majorVPDsibships[vpdSibIndex].length-3:0)
		        			+"\t"+(dominantTransmission?1:0)
		        			+"\t"+(sfh?1:0));
	        	}
            }
	        writer.close();
        } catch (Exception e) {
	        System.err.println("Error writing to "+"FamilyStats.xln");
	        e.printStackTrace();
        }
	}
	
	public static void genderBalances() {
		PrintWriter writer;
		Hashtable<String,String> affectionStatus, genderLookup;
		String[][] majorSibships;
		int[][] counts;
		int gender, affFather, affMother;

		affectionStatus = tools.getBestPDdx();
		genderLookup = tools.getGender();
		majorSibships = tools.getCombonentsOfTheMajorSibships(false, false);

		try {
	        writer = new PrintWriter(new FileWriter(tools.CRF_DIR+"GenderBalance.csv"));
	        writer.println("Family,AffectedFather,AffectedMother,ParentalScore,VPDmales,VPDfemales,VPDscore,VPDweightedScore,AffectedMales,AffectedFemales,AffectedScore,AffectedWeightedScore");
	        for (int i = 0; i<majorSibships.length; i++) {
	        	if (majorSibships[i].length > 1) {
	        		counts = new int[2][3];
	        		affFather = tools.isAffected(affectionStatus, majorSibships[i][0]+"\t"+majorSibships[i][1])?1:0;
	        		affMother = tools.isAffected(affectionStatus, majorSibships[i][0]+"\t"+majorSibships[i][2])?1:0;
	        		for (int j = 3; j<majorSibships[i].length; j++) {
	        			gender = Integer.parseInt(genderLookup.get(majorSibships[i][0]+"\t"+majorSibships[i][j]));
	        			
			        	if (tools.isVPD(affectionStatus.get(majorSibships[i][0]+"\t"+majorSibships[i][j])).equals("1")) {
			        		counts[0][gender]++;
			        	}			        	
			        	if (tools.isAffected(affectionStatus, majorSibships[i][0]+"\t"+majorSibships[i][j])) {
			        		counts[1][gender]++;
			        	}
                    }
		        	writer.println(majorSibships[i][0]+","+affFather+","+affMother+","+(affFather-affMother)+","+parseScores(counts[0])+","+parseScores(counts[1]));
	        	}
            }
	        writer.close();
        } catch (Exception e) {
	        System.err.println("Error writing to "+tools.CRF_DIR+"GenderBalance.csv");
	        e.printStackTrace();
        }
	}

	public static String parseScores(int[] counts) {
		return counts[1]+","+counts[2]+","+(counts[1]-counts[2])+","+ext.formDeci(counts[1]+counts[2]>0?((double)counts[1]/(double)(counts[1]+counts[2])*2-1):0, 3);		
	}
	
	public static void main(String[] args) {
	    int numArgs = args.length;
	    boolean genderBalance = false;

	    String usage = "\n"+
	    "park.DescribeFams requires 0-1 arguments\n"+
	    "   (1) describe fams (default)\n"+
	    "   (2) describe gender balances in major sibships (i.e. -genderBalance (not the default))\n"+
	    "";

	    for (int i = 0; i<args.length; i++) {
		    if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
			    System.err.println(usage);
			    System.exit(1);
		    } else if (args[i].startsWith("-genderBalance")) {
			    genderBalance = true;
			    numArgs--;
		    }
	    }
	    if (numArgs!=0) {
		    System.err.println(usage);
		    System.exit(1);
	    }
	    try {
		    if (genderBalance) {
		    	genderBalances();
		    } else {
		    	describe();
		    }
	    } catch (Exception e) {
		    e.printStackTrace();
	    }
    }
}
