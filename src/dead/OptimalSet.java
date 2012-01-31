package dead;

import java.io.*;
import java.util.*;
import common.*;

public class OptimalSet {
	public static final String DEFAULT_FILENAME = "C:\\Documents and Settings\\npankrat\\My Documents\\nostalgia\\ninjas.txt";
//	public static final double[] THRESHOLDS = {0, 0.05, 0.10, 0.20, 0.53, 1.00};
//	public static final double[] THRESHOLDS = {0.30, 0.35, 0.40, 0.45, 0.55};
	public static final double[] THRESHOLDS = {0.30, 0.45, 0.60, 0.75, 0.90, 1.0};
	public static final int DOLLAR = 2592000;
	public static final int CREW = 271;
//	public static final int CREW = 1000;
//	public static final int CREW = 1000;

//	public static final String DEFAULT_FILENAME = "C:\\Documents and Settings\\npankrat\\My Documents\\nostalgia\\wizards.txt";
//	public static final double[] THRESHOLDS = {0.01, 0.02, 0.04, 0.06, 0.10, 0.20, 1};
//	public static final int DOLLAR = 857620000;
//	public static final int CREW = 331;
	
	
	public static void optimize(String filename) {
        PrintWriter writer;
        String[] line, types;
        Hashtable<String,Vector<String>> hash = new Hashtable<String,Vector<String>>();
        Vector<String> v = new Vector<String>();
        int sumUpkeep;
        long time;
        String[][] names;
        int[][] atk, def, upkeep;
        int[][][] bestSets, bestSetValues;
        int[] limits, counters, sums;
        
        
    	hash = HashVec.loadFileToHashVec(filename, 0, new int[] {1, 2, 3, 4}, "\t", true, false);
        types = HashVec.getKeys(hash);
        names = new String[types.length][];
        atk = new int[types.length][];
        def = new int[types.length][];
        upkeep = new int[types.length][];
        for (int i = 0; i<types.length; i++) {
        	v = hash.get(types[i]);
        	names[i] = new String[v.size()];
        	atk[i] = new int[v.size()];
        	def[i] = new int[v.size()];
        	upkeep[i] = new int[v.size()];
        	for (int j = 0; j<v.size(); j++) {
	        	line = v.elementAt(j).trim().split("\t", -1);
	        	names[i][j] = line[0];
	        	atk[i][j] = Integer.parseInt(line[1]);
	        	def[i][j] = Integer.parseInt(line[2]);
	        	upkeep[i][j] = Integer.parseInt(line[3]);
            }
        }
        
        limits = new int[THRESHOLDS.length];
        for (int i = 0; i<THRESHOLDS.length; i++) {
        	limits[i] = (int)(DOLLAR/CREW*THRESHOLDS[i]);
        }
        bestSets = new int[THRESHOLDS.length][3][types.length];
        bestSetValues = new int[THRESHOLDS.length][3][4];
        counters = Array.intArray(names.length, 0);
        time = new Date().getTime();
        while (counters[0] < names[0].length) {
        	if (counters[0] == 21 && counters[1] == 18 && counters[2] == 23) {
        		System.out.println(names[0][counters[0]]+"/"+names[1][counters[1]]+"/"+names[2][counters[2]]);
        	}
    		sumUpkeep = 0;
    		for (int j = 0; j<counters.length; j++) {
    			sumUpkeep += upkeep[j][counters[j]];
            }
        	for (int i = 0; i<limits.length; i++) {
        		if (sumUpkeep <= limits[i]) {
        			sums = new int[3];
	        		for (int j = 0; j<counters.length; j++) {
	        			sums[0] += atk[j][counters[j]];
	        			sums[1] += def[j][counters[j]];
	        			sums[2] += atk[j][counters[j]]+def[j][counters[j]];
                    }
	            	if (counters[0] == 21 && counters[1] == 18 && counters[2] == 23) {
//	            		System.out.println(sums[0]+"/"+sums[1]+"/"+sumUpkeep);
	            	}
        			for (int k = 0; k<3; k++) {
		        		if (sums[k] > bestSetValues[i][k][0] || (sums[k] == bestSetValues[i][k][0] && sumUpkeep < bestSetValues[i][k][3])) {
		        			bestSetValues[i][k][0] = sums[k];
		        			bestSetValues[i][k][1] = sums[0];
		        			bestSetValues[i][k][2] = sums[1];
		        			bestSetValues[i][k][3] = sumUpkeep;
		        			for (int j = 0; j<counters.length; j++) {
		        				bestSets[i][k][j] = counters[j];
                            }
		        		}      			
        			}
        		}
            }

        	counters[types.length-1]++;
        	for (int i = counters.length-1; i>0; i--) {
        		if (counters[i] == names[i].length) {
        			counters[i] = 0;
        			counters[i-1]++;
        		}
            }
        }
        
        try {
	        writer = new PrintWriter(new FileWriter(ext.rootOf(filename)+"_optimal.out"));
	        for (int i = 0; i<limits.length; i++) {
	        	System.out.println("At "+(int)(THRESHOLDS[i]*100)+"% ("+limits[i]+")");
	        	for (int k = 0; k<3; k++) {
	        		if (k==0) {
	        			System.out.print("Max ATK: ");	
	        		} else if (k==1) {
	        			System.out.print("Max DEF: ");	
	        		} else {
	        			System.out.print("Max PTS: ");	
	        		}
		        	System.out.print(ext.formStr(bestSetValues[i][k][1]*CREW+"", 5)+" atk "+ext.formStr(bestSetValues[i][k][2]*CREW+"", 5)+" def "+ext.formStr(bestSetValues[i][k][3]+"", 5)+" ("+ext.formDeci(bestSetValues[i][k][3]/((double)DOLLAR/CREW)*100, 1)+"%) upkeep each using ");
		        	for (int j = 0; j<types.length; j++) {
		        		System.out.print((j==0?"":"/")+names[j][bestSets[i][k][j]]+"("+atk[j][bestSets[i][k][j]]+","+def[j][bestSets[i][k][j]]+")");
	                }
		        	System.out.println();
                }
	        	System.out.println();
            }
	        writer.close();
        } catch (Exception e) {
	        System.err.println("Error writing to "+ext.rootOf(filename)+"_optimal.out");
	        e.printStackTrace();
        }
        
        
//        String[] test = {"Hardened Tonfas of Earth", "Windmill Block", "Rebirth"};
//      String[] test = {"Hardened Tonfas of Earth", "Shadow Strike", "Rebirth"};
      String[] test = {"Gemini Sais", "Shadow Strike", "Rebirth"};
//        String[] test = {"Hardened Tonfas of Earth", "Shadow Strike", "Stellar Storm"};
//        String[] test = {"Star Sword of Truth", "Shadow Strike", "Stellar Storm"};
//        String[] test = {"Doom Scythe", "Shadow Strike", "Stellar Storm"};
//        String[] test = {"Striking Wind", "Mana Potion", "Magic Ship"};
        counters = new int[3];
		sums = new int[3];
		sumUpkeep = 0;
        for (int i = 0; i<3; i++) {
        	counters[i] = ext.indexOfStr(test[i], names[i]);
        	if (counters[i] == -1) {
        		System.err.println("Error - could not find '"+test[i]+"'");
        	}
			sums[0] += atk[i][counters[i]];
			sums[1] += def[i][counters[i]];
			sumUpkeep += upkeep[i][counters[i]];
        }
        System.out.print("Check  : ");
    	System.out.print(ext.formStr(sums[0]*CREW+"", 5)+" atk "+ext.formStr(sums[1]*CREW+"", 5)+" def "+ext.formStr(sumUpkeep+"", 5)+" ("+ext.formDeci(sumUpkeep/((double)DOLLAR/CREW)*100, 1)+"%) upkeep each using ");
    	for (int j = 0; j<types.length; j++) {
    		System.out.print((j==0?"":"/")+names[j][counters[j]]+"("+atk[j][counters[j]]+","+def[j][counters[j]]+")");
        }
    	System.out.println();
    	System.out.println();
        System.out.println("Finished in "+ext.getTimeElapsed(time));
        
	}

	public static void main(String[] args) {
	    int numArgs = args.length;
	    String filename = DEFAULT_FILENAME;

	    String usage = "\n"+
	    "eglence.OptimalSet requires 0-1 arguments\n"+
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
		    optimize(filename);
	    } catch (Exception e) {
		    e.printStackTrace();
	    }
    }
	
}
