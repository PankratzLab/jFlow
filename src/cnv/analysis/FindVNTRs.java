// -Xms1024M -Xmx1024M
package cnv.analysis;

import java.io.*;
import java.util.*;

//import stats.*;
//import cnv.filesys.*;
import common.*;

public class FindVNTRs {
	public static final String[] EXPECTED_HEADER = {"Marker", "Chr", "Position"};
	
	public static void find(String filename) {
		BufferedReader reader;
        PrintWriter writer;
        String[] line, keys, values;
        int[] counts;
        String trav;
        Hashtable<String,Vector<String>> hash;
        Vector<String> v;
        CountVector cv;
        
        hash = new Hashtable<String,Vector<String>>();
        try {
	        reader = new BufferedReader(new FileReader(filename));
	        line = reader.readLine().trim().split("[\\s]+");
	        ext.checkHeader(line, EXPECTED_HEADER, true);
	        while (reader.ready()) {
	        	line = reader.readLine().trim().split("[\\s]+");
	        	if (line[0].toLowerCase().startsWith("cnv")) {
	        		trav = line[0].toLowerCase().substring(3);
	        		if (trav.indexOf("p") == -1) {
	        			System.out.println(line[0]+" does not have a P");
	        		} else {
	        			trav = trav.substring(0, trav.indexOf("p"));
	        			try {
	        				Integer.parseInt(trav);
//		        			HashVec.addToHashVec(hash, trav, Array.toStr(line), false);
		        			HashVec.addToHashVec(hash, trav, line[0], false);
	        			} catch (NumberFormatException nfe) {
	        				System.out.println("More than a number going on between cnv and P: "+line[0]);
	        			}
	        		}
	        	}
	        	
	        }
	        reader.close();
        } catch (FileNotFoundException fnfe) {
	        System.err.println("Error: file \""+filename+"\" not found in current directory");
	        System.exit(1);
        } catch (IOException ioe) {
	        System.err.println("Error reading file \""+filename+"\"");
	        System.exit(2);
        }
        
        try {
	        writer = new PrintWriter(new FileWriter(filename+"_vntrs.xln"));
	        keys = HashVec.getKeys(hash);
	        cv = new CountVector();
	        for (int i = 0; i<keys.length; i++) {
	        	v = hash.get(keys[i]);
	        	cv.add(v.size()+"");
	        	if (v.size() == 3) {
//	        		for (int j = 0; j<v.size(); j++) {
//	        			writer.println(v.elementAt(j));
//                    }
	        		writer.println(Array.toStr(Array.toStringArray(v)));
	        	}
            }
	        cv.sort(true);
	        counts = cv.getCounts();
	        values = cv.getValues();
	        for (int i = 0; i<values.length; i++) {
	        	System.out.println(values[i]+" (n="+counts[i]+")");
            }
	        writer.close();
        } catch (Exception e) {
	        System.err.println("Error writing to "+filename+"_vntrs.xln");
	        e.printStackTrace();
        }
	}
	
//	public static void workup(Project proj, String filename, String pedfile) {
//		BufferedReader reader;
//        PrintWriter writer;
//        String[] line, files;
//        String trav;
//        Hashtable<String,String> hash;
//        int count, numSamples;
//        MarkerDataCollection collection;
//        MarkerLookup markerLookup;
//        String[] famIndPairs, affStats;
//        SampleData sampleData;
//        String[] samples;
//        boolean[] use;
//        MarkerData[] markerData;
//        float[] lrrs;
//        float[][] data;
//        double[] means;
//        double[] affs;
//        DoubleVector corrs;
//        LogisticRegression reg;
//        
//		files = Files.list(proj.getDir(Project.PLOT_DIRECTORY), ".scat", false);
//		if (files.length != 1) {
//			System.err.println("Error - found "+files.length+" .scat files ");
//			return;
//		}
//
//		sampleData = proj.getSampleData(2, false);
//		famIndPairs = HashVec.loadFileToStringArray(pedfile, false, false, new int[] {0, 1}, false);
//		affStats = HashVec.loadFileToStringArray(pedfile, false, false, new int[] {5}, false);
//		hash = new Hashtable<String,String>();
//		for (int i = 0; i<famIndPairs.length; i++) {
//			hash.put(sampleData.lookup(famIndPairs[i])[0], affStats[i]);
//        }
//		samples = proj.getSamples();
//		use = new boolean[samples.length];
//		for (int i = 0; i<samples.length; i++) {
//			use[i] = hash.containsKey(samples[i]) && Integer.parseInt(hash.get(samples[i])) > 0;
//        }
//		numSamples = Array.booleanArraySum(use);
//		System.out.println("Able to use "+numSamples+" of the "+hash.size()+" individuals found in "+ext.removeDirectoryInfo(pedfile)+" (only using indiviudals with non-missing phenotype and non-missing array data");
//		affs = new double[numSamples];
//		count = 0;
//		for (int i = 0; i<samples.length; i++) {
//			if (use[i]) {
//				affs[count] = Double.parseDouble(hash.get(samples[i]))-1;
//				count++;
//			}
//        }
//        
//        markerLookup = proj.getMarkerLookup();
//        collection = MarkerDataCollection.load(proj.getDir(Project.PLOT_DIRECTORY)+files[0], false);
//        markerData = collection.getCollection();
//        
//        try {
//	        reader = new BufferedReader(new FileReader(filename));
//	        writer = new PrintWriter(new FileWriter(ext.rootOf(filename, false)+"_corr.xln"));
//	        writer.println("N_markers\tMeanCorr\tAffMean\tUnaffMean\tpval\tMarkers&Corrs");
//	        while (reader.ready()) {
//	        	line = reader.readLine().trim().split("[\\s]+");
//	        	data = new float[line.length][];
//	        	for (int i = 0; i<line.length; i++) {
//	        		trav = markerLookup.get(line[i]);
//	        		if (trav == null) {
//	        			System.err.println("Error - '"+line[i]+"' not present in the dataset (at least in the MarkerLookup)");
//	        		} else {
//	        			lrrs = markerData[Integer.parseInt(trav.split("[\\s]+")[1])].getLRRs();
//	        			data[i] = new float[numSamples];
//	        			count = 0;
//	        			for (int j = 0; j<lrrs.length; j++) {
//	        				if (use[j]) {
//	        					data[i][count] = lrrs[j];
//	        					count++;
//	        				}
//                        }
//	        		}
//                }
//	        	writer.print(line.length);
//	        	corrs = new DoubleVector();
//	        	for (int i = 0; i<data.length-1; i++) {
//	        		for (int j = i+1; j<data.length; j++) {
//	        			corrs.add(Correlation.Pearson(Matrix.removeRowsWithNaN(new double[][] {Array.toDoubleArray(data[i]), Array.toDoubleArray(data[j])}))[0]);
//                    }
//                }
//	        	writer.print("\t"+Array.mean(corrs.toArray()));
//	        	
//	        	means = new double[numSamples];
//	        	for (int i = 0; i<numSamples; i++) {
//	        		count = 0;
//	        		for (int j = 0; j<data.length; j++) {
//	        			if (!Float.isNaN(data[j][i])) {
//	        				means[i] += data[j][i];
//	        				count++;
//	        			}
//                    }
//	        		means[i] /= count;
//                }
//	        	writer.print("\t"+Array.meanIf(means, affs, 1.0, true)+"\t"+Array.meanIf(means, affs, 0.0, true));
//	        	reg = new LogisticRegression(affs, means);
//	        	writer.print("\t"+reg.getSigs()[1]);
//	        	
//	        	writer.println("\t"+Array.toStr(line)+"\t"+Array.toStr(corrs.toArray()));
//	        }
//	        reader.close();
//            writer.close();
//        } catch (FileNotFoundException fnfe) {
//	        System.err.println("Error: file \""+filename+"\" not found in current directory");
//	        System.exit(1);
//        } catch (IOException ioe) {
//	        System.err.println("Error reading file \""+filename+"\"");
//	        System.exit(2);
//        }
//	}
	
	public static void main(String[] args) {
	    int numArgs = args.length;
	    String filename = "C:\\Documents and Settings\\npankrat\\My Documents\\CNV_PD\\cnvPositions.txt";
//	    String workup = "C:\\Documents and Settings\\npankrat\\My Documents\\CNV_PD\\cnvPositions.txt_vntrs.xln";
//	    String workup = "C:\\Documents and Settings\\npankrat\\My Documents\\CNV_PD\\rand1000.txt";
//	    String workup = "C:\\Documents and Settings\\npankrat\\My Documents\\CNV_PD\\onePer1730.txt";
//	    String proj = Project.DEFAULT_SCATTER_PROJECT;
//	    String pedfile = "C:\\Documents and Settings\\npankrat\\My Documents\\CNV_PD\\conf1672.fam";

	    String usage = "\n"+
	    "cnv.analysis.FindVNTRs requires 0-1 arguments\n"+
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
//	    	if (!workup.equals("")) {
//	    		workup(new Project(proj, false), workup, pedfile);
//	    	} else {
	    		find(filename);
//	    	}
	    } catch (Exception e) {
		    e.printStackTrace();
	    }
    }
}
