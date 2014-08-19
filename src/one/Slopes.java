package one;

import java.io.*;
import java.util.*;

import parse.GenParser;
import stats.LeastSquares;
import common.*;

public class Slopes {

	public static void parse(String dir, String filename) {
		BufferedReader reader;
        PrintWriter writer;
        String[] line;
        Hashtable<String,String> hash = new Hashtable<String,String>();
        int timepoints;
        Vector<String> allValues, preValues, postValues;
        Vector<double[]> allDates, preDates, postDates;
        double date, value, refDate;
        LeastSquares linear;
        
        try {
	        reader = new BufferedReader(new FileReader(dir+filename));
	        writer = new PrintWriter(new FileWriter(dir+ext.rootOf(filename)+"_slopes.xln"));
	        writer.println("MRN\tSlopeAll\tN\tSlopePre\tN\tSlopePost\tN");
	        line = reader.readLine().trim().split("[\\s]+");
	        if (!line[0].equals("MRN")) {
	        	System.err.println("Error - expecting header with MRN");
	        }
	        while (reader.ready()) {
	        	line = reader.readLine().trim().split("[\\s]+");
	        	refDate = parseDate(line[1]);
	        	timepoints = (line.length-2)/2;
	        	allValues = new Vector<String>();
	        	allDates = new Vector<double[]>();
	        	preValues = new Vector<String>();
	        	preDates = new Vector<double[]>();
	        	postValues = new Vector<String>();
	        	postDates = new Vector<double[]>();
	        	for (int i = 0; i<timepoints; i++) {
	        		if (line[2+i*2+0].equals(".") && line[2+i*2+1].equals(".")) {
	        		} else if (line[2+i*2+0].equals(".")) {
	        			System.err.println("Error - missing value for time point #"+(i+1)+" for "+line[0]);
	        		} else if (line[2+i*2+1].equals(".")) {
	        			System.err.println("Error - missing date for time point #"+(i+1)+" for "+line[0]);
	        		} else {
	        			value = Double.parseDouble(line[2+i*2+0]);
	        			date = parseDate(line[2+i*2+1]);
	        			if (hash.containsKey(date+"")) {
//	        				System.err.println("Error - duplicate dates ("+line[2+i*2+1]+") for "+line[0]);
	        			} else {
		        			allValues.add(value+"");
		        			allDates.add(new double[] {date});
		        			hash.put(date+"", "");
		        			if (date < refDate) {
		        				preValues.add(value+"");
		        				preDates.add(new double[] {date});
		        			}
		        			if (date > refDate - 1) {
		        				postValues.add(value+"");
		        				postDates.add(new double[] {date});
		        			}
	        			}
	        		}
                }
	        	writer.print(line[0]);
	        	if (allValues.size() >= 2) {
	        		linear = new LeastSquares(allValues, allDates);
	        		if (linear.analysisFailed()) {
	        			System.err.println("Error - problem with 'all' data for "+line[0]);
	        		}
	        		writer.print("\t"+linear.getBetas()[1]+"\t"+allValues.size());
	        	} else {
	        		writer.print("\t.\t"+allValues.size());
	        	}
	        	if (preValues.size() >= 2) {
	        		linear = new LeastSquares(preValues, preDates);
	        		if (linear.analysisFailed()) {
	        			System.err.println("Error - problem with 'pre' data for "+line[0]);
	        		}
	        		writer.print("\t"+linear.getBetas()[1]+"\t"+preValues.size());
	        	} else {
	        		writer.print("\t.\t"+preValues.size());
	        	}
	        	if (postValues.size() >= 2) {
	        		linear = new LeastSquares(postValues, postDates);
	        		if (linear.analysisFailed()) {
	        			System.err.println("Error - problem with 'post' data for "+line[0]);
	        		}
	        		writer.print("\t"+linear.getBetas()[1]+"\t"+postValues.size());
	        	} else {
	        		writer.print("\t.\t"+postValues.size());
	        	}
	        	writer.println();
	        }
	        reader.close();
            writer.close();
        } catch (FileNotFoundException fnfe) {
	        System.err.println("Error: file \""+dir+filename+"\" not found in current directory");
	        System.exit(1);
        } catch (IOException ioe) {
	        System.err.println("Error reading file \""+dir+filename+"\"");
	        System.exit(2);
        }
	}
	
	public static double parseDate(String str) {
		String[] line;
		int month, date, year;
		
		line = str.trim().split("/");
		if (line.length != 3) {
			System.err.println("Error - invalid date ("+str+")");
		}
		try {
			month = Integer.parseInt(line[0]);
			date = Integer.parseInt(line[1]);
			year = Integer.parseInt(line[2]);
		} catch (NumberFormatException nfe) {
			nfe.printStackTrace();
			return Double.NaN;
		}
		
		return new GregorianCalendar(year, month, date).getTimeInMillis()/1000/60/60/24/365.25*12;
	}
	
	public static void parseListwise(String dir, String filename, int idIndex, String dateIndices, String valueIndices, String refFile) {
        PrintWriter writer;
        String[] line;
        Hashtable<String,String> refs, dateCheck;
        Hashtable<String,Vector<String>> hash;
        DoubleVector[] values, dates;
        Vector<String> v;
        double date, value, refDate;
        LeastSquares linear;
        String trait;
        String[] ids;
        int[] dateIs, valueIs, allIs;

        if (refFile == null) {
        	refs = new Hashtable<String,String>();
        } else {
        	line = Files.getHeaderOfFile(dir+refFile, "\t", new Logger());
            ext.checkHeader(line, new String[] {"id", "date"}, false);
            refs = HashVec.loadFileToHashString(dir+refFile, true);
        }
        
        dateIs = Array.toIntArray(dateIndices.trim().split(","));
        valueIs = Array.toIntArray(valueIndices.trim().split(","));
        if (dateIs.length != valueIs.length) {
        	System.err.println("Error - number of value indices and date indices does not match up");
        	return;
        }
        allIs = new int[dateIs.length*2];
        for (int i = 0; i<dateIs.length; i++) {
        	allIs[i*2+0] = dateIs[i];
        	allIs[i*2+1] = valueIs[i];
        }
        
    	line = Files.getHeaderOfFile(dir+filename, "\t", new Logger());
    	trait = line[valueIs[0]];
    	if (valueIs.length > 1) {
    		trait = trait.substring(0, trait.lastIndexOf("_"));
    	}
        ext.checkHeader(line, new String[] {"id"}, new int[] {0}, false, new Logger(), true);
        hash = HashVec.loadFileToHashVec(dir+filename, idIndex, allIs, "\t", true, false);
        
        ids = HashVec.getKeys(hash, true, true);
        try {
	        writer = new PrintWriter(new FileWriter(dir+trait+"_slopes"+(refFile==null?"":"_after_"+ext.rootOf(refFile))+".xln"));
	        writer.println("id\t"+trait+"_slope\tN"+(refFile == null?"":"\t"+trait+"_slope_pre_"+ext.rootOf(refFile)+"\tN\t"+trait+"_slope_post_"+ext.rootOf(refFile)+"\tN"));
	        for (int i = 0; i<ids.length; i++) {
	        	v = hash.get(ids[i]);
	        	if (refs.containsKey(ids[i])) {
	        		refDate = Double.parseDouble(refs.get(ids[i]));
	        	} else {
	        		refDate = Double.POSITIVE_INFINITY;
	        	}
	        	values = DoubleVector.newDoubleVectors(3);
	        	dates = DoubleVector.newDoubleVectors(3);
	        	dateCheck = new Hashtable<String,String>();
	        	for (int j = 0; j<v.size(); j++) {
		        	line = v.elementAt(j).split("[\\s]+");
		        	for (int k = 0; k<line.length/2; k++) {
			        	date = GenParser.procDouble(line[k*2+0]);
			        	value = GenParser.procDouble(line[k*2+1]);
			        	
			        	if (!Double.isNaN(date) && !Double.isNaN(value)) {
		        			if (dateCheck.containsKey(date+"")) {
		        				if (!dateCheck.get(date+"").equals(value+"")) {
			        				System.err.println("Error - duplicate dates ("+line[0]+") for "+ids[i]+" and yet different values: "+dateCheck.get(date+"")+" and "+value);
		        				}
		        			} else {
			        			values[0].add(value);
			        			dates[0].add(date);
			        			dateCheck.put(date+"", value+"");
			        			if (date < refDate) {
				        			values[1].add(value);
				        			dates[1].add(date);
			        			}
			        			if (date > refDate - 30) {
				        			values[2].add(value);
				        			dates[2].add(date);
			        			}
		        			}
		        		}
                    }
                }
	        	writer.print(ids[i]);
	        	for (int j = 0; j<(refFile == null?1:3); j++) {
		        	if (values[j].size() >= 2) {
		        		if (Array.variance(values[j].toArray()) == 0) {
		        			System.err.println("Warning - no variation for '"+(j==0?"all":(j==1?"pre":"post"))+"' data for "+ids[i]+"; slope is zero");
			        		writer.print("\t0");
		        		} else {
			        		linear = new LeastSquares(values[j].toArray(), dates[j].toArray());
			        		if (linear.analysisFailed()) {
			        			System.err.println("Error - problem with '"+(j==0?"all":(j==1?"pre":"post"))+"' data for "+ids[i]);
			        		}
			        		writer.print("\t"+linear.getBetas()[1]);
		        		}
		        	} else {
		        		writer.print("\t.");
		        	}
	        		writer.print("\t"+values[j].size());
                }
	        	writer.println();
	        }
            writer.close();
        } catch (FileNotFoundException fnfe) {
	        System.err.println("Error: file \""+dir+filename+"\" not found in current directory");
	        System.exit(1);
        } catch (IOException ioe) {
	        System.err.println("Error reading file \""+dir+filename+"\"");
	        System.exit(2);
        }
	}
	
	public static void checkIfBefore(String dir, String filename, int idIndex, String dateIndices, String valueIndices, String refFile) {
        PrintWriter writer;
        String[] line;
        Hashtable<String,String> refs;
        Hashtable<String,Vector<String>> hash;
        Vector<String> v;
        double date, value, refDate;
        String trait;
        String[] ids;
        int[] dateIs, valueIs, allIs;
        double firstEver, firstBefore, firstAfter;

        if (refFile == null) {
        	refs = new Hashtable<String,String>();
        } else {
        	line = Files.getHeaderOfFile(dir+refFile, "\t", new Logger());
            ext.checkHeader(line, new String[] {"id", "date"}, false);
            refs = HashVec.loadFileToHashString(dir+refFile, true);
        }
        
        dateIs = Array.toIntArray(dateIndices.trim().split(","));
        valueIs = Array.toIntArray(valueIndices.trim().split(","));
        if (dateIs.length != valueIs.length) {
        	System.err.println("Error - number of value indices and date indices does not match up");
        	return;
        }
        allIs = new int[dateIs.length*2];
        for (int i = 0; i<dateIs.length; i++) {
        	allIs[i*2+0] = dateIs[i];
        	allIs[i*2+1] = valueIs[i];
        }
        
    	line = Files.getHeaderOfFile(dir+filename, "\t", new Logger());
    	trait = line[valueIs[0]];
    	if (valueIs.length > 1) {
    		trait = trait.substring(0, trait.lastIndexOf("_"));
    	}
        ext.checkHeader(line, new String[] {"id"}, new int[] {0}, false, new Logger(), true);
        hash = HashVec.loadFileToHashVec(dir+filename, idIndex, allIs, "\t", true, false);
        
        ids = HashVec.getKeys(hash, true, true);
        try {
	        writer = new PrintWriter(new FileWriter(dir+trait+"_first"+(refFile==null?"":"Before_"+ext.rootOf(refFile))+".xln"));
	        writer.println("id\t"+trait+"_ever\t"+trait+"_firstEver\t"+trait+"_before_"+ext.rootOf(refFile)+"\t"+trait+"_firstBefore_"+ext.rootOf(refFile)+"\t"+trait+"_after_"+ext.rootOf(refFile)+"\t"+trait+"_firstAfter_"+ext.rootOf(refFile));
	        for (int i = 0; i<ids.length; i++) {
	        	v = hash.get(ids[i]);
	        	if (refs.containsKey(ids[i])) {
	        		refDate = Double.parseDouble(refs.get(ids[i]));
	        	} else {
	        		refDate = Double.POSITIVE_INFINITY;
	        	}

	        	firstEver = Double.POSITIVE_INFINITY;
	        	firstBefore = Double.POSITIVE_INFINITY;
	        	firstAfter = Double.POSITIVE_INFINITY;
	        	
	        	for (int j = 0; j<v.size(); j++) {
		        	line = v.elementAt(j).split("[\\s]+");
		        	for (int k = 0; k<line.length/2; k++) {
			        	date = GenParser.procDouble(line[k*2+0]);
			        	value = GenParser.procDouble(line[k*2+1]);
			        	
			        	if (!Double.isNaN(date) && !Double.isNaN(value)) {
			        		if (value > 0) {
			        			firstEver = Math.min(date, firstEver);
			        			if (date < refDate) {
			        				firstBefore = Math.min(date, firstBefore);
			        			}
			        			if (date >= refDate) {
			        				firstAfter = Math.min(date, firstAfter);
			        			}
			        		}
		        		}
                    }
                }
	        	writer.println(ids[i]+"\t"+(firstEver == Double.POSITIVE_INFINITY?"0\t.":"1\t"+firstEver)+"\t"+(firstBefore == Double.POSITIVE_INFINITY?"0\t.":"1\t"+firstBefore)+"\t"+(firstAfter == Double.POSITIVE_INFINITY?"0\t.":"1\t"+firstAfter));
	        }
            writer.close();
        } catch (FileNotFoundException fnfe) {
	        System.err.println("Error: file \""+dir+filename+"\" not found in current directory");
	        System.exit(1);
        } catch (IOException ioe) {
	        System.err.println("Error reading file \""+dir+filename+"\"");
	        System.exit(2);
        }
	}

	public static void main(String[] args) {
	    int numArgs = args.length;
//	    String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\ALS\\Sleep_Sheetal\\";
//	    String filename = "ALSFRS.txt";
//	    String filename = "FVC.txt";
//	    String filename = "Bulbar.txt";
//	    String filename = "Resp.txt";
	    boolean list = true;
	    String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\ALS\\NEALS_db\\Celebrex\\";
//	    String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\ALS\\NEALS_db\\Topiramate\\";
	    String refFile = "SeriousPulmonaryEvent.dat";
	    int idIndex = 0;
//	    String filename = "alsfrs.dat";
//	    String dateIndices = "1";
//	    String valueIndices = "5";
//	    String filename = "BMI_calcs2.xln";
//	    String dateIndices = "3,5,7,9,11,13,15,17,19";
//	    String valueIndices ="2,4,6,8,10,12,14,16,18";
//	    String dateIndices = "3,5,7,9,11,13,15,17";
//	    String valueIndices ="2,4,6,8,10,12,14,16";
	    String filename = "fvc_fixed.dat";
	    String dateIndices = "1";
//	    String valueIndices = "2";
	    String ifBefore = "respiratory_rate.dat";
	    String valueIndices = "3";
		
	    String usage = "\n"+"one.Slopes requires 0-1 arguments\n"+"   (1) filename (i.e. file="+filename+" (default))\n"+"";

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
	    	if (ifBefore != null && !ifBefore.equals("")) {
	    		checkIfBefore(dir, ifBefore, idIndex, dateIndices, valueIndices, refFile);
	    	} else if (list) {
	    		parseListwise(dir, filename, idIndex, dateIndices, valueIndices, refFile);
	    	} else {
	    		parse(dir, filename);
	    	}
	    } catch (Exception e) {
		    e.printStackTrace();
	    }
    }
}
