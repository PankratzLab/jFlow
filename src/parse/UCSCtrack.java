package parse;

import java.io.*;
import java.util.*;

import stats.Maths;
import common.*;

public class UCSCtrack {
	public static final int SIG_FIGS = 2;

	public static void describeFromParameters(String filename, Logger log) {
		BufferedReader reader;
		PrintWriter writer;
        String resultsFile, mapFile;
        String[] line;
        int markerIndex, valueIndex, chrIndex, startIndex, stopIndex;
        boolean ignoreCase, ignoreFirstLine, commaDelimited, tabDelimited;
        Hashtable<String,Double> hash;
        String outfile;  
        int start, stop, prevStart;
        double value;
        boolean ignoreWhenAbsentInResultsFile, separateMap, neg_log;
        Vector<String> paramV;

    	String[] filterValues, filterOps;
    	int[] filterCols;
    	Vector<String> filters;
    	boolean passedFilter;
        
		paramV = Files.parseControlFile(filename, "ucsc", new String[] {"plink.assoc 1 8 neg_log chr=0 start=2 stop=2 out=assoc.bed", "#", "#results.xln tab noHeader 0 1 neg_log", "#plink.bim 1 0 start=3 stop=3 ignoreWhenAbsentInResultsFile", "#", "#results.csv , ignorecase 0 1", "#markers.csv , header 1 chr=0 start=3 stop=3"}, log);
		if (paramV == null) {
			log.reportError("Failed...");
			return;
		}

		if (paramV.size() > 2) {
        	log.reportError("Error - can only have one or two lines of paramters (results or results+map), not "+paramV.size()+" lines");
        	return;
        } else if (paramV.size() == 2) {
        	separateMap = true;
        } else if (paramV.size() == 1) {
        	separateMap = false;
        } else {
        	log.reportError("Difficult to parse without any parameters...");
        	return;
        }
        
    	line = paramV.elementAt(0).trim().split("[\\s]+");
    	resultsFile = line[0];
    	markerIndex = valueIndex = chrIndex = startIndex = stopIndex = -1;
    	neg_log = false;
    	ignoreCase = false;
    	commaDelimited = false;
    	tabDelimited = false;
    	ignoreFirstLine = true;
    	filters = new Vector<String>();
    	outfile = ext.rootOf(resultsFile)+".bed";
    	for (int j = 1; j<line.length; j++) {
    		if (line[j].equalsIgnoreCase("ignoreCase")) {
    			ignoreCase = true;
    		} else if (line[j].equalsIgnoreCase("noHeader")) {
    			ignoreFirstLine = false;
    		} else if (line[j].startsWith("out=")) {
    			outfile = line[j].split("=")[1];
    		} else if (line[j].startsWith("chr=")) {
    			chrIndex = Integer.parseInt(line[j].split("=")[1]);
    		} else if (line[j].startsWith("start=")) {
    			startIndex = Integer.parseInt(line[j].split("=")[1]);
    		} else if (line[j].startsWith("stop=")) {
    			stopIndex = Integer.parseInt(line[j].split("=")[1]);
    		} else if (line[j].equals(",")) {
    			commaDelimited = true;
    		} else if (line[j].equals("tab")) {
    			tabDelimited = true;
    		} else if (line[j].equals("neg_log")) {
    			neg_log = true;
    		} else if (line[j].startsWith("!")) {
    			filters.add(line[j].substring(1));
    		} else if (markerIndex == -1){
    			markerIndex = Integer.parseInt(line[j]);
    		} else if (valueIndex == -1){
    			valueIndex = Integer.parseInt(line[j]);
    		} else {
    			log.reportError("Error - don't know what to do with argument: "+line[j]);
    		}
        }
    	
    	if (!separateMap && (chrIndex == -1 || startIndex == -1)) {
    		log.reportError("Error - if a separate mapfile is not provided, then the chr index and start index must be provided");
    	}
    	if (stopIndex == -1) {
    		stopIndex = startIndex;
    	}
    	
    	filterCols = new int[filters.size()];
    	filterOps = new String[filters.size()];
    	filterValues = new String[filters.size()];


    	for (int i = 0; i<filters.size(); i++) {
    		line = Maths.parseOp(filters.elementAt(i));
    		if (line == null) {
    			log.reportError("Error - invalid token in filter('"+filters.elementAt(i)+"')");
    			return;
    		}
    		filterOps[i] = line[1];
			filterCols[i] = Integer.parseInt(line[0]);
			filterValues[i] = line[2];
        }

    	hash = new Hashtable<String,Double>();
    	try {
	        writer = new PrintWriter(new FileWriter(outfile));
			writer.println("track type=wiggle_0 name=\""+ext.rootOf(outfile)+"\"");
	        try {
		        reader = new BufferedReader(new FileReader(resultsFile));
		        if (ignoreFirstLine) {
		        	reader.readLine();
		        }
		        prevStart = -1;
				while (reader.ready()) {
					line = reader.readLine().trim().split(commaDelimited?",":(tabDelimited?"\t":"[\\s]+"));
			    	passedFilter = true;
			    	for (int i = 0; i<filterOps.length && passedFilter; i++) {
			    		if (filterOps[i].equals("=")) {
			    			if (!line[filterCols[i]].equals(filterValues[i])) {
			    				passedFilter = false;
			        		}
			    		} else if (filterOps[i].equals("!")) {
			    			if (line[filterCols[i]].equals(filterValues[i])) {
			    				passedFilter = false;
			    			}
			    		} else if (ext.isMissingValue(line[filterCols[i]]) || !Maths.op(Double.parseDouble(line[filterCols[i]]), Double.parseDouble(filterValues[i]), filterOps[i])) {
			    			passedFilter = false;
			    		}
			        }
			    	if (passedFilter) {
						if (line[valueIndex].equals(".")||line[valueIndex].toUpperCase().startsWith("NA")) {
							value = -1;
						} else {
							value = Double.parseDouble(line[valueIndex]);
							if (neg_log) {
								value = -1*Math.log10(value);
							}
						}
						if (separateMap) {
							hash.put(ignoreCase?line[markerIndex].toLowerCase():line[markerIndex], Double.valueOf(value));
						} else {
							start = Integer.parseInt(line[startIndex]);
							stop = Integer.parseInt(line[stopIndex]);
							if (start == stop) {
								stop += 1;
							}
							if (start == prevStart) {
								start++;
								stop++;
							}
							writer.println("chr"+Positions.chromosomeNumberInverse(Integer.parseInt(line[chrIndex]))+"\t"+start+"\t"+stop+"\t"+ext.formDeci(value, SIG_FIGS));
							prevStart = start;
						}
			    	}
				}
		        
		        reader.close();
	        } catch (FileNotFoundException fnfe) {
	        	log.reportError("Error: file \""+resultsFile+"\" not found in current directory");
	        	log.reportException(fnfe);
		        return;
	        } catch (IOException ioe) {
	        	log.reportError("Error reading file \""+resultsFile+"\"");
	        	log.reportException(ioe);
		        return;
	        }
	        
	        if (separateMap) {
	        	line = paramV.elementAt(1).trim().split("[\\s]+");
	        	mapFile = line[0];
	        	markerIndex = chrIndex = startIndex = stopIndex = -1;
	        	ignoreCase = false;
	        	commaDelimited = false;
	        	tabDelimited = false;
	        	ignoreFirstLine = false;
	        	ignoreWhenAbsentInResultsFile = false;
	        	outfile = ext.rootOf(resultsFile)+"_described.xln";
	        	for (int j = 1; j<line.length; j++) {
	        		if (line[j].equalsIgnoreCase("header")) {
	        			ignoreFirstLine = true;
	        		} else if (line[j].startsWith("chr=")) {
	        			chrIndex = Integer.parseInt(line[j].split("=")[1]);
	        		} else if (line[j].startsWith("start=")) {
	        			startIndex = Integer.parseInt(line[j].split("=")[1]);
	        		} else if (line[j].startsWith("stop=")) {
	        			stopIndex = Integer.parseInt(line[j].split("=")[1]);
	        		} else if (line[j].equals(",")) {
	        			commaDelimited = true;
	        		} else if (line[j].equals("tab")) {
	        			tabDelimited = true;
	        		} else if (line[j].equals("ignoreWhenAbsentInResultsFile")) {
	        			ignoreWhenAbsentInResultsFile = true;
	        		} else if (markerIndex == -1){
	        			markerIndex = Integer.parseInt(line[j]);
	        		} else {
	        			log.reportError("Error - don't know what to do with argument: "+line[j]);
	        		}
	            }
	        	
	        	if (stopIndex == -1) {
	        		stopIndex = startIndex;
	        	}	        	
	        	
		        try {
			        reader = new BufferedReader(new FileReader(mapFile));
			        if (ignoreFirstLine) {
			        	reader.readLine();
			        }
			        prevStart = -1;
					while (reader.ready()) {
						line = reader.readLine().trim().split(commaDelimited?",":(tabDelimited?"\t":"[\\s]+"));
						start = Integer.parseInt(line[startIndex]);
						stop = Integer.parseInt(line[stopIndex]);
						if (start == stop) {
							stop += 1;
						}
						if (start == prevStart) {
							start++;
							stop++;
						}
						if (hash.containsKey(ignoreCase?line[markerIndex].toLowerCase():line[markerIndex])) {
							value = hash.get(ignoreCase?line[markerIndex].toLowerCase():line[markerIndex]).doubleValue();
							writer.println("chr"+Positions.chromosomeNumberInverse(Integer.parseInt(line[chrIndex]))+"\t"+start+"\t"+stop+"\t"+ext.formDeci(value, SIG_FIGS));
						} else if (!ignoreWhenAbsentInResultsFile) {
							writer.println("chr"+Positions.chromosomeNumberInverse(Integer.parseInt(line[chrIndex]))+"\t"+start+"\t"+stop+"\t"+ext.formDeci(-2, SIG_FIGS));
						}
						prevStart = start;
					}
			        
			        reader.close();
		        } catch (FileNotFoundException fnfe) {
		        	log.reportError("Error: file \""+mapFile+"\" not found in current directory");
		        	log.reportException(fnfe);
			        return;
		        } catch (IOException ioe) {
		        	log.reportError("Error reading file \""+mapFile+"\"");
		        	log.reportException(ioe);
			        return;
		        }        
	        }
	        writer.close();
        } catch (Exception e) {
        	log.reportError("Error writing to file \""+filename+"\"");
        	log.reportException(e);
        }
	}
}
