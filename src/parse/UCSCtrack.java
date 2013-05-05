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
        boolean displayMarkerName;
        int unplaced, par, mitochondrial;

        int mapMarkerIndex;
        boolean ignoreFirstLineOfMap, mapCommaDelimited, mapTabDelimited;
        
    	String[] filterValues, filterOps;
    	int[] filterCols;
    	Vector<String> filters;
    	boolean passedFilter;
        
		paramV = Files.parseControlFile(filename, "ucsc", new String[] {
				"plink.assoc 1 8 neg_log chr=0 start=2 stop=2 out=assoc.bed", 
				"results.xln tab noHeader 0 1 neg_log map=plink.bim mapMarkerIndex=1 chr=0 start=3 stop=3 ignoreWhenAbsentInResultsFile", 
				"results.csv , ignorecase 0 1 map=markers.csv mapComma header mapMarkerIndex=1 chr=0 start=3 stop=3"
				}, log);
		if (paramV == null) {
			log.reportError("Failed...");
			return;
		}

		while (paramV.size() > 0) {
	    	line = paramV.remove(0).trim().split("[\\s]+");
	    	resultsFile = line[0];
	    	markerIndex = mapMarkerIndex = valueIndex = chrIndex = startIndex = stopIndex = -1;
	    	neg_log = false;
	    	ignoreCase = false;
	    	commaDelimited = false;
	    	tabDelimited = false;
	    	ignoreFirstLine = true;
	    	displayMarkerName = false;

	    	mapFile = null;
        	separateMap = false;
	    	ignoreCase = false;
	    	mapCommaDelimited = false;
	    	mapTabDelimited = false;
	    	ignoreFirstLineOfMap = false;
	    	ignoreWhenAbsentInResultsFile = false;
	    	outfile = ext.rootOf(resultsFile)+"_described.xln";

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
	    		} else if (line[j].equals("displayMarkerName")) {
	    			displayMarkerName = true;
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
	    		} else if (line[j].equalsIgnoreCase("mapHeader")) {
	    			ignoreFirstLineOfMap = true;
	    		} else if (line[j].startsWith("map=")) {
	    			mapFile = ext.parseStringArg(line[j], null);
	    			separateMap = true;
	    		} else if (line[j].equals("mapComma")) {
	    			mapCommaDelimited = true;
	    		} else if (line[j].equals("mapTab")) {
	    			mapTabDelimited = true;
	    		} else if (line[j].equals("ignoreWhenAbsentInResultsFile")) {
	    			ignoreWhenAbsentInResultsFile = true;
	    		} else if (line[j].startsWith("mapMarkerIndex=")){
	    			mapMarkerIndex = ext.parseIntArg(line[j]);
	    		} else {
	    			log.reportError("Error - don't know what to do with argument: "+line[j]);
	    		}
	        }
	    	
	    	if (separateMap && mapMarkerIndex == -1) {
	    		log.reportError("Error - separate map provided, but not a mapMarkerIndex");
	    		return;
	    	}
	    	if (chrIndex == -1 || startIndex == -1) {
	    		log.reportError("Error - the chr index and start index must be provided, regardless of whether it's in the results file or the map file");
	    		return;
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

	    	par = 0;
	    	unplaced = 0;
	    	mitochondrial = 0;
	    	hash = new Hashtable<String,Double>();
	    	try {
		        writer = Files.getAppropriateWriter(outfile);
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
				    		if (valueIndex == -1) {
								value = 1;
							} else if (line[valueIndex].equals(".")||line[valueIndex].toUpperCase().startsWith("NA")) {
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
								if (line[chrIndex].equals("0")) {
									unplaced++;
								} else if (line[chrIndex].equals("25")) {
									writer.println("chrX\t"+start+"\t"+stop+"\t"+(displayMarkerName?"\t"+line[markerIndex]:"")+"\t"+ext.formDeci(-2, SIG_FIGS));
									writer.println("chrY\t"+start+"\t"+stop+"\t"+(displayMarkerName?"\t"+line[markerIndex]:"")+"\t"+ext.formDeci(-2, SIG_FIGS));
									par++;
								} else if (line[chrIndex].equals("26")) {
									mitochondrial++;
								} else {
									writer.println("chr"+Positions.chromosomeNumberInverse(Integer.parseInt(line[chrIndex]))+"\t"+start+"\t"+stop+(displayMarkerName?"\t"+line[markerIndex]:"")+"\t"+ext.formDeci(value, SIG_FIGS));
								}
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
			        try {
				        reader = new BufferedReader(new FileReader(mapFile));
				        if (ignoreFirstLineOfMap) {
				        	reader.readLine();
				        }
				        prevStart = -1;
						while (reader.ready()) {
							line = reader.readLine().trim().split(mapCommaDelimited?",":(mapTabDelimited?"\t":"[\\s]+"));
							start = Integer.parseInt(line[startIndex]);
							stop = Integer.parseInt(line[stopIndex]);
							if (start == stop) {
								stop += 1;
							}
							if (start == prevStart) {
								start++;
								stop++;
							}
							if (hash.containsKey(ignoreCase?line[mapMarkerIndex].toLowerCase():line[mapMarkerIndex])) {
								value = hash.get(ignoreCase?line[mapMarkerIndex].toLowerCase():line[mapMarkerIndex]).doubleValue();
								if (line[chrIndex].equals("0")) {
									unplaced++;
								} else if (line[chrIndex].equals("25")) {
									writer.println("chrX\t"+start+"\t"+stop+"\t"+(displayMarkerName?"\t"+line[mapMarkerIndex]:"")+"\t"+ext.formDeci(-2, SIG_FIGS));
									writer.println("chrY\t"+start+"\t"+stop+"\t"+(displayMarkerName?"\t"+line[mapMarkerIndex]:"")+"\t"+ext.formDeci(-2, SIG_FIGS));
									par++;
								} else if (line[chrIndex].equals("26")) {
									mitochondrial++;
								} else {
									writer.println("chr"+Positions.chromosomeNumberInverse(Integer.parseInt(line[chrIndex]))+"\t"+start+"\t"+stop+(displayMarkerName?"\t"+line[mapMarkerIndex]:"")+"\t"+ext.formDeci(value, SIG_FIGS));
								}
							} else if (!ignoreWhenAbsentInResultsFile) {
								if (line[chrIndex].equals("0")) {
									unplaced++;
								} else if (line[chrIndex].equals("25")) {
									writer.println("chrX\t"+start+"\t"+stop+"\t"+(displayMarkerName?"\t"+line[mapMarkerIndex]:"")+"\t"+ext.formDeci(-2, SIG_FIGS));
									writer.println("chrY\t"+start+"\t"+stop+"\t"+(displayMarkerName?"\t"+line[mapMarkerIndex]:"")+"\t"+ext.formDeci(-2, SIG_FIGS));
									par++;
								} else if (line[chrIndex].equals("26")) {
									mitochondrial++;
								} else {
									writer.println("chr"+Positions.chromosomeNumberInverse(Integer.parseInt(line[chrIndex]))+"\t"+start+"\t"+stop+"\t"+(displayMarkerName?"\t"+line[mapMarkerIndex]:"")+"\t"+ext.formDeci(-2, SIG_FIGS));
								}
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

	    	if (unplaced > 0) {
	        	log.reportError("Warning - there were "+unplaced+" marker(s) on chromosome zero that were not exported to "+outfile);
	        }
	    	if (mitochondrial > 0) {
	        	log.reportError("Warning - there were "+mitochondrial+" marker(s) from the mitochondria that were not exported to "+outfile);
	        }
	    	if (par > 0) {
	        	log.reportError("Warning - there were "+par+" marker(s) in pseudo-autosomal regions that were duplicated, one for the X and one for the Y inside "+outfile);
	        }
		}    	
	}
}
