package parse;

import java.io.*;
import java.util.*;

import stats.Maths;
import common.*;

public class GenParser {
	public static final int ADD = 0;
	public static final int SUBTRACT = 1;
	public static final int MULTIPLY = 2;
	public static final int DIVIDE = 3;
	public static final String OPERATORS = "+-*/";
	public static final String[][] LINUX_SUBSTITUTIONS = {{"<=", "LTE"}, {"<", "LT"}, {">=", "GTE"}, {">", "GT"},
										{"=>", "replace"}, {"$", "F"}, {"*", "x"}, {"/", "div"}};

	private Logger log;
	private BufferedReader reader;
	private boolean commaDelimited, tabDelimited, serializeFlag;
	private String[] colNames, originalColumnNames, filterValues, filterOps;
	private int[] cols, filterCols;
	private String[][] computes;
	private String[][] replaces, ops;
	private String outfile;
	private String serialFilename;
	private boolean failed;
	private boolean header;
	private String[] failCodes;
	private boolean forceFailCodes;
	private String replaceBlanks;

	public GenParser(String[] line, Logger logger) {
		String file;
		String[] tokens;
		Vector<String> filters, columns, comps;
		Vector<String[]> replacesV;
		String trav;
		int skip;
		String[][] specialOps;
		
		log = logger;
		commaDelimited = false;
		tabDelimited = false;
		serializeFlag = false;
		header = true;
		failed = false;
		replaceBlanks = null;
    	file = line[0];
    	skip = -2;
    	outfile = ext.rootOf(file)+"_parsed.xln";
    	filters = new Vector<String>();
    	columns = new Vector<String>();
    	replacesV = new Vector<String[]>();
    	serialFilename = parseSerialFilename(line);
    	forceFailCodes = false;
    	for (int j = 1; j<line.length; j++) {
    		if (line[j].equalsIgnoreCase("noHeader")) {
    			skip = 0;
    			header = false;
    		} else if (line[j].startsWith("skip=")) {
    			skip = Integer.parseInt(line[j].split("=")[1]);
    		} else if (line[j].startsWith("out=")) {
    			outfile = line[j].split("=")[1];
    		} else if (line[j].startsWith("replace=")) {
    			replaceBlanks = line[j].split("=")[1];
    		} else if (line[j].equals(",")) {
    			commaDelimited = true;
    		} else if (line[j].equalsIgnoreCase("fail")) {
    			forceFailCodes = true;
    		} else if (line[j].equalsIgnoreCase("tab")) {
    			tabDelimited = true;
    		} else if (line[j].equalsIgnoreCase("ser")) {
    			serializeFlag = true;
    		} else if (line[j].startsWith("!")) {
    			filters.add(line[j].substring(1));
    		} else if (line[j].indexOf("=>") > 0) {
    			replacesV.add(new String[] {line[j].substring(0, line[j].indexOf("=>")), line[j].substring(line[j].indexOf("=>")+2)});
    		} else {
    			columns.add(line[j]);
    		}
        }
    	filterCols = new int[filters.size()];
    	filterOps = new String[filters.size()];
    	filterValues = new String[filters.size()];
    	replaces = Matrix.toStringArrays(replacesV);

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

    	cols = new int[columns.size()];
    	colNames = new String[columns.size()];
    	failCodes = Array.stringArray(columns.size(), ".");
    	comps = new Vector<String>();
    	for (int j = 0; j<columns.size(); j++) {
    		trav = columns.elementAt(j);
    		if (trav.indexOf(";") >= 0) {
    			failCodes[j] = trav.split(";", -1)[1];
    			trav = trav.split(";", -1)[0];
    		}
    		tokens = trav.split("=");
    		colNames[j] = tokens.length==1?null:tokens[1];
  			if (colNames[j] != null && colNames[j].equals("")) {
   				colNames[j] = null;
   			}
    		if (tokens[0].startsWith("$")) {
    			comps.add(tokens[0].substring(1));
    			cols[j] = -1*comps.size();
    			if (colNames[j] == null) {
    				colNames[j] = "Function"+comps.size();
    			}
    		} else {
    			try {
    				cols[j] = Integer.parseInt(tokens[0]);
				} catch (NumberFormatException nfe) {
					log.reportError("Error - invalid operator or something... \""+tokens[0]+"\"");
					log.reportException(nfe);
					failed = true;
				}
    		}
        }
    	
    	computes = new String[comps.size()][];
    	ops = new String[comps.size()][];
    	for (int i = 0; i<computes.length; i++) {
    		if (comps.elementAt(i).startsWith("@")) {
    			computes[i] = new String[] {comps.elementAt(i)};
    		} else if (comps.elementAt(i).startsWith("%")) {
        		computes[i] = comps.elementAt(i).split(";", -1);
    		} else {
	    		specialOps = ext.getOperatorsOperatorIndicesAndSplit(comps.elementAt(i), OPERATORS);
	    		ops[i] = specialOps[0];
	    		computes[i] = specialOps[2];
    		}
        }
    	
    	try {
            reader = new BufferedReader(new FileReader(file));
            if (skip == -2) {
            	originalColumnNames = ext.replaceAllWith(reader.readLine().trim(), replaces).split(commaDelimited?",":(tabDelimited?"\\t":"[\\s]+"));
            	for (int j = 0; j<cols.length; j++) {
            		if (colNames[j] == null) {
            			colNames[j] = originalColumnNames[cols[j]];
            		}
                }
            } else {
            	for (int i = 0; i<skip; i++) {
            		reader.readLine();
                }
            }
        } catch (FileNotFoundException fnfe) {
            log.reportError("Error: file \""+file+"\" not found in current directory");
            return;
        } catch (IOException ioe) {
            log.reportError("Error reading file \""+file+"\"");
            return;
        }
	}

	public static String parseSerialFilename(String[] line) {
		String str;

		str = line[0];
		for (int i = 0; i<line.length; i++) {
			str += "_"+ext.replaceAllWith(line[i], LINUX_SUBSTITUTIONS);
        }
		str += ".ser";

		return str;
	}
	
	public String getSerialFilename() {
		return serialFilename;
	}
	
	public boolean toBeSerialized() {
		return serializeFlag;
	}
	
	public int getNumColumns() {
		return cols.length;
	}

	public String[] getColumnNames() {
		return colNames;
	}
	
	public String[] getOriginalColumnNames() {
		return originalColumnNames;
	}
	
	public boolean ready() {
		if (failed) {
			return false;
		}
		try {
			return reader.ready();
		} catch (IOException ioe) {
			log.reportError("Error accessing GenParser file");
			log.reportException(ioe);
			return false;
		}
	}
	
	public String getOutfile() {
		return outfile;
	}
	
	public void close() {
		try {
			reader.close();
		} catch (IOException ioe) {
			log.reportError("Error closing GenParser file");
			log.reportException(ioe);
		}
	}
	
	public boolean hasHeader() {
		return header;
	}
	
	public boolean forcingFailCodes() {
		return forceFailCodes;
	}
	
	public String[] getFailCodes() {
		return failCodes;
	}
	
	// returns null if the line fails criteria
	public String[] nextLine() {
		String[] line, parsed;
		boolean passedFilter;
		int spi;
		double num, d;
		boolean match;
		
		try {
			line = ext.replaceAllWith(commaDelimited||tabDelimited?reader.readLine():reader.readLine().trim(), replaces).split(commaDelimited?",":(tabDelimited?"\\t":"[\\s]+"), -1);
		} catch (IOException ioe) {
			log.reportException(ioe);
			return Array.stringArray(cols.length, "NaN");
		}
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
    	
    	parsed = Array.stringArray(cols.length, "!found or failed");
    	if (forceFailCodes && cols[0] >= 0) {
    		parsed[0] = line[cols[0]];
    	}
    	if (passedFilter) {
    		for (int j = 0; j<cols.length; j++) {
    			if (cols[j] < 0) {
    				spi = -1*cols[j]-1;
    				if (computes[spi].length == 1 && computes[spi][0].startsWith("@")) {
        				parsed[j] = computes[spi][0].substring(1);
    				} else if (computes[spi][0].startsWith("%")) {
            			parsed[j] = line[Integer.parseInt(computes[spi][0].substring(1))];
            			match = false;
            			for (int i = 1; i<computes[spi].length-1; i++) {
            				if (computes[spi][i].startsWith(parsed[j]+":")) {
            					parsed[j] = computes[spi][i].substring(computes[spi][i].indexOf(":")+1);
            					match = true;
            				}
                        }
            			if (!match && !computes[spi][computes[spi].length-1].equals("")) {
            				parsed[j] = computes[spi][computes[spi].length-1];
            			}
    				} else {
	    				num = Double.NaN;
	    				for (int i = 0; i<computes[spi].length; i++) {
	    					if (computes[spi][i].startsWith("#")) {
	    						d = procDouble(computes[spi][i].substring(1), Double.NaN, log); 
	    					} else {
	    						try {
	    							d = procDouble(line[Integer.parseInt(computes[spi][i])], Double.NaN, log);
	    						} catch (NumberFormatException nfe) {
	    							log.reportError("Error - invalid operator or something... \""+computes[spi][i]+"\"");
	    							log.reportException(nfe);
	    							failed = true;
	    							return parsed;
	    						}
	    					}    					
	    					if (i == 0) {
	    						num = d;
	    					} else {
	        					switch (ops[spi][i-1].charAt(0)) {
	        					case '+':
	        						num += d;
	        						break;
	        					case '-':
	        						num -= d;
	        						break;
	        					case '*':
	        						num *= d;
	        						break;
	        					case '/':
	        						num /= d;
	        						break;
	        					}
	    					}
	                    }
	    				parsed[j] = ext.formDeci(num, 10, false);
    				}
    			} else {
    				parsed[j] = line[cols[j]];
    				if (parsed[j].equals("") && replaceBlanks != null) {
        				parsed[j] = replaceBlanks;
    				}
    			}
            }
    	} else if (!forceFailCodes) {
    		parsed = null;
    	}

		for (int j = 0; parsed != null && j<cols.length; j++) {
			if (parsed[j].equals("!found or failed")) {
				parsed[j] = failCodes[j];
			}
		}
    	
    	return parsed;
	}
	
	public static void parse(String[] line, Logger log) {
		PrintWriter writer;
		long time;
		GenParser parser;
		String[] header, trav;
		int nonNull;
		String filename, delimiter;
		
        time = new Date().getTime();
        parser = new GenParser(line, log);
    	
    	try {
    		filename = parser.getOutfile();
    		if (filename.endsWith(".csv")) {
        		delimiter = ",";
    		} else {
    			delimiter = "\t";
    		}
            writer = new PrintWriter(new FileWriter(filename));
            log.report("Parsing '"+line[0]+"'");
            header = parser.getColumnNames();
            nonNull = 0;
            for (int i = 0; i<header.length; i++) {
            	if (header[i] != null) {
            		nonNull++;
            	}
            }
            if (nonNull > 0	&& parser.hasHeader()) {
            	writer.println(Array.toStr(header, delimiter));
            }
            while (parser.ready()) {
            	trav = parser.nextLine();
            	if (trav != null) {
            		writer.println(Array.toStr(trav, delimiter));
            	}
            }
            parser.close();
            writer.close();
        } catch (FileNotFoundException fnfe) {
            log.reportError("Error: file \""+line[0]+"\" not found in current directory");
            return;
        } catch (IOException ioe) {
            log.reportError("Error reading file \""+line[0]+"\"");
            return;
        }	
		log.report(" in "+ext.getTimeElapsed(time));        
	}
	
	public static double procDouble(String str) {
		return procDouble(str, Double.NaN, null);
	}
	
	public static double procDouble(String str, double defaultNull, Logger log) {
		if (ext.isMissingValue(str)) {
			return defaultNull;
		} else {
			try {
				return Double.parseDouble(str);
			} catch (Exception e) {
				if (log != null) {
					log.reportError("Error - '"+str+"' is not a valid double");
				} else {
					System.out.println("Error - '"+str+"' is not a valid double");
				}
				return defaultNull;
			}
		}
	}

	public static void parseFile(String filename, Logger log) {
        String[][] params;

		params = Files.parseControlFile(filename, false, "parse", new String[] {"data.csv , out=parsed_data.xln !11!NA !10>-5 !10<5 !1=2 0 1=chr 2=pos 10 11 $12*13=effN;0"}, log);
		if (params != null) {
        	parse(params[0], log);
		}
	}
	
	public static void main(String[] args) {
	    String usage = "\n"+
	    "parse.GenParser requires 2+ arguments and will take the form of something like:\n"+
	    "    data.csv , out=parsed_data.xln !1=2 !11!NA !10>-5 !10<5 0 1 2 10 11\n"+
	    "\n"+
	    "  where at a minimum you will have the file name and the columns you are intrested in.\n"+
	    "\n"+
	    "  Optional parameters include:\n"+
	    "    ,              file is comma delimited (default is whitespace delimited)\n"+
	    "    tab            file is tab delimited (i.e. can include blanks and spaces)\n"+
	    "    out=file.txt   delineates a specific output file name\n"+
	    "    skip=5         skip the first 5 rows (default is skip=1)\n"+
	    "    noHeader       will start parsing the first row (i.e. skip=0)\n"+
	    "    RS=>rs         replaces all instances of 'RS' in the line with 'rs'\n"+
	    "    1=chr          report index 1 and call the column \"chr\"\n"+
	    "    !              anything starting with ! indicates a filter\n"+
	    "    !1=2           includes only those rows with a 2 in column index 1 (column 2)\n"+
	    "    !11!NA         will filter out any line where column index 11 equals 'NA'\n"+
	    "    !10>-5         includes lines where the number in column index 10 is greater than -5\n"+
	    "    $              anything starting with $ indicates a computation\n"+
	    "    $12*13=effN    multiply column index 13 by column index 12\n"+
	    "    $#12*13=rate   multiply column index 13 by a constant (here 12.0)\n"+
	    "    replace=.		replace empty columns with a period\n"+
	    "    fail			normally row is absent if the row fails criteria, this force fail codes\n"+
	    "    12=Needed;NA	fail code for this column will be set to \"NA\" instead of the default \".\"\n"+
	    "    @1.1 or @fini	set to this value no matter what\n"+
	    "    $%2;3:1;2:1;1:0;.	check value in column index 2 and if it is 3 or 2 then replace with 1, if 1 then replace with 0, otherwise set to missing\n"+
	    "";

//	    args = new String [] {"subset.genome_parsed.xln out=admixedUnrelatedsCheck.txt 0 1 2 3 6 7 8 9 !9>0.2"};
	    if (args.length == 1) {
	    	args = ext.removeQuotes(args[0]).trim().split("[\\s]+");
	    }
	    if (args.length < 2) {
		    System.err.println(usage);
		    System.exit(1);
	    }
	    
	    try {
		    parse(args, new Logger());
	    } catch (Exception e) {
		    e.printStackTrace();
	    }
    }
}
