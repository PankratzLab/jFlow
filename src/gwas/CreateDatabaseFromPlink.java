package gwas;

import java.io.*;
import java.util.*;

import mining.Transformations;
import common.*;
import filesys.*;
import stats.*;

public class CreateDatabaseFromPlink {
	public static final String[] POSSIBLE_MODELS = {"ADD", "DOM", "REC"};
	public static final int ADDITIVE_MODEL = 0;
	public static final int DOMINANT_MODEL = 1;
	public static final int RECESSIVE_MODEL = 2;
	
	public static void createDatabase(String filename, Logger log) {
		createDatabase("", filename, log);
	}
	
	public static void createDatabase(String dir, String filename, Logger log) {
		BufferedReader reader;
        PrintWriter writer;
        String[] line;
        Hashtable<String,String> notes;
//        Vector<String> v = new Vector<String>();
        int count;
        String[] markerNames;
        String pedfile, mapfile, freqfile, outfile;
        String[][] alleles;
        Vector<int[]> vModels;
        int[][] models;
        String allele;
        int markerIndex, modelIndex;
        String[][] params;
        SnpMarkerSet markerSet;
        boolean labelAlleles, listPositions, displayNotes, maskModel, maskFIDs, maskSex, maskAffStat, allMarkers, commaDelimited;
        byte[] chrs;
        int[] positions;
        String delimiter;
        
        params = Files.parseControlFile(filename, false, "plink", new String[] {"plink.ped", "plink.map", "plink.frq", "newDatabase.xln labelAlleles=false listPositions=false displayNotes=false maskModel=true maskFamID=false maskSex=false maskAffectionStatus=false allMarkers=true commaDelimited=false", "# The following are examples if you just want a few markers from the file instead of all", "# rs11550605\tADD", "# rs10808555\tREC ADD note=CriticalMutation", "# rs7837328\tADD DOM REC", "# rs3730881\tADD DOM"}, log);
        
        if (params == null) {
        	return;
        }
        
        labelAlleles = false;
        listPositions = false;
        displayNotes = false;
        maskModel = false;
        maskFIDs = false;
        maskSex = false;
        maskAffStat = false;
        allMarkers = false;
        commaDelimited = false;
    	pedfile = params[0][0];
    	mapfile = params[1][0];
    	markerSet = new SnpMarkerSet(dir+mapfile, false, log);
        markerNames = markerSet.getMarkerNames();
        chrs = markerSet.getChrs();
        positions = markerSet.getPositions();
        notes = new Hashtable<String, String>();
        
    	freqfile = params[2][0];
    	alleles = HashVec.loadFileToStringMatrix(dir+freqfile, true, new int[] {2,3}, false);
    	outfile = params[3][0];    	
    	for (int i = 1; i < params[3].length; i++) {
    		if (params[3][i].startsWith("labelAlleles=")) {
        		labelAlleles = ext.parseBooleanArg(params[3][i], log);
    		} else if (params[3][i].startsWith("listPositions=")) {
    			listPositions = ext.parseBooleanArg(params[3][i], log);
    		} else if (params[3][i].startsWith("displayNotes=")) {
    			displayNotes = ext.parseBooleanArg(params[3][i], log);
    		} else if (params[3][i].startsWith("maskModel=")) {
    			maskModel = ext.parseBooleanArg(params[3][i], log);
    		} else if (params[3][i].startsWith("maskFamID=")) {
    			maskFIDs = ext.parseBooleanArg(params[3][i], log);
    		} else if (params[3][i].startsWith("maskSex=")) {
    			maskSex = ext.parseBooleanArg(params[3][i], log);
    		} else if (params[3][i].startsWith("maskAffectionStatus=")) {
    			maskAffStat = ext.parseBooleanArg(params[3][i], log);    			
    		} else if (params[3][i].startsWith("allMarkers=")) {
    			allMarkers = ext.parseBooleanArg(params[3][i], log);
    		} else if (params[3][i].startsWith("commaDelimited=")) {
    			commaDelimited = ext.parseBooleanArg(params[3][i], log);
    		} else {
    			log.reportError("Error - don't know what to do with argument: "+params[3][i]);
    		}
		}

    	if (!Array.equals(markerNames, HashVec.loadFileToStringArray(dir+freqfile, true, new int[] {1}, false), true)) {
    		log.reportError("Error - the freq file does not match the map file");
    		return;
    	}

        vModels = new Vector<int[]>();
        if (allMarkers) {
        	if (params.length > 4) {
        		System.err.println("Error - The allMarkers flagged was set to true, so all markers will be exported and the following lines will be ignored");
        		for (int i = 4; i < params.length; i++) {
        			System.err.println("'"+Array.toStr(params[i], "  ")+"'");
				}
        	}
        	if (!maskModel) {
        		maskModel = true;
        		System.out.println("Since all markers are being exported, masking the default model ("+POSSIBLE_MODELS[0]+")");
        	}
        	for (int i = 0; i < markerNames.length; i++) {
    			vModels.add(new int[] {i, 0});
			}
        } else {
	        for (int i = 4; i < params.length; i++) {
	        	markerIndex = ext.indexOfStr(params[i][0], markerNames);
	        	if (markerIndex == -1) {
	        		log.reportError("Error - marker '"+params[i][0]+"' was not found in "+mapfile);
	        	} else {
	        		for (int j = 1; j<params[i].length; j++) {
	        			if (params[i][j].toLowerCase().startsWith("note=")) {
	        				notes.put(params[i][0], params[i][j].substring(5));
	        			} else {
		        			modelIndex = ext.indexOfStr(params[i][j], POSSIBLE_MODELS);
			        		if (modelIndex == -1) {
			        			log.reportError("Error - '"+params[i][j]+"' is an invalid model");
			        		} else {
			        			vModels.add(new int[] {markerIndex, modelIndex});
			        		}
	        			}
	                }
	        	}
			}
        }        
        models = Matrix.toMatrix(vModels);
        
        if (commaDelimited) {
        	delimiter = ",";
        } else {
        	delimiter = "\t";
        }
        try {
            reader = new BufferedReader(new FileReader(dir+pedfile));
            writer = new PrintWriter(new FileWriter(dir+outfile));
            if (listPositions) {
                writer.print((maskFIDs?"":delimiter)+(maskSex?"":delimiter)+(maskAffStat?"":delimiter)+"Chr:");
                for (int i = 0; i<models.length; i++) {
                	writer.print(delimiter+chrs[models[i][0]]);
                }
                writer.println();
                writer.print((maskFIDs?"":delimiter)+(maskSex?"":delimiter)+(maskAffStat?"":delimiter)+"Position:");
                for (int i = 0; i<models.length; i++) {
                	writer.print(delimiter+positions[models[i][0]]);
                }
                writer.println();
            }
            if (labelAlleles) {
                writer.print((maskFIDs?"":delimiter)+(maskSex?"":delimiter)+(maskAffStat?"":delimiter)+"Minor allele, counted:");
                for (int i = 0; i<models.length; i++) {
                	writer.print(delimiter+alleles[models[i][0]][0]);
                }
                writer.println();
                writer.print((maskFIDs?"":delimiter)+(maskSex?"":delimiter)+(maskAffStat?"":delimiter)+"Other allele:");
                for (int i = 0; i<models.length; i++) {
                	writer.print(delimiter+alleles[models[i][0]][1]);
                }
                writer.println();
            }
            if (displayNotes) {
                writer.print((maskFIDs?"":delimiter)+(maskSex?"":delimiter)+(maskAffStat?"":delimiter)+"Note:");
                for (int i = 0; i<models.length; i++) {
                	writer.print(delimiter+(notes.containsKey(markerNames[models[i][0]])?notes.get(markerNames[models[i][0]]):"."));
                }
                writer.println();
            }
            writer.print((maskFIDs?"":"FID"+delimiter)+"IID"+(maskSex?"":delimiter+"Gender")+(maskAffStat?"":delimiter+"Affected"));
            for (int i = 0; i<models.length; i++) {
            	writer.print(delimiter+markerNames[models[i][0]]+(maskModel?"":"_"+POSSIBLE_MODELS[models[i][1]]));
            }
            writer.println();
            while (reader.ready()) {
            	line = reader.readLine().trim().split("[\\s]+");
            	writer.print((maskFIDs?"":line[0]+delimiter)+line[1]+(maskSex?"":delimiter+(line[4].equals("1")?"1":(line[4].equals("2")?"0":".")))+(maskAffStat?"":delimiter+(line[5].equals("2")?"1":(line[5].equals("1")?"0":"."))));
            	for (int i = 0; i<models.length; i++) {
                	count = 0;
            		for (int j = 0; j<2; j++) {
            			allele = line[6+models[i][0]*2+j];
            			if (allele.equals("0")) {
            				count = -9;
            			} else if (allele.equals(alleles[models[i][0]][0])) {
            				count++;
            			} else if (!allele.equals(alleles[models[i][0]][1])) {
            				log.reportError("Error - mismatched alleles, expecting "+alleles[models[i][0]][0]+" or "+alleles[models[i][0]][1]+" and got '"+allele+"'");
            			}

                    }
            		if (count < 0) {
            			writer.print(delimiter+".");
            		} else {
            			switch (models[i][1]) {
            			case ADDITIVE_MODEL:
            				writer.print(delimiter+count);
            				break;
            			case DOMINANT_MODEL:
            				writer.print(delimiter+(count>=1?1:0));
            				break;
            			case RECESSIVE_MODEL:
            				writer.print(delimiter+(count>=2?1:0));
            				break;
            			default:
            				log.reportError("Error - model not translated");
            				writer.print(delimiter+"."+count);
            				break;
            			}
            		}
                }
            	writer.println();
            }
            writer.close();
            reader.close();
        } catch (FileNotFoundException fnfe) {
        	log.reportError("Error: file \""+pedfile+"\" not found in current directory");
        	log.reportException(fnfe);
            return;
        } catch (IOException ioe) {
        	log.reportError("Error reading file \""+pedfile+"\"");
        	log.reportException(ioe);
        	return;
        }
	}

	public static void createScore(String filename, Logger log) {
		BufferedReader reader;
        PrintWriter writer;
        String[] line, header;
        Vector<String[]> v = new Vector<String[]>();
        String dbfile, outfile;
//        double[] freqs;
        int[] varIndices;
        double[] weights;
        double constant;
        String[] ids;
        double[] scores, normalized;
        int numInds;
        boolean normalize;
        double[][] data;
        String[][] params;

		params = Files.parseControlFile(filename, false, "score", new String[] { "database.xln normalize", "newScore.xln", "Gender  1.006", "rs3730881_ADD  -2.297", "rs10808555_REC  0.800", "rs3087367_DOM  -0.498", "Constant  -0.360" }, log);
		if (params != null) {
        	dbfile = params[0][0];
        	normalize = false;
        	for (int i = 1; i<params[0].length; i++) {
        		if (params[0][i].equals("normalize")) {
        			normalize = true;
        		} else {
        			log.reportError("Error - unknown argument after database filename: "+params[0][i]);
        		}
            }
        	outfile = params[1][0];

        	constant = 0;
        	for (int i = 2; i < params.length; i++) {
	        	if (params[i][0].toLowerCase().startsWith("constant")) {
	        		constant = Double.parseDouble(params[i][1]);
	        	} else {
	        		v.add(params[i]);
	        	}
			}

        	varIndices = new int[v.size()];
	        weights = new double[v.size()];
	        numInds = Files.countLines(dbfile, 1);
	        try {
                reader = new BufferedReader(new FileReader(dbfile));
                header = reader.readLine().trim().split("[\\s]+");
		        for (int i = 0; i<v.size(); i++) {
		        	line = v.elementAt(i);

		        	varIndices[i] = ext.indexOfStr(line[0], header);
		        	if (varIndices[i] == -1) {
		        		log.reportError("Error - variable '"+line[0]+"' was not found in "+dbfile);
		        	} else {
		        		weights[i] = Double.parseDouble(line[1]);
		        	}
		        }
                
		        ids = new String[numInds];
		        scores = new double[numInds];
		        
		        if (normalize) {
		        	data = new double[varIndices.length][numInds];
			        for (int i = 0; i<numInds; i++) {
	                	line = reader.readLine().trim().split("[\\s]+");
	                	ids[i] = line[0]+"\t"+line[1];
	                	for (int j = 0; j<varIndices.length; j++) {
	                		data[j][i] = Double.parseDouble(line[varIndices[j]]);
	                	}
                    }
                	for (int j = 0; j<varIndices.length; j++) {
                		data[j] = Transformations.transform(data[j], Transformations.NORMALIZE);
                	}
			        for (int i = 0; i<numInds; i++) {
	                	scores[i] = constant;
	                	for (int j = 0; j<weights.length; j++) {
	                		scores[i] += weights[j]*data[j][i];
	                	}
                    }
		        } else {
			        for (int i = 0; i<numInds; i++) {
	                	line = reader.readLine().trim().split("[\\s]+");
	                	ids[i] = line[0]+"\t"+line[1];
	                	scores[i] = constant;
	                	for (int j = 0; j<weights.length; j++) {
	                		scores[i] += weights[j]*Double.parseDouble(line[varIndices[j]]);
	                	}
                    }
		        }
                reader.close();
                
                normalized = Transformations.transform(scores, Transformations.NORMALIZE);
                try {
	                writer = new PrintWriter(new FileWriter(outfile));
	                writer.println(header[0]+"\t"+header[1]+"\tscore\tnormalized\tlogit");
	                for (int i = 0; i<ids.length; i++) {
	                	writer.println(ids[i]+"\t"+scores[i]+"\t"+normalized[i]+"\t"+Maths.logit(scores[i]));
                    }
	                writer.close();
                } catch (Exception e) {
                	log.reportError("Error writing to "+outfile);
                	log.reportException(e);
                }
            } catch (FileNotFoundException fnfe) {
            	log.reportError("Error: file \""+dbfile+"\" not found in current directory");
                return;
            } catch (IOException ioe) {
            	log.reportError("Error reading file \""+dbfile+"\"");
                return;
            }
		}
	}
	
	public static void createCountsMatrix(String filename, Logger log) {
		BufferedReader reader;
        PrintWriter writer;
        String[] line;
        int count, index;
        String[] markerNames;
        String tpedfile, freqfile, outfile, temp;
        String[][] alleles;
        String allele;
        String[][] params;
        boolean anyMissing = false;
        int numIndividuals = -1;
        Hashtable<String,String> hash;
        IntVector monomorphs;
        
        params = Files.parseControlFile(filename, false, "simpleM", new String[] {"plink.tped", "plink.frq", "simpleM.txt", "# in order to generate the tped file, use the following command:", "# plink --bfile plink --recode --transpose"}, log);
        if (params == null) {
        	return;
        }
        
    	tpedfile = params[0][0];
    	freqfile = params[1][0];
        markerNames = HashVec.loadFileToStringArray(freqfile, true, new int[] {1}, false);
    	alleles = HashVec.loadFileToStringMatrix(freqfile, true, new int[] {2,3}, false);
    	outfile = params[2][0];    	

        monomorphs = new IntVector();
        try {
            reader = new BufferedReader(new FileReader(tpedfile));
            writer = new PrintWriter(new FileWriter(outfile));
        	for (int i = 0; i<markerNames.length; i++) {
        		hash = new Hashtable<String, String>();
            	line = reader.readLine().trim().split("[\\s]+");
            	if (!line[1].equals(markerNames[i])) {
            		log.reportError("Error - the freq file does not match the map file at line "+(i+1)+"; expecting "+markerNames[i]+", found "+line[1]);
            		reader.close();
            		writer.close();
            		return;
            	}
            	if (i==0) {
            		numIndividuals = (line.length - 4) / 2;
            	} else if ((line.length - 4) / 2 != numIndividuals) {
            		log.reportError("Error - mismatched number of indiviudals in line "+(i+1)+"; due to previous lines, expecting "+numIndividuals+"*2+4="+(numIndividuals*2+4)+", found "+line.length); 
            	}
        		for (int j = 0; j<numIndividuals; j++) {
                   	count = 0;
               		for (int k = 0; k<2; k++) {
            			allele = line[4+j*2+k];
            			if (allele.equals("0")) {
            				count = -9;
            			} else if (allele.equals(alleles[i][0])) {
            				count++;
            			} else if (!allele.equals(alleles[i][1])) {
            				log.reportError("Error - mismatched alleles, expecting "+alleles[i][0]+" or "+alleles[i][1]+" and got '"+allele+"'");
            			}
                    }
            		if (count < 0) {
            			if (!anyMissing) {
            				log.reportError("Warning - there are missing values in the genotype data; setting value to zero; impute if you want better accuracy");
            				anyMissing = true;            			
            			}
            			count = 0;
            		}
               		hash.put(count+"", "");
            		writer.print((j==0?"":"\t")+count);
                }
            	writer.println();
            	if (hash.size() < 2) {
        			if (monomorphs.size() == 0) {
        				log.reportError("Warning - there are monomorphic markers in the genotype data; these will be filtered out since it would otherwise cause the simpleM algorithm to fail");
        			}
        			monomorphs.add(i);
            	}
            }
            writer.close();
            reader.close();
        } catch (FileNotFoundException fnfe) {
        	log.reportError("Error: file \""+tpedfile+"\" not found in current directory");
        	log.reportException(fnfe);
            return;
        } catch (IOException ioe) {
        	log.reportError("Error reading file \""+tpedfile+"\"");
        	log.reportException(ioe);
        	return;
        }
        
        if (monomorphs.size() > 0) {
        	log.report("There were "+monomorphs.size()+" monomorphic markers that will be filtered out");
        	Files.writeList(Array.subArray(markerNames, monomorphs.toArray()), ext.parseDirectoryOfFile(filename)+"monomorphs.dat");
    		temp = outfile+".bak";
    		new File(outfile).renameTo(new File(temp));
        	try {
				reader = new BufferedReader(new FileReader(temp));
				writer = new PrintWriter(new FileWriter(outfile));
				count = 0;
				index = 0;
				while (reader.ready()) {
					temp = reader.readLine();
					if (index < monomorphs.size() && count == monomorphs.elementAt(index)) {
						index++;
					} else {
						writer.println(temp);
					}
					count++;
				}
				writer.close();
				reader.close();
			} catch (FileNotFoundException fnfe) {
				System.err.println("Error: file \"" + temp + "\" not found in current directory");
				System.exit(1);
			} catch (IOException ioe) {
				System.err.println("Error reading file \"" + temp + "\"");
				System.exit(2);
			}
			new File(temp).delete();
        }
        
        try {
			reader = new BufferedReader(new FileReader("/home/npankrat/bin/simpleM_Exec.R"));
			writer = new PrintWriter(new FileWriter(ext.parseDirectoryOfFile(filename)+"simpleM.R"));
			while (reader.ready()) {
				temp = reader.readLine();
				if (temp.startsWith("fn_In <- ")) {
					writer.println("fn_In <- \""+ext.parseDirectoryOfFile(filename)+"simpleM.txt"+"\"");
				} else {
					writer.println(temp);
				}
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \"" + "/home/npankrat/bin/simpleM_Exec.R" + "\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \"" + "/home/npankrat/bin/simpleM_Exec.R" + "\"");
			System.exit(2);
		}
	}
	
	public static void toGWAF(String pedfile, String mapfile, String freqfile, String rlinker, String outfile) {
		BufferedReader reader;
        PrintWriter writer;
        String[] line;
        int count;
        String[] markerNames;
        String[][] alleles;
        String allele;
        SnpMarkerSet markerSet;
        Logger log;
		Hashtable<String, String> hash;
        
        log = new Logger();
        
    	markerSet = new SnpMarkerSet(mapfile, false, log);
        markerNames = markerSet.getMarkerNames();
    	alleles = HashVec.loadFileToStringMatrix(freqfile, true, new int[] {0,2,3}, false);
    	if (markerNames.length != alleles.length) {
    		log.reportError("Error - different number of markers in map file and frq file");
			System.exit(1);
    	}
    	for (int i = 0; i < markerNames.length; i++) {
    		if (markerNames[i].equals(alleles[i][0])) {
    			log.reportError("Error - marker name mismatch at index "+i+": found "+markerNames[i]+" in map file and "+alleles[i][0]+" in frq file");
    			System.exit(1);
    		}
		}
    	
    	if (rlinker == null) {
    		hash = new Hashtable<String, String>();
    	} else {
    		hash = HashVec.loadFileToHashString(rlinker, new int[] {0}, new int[] {2}, false, "", false, false, false);
    	}

        try {
            reader = new BufferedReader(new FileReader(pedfile));
            writer = new PrintWriter(new FileWriter(outfile));

            writer.print("id");
            for (int i = 0; i<markerNames.length; i++) {
            	writer.print(","+markerNames[i]);
            }
            writer.println();

            while (reader.ready()) {
            	line = reader.readLine().trim().split("[\\s]+");
            	if (rlinker != null) {
            		if (hash.containsKey(line[1])) {
            			line[1] = hash.get(line[1]);
            		} else {
            			System.err.println("Error - failed to lookup "+line[1]+" in the Rlinker file: "+rlinker);
            		}            		
            	}
            	writer.print(line[1]);
            	for (int i = 0; i<markerNames.length; i++) {
                	count = 0;
            		for (int j = 0; j<2; j++) {
            			allele = line[6+i*2+j];
            			if (allele.equals("0")) {
            				count = -9;
            			} else if (allele.equals(alleles[i][1])) {
            				count++;
            			} else if (!allele.equals(alleles[i][2])) {
            				System.err.println("Error - mismatched alleles, expecting "+alleles[i][1]+" or "+alleles[i][2]+" and got '"+allele+"'");
            			}

                    }
       				writer.print(","+(count < 0?"":count));
                }
            	writer.println();
            }
            writer.close();
            reader.close();
        } catch (FileNotFoundException fnfe) {
        	log.reportError("Error: file \""+pedfile+"\" not found in current directory");
        	log.reportException(fnfe);
            return;
        } catch (IOException ioe) {
        	log.reportError("Error reading file \""+pedfile+"\"");
        	log.reportException(ioe);
        	return;
        }
	}

	public static void gwafFromParamters(String filename, Logger log) {
		Vector<String> params;
		String pedfile, mapfile, freqfile, rlinker, outfile;
		
		params = Files.parseControlFile(filename, "gwaf", new String[] {"plink.ped", "plink.map", "plink.frq", "leslie_lange.FHS.IBC.CEU.Rlinker", "gwaf.csv", "# set Rlinker file to null if none is needed"}, log);
		if (params != null) {
    		pedfile = params.remove(0);
    		mapfile = params.remove(0);
    		freqfile = params.remove(0);
    		rlinker = params.remove(0);
    		if (rlinker.toLowerCase().equals("null")) {
    			rlinker = null;
    		}
    		outfile = params.remove(0);
//    		for (int i = 0; i < files.length; i++) {
//    			line = params.elementAt(i).trim().split("[\\s]+");
//    			files[i] = line[0];
//    			for (int j = 1; j < line.length; j++) {
//        			if (line[j].startsWith("skip=") || line[j].startsWith("skips=")) {
//            			skips[i] = Integer.parseInt(line[j].split("=")[1]);
//        			} else if (line[j].equals(",")) {
//            			delimiters[i] = ",";
//        			} else if (line[j].equalsIgnoreCase("tab")) {
//            			delimiters[i] = "\t";
//        			} else {
//        				System.err.println("Error - invalid argument (\""+line[j]+"\") ");
//        			}
//				}
//			}
    		toGWAF(pedfile, mapfile, freqfile, rlinker, outfile);
		}
	}
	
	public static void addHeaderToPed(String root, String output, boolean split) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, batch_dna;
		Vector<String> v = new Vector<String>();
		int count;
		String[] markerNames;
		
		try {
			reader = Files.getAppropriateReader(root+".map");
			while (reader.ready()) {
				v.add(reader.readLine().trim().split("[\\s]+")[1]);
			}
			reader.close();
			
			markerNames = new String[v.size()];
			for (int i = 0; i < markerNames.length; i++) {
				markerNames[i] = v.elementAt(i);
			}

			reader = Files.getAppropriateReader(root+".ped");
			writer = Files.getAppropriateWriter(output);
			writer.print((split?"BATCH\t":"")+"FID\tIID\tFATHER\tMOTHER\tSEX\tAFFECTED");
			for (int i = 0; i < markerNames.length; i++) {
				writer.print("\t"+markerNames[i]+"A"+"\t"+markerNames[i]+"B");
			}
			writer.println();
			count = 0;
			while (reader.ready()) {
				count++;
				line = reader.readLine().trim().split("[\\s]+");
				if (line.length != markerNames.length*2+6) {
					System.err.println("Error - line "+count+" has "+line.length+" columns instead of the expected "+markerNames.length*2+6+" (#markers*2+6)");
					try {
						writer = new PrintWriter(new FileWriter("ERROR_NUMBER_OF_COLUMNS_DON'T_MATCHUP_FOR_THESE_LINES.TXT", true));
						writer.println(count);
						writer.close();
					} catch (Exception e) {
						System.err.println("Error writing to " + "ERROR_COULD_NOT_OPEN_"+root+".map_FILE.TXT");
						e.printStackTrace();
					}
				}
				for (int i = 0; i < line.length; i++) {
					if (i==0) {
						if (split) {
							batch_dna = line[i].split("_");
							if (batch_dna.length != 2) {
								try {
									writer = new PrintWriter(new FileWriter("ERROR_THESE_FIDS_AREN'T_CONCATENATED_BY_AN_UNDERSCORE.TXT", true));
									writer.println("line"+count+"\t"+line[i]);
									writer.close();
								} catch (Exception e) {
									System.err.println("Error writing to " + "ERROR_COULD_NOT_OPEN_"+root+".map_FILE.TXT");
									e.printStackTrace();
								}
								writer.print(line[i]+"\t"+line[i]);
							} else {
								writer.print(batch_dna[0]+"\t"+batch_dna[1]);
							}
						} else {
							writer.print(line[i]);
						}
					} else {
						writer.print("\t"+line[i]);
					}
				}
				writer.println();
			}
			reader.close();
			writer.close();			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String pedfile = "plink.ped";
		String mapfile = "plink.map";
		String freqfile = "plink.frq";
		String rlinker = null;
		String outfile = "gwaf.csv";

		String root = "plink";
		String output = "plink.txt";
		boolean split = false;

		String usage = "\n" +
		"gwas.CreateDatabaseFromPlink\n" +
		" *** see Launch for createDatabase, createScore, and createCountsMatrix (for simpleM)***\n" + 
		"  For convert toGWAF, use the following arguments:\n" + 
		"   (1) .ped filename (i.e. ped=" + pedfile + " (default))\n" + 
		"   (2) .map filename (i.e. map=" + mapfile + " (default))\n" + 
		"   (3) .frq filename (i.e. frq=" + freqfile + " (default))\n" + 
		"   (4) String->int Rlinker filename with CARe_ID FID IID (i.e. rlinker=" + rlinker + " (default))\n" + 
		"   (5) output filename (i.e. gwaf=" + outfile + " (default))\n" +
		" OR\n"+
		"   (1) root of plink file (i.e. root=" + root + " (default))\n" + 
		"   (2) name of output file (i.e. output=" + output + " (default))\n" + 
		"   (3) split the FID (i.e. -split (not the default))\n" + 
		"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("ped=")) {
				pedfile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("map=")) {
				mapfile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("frq=")) {
				freqfile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("rlinker=")) {
				rlinker = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("gwaf=")) {
				outfile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("root=")) {
				root = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("output=")) {
				output = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("-split")) {
				split = true;
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		
//		String dir = "D:\\BOSS\\GWAF\\";
//		pedfile = dir+"plink.ped";
//		mapfile = dir+"plink.map";
//		freqfile = dir+"plink.frq";
//		rlinker = null;
//		outfile = dir+"gwaf.csv";

		String dir = "D:/CARe/CARe_geno_data_and_misc/IBC/FHS/iSELECT/gwaf/";
		pedfile = dir+"plink.ped";
		mapfile = dir+"plink.map";
		freqfile = dir+"plink.frq";
		rlinker = dir+"leslie_lange.FHS.IBC.CEU.Rlinker";
		outfile = dir+"gwaf.csv";
		
		try {
//			createCountsMatrix("simpleM.crf", new Logger());
			if (pedfile != null && Files.exists(pedfile)) {
				toGWAF(pedfile, mapfile, freqfile, rlinker, outfile);
			} else if (root != null && Files.exists(root+".ped")) {
				addHeaderToPed(root, output, split);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
