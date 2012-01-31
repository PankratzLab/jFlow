// -Xms1536M -Xmx1536M
//
// used to filter out markers because of hwe, but since Haploview doesn't limit itself to controls, I truned this off (-hwcutoff 0)
// always filter on hwe before using as a reference
//
// Might need to set a minimum distance to previous tag (default 32 bp, but might consider 1 Mb for unlinked markers
//
// should the chromosomal positions be stored in LDdatabase and serialized to be used in subsequent calculations?
//
// still doesn't take into consideration distance (i.e. can't have two probes within 32 bp of each other
//
package gwas;

import java.io.*;
import java.util.*;

import common.*;
import filesys.*;

public class IndependentSNPs {
//	public static final float DEFAULT_PVAL_THRESHOLD = 5E-4;
	public static final float DEFAULT_PVAL_THRESHOLD = 1;
	public static final float DEFAULT_R2_THRESHOLD = 0.80f;
//	public static final float DEFAULT_R2_THRESHOLD = 0.05f;
	public static final String DEFAULT_DIR_ILLUMINA_SCORES = "C:\\Documents and Settings\\npankrat\\My Documents\\tWork\\Consortium\\selection_db\\";
//	public static final String DEFAULT_DIR_ILLUMINA_SCORES = null;
//	public static final String DEFAULT_FORCED_DATASET = "/home/npankrat/NCBI/HapMap/CEU_660/CEU_660_founders.bim";
	public static final String DEFAULT_FORCED_DATASET = null;
//	public static final String DEFAULT_FORCE_INCLUDE_BEFORE_TAGGING_IF_MEETS_CRITERIA = "preferredTags.dat";
	public static final String DEFAULT_FORCE_INCLUDE_BEFORE_TAGGING_IF_MEETS_CRITERIA = null;
	public static final String DEFAULT_FORCE_INCLUDE_AFTER_TAGGING_IF_MEETS_CRITERIA = "common_nonsynonymous+UTR.dat";
	public static final String DEFAULT_FORCE_INCLUDE_REGARDLESS = "missingGenotypes.dat";
	public static final String[] ILLUMINA_TARGET_COLUMNS = {"Locus_Name", "Failure_Message", "Final_Score"};
//	public static final float DEFAULT_SCORE_THRESHOLD = 0;
	public static final String DEFAULT_SCORE_THRESHOLDS = "0.4<1E-5;0.6";
	public static final float DEFAULT_SCORE_DIFF_THRESHOLD = 0.2f;
	public static final float DEFAULT_SCORE_CLASS_BUMP = 0;
	public static final String DEFAULT_LD_ROOT = LDdatabase.MASTER_HAPMAP_ROOT;
	public static final int DEFAULT_POSITION_REPORT_TYPE = LDdatabase.REPORT_LAST;
	public static final int DEFAULT_LD_COMPARISON_TYPE = LDdatabase.COMP_FIRST;
	
	public static final int NUM_CHROMOSOMES = 27;
	public static final double[] SCORE_BIN_FAIL_RATES = {1.0, 0.55, 0.45, 0.45, 0.30, 0.12, 0.10};

	public static void findOptimalSet(String dir, String filename, String outputRoot, int numSNPs, float pval_threshold, String ldRoots, int position_reportType, int r2_compType, float r2_threshold, String directoryOfIlluminaScores, String scoreThresholds, float scoreDiffThreshold, float scoreClassBump, String filteringDataset, String forceBeforeFile, String forceAfterFile, String forceRegardlessFile, Logger log) {
		BufferedReader reader;
        PrintWriter writer;
        Vector<String> tags, untaggedTags, checkTags;
        long time;
        ResultSet results;
        IntVector[] chrIVs;
        IntVector iv, beforeIndicesVector;
        String[] markerNames, superset, subset, forcedKeys;
        float[] pvalues, pvals;
        LDdatabase[] lddbs;
        int[] positions, keys;
        Hashtable<String,Hashtable<String,String>> chrHashes;
        Hashtable<String,String> forceBefore, forceAfter, forceRegardless;
        String trav;
        int index, offset;
//        StringLDdb chrLDdb;
        LongLDdb[] chrLDdbs;
        float r2, pval;
        String[] files, line, dirs;
        Hashtable<String,Float> scores;
        int[] indices;
        Vector<String> missingMarkers, missingIlluminaValues;
        Hashtable<String, String> chrHash, merges, allMissingMarkers;
        float score, travScore;
        int bestIndex;
        double bestDiff, diff;
        int beforeIndex;
        int[] beforeIndices;
        float[] pThreshs, scoreThreshs;
        boolean passes;
        String format;
        
        if (filename.indexOf(";") > 0) {
        	format = filename.substring(filename.indexOf(";")+1);
        	filename = filename.substring(0, filename.indexOf(";"));
        } else {
        	format = "pvalues";
        }
        if (!new File(dir+filename+".rset").exists()) {
    		prep(dir+filename, format, log);
    	}
    	
        log.report("Reading results file: "+filename, false, true);
        time = new Date().getTime();
		results = ResultSet.load(dir+filename+".rset", false, false);
		log.report("...finished in "+ext.getTimeElapsed(time));
		
        scores = new Hashtable<String,Float>();
        merges = new Hashtable<String,String>();
        if (directoryOfIlluminaScores != null) {
        	dirs = directoryOfIlluminaScores.split(";");
        	for (int j = 0; j<dirs.length; j++) {
        		dirs[j] = ext.verifyDirFormat(dirs[j]);
            }
        	for (int j = 0; j<dirs.length; j++) {
            	files = Files.list(dirs[j], ".csv", false);
            	for (int i = 0; i<files.length; i++) {
            		try {
    	                reader = new BufferedReader(new FileReader(dirs[j]+files[i]));
    	                line = new String[] {""};
    	                while (reader.ready() && !line[0].startsWith("Locus_Name")) {
    	                	line = reader.readLine().trim().split(",");
    	                }
    	                indices = ext.indexFactors(ILLUMINA_TARGET_COLUMNS, line, false, null, false, false);
    	                if (!reader.ready()) {
    	                	log.reportError("Error - failed to parse '"+files[i]+"'; no line started with 'Locus_Name'");
    	                } else if (indices[1] == -1 && indices[2] == -1) {
    		                log.reportError("Error - Illumina files must have one of the following columns: "+Array.toStr(Array.subArray(ILLUMINA_TARGET_COLUMNS, 1), ", "));
    		                log.reportError("      - only found: "+Array.toStr(line, ", "));
    	                } else {
    		                while (reader.ready()) {
    		                	line = reader.readLine().trim().split(",", -1);
    		                	if (indices[1] >= 0) {
    		                		if (line[indices[1]].startsWith("merged to ")) {
    		                			merges.put(line[indices[0]], line[indices[1]].substring(10));
    		                		} else if (line[indices[1]].endsWith("not supported") || line[indices[1]].endsWith("assembly") || line[indices[1]].endsWith("alternative")) {
    		                			scores.put(line[indices[0]], Float.parseFloat("-8"));
    		                		} else {
    		                			log.reportError("Error - unknown Failure_Message (assuming it's bad): "+line[indices[1]]);
    		                			scores.put(line[indices[0]], Float.parseFloat("-7"));
    		                		}
    		                	} else {
    		                		if (line[indices[2]].equals("")) {
    		                			scores.put(line[indices[0]], Float.parseFloat("0.05"));
    		                		} else {
    		                			if (scores.containsKey(line[indices[0]])) {
    		                				if (Math.min(scores.get(line[indices[0]]).floatValue(), Float.parseFloat(line[indices[2]])) < -900) {
        		                				scores.put(line[indices[0]], Float.parseFloat("-999"));
    		                				} else {
    		                					scores.put(line[indices[0]], Math.max(scores.get(line[indices[0]]).floatValue(), Float.parseFloat(line[indices[2]])));
    		                				}
    		                			} else {
    		                				scores.put(line[indices[0]], Float.parseFloat(line[indices[2]]));
    		                			}
    		                		}
    		                	}
    		                }
    	                }
    	                reader.close();
                    } catch (FileNotFoundException fnfe) {
    	                log.reportError("Error: file \""+dirs[j]+files[i]+"\" not found in current directory");
    	                return;
                    } catch (IOException ioe) {
    	                log.reportError("Error reading file \""+dirs[j]+files[i]+"\"");
    	                return;
                    }
                }
            	if (files.length == 0) {
            		log.reportError("Error - no illumina scores files found in "+dirs[j]+"; no warnings will be provided");
            		dirs[j] = null;
            	}
            }
        }
		
		lddbs = LDdatabase.multiLDdbs(ldRoots.split(";"), log);
		
		iv = new IntVector();
		
		markerNames = results.getMarkerNames();
		pvalues = results.getPvals();
		
		for (int i = 0; i<pvalues.length; i++) {
			if (pvalues[i] < pval_threshold) {
				iv.add(i);
			}
        }

		superset = new String[iv.size()];
        for (int j = 0; j<iv.size(); j++) {
        	superset[j] = markerNames[iv.elementAt(j)];	
        }
        LDdatabase.multiUpdate(lddbs, superset, ext.replaceDirectoryCharsWithUnderscore(dir+filename, 2));
        
		if (forceBeforeFile != null && new File(dir+forceBeforeFile).exists()) {
			log.report("Forcing those SNPs in '"+forceBeforeFile+"' to be included as tags if p-value met");
			forceBefore = HashVec.loadFileToHashNull(dir+forceBeforeFile, false);
		} else {
			log.report((forceBeforeFile==null?"No ":"No file named '"+forceBeforeFile+"'; no ")+"SNPs will be forced to be included as tags");
			forceBefore = new Hashtable<String,String>();
		}
		if (forceAfterFile != null && new File(dir+forceAfterFile).exists()) {
			log.report("Forcing those SNPs in '"+forceAfterFile+"' to be included as tags if p-value met");
			forceAfter = HashVec.loadFileToHashNull(dir+forceAfterFile, false);
		} else {
			log.report((forceAfterFile==null?"No ":"No file named '"+forceAfterFile+"'; no ")+"SNPs will be forced to be included as tags after tagging");
			forceAfter = new Hashtable<String,String>();
		}
		if (forceRegardlessFile != null && new File(dir+forceRegardlessFile).exists()) {
			log.report("Forcing those SNPs in '"+forceRegardlessFile+"' to be included as tags regardless of p-value");
			forceRegardless = HashVec.loadFileToHashNull(dir+forceRegardlessFile, false);
			forcedKeys = HashVec.getKeys(forceRegardless);
		} else {
			log.report((forceRegardlessFile==null?"No ":"No file named '"+forceRegardlessFile+"'; no ")+"additional SNPs will be forced regardless of p-value");
			forceRegardless = new Hashtable<String,String>();
			forcedKeys = new String[0];
		}
		
		line = scoreThresholds.trim().split(";");
		scoreThreshs = new float[line.length];
		pThreshs = new float[line.length];
		for (int i = 0; i<line.length; i++) {
			scoreThreshs[i] = Float.parseFloat(line[i].split("<")[0]);
			if (line[i].split("<").length > 1) {
				pThreshs[i] = Float.parseFloat(line[i].split("<")[1]);
			} else {
				pThreshs[i] = 1;
			}
        }
		
		chrHashes = LDdatabase.getChrHashes(lddbs);
		tags = new Vector<String>();
		missingIlluminaValues = new Vector<String>();
		allMissingMarkers = new Hashtable<String,String>();
    	missingMarkers = new Vector<String>();
    	chrIVs = IntVector.newIntVectors(27);
        for (int j = 0; j<iv.size(); j++) {
			score = getScore(superset[j], scores, merges, missingIlluminaValues, log);
			passes = false;
			for (int i = 0; i<scoreThreshs.length; i++) {
				if (score >= scoreThreshs[i] && pvalues[iv.elementAt(j)] < pThreshs[i]) {
					passes = true;
				}
            }
//			if (directoryOfIlluminaScores == null || score > scoreThreshold) {
//			if (scores.size() == 0 || score > scoreThreshold || (pvalues[iv.elementAt(j)] < 1E-4 && score > 0.4)) {
			if (scores.size() == 0 || passes) {
				trav = LDdatabase.multiHashCheck(chrHashes, superset[j], position_reportType, log);
				if (trav == null) {
//					log.reportError("Error - the chrHash derived from the marker file does not contain marker "+markerNames[ivs[i].elementAt(j)]);
					HashVec.addIfAbsent(superset[j], missingMarkers);
					allMissingMarkers.put(superset[j], "");
	    			chrIVs[0].add(iv.elementAt(j));
				} else {
	    			chrIVs[Byte.parseByte(trav.split("[\\s]+")[0])].add(iv.elementAt(j));
				}
			}
        }

        for (int chr = 0; chr<NUM_CHROMOSOMES; chr++) {
        	subset = new String[chrIVs[chr].size()];
        	pvals = new float[chrIVs[chr].size()];
        	positions = new int[chrIVs[chr].size()];
            beforeIndicesVector = new IntVector();
            for (int j = 0; j<positions.length; j++) {
            	subset[j] = markerNames[chrIVs[chr].elementAt(j)];
            	pvals[j] = pvalues[chrIVs[chr].elementAt(j)];
				trav = LDdatabase.multiHashCheck(chrHashes, markerNames[chrIVs[chr].elementAt(j)], position_reportType, log);
    			if (trav == null) {
    				positions[j] = j;
    			} else {
    				positions[j] = Integer.parseInt(trav.split("[\\s]+")[1]);
    			}
    			if (forceBefore.containsKey(subset[j])) {
    				forceBefore.put(subset[j], pvals[j]+"");
    				beforeIndicesVector.add(j);
    			}
    			if (forceAfter.containsKey(subset[j])) {
    				forceAfter.put(subset[j], pvals[j]+"");
    			}
    			if (forceRegardless.containsKey(subset[j])) {
    				forceRegardless.put(subset[j], pvals[j]+"");
    			}
            }
            
            keys = Sort.quicksort(positions);
            subset = Sort.putInOrder(subset, keys);
            positions = Sort.putInOrder(positions, keys);
            pvals = Sort.putInOrder(pvals, keys);

            beforeIndices = beforeIndicesVector.toArray();
            beforeIndex = 0;            
            chrLDdbs = LDdatabase.getChrLDdbs(lddbs, chr);
            while (pvals.length > 0 && Array.min(pvals) < 2) {
            	index = -1;
            	pval = -1;
            	
            	if (beforeIndex < beforeIndices.length) {
            		index = beforeIndices[beforeIndex];
        			pval = pvals[index]; 
        			pvals[index] = 3;
            		beforeIndex++;
            	} else {
            		index = Array.indexOfMin(pvals);
        			pval = pvals[index]; 
                	score = getScore(subset[index], scores, merges, missingIlluminaValues, log);
                	bestDiff = scoreDiffThreshold;
                	bestIndex = -1;
                	for (int j = 0; j<2; j++) {
                    	offset = j==0?-1:1;
                    	for (int k = 1; index+k*offset >= 0 && 
                    					index+k*offset < pvals.length && 
                    					Math.abs(positions[index] - positions[index+k*offset]) < LDdatabase.BP_LIMIT; k++) {
                    		if (pvals[index+k*offset] < 2) {
                    			r2 = LDdatabase.multiCheck(chrLDdbs, subset[index], subset[index+k*offset], r2_compType, log);
    	                		if (r2 == LDdatabase.MISSING_INFO) {
    	                			if (!allMissingMarkers.containsKey(subset[index]) || !allMissingMarkers.containsKey(subset[index+k*offset])) {
    	                				log.reportError("Error - missing LD info for "+subset[index]+"/"+subset[index+k*offset]+" pair");
    	                			}
    	                		} else if (r2 > r2_threshold){
    	                			travScore = getScore(subset[index+k*offset], scores, merges, missingIlluminaValues, log);
    	                			diff = (travScore - score) / (-1*Math.log10(pval) - -1*Math.log10(pvals[index+k*offset])) + scoreClassBump*(getScoreClass(travScore)-getScoreClass(score));
    	                			if (travScore > score && diff  > bestDiff) {
    	                				bestDiff = diff;
    	                				bestIndex = index+k*offset;
    	                			}
    	                		}
                    		}
                        }
                	}
                	if (bestIndex >= 0) {
//                		log.report("From "+subset[index]+" to "+subset[bestIndex]);
                		index = bestIndex;
                	}
        			pval = pvals[index]; 
        			pvals[index] = 2;
            	}
            	tags.add(subset[index]+"\t"+chr+"\t"+positions[index]+"\t"+pval+"\t"+getScore(subset[index], scores, merges, missingIlluminaValues, log));
            	
            	for (int j = 0; j<2; j++) {
                	offset = j==0?-1:1;
                	for (int k = 1; index+k*offset >= 0 && 
                					index+k*offset < pvals.length && 
                					Math.abs(positions[index] - positions[index+k*offset]) < LDdatabase.BP_LIMIT; k++) {
                		if (pvals[index+k*offset] < 2) {
                			r2 = LDdatabase.multiCheck(chrLDdbs, subset[index], subset[index+k*offset], r2_compType, log);
	                		if (r2 == LDdatabase.MISSING_INFO) {
	                			if (!allMissingMarkers.containsKey(subset[index]) || !allMissingMarkers.containsKey(subset[index+k*offset])) {
	                				log.reportError("Error - missing LD info for "+subset[index]+"/"+subset[index+k*offset]+" pair");
	                			}
	                		} else if (r2 > r2_threshold){
	                			pvals[index+k*offset] = 2;
	                		}
                		}
                    }
            	}
            }
            for (int j = 0; j<subset.length; j++) {
            	if (forceAfter.containsKey(subset[j]) && Double.parseDouble(forceAfter.get(subset[j])) < pval_threshold) {
            		HashVec.addIfAbsent(subset[j]+"\t"+chr+"\t"+positions[j]+"\t"+forceAfter.get(subset[j])+"\t"+getScore(subset[j], scores, merges, missingIlluminaValues, log), tags);
            	}
            }
        }
        
        for (int j = 0; j<forcedKeys.length; j++) {
			trav = LDdatabase.multiHashCheck(chrHashes, forcedKeys[j], position_reportType, log);
        	HashVec.addIfAbsent(forcedKeys[j]+"\t"+(trav == null?".\t.":trav)+"\t"+forceRegardless.get(forcedKeys[j])+"\t"+getScore(forcedKeys[j], scores, merges, missingIlluminaValues, log), tags);
        }
        if (missingMarkers.size() > 0) {
        	log.reportError("Error - Missing "+missingMarkers.size()+" markers for threshold "+pval_threshold);
        	Files.writeList(Array.toStringArray(missingMarkers), dir+(outputRoot == null?pval_threshold:outputRoot)+"_missingValues.txt");
        }
        if (missingIlluminaValues.size() > 0) {
        	log.reportError("Error - Missing Illumina design scores for "+missingIlluminaValues.size()+" markers");
        	missingIlluminaValues.insertElementAt("Locus_Name", 0);
        	Files.writeList(Array.toStringArray(missingIlluminaValues), dir+(outputRoot == null?pval_threshold:outputRoot)+"_ill.txt");
        }
        
		untaggedTags = new Vector<String>();
		checkTags = new Vector<String>();
        if (filteringDataset != null) {
        	log.report("Loading which markers are present in the "+filteringDataset+" dataset...", false, true);
			chrHash = SnpMarkerSet.loadSnpMarkerSetToChrHash(filteringDataset);
			log.report("done");
            for (int j = 0; j<tags.size(); j++) {
            	trav = tags.elementAt(j).split("[\\s]+")[0];
            	if (!chrHash.containsKey(trav)) {
            		untaggedTags.add(trav+"\t1");
            	} else {
            		checkTags.add(trav+"\t1");
            	}
            }
        }
        
		try {
	        writer = new PrintWriter(new FileWriter(dir+(outputRoot == null?pval_threshold:outputRoot)+"_summary.xln"));
	        writer.println("\tp<"+pval_threshold);
	        writer.println("Total SNPs meeting threshold"+"\t"+iv.size());
	        writer.println("Independent index SNPs meeting threshold"+"\t"+tags.size());

	        if (filteringDataset != null) {
		        writer.println("Independent index SNPs meeting threshold on filtering array"+"\t"+checkTags.size());
		        writer.println("Independent index SNPs meeting threshold not on filtering array"+"\t"+untaggedTags.size());
			    log.report(untaggedTags.size()+" SNPs");
	        }
	        
	        writer.close();
        } catch (Exception e) {
	        log.reportError("Error writing to "+(outputRoot == null?pval_threshold:outputRoot)+"_summary.xln");
	        e.printStackTrace();
        }
        
        Files.writeList(trimList(tags, numSNPs), dir+(outputRoot == null?pval_threshold:outputRoot)+"_tags.xln");
        if (filteringDataset != null) {
            Files.writeList(Array.toStringArray(checkTags), dir+(outputRoot == null?pval_threshold:outputRoot)+"_checkArray_tags.xln");
            Files.writeList(Array.toStringArray(untaggedTags), dir+(outputRoot == null?pval_threshold:outputRoot)+"_untagged_tags.xln");
        }
	}

	public static String[] trimList(Vector<String> tags, int numSNPs) {
		String[] finalList;
		double[] values;
		int[] order;
		
		if (tags.size() < numSNPs) {
			return Array.toStringArray(tags);
		}
		
		values = new double[tags.size()];
		for (int i = 0; i<values.length; i++) {
			values[i] = Double.parseDouble(tags.elementAt(i).split("[\\s]+")[2]);
        }
		order = Sort.quicksort(values);
		
		finalList = new String[numSNPs];
		for (int i = 0; i<numSNPs; i++) {
			finalList[i] = tags.elementAt(order[i]);
        }
		
		return finalList;
	}
	
	public static int getScoreClass(float score) {
		if (score > 1.0) {
			return 6;
		} else if (score > 0.8) {
			return 5;
		} else if (score > 0.6) {
			return 4;
		} else if (score > 0.4) {
			return 3;
		} else if (score > 0.2) {
			return 2;
		} else if (score > 0) {
			return 1;
		} else {
			return 0;
		}
	}
	
	public static float getScore(String element, Hashtable<String,Float> scores, Hashtable<String,String> merges, Vector<String> missingValues, Logger log) {
		if (scores.containsKey(element)) {
			return scores.get(element).floatValue();
		} else if (merges.containsKey(element)) {
			if (scores.containsKey(merges.get(element))) {
				return scores.get(merges.get(element)).floatValue();
			} else {
				HashVec.addIfAbsent(merges.get(element), missingValues);
				log.reportError("Error - "+element+" merged with "+merges.get(element)+"; so we'll need to request info for it");
				return -9;
			}
		} else {
			HashVec.addIfAbsent(merges.get(element), missingValues);
			return -9;
		}
	}
	
	public static void prep(String filename, String format, Logger log) {
		log.report("Prepping file using "+format+" format...");
        long time = new Date().getTime();
        if (format.equals("pvalues")) {
    		new ResultSet(filename, ResultSet.PVALUES_ONLY).serialize(filename+".rset");
        } else if (format.equals("metal")) {
    		new ResultSet(filename, ResultSet.METAL_TBL_FORMAT).serialize(filename+".rset");
        } else if (format.equals("generic")) {
    		new ResultSet(filename, ResultSet.GENERIC_FORMAT).serialize(filename+".rset");
        } else {
        	log.reportError("Error - '"+format+"' is not a recognized file format for a set of results (expecting \"pvalues\", \"metal\", or \"generic\")");
        }
		log.report("Finished prepping in "+ext.getTimeElapsed(time));
	}
	
	public static String[] findUnique(String[] inThese, String[] butNotInThese) {
        Vector<String> v = new Vector<String>();
        
        for (int i = 0; i<inThese.length; i++) {
        	if (ext.indexOfStr(inThese[i], butNotInThese) == -1) {
        		v.add(inThese[i]);
        	}
        }
		
		return Array.toStringArray(v);
	}
	
	public static Vector<String> checkCoverage(String hitTags, String arrayTags, Logger log) {
		BufferedReader reader;
        String[] line;
        Hashtable<String,String> hash;
        Vector<String> v = new Vector<String>();
        int numAlleles;
        
        hash = HashVec.loadFileToHashString(hitTags, 0, new int[] {0}, "\t", false);
        try {
	        reader = new BufferedReader(new FileReader(arrayTags));
        	line = reader.readLine().trim().split("[\\s]+");
        	if (!line[0].equals("#captured") || !line[4].equals("alleles")) {
        		log.reportError("Error - Haploview .TAGS file has changed format and needs to be addressed");
        		System.exit(1);
        	}
        	numAlleles = Integer.parseInt(line[3]);
        	if (!reader.readLine().startsWith("#captured") || !reader.readLine().startsWith("#using")) {
        		log.reportError("Error - Haploview .TAGS file has changed format and needs to be addressed");
        		System.exit(1);
        	}
        	for (int i = 0; i<numAlleles; i++) {
	        	line = reader.readLine().trim().split("[\\s]+");
	        	if (hash.containsKey(line[0]) && line.length == 1) {
	        		v.add(line[0]);
	        	}
	        }
	        reader.close();
        } catch (FileNotFoundException fnfe) {
	        log.reportError("Error: file \""+arrayTags+"\" not found in current directory");
	        System.exit(1);
        } catch (IOException ioe) {
	        log.reportError("Error reading file \""+arrayTags+"\"");
	        System.exit(2);
        }
        
		return v;
	}
	
	public static void compareSelectionParameters(int numSNPs, String dir, String filename, float pval_threshold, String ldRoots, int position_reportType, int r2_compType, float r2_threshold, String directoryOfIlluminaScores, String filteringDataset, String forceBeforeFile, String forceAfterFile, String forceRegardlessFile, Logger log) {
        PrintWriter writer;
        String[] line;
        
		// scoreThreshold, scoreDiffThreshold, scoreClassBump
		float[][] params = new float[][] {
//				{Float.NEGATIVE_INFINITY, Float.POSITIVE_INFINITY, 0},
//				{0, Float.POSITIVE_INFINITY, 0},
//				{0, 1, 0},
//				{0, 0.5f, 0},
//				{0, 0.2f, 0},
//				{0, 1, 0.1f},
//				{0, 0.5f, 0.1f},
//				{0, 0.2f, 0.1f},
				{0.4f, 0.5f, 0.1f},
				{0.6f, 0.5f, 0.1f},
//				{0.8f, 0.5f, 0.1f},
//				{1, 0.5f, 0.1f},
//				{0.4f, Float.POSITIVE_INFINITY, 0},
//				{0.6f, Float.POSITIVE_INFINITY, 0},
//				{0.8f, Float.POSITIVE_INFINITY, 0},
//				{1, Float.POSITIVE_INFINITY, 0},
				};
		
		double sumPval;
		double effSumPval;
		int[] bins;
		double failRate;
		int bin;
		String[][] data;
		int[] keys;
		
		try {
	        writer = new PrintWriter(new FileWriter(dir+"comparison.xln"));
	        writer.println("scoreThreshold\tscoreDiffThreshold\tscoreClassBump\t\tPredictedFailRate\tSumLogPvals\tEffectiveSum\t%Tagged\t#failDesign\t#0<score<0.2\t#0.2<=score<0.4\t#0.4<=score<0.6\t#0.6<=score<0.8\t#0.8<=score<1.0\t#score=1.1");
	        for (int i = 0; i<params.length; i++) {
    			findOptimalSet(dir, filename, params[i][0]+"_"+params[i][1]+"_"+params[i][2], numSNPs, pval_threshold, ldRoots, position_reportType, r2_compType, r2_threshold, directoryOfIlluminaScores, params[i][0]+"", params[i][1], params[i][2], filteringDataset, forceBeforeFile, forceAfterFile, forceRegardlessFile, log);
	        	data = HashVec.loadFileToStringMatrix(dir+params[i][0]+"_"+params[i][1]+"_"+params[i][2]+"_tags.xln", false, new int[] {0, 1, 2, 3}, false);
	        	keys = Sort.quicksort(Array.toDoubleArray(Matrix.extractColumn(data, 2)));

        		bins = new int[7];
                failRate = -1;
        		sumPval = effSumPval = 0;
	        	for (int j = 0; j<keys.length && j<numSNPs; j++) {
                	line = data[keys[j]];
                	bin = getScoreClass(Float.parseFloat(line[3]));
                	bins[bin]++;
                	sumPval += bin==0?0:(-1*Math.log10(Double.parseDouble(line[2]))); 
                	effSumPval += -1*Math.log10(Double.parseDouble(line[2]))*SCORE_BIN_FAIL_RATES[bin];
                }
                failRate = 0;
                for (int j = 0; j<bins.length; j++) {
                	failRate += SCORE_BIN_FAIL_RATES[j]*bins[j];
                }
                failRate /= Array.sum(bins);
            	writer.println(Array.toStr(Array.toDoubleArray(params[i]), 2, 2, "\t")+"\t\t"+failRate+"\t"+sumPval+"\t"+effSumPval+"\t"+"."+"\t"+Array.toStr(bins));
            }
	        writer.close();
        } catch (Exception e) {
	        log.reportError("Error writing to "+dir+"comparison.xln");
	        e.printStackTrace();
        }
		
	}
	
	public static void selectFromParameters(String filename, Logger log) {
        String trav;
        String resultsFile = "hits.txt", outputRoot = "tags", ldRoots = DEFAULT_LD_ROOT;
        float pval_threshold = DEFAULT_PVAL_THRESHOLD, r2_threshold = DEFAULT_R2_THRESHOLD;
        String directoryOfIlluminaScores = null;
        String scoreThresholds = "0";
        float scoreDiffThreshold = 0.5f, scoreClassBump = 0.1f;
        String filteringDataset = null, forceBeforeFile = null, forceAfterFile = null, forceRegardlessFile = null;
        int numSNPs = Integer.MAX_VALUE;
	    int position_reportType = DEFAULT_POSITION_REPORT_TYPE;
	    int r2_compType = DEFAULT_LD_COMPARISON_TYPE;
        Vector<String> params;

		params = Files.parseControlFile(filename, "indep", new String[] {"resultsFile=hits.txt", "outputRoot=tags", "pval_threshold=0.0001", "# ldRoots should be separated by a semicolon", "ldRoots=/home/npankrat/NCBI/HapMap/CEU_founders/CEU_founders", "# reportType options: 1=first valid position in order of ldRoots, 2=last valid position, 3=requires unanimous, otherwise missing, 4=same as 3 but reports error", "reportType=2", "# compType options: 1=first valid r2 value, 2=minimum r2 value, 3=maximum r2 value", "compType=1", "r2_threshold=0.80", "directoryOfIlluminaScores=null", "scoreThresholds=0.4<1E-5;0.6", "scoreDiffThreshold=0.50", "scoreClassBump=0.10", "filteringDataset=null", "forceBeforeFile=null", "forceAfterFile=null", "forceRegardlessFile=null", "#numSNPs=100"}, log);
		if (params != null) {
			for (int i = 0; i < params.size(); i++) {
	        	trav = params.elementAt(i).trim();
	        	if (trav.startsWith("resultsFile=")) {
	        		resultsFile = ext.parseStringArg(trav, null);
    		    } else if (trav.startsWith("outputRoot=")) {
    			    outputRoot = ext.parseStringArg(trav, null);
    		    } else if (trav.startsWith("pval_threshold=")) {
    		    	pval_threshold = ext.parseFloatArg(trav);
    		    } else if (trav.startsWith("ldRoots=")) {
    		    	ldRoots = ext.parseStringArg(trav, null);
    		    } else if (trav.startsWith("reportType=")) {
    		    	position_reportType = ext.parseIntArg(trav);
    		    } else if (trav.startsWith("compType=")) {
    		    	r2_compType = ext.parseIntArg(trav);
    		    } else if (trav.startsWith("r2_threshold=")) {
    			    r2_threshold = ext.parseFloatArg(trav);
    		    } else if (trav.startsWith("directoryOfIlluminaScores=")) {
    		    	directoryOfIlluminaScores = ext.parseStringArg(trav, null);
    		    } else if (trav.startsWith("scoreThresholds=")) {
    		    	scoreThresholds = ext.parseStringArg(trav, "0");
    		    } else if (trav.startsWith("scoreDiffThreshold=")) {
    		    	scoreDiffThreshold = ext.parseFloatArg(trav);
    		    } else if (trav.startsWith("scoreClassBump=")) {
    		    	scoreClassBump = ext.parseFloatArg(trav);
    		    } else if (trav.startsWith("filteringDataset=")) {
    			    filteringDataset = ext.parseStringArg(trav, null);
    		    } else if (trav.startsWith("forceBeforeFile=")) {
    			    forceBeforeFile = ext.parseStringArg(trav, null);
    		    } else if (trav.startsWith("forceAfterFile=")) {
    		    	forceAfterFile = ext.parseStringArg(trav, null);
    		    } else if (trav.startsWith("forceRegardlessFile=")) {
    			    forceRegardlessFile = ext.parseStringArg(trav, null);
    		    } else if (trav.startsWith("numSNPs=")) {
    			    numSNPs = ext.parseIntArg(trav);
    		    } else if (!trav.startsWith("#")){
    		    	log.reportError("Error - don't know what to do with argument: "+trav);
    		    }
	        }

	        findOptimalSet(ext.parseDirectoryOfFile(filename), resultsFile, outputRoot, numSNPs, pval_threshold, ldRoots, position_reportType, r2_compType, r2_threshold, directoryOfIlluminaScores, scoreThresholds, scoreDiffThreshold, scoreClassBump, filteringDataset, forceBeforeFile, forceAfterFile, forceRegardlessFile, log);
//	        findOptimalSet("", resultsFile, outputRoot, numSNPs, pval_threshold, ldRoots, position_reportType, r2_compType, r2_threshold, directoryOfIlluminaScores, scoreThresholds, scoreDiffThreshold, scoreClassBump, filteringDataset, forceBeforeFile, forceAfterFile, forceRegardlessFile);
		}
        
	}
	
	public static void main(String[] args) {
	    int numArgs = args.length;
	    float r2_threshold = DEFAULT_R2_THRESHOLD;
	    String ldRoots = DEFAULT_LD_ROOT;
	    float pval_threshold = DEFAULT_PVAL_THRESHOLD;
	    String outputRoot = null;
	    String dirIlluminaScores = DEFAULT_DIR_ILLUMINA_SCORES;
	    String filteringDataset = DEFAULT_FORCED_DATASET;
	    String forceBeforeFile = DEFAULT_FORCE_INCLUDE_BEFORE_TAGGING_IF_MEETS_CRITERIA;
	    String forceAfterFile = DEFAULT_FORCE_INCLUDE_AFTER_TAGGING_IF_MEETS_CRITERIA;
	    String forceRegardlessFile = DEFAULT_FORCE_INCLUDE_REGARDLESS;
	    String scoreThresholds = DEFAULT_SCORE_THRESHOLDS;
	    float scoreDiffThreshold = DEFAULT_SCORE_DIFF_THRESHOLD;
	    float scoreClassBump = DEFAULT_SCORE_CLASS_BUMP;
	    boolean compare = false;
	    int numSNPs = Integer.MAX_VALUE;
	    int position_reportType = DEFAULT_POSITION_REPORT_TYPE;
	    int r2_compType = DEFAULT_LD_COMPARISON_TYPE;
	    
//	    String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\UMN\\Folson\\VTE_meta_analysis\\finalAnalysis\\06 tags\\";
//	    String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\UMN\\Folson\\VTE_meta_analysis\\finalAnalysis\\17 after dealing with CHS\\";
//	    String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\UMN\\Folson\\VTE_meta_analysis\\finalAnalysis\\11 X chromosome\\";
//	    String filename = "META_ANALYSIS_beta_se1.tbl";
//	    String filename = "META_ANALYSIS_beta_se1_parsed_noMissing.xln";
//	    String filename = "META_ANALYSIS_beta_se1_parsed_noMissing_noMAF_LTE_0.278.xln";
//	    String filename = "hits_wMAF_GTE0.024_described.xln";
	    
	    String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\tWork\\Consortium\\analysisOfImputation\\Aff_AAO_combo\\";
	    String filename = "minComputed.xln";

//	    String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\tWork\\Consortium\\analyses\\NGRC\\IlluminaAshk\\";
//	    String filename = "NGRC_Ashk_pvals.xln";

//	    String dir = "C:\\Documents and Settings\\npankrat\\My Documents\\tWork\\Consortium\\Fst\\Discovery\\TopFst\\";
//	    String filename = "TopFst.txt";

	    String usage = "\n"+
	    "gwas.IndependentSNPs requires 0-1 arguments\n"+
	    "   (1) directory (i.e. dir="+dir+" (default))\n"+
	    "   (2) filename (i.e. file="+filename+" (default))\n"+
	    "   (3) roots of PLINK files to test for LD (i.e. ldRoots="+ldRoots+" (default; separate roots with semicolon))\n"+
	    "   (4) r^2 threshold (i.e. r2="+r2_threshold+" (default))\n"+
	    "   (5) (optional) directory with Illumina score .csv files (i.e. scores="+dirIlluminaScores+" (default))\n"+
	    "   (6) (optional) thresholds for Illumina scores (i.e. scoreThresholds="+scoreThresholds+" (default))\n"+
	    "   (7) (optional) difference in score_diff/p_diff before swap (i.e. scoreDiffThreshold="+scoreDiffThreshold+" (default))\n"+
	    "   (8) (optional) bump from a categorical class difference (i.e. scoreClassBump="+scoreClassBump+" (default))\n"+
	    "   (9) (optional) dataset to filter tags (i.e. filter="+filteringDataset+" (default))\n"+
	    "   (10) (optional) file of SNPs to force before tagging [if they meet minimum p-value] (i.e. before="+forceBeforeFile+" (default))\n"+
	    "   (11) (optional) file of SNPs to force after tagging [if they meet minimum p-value] (i.e. after="+forceAfterFile+" (default))\n"+
	    "   (12) (optional) file of SNPs to force regardless of p-value (i.e. regardless="+forceRegardlessFile+" (default))\n"+
	    "   (13) (optional) compare multiple models (i.e. -compare (not the default))\n"+
	    "   (14) (optional) number of SNPs (i.e. snps="+numSNPs+" (default))\n"+
	    "   (15) (optional) marker position report type (i.e. reportType="+position_reportType+" (default; options: ))\n"+
	    "   (16) (optional) LD comparison type (i.e. compType="+r2_compType+" (default; options: ))\n"+
	    "";

	    for (int i = 0; i<args.length; i++) {
		    if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
		    	System.err.println(usage);
			    System.exit(1);
		    } else if (args[i].startsWith("dir=")) {
			    filename = ext.parseStringArg(args[i], "");
			    numArgs--;
		    } else if (args[i].startsWith("file=")) {
			    filename = args[i].split("=")[1];
			    numArgs--;
		    } else if (args[i].startsWith("ldRoots=")) {
		    	ldRoots = ext.parseStringArg(args[i], null);
			    numArgs--;
		    } else if (args[i].startsWith("reportType=")) {
		    	position_reportType = ext.parseIntArg(args[i]);
			    numArgs--;
		    } else if (args[i].startsWith("compType=")) {
		    	r2_compType = ext.parseIntArg(args[i]);
			    numArgs--;
		    } else if (args[i].startsWith("r2=")) {
			    r2_threshold = Float.parseFloat(args[i].split("=")[1]);
			    numArgs--;
		    } else if (args[i].startsWith("scores=")) {
		    	dirIlluminaScores = ext.parseStringArg(args[i], null);
			    numArgs--;
		    } else if (args[i].startsWith("filter=")) {
			    filteringDataset = ext.parseStringArg(args[i], null);
			    numArgs--;
		    } else if (args[i].startsWith("scoreThresholds=")) {
		    	scoreThresholds = ext.parseStringArg(args[i], "0");
			    numArgs--;
		    } else if (args[i].startsWith("scoreDiffThreshold=")) {
		    	scoreDiffThreshold = ext.parseFloatArg(args[i]);
			    numArgs--;
		    } else if (args[i].startsWith("scoreClassBump=")) {
		    	scoreClassBump = ext.parseFloatArg(args[i]);
			    numArgs--;
		    } else if (args[i].startsWith("before=")) {
			    forceBeforeFile = ext.parseStringArg(args[i], null);
			    numArgs--;
		    } else if (args[i].startsWith("after=")) {
			    forceAfterFile = ext.parseStringArg(args[i], null);
			    numArgs--;
		    } else if (args[i].startsWith("regardless=")) {
			    forceRegardlessFile = ext.parseStringArg(args[i], null);
			    numArgs--;
		    } else if (args[i].startsWith("-compare")) {
			    compare = true;
			    numArgs--;
		    } else if (args[i].startsWith("snps=")) {
			    numSNPs = ext.parseIntArg(args[i]);
			    numArgs--;
		    } else {
		    	System.err.println("Error - don't know what to do with argument: "+args[i]);
		    }
	    }
	    if (numArgs!=0) {
	    	System.err.println(usage);
		    System.exit(1);
	    }

//	    compare = true;
//	    numSNPs = 768;
	    
	    filename = "indepHapMap.crf";
	    selectFromParameters(filename, new Logger(null));
	    System.exit(1);
	    
	    
	    try {
	    	if (compare) {
	    		compareSelectionParameters(numSNPs, dir, filename, pval_threshold, ldRoots, position_reportType, r2_compType, r2_threshold, dirIlluminaScores, filteringDataset, forceBeforeFile, forceAfterFile, forceRegardlessFile, new Logger(null));
	    	} else {
	    		findOptimalSet(dir, filename, outputRoot, numSNPs, pval_threshold, ldRoots, position_reportType, r2_compType, r2_threshold, dirIlluminaScores, scoreThresholds, scoreDiffThreshold, scoreClassBump, filteringDataset, forceBeforeFile, forceAfterFile, forceRegardlessFile, new Logger(null));
	    	}
	    } catch (Exception e) {
		    e.printStackTrace();
	    }
    }
}
