// -Xms1536M -Xmx1536M
//
// used to filter out markers because of hwe, but since Haploview doesn't limit itself to controls,
// I truned this off (-hwcutoff 0)
// always filter on hwe before using as a reference
//
// should the chromosomal positions be stored in LDdatabase and serialized to be used in subsequent
// calculations?
package org.genvisis.dead;

// import java.io.*;
// import java.util.*;
//
// import common.*;
// import filesys.*;

public class CopyAlgorithm {
	// public static final double DEFAULT_INDEX_THRESHOLD = 0.001;
	// public static final double DEFAULT_INCLUSION_THRESHOLD = 0.01;
	// public static final int DEFAULT_WINDOW = 150000;
	//
	// public static final int NUM_CHROMOSOMES = 27;
	//
	// public static void findOptimalSet(String dir, String filename, String outputRoot, int numSNPs,
	// float pval_threshold, String ldRoot, float r2_threshold, String directoryOfIlluminaScores,
	// float scoreThreshold, float scoreDiffThreshold, float scoreClassBump, String filteringDataset,
	// String forceBeforeFile, String forceAfterFile, String forceRegardlessFile) {
	// BufferedReader reader;
	// PrintWriter writer;
	// Vector<String> tags, untaggedTags, checkTags;
	// long time;
	// ResultSet results;
	// IntVector[] chrIVs;
	// IntVector iv, beforeIndicesVector;
	// String[] markerNames, superset, subset, forcedKeys;
	// float[] pvalues, pvals;
	// LDdatabase lddb;
	// int[] positions, keys;
	// Hashtable<String,String> chrHash, forceBefore, forceAfter, forceRegardless;
	// String trav;
	// int index, offset;
	//// StringLDdb chrLDdb;
	// LongLDdb chrLDdb;
	// float r2, pval;
	// String[] files, line;
	// Hashtable<String,Float> scores;
	// int[] indices;
	// Vector<String> missingMarkers, missingIlluminaValues;
	// Hashtable<String, String> merges, allMissingMarkers;
	// float score, travScore;
	// int bestIndex;
	// double bestDiff, diff;
	// int beforeIndex;
	// int[] beforeIndices;
	//
	// if (!new File(dir+filename+".rset").exists()) {
	// prep(dir+filename);
	// }
	// System.out.print("Reading results file...");
	// time = new Date().getTime();
	// results = ResultSet.load(dir+filename+".rset", false, false);
	// System.out.println("...finished in "+ext.getTimeElapsed(time));
	//
	// scores = new Hashtable<String,Float>();
	// merges = new Hashtable<String,String>();
	// if (directoryOfIlluminaScores != null) {
	// files = Files.list(directoryOfIlluminaScores, ".csv", false);
	// for (int i = 0; i<files.length; i++) {
	// try {
	// reader = new BufferedReader(new FileReader(directoryOfIlluminaScores+files[i]));
	// line = new String[] {""};
	// while (reader.ready() && !line[0].startsWith("Locus_Name")) {
	// line = reader.readLine().trim().split(",");
	// }
	// indices = ext.indexFactors(ILLUMINA_TARGET_COLUMNS, line, false, null, false, false);
	// if (!reader.ready()) {
	// System.err.println("Error - failed to parse '"+files[i]+"'; no line started with
	// 'Locus_Name'");
	// } else if (indices[1] == -1 && indices[2] == -1) {
	// System.err.println("Error - Illumina files must have one of the following columns:
	// "+Array.toStr(Array.subArray(ILLUMINA_TARGET_COLUMNS, 1), ", "));
	// System.err.println(" - only found: "+Array.toStr(line, ", "));
	// } else {
	// while (reader.ready()) {
	// line = reader.readLine().trim().split(",", -1);
	// if (indices[1] >= 0) {
	// if (line[indices[1]].startsWith("merged to ")) {
	// merges.put(line[indices[0]], line[indices[1]].substring(10));
	// } else if (line[indices[1]].endsWith("not supported")) {
	// scores.put(line[indices[0]], Float.parseFloat("-8"));
	// } else {
	// System.err.println("Error - unknown Failure_Message (assuming it's bad): "+line[indices[1]]);
	// scores.put(line[indices[0]], Float.parseFloat("-7"));
	// }
	// } else {
	// if (line[indices[2]].equals("")) {
	// scores.put(line[indices[0]], Float.parseFloat("0.05"));
	// } else {
	// scores.put(line[indices[0]], Float.parseFloat(line[indices[2]]));
	// }
	// }
	// }
	// }
	// reader.close();
	// } catch (FileNotFoundException fnfe) {
	// System.err.println("Error: file \""+directoryOfIlluminaScores+files[i]+"\" not found in current
	// directory");
	// System.exit(1);
	// } catch (IOException ioe) {
	// System.err.println("Error reading file \""+directoryOfIlluminaScores+files[i]+"\"");
	// System.exit(2);
	// }
	// }
	// if (files.length == 0) {
	// System.err.println("Error - no illumina scores files found in score directory; no warnings will
	// be provided");
	// directoryOfIlluminaScores = null;
	// }
	// }
	//
	// lddb = new LDdatabase(ldRoot, LDdatabase.TYPE_LONG);
	//
	// iv = new IntVector();
	//
	// markerNames = results.getMarkerNames();
	// pvalues = results.getPvals();
	//
	// for (int i = 0; i<pvalues.length; i++) {
	// if (pvalues[i] < pval_threshold) {
	// iv.add(i);
	// }
	// }
	//
	// superset = new String[iv.size()];
	// for (int j = 0; j<iv.size(); j++) {
	// superset[j] = markerNames[iv.elementAt(j)];
	// }
	// lddb.updateWithTheseMarkers(superset, ext.replaceDirectoryCharsWithUnderscore(dir+filename,
	// 2));
	//
	// if (forceBeforeFile != null && new File(dir+forceBeforeFile).exists()) {
	// System.out.println("Forcing those SNPs in '"+forceBeforeFile+"' to be included as tags if
	// p-value met");
	// forceBefore = HashVec.loadFileToHashNull(dir+forceBeforeFile, false);
	// } else {
	// System.out.println((forceBeforeFile==null?"No ":"No file named '"+forceBeforeFile+"'; no
	// ")+"SNPs will be forced to be included as tags");
	// forceBefore = new Hashtable<String,String>();
	// }
	// if (forceAfterFile != null && new File(dir+forceAfterFile).exists()) {
	// System.out.println("Forcing those SNPs in '"+forceAfterFile+"' to be included as tags if
	// p-value met");
	// forceAfter = HashVec.loadFileToHashNull(dir+forceAfterFile, false);
	// } else {
	// System.out.println((forceAfterFile==null?"No ":"No file named '"+forceAfterFile+"'; no ")+"SNPs
	// will be forced to be included as tags after tagging");
	// forceAfter = new Hashtable<String,String>();
	// }
	// if (forceRegardlessFile != null && new File(dir+forceRegardlessFile).exists()) {
	// System.out.println("Forcing those SNPs in '"+forceRegardlessFile+"' to be included as tags
	// regardless of p-value");
	// forceRegardless = HashVec.loadFileToHashNull(dir+forceRegardlessFile, false);
	// forcedKeys = HashVec.getKeys(forceRegardless);
	// } else {
	// System.out.println((forceRegardlessFile==null?"No ":"No file named '"+forceRegardlessFile+"';
	// no ")+"additional SNPs will be forced regardless of p-value");
	// forceRegardless = new Hashtable<String,String>();
	// forcedKeys = new String[0];
	// }
	//
	// chrHash = lddb.getChrHash();
	// tags = new Vector<String>();
	// missingIlluminaValues = new Vector<String>();
	// allMissingMarkers = new Hashtable<String,String>();
	// missingMarkers = new Vector<String>();
	// chrIVs = IntVector.newIntVectors(27);
	// for (int j = 0; j<iv.size(); j++) {
	// score = getScore(superset[j], scores, merges, missingIlluminaValues);
	//// if (directoryOfIlluminaScores == null || score > scoreThreshold) {
	// if (directoryOfIlluminaScores == null || score > scoreThreshold || (pvalues[iv.elementAt(j)] <
	// 1E-4 && score > 0.4)) {
	// trav = chrHash.get(superset[j]);
	// if (trav == null) {
	//// System.err.println("Error - the chrHash derived from the marker file does not contain marker
	// "+markerNames[ivs[i].elementAt(j)]);
	// HashVec.addIfAbsent(superset[j], missingMarkers);
	// allMissingMarkers.put(superset[j], "");
	// chrIVs[0].add(iv.elementAt(j));
	// } else {
	// chrIVs[Byte.parseByte(trav.split("[\\s]+")[0])].add(iv.elementAt(j));
	// }
	// }
	// }
	//
	// for (int chr = 0; chr<NUM_CHROMOSOMES; chr++) {
	// subset = new String[chrIVs[chr].size()];
	// pvals = new float[chrIVs[chr].size()];
	// positions = new int[chrIVs[chr].size()];
	// beforeIndicesVector = new IntVector();
	// for (int j = 0; j<positions.length; j++) {
	// subset[j] = markerNames[chrIVs[chr].elementAt(j)];
	// pvals[j] = pvalues[chrIVs[chr].elementAt(j)];
	// trav = chrHash.get(subset[j]);
	// if (trav == null) {
	// positions[j] = j;
	// } else {
	// positions[j] = Integer.parseInt(trav.split("[\\s]+")[1]);
	// }
	// if (forceBefore.containsKey(subset[j])) {
	// forceBefore.put(subset[j], pvals[j]+"");
	// beforeIndicesVector.add(j);
	// }
	// if (forceAfter.containsKey(subset[j])) {
	// forceAfter.put(subset[j], pvals[j]+"");
	// }
	// if (forceRegardless.containsKey(subset[j])) {
	// forceRegardless.put(subset[j], pvals[j]+"");
	// }
	// }
	// keys = Sort.quicksort(positions);
	// subset = Sort.putInOrder(subset, keys);
	// positions = Sort.putInOrder(positions, keys);
	// pvals = Sort.putInOrder(pvals, keys);
	//
	// beforeIndices = beforeIndicesVector.toArray();
	// beforeIndex = 0;
	// chrLDdb = lddb.getLongChrLDdb(chr);
	// while (pvals.length > 0 && Array.min(pvals) < 2) {
	// index = -1;
	// pval = -1;
	//
	// if (beforeIndex < beforeIndices.length) {
	// index = beforeIndices[beforeIndex];
	// pvals[index] = 3;
	// pval = pvals[index];
	// beforeIndex++;
	// } else {
	// index = Array.indexOfMin(pvals);
	// pval = pvals[index];
	// score = getScore(subset[index], scores, merges, missingIlluminaValues);
	// bestDiff = scoreDiffThreshold;
	// bestIndex = -1;
	// for (int j = 0; j<2; j++) {
	// offset = j==0?-1:1;
	// for (int k = 1; index+k*offset >= 0 &&
	// index+k*offset < pvals.length &&
	// Math.abs(positions[index] - positions[index+k*offset]) < LDdatabase.BP_LIMIT; k++) {
	// if (pvals[index+k*offset] < 2) {
	// r2 = chrLDdb.get(subset[index], subset[index+k*offset]);
	// if (r2 == LDdatabase.MISSING_INFO) {
	// if (!allMissingMarkers.containsKey(subset[index]) ||
	// !allMissingMarkers.containsKey(subset[index+k*offset])) {
	// System.err.println("Error - missing LD info for "+subset[index]+"/"+subset[index+k*offset]+"
	// pair");
	// }
	// } else if (r2 > r2_threshold){
	// travScore = getScore(subset[index+k*offset], scores, merges, missingIlluminaValues);
	// diff = (travScore - score) / (-1*Math.log10(pval) - -1*Math.log10(pvals[index+k*offset])) +
	// scoreClassBump*(getScoreClass(travScore)-getScoreClass(score));
	// if (travScore > score && diff > bestDiff) {
	// bestDiff = diff;
	// bestIndex = index+k*offset;
	// }
	// }
	// }
	// }
	// }
	// if (bestIndex >= 0) {
	//// System.out.println("From "+subset[index]+" to "+subset[bestIndex]);
	// index = bestIndex;
	// }
	// pval = pvals[index];
	// pvals[index] = 2;
	// }
	// tags.add(subset[index]+"\t"+chr+"\t"+pval+"\t"+getScore(subset[index], scores, merges,
	// missingIlluminaValues));
	//
	// for (int j = 0; j<2; j++) {
	// offset = j==0?-1:1;
	// for (int k = 1; index+k*offset >= 0 &&
	// index+k*offset < pvals.length &&
	// Math.abs(positions[index] - positions[index+k*offset]) < LDdatabase.BP_LIMIT; k++) {
	// if (pvals[index+k*offset] < 2) {
	// r2 = chrLDdb.get(subset[index], subset[index+k*offset]);
	// if (r2 == LDdatabase.MISSING_INFO) {
	// if (!allMissingMarkers.containsKey(subset[index]) ||
	// !allMissingMarkers.containsKey(subset[index+k*offset])) {
	// System.err.println("Error - missing LD info for "+subset[index]+"/"+subset[index+k*offset]+"
	// pair");
	// }
	// } else if (r2 > r2_threshold){
	// pvals[index+k*offset] = 2;
	// }
	// }
	// }
	// }
	// }
	// for (int j = 0; j<subset.length; j++) {
	// if (forceAfter.containsKey(subset[j]) && Double.parseDouble(forceAfter.get(subset[j])) <
	// pval_threshold) {
	// HashVec.addIfAbsent(subset[j]+"\t"+forceAfter.get(subset[j])+"\t"+getScore(subset[j], scores,
	// merges, missingIlluminaValues), tags);
	// }
	// }
	// }
	//
	// for (int j = 0; j<forcedKeys.length; j++) {
	// HashVec.addIfAbsent(forcedKeys[j]+"\t"+forceRegardless.get(forcedKeys[j])+"\t"+getScore(forcedKeys[j],
	// scores, merges, missingIlluminaValues), tags);
	// }
	// if (missingMarkers.size() > 0) {
	// System.err.println("Error - Missing "+missingMarkers.size()+" markers for threshold
	// "+pval_threshold);
	// Files.writeList(Array.toStringArray(missingMarkers), dir+pval_threshold+"_missingValues.txt");
	// }
	//
	// untaggedTags = new Vector<String>();
	// checkTags = new Vector<String>();
	// if (filteringDataset != null) {
	// chrHash.clear();
	// System.out.print("Loading which markers are present in the "+filteringDataset+" dataset...");
	// chrHash = SnpMarkerSet.loadSnpMarkerSetToChrHash(filteringDataset);
	// System.out.println("done");
	// for (int j = 0; j<tags.size(); j++) {
	// trav = tags.elementAt(j).split("[\\s]+")[0];
	// if (!chrHash.containsKey(trav)) {
	// untaggedTags.add(trav+"\t1");
	// } else {
	// checkTags.add(trav+"\t1");
	// }
	// }
	// }
	//
	// try {
	// writer = new PrintWriter(new FileWriter(dir+"Hits_summary.xln"));
	// writer.println("\tp<"+pval_threshold);
	// writer.println("Total SNPs meeting threshold"+"\t"+iv.size());
	// writer.println("Independent index SNPs meeting threshold"+"\t"+tags.size());
	//
	// if (filteringDataset != null) {
	// writer.println("Independent index SNPs meeting threshold on filtering
	// array"+"\t"+checkTags.size());
	// writer.println("Independent index SNPs meeting threshold not on filtering
	// array"+"\t"+untaggedTags.size());
	// System.out.println(untaggedTags.size()+" SNPs");
	// }
	//
	// writer.close();
	// } catch (Exception e) {
	// System.err.println("Error writing to "+"Hits_summary.xln");
	// e.printStackTrace();
	// }
	//
	// Files.writeList(trimList(tags, numSNPs), dir+(outputRoot ==
	// null?pval_threshold:outputRoot)+"_tags.xln");
	// if (filteringDataset != null) {
	// Files.writeList(Array.toStringArray(checkTags), dir+(outputRoot ==
	// null?pval_threshold:outputRoot)+"_checkArray_tags.xln");
	// Files.writeList(Array.toStringArray(untaggedTags), dir+(outputRoot ==
	// null?pval_threshold:outputRoot)+"_untagged_tags.xln");
	// }
	// }
	//
	// public static String[] trimList(Vector<String> tags, int numSNPs) {
	// String[] finalList;
	// double[] values;
	// int[] order;
	//
	// if (tags.size() < numSNPs) {
	// return Array.toStringArray(tags);
	// }
	//
	// values = new double[tags.size()];
	// for (int i = 0; i<values.length; i++) {
	// values[i] = Double.parseDouble(tags.elementAt(i).split("[\\s]+")[2]);
	// }
	// order = Sort.quicksort(values);
	//
	// finalList = new String[numSNPs];
	// for (int i = 0; i<numSNPs; i++) {
	// finalList[i] = tags.elementAt(order[i]);
	// }
	//
	// return finalList;
	// }
	//
	// public static int getScoreClass(float score) {
	// if (score > 1.0) {
	// return 6;
	// } else if (score > 0.8) {
	// return 5;
	// } else if (score > 0.6) {
	// return 4;
	// } else if (score > 0.4) {
	// return 3;
	// } else if (score > 0.2) {
	// return 2;
	// } else if (score > 0) {
	// return 1;
	// } else {
	// return 0;
	// }
	// }
	//
	// public static float getScore(String element, Hashtable<String,Float> scores,
	// Hashtable<String,String> merges, Vector<String> missingValues) {
	// if (scores.containsKey(element)) {
	// return scores.get(element).floatValue();
	// } else if (merges.containsKey(element)) {
	// if (scores.containsKey(merges.get(element))) {
	// return scores.get(merges.get(element)).floatValue();
	// } else {
	// HashVec.addIfAbsent(merges.get(element), missingValues);
	// System.err.println("Error - "+element+" merged with "+merges.get(element)+"; so we'll need to
	// request info for it");
	// return -9;
	// }
	// } else {
	// HashVec.addIfAbsent(merges.get(element), missingValues);
	// return -9;
	// }
	// }
	//
	// public static void prep(String filename) {
	// System.out.println("Prepping file...");
	// long time = new Date().getTime();
	//// new ResultSet(filename, ResultSet.METAL_TBL_FORMAT).serialize(filename+".rset");
	// new ResultSet(filename, ResultSet.PVALUES_ONLY).serialize(filename+".rset");
	// System.out.println("Finished prepping in "+ext.getTimeElapsed(time));
	// }
	//
	// public static String[] findUnique(String[] inThese, String[] butNotInThese) {
	// Vector<String> v = new Vector<String>();
	//
	// for (int i = 0; i<inThese.length; i++) {
	// if (ext.indexOfStr(inThese[i], butNotInThese) == -1) {
	// v.add(inThese[i]);
	// }
	// }
	//
	// return Array.toStringArray(v);
	// }
	//
	// public static Vector<String> checkCoverage(String hitTags, String arrayTags) {
	// BufferedReader reader;
	// String[] line;
	// Hashtable<String,String> hash;
	// Vector<String> v = new Vector<String>();
	// int numAlleles;
	//
	// hash = HashVec.loadFileToHashString(hitTags, 0, new int[] {0}, "\t", false);
	// try {
	// reader = new BufferedReader(new FileReader(arrayTags));
	// line = reader.readLine().trim().split("[\\s]+");
	// if (!line[0].equals("#captured") || !line[4].equals("alleles")) {
	// System.err.println("Error - Haploview .TAGS file has changed format and needs to be
	// addressed");
	// System.exit(1);
	// }
	// numAlleles = Integer.parseInt(line[3]);
	// if (!reader.readLine().startsWith("#captured") || !reader.readLine().startsWith("#using")) {
	// System.err.println("Error - Haploview .TAGS file has changed format and needs to be
	// addressed");
	// System.exit(1);
	// }
	// for (int i = 0; i<numAlleles; i++) {
	// line = reader.readLine().trim().split("[\\s]+");
	// if (hash.containsKey(line[0]) && line.length == 1) {
	// v.add(line[0]);
	// }
	// }
	// reader.close();
	// } catch (FileNotFoundException fnfe) {
	// System.err.println("Error: file \""+arrayTags+"\" not found in current directory");
	// System.exit(1);
	// } catch (IOException ioe) {
	// System.err.println("Error reading file \""+arrayTags+"\"");
	// System.exit(2);
	// }
	//
	// return v;
	// }
	//
	// public static void compareSelectionParameters(int numSNPs, String dir, String filename, float
	// pval_threshold, String ldRoot, float r2_threshold, String directoryOfIlluminaScores, String
	// filteringDataset, String forceBeforeFile, String forceAfterFile, String forceRegardlessFile) {
	// PrintWriter writer;
	// String[] line;
	//
	// // scoreThreshold, scoreDiffThreshold, scoreClassBump
	// float[][] params = new float[][] {
	//// {Float.NEGATIVE_INFINITY, Float.POSITIVE_INFINITY, 0},
	//// {0, Float.POSITIVE_INFINITY, 0},
	//// {0, 1, 0},
	//// {0, 0.5f, 0},
	//// {0, 0.2f, 0},
	//// {0, 1, 0.1f},
	//// {0, 0.5f, 0.1f},
	//// {0, 0.2f, 0.1f},
	// {0.4f, 0.5f, 0.1f},
	// {0.6f, 0.5f, 0.1f},
	//// {0.8f, 0.5f, 0.1f},
	//// {1, 0.5f, 0.1f},
	//// {0.4f, Float.POSITIVE_INFINITY, 0},
	//// {0.6f, Float.POSITIVE_INFINITY, 0},
	//// {0.8f, Float.POSITIVE_INFINITY, 0},
	//// {1, Float.POSITIVE_INFINITY, 0},
	// };
	//
	// double sumPval;
	// double effSumPval;
	// int[] bins;
	// double failRate;
	// int bin;
	// String[][] data;
	// int[] keys;
	//
	// try {
	// writer = new PrintWriter(new FileWriter(dir+"comparison.xln"));
	// writer.println("scoreThreshold\tscoreDiffThreshold\tscoreClassBump\t\tPredictedFailRate\tSumLogPvals\tEffectiveSum\t%Tagged\t#failDesign\t#0<score<0.2\t#0.2<=score<0.4\t#0.4<=score<0.6\t#0.6<=score<0.8\t#0.8<=score<1.0\t#score=1.1");
	// for (int i = 0; i<params.length; i++) {
	// findOptimalSet(dir, filename, params[i][0]+"_"+params[i][1]+"_"+params[i][2], numSNPs,
	// pval_threshold, ldRoot, r2_threshold, directoryOfIlluminaScores, params[i][0], params[i][1],
	// params[i][2], filteringDataset, forceBeforeFile, forceAfterFile, forceRegardlessFile);
	// data =
	// HashVec.loadFileToStringMatrix(dir+params[i][0]+"_"+params[i][1]+"_"+params[i][2]+"_tags.xln",
	// false, new int[] {0, 1, 2, 3}, false);
	// keys = Sort.quicksort(Array.toDoubleArray(Matrix.extractColumn(data, 2)));
	//
	// bins = new int[7];
	// failRate = -1;
	// sumPval = effSumPval = 0;
	// for (int j = 0; j<keys.length && j<numSNPs; j++) {
	// line = data[keys[j]];
	// bin = getScoreClass(Float.parseFloat(line[3]));
	// bins[bin]++;
	// sumPval += bin==0?0:(-1*Math.log10(Double.parseDouble(line[2])));
	// effSumPval += -1*Math.log10(Double.parseDouble(line[2]))*SCORE_BIN_FAIL_RATES[bin];
	// }
	// failRate = 0;
	// for (int j = 0; j<bins.length; j++) {
	// failRate += SCORE_BIN_FAIL_RATES[j]*bins[j];
	// }
	// failRate /= Array.sum(bins);
	// writer.println(Array.toStr(Array.toDoubleArray(params[i]), 2, 2,
	// "\t")+"\t\t"+failRate+"\t"+sumPval+"\t"+effSumPval+"\t"+"."+"\t"+Array.toStr(bins));
	// }
	// writer.close();
	// } catch (Exception e) {
	// System.err.println("Error writing to "+dir+"comparison.xln");
	// e.printStackTrace();
	// }
	//
	// }
	//
	// public static void selectFromParameters(String filename, Logger log) {
	// BufferedReader reader;
	// String[] line;
	// String trav;
	// String resultsFile = "hits.txt", outputRoot = "tags", ldRoot = DEFAULT_LD_ROOT;
	// float pval_threshold = DEFAULT_PVAL_THRESHOLD, r2_threshold = DEFAULT_R2_THRESHOLD;
	// String directoryOfIlluminaScores = null;
	// float scoreThreshold = 0, scoreDiffThreshold = 0.5f, scoreClassBump = 0.1f;
	// String filteringDataset = null, forceBeforeFile = null, forceAfterFile = null,
	// forceRegardlessFile = null;
	// int numSNPs = Integer.MAX_VALUE;
	//
	// try {
	// reader = new BufferedReader(new FileReader(filename));
	// line = reader.readLine().trim().split("[\\s]+");
	// if (!line[0].equals("indep")) {
	// log.reportError("Error - file must start with the line 'indep'");
	// return;
	// }
	//
	// if (!reader.ready()) {
	// reader.close();
	// generateDefaultParameters(filename);
	// } else {
	// while (reader.ready()) {
	// trav = reader.readLine().trim();
	// if (trav.startsWith("resultsFile=")) {
	// resultsFile = ext.parseStringArg(trav, null);
	// } else if (trav.startsWith("outputRoot=")) {
	// outputRoot = ext.parseStringArg(trav, null);
	// } else if (trav.startsWith("pval_threshold=")) {
	// pval_threshold = ext.parseFloatArg(trav);
	// } else if (trav.startsWith("ldRoot=")) {
	// ldRoot = ext.parseStringArg(trav, null);
	// } else if (trav.startsWith("r2_threshold=")) {
	// r2_threshold = ext.parseFloatArg(trav);
	// } else if (trav.startsWith("directoryOfIlluminaScores=")) {
	// directoryOfIlluminaScores = ext.parseStringArg(trav, null);
	// } else if (trav.startsWith("scoreThreshold=")) {
	// scoreThreshold = ext.parseFloatArg(trav);
	// } else if (trav.startsWith("scoreDiffThreshold=")) {
	// scoreDiffThreshold = ext.parseFloatArg(trav);
	// } else if (trav.startsWith("scoreClassBump=")) {
	// scoreClassBump = ext.parseFloatArg(trav);
	// } else if (trav.startsWith("filteringDataset=")) {
	// filteringDataset = ext.parseStringArg(trav, null);
	// } else if (trav.startsWith("forceBeforeFile=")) {
	// forceBeforeFile = ext.parseStringArg(trav, null);
	// } else if (trav.startsWith("forceAfterFile=")) {
	// forceAfterFile = ext.parseStringArg(trav, null);
	// } else if (trav.startsWith("forceRegardlessFile=")) {
	// forceRegardlessFile = ext.parseStringArg(trav, null);
	// } else if (trav.startsWith("numSNPs=")) {
	// numSNPs = ext.parseIntArg(trav);
	// } else if (!trav.startsWith("#")){
	// System.err.println("Error - don't know what to do with argument: "+trav);
	// }
	// }
	//
	// findOptimalSet("", resultsFile, outputRoot, numSNPs, pval_threshold, ldRoot, r2_threshold,
	// directoryOfIlluminaScores, scoreThreshold, scoreDiffThreshold, scoreClassBump,
	// filteringDataset, forceBeforeFile, forceAfterFile, forceRegardlessFile);
	// }
	// reader.close();
	// } catch (FileNotFoundException fnfe) {
	// log.reportError("Error: file \""+filename+"\" not found in current directory");
	// log.reportException(fnfe);
	// System.exit(1);
	// } catch (IOException ioe) {
	// log.reportError("Error reading file \""+filename+"\"");
	// log.reportException(ioe);
	// System.exit(2);
	// }
	//
	// }
	//
	// public static void generateDefaultParameters(String filename) {
	// PrintWriter writer;
	//
	// try {
	// writer = new PrintWriter(new FileWriter(filename));
	// writer.println("indep");
	// writer.println("resultsFile=hits.txt");
	// writer.println("outputRoot=tags");
	// writer.println("pval_threshold=0.0001");
	// writer.println("ldRoot=/home/npankrat/NCBI/HapMap/CEU_founders/CEU_founders");
	// writer.println("r2_threshold=0.80");
	// writer.println("directoryOfIlluminaScores=null;");
	// writer.println("scoreThreshold=0");
	// writer.println("scoreDiffThreshold=0.50");
	// writer.println("scoreClassBump=0.10");
	// writer.println("filteringDataset=null");
	// writer.println("forceBeforeFile=null");
	// writer.println("forceAfterFile=null");
	// writer.println("forceRegardlessFile=null");
	// writer.println("#numSNPs=100");
	// writer.close();
	// } catch (Exception e) {
	// System.err.println("Error writing to "+filename);
	// e.printStackTrace();
	// }
	// }
	//
	// public static void main(String[] args) {
	// int numArgs = args.length;
	// float r2_threshold = DEFAULT_R2_THRESHOLD;
	// String ldRoot = DEFAULT_LD_ROOT;
	// float pval_threshold = DEFAULT_PVAL_THRESHOLD;
	// String outputRoot = null;
	// String dirIlluminaScores = DEFAULT_DIR_ILLUMINA_SCORES;
	// String filteringDataset = DEFAULT_FORCED_DATASET;
	// String forceBeforeFile = DEFAULT_FORCE_INCLUDE_BEFORE_TAGGING_IF_MEETS_CRITERIA;
	// String forceAfterFile = DEFAULT_FORCE_INCLUDE_AFTER_TAGGING_IF_MEETS_CRITERIA;
	// String forceRegardlessFile = DEFAULT_FORCE_INCLUDE_REGARDLESS;
	// float scoreThreshold = DEFAULT_SCORE_THRESHOLD;
	// float scoreDiffThreshold = DEFAULT_SCORE_DIFF_THRESHOLD;
	// float scoreClassBump = DEFAULT_SCORE_CLASS_BUMP;
	// boolean compare = false;
	// int numSNPs = Integer.MAX_VALUE;
	//
	//// String dir = "C:\\Documents and Settings\\npankrat\\My
	// Documents\\UMN\\Folson\\VTE_meta_analysis\\finalAnalysis\\06 tags\\";
	//// String dir = "C:\\Documents and Settings\\npankrat\\My
	// Documents\\UMN\\Folson\\VTE_meta_analysis\\finalAnalysis\\17 after dealing with CHS\\";
	//// String dir = "C:\\Documents and Settings\\npankrat\\My
	// Documents\\UMN\\Folson\\VTE_meta_analysis\\finalAnalysis\\11 X chromosome\\";
	//// String filename = "META_ANALYSIS_beta_se1.tbl";
	//// String filename = "META_ANALYSIS_beta_se1_parsed_noMissing.xln";
	//// String filename = "META_ANALYSIS_beta_se1_parsed_noMissing_noMAF_LTE_0.278.xln";
	//// String filename = "hits_wMAF_GTE0.024_described.xln";
	//
	// String dir = "C:\\Documents and Settings\\npankrat\\My
	// Documents\\tWork\\Consortium\\analysisOfImputation\\Aff_AAO_combo\\";
	// String filename = "minComputed.xln";
	//
	//// String dir = "C:\\Documents and Settings\\npankrat\\My
	// Documents\\tWork\\Consortium\\analyses\\NGRC\\IlluminaAshk\\";
	//// String filename = "NGRC_Ashk_pvals.xln";
	//
	//// String dir = "C:\\Documents and Settings\\npankrat\\My
	// Documents\\tWork\\Consortium\\Fst\\Discovery\\TopFst\\";
	//// String filename = "TopFst.txt";
	//
	// String usage = "\n"+
	// "gwas.IndependentSNPs requires 0-1 arguments\n"+
	// " (1) directory (i.e. dir="+dir+" (default))\n"+
	// " (2) filename (i.e. file="+filename+" (default))\n"+
	// " (3) root of PLINK file to test for LD (i.e. ldRoot="+ldRoot+" (default))\n"+
	// " (4) r^2 threshold (i.e. r2="+r2_threshold+" (default))\n"+
	// " (5) (optional) directory with Illumina score .csv files (i.e. scores="+dirIlluminaScores+"
	// (default))\n"+
	// " (6) (optional) dataset to filter tags (i.e. filter="+filteringDataset+" (default))\n"+
	// " (7) (optional) file of SNPs to force before tagging [if they meet minimum p-value] (i.e.
	// before="+forceBeforeFile+" (default))\n"+
	// " (8) (optional) file of SNPs to force after tagging [if they meet minimum p-value] (i.e.
	// after="+forceAfterFile+" (default))\n"+
	// " (9) (optional) file of SNPs to force regardless of p-value (i.e.
	// regardless="+forceRegardlessFile+" (default))\n"+
	// " (10) (optional) compare multiple models (i.e. -compare (not the default))\n"+
	// " (11) (optional) number of SNPs (i.e. snps="+numSNPs+" (default))\n"+
	// "";
	//
	// for (int i = 0; i<args.length; i++) {
	// if
	// (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help"))
	// {
	// System.err.println(usage);
	// System.exit(1);
	// } else if (args[i].startsWith("dir=")) {
	// filename = ext.parseStringArg(args[i], "");
	// numArgs--;
	// } else if (args[i].startsWith("file=")) {
	// filename = args[i].split("=")[1];
	// numArgs--;
	// } else if (args[i].startsWith("ldRoot=")) {
	// ldRoot = ext.parseStringArg(args[i], null);
	// numArgs--;
	// } else if (args[i].startsWith("r2=")) {
	// r2_threshold = Float.parseFloat(args[i].split("=")[1]);
	// numArgs--;
	// } else if (args[i].startsWith("scores=")) {
	// dirIlluminaScores = ext.parseStringArg(args[i], null);
	// numArgs--;
	// } else if (args[i].startsWith("filter=")) {
	// filteringDataset = ext.parseStringArg(args[i], null);
	// numArgs--;
	// } else if (args[i].startsWith("before=")) {
	// forceBeforeFile = ext.parseStringArg(args[i], null);
	// numArgs--;
	// } else if (args[i].startsWith("after=")) {
	// forceAfterFile = ext.parseStringArg(args[i], null);
	// numArgs--;
	// } else if (args[i].startsWith("regardless=")) {
	// forceRegardlessFile = ext.parseStringArg(args[i], null);
	// numArgs--;
	// } else if (args[i].startsWith("-compare")) {
	// compare = true;
	// numArgs--;
	// } else if (args[i].startsWith("snps=")) {
	// numSNPs = ext.parseIntArg(args[i]);
	// numArgs--;
	// } else {
	// System.err.println("Error - don't know what to do with argument: "+args[i]);
	// }
	// }
	// if (numArgs!=0) {
	// System.err.println(usage);
	// System.exit(1);
	// }
	//
	// compare = true;
	// numSNPs = 768;
	//
	// try {
	// if (compare) {
	// compareSelectionParameters(numSNPs, dir, filename, pval_threshold, ldRoot, r2_threshold,
	// dirIlluminaScores, filteringDataset, forceBeforeFile, forceAfterFile, forceRegardlessFile);
	// } else {
	// findOptimalSet(dir, filename, outputRoot, numSNPs, pval_threshold, ldRoot, r2_threshold,
	// dirIlluminaScores, scoreThreshold, scoreDiffThreshold, scoreClassBump, filteringDataset,
	// forceBeforeFile, forceAfterFile, forceRegardlessFile);
	// }
	// } catch (Exception e) {
	// e.printStackTrace();
	// }
	// }
}
