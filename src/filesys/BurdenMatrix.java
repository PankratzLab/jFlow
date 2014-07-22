// expand to dynamically load/save a certain chunk of markers at a time  
package filesys;

import java.io.*;
import java.util.*;

import common.*;
import stats.*;

public class BurdenMatrix implements Serializable {
	public static final long serialVersionUID = 1L;

	public static final String[] ALL_KNOWN_ANNOTATIONS = {"ncRNA", "stoploss_SNV", "synonymous_SNV", "intronic", "UTR3", "upstream", "nonsynonymous_SNV", "splicing", "stopgain_SNV", "intergenic", "stoploss_SNV", "downstream", "stopgain_SNV", "synonymous_SNV", "unknown", "nonsynonymous_SNV", "UTR5", "exonic_splicing"};
	public static final String[] DEFAULT_ANNOTATIONS_TO_INCLUDE = {"TRUE", "1", "stopgain_SNV", "stoploss_SNV", "nonsynonymous_SNV", "splicing", "exonic_splicing"};
	public static final String[] ANNOTATION_DELIMITERS = {"/"};

	private SnpMarkerSet markerSet;
	private String[] ids;
	private boolean[] useForImputation;
	private double[][] values;
	private String[] geneNames;
	private int[][] startAndStopPositions;
	private int[] numberOfVariants;
	private int[][] numberOfAlleles;
	private double[] observedFreqs;
	
	public BurdenMatrix(GenotypeMatrix gens, double mafThreshold, String[] annotationsToInclude, String[] allKnownAnnotations, boolean additiveVariants, double mafThresholdToStartImputing, String[] imputeUsingDataFreqFromTheseIDsNotAnnotationFreq, String variantWeightsFile, Logger log) {
		String[] line;
		String[] markerNames;
		byte[][] counts;
		String[][] annotation;
		int geneAnnotation, mafAnnotation, functionAnnotation;
		Hashtable<String, IntVector> geneMappingHash;
		IntVector indicesOfVariantsInGene;
		boolean include;
		double[] alleleFreqs;
		String function;
		CountVector unknownFunctions;
		String gene;
		byte[] chrs;
		int[] positions;
		int index, position, chr, start, stop;
		int[] indices;
		double trav;
		int variantIndex;
		HashSet<String> idsToUse;
		int countObserved;
		String[][] importedWeights;
		Hashtable<String, Double> weightsHash;
		double[] weights;
		int failedForMAF, failedForFunction, included;
		
		ids = gens.getIds();
		markerSet = gens.getMarkerSet();
		counts = gens.getGenotypeCounts();
		annotation = markerSet.getAnnotation();
		markerNames = markerSet.getMarkerNames();

		if (variantWeightsFile != null) {
			importedWeights = HashVec.loadFileToStringMatrix(variantWeightsFile, false, new int [] {0,1}, false);
			weightsHash = new Hashtable<String, Double>();
			for (int i = 0; i < importedWeights.length; i++) {
				try {
					weightsHash.put(importedWeights[i][0], new Double(importedWeights[i][1]));
				} catch (Exception e) {
					if (i>0) {
						log.reportError("Error - '"+importedWeights[i][1]+"' is an invalid weight for marker "+importedWeights[i][0]);
					}
				}
			}

			weights = new double[markerNames.length];
			for (int i = 0; i < markerNames.length; i++) {
				try {
					weights[i] = weightsHash.get(markerNames[i]).doubleValue();
				} catch (Exception e) {
					log.reportError("Error - no valid weight could be found for marker "+markerNames[i]);
				} 
			}
		} else {
			weights = null;
		}
		
		geneAnnotation = 0;
		mafAnnotation = -1;
		functionAnnotation = -1;
		// checks to see if we have enough annotation to build the specified model
		if (annotation == null || annotation[0].length == 0) {
			log.reportError("Error - no annotation available, need at least gene/group name to determine membership");
			return;
		} else if (mafThreshold < 0 && annotationsToInclude == null) {
			log.report("No annotation beyond gene/group name will be included in the burden scores");
		} else if (annotation[0].length >= 3) {
			log.report("Making the assumption that the order of the annotation is gene/group name, then MAF, then functional annotation");
			mafAnnotation = 1;
			functionAnnotation = 2;
			if (annotation[0].length > 3) {
				log.reportError("Warning - found more annotation than is required; assuming everything else is coming after the three expected columns; if this is not correct, there could be problems");
			}
		} else if (annotation[0].length == 2) {
			log.report("Found only one bit of annotation besides gene/group name");
			if (mafThreshold < 0) {
				log.report("Filtering on functional annotation; all frequencies will be included");
				mafAnnotation = -1;
				functionAnnotation = 1;
			} else if (annotationsToInclude == null) {
				log.report("Only filtering based on MAF, not functional annotation");
				mafAnnotation = 1;
				functionAnnotation = -1;
			} else {
				log.reportError("Error - not enough annotation available to implement desired model");
				return;
			}
		} else if (annotation[0].length == 1 && (mafThreshold > 0 || annotationsToInclude != null)) {
			log.reportError("Error - not enough annotation available to implement desired model");
			return;
		}
		
		included = failedForMAF = failedForFunction = 0;
		geneMappingHash = new Hashtable<String, IntVector>();
		unknownFunctions = new CountVector();
		alleleFreqs = new double[markerNames.length];
		for (int i = 0; i < markerNames.length; i++) {
			include = true;
			
			// if filtering on minor allele frequency (maf), then check to see if the maf is less than or equal to the threshold 
			if (mafAnnotation >= 0) {
				if (ext.containsAny(annotation[i][mafAnnotation], ANNOTATION_DELIMITERS)) {
					log.reportError("Error - don't kow what to do with multiple allele frequencies ('"+annotation[i][mafAnnotation]+"' for marker "+markerNames[i]+"); aborting");
					return;
				}
				try {
					if (annotation[i][mafAnnotation].equals("NA")) {
						alleleFreqs[i] = 0;
					} else {
						alleleFreqs[i] = Double.parseDouble(annotation[i][mafAnnotation]);
					}
				} catch (NumberFormatException nfe) {
					log.reportError("Error - invalid minor allele frequency '"+annotation[i][mafAnnotation]+"' for marker "+markerNames[i]+"; aborting");
					return;
				}
				// less than threshold (e.g., maf<0.05), not less than or equal to threshold
				if (alleleFreqs[i] >= mafThreshold && alleleFreqs[i] <= 1-mafThreshold) {

				// less than or equal to threshold (e.g., maf<=0.05), not less than threshold
//				if (alleleFreqs[i] > mafThreshold && alleleFreqs[i] < 1-mafThreshold) {
					include = false;
					failedForMAF++;
				}
			}

			// if filtering on function, then check to see if the annotation is one of the included functions 
			if (functionAnnotation >= 0) {
				function = annotation[i][functionAnnotation];
				for (int j = 0; j < ANNOTATION_DELIMITERS.length; j++) {
					function = ext.replaceAllWith(function, ANNOTATION_DELIMITERS[j], "\t");
				}
				line = function.split("\t");
				for (int j = 0; j < line.length; j++) {
					if (ext.indexOfStr(line[j], annotationsToInclude) == -1) {
						failedForFunction++;
						include = false;
						// if not a known function, then save counts and report later
						if (ext.indexOfStr(line[j], allKnownAnnotations) == -1) {
							unknownFunctions.add(line[j]);
						}					
					}
				}
			}
			
			// if still valid, then map index to the gene
			if (include) {
				gene = annotation[i][geneAnnotation];
				for (int j = 0; j < ANNOTATION_DELIMITERS.length; j++) {
					gene = ext.replaceAllWith(gene, ANNOTATION_DELIMITERS[j], "\t");
				}
				line = gene.split("\t");
				if (line.length > 1) {
					log.report("Multiple mappings for marker "+markerNames[i]+" (n="+line.length+")");
				}
				for (int j = 0; j < line.length; j++) {
					if (geneMappingHash.containsKey(line[j])) {
						indicesOfVariantsInGene = geneMappingHash.get(line[j]);
					} else {
						geneMappingHash.put(line[j], indicesOfVariantsInGene = new IntVector());
					}
					indicesOfVariantsInGene.add(i);
					included++;
				}
			}
		}
		
		geneNames = HashVec.getKeys(geneMappingHash, true, false);
		System.out.println("Included "+included+" variant to gene mappings in "+geneNames.length+" genes ("+markerNames.length+" markers attempted; "+failedForMAF+" failed for maf; "+failedForFunction+" failed for function)");
		chrs = markerSet.getChrs();
		positions = markerSet.getPositions();
		numberOfVariants = new int[geneNames.length];
		numberOfAlleles = new int[geneNames.length][ids.length];
		startAndStopPositions = new int[geneNames.length][];
		for (int i = 0; i < geneNames.length; i++) {
			indicesOfVariantsInGene = geneMappingHash.get(geneNames[i]);
			numberOfVariants[i] = indicesOfVariantsInGene.size();
			chr = Integer.MIN_VALUE;
			start = Integer.MAX_VALUE;
			stop = Integer.MIN_VALUE;
			for (int j = 0; j < indicesOfVariantsInGene.size(); j++) {
				index = indicesOfVariantsInGene.elementAt(j);
				position = positions[index];
				if (chr == Integer.MIN_VALUE) {
					chr = chrs[index];
				} else if (chr != chrs[index] && !geneNames[i].equals("NA")) {
					log.reportError("Error - variant "+markerNames[index]+" is on a different chromosome (chr"+chrs[index]+") than previous variants in gene "+geneNames[i]+" (chr"+chr+")");
				}
				if (position < start) {
					start = position;
				}
				if (position > stop) {
					stop = position;
				}
			}
			startAndStopPositions[i] = new int[] {chr, start, stop};
		}
		
		if (imputeUsingDataFreqFromTheseIDsNotAnnotationFreq != null) {
			if (imputeUsingDataFreqFromTheseIDsNotAnnotationFreq.length == 0) {
				useForImputation = Array.booleanArray(ids.length, true);
			} else {
				useForImputation = new boolean[ids.length];
				idsToUse = HashVec.loadToHashSet(imputeUsingDataFreqFromTheseIDsNotAnnotationFreq);
				for (int i = 0; i < ids.length; i++) {
					useForImputation[i] = idsToUse.contains(ids[i]);
				}
			}
			observedFreqs = Array.doubleArray(markerNames.length, Double.MAX_VALUE);
		} else {
			useForImputation = null;
			observedFreqs = null;
		}
		
		values = new double[geneNames.length][ids.length];
		for (int geneIndex = 0; geneIndex < geneNames.length; geneIndex++) {
			indices = geneMappingHash.get(geneNames[geneIndex]).toArray();
			for (int idIndex = 0; idIndex < ids.length; idIndex++) {
				for (int i = 0; i < indices.length; i++) {
					trav = 0;
					variantIndex = indices[i];
					if (counts[variantIndex][idIndex] > 0) {
						trav = counts[variantIndex][idIndex];
					} else if (counts[variantIndex][idIndex] < 0) {
						if (imputeUsingDataFreqFromTheseIDsNotAnnotationFreq != null) {
							if (alleleFreqs[variantIndex] >= mafThresholdToStartImputing) {
								if (observedFreqs[variantIndex] == Double.MAX_VALUE) {
									observedFreqs[variantIndex] = 0;
									countObserved = 0;
									for (int j = 0; j < ids.length; j++) {
										if (useForImputation[j] && counts[variantIndex][j] >= 0) {
											observedFreqs[variantIndex] += counts[variantIndex][j];
											countObserved++;
										}
									}
									observedFreqs[variantIndex] = observedFreqs[variantIndex] / (double)countObserved;									
								}
								trav = observedFreqs[variantIndex];
							}							
						} else {
							if (alleleFreqs[variantIndex] >= mafThresholdToStartImputing) {
								trav = alleleFreqs[variantIndex];
							}
						}
					}
					if (alleleFreqs[variantIndex] >= 1-mafThreshold) {
						trav = 2-trav;
					}
					if (!additiveVariants) {
						trav = Math.min(trav, 1);
					}
					numberOfAlleles[geneIndex][idIndex] += trav;
					if (weights != null) {
						trav *= weights[variantIndex];
					}
					
					values[geneIndex][idIndex] += trav;
				}
//				if (maxValueOfOne && value > 0) {
//					value = 1;
//				}
			}
		}

	}

	public void writeToFile(String newMatrixFile, String phenotypeFileToMergeWith, String newMapFile, Logger log) {
		PrintWriter writer;
		int count;
		String head, delimiter;
		String[][] data;
		boolean isValid;
		Hashtable<String, String> hash;
		
		if (phenotypeFileToMergeWith != null) {
			head = Files.getFirstNLinesOfFile(phenotypeFileToMergeWith, 1, log)[0];
			delimiter = ext.determineDelimiter(head);
			if (ext.indexOfStr(head.split(delimiter)[0], ext.COMMON_IDS, false, true) == -1) {
				System.err.println("Warning - did not recognize the header; need an id in the first column; looking for complete data after that for inclusion");
				System.err.println("Found :  "+head);
			}
			
			hash = new Hashtable<String, String>();
			data = HashVec.loadFileToStringMatrix(phenotypeFileToMergeWith, false, null, delimiter, false, 1000, false);
			for (int i = 0; i < data.length; i++) {
				isValid = true;				
				for (int j = 1; j < data[i].length; j++) {
					try {
						Double.parseDouble(data[i][j]);
					} catch (Exception e) {
						isValid = false;
					}
				}
				if (isValid) {
					hash.put(data[i][0], "");
				}
			}
		} else {
			hash = null;
		}
		
		
		count = 0;
		try {
			writer = new PrintWriter(new FileWriter(newMatrixFile));
//			writer.print("gene");
//			for (int i = 0; i < ids.length; i++) {
//				if (hash == null || hash.containsKey(ids[i])) {
//					writer.print(","+ids[i]);
//				}
//			}

			for (int i = 0; i < ids.length; i++) {
				if (hash == null || hash.containsKey(ids[i])) {
					writer.print((i==0?"\"":"\",\"")+ids[i]);
				}
			}
			
			writer.println("\"");
			for (int geneIndex = 0; geneIndex < geneNames.length; geneIndex++) {
//				writer.print(geneNames[geneIndex]);
				writer.print("\""+geneNames[geneIndex]+"\"");
				
				for (int idIndex = 0; idIndex < ids.length; idIndex++) {
					if (hash == null || hash.containsKey(ids[idIndex])) {
						writer.print(","+values[geneIndex][idIndex]);
					}
				}
				writer.println();
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + newMatrixFile);
			e.printStackTrace();
		}
		if (count > 0) {
			log.report("There were "+count+" byte outliers in the dataset");
		}

		try {
			writer = new PrintWriter(new FileWriter(newMapFile));
			writer.println("AnalysisUnit,CHROM,POS,POS_STOP,REF,ALT");
			for (int i = 0; i < geneNames.length; i++) {
				writer.println(geneNames[i]+","+Array.toStr(startAndStopPositions[i], ",")+",N,B");
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + newMapFile);
			e.printStackTrace();
		}
	}
	
	public SnpMarkerSet getMarkerSet() {
		return markerSet;
	}

	public String[] getIds() {
		return ids;
	}

	public double[][] getBurdenScores() {
		return values;
	}

	public void analyze(String phenoFile, String phenoMissingValue, String geneListSubset, String outputFile, boolean verbose, Logger log) {
		PrintWriter writer;
		String[] line;
		Hashtable<String, String> hash;
		HashSet<String> genes;
		int countSamplesUsed;
		String[] traits;
		boolean[] use, analyze;
		double[] deps;
		double[][] indeps;
		boolean logistic;
		RegressionModel model;
		double[] betas, stderrs, pvals;
		String[] names, namesUsed;
		double[] snpCounts;
		int indexOfBeta;
		String delimiter;

		traits = Files.getHeaderOfFile(phenoFile, Files.determineDelimiter(phenoFile, log), log);
		names = Array.subArray(traits, 1);
		hash = HashVec.loadFileToHashString(phenoFile, new int[] {0}, Array.subArray(Array.intArray(traits.length), 1, traits.length), phenoFile.endsWith(".csv"), "\t", true, false, false);
		traits = Array.subArray(traits, 1);
		log.report("Missing phenotype is set to '"+phenoMissingValue+"'");

//		markerNames = markerSet.getMarkerNames();
//		chrs = markerSet.getChrs();
//		positions = markerSet.getPositions();
//		alleles = markerSet.getAlleles();
		analyze = Array.booleanArray(geneNames.length, true);
		if (geneListSubset != null) {
			genes = HashVec.loadFileToHashSet(geneListSubset, false);
			for (int i = 0; i < geneNames.length; i++) {
				if (!genes.contains(geneNames[i]) || geneNames[i].equals("NA")) {
					analyze[i] = false;
				}
			}
		}

		for (int i = 0; i < geneNames.length; i++) {
			if (geneNames[i].equals("NA")) {
				analyze[i] = false;
			}
		}

		use = Array.booleanArray(ids.length, true);
		for (int i = 0; i < ids.length; i++) {
			if (hash.containsKey(ids[i])) {
				line = hash.get(ids[i]).split("[\\s]+");
				for (int j = 0; j < line.length; j++) {
					if (!ext.isValidDouble(line[j]) || line[j].equals(phenoMissingValue)) {
						use[i] = false;
					}
				}
			} else {
				use[i] = false;
			}
		}
		
		deps = new double[Array.booleanArraySum(use)];
		indeps = new double[deps.length][traits.length];
		log.report("There are "+deps.length+" samples with complete data", true, verbose);
		countSamplesUsed = 0;
		for (int i = 0; i < ids.length; i++) {
			if (use[i]) {
				line = hash.get(ids[i]).split("[\\s]+");
				deps[countSamplesUsed] = Double.parseDouble(line[0]);
				for (int j = 1; j < traits.length; j++) {
					indeps[countSamplesUsed][j] = Double.parseDouble(line[j]);
				}
				countSamplesUsed++;
			}
		}
		logistic = RegressionModel.isBinaryTrait(Array.toStr(deps).split("[\\s]+"), log);
		log.report("Running a "+(logistic?"logistic":"linear")+" model for trait '"+traits[0]+"'", true, verbose);
		try {
			writer = new PrintWriter(new FileWriter(outputFile));
			delimiter = outputFile.endsWith(".csv")?",":"\t";
			writer.println("gene"+delimiter+"region"+delimiter+"chr"+delimiter+"start"+delimiter+"stop"+delimiter+"NSNPS"+delimiter+"NALLELES"+delimiter+"SUMOFSCORES"+delimiter+"N"+delimiter+"mean"+delimiter+"var"+delimiter+"beta"+delimiter+"se"+delimiter+"pval"); // \tNALLELES
			for (int i = 0; i < geneNames.length; i++) {
				if (analyze[i]) {
					countSamplesUsed = 0;
					for (int j = 0; j < ids.length; j++) {
						if (use[j]) {
							indeps[countSamplesUsed][0] = (double)values[i][j];
							countSamplesUsed++;
						}
					}
					names[0] = "gene";
					model = logistic?new LogisticRegression(deps, indeps, names, false, false):new LeastSquares(deps, indeps, names, false, false);
					betas = model.getBetas();
					stderrs = model.getSEofBs();
					pvals = model.getSigs();
//					stats = model.getStats();
//					int sigfig = 6;
//					System.err.println(betas.length+"\t"+traits.length);
					namesUsed = model.getVarNames();
					if (model.getFinalDependentVariables().length != countSamplesUsed) {
						log.reportError("Mismatched number of missing values for gene/region '"+geneNames[i]+"' ("+model.getFinalDependentVariables().length+", expected "+countSamplesUsed+"); might want to check missing value codes: "+model.getFinalIndependentVariables()[2][0]);
					}
					indexOfBeta = ext.indexOfStr("gene", namesUsed);
					if (indexOfBeta == -1) {
						// skip if no phenotyped indiviudal had the variant
//					writer.print(geneNames[i]+delimiter+"chr"+startAndStopPositions[i][0]+":"+startAndStopPositions[i][1]+"_"+startAndStopPositions[i][2]+delimiter+startAndStopPositions[i][0]+delimiter+startAndStopPositions[i][1]+delimiter+startAndStopPositions[i][2]);
//					writer.println(delimiter+"0"+delimiter+"0"+delimiter+"NA"+delimiter+"NA"+delimiter+"NA"+delimiter+"NA"+delimiter+"NA");
					} else {
						if (geneNames[i].equals("NA")) {
							System.out.println("'"+geneNames[i]+"'\t"+i+"\t"+analyze[i]);
						}
						writer.print(geneNames[i]+delimiter+"chr"+startAndStopPositions[i][0]+":"+startAndStopPositions[i][1]+"_"+startAndStopPositions[i][2]+delimiter+startAndStopPositions[i][0]+delimiter+startAndStopPositions[i][1]+delimiter+startAndStopPositions[i][2]);
						snpCounts = Matrix.extractColumn(indeps, 0);
//						System.out.println(snpCounts.length);
//						writer.println(delimiter+numberOfVariants[i]+delimiter+Array.sum(snpCounts)+delimiter+countSamplesUsed+delimiter+Array.mean(snpCounts)+delimiter+Array.variance(snpCounts)+delimiter+ext.formDeci(betas[indexOfBeta], sigfig, true)+delimiter+ext.formDeci(stderrs[indexOfBeta], sigfig, true)+delimiter+ext.prettyP(pvals[indexOfBeta], sigfig, 4, 3, true));
						writer.println(delimiter+numberOfVariants[i]+delimiter+Array.sum(Array.subArray(numberOfAlleles[i], use))+delimiter+Array.sum(snpCounts)+delimiter+countSamplesUsed+delimiter+Array.mean(snpCounts)+delimiter+Array.variance(snpCounts)+delimiter+betas[indexOfBeta]+delimiter+stderrs[indexOfBeta]+delimiter+pvals[indexOfBeta]); // delimiter+Array.sum(snpCounts)+
					}
					writer.flush();
					
				}
			}
			writer.close();
		} catch (Exception e) {
			System.err.println("Error writing to " + outputFile);
			e.printStackTrace();
		}
	}
	
	public void serialize(String filename) {
		Files.writeSerial(this, filename);
	}

	public static BurdenMatrix load(String filename, boolean jar) {
		return (BurdenMatrix)Files.readSerial(filename, jar, true);
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String genoFile = "SeqGenotypes.dat";
		String annotationFile = null;
		String phenoFile = null;
		String logfile = null;
		GenotypeMatrix gens;
		BurdenMatrix burden;
		Logger log;
		String dir;
		double mafThreshold = 0.05;
		String[] annotationsToInclude = DEFAULT_ANNOTATIONS_TO_INCLUDE;
		String[] allKnownAnnotations = ALL_KNOWN_ANNOTATIONS;
		boolean additiveVariants = false;
//		boolean maxValueOfOne = false;
		double mafThresholdToStartImputing = 999;
		String[] imputeUsingDataFreqFromTheseIDsNotAnnotationFreq = null;
		String variantWeightsFile = null;

		String usage = "\n" + 
		"filesys.GenotypeMatrix requires 0-1 arguments\n" + 
		"   (1) filename (i.e. file=" + genoFile + " (default))\n" + 
		"";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				genoFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			log = new Logger(logfile);
			dir = "D:/LITE/CHARGE-S/ARIC_CHAGE_S_Freeze2/";
			if (!Files.exists(dir)) {
				dir = "C:/LITE/CHARGE-S/ARIC_CHAGE_S_Freeze2/";
			}
			genoFile = dir+"aric_genotypes_frz2_final.csv";
//			genoFile = dir+"slim.csv";
//			annotationFile = dir+"SNPInfo_ExomeFreeze2_120810.csv";
//			annotationFile = dir+"slimMarkers2.csv";
			annotationFile = "SNPInfo_ExomeFreeze2_120810_aafSlim.csv";
//			phenoFile = dir+"phenoCensored2.csv";

//			genoFile = dir+"chr1pter.csv";
//			annotationFile = dir+"chr1pter_SNPInfo.csv";
//			phenoFile = dir+"lnFibrinogen_alanna.csv";

			
			
//			genoFile = dir+"chr1pter.csv";
//			genoFile = dir+"samd.csv";
//			annotationFile = "SNPInfo_ExomeFreeze2_120810_min.csv";
//			phenoFile = dir+"lnFibrinogen_alanna.csv";
//			phenoFile = dir+"lnFibrinogen1b.csv";
//			phenoFile = dir+"lnFibrinogen_vanilla.csv";
			
			phenoFile = dir+"lnFibrinogen.csv";
//			phenoFile = dir+"FVII.csv";
//			phenoFile = dir+"FVIII.csv";
//			phenoFile = dir+"vWF.csv";
			
			
//			maxValueOfOne = false;
			mafThresholdToStartImputing = 0.01;
			imputeUsingDataFreqFromTheseIDsNotAnnotationFreq = RegressionModel.getIDsWithCompleteData(dir+"phenall.csv", false, log);
			

			if (Files.exists(genoFile+".ser")) {
				System.out.println("Loading serialized version: "+genoFile+".ser");
				gens = GenotypeMatrix.load(genoFile+".ser", false);
			} else {
				System.out.println("Loading: "+genoFile);
				gens = new GenotypeMatrix(genoFile, null, annotationFile, log);
				gens.serialize(genoFile+".ser");
			}

			variantWeightsFile = null;
			
			if (Files.exists(genoFile+".T5_burden.ser")) {
				System.out.println("Loading serialized version: "+genoFile+".T5_burden.ser");
				burden = load(genoFile+".T5_burden.ser", false);
			} else {
				System.out.println("Generating: "+genoFile+".T5_burden");
				if (imputeUsingDataFreqFromTheseIDsNotAnnotationFreq != null) {
					System.out.println("Missing values will be replaced with the mean value across the "+imputeUsingDataFreqFromTheseIDsNotAnnotationFreq.length+" samples with complete phenotypic data");
				}
				burden = new BurdenMatrix(gens, mafThreshold, annotationsToInclude, allKnownAnnotations, additiveVariants, mafThresholdToStartImputing, imputeUsingDataFreqFromTheseIDsNotAnnotationFreq, variantWeightsFile, log);
				burden.writeToFile(genoFile+".T5_burden.csv", phenoFile, genoFile+".T5_info.csv", log);
				burden.serialize(genoFile+".T5_burden.ser");
			}
//			gens.writeToPlinkFiles(ext.rootOf(genoFile, false));
			System.out.println("Analyzing: "+phenoFile);
//			burden.analyze(phenoFile, "NA", null, ext.rootOf(phenoFile, false)+".burden.se.metal", true, log);

			phenoFile = dir+"lnFibrinogen.csv";
			burden.analyze(phenoFile, "NA", null, ext.parseDirectoryOfFile(phenoFile)+"ARIC."+ext.rootOf(phenoFile)+".T5.EA."+ext.getDate(new Date(),"")+".csv", true, log);

			phenoFile = dir+"FVII.csv";
			burden.analyze(phenoFile, "NA", null, ext.parseDirectoryOfFile(phenoFile)+"ARIC."+ext.rootOf(phenoFile)+".T5.EA."+ext.getDate(new Date(),"")+".csv", true, log);

			phenoFile = dir+"FVIII.csv";
			burden.analyze(phenoFile, "NA", null, ext.parseDirectoryOfFile(phenoFile)+"ARIC."+ext.rootOf(phenoFile)+".T5.EA."+ext.getDate(new Date(),"")+".csv", true, log);
			
			phenoFile = dir+"vWF.csv";
			burden.analyze(phenoFile, "NA", null, ext.parseDirectoryOfFile(phenoFile)+"ARIC."+ext.rootOf(phenoFile)+".T5.EA."+ext.getDate(new Date(),"")+".csv", true, log);
			
			
			variantWeightsFile = dir+"MadsenBrowningWeights.dat";
			
			if (Files.exists(genoFile+".MadBr_burden.ser")) {
				System.out.println("Loading serialized version: "+genoFile+".MadBr_burden.ser");
				burden = load(genoFile+".MadBr_burden.ser", false);
			} else {
				System.out.println("Generating: "+genoFile+".MadBr_burden");
				if (imputeUsingDataFreqFromTheseIDsNotAnnotationFreq != null) {
					System.out.println("Missing values will be replaced with the mean value across the "+imputeUsingDataFreqFromTheseIDsNotAnnotationFreq.length+" samples with complete phenotypic data");
				}
				burden = new BurdenMatrix(gens, mafThreshold, annotationsToInclude, allKnownAnnotations, additiveVariants, mafThresholdToStartImputing, imputeUsingDataFreqFromTheseIDsNotAnnotationFreq, variantWeightsFile, log);
				burden.writeToFile(genoFile+".MadBr_burden.csv", phenoFile, genoFile+".MadBr_info.csv", log);
				burden.serialize(genoFile+".MadBr_burden.ser");
			}
//			gens.writeToPlinkFiles(ext.rootOf(genoFile, false));
			System.out.println("Analyzing: "+phenoFile);
//			burden.analyze(phenoFile, "NA", null, ext.rootOf(phenoFile, false)+".burden.se.metal", true, log);

			phenoFile = dir+"lnFibrinogen.csv";
			burden.analyze(phenoFile, "NA", null, ext.parseDirectoryOfFile(phenoFile)+"ARIC."+ext.rootOf(phenoFile)+".T5MB.EA."+ext.getDate(new Date(),"")+".csv", true, log);

			phenoFile = dir+"FVII.csv";
			burden.analyze(phenoFile, "NA", null, ext.parseDirectoryOfFile(phenoFile)+"ARIC."+ext.rootOf(phenoFile)+".T5MB.EA."+ext.getDate(new Date(),"")+".csv", true, log);

			phenoFile = dir+"FVIII.csv";
			burden.analyze(phenoFile, "NA", null, ext.parseDirectoryOfFile(phenoFile)+"ARIC."+ext.rootOf(phenoFile)+".T5MB.EA."+ext.getDate(new Date(),"")+".csv", true, log);
			
			phenoFile = dir+"vWF.csv";
			burden.analyze(phenoFile, "NA", null, ext.parseDirectoryOfFile(phenoFile)+"ARIC."+ext.rootOf(phenoFile)+".T5MB.EA."+ext.getDate(new Date(),"")+".csv", true, log);

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

//	main
//	" \n" + 
//	"       Type Extension Description\n" + 
//	"        -1   [any]    auto-detect from extension\n" + 
//	"        0   .csv      CHARGE-S Houston style .csv file\n" + 
//	"";
}
