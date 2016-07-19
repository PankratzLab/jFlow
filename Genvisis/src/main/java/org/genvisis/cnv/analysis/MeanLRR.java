package org.genvisis.cnv.analysis;

import java.io.*;
import java.util.*;

import org.genvisis.cnv.filesys.*;
import org.genvisis.cnv.manage.MarkerDataLoader;
import org.genvisis.cnv.manage.Transforms;
import org.genvisis.common.*;
import org.genvisis.filesys.CNVariant;
import org.genvisis.filesys.Segment;
import org.genvisis.stats.*;

public class MeanLRR {
	public static void createFilesFromFullSample(Project proj, String regionsFile) {
		PrintWriter writer;
		MarkerSet markerSet;
		SampleList sampleList;
		Sample samp;
		String[] samples, markerNames;
		byte[] chrs;
		int[] positions;
		Segment[] regions;
		IntVector[] components;
		int[][] indices;
		float[] lrrs;
		float[][][] data;
		double sum;
		int count;
		long time;
		int[] numberOfMarkers;
		Hashtable<String, String> hash;
		int[][] markerChrIndices;
		boolean[] transChrs;
		Logger log;
		
		log = proj.getLog();
		log.report("Computing LRR values from FullSample. While this is slower, it allows for additional columns containing normalized values.");
		
		markerSet = proj.getMarkerSet();
		markerNames = markerSet.getMarkerNames();
		chrs = markerSet.getChrs();
		positions = markerSet.getPositions();

		hash = proj.getFilteredHash();

		regions = CNVariant.loadUCSCregions(proj.PROJECT_DIRECTORY.getValue()+regionsFile, false);
		components = IntVector.newIntVectors(regions.length);		
		for (int i = 0; i<positions.length; i++) {
			for (int j = 0; j<regions.length; j++) {
				if (chrs[i] == regions[j].getChr() && positions[i] >= regions[j].getStart() && positions[i] <= regions[j].getStop()) {
					if (hash.containsKey(markerNames[i])) {
						log.report(markerNames[i]+" was filtered out");
					} else {
						components[j].add(i);
					}
				}
            }
        }
		indices = IntVector.toIntMatrix(components);
		numberOfMarkers = new int[regions.length];
		for (int i = 0; i<regions.length; i++) {
			numberOfMarkers[i] = indices[i].length;
        }
		
		transChrs = Array.booleanArray(27, false);
		try {
			writer = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue()+"MarkersIn_"+ext.rootOf(regionsFile)+".xln"));
			for (int i = 0; i<regions.length; i++) {
				writer.print(regions[i].getUCSClocation()+"\t"+numberOfMarkers[i]+"\t"+regions[i].getChr()+"\t"+regions[i].getStart()+"\t"+regions[i].getStop());
				transChrs[regions[i].getChr()] = true;
				for (int j = 0; j<indices[i].length; j++) {
					writer.print("\t"+markerNames[indices[i][j]]);
                }
				writer.println();
	        }
            writer.close();
        } catch (Exception e) {
	        log.reportError("Error writing the list of marker names within the regions");
	        e.printStackTrace();
        }
		
		markerChrIndices = markerSet.getIndicesByChr();
		sampleList = proj.getSampleList(); 
		samples = sampleList.getSamples();
		data = new float[regions.length][samples.length][Transforms.TRANSFORMATION_TYPES.length];
		log.report("Computing mean Log R ratios for:");
        time = new Date().getTime();
		for (int i = 0; i<samples.length; i++) {
			if (i%100 == 0) {
				log.report((i+1)+" of "+samples.length+" ("+ext.getTimeElapsed(time)+")");
		        time = new Date().getTime();
			}
			samp = proj.getPartialSampleFromRandomAccessFile(samples[i]);
			for (int trans = 0; trans < Transforms.TRANSFORMATION_TYPES.length; trans++) {
				lrrs = samp.getLRRs();
				if (trans > 0) {
					lrrs = Transforms.transform(lrrs, trans, markerChrIndices, transChrs);
				}				
				for (int j = 0; j<regions.length; j++) {
					sum = 0;
					count = 0;
					for (int k = 0; k<indices[j].length; k++) {
						if (!Double.isNaN(lrrs[indices[j][k]])) {
							sum += lrrs[indices[j][k]];
							count++;
						}
	                }
					data[j][i][trans] = (float)(sum/(double)count);
	            }			
			}
        }
		
		new MeanLRRset(sampleList.getFingerprint(), regions, numberOfMarkers, data, Transforms.TRANFORMATIONS).serialize(proj.PROJECT_DIRECTORY.getValue()+ext.rootOf(regionsFile)+".mlrr");
	}
	
	public static void createFilesFromMarkerData(Project proj, String regionsFile) {
		PrintWriter writer;
		MarkerSet markerSet;
		SampleList sampleList;
		String[] samples, markerNames;
		byte[] chrs;
		int[] positions;
		Segment[] regions;
		Hashtable<String,Vector<String>> components;
		float[] lrrs;
		float[][][] data;
		int[] counts;
		int[] numberOfMarkers;
		Hashtable<String, String> hash;
		MarkerDataLoader markerDataLoader;
		MarkerData markerData;
		Logger log;
		
		log = proj.getLog();		
		log.report("Computing LRR values from MarkerData. This is faster, but does not allow for additional columns containing normalized values.");
		
		markerSet = proj.getMarkerSet();
		markerNames = markerSet.getMarkerNames();
		chrs = markerSet.getChrs();
		positions = markerSet.getPositions();

		hash = proj.getFilteredHash();

		regions = CNVariant.loadUCSCregions(proj.PROJECT_DIRECTORY.getValue()+regionsFile, false);
		if (regions == null) {
			return;
		}
		components = new Hashtable<String, Vector<String>>();		
		for (int i = 0; i<positions.length; i++) {
			for (int j = 0; j<regions.length; j++) {
				if (chrs[i] == regions[j].getChr() && positions[i] >= regions[j].getStart() && positions[i] <= regions[j].getStop()) {
					if (hash.containsKey(markerNames[i])) {
						log.report(markerNames[i]+" was filtered out");
					} else {
						HashVec.addToHashVec(components, j+"", markerNames[i], false);
					}
				}
            }
        }
		
		sampleList = proj.getSampleList(); 
		samples = sampleList.getSamples();
		numberOfMarkers = new int[regions.length];
		data = new float[regions.length][samples.length][1]; // only mean will be computed; no normalization
		try {
			writer = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue()+"MarkersIn_"+ext.rootOf(regionsFile)+".xln"));
			for (int i = 0; i<regions.length; i++) {
				markerNames = Array.toStringArray(components.get(i+""));
				numberOfMarkers[i] = markerNames.length;
				writer.println(regions[i].getUCSClocation()+"\t"+markerNames.length+"\t"+regions[i].getChr()+"\t"+regions[i].getStart()+"\t"+regions[i].getStop()+"\t"+Array.toStr(markerNames));
				log.report("Computing mean Log R ratios for: "+regions[i].getUCSClocation());
				markerDataLoader = MarkerDataLoader.loadMarkerDataFromListInSeparateThread(proj, markerNames);
				counts = new int[samples.length];
				for (int j = 0; j < markerNames.length; j++) {
					markerData = markerDataLoader.requestMarkerData(j);
					if (markerData.getFingerprint() != sampleList.getFingerprint()) {
						log.reportError("Error - mismatched fingerprint for "+markerData.getMarkerName());
					}
					
					lrrs = markerData.getLRRs();
					for (int k = 0; k < samples.length; k++) {
						if (!ext.isMissingValue(lrrs[k]+"")) {
							data[i][k][0] += lrrs[k];
							counts[k]++;
						}
					}
				}
				for (int j = 0; j < samples.length; j++) {
					data[i][j][0] /= (float)counts[j];
				}
			}
            writer.close();
        } catch (Exception e) {
	        log.reportError("Error writing the list of marker names within the regions");
	        e.printStackTrace();
        }
		
		new MeanLRRset(sampleList.getFingerprint(), regions, numberOfMarkers, data, new String[] {Transforms.TRANFORMATIONS[0]}).serialize(proj.PROJECT_DIRECTORY.getValue()+ext.rootOf(regionsFile)+".mlrr");
	}

	public static void analyze(Project proj, String phenotype, String mlrrSetFile) {
        PrintWriter writer;
		MeanLRRset mlrrSet;
		float[][][] data;
		SampleList sampleList;
		String[] samples;
		Segment[] regions;
		int[] numberOfMarkers;
		RegressionModel model;
		Hashtable<String, String> hash;
		double[] pheno;
		String phen;
		Logger log;
		
		log = proj.getLog();		
//		hash = HashVec.loadFileToHashString(proj.getFilename(proj.SAMPLE_DATA_FILENAME), "DNA", new String[] {phenotype}, "");
		hash = HashVec.loadFileToHashString(proj.SAMPLE_DATA_FILENAME.getValue(), "DNA", new String[] {phenotype}, "");
		sampleList = proj.getSampleList(); 
		samples = sampleList.getSamples();
		pheno = new double[samples.length];
		for (int i = 0; i<samples.length; i++) {
			if (hash.containsKey(samples[i])) {
				phen = hash.get(samples[i]);
				if (ext.isMissingValue(phen)) {
					pheno[i] = Double.NaN;
				} else {
					pheno[i] = Double.parseDouble(phen);
				}
			} else {
				log.reportError("Error - no phenotypic data for sample: "+samples[i]);
				pheno[i] = Double.NaN;
			}
        }	
		
		mlrrSet = MeanLRRset.load(mlrrSetFile, false);
		data = mlrrSet.getData();
		regions = mlrrSet.getRegions();
		numberOfMarkers = mlrrSet.getNumerOfMarkersPerRegion();
		if (mlrrSet.getSampleFingerprint() != sampleList.getFingerprint()) {
			log.reportError("Error - the SampleList fingerprint for the MeanLRRset ("+mlrrSet.getSampleFingerprint()+") does not match the Project's SampleList fingerprint ("+sampleList.getFingerprint()+")");
			return;
		}
		
		try {
			writer = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue()+ext.rootOf(mlrrSetFile)+".xln"));
			writer.println("Region\tNumberOfMarkersInRegion\tChr\tStart\tStop\tBeta\tOR\tStat\tp-value");
			for (int i = 0; i<regions.length; i++) {
				writer.print(regions[i].getUCSClocation()+"\t"+numberOfMarkers[i]+"\t"+regions[i].getChr()+"\t"+regions[i].getStart()+"\t"+regions[i].getStop());
				if (numberOfMarkers[i] > 0) {
					model = RegressionModel.determineAppropriate(pheno, FloatVector.toDoubleArray(Matrix.extractColumn(data[i], 0)), false, false);
					if (model.analysisFailed()) {
						writer.println("\t.\t.\t.\t.");
					} else {
						writer.println("\t"+ext.formDeci(model.getBetas()[1], 5)+"\t"+(model.isLogistic()?ext.formDeci(Math.exp(model.getBetas()[1]), 3):".")+"\t"+ext.formDeci(model.getStats()[1], 3)+"\t"+ext.prettyP(model.getSigs()[1]));
					}
				} else {
					writer.println("\t.\t.\t.\t.");
				}
				writer.flush();
	        }
            writer.close();
        } catch (Exception e) {
	        log.reportError("Error writing to "+ext.rootOf(mlrrSetFile)+".xln");
	        e.printStackTrace();
        }
	}

	public static void dump(Project proj, String[] phenotypes, String mlrrSetFile, String regionToDumpOrNullForAll, int transformationToUse) {
        PrintWriter writer;
		MeanLRRset mlrrSet;
		float[][][] data;
		SampleList sampleList;
		String[] samples;
		Segment[] regions;
		Hashtable<String, String> hash;
		int index;
		String[] transformations;
		Logger log;
		
		log = proj.getLog();		
		if (!Files.exists(mlrrSetFile)) {
			log.report("Error - mlrr dataset file '"+mlrrSetFile+"' was never created");
			return;
		}
		
//		hash = HashVec.loadFileToHashString(proj.getFilename(proj.SAMPLE_DATA_FILENAME), "DNA", phenotypes, "\t");
		hash = HashVec.loadFileToHashString(proj.SAMPLE_DATA_FILENAME.getValue(), "DNA", phenotypes, "\t");
		sampleList = proj.getSampleList(); 
		samples = sampleList.getSamples();
		
		mlrrSet = MeanLRRset.load(mlrrSetFile, false);
		data = mlrrSet.getData();
		regions = mlrrSet.getRegions();
		if (mlrrSet.getSampleFingerprint() != sampleList.getFingerprint()) {
			log.reportError("Error - the SampleList fingerprint for the MeanLRRset ("+mlrrSet.getSampleFingerprint()+") does not match the Project's SampleList fingerprint ("+sampleList.getFingerprint()+")");
			return;
		}
		
		if (regionToDumpOrNullForAll != null) {
			regionToDumpOrNullForAll = ext.replaceAllWith(regionToDumpOrNullForAll, new String[][] {{",", ""}});
		}

		try {
			index = -1;
			if (regionToDumpOrNullForAll == null) {
				writer = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue()+ext.rootOf(mlrrSetFile)+"_dump.xln"));
			} else {
				writer = new PrintWriter(new FileWriter(proj.PROJECT_DIRECTORY.getValue()+ext.replaceAllWith(regionToDumpOrNullForAll, ":", "_") +".xln"));
				
				for (int i = 0; i<regions.length; i++) {
					if (regionToDumpOrNullForAll.equals(regions[i].getUCSClocation())) {
						index = i;
					}
		        }
				if (index == -1) {
					log.reportError("Error - Region flagged to be dumped ("+regionToDumpOrNullForAll+") was not found in the MeanLRRset's regions list");
					return;
				}

			}
			writer.print("Sample");
			for (int i = 0; i < phenotypes.length; i++) {
				writer.print("\t"+phenotypes[i].substring(phenotypes[i].lastIndexOf("=")+1));
			}
			transformations = mlrrSet.getTransformations();
			for (int i = 0; i < regions.length; i++) {
				if (transformationToUse == -1) {
					for (int j = 0; j < transformations.length; j++) {
						writer.print("\t"+regions[i].getUCSClocation()+"_"+ext.replaceWithLinuxSafeCharacters(transformations[j], false));
					}
				} else {
					writer.print("\t"+regions[i].getUCSClocation()+"_"+ext.replaceWithLinuxSafeCharacters(transformations[transformationToUse], false));
				}
			}
			writer.println();
			for (int i = 0; i<samples.length; i++) {
				writer.print(samples[i]+"\t"+(hash.containsKey(samples[i])?hash.get(samples[i]):Array.stringArray(phenotypes.length, ".")));
				if (regionToDumpOrNullForAll == null) {
					for (int j = 0; j < regions.length; j++) {
						if (transformationToUse == -1) {
							writer.print("\t"+Array.toStr(data[j][i]));
						} else {
							writer.print("\t"+data[j][i][transformationToUse]);
						}
					}
				} else {
					writer.print("\t"+Array.toStr(data[index][i]));
				}
				writer.println();
	        }
            writer.close();
        } catch (Exception e) {
	        log.reportError("Error writing to "+ext.rootOf(mlrrSetFile)+".xln");
	        e.printStackTrace();
        }		
	}
	
	public static void fromParameters(String filename, Logger log) {
		Vector<String> params;

		params = Files.parseControlFile(filename, "MeanLRR", new String[] {
				"proj=/home/npankrat/projects/GEDI.properties",
				"# list of regions for which to compute the mean LRR: use the format chr8:25129632-25130278 using one region per line",
				"regions=listOfCNPs.txt",
				"# phenotype in SampleData.txt to use; delimit with a comma if several are desired in an export file, only first will be analyzed",
				"pheno=CLASS=Use32Phenotype,CLASS=Sex,Filter=Age",
				"# compute transforms as well (takes much much longer, as in ~20 minutes versus 1 second, since it loads all project data into)",
				"transform=false",
				"analyze=false",
				"dumpAll=true",
				"# If dumping all and only one transform is desired, then delineate the index of the transformation to export, -1 will export all transforms",
				"transIndex=-1",
				"# Alternatively, if all transformations for only a particular region are desired, then uncomment and list here (must be in the regions file)",
				"#dump=chr8:25129632-25130278"
		}, log);

		if (params != null) {
			params.add("log=" + log.getFilename());
			main(Array.toStringArray(params));
		}
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		String regions = "cnps.txt";
		boolean dumpAll = false;
		String dump = null;
		String[] phenotypes = new String[] {"CLASS=Used"};
		int transIndex = -1;
		boolean transform = false;
		boolean analyze = false;
		long time;
		String logfile = null;
		Logger log;
		
//		filename = "/home/npankrat/projects/GEDI_exome.properties";
//		phenotype = "Class=SexReversal";
//		regions = "PAR2.txt";
//		dump = "chr23:154879620-155227607"; // Marker 6

//		filename = "/home/npankrat/projects/GEDI.properties";
//		phenotypes = new String[] {"Class=UsedInCNVAnlayses", "Final_LRR_SD", "Final_conf15usedCount"};
//		regions = "firstCNP.txt";
//		dumpAll = true;
		
		// warning cannot currently accommodate commas
//		dump = "chr1:104187899-104303967"; // AMY1A (processes starch)
//		AMY1A	chr1:104095896-104122149
//		AMY2A	chr1:104159999-104168400
//		AMY2B	chr1:104095896-104122149

//		filename = "/home/npankrat/projects/GEDI.properties";
//		phenotypes = new String[] {"CLASS=Suitable for CNV", "Final_LRR_SD", "Final_conf15usedCount"};

		
//		filename = "/home/npankrat/projects/SingaporeReplication.properties";
//		phenotypes = new String[] {"Class=LQ LRR_SD;2=Bad;3=KindaBad", "Class=BAF_outliers", "Class=Estimated Sex;1=Male;2=Female;3=Klinefelter;4=Mosaic Klinefelter;5=Triple X;6=Turner;7=Mosaic Turner", "Class=Exclude"};
//		dumpAll = true;
//		transform = true;
//		
//		phenotypes = new String[] {"Class=Exclude"};
//		dump = "chr20:58,433,777-58,497,422";
//		dumpAll = true;

		Project proj;

		String usage = "\n"+
				"cnv.analysis.MeanLRR requires 0-1 arguments\n"+"" +
				"   (1) project properties filename (i.e. proj="+org.genvisis.cnv.Launch.getDefaultDebugProjectFile(false)+" (default))\n"+
				"   (2) filename of the regions in UCSC format (chr8:25129632-25130278) (i.e. regions="+regions+" (default))\n"+
				"   (3) phenotype in SampleData.txt; delimit with a comma for export, only first will be analyzed (i.e. pheno="+Array.toStr(phenotypes, ",")+" (default))\n"+
				"   (4) compute transforms as well (takes much much longer) (i.e. transform=false (default))\n"+
				"   (5) run a regression model using the first phenotype (i.e. analyze=false (default))\n"+
				"  ADD the following if you want to dump the data to a text file\n"+
				"   (6) dump all regions for a particular transformation (i.e. dumpAll="+dumpAll+" (default))\n"+
				"   (7) the index of the transformation to export or -1 for all (i.e. transIndex="+transIndex+" (default))\n"+
				"  OR:\n"+
				"   (6) dump all transformations for a particular region (i.e. dump=chr8:25129632-25130278 (not the default))\n"+
				"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				return;
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("regions=")) {
				regions = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("pheno=")) {
				phenotypes = ext.parseStringArg(args[i], null).split(",");
				numArgs--;
			} else if (args[i].startsWith("dump=")) {
				dump = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("dumpAll=")) {
				dumpAll = ext.parseBooleanArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("transform=")) {
				transform = ext.parseBooleanArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("transIndex=")) {
				transIndex = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("analyze=")) {
				analyze = ext.parseBooleanArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = ext.parseStringArg(args[i], null);
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			return;
		}
		try {
			time = new Date().getTime();
			proj = new Project(filename, logfile, false);
			log = proj.getLog();
			if (!new File(proj.PROJECT_DIRECTORY.getValue()+ext.rootOf(regions)+".mlrr").exists()) {
				if (transform) {
					createFilesFromFullSample(proj, regions);
				} else {
					createFilesFromMarkerData(proj, regions);
				}
			}
			if (dumpAll) {
				log.report("Dumping all regions");
				dump(proj, phenotypes, proj.PROJECT_DIRECTORY.getValue()+ext.rootOf(regions)+".mlrr", null, transIndex);
			} else if (dump != null) {
				log.report("Dumping "+dump);
				dump(proj, phenotypes, proj.PROJECT_DIRECTORY.getValue()+ext.rootOf(regions)+".mlrr", dump, -1);
			}
			if (analyze) {
				log.report("Analyzing "+regions+" using "+phenotypes[0]);
				analyze(proj, phenotypes[0], proj.PROJECT_DIRECTORY.getValue()+ext.rootOf(regions)+".mlrr");
			}
			log.report("Finished in " + ext.getTimeElapsed(time));
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
