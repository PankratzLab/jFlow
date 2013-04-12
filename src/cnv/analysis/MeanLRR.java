// -Xms1024M -Xmx1024M is probably not needed
package cnv.analysis;

import java.io.*;
import java.util.*;

import cnv.filesys.*;
import cnv.manage.Transforms;
import cnv.var.CNVariant;
import common.*;
import filesys.Segment;
import stats.*;

public class MeanLRR {
	public static void createFilesFromFullSample(Project proj, String cnpFile) {
		PrintWriter writer;
		MarkerSet markerSet;
		SampleList sampleList;
		Sample samp;
		String[] samples, markerNames;
		byte[] chrs;
		int[] positions;
		Segment[] cnps;
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
		
		System.out.println("Computing LRR values from FullSample. While this is slower, it allows for additional columns containing normalized values.");
		
		markerSet = proj.getMarkerSet();
		markerNames = markerSet.getMarkerNames();
		chrs = markerSet.getChrs();
		positions = markerSet.getPositions();

		hash = proj.getFilteredHash();

		cnps = CNVariant.loadUCSCregions(proj.getProjectDir()+cnpFile, false);
		components = IntVector.newIntVectors(cnps.length);		
		for (int i = 0; i<positions.length; i++) {
			for (int j = 0; j<cnps.length; j++) {
				if (chrs[i] == cnps[j].getChr() && positions[i] >= cnps[j].getStart() && positions[i] <= cnps[j].getStop()) {
					if (hash.containsKey(markerNames[i])) {
						System.out.println(markerNames[i]+" was filtered out");
					} else {
						components[j].add(i);
					}
				}
            }
        }
		indices = IntVector.toIntMatrix(components);
		numberOfMarkers = new int[cnps.length];
		for (int i = 0; i<cnps.length; i++) {
			numberOfMarkers[i] = indices[i].length;
        }
		
		transChrs = Array.booleanArray(27, false);
		try {
			writer = new PrintWriter(new FileWriter(proj.getProjectDir()+"MarkersIn_"+ext.rootOf(cnpFile)+".xln"));
			for (int i = 0; i<cnps.length; i++) {
				writer.print(cnps[i].getUCSClocation()+"\t"+numberOfMarkers[i]+"\t"+cnps[i].getChr()+"\t"+cnps[i].getStart()+"\t"+cnps[i].getStop());
				transChrs[cnps[i].getChr()] = true;
				for (int j = 0; j<indices[i].length; j++) {
					writer.print("\t"+markerNames[indices[i][j]]);
                }
				writer.println();
	        }
            writer.close();
        } catch (Exception e) {
	        System.err.println("Error writing marker Names in CNPs");
	        e.printStackTrace();
        }
		
		markerChrIndices = markerSet.getIndicesByChr();
		sampleList = proj.getSampleList(); 
		samples = sampleList.getSamples();
		data = new float[cnps.length][samples.length][Transforms.TRANSFORMATION_TYPES.length];
		System.out.println("Computing mean Log R ratios for:");
        time = new Date().getTime();
		for (int i = 0; i<samples.length; i++) {
			if (i%100 == 0) {
				System.out.println((i+1)+" of "+samples.length+" ("+ext.getTimeElapsed(time)+")");
		        time = new Date().getTime();
			}
			samp = proj.getSample(samples[i]);
			for (int trans = 0; trans < Transforms.TRANSFORMATION_TYPES.length; trans++) {
				lrrs = samp.getLRRs();
				if (trans > 0) {
					lrrs = Transforms.transform(lrrs, trans, markerChrIndices, transChrs);
				}				
				for (int j = 0; j<cnps.length; j++) {
					sum = 0;
					count = 0;
					for (int k = 0; k<indices[j].length; k++) {
						if (!Double.isNaN(lrrs[indices[j][k]])) {
							sum += lrrs[indices[j][k]];
							count++;
						}
	                }
					data[j][i][trans] = (float)(sum/count);
	            }			
			}
        }
		
		new MeanLRRset(sampleList.getFingerprint(), cnps, numberOfMarkers, data).serialize(proj.getProjectDir()+ext.rootOf(cnpFile)+".mlrr");
	}
	
	public static void createFilesFromMarkerData(Project proj, String cnpFile) {
		PrintWriter writer;
		MarkerSet markerSet;
		SampleList sampleList;
//		Sample samp;
		String[] samples, markerNames;
		byte[] chrs;
		int[] positions;
		Segment[] cnps;
		IntVector[] components;
		int[][] indices;
//		float[] lrrs;
		float[][][] data;
//		double sum;
//		int count;
//		long time;
		int[] numberOfMarkers;
		Hashtable<String, String> hash;
//		int[][] markerChrIndices;
		boolean[] transChrs;
		
		System.out.println("Computing LRR values from MarkerData. This is faster, but does not allow for additional columns containing normalized values.");
		
		markerSet = proj.getMarkerSet();
		markerNames = markerSet.getMarkerNames();
		chrs = markerSet.getChrs();
		positions = markerSet.getPositions();

		hash = proj.getFilteredHash();

		cnps = CNVariant.loadUCSCregions(proj.getProjectDir()+cnpFile, false);
		components = IntVector.newIntVectors(cnps.length);		
		for (int i = 0; i<positions.length; i++) {
			for (int j = 0; j<cnps.length; j++) {
				if (chrs[i] == cnps[j].getChr() && positions[i] >= cnps[j].getStart() && positions[i] <= cnps[j].getStop()) {
					if (hash.containsKey(markerNames[i])) {
						System.out.println(markerNames[i]+" was filtered out");
					} else {
						components[j].add(i);
					}
				}
            }
        }
		indices = IntVector.toIntMatrix(components);
		numberOfMarkers = new int[cnps.length];
		for (int i = 0; i<cnps.length; i++) {
			numberOfMarkers[i] = indices[i].length;
        }
		
		transChrs = Array.booleanArray(27, false);
		try {
			writer = new PrintWriter(new FileWriter(proj.getProjectDir()+"MarkersIn_"+ext.rootOf(cnpFile)+".xln"));
			for (int i = 0; i<cnps.length; i++) {
				writer.print(cnps[i].getUCSClocation()+"\t"+numberOfMarkers[i]+"\t"+cnps[i].getChr()+"\t"+cnps[i].getStart()+"\t"+cnps[i].getStop());
				transChrs[cnps[i].getChr()] = true;
				for (int j = 0; j<indices[i].length; j++) {
					writer.print("\t"+markerNames[indices[i][j]]);
                }
				writer.println();
	        }
            writer.close();
        } catch (Exception e) {
	        System.err.println("Error writing marker Names in CNPs");
	        e.printStackTrace();
        }
		
//		markerChrIndices = markerSet.getIndicesByChr();
		sampleList = proj.getSampleList(); 
		samples = sampleList.getSamples();
		data = new float[cnps.length][samples.length][1]; // only mean will be computed; no normalization
		System.out.println("Computing mean Log R ratios for:");
//        time = new Date().getTime();
//		for (int i = 0; i<samples.length; i++) {	// to be replaced at a later date with marker data
//			if (i%100 == 0) {
//				System.out.println((i+1)+" of "+samples.length+" ("+ext.getTimeElapsed(time)+")");
//		        time = new Date().getTime();
//			}
//			samp = proj.getSample(samples[i]);
//			for (int trans = 0; trans < Transforms.TRANSFORMATION_TYPES.length; trans++) {
//				lrrs = samp.getLRRs();
//				if (trans > 0) {
//					lrrs = Transforms.transform(lrrs, trans, markerChrIndices, transChrs);
//				}				
//				for (int j = 0; j<cnps.length; j++) {
//					sum = 0;
//					count = 0;
//					for (int k = 0; k<indices[j].length; k++) {
//						if (!Double.isNaN(lrrs[indices[j][k]])) {
//							sum += lrrs[indices[j][k]];
//							count++;
//						}
//	                }
//					data[j][i][trans] = (float)(sum/count);
//	            }			
//			}
//        }
		
		new MeanLRRset(sampleList.getFingerprint(), cnps, numberOfMarkers, data).serialize(proj.getProjectDir()+ext.rootOf(cnpFile)+".mlrr");
	}

	public static void analyze(Project proj, String phenotype, String mlrrSetFile) {
        PrintWriter writer;
		MeanLRRset mlrrSet;
		float[][][] data;
		SampleList sampleList;
		String[] samples;
		Segment[] cnps;
		int[] numberOfMarkers;
		RegressionModel model;
		Hashtable<String, String> hash;
		double[] pheno;
		String phen;
		
		hash = HashVec.loadFileToHashString(proj.getFilename(Project.SAMPLE_DATA_FILENAME), "DNA", new String[] {phenotype}, "");
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
				System.err.println("Error - no phenotypic data for sample: "+samples[i]);
				pheno[i] = Double.NaN;
			}
        }	
		
		mlrrSet = MeanLRRset.load(mlrrSetFile, false);
		data = mlrrSet.getData();
		cnps = mlrrSet.getCnps();
		numberOfMarkers = mlrrSet.getNumerOfMarkersPerCNV();
		if (mlrrSet.getSampleFingerprint() != sampleList.getFingerprint()) {
			System.err.println("Error - the SampleList fingerprint for the MeanLRRset ("+mlrrSet.getSampleFingerprint()+") does not match the Project's SampleList fingerprint ("+sampleList.getFingerprint()+")");
			System.exit(1);
		}
		
		try {
			writer = new PrintWriter(new FileWriter(proj.getProjectDir()+ext.rootOf(mlrrSetFile)+".xln"));
			writer.println("CNV\tNumberOfMarkersInCNV\tChr\tStart\tStop\tBeta\tOR\tStat\tp-value");
			for (int i = 0; i<cnps.length; i++) {
//				System.out.println((i+1)+" of "+cnps.length);
				writer.print(cnps[i].getUCSClocation()+"\t"+numberOfMarkers[i]+"\t"+cnps[i].getChr()+"\t"+cnps[i].getStart()+"\t"+cnps[i].getStop());
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
	        System.err.println("Error writing to "+ext.rootOf(mlrrSetFile)+".xln");
	        e.printStackTrace();
        }
	}

	public static void dump(Project proj, String phenotype, String mlrrSetFile, String cnvToDump) {
        PrintWriter writer;
		MeanLRRset mlrrSet;
		float[][][] data;
		SampleList sampleList;
		String[] samples;
		Segment[] cnps;
		Hashtable<String, String> hash;
		double[] pheno;
		String phen;
		int index;
		
		hash = HashVec.loadFileToHashString(proj.getFilename(Project.SAMPLE_DATA_FILENAME), "DNA", new String[] {phenotype}, "");
		sampleList = proj.getSampleList(); 
		samples = sampleList.getSamples();
		pheno = new double[samples.length];
		for (int i = 0; i<samples.length; i++) {
			if (hash.containsKey(samples[i])) {
				phen = hash.get(samples[i]);
				if (phen.equals(".") || phen.equals("NaN")) {
					pheno[i] = Double.NaN;
				} else {
					pheno[i] = Double.parseDouble(phen);
				}
			} else {
				System.err.println("Error - no phenotypic data for sample: "+samples[i]);
				pheno[i] = Double.NaN;
			}
        }	
		
		mlrrSet = MeanLRRset.load(mlrrSetFile, false);
		data = mlrrSet.getData();
		cnps = mlrrSet.getCnps();
		if (mlrrSet.getSampleFingerprint() != sampleList.getFingerprint()) {
			System.err.println("Error - the SampleList fingerprint for the MeanLRRset ("+mlrrSet.getSampleFingerprint()+") does not match the Project's SampleList fingerprint ("+sampleList.getFingerprint()+")");
			System.exit(1);
		}
		
		index = -1;
		for (int i = 0; i<cnps.length; i++) {
			if (cnvToDump.equals(cnps[i].getUCSClocation())) {
				index = i;
			}
        }
		if (index == -1) {
			System.err.println("Error - CNV flagged to be dumped ("+cnvToDump+") was not found in the MeanLRRset's CNP list");
			System.exit(1);
		}
		
		try {
			writer = new PrintWriter(new FileWriter(proj.getProjectDir()+ext.replaceAllWith(cnvToDump, ":", "_") +".xln"));
			writer.println("Sample\tphenotype\t"+Array.toStr(Transforms.TRANFORMATIONS));
			for (int i = 0; i<samples.length; i++) {
				writer.println(samples[i]+"\t"+pheno[i]+"\t"+Array.toStr(data[index][i]));
	        }
            writer.close();
        } catch (Exception e) {
	        System.err.println("Error writing to "+ext.rootOf(mlrrSetFile)+".xln");
	        e.printStackTrace();
        }		
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = Project.DEFAULT_PROJECT;
//		String cnps = "cnps.cnv";
//		String cnps = "C:\\Documents and Settings\\npankrat\\My Documents\\CNV_PD\\all_CNPs.txt";
//		String cnps = "chr5_CNP.txt";
//		String cnps = "USP32.txt";
//		String cnps = "DOCK5_CNPs.txt";
		String cnps = "Colins_SLC2A3.txt";
		
		String phenotype = "CLASS=Use32Phenotype";
//		String phenotype = "CLASS=Gender";
//		String dump = "";
//		String dump = "chr6:31467630-31559451";
//		String dump = "chr9:38906782-39964796";
//		String dump = "chr5:151497267-151499003";
//		String dump = "chr8:25129632-25130278";
//		String dump = "chr8:25129632-25129964"; // first 3
//		String dump = "chr8:25130171-25130171"; // problem 4th
//		String dump = "chr8:25130231-25130278"; // last 2
//		String dump = "chr8:25129632-25129632"; // Marker 1
//		String dump = "chr8:25129709-25129709"; // Marker 2
//		String dump = "chr8:25129964-25129964"; // Marker 3
//		String dump = "chr8:25130171-25130171"; // Marker 4
//		String dump = "chr8:25130231-25130231"; // Marker 5
//		String dump = "chr8:25130278-25130278"; // Marker 6
		String dump = "chr12:7884583-8017012"; // Marker 6
		
		filename = "/home/npankrat/projects/GEDI_exome.properties";
		phenotype = "Class=SexReversal";
		cnps = "PAR2.txt";
		dump = "chr23:154879620-155227607"; // Marker 6

		
		Project proj;

		String usage = "\n"+
				"cnv.analysis.MeanLRR requires 0-1 arguments\n"+"" +
				"   (1) project file (i.e. proj="+filename+" (default))\n"+
				"   (2) filename of common CNVs (i.e. cnps="+cnps+" (default))\n"+
				"   (3) phenotype in SampleData.txt (i.e. pheno="+phenotype+" (default))\n"+
				"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("cnps=")) {
				cnps = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("pheno=")) {
				phenotype = args[i].split("=")[1];
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			proj = new Project(filename, false);
			if (!new File(proj.getProjectDir()+ext.rootOf(cnps)+".mlrr").exists()) {
//				createFilesFromFullSample(proj, cnps);
				createFilesFromMarkerData(proj, cnps);
			}
			if (!dump.equals("")) {
				dump(proj, phenotype, proj.getProjectDir()+ext.rootOf(cnps)+".mlrr", dump);
			} else {
				analyze(proj, phenotype, proj.getProjectDir()+ext.rootOf(cnps)+".mlrr");
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
