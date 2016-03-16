package one.ben;
import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.Toolkit;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.Transferable;
import java.awt.datatransfer.UnsupportedFlavorException;
import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.MalformedURLException;
import java.net.URISyntaxException;
import java.net.URL;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Random;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.WindowConstants;

import one.ben.ConditionalAnalysisPipeline.ConditionalAnalysisToolset_FAST;
import common.Array;
import common.Files;
import common.HashVec;
import common.Logger;
import common.Positions;
import common.ext;
import cnv.filesys.Centroids;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.manage.UCSCtrack;
import cnv.qc.CNVFilter;
import cnv.qc.CNVFilter.CNVFilterPass;
import cnv.var.CNVariant;
import cnv.var.SampleData;
import filesys.DosageData;
import filesys.MergeExtractPipeline;
import filesys.Segment;
import filesys.SnpMarkerSet;

public class lab {
	
    
    private static void createRandomSelectionFile() {
        String sourceFile = "F:/testing/all.txt";
        int[] sizes = {1000, 10000, 50000, 100000, 500000, 1000000};
        String[] data = HashVec.loadFileToStringArray(sourceFile, false, new int[]{0}, false);
        
        for (int size : sizes) {
            System.out.println("Printing " + size);
            String outFile = ext.rootOf(sourceFile, false) + "test_" + size + ".txt";
            PrintWriter writer = Files.getAppropriateWriter(outFile);
            Random rand = new Random();
            for (int i = 0; i < size; i++) {
                writer.println(data[rand.nextInt(data.length)]);
            }
            writer.flush();
            writer.close();
        }
        
    }
    
    private static void trimAlleles() {
        
        String dir = "/scratch.global/cole0482/CARDIA_CHS/";
        String dir2 = "/scratch.global/cole0482/CARDIA_CHS/oneAllele/";
        String[] fileRoots = new String[]{
            "chr1_P1_V3",  
            "chr2_P1_V3",
            "chr3_P1_V3",
            "chr4_P1_V3",
            "chr5_P1_V3",
            "chr6_P1_V3",
            "chr7_P1_V3",
            "chr8_P1_V3",
            "chr9_P1_V3",
            "chr10_P1_V3",
            "chr11_P1_V3",
            "chr12_P1_V3",
            "chr13_P1_V3",
            "chr14_P1_V3",
            "chr15_P1_V3",
            "chr16_P1_V3",
            "chr17_P1_V3",
            "chr18_P1_V3",
            "chr19_P1_V3",  
            "chr20_P1_V3",
            "chr21_P1_V3",
            "chr22_P1_V3",
            "chr23_FEMALE_no.auto.P1_V3",  
            "chr23_MALE_no.auto.P1_V3"    
        };
        
        for (String file : fileRoots) {
            System.out.println("Processing file : " + file);
//            String dataFile = dir + file + ".dose";
            String infoFile = dir + file + ".info";
            String newInfoFile = dir2 + file + ".info";
            
            try {
                BufferedReader reader = Files.getAppropriateReader(infoFile);
                PrintWriter writer = Files.getAppropriateWriter(newInfoFile);
                
                writer.println(reader.readLine()); // header
                String line = null;
                while ((line = reader.readLine()) != null) {
                    
                    String[] parts = line.split("[\\s+]");
                    String[] snpParts = parts[0].split(":");
                    if (snpParts.length == 3) {
                        StringBuilder newLine = new StringBuilder();
                        newLine.append(snpParts[0]).append(":").append(snpParts[1]).append(":").append(snpParts[2].charAt(0));
                        for (int i = 1; i < parts.length; i++) {
                            newLine.append(" ").append(parts[i]);
                        }
                        writer.println(newLine.toString());
                    } else {
                        writer.println(line);
                    }
                    
                }
                writer.flush();
                writer.close();
                reader.close();
            } catch (IOException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
            
        }
        
    }
    
    
    private static void countCNVsForIndividuals(String indivFile, String cnvFile, String outFile) throws IOException {
		Hashtable<String, String> sampleKeyHash = new Hashtable<String, String>();
		BufferedReader reader = new BufferedReader(new FileReader(indivFile));
		String line = null;
		while ((line = reader.readLine()) != null) {
			String[] tmp = line.split("\t");
			sampleKeyHash.put(tmp[0] + "\t" + tmp[1], "");
		}
		reader.close();
		
		Vector<CNVariant> cnvs = CNVariant.loadPlinkFile(cnvFile, sampleKeyHash, true, false);
		
		PrintWriter writer = Files.getAppropriateWriter(outFile);
		for (CNVariant cnv : cnvs) {
			writer.println(cnv.toPlinkFormat());
		}
		writer.flush();
		writer.close();
	}
	
	private static void filterForMarkers(String dataFile, String markerFile) throws IOException {
		BufferedReader reader = Files.getAppropriateReader(markerFile);
		HashMap<String, String[]> mkrsBim = new HashMap<String, String[]>();
		
		String line = reader.readLine();
		do {
			String[] parts = line.split("[\\s]+");
//			String chr = parts[0];
			String mkr = parts[1];
			mkrsBim.put(mkr, parts);
		} while ((line = reader.readLine()) != null);
		
		reader.close();
		
		reader = Files.getAppropriateReader(dataFile);
		PrintWriter writer = new PrintWriter(ext.rootOf(dataFile, false) + ".trimmed");
		PrintWriter writerBim = new PrintWriter(ext.rootOf(dataFile, false) + ".trimmed.bim");
		
		reader.readLine();
		writer.println("MarkerName\tChr\tPosition\tAllele1\tAllele2\tFreq.Allele1.HapMapCEU\tb\tSE\tp\tN");
		while((line = reader.readLine()) != null) {
			String[] parts = line.split("[\\s]+");
			if (mkrsBim.containsKey(parts[0])) {
				String[] bimParts = mkrsBim.get(parts[0]);
				writer.print(parts[0]);
				writer.print("\t");
				writer.print(bimParts[0]);
				writer.print("\t");
				writer.print(bimParts[3]);
				writer.print("\t");
				writer.println(Array.toStr(Array.subArray(parts, 1)));
				writerBim.println(Array.toStr(bimParts, "\t"));
			}
		}
		writer.flush();
		writer.close();
		writerBim.flush();
		writerBim.close();
		reader.close();
		
	}
	
	private static void writeRegions() {
		String filePrefix = "D:/data/gedi_gwas/data/corrected/LRR_MEDIAN_Median_LRR_XY_Raw_";
		for (int i = 1; i < 5; i++) {
			String file = filePrefix + i + ".xln";
			String[] header = Files.getHeaderOfFile(file, null);
			for (String hdr : header) {
				if (hdr.contains("MEDIAN")) {
					System.out.println("chr" + hdr.split("chr")[1]);
				}
			}
		}
		
		
		
	}
	
	private static void stripMarkerNames() {
		String bimFileTemplate = "D:/LD qsub/chr#.bim";
		String outFileTemplate = "D:/LD qsub/chr#_noMarkers.bim";
		
		for (int i = 1; i < 23; i++) {
			BufferedReader reader;
			PrintWriter writer;
			try {
				reader = Files.getAppropriateReader(bimFileTemplate.replace("#", "" + i));
				writer = Files.getAppropriateWriter(outFileTemplate.replace("#", "" + i));
				
				String line = reader.readLine();
				do {
					String[] splitLine = line.split("\t"); 
					if (splitLine[1].startsWith("rs")) {
						String[] mkr = splitLine[1].split(":");
						mkr[0] = splitLine[0];
						splitLine[1] = Array.toStr(mkr, ":");
						writer.println(Array.toStr(splitLine, "\t"));
					} else {
						writer.println(line);
					}
					
				} while((line = reader.readLine()) != null);
				
				reader.close();
				writer.flush();
				writer.close();
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
		}
		
	}
	
	private static void compareMarkers() {
		String markerFile = "D:/LD qsub/markers.txt";
		String bimFileTemplate = "D:/LD qsub/chr#_eu_unrel.bim";
		String markerOut1 = "D:/LD qsub/markersFound.txt";
		String markerOut2 = "D:/LD qsub/markersMissing.txt";
		
		System.out.println("Reading markers...");
		HashSet<String> markers = new HashSet<String>();
		try {
			BufferedReader markerReader = new BufferedReader(new FileReader(markerFile));
			markerReader.readLine();
			String line = null;
			while((line = markerReader.readLine()) != null) {
				markers.add(line.trim());
			}
			markerReader.close();
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		System.out.println("..." + markers.size() + " marker names read");
		
		HashSet<String> found = new HashSet<String>();
		
		PrintWriter writer1;
		PrintWriter writer2;
		writer1 = Files.getAppropriateWriter(markerOut1);
		writer2 = Files.getAppropriateWriter(markerOut2);
		for (int i = 1; i < 23; i++) {
			BufferedReader reader;
			try {
				reader = Files.getAppropriateReader(bimFileTemplate.replace("#", "" + i));
				
				String line = reader.readLine();
				while((line = reader.readLine()) != null) {
					String[] splitLine = line.split("\t"); 
					if (splitLine[1].startsWith("rs")) {
						String marker = splitLine[1].split(":")[0];
						
						if (markers.contains(marker)) {
							if (!found.contains(marker)) {
								writer1.println(marker);
								found.add(marker);
								markers.remove(marker);
							}
						}
						
					}
				}
				
				reader.close();
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
		}
		
		for (String leftovers : markers) {
			writer2.println(leftovers);
		}
		
		writer1.flush();
		writer2.flush();
		writer1.close();
		writer2.close();
		
	}
	
	private static void formatLDResults() {
		String fileTemplate = "D:/LD qsub/chr#_eu_markers.ld";
		String markerFile = "D:/LD qsub/markerData.txt";
		String outfile1name = "D:/LD qsub/eu_snps.ld";
		String outfile2name = "D:/LD qsub/eu_markers.ld";
		
		System.out.println(ext.getTime()+"]\tReading markers...");
		HashSet<String> markers = new HashSet<String>();
		int cnt = 0;
		try {
			BufferedReader markerReader = new BufferedReader(new FileReader(markerFile));
			String line = markerReader.readLine();
			while((line = markerReader.readLine()) != null) {
				String[] temp = line.split("\t");
				markers.add(temp[0]);
				markers.add("chr" + temp[1] + ":" + temp[2]);
				cnt++;
			}
			markerReader.close();
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		System.out.println(ext.getTime()+"]\t..." + cnt + " marker names read");

		String header1 = "SNP_A\tSNP_B\tR2";
		String header2 = "SNP\t#Variants\tLD_SNPS";
		HashMap<String, ArrayList<String>> markerResults = new HashMap<String, ArrayList<String>>();
		
		try {
			PrintWriter writer;
			
			writer = Files.getAppropriateWriter(outfile1name);
			writer.println(header1);
			
			for (int i = 1; i < 23; i++) {
				BufferedReader reader;
				reader = Files.getAppropriateReader(fileTemplate.replace("#", i+""));
				reader.readLine();
				
				String line = null;
				while((line = reader.readLine()) != null) {
					String[] lineParts = line.split("\t");
					
					String snp1 = lineParts[0];
					String snp2 = lineParts[1];
					String r2str = lineParts[2];
					
					if (!markers.contains(snp1)) {
						writer.println(snp2 + "\t" + snp1 + "\t" + r2str);
					}
					
					if (markers.contains(snp1)) {
						ArrayList<String> snps = markerResults.get(snp1);
						if (null == snps) {
							snps = new ArrayList<String>();
							markerResults.put(snp1, snps);
						}
						snps.add(snp2);
					} else {
						ArrayList<String> snps = markerResults.get(snp2);
						if (null == snps) {
							snps = new ArrayList<String>();
							markerResults.put(snp2, snps);
						}
						snps.add(snp1);
					}
					
					
				}
				reader.close();
			}
			writer.flush();
			writer.close();
			
			writer = Files.getAppropriateWriter(outfile2name);
			writer.println(header2);
			for (java.util.Map.Entry<String, ArrayList<String>> entry : markerResults.entrySet()) {
				writer.print(entry.getKey());
				writer.print("\t");
				writer.print(entry.getValue().size());
				writer.print("\t");
				for (String str : entry.getValue()) {
					writer.print(str);
					writer.print(";");
				}
				writer.println();
				
				markers.remove(entry.getKey());
			}
			
			for (String marker : markers) {
				writer.println(marker + "\t0\t");
			}
			
			writer.flush();
			writer.close();
			
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		System.out.println(ext.getTime()+"]\t...Done!");
		
	}

	private static void filterLDFiles(double pct) {
		String markerFile = "D:/LD qsub/markerData.txt";
		String gzFilenameTemplate = "D:/LD qsub/chr#_eu.ld.gz";
		String outfileTemplate = "D:/LD qsub/chr#_eu_markers.ld";
		
		System.out.println(ext.getTime()+"]\tReading markers...");
		HashSet<String> markers = new HashSet<String>();
		int cnt = 0;
		try {
			BufferedReader markerReader = new BufferedReader(new FileReader(markerFile));
			String line = markerReader.readLine();
			while((line = markerReader.readLine()) != null) {
				String[] temp = line.split("\t");
				markers.add(temp[0]);
				markers.add("chr" + temp[1] + ":" + temp[2]);
				cnt++;
			}
			markerReader.close();
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		System.out.println(ext.getTime()+"]\t..." + (cnt + 1) + " marker names read");
		
		System.out.println(ext.getTime()+"]\tWriting LD Files...");
//		for (int i = 1; i < 23; i++) {
		for (int i = 6; i < 9; i++) {
//		for (int i = 17; i < 18; i++) {
			System.out.println(ext.getTime()+"]\t...chr " + i);
			BufferedReader gzReader;
			try {
				gzReader = Files.getAppropriateReader(gzFilenameTemplate.replace("#", i+""));
				PrintWriter writer = Files.getAppropriateWriter(outfileTemplate.replace("#", i+""));
				
				// skip header
				gzReader.readLine();
				writer.println("SNP_A\tSNP_B\tR2");
				cnt = 0;
				String temp = "";
				String[] line;
				while((temp = gzReader.readLine()) != null) {
//					if (cnt % 100 == 0) System.out.print(".");
					line = temp.split("[\\s]+");
					if (line.length < 8) {
						System.out.println(Array.toStr(line));
						continue;
					}
					if (ext.isValidDouble(line[7]) && Double.parseDouble(line[7]) >= pct) {
						String rs1 = line[3].split(":")[0];
						if (!rs1.startsWith("rs")) {
							rs1 = "chr" + rs1 + ":" + line[3].split(":")[1];
						}
						
						String rs2 = line[6].split(":")[0];
						if (!rs2.startsWith("rs")) {
							rs2 = "chr" + rs2 + ":" + line[6].split(":")[1];
						}

						if (markers.contains(rs1) || markers.contains(rs2)) {
							writer.println(rs1 + "\t" + rs2 + "\t" + line[7]);
							cnt++;
						}
					}
				}
				System.out.println();
				System.out.println(ext.getTime()+"]\t..." + cnt + " markers found");
				writer.flush();
				writer.close();
				gzReader.close();
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		System.out.println(ext.getTime()+"]\t... Done!");
	}
	
	private static void filterFAMFile() {
		String famFile = "D:/SIDS and IQ/penncnv.fam";

		String filterFile1 = "D:/SIDS and IQ/IQ/excludeList.txt";
		String filterFile2 = "D:/SIDS and IQ/ids_keep.dat";
		String outputFile1 = "D:/SIDS and IQ/IQ/use20Blood.fam";
		String outputFile2 = "D:/SIDS and IQ/IQ/unrelated.fam";
		
		HashSet<String> indivList1 = HashVec.convertHashNullToHashSet(HashVec.loadFileToHashString(filterFile1, new int[]{0, 1}, null, false, "\t", false, false, false));
		HashSet<String> indivList2 = HashVec.convertHashNullToHashSet(HashVec.loadFileToHashString(filterFile2, new int[]{0, 1}, null, false, "\t", false, false, false));
		
		try {
			BufferedReader reader = new BufferedReader(new FileReader(famFile));
			PrintWriter writer1 = new PrintWriter(new FileWriter(outputFile1));
			PrintWriter writer2 = new PrintWriter(new FileWriter(outputFile2));
			
			String temp = "";
			String[] line;
			while(reader.ready()) {
				temp = reader.readLine();
				line = temp.split("\t");
				String indiv = line[0] + "\t" + line[1];
				
				if (indivList1.contains(indiv)) {
					writer1.println(temp);
				}
				if (indivList2.contains(indiv)) {
					writer2.println(temp);
				}
			}
			
			reader.close();
			writer1.flush();
			writer1.close();
			writer2.flush();
			writer2.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private static void mockupGUI() {
		JFrame frame = new JFrame();
		frame.setLayout(new FlowLayout());
		JCheckBox chkBx = new JCheckBox("Derive from Centroids:");
		chkBx.setFont(new Font("Arial", 0, 14));
		JComboBox<String> comboBx = new JComboBox<String>(new String[]{"sexSpecific_Male", "sexSpecific_Female"});
		comboBx.setFont(new Font("Arial", 0, 14));
		frame.add(chkBx);
		frame.add(comboBx);
		frame.add(new JButton());
		frame.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
		frame.setVisible(true);
	}
	
	
	private static void fakeAutoCentFile(Project proj, String centFile) {
		byte[] markerChrs;
		String[] allMarkers;
		MarkerSet ms;
		PrintWriter writer;
		Vector<String> markerList;

		ms = proj.getMarkerSet();
		allMarkers = ms.getMarkerNames();
		markerChrs = ms.getChrs();
		markerList = new Vector<String>();
		
		Centroids centroid = Centroids.load(centFile, false);
		float[][][] newCentroids = new float[allMarkers.length][][];
		
		int breakInd = 0;
		for (int i = 0; i < markerChrs.length; i++) {
			switch(markerChrs[i]) {
				case 23:
				case 24:
				case 25:
				case 26:
					if (breakInd == 0) breakInd = i;
					break;
				default:
					markerList.add(allMarkers[i]);
					break;
			}
		}
		
		for (int i = 0; i < breakInd; i++) {
			newCentroids[i] = new float[][]{{Float.NaN, Float.NaN}, {Float.NaN, Float.NaN}, {Float.NaN, Float.NaN}};
		}
		for (int i = breakInd; i < newCentroids.length; i++) {
			newCentroids[i] = centroid.getCentroids()[i - breakInd];
		}
		
		Centroids newCentObj = new Centroids(newCentroids, ms.getFingerprint());
		newCentObj.serialize(centFile + ".faked");
		String dir = proj.PROJECT_DIRECTORY.getValue();
		String relativeFile = centFile.substring(dir.length());
		Centroids.exportToText(proj, relativeFile + ".faked", relativeFile + ".faked.txt");
		
	}
	
	private static void idSwap(Project proj, String fileIn) throws IOException {
		BufferedReader reader = new BufferedReader(new FileReader(fileIn));
		String outFile = ext.rootOf(fileIn, false) + ".ids";
		PrintWriter writer = new PrintWriter(new FileWriter(outFile));
		
		SampleData sampleData = proj.getSampleData(0, false);
		
		while(reader.ready()) {
			String line = reader.readLine();
			String[] keyVal = sampleData.lookup(line);
			writer.println(Array.toStr(keyVal, "\t"));
		}
		writer.flush();
		writer.close();
		reader.close();
	}
	
	private static void writeExclusions() {
		String projectFile = "D:/projects/gedi_gwas.properties";
		String outFile = "D:/data/gedi_gwas/data/indivList.txt";
		Project proj = new Project(projectFile, false);
		String[] samps = proj.getSamples();
		boolean[] excl = proj.getSamplesToExclude();
		PrintWriter writer = Files.getAppropriateWriter(outFile);
		for (int i = 0; i < samps.length; i++) {
			writer.println(samps[i] + "\t" + (excl[i] ? "1" : "0"));
		}
		writer.flush();
		writer.close();
	}
	
	
	private static void filter() throws IOException {
		CNVariant[] cnvList1 = CNVariant.loadPlinkFile("D:/SIDS and IQ/IQ/filtertest1.cnv", false);
		CNVariant[] cnvList2 = CNVariant.loadPlinkFile("D:/SIDS and IQ/IQ/filtertest2.cnv", false);
		
		
		HashSet<CNVariant> list1 = new HashSet<CNVariant>();
		for (CNVariant cnv : cnvList1) {
			list1.add(cnv);
		}
		HashSet<CNVariant> list2 = new HashSet<CNVariant>();
		for (CNVariant cnv : cnvList2) {
			list2.add(cnv);
		}
		
		ArrayList<CNVariant> retain1 = new ArrayList<CNVariant>();
		for (CNVariant cnv : list1) {
			if (list2.contains(cnv)) {
				continue;
			}
			retain1.add(cnv);
		}
		
		ArrayList<CNVariant> retain2 = new ArrayList<CNVariant>();
		for (CNVariant cnv : list2) {
			if (list1.contains(cnv)) {
				continue;
			}
			retain2.add(cnv);
		}
		
		PrintWriter writer1 = new PrintWriter(new FileWriter("D:/SIDS and IQ/IQ/filtercheck1.cnv"));
		PrintWriter writer2 = new PrintWriter(new FileWriter("D:/SIDS and IQ/IQ/filtercheck2.cnv"));
		writer1.println(Array.toStr(CNVariant.PLINK_CNV_HEADER, "\t"));
		writer2.println(Array.toStr(CNVariant.PLINK_CNV_HEADER, "\t"));
		
		for (CNVariant cnv : retain1) {
			writer2.println(cnv.toPlinkFormat());
		}
		for (CNVariant cnv : retain2) {
			writer2.println(cnv.toPlinkFormat());
		}

		writer1.flush();
		writer1.close();
		writer2.flush();
		writer2.close();
		
	}
	

	public static void filterCentromeric(String dir, String in, String out, String markerSetFilenameToBreakUpCentromeres, int build, Logger log) {
		BufferedReader reader;
		PrintWriter writer;
		String[] line;
		CNVariant cnv;
		Segment[] centromereMidpoints;
		int[][] centromereBoundaries;
		
		centromereBoundaries = Positions.determineCentromereBoundariesFromMarkerSet(markerSetFilenameToBreakUpCentromeres, build, log);
		centromereMidpoints = Positions.computeCentromereMidpoints(centromereBoundaries);
		
		try {
			reader = new BufferedReader(new FileReader(dir+in));
			writer = new PrintWriter(new FileWriter(dir+out));
			writer.println(reader.readLine());
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				cnv = new CNVariant(line);
				if (cnv.overlaps(centromereMidpoints[cnv.getChr()])) {
					writer.println(Array.toStr(line));
				}
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+dir+in+"\" not found in current directory");
			return;
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+dir+in+"\"");
			return;
		}
	}
	
	public static void breakCentromeric() throws IOException {
		CNVFilter filter = new CNVFilter(null);
		filter.setBreakupCentromeres(true);
		filter.setCentromereBoundariesFromFile("D:/data/gedi_gwas/data/markers.bim");
		filter.computeCentromereMidPoints();
		
		CNVariant[] centromeric = CNVariant.loadPlinkFile("D:/SIDS and IQ/IQ/merged.cnv", false);
		
		PrintWriter writer = new PrintWriter(new FileWriter("D:/SIDS and IQ/IQ/merged_split.cnv"));
		writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER, "\t"));
		
		for (CNVariant cnv : centromeric) {
			CNVFilterPass fp = filter.getCNVFilterPass(cnv);
			if (fp.isCentromeric()) {
				CNVariant[] broken = filter.breakUpCentromere(fp, cnv);
				for (CNVariant newcnv : broken) {
					writer.println(newcnv.toPlinkFormat());
				}
			} else {
				writer.println(cnv.toPlinkFormat());
			}
		}
		writer.flush();
		writer.close();
		
	}
	
	
	private static void filterWrong() throws IOException {
		CNVariant[] cnvList = CNVariant.loadPlinkFile("D:/SIDS and IQ/IQ/c10_p3,15.cnv", false);

		PrintWriter writer = new PrintWriter(new FileWriter("D:/CNVcheck.cnv"));
		writer.println(Array.toStr(CNVariant.PLINK_CNV_HEADER, "\t"));
		
		for (CNVariant cnv : cnvList) {
			if (cnv.getStart() > cnv.getStop()) {
				writer.println(cnv.toPlinkFormat());
			}
		}
		writer.flush();
		writer.close();
		
	}
	
	private static void processSNPPositions(String file, String out) throws IOException {
	    BufferedReader reader = Files.getAppropriateReader(file);
        PrintWriter writer = Files.getAppropriateWriter(out);
        
        String line = reader.readLine();
        writer.println(line);
        while((line = reader.readLine()) != null) {
            String[] parts = line.split("[\\s]+");
            
            String rs = parts[0];
            
            if (rs.indexOf(":") != -1) {
                String[] sub = rs.split(":");
                writer.println(rs + "\t" + sub[0] + "\t" + sub[1]);
            } else {
                writer.println(line);
            }
            
        }
        
        writer.flush();
        writer.close();
	}
	
	private static void orderAlleles(String file, String out) throws IOException {
	    BufferedReader reader = Files.getAppropriateReader(file);
	    PrintWriter writer = Files.getAppropriateWriter(out);
	    
	    String line = null;
	    while((line = reader.readLine()) != null) {
	        String[] parts = line.split("[\\s]+");
	        
	        String rs = parts[0];
	        String a1 = parts[1];
	        String a2 = parts[2];
	        String ma = parts[3];
	        
	        String oa = ma.equalsIgnoreCase(a1) ? a2 : a1;
	        writer.println(rs + " " + ma + " " + oa);
	    }
	    
	    writer.flush();
	    writer.close();
	}
	
	public enum TEST {
	    TEST1,
	    TEST2,
	    TEST3;
	}
	
	public static <T extends Enum<T>> void enumValues(Class<T> enumType) {
        for (T c : enumType.getEnumConstants()) {
             System.out.println(c.name());
        }
	}
	
	public static void main(String[] args) {
		int numArgs = args.length;
		Project proj;
		String filename = "lab.dat";
		String logfile = null;
		String file = null;
		
		boolean test = true;
		if (test) {

		    String markersFile = "/scratch.global/cole0482/merge/mkrs.txt";
		    String regionsFile = "/scratch.global/cole0482/merge/rgns.txt";
		    String outFile = "/scratch.global/cole0482/test.db.xln.gz";
		    String mapOutFile = "/scratch.global/cole0482/mapOut.xln";

	        
		    MergeExtractPipeline pipeline = new MergeExtractPipeline();
//		    pipeline.setMarkers(markersFile);		    
		    pipeline.setRunDirectory("/scratch.global/cole0482/merge/", true);
		    pipeline.setOutputFormat(DosageData.DATABASE_DOSE_FORMAT);
		    pipeline.setOutputFiles(outFile, mapOutFile);
		    pipeline.setRenamePlinkMarkers(true);
//		    pipeline.addDataSource("/scratch.global/cole0482/merge/blacks/", "gwas.bed", "gwas.bim", "gwas.fam");
		    pipeline.addDataSource("/scratch.global/cole0482/merge/blacks/", "exome.bed", "exome.bim", "exome.fam");
		    pipeline.addDataSource("/scratch.global/cole0482/merge/blacks/", "metab.bed", "metab.bim", "metab.fam");
		    // add more;
		    pipeline.run();
		    
		    
//		    String doseFile1 = "/home/pankarne/shared/ARIC_Genomics_Data/GWAS_Chip/1000G/ARIC.whites.impute2/chr3.90069244.95069244.impute2.gz";
//		    String mapFile1 = "/home/pankarne/shared/ARIC_Genomics_Data/GWAS_Chip/1000G/ARIC.whites.impute2/chr3.90069244.95069244.impute2_info";
//		    
//		    String doseFile2 = "/home/pankarne/shared/ARIC_Genomics_Data/GWAS_Chip/1000G/ARIC.whites.impute2/chr3.95069244.100069244.impute2.gz";
//		    String mapFile2 = "/home/pankarne/shared/ARIC_Genomics_Data/GWAS_Chip/1000G/ARIC.whites.impute2/chr3.95069244.100069244.impute2_info";
//		    
//		    String idFile = "/home/pankarne/cole0482/EA.indiv.dup";
//		    String outFile = "/scratch.global/cole0482/test.db.xln.gz";
//		    String mapOutFile = "/scratch.global/cole0482/mapOut.xln";
//		    
//		    DosageData dd1 = new DosageData(doseFile1, idFile, mapFile1, DosageData.IMPUTE2_DOSE_FORMAT, null, true, null);
//		    DosageData dd2 = new DosageData(doseFile2, idFile, mapFile2, DosageData.IMPUTE2_DOSE_FORMAT, null, true, null);
//		    DosageData dd3 = DosageData.combine(dd1, dd2);
//		    dd1 = null;
//		    dd2 = null;
//		    dd1 = DosageData.loadPlinkBinary(dir, plinkRoot);
//		    dd2 = DosageData.combine(dd3, dd1);
//		    dd1 = null;
//		    dd3 = null;
//		    dd1 = DosageData.loadPlinkBinary(dir2, plinkRoot2);
//		    dd3 = DosageData.combine(dd2, dd1);
//		    dd3.writeToFile(outFile, mapOutFile, null, DosageData.DATABASE_DOSE_FORMAT, null);
//		    System.out.println("complete!");
		    
//		    createRandomSelectionFile();
//		    try {
//    		    String checkfile = "D:/data/gedi_gwas/data/cluster.genome.gz";
//    		    BufferedReader reader = Files.getAppropriateReader(checkfile);
//    		    String line = reader.readLine();
//    		    int cnt = 1;
//    		    while ((line = reader.readLine()) != null) {
//    		        cnt++;
//    		    }
//    		    reader.close();
//    		    System.out.println("Read " + cnt + " lines");
		    
		    
		    
//		    String s = "10000000000";
//		    long l = Long.parseLong(s);
//		    System.out.println(l > Integer.MAX_VALUE);
//		    int d = Integer.parseInt(s);
		    
//		    
//                Process p = Runtime.getRuntime().exec("where notepad.exe");
//                int waitCode = p.waitFor();
//                int outCode = p.exitValue();
//                System.out.println("wait: " + waitCode + "| out: " + outCode);
//		        processData();
//		        concatData();
//		        ConditionalAnalysisPipeline.processOnly(args[0], args[1], args[2]);
//            } catch (IOException e) {
//                // TODO Auto-generated catch block
//                e.printStackTrace();
//            } /*catch (InterruptedException e) {
                // TODO Auto-generated catch block
//                e.printStackTrace();
//            }*/
//		    
//            String dir = "F:/CARDIA processing/";
//
//            String[] files = (new File(dir)).list(new FilenameFilter() {
//                @Override
//                public boolean accept(File dir, String name) {
//                    return name.endsWith("positions.xln");
//                }
//            });
//
//            for (String file2 : files) {
//                try {
//                    processSNPPositions(dir + file2, dir + ext.rootOf(file2) + ".proc.xln");
//                } catch (IOException e) {
                    // TODO Auto-generated catch block
//                    e.printStackTrace();
//                }
//            }

		    
		    return;
		}
		
		
		String usage = "";
		
		if (numArgs == 0) {
			try {
			countCNVsForIndividuals("D:/data/ny_registry/new_york/stats/puv_ids.txt", "D:/data/ny_registry/new_york/stats/recodedM.cnv", "D:/data/ny_registry/new_york/stats/puv_cnvs.cnv");
//			testClipboard();
//				BufferedReader reader = new BufferedReader(new FileReader("D:/ForestPlot/Hb_SingleSNP.csv"));
//				String line = reader.readLine();
//				do {
//					System.out.println(line);
//				} while((line = reader.readLine()) != null);
				
//			filterForMarkers("D:/height/scratch/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq_parsed.xln", "D:/height/scratch/samples_logan.bim");
//			writeRegions();
//			stripMarkerNames();
//				writeExclusions();
//			concatFiles();
//			splitFile();
//				idSwap(new Project("D:/projects/gedi_gwas.properties", false), "D:/data/gedi_gwas/overlap_ids.txt");
//			compareMarkers();
//				filterLDFiles(0.5);
//				formatLDResults();
//				filter();
//				breakCentromeric();
//				filterWrong();
			} catch (IOException e) {
//				 TODO Auto-generated catch block
				e.printStackTrace();
			}
//			mockupGUI();
			return;
		}
		
		
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("file=")) {
				file = args[i].split("=")[1];
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
			proj = new Project(filename, logfile, false);
			if (file != null) {
				idSwap(proj, file);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private static void concatData() throws IOException {
        String file = "C:/workspace/temp/temp2/region1_snpInfo";
	    String indivFile = "C:/workspace/temp/temp2/EA.indiv";
        String traitFile = "C:/workspace/temp/temp2/ARIC_EA_ProteinC.trait";
        
        String[] iids = HashVec.loadFileToStringArray(indivFile, false, new int[]{0}, false);
        HashMap<String, Integer> indexMap = new HashMap<String, Integer>();
        for (int i = 0; i < iids.length; i++) {
            indexMap.put(iids[i], i);
        }
        
        String[] genoData = HashVec.loadFileToStringArray(file + "_collapsed.xln", false, null, false);
        
        String[] traitData = HashVec.loadFileToStringArray(traitFile, false, null, false);
        
        PrintWriter writer = Files.getAppropriateWriter(file + "_withTrait.xln");
        
        writer.println(traitData[0] + "\t" + genoData[0]);
        for (int i = 1; i < traitData.length; i++) {
            String iid = traitData[i].split("[\\s]+")[1];
            Integer intgr = indexMap.get(iid); 
            if (intgr == null) continue;
            int index = intgr.intValue();
            writer.println(traitData[i] + "\t" + genoData[6 + index]);
            
            
        }
        writer.flush();
        writer.close();
        
	}
	
	private static void processData() throws IOException {
        String file = "C:/workspace/temp/temp2/region1_snpInfo";
        
        BufferedReader reader = Files.getAppropriateReader(file + ".txt");
        PrintWriter writer = Files.getAppropriateWriter(file + "_collapsed.xln");
        
        for (int i = 0; i < 6; i++) {
            String line = reader.readLine();
            writer.println(line);
        }
        String line = null;
        String line2, line3 = null;
        String[] parts1, parts2, parts3;
        while((line = reader.readLine()) != null) {
            line2 = reader.readLine();
            line3 = reader.readLine();
            
            parts1 = line.split("[\\s]+");
            parts2 = line2.split("[\\s]+");
            parts3 = line3.split("[\\s]+");
            
            String[] newVals = new String[parts1.length];
            for (int i = 0; i < parts1.length; i++) {
                double d2 = Double.parseDouble(parts2[i]);
                double d3 = Double.parseDouble(parts3[i]);
                newVals[i] = "" + (d2 + (2 * d3));
            }
            writer.println(Array.toStr(newVals));
            
        }
        writer.flush();
        writer.close();
        
    }

    private static void testClipboard() {
		Clipboard systemClipboard;
		Transferable contents;
		systemClipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
		contents = systemClipboard.getContents(null);
		if (contents==null) {
			System.out.println("Clipboard is empty");
		} else {
			try {
//				DataFlavor[] systemFlavors = systemClipboard.getAvailableDataFlavors();
//				DataFlavor[] transFlavors = contents.getTransferDataFlavors();
				DataFlavor RTF_FLAVOR = new DataFlavor("text/rtf", "Rich Formatted Text");

				if (contents.isDataFlavorSupported(RTF_FLAVOR)) {
					ByteArrayInputStream obj = (ByteArrayInputStream) contents.getTransferData(RTF_FLAVOR);
					int n = obj.available();
					byte[] bytes = new byte[n];
					obj.read(bytes, 0, n);
					String s = new String(bytes, StandardCharsets.UTF_8);
					System.out.println(s);
				}
				System.out.println();
//				if (contents.isDataFlavorSupported(DataFlavor.stringFlavor)) {
//					return (String)contents.getTransferData(DataFlavor.stringFlavor);
//				}
			} catch (UnsupportedFlavorException ufe) {
				ufe.printStackTrace();
			} catch (IOException ioe) {
				ioe.printStackTrace();
			}
		}
	}
	
}



